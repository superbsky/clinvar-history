"""Reusable offline analyses over compact ClinVar event histories."""

import calendar
from collections import Counter
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from datetime import date
from enum import StrEnum
from pathlib import Path

import pyarrow.parquet as pq

from .attribution import AggregateAttribution, AttributionMechanism
from .classification import ClassificationState
from .diff import EventType, VcvEvent
from .scv_diff import ScvEvent, ScvEventType

_BENIGN_STATES = {
    ClassificationState.BENIGN,
    ClassificationState.LIKELY_BENIGN,
}
_PATHOGENIC_STATES = {
    ClassificationState.PATHOGENIC,
    ClassificationState.LIKELY_PATHOGENIC,
}
_RESOLVED_STATES = _BENIGN_STATES | _PATHOGENIC_STATES


@dataclass(frozen=True, slots=True)
class SamplingBias:
    true_transitions: int
    observed_transitions: int
    event_recall: float | None
    missed_recurrent_transitions: int
    missed_conflict_episodes: int
    mean_event_date_error_days: float | None


@dataclass(frozen=True, slots=True)
class CohortVusMetrics:
    period_start: date
    period_end: date
    vus_at_start: int
    vus_resolved: int
    new_vus: int
    reclassification_rate: float | None
    discovery_rate: float | None


@dataclass(frozen=True, slots=True)
class RenovoQuarterlyMetrics:
    period_start: date
    period_end: date
    vus_previous: int
    vus_current: int
    reclassified_vus: int
    net_vus_change: int
    reclassification_rate: float | None
    discovery_rate: float | None
    discovery_reclassification_ratio: float | None


@dataclass(frozen=True, slots=True)
class TransitionCount:
    old_classification: ClassificationState
    new_classification: ClassificationState
    count: int


@dataclass(frozen=True, slots=True)
class StateEpisode:
    variation_id: int
    classification: ClassificationState
    started_at: date
    ended_at: date
    duration_days: int
    censored: bool


class CompetingOutcome(StrEnum):
    BENIGN = "benign"
    PATHOGENIC = "pathogenic"
    CONFLICT = "conflict"
    OTHER_EXIT = "other_exit"
    CENSORED = "censored"


@dataclass(frozen=True, slots=True)
class VusResolution:
    variation_id: int
    started_at: date
    ended_at: date
    duration_days: int
    outcome: CompetingOutcome
    review_status: str | None = None
    gene_info: str | None = None
    conditions: str | None = None
    variant_type: str | None = None


@dataclass(frozen=True, slots=True)
class ScvVusResolution:
    scv_accession: str
    variation_id: int
    started_at: date
    ended_at: date
    duration_days: int
    outcome: CompetingOutcome
    submitter: str | None = None
    collection_method: str | None = None
    review_status: str | None = None
    conditions: str | None = None


def _events_by_variation(
    events: Iterable[VcvEvent],
) -> dict[int, list[VcvEvent]]:
    grouped: dict[int, list[VcvEvent]] = {}
    for event in events:
        grouped.setdefault(event.variation_id, []).append(event)
    for history in grouped.values():
        history.sort(key=lambda event: (event.observed_at, event.event_type.value))
    return grouped


def _add_months(value: date, months: int) -> date:
    month_index = value.year * 12 + value.month - 1 + months
    year = month_index // 12
    month = month_index % 12 + 1
    day = min(value.day, calendar.monthrange(year, month)[1])
    return date(year, month, day)


def month_grid(start: date, end: date, months: int) -> list[date]:
    """Build a calendar grid without silently moving releases to day one."""
    if months < 1:
        raise ValueError("months must be positive")
    dates = [start]
    current = start
    while current < end:
        current = min(_add_months(current, months), end)
        dates.append(current)
    return dates


def _state_at(
    history: Sequence[VcvEvent],
    observed_at: date,
) -> ClassificationState | None:
    state = None
    for event in history:
        if event.observed_at > observed_at:
            break
        if event.event_type is EventType.VARIANT_REMOVED:
            state = None
        elif event.event_type in {
            EventType.VARIANT_ADDED,
            EventType.CLASSIFICATION_CHANGED,
        }:
            state = event.new_classification
    return state


def _snapshot(
    grouped: dict[int, list[VcvEvent]],
    observed_at: date,
) -> dict[int, ClassificationState]:
    states = {}
    for variation_id, history in grouped.items():
        state = _state_at(history, observed_at)
        if state is not None:
            states[variation_id] = state
    return states


def sampling_bias(
    events: Iterable[VcvEvent],
    sample_dates: Sequence[date],
) -> SamplingBias:
    """Compare all observed classification events with down-sampled snapshots.

    Sweep the chronological event log once instead of rebuilding a complete
    variant snapshot for every sampling date.
    """
    if len(sample_dates) < 2:
        raise ValueError("At least two sample dates are required")
    dates = sorted(set(sample_dates))
    if len(dates) != len(sample_dates):
        raise ValueError("Sample dates must be unique")

    event_list = list(events)
    event_list.sort(
        key=lambda event: (
            event.observed_at,
            event.variation_id,
            event.event_type.value,
        )
    )
    first, last = dates[0], dates[-1]
    true_transitions = sum(
        1
        for event in event_list
        if event.event_type is EventType.CLASSIFICATION_CHANGED
        and event.old_classification != event.new_classification
        and first < event.observed_at <= last
    )
    true_conflict_episodes = sum(
        1
        for event in event_list
        if event.event_type is EventType.CLASSIFICATION_CHANGED
        and event.old_classification != event.new_classification
        and event.new_classification is ClassificationState.CONFLICT
        and first < event.observed_at <= last
    )

    observed_transitions = 0
    observed_conflict_episodes = 0
    date_errors = []
    states: dict[int, ClassificationState] = {}
    event_index = 0

    def apply_state_event(event: VcvEvent) -> None:
        if event.event_type is EventType.VARIANT_REMOVED:
            states.pop(event.variation_id, None)
        elif event.event_type in {
            EventType.VARIANT_ADDED,
            EventType.CLASSIFICATION_CHANGED,
        }:
            if event.new_classification is None:
                states.pop(event.variation_id, None)
            else:
                states[event.variation_id] = event.new_classification

    while (
        event_index < len(event_list) and event_list[event_index].observed_at <= first
    ):
        apply_state_event(event_list[event_index])
        event_index += 1

    for observed_at in dates[1:]:
        interval_start_states: dict[int, ClassificationState | None] = {}
        latest_transition_dates: dict[tuple[int, ClassificationState | None], date] = {}
        while (
            event_index < len(event_list)
            and event_list[event_index].observed_at <= observed_at
        ):
            event = event_list[event_index]
            if event.event_type in {
                EventType.VARIANT_ADDED,
                EventType.VARIANT_REMOVED,
                EventType.CLASSIFICATION_CHANGED,
            }:
                interval_start_states.setdefault(
                    event.variation_id,
                    states.get(event.variation_id),
                )
                if (
                    event.event_type is EventType.CLASSIFICATION_CHANGED
                    and event.old_classification != event.new_classification
                ):
                    latest_transition_dates[
                        (event.variation_id, event.new_classification)
                    ] = event.observed_at
                apply_state_event(event)
            event_index += 1

        changed_ids = {
            variation_id
            for variation_id, previous_state in interval_start_states.items()
            if previous_state is not None
            and (current_state := states.get(variation_id)) is not None
            and previous_state != current_state
        }
        observed_transitions += len(changed_ids)
        observed_conflict_episodes += sum(
            states[variation_id] is ClassificationState.CONFLICT
            for variation_id in changed_ids
        )
        for variation_id in changed_ids:
            actual = latest_transition_dates.get((variation_id, states[variation_id]))
            if actual is not None:
                date_errors.append((observed_at - actual).days)

    recall = observed_transitions / true_transitions if true_transitions else None
    return SamplingBias(
        true_transitions=true_transitions,
        observed_transitions=observed_transitions,
        event_recall=recall,
        missed_recurrent_transitions=max(
            true_transitions - observed_transitions,
            0,
        ),
        missed_conflict_episodes=max(
            true_conflict_episodes - observed_conflict_episodes,
            0,
        ),
        mean_event_date_error_days=(
            sum(date_errors) / len(date_errors) if date_errors else None
        ),
    )


def cohort_vus_metrics(
    events: Iterable[VcvEvent],
    sample_dates: Sequence[date],
) -> list[CohortVusMetrics]:
    """Measure resolution of a starting VUS cohort and newly entering VUS."""
    if len(sample_dates) < 2:
        raise ValueError("At least two sample dates are required")
    dates = sorted(set(sample_dates))
    grouped = _events_by_variation(events)
    snapshots = [_snapshot(grouped, observed_at) for observed_at in dates]
    results = []

    for period_start, period_end, start, end in zip(
        dates[:-1],
        dates[1:],
        snapshots[:-1],
        snapshots[1:],
        strict=True,
    ):
        starting_vus = {
            variation_id
            for variation_id, state in start.items()
            if state is ClassificationState.VUS
        }
        resolved = sum(
            end.get(variation_id) in _RESOLVED_STATES for variation_id in starting_vus
        )
        new_vus = sum(
            state is ClassificationState.VUS and variation_id not in start
            for variation_id, state in end.items()
        )
        denominator = len(starting_vus)
        results.append(
            CohortVusMetrics(
                period_start=period_start,
                period_end=period_end,
                vus_at_start=denominator,
                vus_resolved=resolved,
                new_vus=new_vus,
                reclassification_rate=(resolved / denominator if denominator else None),
                discovery_rate=(new_vus / denominator if denominator else None),
            )
        )
    return results


def renovo_quarterly_metrics(
    events: Iterable[VcvEvent],
    sample_dates: Sequence[date],
) -> list[RenovoQuarterlyMetrics]:
    """Reproduce Bonetti et al.'s consecutive-snapshot rate formulas."""
    if len(sample_dates) < 2:
        raise ValueError("At least two sample dates are required")
    dates = sorted(set(sample_dates))
    grouped = _events_by_variation(events)
    snapshots = [_snapshot(grouped, observed_at) for observed_at in dates]
    results = []

    for period_start, period_end, previous, current in zip(
        dates[:-1],
        dates[1:],
        snapshots[:-1],
        snapshots[1:],
        strict=True,
    ):
        previous_vus = {
            variation_id
            for variation_id, state in previous.items()
            if state is ClassificationState.VUS
        }
        current_vus = {
            variation_id
            for variation_id, state in current.items()
            if state is ClassificationState.VUS
        }
        reclassified = sum(
            current.get(variation_id) in _RESOLVED_STATES
            for variation_id in previous_vus
        )
        denominator = len(current_vus)
        net_change = denominator - len(previous_vus)
        reclassification_rate = reclassified / denominator if denominator else None
        discovery_rate = net_change / denominator if denominator else None
        ratio = (
            discovery_rate / reclassification_rate
            if reclassification_rate not in {None, 0} and discovery_rate is not None
            else None
        )
        results.append(
            RenovoQuarterlyMetrics(
                period_start=period_start,
                period_end=period_end,
                vus_previous=len(previous_vus),
                vus_current=denominator,
                reclassified_vus=reclassified,
                net_vus_change=net_change,
                reclassification_rate=reclassification_rate,
                discovery_rate=discovery_rate,
                discovery_reclassification_ratio=ratio,
            )
        )
    return results


def transition_counts(events: Iterable[VcvEvent]) -> list[TransitionCount]:
    """Count every aggregate classification transition, including recurrence."""
    counts = Counter(
        (event.old_classification, event.new_classification)
        for event in events
        if event.event_type is EventType.CLASSIFICATION_CHANGED
        and event.old_classification is not None
        and event.new_classification is not None
        and event.old_classification != event.new_classification
    )
    return [
        TransitionCount(old, new, count)
        for (old, new), count in sorted(
            counts.items(),
            key=lambda item: (item[0][0].value, item[0][1].value),
        )
    ]


def state_episodes(
    events: Iterable[VcvEvent],
    *,
    study_end: date,
) -> list[StateEpisode]:
    """Construct half-open state intervals and right-censor active states."""
    episodes = []
    for variation_id, history in _events_by_variation(events).items():
        state: ClassificationState | None = None
        started_at: date | None = None
        for event in history:
            if event.observed_at > study_end:
                break
            if event.event_type not in {
                EventType.VARIANT_ADDED,
                EventType.VARIANT_REMOVED,
                EventType.CLASSIFICATION_CHANGED,
            }:
                continue
            if (
                event.event_type is EventType.CLASSIFICATION_CHANGED
                and event.old_classification == event.new_classification
            ):
                continue
            if state is not None and started_at is not None:
                episodes.append(
                    StateEpisode(
                        variation_id=variation_id,
                        classification=state,
                        started_at=started_at,
                        ended_at=event.observed_at,
                        duration_days=(event.observed_at - started_at).days,
                        censored=False,
                    )
                )
            state = (
                None
                if event.event_type is EventType.VARIANT_REMOVED
                else event.new_classification
            )
            started_at = event.observed_at if state is not None else None

        if state is not None and started_at is not None and started_at <= study_end:
            episodes.append(
                StateEpisode(
                    variation_id=variation_id,
                    classification=state,
                    started_at=started_at,
                    ended_at=study_end,
                    duration_days=(study_end - started_at).days,
                    censored=True,
                )
            )
    return sorted(
        episodes,
        key=lambda episode: (episode.variation_id, episode.started_at),
    )


def vus_resolutions(
    events: Iterable[VcvEvent],
    *,
    study_end: date,
    include_left_censored: bool = False,
) -> list[VusResolution]:
    """Measure each VUS episode to benign/pathogenic resolution or censoring."""
    event_list = list(events)
    first_snapshot = min(
        (event.observed_at for event in event_list),
        default=None,
    )
    results = []
    for variation_id, history in _events_by_variation(event_list).items():
        for index, event in enumerate(history):
            if event.observed_at > study_end:
                break
            if event.event_type not in {
                EventType.VARIANT_ADDED,
                EventType.CLASSIFICATION_CHANGED,
            }:
                continue
            if event.new_classification is not ClassificationState.VUS:
                continue
            if event.old_classification is ClassificationState.VUS:
                continue
            if (
                not include_left_censored
                and event.event_type is EventType.VARIANT_ADDED
                and event.observed_at == first_snapshot
            ):
                continue

            outcome = CompetingOutcome.CENSORED
            ended_at = study_end
            for later in history[index + 1 :]:
                if later.observed_at > study_end:
                    break
                if later.event_type not in {
                    EventType.VARIANT_REMOVED,
                    EventType.CLASSIFICATION_CHANGED,
                }:
                    continue
                if later.event_type is EventType.VARIANT_REMOVED:
                    ended_at = later.observed_at
                    break
                if later.new_classification in _BENIGN_STATES:
                    outcome = CompetingOutcome.BENIGN
                    ended_at = later.observed_at
                    break
                if later.new_classification in _PATHOGENIC_STATES:
                    outcome = CompetingOutcome.PATHOGENIC
                    ended_at = later.observed_at
                    break
                if later.new_classification is ClassificationState.CONFLICT:
                    outcome = CompetingOutcome.CONFLICT
                    ended_at = later.observed_at
                    break
                if later.new_classification is not ClassificationState.VUS:
                    outcome = CompetingOutcome.OTHER_EXIT
                    ended_at = later.observed_at
                    break

            results.append(
                VusResolution(
                    variation_id=variation_id,
                    started_at=event.observed_at,
                    ended_at=ended_at,
                    duration_days=(ended_at - event.observed_at).days,
                    outcome=outcome,
                    review_status=event.new_review_status,
                    gene_info=event.new_gene_info,
                    conditions=event.new_conditions,
                    variant_type=event.new_variant_type,
                )
            )
    return sorted(
        results,
        key=lambda result: (result.variation_id, result.started_at),
    )


def scv_vus_resolutions(
    events: Iterable[ScvEvent],
    *,
    study_end: date,
    include_left_censored: bool = False,
) -> list[ScvVusResolution]:
    """Measure submitter-level VUS episodes and competing resolutions."""
    event_list = list(events)
    first_snapshot = min(
        (event.observed_at for event in event_list),
        default=None,
    )
    grouped: dict[str, list[ScvEvent]] = {}
    for event in event_list:
        grouped.setdefault(event.scv_accession, []).append(event)
    for history in grouped.values():
        history.sort(key=lambda event: (event.observed_at, event.event_type.value))

    results = []
    for accession, history in grouped.items():
        for index, event in enumerate(history):
            if event.observed_at > study_end:
                break
            if event.event_type not in {
                ScvEventType.SCV_ADDED,
                ScvEventType.SCV_RECLASSIFIED,
            }:
                continue
            if event.new_classification is not ClassificationState.VUS:
                continue
            if event.old_classification is ClassificationState.VUS:
                continue
            if (
                not include_left_censored
                and event.event_type is ScvEventType.SCV_ADDED
                and event.observed_at == first_snapshot
            ):
                continue

            outcome = CompetingOutcome.CENSORED
            ended_at = study_end
            for later in history[index + 1 :]:
                if later.observed_at > study_end:
                    break
                if later.event_type is ScvEventType.SCV_REMOVED:
                    ended_at = later.observed_at
                    break
                if later.event_type is not ScvEventType.SCV_RECLASSIFIED:
                    continue
                if later.new_classification in _BENIGN_STATES:
                    outcome = CompetingOutcome.BENIGN
                    ended_at = later.observed_at
                    break
                if later.new_classification in _PATHOGENIC_STATES:
                    outcome = CompetingOutcome.PATHOGENIC
                    ended_at = later.observed_at
                    break
                if later.new_classification is ClassificationState.CONFLICT:
                    outcome = CompetingOutcome.CONFLICT
                    ended_at = later.observed_at
                    break
                if later.new_classification is not ClassificationState.VUS:
                    outcome = CompetingOutcome.OTHER_EXIT
                    ended_at = later.observed_at
                    break

            results.append(
                ScvVusResolution(
                    scv_accession=accession,
                    variation_id=event.variation_id,
                    started_at=event.observed_at,
                    ended_at=ended_at,
                    duration_days=(ended_at - event.observed_at).days,
                    outcome=outcome,
                    submitter=event.new_submitter,
                    collection_method=event.new_collection_method,
                    review_status=event.new_review_status,
                    conditions=event.new_conditions,
                )
            )
    return sorted(
        results,
        key=lambda result: (result.scv_accession, result.started_at),
    )


def mechanism_counts(
    attributions: Iterable[AggregateAttribution],
) -> dict[AttributionMechanism, int]:
    """Count conservative primary attribution labels."""
    return dict(Counter(item.primary_mechanism for item in attributions))


def _optional_classification(value: str | None) -> ClassificationState | None:
    return ClassificationState(value) if value is not None else None


def read_vcv_event_log(root: Path) -> list[VcvEvent]:
    """Read the partitioned aggregate event log without network access."""
    events = []
    for path in sorted((root / "events" / "vcv").rglob("*.parquet")):
        for row in pq.read_table(path).to_pylist():
            events.append(
                VcvEvent(
                    variation_id=row["variation_id"],
                    event_type=EventType(row["event_type"]),
                    observed_at=row["observed_at"],
                    old_classification=_optional_classification(
                        row["old_classification"]
                    ),
                    new_classification=_optional_classification(
                        row["new_classification"]
                    ),
                    old_raw_classification=row["old_raw_classification"],
                    new_raw_classification=row["new_raw_classification"],
                    old_review_status=row["old_review_status"],
                    new_review_status=row["new_review_status"],
                    old_conditions=row["old_conditions"],
                    new_conditions=row["new_conditions"],
                    old_gene_info=row.get("old_gene_info"),
                    new_gene_info=row.get("new_gene_info"),
                    old_variant_type=row.get("old_variant_type"),
                    new_variant_type=row.get("new_variant_type"),
                )
            )
    return events


def read_scv_event_log(root: Path) -> list[ScvEvent]:
    """Read the partitioned submission event log without network access."""
    events = []
    for path in sorted((root / "events" / "scv").rglob("*.parquet")):
        for row in pq.read_table(path).to_pylist():
            events.append(
                ScvEvent(
                    scv_accession=row["scv_accession"],
                    variation_id=row["variation_id"],
                    event_type=ScvEventType(row["event_type"]),
                    observed_at=row["observed_at"],
                    old_version=row["old_version"],
                    new_version=row["new_version"],
                    old_classification=_optional_classification(
                        row["old_classification"]
                    ),
                    new_classification=_optional_classification(
                        row["new_classification"]
                    ),
                    old_raw_classification=row["old_raw_classification"],
                    new_raw_classification=row["new_raw_classification"],
                    old_review_status=row["old_review_status"],
                    new_review_status=row["new_review_status"],
                    old_conditions=row["old_conditions"],
                    new_conditions=row["new_conditions"],
                    old_submitter=row.get("old_submitter"),
                    new_submitter=row.get("new_submitter"),
                    old_collection_method=row.get("old_collection_method"),
                    new_collection_method=row.get("new_collection_method"),
                )
            )
    return events


def read_attributions(path: Path) -> list[AggregateAttribution]:
    """Read conservative attribution results from Parquet."""
    return [
        AggregateAttribution(
            variation_id=row["variation_id"],
            observed_at=row["observed_at"],
            primary_mechanism=AttributionMechanism(row["primary_mechanism"]),
            candidate_mechanisms=tuple(
                AttributionMechanism(item) for item in row["candidate_mechanisms"]
            ),
            supporting_event_ids=tuple(row["supporting_event_ids"]),
        )
        for row in pq.read_table(path).to_pylist()
    ]
