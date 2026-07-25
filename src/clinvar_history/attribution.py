"""Conservative candidate attribution for aggregate ClinVar changes."""

from collections.abc import Iterable, Iterator, Sequence
from dataclasses import dataclass
from datetime import date
from enum import StrEnum
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from .classification import ClassificationState
from .diff import EventType, VcvEvent
from .scv_diff import ScvEvent, ScvEventType


class AttributionMechanism(StrEnum):
    EXACT_MATCHING_SCV_RECLASSIFICATION = "exact_matching_scv_reclassification"
    EXPERT_PANEL_ADDED = "expert_panel_added"
    SCV_ADDED = "scv_added"
    SCV_REMOVED = "scv_removed"
    REVIEW_PRECEDENCE_CHANGED = "review_precedence_changed"
    CONDITION_CHANGED = "condition_changed"
    MULTIPLE_CANDIDATES = "multiple_candidates"
    UNEXPLAINED = "unexplained"


_MECHANISM_ORDER = (
    AttributionMechanism.EXACT_MATCHING_SCV_RECLASSIFICATION,
    AttributionMechanism.EXPERT_PANEL_ADDED,
    AttributionMechanism.SCV_ADDED,
    AttributionMechanism.SCV_REMOVED,
    AttributionMechanism.REVIEW_PRECEDENCE_CHANGED,
    AttributionMechanism.CONDITION_CHANGED,
)


@dataclass(frozen=True, slots=True)
class AggregateAttribution:
    variation_id: int
    observed_at: date
    primary_mechanism: AttributionMechanism
    candidate_mechanisms: tuple[AttributionMechanism, ...]
    supporting_event_ids: tuple[str, ...]


def scv_event_id(event: ScvEvent) -> str:
    """Return a stable identifier for an SCV event partition row."""
    return f"{event.scv_accession}:{event.observed_at}:{event.event_type.value}"


def _candidate(
    event: ScvEvent,
    aggregate_event: VcvEvent,
) -> AttributionMechanism | None:
    if event.event_type is ScvEventType.SCV_RECLASSIFIED:
        if (
            event.old_classification != aggregate_event.old_classification
            or event.new_classification != aggregate_event.new_classification
            or event.old_classification == event.new_classification
        ):
            return None
        return AttributionMechanism.EXACT_MATCHING_SCV_RECLASSIFICATION
    if event.event_type is ScvEventType.SCV_ADDED:
        if "expert panel" in (event.new_review_status or "").lower():
            return AttributionMechanism.EXPERT_PANEL_ADDED
        return AttributionMechanism.SCV_ADDED
    if event.event_type is ScvEventType.SCV_REMOVED:
        return AttributionMechanism.SCV_REMOVED
    if event.event_type is ScvEventType.REVIEW_STATUS_CHANGED:
        return AttributionMechanism.REVIEW_PRECEDENCE_CHANGED
    if event.event_type is ScvEventType.CONDITION_SET_CHANGED:
        return AttributionMechanism.CONDITION_CHANGED
    return None


def attribute_aggregate_event(
    aggregate_event: VcvEvent,
    scv_events: Iterable[ScvEvent],
) -> AggregateAttribution:
    """Report candidate mechanisms, avoiding a causal label when they compete."""
    supporting: list[tuple[AttributionMechanism, str]] = []
    for event in scv_events:
        if event.variation_id != aggregate_event.variation_id:
            continue
        mechanism = _candidate(event, aggregate_event)
        if mechanism is not None:
            supporting.append((mechanism, scv_event_id(event)))

    present = {mechanism for mechanism, _ in supporting}
    candidates = tuple(
        mechanism for mechanism in _MECHANISM_ORDER if mechanism in present
    )
    if not candidates:
        primary = AttributionMechanism.UNEXPLAINED
    elif len(candidates) == 1:
        primary = candidates[0]
    else:
        primary = AttributionMechanism.MULTIPLE_CANDIDATES

    return AggregateAttribution(
        variation_id=aggregate_event.variation_id,
        observed_at=aggregate_event.observed_at,
        primary_mechanism=primary,
        candidate_mechanisms=candidates,
        supporting_event_ids=tuple(sorted({event_id for _, event_id in supporting})),
    )


def exclude_bootstrap_additions(events: Sequence[ScvEvent]) -> list[ScvEvent]:
    """Exclude first-snapshot inventory additions from mechanism candidates."""
    if not events:
        return []
    first_snapshot = min(event.observed_at for event in events)
    return [
        event
        for event in events
        if not (
            event.observed_at == first_snapshot
            and event.event_type is ScvEventType.SCV_ADDED
        )
    ]


def _optional_classification(value: str | None) -> ClassificationState | None:
    return ClassificationState(value) if value is not None else None


def _load_vcv_events(root: Path) -> Iterator[VcvEvent]:
    for path in sorted((root / "events" / "vcv").rglob("*.parquet")):
        for row in pq.read_table(path).to_pylist():
            yield VcvEvent(
                variation_id=row["variation_id"],
                event_type=EventType(row["event_type"]),
                observed_at=row["observed_at"],
                old_classification=_optional_classification(row["old_classification"]),
                new_classification=_optional_classification(row["new_classification"]),
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


def _load_scv_events(root: Path) -> Iterator[ScvEvent]:
    for path in sorted((root / "events" / "scv").rglob("*.parquet")):
        for row in pq.read_table(path).to_pylist():
            yield ScvEvent(
                scv_accession=row["scv_accession"],
                variation_id=row["variation_id"],
                event_type=ScvEventType(row["event_type"]),
                observed_at=row["observed_at"],
                old_version=row["old_version"],
                new_version=row["new_version"],
                old_classification=_optional_classification(row["old_classification"]),
                new_classification=_optional_classification(row["new_classification"]),
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


def _snapshot_predecessors(root: Path, layer: str) -> dict[date, date]:
    dates = sorted(
        {
            date.fromisoformat(path.parent.name.removeprefix("release_date="))
            for path in (root / "events" / layer).rglob("*.parquet")
        }
    )
    return dict(zip(dates[1:], dates[:-1], strict=True))


def _observation_intervals_overlap(
    first_start: date,
    first_end: date,
    second_start: date,
    second_end: date,
) -> bool:
    """Return whether two interval-censored windows (start, end] overlap."""
    return max(first_start, second_start) < min(first_end, second_end)


def attribute_event_tables(
    root: Path,
    output: Path,
    *,
    max_lag_days: int = 35,
) -> list[AggregateAttribution]:
    """Attribute VCV changes using overlapping VCV/SCV observation intervals."""
    vcv_predecessors = _snapshot_predecessors(root, "vcv")
    scv_predecessors = _snapshot_predecessors(root, "scv")
    vcv_events = [
        event
        for event in _load_vcv_events(root)
        if event.event_type is EventType.CLASSIFICATION_CHANGED
        and event.old_classification != event.new_classification
    ]
    scv_by_variation: dict[int, list[ScvEvent]] = {}
    first_scv_snapshot: date | None = None
    for event in _load_scv_events(root):
        if first_scv_snapshot is None:
            first_scv_snapshot = event.observed_at
        if (
            event.observed_at == first_scv_snapshot
            and event.event_type is ScvEventType.SCV_ADDED
        ):
            continue
        scv_by_variation.setdefault(event.variation_id, []).append(event)

    attributions = []
    for event in vcv_events:
        vcv_interval_start = vcv_predecessors.get(event.observed_at)
        nearby = [
            candidate
            for candidate in scv_by_variation.get(event.variation_id, [])
            if abs((candidate.observed_at - event.observed_at).days) <= max_lag_days
            and vcv_interval_start is not None
            and (scv_interval_start := scv_predecessors.get(candidate.observed_at))
            is not None
            and _observation_intervals_overlap(
                vcv_interval_start,
                event.observed_at,
                scv_interval_start,
                candidate.observed_at,
            )
        ]
        attributions.append(attribute_aggregate_event(event, nearby))

    _write_attributions(attributions, output)
    return attributions


def _write_attributions(
    attributions: Sequence[AggregateAttribution],
    output: Path,
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    table = pa.Table.from_pylist(
        [
            {
                "variation_id": item.variation_id,
                "observed_at": item.observed_at,
                "primary_mechanism": item.primary_mechanism.value,
                "candidate_mechanisms": [
                    mechanism.value for mechanism in item.candidate_mechanisms
                ],
                "supporting_event_ids": list(item.supporting_event_ids),
            }
            for item in attributions
        ],
        schema=pa.schema(
            [
                ("variation_id", pa.int64()),
                ("observed_at", pa.date32()),
                ("primary_mechanism", pa.string()),
                ("candidate_mechanisms", pa.list_(pa.string())),
                ("supporting_event_ids", pa.list_(pa.string())),
            ]
        ),
    )
    pq.write_table(table, output, compression="zstd")
