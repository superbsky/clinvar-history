from dataclasses import replace
from datetime import date

import pytest
from click.testing import CliRunner

from clinvar_history.attribution import (
    AttributionMechanism,
    attribute_aggregate_event,
    attribute_event_tables,
    exclude_bootstrap_additions,
)
from clinvar_history.classification import ClassificationState
from clinvar_history.cli import main
from clinvar_history.diff import EventType, VcvEvent
from clinvar_history.scv_diff import ScvEvent, ScvEventType
from clinvar_history.storage import write_scv_events, write_vcv_events

OBSERVED_AT = date(2024, 2, 1)


@pytest.fixture
def aggregate_event():
    return VcvEvent(
        variation_id=123,
        event_type=EventType.CLASSIFICATION_CHANGED,
        observed_at=OBSERVED_AT,
        old_classification=ClassificationState.VUS,
        new_classification=ClassificationState.LIKELY_PATHOGENIC,
        old_raw_classification="Uncertain significance",
        new_raw_classification="Likely pathogenic",
        old_review_status="criteria provided, single submitter",
        new_review_status="criteria provided, multiple submitters, no conflicts",
        old_conditions="MONDO:1",
        new_conditions="MONDO:1",
    )


def _scv_event(
    event_type: ScvEventType,
    *,
    accession: str = "SCV000000001",
    review_status: str | None = "criteria provided, single submitter",
) -> ScvEvent:
    return ScvEvent(
        scv_accession=accession,
        variation_id=123,
        event_type=event_type,
        observed_at=OBSERVED_AT,
        old_version=1,
        new_version=2,
        old_classification=ClassificationState.VUS,
        new_classification=ClassificationState.LIKELY_PATHOGENIC,
        old_raw_classification="Uncertain significance",
        new_raw_classification="Likely pathogenic",
        old_review_status="criteria provided, single submitter",
        new_review_status=review_status,
        old_conditions="MONDO:1",
        new_conditions="MONDO:1",
    )


@pytest.mark.parametrize(
    ("event", "mechanism"),
    [
        (
            _scv_event(ScvEventType.SCV_RECLASSIFIED),
            AttributionMechanism.EXACT_MATCHING_SCV_RECLASSIFICATION,
        ),
        (
            _scv_event(ScvEventType.SCV_ADDED),
            AttributionMechanism.SCV_ADDED,
        ),
        (
            _scv_event(ScvEventType.SCV_REMOVED),
            AttributionMechanism.SCV_REMOVED,
        ),
        (
            _scv_event(ScvEventType.REVIEW_STATUS_CHANGED),
            AttributionMechanism.REVIEW_PRECEDENCE_CHANGED,
        ),
        (
            _scv_event(ScvEventType.CONDITION_SET_CHANGED),
            AttributionMechanism.CONDITION_CHANGED,
        ),
    ],
)
def test_single_candidate_gets_specific_primary(
    aggregate_event,
    event,
    mechanism,
):
    result = attribute_aggregate_event(aggregate_event, [event])

    assert result.primary_mechanism is mechanism
    assert result.candidate_mechanisms == (mechanism,)
    assert result.supporting_event_ids == (
        f"{event.scv_accession}:{event.observed_at}:{event.event_type.value}",
    )


def test_expert_panel_addition_is_more_specific_than_generic_addition(
    aggregate_event,
):
    event = _scv_event(
        ScvEventType.SCV_ADDED,
        review_status="reviewed by expert panel",
    )

    result = attribute_aggregate_event(aggregate_event, [event])

    assert result.primary_mechanism is AttributionMechanism.EXPERT_PANEL_ADDED
    assert result.candidate_mechanisms == (AttributionMechanism.EXPERT_PANEL_ADDED,)


def test_multiple_candidate_types_are_not_presented_as_causal(
    aggregate_event,
):
    events = [
        _scv_event(ScvEventType.SCV_RECLASSIFIED),
        _scv_event(
            ScvEventType.SCV_ADDED,
            accession="SCV000000002",
        ),
    ]

    result = attribute_aggregate_event(aggregate_event, events)

    assert result.primary_mechanism is AttributionMechanism.MULTIPLE_CANDIDATES
    assert result.candidate_mechanisms == (
        AttributionMechanism.EXACT_MATCHING_SCV_RECLASSIFICATION,
        AttributionMechanism.SCV_ADDED,
    )
    assert len(result.supporting_event_ids) == 2


def test_unexplained_when_only_neutral_or_unrelated_events_exist(
    aggregate_event,
):
    neutral = _scv_event(ScvEventType.SCV_VERSION_CHANGED)
    unrelated = replace(
        _scv_event(ScvEventType.SCV_RECLASSIFIED),
        variation_id=999,
    )

    result = attribute_aggregate_event(aggregate_event, [neutral, unrelated])

    assert result.primary_mechanism is AttributionMechanism.UNEXPLAINED
    assert result.candidate_mechanisms == ()
    assert result.supporting_event_ids == ()


def test_direction_incompatible_scv_reclassification_is_not_supporting(
    aggregate_event,
):
    incompatible = replace(
        _scv_event(ScvEventType.SCV_RECLASSIFIED),
        old_classification=ClassificationState.BENIGN,
        new_classification=ClassificationState.LIKELY_BENIGN,
        old_raw_classification="Benign",
        new_raw_classification="Likely benign",
    )

    result = attribute_aggregate_event(aggregate_event, [incompatible])

    assert result.primary_mechanism is AttributionMechanism.UNEXPLAINED
    assert result.supporting_event_ids == ()


def test_first_snapshot_additions_are_not_attribution_candidates():
    bootstrap = replace(
        _scv_event(ScvEventType.SCV_ADDED),
        observed_at=date(2024, 1, 1),
    )
    later_addition = replace(
        _scv_event(ScvEventType.SCV_ADDED, accession="SCV000000002"),
        observed_at=date(2024, 2, 1),
    )

    filtered = exclude_bootstrap_additions([bootstrap, later_addition])

    assert filtered == [later_addition]


def test_attribution_cli_is_available():
    result = CliRunner().invoke(main, ["attribute-events", "--help"])

    assert result.exit_code == 0
    assert "--root" in result.output
    assert "--output" in result.output


def test_table_attribution_requires_overlapping_observation_intervals(
    tmp_path,
    aggregate_event,
):
    bootstrap_vcv = replace(
        aggregate_event,
        event_type=EventType.VARIANT_ADDED,
        observed_at=date(2024, 1, 24),
        old_classification=None,
        old_raw_classification=None,
    )
    early_vcv = replace(aggregate_event, observed_at=date(2024, 1, 31))
    overlapping_vcv = replace(aggregate_event, observed_at=date(2024, 2, 24))
    bootstrap_scv = replace(
        _scv_event(ScvEventType.SCV_ADDED),
        observed_at=date(2024, 2, 1),
        old_version=None,
        old_classification=None,
        old_raw_classification=None,
    )
    later_scv = replace(
        _scv_event(ScvEventType.SCV_RECLASSIFIED),
        observed_at=date(2024, 3, 1),
    )

    for observed_at, events in [
        (date(2024, 1, 24), [bootstrap_vcv]),
        (date(2024, 1, 31), [early_vcv]),
        (date(2024, 2, 24), [overlapping_vcv]),
    ]:
        write_vcv_events(
            events,
            tmp_path
            / "events"
            / "vcv"
            / f"release_date={observed_at}"
            / "events.parquet",
        )
    for observed_at, events in [
        (date(2024, 2, 1), [bootstrap_scv]),
        (date(2024, 3, 1), [later_scv]),
    ]:
        write_scv_events(
            events,
            tmp_path
            / "events"
            / "scv"
            / f"release_date={observed_at}"
            / "events.parquet",
        )

    results = attribute_event_tables(
        tmp_path,
        tmp_path / "attributions.parquet",
        max_lag_days=35,
    )

    assert [row.primary_mechanism for row in results] == [
        AttributionMechanism.UNEXPLAINED,
        AttributionMechanism.EXACT_MATCHING_SCV_RECLASSIFICATION,
    ]
