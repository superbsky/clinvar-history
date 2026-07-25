from dataclasses import replace
from datetime import date
from pathlib import Path

import pytest

from clinvar_history import analysis_core
from clinvar_history.analysis_core import (
    CompetingOutcome,
    cohort_vus_metrics,
    mechanism_counts,
    month_grid,
    renovo_quarterly_metrics,
    sampling_bias,
    scv_vus_resolutions,
    state_episodes,
    transition_counts,
    vus_resolutions,
)
from clinvar_history.attribution import (
    AggregateAttribution,
    AttributionMechanism,
)
from clinvar_history.classification import ClassificationState
from clinvar_history.diff import EventType, VcvEvent
from clinvar_history.scv_diff import ScvEvent, ScvEventType


def _event(
    variation_id: int,
    observed_at: date,
    old: ClassificationState | None,
    new: ClassificationState | None,
) -> VcvEvent:
    if old is None:
        event_type = EventType.VARIANT_ADDED
    elif new is None:
        event_type = EventType.VARIANT_REMOVED
    else:
        event_type = EventType.CLASSIFICATION_CHANGED
    return VcvEvent(
        variation_id=variation_id,
        event_type=event_type,
        observed_at=observed_at,
        old_classification=old,
        new_classification=new,
        old_raw_classification=old.value if old else None,
        new_raw_classification=new.value if new else None,
        old_review_status=None,
        new_review_status=None,
        old_conditions=None,
        new_conditions=None,
    )


@pytest.fixture
def history():
    vus = ClassificationState.VUS
    benign = ClassificationState.BENIGN
    pathogenic = ClassificationState.PATHOGENIC
    likely_pathogenic = ClassificationState.LIKELY_PATHOGENIC
    return [
        _event(1, date(2024, 1, 1), None, vus),
        _event(2, date(2024, 1, 1), None, vus),
        _event(3, date(2024, 1, 1), None, benign),
        _event(4, date(2024, 1, 1), None, vus),
        _event(1, date(2024, 2, 1), vus, likely_pathogenic),
        _event(2, date(2024, 2, 1), vus, benign),
        _event(3, date(2024, 2, 15), benign, pathogenic),
        _event(2, date(2024, 3, 1), benign, vus),
        _event(3, date(2024, 3, 15), pathogenic, benign),
    ]


def test_higher_cadence_recovers_transient_and_repeated_changes(history):
    quarterly = sampling_bias(
        history,
        [date(2024, 1, 1), date(2024, 4, 1)],
    )
    monthly = sampling_bias(
        history,
        [
            date(2024, 1, 1),
            date(2024, 2, 1),
            date(2024, 3, 1),
            date(2024, 4, 1),
        ],
    )

    assert quarterly.true_transitions == 5
    assert quarterly.observed_transitions == 1
    assert quarterly.event_recall == pytest.approx(0.2)
    assert quarterly.missed_recurrent_transitions == 4
    assert quarterly.missed_conflict_episodes == 0
    assert quarterly.mean_event_date_error_days == 60
    assert monthly.true_transitions == 5
    assert monthly.observed_transitions == 5
    assert monthly.event_recall == 1


def test_sampling_bias_counts_transient_conflict_episodes():
    vus = ClassificationState.VUS
    conflict = ClassificationState.CONFLICT
    events = [
        _event(10, date(2024, 1, 1), None, vus),
        _event(10, date(2024, 2, 1), vus, conflict),
        _event(10, date(2024, 3, 1), conflict, vus),
    ]

    result = sampling_bias(
        events,
        [date(2024, 1, 1), date(2024, 4, 1)],
    )

    assert result.missed_conflict_episodes == 1


def test_sampling_bias_sweeps_events_without_rebuilding_full_snapshots(
    monkeypatch, history
):
    def fail_if_called(*_args, **_kwargs):
        raise AssertionError("sampling_bias must not rebuild full snapshots")

    monkeypatch.setattr(analysis_core, "_snapshot", fail_if_called)

    result = sampling_bias(
        history,
        [
            date(2024, 1, 1),
            date(2024, 2, 1),
            date(2024, 3, 1),
            date(2024, 4, 1),
        ],
    )

    assert result.true_transitions == 5
    assert result.observed_transitions == 5
    assert result.event_recall == 1


def test_cohort_vus_metric_uses_starting_risk_set(history):
    metrics = cohort_vus_metrics(
        history,
        [date(2024, 1, 1), date(2024, 4, 1)],
    )

    assert len(metrics) == 1
    quarter = metrics[0]
    assert quarter.vus_at_start == 3
    assert quarter.vus_resolved == 1
    assert quarter.new_vus == 0
    assert quarter.reclassification_rate == pytest.approx(1 / 3)
    assert quarter.discovery_rate == 0


def test_renovo_baseline_uses_snapshot_formula_exactly(history):
    metrics = renovo_quarterly_metrics(
        history,
        [date(2024, 1, 1), date(2024, 4, 1)],
    )

    quarter = metrics[0]
    assert quarter.vus_previous == 3
    assert quarter.vus_current == 2
    assert quarter.reclassified_vus == 1
    assert quarter.net_vus_change == -1
    assert quarter.reclassification_rate == pytest.approx(1 / 2)
    assert quarter.discovery_rate == pytest.approx(-1 / 2)
    assert quarter.discovery_reclassification_ratio == pytest.approx(-1)


def test_transition_matrix_counts_every_transition_including_recurrence(history):
    counts = {
        (row.old_classification, row.new_classification): row.count
        for row in transition_counts(history)
    }

    assert counts[(ClassificationState.VUS, ClassificationState.LIKELY_PATHOGENIC)] == 1
    assert counts[(ClassificationState.VUS, ClassificationState.BENIGN)] == 1
    assert counts[(ClassificationState.BENIGN, ClassificationState.VUS)] == 1
    assert counts[(ClassificationState.BENIGN, ClassificationState.PATHOGENIC)] == 1
    assert counts[(ClassificationState.PATHOGENIC, ClassificationState.BENIGN)] == 1


def test_canonical_analyses_ignore_raw_only_label_changes(history):
    raw_only = _event(
        4,
        date(2024, 2, 10),
        ClassificationState.VUS,
        ClassificationState.VUS,
    )
    raw_only = replace(
        raw_only,
        old_raw_classification="VUS",
        new_raw_classification="Uncertain significance",
    )
    augmented = [*history, raw_only]

    bias = sampling_bias(
        augmented,
        [
            date(2024, 1, 1),
            date(2024, 2, 1),
            date(2024, 3, 1),
            date(2024, 4, 1),
        ],
    )
    counts = transition_counts(augmented)
    vus_episodes = [
        episode
        for episode in state_episodes(augmented, study_end=date(2024, 4, 1))
        if episode.variation_id == 4
    ]

    assert bias.true_transitions == 5
    assert all(row.old_classification != row.new_classification for row in counts)
    assert len(vus_episodes) == 1
    assert vus_episodes[0].duration_days == 91


def test_dwell_times_and_competing_vus_outcomes(history):
    episodes = state_episodes(history, study_end=date(2024, 4, 1))
    first_vus = next(
        episode
        for episode in episodes
        if episode.variation_id == 1
        and episode.classification is ClassificationState.VUS
    )
    resolutions = vus_resolutions(
        history,
        study_end=date(2024, 4, 1),
        include_left_censored=True,
    )
    outcomes = {(row.variation_id, row.started_at): row for row in resolutions}

    assert first_vus.duration_days == 31
    assert first_vus.censored is False
    assert outcomes[(1, date(2024, 1, 1))].outcome is CompetingOutcome.PATHOGENIC
    assert outcomes[(1, date(2024, 1, 1))].duration_days == 31
    assert outcomes[(2, date(2024, 1, 1))].outcome is CompetingOutcome.BENIGN
    assert outcomes[(2, date(2024, 3, 1))].outcome is CompetingOutcome.CENSORED


def test_vus_resolution_excludes_bootstrap_left_censored_episode(history):
    rows = vus_resolutions(history, study_end=date(2024, 4, 1))

    assert {(row.variation_id, row.started_at) for row in rows} == {
        (2, date(2024, 3, 1)),
    }


def test_vus_resolution_records_conflict_as_observed_exit():
    benign = ClassificationState.BENIGN
    vus = ClassificationState.VUS
    conflict = ClassificationState.CONFLICT
    events = [
        _event(10, date(2024, 1, 1), None, benign),
        _event(10, date(2024, 2, 1), benign, vus),
        _event(10, date(2024, 3, 1), vus, conflict),
    ]

    rows = vus_resolutions(events, study_end=date(2024, 4, 1))

    assert len(rows) == 1
    assert rows[0].outcome is CompetingOutcome.CONFLICT
    assert rows[0].ended_at == date(2024, 3, 1)


def test_vus_resolution_ignores_parallel_review_event():
    benign = ClassificationState.BENIGN
    vus = ClassificationState.VUS
    added = _event(9, date(2024, 1, 1), None, benign)
    changed = replace(
        _event(9, date(2024, 2, 1), benign, vus),
        new_review_status="reviewed by expert panel",
        new_gene_info="GENE9:9",
        new_conditions="MONDO:9",
        new_variant_type="single nucleotide variant",
    )
    parallel_review = VcvEvent(
        variation_id=9,
        event_type=EventType.REVIEW_STATUS_CHANGED,
        observed_at=date(2024, 2, 1),
        old_classification=benign,
        new_classification=vus,
        old_raw_classification="Benign",
        new_raw_classification="Uncertain significance",
        old_review_status="old",
        new_review_status="new",
        old_conditions=None,
        new_conditions=None,
    )

    rows = vus_resolutions(
        [added, changed, parallel_review],
        study_end=date(2024, 4, 1),
    )

    assert len(rows) == 1
    assert rows[0].started_at == date(2024, 2, 1)
    assert rows[0].review_status == "reviewed by expert panel"
    assert rows[0].gene_info == "GENE9:9"
    assert rows[0].conditions == "MONDO:9"
    assert rows[0].variant_type == "single nucleotide variant"


def test_scv_vus_resolution_retains_submitter_and_censors_removal():
    vus = ClassificationState.VUS
    pathogenic = ClassificationState.PATHOGENIC
    events = [
        ScvEvent(
            scv_accession="SCV1",
            variation_id=9,
            event_type=ScvEventType.SCV_ADDED,
            observed_at=date(2023, 12, 1),
            old_version=None,
            new_version=1,
            old_classification=None,
            new_classification=pathogenic,
            old_raw_classification=None,
            new_raw_classification="Pathogenic",
            old_review_status=None,
            new_review_status="criteria provided",
            old_conditions=None,
            new_conditions="MONDO:9",
            new_submitter="Lab A",
            new_collection_method="clinical testing",
        ),
        ScvEvent(
            scv_accession="SCV1",
            variation_id=9,
            event_type=ScvEventType.SCV_RECLASSIFIED,
            observed_at=date(2024, 1, 1),
            old_version=1,
            new_version=2,
            old_classification=pathogenic,
            new_classification=vus,
            old_raw_classification="Pathogenic",
            new_raw_classification="Uncertain significance",
            old_review_status="criteria provided",
            new_review_status="criteria provided",
            old_conditions="MONDO:9",
            new_conditions="MONDO:9",
            old_submitter="Lab A",
            new_submitter="Lab A",
            old_collection_method="clinical testing",
            new_collection_method="clinical testing",
        ),
        ScvEvent(
            scv_accession="SCV1",
            variation_id=9,
            event_type=ScvEventType.SCV_RECLASSIFIED,
            observed_at=date(2024, 2, 1),
            old_version=2,
            new_version=3,
            old_classification=vus,
            new_classification=pathogenic,
            old_raw_classification="Uncertain significance",
            new_raw_classification="Pathogenic",
            old_review_status="criteria provided",
            new_review_status="criteria provided",
            old_conditions="MONDO:9",
            new_conditions="MONDO:9",
            old_submitter="Lab A",
            new_submitter="Lab A",
            old_collection_method="clinical testing",
            new_collection_method="clinical testing",
        ),
        ScvEvent(
            scv_accession="SCV1",
            variation_id=9,
            event_type=ScvEventType.SCV_RECLASSIFIED,
            observed_at=date(2024, 3, 1),
            old_version=3,
            new_version=4,
            old_classification=pathogenic,
            new_classification=vus,
            old_raw_classification="Pathogenic",
            new_raw_classification="Uncertain significance",
            old_review_status="criteria provided",
            new_review_status="criteria provided",
            old_conditions="MONDO:9",
            new_conditions="MONDO:9",
            old_submitter="Lab A",
            new_submitter="Lab B",
            old_collection_method="clinical testing",
            new_collection_method="research",
        ),
        ScvEvent(
            scv_accession="SCV1",
            variation_id=9,
            event_type=ScvEventType.SCV_REMOVED,
            observed_at=date(2024, 4, 1),
            old_version=4,
            new_version=None,
            old_classification=vus,
            new_classification=None,
            old_raw_classification="Uncertain significance",
            new_raw_classification=None,
            old_review_status="criteria provided",
            new_review_status=None,
            old_conditions="MONDO:9",
            new_conditions=None,
            old_submitter="Lab B",
            new_submitter=None,
            old_collection_method="research",
            new_collection_method=None,
        ),
    ]
    events.extend(
        [
            replace(
                events[1],
                scv_accession="SCV2",
                variation_id=10,
                observed_at=date(2024, 1, 10),
            ),
            replace(
                events[2],
                scv_accession="SCV2",
                variation_id=10,
                observed_at=date(2024, 2, 10),
                new_classification=ClassificationState.CONFLICT,
                new_raw_classification="Conflicting classifications",
            ),
        ]
    )

    rows = scv_vus_resolutions(events, study_end=date(2024, 5, 1))

    assert [
        (row.scv_accession, row.started_at, row.outcome, row.ended_at) for row in rows
    ] == [
        (
            "SCV1",
            date(2024, 1, 1),
            CompetingOutcome.PATHOGENIC,
            date(2024, 2, 1),
        ),
        (
            "SCV1",
            date(2024, 3, 1),
            CompetingOutcome.CENSORED,
            date(2024, 4, 1),
        ),
        (
            "SCV2",
            date(2024, 1, 10),
            CompetingOutcome.CONFLICT,
            date(2024, 2, 10),
        ),
    ]
    assert rows[0].submitter == "Lab A"
    assert rows[0].collection_method == "clinical testing"
    assert rows[1].submitter == "Lab B"
    assert rows[1].collection_method == "research"


def test_mechanism_summary_counts_primary_labels():
    rows = [
        AggregateAttribution(
            variation_id=1,
            observed_at=date(2024, 1, 1),
            primary_mechanism=(
                AttributionMechanism.EXACT_MATCHING_SCV_RECLASSIFICATION
            ),
            candidate_mechanisms=(
                AttributionMechanism.EXACT_MATCHING_SCV_RECLASSIFICATION,
            ),
            supporting_event_ids=("a",),
        ),
        AggregateAttribution(
            variation_id=2,
            observed_at=date(2024, 1, 1),
            primary_mechanism=AttributionMechanism.UNEXPLAINED,
            candidate_mechanisms=(),
            supporting_event_ids=(),
        ),
        AggregateAttribution(
            variation_id=3,
            observed_at=date(2024, 1, 1),
            primary_mechanism=AttributionMechanism.UNEXPLAINED,
            candidate_mechanisms=(),
            supporting_event_ids=(),
        ),
    ]

    assert mechanism_counts(rows) == {
        AttributionMechanism.EXACT_MATCHING_SCV_RECLASSIFICATION: 1,
        AttributionMechanism.UNEXPLAINED: 2,
    }


def test_analysis_scripts_are_thin_offline_consumers():
    analysis_dir = Path(__file__).parents[1] / "analysis"
    expected = {
        "01_reproduce_renovo.py",
        "02_sampling_bias.py",
        "03_transition_matrix.py",
        "04_time_to_resolution.py",
        "05_event_mechanisms.py",
        "06_scv_time_to_resolution.py",
    }

    scripts = {path.name for path in analysis_dir.glob("*.py")}
    assert expected <= scripts
    for name in expected:
        source = (analysis_dir / name).read_text()
        assert "httpx" not in source
        assert "ftp.ncbi" not in source


def test_month_grid_preserves_release_day_when_calendar_allows():
    assert month_grid(
        date(2020, 1, 6),
        date(2020, 7, 6),
        3,
    ) == [
        date(2020, 1, 6),
        date(2020, 4, 6),
        date(2020, 7, 6),
    ]
