import pytest

from clinvar_history.survival import (
    Observation,
    bootstrap_incidence,
    cumulative_incidence,
)


def test_competing_risk_hand_calculation_and_tied_censoring():
    # At day 1: one benign event, one censor, four at risk. At day 2:
    # one pathogenic and one censor, two at risk. S=3/8, F_B=1/4, F_P=3/8.
    rows = [
        Observation("a", 1, "benign"),
        Observation("b", 1, "censored"),
        Observation("c", 2, "pathogenic"),
        Observation("d", 2, "censored"),
    ]
    points = cumulative_incidence(rows, [0, 1, 2, 3])
    assert points[1]["benign"] == 0.25
    assert points[2]["pathogenic"] == 0.375
    assert points[2]["survival"] == 0.375
    assert points[3]["benign"] is None  # no supported extrapolation
    assert points[2]["at_risk"] == 2


def test_all_censored_and_all_absorbing():
    assert (
        cumulative_incidence([Observation("a", 10, "censored")], [5])[0]["benign"] == 0
    )
    point = cumulative_incidence([Observation("a", 1, "benign")], [10])[0]
    assert point["benign"] == 1
    assert point["survival"] == 0


def test_bootstrap_is_reproducible_and_resamples_whole_variant_clusters():
    rows = [Observation("same", 1, "benign"), Observation("same", 1, "benign")]
    result = bootstrap_incidence(rows, [1], replicates=30, seed=9)
    assert result == bootstrap_incidence(rows, [1], replicates=30, seed=9)
    assert result[0]["benign_lower"] == result[0]["benign_upper"] == 1


@pytest.mark.parametrize(
    "observation", [Observation("a", -1, "benign"), Observation("a", 1, "unknown")]
)
def test_invalid_observations_rejected(observation):
    with pytest.raises(ValueError):
        cumulative_incidence([observation], [1])


def test_first_incident_selection_precedes_cohort_filter_and_clips_followup(tmp_path):
    from datetime import date

    import pyarrow as pa
    import pyarrow.parquet as pq

    from clinvar_history.survival import load_first_episodes

    path = tmp_path / "episodes.parquet"
    rows = [
        {
            "variation_id": 1,
            "started_at": date(2020, 1, 1),
            "ended_at": date(2020, 1, 2),
            "duration_days": 1,
            "outcome": "benign",
            "review_status": "single",
        },
        {
            "variation_id": 1,
            "started_at": date(2020, 2, 1),
            "ended_at": date(2020, 2, 2),
            "duration_days": 1,
            "outcome": "benign",
            "review_status": "single",
        },
        {
            "variation_id": 2,
            "started_at": date(2020, 2, 1),
            "ended_at": date(2020, 3, 1),
            "duration_days": 29,
            "outcome": "pathogenic",
            "review_status": "single",
        },
    ]
    pq.write_table(pa.Table.from_pylist(rows), path)
    result = load_first_episodes(
        path,
        layer="vcv",
        cohort_start=date(2020, 2, 1),
        cohort_end=date(2020, 2, 10),
        study_end=date(2020, 2, 20),
        max_clusters=10,
        seed=3,
    )
    assert [r["variation_id"] for r in result] == [2]
    assert result[0]["ended_at"] == date(2020, 2, 20)
    assert result[0]["outcome"] == "censored"


def test_monthly_snapshots_lose_a_transient_vus_episode(tmp_path):
    from datetime import date

    import pyarrow as pa
    import pyarrow.parquet as pq

    from clinvar_history.survival import snapshot_sensitivity

    path = tmp_path / "events/vcv/release_date=2020-01-01/events.parquet"
    path.parent.mkdir(parents=True)
    rows = [
        {
            "variation_id": vid,
            "event_type": kind,
            "observed_at": date.fromisoformat(day),
            "new_classification": state,
        }
        for vid, kind, day, state in [
            (1, "variant_added", "2020-01-01", "benign"),
            (1, "classification_changed", "2020-01-10", "vus"),
            (1, "classification_changed", "2020-01-20", "benign"),
            (2, "variant_added", "2020-01-01", "benign"),
            (2, "classification_changed", "2020-01-10", "vus"),
            (2, "classification_changed", "2020-02-10", "pathogenic"),
        ]
    ]
    pq.write_table(pa.Table.from_pylist(rows), path)
    result = snapshot_sensitivity(
        tmp_path,
        [{"variation_id": 1, "stratum": "all"}, {"variation_id": 2, "stratum": "all"}],
        release_dates=[date(2020, 1, 1)],
        study_end=date(2020, 3, 1),
        cohort_start=date(2020, 1, 2),
        cohort_end=date(2020, 2, 1),
        months=1,
    )
    assert result == [Observation("2", 29, "pathogenic")]
