from datetime import date

import pyarrow as pa
import pyarrow.parquet as pq

from clinvar_history.manifest import ReleaseFile, ReleaseKind, write_manifest
from clinvar_history.research import reconstruct_episodes


def _dataset(root):
    releases = [
        ReleaseFile(date(2020, 1, d), ReleaseKind.VCF, f"fixture://{d}", 1)
        for d in [1, 8, 15, 22]
    ]
    write_manifest(releases, root / "manifest/releases.parquet")
    records = [
        (1, 1, "variant_added", None, "vus", None, "Uncertain significance"),
        (2, 1, "variant_added", None, "benign", None, "Benign"),
        (
            2,
            8,
            "classification_changed",
            "benign",
            "vus",
            "Benign",
            "Uncertain significance",
        ),
        (
            2,
            15,
            "classification_changed",
            "vus",
            "benign",
            "Uncertain significance",
            "Benign",
        ),
        (
            2,
            22,
            "classification_changed",
            "benign",
            "vus",
            "Benign",
            "Uncertain significance",
        ),
        (3, 8, "variant_added", None, "vus", None, "Uncertain significance"),
        (
            3,
            15,
            "classification_changed",
            "vus",
            "other",
            "Uncertain significance",
            "Benign/Likely benign",
        ),
    ]
    rows = [
        {
            "variation_id": v,
            "event_type": t,
            "observed_at": date(2020, 1, d),
            "old_classification": old,
            "new_classification": new,
            "old_raw_classification": ro,
            "new_raw_classification": rn,
            "new_review_status": "single",
            "new_gene_info": "GENE:1",
            "new_variant_type": "SNV",
            "new_conditions": "condition",
        }
        for v, d, t, old, new, ro, rn in records
    ]
    path = root / "events/vcv/release_date=2020-01-01/events.parquet"
    path.parent.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist(rows), path)


def test_independent_episode_reconstruction_preserves_recurrence_and_censors_at_source_end(
    tmp_path,
):
    _dataset(tmp_path)
    output = tmp_path / "episodes.parquet"
    reconstruct_episodes(tmp_path, output, layer="vcv", policy="strict")
    rows = pq.read_table(output).to_pylist()
    assert not any(r["unit"] == "1" for r in rows)  # bootstrap inventory
    second = [r for r in rows if r["unit"] == "2"]
    assert [(r["started_at"], r["ended_at"], r["outcome"]) for r in second] == [
        (date(2020, 1, 8), date(2020, 1, 15), "benign"),
        (date(2020, 1, 22), date(2020, 1, 22), "censored"),
    ]
    assert second[0]["entry_lower"] == date(2020, 1, 1)
    assert second[0]["exit_lower"] == date(2020, 1, 8)
    assert second[1]["censor_reason"] == "administrative"


def test_endpoint_family_sensitivity_does_not_force_composites_into_an_ordinal_state(
    tmp_path,
):
    _dataset(tmp_path)
    strict, family = tmp_path / "strict.parquet", tmp_path / "family.parquet"
    reconstruct_episodes(tmp_path, strict, layer="vcv", policy="strict")
    reconstruct_episodes(tmp_path, family, layer="vcv", policy="family")
    assert (
        next(r for r in pq.read_table(strict).to_pylist() if r["unit"] == "3")[
            "outcome"
        ]
        == "other_exit"
    )
    assert (
        next(r for r in pq.read_table(family).to_pylist() if r["unit"] == "3")[
            "outcome"
        ]
        == "benign"
    )


def test_excluding_an_observation_collapses_the_path_without_dropping_state_updates(
    tmp_path,
):
    _dataset(tmp_path)
    output = tmp_path / "excluded.parquet"
    reconstruct_episodes(
        tmp_path,
        output,
        layer="vcv",
        policy="family",
        excluded_dates=(date(2020, 1, 8),),
    )
    rows = pq.read_table(output).to_pylist()
    assert [(r["unit"], r["started_at"]) for r in rows] == [("2", date(2020, 1, 22))]
