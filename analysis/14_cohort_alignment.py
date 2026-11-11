"""Align the same first episodes across cadences and describe recurrent VUS entry."""

import argparse
from pathlib import Path

import pyarrow as pa

from clinvar_history.checkpoints import atomic_json, file_digest, publish_checkpoint
from clinvar_history.research import _copy, research_connection
from clinvar_history.research_statistics import estimate_counts, load_cluster_counts


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path("results/full-research"))
    args = p.parse_args()
    output = args.root / "alignment"
    output.mkdir(parents=True, exist_ok=True)
    summaries = []
    for layer in ["vcv", "scv"]:
        native = args.root / f"episodes_{layer}_family_0.parquet"
        monthly = args.root / f"episodes_{layer}_family_1.parquet"
        quarterly = args.root / f"episodes_{layer}_family_3.parquet"
        aligned = output / f"{layer}_aligned_units.parquet"
        with research_connection() as c:
            _copy(
                c,
                """SELECT n.unit FROM read_parquet(?) n
                JOIN read_parquet(?) m USING(unit) JOIN read_parquet(?) q USING(unit)
                WHERE n.incident_rank=1 AND m.incident_rank=1 AND q.incident_rank=1
                  AND NOT n.identity_migrated AND NOT m.identity_migrated AND NOT q.identity_migrated
                  AND m.started_at>=n.started_at AND (m.started_at<n.ended_at OR (n.outcome='censored' AND m.started_at=n.ended_at))
                  AND q.started_at>=n.started_at AND (q.started_at<n.ended_at OR (n.outcome='censored' AND q.started_at=n.ended_at))""",
                aligned,
                [str(native), str(monthly), str(quarterly)],
            )
            _copy(
                c,
                """SELECT year(started_at) AS entry_year,entry_origin,outcome,
                count(*) AS episodes,count(DISTINCT unit) AS units,
                count(*) FILTER(WHERE incident_rank>1) AS recurrent_episodes,
                count(*) FILTER(WHERE censor_reason='removal') AS removals
                FROM read_parquet(?) GROUP BY ALL ORDER BY entry_year,entry_origin,outcome""",
                output / f"{layer}_recurrence.parquet",
                [str(native)],
            )
            n_aligned = c.execute(
                "SELECT count(*) FROM read_parquet(?)", [str(aligned)]
            ).fetchone()[0]
            total, units, recurrent, removed, migrated = c.execute(
                """SELECT count(*),count(DISTINCT unit),count(*) FILTER(WHERE incident_rank>1),count(*) FILTER(WHERE censor_reason='removal'),count(*) FILTER(WHERE identity_migrated) FROM read_parquet(?)""",
                [str(native)],
            ).fetchone()
            summaries.append(
                {
                    "layer": layer,
                    "all_episodes": total,
                    "units": units,
                    "recurrent_episodes": recurrent,
                    "removal_censors": removed,
                    "identity_migrated_episodes": migrated,
                    "aligned_first_episodes": n_aligned,
                }
            )
        rows = []
        curves = []
        for name, source, where, params in [
            (
                "aligned_native",
                native,
                "incident_rank=1 AND unit IN (SELECT unit FROM read_parquet(?))",
                [str(aligned)],
            ),
            (
                "aligned_monthly",
                monthly,
                "incident_rank=1 AND unit IN (SELECT unit FROM read_parquet(?))",
                [str(aligned)],
            ),
            (
                "aligned_quarterly",
                quarterly,
                "incident_rank=1 AND unit IN (SELECT unit FROM read_parquet(?))",
                [str(aligned)],
            ),
            ("all_episodes_clustered", native, "NOT identity_migrated", []),
        ]:
            print(layer, name, flush=True)
            data = load_cluster_counts(source, where=where, parameters=params)
            point, curve, _ = estimate_counts(
                data, [365, 730, 1095, 1825], replicates=1000, seed=20261111
            )
            for r in point:
                r.update({"layer": layer, "analysis": name})
            rows.extend(point)
            curves.append(
                pa.table(
                    {
                        **curve,
                        "layer": [layer] * len(curve["day"]),
                        "analysis": [name] * len(curve["day"]),
                    }
                )
            )
        identity = {
            "stage": "cohort-alignment",
            "source": {str(p): file_digest(p) for p in [native, monthly, quarterly]},
            "implementation": {
                str(p): file_digest(p)
                for p in [
                    Path(__file__),
                    Path("src/clinvar_history/research_statistics.py"),
                ]
            },
            "replicates": 1000,
            "seed": 20261111,
        }
        publish_checkpoint(
            pa.Table.from_pylist(rows), output / f"{layer}_incidence.parquet", identity
        )
        publish_checkpoint(
            pa.concat_tables(curves),
            output / f"{layer}_curves.parquet",
            {**identity, "artifact": "curves"},
        )
    atomic_json(summaries, output / "summary.json")


if __name__ == "__main__":
    main()
