"""Full-population first-incident cumulative incidence with exact cluster bootstrap."""

import argparse
import hashlib
from datetime import timedelta
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from clinvar_history.checkpoints import (
    atomic_json,
    checkpoint_ready,
    digest,
    file_digest,
    publish_checkpoint,
)
from clinvar_history.research import research_connection
from clinvar_history.research_statistics import estimate_counts, load_cluster_counts


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("results/full-research"))
    parser.add_argument("--replicates", type=int, default=1000)
    parser.add_argument("--seed", type=int, default=20261111)
    parser.add_argument("--primary-only", action="store_true")
    args = parser.parse_args()
    destination = args.root / "incidence"
    destination.mkdir(parents=True, exist_ok=True)
    implementation = digest(
        {
            str(p): file_digest(p)
            for p in [
                Path(__file__),
                Path("src/clinvar_history/research_statistics.py"),
            ]
        }
    )
    jobs = []

    def add(
        layer,
        name,
        *,
        policy="family",
        cadence=0,
        where="incident_rank=1 AND NOT identity_migrated",
        parameters=(),
        timing="observed",
        removal=False,
        source=None,
        stratum="all",
        value="all",
    ):
        jobs.append(
            {
                "layer": layer,
                "name": name,
                "source": str(
                    source or args.root / f"episodes_{layer}_{policy}_{cadence}.parquet"
                ),
                "policy": policy,
                "cadence": cadence,
                "where": where,
                "parameters": list(parameters),
                "timing": timing,
                "removal": removal,
                "stratum": stratum,
                "value": str(value),
            }
        )

    for layer in ["vcv", "scv"]:
        add(layer, "primary")
        if args.primary_only:
            continue
        add(layer, "strict_labels", policy="strict")
        add(layer, "monthly", cadence=1)
        add(layer, "quarterly", cadence=3)
        for timing in ["midpoint", "lower", "upper"]:
            add(layer, "timing_" + timing, timing=timing)
        add(layer, "removal_competing", removal=True)
        if layer == "scv":
            add(layer, "include_identity_migration", where="incident_rank=1")
        if layer == "vcv":
            add(
                layer,
                "exclude_20220812",
                source=args.root / "episodes_vcv_family_no20220812.parquet",
            )
        path = args.root / f"episodes_{layer}_family_0.parquet"
        with research_connection() as c:
            end = c.execute(
                "SELECT max(ended_at) FROM read_parquet(?)", [str(path)]
            ).fetchone()[0]
            for horizon in [365, 1095, 1825]:
                add(
                    layer,
                    f"mature_{horizon}",
                    where="incident_rank=1 AND NOT identity_migrated AND started_at<=?::DATE",
                    parameters=[str(end - timedelta(days=horizon))],
                    stratum="minimum_administrative_followup",
                    value=horizon,
                )
            for year in range(2020, 2027):
                add(
                    layer,
                    f"entry_{year}",
                    where="incident_rank=1 AND NOT identity_migrated AND year(started_at)=?",
                    parameters=[year],
                    stratum="entry_year",
                    value=year,
                )
            fields = (
                [
                    "entry_origin",
                    "review_status",
                    "variant_type",
                    "gene_info",
                    "conditions",
                ]
                if layer == "vcv"
                else [
                    "entry_origin",
                    "review_status",
                    "collection_method",
                    "submitter",
                    "conditions",
                ]
            )
            for field in fields:
                limit = 10 if field in {"gene_info", "submitter"} else 6
                rows = c.execute(
                    f"SELECT {field},count(*) AS n FROM read_parquet(?) WHERE incident_rank=1 AND NOT identity_migrated AND {field} IS NOT NULL GROUP BY {field} HAVING count(DISTINCT variation_id)>=1000 ORDER BY n DESC LIMIT ?",
                    [str(path), limit],
                ).fetchall()
                for value, n in rows:
                    add(
                        layer,
                        f"{field}_{digest(value)[:8]}",
                        where=f"incident_rank=1 AND NOT identity_migrated AND {field}=?",
                        parameters=[value],
                        stratum=field,
                        value=value,
                    )
    completed = []
    hashes = {}
    for i, job in enumerate(jobs):
        source = Path(job["source"])
        if str(source) not in hashes:
            hashes[str(source)] = file_digest(source)
        identity = {
            "stage": "full-cohort-incidence",
            "format": 1,
            "input": hashes[str(source)],
            "implementation": implementation,
            "parameters": job,
            "replicates": args.replicates,
            "seed": args.seed,
        }
        output = destination / f"{job['layer']}_{job['name']}.parquet"
        curve_path = output.with_name(output.stem + "_curve.parquet")
        if checkpoint_ready(output, identity) and checkpoint_ready(
            curve_path, {**identity, "artifact": "curve"}
        ):
            print(f"Reusing {i + 1}/{len(jobs)} {output.stem}", flush=True)
        else:
            print(f"Estimating {i + 1}/{len(jobs)} {output.stem}", flush=True)
            data = load_cluster_counts(
                source,
                where=job["where"],
                parameters=job["parameters"],
                timing=job["timing"],
                removal_as_exit=job["removal"],
            )
            seed = args.seed + int(
                hashlib.sha256((job["layer"] + job["name"]).encode()).hexdigest()[:6],
                16,
            )
            rows, curve, _ = estimate_counts(
                data, [365, 730, 1095, 1825], replicates=args.replicates, seed=seed
            )
            for row in rows:
                row.update(
                    {
                        "layer": job["layer"],
                        "analysis": job["name"],
                        "stratum": job["stratum"],
                        "stratum_value": job["value"],
                        "timing": job["timing"],
                        "endpoint_policy": job["policy"],
                        "cadence_months": job["cadence"],
                        "removal_as_exit": job["removal"],
                    }
                )
            publish_checkpoint(pa.Table.from_pylist(rows), output, identity)
            publish_checkpoint(
                pa.table(curve), curve_path, {**identity, "artifact": "curve"}
            )
        completed.append(output)
        atomic_json(
            {
                "jobs": jobs,
                "completed": [str(p) for p in completed],
                "implementation": implementation,
            },
            destination / "progress.json",
        )
    table = pa.concat_tables(
        [pq.read_table(path) for path in completed], promote_options="permissive"
    )
    publish_checkpoint(
        table,
        args.root / "full_cohort_incidence.parquet",
        {
            "stage": "full-cohort-summary",
            "implementation": implementation,
            "jobs": {str(p): file_digest(p) for p in completed},
        },
    )
    print(f"Completed {len(jobs)} full-population analyses", flush=True)


if __name__ == "__main__":
    main()
