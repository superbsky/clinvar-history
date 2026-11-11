"""Defined 100-event audit of the unexplained August 12 aggregate transitions."""

import argparse
from pathlib import Path

from clinvar_history.checkpoints import atomic_json, file_digest
from clinvar_history.research import _copy, research_connection
from clinvar_history.research_archives import extract_archived_records


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path("data/full"))
    p.add_argument(
        "--attributions",
        type=Path,
        default=Path("results/full-v1/attributions.parquet"),
    )
    p.add_argument("--archives", type=Path, default=Path("artifacts/research-xml"))
    p.add_argument(
        "--output", type=Path, default=Path("results/full-research/xml-historical")
    )
    p.add_argument("--seed", type=int, default=20261111)
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    with research_connection() as c:
        _copy(
            c,
            """WITH candidates AS (
          SELECT e.variation_id,e.observed_at,e.old_raw_classification,e.new_raw_classification,
           e.old_raw_classification||' -> '||e.new_raw_classification AS stratum
          FROM read_parquet(?) a JOIN read_parquet(?) e USING(variation_id,observed_at)
          WHERE a.primary_mechanism='unexplained' AND e.observed_at='2022-08-12'::DATE AND e.event_type='classification_changed' AND e.old_classification!=e.new_classification
        ), ranked AS (
          SELECT *,count(*) OVER(PARTITION BY stratum) AS population_in_stratum,
           row_number() OVER(PARTITION BY stratum ORDER BY sha256(CAST(variation_id AS VARCHAR)||?)) AS rank
          FROM candidates
        ), selected AS (
          SELECT * FROM ranked ORDER BY rank,sha256(stratum||?) LIMIT 100
        ) SELECT *,count(*) OVER(PARTITION BY stratum) AS sampled_in_stratum FROM selected""",
            args.output / "sample.parquet",
            [
                str(args.attributions),
                str(args.root / "events/vcv/release_date=*/events.parquet"),
                str(args.seed),
                str(args.seed),
            ],
        )
        ids = [
            r[0]
            for r in c.execute(
                "SELECT DISTINCT variation_id FROM read_parquet(?)",
                [str(args.output / "sample.parquet")],
            ).fetchall()
        ]
    atomic_json(
        {
            "design": "Round-robin allocation across raw direction strata, SHA256 seeded selection within strata; one dated event per selected variant. Population is unexplained canonical transitions observed on 2022-08-12, not all unexplained history.",
            "seed": args.seed,
            "sample_sha256": file_digest(args.output / "sample.parquet"),
            "requested_variants": len(ids),
        },
        args.output / "design.json",
    )
    for month in ["08", "09"]:
        path = args.archives / f"ClinVarVariationRelease_2022-{month}.xml.gz"
        print(f"Extracting {path.name}", flush=True)
        extract_archived_records(path, ids, args.output / f"2022-{month}")


if __name__ == "__main__":
    main()
