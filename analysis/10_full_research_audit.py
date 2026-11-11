"""Independently reconstruct full VUS cohorts and audit research assumptions."""

import argparse
from datetime import date
from itertools import pairwise
from pathlib import Path

from clinvar_history.checkpoints import (
    atomic_json,
    dataset_snapshot,
    implementation_identity,
)
from clinvar_history.manifest import ReleaseKind, read_manifest
from clinvar_history.research import (
    _copy,
    compare_legacy,
    endpoint_sql,
    reconstruct_episodes,
    research_connection,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path("data/full"))
    parser.add_argument("--legacy-results", type=Path, default=Path("results/full-v1"))
    parser.add_argument("--output", type=Path, default=Path("results/full-research"))
    parser.add_argument(
        "--legacy-end", type=date.fromisoformat, default=date(2026, 7, 24)
    )
    parser.add_argument(
        "--baseline-root",
        type=Path,
        help="Original retained root for exact legacy comparison",
    )
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    releases = read_manifest(args.root / "manifest/releases.parquet")
    snapshot = dataset_snapshot(args.root, releases)
    summary = {
        "input": snapshot,
        "implementation": implementation_identity(),
        "layers": {},
        "legacy_end": str(args.legacy_end),
    }
    with research_connection(args.output / ".scratch") as c:
        for layer, kind, legacy_name in [
            ("vcv", ReleaseKind.VCF, "vus_time_to_resolution"),
            ("scv", ReleaseKind.SUBMISSION_SUMMARY, "scv_vus_time_to_resolution"),
        ]:
            days = sorted(r.release_date for r in releases if r.kind is kind)
            gaps = [(b - a).days for a, b in pairwise(days)]
            report = {
                "source_first": str(days[0]),
                "source_last": str(days[-1]),
                "source_releases": len(days),
                "gap_min_days": min(gaps),
                "gap_max_days": max(gaps),
                "unobserved_days_to_legacy_end": (args.legacy_end - days[-1]).days,
                "reconstructions": {},
            }
            print(f"{layer}: independently checking every legacy episode", flush=True)
            legacy_path = args.output / f"legacy_check_{layer}.parquet"
            reconstruct_episodes(
                args.baseline_root or args.root,
                legacy_path,
                layer=layer,
                policy="strict",
                administrative_end=args.legacy_end,
            )
            comparison = compare_legacy(
                c, legacy_path, args.legacy_results / f"{legacy_name}.parquet", layer
            )
            report["legacy_episode_comparison"] = comparison
            if not comparison["exact"]:
                atomic_json(summary, args.output / "audit_incomplete.json")
                raise ValueError(
                    f"Independent episode reconstruction disagrees: {layer} {comparison}"
                )
            for policy in ["strict", "family"]:
                for cadence in [0, 1, 3]:
                    path = args.output / f"episodes_{layer}_{policy}_{cadence}.parquet"
                    print(
                        f"Reconstructing {layer} {policy} cadence={cadence}", flush=True
                    )
                    episodes, units = reconstruct_episodes(
                        args.root,
                        path,
                        layer=layer,
                        policy=policy,
                        cadence_months=cadence,
                    )
                    report["reconstructions"][f"{policy}_{cadence}"] = {
                        "episodes": episodes,
                        "first_incident_units": units,
                        "path": str(path),
                    }
                    _copy(
                        c,
                        "SELECT year(started_at) AS entry_year, outcome, censor_reason, identity_migrated, count(*) AS episodes, count(*) FILTER (WHERE incident_rank=1) AS first_incident_units FROM read_parquet(?) GROUP BY ALL ORDER BY entry_year,outcome",
                        args.output
                        / f"cohort_inventory_{layer}_{policy}_{cadence}.parquet",
                        [str(path)],
                    )
            event_type = (
                "classification_changed" if layer == "vcv" else "scv_reclassified"
            )
            source = str(args.root / f"events/{layer}/release_date=*/events.parquet")
            old = endpoint_sql("old_raw_classification", "old_classification", "family")
            new = endpoint_sql("new_raw_classification", "new_classification", "family")
            _copy(
                c,
                f"""WITH t AS (SELECT *,{old} AS old_family,{new} AS new_family FROM read_parquet(?) WHERE event_type=?)
                SELECT observed_at,old_classification,new_classification,old_raw_classification,new_raw_classification,
                       old_family,new_family,count(*) AS events,
                       CASE WHEN old_family=new_family THEN 'within_family_or_synonym'
                            WHEN old_family IN ('benign','pathogenic') AND new_family IN ('benign','pathogenic') THEN 'benign_pathogenic_switch'
                            WHEN old_family IN ('pathogenic','vus') AND new_family IN ('pathogenic','vus') THEN 'pathogenic_vus_switch'
                            WHEN old_family IN ('benign','vus') AND new_family IN ('benign','vus') THEN 'benign_vus_switch'
                            WHEN old_family='conflict' OR new_family='conflict' THEN 'conflict_entry_exit'
                            ELSE 'other_or_missing' END AS transition_group
                FROM t GROUP BY ALL ORDER BY observed_at""",
                args.output / f"raw_transition_audit_{layer}.parquet",
                [source, event_type],
            )
            key = "variation_id" if layer == "vcv" else "scv_accession"
            report["identity_migrated_units"] = c.execute(
                f"SELECT count(*) FROM (SELECT {key} FROM read_parquet(?) GROUP BY {key} HAVING count(DISTINCT variation_id)>1)",
                [source],
            ).fetchone()[0]
            summary["layers"][layer] = report
            atomic_json(summary, args.output / "audit.json")
    if date(2022, 8, 12) in {
        r.release_date for r in releases if r.kind is ReleaseKind.VCF
    }:
        reconstruct_episodes(
            args.root,
            args.output / "episodes_vcv_family_no20220812.parquet",
            layer="vcv",
            excluded_dates=(date(2022, 8, 12),),
        )
    print("Full independent reconstruction audit complete", flush=True)


if __name__ == "__main__":
    main()
