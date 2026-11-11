"""Bounded first-incident VUS competing-risk pilot with variant-cluster intervals."""

import argparse
from collections import Counter
from datetime import date
from pathlib import Path

import pyarrow as pa

from clinvar_history.checkpoints import (
    VERSIONS,
    digest,
    file_digest,
    implementation_identity,
    publish_checkpoint,
)
from clinvar_history.manifest import ReleaseKind, read_manifest
from clinvar_history.survival import (
    bootstrap_incidence,
    load_first_episodes,
    observations,
    snapshot_sensitivity,
)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--root",
        type=Path,
        help="VCV event root for actual monthly/quarterly sensitivity",
    )
    parser.add_argument("--episodes", type=Path, required=True)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--layer", choices=["vcv", "scv"], required=True)
    parser.add_argument("--cohort-start", type=date.fromisoformat, required=True)
    parser.add_argument("--cohort-end", type=date.fromisoformat, required=True)
    parser.add_argument("--study-end", type=date.fromisoformat, required=True)
    parser.add_argument("--max-clusters", type=int, default=5000)
    parser.add_argument("--replicates", type=int, default=200)
    parser.add_argument("--seed", type=int, default=20261111)
    parser.add_argument(
        "--stratum",
        choices=["review_status", "variant_type", "collection_method"],
        default="review_status",
    )
    parser.add_argument("--horizons", type=int, nargs="+", default=[365, 730, 1095])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    episodes = load_first_episodes(
        args.episodes,
        layer=args.layer,
        cohort_start=args.cohort_start,
        cohort_end=args.cohort_end,
        study_end=args.study_end,
        max_clusters=args.max_clusters,
        seed=args.seed,
        stratum=args.stratum,
    )
    if not episodes:
        raise SystemExit("No first incident episodes in the requested cohort")
    kind = ReleaseKind.VCF if args.layer == "vcv" else ReleaseKind.SUBMISSION_SUMMARY
    release_dates = [
        r.release_date
        for r in read_manifest(args.manifest)
        if r.kind is kind and r.release_date <= args.study_end
    ]
    rows = []
    for timing, removal in [
        ("observed", False),
        ("midpoint", False),
        ("observed", True),
    ]:
        data = observations(
            episodes,
            release_dates=release_dates,
            timing=timing,
            removal_as_exit=removal,
            study_end=args.study_end,
        )
        groups = {"all": data}
        if timing == "observed" and not removal:
            for stratum in sorted({r.stratum for r in data}):
                selected = [r for r in data if r.stratum == stratum]
                if len({r.cluster for r in selected}) >= 30:
                    groups[f"{args.stratum}={stratum}"] = selected
        for stratum, selected in groups.items():
            estimates = bootstrap_incidence(
                selected, args.horizons, replicates=args.replicates, seed=args.seed
            )
            counts = Counter(r.outcome for r in selected)
            for estimate in estimates:
                rows.append(
                    {
                        "layer": args.layer,
                        "timing": timing,
                        "removal_as_exit": removal,
                        "stratum": stratum,
                        "episodes": len(selected),
                        "clusters": len({r.cluster for r in selected}),
                        "censored": counts["censored"],
                        **estimate,
                    }
                )
    if args.root is not None:
        if args.layer != "vcv":
            raise SystemExit("Snapshot sensitivity currently supports the VCV layer")
        for months in (1, 3):
            selected = snapshot_sensitivity(
                args.root,
                episodes,
                release_dates=release_dates,
                study_end=args.study_end,
                cohort_start=args.cohort_start,
                cohort_end=args.cohort_end,
                months=months,
            )
            if not selected:
                print(f"No eligible episodes at {months}-month cadence")
                continue
            for estimate in bootstrap_incidence(
                selected, args.horizons, replicates=args.replicates, seed=args.seed
            ):
                rows.append(
                    {
                        "layer": args.layer,
                        "timing": f"snapshots_{months}_months",
                        "removal_as_exit": False,
                        "stratum": "all",
                        "episodes": len(selected),
                        "clusters": len({r.cluster for r in selected}),
                        "censored": sum(r.outcome == "censored" for r in selected),
                        **estimate,
                    }
                )
    identity = {
        "stage": "competing-risk-pilot",
        "versions": VERSIONS,
        "implementation": implementation_identity(),
        "episodes_sha256": file_digest(args.episodes),
        "manifest_sha256": file_digest(args.manifest),
        "snapshot_events_sha256": digest(
            {
                str(p.relative_to(args.root)): file_digest(p)
                for p in sorted((args.root / "events/vcv").rglob("*.parquet"))
            }
        )
        if args.root
        else None,
        "parameters": {
            k: str(v) if isinstance(v, (date, Path)) else v
            for k, v in vars(args).items()
            if k != "output"
        },
    }
    publish_checkpoint(pa.Table.from_pylist(rows), args.output, identity)
    print(
        f"Wrote {len(rows)} estimates from {len(episodes)} {args.layer} episodes to {args.output}"
    )


if __name__ == "__main__":
    main()
