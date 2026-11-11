"""Clone retained data safely and append newly discoverable releases in the same horizon."""

import argparse
import os
import shutil
from datetime import date
from pathlib import Path

import httpx

from clinvar_history.checkpoints import atomic_json, stamp_checkpoint
from clinvar_history.manifest import (
    ReleaseKind,
    default_index_urls,
    discover_manifest,
    read_manifest,
    write_manifest,
)
from clinvar_history.pipeline import (
    http_vcf_source,
    ingest_vcf_releases,
    resolve_release_checksums,
)
from clinvar_history.runner import manifest_request


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source", type=Path, default=Path("data/full"))
    parser.add_argument("--target", type=Path, default=Path("data/research"))
    parser.add_argument("--end", type=date.fromisoformat, default=date(2026, 7, 24))
    args = parser.parse_args()
    if args.source.resolve() == args.target.resolve():
        raise ValueError("Research clone must use a separate root")
    if not args.target.exists():
        args.target.mkdir(parents=True)
        for path in args.source.rglob("*"):
            rel = path.relative_to(args.source)
            if path.is_dir():
                (args.target / rel).mkdir(exist_ok=True)
            elif path.suffix == ".parquet":
                os.link(path, args.target / rel)
            else:
                shutil.copy2(path, args.target / rel)
    manifest = args.target / "manifest/releases.parquet"
    retained = read_manifest(manifest)
    start = min(r.release_date for r in retained)
    with httpx.Client(timeout=60, follow_redirects=True) as client:
        discovered = discover_manifest(
            default_index_urls(start, args.end),
            start=start,
            end=args.end,
            client=client,
        )
        known = {r.url for r in retained}
        new = [r for r in discovered if r.url not in known]
        for release in new:
            previous = max(r.release_date for r in retained if r.kind is release.kind)
            if release.release_date <= previous:
                raise ValueError(
                    "New source inserts historical releases; requires an independent backfill"
                )
            if release.kind is not ReleaseKind.VCF:
                raise ValueError(
                    "Unexpected new SCV release; review the source horizon"
                )
        new = resolve_release_checksums(new, client=client)
        releases = sorted(
            [*retained, *new], key=lambda r: (r.release_date, r.kind.value, r.url)
        )
        write_manifest(releases, manifest)
        stamp_checkpoint(manifest, manifest_request(start, args.end))

        def source(release):
            return http_vcf_source(release, client=client, expected_md5=release.md5)

        summary = ingest_vcf_releases(releases, root=args.target, source=source)
    report = {
        "retained_releases": len(retained),
        "research_releases": len(releases),
        "new_urls": [r.url for r in new],
        "new_source_bytes": sum(r.size_bytes for r in new),
        "processed": summary.processed,
        "skipped": summary.skipped,
        "data_root": str(args.target),
        "last_vcf": str(
            max(r.release_date for r in releases if r.kind is ReleaseKind.VCF)
        ),
    }
    atomic_json(report, Path("results/full-research/source_update.json"))
    print(report, flush=True)


if __name__ == "__main__":
    main()
