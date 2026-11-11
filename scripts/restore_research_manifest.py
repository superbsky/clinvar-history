"""Restore a pinned public-source CSV manifest for a new, empty data root."""

import argparse
import csv
from datetime import date
from pathlib import Path

from clinvar_history.checkpoints import stamp_checkpoint
from clinvar_history.manifest import ReleaseFile, ReleaseKind, write_manifest
from clinvar_history.runner import manifest_request


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--csv", type=Path, default=Path("docs/research-tables/baseline-manifest.csv")
    )
    p.add_argument(
        "--output", type=Path, default=Path("data/full/manifest/releases.parquet")
    )
    args = p.parse_args()
    if args.output.exists():
        raise ValueError("Refuse to replace an existing manifest")
    with args.csv.open(newline="") as stream:
        rows = list(csv.DictReader(stream))
    releases = [
        ReleaseFile(
            release_date=date.fromisoformat(r["release_date"]),
            kind=ReleaseKind(r["kind"]),
            url=r["url"],
            size_bytes=int(r["size_bytes"]),
            md5=r["md5"] or None,
            source_month=r["source_month"] or None,
        )
        for r in rows
    ]
    write_manifest(releases, args.output)
    stamp_checkpoint(args.output, manifest_request(date(2020, 1, 1), date(2026, 7, 24)))
    print(
        f"Restored {len(releases)} pinned sources; source bytes must still pass recorded checksums"
    )


if __name__ == "__main__":
    main()
