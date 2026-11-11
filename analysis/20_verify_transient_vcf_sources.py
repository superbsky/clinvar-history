"""Verify the 117-variant July 2023 round trip against three complete source streams."""

import argparse
import csv
from datetime import UTC, date, datetime
from pathlib import Path

import httpx
import pyarrow as pa
from pyarrow import csv as pcsv

from clinvar_history.checkpoints import atomic_json, file_digest
from clinvar_history.manifest import ReleaseKind, read_manifest
from clinvar_history.pipeline import GzipLineStream, fetch_release_md5
from clinvar_history.vcf_parser import parse_info


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--offline", action="store_true")
    args = p.parse_args()
    root = Path("docs/research-tables/liang-crosscheck")
    cache = Path("artifacts/liang-crosscheck/source-vcf")
    cache.mkdir(parents=True, exist_ok=True)
    histories = list(
        csv.DictReader((root / "native_only_primary_policy_check.csv").open())
    )
    selected = [
        r for r in histories if r["observed_at"] in ["2023-07-10", "2023-07-17"]
    ]
    ids = {int(r["variation_id"]) for r in selected}
    assert len(ids) == 117 and len(selected) == 234
    releases = [
        r
        for r in read_manifest(Path("data/research/manifest/releases.parquet"))
        if r.kind is ReleaseKind.VCF
        and r.release_date in {date(2023, 7, 2), date(2023, 7, 10), date(2023, 7, 17)}
    ]
    assert len(releases) == 3
    all_rows = []
    receipts = []
    total = 0
    with httpx.Client(timeout=90, follow_redirects=True) as client:
        for release in sorted(releases, key=lambda r: r.release_date):
            output = cache / (str(release.release_date) + ".json")
            if args.offline:
                import json

                payload = json.loads(output.read_text())
                if (
                    file_digest(cache / (str(release.release_date) + ".tsv"))
                    != payload["raw_extract_sha256"]
                ):
                    raise ValueError("Raw source-extract cache changed")
                if payload["source_md5"] != release.md5:
                    raise ValueError("Source identity changed")
                receipts.append(payload)
                all_rows.extend(payload["selected_rows"])
                continue
            official = fetch_release_md5(release, client=client)
            if official != release.md5:
                raise ValueError(
                    "Official source checksum changed since the pinned study"
                )
            print(f"Verifying full source {release.release_date}", flush=True)
            count = 0
            selected_rows = []
            raw_rows = []
            header_date = None
            source_bytes = 0
            with client.stream("GET", release.url) as response:
                response.raise_for_status()

                def chunks():
                    nonlocal source_bytes, total
                    for chunk in response.iter_bytes(chunk_size=65536):
                        source_bytes += len(chunk)
                        total += len(chunk)
                        if source_bytes > 200_000_000 or total > 600_000_000:
                            raise ValueError("Source transfer budget exceeded")
                        yield chunk

                stream = GzipLineStream(chunks(), expected_md5=release.md5)
                for line in stream:
                    if line.startswith(b"#"):
                        if line.startswith(b"##fileDate="):
                            header_date = line.decode().strip().split("=", 1)[1]
                        continue
                    count += 1
                    early = line.split(b"\t", 3)
                    if int(early[2]) not in ids:
                        continue
                    text = line.decode().rstrip("\r\n")
                    columns = text.split("\t")
                    info = parse_info(columns[7])
                    raw_rows.append(text)
                    selected_rows.append(
                        {
                            "release_date": str(release.release_date),
                            "variation_id": int(columns[2]),
                            "raw_classification": info.get("CLNSIG"),
                            "review_status": info.get("CLNREVSTAT"),
                            "gene_info": info.get("GENEINFO"),
                        }
                    )
                if stream.calculated_md5 != release.md5:
                    raise ValueError("Source was not fully consumed/verified")
                modified = response.headers.get("Last-Modified")
            if (
                len(selected_rows) != 117
                or len({r["variation_id"] for r in selected_rows}) != 117
            ):
                raise ValueError("Missing or duplicated target variant in source")
            extract = cache / (str(release.release_date) + ".tsv")
            extract.write_text("\n".join(raw_rows) + "\n")
            payload = {
                "url": release.url,
                "source_md5": stream.calculated_md5,
                "official_md5": official,
                "compressed_bytes": source_bytes,
                "all_data_rows_scanned": count,
                "selected_records": len(selected_rows),
                "header_file_date": header_date,
                "http_last_modified": modified,
                "retrieved_at": datetime.now(UTC).isoformat(),
                "raw_extract_sha256": file_digest(extract),
                "selected_rows": selected_rows,
            }
            atomic_json(payload, output)
            receipts.append(payload)
            all_rows.extend(selected_rows)
    actual = {
        (r["release_date"], r["variation_id"]): r["raw_classification"]
        for r in all_rows
    }
    for event in selected:
        vid = int(event["variation_id"])
        assert actual[event["observed_at"], vid] == event["new_raw_classification"]
        if event["observed_at"] == "2023-07-10":
            assert actual["2023-07-02", vid] == event["old_raw_classification"]
    for vid in ids:
        assert actual["2023-07-02", vid] == actual["2023-07-17", vid]
    pcsv.write_csv(
        pa.Table.from_pylist(all_rows), root / "verified_source_vcf_rows.csv"
    )
    atomic_json(
        {
            "variants": len(ids),
            "implementation_sha256": file_digest(Path(__file__)),
            "target_table_sha256": file_digest(
                root / "native_only_primary_policy_check.csv"
            ),
            "manifest_sha256": file_digest(
                Path("data/research/manifest/releases.parquet")
            ),
            "selected_source_rows": len(all_rows),
            "raw_label_round_trips_verified": len(ids),
            "compressed_source_bytes": sum(r["compressed_bytes"] for r in receipts),
            "verification": "Every compressed source fully consumed and matched both the pinned and freshly fetched official MD5; all target IDs unique; source labels match each saved event and return to pre-window label.",
            "source_receipts": [
                {k: v for k, v in r.items() if k != "selected_rows"} for r in receipts
            ],
            "raw_archives_retained": 0,
            "mechanism": "Not established by VCF labels; clinical reinterpretation is not inferred.",
        },
        root / "source_vcf_verification.json",
    )
    print(
        f"Confirmed {len(ids)} exact source-label round trips in {len(all_rows)} selected source rows",
        flush=True,
    )


if __name__ == "__main__":
    main()
