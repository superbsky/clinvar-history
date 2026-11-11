"""Compare dated XML records and retrieve a bounded intervening-version series."""

import argparse
import json
from collections import Counter
from pathlib import Path

import httpx
import pyarrow as pa
import pyarrow.parquet as pq

from clinvar_history.checkpoints import atomic_json, file_digest
from clinvar_history.xml_review import compare_versions


def assertion_signature(record, fields):
    return sorted(
        tuple(
            a.get(f) if not isinstance(a.get(f), list) else tuple(a[f]) for f in fields
        )
        for a in record["assertions"]
    )


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--fetch-versions", action="store_true")
    p.add_argument("--offline", action="store_true")
    args = p.parse_args()
    root = Path("results/full-research")
    hist = root / "xml-historical"
    design = pq.read_table(hist / "sample.parquet").to_pylist()
    before = json.loads((hist / "2022-08/extraction.json").read_text())
    after = json.loads((hist / "2022-09/extraction.json").read_text())
    old = {r["variation_id"]: r for r in before["records"]}
    new = {r["variation_id"]: r for r in after["records"]}
    comparisons = []
    for event in design:
        vid = event["variation_id"]
        a, b = old.get(vid), new.get(vid)
        row = dict(event)
        row.update(
            {"before_available": a is not None, "after_available": b is not None}
        )
        if a and b and "error" not in a and "error" not in b:
            row.update(
                {
                    "before_version": int(a["version"]),
                    "after_version": int(b["version"]),
                    "aggregate_unchanged": a.get("aggregate_classifications")
                    == b.get("aggregate_classifications"),
                    "before_labels": a.get("aggregate_classifications", []),
                    "after_labels": b.get("aggregate_classifications", []),
                    "assertion_labels_versions_unchanged": assertion_signature(
                        a, ["accession", "version", "classification"]
                    )
                    == assertion_signature(
                        b, ["accession", "version", "classification"]
                    ),
                    "assertion_metadata_unchanged": assertion_signature(
                        a,
                        [
                            "accession",
                            "version",
                            "classification",
                            "review_status",
                            "conditions",
                            "contributes",
                        ],
                    )
                    == assertion_signature(
                        b,
                        [
                            "accession",
                            "version",
                            "classification",
                            "review_status",
                            "conditions",
                            "contributes",
                        ],
                    ),
                    "aggregate_review_unchanged": a.get("review_status")
                    == b.get("review_status"),
                    "before_aggregate_matches_vcf_pre_window": ";".join(
                        a.get("aggregate_classifications", [])
                    )
                    .strip()
                    .lower()
                    == event["old_raw_classification"].strip().lower(),
                    "both_interpreted": a["record_type"] == "interpreted"
                    and b["record_type"] == "interpreted",
                }
            )
        comparisons.append(row)
    pq.write_table(
        pa.Table.from_pylist(comparisons),
        hist / "comparisons.parquet",
        compression="zstd",
    )
    counts = {
        field: sum(bool(r.get(field)) for r in comparisons)
        for field in [
            "before_available",
            "after_available",
            "aggregate_unchanged",
            "assertion_labels_versions_unchanged",
            "assertion_metadata_unchanged",
            "aggregate_review_unchanged",
            "before_aggregate_matches_vcf_pre_window",
            "both_interpreted",
        ]
    }
    counts.update(
        {
            "selected_events": len(comparisons),
            "raw_direction_strata": len({r["stratum"] for r in comparisons}),
            "vcv_version_deltas": dict(
                Counter(
                    str(r["after_version"] - r["before_version"])
                    for r in comparisons
                    if "after_version" in r
                )
            ),
            "archive_vintage_warning": "Both dated source files were last modified in March 2026; root ReleaseDate is 2022-08-01/2022-09-01. They are the currently available dated archives, not a verified immutable 2022 byte capture.",
            "historical_mechanism": "Endpoint consistency and stable SCV evidence do not identify the intervening process or prove a clinical cause.",
        }
    )
    atomic_json(counts, hist / "summary.json")
    current = json.loads((root / "xml-current/review.json").read_text())
    current_counts = {
        field: sum(
            bool(e["review"].get("extra_fields", {}).get(field))
            for e in current["events"]
        )
        for field in [
            "aggregate_explanation",
            "assertion_methods",
            "contribution_flags",
            "submission_comments",
            "submission_citations",
        ]
    }
    current_counts.update(
        {
            "events": len(current["events"]),
            "statuses": current["status_counts"],
            "bytes": current["total_xml_bytes"],
            "records_with_noncontributing_scv": sum(
                any(
                    a.get("contributes") == "false"
                    for a in e["review"].get("assertions", [])
                )
                for e in current["events"]
            ),
            "historically_resolved_from_current_xml_only": 0,
        }
    )
    atomic_json(current_counts, root / "xml-current/summary.json")
    print({"historical": counts, "current": current_counts}, flush=True)
    if not args.fetch_versions:
        return
    # Ten largest raw-direction strata; deterministic first sampled case with an
    # exactly two-step VCV version increase. This is targeted, not random validation.
    selected = {}
    for r in sorted(
        comparisons,
        key=lambda r: (-r["population_in_stratum"], r["rank"], r["variation_id"]),
    ):
        if (
            r.get("after_version", 0) - r.get("before_version", 0) == 2
            and r["stratum"] not in selected
        ):
            selected[r["stratum"]] = r
        if len(selected) == 10:
            break
    series = []
    for event in selected.values():
        vid = event["variation_id"]
        version = event["before_version"]
        print(
            f"Retrieving targeted version series {vid}: {version}/{version + 1}/{version + 2}",
            flush=True,
        )
        row = {
            "event": event,
            "selection": "First seeded sampled case with VCV delta=2 in each of the ten largest raw-direction strata",
            "max_record_bytes": 2000000,
            "max_total_requested_records": 30,
        }
        try:
            first = compare_versions(
                vid,
                version,
                version + 1,
                root / f"xml-versions/{vid}",
                offline=args.offline,
            )
            second = compare_versions(
                vid,
                version + 1,
                version + 2,
                root / f"xml-versions/{vid}",
                offline=args.offline,
            )
            row.update(
                {
                    "before": first["old"],
                    "middle": first["new"],
                    "after": second["new"],
                    "first_scv_changes": first["scv_changes"],
                    "second_scv_changes": second["scv_changes"],
                    "status": "retrieved",
                }
            )
        except (httpx.HTTPError, ValueError, OSError) as exc:
            row.update(
                {"status": "unavailable", "error": f"{type(exc).__name__}: {exc}"}
            )
        series.append(row)
        atomic_json(
            {
                "cases": series,
                "sample_sha256": file_digest(hist / "sample.parquet"),
                "purpose": "Assess whether explicitly versioned records retain the intervening aggregate change; version retrieval may include later metadata edits.",
            },
            root / "xml-versions/series.json",
        )


if __name__ == "__main__":
    main()
