"""Bounded current-XML evidence review; never infer historical causes."""

import hashlib
import json
import time
import xml.etree.ElementTree as ET
from collections import Counter
from datetime import UTC, datetime
from pathlib import Path

import duckdb
import httpx

from .checkpoints import atomic_json, digest, file_digest

EFETCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"


def parse_vcv_xml(content: bytes, variation_id: int) -> dict:
    if b"<!DOCTYPE" in content.upper() or b"<!ENTITY" in content.upper():
        raise ValueError("DTD/entity XML is unsupported")
    root = ET.fromstring(content)
    for node in root.iter():
        node.tag = node.tag.rsplit("}", 1)[-1]
    records = [
        r
        for r in root.iter("VariationArchive")
        if r.get("VariationID") == str(variation_id)
    ]
    if len(records) != 1:
        raise ValueError("XML does not contain exactly the requested VariationID")
    record = records[0]
    gc = record.find("./ClassifiedRecord/Classifications/GermlineClassification")
    assertions = []
    for assertion in record.findall(
        "./ClassifiedRecord/ClinicalAssertionList/ClinicalAssertion"
    ):
        accession = assertion.find("ClinVarAccession")
        if accession is None:
            continue
        classification = assertion.find("Classification")
        if classification is None:
            continue
        germline = classification.findtext("GermlineClassification")
        assertions.append(
            {
                "accession": accession.get("Accession"),
                "version": accession.get("Version"),
                "date_created": accession.get("DateCreated"),
                "date_updated": accession.get("DateUpdated"),
                "submission_date": assertion.get("SubmissionDate"),
                "submitter": accession.get("SubmitterName"),
                "classification": germline,
                "classification_type": "germline" if germline is not None else "other",
                "review_status": classification.findtext("ReviewStatus"),
                "last_evaluated": classification.get("DateLastEvaluated"),
                "contributes": assertion.get("ContributesToAggregateClassification"),
                "comments": [n.text for n in assertion.iter("Comment") if n.text],
                "assertion_methods": [
                    n.text
                    for n in assertion.iter("Attribute")
                    if n.get("Type") == "AssertionMethod" and n.text
                ],
                "pubmed_ids": sorted(
                    {
                        n.text
                        for n in assertion.iter("ID")
                        if n.get("Source") == "PubMed" and n.text
                    }
                ),
            }
        )
    return {
        "variation_id": variation_id,
        "accession": record.get("Accession"),
        "version": record.get("Version"),
        "date_last_updated": record.get("DateLastUpdated"),
        "most_recent_submission": record.get("MostRecentSubmission"),
        "record_status": record.findtext("RecordStatus"),
        "current_germline_classification": gc.findtext("Description")
        if gc is not None
        else None,
        "aggregate_explanations": [n.text for n in gc.findall("Explanation") if n.text]
        if gc is not None
        else [],
        "assertions": assertions,
    }


def sample_unexplained(
    root: Path, attributions: Path, *, limit=12, seed=20261111
) -> list[dict]:
    if not 1 <= limit <= 100:
        raise ValueError("XML sample size must be between 1 and 100")
    with duckdb.connect() as connection:
        query = """
            WITH candidates AS (
                SELECT a.variation_id, a.observed_at, e.old_classification, e.new_classification,
                       CAST(year(a.observed_at) AS VARCHAR) || ':' || e.old_classification || ':' || e.new_classification AS stratum
                FROM read_parquet(?) a JOIN read_parquet(?) e USING (variation_id, observed_at)
                WHERE a.primary_mechanism = 'unexplained' AND e.event_type = 'classification_changed'
                  AND e.old_classification != e.new_classification
            ), ranked AS (
                SELECT *, count(*) OVER (PARTITION BY stratum) AS population_in_stratum,
                       row_number() OVER (PARTITION BY stratum ORDER BY sha256(CAST(variation_id AS VARCHAR) || CAST(observed_at AS VARCHAR) || ?)) AS rank_in_stratum
                FROM candidates
            )
            SELECT * FROM ranked ORDER BY rank_in_stratum, sha256(stratum || ?) LIMIT ?
        """
        cursor = connection.execute(
            query,
            [
                str(attributions),
                str(root / "events/vcv/release_date=*/events.parquet"),
                str(seed),
                str(seed),
                limit,
            ],
        )
        names = [r[0] for r in cursor.description]
        return [dict(zip(names, row, strict=True)) for row in cursor.fetchall()]


def fetch_xml(
    variation_id: int,
    *,
    client: httpx.Client,
    max_bytes: int,
    version: int | None = None,
) -> tuple[bytes, str]:
    content = bytearray()
    try:
        with client.stream(
            "GET",
            EFETCH,
            params={
                "db": "clinvar",
                "rettype": "vcv",
                **({"is_variationid": "true"} if version is None else {}),
                "id": str(variation_id)
                if version is None
                else f"VCV{variation_id:09d}.{version}",
                "tool": "clinvar-history",
            },
        ) as response:
            response.raise_for_status()
            deadline = time.monotonic() + 30
            for chunk in response.iter_bytes(chunk_size=min(65536, max_bytes + 1)):
                if time.monotonic() > deadline:
                    raise ValueError("XML response exceeds time budget")
                if len(content) + len(chunk) > max_bytes:
                    raise ValueError("XML response exceeds byte budget")
                content.extend(chunk)
            url = str(response.url)
    finally:
        time.sleep(0.4)  # throttle successful and failed HTTP requests alike
    return bytes(content), url


def review_sample(
    root: Path,
    attributions: Path,
    output: Path,
    *,
    limit=12,
    seed=20261111,
    max_record_bytes=2_000_000,
    max_total_bytes=20_000_000,
    offline=False,
) -> dict:
    if max_record_bytes < 1 or max_total_bytes < 1:
        raise ValueError("XML byte budgets must be positive")
    sample = sample_unexplained(root, attributions, limit=limit, seed=seed)
    if not sample:
        raise ValueError("No unexplained transitions to review")
    output.mkdir(parents=True, exist_ok=True)
    ids = sorted({r["variation_id"] for r in sample})
    with duckdb.connect() as connection:
        # Preserved local history is compared with current accessions, not treated as
        # a historical XML snapshot. Read the bounded set of variants in one pass.
        local = connection.execute(
            "SELECT DISTINCT variation_id, scv_accession FROM read_parquet(?) WHERE variation_id IN (SELECT unnest(?))",
            [str(root / "events/scv/release_date=*/events.parquet"), ids],
        ).fetchall()
    local_accessions = {vid: {scv for v, scv in local if v == vid} for vid in ids}
    reviewed = {}
    total_bytes = 0
    charged_bytes = 0
    with httpx.Client(timeout=30, follow_redirects=True) as client:
        for vid in ids:
            cache = output / f"{vid}.xml"
            provenance = output / f"{vid}.json"
            try:
                if cache.exists() and provenance.exists():
                    saved = json.loads(provenance.read_text())
                    if (
                        cache.stat().st_size > max_record_bytes
                        or charged_bytes + cache.stat().st_size > max_total_bytes
                        or file_digest(cache) != saved["sha256"]
                    ):
                        raise ValueError("Cached XML fails integrity or size check")
                    content = cache.read_bytes()
                    charged_bytes += len(content)
                else:
                    if offline:
                        raise ValueError("No cached XML in offline mode")
                    remaining = min(max_record_bytes, max_total_bytes - charged_bytes)
                    if remaining <= 0:
                        raise ValueError("Total XML byte budget exhausted")
                    charged_bytes += remaining  # reserve failed requests conservatively
                    content, url = fetch_xml(vid, client=client, max_bytes=remaining)
                    charged_bytes -= remaining - len(content)
                    parse_vcv_xml(content, vid)
                    cache.write_bytes(content)
                    saved = {
                        "url": url,
                        "fetched_at": datetime.now(UTC).isoformat(),
                        "sha256": hashlib.sha256(content).hexdigest(),
                        "bytes": len(content),
                    }
                    atomic_json(saved, provenance)
                total_bytes += len(content)
                if total_bytes > max_total_bytes:
                    raise ValueError("Total XML byte budget exhausted")
                parsed = parse_vcv_xml(content, vid)
                current = {a["accession"] for a in parsed["assertions"]}
                extras = {
                    "aggregate_explanation": bool(parsed["aggregate_explanations"]),
                    "submission_comments": any(
                        a["comments"] for a in parsed["assertions"]
                    ),
                    "assertion_methods": any(
                        a["assertion_methods"] for a in parsed["assertions"]
                    ),
                    "submission_citations": any(
                        a["pubmed_ids"] for a in parsed["assertions"]
                    ),
                    "contribution_flags": any(
                        a["contributes"] is not None for a in parsed["assertions"]
                    ),
                }
                reviewed[vid] = {
                    "status": "reviewed_current_record",
                    "provenance": saved,
                    "extra_fields": extras,
                    "local_scv_accessions": len(local_accessions[vid]),
                    "local_accessions_absent_from_current_xml": sorted(
                        local_accessions[vid] - current
                    ),
                    "historical_mechanism": "unresolved_current_xml_only",
                    **parsed,
                }
            except (
                httpx.HTTPError,
                ValueError,
                ET.ParseError,
                OSError,
                KeyError,
                TypeError,
            ) as error:
                reviewed[vid] = {
                    "status": "unavailable",
                    "error": str(error),
                    "historical_mechanism": "unresolved_current_xml_only",
                }
    result = {
        "format": 1,
        "sample_seed": seed,
        "sample_limit": limit,
        "attributions_sha256": file_digest(attributions),
        "sampling": "deterministic round robin across year and transition direction strata",
        "total_xml_bytes": total_bytes,
        "budget_charged_bytes": charged_bytes,
        "local_scv_events_sha256": digest(
            {
                str(p.relative_to(root)): file_digest(p)
                for p in sorted((root / "events/scv").rglob("*.parquet"))
            }
        ),
        "budgets": {"record": max_record_bytes, "total": max_total_bytes},
        "events": [
            {**event, "review": reviewed[event["variation_id"]]} for event in sample
        ],
        "status_counts": dict(
            Counter(reviewed[e["variation_id"]]["status"] for e in sample)
        ),
    }
    result["sample_identity"] = digest(sample)
    atomic_json(result, output / "review.json")
    return result


def compare_versions(
    variation_id: int,
    old_version: int,
    new_version: int,
    output: Path,
    *,
    offline=False,
    max_record_bytes=2_000_000,
):
    """Retrieve two explicitly versioned VCVs; preserve uncertainty about dates."""
    if not 0 < old_version < new_version or variation_id < 1 or max_record_bytes < 1:
        raise ValueError("Invalid variation, version pair or byte budget")
    output.mkdir(parents=True, exist_ok=True)
    records = []
    with httpx.Client(timeout=30, follow_redirects=True) as client:
        for version in (old_version, new_version):
            cache = output / f"VCV{variation_id:09d}.{version}.xml"
            meta = cache.with_suffix(".json")
            if cache.exists() and meta.exists():
                provenance = json.loads(meta.read_text())
                if (
                    cache.stat().st_size > max_record_bytes
                    or file_digest(cache) != provenance["sha256"]
                ):
                    raise ValueError("Versioned XML cache is corrupt or exceeds budget")
                content = cache.read_bytes()
            else:
                if offline:
                    raise ValueError("Versioned XML cache missing in offline mode")
                content, url = fetch_xml(
                    variation_id,
                    client=client,
                    max_bytes=max_record_bytes,
                    version=version,
                )
                parsed = parse_vcv_xml(content, variation_id)
                if parsed["version"] != str(version):
                    raise ValueError("NCBI returned a different VCV version")
                cache.write_bytes(content)
                provenance = {
                    "url": url,
                    "fetched_at": datetime.now(UTC).isoformat(),
                    "sha256": hashlib.sha256(content).hexdigest(),
                    "bytes": len(content),
                }
                atomic_json(provenance, meta)
            parsed = parse_vcv_xml(content, variation_id)
            if parsed["version"] != str(version):
                raise ValueError(
                    "Cached VCV version does not match the requested version"
                )
            records.append({"provenance": provenance, **parsed})
    before, after = ({a["accession"]: a for a in r["assertions"]} for r in records)
    changes = []
    for accession in sorted(before.keys() | after.keys()):
        old, new = before.get(accession), after.get(accession)
        changed_fields = [
            k
            for k in ("classification", "review_status", "contributes", "version")
            if old and new and old[k] != new[k]
        ]
        if old is None or new is None or changed_fields:
            changes.append(
                {
                    "accession": accession,
                    "change": "added"
                    if old is None
                    else "removed"
                    if new is None
                    else "updated",
                    "changed_fields": changed_fields,
                    "before": old,
                    "after": new,
                }
            )
    comparison = {
        "variation_id": variation_id,
        "old": records[0],
        "new": records[1],
        "scv_changes": changes,
        "interpretation": "Differences between explicitly requested VCV versions; not immutable release snapshots or proven historical causes",
        "max_total_bytes": 2 * max_record_bytes,
    }
    atomic_json(comparison, output / "comparison.json")
    return comparison
