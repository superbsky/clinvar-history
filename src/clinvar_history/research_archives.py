"""Bounded, checksum-verified historical archive acquisition for research."""

import hashlib
import re
import time
from datetime import UTC, datetime
from pathlib import Path
from urllib.parse import urlsplit

import httpx

from .checkpoints import atomic_json


def download_archive(
    url: str, output: Path, *, expected_bytes: int, max_bytes=3_000_000_000
):
    if (
        urlsplit(url).hostname != "ftp.ncbi.nlm.nih.gov"
        or not 0 < expected_bytes <= max_bytes
    ):
        raise ValueError("Archive host or byte budget is invalid")
    output.parent.mkdir(parents=True, exist_ok=True)
    with httpx.Client(timeout=120, follow_redirects=True) as client:
        response = client.get(url + ".md5")
        response.raise_for_status()
        expected = response.text.split()[0].lower()
        if re.fullmatch("[0-9a-f]{32}", expected) is None:
            raise ValueError("Invalid MD5 sidecar")
        if output.exists():
            with output.open("rb") as stream:
                actual = hashlib.file_digest(stream, "md5").hexdigest()
            if output.stat().st_size != expected_bytes or actual != expected:
                raise ValueError("Existing archive differs from official size/checksum")
            return {"url": url, "bytes": expected_bytes, "md5": actual, "reused": True}
        pending = output.with_suffix(output.suffix + ".part")
        offset = pending.stat().st_size if pending.exists() else 0
        if offset > expected_bytes:
            raise ValueError("Partial archive exceeds declared size")
        headers = {"Range": f"bytes={offset}-"} if offset else {}
        digest = hashlib.md5(usedforsecurity=False)
        if offset:
            with pending.open("rb") as stream:
                for chunk in iter(lambda: stream.read(4 * 1024**2), b""):
                    digest.update(chunk)
        with client.stream("GET", url, headers=headers) as response:
            response.raise_for_status()
            if offset and response.status_code == 200:
                offset = 0
                digest = hashlib.md5(usedforsecurity=False)
            elif offset and (
                response.status_code != 206
                or not response.headers.get("Content-Range", "").startswith(
                    f"bytes {offset}-"
                )
            ):
                raise ValueError("Server did not honor the requested resume range")
            total = offset
            last = time.monotonic()
            with pending.open("ab" if offset else "wb") as stream:
                for chunk in response.iter_bytes(chunk_size=4 * 1024**2):
                    total += len(chunk)
                    if total > expected_bytes:
                        raise ValueError("Archive exceeds declared byte budget")
                    stream.write(chunk)
                    digest.update(chunk)
                    if time.monotonic() - last > 30:
                        print(
                            f"{output.name}: {total:,}/{expected_bytes:,} bytes",
                            flush=True,
                        )
                        last = time.monotonic()
            if total != expected_bytes or digest.hexdigest() != expected:
                raise ValueError(
                    "Downloaded archive fails size or official checksum validation"
                )
            receipt = {
                "url": url,
                "bytes": total,
                "md5": expected,
                "http_last_modified": response.headers.get("Last-Modified"),
                "verified_at": datetime.now(UTC).isoformat(),
                "reused": False,
            }
        pending.replace(output)
        atomic_json(receipt, output.with_suffix(output.suffix + ".json"))
        return receipt


def parse_archived_record(record):
    """Parse the old VCV archive schema without interpreting changes as causes."""
    interpreted = record.find("InterpretedRecord")
    if interpreted is None:
        return {
            "variation_id": int(record.get("VariationID")),
            "record_type": record.get("RecordType"),
            "version": record.get("Version"),
            "assertions": [],
        }
    assertions = []
    for assertion in interpreted.findall("./ClinicalAssertionList/ClinicalAssertion"):
        accession = assertion.find("ClinVarAccession")
        if accession is None:
            continue
        interpretation = assertion.find("Interpretation")
        assertions.append(
            {
                "accession": accession.get("Accession"),
                "version": accession.get("Version"),
                "date_created": accession.get("DateCreated"),
                "date_updated": accession.get("DateUpdated"),
                "submission_date": assertion.get("SubmissionDate"),
                "submitter": accession.get("SubmitterName"),
                "classification": interpretation.findtext("Description")
                if interpretation is not None
                else None,
                "review_status": assertion.findtext("ReviewStatus"),
                "last_evaluated": interpretation.get("DateLastEvaluated")
                if interpretation is not None
                else None,
                "record_status": assertion.findtext("RecordStatus"),
                "contributes": assertion.get("ContributesToAggregateClassification"),
                "conditions": sorted(
                    {
                        n.text
                        for n in assertion.findall("./TraitSet/Trait/Name/ElementValue")
                        if n.text
                    }
                ),
                "comments": [n.text for n in assertion.iter("Comment") if n.text],
                "pubmed_ids": sorted(
                    {
                        n.text
                        for n in assertion.iter("ID")
                        if n.get("Source") == "PubMed" and n.text
                    }
                ),
                "assertion_methods": [
                    n.text
                    for n in assertion.iter("Attribute")
                    if n.get("Type") == "AssertionMethod" and n.text
                ],
            }
        )
    return {
        "variation_id": int(record.get("VariationID")),
        "record_type": record.get("RecordType"),
        "accession": record.get("Accession"),
        "version": record.get("Version"),
        "date_last_updated": record.get("DateLastUpdated"),
        "review_status": interpreted.findtext("ReviewStatus"),
        "aggregate_classifications": [
            n.text
            for n in interpreted.findall("./Interpretations/Interpretation/Description")
            if n.text
        ],
        "aggregate_explanations": [
            n.text
            for n in interpreted.findall("./Interpretations/Interpretation/Explanation")
            if n.text
        ],
        "assertions": assertions,
    }


def extract_archived_records(
    archive: Path, variation_ids, output: Path, *, max_record_bytes=8_000_000
):
    """Stream a complete gzip XML, clear siblings, and retain only selected records."""
    import gzip
    import json

    from lxml import etree

    from .checkpoints import file_digest

    targets = {int(v) for v in variation_ids}
    source = json.loads(archive.with_suffix(archive.suffix + ".json").read_text())
    output.mkdir(parents=True, exist_ok=True)
    found = {}
    scanned = 0
    last = time.monotonic()
    with gzip.open(archive, "rb") as stream:
        context = etree.iterparse(
            stream,
            events=("start", "end"),
            resolve_entities=False,
            load_dtd=False,
            no_network=True,
        )
        root_attributes = None
        for event, record in context:
            if root_attributes is None and event == "start":
                root_attributes = dict(record.attrib)
            if event != "end" or record.tag != "VariationArchive":
                continue
            scanned += 1
            vid = int(record.get("VariationID"))
            if vid in targets:
                if vid in found:
                    raise ValueError("Duplicate VariationID in archive")
                content = etree.tostring(record)
                if len(content) > max_record_bytes:
                    found[vid] = {
                        "variation_id": vid,
                        "error": "Selected record exceeds byte limit",
                        "bytes": len(content),
                    }
                else:
                    parsed = parse_archived_record(record)
                    path = output / f"{vid}.xml"
                    path.write_bytes(content)
                    parsed.update(
                        {
                            "bytes": len(content),
                            "sha256": file_digest(path),
                            "source_release_date": root_attributes.get("ReleaseDate"),
                        }
                    )
                    atomic_json(parsed, output / f"{vid}.json")
                    found[vid] = parsed
            record.clear()
            while record.getprevious() is not None:
                del record.getparent()[0]
            if time.monotonic() - last > 30:
                print(
                    f"{archive.name}: scanned {scanned:,}, found {len(found)}/{len(targets)}",
                    flush=True,
                )
                last = time.monotonic()
    result = {
        "archive": str(archive),
        "source": source,
        "root_attributes": root_attributes,
        "scanned_records": scanned,
        "requested": len(targets),
        "found": len(found),
        "missing": sorted(targets - found.keys()),
        "records": [found[v] for v in sorted(found)],
    }
    atomic_json(result, output / "extraction.json")
    return result
