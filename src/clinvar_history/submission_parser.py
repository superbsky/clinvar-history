"""Streaming parser for ClinVar's submission_summary tab-delimited files."""

import csv
import re
import sys
import time
from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from datetime import date

from .classification import ClassificationState, normalize


@dataclass(frozen=True, slots=True)
class SubmissionRecord:
    variation_id: int
    scv_accession: str
    scv_version: int | None
    raw_classification: str | None
    classification: ClassificationState
    classification_ambiguous: bool
    last_evaluated: date | None
    review_status: str | None
    conditions: str | None
    submitter: str | None
    collection_method: str | None


_SCV_PATTERN = re.compile(r"^(SCV\d+)\s*(?:\.\s*(\d+))?$", re.IGNORECASE)
_REQUIRED_COLUMNS = {
    "VariationID",
    "ClinicalSignificance",
    "ReviewStatus",
    "SCV",
}


def _as_text(line: str | bytes) -> str:
    return line.decode("utf-8") if isinstance(line, bytes) else line


def _optional(value: str | None) -> str | None:
    if value is None:
        return None
    stripped = value.strip()
    return None if stripped in {"", "-", "na"} else stripped


def _optional_date(value: str | None) -> date | None:
    cleaned = _optional(value)
    if cleaned is None:
        return None
    try:
        return date.fromisoformat(cleaned)
    except ValueError:
        parsed = time.strptime(cleaned, "%b %d, %Y")
        return date(parsed.tm_year, parsed.tm_mon, parsed.tm_mday)


def parse_scv_accession(raw: str) -> tuple[str, int | None]:
    """Split tolerant historical SCV accession/version formatting."""
    match = _SCV_PATTERN.fullmatch(raw.strip())
    if match is None:
        raise ValueError(f"Invalid SCV accession: {raw!r}")
    version = int(match.group(2)) if match.group(2) is not None else None
    return match.group(1).upper(), version


def parse_submissions(
    stream: Iterable[str | bytes],
) -> Iterator[SubmissionRecord]:
    """Yield compact current SCV states using column names across schema versions."""
    csv.field_size_limit(sys.maxsize)
    lines = iter(stream)
    header: str | None = None
    for raw_line in lines:
        line = _as_text(raw_line).rstrip("\r\n")
        candidate = line.removeprefix("#")
        if candidate.startswith("VariationID\t"):
            header = candidate
            break

    if header is None:
        raise ValueError("submission_summary header was not found")

    columns = header.split("\t")
    missing = sorted(_REQUIRED_COLUMNS - set(columns))
    if missing:
        raise ValueError(f"Missing required columns: {', '.join(missing)}")

    reader = csv.DictReader(
        (_as_text(line) for line in lines),
        fieldnames=columns,
        delimiter="\t",
        quoting=csv.QUOTE_NONE,
    )
    for row_number, row in enumerate(reader, start=2):
        if not row or not row.get("VariationID"):
            continue
        if row["VariationID"].startswith("#"):
            continue

        if row.get("SCV") is None:
            raise ValueError(
                "Missing SCV at submission data row "
                f"{row_number} for VariationID={row.get('VariationID')!r}; "
                f"populated_columns={sum(value is not None for value in row.values())}"
            )
        accession, version = parse_scv_accession(row["SCV"])
        raw_classification = _optional(row.get("ClinicalSignificance"))
        normalized = normalize(raw_classification)
        conditions = _optional(row.get("ReportedPhenotypeInfo")) or _optional(
            row.get("SubmittedPhenotypeInfo")
        )

        yield SubmissionRecord(
            variation_id=int(row["VariationID"]),
            scv_accession=accession,
            scv_version=version,
            raw_classification=raw_classification,
            classification=normalized.state,
            classification_ambiguous=normalized.ambiguous,
            last_evaluated=_optional_date(row.get("DateLastEvaluated")),
            review_status=_optional(row.get("ReviewStatus")),
            conditions=conditions,
            submitter=_optional(row.get("Submitter")),
            collection_method=_optional(row.get("CollectionMethod")),
        )
