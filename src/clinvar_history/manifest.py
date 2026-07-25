"""Discover immutable ClinVar release files from official FTP indexes."""

import re
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from datetime import date
from enum import StrEnum
from pathlib import Path
from urllib.parse import urljoin

import httpx
import pyarrow as pa
import pyarrow.parquet as pq


class ReleaseKind(StrEnum):
    VCF = "vcf"
    SUBMISSION_SUMMARY = "submission_summary"


@dataclass(frozen=True, slots=True)
class ReleaseFile:
    release_date: date
    kind: ReleaseKind
    url: str
    size_bytes: int
    source_month: str | None = None
    md5: str | None = None


_SIZE_PATTERN = re.compile(r"^(?P<number>[\d,.]+)\s*(?P<unit>[KMG]?)$")
_VCF_PATTERN = re.compile(
    r'href="(?P<name>clinvar_(?P<date>\d{8})\.vcf\.gz)"'
    r".*?\s(?P<size>[\d.]+[KMG]?)\s*(?:</pre>)?$",
    re.MULTILINE,
)
_SUBMISSION_PATTERN = re.compile(
    r'href="(?P<name>submission_summary_(?P<date>\d{4}-\d{2})\.txt\.gz)"'
    r"[^>]*>[^<]*</a></td><td class=\"rj\">(?P<size>[\d,]+)</td>"
    r"<td>\d{4}-\d{2}-\d{2}(?:\s+\d{2}:\d{2}:\d{2})?</td>",
)


def parse_size(value: str) -> int:
    """Parse Apache FTP byte counts and binary K/M/G suffixes."""
    match = _SIZE_PATTERN.fullmatch(value.strip())
    if match is None:
        raise ValueError(f"Unsupported FTP size: {value!r}")

    number = float(match.group("number").replace(",", ""))
    multiplier = {
        "": 1,
        "K": 1024,
        "M": 1024**2,
        "G": 1024**3,
    }[match.group("unit")]
    return int(number * multiplier)


def _unique_sorted(releases: Iterable[ReleaseFile]) -> list[ReleaseFile]:
    unique = {(item.kind, item.url): item for item in releases}
    return sorted(
        unique.values(),
        key=lambda item: (item.release_date, item.kind.value, item.url),
    )


def parse_vcf_archive_index(html: str, base_url: str) -> list[ReleaseFile]:
    """Parse complete non-PAPU VCF releases from an Apache index."""
    releases = []
    for match in _VCF_PATTERN.finditer(html):
        releases.append(
            ReleaseFile(
                release_date=date.fromisoformat(
                    f"{match['date'][:4]}-{match['date'][4:6]}-{match['date'][6:]}"
                ),
                kind=ReleaseKind.VCF,
                url=urljoin(base_url, match["name"]),
                size_bytes=parse_size(match["size"]),
            )
        )
    return _unique_sorted(releases)


def parse_submission_archive_index(html: str, base_url: str) -> list[ReleaseFile]:
    """Parse monthly submission-summary releases from an FTP table index."""
    releases = []
    for match in _SUBMISSION_PATTERN.finditer(html):
        releases.append(
            ReleaseFile(
                # The FTP "Released" field can reflect a later archive
                # republishing date. The filename's source month is the stable
                # longitudinal key for these monthly snapshots.
                release_date=date.fromisoformat(f"{match['date']}-01"),
                kind=ReleaseKind.SUBMISSION_SUMMARY,
                url=urljoin(base_url, match["name"]),
                size_bytes=parse_size(match["size"]),
                source_month=match["date"],
            )
        )
    return _unique_sorted(releases)


def discover_manifest(
    index_urls: Sequence[tuple[ReleaseKind, str]],
    *,
    start: date,
    end: date,
    client: httpx.Client,
) -> list[ReleaseFile]:
    """Fetch release indexes and return a bounded deterministic manifest."""
    releases: list[ReleaseFile] = []
    parsers = {
        ReleaseKind.VCF: parse_vcf_archive_index,
        ReleaseKind.SUBMISSION_SUMMARY: parse_submission_archive_index,
    }
    for kind, url in index_urls:
        response = client.get(url)
        response.raise_for_status()
        releases.extend(parsers[kind](response.text, url))

    return [
        item for item in _unique_sorted(releases) if start <= item.release_date <= end
    ]


def write_manifest(releases: Sequence[ReleaseFile], output: Path) -> None:
    """Write the release manifest as a compact Parquet table."""
    output.parent.mkdir(parents=True, exist_ok=True)
    table = pa.table(
        {
            "release_date": pa.array(
                [item.release_date for item in releases], type=pa.date32()
            ),
            "source_month": [item.source_month for item in releases],
            "kind": [item.kind.value for item in releases],
            "url": [item.url for item in releases],
            "size_bytes": [item.size_bytes for item in releases],
            "md5": [item.md5 for item in releases],
        }
    )
    pq.write_table(table, output, compression="zstd")


def read_manifest(path: Path) -> list[ReleaseFile]:
    """Restore release metadata from a Parquet manifest."""
    return [
        ReleaseFile(
            release_date=row["release_date"],
            kind=ReleaseKind(row["kind"]),
            url=row["url"],
            size_bytes=row["size_bytes"],
            source_month=row.get("source_month"),
            md5=row.get("md5"),
        )
        for row in pq.read_table(path).to_pylist()
    ]


def default_index_urls(start: date, end: date) -> list[tuple[ReleaseKind, str]]:
    """Return official index URLs that can contain the requested date range."""
    urls = [
        (
            ReleaseKind.VCF,
            f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/vcf_GRCh38/archive_2.0/{year}/",
        )
        for year in range(start.year, end.year + 1)
    ]
    urls.extend(
        (
            ReleaseKind.SUBMISSION_SUMMARY,
            f"https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/{year}/",
        )
        for year in range(start.year, min(end.year, 2024) + 1)
    )
    if end.year >= 2025:
        urls.append(
            (
                ReleaseKind.SUBMISSION_SUMMARY,
                "https://ftp.ncbi.nlm.nih.gov/pub/clinvar/tab_delimited/archive/",
            )
        )
    return urls
