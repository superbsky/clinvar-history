from datetime import date
from pathlib import Path

import httpx
import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

from clinvar_history.cli import main
from clinvar_history.manifest import (
    ReleaseKind,
    discover_manifest,
    parse_size,
    parse_submission_archive_index,
    parse_vcf_archive_index,
    read_manifest,
    write_manifest,
)

FIXTURES = Path(__file__).parent / "fixtures"


@pytest.mark.parametrize(
    ("value", "expected"),
    [
        ("132", 132),
        ("500K", 500 * 1024),
        ("92M", 92 * 1024**2),
        ("1.5G", int(1.5 * 1024**3)),
        ("201,500,000", 201_500_000),
    ],
)
def test_parse_size_supports_ftp_units_and_exact_bytes(value, expected):
    assert parse_size(value) == expected


def test_parse_vcf_archive_excludes_sidecars_papu_duplicates_and_sorts():
    html = (FIXTURES / "vcf_archive_index.html").read_text()

    releases = parse_vcf_archive_index(
        html, "https://ftp.example.test/vcf/archive/2024/"
    )

    assert [(item.release_date, item.size_bytes) for item in releases] == [
        (date(2024, 1, 1), 91 * 1024**2),
        (date(2024, 1, 8), 92 * 1024**2),
    ]
    assert all(item.kind is ReleaseKind.VCF for item in releases)
    assert releases[0].url.endswith("/clinvar_20240101.vcf.gz")


def test_parse_submission_archive_selects_only_submission_summaries():
    html = (FIXTURES / "submission_archive_index.html").read_text()

    releases = parse_submission_archive_index(
        html, "https://ftp.example.test/tab/archive/"
    )

    assert [(item.release_date, item.size_bytes) for item in releases] == [
        (date(2024, 1, 1), 200_000_000),
        (date(2024, 2, 1), 201_500_000),
    ]
    assert [item.source_month for item in releases] == ["2024-01", "2024-02"]
    assert all(item.kind is ReleaseKind.SUBMISSION_SUMMARY for item in releases)


def test_discover_manifest_fetches_filters_and_deduplicates_indexes():
    vcf_html = (FIXTURES / "vcf_archive_index.html").read_text()
    submission_html = (FIXTURES / "submission_archive_index.html").read_text()

    def handler(request):
        body = vcf_html if "vcf" in str(request.url) else submission_html
        return httpx.Response(200, text=body)

    with httpx.Client(transport=httpx.MockTransport(handler)) as client:
        releases = discover_manifest(
            [
                (ReleaseKind.VCF, "https://example.test/vcf"),
                (ReleaseKind.SUBMISSION_SUMMARY, "https://example.test/submissions"),
            ],
            start=date(2024, 1, 8),
            end=date(2024, 2, 8),
            client=client,
        )

    assert [(item.kind, item.release_date) for item in releases] == [
        (ReleaseKind.VCF, date(2024, 1, 8)),
        (ReleaseKind.SUBMISSION_SUMMARY, date(2024, 2, 1)),
    ]


def test_write_manifest_creates_queryable_parquet(tmp_path):
    html = (FIXTURES / "vcf_archive_index.html").read_text()
    releases = parse_vcf_archive_index(html, "https://example.test/vcf/")
    output = tmp_path / "manifest.parquet"

    write_manifest(releases, output)

    table = pq.read_table(output)
    assert table.num_rows == 2
    assert table.column_names == [
        "release_date",
        "source_month",
        "kind",
        "url",
        "size_bytes",
        "md5",
    ]
    assert read_manifest(output) == releases


def test_build_manifest_command_exposes_bounded_date_options():
    result = CliRunner().invoke(main, ["build-manifest", "--help"])

    assert result.exit_code == 0
    assert "--start" in result.output
    assert "--end" in result.output
    assert "--output" in result.output
