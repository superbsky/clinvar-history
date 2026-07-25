import gzip
import hashlib
from collections.abc import Iterator
from contextlib import contextmanager
from dataclasses import replace
from datetime import date
from pathlib import Path

import duckdb
import httpx
import pyarrow.parquet as pq
import pytest
from click.testing import CliRunner

import clinvar_history.pipeline as pipeline_module
from clinvar_history.cli import main
from clinvar_history.manifest import ReleaseFile, ReleaseKind
from clinvar_history.pipeline import (
    GzipLineStream,
    IngestionSummary,
    fetch_release_md5,
    ingest_vcf_releases,
    iter_gzip_lines,
    resolve_release_checksums,
)
from clinvar_history.storage import read_vcf_state
from clinvar_history.vcf_parser import parse_vcf

FIXTURES = Path(__file__).parent / "fixtures"


def _release(day: int, fixture: str) -> ReleaseFile:
    return ReleaseFile(
        release_date=date(2024, 1, day),
        kind=ReleaseKind.VCF,
        url=f"fixture://{fixture}",
        size_bytes=(FIXTURES / fixture).stat().st_size,
    )


@contextmanager
def _fixture_source(
    release: ReleaseFile,
) -> Iterator[Iterator[bytes]]:
    fixture = release.url.removeprefix("fixture://")
    with (FIXTURES / fixture).open("rb") as stream:
        yield iter(stream)


def _ledger_rows(root: Path) -> list[tuple]:
    with duckdb.connect(str(root / "clinvar.duckdb"), read_only=True) as connection:
        return connection.execute(
            """
            SELECT release_date, url, status, error
            FROM pipeline_runs
            ORDER BY release_date
            """
        ).fetchall()


def test_ingestion_sorts_releases_and_replaces_latest_state(tmp_path):
    older = _release(1, "release_a.vcf")
    newer = _release(8, "release_b.vcf")
    calls = []

    @contextmanager
    def source(release):
        calls.append(release.release_date)
        with _fixture_source(release) as lines:
            yield lines

    summary = ingest_vcf_releases([newer, older], root=tmp_path, source=source)

    assert summary == IngestionSummary(processed=2, skipped=0)
    assert calls == [date(2024, 1, 1), date(2024, 1, 8)]
    with (FIXTURES / "release_b.vcf").open() as stream:
        expected = list(parse_vcf(stream))
    assert read_vcf_state(tmp_path / "state" / "vcv_latest.parquet") == expected
    event_files = sorted((tmp_path / "events" / "vcv").rglob("*.parquet"))
    assert len(event_files) == 2
    assert sum(pq.read_table(path).num_rows for path in event_files) == 6
    assert [row[2] for row in _ledger_rows(tmp_path)] == ["completed", "completed"]


def test_completed_releases_are_skipped_without_duplicate_events(tmp_path):
    releases = [_release(1, "release_a.vcf"), _release(8, "release_b.vcf")]
    first = ingest_vcf_releases(releases, root=tmp_path, source=_fixture_source)
    rows_before = sum(
        pq.read_table(path).num_rows
        for path in (tmp_path / "events" / "vcv").rglob("*.parquet")
    )

    def unexpected_source(_release):
        raise AssertionError("completed releases must not be downloaded again")

    second = ingest_vcf_releases(
        releases,
        root=tmp_path,
        source=unexpected_source,
    )
    rows_after = sum(
        pq.read_table(path).num_rows
        for path in (tmp_path / "events" / "vcv").rglob("*.parquet")
    )

    assert first == IngestionSummary(processed=2, skipped=0)
    assert second == IngestionSummary(processed=0, skipped=2)
    assert rows_after == rows_before
    assert len(_ledger_rows(tmp_path)) == 2


def test_failed_release_is_recorded_and_can_be_retried(tmp_path):
    older = _release(1, "release_a.vcf")
    newer = _release(8, "release_b.vcf")

    @contextmanager
    def failing_source(release):
        if release == newer:
            raise OSError("temporary download failure")
        with _fixture_source(release) as lines:
            yield lines

    with pytest.raises(OSError, match="temporary download failure"):
        ingest_vcf_releases([older, newer], root=tmp_path, source=failing_source)

    failed_rows = _ledger_rows(tmp_path)
    assert [row[2] for row in failed_rows] == ["completed", "failed"]
    assert "temporary download failure" in failed_rows[1][3]

    summary = ingest_vcf_releases(
        [older, newer],
        root=tmp_path,
        source=_fixture_source,
    )

    assert summary == IngestionSummary(processed=1, skipped=1)
    assert [row[2] for row in _ledger_rows(tmp_path)] == ["completed", "completed"]
    with (FIXTURES / "release_b.vcf").open() as stream:
        assert read_vcf_state(tmp_path / "state" / "vcv_latest.parquet") == list(
            parse_vcf(stream)
        )


def test_retry_after_state_write_rebuilds_events_from_last_completed_checkpoint(
    tmp_path,
    monkeypatch,
):
    older = _release(1, "release_a.vcf")
    newer = _release(8, "release_b.vcf")
    ingest_vcf_releases([older], root=tmp_path, source=_fixture_source)
    real_mark_completed = pipeline_module._mark_completed

    def fail_after_artifact_writes(_connection, _release, _source_md5):
        raise OSError("crash before ledger commit")

    monkeypatch.setattr(
        pipeline_module,
        "_mark_completed",
        fail_after_artifact_writes,
    )
    with pytest.raises(OSError, match="crash before ledger commit"):
        ingest_vcf_releases([older, newer], root=tmp_path, source=_fixture_source)

    monkeypatch.setattr(pipeline_module, "_mark_completed", real_mark_completed)
    ingest_vcf_releases([older, newer], root=tmp_path, source=_fixture_source)

    newer_events = (
        tmp_path / "events" / "vcv" / "release_date=2024-01-08" / "events.parquet"
    )
    assert pq.read_table(newer_events).num_rows == 3


def test_expanded_manifest_cannot_backfill_before_completed_release(tmp_path):
    older = _release(1, "release_a.vcf")
    newer = _release(8, "release_b.vcf")
    ingest_vcf_releases([newer], root=tmp_path, source=_fixture_source)

    with pytest.raises(ValueError, match="completed releases must form"):
        ingest_vcf_releases(
            [older, newer],
            root=tmp_path,
            source=_fixture_source,
        )

    with (FIXTURES / "release_b.vcf").open() as stream:
        assert read_vcf_state(tmp_path / "state" / "vcv_latest.parquet") == list(
            parse_vcf(stream)
        )


def test_gzip_streaming_validates_compressed_checksum():
    compressed = gzip.compress(b"first\nsecond\n")
    expected_md5 = hashlib.md5(compressed, usedforsecurity=False).hexdigest()

    lines = list(
        iter_gzip_lines(
            [compressed[:3], compressed[3:11], compressed[11:]],
            expected_md5=expected_md5,
        )
    )

    assert lines == [b"first\n", b"second\n"]
    with pytest.raises(ValueError, match="checksum"):
        list(iter_gzip_lines([compressed], expected_md5="0" * 32))


def test_gzip_streaming_reads_all_concatenated_bgzf_members():
    compressed = gzip.compress(b"first\n") + gzip.compress(b"second\n")

    lines = list(
        iter_gzip_lines(
            [compressed[:9], compressed[9:25], compressed[25:]],
        )
    )

    assert lines == [b"first\n", b"second\n"]


def test_gzip_stream_exposes_calculated_checksum_without_sidecar():
    compressed = gzip.compress(b"content\n")
    stream = GzipLineStream([compressed])

    assert list(stream) == [b"content\n"]
    assert (
        stream.calculated_md5
        == hashlib.md5(
            compressed,
            usedforsecurity=False,
        ).hexdigest()
    )


def test_checksum_is_discovered_when_sidecar_exists():
    checksum = "0123456789abcdef0123456789abcdef"

    def handler(request):
        assert str(request.url).endswith(".md5")
        return httpx.Response(200, text=f"{checksum}  release.vcf.gz\n")

    release = _release(1, "release_a.vcf")
    with httpx.Client(transport=httpx.MockTransport(handler)) as client:
        assert fetch_release_md5(release, client=client) == checksum


def test_missing_checksum_sidecar_is_optional():
    def handler(_request):
        return httpx.Response(404)

    release = _release(1, "release_a.vcf")
    with httpx.Client(transport=httpx.MockTransport(handler)) as client:
        assert fetch_release_md5(release, client=client) is None


def test_resolved_checksum_is_persisted_in_ledger(tmp_path):
    checksum = "0123456789abcdef0123456789abcdef"

    def handler(_request):
        return httpx.Response(200, text=f"{checksum}  release.vcf.gz\n")

    release = _release(1, "release_a.vcf")
    with httpx.Client(transport=httpx.MockTransport(handler)) as client:
        resolved = resolve_release_checksums([release], client=client)

    ingest_vcf_releases(resolved, root=tmp_path, source=_fixture_source)
    with duckdb.connect(
        str(tmp_path / "clinvar.duckdb"),
        read_only=True,
    ) as connection:
        stored = connection.execute("SELECT source_md5 FROM pipeline_runs").fetchone()[
            0
        ]

    assert resolved[0].md5 == checksum
    assert stored == checksum


def test_completed_release_rejects_changed_official_checksum(tmp_path):
    first_checksum = "0123456789abcdef0123456789abcdef"
    changed_checksum = "abcdef0123456789abcdef0123456789"
    release = replace(_release(1, "release_a.vcf"), md5=first_checksum)
    ingest_vcf_releases([release], root=tmp_path, source=_fixture_source)

    with pytest.raises(ValueError, match="checksum changed"):
        ingest_vcf_releases(
            [replace(release, md5=changed_checksum)],
            root=tmp_path,
            source=_fixture_source,
        )


def test_calculated_checksum_is_recorded_without_official_sidecar(tmp_path):
    release = _release(1, "release_a.vcf")
    compressed = gzip.compress((FIXTURES / "release_a.vcf").read_bytes())
    expected = hashlib.md5(compressed, usedforsecurity=False).hexdigest()

    @contextmanager
    def source(_release):
        yield GzipLineStream([compressed])

    ingest_vcf_releases([release], root=tmp_path, source=source)
    with duckdb.connect(
        str(tmp_path / "clinvar.duckdb"),
        read_only=True,
    ) as connection:
        stored = connection.execute("SELECT source_md5 FROM pipeline_runs").fetchone()[
            0
        ]

    assert stored == expected


def test_ingest_vcf_cli_exposes_bounded_storage_options():
    result = CliRunner().invoke(main, ["ingest-vcf", "--help"])

    assert result.exit_code == 0
    assert "--manifest" in result.output
    assert "--root" in result.output
    assert "--keep-raw" in result.output
