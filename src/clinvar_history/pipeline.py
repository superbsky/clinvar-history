"""Resumable, disk-bounded ingestion orchestration."""

import hashlib
import json
import os
import re
import tempfile
import zlib
from collections.abc import Callable, Iterable, Iterator, Sequence
from contextlib import AbstractContextManager, contextmanager
from dataclasses import dataclass, replace
from datetime import UTC, datetime
from pathlib import Path
from typing import BinaryIO, Protocol
from urllib.parse import urlparse

import duckdb
import httpx

from .diff import diff_vcf_states
from .manifest import ReleaseFile, ReleaseKind
from .scv_diff import consolidate_scv_records, diff_scv_states
from .storage import (
    read_scv_state,
    read_vcf_state,
    write_scv_events,
    write_scv_state,
    write_vcf_state,
    write_vcv_events,
)
from .submission_parser import parse_submissions
from .vcf_parser import parse_vcf


class ReleaseSource(Protocol):
    """Open a release as an iterable of decompressed VCF lines."""

    def __call__(
        self, release: ReleaseFile
    ) -> AbstractContextManager[Iterable[str | bytes]]: ...


@dataclass(frozen=True, slots=True)
class IngestionSummary:
    processed: int
    skipped: int


@dataclass(slots=True)
class GzipLineStream:
    """Iterable decompressed lines with a checksum populated after exhaustion."""

    chunks: Iterable[bytes]
    expected_md5: str | None = None
    calculated_md5: str | None = None

    def __iter__(self) -> Iterator[bytes]:
        def capture(checksum: str) -> None:
            self.calculated_md5 = checksum

        yield from iter_gzip_lines(
            self.chunks,
            expected_md5=self.expected_md5,
            checksum_sink=capture,
        )


def fetch_release_md5(
    release: ReleaseFile,
    *,
    client: httpx.Client,
) -> str | None:
    """Return an official sidecar checksum when the release publishes one."""
    response = client.get(f"{release.url}.md5")
    if response.status_code == httpx.codes.NOT_FOUND:
        return None
    response.raise_for_status()
    fields = response.text.split()
    if not fields or re.fullmatch(r"[0-9a-fA-F]{32}", fields[0]) is None:
        raise ValueError(f"Invalid MD5 sidecar for {release.url}")
    return fields[0].lower()


def resolve_release_checksums(
    releases: Sequence[ReleaseFile],
    *,
    client: httpx.Client,
) -> list[ReleaseFile]:
    """Populate missing official checksums while preserving manifest order."""
    return [
        release
        if release.md5 is not None
        else replace(release, md5=fetch_release_md5(release, client=client))
        for release in releases
    ]


def iter_gzip_lines(
    chunks: Iterable[bytes],
    *,
    expected_md5: str | None = None,
    checksum_sink: Callable[[str], None] | None = None,
) -> Iterator[bytes]:
    """Incrementally decompress gzip chunks and optionally validate their MD5."""
    digest = hashlib.md5(usedforsecurity=False)
    decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
    pending = b""

    for chunk in chunks:
        digest.update(chunk)
        compressed = chunk
        while compressed:
            if decompressor.eof:
                decompressor = zlib.decompressobj(16 + zlib.MAX_WBITS)
            pending += decompressor.decompress(compressed)
            compressed = decompressor.unused_data
            while b"\n" in pending:
                line, pending = pending.split(b"\n", 1)
                yield line + b"\n"

    pending += decompressor.flush()
    if pending:
        yield pending
    if not decompressor.eof:
        raise ValueError("Incomplete gzip stream")
    calculated_md5 = digest.hexdigest().lower()
    if checksum_sink is not None:
        checksum_sink(calculated_md5)
    if expected_md5 is not None and calculated_md5 != expected_md5.lower():
        raise ValueError("Compressed release checksum does not match expected MD5")


@contextmanager
def http_vcf_source(
    release: ReleaseFile,
    *,
    client: httpx.Client,
    keep_raw: bool = False,
    raw_dir: Path | None = None,
    expected_md5: str | None = None,
) -> Iterator[Iterable[bytes]]:
    """Stream and decompress one HTTP release without retaining it by default."""
    raw_stream: BinaryIO | None = None
    if keep_raw:
        if raw_dir is None:
            raise ValueError("raw_dir is required when keep_raw is enabled")
        raw_dir.mkdir(parents=True, exist_ok=True)
        filename = Path(urlparse(release.url).path).name
        raw_stream = (raw_dir / filename).open("wb")

    try:
        with client.stream("GET", release.url) as response:
            response.raise_for_status()

            def compressed_chunks() -> Iterator[bytes]:
                for chunk in response.iter_bytes():
                    if raw_stream is not None:
                        raw_stream.write(chunk)
                    yield chunk

            yield GzipLineStream(
                compressed_chunks(),
                expected_md5,
            )
    finally:
        if raw_stream is not None:
            raw_stream.close()


def _create_ledger(connection: duckdb.DuckDBPyConnection) -> None:
    connection.execute(
        """
        CREATE TABLE IF NOT EXISTS pipeline_runs (
            kind VARCHAR NOT NULL,
            url VARCHAR NOT NULL,
            release_date DATE NOT NULL,
            status VARCHAR NOT NULL,
            started_at TIMESTAMPTZ NOT NULL,
            completed_at TIMESTAMPTZ,
            error VARCHAR,
            source_md5 VARCHAR,
            official_md5 VARCHAR,
            PRIMARY KEY (kind, url)
        )
        """
    )
    connection.execute(
        "ALTER TABLE pipeline_runs ADD COLUMN IF NOT EXISTS source_md5 VARCHAR"
    )
    connection.execute(
        "ALTER TABLE pipeline_runs ADD COLUMN IF NOT EXISTS official_md5 VARCHAR"
    )
    connection.execute(
        "ALTER TABLE pipeline_runs ADD COLUMN IF NOT EXISTS release_identity VARCHAR"
    )
    connection.execute(
        "ALTER TABLE pipeline_runs ADD COLUMN IF NOT EXISTS processing_versions VARCHAR"
    )


def _validate_completed_identity(connection, releases, kind) -> None:
    from .checkpoints import processing_versions, release_identity, validate_manifest

    if releases:
        validate_manifest(releases)
    columns = {r[0] for r in connection.execute("DESCRIBE pipeline_runs").fetchall()}
    if not {"release_identity", "processing_versions"} <= columns:
        raise ValueError("Legacy dataset needs adopt-legacy before checkpoint reuse")
    stored = {
        url: (identity, versions, release_date)
        for url, identity, versions, release_date in connection.execute(
            "SELECT url, release_identity, processing_versions, release_date FROM pipeline_runs "
            "WHERE kind = ? AND status = 'completed'",
            [kind.value],
        ).fetchall()
    }
    ingestion_versions = processing_versions()
    for release in releases:
        if release.url not in stored:
            continue
        identity, versions, stored_date = stored[release.url]
        if identity is None or versions is None:
            raise ValueError(
                "Legacy dataset needs adopt-legacy before checkpoint reuse"
            )
        if json.loads(versions) != ingestion_versions:
            raise ValueError(
                "Dataset processing versions changed; rebuild in a new root"
            )
        if identity != release_identity(release) or stored_date != release.release_date:
            raise ValueError("Completed release identity changed; use a new root")


def _completed_urls(
    connection: duckdb.DuckDBPyConnection,
    kind: ReleaseKind,
) -> set[str]:
    return {
        row[0]
        for row in connection.execute(
            """
            SELECT url
            FROM pipeline_runs
            WHERE kind = ? AND status = 'completed'
            """,
            [kind.value],
        ).fetchall()
    }


def _validate_completed_checksums(
    connection: duckdb.DuckDBPyConnection,
    releases: Sequence[ReleaseFile],
    kind: ReleaseKind,
) -> None:
    stored = dict(
        connection.execute(
            """
            SELECT url, source_md5
            FROM pipeline_runs
            WHERE kind = ? AND status = 'completed'
            """,
            [kind.value],
        ).fetchall()
    )
    for release in releases:
        if (
            release.url in stored
            and release.md5 is not None
            and stored[release.url] != release.md5
        ):
            raise ValueError(
                f"Official checksum changed for completed release {release.url}; "
                "use a new root to rebuild its artifacts"
            )


def _validate_completed_prefix(
    releases: Sequence[ReleaseFile],
    completed: set[str],
) -> None:
    manifest_urls = {release.url for release in releases}
    if completed - manifest_urls:
        raise ValueError(
            "Manifest must include the completed release history for this root"
        )
    incomplete_seen = False
    for release in releases:
        if release.url not in completed:
            incomplete_seen = True
        elif incomplete_seen:
            raise ValueError(
                "completed releases must form a chronological manifest prefix; "
                "use a new root to backfill earlier releases"
            )


def _mark_started(
    connection: duckdb.DuckDBPyConnection,
    release: ReleaseFile,
) -> None:
    from .checkpoints import canonical, processing_versions, release_identity

    connection.execute(
        """
        INSERT INTO pipeline_runs
            (
                kind, url, release_date, status, started_at,
                completed_at, error, source_md5, official_md5,
                release_identity, processing_versions
            )
        VALUES (?, ?, ?, 'processing', ?, NULL, NULL, NULL, ?, ?, ?)
        ON CONFLICT (kind, url) DO UPDATE SET
            release_date = excluded.release_date,
            status = excluded.status,
            started_at = excluded.started_at,
            completed_at = NULL,
            error = NULL,
            source_md5 = NULL,
            official_md5 = excluded.official_md5,
            release_identity = excluded.release_identity,
            processing_versions = excluded.processing_versions
        """,
        [
            release.kind.value,
            release.url,
            release.release_date,
            datetime.now(UTC),
            release.md5,
            release_identity(release),
            canonical(processing_versions()),
        ],
    )


def _mark_completed(
    connection: duckdb.DuckDBPyConnection,
    release: ReleaseFile,
    source_md5: str | None,
) -> None:
    connection.execute(
        """
        UPDATE pipeline_runs
        SET
            status = 'completed',
            completed_at = ?,
            error = NULL,
            source_md5 = ?
        WHERE kind = ? AND url = ?
        """,
        [datetime.now(UTC), source_md5, release.kind.value, release.url],
    )


def _mark_failed(
    connection: duckdb.DuckDBPyConnection,
    release: ReleaseFile,
    error: Exception,
) -> None:
    connection.execute(
        """
        UPDATE pipeline_runs
        SET status = 'failed', completed_at = ?, error = ?
        WHERE kind = ? AND url = ?
        """,
        [datetime.now(UTC), str(error), release.kind.value, release.url],
    )


def _vcv_event_path(root: Path, release: ReleaseFile) -> Path:
    return (
        root
        / "events"
        / "vcv"
        / f"release_date={release.release_date.isoformat()}"
        / "events.parquet"
    )


def _checkpoint_path(
    root: Path,
    kind: ReleaseKind,
    release: ReleaseFile,
) -> Path:
    url_hash = hashlib.sha256(release.url.encode()).hexdigest()[:12]
    return (
        root
        / "state"
        / "checkpoints"
        / kind.value
        / f"{release.release_date.isoformat()}-{url_hash}.parquet"
    )


def _publish_latest(checkpoint: Path, latest: Path) -> None:
    """Atomically point the public latest path at a completed checkpoint inode."""
    latest.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=latest.parent,
            prefix=f".{latest.name}.",
            suffix=".link",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
        temporary_path.unlink()
        os.link(checkpoint, temporary_path)
        os.replace(temporary_path, latest)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def _prune_checkpoints(root: Path, kind: ReleaseKind, keep: Path) -> None:
    checkpoint_dir = root / "state" / "checkpoints" / kind.value
    for path in checkpoint_dir.glob("*.parquet"):
        if path != keep:
            path.unlink()


def ingest_vcf_releases(
    releases: Sequence[ReleaseFile],
    *,
    root: Path,
    source: ReleaseSource,
) -> IngestionSummary:
    """Ingest VCF releases chronologically and resume from a DuckDB ledger."""
    root.mkdir(parents=True, exist_ok=True)
    latest_state = root / "state" / "vcv_latest.parquet"
    vcf_releases = sorted(
        (release for release in releases if release.kind is ReleaseKind.VCF),
        key=lambda release: (release.release_date, release.url),
    )
    processed = 0
    skipped = 0

    with duckdb.connect(str(root / "clinvar.duckdb")) as connection:
        _create_ledger(connection)
        completed = _completed_urls(connection, ReleaseKind.VCF)
        _validate_completed_checksums(
            connection,
            vcf_releases,
            ReleaseKind.VCF,
        )
        _validate_completed_prefix(vcf_releases, completed)
        _validate_completed_identity(connection, vcf_releases, ReleaseKind.VCF)
        completed_releases = [
            release for release in vcf_releases if release.url in completed
        ]
        if completed_releases:
            last_completed = max(
                completed_releases,
                key=lambda release: (release.release_date, release.url),
            )
            checkpoint = _checkpoint_path(root, ReleaseKind.VCF, last_completed)
            previous = read_vcf_state(
                checkpoint if checkpoint.exists() else latest_state
            )
        else:
            previous = []

        for release in vcf_releases:
            if release.url in completed:
                skipped += 1
                continue

            _mark_started(connection, release)
            try:
                with source(release) as lines:
                    current = list(parse_vcf(lines))
                source_md5 = getattr(lines, "calculated_md5", None) or release.md5
                events = list(
                    diff_vcf_states(
                        previous,
                        current,
                        observed_at=release.release_date,
                    )
                )
                checkpoint = _checkpoint_path(root, ReleaseKind.VCF, release)
                write_vcv_events(events, _vcv_event_path(root, release))
                write_vcf_state(current, checkpoint)
                _publish_latest(checkpoint, latest_state)
                _mark_completed(connection, release, source_md5)
            except Exception as error:
                _mark_failed(connection, release, error)
                raise

            completed.add(release.url)
            previous = current
            processed += 1
            _prune_checkpoints(root, ReleaseKind.VCF, checkpoint)

    return IngestionSummary(processed=processed, skipped=skipped)


def _scv_event_path(root: Path, release: ReleaseFile) -> Path:
    return (
        root
        / "events"
        / "scv"
        / f"release_date={release.release_date.isoformat()}"
        / "events.parquet"
    )


def ingest_submission_releases(
    releases: Sequence[ReleaseFile],
    *,
    root: Path,
    source: ReleaseSource,
) -> IngestionSummary:
    """Ingest submission summaries chronologically with the shared ledger."""
    root.mkdir(parents=True, exist_ok=True)
    latest_state = root / "state" / "scv_latest.parquet"
    submission_releases = sorted(
        (
            release
            for release in releases
            if release.kind is ReleaseKind.SUBMISSION_SUMMARY
        ),
        key=lambda release: (release.release_date, release.url),
    )
    processed = 0
    skipped = 0

    with duckdb.connect(str(root / "clinvar.duckdb")) as connection:
        _create_ledger(connection)
        completed = _completed_urls(connection, ReleaseKind.SUBMISSION_SUMMARY)
        _validate_completed_checksums(
            connection,
            submission_releases,
            ReleaseKind.SUBMISSION_SUMMARY,
        )
        _validate_completed_prefix(submission_releases, completed)
        _validate_completed_identity(
            connection, submission_releases, ReleaseKind.SUBMISSION_SUMMARY
        )
        completed_releases = [
            release for release in submission_releases if release.url in completed
        ]
        if completed_releases:
            last_completed = max(
                completed_releases,
                key=lambda release: (release.release_date, release.url),
            )
            checkpoint = _checkpoint_path(
                root,
                ReleaseKind.SUBMISSION_SUMMARY,
                last_completed,
            )
            previous = read_scv_state(
                checkpoint if checkpoint.exists() else latest_state
            )
        else:
            previous = []

        for release in submission_releases:
            if release.url in completed:
                skipped += 1
                continue

            _mark_started(connection, release)
            try:
                with source(release) as lines:
                    current = consolidate_scv_records(parse_submissions(lines))
                source_md5 = getattr(lines, "calculated_md5", None) or release.md5
                events = list(
                    diff_scv_states(
                        previous,
                        current,
                        observed_at=release.release_date,
                    )
                )
                checkpoint = _checkpoint_path(
                    root,
                    ReleaseKind.SUBMISSION_SUMMARY,
                    release,
                )
                write_scv_events(events, _scv_event_path(root, release))
                write_scv_state(current, checkpoint)
                _publish_latest(checkpoint, latest_state)
                _mark_completed(connection, release, source_md5)
            except Exception as error:
                _mark_failed(connection, release, error)
                raise

            completed.add(release.url)
            previous = current
            processed += 1
            _prune_checkpoints(
                root,
                ReleaseKind.SUBMISSION_SUMMARY,
                checkpoint,
            )

    return IngestionSummary(processed=processed, skipped=skipped)
