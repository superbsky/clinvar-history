"""Managed full-study runner with semantic, content-verified checkpoints."""

import os
import shutil
import subprocess
import sys
import tempfile
import time
from datetime import date
from pathlib import Path

import click
import duckdb

from .checkpoints import (
    canonical,
    checkpoint_ready,
    dataset_snapshot,
    file_digest,
    implementation_identity,
    metadata_path,
    processing_versions,
    release_identity,
    stamp_checkpoint,
    validate_manifest,
)
from .manifest import ReleaseKind, read_manifest


def manifest_request(start: date, end: date) -> dict:
    return {
        "stage": "manifest",
        "start": start.isoformat(),
        "end": end.isoformat(),
        "discovery": file_digest(Path(__file__).with_name("manifest.py")),
        "format": 1,
    }


def _run(command: list[str], attempts: int = 1, retry_delay: float = 0) -> None:
    for attempt in range(1, attempts + 1):
        result = subprocess.run(command, check=False)
        if result.returncode == 0:
            return
        if attempt == attempts:
            raise click.ClickException(
                f"Stage failed after {attempts} attempt(s): {command}"
            )
        click.echo(
            f"Stage failed ({result.returncode}); retry {attempt + 1}/{attempts}"
        )
        time.sleep(retry_delay)


def _validate_extension(root: Path, releases) -> None:
    from .pipeline import (
        _completed_urls,
        _validate_completed_checksums,
        _validate_completed_identity,
        _validate_completed_prefix,
    )

    if not (root / "clinvar.duckdb").exists():
        return
    with duckdb.connect(str(root / "clinvar.duckdb"), read_only=True) as connection:
        for kind in ReleaseKind:
            subset = sorted(
                (r for r in releases if r.kind is kind),
                key=lambda r: (r.release_date, r.url),
            )
            _validate_completed_prefix(subset, _completed_urls(connection, kind))
            _validate_completed_checksums(connection, subset, kind)
            _validate_completed_identity(connection, subset, kind)


def run_study(
    root: Path,
    results: Path,
    start: date,
    end: date,
    *,
    offline: bool,
    refresh_manifest: bool,
    max_lag_days: int,
    attempts: int,
    retry_delay: float,
) -> None:
    if start >= end:
        raise ValueError("Study start must precede end")
    root.mkdir(parents=True, exist_ok=True)
    results.mkdir(parents=True, exist_ok=True)
    manifest = root / "manifest" / "releases.parquet"
    request = manifest_request(start, end)
    if refresh_manifest or not checkpoint_ready(manifest, request):
        if offline:
            raise ValueError(
                "Offline manifest is incompatible or unverified; build-manifest or adopt-legacy first"
            )
        manifest.parent.mkdir(parents=True, exist_ok=True)
        temporary = manifest.with_name(".releases.pending.parquet")
        try:
            _run(
                [
                    sys.executable,
                    "-m",
                    "clinvar_history.cli",
                    "build-manifest",
                    "--start",
                    str(start),
                    "--end",
                    str(end),
                    "--output",
                    str(temporary),
                ],
                attempts,
                retry_delay,
            )
            _validate_extension(root, read_manifest(temporary))
            os.replace(temporary, manifest)
            os.replace(metadata_path(temporary), metadata_path(manifest))
        finally:
            temporary.unlink(missing_ok=True)
            metadata_path(temporary).unlink(missing_ok=True)
    releases = read_manifest(manifest)
    validate_manifest(releases)
    if sum(r.kind is ReleaseKind.VCF for r in releases) < 2 or not any(
        r.kind is ReleaseKind.SUBMISSION_SUMMARY for r in releases
    ):
        raise ValueError("Study requires at least two VCF releases and an SCV release")
    _validate_extension(root, releases)
    if not offline:
        # Sources stream through memory. Reserve a conservative compressed-size
        # multiplier plus two latest-state copies for atomic replacements.
        state_bytes = sum(
            p.stat().st_size for p in (root / "state").glob("*_latest.parquet")
        )
        required = max(
            512 * 1024**2, 4 * max(r.size_bytes for r in releases) + 2 * state_bytes
        )
        free = shutil.disk_usage(root).free
        if free < required:
            raise ValueError(
                f"Insufficient free disk: {free} bytes available, {required} required"
            )
        click.echo(f"Disk preflight: {free} bytes free; reserve {required}")
        for command in ("ingest-vcf", "ingest-submissions"):
            _run(
                [
                    sys.executable,
                    "-m",
                    "clinvar_history.cli",
                    command,
                    "--manifest",
                    str(manifest),
                    "--root",
                    str(root),
                ],
                attempts,
                retry_delay,
            )
    with duckdb.connect(str(root / "clinvar.duckdb"), read_only=True):
        snapshot = dataset_snapshot(root, releases)
        vcf_dates = [r.release_date for r in releases if r.kind is ReleaseKind.VCF]
        if len(vcf_dates) < 2 or not any(
            r.kind is ReleaseKind.SUBMISSION_SUMMARY for r in releases
        ):
            raise ValueError(
                "Study requires at least two VCF releases and an SCV release"
            )
        analysis_start = min(vcf_dates)
        implementation = implementation_identity()
        common = {
            "dataset": snapshot,
            "start": str(start),
            "end": str(end),
            "analysis_start": str(analysis_start),
            "implementation": implementation,
        }
        stages = [
            (
                "attributions",
                ["-m", "clinvar_history.cli", "attribute-events"],
                ["--max-lag-days", str(max_lag_days)],
            ),
            (
                "renovo_quarterly",
                ["-m", "clinvar_history.stages", "--stage", "renovo_quarterly"],
                ["--start", str(analysis_start), "--end", str(end)],
            ),
            (
                "sampling_bias",
                ["-m", "clinvar_history.stages", "--stage", "sampling_bias"],
                [],
            ),
            (
                "transition_matrix",
                ["-m", "clinvar_history.stages", "--stage", "transition_matrix"],
                [],
            ),
            (
                "vus_time_to_resolution",
                ["-m", "clinvar_history.stages", "--stage", "vus_time_to_resolution"],
                ["--study-end", str(end)],
            ),
            (
                "event_mechanisms",
                ["-m", "clinvar_history.stages", "--stage", "event_mechanisms"],
                ["--attributions", str(results / "attributions.parquet")],
            ),
            (
                "scv_vus_time_to_resolution",
                [
                    "-m",
                    "clinvar_history.stages",
                    "--stage",
                    "scv_vus_time_to_resolution",
                ],
                ["--study-end", str(end)],
            ),
        ]
        for name, command, parameters in stages:
            identity = {**common, "stage": name, "parameters": parameters}
            if name == "event_mechanisms":
                identity["attributions_sha256"] = file_digest(
                    results / "attributions.parquet"
                )
            output = results / f"{name}.parquet"
            if checkpoint_ready(output, identity):
                click.echo(f"Reusing verified {name}")
                continue
            click.echo(f"Computing {name}: checkpoint missing or incompatible")
            with tempfile.TemporaryDirectory(
                prefix=f".{name}.", dir=results
            ) as scratch:
                pending = Path(scratch) / "result.parquet"
                _run(
                    [
                        sys.executable,
                        *command,
                        "--root",
                        str(root),
                        "--output",
                        str(pending),
                        *parameters,
                    ]
                )
                if implementation_identity() != implementation:
                    raise ValueError(
                        "Analysis source changed during calculation; rerun the study"
                    )
                stamp_checkpoint(pending, identity)
                os.replace(pending, output)
                os.replace(metadata_path(pending), metadata_path(output))
        click.echo("Study completed with verified checkpoints")


def adopt_legacy(root: Path, results: Path, start: date, end: date) -> None:
    """Explicitly register audited v0.1 data; never stamp legacy analysis outputs."""
    from .pipeline import _create_ledger
    from .reporting import audit_full_backfill

    manifest = root / "manifest" / "releases.parquet"
    releases = read_manifest(manifest)
    with duckdb.connect(str(root / "clinvar.duckdb"), read_only=True) as connection:
        columns = {
            r[0] for r in connection.execute("DESCRIBE pipeline_runs").fetchall()
        }
        if {"release_identity", "processing_versions"} <= columns:
            existing = connection.execute(
                "SELECT release_identity, processing_versions FROM pipeline_runs"
            ).fetchall()
            if any(
                identity is not None or versions is not None
                for identity, versions in existing
            ):
                _validate_extension(root, releases)
                raise ValueError(
                    "Dataset already has provenance; adoption cannot overwrite it"
                )
    if start > min(r.release_date for r in releases) or end < max(
        r.release_date for r in releases
    ):
        raise ValueError("Requested dates exclude retained releases")
    audit = audit_full_backfill(root, results)
    if audit.study_end != end or start > audit.study_start:
        raise ValueError("Requested horizon differs from audited study")
    with duckdb.connect(str(root / "clinvar.duckdb")) as connection:
        _create_ledger(connection)
        connection.begin()
        try:
            for release in releases:
                connection.execute(
                    "UPDATE pipeline_runs SET release_identity = ?, processing_versions = ? WHERE kind = ? AND url = ?",
                    [
                        release_identity(release),
                        canonical(processing_versions()),
                        release.kind.value,
                        release.url,
                    ],
                )
            connection.commit()
        except Exception:
            connection.rollback()
            raise
    stamp_checkpoint(manifest, manifest_request(start, end))
    click.echo(
        f"Adopted {len(releases)} audited v0.1 releases; legacy analysis outputs remain unverified"
    )
