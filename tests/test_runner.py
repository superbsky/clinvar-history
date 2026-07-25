import os
import subprocess
from pathlib import Path


def test_full_runner_reuses_manifest_and_retries_ingestion(tmp_path):
    project_root = Path(__file__).parents[1]
    data_root = tmp_path / "data"
    results_root = tmp_path / "results"
    manifest = data_root / "manifest" / "releases.parquet"
    manifest.parent.mkdir(parents=True)
    manifest.touch()

    calls = tmp_path / "calls.log"
    counter = tmp_path / "counter"
    fake_uv = tmp_path / "fake-uv"
    fake_uv.write_text(
        """#!/usr/bin/env bash
set -eu
printf '%s\n' "$*" >> "${FAKE_UV_CALLS}"
if [[ "$*" == *"ingest-vcf"* ]] && [[ ! -e "${FAKE_UV_COUNTER}" ]]; then
  touch "${FAKE_UV_COUNTER}"
  exit 1
fi
"""
    )
    fake_uv.chmod(0o755)

    result = subprocess.run(
        [project_root / "scripts" / "run_full_backfill.sh"],
        cwd=project_root,
        env={
            **os.environ,
            "CLINVAR_DATA_ROOT": str(data_root),
            "CLINVAR_RESULTS_ROOT": str(results_root),
            "CLINVAR_UV_BIN": str(fake_uv),
            "CLINVAR_RETRY_DELAY_SECONDS": "0",
            "FAKE_UV_CALLS": str(calls),
            "FAKE_UV_COUNTER": str(counter),
        },
        capture_output=True,
        check=False,
        text=True,
        timeout=10,
    )

    invoked = calls.read_text()
    assert result.returncode == 0, result.stderr
    assert "build-manifest" not in invoked
    assert invoked.count("ingest-vcf") == 2
    assert invoked.count("ingest-submissions") == 1


def test_full_runner_reuses_valid_existing_analysis_outputs(tmp_path):
    project_root = Path(__file__).parents[1]
    data_root = tmp_path / "data"
    results_root = tmp_path / "results"
    manifest = data_root / "manifest" / "releases.parquet"
    manifest.parent.mkdir(parents=True)
    manifest.touch()
    results_root.mkdir()
    for name in (
        "attributions.parquet",
        "renovo_quarterly.parquet",
        "sampling_bias.parquet",
        "transition_matrix.parquet",
        "vus_time_to_resolution.parquet",
        "event_mechanisms.parquet",
        "scv_vus_time_to_resolution.parquet",
    ):
        (results_root / name).write_text("existing")

    calls = tmp_path / "calls.log"
    fake_uv = tmp_path / "fake-uv"
    fake_uv.write_text(
        """#!/usr/bin/env bash
set -eu
printf '%s\n' "$*" >> "${FAKE_UV_CALLS}"
"""
    )
    fake_uv.chmod(0o755)

    result = subprocess.run(
        [project_root / "scripts" / "run_full_backfill.sh"],
        cwd=project_root,
        env={
            **os.environ,
            "CLINVAR_DATA_ROOT": str(data_root),
            "CLINVAR_RESULTS_ROOT": str(results_root),
            "CLINVAR_UV_BIN": str(fake_uv),
            "FAKE_UV_CALLS": str(calls),
        },
        capture_output=True,
        check=False,
        text=True,
        timeout=10,
    )

    invoked = calls.read_text()
    assert result.returncode == 0, result.stderr
    assert "attribute-events" not in invoked
    assert "analysis/01_reproduce_renovo.py" not in invoked
    assert "analysis/02_sampling_bias.py" not in invoked
    assert "analysis/03_transition_matrix.py" not in invoked
    assert "analysis/04_time_to_resolution.py" not in invoked
    assert "analysis/05_event_mechanisms.py" not in invoked
    assert "analysis/06_scv_time_to_resolution.py" not in invoked
