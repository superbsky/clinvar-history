import json
import sys
from contextlib import contextmanager
from datetime import date
from pathlib import Path

import pytest
from click.testing import CliRunner

from clinvar_history.checkpoints import metadata_path, stamp_checkpoint
from clinvar_history.cli import main
from clinvar_history.manifest import ReleaseFile, ReleaseKind, write_manifest
from clinvar_history.pipeline import ingest_submission_releases, ingest_vcf_releases
from clinvar_history.runner import _run, manifest_request

FIXTURES = Path(__file__).parent / "fixtures"


def _dataset(root):
    releases = [
        ReleaseFile(date(2024, 1, day), kind, f"fixture://{name}", 1)
        for day, kind, name in [
            (1, ReleaseKind.VCF, "release_a.vcf"),
            (8, ReleaseKind.VCF, "release_b.vcf"),
            (1, ReleaseKind.SUBMISSION_SUMMARY, "submission_a.txt"),
            (8, ReleaseKind.SUBMISSION_SUMMARY, "submission_b.txt"),
        ]
    ]

    @contextmanager
    def source(release):
        with (FIXTURES / release.url.removeprefix("fixture://")).open() as stream:
            yield stream

    ingest_vcf_releases(releases, root=root, source=source)
    ingest_submission_releases(releases, root=root, source=source)
    manifest = root / "manifest/releases.parquet"
    write_manifest(releases, manifest)
    stamp_checkpoint(manifest, manifest_request(date(2024, 1, 1), date(2024, 1, 8)))


def test_offline_runner_computes_then_reuses_and_invalidates_parameters(tmp_path):
    root, results = tmp_path / "data", tmp_path / "results"
    _dataset(root)
    args = [
        "run-study",
        "--root",
        str(root),
        "--results",
        str(results),
        "--start",
        "2024-01-01",
        "--end",
        "2024-01-08",
        "--offline",
    ]
    runner = CliRunner()
    first = runner.invoke(main, args)
    assert first.exit_code == 0, first.output
    assert first.output.count("Computing") == 7
    import pyarrow.parquet as pq

    assert {"scv_accession", "started_at", "outcome"} <= set(
        pq.read_schema(results / "scv_vus_time_to_resolution.parquet").names
    )
    before = {p: p.stat().st_mtime_ns for p in results.glob("*.parquet")}
    second = runner.invoke(main, args)
    assert second.exit_code == 0, second.output
    assert second.output.count("Reusing verified") == 7
    assert before == {p: p.stat().st_mtime_ns for p in results.glob("*.parquet")}
    third = runner.invoke(main, [*args, "--max-lag-days", "7"])
    assert third.exit_code == 0, third.output
    assert "Computing attributions" in third.output
    assert "Reusing verified transition_matrix" in third.output
    changed_dates = runner.invoke(main, [*args, "--end", "2024-01-09"])
    assert changed_dates.exit_code != 0
    assert "incompatible" in changed_dates.output
    path = metadata_path(results / "transition_matrix.parquet")
    record = json.loads(path.read_text())
    record["identity"]["implementation"] = "old"
    path.write_text(json.dumps(record))
    fourth = runner.invoke(main, args)
    assert fourth.exit_code == 0, fourth.output
    assert "Computing transition_matrix" in fourth.output


def test_failed_stage_has_bounded_retries(tmp_path):
    counter = tmp_path / "attempts"
    code = "from pathlib import Path; p=Path(__import__('sys').argv[1]); p.write_text(p.read_text()+'x' if p.exists() else 'x'); raise SystemExit(7)"
    with pytest.raises(Exception, match="after 3 attempt"):
        _run([sys.executable, "-c", code, str(counter)], attempts=3)
    assert counter.read_text() == "xxx"


def test_interrupted_recalculation_preserves_published_result(tmp_path, monkeypatch):
    import click

    import clinvar_history.runner as managed

    root, results = tmp_path / "data", tmp_path / "results"
    _dataset(root)
    args = [
        "run-study",
        "--root",
        str(root),
        "--results",
        str(results),
        "--start",
        "2024-01-01",
        "--end",
        "2024-01-08",
        "--offline",
    ]
    runner = CliRunner()
    assert runner.invoke(main, args).exit_code == 0
    output = results / "transition_matrix.parquet"
    before = output.read_bytes()
    metadata_path(output).unlink()

    def failing_stage(command, *args, **kwargs):
        pending = Path(command[command.index("--output") + 1])
        pending.write_bytes(b"partial output")
        raise click.ClickException("simulated interruption")

    monkeypatch.setattr(managed, "_run", failing_stage)
    result = runner.invoke(main, args)
    assert result.exit_code != 0
    assert output.read_bytes() == before
    assert not metadata_path(output).exists()


def test_adoption_cannot_relabel_an_already_versioned_dataset(tmp_path):
    root, results = tmp_path / "data", tmp_path / "results"
    _dataset(root)
    results.mkdir()
    result = CliRunner().invoke(
        main,
        [
            "adopt-legacy",
            "--root",
            str(root),
            "--results",
            str(results),
            "--start",
            "2024-01-01",
            "--end",
            "2024-01-08",
        ],
    )
    assert result.exit_code != 0
    assert "cannot overwrite" in result.output
