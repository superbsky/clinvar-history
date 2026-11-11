"""Content identities for managed manifests, datasets and result checkpoints."""

import hashlib
import json
import os
import sys
import tempfile
from dataclasses import asdict
from pathlib import Path

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from .manifest import ReleaseFile, ReleaseKind
from .storage import _atomic_write_table

# Bump when persisted semantics change, including changes not visible in a schema.
VERSIONS = {"schema": 1, "normalization": 1, "ingestion": 1, "analysis": 1}


def processing_versions() -> dict:
    package = Path(__file__).parent
    sources = (
        "classification.py",
        "vcf_parser.py",
        "submission_parser.py",
        "diff.py",
        "scv_diff.py",
        "storage.py",
        "pipeline.py",
    )
    return {
        **{k: v for k, v in VERSIONS.items() if k != "analysis"},
        "implementation": digest(
            {name: file_digest(package / name) for name in sources}
        ),
    }


def canonical(value) -> str:
    return json.dumps(value, sort_keys=True, separators=(",", ":"), default=str)


def digest(value) -> str:
    return hashlib.sha256(canonical(value).encode()).hexdigest()


def file_digest(path: Path) -> str:
    with path.open("rb") as stream:
        return hashlib.file_digest(stream, "sha256").hexdigest()


def release_identity(release: ReleaseFile) -> str:
    return digest(asdict(release))


def manifest_identity(releases: list[ReleaseFile]) -> str:
    validate_manifest(releases)
    return digest(sorted(release_identity(release) for release in releases))


def validate_manifest(releases) -> None:
    if not releases:
        raise ValueError("Manifest is empty")
    keys = [(r.kind, r.release_date) for r in releases]
    if len(keys) != len(set(keys)) or len({r.url for r in releases}) != len(releases):
        raise ValueError("Manifest has duplicate source dates or URLs")


def metadata_path(path: Path) -> Path:
    return path.with_suffix(path.suffix + ".checkpoint.json")


def atomic_json(value, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", dir=output.parent, delete=False
        ) as f:
            temporary = Path(f.name)
            f.write(canonical(value) + "\n")
        os.replace(temporary, output)
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


def stamp_checkpoint(path: Path, identity: dict) -> None:
    # Decode every batch before accepting an output, including all row groups.
    for _ in pq.ParquetFile(path).iter_batches(batch_size=65536):
        pass
    atomic_json(
        {"format": 1, "identity": identity, "sha256": file_digest(path)},
        metadata_path(path),
    )


def checkpoint_ready(path: Path, identity: dict) -> bool:
    try:
        record = json.loads(metadata_path(path).read_text())
        return (
            record["format"] == 1
            and canonical(record["identity"]) == canonical(identity)
            and record["sha256"] == file_digest(path)
            and pq.ParquetFile(path).metadata is not None
        )
    except (OSError, ValueError, KeyError, TypeError, pa.ArrowException):
        return False


def publish_checkpoint(table: pa.Table, output: Path, identity: dict) -> None:
    _atomic_write_table(table, output)
    stamp_checkpoint(output, identity)


def implementation_identity() -> str:
    """Conservatively invalidate analysis caches after source or lockfile edits."""
    package = Path(__file__).parent
    project = package.parents[1]
    paths = [*package.glob("*.py"), *(project / "analysis").glob("*.py")]
    paths.extend(p for p in [project / "uv.lock"] if p.exists())
    return digest(
        {
            "sources": {
                str(p.relative_to(project)): file_digest(p) for p in sorted(paths)
            },
            "runtime": {
                "python": sys.version,
                "duckdb": duckdb.__version__,
                "pyarrow": pa.__version__,
            },
        }
    )


def dataset_snapshot(root: Path, releases: list[ReleaseFile]) -> dict:
    """Require completed matching history and fingerprint actual analysis inputs."""
    from .pipeline import _validate_completed_checksums, _validate_completed_identity

    validate_manifest(releases)
    with duckdb.connect(str(root / "clinvar.duckdb"), read_only=True) as connection:
        rows = connection.execute(
            "SELECT kind, url, release_date, status, source_md5 FROM pipeline_runs"
        ).fetchall()
        if {(r.kind.value, r.url, r.release_date) for r in releases} != {
            row[:3] for row in rows
        } or any(row[3] != "completed" for row in rows):
            raise ValueError("Dataset ledger must exactly match the completed manifest")
        for kind in ReleaseKind:
            subset = [r for r in releases if r.kind is kind]
            _validate_completed_identity(connection, subset, kind)
            _validate_completed_checksums(connection, subset, kind)
    expected = {
        root
        / "events"
        / ("vcv" if r.kind is ReleaseKind.VCF else "scv")
        / f"release_date={r.release_date}"
        / "events.parquet"
        for r in releases
    }
    actual = set((root / "events").rglob("*.parquet"))
    if expected != actual:
        raise ValueError("Dataset event partitions do not match the manifest")
    return {
        "manifest": manifest_identity(releases),
        "versions": VERSIONS,
        "ledger": digest(sorted(rows)),
        "events": digest(
            {str(p.relative_to(root)): file_digest(p) for p in sorted(actual)}
        ),
    }
