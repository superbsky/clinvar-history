from datetime import date

import pyarrow as pa
import pyarrow.parquet as pq

from clinvar_history.checkpoints import (
    checkpoint_ready,
    manifest_identity,
    publish_checkpoint,
)
from clinvar_history.manifest import ReleaseFile, ReleaseKind


def test_manifest_identity_ignores_order_but_tracks_source_and_dates():
    a = ReleaseFile(date(2024, 1, 1), ReleaseKind.VCF, "https://a", 10)
    b = ReleaseFile(date(2024, 1, 8), ReleaseKind.VCF, "https://b", 11)
    assert manifest_identity([a, b]) == manifest_identity([b, a])
    assert manifest_identity([a]) != manifest_identity([b])


def test_readable_legacy_output_is_not_a_checkpoint(tmp_path):
    path = tmp_path / "output.parquet"
    pq.write_table(pa.table({"value": [1]}), path)
    assert not checkpoint_ready(path, {"study_end": "2024-01-01"})


def test_checkpoint_tracks_parameters_versions_and_output_integrity(tmp_path):
    path = tmp_path / "output.parquet"
    identity = {"study_end": "2024-01-01", "normalization": 1, "lag": 35}
    publish_checkpoint(pa.table({"value": [1]}), path, identity)
    assert checkpoint_ready(path, identity)
    for key, value in [("study_end", "2024-02-01"), ("normalization", 2), ("lag", 7)]:
        assert not checkpoint_ready(path, {**identity, key: value})
    pq.write_table(pa.table({"value": [2]}), path)
    assert not checkpoint_ready(path, identity)


def test_interrupted_or_corrupt_metadata_is_not_reused(tmp_path):
    path = tmp_path / "output.parquet"
    publish_checkpoint(pa.table({"value": [1]}), path, {"stage": "test"})
    path.with_suffix(".parquet.checkpoint.json").write_text("{")
    assert not checkpoint_ready(path, {"stage": "test"})
