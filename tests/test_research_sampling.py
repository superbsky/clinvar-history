from datetime import date

import pyarrow as pa
import pyarrow.parquet as pq
import pytest

from clinvar_history.manifest import ReleaseFile, ReleaseKind, write_manifest
from clinvar_history.research_sampling import audit_sampling


def test_net_endpoint_edge_need_not_recover_an_observed_edge(tmp_path):
    dates = [date(2020, 1, d) for d in [1, 8, 15, 22]]
    write_manifest(
        [ReleaseFile(d, ReleaseKind.VCF, f"fixture://{d}", 1) for d in dates],
        tmp_path / "manifest/releases.parquet",
    )
    labels = ["benign", "vus", "pathogenic"]
    rows = [
        {
            "variation_id": 1,
            "observed_at": d,
            "event_type": "variant_added" if i == 0 else "classification_changed",
            "new_classification": labels[i],
            "new_raw_classification": labels[i],
        }
        for i, d in enumerate(dates[:3])
    ]
    path = tmp_path / "events/vcv/release_date=2020-01-01/events.parquet"
    path.parent.mkdir(parents=True)
    pq.write_table(pa.Table.from_pylist(rows), path)
    native = audit_sampling(tmp_path, tmp_path / "native.parquet")
    assert native["edge_recall"] == 1
    coarse = audit_sampling(tmp_path, tmp_path / "monthly.parquet", months=1)
    assert coarse["reference_edges"] == 2
    assert coarse["count_recovery"] == 0.5
    assert coarse["edge_recall"] == 0
    assert coarse["synthetic_edges"] == 1
    omitted = audit_sampling(
        tmp_path, tmp_path / "omit.parquet", excluded_dates=(dates[1],)
    )
    assert omitted["reference_edges"] == 1
    assert omitted["edge_recall"] == 1
    with pytest.raises(ValueError, match="boundary"):
        audit_sampling(tmp_path, tmp_path / "bad.parquet", excluded_dates=(dates[0],))
