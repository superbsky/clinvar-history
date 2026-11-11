"""Sensitivity to observing the reversible 12–23 August 2022 VCF state window."""

from datetime import date
from pathlib import Path

import pyarrow as pa

from clinvar_history.checkpoints import atomic_json, file_digest, publish_checkpoint
from clinvar_history.research import reconstruct_episodes
from clinvar_history.research_sampling import audit_sampling
from clinvar_history.research_statistics import estimate_counts, load_cluster_counts


def main():
    root = Path("results/full-research")
    output = root / "reversible-batch"
    output.mkdir(parents=True, exist_ok=True)
    excluded = (date(2022, 8, 12), date(2022, 8, 16))
    source = output / "vcv_family_episodes.parquet"
    print(
        "Reconstructing history without observing either transient-window release",
        flush=True,
    )
    reconstruct_episodes(
        Path("data/research"), source, layer="vcv", excluded_dates=excluded
    )
    data = load_cluster_counts(source)
    rows, curve, _ = estimate_counts(
        data, [365, 730, 1095, 1825], replicates=1000, seed=20261111
    )
    identity = {
        "stage": "reversible-window-sensitivity",
        "source_sha256": file_digest(source),
        "excluded_dates": [str(d) for d in excluded],
        "implementation": {
            str(p): file_digest(p)
            for p in [
                Path(__file__),
                Path("src/clinvar_history/research_sampling.py"),
                Path("src/clinvar_history/research_statistics.py"),
                Path("src/clinvar_history/research.py"),
            ]
        },
    }
    publish_checkpoint(
        pa.Table.from_pylist(rows), output / "incidence.parquet", identity
    )
    publish_checkpoint(
        pa.table(curve), output / "curve.parquet", {**identity, "artifact": "curve"}
    )
    sampling = []
    for policy in ["strict", "family"]:
        for name, months, days in [
            ("native", 0, 0),
            ("weekly", 0, 7),
            ("monthly", 1, 0),
            ("quarterly", 3, 0),
        ]:
            print(policy, name, flush=True)
            sampling.append(
                audit_sampling(
                    Path("data/research"),
                    output / f"{policy}_{name}.parquet",
                    policy=policy,
                    months=months,
                    days=days,
                    excluded_dates=excluded,
                )
            )
    atomic_json(sampling, output / "sampling.json")


if __name__ == "__main__":
    main()
