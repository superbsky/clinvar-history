"""Scalable exact cluster-signature bootstrap for full-cohort cumulative incidence."""

from dataclasses import dataclass
from pathlib import Path

import numpy as np
from scipy import sparse

from .research import research_connection

CAUSES = ("benign", "pathogenic", "conflict", "other_exit")


@dataclass
class ClusterCounts:
    codes: np.ndarray
    matrix: sparse.csr_matrix
    frequencies: np.ndarray

    @property
    def clusters(self):
        return int(self.frequencies.sum())


def compress_signatures(signatures, frequencies):
    """Identical cluster outcome vectors are exchangeable in an exact bootstrap."""
    codes = np.asarray(
        sorted({v for signature in signatures for v in signature}), dtype=np.int64
    )
    if codes.size == 0:
        raise ValueError("No observations")
    positions = {int(v): i for i, v in enumerate(codes)}
    rows = []
    cols = []
    for j, signature in enumerate(signatures):
        rows.extend(positions[int(v)] for v in signature)
        cols.extend([j] * len(signature))
    matrix = sparse.coo_matrix(
        (np.ones(len(rows)), (rows, cols)), shape=(len(codes), len(signatures))
    ).tocsr()
    return ClusterCounts(codes, matrix, np.asarray(frequencies, dtype=np.int64))


def load_cluster_counts(
    path: Path,
    *,
    where="incident_rank=1 AND NOT identity_migrated",
    parameters=(),
    timing="observed",
    removal_as_exit=False,
):
    if timing not in {"observed", "midpoint", "lower", "upper"}:
        raise ValueError("Invalid timing assumption")
    duration = {
        "observed": "2*duration_days",
        "midpoint": "2*duration_days + date_diff('day',entry_lower,started_at) - CASE WHEN outcome!='censored' THEN date_diff('day',exit_lower,ended_at) ELSE 0 END",
        "lower": "CASE WHEN outcome='censored' THEN 2*duration_days ELSE 2*duration_lower_days END",
        "upper": "2*duration_upper_days",
    }[timing]
    cause = (
        "CASE WHEN censor_reason='removal' THEN 'other_exit' ELSE outcome END"
        if removal_as_exit
        else "outcome"
    )
    with research_connection() as c:
        query = f"""WITH observations AS (
            SELECT variation_id, greatest(0,({duration}))::BIGINT*5 +
                CASE ({cause}) WHEN 'censored' THEN 0 WHEN 'benign' THEN 1 WHEN 'pathogenic' THEN 2
                               WHEN 'conflict' THEN 3 WHEN 'other_exit' THEN 4 ELSE NULL END AS code
            FROM read_parquet(?) WHERE {where}
        ), clusters AS (
            SELECT variation_id, list(code ORDER BY code) AS signature FROM observations GROUP BY variation_id
        ) SELECT signature,count(*) AS frequency FROM clusters GROUP BY signature ORDER BY signature"""
        rows = c.execute(query, [str(path), *parameters]).fetchall()
    return compress_signatures([r[0] for r in rows], [r[1] for r in rows])


def curves_from_counts(codes, counts):
    """Vectorized tied-time Aalen–Johansen for one or a batch of count vectors."""
    counts = np.asarray(counts, dtype=float)
    if counts.ndim == 1:
        counts = counts[:, None]
    size = (int(np.max(codes)) // 5 + 1) * 5
    dense = np.zeros((size, counts.shape[1]))
    dense[codes] = counts
    cells = dense.reshape(-1, 5, counts.shape[1])
    exited = cells.sum(axis=1)
    at_risk = exited.sum(axis=0)[None, :] - np.vstack(
        [np.zeros((1, counts.shape[1])), np.cumsum(exited, axis=0)[:-1]]
    )
    all_events = cells[:, 1:, :].sum(axis=1)
    hazard = np.divide(
        all_events, at_risk, out=np.zeros_like(all_events), where=at_risk > 0
    )
    survival = np.cumprod(1 - hazard, axis=0)
    left = np.vstack([np.ones((1, counts.shape[1])), survival[:-1]])
    hazards = np.divide(
        cells[:, 1:, :],
        at_risk[:, None, :],
        out=np.zeros_like(cells[:, 1:, :]),
        where=at_risk[:, None, :] > 0,
    )
    incidence = np.cumsum(left[:, None, :] * hazards, axis=0)
    max_time = np.max(np.where(exited > 0, np.arange(len(exited))[:, None], -1), axis=0)
    return incidence, survival, at_risk, max_time


def estimate_counts(data: ClusterCounts, horizons, *, replicates=1000, seed=20261111):
    if replicates < 100 or any(h < 0 for h in horizons):
        raise ValueError("Invalid bootstrap or horizon")
    counts = np.asarray(data.matrix @ data.frequencies).ravel()
    incidence, survival, at_risk, max_time = curves_from_counts(data.codes, counts)
    n = int(counts.sum())
    horizon_ticks = np.asarray(horizons, dtype=int) * 2
    indices = np.minimum(horizon_ticks, len(incidence) - 1)
    points = incidence[indices, :, 0]
    supported = (horizon_ticks <= max_time[0]) | (survival[-1, 0] < 1e-12)
    probabilities = data.frequencies / data.clusters
    rng = np.random.default_rng(seed)
    draws = np.full((replicates, len(horizons), 4), np.nan)
    batch = max(1, min(50, 16_000_000 // len(data.frequencies)))
    for first in range(0, replicates, batch):
        amount = min(batch, replicates - first)
        weights = rng.multinomial(data.clusters, probabilities, size=amount).T
        result = np.asarray(data.matrix @ weights)
        ci, s, _, last = curves_from_counts(data.codes, result)
        keep = (horizon_ticks[:, None] <= last[None, :]) | (s[-1][None, :] < 1e-12)
        values = np.transpose(ci[indices], (2, 0, 1))
        values[~keep.T] = np.nan
        draws[first : first + amount] = values
    result = []
    for i, horizon in enumerate(horizons):
        row = {
            "day": int(horizon),
            "episodes": n,
            "clusters": data.clusters,
            "at_risk": int(at_risk[indices[i], 0])
            if horizon_ticks[i] < len(at_risk)
            else 0,
            "supported": bool(supported[i]),
            "bootstrap_replicates": replicates,
        }
        for k, cause in enumerate(CAUSES):
            values = draws[:, i, k]
            values = values[np.isfinite(values)]
            enough = supported[i] and len(values) >= 0.95 * replicates
            lo, hi = np.quantile(values, [0.025, 0.975]) if enough else (None, None)
            row[cause] = float(points[i, k]) if supported[i] else None
            row[cause + "_lower"] = float(lo) if lo is not None else None
            row[cause + "_upper"] = float(hi) if hi is not None else None
            row[cause + "_bootstrap_supported"] = len(values)
            half = draws[: replicates // 2, i, k]
            half = half[np.isfinite(half)]
            row[cause + "_ci_half_run_max_change"] = (
                float(np.max(np.abs(np.quantile(half, [0.025, 0.975]) - [lo, hi])))
                if enough and len(half) >= 0.475 * replicates
                else None
            )
        result.append(row)
    curve = {
        "day": np.arange(len(incidence)) / 2,
        "survival": survival[:, 0],
        "at_risk": at_risk[:, 0],
    }
    curve.update({cause: incidence[:, k, 0] for k, cause in enumerate(CAUSES)})
    return result, curve, draws
