import pytest

np = pytest.importorskip("numpy")
pytest.importorskip("scipy")
pytest.importorskip("statsmodels")
from clinvar_history.research_statistics import (
    compress_signatures,
    curves_from_counts,
    estimate_counts,
)


def test_exact_tied_cumulative_incidence_matches_independent_statsmodels():
    from statsmodels.duration.survfunc import CumIncidenceRight

    rng = np.random.default_rng(9)
    times = rng.integers(0, 30, 500)
    status = rng.integers(0, 5, 500)
    codes = times * 2 * 5 + status
    unique, counts = np.unique(codes, return_counts=True)
    ours, _, _, _ = curves_from_counts(unique, counts)
    independent = CumIncidenceRight(times, status, freq_weights=np.ones(len(times)))
    for k in range(4):
        assert np.allclose(
            ours[np.asarray(independent.times, dtype=int) * 2, k, 0],
            independent.cinc[k],
            rtol=0,
            atol=1e-12,
        )


def test_signature_bootstrap_preserves_whole_scv_variant_clusters():
    # Two variant clusters, two SCVs each: resampling clusters yields only
    # proportions 0, 1/2 or 1, never the quarter fractions of an SCV bootstrap.
    data = compress_signatures([[11, 11], [12, 12]], [1, 1])
    result, _, draws = estimate_counts(data, [1], replicates=1000, seed=3)
    assert data.clusters == 2
    assert result[0]["episodes"] == 4
    assert result[0]["benign"] == 0.5
    assert set(draws[:, 0, 0]) <= {0, 0.5, 1}
