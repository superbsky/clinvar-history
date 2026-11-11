"""Independent SQL reconstruction and audits for full-cohort ClinVar research."""

import bisect
import os
from datetime import date
from pathlib import Path

import duckdb
import pyarrow as pa

from .analysis_core import month_grid
from .manifest import ReleaseKind, read_manifest


def endpoint_sql(raw: str, canonical: str, policy: str) -> str:
    """Analysis-only family bins; preserve the original persisted normalization."""
    if policy == "strict":
        return f"coalesce({canonical}, 'missing')"
    if policy != "family":
        raise ValueError("Endpoint policy must be strict or family")
    return f"""CASE
        WHEN {canonical} IN ('benign', 'likely_benign')
          OR lower(trim({raw})) IN ('benign/likely benign','likely benign/benign') THEN 'benign'
        WHEN {canonical} IN ('pathogenic', 'likely_pathogenic')
          OR lower(trim({raw})) IN ('pathogenic/likely pathogenic','likely pathogenic/pathogenic') THEN 'pathogenic'
        ELSE coalesce({canonical}, 'missing') END"""


def research_connection(scratch: Path | None = None):
    connection = duckdb.connect()
    connection.execute("SET memory_limit='4GB'")
    connection.execute("SET threads=4")
    if scratch is not None:
        scratch.mkdir(parents=True, exist_ok=True)
        connection.execute("SET temp_directory=?", [str(scratch)])
        connection.execute("SET max_temp_directory_size='16GB'")
    return connection


def _copy(connection, query: str, output: Path, parameters=None):
    output.parent.mkdir(parents=True, exist_ok=True)
    pending = output.with_name("." + output.name + ".pending")
    try:
        connection.execute(
            "COPY ("
            + query
            + ") TO '"
            + str(pending).replace("'", "''")
            + "' (FORMAT PARQUET, COMPRESSION ZSTD)",
            parameters or [],
        )
        os.replace(pending, output)
    finally:
        pending.unlink(missing_ok=True)


def reconstruct_episodes(
    root: Path,
    output: Path,
    *,
    layer: str,
    policy="family",
    cadence_months=0,
    administrative_end: date | None = None,
    excluded_dates: tuple[date, ...] = (),
):
    """Reconstruct all incident VUS episodes with lead/lag over observed states.

    All variants/accessions are retained. Coarser observations collapse source
    changes within each sample interval before episodes are identified. An
    administrative_end override is only for comparison with legacy exports.
    """
    if layer not in {"vcv", "scv"} or cadence_months not in {0, 1, 3}:
        raise ValueError("Unsupported layer or cadence")
    kind = ReleaseKind.VCF if layer == "vcv" else ReleaseKind.SUBMISSION_SUMMARY
    releases = sorted(
        {
            r.release_date
            for r in read_manifest(root / "manifest/releases.parquet")
            if r.kind is kind
        }
    )
    if not releases:
        raise ValueError("No source releases")
    cutoff = administrative_end or releases[-1]
    if cutoff < releases[-1]:
        raise ValueError("Administrative override cannot precede source end")
    grid = (
        releases
        if not cadence_months
        else month_grid(releases[0], releases[-1], cadence_months)
    )
    if excluded_dates:
        if (
            cadence_months
            or releases[0] in excluded_dates
            or releases[-1] in excluded_dates
        ):
            raise ValueError("Exclude only internal release dates at native cadence")
        if not set(excluded_dates) <= set(releases):
            raise ValueError("Excluded date is not a source release")
        grid = [d for d in grid if d not in excluded_dates]

    source = str(root / f"events/{layer}/release_date=*/events.parquet")
    key = "variation_id" if layer == "vcv" else "scv_accession"
    added = "variant_added" if layer == "vcv" else "scv_added"
    removed = "variant_removed" if layer == "vcv" else "scv_removed"
    changed = "classification_changed" if layer == "vcv" else "scv_reclassified"
    new = endpoint_sql("new_raw_classification", "new_classification", policy)
    metadata = (
        "new_gene_info AS gene_info, new_variant_type AS variant_type, NULL::VARCHAR AS submitter, NULL::VARCHAR AS collection_method"
        if layer == "vcv"
        else "NULL::VARCHAR AS gene_info, NULL::VARCHAR AS variant_type, new_submitter AS submitter, new_collection_method AS collection_method"
    )
    with research_connection(output.parent / ".scratch") as c:
        c.register(
            "release_map",
            pa.table(
                {
                    "observed_at": pa.array(releases, type=pa.date32()),
                    "sample_at": pa.array(
                        [grid[bisect.bisect_left(grid, d)] for d in releases],
                        type=pa.date32(),
                    ),
                }
            ),
        )
        c.register(
            "sample_grid",
            pa.table(
                {
                    "sample_at": pa.array(grid, type=pa.date32()),
                    "previous_sample": pa.array([None, *grid[:-1]], type=pa.date32()),
                }
            ),
        )
        c.execute(f"""CREATE TEMP VIEW unit_identity AS
            SELECT CAST({key} AS VARCHAR) AS unit, min(observed_at) AS first_seen,
                   count(DISTINCT variation_id)>1 AS identity_migrated
            FROM read_parquet('{source.replace("'", "''")}') GROUP BY {key}""")
        query = f"""
            WITH raw AS (
                SELECT CAST({key} AS VARCHAR) AS unit, variation_id, event_type, observed_at,
                       CASE WHEN event_type='{removed}' THEN 'absent' ELSE {new} END AS state,
                       new_raw_classification AS raw_label, new_review_status AS review_status,
                       new_conditions AS conditions, {metadata}
                FROM read_parquet(?) WHERE event_type IN ('{added}', '{removed}', '{changed}')
            ), sampled AS (
                SELECT raw.*, release_map.sample_at,
                       row_number() OVER (PARTITION BY unit, sample_at ORDER BY raw.observed_at DESC, event_type DESC) AS sample_rank
                FROM raw JOIN release_map USING (observed_at)
            ), lagged AS (
                SELECT *, lag(state) OVER (PARTITION BY unit ORDER BY sample_at) AS previous_state
                FROM sampled WHERE sample_rank=1
            ), changes AS (
                SELECT * FROM lagged WHERE previous_state IS NULL OR state IS DISTINCT FROM previous_state
            ), linked AS (
                SELECT *, lead(sample_at) OVER (PARTITION BY unit ORDER BY sample_at) AS next_at,
                       lead(state) OVER (PARTITION BY unit ORDER BY sample_at) AS next_state,
                       lead(raw_label) OVER (PARTITION BY unit ORDER BY sample_at) AS exit_raw_label
                FROM changes
            ), episodes AS (
                SELECT unit, variation_id, l.sample_at AS started_at, coalesce(next_at, ?::DATE) AS ended_at,
                       g.previous_sample AS entry_lower,
                       CASE WHEN next_at IS NULL THEN ?::DATE ELSE h.previous_sample END AS exit_lower,
                       CASE WHEN next_at IS NULL OR next_state='absent' THEN 'censored'
                            WHEN next_state IN ('benign','likely_benign') THEN 'benign'
                            WHEN next_state IN ('pathogenic','likely_pathogenic') THEN 'pathogenic'
                            WHEN next_state='conflict' THEN 'conflict' ELSE 'other_exit' END AS outcome,
                       CASE WHEN next_at IS NULL THEN 'administrative'
                            WHEN next_state='absent' THEN 'removal' ELSE NULL END AS censor_reason,
                       CASE WHEN previous_state IS NULL THEN 'new_presence'
                            WHEN previous_state='absent' THEN 'reappearance' ELSE 'reclassification' END AS entry_origin,
                       review_status, gene_info, variant_type, submitter, collection_method, conditions,
                       raw_label AS entry_raw_label, exit_raw_label
                FROM linked l JOIN sample_grid g ON l.sample_at=g.sample_at
                LEFT JOIN sample_grid h ON l.next_at=h.sample_at
                WHERE state='vus' AND l.sample_at > ?::DATE
            )
            SELECT e.*, date_diff('day', started_at, ended_at)::INTEGER AS duration_days,
                   greatest(0,date_diff('day', started_at, exit_lower))::INTEGER AS duration_lower_days,
                   date_diff('day', entry_lower, ended_at)::INTEGER AS duration_upper_days,
                   row_number() OVER (PARTITION BY e.unit ORDER BY started_at) AS incident_rank,
                   i.identity_migrated, i.first_seen,
                   '{layer}' AS layer, '{policy}' AS endpoint_policy, {cadence_months} AS cadence_months
            FROM episodes e JOIN unit_identity i USING (unit) ORDER BY unit, started_at
        """
        _copy(c, query, output, [source, cutoff, cutoff, grid[0]])
        return c.execute(
            "SELECT count(*), count(DISTINCT unit) FROM read_parquet(?)", [str(output)]
        ).fetchone()


def compare_legacy(connection, reconstructed: Path, legacy: Path, layer: str):
    key = "variation_id" if layer == "vcv" else "scv_accession"
    query = f"""
        WITH expected AS (
            SELECT CAST({key} AS VARCHAR) AS unit, started_at, ended_at, outcome, duration_days FROM read_parquet(?)
        ), actual AS (
            SELECT unit, started_at, ended_at, outcome, duration_days FROM read_parquet(?)
        )
        SELECT (SELECT count(*) FROM (SELECT * FROM expected EXCEPT ALL SELECT * FROM actual)),
               (SELECT count(*) FROM (SELECT * FROM actual EXCEPT ALL SELECT * FROM expected))
    """
    missing, extra = connection.execute(
        query, [str(legacy), str(reconstructed)]
    ).fetchone()
    return {
        "missing_from_reconstruction": missing,
        "extra_in_reconstruction": extra,
        "exact": missing == extra == 0,
    }
