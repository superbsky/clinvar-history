"""Distinguish net endpoint counts from recovery of recorded transition edges."""

import bisect
from datetime import date, timedelta
from pathlib import Path

import pyarrow as pa

from .analysis_core import month_grid
from .manifest import ReleaseKind, read_manifest
from .research import _copy, endpoint_sql, research_connection


def audit_sampling(
    root: Path,
    output: Path,
    *,
    policy="family",
    months=0,
    days=0,
    excluded_dates: tuple[date, ...] = (),
):
    if months < 0 or days < 0 or (months and days):
        raise ValueError("Choose one nonnegative sampling interval")
    releases = sorted(
        r.release_date
        for r in read_manifest(root / "manifest/releases.parquet")
        if r.kind is ReleaseKind.VCF
    )
    if not releases or not set(excluded_dates) <= set(releases):
        raise ValueError("Invalid source dates")
    if releases[0] in excluded_dates or releases[-1] in excluded_dates:
        raise ValueError("Cannot exclude either observation boundary")
    retained = [d for d in releases if d not in excluded_dates]
    grid = month_grid(retained[0], retained[-1], months) if months else retained
    if days:
        grid = [retained[0]]
        while grid[-1] < retained[-1]:
            grid.append(min(grid[-1] + timedelta(days=days), retained[-1]))
    source = str(root / "events/vcv/release_date=*/events.parquet")
    state = endpoint_sql("new_raw_classification", "new_classification", policy)
    with research_connection(output.parent / ".scratch") as c:
        c.register(
            "maps",
            pa.table(
                {
                    "observed_at": pa.array(releases, type=pa.date32()),
                    "reference_at": pa.array(
                        [retained[bisect.bisect_left(retained, d)] for d in releases],
                        type=pa.date32(),
                    ),
                }
            ),
        )
        c.register(
            "sample_map",
            pa.table(
                {
                    "reference_at": pa.array(retained, type=pa.date32()),
                    "sample_at": pa.array(
                        [grid[bisect.bisect_left(grid, d)] for d in retained],
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
        c.execute(
            f"""CREATE TEMP TABLE reference_states AS WITH raw AS (
                SELECT variation_id, observed_at, reference_at,
                    CASE WHEN event_type='variant_removed' THEN 'absent' ELSE {state} END AS state,
                    row_number() OVER (PARTITION BY variation_id,reference_at ORDER BY observed_at DESC,event_type DESC) AS rn
                FROM read_parquet(?) JOIN maps USING(observed_at)
                WHERE event_type IN ('variant_added','variant_removed','classification_changed')
            ), lagged AS (
                SELECT *,lag(state) OVER(PARTITION BY variation_id ORDER BY reference_at) AS previous
                FROM raw WHERE rn=1
            ) SELECT variation_id,reference_at,state FROM lagged WHERE previous IS NULL OR state!=previous""",
            [source],
        )
        c.execute("""CREATE TEMP TABLE reference_edges AS WITH linked AS (
            SELECT *,lag(state) OVER(PARTITION BY variation_id ORDER BY reference_at) AS old_state FROM reference_states
        ) SELECT variation_id,reference_at,old_state,state AS new_state FROM linked
          WHERE old_state IS NOT NULL AND old_state!='absent' AND state!='absent' AND state!=old_state""")
        query = """WITH ranked AS (
            SELECT r.*,m.sample_at,row_number() OVER(PARTITION BY variation_id,sample_at ORDER BY reference_at DESC) AS rn
            FROM reference_states r JOIN sample_map m USING(reference_at)
        ), linked AS (
            SELECT *,lag(state) OVER(PARTITION BY variation_id ORDER BY sample_at) AS old_state
            FROM ranked WHERE rn=1
        ), edges AS (
            SELECT variation_id,sample_at,previous_sample,old_state,state AS new_state
            FROM linked JOIN sample_grid USING(sample_at)
            WHERE old_state IS NOT NULL AND state!=old_state AND old_state!='absent' AND state!='absent'
        ) SELECT e.*,max(r.reference_at) AS matched_at,
            date_diff('day',max(r.reference_at),e.sample_at) AS detection_lag_days,
            count(r.reference_at) AS matching_reference_edges
        FROM edges e LEFT JOIN reference_edges r
          ON e.variation_id=r.variation_id AND e.old_state=r.old_state AND e.new_state=r.new_state
          AND r.reference_at>e.previous_sample AND r.reference_at<=e.sample_at
        GROUP BY e.variation_id,e.sample_at,e.previous_sample,e.old_state,e.new_state ORDER BY variation_id,sample_at"""
        _copy(c, query, output)
        ref, conflicts = c.execute(
            "SELECT count(*),count(*) FILTER(WHERE new_state='conflict') FROM reference_edges"
        ).fetchone()
        n, matched, lag, synthetic, conflict_net, conflict_match = c.execute(
            """SELECT count(*),count(*) FILTER(WHERE matched_at IS NOT NULL),avg(detection_lag_days),
                count(*) FILTER(WHERE matched_at IS NULL),count(*) FILTER(WHERE new_state='conflict'),
                count(*) FILTER(WHERE new_state='conflict' AND matched_at IS NOT NULL)
            FROM read_parquet(?)""",
            [str(output)],
        ).fetchone()
    return {
        "endpoint_policy": policy,
        "months": months,
        "days": days,
        "excluded_dates": [str(d) for d in excluded_dates],
        "source_releases": len(releases),
        "reference_observations": len(retained),
        "sample_observations": len(grid),
        "reference_edges": ref,
        "sampled_net_edges": n,
        "count_recovery": n / ref if ref else None,
        "matched_edges": matched,
        "edge_recall": matched / ref if ref else None,
        "synthetic_edges": synthetic,
        "synthetic_fraction": synthetic / n if n else None,
        "mean_matched_lag_days": lag,
        "reference_conflict_entries": conflicts,
        "sampled_conflict_entries": conflict_net,
        "matched_conflict_entries": conflict_match,
    }
