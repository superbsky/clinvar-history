"""Independently reproduce conservative attribution and expose SCV coverage limits."""

from pathlib import Path

from clinvar_history.checkpoints import atomic_json
from clinvar_history.manifest import ReleaseKind, read_manifest
from clinvar_history.research import _copy, endpoint_sql, research_connection


def main():
    root = Path("data/research")
    output = Path("results/full-research")
    releases = read_manifest(root / "manifest/releases.parquet")
    scv_end = max(
        r.release_date for r in releases if r.kind is ReleaseKind.SUBMISSION_SUMMARY
    )
    import pyarrow as pa

    with research_connection(output / ".scratch") as c:
        for layer, kind in [
            ("vcv", ReleaseKind.VCF),
            ("scv", ReleaseKind.SUBMISSION_SUMMARY),
        ]:
            dates = sorted(r.release_date for r in releases if r.kind is kind)
            c.register(
                layer + "_dates",
                pa.table(
                    {
                        "observed_at": pa.array(dates, type=pa.date32()),
                        "lower": pa.array([None, *dates[:-1]], type=pa.date32()),
                    }
                ),
            )
        old = endpoint_sql("old_raw_classification", "old_classification", "family")
        new = endpoint_sql("new_raw_classification", "new_classification", "family")
        c.execute(
            f"""CREATE TEMP TABLE aggregate_changes AS SELECT e.*,d.lower, {old} AS old_family,{new} AS new_family
            FROM read_parquet(?) e JOIN vcv_dates d USING(observed_at)
            WHERE event_type='classification_changed' AND old_classification!=new_classification""",
            [str(root / "events/vcv/release_date=*/events.parquet")],
        )
        c.execute(
            """CREATE TEMP TABLE submission_changes AS SELECT e.*,d.lower FROM read_parquet(?) e JOIN scv_dates d USING(observed_at)
          WHERE d.lower IS NOT NULL AND event_type IN ('scv_reclassified','scv_added','scv_removed','review_status_changed','condition_set_changed')""",
            [str(root / "events/scv/release_date=*/events.parquet")],
        )
        print("Joining temporally compatible events", flush=True)
        c.execute("""CREATE TEMP TABLE candidates AS SELECT v.variation_id,v.observed_at,
            CASE WHEN s.event_type='scv_reclassified' AND s.old_classification=v.old_classification AND s.new_classification=v.new_classification AND s.old_classification!=s.new_classification THEN 'exact_matching_scv_reclassification'
                 WHEN s.event_type='scv_added' AND lower(coalesce(s.new_review_status,'')) LIKE '%expert panel%' THEN 'expert_panel_added'
                 WHEN s.event_type='scv_added' THEN 'scv_added'
                 WHEN s.event_type='scv_removed' THEN 'scv_removed'
                 WHEN s.event_type='review_status_changed' THEN 'review_precedence_changed'
                 WHEN s.event_type='condition_set_changed' THEN 'condition_changed' END AS mechanism,
            s.scv_accession||':'||CAST(s.observed_at AS VARCHAR)||':'||s.event_type AS supporting_id
          FROM aggregate_changes v JOIN submission_changes s ON v.variation_id=s.variation_id
            AND abs(date_diff('day',v.observed_at,s.observed_at))<=35
            AND greatest(v.lower,s.lower)<least(v.observed_at,s.observed_at)""")
        _copy(
            c,
            """WITH grouped AS (SELECT variation_id,observed_at,list_sort(list(DISTINCT mechanism) FILTER(WHERE mechanism IS NOT NULL)) AS candidate_mechanisms,
          list_sort(list(DISTINCT supporting_id) FILTER(WHERE mechanism IS NOT NULL)) AS supporting_event_ids
          FROM candidates GROUP BY variation_id,observed_at)
          SELECT v.variation_id,v.observed_at,v.old_classification,v.new_classification,v.old_family,v.new_family,
            CASE WHEN coalesce(len(g.candidate_mechanisms),0)=0 THEN 'unexplained' WHEN len(g.candidate_mechanisms)=1 THEN g.candidate_mechanisms[1] ELSE 'multiple_candidates' END AS primary_mechanism,
            coalesce(g.candidate_mechanisms,[]::VARCHAR[]) AS candidate_mechanisms,coalesce(g.supporting_event_ids,[]::VARCHAR[]) AS supporting_event_ids,
            CASE WHEN v.lower>=?::DATE THEN 'no_scv_interval_overlap_after_source_end' WHEN v.observed_at>?::DATE THEN 'partial_scv_interval_coverage' ELSE 'within_scv_observation_horizon' END AS scv_observability
          FROM aggregate_changes v LEFT JOIN grouped g USING(variation_id,observed_at) ORDER BY v.variation_id,v.observed_at""",
            output / "attributions_reconstructed.parquet",
            [str(scv_end), str(scv_end)],
        )
        discrepancies = c.execute(
            """SELECT count(*) FROM read_parquet(?) baseline FULL JOIN (SELECT * FROM read_parquet(?) WHERE observed_at<='2026-07-15'::DATE) current USING(variation_id,observed_at)
          WHERE baseline.primary_mechanism IS DISTINCT FROM current.primary_mechanism OR list_sort(baseline.candidate_mechanisms) IS DISTINCT FROM current.candidate_mechanisms OR baseline.supporting_event_ids IS DISTINCT FROM current.supporting_event_ids""",
            [
                "results/full-v1/attributions.parquet",
                str(output / "attributions_reconstructed.parquet"),
            ],
        ).fetchone()[0]
        _copy(
            c,
            """SELECT primary_mechanism,scv_observability,old_family=new_family AS within_same_family,count(*) AS transitions FROM read_parquet(?) GROUP BY ALL ORDER BY transitions DESC""",
            output / "mechanism_coverage.parquet",
            [str(output / "attributions_reconstructed.parquet")],
        )
        result = {
            "baseline_rows_compared": 479525,
            "discrepancies": discrepancies,
            "scv_last_nominal_source_date": str(scv_end),
            "max_lag_days": 35,
            "interval_overlap": "max(lower_vcv,lower_scv)<min(upper_vcv,upper_scv)",
            "causal_validation": False,
        }
        atomic_json(result, output / "attribution_audit.json")
        assert discrepancies == 0, result
        print(result, flush=True)


if __name__ == "__main__":
    main()
