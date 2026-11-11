"""Cross-check pinned published case histories and coverage; keep estimands distinct."""

import ast
import bisect
import calendar
import csv
import json
import subprocess
from datetime import date
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from clinvar_history.checkpoints import atomic_json, dataset_snapshot, file_digest
from clinvar_history.cross_study import family, liang_category, parse_history
from clinvar_history.manifest import ReleaseKind, read_manifest
from clinvar_history.research import _copy, endpoint_sql, research_connection
from clinvar_history.research_statistics import curves_from_counts, load_cluster_counts

THEIRS = Path("artifacts/references/clinvar-longitudinal")
ROOT = Path("results/liang-crosscheck")
DOCS = Path("docs/research-tables/liang-crosscheck")
DATA = Path("data/research")
THEIR_COMMIT = "dae628001d2b335068cef74e262e4431b288eb00"


def isolated_function(path, name, globals_):
    """Compile only the inspected pure function, never import repository startup code."""
    tree = ast.parse(path.read_text())
    node = next(
        n for n in tree.body if isinstance(n, ast.FunctionDef) and n.name == name
    )
    exec(compile(ast.Module(body=[node], type_ignores=[]), str(path), "exec"), globals_)  # noqa: S102 - reviewed pure function from pinned source
    return globals_[name]


def main():
    actual_commit = subprocess.check_output(
        ["git", "-C", str(THEIRS), "rev-parse", "HEAD"], text=True
    ).strip()
    if actual_commit != THEIR_COMMIT:
        raise ValueError("Reference checkout is not the pinned comparison version")
    subprocess.run(
        [
            "git",
            "-C",
            str(THEIRS),
            "diff",
            "--quiet",
            THEIR_COMMIT,
            "--",
            "scripts/clinvar_ftp_only_v5.py",
            "revision_analyses/29_km_r3.py",
        ],
        check=True,
    )
    ROOT.mkdir(parents=True, exist_ok=True)
    DOCS.mkdir(parents=True, exist_ok=True)
    source = str(DATA / "events/vcv/release_date=*/events.parquet")
    case_path = THEIRS / "data/plp_blb_reversal_variants.csv"
    case_rows = list(csv.DictReader(case_path.open()))
    assert len(case_rows) == 6340
    histories = []
    raw_history_count = 0
    for row in case_rows:
        events = parse_history(row["Event_History"])
        assert len(events) == int(row["Num_Reclass_Events"])
        raw_history_count += len(events)
        for event in events:
            old, new = family(event["old_category"]), family(event["new_category"])
            if (
                {old, new} == {"benign", "pathogenic"}
                and event["previous_month"] >= date(2020, 1, 1)
                and event["current_month"] <= date(2025, 7, 1)
            ):
                end = date(
                    event["current_month"].year,
                    event["current_month"].month,
                    calendar.monthrange(
                        event["current_month"].year, event["current_month"].month
                    )[1],
                )
                histories.append(
                    {
                        "variation_id": int(row["Variation_ID"]),
                        **event,
                        "old_family": old,
                        "new_family": new,
                        "nominal_window_end": end,
                    }
                )
    atomic_json(
        {
            "pinned_commit": THEIR_COMMIT,
            "published_rows": len(case_rows),
            "parsed_history_entries": raw_history_count,
            "eligible_reversal_edges": len(histories),
            "selection": "Reversal-case table only; previous month >=2020-01, current month <=2025-07. This is not the complete event table or a population recall denominator.",
        },
        ROOT / "design.json",
    )
    cursor_root = read_manifest(DATA / "manifest/releases.parquet")
    last_vcf = max(
        r.release_date
        for r in cursor_root
        if r.kind is ReleaseKind.VCF and r.release_date <= date(2025, 7, 1)
    )
    with research_connection(ROOT / ".scratch") as c:
        print(
            "Reconstructing all VCF states at the shared nominal July 2025 boundary",
            flush=True,
        )
        c.execute(
            """CREATE TEMP TABLE latest AS WITH ranked AS (
            SELECT variation_id,event_type,new_raw_classification,new_classification,new_gene_info,new_variant_type,
                observed_at,row_number() OVER(PARTITION BY variation_id ORDER BY observed_at DESC,event_type DESC) rn
            FROM read_parquet(?) WHERE observed_at<='2025-07-01'::DATE
                AND event_type IN ('variant_added','variant_removed','classification_changed')
          ) SELECT * EXCLUDE(rn) FROM ranked WHERE rn=1 AND event_type!='variant_removed'""",
            [source],
        )
        # Map every distinct raw value in Python; verify our mirror against the
        # exact pinned source function before any numeric comparison.
        labels = [
            r[0]
            for r in c.execute(
                "SELECT DISTINCT new_raw_classification FROM latest"
            ).fetchall()
        ]
        labels.extend(row["Current_Classification"] for row in case_rows)
        labels = list(dict.fromkeys(labels))
        original = isolated_function(
            THEIRS / "scripts/clinvar_ftp_only_v5.py", "categorize_classification", {}
        )
        assert all(liang_category(v) == original(v) for v in labels)
        c.register(
            "labels",
            pa.Table.from_pylist(
                [
                    {
                        "raw_label": v,
                        "paper_category": liang_category(v),
                        "paper_family": family(liang_category(v)),
                    }
                    for v in labels
                ]
            ),
        )
        c.execute(
            "CREATE TEMP TABLE mapped AS SELECT l.*,p.paper_category,p.paper_family FROM latest l JOIN labels p ON l.new_raw_classification IS NOT DISTINCT FROM p.raw_label"
        )
        c.execute(
            "CREATE TEMP TABLE paper_presence AS SELECT * FROM read_csv(?,header=true,types={'variation_id':'BIGINT'})",
            [str(THEIRS / "revision_analyses/outputs/variant_origin_map.csv.gz")],
        )
        coverage = c.execute("""SELECT (SELECT count(*) FROM mapped) AS our_vcf_units,
            (SELECT count(*) FROM paper_presence) AS published_xml_units,
            count(*) AS shared_units,count(*) FILTER(WHERE has_somatic) AS shared_with_somatic_origin
            FROM mapped JOIN paper_presence USING(variation_id)""").fetchone()
        c.register(
            "cases",
            pa.Table.from_pylist(
                [
                    {
                        "variation_id": int(row["Variation_ID"]),
                        "paper_raw": row["Current_Classification"],
                        "reported_category": row["Current_Category"],
                        "paper_first_month": row["First_Date"],
                        "paper_family": family(
                            liang_category(row["Current_Classification"])
                        ),
                    }
                    for row in case_rows
                ]
            ),
        )
        _copy(
            c,
            """SELECT p.*,m.observed_at AS vcf_latest_state_change,m.new_raw_classification AS vcf_raw,
            m.paper_category AS vcf_paper_category,m.paper_family AS vcf_paper_family,
            m.new_classification AS vcf_strict_classification,
            m.variation_id IS NOT NULL AS vcf_present,
            EXISTS(SELECT 1 FROM paper_presence f WHERE f.variation_id=p.variation_id) AS paper_final_present,
            m.paper_category=p.reported_category AS exact_category_agreement,
            m.paper_family=p.paper_family AS family_agreement
          FROM cases p LEFT JOIN mapped m USING(variation_id) ORDER BY variation_id""",
            ROOT / "case_snapshot_comparison.parquet",
        )
        snapshot_stats = c.execute(
            """SELECT count(*) AS cases,count(*) FILTER(WHERE vcf_present) AS present,
            count(*) FILTER(WHERE exact_category_agreement) AS exact_category_agreement,
            count(*) FILTER(WHERE family_agreement) AS family_agreement,
            count(*) FILTER(WHERE paper_final_present) AS paper_final_present,
            count(*) FILTER(WHERE paper_final_present AND vcf_present) AS shared_final_present FROM read_parquet(?)""",
            [str(ROOT / "case_snapshot_comparison.parquet")],
        ).fetchone()
        # Clinical-family edges under the paper's mapping, across ALL source
        # changes, including raw other->other changes missed by canonical counts.
        print("Comparing all 361 eligible published reversal edges", flush=True)
        c.execute(
            """CREATE TEMP TABLE raw_labels AS SELECT DISTINCT new_raw_classification AS raw_label
           FROM read_parquet(?) WHERE observed_at<='2025-08-01'::DATE""",
            [source],
        )
        all_labels = [
            r[0] for r in c.execute("SELECT raw_label FROM raw_labels").fetchall()
        ]
        assert all(liang_category(v) == original(v) for v in all_labels)
        c.register(
            "all_labels",
            pa.Table.from_pylist(
                [
                    {"raw_label": v, "paper_family": family(liang_category(v))}
                    for v in all_labels
                ]
            ),
        )
        c.execute(
            """CREATE TEMP TABLE native_states AS WITH raw AS (
            SELECT e.variation_id,e.observed_at,CASE WHEN event_type='variant_removed' THEN 'absent' ELSE a.paper_family END AS state
            FROM read_parquet(?) e LEFT JOIN all_labels a ON e.new_raw_classification IS NOT DISTINCT FROM a.raw_label
            WHERE e.observed_at<='2025-08-01'::DATE AND event_type IN ('variant_added','variant_removed','classification_changed')
          ), linked AS (SELECT *,lag(state) OVER(PARTITION BY variation_id ORDER BY observed_at) AS old_family FROM raw)
          SELECT * FROM linked WHERE old_family IS NULL OR state!=old_family""",
            [source],
        )
        c.execute("""CREATE TEMP TABLE native_reversals AS WITH linked AS (
            SELECT variation_id,observed_at,state AS new_family,lag(state) OVER(PARTITION BY variation_id ORDER BY observed_at) AS old_family FROM native_states
          ) SELECT * FROM linked WHERE (old_family='benign' AND new_family='pathogenic') OR (old_family='pathogenic' AND new_family='benign')""")
        c.register("published_edges", pa.Table.from_pylist(histories))
        _copy(
            c,
            """SELECT p.*,count(n.observed_at) AS candidate_reference_edges,
            min(n.observed_at) AS first_candidate_date,
            count(n.observed_at) FILTER(WHERE date_trunc('month',n.observed_at)=p.current_month) AS same_month_edges,
            count(n.observed_at) FILTER(WHERE abs(date_diff('day',n.observed_at,p.current_month))<=31) AS within_31_days_of_nominal_boundary,
            count(n.observed_at) FILTER(WHERE n.observed_at>p.previous_month AND n.observed_at<=p.current_month) AS inside_nominal_start_to_start_interval
          FROM published_edges p LEFT JOIN native_reversals n
            ON p.variation_id=n.variation_id AND p.old_family=n.old_family AND p.new_family=n.new_family
            AND n.observed_at>p.previous_month AND n.observed_at<=least(p.nominal_window_end,DATE '2025-07-01')
          GROUP BY p.variation_id,p.previous_month,p.old_category,p.current_month,p.new_category,p.mechanism,p.old_family,p.new_family,p.nominal_window_end ORDER BY variation_id,current_month""",
            ROOT / "published_reversal_edge_comparison.parquet",
        )
        c.execute("""CREATE TEMP TABLE edge_endpoints AS SELECT p.variation_id,p.previous_month,p.current_month,
            arg_max(n.state,n.observed_at) FILTER(WHERE n.observed_at<=p.previous_month) AS vcf_start_family,
            arg_max(n.state,n.observed_at) FILTER(WHERE n.observed_at<=p.current_month) AS vcf_end_family
            FROM published_edges p LEFT JOIN native_states n ON p.variation_id=n.variation_id AND n.observed_at<=p.current_month GROUP BY p.variation_id,p.previous_month,p.current_month""")
        _copy(
            c,
            """SELECT e.*,b.vcf_start_family,b.vcf_end_family,
            b.vcf_start_family=e.old_family AND b.vcf_end_family=e.new_family AS nominal_endpoints_compatible
            FROM read_parquet(?) e LEFT JOIN edge_endpoints b USING(variation_id,previous_month,current_month)""",
            ROOT / "reversal_endpoint_comparison.parquet",
            [str(ROOT / "published_reversal_edge_comparison.parquet")],
        )
        source_dates = sorted(
            r.release_date
            for r in cursor_root
            if r.kind is ReleaseKind.VCF and r.release_date <= date(2025, 7, 1)
        )
        grid = [
            date(2020, 1, 6),
            *[
                date(y, m, 1)
                for y in range(2020, 2026)
                for m in range(1, 13)
                if date(2020, 2, 1) <= date(y, m, 1) <= date(2025, 7, 1)
            ],
        ]
        c.register(
            "monthly_map",
            pa.table(
                {
                    "observed_at": pa.array(source_dates, type=pa.date32()),
                    "sample_at": pa.array(
                        [grid[bisect.bisect_left(grid, d)] for d in source_dates],
                        type=pa.date32(),
                    ),
                }
            ),
        )
        c.execute("""CREATE TEMP TABLE monthly_edges AS WITH ranked AS (
            SELECT n.*,m.sample_at,row_number() OVER(PARTITION BY variation_id,sample_at ORDER BY n.observed_at DESC) AS rn FROM native_states n JOIN monthly_map m USING(observed_at)
          ), linked AS (SELECT variation_id,sample_at,state AS new_family,lag(state) OVER(PARTITION BY variation_id ORDER BY sample_at) AS old_family FROM ranked WHERE rn=1)
          SELECT * FROM linked WHERE (old_family='benign' AND new_family='pathogenic') OR (old_family='pathogenic' AND new_family='benign')""")
        _copy(
            c,
            """SELECT p.variation_id,p.previous_month,p.current_month,p.old_family,p.new_family,
            e.observed_at,e.old_raw_classification,e.new_raw_classification,
            e.old_classification AS vcf_old_strict,e.new_classification AS vcf_new_strict
            FROM read_parquet(?) p JOIN read_parquet(?) e USING(variation_id)
            WHERE p.candidate_reference_edges=0 AND e.event_type='classification_changed'
                AND e.observed_at>p.previous_month AND e.observed_at<=p.current_month
            ORDER BY p.variation_id,p.current_month,e.observed_at""",
            ROOT / "unmatched_published_edge_paths.parquet",
            [str(ROOT / "reversal_endpoint_comparison.parquet"), source],
        )
        endpoint_stats = c.execute(
            """SELECT count(*) FILTER(WHERE nominal_endpoints_compatible),count(*) FILTER(WHERE candidate_reference_edges=0 AND nominal_endpoints_compatible) FROM read_parquet(?)""",
            [str(ROOT / "reversal_endpoint_comparison.parquet")],
        ).fetchone()
        monthly_stats = c.execute(
            """SELECT count(*),count(DISTINCT m.variation_id),count(*) FILTER(WHERE m.variation_id IN(SELECT variation_id FROM cases)) FROM monthly_edges m JOIN paper_presence p USING(variation_id)"""
        ).fetchone()
        _copy(
            c,
            """SELECT n.*,
            EXISTS(SELECT 1 FROM monthly_edges m WHERE m.variation_id=n.variation_id) AS variant_has_monthly_reversal,
            EXISTS(SELECT 1 FROM cases k WHERE k.variation_id=n.variation_id) AS variant_in_published_reversal_table
            FROM native_reversals n JOIN paper_presence p USING(variation_id)
            WHERE n.observed_at>'2020-01-06'::DATE AND n.observed_at<='2025-07-01'::DATE""",
            ROOT / "native_shared_reversals.parquet",
        )
        invisible_variants = c.execute(
            "SELECT count(DISTINCT variation_id) FROM read_parquet(?) WHERE NOT variant_has_monthly_reversal",
            [str(ROOT / "native_shared_reversals.parquet")],
        ).fetchone()[0]
        _copy(
            c,
            """SELECT m.*,count(n.observed_at) AS matching_native_edges
            FROM monthly_edges m JOIN paper_presence p USING(variation_id) LEFT JOIN native_reversals n
              ON m.variation_id=n.variation_id AND m.old_family=n.old_family AND m.new_family=n.new_family
              AND n.observed_at>CASE WHEN m.sample_at='2020-02-01'::DATE THEN '2020-01-06'::DATE ELSE m.sample_at-INTERVAL 1 MONTH END
              AND n.observed_at<=m.sample_at
            GROUP BY m.variation_id,m.sample_at,m.new_family,m.old_family""",
            ROOT / "monthly_shared_reversals.parquet",
        )
        matched_monthly = c.execute(
            "SELECT count(*) FILTER(WHERE matching_native_edges>0) FROM read_parquet(?)",
            [str(ROOT / "monthly_shared_reversals.parquet")],
        ).fetchone()[0]
        monthly_published_agreement = c.execute(
            """SELECT count(*) FROM published_edges p WHERE EXISTS(SELECT 1 FROM monthly_edges m WHERE p.variation_id=m.variation_id AND p.current_month=m.sample_at AND p.old_family=m.old_family AND p.new_family=m.new_family)"""
        ).fetchone()[0]
        edge_stats = c.execute(
            """SELECT count(*),count(*) FILTER(WHERE candidate_reference_edges>0),count(*) FILTER(WHERE same_month_edges>0),count(*) FILTER(WHERE within_31_days_of_nominal_boundary>0),count(*) FILTER(WHERE inside_nominal_start_to_start_interval>0) FROM read_parquet(?)""",
            [str(ROOT / "published_reversal_edge_comparison.parquet")],
        ).fetchone()
        _copy(
            c,
            """SELECT paper_raw,reported_category,vcf_raw,vcf_paper_category,count(*) AS cases FROM read_parquet(?) WHERE vcf_present AND NOT exact_category_agreement GROUP BY ALL ORDER BY cases DESC""",
            ROOT / "snapshot_disagreement_labels.parquet",
            [str(ROOT / "case_snapshot_comparison.parquet")],
        )
        # Eligible native direct reversals in shared FINAL IDs. This is a second
        # denominator; it cannot establish paper sensitivity without all histories.
        our_reverse = c.execute("""SELECT count(*),count(DISTINCT n.variation_id),count(*) FILTER(WHERE n.variation_id IN (SELECT variation_id FROM cases))
            FROM native_reversals n JOIN paper_presence p USING(variation_id)
            WHERE n.observed_at>'2020-01-06'::DATE AND n.observed_at<='2025-07-01'::DATE""").fetchone()
        old_primary = endpoint_sql(
            "e.old_raw_classification", "e.old_classification", "family"
        )
        new_primary = endpoint_sql(
            "e.new_raw_classification", "e.new_classification", "family"
        )
        _copy(
            c,
            f"""SELECT n.*,e.old_raw_classification,e.new_raw_classification,
            {old_primary} AS primary_old_family,{new_primary} AS primary_new_family,
            (({old_primary}='benign' AND {new_primary}='pathogenic') OR ({old_primary}='pathogenic' AND {new_primary}='benign')) AS primary_policy_direct_reversal
            FROM read_parquet(?) n JOIN read_parquet(?) e USING(variation_id,observed_at)
            WHERE NOT n.variant_has_monthly_reversal AND e.event_type='classification_changed' ORDER BY variation_id,observed_at""",
            ROOT / "native_only_primary_policy_check.parquet",
            [str(ROOT / "native_shared_reversals.parquet"), source],
        )
        native_only_stats = c.execute(
            "SELECT count(*),count(DISTINCT variation_id),count(*) FILTER(WHERE primary_policy_direct_reversal),count(DISTINCT variation_id) FILTER(WHERE primary_policy_direct_reversal) FROM read_parquet(?)",
            [str(ROOT / "native_only_primary_policy_check.parquet")],
        ).fetchone()
        _copy(
            c,
            """SELECT variation_id,count(*) AS native_reversal_edges,min(observed_at) AS first_observed,max(observed_at) AS last_observed,bool_or(primary_policy_direct_reversal) AS has_primary_policy_direct_reversal FROM read_parquet(?) GROUP BY variation_id ORDER BY variation_id""",
            ROOT / "native_only_variants.parquet",
            [str(ROOT / "native_only_primary_policy_check.parquet")],
        )
        _copy(
            c,
            "SELECT observed_at,count(*) AS edges FROM read_parquet(?) GROUP BY observed_at ORDER BY edges DESC,observed_at",
            ROOT / "native_only_date_counts.parquet",
            [str(ROOT / "native_only_primary_policy_check.parquet")],
        )
        batch_2023 = c.execute(
            """SELECT count(*) FROM (SELECT variation_id FROM read_parquet(?)
            WHERE observed_at BETWEEN '2023-07-10'::DATE AND '2023-07-17'::DATE GROUP BY variation_id HAVING count(*)=2
              AND arg_min(old_raw_classification,observed_at)=arg_max(new_raw_classification,observed_at)
              AND arg_min(new_raw_classification,observed_at)=arg_max(old_raw_classification,observed_at))""",
            [str(ROOT / "native_only_primary_policy_check.parquet")],
        ).fetchone()[0]
    point_checks = []
    their_aj = isolated_function(
        THEIRS / "revision_analyses/29_km_r3.py",
        "aj_from_cells",
        {"np": np, "CAUSES": ["benign", "pathogenic", "conflict", "other_exit"]},
    )
    for layer in ["vcv", "scv"]:
        data = load_cluster_counts(
            Path(f"results/full-research/episodes_{layer}_family_0.parquet")
        )
        counts = np.asarray(data.matrix @ data.frequencies).ravel()
        cells = {}
        for code, n in zip(data.codes, counts, strict=True):
            cell = cells.setdefault(int(code) // 5, np.zeros(5))
            status = int(code) % 5
            cell[status - 1 if status else 4] += n
        times, other = their_aj(cells, cells.keys())
        ours, *_ = curves_from_counts(data.codes, counts)
        err = max(
            float(np.max(np.abs(ours[times, k, 0] - other[cause])))
            for k, cause in enumerate(
                ["benign", "pathogenic", "conflict", "other_exit"]
            )
        )
        assert err < 1e-10
        point_checks.append(
            {
                "layer": layer,
                "maximum_absolute_difference": err,
                "input": "Our complete incident cohort, identical durations and four causes in both estimators",
                "limitation": "This checks the shared AJ formula, not reproduction of their published all-entry cohort or six cause definitions.",
            }
        )
    summary = {
        "their_commit": THEIR_COMMIT,
        "our_reference_commit": "b63639dd75da155b6b5084ba38ddd7760e0c4b89",
        "comparison_date": "2026-11-11",
        "nominal_snapshot": "2025-07-01",
        "actual_last_vcf": str(last_vcf),
        "mapping_source_equivalence_distinct_labels": len(all_labels),
        "published_history_entries_parsed": raw_history_count,
        "coverage": dict(
            zip(
                [
                    "our_vcf_units",
                    "published_xml_units",
                    "shared_units",
                    "shared_with_somatic_origin",
                ],
                coverage,
                strict=True,
            )
        ),
        "case_snapshot": dict(
            zip(
                [
                    "selected_cases",
                    "vcf_present",
                    "exact_category_agreement",
                    "family_agreement",
                    "paper_final_present",
                    "shared_final_present",
                ],
                snapshot_stats,
                strict=True,
            )
        ),
        "published_reversal_edges": dict(
            zip(
                [
                    "eligible_edges",
                    "matched_in_expanded_nominal_window",
                    "matched_same_calendar_month",
                    "matched_within_31_days",
                    "matched_nominal_start_to_start",
                ],
                edge_stats,
                strict=True,
            )
        ),
        "native_reversals_on_shared_final_ids": dict(
            zip(
                ["edges", "variants", "edges_in_published_selected_case_ids"],
                our_reverse,
                strict=True,
            )
        ),
        "coarse_endpoint_checks": {
            "published_edges_with_matching_nominal_endpoints": endpoint_stats[0],
            "without_direct_native_edge_but_matching_endpoints": endpoint_stats[1],
            "monthly_edges_on_shared_final_ids": monthly_stats[0],
            "monthly_edges_with_a_matching_native_edge": matched_monthly,
            "monthly_synthesized_edges": monthly_stats[0] - matched_monthly,
            "native_variants_without_a_monthly_reversal": invisible_variants,
            "monthly_variants": monthly_stats[1],
            "monthly_edges_in_published_case_ids": monthly_stats[2],
            "published_edges_matching_monthly_direction_and_boundary": monthly_published_agreement,
        },
        "native_only_primary_policy": dict(
            zip(
                [
                    "edges",
                    "variants",
                    "primary_policy_edges",
                    "primary_policy_variants",
                ],
                native_only_stats,
                strict=True,
            )
        ),
        "our_dataset": dataset_snapshot(DATA, cursor_root),
        "implementation_sha256": {
            str(p): file_digest(p)
            for p in [
                Path("analysis/19_liang_crosscheck.py"),
                Path("src/clinvar_history/cross_study.py"),
            ]
        },
        "native_only_temporal_concentration": {
            "exact_raw_round_trip_variants_2023_07_10_to_17": batch_2023,
            "interpretation": "Post hoc concentration in one exact seven-day raw-label round trip; not independent clinical reinterpretations or an established cause.",
        },
        "same_input_estimator_checks": point_checks,
        "input_sha256": {
            str(p): file_digest(p)
            for p in [
                case_path,
                THEIRS / "revision_analyses/outputs/variant_origin_map.csv.gz",
                THEIRS / "revision_analyses/outputs/km_r3.json",
                THEIRS / "scripts/clinvar_ftp_only_v5.py",
                THEIRS / "revision_analyses/29_km_r3.py",
            ]
        },
        "limits": [
            "Published reversal-case table is selected, not the full event population.",
            "Monthly nominal labels and VCF release dates are not identical observations. Expanded-window agreement is candidate compatibility, not exact replication.",
            "Full per-variant summary, event table and checkpoints are not redistributed in the pinned public repository.",
            "New November findings must not be retroactively assigned priority to the July public release.",
        ],
    }
    atomic_json(summary, ROOT / "summary.json")
    atomic_json(summary, DOCS / "summary.json")
    for name in [
        "case_snapshot_comparison",
        "published_reversal_edge_comparison",
        "snapshot_disagreement_labels",
        "reversal_endpoint_comparison",
        "native_shared_reversals",
        "monthly_shared_reversals",
        "unmatched_published_edge_paths",
        "native_only_primary_policy_check",
        "native_only_variants",
        "native_only_date_counts",
    ]:
        from pyarrow import csv as pcsv

        pcsv.write_csv(
            pq.read_table(ROOT / (name + ".parquet")), DOCS / (name + ".csv")
        )
    print(json.dumps(summary, indent=2), flush=True)


if __name__ == "__main__":
    main()
