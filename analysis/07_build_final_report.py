"""Audit the completed backfill and regenerate the final report figures."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

from clinvar_history.reporting import build_full_report_summary, render_report_figures


def _parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--data-root", type=Path, default=Path("data/full"))
    parser.add_argument("--results-root", type=Path, default=Path("results/full-v1"))
    parser.add_argument(
        "--figures-dir",
        type=Path,
        default=Path("docs/figures/full-backfill"),
    )
    parser.add_argument("--expected-vcf-releases", type=int, default=338)
    parser.add_argument("--expected-scv-releases", type=int, default=79)
    return parser


def _outcome_payload(outcome) -> dict:
    return {
        "total": outcome.total,
        "counts": outcome.counts,
        "shares": outcome.shares,
        "median_days": outcome.median_days,
    }


def _payload(summary, figures: dict[str, Path]) -> dict:
    audit = summary.audit
    quarterly = summary.quarterly
    return {
        "audit": {
            "ok": audit.ok,
            "vcf_releases": audit.vcf_releases,
            "scv_releases": audit.scv_releases,
            "completed_ledger_rows": audit.completed_ledger_rows,
            "available_official_checksums": audit.available_official_checksums,
            "matching_official_checksums": audit.matching_official_checksums,
            "cumulative_source_bytes": audit.cumulative_source_bytes,
            "retained_data_bytes": audit.retained_data_bytes,
            "results_bytes": audit.results_bytes,
            "vcv_state_rows": audit.vcv_state_rows,
            "scv_state_rows": audit.scv_state_rows,
            "vcv_event_rows": audit.vcv_event_rows,
            "scv_event_rows": audit.scv_event_rows,
            "result_rows": audit.result_rows,
            "study_start": audit.study_start.isoformat(),
            "study_end": audit.study_end.isoformat(),
        },
        "headline": {
            "initial_vus": quarterly.initial_vus,
            "final_vus": quarterly.final_vus,
            "net_vus_growth": quarterly.net_vus_growth,
            "fold_change": quarterly.fold_change,
            "total_reclassified": quarterly.total_reclassified,
            "complete_quarters": quarterly.complete_quarters,
            "partial_intervals": quarterly.partial_intervals,
            "cross_table_transition_total": audit.cross_table_transition_total,
        },
        "periods": {
            "all_complete": {
                "quarters": quarterly.complete.quarters,
                "mean_reclassification": quarterly.complete.mean_reclassification,
                "median_reclassification": quarterly.complete.median_reclassification,
                "mean_discovery": quarterly.complete.mean_discovery,
                "net_discovery_reclassification_ratio": (
                    quarterly.complete.net_discovery_reclassification_ratio
                ),
            },
            "before_split": {
                "end": quarterly.split_date.isoformat(),
                "quarters": quarterly.before_split.quarters,
                "mean_reclassification": quarterly.before_split.mean_reclassification,
                "median_reclassification": quarterly.before_split.median_reclassification,
                "mean_discovery": quarterly.before_split.mean_discovery,
                "net_discovery_reclassification_ratio": (
                    quarterly.before_split.net_discovery_reclassification_ratio
                ),
            },
            "after_split": {
                "start": quarterly.split_date.isoformat(),
                "quarters": quarterly.after_split.quarters,
                "mean_reclassification": quarterly.after_split.mean_reclassification,
                "median_reclassification": quarterly.after_split.median_reclassification,
                "mean_discovery": quarterly.after_split.mean_discovery,
                "net_discovery_reclassification_ratio": (
                    quarterly.after_split.net_discovery_reclassification_ratio
                ),
            },
        },
        "sampling": {
            cadence: {
                "true_transitions": row.true_transitions,
                "observed_transitions": row.observed_transitions,
                "recall": row.recall,
                "missed_transitions": row.missed_transitions,
                "missed_conflict_episodes": row.missed_conflict_episodes,
                "mean_event_date_error_days": row.mean_event_date_error_days,
            }
            for cadence, row in summary.sampling.items()
        },
        "mechanisms": [
            {
                "mechanism": row.mechanism,
                "count": row.count,
                "share": row.share,
            }
            for row in summary.mechanisms
        ],
        "outcomes": {
            "vcv": _outcome_payload(summary.vcv_outcomes),
            "scv": _outcome_payload(summary.scv_outcomes),
        },
        "figures": {name: str(path) for name, path in figures.items()},
    }


def main() -> None:
    args = _parser().parse_args()
    summary = build_full_report_summary(
        args.data_root,
        args.results_root,
        expected_vcf_releases=args.expected_vcf_releases,
        expected_scv_releases=args.expected_scv_releases,
    )
    figures = render_report_figures(summary, args.figures_dir)
    print(json.dumps(_payload(summary, figures), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
