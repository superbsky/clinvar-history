"""Offline managed analysis stages, also usable from an installed wheel."""

import argparse
from dataclasses import asdict, is_dataclass
from datetime import date, timedelta
from enum import Enum
from pathlib import Path

import pyarrow as pa

from . import analysis_core as analysis
from .storage import _atomic_write_table

_EPISODE = [
    ("variation_id", pa.int64()),
    ("started_at", pa.date32()),
    ("ended_at", pa.date32()),
    ("duration_days", pa.int64()),
    ("outcome", pa.string()),
]
_SCHEMAS = {
    "renovo_quarterly": pa.schema(
        [(n, pa.date32()) for n in ("period_start", "period_end")]
        + [
            (n, pa.int64())
            for n in (
                "vus_previous",
                "vus_current",
                "reclassified_vus",
                "net_vus_change",
            )
        ]
        + [
            (n, pa.float64())
            for n in (
                "reclassification_rate",
                "discovery_rate",
                "discovery_reclassification_ratio",
            )
        ]
    ),
    "sampling_bias": pa.schema(
        [
            ("cadence", pa.string()),
            ("true_transitions", pa.int64()),
            ("observed_transitions", pa.int64()),
            ("event_recall", pa.float64()),
            ("missed_recurrent_transitions", pa.int64()),
            ("missed_conflict_episodes", pa.int64()),
            ("mean_event_date_error_days", pa.float64()),
        ]
    ),
    "transition_matrix": pa.schema(
        [
            ("old_classification", pa.string()),
            ("new_classification", pa.string()),
            ("count", pa.int64()),
        ]
    ),
    "vus_time_to_resolution": pa.schema(
        _EPISODE
        + [
            (n, pa.string())
            for n in ("review_status", "gene_info", "conditions", "variant_type")
        ]
    ),
    "scv_vus_time_to_resolution": pa.schema(
        [("scv_accession", pa.string())]
        + _EPISODE
        + [
            (n, pa.string())
            for n in ("submitter", "collection_method", "review_status", "conditions")
        ]
    ),
    "event_mechanisms": pa.schema(
        [("primary_mechanism", pa.string()), ("count", pa.int64())]
    ),
}


def _plain(value):
    if isinstance(value, Enum):
        return value.value
    if is_dataclass(value):
        return _plain(asdict(value))
    if isinstance(value, dict):
        return {str(_plain(k)): _plain(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [_plain(v) for v in value]
    return value


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--stage",
        required=True,
        choices=[
            "renovo_quarterly",
            "sampling_bias",
            "transition_matrix",
            "vus_time_to_resolution",
            "event_mechanisms",
            "scv_vus_time_to_resolution",
        ],
    )
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--start", type=date.fromisoformat)
    parser.add_argument("--end", type=date.fromisoformat)
    parser.add_argument("--study-end", type=date.fromisoformat)
    parser.add_argument("--attributions", type=Path)
    args = parser.parse_args()
    if args.stage == "event_mechanisms":
        rows = [
            {"primary_mechanism": k.value, "count": v}
            for k, v in analysis.mechanism_counts(
                analysis.read_attributions(args.attributions)
            ).items()
        ]
    elif args.stage == "scv_vus_time_to_resolution":
        rows = analysis.scv_vus_resolutions(
            analysis.read_scv_event_log(args.root), study_end=args.study_end
        )
    else:
        events = analysis.read_vcv_event_log(args.root)
        if not events:
            raise SystemExit("No VCV events found")
        if args.stage == "renovo_quarterly":
            rows = analysis.renovo_quarterly_metrics(
                events, analysis.month_grid(args.start, args.end, 3)
            )
        elif args.stage == "transition_matrix":
            rows = analysis.transition_counts(events)
        elif args.stage == "vus_time_to_resolution":
            rows = analysis.vus_resolutions(events, study_end=args.study_end)
        else:
            start, end = (
                min(e.observed_at for e in events),
                max(e.observed_at for e in events),
            )
            weekly = [start]
            while weekly[-1] < end:
                weekly.append(min(weekly[-1] + timedelta(days=7), end))
            grids = {
                "weekly": weekly,
                "monthly": analysis.month_grid(start, end, 1),
                "quarterly": analysis.month_grid(start, end, 3),
            }
            rows = [
                {"cadence": cadence, **asdict(analysis.sampling_bias(events, dates))}
                for cadence, dates in grids.items()
            ]
    _atomic_write_table(
        pa.Table.from_pylist(
            [_plain(row) for row in rows], schema=_SCHEMAS[args.stage]
        ),
        args.output,
    )


if __name__ == "__main__":
    main()
