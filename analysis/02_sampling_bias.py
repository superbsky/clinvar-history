"""Compare event recall at weekly, monthly, and quarterly cadence."""

from dataclasses import asdict
from datetime import timedelta

from _common import base_parser, write_rows

from clinvar_history.analysis_core import (
    month_grid,
    read_vcv_event_log,
    sampling_bias,
)


def _weekly_grid(start, end):
    dates = [start]
    while dates[-1] < end:
        dates.append(min(dates[-1] + timedelta(days=7), end))
    return dates


def main() -> None:
    parser = base_parser(
        __doc__,
        "results/v1/sampling_bias.parquet",
    )
    args = parser.parse_args()
    events = read_vcv_event_log(args.root)
    event_dates = [event.observed_at for event in events]
    if not event_dates:
        raise SystemExit("No VCV events found")
    start, end = min(event_dates), max(event_dates)

    grids = {
        "weekly": _weekly_grid(start, end),
        "monthly": month_grid(start, end, 1),
        "quarterly": month_grid(start, end, 3),
    }
    rows = [
        {"cadence": cadence, **asdict(sampling_bias(events, dates))}
        for cadence, dates in grids.items()
    ]
    write_rows(rows, args.output)


if __name__ == "__main__":
    main()
