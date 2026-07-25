"""Reproduce quarterly snapshot-based VUS rates from local event Parquet."""

from datetime import date

from _common import base_parser, write_rows

from clinvar_history.analysis_core import (
    month_grid,
    read_vcv_event_log,
    renovo_quarterly_metrics,
)


def main() -> None:
    parser = base_parser(
        __doc__,
        "results/v1/renovo_quarterly.parquet",
    )
    parser.add_argument("--start", type=date.fromisoformat, required=True)
    parser.add_argument("--end", type=date.fromisoformat, required=True)
    args = parser.parse_args()

    events = read_vcv_event_log(args.root)
    rows = renovo_quarterly_metrics(
        events,
        month_grid(args.start, args.end, 3),
    )
    write_rows(rows, args.output)


if __name__ == "__main__":
    main()
