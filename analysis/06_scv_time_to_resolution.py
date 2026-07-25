"""Measure submitter-level VUS time to resolution or another observed exit."""

from datetime import date

from _common import base_parser, write_rows

from clinvar_history.analysis_core import (
    read_scv_event_log,
    scv_vus_resolutions,
)


def main() -> None:
    parser = base_parser(
        __doc__,
        "results/v1/scv_vus_time_to_resolution.parquet",
    )
    parser.add_argument("--study-end", type=date.fromisoformat)
    args = parser.parse_args()
    events = read_scv_event_log(args.root)
    event_dates = [event.observed_at for event in events]
    if not event_dates:
        raise SystemExit("No SCV events found")
    study_end = args.study_end or max(event_dates)
    write_rows(
        scv_vus_resolutions(events, study_end=study_end),
        args.output,
    )


if __name__ == "__main__":
    main()
