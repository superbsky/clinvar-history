"""Measure aggregate VUS time to resolution or another observed exit."""

from datetime import date

from _common import base_parser, write_rows

from clinvar_history.analysis_core import read_vcv_event_log, vus_resolutions


def main() -> None:
    parser = base_parser(
        __doc__,
        "results/v1/vus_time_to_resolution.parquet",
    )
    parser.add_argument("--study-end", type=date.fromisoformat)
    args = parser.parse_args()
    events = read_vcv_event_log(args.root)
    event_dates = [event.observed_at for event in events]
    if not event_dates:
        raise SystemExit("No VCV events found")
    study_end = args.study_end or max(event_dates)
    write_rows(
        vus_resolutions(events, study_end=study_end),
        args.output,
    )


if __name__ == "__main__":
    main()
