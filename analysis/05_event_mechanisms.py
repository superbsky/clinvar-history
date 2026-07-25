"""Summarize conservative aggregate-event attribution mechanisms."""

from pathlib import Path

from _common import base_parser, write_rows

from clinvar_history.analysis_core import mechanism_counts, read_attributions


def main() -> None:
    parser = base_parser(
        __doc__,
        "results/v1/event_mechanisms.parquet",
    )
    parser.add_argument(
        "--attributions",
        type=Path,
        default="results/attributions.parquet",
    )
    args = parser.parse_args()
    counts = mechanism_counts(read_attributions(args.attributions))
    rows = [
        {"primary_mechanism": mechanism.value, "count": count}
        for mechanism, count in counts.items()
    ]
    write_rows(rows, args.output)


if __name__ == "__main__":
    main()
