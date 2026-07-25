"""Count all aggregate classification transitions from local event Parquet."""

from _common import base_parser, write_rows

from clinvar_history.analysis_core import read_vcv_event_log, transition_counts


def main() -> None:
    parser = base_parser(
        __doc__,
        "results/v1/transition_matrix.parquet",
    )
    args = parser.parse_args()
    rows = transition_counts(read_vcv_event_log(args.root))
    write_rows(rows, args.output)


if __name__ == "__main__":
    main()
