from datetime import date
from pathlib import Path

import pyarrow.parquet as pq

from clinvar_history.diff import diff_vcf_states
from clinvar_history.storage import (
    append_vcv_events,
    read_vcf_state,
    write_vcf_state,
)
from clinvar_history.vcf_parser import parse_vcf

FIXTURES = Path(__file__).parent / "fixtures"


def _records(name):
    with (FIXTURES / name).open() as stream:
        return list(parse_vcf(stream))


def test_vcf_state_round_trip_preserves_compact_records(tmp_path):
    output = tmp_path / "state" / "latest.parquet"
    records = _records("release_a.vcf")

    write_vcf_state(records, output)

    restored = read_vcf_state(output)
    assert restored == records
    assert pq.read_table(output).num_rows == 3
    assert not list(output.parent.glob("*.tmp"))


def test_state_replacement_is_atomic_and_does_not_append(tmp_path):
    output = tmp_path / "latest.parquet"
    write_vcf_state(_records("release_a.vcf"), output)

    write_vcf_state(_records("release_b.vcf"), output)

    assert read_vcf_state(output) == _records("release_b.vcf")
    assert not list(tmp_path.glob("*.tmp"))


def test_append_events_preserves_existing_event_rows(tmp_path):
    previous = _records("release_a.vcf")
    current = _records("release_b.vcf")
    first = list(diff_vcf_states(previous, current, observed_at=date(2024, 1, 8)))
    second = list(diff_vcf_states(current, previous, observed_at=date(2024, 1, 15)))
    output = tmp_path / "events.parquet"

    append_vcv_events(first, output)
    append_vcv_events(second, output)

    table = pq.read_table(output)
    assert table.num_rows == len(first) + len(second)
    assert table.column_names == [
        "variation_id",
        "event_type",
        "observed_at",
        "old_classification",
        "new_classification",
        "old_raw_classification",
        "new_raw_classification",
        "old_review_status",
        "new_review_status",
        "old_conditions",
        "new_conditions",
        "old_gene_info",
        "new_gene_info",
        "old_variant_type",
        "new_variant_type",
    ]
