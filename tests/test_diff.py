from dataclasses import replace
from datetime import date
from pathlib import Path

from clinvar_history.classification import ClassificationState
from clinvar_history.diff import EventType, diff_vcf_states
from clinvar_history.vcf_parser import parse_vcf

FIXTURES = Path(__file__).parent / "fixtures"


def _records(name):
    with (FIXTURES / name).open() as stream:
        return list(parse_vcf(stream))


def test_diff_reports_add_remove_and_classification_change_only():
    events = list(
        diff_vcf_states(
            _records("release_a.vcf"),
            _records("release_b.vcf"),
            observed_at=date(2024, 1, 8),
        )
    )

    assert [(event.variation_id, event.event_type) for event in events] == [
        (123, EventType.CLASSIFICATION_CHANGED),
        (125, EventType.VARIANT_REMOVED),
        (126, EventType.VARIANT_ADDED),
    ]
    changed = events[0]
    assert changed.old_classification is ClassificationState.VUS
    assert changed.new_classification is ClassificationState.LIKELY_PATHOGENIC
    assert changed.observed_at == date(2024, 1, 8)
    assert changed.old_gene_info == "GENE1:1"
    assert changed.new_gene_info == "GENE1:1"
    assert changed.old_variant_type == "single nucleotide variant"
    assert changed.new_variant_type == "single nucleotide variant"


def test_diff_emits_separate_review_and_condition_events():
    previous = _records("release_a.vcf")[0]
    current = replace(
        previous,
        review_status="reviewed by expert panel",
        conditions="MONDO:9999",
    )

    events = list(diff_vcf_states([previous], [current], observed_at=date(2024, 2, 1)))

    assert [event.event_type for event in events] == [
        EventType.REVIEW_STATUS_CHANGED,
        EventType.CONDITION_SET_CHANGED,
    ]


def test_diff_rejects_duplicate_variation_ids():
    record = _records("release_a.vcf")[0]

    try:
        list(diff_vcf_states([record, record], [record], observed_at=date(2024, 1, 8)))
    except ValueError as error:
        assert "Duplicate VariationID 123" in str(error)
    else:
        raise AssertionError("duplicate state key was accepted")
