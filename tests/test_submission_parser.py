from datetime import date
from io import StringIO
from pathlib import Path

import pytest

from clinvar_history.classification import ClassificationState
from clinvar_history.submission_parser import parse_scv_accession, parse_submissions

FIXTURES = Path(__file__).parent / "fixtures"


def _records(name):
    with (FIXTURES / name).open() as stream:
        return list(parse_submissions(stream))


def test_parser_skips_comment_preamble_and_reads_old_schema():
    records = _records("submission_a.txt")

    assert len(records) == 4
    first = records[0]
    assert first.variation_id == 123
    assert first.scv_accession == "SCV000000001"
    assert first.scv_version == 1
    assert first.classification is ClassificationState.VUS
    assert first.raw_classification == "Uncertain significance"
    assert first.last_evaluated == date(2023, 12, 1)
    assert first.review_status == "criteria provided, single submitter"
    assert first.conditions == "MONDO:0001, Disease one"
    assert first.submitter == "Lab A"
    assert first.collection_method == "clinical testing"


def test_parser_accepts_appended_columns_and_spaced_scv_version():
    old_records = _records("submission_a.txt")
    current_records = _records("submission_b.txt")

    assert old_records[1].scv_accession == "SCV000000002"
    assert old_records[1].scv_version == 1
    assert current_records[0].scv_version == 2
    assert current_records[-1].classification is ClassificationState.VUS


def test_parser_accepts_historical_human_readable_evaluation_date():
    source = StringIO(
        "#VariationID\tClinicalSignificance\tDateLastEvaluated\t"
        "SubmittedPhenotypeInfo\tReportedPhenotypeInfo\tReviewStatus\t"
        "CollectionMethod\tSubmitter\tSCV\n"
        "123\tBenign\tJun 29, 2010\tDisease\tMONDO:1\tcriteria provided\t"
        "clinical testing\tLab A\tSCV000000001.1\n"
    )

    record = next(parse_submissions(source))

    assert record.last_evaluated == date(2010, 6, 29)


def test_parser_accepts_very_large_phenotype_field():
    phenotype = "MONDO:1," + "x" * 150_000
    source = StringIO(
        "#VariationID\tClinicalSignificance\tDateLastEvaluated\t"
        "SubmittedPhenotypeInfo\tReportedPhenotypeInfo\tReviewStatus\t"
        "CollectionMethod\tSubmitter\tSCV\n"
        f"123\tBenign\t2020-01-01\tDisease\t{phenotype}\tcriteria provided\t"
        "clinical testing\tLab A\tSCV000000001.1\n"
    )

    record = next(parse_submissions(source))

    assert record.conditions == phenotype


def test_parser_treats_unmatched_description_quote_as_literal_tsv_text():
    source = StringIO(
        "#VariationID\tClinicalSignificance\tDateLastEvaluated\tDescription\t"
        "SubmittedPhenotypeInfo\tReportedPhenotypeInfo\tReviewStatus\t"
        "CollectionMethod\tOriginCounts\tSubmitter\tSCV\n"
        '123\tLikely pathogenic\t2025-08-01\t"Unclosed source quote\tDisease\t'
        "MONDO:1\tcriteria provided\tclinical testing\tgermline:1\tLab A\t"
        "SCV000000001.1\n"
        "124\tBenign\t2025-08-02\t-\tDisease\tMONDO:2\tcriteria provided\t"
        "clinical testing\tgermline:1\tLab B\tSCV000000002.1\n"
    )

    records = list(parse_submissions(source))

    assert [record.scv_accession for record in records] == [
        "SCV000000001",
        "SCV000000002",
    ]


@pytest.mark.parametrize(
    ("raw", "expected"),
    [
        ("SCV000041485.7", ("SCV000041485", 7)),
        ("SCV000041485 .1", ("SCV000041485", 1)),
        ("SCV000041485", ("SCV000041485", None)),
    ],
)
def test_scv_accession_variants(raw, expected):
    assert parse_scv_accession(raw) == expected


def test_parser_rejects_missing_required_header():
    stream = StringIO("#VariationID\tClinicalSignificance\n123\tBenign\n")

    with pytest.raises(ValueError, match="SCV"):
        list(parse_submissions(stream))
