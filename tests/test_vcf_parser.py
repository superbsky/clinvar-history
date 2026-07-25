from io import BytesIO, StringIO
from pathlib import Path

import pytest

from clinvar_history.classification import ClassificationState
from clinvar_history.vcf_parser import parse_info, parse_vcf

FIXTURES = Path(__file__).parent / "fixtures"


def test_parse_info_decodes_vcf_underscores_and_flags():
    result = parse_info(
        "CLNSIG=Likely_benign;CLNREVSTAT=criteria_provided,_single_submitter;FLAG"
    )

    assert result == {
        "CLNSIG": "Likely benign",
        "CLNREVSTAT": "criteria provided, single submitter",
        "FLAG": None,
    }


def test_parse_vcf_extracts_compact_variant_state_from_text_stream():
    with (FIXTURES / "release_a.vcf").open() as stream:
        records = list(parse_vcf(stream))

    assert len(records) == 3
    first = records[0]
    assert first.variation_id == 123
    assert first.allele_id == 456
    assert first.chromosome == "1"
    assert first.position == 100
    assert first.reference == "A"
    assert first.alternate == "G"
    assert first.raw_classification == "Uncertain significance"
    assert first.classification is ClassificationState.VUS
    assert first.classification_ambiguous is False
    assert first.review_status == "criteria provided, single submitter"
    assert first.gene_info == "GENE1:1"
    assert first.conditions == "MONDO:0001|MedGen:C0001"
    assert first.variant_type == "single nucleotide variant"


def test_parse_vcf_accepts_binary_streams_and_normalizes_single_vcf_terms():
    payload = (FIXTURES / "release_a.vcf").read_bytes()

    records = list(parse_vcf(BytesIO(payload)))

    assert records[2].raw_classification == "Likely benign"
    assert records[2].classification is ClassificationState.LIKELY_BENIGN


def test_parse_vcf_keeps_missing_optional_info_fields_empty():
    with (FIXTURES / "release_a.vcf").open() as stream:
        record = list(parse_vcf(stream))[1]

    assert record.allele_id is None
    assert record.review_status is None
    assert record.gene_info is None
    assert record.conditions is None
    assert record.variant_type is None


def test_parse_vcf_rejects_malformed_data_rows():
    stream = StringIO("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n1\t2\t3\n")

    with pytest.raises(ValueError, match="Expected at least 8 VCF columns"):
        list(parse_vcf(stream))
