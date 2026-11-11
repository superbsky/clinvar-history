from datetime import date

import pytest

from clinvar_history.cross_study import family, liang_category, parse_history


def test_comparison_policy_exposes_substring_differences_without_changing_ingestion():
    from clinvar_history.classification import normalize

    assert liang_category("Uncertain risk allele") == "VUS"
    assert normalize("Uncertain risk allele").state == "other"
    assert liang_category("Pathogenic; drug response") == "P"
    assert normalize("Pathogenic; drug response").state == "other"
    assert family(liang_category("Benign/Likely benign")) == "benign"
    assert (
        liang_category("Conflicting interpretations of pathogenicity; other")
        == "Conflicting"
    )


def test_published_history_parsing_preserves_month_intervals_and_direction():
    rows = parse_history(
        "2022-07:P/LP→2022-08:LB[reinterpretation] | 2022-08:LB→2022-09:VUS[new_submitter]"
    )
    assert rows[0]["previous_month"] == date(2022, 7, 1)
    assert rows[0]["current_month"] == date(2022, 8, 1)
    assert family(rows[0]["old_category"]) == "pathogenic"
    assert family(rows[0]["new_category"]) == "benign"
    with pytest.raises(ValueError):
        parse_history("not an event")
