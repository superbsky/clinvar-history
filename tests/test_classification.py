import pytest

from clinvar_history.classification import ClassificationState, normalize


@pytest.mark.parametrize(
    ("raw", "state"),
    [
        ("Benign", ClassificationState.BENIGN),
        ("Likely benign", ClassificationState.LIKELY_BENIGN),
        ("Uncertain significance", ClassificationState.VUS),
        ("VUS-high", ClassificationState.VUS),
        ("Likely pathogenic", ClassificationState.LIKELY_PATHOGENIC),
        ("Pathogenic", ClassificationState.PATHOGENIC),
        (
            "Conflicting classifications of pathogenicity",
            ClassificationState.CONFLICT,
        ),
        ("Conflicting interpretations of pathogenicity", ClassificationState.CONFLICT),
    ],
)
def test_normalizes_recognized_germline_classifications(raw, state):
    result = normalize(raw)

    assert result.raw == raw
    assert result.state is state
    assert result.ambiguous is False


@pytest.mark.parametrize("raw", [None, "", "-", "not provided"])
def test_normalizes_empty_classifications_as_missing(raw):
    result = normalize(raw)

    assert result.raw == raw
    assert result.state is ClassificationState.MISSING
    assert result.ambiguous is False


@pytest.mark.parametrize(
    "raw",
    [
        "Pathogenic; drug response",
        "Benign/Likely benign",
        "Pathogenic/Likely pathogenic",
        "risk factor",
        "a future ClinVar term",
    ],
)
def test_preserves_unknown_or_composite_values_without_forcing_an_order(raw):
    result = normalize(raw)

    assert result.raw == raw
    assert result.state is ClassificationState.OTHER
    assert result.ambiguous is True


def test_normalization_result_is_immutable():
    result = normalize("Benign")

    with pytest.raises((AttributeError, TypeError)):
        result.state = ClassificationState.PATHOGENIC
