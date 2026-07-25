"""Normalize ClinVar germline classification labels without losing source data."""

from dataclasses import dataclass
from enum import StrEnum


class ClassificationState(StrEnum):
    BENIGN = "benign"
    LIKELY_BENIGN = "likely_benign"
    VUS = "vus"
    LIKELY_PATHOGENIC = "likely_pathogenic"
    PATHOGENIC = "pathogenic"
    CONFLICT = "conflict"
    OTHER = "other"
    MISSING = "missing"


@dataclass(frozen=True, slots=True)
class NormalizedClassification:
    raw: str | None
    state: ClassificationState
    ambiguous: bool


_EXACT_STATES = {
    "benign": ClassificationState.BENIGN,
    "likely benign": ClassificationState.LIKELY_BENIGN,
    "uncertain significance": ClassificationState.VUS,
    "vus-high": ClassificationState.VUS,
    "vus-mid": ClassificationState.VUS,
    "vus-low": ClassificationState.VUS,
    "likely pathogenic": ClassificationState.LIKELY_PATHOGENIC,
    "pathogenic": ClassificationState.PATHOGENIC,
    "conflicting classifications of pathogenicity": ClassificationState.CONFLICT,
    "conflicting interpretations of pathogenicity": ClassificationState.CONFLICT,
}
_MISSING_VALUES = {"", "-", "not provided"}


def normalize(raw: str | None) -> NormalizedClassification:
    """Return a canonical state while preserving the exact source value."""
    if raw is None:
        return NormalizedClassification(raw, ClassificationState.MISSING, False)

    key = raw.strip().lower()
    if key in _MISSING_VALUES:
        return NormalizedClassification(raw, ClassificationState.MISSING, False)

    state = _EXACT_STATES.get(key)
    if state is not None:
        return NormalizedClassification(raw, state, False)

    return NormalizedClassification(raw, ClassificationState.OTHER, True)
