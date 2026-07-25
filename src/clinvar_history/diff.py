"""Compare compact ClinVar states and emit only observable changes."""

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from datetime import date
from enum import StrEnum

from .classification import ClassificationState
from .vcf_parser import VcfRecord


class EventType(StrEnum):
    VARIANT_ADDED = "variant_added"
    VARIANT_REMOVED = "variant_removed"
    CLASSIFICATION_CHANGED = "classification_changed"
    REVIEW_STATUS_CHANGED = "review_status_changed"
    CONDITION_SET_CHANGED = "condition_set_changed"


@dataclass(frozen=True, slots=True)
class VcvEvent:
    variation_id: int
    event_type: EventType
    observed_at: date
    old_classification: ClassificationState | None
    new_classification: ClassificationState | None
    old_raw_classification: str | None
    new_raw_classification: str | None
    old_review_status: str | None
    new_review_status: str | None
    old_conditions: str | None
    new_conditions: str | None
    old_gene_info: str | None = None
    new_gene_info: str | None = None
    old_variant_type: str | None = None
    new_variant_type: str | None = None


def _index(records: Iterable[VcfRecord]) -> dict[int, VcfRecord]:
    indexed: dict[int, VcfRecord] = {}
    for record in records:
        if record.variation_id in indexed:
            raise ValueError(f"Duplicate VariationID {record.variation_id}")
        indexed[record.variation_id] = record
    return indexed


def _event(
    variation_id: int,
    event_type: EventType,
    observed_at: date,
    old: VcfRecord | None,
    new: VcfRecord | None,
) -> VcvEvent:
    return VcvEvent(
        variation_id=variation_id,
        event_type=event_type,
        observed_at=observed_at,
        old_classification=old.classification if old else None,
        new_classification=new.classification if new else None,
        old_raw_classification=old.raw_classification if old else None,
        new_raw_classification=new.raw_classification if new else None,
        old_review_status=old.review_status if old else None,
        new_review_status=new.review_status if new else None,
        old_conditions=old.conditions if old else None,
        new_conditions=new.conditions if new else None,
        old_gene_info=old.gene_info if old else None,
        new_gene_info=new.gene_info if new else None,
        old_variant_type=old.variant_type if old else None,
        new_variant_type=new.variant_type if new else None,
    )


def diff_vcf_states(
    previous: Iterable[VcfRecord],
    current: Iterable[VcfRecord],
    *,
    observed_at: date,
) -> Iterator[VcvEvent]:
    """Yield deterministic events between two complete VCF states."""
    old_by_id = _index(previous)
    new_by_id = _index(current)

    for variation_id in sorted(old_by_id.keys() | new_by_id.keys()):
        old = old_by_id.get(variation_id)
        new = new_by_id.get(variation_id)
        if old is None:
            yield _event(variation_id, EventType.VARIANT_ADDED, observed_at, old, new)
            continue
        if new is None:
            yield _event(variation_id, EventType.VARIANT_REMOVED, observed_at, old, new)
            continue

        if (
            old.classification != new.classification
            or old.raw_classification != new.raw_classification
        ):
            yield _event(
                variation_id,
                EventType.CLASSIFICATION_CHANGED,
                observed_at,
                old,
                new,
            )
        if old.review_status != new.review_status:
            yield _event(
                variation_id,
                EventType.REVIEW_STATUS_CHANGED,
                observed_at,
                old,
                new,
            )
        if old.conditions != new.conditions:
            yield _event(
                variation_id,
                EventType.CONDITION_SET_CHANGED,
                observed_at,
                old,
                new,
            )
