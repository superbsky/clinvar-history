"""Compare submission-level ClinVar states without conflating update mechanisms."""

import json
from collections.abc import Iterable, Iterator
from dataclasses import dataclass, replace
from datetime import date
from enum import StrEnum

from .classification import ClassificationState
from .submission_parser import SubmissionRecord

_DUPLICATE_CONFLICT_PREFIX = "duplicate_conflict:"


class ScvEventType(StrEnum):
    SCV_ADDED = "scv_added"
    SCV_REMOVED = "scv_removed"
    SCV_VERSION_CHANGED = "scv_version_changed"
    SCV_RECLASSIFIED = "scv_reclassified"
    REVIEW_STATUS_CHANGED = "review_status_changed"
    CONDITION_SET_CHANGED = "condition_set_changed"


@dataclass(frozen=True, slots=True)
class ScvEvent:
    scv_accession: str
    variation_id: int
    event_type: ScvEventType
    observed_at: date
    old_version: int | None
    new_version: int | None
    old_classification: ClassificationState | None
    new_classification: ClassificationState | None
    old_raw_classification: str | None
    new_raw_classification: str | None
    old_review_status: str | None
    new_review_status: str | None
    old_conditions: str | None
    new_conditions: str | None
    old_submitter: str | None = None
    new_submitter: str | None = None
    old_collection_method: str | None = None
    new_collection_method: str | None = None


def _duplicate_classification_labels(record: SubmissionRecord) -> set[str]:
    raw = record.raw_classification
    if raw is None:
        return set()
    if raw.startswith(_DUPLICATE_CONFLICT_PREFIX):
        decoded = json.loads(raw.removeprefix(_DUPLICATE_CONFLICT_PREFIX))
        if isinstance(decoded, list) and all(isinstance(item, str) for item in decoded):
            return set(decoded)
    return {raw}


def _merge_tied_classification_conflict(
    first: SubmissionRecord,
    second: SubmissionRecord,
) -> SubmissionRecord:
    labels = sorted(
        _duplicate_classification_labels(first)
        | _duplicate_classification_labels(second)
    )
    base = min(
        (first, second),
        key=lambda record: (
            record.raw_classification or "",
            record.review_status or "",
            record.conditions or "",
            record.submitter or "",
            record.collection_method or "",
        ),
    )
    return replace(
        base,
        raw_classification=(
            _DUPLICATE_CONFLICT_PREFIX + json.dumps(labels, ensure_ascii=False)
        ),
        classification=ClassificationState.CONFLICT,
        classification_ambiguous=False,
    )


def _index(records: Iterable[SubmissionRecord]) -> dict[str, SubmissionRecord]:
    indexed: dict[str, SubmissionRecord] = {}
    for record in records:
        existing = indexed.get(record.scv_accession)
        if existing is not None:
            same_identity = (
                existing.variation_id,
                existing.scv_version,
            ) == (
                record.variation_id,
                record.scv_version,
            )
            if not same_identity:
                raise ValueError(f"Duplicate SCV accession {record.scv_accession}")
            if (
                existing.classification is ClassificationState.MISSING
                and record.classification is not ClassificationState.MISSING
            ):
                indexed[record.scv_accession] = record
                continue
            if (
                record.classification is ClassificationState.MISSING
                and existing.classification is not ClassificationState.MISSING
            ):
                continue
            same_classification = record.classification is existing.classification and (
                record.classification is ClassificationState.MISSING
                or record.raw_classification == existing.raw_classification
            )
            existing_date = existing.last_evaluated or date.min
            record_date = record.last_evaluated or date.min
            if not same_classification:
                if existing_date == record_date:
                    indexed[record.scv_accession] = _merge_tied_classification_conflict(
                        existing, record
                    )
                    continue
                if record_date > existing_date:
                    indexed[record.scv_accession] = record
                continue
            if record_date > existing_date:
                indexed[record.scv_accession] = record
            continue
        indexed[record.scv_accession] = record
    return indexed


def consolidate_scv_records(
    records: Iterable[SubmissionRecord],
) -> list[SubmissionRecord]:
    """Return one deterministic current row per SCV accession."""
    return sorted(
        _index(records).values(),
        key=lambda record: record.scv_accession,
    )


def _event(
    event_type: ScvEventType,
    observed_at: date,
    old: SubmissionRecord | None,
    new: SubmissionRecord | None,
) -> ScvEvent:
    record = new or old
    assert record is not None
    return ScvEvent(
        scv_accession=record.scv_accession,
        variation_id=record.variation_id,
        event_type=event_type,
        observed_at=observed_at,
        old_version=old.scv_version if old else None,
        new_version=new.scv_version if new else None,
        old_classification=old.classification if old else None,
        new_classification=new.classification if new else None,
        old_raw_classification=old.raw_classification if old else None,
        new_raw_classification=new.raw_classification if new else None,
        old_review_status=old.review_status if old else None,
        new_review_status=new.review_status if new else None,
        old_conditions=old.conditions if old else None,
        new_conditions=new.conditions if new else None,
        old_submitter=old.submitter if old else None,
        new_submitter=new.submitter if new else None,
        old_collection_method=old.collection_method if old else None,
        new_collection_method=new.collection_method if new else None,
    )


def diff_scv_states(
    previous: Iterable[SubmissionRecord],
    current: Iterable[SubmissionRecord],
    *,
    observed_at: date,
) -> Iterator[ScvEvent]:
    """Yield distinct SCV lifecycle and field-change events."""
    old_by_accession = _index(previous)
    new_by_accession = _index(current)

    for accession in sorted(old_by_accession.keys() | new_by_accession.keys()):
        old = old_by_accession.get(accession)
        new = new_by_accession.get(accession)
        if old is None:
            yield _event(ScvEventType.SCV_ADDED, observed_at, old, new)
            continue
        if new is None:
            yield _event(ScvEventType.SCV_REMOVED, observed_at, old, new)
            continue

        if old.scv_version != new.scv_version:
            yield _event(
                ScvEventType.SCV_VERSION_CHANGED,
                observed_at,
                old,
                new,
            )
        if (
            old.classification != new.classification
            or old.raw_classification != new.raw_classification
        ):
            yield _event(ScvEventType.SCV_RECLASSIFIED, observed_at, old, new)
        if old.review_status != new.review_status:
            yield _event(
                ScvEventType.REVIEW_STATUS_CHANGED,
                observed_at,
                old,
                new,
            )
        if old.conditions != new.conditions:
            yield _event(
                ScvEventType.CONDITION_SET_CHANGED,
                observed_at,
                old,
                new,
            )
