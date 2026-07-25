"""Atomic Parquet persistence for compact states and event logs."""

import os
import tempfile
from collections.abc import Sequence
from pathlib import Path

import pyarrow as pa
import pyarrow.parquet as pq

from .classification import ClassificationState
from .diff import VcvEvent
from .scv_diff import ScvEvent
from .submission_parser import SubmissionRecord
from .vcf_parser import VcfRecord

_VCF_STATE_SCHEMA = pa.schema(
    [
        ("variation_id", pa.int64()),
        ("allele_id", pa.int64()),
        ("chromosome", pa.string()),
        ("position", pa.int64()),
        ("reference", pa.string()),
        ("alternate", pa.string()),
        ("raw_classification", pa.string()),
        ("classification", pa.string()),
        ("classification_ambiguous", pa.bool_()),
        ("review_status", pa.string()),
        ("gene_info", pa.string()),
        ("conditions", pa.string()),
        ("variant_type", pa.string()),
    ]
)
_VCV_EVENT_SCHEMA = pa.schema(
    [
        ("variation_id", pa.int64()),
        ("event_type", pa.string()),
        ("observed_at", pa.date32()),
        ("old_classification", pa.string()),
        ("new_classification", pa.string()),
        ("old_raw_classification", pa.string()),
        ("new_raw_classification", pa.string()),
        ("old_review_status", pa.string()),
        ("new_review_status", pa.string()),
        ("old_conditions", pa.string()),
        ("new_conditions", pa.string()),
        ("old_gene_info", pa.string()),
        ("new_gene_info", pa.string()),
        ("old_variant_type", pa.string()),
        ("new_variant_type", pa.string()),
    ]
)
_SCV_STATE_SCHEMA = pa.schema(
    [
        ("variation_id", pa.int64()),
        ("scv_accession", pa.string()),
        ("scv_version", pa.int64()),
        ("raw_classification", pa.string()),
        ("classification", pa.string()),
        ("classification_ambiguous", pa.bool_()),
        ("last_evaluated", pa.date32()),
        ("review_status", pa.string()),
        ("conditions", pa.string()),
        ("submitter", pa.string()),
        ("collection_method", pa.string()),
    ]
)
_SCV_EVENT_SCHEMA = pa.schema(
    [
        ("scv_accession", pa.string()),
        ("variation_id", pa.int64()),
        ("event_type", pa.string()),
        ("observed_at", pa.date32()),
        ("old_version", pa.int64()),
        ("new_version", pa.int64()),
        ("old_classification", pa.string()),
        ("new_classification", pa.string()),
        ("old_raw_classification", pa.string()),
        ("new_raw_classification", pa.string()),
        ("old_review_status", pa.string()),
        ("new_review_status", pa.string()),
        ("old_conditions", pa.string()),
        ("new_conditions", pa.string()),
        ("old_submitter", pa.string()),
        ("new_submitter", pa.string()),
        ("old_collection_method", pa.string()),
        ("new_collection_method", pa.string()),
    ]
)


def _atomic_write_table(table: pa.Table, output: Path) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)
    temporary_path: Path | None = None
    try:
        with tempfile.NamedTemporaryFile(
            dir=output.parent,
            prefix=f".{output.name}.",
            suffix=".tmp",
            delete=False,
        ) as temporary:
            temporary_path = Path(temporary.name)
        pq.write_table(table, temporary_path, compression="zstd")
        os.replace(temporary_path, output)
    finally:
        if temporary_path is not None and temporary_path.exists():
            temporary_path.unlink()


def _vcf_state_table(records: Sequence[VcfRecord]) -> pa.Table:
    rows = [
        {
            "variation_id": record.variation_id,
            "allele_id": record.allele_id,
            "chromosome": record.chromosome,
            "position": record.position,
            "reference": record.reference,
            "alternate": record.alternate,
            "raw_classification": record.raw_classification,
            "classification": record.classification.value,
            "classification_ambiguous": record.classification_ambiguous,
            "review_status": record.review_status,
            "gene_info": record.gene_info,
            "conditions": record.conditions,
            "variant_type": record.variant_type,
        }
        for record in records
    ]
    return pa.Table.from_pylist(rows, schema=_VCF_STATE_SCHEMA)


def write_vcf_state(records: Sequence[VcfRecord], output: Path) -> None:
    """Atomically replace the compact latest-state Parquet file."""
    _atomic_write_table(_vcf_state_table(records), output)


def read_vcf_state(path: Path) -> list[VcfRecord]:
    """Restore compact VCF records from Parquet."""
    rows = pq.read_table(path).to_pylist()
    return [
        VcfRecord(
            variation_id=row["variation_id"],
            allele_id=row["allele_id"],
            chromosome=row["chromosome"],
            position=row["position"],
            reference=row["reference"],
            alternate=row["alternate"],
            raw_classification=row["raw_classification"],
            classification=ClassificationState(row["classification"]),
            classification_ambiguous=row["classification_ambiguous"],
            review_status=row["review_status"],
            gene_info=row["gene_info"],
            conditions=row["conditions"],
            variant_type=row["variant_type"],
        )
        for row in rows
    ]


def _vcv_event_table(events: Sequence[VcvEvent]) -> pa.Table:
    rows = [
        {
            "variation_id": event.variation_id,
            "event_type": event.event_type.value,
            "observed_at": event.observed_at,
            "old_classification": (
                event.old_classification.value if event.old_classification else None
            ),
            "new_classification": (
                event.new_classification.value if event.new_classification else None
            ),
            "old_raw_classification": event.old_raw_classification,
            "new_raw_classification": event.new_raw_classification,
            "old_review_status": event.old_review_status,
            "new_review_status": event.new_review_status,
            "old_conditions": event.old_conditions,
            "new_conditions": event.new_conditions,
            "old_gene_info": event.old_gene_info,
            "new_gene_info": event.new_gene_info,
            "old_variant_type": event.old_variant_type,
            "new_variant_type": event.new_variant_type,
        }
        for event in events
    ]
    return pa.Table.from_pylist(rows, schema=_VCV_EVENT_SCHEMA)


def append_vcv_events(events: Sequence[VcvEvent], output: Path) -> None:
    """Atomically append aggregate events to a Parquet log."""
    incoming = _vcv_event_table(events)
    if output.exists():
        existing = pq.read_table(output, schema=_VCV_EVENT_SCHEMA)
        incoming = pa.concat_tables([existing, incoming])
    _atomic_write_table(incoming, output)


def write_vcv_events(events: Sequence[VcvEvent], output: Path) -> None:
    """Atomically replace one release's aggregate-event partition."""
    _atomic_write_table(_vcv_event_table(events), output)


def _scv_state_table(records: Sequence[SubmissionRecord]) -> pa.Table:
    return pa.Table.from_pylist(
        [
            {
                "variation_id": record.variation_id,
                "scv_accession": record.scv_accession,
                "scv_version": record.scv_version,
                "raw_classification": record.raw_classification,
                "classification": record.classification.value,
                "classification_ambiguous": record.classification_ambiguous,
                "last_evaluated": record.last_evaluated,
                "review_status": record.review_status,
                "conditions": record.conditions,
                "submitter": record.submitter,
                "collection_method": record.collection_method,
            }
            for record in records
        ],
        schema=_SCV_STATE_SCHEMA,
    )


def write_scv_state(records: Sequence[SubmissionRecord], output: Path) -> None:
    """Atomically replace the compact current submission state."""
    _atomic_write_table(_scv_state_table(records), output)


def read_scv_state(path: Path) -> list[SubmissionRecord]:
    """Restore compact submission records from Parquet."""
    return [
        SubmissionRecord(
            variation_id=row["variation_id"],
            scv_accession=row["scv_accession"],
            scv_version=row["scv_version"],
            raw_classification=row["raw_classification"],
            classification=ClassificationState(row["classification"]),
            classification_ambiguous=row["classification_ambiguous"],
            last_evaluated=row["last_evaluated"],
            review_status=row["review_status"],
            conditions=row["conditions"],
            submitter=row["submitter"],
            collection_method=row["collection_method"],
        )
        for row in pq.read_table(path).to_pylist()
    ]


def write_scv_events(events: Sequence[ScvEvent], output: Path) -> None:
    """Atomically replace one release's submission-event partition."""
    rows = [
        {
            "scv_accession": event.scv_accession,
            "variation_id": event.variation_id,
            "event_type": event.event_type.value,
            "observed_at": event.observed_at,
            "old_version": event.old_version,
            "new_version": event.new_version,
            "old_classification": (
                event.old_classification.value if event.old_classification else None
            ),
            "new_classification": (
                event.new_classification.value if event.new_classification else None
            ),
            "old_raw_classification": event.old_raw_classification,
            "new_raw_classification": event.new_raw_classification,
            "old_review_status": event.old_review_status,
            "new_review_status": event.new_review_status,
            "old_conditions": event.old_conditions,
            "new_conditions": event.new_conditions,
            "old_submitter": event.old_submitter,
            "new_submitter": event.new_submitter,
            "old_collection_method": event.old_collection_method,
            "new_collection_method": event.new_collection_method,
        }
        for event in events
    ]
    _atomic_write_table(
        pa.Table.from_pylist(rows, schema=_SCV_EVENT_SCHEMA),
        output,
    )
