"""Streaming parser for the compact fields used from ClinVar VCF releases."""

from collections.abc import Iterable, Iterator
from dataclasses import dataclass
from urllib.parse import unquote

from .classification import ClassificationState, normalize


@dataclass(frozen=True, slots=True)
class VcfRecord:
    variation_id: int
    allele_id: int | None
    chromosome: str
    position: int
    reference: str
    alternate: str
    raw_classification: str | None
    classification: ClassificationState
    classification_ambiguous: bool
    review_status: str | None
    gene_info: str | None
    conditions: str | None
    variant_type: str | None


def _decode_info_value(value: str) -> str:
    return unquote(value).replace("_", " ")


def parse_info(raw: str) -> dict[str, str | None]:
    """Parse and decode a VCF INFO column."""
    parsed: dict[str, str | None] = {}
    for item in raw.split(";"):
        if not item:
            continue
        if "=" not in item:
            parsed[item] = None
            continue
        key, value = item.split("=", 1)
        parsed[key] = _decode_info_value(value)
    return parsed


def _as_text(line: str | bytes) -> str:
    if isinstance(line, bytes):
        return line.decode("utf-8")
    return line


def _optional_int(value: str | None) -> int | None:
    if value in {None, "", "."}:
        return None
    return int(value)


def parse_vcf(stream: Iterable[str | bytes]) -> Iterator[VcfRecord]:
    """Yield compact variant states from a text or binary VCF stream."""
    for line_number, raw_line in enumerate(stream, start=1):
        line = _as_text(raw_line).rstrip("\r\n")
        if not line or line.startswith("#"):
            continue

        columns = line.split("\t")
        if len(columns) < 8:
            raise ValueError(
                f"Expected at least 8 VCF columns at line {line_number}, "
                f"found {len(columns)}"
            )

        chromosome, position, variation_id, reference, alternate = columns[:5]
        info = parse_info(columns[7])
        raw_classification = info.get("CLNSIG")
        normalized = normalize(raw_classification)
        conditions = info.get("CLNDISDB") or info.get("CLNDN")

        yield VcfRecord(
            variation_id=int(variation_id),
            allele_id=_optional_int(info.get("ALLELEID")),
            chromosome=chromosome,
            position=int(position),
            reference=reference,
            alternate=alternate,
            raw_classification=raw_classification,
            classification=normalized.state,
            classification_ambiguous=normalized.ambiguous,
            review_status=info.get("CLNREVSTAT"),
            gene_info=info.get("GENEINFO"),
            conditions=conditions,
            variant_type=info.get("CLNVC"),
        )
