"""Validated summaries and dependency-free figures for the final report."""

from __future__ import annotations

import math
from dataclasses import dataclass
from datetime import date
from pathlib import Path
from statistics import mean, median
from typing import Any

import duckdb
import pyarrow as pa
import pyarrow.parquet as pq

from .classification import ClassificationState


class ReportValidationError(ValueError):
    """Raised when a report input violates a required invariant."""


@dataclass(frozen=True, slots=True)
class QuarterlyPoint:
    period_start: date
    period_end: date
    vus_previous: int
    vus_current: int
    reclassified_vus: int
    net_vus_change: int
    reclassification_rate: float
    discovery_rate: float
    discovery_reclassification_ratio: float | None
    complete: bool


@dataclass(frozen=True, slots=True)
class PeriodSummary:
    quarters: int
    mean_reclassification: float
    median_reclassification: float
    mean_discovery: float
    net_discovery_reclassification_ratio: float


@dataclass(frozen=True, slots=True)
class QuarterlySummary:
    points: tuple[QuarterlyPoint, ...]
    initial_vus: int
    final_vus: int
    net_vus_growth: int
    fold_change: float
    total_reclassified: int
    complete_quarters: int
    partial_intervals: int
    complete: PeriodSummary
    before_split: PeriodSummary
    after_split: PeriodSummary
    split_date: date


@dataclass(frozen=True, slots=True)
class CadenceSummary:
    cadence: str
    true_transitions: int
    observed_transitions: int
    recall: float
    missed_transitions: int
    missed_conflict_episodes: int
    mean_event_date_error_days: float


@dataclass(frozen=True, slots=True)
class MechanismSummary:
    mechanism: str
    count: int
    share: float


@dataclass(frozen=True, slots=True)
class OutcomeSummary:
    total: int
    counts: dict[str, int]
    shares: dict[str, float]
    median_days: dict[str, float]


@dataclass(frozen=True, slots=True)
class TransitionSummary:
    states: tuple[str, ...]
    counts: dict[tuple[str, str], int]
    total: int


@dataclass(frozen=True, slots=True)
class FullBackfillAudit:
    ok: bool
    vcf_releases: int
    scv_releases: int
    completed_ledger_rows: int
    available_official_checksums: int
    matching_official_checksums: int
    cumulative_source_bytes: int
    retained_data_bytes: int
    results_bytes: int
    vcv_state_rows: int
    scv_state_rows: int
    vcv_event_rows: int
    scv_event_rows: int
    result_rows: dict[str, int]
    cross_table_transition_total: int
    study_start: date
    study_end: date


@dataclass(frozen=True, slots=True)
class FullReportSummary:
    audit: FullBackfillAudit
    quarterly: QuarterlySummary
    sampling: dict[str, CadenceSummary]
    mechanisms: tuple[MechanismSummary, ...]
    transitions: TransitionSummary
    vcv_outcomes: OutcomeSummary
    scv_outcomes: OutcomeSummary


_OUTCOMES = ("censored", "conflict", "benign", "pathogenic", "other_exit")
_STATE_ORDER = tuple(state.value for state in ClassificationState)
_RESULT_COLUMNS = {
    "attributions.parquet": {
        "variation_id",
        "observed_at",
        "primary_mechanism",
        "candidate_mechanisms",
        "supporting_event_ids",
    },
    "event_mechanisms.parquet": {"primary_mechanism", "count"},
    "renovo_quarterly.parquet": {
        "period_start",
        "period_end",
        "vus_previous",
        "vus_current",
        "reclassified_vus",
        "net_vus_change",
        "reclassification_rate",
        "discovery_rate",
        "discovery_reclassification_ratio",
    },
    "sampling_bias.parquet": {
        "cadence",
        "true_transitions",
        "observed_transitions",
        "event_recall",
        "missed_recurrent_transitions",
        "missed_conflict_episodes",
        "mean_event_date_error_days",
    },
    "scv_vus_time_to_resolution.parquet": {
        "scv_accession",
        "variation_id",
        "started_at",
        "ended_at",
        "duration_days",
        "outcome",
    },
    "transition_matrix.parquet": {
        "old_classification",
        "new_classification",
        "count",
    },
    "vus_time_to_resolution.parquet": {
        "variation_id",
        "started_at",
        "ended_at",
        "duration_days",
        "outcome",
    },
}


def _rows(table: pa.Table, required: set[str]) -> list[dict[str, Any]]:
    missing = required.difference(table.column_names)
    if missing:
        names = ", ".join(sorted(missing))
        raise ReportValidationError(f"Missing required columns: {names}")
    rows = table.to_pylist()
    if not rows:
        raise ReportValidationError("Report input table is empty")
    return rows


def _finite(value: float, label: str) -> float:
    result = float(value)
    if not math.isfinite(result):
        raise ReportValidationError(f"Non-finite {label}")
    return result


def _require_close(actual: float, expected: float, label: str) -> None:
    if not math.isclose(actual, expected, rel_tol=1e-9, abs_tol=1e-12):
        raise ReportValidationError(
            f"Declared {label} does not match recomputed value: "
            f"{actual!r} != {expected!r}"
        )


def _period_summary(points: list[QuarterlyPoint]) -> PeriodSummary:
    if not points:
        raise ReportValidationError("Quarterly comparison period is empty")
    reclassified = sum(point.reclassified_vus for point in points)
    net_discovery = sum(point.net_vus_change for point in points)
    if reclassified == 0:
        raise ReportValidationError("Quarterly comparison has zero reclassifications")
    return PeriodSummary(
        quarters=len(points),
        mean_reclassification=mean(point.reclassification_rate for point in points),
        median_reclassification=median(point.reclassification_rate for point in points),
        mean_discovery=mean(point.discovery_rate for point in points),
        net_discovery_reclassification_ratio=net_discovery / reclassified,
    )


def summarize_quarters(table: pa.Table, split_date: date) -> QuarterlySummary:
    """Validate and summarize consecutive-snapshot quarterly VUS metrics."""
    required = {
        "period_start",
        "period_end",
        "vus_previous",
        "vus_current",
        "reclassified_vus",
        "net_vus_change",
        "reclassification_rate",
        "discovery_rate",
        "discovery_reclassification_ratio",
    }
    source_rows = sorted(_rows(table, required), key=lambda row: row["period_start"])
    points: list[QuarterlyPoint] = []
    previous_current: int | None = None
    for row in source_rows:
        period_start = row["period_start"]
        period_end = row["period_end"]
        vus_previous = int(row["vus_previous"])
        vus_current = int(row["vus_current"])
        reclassified = int(row["reclassified_vus"])
        net_change = int(row["net_vus_change"])
        if period_end <= period_start:
            raise ReportValidationError("Quarterly interval is not positive")
        if previous_current is not None and vus_previous != previous_current:
            raise ReportValidationError("Quarterly VUS inventory chain is inconsistent")
        if vus_current <= 0 or reclassified < 0:
            raise ReportValidationError("Quarterly counts must be non-negative")
        if net_change != vus_current - vus_previous:
            raise ReportValidationError("Quarterly net VUS change is inconsistent")

        reclassification_rate = reclassified / vus_current
        discovery_rate = net_change / vus_current
        declared_reclassification = _finite(
            row["reclassification_rate"], "quarterly reclassification rate"
        )
        declared_discovery = _finite(row["discovery_rate"], "quarterly discovery rate")
        _require_close(
            declared_reclassification,
            reclassification_rate,
            "quarterly reclassification rate",
        )
        _require_close(declared_discovery, discovery_rate, "quarterly discovery rate")

        declared_ratio = row["discovery_reclassification_ratio"]
        ratio = None if reclassified == 0 else net_change / reclassified
        if ratio is not None:
            _require_close(
                _finite(declared_ratio, "discovery/reclassification ratio"),
                ratio,
                "discovery/reclassification ratio",
            )
        elif declared_ratio is not None:
            raise ReportValidationError(
                "Discovery/reclassification ratio must be null when denominator is zero"
            )

        points.append(
            QuarterlyPoint(
                period_start=period_start,
                period_end=period_end,
                vus_previous=vus_previous,
                vus_current=vus_current,
                reclassified_vus=reclassified,
                net_vus_change=net_change,
                reclassification_rate=reclassification_rate,
                discovery_rate=discovery_rate,
                discovery_reclassification_ratio=ratio,
                complete=(period_end - period_start).days >= 80,
            )
        )
        previous_current = vus_current

    complete_points = [point for point in points if point.complete]
    before = [point for point in complete_points if point.period_end <= split_date]
    after = [point for point in complete_points if point.period_start >= split_date]
    initial_vus = points[0].vus_previous
    final_vus = points[-1].vus_current
    return QuarterlySummary(
        points=tuple(points),
        initial_vus=initial_vus,
        final_vus=final_vus,
        net_vus_growth=final_vus - initial_vus,
        fold_change=final_vus / initial_vus,
        total_reclassified=sum(point.reclassified_vus for point in points),
        complete_quarters=len(complete_points),
        partial_intervals=len(points) - len(complete_points),
        complete=_period_summary(complete_points),
        before_split=_period_summary(before),
        after_split=_period_summary(after),
        split_date=split_date,
    )


def summarize_sampling(table: pa.Table) -> dict[str, CadenceSummary]:
    """Validate cadence rows and recompute event recall."""
    required = {
        "cadence",
        "true_transitions",
        "observed_transitions",
        "event_recall",
        "missed_recurrent_transitions",
        "missed_conflict_episodes",
        "mean_event_date_error_days",
    }
    result: dict[str, CadenceSummary] = {}
    for row in _rows(table, required):
        cadence = str(row["cadence"])
        if cadence in result:
            raise ReportValidationError(f"Duplicate sampling cadence: {cadence}")
        true = int(row["true_transitions"])
        observed = int(row["observed_transitions"])
        missed = true - observed
        declared_missed = int(row["missed_recurrent_transitions"])
        if true <= 0 or not 0 <= observed <= true:
            raise ReportValidationError("Invalid sampling transition counts")
        if declared_missed != missed:
            raise ReportValidationError("Declared missed transitions are inconsistent")
        recall = observed / true
        _require_close(
            _finite(row["event_recall"], "sampling recall"),
            recall,
            "sampling recall",
        )
        result[cadence] = CadenceSummary(
            cadence=cadence,
            true_transitions=true,
            observed_transitions=observed,
            recall=recall,
            missed_transitions=missed,
            missed_conflict_episodes=int(row["missed_conflict_episodes"]),
            mean_event_date_error_days=_finite(
                row["mean_event_date_error_days"], "mean event-date error"
            ),
        )
    return result


def summarize_mechanisms(table: pa.Table) -> tuple[MechanismSummary, ...]:
    """Return attribution mechanisms sorted by count with recomputed shares."""
    rows = _rows(table, {"primary_mechanism", "count"})
    counts: dict[str, int] = {}
    for row in rows:
        mechanism = str(row["primary_mechanism"])
        if mechanism in counts:
            raise ReportValidationError(f"Duplicate mechanism: {mechanism}")
        count = int(row["count"])
        if count < 0:
            raise ReportValidationError("Mechanism count must be non-negative")
        counts[mechanism] = count
    total = sum(counts.values())
    if total == 0:
        raise ReportValidationError("Mechanism total is zero")
    return tuple(
        MechanismSummary(mechanism=name, count=count, share=count / total)
        for name, count in sorted(counts.items(), key=lambda item: (-item[1], item[0]))
    )


def summarize_outcomes(table: pa.Table) -> OutcomeSummary:
    """Summarize raw episode outcomes and outcome-specific observed durations."""
    rows = _rows(table, {"duration_days", "outcome"})
    counts = {outcome: 0 for outcome in _OUTCOMES}
    durations: dict[str, list[int]] = {outcome: [] for outcome in _OUTCOMES}
    for row in rows:
        outcome = str(row["outcome"])
        if outcome not in counts:
            raise ReportValidationError(f"Unknown episode outcome: {outcome}")
        duration = int(row["duration_days"])
        if duration < 0:
            raise ReportValidationError("Episode duration must be non-negative")
        counts[outcome] += 1
        durations[outcome].append(duration)
    total = len(rows)
    return OutcomeSummary(
        total=total,
        counts=counts,
        shares={outcome: count / total for outcome, count in counts.items()},
        median_days={
            outcome: median(values) for outcome, values in durations.items() if values
        },
    )


def summarize_transitions(table: pa.Table) -> TransitionSummary:
    """Validate a directed canonical transition matrix."""
    rows = _rows(table, {"old_classification", "new_classification", "count"})
    counts: dict[tuple[str, str], int] = {}
    valid_states = set(_STATE_ORDER)
    for row in rows:
        pair = (str(row["old_classification"]), str(row["new_classification"]))
        if pair[0] not in valid_states or pair[1] not in valid_states:
            raise ReportValidationError(f"Unknown canonical transition state: {pair}")
        if pair in counts:
            raise ReportValidationError(f"Duplicate transition pair: {pair}")
        count = int(row["count"])
        if count < 0:
            raise ReportValidationError("Transition count must be non-negative")
        counts[pair] = count
    return TransitionSummary(
        states=_STATE_ORDER,
        counts=counts,
        total=sum(counts.values()),
    )


def _file_bytes(root: Path) -> int:
    total = 0
    seen: set[tuple[int, int]] = set()
    for path in root.rglob("*"):
        if not path.is_file():
            continue
        stat = path.stat()
        identity = (stat.st_dev, stat.st_ino)
        if identity in seen:
            continue
        seen.add(identity)
        total += getattr(stat, "st_blocks", (stat.st_size + 511) // 512) * 512
    return total


def _parquet_rows(path: Path) -> int:
    try:
        return pq.ParquetFile(path).metadata.num_rows
    except (FileNotFoundError, OSError, pa.ArrowException) as error:
        raise ReportValidationError(f"Unreadable Parquet artifact: {path}") from error


def _validate_result_files(results_root: Path) -> dict[str, int]:
    row_counts: dict[str, int] = {}
    for filename, required in _RESULT_COLUMNS.items():
        path = results_root / filename
        rows = _parquet_rows(path)
        columns = set(pq.ParquetFile(path).schema_arrow.names)
        missing = required.difference(columns)
        if missing:
            names = ", ".join(sorted(missing))
            raise ReportValidationError(f"{filename} is missing columns: {names}")
        row_counts[filename] = rows
    return row_counts


def _validate_latest_state(
    connection: duckdb.DuckDBPyConnection,
    path: Path,
    key: str,
) -> int:
    states = list(_STATE_ORDER)
    placeholders = ", ".join("?" for _ in states)
    query = f"""
        SELECT
            count(*) AS rows,
            count(DISTINCT {key}) AS unique_keys,
            count(*) FILTER (
                WHERE classification IS NULL
                   OR classification NOT IN ({placeholders})
            ) AS invalid_states
        FROM read_parquet(?)
    """
    rows, unique_keys, invalid_states = connection.execute(
        query,
        [*states, str(path)],
    ).fetchone()
    if rows != unique_keys:
        raise ReportValidationError(f"Latest-state primary key is not unique: {path}")
    if invalid_states:
        raise ReportValidationError(
            f"Latest state contains invalid classification: {path}"
        )
    return rows


def _validate_event_store(
    connection: duckdb.DuckDBPyConnection,
    root: Path,
    expected_partitions: int,
) -> int:
    paths = sorted(root.glob("release_date=*/events.parquet"))
    if len(paths) != expected_partitions:
        raise ReportValidationError(
            f"Event partition count mismatch for {root}: "
            f"{len(paths)} != {expected_partitions}"
        )
    rows = sum(_parquet_rows(path) for path in paths)
    mismatch = connection.execute(
        """
        SELECT count(*)
        FROM read_parquet(?, filename = true)
        WHERE observed_at != CAST(
            regexp_extract(
                filename,
                'release_date=([0-9]{4}-[0-9]{2}-[0-9]{2})',
                1
            ) AS DATE
        )
        """,
        [str(root / "release_date=*" / "events.parquet")],
    ).fetchone()[0]
    if mismatch:
        raise ReportValidationError(
            f"Event observed_at does not match partition date in {root}"
        )
    return rows


def _validate_episodes(
    connection: duckdb.DuckDBPyConnection,
    path: Path,
    key_columns: tuple[str, ...],
    study_end: date,
) -> None:
    keys = ", ".join(key_columns)
    duplicate_groups = connection.execute(
        f"""
        SELECT count(*)
        FROM (
            SELECT {keys}, count(*) AS occurrences
            FROM read_parquet(?)
            GROUP BY {keys}
            HAVING occurrences > 1
        )
        """,
        [str(path)],
    ).fetchone()[0]
    if duplicate_groups:
        raise ReportValidationError(f"Episode key is not unique: {path}")

    placeholders = ", ".join("?" for _ in _OUTCOMES)
    invalid_duration, invalid_outcome, beyond_study_end = connection.execute(
        f"""
        SELECT
            count(*) FILTER (
                WHERE started_at IS NULL
                   OR ended_at IS NULL
                   OR ended_at < started_at
                   OR duration_days != date_diff('day', started_at, ended_at)
            ) AS invalid_duration,
            count(*) FILTER (
                WHERE outcome IS NULL OR outcome NOT IN ({placeholders})
            ) AS invalid_outcome,
            count(*) FILTER (WHERE ended_at > ?) AS beyond_study_end
        FROM read_parquet(?)
        """,
        [*_OUTCOMES, study_end, str(path)],
    ).fetchone()
    if invalid_duration:
        raise ReportValidationError(f"Invalid episode duration in {path}")
    if invalid_outcome:
        raise ReportValidationError(f"Invalid episode outcome in {path}")
    if beyond_study_end:
        raise ReportValidationError(f"Episode exceeds study end in {path}")


def _comparison_split_date(table: pa.Table) -> date:
    rows = sorted(
        table.select(["period_start", "period_end"]).to_pylist(),
        key=lambda row: row["period_start"],
    )
    if not rows:
        raise ReportValidationError("Quarterly metrics table is empty")
    complete = [
        row for row in rows if (row["period_end"] - row["period_start"]).days >= 80
    ]
    if len(complete) < 2:
        raise ReportValidationError("At least two complete quarters are required")
    preferred = date(2024, 4, 6)
    has_before = any(row["period_end"] <= preferred for row in complete)
    has_after = any(row["period_start"] >= preferred for row in complete)
    if has_before and has_after:
        return preferred
    return complete[len(complete) // 2]["period_start"]


def audit_full_backfill(
    data_root: Path,
    results_root: Path,
    *,
    expected_vcf_releases: int | None = None,
    expected_scv_releases: int | None = None,
) -> FullBackfillAudit:
    """Independently validate the retained full-backfill artifacts."""
    data_root = Path(data_root)
    results_root = Path(results_root)
    manifest_path = data_root / "manifest" / "releases.parquet"
    ledger_path = data_root / "clinvar.duckdb"
    try:
        manifest = pq.read_table(manifest_path)
    except (FileNotFoundError, OSError, pa.ArrowException) as error:
        raise ReportValidationError(
            f"Unreadable release manifest: {manifest_path}"
        ) from error
    manifest_rows = _rows(
        manifest,
        {"release_date", "kind", "url", "size_bytes", "md5"},
    )
    urls = [str(row["url"]) for row in manifest_rows]
    if len(urls) != len(set(urls)):
        raise ReportValidationError("Manifest URLs are not unique")
    kinds = {"vcf", "submission_summary"}
    if {str(row["kind"]) for row in manifest_rows} != kinds:
        raise ReportValidationError("Manifest source kinds are invalid")
    vcf_releases = sum(row["kind"] == "vcf" for row in manifest_rows)
    scv_releases = sum(row["kind"] == "submission_summary" for row in manifest_rows)
    if expected_vcf_releases is not None and vcf_releases != expected_vcf_releases:
        raise ReportValidationError(
            f"VCF release count mismatch: {vcf_releases} != {expected_vcf_releases}"
        )
    if expected_scv_releases is not None and scv_releases != expected_scv_releases:
        raise ReportValidationError(
            f"SCV release count mismatch: {scv_releases} != {expected_scv_releases}"
        )

    try:
        connection = duckdb.connect(str(ledger_path), read_only=True)
    except duckdb.Error as error:
        raise ReportValidationError(
            f"Unreadable processing ledger: {ledger_path}"
        ) from error
    try:
        ledger_rows = connection.execute(
            """
            SELECT kind, url, release_date, status, source_md5, official_md5
            FROM pipeline_runs
            """
        ).fetchall()
        manifest_keys = {
            (str(row["kind"]), str(row["url"]), row["release_date"])
            for row in manifest_rows
        }
        ledger_keys = {(row[0], row[1], row[2]) for row in ledger_rows}
        if manifest_keys != ledger_keys:
            raise ReportValidationError("manifest/ledger release keys do not match")
        if any(row[3] != "completed" for row in ledger_rows):
            raise ReportValidationError(
                "Processing ledger contains incomplete releases"
            )
        if any(row[4] is None for row in ledger_rows):
            raise ReportValidationError(
                "Processing ledger is missing a computed checksum"
            )
        with_official = [row for row in ledger_rows if row[5] is not None]
        checksum_matches = sum(row[4] == row[5] for row in with_official)
        if checksum_matches != len(with_official):
            raise ReportValidationError(
                "Computed checksum does not match official checksum"
            )

        vcv_state_rows = _validate_latest_state(
            connection,
            data_root / "state" / "vcv_latest.parquet",
            "variation_id",
        )
        scv_state_rows = _validate_latest_state(
            connection,
            data_root / "state" / "scv_latest.parquet",
            "scv_accession",
        )
        vcv_event_rows = _validate_event_store(
            connection,
            data_root / "events" / "vcv",
            vcf_releases,
        )
        scv_event_rows = _validate_event_store(
            connection,
            data_root / "events" / "scv",
            scv_releases,
        )

        result_rows = _validate_result_files(results_root)
        mechanism_summary = summarize_mechanisms(
            pq.read_table(results_root / "event_mechanisms.parquet")
        )
        transition_summary = summarize_transitions(
            pq.read_table(results_root / "transition_matrix.parquet")
        )
        sampling_summary = summarize_sampling(
            pq.read_table(results_root / "sampling_bias.parquet")
        )
        attribution_total = result_rows["attributions.parquet"]
        totals = {
            attribution_total,
            sum(row.count for row in mechanism_summary),
            transition_summary.total,
            *(row.true_transitions for row in sampling_summary.values()),
        }
        if len(totals) != 1:
            raise ReportValidationError(
                f"Cross-table transition totals do not match: {sorted(totals)}"
            )

        quarterly_table = pq.read_table(results_root / "renovo_quarterly.parquet")
        quarterly = summarize_quarters(
            quarterly_table,
            split_date=_comparison_split_date(quarterly_table),
        )
        study_end = quarterly.points[-1].period_end
        _validate_episodes(
            connection,
            results_root / "vus_time_to_resolution.parquet",
            ("variation_id", "started_at"),
            study_end,
        )
        _validate_episodes(
            connection,
            results_root / "scv_vus_time_to_resolution.parquet",
            ("scv_accession", "started_at"),
            study_end,
        )
    finally:
        connection.close()

    return FullBackfillAudit(
        ok=True,
        vcf_releases=vcf_releases,
        scv_releases=scv_releases,
        completed_ledger_rows=len(ledger_rows),
        available_official_checksums=len(with_official),
        matching_official_checksums=checksum_matches,
        cumulative_source_bytes=sum(
            int(row["size_bytes"] or 0) for row in manifest_rows
        ),
        retained_data_bytes=_file_bytes(data_root),
        results_bytes=_file_bytes(results_root),
        vcv_state_rows=vcv_state_rows,
        scv_state_rows=scv_state_rows,
        vcv_event_rows=vcv_event_rows,
        scv_event_rows=scv_event_rows,
        result_rows=result_rows,
        cross_table_transition_total=attribution_total,
        study_start=quarterly.points[0].period_start,
        study_end=study_end,
    )


def _summarize_outcome_path(path: Path) -> OutcomeSummary:
    connection = duckdb.connect()
    try:
        rows = connection.execute(
            """
            SELECT outcome, count(*) AS episodes, median(duration_days) AS median_days
            FROM read_parquet(?)
            GROUP BY outcome
            """,
            [str(path)],
        ).fetchall()
    finally:
        connection.close()
    counts = {outcome: 0 for outcome in _OUTCOMES}
    medians: dict[str, float] = {}
    for outcome, count, median_days in rows:
        if outcome not in counts:
            raise ReportValidationError(f"Unknown episode outcome: {outcome}")
        counts[outcome] = int(count)
        medians[outcome] = float(median_days)
    total = sum(counts.values())
    if total == 0:
        raise ReportValidationError(f"Episode table is empty: {path}")
    return OutcomeSummary(
        total=total,
        counts=counts,
        shares={outcome: count / total for outcome, count in counts.items()},
        median_days=medians,
    )


def build_full_report_summary(
    data_root: Path,
    results_root: Path,
    *,
    expected_vcf_releases: int | None = None,
    expected_scv_releases: int | None = None,
) -> FullReportSummary:
    """Audit the run and derive the single summary used by report figures."""
    data_root = Path(data_root)
    results_root = Path(results_root)
    audit = audit_full_backfill(
        data_root,
        results_root,
        expected_vcf_releases=expected_vcf_releases,
        expected_scv_releases=expected_scv_releases,
    )
    quarterly_table = pq.read_table(results_root / "renovo_quarterly.parquet")
    return FullReportSummary(
        audit=audit,
        quarterly=summarize_quarters(
            quarterly_table,
            split_date=_comparison_split_date(quarterly_table),
        ),
        sampling=summarize_sampling(
            pq.read_table(results_root / "sampling_bias.parquet")
        ),
        mechanisms=summarize_mechanisms(
            pq.read_table(results_root / "event_mechanisms.parquet")
        ),
        transitions=summarize_transitions(
            pq.read_table(results_root / "transition_matrix.parquet")
        ),
        vcv_outcomes=_summarize_outcome_path(
            results_root / "vus_time_to_resolution.parquet"
        ),
        scv_outcomes=_summarize_outcome_path(
            results_root / "scv_vus_time_to_resolution.parquet"
        ),
    )


def render_report_figures(
    summary: FullReportSummary,
    output_dir: Path,
) -> dict[str, Path]:
    """Render the five canonical report figures as deterministic SVG files."""
    from .report_figures import render_report_figures as render

    return render(summary, Path(output_dir))
