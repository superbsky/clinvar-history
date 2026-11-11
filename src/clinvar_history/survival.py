"""First-episode Aalen–Johansen estimates on observed release-time coordinates."""

import bisect
import math
import random
from collections import Counter, defaultdict
from dataclasses import dataclass
from datetime import date
from pathlib import Path

import duckdb

CAUSES = ("benign", "pathogenic", "conflict", "other_exit")


@dataclass(frozen=True, slots=True)
class Observation:
    cluster: str
    duration: float
    outcome: str
    stratum: str = "all"


def _counts(rows, weights=None):
    counts = defaultdict(Counter)
    for row in rows:
        if (
            not math.isfinite(row.duration)
            or row.duration < 0
            or row.outcome not in (*CAUSES, "censored")
        ):
            raise ValueError("Invalid duration or competing outcome")
        count = 1 if weights is None else weights[row.cluster]
        if count:
            counts[row.duration][row.outcome] += count
    return counts


def _estimate(counts, horizons):
    at_risk = sum(sum(c.values()) for c in counts.values())
    if not at_risk:
        raise ValueError("No observations")
    survival = 1.0
    incidence = dict.fromkeys(CAUSES, 0.0)
    times = sorted(counts)
    index = 0
    result = []
    for horizon in sorted(set(horizons)):
        risk_at_horizon = at_risk
        while index < len(times) and times[index] <= horizon:
            time = times[index]
            group = counts[time]
            risk_at_horizon = at_risk
            events = sum(group[cause] for cause in CAUSES)
            for cause in CAUSES:
                incidence[cause] += survival * group[cause] / at_risk
            survival *= 1 - events / at_risk
            # Events occur before censors at tied observation times.
            at_risk -= sum(group.values())
            index += 1
            if time < horizon:
                risk_at_horizon = at_risk
        supported = horizon <= times[-1] or survival <= 1e-12
        result.append(
            {
                "day": horizon,
                "at_risk": risk_at_horizon,
                "supported": supported,
                "survival": survival if supported else None,
                **{c: incidence[c] if supported else None for c in CAUSES},
            }
        )
    return result


def cumulative_incidence(rows, horizons):
    """Discrete tied-time AJ; retain censoring in the risk set at an event tie."""
    if any(not math.isfinite(h) or h < 0 for h in horizons):
        raise ValueError("Horizons must be finite and non-negative")
    return _estimate(_counts(rows), horizons)


def _quantile(values, probability):
    values = sorted(values)
    position = (len(values) - 1) * probability
    lower = int(position)
    upper = min(lower + 1, len(values) - 1)
    return values[lower] + (values[upper] - values[lower]) * (position - lower)


def bootstrap_incidence(rows, horizons, *, replicates=200, seed=20261111):
    """Pointwise percentile intervals from independent variant-cluster resampling."""
    if replicates < 20:
        raise ValueError("Use at least 20 bootstrap replicates")
    point = cumulative_incidence(rows, horizons)
    clusters = sorted({row.cluster for row in rows})
    rng = random.Random(seed)
    draws = []
    for _ in range(replicates):
        weights = Counter(rng.choices(clusters, k=len(clusters)))
        draws.append(_estimate(_counts(rows, weights), horizons))
    for i, estimate in enumerate(point):
        for cause in CAUSES:
            values = [draw[i][cause] for draw in draws if draw[i][cause] is not None]
            estimate[f"{cause}_bootstrap_supported"] = len(values)
            enough = estimate[cause] is not None and len(values) >= math.ceil(
                replicates * 0.95
            )
            estimate[f"{cause}_lower"] = _quantile(values, 0.025) if enough else None
            estimate[f"{cause}_upper"] = _quantile(values, 0.975) if enough else None
    return point


def load_first_episodes(
    path: Path,
    *,
    layer: str,
    cohort_start: date,
    cohort_end: date,
    study_end: date,
    max_clusters: int,
    seed: int,
    stratum: str = "review_status",
) -> list[dict]:
    """Select first incident episodes before cohort filtering, with stable sampling."""
    if layer not in {"vcv", "scv"} or stratum not in {
        "review_status",
        "variant_type",
        "collection_method",
    }:
        raise ValueError("Unsupported layer or stratum")
    if not cohort_start <= cohort_end <= study_end or max_clusters < 1:
        raise ValueError("Invalid cohort dates or cluster budget")
    key = "variation_id" if layer == "vcv" else "scv_accession"
    with duckdb.connect() as connection:
        columns = {
            r[0]
            for r in connection.execute(
                "DESCRIBE SELECT * FROM read_parquet(?)", [str(path)]
            ).fetchall()
        }
        if stratum not in columns:
            raise ValueError(f"Stratum {stratum} is absent")
        invalid = connection.execute(
            "SELECT count(*) FROM read_parquet(?) WHERE started_at IS NULL OR ended_at IS NULL OR ended_at < started_at OR duration_days != date_diff('day', started_at, ended_at) OR outcome NOT IN ('benign','pathogenic','conflict','other_exit','censored') OR outcome IS NULL",
            [str(path)],
        ).fetchone()[0]
        if invalid:
            raise ValueError("Invalid input episodes")
        duplicates = connection.execute(
            f"SELECT count(*) FROM (SELECT {key}, started_at FROM read_parquet(?) GROUP BY {key}, started_at HAVING count(*) > 1)",
            [str(path)],
        ).fetchone()[0]
        if duplicates:
            raise ValueError("Duplicate episode keys")
        query = f"""
            WITH first AS (
                SELECT *, row_number() OVER (PARTITION BY {key} ORDER BY started_at) AS episode_rank
                FROM read_parquet(?)
            ), cohort AS (
                SELECT * FROM first WHERE episode_rank = 1 AND started_at BETWEEN ? AND ?
            ), clusters AS (
                SELECT DISTINCT variation_id FROM cohort
                ORDER BY sha256(CAST(variation_id AS VARCHAR) || ?) LIMIT ?
            )
            SELECT CAST({key} AS VARCHAR) AS unit, variation_id, started_at,
                   LEAST(ended_at, ?::DATE) AS ended_at,
                   CASE WHEN ended_at > ?::DATE THEN 'censored' ELSE outcome END AS outcome,
                   COALESCE({stratum}, 'missing') AS stratum
            FROM cohort JOIN clusters USING (variation_id) ORDER BY variation_id, unit
        """
        cursor = connection.execute(
            query,
            [
                str(path),
                cohort_start,
                cohort_end,
                str(seed),
                max_clusters,
                study_end,
                study_end,
            ],
        )
        names = [c[0] for c in cursor.description]
        return [dict(zip(names, row, strict=True)) for row in cursor.fetchall()]


def observations(
    episodes, *, release_dates, timing="observed", removal_as_exit=False, study_end=None
):
    """Midpoint sensitivity changes observation coordinates, not biological truth."""
    dates = sorted(set(release_dates))
    if timing not in {"observed", "midpoint"} or not dates:
        raise ValueError("Invalid timing model or release grid")

    def midpoint(value):
        i = bisect.bisect_left(dates, value)
        if i == len(dates) or dates[i] != value or i == 0:
            raise ValueError("Incident boundary not covered by release grid")
        return (value.toordinal() + dates[i - 1].toordinal()) / 2

    result = []
    study_end = study_end or dates[-1]
    for episode in episodes:
        start, end = episode["started_at"], episode["ended_at"]
        outcome = episode["outcome"]
        if removal_as_exit and outcome == "censored" and end < study_end:
            outcome = "other_exit"
        start_coordinate = (
            midpoint(start) if timing == "midpoint" else start.toordinal()
        )
        end_coordinate = (
            midpoint(end)
            if timing == "midpoint" and outcome != "censored"
            else end.toordinal()
        )
        result.append(
            Observation(
                str(episode["variation_id"]),
                max(0, end_coordinate - start_coordinate),
                outcome,
                episode["stratum"],
            )
        )
    return result


def snapshot_sensitivity(
    root, episodes, *, release_dates, study_end, cohort_start, cohort_end, months
):
    """Actually reconstruct coarser VCV snapshots for the paired pilot variants."""
    from .analysis_core import month_grid

    ids = sorted({e["variation_id"] for e in episodes})
    strata = {e["variation_id"]: e["stratum"] for e in episodes}
    with duckdb.connect() as connection:
        source = connection.execute(
            "SELECT variation_id, event_type, observed_at, new_classification FROM read_parquet(?) "
            "WHERE variation_id IN (SELECT unnest(?)) AND event_type IN ('variant_added', 'variant_removed', 'classification_changed') "
            "ORDER BY variation_id, observed_at, event_type",
            [str(root / "events/vcv/release_date=*/events.parquet"), ids],
        ).fetchall()
    histories = defaultdict(list)
    for vid, event_type, observed_at, classification in source:
        histories[vid].append(
            (observed_at, None if event_type == "variant_removed" else classification)
        )
    grid = month_grid(min(release_dates), study_end, months)
    result = []
    for vid in ids:
        history = histories[vid]
        index, state, previous, started = 0, None, None, None
        for i, observed_at in enumerate(grid):
            while index < len(history) and history[index][0] <= observed_at:
                state = history[index][1]
                index += 1
            if i and previous != "vus" and state == "vus":
                started = observed_at
            if started is not None and state != "vus":
                if cohort_start <= started <= cohort_end:
                    outcome = (
                        "benign"
                        if state in {"benign", "likely_benign"}
                        else "pathogenic"
                        if state in {"pathogenic", "likely_pathogenic"}
                        else "conflict"
                        if state == "conflict"
                        else "censored"
                        if state is None
                        else "other_exit"
                    )
                    result.append(
                        Observation(
                            str(vid), (observed_at - started).days, outcome, strata[vid]
                        )
                    )
                break  # first incident observed episode, even if outside the cohort
            previous = state
        else:
            if started is not None and cohort_start <= started <= cohort_end:
                result.append(
                    Observation(
                        str(vid), (study_end - started).days, "censored", strata[vid]
                    )
                )
    return result
