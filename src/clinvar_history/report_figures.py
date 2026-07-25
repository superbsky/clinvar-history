"""Deterministic, dependency-free SVG figures for the full-backfill report."""

from __future__ import annotations

import math
from collections.abc import Iterable
from datetime import date
from html import escape
from pathlib import Path
from typing import TYPE_CHECKING

if TYPE_CHECKING:
    from .reporting import FullReportSummary, OutcomeSummary


_BLUE = "#0072B2"
_SKY = "#56B4E9"
_GREEN = "#009E73"
_ORANGE = "#E69F00"
_VERMILION = "#D55E00"
_PURPLE = "#CC79A7"
_YELLOW = "#F0E442"
_INK = "#17212B"
_MUTED = "#5B6775"
_GRID = "#D8DEE6"
_LIGHT = "#F4F7FA"
_WHITE = "#FFFFFF"

_STYLE = """
text { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Helvetica, Arial, sans-serif; }
.title { font-size: 30px; font-weight: 700; letter-spacing: -0.3px; }
.subtitle { font-size: 15px; fill: #5B6775; }
.panel-title { font-size: 17px; font-weight: 700; }
.axis { font-size: 12px; fill: #5B6775; }
.label { font-size: 14px; }
.value { font-size: 13px; font-weight: 650; }
.small { font-size: 11px; fill: #5B6775; }
.note { font-size: 13px; fill: #5B6775; }
.matrix-column-label { font-size: 12px; fill: #5B6775; }
.grid { stroke: #D8DEE6; stroke-width: 1; }
.axis-line { stroke: #8A96A3; stroke-width: 1; }
""".strip()


def _number(value: float) -> str:
    if isinstance(value, int):
        return str(value)
    rounded = f"{value:.2f}".rstrip("0").rstrip(".")
    return "0" if rounded == "-0" else rounded


def _attributes(attributes: dict[str, object]) -> str:
    rendered = []
    for name, value in attributes.items():
        if value is None:
            continue
        name = "class" if name == "class_" else name.replace("_", "-")
        rendered.append(
            f'{name}="{escape(_number(value) if isinstance(value, (int, float)) else str(value), quote=True)}"'
        )
    return " " + " ".join(rendered) if rendered else ""


class _Svg:
    def __init__(self, width: int, height: int, slug: str, title: str, desc: str):
        self.width = width
        self.height = height
        self.parts = [
            '<?xml version="1.0" encoding="UTF-8"?>',
            (
                f'<svg xmlns="http://www.w3.org/2000/svg" '
                f'viewBox="0 0 {width} {height}" role="img" '
                f'id="{slug}" aria-labelledby="{slug}-title {slug}-desc">'
            ),
            f'<title id="{slug}-title">{escape(title)}</title>',
            f'<desc id="{slug}-desc">{escape(desc)}</desc>',
            f"<style>{_STYLE}</style>",
            f'<rect width="{width}" height="{height}" fill="{_WHITE}"/>',
        ]

    def element(self, tag: str, **attributes: object) -> None:
        self.parts.append(f"<{tag}{_attributes(attributes)}/>")

    def text(self, x: float, y: float, value: object, **attributes: object) -> None:
        attrs = {"x": x, "y": y, "fill": _INK, **attributes}
        self.parts.append(f"<text{_attributes(attrs)}>{escape(str(value))}</text>")

    def finish(self) -> str:
        return "\n".join([*self.parts, "</svg>", ""])


def _header(svg: _Svg, title: str, subtitle: str) -> None:
    svg.text(56, 48, title, class_="title")
    svg.text(56, 75, subtitle, class_="subtitle")


def _footer(svg: _Svg, text: str) -> None:
    svg.text(56, svg.height - 24, text, class_="small")


def _linear(value: float, low: float, high: float, start: float, end: float) -> float:
    if math.isclose(low, high):
        return (start + end) / 2
    return start + (value - low) / (high - low) * (end - start)


def _nice_max(value: float) -> float:
    if value <= 0:
        return 1.0
    exponent = math.floor(math.log10(value))
    fraction = value / 10**exponent
    for step in (1, 2, 2.5, 5, 10):
        if fraction <= step:
            return step * 10**exponent
    return 10 ** (exponent + 1)


def _compact_count(value: float) -> str:
    if abs(value) >= 1_000_000:
        return f"{value / 1_000_000:.1f}M".replace(".0M", "M")
    if abs(value) >= 1_000:
        return f"{value / 1_000:.0f}k"
    return f"{value:.0f}"


def _count(value: int) -> str:
    return f"{value:,}"


def _percent(value: float, digits: int = 2) -> str:
    return f"{value * 100:.{digits}f}%"


def _draw_y_grid(
    svg: _Svg,
    x_left: float,
    x_right: float,
    y_top: float,
    y_bottom: float,
    maximum: float,
    formatter,
    ticks: int = 4,
) -> None:
    for index in range(ticks + 1):
        value = maximum * index / ticks
        y = _linear(value, 0, maximum, y_bottom, y_top)
        svg.element("line", x1=x_left, y1=y, x2=x_right, y2=y, class_="grid")
        svg.text(x_left - 12, y + 4, formatter(value), class_="axis", text_anchor="end")


def _date_position(
    value: date, start: date, end: date, x_left: float, x_right: float
) -> float:
    return _linear(
        value.toordinal(), start.toordinal(), end.toordinal(), x_left, x_right
    )


def _draw_year_axis(
    svg: _Svg,
    start: date,
    end: date,
    x_left: float,
    x_right: float,
    y: float,
) -> None:
    svg.element("line", x1=x_left, y1=y, x2=x_right, y2=y, class_="axis-line")
    for year in range(start.year, end.year + 1):
        value = max(start, min(end, date(year, 1, 1)))
        if value < start or value > end:
            continue
        x = _date_position(value, start, end, x_left, x_right)
        svg.element("line", x1=x, y1=y, x2=x, y2=y + 6, class_="axis-line")
        svg.text(x, y + 22, year, class_="axis", text_anchor="middle")


def _draw_series(
    svg: _Svg,
    coordinates: Iterable[tuple[float, float]],
    color: str,
    *,
    width: float = 3,
) -> list[tuple[float, float]]:
    points = list(coordinates)
    svg.element(
        "polyline",
        points=" ".join(f"{_number(x)},{_number(y)}" for x, y in points),
        fill="none",
        stroke=color,
        stroke_width=width,
        stroke_linejoin="round",
        stroke_linecap="round",
    )
    for x, y in points:
        svg.element(
            "circle", cx=x, cy=y, r=3.5, fill=color, stroke=_WHITE, stroke_width=1
        )
    return points


def _trend_figure(summary: FullReportSummary) -> str:
    quarterly = summary.quarterly
    svg = _Svg(
        1200,
        920,
        "vus-quarterly-trend",
        "VUS inventory and quarterly dynamics",
        (
            "Three aligned time-series panels show aggregate VUS inventory, "
            "quarterly discovery rate, and quarterly reclassification rate from "
            f"{quarterly.points[0].period_start} through {quarterly.points[-1].period_end}."
        ),
    )
    _header(
        svg,
        "VUS inventory and quarterly dynamics",
        "Inventory rose steadily while discovery remained much larger than reclassification.",
    )
    x_left, x_right = 115.0, 1140.0
    start = quarterly.points[0].period_start
    end = quarterly.points[-1].period_end
    panels = [
        (110.0, 340.0, "A  Aggregate VUS inventory", _BLUE),
        (405.0, 600.0, "B  Quarterly discovery rate", _ORANGE),
        (665.0, 820.0, "C  Quarterly reclassification rate", _PURPLE),
    ]

    inventory_values = [
        quarterly.initial_vus,
        *(point.vus_current for point in quarterly.points),
    ]
    inventory_dates = [start, *(point.period_end for point in quarterly.points)]
    inventory_max = _nice_max(max(inventory_values) * 1.02)
    top, bottom, title, color = panels[0]
    svg.text(56, top - 18, title, class_="panel-title")
    _draw_y_grid(svg, x_left, x_right, top, bottom, inventory_max, _compact_count)
    inventory_points = _draw_series(
        svg,
        (
            (
                _date_position(value_date, start, end, x_left, x_right),
                _linear(value, 0, inventory_max, bottom, top),
            )
            for value_date, value in zip(inventory_dates, inventory_values, strict=True)
        ),
        color,
    )
    svg.text(
        inventory_points[0][0] + 7,
        inventory_points[0][1] - 9,
        _count(quarterly.initial_vus),
        class_="value",
        fill=color,
    )
    svg.text(
        inventory_points[-1][0] - 5,
        inventory_points[-1][1] - 11,
        _count(quarterly.final_vus),
        class_="value",
        fill=color,
        text_anchor="end",
    )

    for panel_index, (attribute, label) in enumerate(
        (("discovery_rate", panels[1]), ("reclassification_rate", panels[2])),
        start=1,
    ):
        top, bottom, title, color = label
        values = [getattr(point, attribute) for point in quarterly.points]
        maximum = _nice_max(max(values) * 1.05)
        svg.text(56, top - 18, title, class_="panel-title")
        _draw_y_grid(
            svg,
            x_left,
            x_right,
            top,
            bottom,
            maximum,
            lambda value: f"{value * 100:.1f}%",
        )
        coordinates = [
            (
                _date_position(point.period_end, start, end, x_left, x_right),
                _linear(getattr(point, attribute), 0, maximum, bottom, top),
            )
            for point in quarterly.points
        ]
        _draw_series(svg, coordinates, color, width=2.5)
        partial = quarterly.points[-1]
        if not partial.complete:
            x, y = coordinates[-1]
            svg.element(
                "circle", cx=x, cy=y, r=5, fill=_WHITE, stroke=color, stroke_width=2
            )
        if panel_index == 2:
            svg.text(
                x_right,
                top - 18,
                "Open circle: partial final interval",
                class_="small",
                text_anchor="end",
            )

    split_x = _date_position(quarterly.split_date, start, end, x_left, x_right)
    for top, bottom, _, _ in panels:
        svg.element(
            "line",
            x1=split_x,
            y1=top,
            x2=split_x,
            y2=bottom,
            stroke=_MUTED,
            stroke_width=1.5,
            stroke_dasharray="5 5",
        )
    svg.text(
        split_x + 7,
        panels[0][0] + 18,
        quarterly.split_date.isoformat(),
        class_="small",
    )
    _draw_year_axis(svg, start, end, x_left, x_right, panels[-1][1])
    _footer(
        svg,
        "Source: results/full-v1/renovo_quarterly.parquet · Final interval is partial.",
    )
    return svg.finish()


def _temporal_resolution_figure(summary: FullReportSummary) -> str:
    svg = _Svg(
        1200,
        660,
        "temporal-resolution",
        "Information loss at coarser temporal resolution",
        "Recall, missed canonical transitions, missed conflict episodes, and mean date error are compared for 7-day, monthly, and quarterly endpoint grids.",
    )
    _header(
        svg,
        "Information loss at coarser temporal resolution",
        "Release-level history preserves transient and recurrent episodes hidden between endpoints.",
    )
    order = [
        name for name in ("weekly", "monthly", "quarterly") if name in summary.sampling
    ]
    labels = {"weekly": "7-day grid", "monthly": "Monthly", "quarterly": "Quarterly"}
    maximum_missed = _nice_max(
        max(summary.sampling[name].missed_transitions for name in order) * 1.03
    )
    svg.text(250, 112, "Transition recall", class_="panel-title", text_anchor="middle")
    svg.text(
        870,
        112,
        "Episodes missed between endpoints",
        class_="panel-title",
        text_anchor="middle",
    )
    for tick in range(5):
        share = tick / 4
        x = _linear(share, 0, 1, 160, 520)
        svg.element("line", x1=x, y1=130, x2=x, y2=540, class_="grid")
        svg.text(x, 562, f"{share * 100:.0f}%", class_="axis", text_anchor="middle")
    for tick in range(5):
        value = maximum_missed * tick / 4
        x = _linear(value, 0, maximum_missed, 710, 1130)
        svg.element("line", x1=x, y1=130, x2=x, y2=540, class_="grid")
        svg.text(x, 562, _compact_count(value), class_="axis", text_anchor="middle")

    for index, cadence in enumerate(order):
        row = summary.sampling[cadence]
        y = 185 + index * 135
        svg.text(56, y + 6, labels[cadence], class_="label")
        recall_width = 360 * row.recall
        svg.element("rect", x=160, y=y - 16, width=360, height=28, rx=4, fill=_LIGHT)
        svg.element(
            "rect", x=160, y=y - 16, width=recall_width, height=28, rx=4, fill=_BLUE
        )
        svg.text(
            160 + recall_width - 7,
            y + 4,
            _percent(row.recall),
            class_="value",
            fill=_WHITE,
            text_anchor="end",
        )

        transition_width = _linear(row.missed_transitions, 0, maximum_missed, 0, 420)
        conflict_width = _linear(
            row.missed_conflict_episodes, 0, maximum_missed, 0, 420
        )
        svg.element(
            "rect",
            x=710,
            y=y - 24,
            width=transition_width,
            height=20,
            rx=3,
            fill=_ORANGE,
        )
        svg.element(
            "rect", x=710, y=y + 4, width=conflict_width, height=20, rx=3, fill=_PURPLE
        )
        svg.text(
            716 + transition_width,
            y - 9,
            _count(row.missed_transitions),
            class_="value",
        )
        svg.text(
            716 + conflict_width,
            y + 20,
            _count(row.missed_conflict_episodes),
            class_="value",
        )
        svg.text(
            710,
            y + 48,
            f"Mean event-date error: {row.mean_event_date_error_days:.2f} days",
            class_="small",
        )

    svg.element("rect", x=710, y=590, width=18, height=12, fill=_ORANGE)
    svg.text(736, 601, "Missed transitions", class_="small")
    svg.element("rect", x=870, y=590, width=18, height=12, fill=_PURPLE)
    svg.text(896, 601, "Missed conflict episodes", class_="small")
    _footer(
        svg,
        "Source: results/full-v1/sampling_bias.parquet · Recall denominator: all canonical release-level transitions.",
    )
    return svg.finish()


_STATE_LABELS = {
    "benign": "Benign",
    "likely_benign": "Likely benign",
    "vus": "VUS",
    "likely_pathogenic": "Likely pathogenic",
    "pathogenic": "Pathogenic",
    "conflict": "Conflict",
    "other": "Other",
    "missing": "Missing",
}


def _blend(low: str, high: str, amount: float) -> str:
    low_rgb = tuple(int(low[index : index + 2], 16) for index in (1, 3, 5))
    high_rgb = tuple(int(high[index : index + 2], 16) for index in (1, 3, 5))
    rgb = tuple(
        round(start + (end - start) * amount)
        for start, end in zip(low_rgb, high_rgb, strict=True)
    )
    return "#" + "".join(f"{channel:02X}" for channel in rgb)


def _transition_matrix_figure(summary: FullReportSummary) -> str:
    transition = summary.transitions
    svg = _Svg(
        1120,
        1050,
        "transition-matrix",
        "Directed canonical classification transitions",
        "An eight-by-eight matrix shows counts for each old-to-new canonical classification pair; color intensity uses a logarithmic scale.",
    )
    _header(
        svg,
        "Directed canonical classification transitions",
        f"All {_count(transition.total)} transitions; rows are old states and columns are new states.",
    )
    states = transition.states
    cell = 82
    x0, y0 = 330, 230
    maximum = max(transition.counts.values(), default=1)
    svg.text(
        56,
        y0 + cell * 4,
        "Old classification",
        class_="panel-title",
        transform=f"rotate(-90 56 {y0 + cell * 4})",
        text_anchor="middle",
    )
    svg.text(
        x0 + cell * 4,
        108,
        "New classification",
        class_="panel-title",
        text_anchor="middle",
    )
    for row_index, old_state in enumerate(states):
        y = y0 + row_index * cell
        svg.text(
            x0 - 18,
            y + cell / 2 + 4,
            _STATE_LABELS[old_state],
            class_="axis",
            text_anchor="end",
        )
        for column, new_state in enumerate(states):
            count = transition.counts.get((old_state, new_state), 0)
            intensity = math.log1p(count) / math.log1p(maximum) if maximum else 0
            fill = _blend("#F4F7FA", _BLUE, intensity)
            x = x0 + column * cell
            svg.element(
                "rect",
                x=x,
                y=y,
                width=cell - 2,
                height=cell - 2,
                rx=3,
                fill=fill,
                stroke=_WHITE,
                stroke_width=1,
                class_="matrix-cell",
            )
            if count:
                text_color = _WHITE if intensity > 0.58 else _INK
                svg.text(
                    x + cell / 2 - 1,
                    y + cell / 2 + 4,
                    _count(count),
                    class_="value",
                    fill=text_color,
                    text_anchor="middle",
                )

    for column, state in enumerate(states):
        x = x0 + column * cell + 15
        y = y0 - 13
        svg.text(
            x,
            y,
            _STATE_LABELS[state],
            class_="matrix-column-label",
            text_anchor="start",
            transform=f"rotate(-48 {x} {y})",
        )

    legend_y = y0 + cell * len(states) + 38
    svg.text(x0, legend_y, "Log-scaled count", class_="small")
    for index in range(8):
        amount = index / 7
        svg.element(
            "rect",
            x=x0 + 115 + index * 34,
            y=legend_y - 12,
            width=34,
            height=14,
            fill=_blend("#F4F7FA", _BLUE, amount),
        )
    svg.text(x0 + 105, legend_y + 22, "0", class_="small", text_anchor="middle")
    svg.text(
        x0 + 115 + 8 * 34,
        legend_y + 22,
        _count(maximum),
        class_="small",
        text_anchor="end",
    )
    _footer(
        svg,
        "Source: results/full-v1/transition_matrix.parquet · Cell labels are exact counts; only color is log-scaled.",
    )
    return svg.finish()


_MECHANISM_LABELS = {
    "scv_added": "SCV added",
    "unexplained": "Unexplained",
    "exact_matching_scv_reclassification": "Exact matching SCV reclassification",
    "multiple_candidates": "Multiple candidates",
    "review_precedence_changed": "Review precedence changed",
    "condition_changed": "Condition changed",
    "scv_removed": "SCV removed",
    "expert_panel_added": "Expert panel added",
}


def _mechanisms_figure(summary: FullReportSummary) -> str:
    svg = _Svg(
        1200,
        720,
        "candidate-mechanisms",
        "Conservative candidate mechanisms for aggregate transitions",
        "A descending horizontal bar chart shows counts and shares for eight interval-overlap attribution categories.",
    )
    _header(
        svg,
        "Conservative candidate mechanisms",
        "Categories are compatible temporal evidence, not proven causes.",
    )
    rows = summary.mechanisms
    maximum = _nice_max(max(row.count for row in rows) * 1.02)
    x0, x1 = 355.0, 1030.0
    for tick in range(6):
        value = maximum * tick / 5
        x = _linear(value, 0, maximum, x0, x1)
        svg.element("line", x1=x, y1=115, x2=x, y2=625, class_="grid")
        svg.text(x, 646, _compact_count(value), class_="axis", text_anchor="middle")
    colors = (_BLUE, _MUTED, _GREEN, _PURPLE, _SKY, _ORANGE, _VERMILION, _YELLOW)
    for index, row in enumerate(rows):
        y = 145 + index * 61
        label = _MECHANISM_LABELS.get(
            row.mechanism, row.mechanism.replace("_", " ").title()
        )
        width = _linear(row.count, 0, maximum, 0, x1 - x0)
        svg.text(x0 - 18, y + 17, label, class_="label", text_anchor="end")
        svg.element(
            "rect",
            x=x0,
            y=y,
            width=width,
            height=28,
            rx=4,
            fill=colors[index % len(colors)],
        )
        svg.text(
            x0 + width + 9,
            y + 19,
            f"{_count(row.count)}  ·  {_percent(row.share)}",
            class_="value",
        )
    _footer(
        svg,
        "Source: results/full-v1/event_mechanisms.parquet · Total = 479,525 canonical transitions.",
    )
    return svg.finish()


_OUTCOME_LABELS = {
    "censored": "Censored",
    "conflict": "Conflict",
    "benign": "Benign / likely benign",
    "pathogenic": "Pathogenic / likely pathogenic",
    "other_exit": "Other exit",
}
_OUTCOME_COLORS = {
    "censored": "#A7B0BA",
    "conflict": _ORANGE,
    "benign": _GREEN,
    "pathogenic": _VERMILION,
    "other_exit": _PURPLE,
}


def _stacked_bar(
    svg: _Svg,
    outcome: OutcomeSummary,
    order: tuple[str, ...],
    x: float,
    y: float,
    width: float,
    height: float,
    denominator: int,
) -> None:
    cursor = x
    if denominator == 0:
        svg.element("rect", x=x, y=y, width=width, height=height, rx=4, fill=_LIGHT)
        return
    for name in order:
        share = outcome.counts[name] / denominator
        segment = width * share
        svg.element(
            "rect",
            x=cursor,
            y=y,
            width=segment,
            height=height,
            fill=_OUTCOME_COLORS[name],
        )
        if share >= 0.075:
            svg.text(
                cursor + segment / 2,
                y + height / 2 + 5,
                _percent(share, 1),
                class_="value",
                fill=_WHITE,
                text_anchor="middle",
            )
        cursor += segment
    svg.element(
        "rect",
        x=x,
        y=y,
        width=width,
        height=height,
        rx=4,
        fill="none",
        stroke=_MUTED,
        stroke_width=1,
    )


def _median_text(outcome: OutcomeSummary) -> str:
    benign = outcome.median_days.get("benign")
    pathogenic = outcome.median_days.get("pathogenic")
    benign_text = "—" if benign is None else f"{benign:.0f} d"
    pathogenic_text = "—" if pathogenic is None else f"{pathogenic:.0f} d"
    return f"Median observed duration · benign {benign_text} · pathogenic {pathogenic_text}"


def _outcomes_figure(summary: FullReportSummary) -> str:
    svg = _Svg(
        1200,
        850,
        "vus-episode-outcomes",
        "Observed outcomes of post-bootstrap VUS episodes",
        "Aggregate VCV and submitter-level SCV VUS episodes are shown as all-episode compositions and as observed-exit compositions excluding right-censored episodes.",
    )
    _header(
        svg,
        "Observed outcomes of post-bootstrap VUS episodes",
        "The lower panel expands observed exits; raw fractions are not cumulative-incidence estimates.",
    )
    datasets = (
        ("Aggregate VCV", summary.vcv_outcomes),
        ("Submitter SCV", summary.scv_outcomes),
    )
    all_order = ("censored", "conflict", "benign", "pathogenic", "other_exit")
    exit_order = ("conflict", "benign", "pathogenic", "other_exit")
    x, width = 270.0, 830.0
    svg.text(56, 120, "A  All episodes", class_="panel-title")
    for index, (label, outcome) in enumerate(datasets):
        y = 160 + index * 115
        svg.text(x - 18, y + 25, label, class_="label", text_anchor="end")
        _stacked_bar(svg, outcome, all_order, x, y, width, 42, outcome.total)
        svg.text(
            x,
            y + 65,
            f"n = {_count(outcome.total)} · censored {_percent(outcome.shares['censored'])}",
            class_="small",
        )

    svg.text(56, 430, "B  Observed exits only", class_="panel-title")
    for index, (label, outcome) in enumerate(datasets):
        y = 470 + index * 125
        observed = outcome.total - outcome.counts["censored"]
        svg.text(x - 18, y + 25, label, class_="label", text_anchor="end")
        _stacked_bar(svg, outcome, exit_order, x, y, width, 42, observed)
        svg.text(
            x,
            y + 64,
            f"Observed exits n = {_count(observed)} · {_median_text(outcome)}",
            class_="small",
        )

    legend_y = 735
    legend_x = 155
    for index, name in enumerate(all_order):
        x_position = legend_x + index * 205
        svg.element(
            "rect",
            x=x_position,
            y=legend_y,
            width=18,
            height=14,
            fill=_OUTCOME_COLORS[name],
        )
        svg.text(x_position + 26, legend_y + 12, _OUTCOME_LABELS[name], class_="small")
    _footer(
        svg,
        "Source: results/full-v1/{vus_time_to_resolution,scv_vus_time_to_resolution}.parquet.",
    )
    return svg.finish()


def render_report_figures(
    summary: FullReportSummary,
    output_dir: Path,
) -> dict[str, Path]:
    """Render and atomically write every canonical report figure."""
    figures = {
        "vus-quarterly-trend.svg": _trend_figure(summary),
        "temporal-resolution.svg": _temporal_resolution_figure(summary),
        "transition-matrix.svg": _transition_matrix_figure(summary),
        "candidate-mechanisms.svg": _mechanisms_figure(summary),
        "vus-episode-outcomes.svg": _outcomes_figure(summary),
    }
    output_dir.mkdir(parents=True, exist_ok=True)
    paths: dict[str, Path] = {}
    for name, content in figures.items():
        path = output_dir / name
        temporary = path.with_suffix(path.suffix + ".tmp")
        temporary.write_text(content, encoding="utf-8")
        temporary.replace(path)
        paths[name] = path
    return paths
