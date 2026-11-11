"""Explicit comparison policies for the pinned Liang et al. public results."""

import re
from datetime import date

# Mirrors categorize_classification in the MIT-licensed source repository,
# commit dae628001d2b335068cef74e262e4431b288eb00. This is a comparison policy,
# not a recommendation to use substring classification in ingestion.
# Copyright (c) 2026 Yilun Liang and co-authors; MIT license reproduced in
# docs/third-party-liang-license.txt.


def liang_category(label):
    if not label:
        return "Other"
    text = label.lower().strip()
    if "conflicting" in text:
        return "Conflicting"
    if "uncertain" in text or "vus" in text:
        return "VUS"
    if "pathogenic/likely" in text or "pathogenic, likely" in text:
        return "P/LP"
    if "benign/likely" in text or "benign, likely" in text:
        return "B/LB"
    if "likely pathogenic" in text:
        return "LP"
    if "pathogenic" in text:
        return "P"
    if "likely benign" in text:
        return "LB"
    if "benign" in text:
        return "B"
    return "Other"


def family(category):
    if category in {"P", "LP", "P/LP"}:
        return "pathogenic"
    if category in {"B", "LB", "B/LB"}:
        return "benign"
    return {"VUS": "vus", "Conflicting": "conflict", "Other": "other"}.get(
        category, category
    )


def parse_history(value):
    rows = []
    for text in value.split(" | ") if value else []:
        m = re.fullmatch(r"(\d{4}-\d{2}):(.+?)→(\d{4}-\d{2}):(.+?)\[([^\]]+)\]", text)
        if m is None:
            raise ValueError(f"Unexpected history entry: {text}")
        rows.append(
            {
                "previous_month": date.fromisoformat(m[1] + "-01"),
                "old_category": m[2],
                "current_month": date.fromisoformat(m[3] + "-01"),
                "new_category": m[4],
                "mechanism": m[5],
            }
        )
    return rows
