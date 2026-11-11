"""Audit count recovery, matched edges and the reversible August 2022 batch."""

import argparse
from datetime import date
from pathlib import Path

from clinvar_history.checkpoints import atomic_json
from clinvar_history.research import _copy, research_connection
from clinvar_history.research_sampling import audit_sampling


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--root", type=Path, default=Path("data/research"))
    p.add_argument(
        "--output", type=Path, default=Path("results/full-research/sampling")
    )
    args = p.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)
    results = []
    for policy in ["strict", "family"]:
        for excluded in [(), (date(2022, 8, 12),)]:
            for name, months, days in [
                ("native", 0, 0),
                ("weekly", 0, 7),
                ("monthly", 1, 0),
                ("quarterly", 3, 0),
            ]:
                key = f"{policy}_{name}" + ("_omit20220812" if excluded else "")
                print(key, flush=True)
                results.append(
                    audit_sampling(
                        args.root,
                        args.output / (key + ".parquet"),
                        policy=policy,
                        months=months,
                        days=days,
                        excluded_dates=excluded,
                    )
                )
                atomic_json(results, args.output / "summary.json")
    with research_connection() as c:
        source = str(args.root / "events/vcv/release_date=*/events.parquet")
        _copy(
            c,
            """WITH changes AS (
            SELECT *,lead(observed_at) OVER(PARTITION BY variation_id ORDER BY observed_at) AS next_at,
              lead(old_raw_classification) OVER(PARTITION BY variation_id ORDER BY observed_at) AS next_old,
              lead(new_raw_classification) OVER(PARTITION BY variation_id ORDER BY observed_at) AS next_new
            FROM read_parquet(?) WHERE event_type='classification_changed'
        ) SELECT observed_at,old_raw_classification,new_raw_classification,count(*) AS transitions,
            count(*) FILTER(WHERE next_old=new_raw_classification AND next_new=old_raw_classification AND date_diff('day',observed_at,next_at)<=21) AS exact_reversal_within_21_days
        FROM changes WHERE old_classification!=new_classification GROUP BY ALL ORDER BY transitions DESC""",
            args.output / "reversals.parquet",
            [source],
        )


if __name__ == "__main__":
    main()
