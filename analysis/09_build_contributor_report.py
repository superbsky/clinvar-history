"""Write a compact reproducible report from the contributor pilot artifacts."""

import argparse
import json
from pathlib import Path

import pyarrow.parquet as pq

from clinvar_history.checkpoints import file_digest


def percent(value):
    return "unsupported" if value is None else f"{100 * value:.2f}%"


def interval(row, cause):
    return f"{percent(row[cause])} ({percent(row[cause + '_lower'])}–{percent(row[cause + '_upper'])})"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--pilot-root", type=Path, default=Path("results/contributor-pilot")
    )
    parser.add_argument(
        "--output", type=Path, default=Path("results/contributor-pilot/report.md")
    )
    args = parser.parse_args()
    review = json.loads((args.pilot_root / "xml/review.json").read_text())
    fetched_dates = sorted(
        {
            e["review"]["provenance"]["fetched_at"][:10]
            for e in review["events"]
            if "provenance" in e["review"]
        }
    )
    fetched_label = ", ".join(fetched_dates) if fetched_dates else "unavailable"
    lines = [
        "# Contributor pilot results",
        "",
        "This exploratory analysis uses",
        f"the retained study ending 24 July 2026 and current XML fetched on {fetched_label}.",
        "It does not update the historical study horizon. The",
        "research protocol in docs/research-protocol.md defines the estimands and limitations.",
        "",
        "## First incident VUS episodes",
        "",
        "The VCV cohort includes first observed incident episodes starting 7 January–31",
        "December 2020; the SCV cohort starts 1 February–31 December 2020. Bootstrap",
        "inventory episodes were excluded by the upstream episode exporter. The first",
        "incident episode is selected before filtering cohort dates. Each layer samples",
        "5,000 distinct VariationIDs by SHA-256 ordering with seed 20261111. All eligible",
        "first SCV episodes within sampled variants are retained.",
        "",
        "The Aalen–Johansen estimator treats benign, pathogenic, conflict and other",
        "first exits as competing outcomes. Administrative end and removal are censored",
        "in the primary analysis. The following are cumulative incidence estimates",
        "before a competing first exit, not eventual resolution through later states.",
        "Parentheses show pointwise 95% percentile intervals from 200 variant-cluster",
        "bootstrap replicates. These are descriptive pilot intervals, not clinical risks.",
        "",
        "| Layer | Days | Episodes | At risk immediately before horizon | Benign or likely benign | Pathogenic or likely pathogenic |",
        "|---|---:|---:|---:|---|---|",
    ]
    tables = {}
    for layer in ["vcv", "scv"]:
        path = args.pilot_root / f"{layer}_competing_risks.parquet"
        metadata = json.loads(path.with_suffix(".parquet.checkpoint.json").read_text())
        expected = {
            "cohort_start": "2020-01-07" if layer == "vcv" else "2020-02-01",
            "cohort_end": "2020-12-31",
            "study_end": "2026-07-24",
            "max_clusters": 5000,
            "replicates": 200,
            "seed": 20261111,
        }
        if metadata["sha256"] != file_digest(path) or any(
            metadata["identity"]["parameters"][k] != v for k, v in expected.items()
        ):
            raise ValueError(
                "Report requires verified artifacts with the documented pilot parameters"
            )
        tables[layer] = pq.read_table(path).to_pylist()
        for row in tables[layer]:
            if (
                row["timing"] == "observed"
                and not row["removal_as_exit"]
                and row["stratum"] == "all"
            ):
                lines.append(
                    f"| {layer.upper()} | {row['day']} | {row['episodes']:,} | {row['at_risk']:,} | {interval(row, 'benign')} | {interval(row, 'pathogenic')} |"
                )
    lines += [
        "",
        "The VCV and SCV cohorts have different observation units and entry windows;",
        "their estimates should not be interpreted as a paired aggregate-versus-laboratory",
        "effect. SCV uncertainty resamples whole VariationIDs. Review-status strata",
        "with at least 30 sampled variant clusters are included in the Parquet outputs.",
        "",
        "## Sensitivity checks",
        "",
        "VCV monthly and quarterly results reconstruct actual endpoint snapshots for",
        "the same 5,000 source-cohort variants. Coarser observation can lose an episode",
        "or move entry outside the cohort. The resulting populations differ, so these",
        "are descriptive sensitivity estimates rather than a causal cadence comparison.",
        "",
        "| Layer | Timing or removal rule | Eligible episodes | Benign at 1,095 days | Pathogenic at 1,095 days |",
        "|---|---|---:|---:|---:|",
    ]
    for layer, rows in tables.items():
        for row in rows:
            if row["day"] == 1095 and row["stratum"] == "all":
                label = (
                    "removal as competing exit"
                    if row["removal_as_exit"]
                    else row["timing"]
                )
                lines.append(
                    f"| {layer.upper()} | {label} | {row['episodes']:,} | {percent(row['benign'])} | {percent(row['pathogenic'])} |"
                )
    lines += [
        "",
        "Midpoint sensitivity substitutes the midpoint of each adjacent-release interval",
        "for incident entry and observed exits. It is an assumption about observation",
        "coordinates, not an interval-censored likelihood model. Removal sensitivity",
        "treats removal as an additional absorbing other exit. Neither check establishes",
        "independent censoring. No estimate is extrapolated beyond observed support",
        "unless all subjects have already experienced an absorbing outcome.",
        "",
        "## Targeted XML evidence",
        "",
    ]
    review = json.loads((args.pilot_root / "xml/review.json").read_text())
    events = review["events"]
    successful = [
        e for e in events if e["review"]["status"] == "reviewed_current_record"
    ]
    lines += [
        f"All {len(successful)} of {len(events)} sampled unexplained transitions had retrievable current XML.",
        f"The sample used {review['total_xml_bytes']:,} XML bytes ({review['total_xml_bytes'] / 1024:.1f} KiB),",
        "with limits of 12 events, 2 MB per record and 20 MB total. Selection rotates",
        "deterministically across year and transition-direction strata. This allocation",
        "is deliberately balanced and is not a population-proportional random sample.",
        "",
        "| Extra XML evidence absent from compact event fields | Sampled events with field |",
        "|---|---:|",
    ]
    for field in [
        "aggregate_explanation",
        "submission_comments",
        "assertion_methods",
        "submission_citations",
        "contribution_flags",
    ]:
        lines.append(
            f"| {field.replace('_', ' ')} | {sum(e['review']['extra_fields'][field] for e in successful)}/{len(events)} |"
        )
    excluded = sum(
        any(a["contributes"] == "false" for a in e["review"]["assertions"])
        for e in successful
    )
    lines += [
        "",
        f"{excluded}/{len(events)} records contain at least one SCV explicitly excluded from the current",
        "aggregate classification. These flags motivate historical follow-up; they do",
        "not establish why an earlier aggregate transition occurred. Every sampled",
        "historical mechanism remains unresolved because only current XML was fetched.",
        "",
        "| VariationID | Historical date | Recorded transition | Current VCV version | Local historical SCVs / current XML SCVs |",
        "|---|---|---|---:|---:|",
    ]
    for event in events:
        r = event["review"]
        vid = event["variation_id"]
        lines.append(
            f"| [{vid}](https://www.ncbi.nlm.nih.gov/clinvar/variation/{vid}/) | {event['observed_at']} | {event['old_classification']} → {event['new_classification']} | {r.get('version', 'unavailable')} | {r.get('local_scv_accessions', 0)} / {len(r.get('assertions', []))} |"
        )
    lines += [
        "",
        "Local historical SCVs count distinct accessions observed anywhere in retained",
        "SCV history for the sampled variant; they are not submissions active at the",
        "transition date. Current XML may include later submissions or metadata edits.",
        "No historical attribution labels were changed by this review.",
    ]
    pair_path = args.pilot_root / "xml-version-pair/comparison.json"
    if pair_path.exists():
        pair = json.loads(pair_path.read_text())
        if (
            pair["variation_id"] != 978448
            or pair["old"]["version"] != "1"
            or pair["new"]["version"] != "2"
        ):
            raise ValueError(
                "Versioned case report expects VariationID 978448 versions 1 and 2"
            )
        lines += [
            "",
            "## Versioned XML follow-up",
            "",
            "A targeted follow-up selected VariationID 978448 because the current record",
            "had only two VCV versions and an excluded older SCV. This is a feasibility",
            "case chosen after the balanced sample review, not a population estimate.",
            "",
            "The [VCV version 1](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&id=VCV000978448.1)",
            "contains one likely-pathogenic SCV, SCV001433639. The",
            "[VCV version 2](https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=clinvar&rettype=vcv&id=VCV000978448.2)",
            "retains that SCV with the same label, marks it noncontributing, and adds",
            "SCV001983588 with a VUS label, supplied criteria and contribution enabled.",
            "The added SCV was submitted on 22 October 2021 and created in ClinVar on",
            "30 October, the VCV event date. Its first appearance in retained monthly SCV",
            "summaries is 1 December 2021, 32 days later. The VCV observation interval",
            "(25 October, 30 October] does not overlap the SCV interval (1 November,",
            "1 December].",
            "",
            "This is additional evidence compatible with submission composition and",
            "review precedence, and a concrete example of evidence missed by conservative",
            "interval-overlap attribution. It does not justify relaxing that rule globally.",
            "Version-specific XML can receive metadata edits without a VCV version bump;",
            "these fetched records are not immutable release-date snapshots. The historical",
            "attribution remains unchanged pending curator validation. The two XML records",
            f"used {sum(pair[k]['provenance']['bytes'] for k in ('old', 'new')):,} bytes in total.",
        ]
    lines += [
        "",
        "## Reproduce",
        "",
        "Run the commands in the README section “Contributor extensions”.",
        "Generated Parquets, cached XML and JSON provenance remain outside Git. Each",
        "XML response retains its retrieval URL, timestamp and SHA-256. Offline reruns",
        "verify cached XML hashes. Pilot Parquets record source hashes, parameters and",
        "implementation identity in adjacent checkpoint JSON files.",
        "",
        "| Local pilot input | SHA-256 |",
        "|---|---|",
    ]
    for path in [
        args.pilot_root / "vcv_competing_risks.parquet",
        args.pilot_root / "scv_competing_risks.parquet",
        args.pilot_root / "xml/review.json",
        *([pair_path] if pair_path.exists() else []),
    ]:
        lines.append(f"| {path.name} | `{file_digest(path)}` |")
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text("\n".join(lines) + "\n")
    print(f"Wrote {args.output}")


if __name__ == "__main__":
    main()
