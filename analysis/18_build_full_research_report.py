"""Build auditable full-research tables, figures and a results narrative."""

import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pyarrow.parquet as pq
from pyarrow import csv
from statsmodels.duration.survfunc import CumIncidenceRight

from clinvar_history.checkpoints import atomic_json, dataset_snapshot, file_digest
from clinvar_history.classification import ClassificationState, normalize
from clinvar_history.manifest import read_manifest
from clinvar_history.research import research_connection
from clinvar_history.research_statistics import curves_from_counts, load_cluster_counts

ROOT = Path("results/full-research")
DOCS = Path("docs")
FIG = DOCS / "figures/full-research"
TABLES = DOCS / "research-tables"
COLORS = {"benign": "#2171b5", "pathogenic": "#c44e52", "conflict": "#b07a12"}


def read(path):
    return pq.read_table(path).to_pylist()


def key_label(text):
    n = normalize(text)
    return (
        n.state.value
        if n.state is not ClassificationState.OTHER
        else text.casefold().strip()
    )


def export(table, name):
    # CSV is an interoperable research export, not a formatted spreadsheet.
    csv.write_csv(table, TABLES / (name + ".csv"))


def save(fig, name):
    svg = FIG / (name + ".svg")
    fig.savefig(svg, bbox_inches="tight", metadata={"Date": None})
    svg.write_text(
        "\n".join(line.rstrip() for line in svg.read_text().splitlines()) + "\n"
    )
    fig.savefig(FIG / (name + ".png"), bbox_inches="tight", dpi=160)
    plt.close(fig)


def md_table(columns, rows):
    return "\n".join(
        [
            "| " + " | ".join(columns) + " |",
            "| " + " | ".join(["---"] * len(columns)) + " |",
            *["| " + " | ".join(str(v) for v in row) + " |" for row in rows],
        ]
    )


def pct(x):
    return "—" if x is None else f"{100 * x:.3f}"


def interval(row, cause):
    return (
        f"{pct(row[cause])} ({pct(row[cause + '_lower'])}–{pct(row[cause + '_upper'])})"
    )


def main():
    FIG.mkdir(parents=True, exist_ok=True)
    TABLES.mkdir(parents=True, exist_ok=True)
    plt.rcParams.update(
        {
            "font.size": 10,
            "axes.spines.top": False,
            "axes.spines.right": False,
            "svg.fonttype": "none",
            "svg.hashsalt": "clinvar-history-full-research-20261111",
            "axes.titleweight": "bold",
        }
    )
    source_root = Path("data/research")
    snapshot = dataset_snapshot(
        source_root, read_manifest(source_root / "manifest/releases.parquet")
    )
    assert snapshot == json.loads((ROOT / "audit.json").read_text())["input"], (
        "Research inputs changed since reconstruction"
    )
    assert (
        snapshot["manifest"]
        == "95bc82bf4724995d7e7191c9b5bbd24258d60e1e626cced5392bc33e8ddfdd4d"
    ), (
        "Narrative is pinned to the November source manifest; revise it for a different study"
    )
    incidence = pq.read_table(ROOT / "full_cohort_incidence.parquet")
    rows = incidence.to_pylist()
    export(incidence, "full-cohort-incidence")
    sampling = json.loads((ROOT / "sampling/summary.json").read_text())
    window = json.loads((ROOT / "reversible-batch/sampling.json").read_text())
    import pyarrow as pa

    scalar_sampling = [
        {k: v for k, v in r.items() if k != "excluded_dates"}
        | {"excluded_dates": ";".join(r["excluded_dates"])}
        for r in sampling + window
    ]
    export(pa.Table.from_pylist(scalar_sampling), "sampling-edge-audit")
    export(pq.read_table(ROOT / "mechanism_coverage.parquet"), "mechanism-coverage")
    export(
        pq.read_table(ROOT / "xml-historical/comparisons.parquet").drop(
            ["before_labels", "after_labels"]
        ),
        "historical-xml-cases",
    )
    checks = []
    for layer in ["vcv", "scv"]:
        data = load_cluster_counts(ROOT / f"episodes_{layer}_family_0.parquet")
        weights = np.asarray(data.matrix @ data.frequencies).ravel()
        other = CumIncidenceRight(
            data.codes // 5 / 2, data.codes % 5, freq_weights=weights
        )
        ours, *_ = curves_from_counts(data.codes, weights)
        err = float(
            np.max(
                np.abs(
                    ours[np.asarray(other.times * 2, dtype=int), :, 0]
                    - np.asarray(other.cinc).T
                )
            )
        )
        assert err < 1e-10
        checks.append(
            {
                "layer": layer,
                "episodes": int(weights.sum()),
                "clusters": data.clusters,
                "independent_implementation": "statsmodels.duration.survfunc.CumIncidenceRight",
                "max_absolute_difference": err,
            }
        )
        export(
            pq.read_table(ROOT / f"alignment/{layer}_incidence.parquet"),
            f"{layer}-aligned-incidence",
        )
        export(
            pq.read_table(ROOT / f"alignment/{layer}_recurrence.parquet"),
            f"{layer}-recurrence",
        )
    atomic_json(checks, ROOT / "independent_estimator_check.json")
    series = json.loads((ROOT / "xml-versions/series.json").read_text())["cases"]
    series_summary = []
    for case in series:
        row = {
            "variation_id": case["event"]["variation_id"],
            "status": case["status"],
            "stratum": case["event"]["stratum"],
        }
        if case["status"] == "retrieved":
            records = [case[k] for k in ["before", "middle", "after"]]
            labels = [r["current_germline_classification"] for r in records]
            row.update(
                {
                    "before_label": labels[0],
                    "middle_label": labels[1],
                    "after_label": labels[2],
                    "matches_vcf_round_trip": key_label(labels[0])
                    == key_label(case["event"]["old_raw_classification"])
                    and key_label(labels[1])
                    == key_label(case["event"]["new_raw_classification"])
                    and key_label(labels[2]) == key_label(labels[0]),
                    "assertion_identity_labels_versions_stable": len(
                        {
                            json.dumps(
                                sorted(
                                    (a["accession"], a["version"], a["classification"])
                                    for a in r["assertions"]
                                )
                            )
                            for r in records
                        }
                    )
                    == 1,
                    "changed_fields_first": ";".join(
                        sorted(
                            {
                                field
                                for change in case["first_scv_changes"]
                                for field in change["changed_fields"]
                            }
                        )
                    ),
                    "changed_fields_second": ";".join(
                        sorted(
                            {
                                field
                                for change in case["second_scv_changes"]
                                for field in change["changed_fields"]
                            }
                        )
                    ),
                    "review_status_round_trip": all(
                        change["change"] == "updated"
                        and change["changed_fields"] == ["review_status"]
                        and change["before"]["review_status"]
                        == "criteria provided, single submitter"
                        and change["after"]["review_status"]
                        == "no assertion criteria provided"
                        for change in case["first_scv_changes"]
                    )
                    and bool(case["first_scv_changes"])
                    and all(
                        change["change"] == "updated"
                        and change["changed_fields"] == ["review_status"]
                        and change["before"]["review_status"]
                        == "no assertion criteria provided"
                        and change["after"]["review_status"]
                        == "criteria provided, single submitter"
                        for change in case["second_scv_changes"]
                    )
                    and bool(case["second_scv_changes"]),
                }
            )
        series_summary.append(row)
    export(pa.Table.from_pylist(series_summary), "targeted-version-series")
    atomic_json(series_summary, ROOT / "xml-versions/summary.json")
    # Figure 1: every canonical release change, no smoothing.
    with research_connection() as c:
        changes = c.execute(
            "SELECT observed_at,sum(events) FROM read_parquet(?) WHERE old_classification!=new_classification GROUP BY observed_at ORDER BY observed_at",
            [str(ROOT / "raw_transition_audit_vcv.parquet")],
        ).fetchall()
    fig, axes = plt.subplots(
        1, 2, figsize=(12, 3.7), gridspec_kw={"width_ratios": [2, 1]}
    )
    axes[0].plot(
        [r[0] for r in changes], [r[1] for r in changes], lw=1, color="#38566b"
    )
    axes[0].set(
        title="Canonical VCV changes by release",
        ylabel="Transitions (log scale)",
        yscale="log",
        xlabel="Observation date",
    )
    selected = [r for r in changes if str(r[0]).startswith("2022-08")]
    axes[1].bar(
        [str(r[0])[5:] for r in selected],
        [r[1] for r in selected],
        color=["#b07a12", "#b07a12", "#38566b", "#38566b"],
    )
    axes[1].set(
        title="August 2022", ylabel="Transitions", xlabel="Release date (MM-DD)"
    )
    axes[1].ticklabel_format(axis="y", style="plain")
    fig.tight_layout()
    save(fig, "01-release-transitions")
    # Figure 2: exact matched reference edges, plus clearly labeled window omission.
    fig, axes = plt.subplots(1, 2, figsize=(11, 3.8), sharey=True)
    for ax, policy, title in zip(
        axes,
        ["strict", "family"],
        ["Strict canonical labels", "Clinical family bins"],
        strict=True,
    ):
        for dataset, label, color, offset in [
            (sampling, "All observations", "#38566b", -0.18),
            (window, "Omit Aug 12 & 16", "#b07a12", 0.18),
        ]:
            vals = [
                next(
                    r
                    for r in dataset
                    if r["endpoint_policy"] == policy
                    and r["months"] == m
                    and r["days"] == d
                    and (dataset is window or not r["excluded_dates"])
                )
                for m, d in [(0, 0), (0, 7), (1, 0), (3, 0)]
            ]
            ax.bar(
                np.arange(4) + offset,
                [100 * r["edge_recall"] for r in vals],
                width=0.35,
                label=label,
                color=color,
            )
        ax.set(
            title=title,
            xticks=np.arange(4),
            xticklabels=["Native", "Weekly", "Monthly", "Quarterly"],
            ylim=(0, 105),
            ylabel="Matched-edge recall (%)",
        )
        ax.legend(
            frameon=True,
            framealpha=1,
            facecolor="white",
            edgecolor="white",
            fontsize=9,
            loc="lower left",
        )
    fig.tight_layout()
    save(fig, "02-sampling-edge-recall")
    # Figure 3: first-exit estimands, separate units and risk sets.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4), sharex=True, sharey=True)
    for ax, layer, title in zip(
        axes,
        ["vcv", "scv"],
        [
            "VCV: first incident episode per variant",
            "SCV: first incident episode per accession",
        ],
        strict=True,
    ):
        curve = read(ROOT / f"incidence/{layer}_primary_curve.parquet")
        for cause, color in COLORS.items():
            ax.step(
                [r["day"] / 365 for r in curve],
                [100 * r[cause] for r in curve],
                where="post",
                color=color,
                label=cause.capitalize(),
            )
        ax.set(
            title=title,
            xlabel="Years since first observed incident VUS",
            ylabel="First-exit cumulative incidence (%)",
            xlim=(0, 5),
            ylim=(0, 8),
        )
        ax.legend(frameon=False)
    fig.tight_layout()
    save(fig, "03-first-exit-curves")
    # Figure 4: entry counts and calendar-cohort heterogeneity.
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for layer, color, offset in [("vcv", "#38566b", -0.18), ("scv", "#b07a12", 0.18)]:
        years = [
            next(
                r
                for r in rows
                if r["layer"] == layer
                and r["analysis"] == f"entry_{year}"
                and r["day"] == 365
            )
            for year in range(2020, 2027)
        ]
        axes[0].bar(
            np.arange(7) + offset,
            [r["episodes"] / 1000 for r in years],
            width=0.35,
            color=color,
            label=layer.upper(),
        )
    axes[0].set(
        title="First incident VUS entry cohorts",
        xticks=np.arange(7),
        xticklabels=range(2020, 2027),
        ylabel="Units (thousands)",
        xlabel="Entry year; 2026 is partial",
    )
    axes[0].legend(frameon=False)
    for layer, style in [("vcv", "-o"), ("scv", "--s")]:
        values = [
            next(
                r
                for r in rows
                if r["layer"] == layer
                and r["analysis"] == f"entry_{year}"
                and r["day"] == 1095
            )
            for year in range(2020, 2024)
        ]
        for cause in ["benign", "pathogenic"]:
            axes[1].plot(
                range(2020, 2024),
                [100 * r[cause] for r in values],
                style,
                color=COLORS[cause],
                label=f"{layer.upper()} {cause}",
            )
    axes[1].set(
        title="Three-year estimates by entry year",
        xticks=range(2020, 2024),
        xlabel="Entry year; 2023 has partial 3y follow-up",
        ylabel="Cumulative incidence (%)",
        ylim=(0, 6),
    )
    axes[1].legend(frameon=False, fontsize=8)
    fig.tight_layout()
    save(fig, "04-entry-cohorts")
    # Figure 5: separate eligibility from cadence among matched first episodes.
    aligned = read(ROOT / "alignment/vcv_incidence.parquet")
    without = next(
        r for r in read(ROOT / "reversible-batch/incidence.parquet") if r["day"] == 1095
    )
    compare = [
        next(
            r
            for r in rows
            if r["layer"] == "vcv" and r["analysis"] == "primary" and r["day"] == 1095
        ),
        *[
            next(r for r in aligned if r["analysis"] == a and r["day"] == 1095)
            for a in ["aligned_native", "aligned_monthly", "aligned_quarterly"]
        ],
        without,
    ]
    fig, axes = plt.subplots(1, 3, figsize=(12, 4.4))
    labels = [
        "All native",
        "Aligned native",
        "Aligned monthly",
        "Aligned quarterly",
        "Window omitted",
    ]
    for ax, cause in zip(axes, COLORS, strict=True):
        values = np.array([100 * r[cause] for r in compare])
        low = np.array([100 * r[cause + "_lower"] for r in compare])
        high = np.array([100 * r[cause + "_upper"] for r in compare])
        ax.barh(
            np.arange(5),
            values,
            color=COLORS[cause],
            alpha=0.85,
            xerr=[values - low, high - values],
            capsize=2,
        )
        ax.set(
            title=cause.capitalize(),
            yticks=np.arange(5),
            yticklabels=labels,
            xlabel="Three-year incidence (%)",
        )
        ax.invert_yaxis()
    fig.tight_layout()
    save(fig, "05-cohort-alignment")
    # Compact evidence bundle pins all report inputs without distributing raw archives.
    paths = sorted(
        {
            p
            for p in ROOT.rglob("*")
            if p.is_file()
            and p.suffix in {".json", ".parquet"}
            and ".scratch" not in p.parts
            and p.name != "report-input-provenance.json"
        }
    )
    atomic_json(
        {
            "format": 1,
            "files": {str(p): file_digest(p) for p in paths},
            "sources": {
                str(p): file_digest(p)
                for p in [
                    *Path("analysis").glob("1*.py"),
                    *Path("src/clinvar_history").glob("*.py"),
                    *Path("scripts").glob("*research*.py"),
                    Path("pyproject.toml"),
                    Path("uv.lock"),
                ]
            },
            "dataset": snapshot,
            "estimator_checks": checks,
        },
        ROOT / "report-input-provenance.json",
    )
    primary = [r for r in rows if r["analysis"] == "primary"]
    incidence_table = md_table(
        [
            "Layer",
            "Days",
            "Units",
            "At risk",
            "Benign % (95% interval)",
            "Pathogenic % (95% interval)",
            "Conflict % (95% interval)",
        ],
        [
            [
                r["layer"].upper(),
                r["day"],
                f"{r['episodes']:,}",
                f"{r['at_risk']:,}",
                interval(r, "benign"),
                interval(r, "pathogenic"),
                interval(r, "conflict"),
            ]
            for r in primary
        ],
    )
    sampling_table = md_table(
        [
            "Policy",
            "Observation",
            "Reference edges",
            "Net edges",
            "Count recovery %",
            "Matched recall %",
            "Synthesized edges",
        ],
        [
            [
                r["endpoint_policy"],
                ("Window omitted; " if r in window else "")
                + (
                    "quarterly"
                    if r["months"] == 3
                    else "monthly"
                    if r["months"] == 1
                    else "weekly"
                    if r["days"]
                    else "native"
                ),
                f"{r['reference_edges']:,}",
                f"{r['sampled_net_edges']:,}",
                f"{100 * r['count_recovery']:.3f}",
                f"{100 * r['edge_recall']:.3f}",
                r["synthetic_edges"],
            ]
            for r in sampling + window
            if (not r["excluded_dates"] or r in window) and r["months"]
        ],
    )
    alignment_table = md_table(
        ["VCV analysis", "Units", "Benign %", "Pathogenic %", "Conflict %"],
        [
            [
                name,
                f"{r['episodes']:,}",
                pct(r["benign"]),
                pct(r["pathogenic"]),
                pct(r["conflict"]),
            ]
            for name, r in zip(labels, compare, strict=True)
        ],
    )
    hist = json.loads((ROOT / "xml-historical/summary.json").read_text())
    cur = json.loads((ROOT / "xml-current/summary.json").read_text())
    nseries = len(series_summary)
    supported = sum(
        r.get("review_status_round_trip", False)
        and r.get("assertion_identity_labels_versions_stable", False)
        for r in series_summary
    )
    text = f"""# Full ClinVar research extension: observation, endpoints and first VUS exits

Research run: 11 November 2026. The [July baseline report](full-backfill-results.md) provides the earlier source population and net-count analysis. The definitions below distinguish count recovery from matched-edge recall and use the extended source population.

The main finding is that the previously reported approximately 20% reduction in monthly transition counts is dominated by a short, highly reversible August 2022 batch. This changes the interpretation of the original result. It does not justify deleting those observations or declaring a database error. Full incident cohorts, cadence-aligned cohorts, dated XML and explicitly versioned XML provide complementary evidence.

## Analytical scope

| Analysis | Evidence | Limit |
|---|---|---|
| VUS duration and competing risks | All first incident VCV/SCV units; 92 prespecified executable analysis jobs; 1,000 variant-cluster bootstrap draws each; timing/removal/label/calendar sensitivities; aligned-cohort and recurrence analyses | Release-time descriptive estimates, with independent-censoring assumptions; no fitted biological interval-censoring likelihood |
| Unexplained transitions | 100 current records; 100 cases in each of two dated archives; {nseries} targeted three-version series; explicit population/selection/provenance | Automated evidence audit, not blinded curator validation or a population causal explanation rate |
| Safe updates | Managed semantic checkpoints and interruption/reuse tests; bounded real update; independently discovered July 20 release ingested into a separate research copy | The original backfill was audited, not downloaded again in full |
| Reproducibility | [Protocol](research-protocol.md), [literature matrix](literature-evidence.md), analysis tables, figures and source manifests | Observational research with explicit source, sampling and censoring limitations |

## Sources and independent reconstruction

The frozen baseline contains 417 releases: 338 VCF and 79 monthly submission summaries. Independent SQL reconstruction exactly matched all 2,276,393 VCV and 2,997,859 SCV exported episodes in that baseline, comparing unit, start, end, outcome and duration with multiset differences in both directions. This is an independent implementation check of the retained data, not an independent acquisition of every source.

Re-discovery within the same configured horizon found one additional available VCF release, 20 July 2026. The research copy therefore contains **418 sources: 339 VCF and 79 SCV summaries**. Its actual observation boundaries are 6 January 2020–20 July 2026 for VCV and nominal 1 January 2020–1 July 2026 for SCV. Administrative censoring uses each layer's actual last observation. The configured 24 July horizon is not treated as an observed release. The old exports had extended censoring nine VCV days and 23 nominal SCV days beyond their last source.

The research copy preserves the baseline and independently reconstructs 2,286,991 VCV and 2,997,859 SCV incident episodes. There are 2,280,091 distinct VCV first-incident units and 2,997,087 SCV accessions. Across all SCV events, 843 accessions have more than one VariationID. Excluding such identities removes 60 first-incident SCVs (63 episodes); the main SCV analysis has 2,997,027 accessions within 2,454,492 variant clusters. Including them is a sensitivity analysis. The VCF subset requires precise GRCh38 coordinates and does not cover the whole ClinVar archive.

Time zero is the first *observed incident* VUS episode, after excluding VUS inventory in the initial snapshot. A later return to VUS in a baseline-prevalent unit can enter this cohort. The estimand is therefore not first lifetime uncertainty. Entry metadata are the fields present on the structural state event, not continuously updated time-varying covariates.

## Endpoint policy changes what counts as a transition

Original normalization retains exact raw labels, maps the five ACMG classes individually and leaves unrecognized or mixed labels as `other`. The research family policy combines B/LB and P/LP, including the exact within-family composite labels, while retaining risk terms and mixed/unknown labels separately. It does not modify ingestion or assign an ordinal level to a composite label.

The extended history contains **481,556 canonical changes**. Among raw classification changes, 166,290 canonical changes stay in the same clinical family; conversely, 28 family changes were hidden inside canonical `other` to `other`. There are 315,294 family-level reference edges. Neither count alone measures a change in a laboratory's biological interpretation.

## Cadence: net counts versus matched recorded edges

A sampled B → P endpoint can summarize a recorded B → VUS → P path. It contributes one net edge but recovers neither original edge. We therefore report count recovery and matched-edge recall separately. A match requires the same variant, old state and new state within that sampling interval; at most one reference edge is credited per sampled edge. Unmatched sampled edges are synthesized summaries. Detection lag is computed only for matched edges. Absence is a separate state and additions/removals do not count as classification edges.

{sampling_table}

![Matched-edge recall](figures/full-research/02-sampling-edge-recall.png)

These are observation-grid reconstructions from the same underlying event history. They describe recovery of recorded states, not accuracy against unobserved clinical decisions. Weekly/monthly/quarterly grids are anchored on the first VCF date and include the final source date; SCV monthly coordinates already have monthly cadence.

## A short reversible batch dominates the original count loss

The 12 August 2022 VCF observation contains 43,584 canonical changes; 43,303 (99.36%) have an exactly reversed next raw-label transition within 21 days. A second observation on 16 August preserves most of the changed state until the 24 August release. Omitting 12 August alone is consequently a weak negative control: unchanged state updates move into the next observation. The substantive sensitivity omits *both* 12 and 16 August observation boundaries, folds their updates forward and compares the final state on 24 August with the last retained state. It does not simply delete event rows and leave an inconsistent history.

Strict reference edges fall from 481,556 to 394,880; monthly net edges remain 388,425. Thus 86,676 of the original 93,131 missing net edges (93.07%) are associated with this window. Monthly count recovery becomes 98.365%, with matched recall 98.176%. Family-level monthly matched recall becomes 97.634%. The sensitivity is post hoc and localizes the result; it is not a corrected ground truth or a recommendation to discard real transitions.

![Changes by release](figures/full-research/01-release-transitions.png)

## Full first-incident cumulative incidence

Aalen–Johansen estimation treats first benign, pathogenic, conflict and other exits as competing outcomes. Removal is censored at observed disappearance in the primary model. Events precede censored units at tied durations. Intervals below are pointwise 95% percentile intervals from 1,000 whole-VariationID bootstrap draws. Identical within-variant duration/outcome signatures are grouped; multinomial resampling of their frequencies is exactly equivalent to ordinary cluster resampling, retaining all SCVs within each sampled variant. It is not a sample of only 5,000 variants.

{incidence_table}

![First-exit curves](figures/full-research/03-first-exit-curves.png)

The independently implemented statsmodels cumulative-incidence estimator matched the complete VCV and SCV point curves to maximum absolute differences below 7e-17. Analytic/tied-time tests and a two-SCV-per-variant bootstrap control separately verify the implementation. These comparisons verify the estimator on the selected observations, not the independent-censoring assumption.

Pooled five-year estimates rely on early cohorts: the VCV risk set declines to 164,665, and SCV to 269,214. Restricting entry to at least 1,825 days before the layer's source end gives five-year benign/pathogenic estimates of 4.024%/1.221% for VCV and 5.599%/1.699% for SCV. The corresponding pooled values are 2.417%/0.660% and 3.761%/1.141%. Calendar composition and unequal follow-up matter; the pooled curve is not a forecast for a variant entering in 2026. Recent entry strata with unsupported horizons are left missing.

![Entry cohorts](figures/full-research/04-entry-cohorts.png)

## Eligibility, timing, recurrence and other sensitivities

We also align the *same first native episode* across all three cadences: a unit must be eligible at every cadence and its coarser entry must occur within the first native VUS episode. This retains 2,260,601 VCV and 2,996,730 SCV units. Simple unit intersection without episode overlap could pair unrelated recurrences, so it is not used.

{alignment_table}

![VCV alignment](figures/full-research/05-cohort-alignment.png)

At three years, full VCV conflict incidence is 3.555%, but it is 2.736% in the aligned native cohort. The 0.819 percentage-point difference already exists before changing timing within that matched cohort. Monthly aligned conflict incidence is 2.773%. Omitting the reversible window gives 2.775% in its own reconstructed cohort. These are descriptive comparisons with changed selection; no causal cadence effect is claimed. Quarterly aligned benign incidence falls from 1.591% to 1.348%, showing that the small effect in the earlier bounded pilot does not generalize to every full-cohort comparison.

At three years, midpoint timing gives benign/pathogenic estimates 1.578%/0.399% for VCV and 2.170%/0.607% for SCV, versus 1.581%/0.400% and 2.192%/0.616% at recorded boundaries. Extreme entry/exit interval placements are also exported; they are timing stress assumptions, **not identified cumulative-incidence bounds**. Treating removal as a competing exit and including migrated SCV identities changes these pooled estimates little. This does not exclude informative censoring.

There are 6,900 recurrent VCV episodes and 772 recurrent SCV episodes beyond the first observed incident episode. Supplementary all-episode estimates cluster by VariationID and weight episodes; they do not estimate the first-episode probability or eventual resolution after any conflict. The primary analysis makes conflict absorbing only for its first-exit estimand. A fitted recurrent multistate or interval-censored likelihood is not part of this analysis.

Exploratory entry-origin, review-status, variant-type, gene-string, condition-string, submitter and collection-method strata are exported when at least 1,000 variant clusters support selection. They describe raw source fields; combined genes or conditions are not exploded into independent observations. No multiplicity-adjusted association tests or causal rankings are made. Some five-year strata have small risk sets and visibly less stable bootstrap endpoints; numbers at risk and half-run versus full-run CI changes are retained for every result. Narrow pooled intervals do not remove structural bias or shared gene/submitter dependence.

## Attribution coverage and historical XML

A second SQL implementation reproduces all 479,525 baseline attribution labels, candidate sets and supporting-event identifiers with zero discrepancies. The extended data retain the same conservative rule: a candidate SCV event must share VariationID, be within 35 days and have strictly overlapping observation intervals. The new export exposes coverage separately: 2,110 transitions have no SCV interval overlap after the last nominal SCV source, and 187 have partial end-of-source coverage. These are observation limitations; they must not automatically be interpreted as intrinsically inexplicable changes.

The current-XML review selected 100 unexplained events across year/direction strata and retrieved all 100 records using {cur["bytes"]:,} bytes. It found aggregate explanations in {cur["aggregate_explanation"]}, assertion methods in {cur["assertion_methods"]}, contribution flags in {cur["contribution_flags"]} and an explicitly noncontributing SCV in {cur["records_with_noncontributing_scv"]}. Current XML alone established no historical cause. This balanced sample is a field-availability audit, not a population explanation rate.

A separate defined sample of 100 unexplained 12 August events spans {hist["raw_direction_strata"]} raw-direction strata. All 100 appear in both dated XML archives. Aggregate labels agree before and after the window in {hist["aggregate_unchanged"]}; SCV identities, versions and labels agree in {hist["assertion_labels_versions_unchanged"]}; the pre-window aggregate exactly matches the pre-window VCF raw label in {hist["before_aggregate_matches_vcf_pre_window"]}. Five cases with changed SCV evidence and two with changed aggregate endpoints are retained explicitly. Stratum population and sample counts accompany every selected event. Sparse allocation in common strata precludes treating these unweighted fractions as population estimates.

The compressed archives total 4,033,553,961 bytes, verified against official MD5 files. Root dates are 1 August and 1 September 2022; the catalog describes later release dates. Both currently available files were last modified in March 2026. They were streamed without retaining full decompressed XML. This is stronger evidence than current XML alone, but not a verified immutable capture of the originally distributed 2022 bytes; neither archive directly observes the transient window.

We therefore retrieved the intervening VCV version for one selected case in each of the ten largest direction strata with a two-version increase. All {nseries} targeted three-version series were retrieved. In {supported}, aggregate labels reproduce the VCF round trip while SCV accession/version/classification stay stable and SCV review status moves from criteria-provided to no-criteria and back. These are targeted examples, not a random validation sample. Version retrieval can include later metadata edits. The evidence is compatible with a temporary review-metadata/aggregation change; it does not identify who changed metadata, why, or whether it was erroneous. Existing candidate attribution labels were not silently rewritten.

## Literature, contribution and interpretation

The [focused evidence matrix](literature-evidence.md) documents searches and the closest overlapping paper and code. General longitudinal ClinVar analysis, survival estimation and heuristic mechanism assessment are already present in the literature. The defensible result here is the separation of endpoint counts, recorded-edge recovery, reversible batch behavior and incident-cohort selection, together with individually tracked SCV observations and reproducible evidence limits. A longer time series alone is not a novelty claim.

The main practical consequence is to inspect concentrated reversible series before generalizing a cadence-loss percentage. Conversely, rare transient records may matter to a researcher studying conflict even when they hardly affect benign/pathogenic incidence. A release history preserves those observations; a coarse endpoint comparison answers a different question. Neither provides biological ground truth or validates clinical decisions.

## Reproduction and remaining review

Use the commands in the [protocol](research-protocol.md). Machine-readable tables are in [research-tables](research-tables/); figures above are generated by `analysis/18_build_full_research_report.py`. Local inputs/results remain under ignored `data/research`, `artifacts/research-xml` and `results/full-research`. `report-input-provenance.json` fingerprints retained evidence; archive receipts and version caches retain retrieval provenance. No raw multi-gigabyte archive is committed to Git.

The estimates describe recorded database outcomes and have not been clinically validated. Calendar-dependent censoring, interval-time assumptions, targeted XML selection and discordant historical cases limit interpretation. Automated evidence checks do not replace independent statistical or ClinVar curation review.
"""
    (DOCS / "full-research-results.md").write_text(text)
    print(
        "Built full research report, 5 figure pairs and 9 research tables", flush=True
    )


if __name__ == "__main__":
    main()
