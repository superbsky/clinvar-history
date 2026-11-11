# Full-population ClinVar observation and first-exit protocol

This protocol describes the November 2026 computational extension. It replaces the earlier 5,000-variant pilot protocol. The [full report](full-research-results.md) gives results, the [literature matrix](literature-evidence.md) records the novelty assessment. The protocol was developed after inspection of the retained dataset; it is not preregistered. The reversible-batch sensitivity and targeted intervening-version selection are explicitly post hoc.

## Research questions

1. Does an independently written reconstruction reproduce every retained incident episode and attribution result? Are configured dates actually observed?
2. How do canonical labels, clinical family bins, net endpoint counts and recovered reference edges differ?
3. Is cadence loss spread through time or concentrated in reversible series? Does omitting the relevant observation window change the interpretation?
4. What are first observed VUS-exit cumulative incidences in complete incident cohorts, with stated follow-up and censoring assumptions? How much does changed cohort eligibility contribute to cadence comparisons?
5. Can richer current, dated and explicitly versioned XML provide inspectable evidence absent from monthly SCV summaries, and what remains uncertain?

Computations report full-population denominators, independent numerical checks, documented samples and source provenance. The [literature matrix](literature-evidence.md) distinguishes existing methods from the questions addressed here.

## Data, units and endpoints

The frozen baseline is 417 public sources. The extended research dataset appends the newly discoverable 20 July 2026 VCF, yielding 339 VCF and 79 SCV sources. Administrative follow-up ends at the actual last source of each layer: 20 July for VCV, nominal 1 July for SCV. Monthly SCV dates are source-month coordinates, not exact submission observation dates. Coverage is the precise-coordinate GRCh38 VCF subset and retained SCV summaries, not all ClinVar records or patient testing.

Independent SQL reconstruction reads structural state events, compresses adjacent identical states and derives episode boundaries with lead/lag. All noninitial VUS entries are eligible. Select each unit's first observed incident episode, rather than its first episode in a filtered year. Initial VUS inventory is excluded; a later recurrence of a baseline-prevalent unit may still enter. VCV uses VariationID; SCV uses accession without version. Flag accessions associated with multiple VariationIDs across history and exclude these in the primary SCV cohort, retaining inclusion as a sensitivity.

The strict endpoint policy preserves the original canonical mapping. The family policy combines B/LB and P/LP, including exact corresponding composite labels. It leaves other mixed/risk/unknown labels uncollapsed. Raw labels and original ingestion events remain unchanged. First exits are benign-family, pathogenic-family, conflict or another non-VUS state. Conflict is absorbing only for the first-exit estimand. Removal is censoring in the primary model and another competing exit in a sensitivity. These are recorded database outcomes, not patient-level risk or biological truth.

## Cadence and cohort alignment

Reconstruct native, regular seven-day, monthly and quarterly observations from the same state history; anchor regular grids on the first VCF date and include the last source date. Coarse state is the last state at the observation boundary. Presence/absence is retained.

For transition audits, exclude additions/removals from classification-edge counts. Report net edge count divided by the number of native reference edges as **count recovery**. A **matched edge** must have the same unit, old state and new state within the sample interval. Credit at most one reference edge per sampled edge, using its latest matching date for lag. Report unmatched synthesized edges separately. Thus B → VUS → P sampled as B → P contributes one count but matches neither edge. This matching does not recover an exact decision time or distinguish every latent occurrence of a repeated identical edge.

For first-episode analyses, reidentify eligibility at each cadence. Then additionally align common first episodes: coarse entry must fall within the first native VUS episode, with administrative equality permitted for censored episodes. Analyze the same aligned units at every cadence, and report the loss from the complete native cohort. These descriptive comparisons separate selection from within-aligned-cohort observation changes without claiming a randomized or causal cadence effect.

A post hoc sensitivity omits native observations on 12 and 16 August 2022, preserving their updates until the next retained observation. Omitting only 12 August is retained as a negative control. Neither analysis changes the persisted event log or certifies an erroneous source release. Quantify exact next-edge raw reversals within 21 days separately from this omission analysis.

## Estimation and uncertainty

Use tied-time Aalen–Johansen cumulative incidence, decrementing survival by all observed competing exits. Count events before censoring at tied times. Report 365-, 730-, 1,095- and 1,825-day estimates, numbers at risk and supported horizons. Estimate complete first-incident cohorts, without variant subsampling.

Use 1,000 seed-fixed nonparametric bootstrap replicates, resampling whole VariationIDs. Multiple SCVs from a variant remain together. Compress identical cluster vectors of duration/outcome codes and resample their frequencies multinomially; grouping exchangeable signatures is exactly equivalent to ordinary cluster resampling. Report pointwise percentile 95% intervals only where at least 95% of replicates support the horizon. Retain maximum half-run versus full-run CI endpoint changes as a Monte Carlo diagnostic. These are not simultaneous bands and do not account for additional shared submitter/gene dependence or structural data uncertainty.

Verify point curves against the independently implemented statsmodels `CumIncidenceRight` on the complete cohorts. Tests cover tied censoring, a cluster-bootstrap control and independently reconstructed paths. Agreement verifies computation, not independent censoring. Calendar-dependent entry and censoring may violate that assumption; report annual entry strata and cohorts with at least one, three or five years of possible administrative follow-up. Unsupported recent-year horizons remain missing; partial follow-up within a supported year is still described explicitly.

Timing sensitivities use adjacent-observation midpoints and extreme feasible entry/exit placements. Censoring and the unit's duration bounds remain explicit. These are stress assumptions, not fitted interval-censored likelihoods or identified upper/lower incidence bounds. Removal-as-exit and inclusion of migrated identities are separate sensitivities. All-episode analysis clusters by variant but estimates an episode-weighted quantity; it is supplementary to first-exit estimates.

Exploratory strata include entry origin, source review status, variant type, gene/condition strings, collection method and submitter where at least 1,000 variant clusters are available; retain the largest six categories, or ten gene/submitter categories. Do not explode multigene/condition strings into independent observations. Do not perform post hoc significance ranking or imply causal differences. Inspect risk sets and Monte Carlo diagnostics before interpreting long-horizon small strata.

## Attribution and XML evidence

Independently reproduce conservative canonical attribution with the same 35-day limit and strict overlap of observation intervals. Compare every baseline label, candidate set and supporting-event identifier. Add source-end observability to the extended export: full, partial or absent SCV interval coverage. Do not label a source-follow-up limitation as a demonstrated absence of a mechanism.

Current XML: select 100 unexplained events by seeded balanced allocation over year/canonical-direction strata; retain stratum population counts and status of every request. Record versions, timestamps, URLs, SHA-256 and explicit byte budgets. Measure availability of explanations, methods, comments, citations and contribution flags. Current records alone cannot resolve historical causes.

Historical XML: select 100 unexplained 12 August events, balanced across raw-direction strata, with seeded SHA-256 order and population/sample counts. Stream the currently available August/September 2022 old-format VCV archives; retain only selected records, official compressed MD5, sizes, HTTP Last-Modified, root date and extraction counts. Report missing/oversize cases. These archives predate/postdate the transient window but were modified in 2026; do not call them immutable original 2022 bytes. Compare aggregate endpoints and SCV identity/version/label/metadata signatures. Unweighted balanced-sample fractions are not population mechanism rates.

For ten largest direction strata with a two-version increase, retrieve the first selected case's before/intermediate/after VCV versions, maximum 30 distinct records and 2 MB each. Record failed or fallback requests explicitly; require exact returned versions. This is targeted case investigation. Compare review metadata and labels separately. Version-specific XML may include later metadata changes. Retain the previous 978448 case as a separate feasibility example; do not use it to infer prevalence. No automated finding is blinded curator adjudication or proof of cause.

## Reproduce

The compact datasets are generated locally from pinned public sources. Pinned public-source manifests are committed as CSV. If starting from an empty checkout, restore the baseline manifest before the managed backfill:

```bash
uv sync --locked --all-groups
uv run python scripts/restore_research_manifest.py
scripts/run_full_backfill.sh
```

This can transfer approximately 41 GiB of compressed input while streaming it. Exact reconstruction depends on the recorded sources remaining available with their recorded checksums. A checksum mismatch must be investigated, not silently accepted as the original dataset. The November run reused and independently audited the retained baseline rather than performing this complete network rebuild.

With the frozen baseline available, the actual extension sequence is:

```bash
uv run python scripts/prepare_full_research.py
uv run --group research python analysis/10_full_research_audit.py \
  --root data/research --baseline-root data/full
uv run --group research python analysis/11_full_cohort_incidence.py
uv run --group research python analysis/12_sampling_edge_audit.py
uv run --group research python analysis/14_cohort_alignment.py
uv run --group research python analysis/15_reversible_batch_sensitivity.py
uv run --group research python analysis/16_attribution_coverage_audit.py
uv run clinvar-history review-xml --root data/full \
  --attributions results/full-v1/attributions.parquet \
  --output results/full-research/xml-current --limit 100 \
  --max-total-bytes 200000000
uv run python scripts/fetch_research_archives.py
uv run --group research python analysis/13_historical_xml_audit.py
uv run --group research python analysis/17_xml_evidence_summary.py --fetch-versions
MPLCONFIGDIR=/tmp/clinvar-matplotlib uv run --group research \
  python analysis/18_build_full_research_report.py
```

Archives require about 4.03 GB compressed; full decompressed archives are never written. Allow additional space for compact data, research results and DuckDB spill. Native ingestion/results publish atomically; the incidence jobs validate per-input and implementation identities on reuse. Archive transfers support verified resumption. Other research scripts are repeatable batch analyses; they do not all provide per-stage automatic checkpoint reuse. `17_xml_evidence_summary.py --fetch-versions --offline` verifies the retained version cache without another request. Current review also supports `--offline`.

For validation, run `uv run --group research pytest -q`, `uv run ruff check .`, `uv run ruff format --check .` and `uv build`. Use `TZ=UTC` for a consistent time-zone environment.

## Interpretation limits

Calendar-dependent censoring, observation-time assumptions and the pooled estimand limit interpretation. Versioned review-status round trips and discordant archival cases require curation expertise. The analysis does not include a fitted interval-censored recurrent multistate model, complete historical XML ingestion or blinded representative curator adjudication.
