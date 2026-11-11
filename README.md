# ClinVar classification history

[![CI](https://github.com/superbsky/clinvar-history/actions/workflows/ci.yml/badge.svg)](https://github.com/superbsky/clinvar-history/actions/workflows/ci.yml)
[![Python 3.12+](https://img.shields.io/badge/Python-3.12%2B-3776AB.svg)](https://github.com/superbsky/clinvar-history)
[![MIT license](https://img.shields.io/badge/license-MIT-blue.svg)](LICENSE)

> **Research use only:** This is a research pipeline, not a clinical
> interpretation tool.

Disk-bounded longitudinal processing of ClinVar aggregate variant (VCV) and
submission-level (SCV) classifications. The pipeline is designed to answer a
question that quarterly snapshots cannot: which intermediate and recurrent
changes occurred, and which aggregate changes have evidence of an underlying
SCV reclassification rather than a change in submission composition?

## Full research extension (November 2026)

The [full research report](docs/full-research-results.md) analyzes complete
incident cohorts, separates net endpoint counts from matched transition-edge
recall, and investigates a short reversible
August 2022 batch with dated and explicitly versioned XML. That window accounts
for 93.07% of the monthly net-count loss in the extended dataset; omitting its
observations is a sensitivity analysis, not a correction of clinical truth.

The [Liang et al. comparison](docs/liang-crosscheck.md) cross-checks
released case histories and the additional content of submonthly observations.

See the [protocol and reproduction commands](docs/research-protocol.md),
[literature assessment](docs/literature-evidence.md), tables and figures in the
report. Research dependencies are optional:
`uv sync --locked --all-groups`.

## Full-backfill results

The frozen July analysis, configured for 6 January 2020–24 July 2026,
produced these audited headlines:

- 417 releases total (338 VCF + 79 SCV), all completed;
- 479,525 canonical VCV transitions;
- aggregate VUS count increased from 217,177 to 2,315,022 (10.66x);
- against the release-level event log, net transition count recovery was 99.91% for the
  regular 7-day grid, 80.58% monthly, and 79.35% quarterly.

The [full English report](docs/full-backfill-results.md) documents the methods,
figures, audit evidence, and limitations.

## What it builds

- weekly GRCh38 VCF compact state and VCV event partitions;
- monthly `submission_summary` compact state and SCV event partitions;
- a DuckDB resume ledger;
- conservative candidate attribution of aggregate changes;
- quarterly baseline, sampling-bias, transition, dwell-time, competing-risk,
  and mechanism result tables.

Full ClinVar XML ingestion remains outside the production pipeline; a separate
bounded research workflow streams two dated XML archives and keeps selected records. Compressed source releases are
streamed, checked against an official `.md5` sidecar when one is available, and
discarded unless `--keep-raw` is explicitly set.

## Install and test

```bash
uv sync --locked
uv run pytest -q
uv run ruff check .
uv run ruff format --check .
```

Python 3.12 or newer is required.

## Run

The canonical full backfill is the restartable runner:

```bash
scripts/run_full_backfill.sh
```

By default, it covers 1 January 2020–24 July 2026, stores the manifest and
retained dataset under `data/full`, and writes all seven result Parquets under
`results/full-v1`. Networked stages retry, completed ingestion is skipped, and
valid result checkpoints are reused.

The managed runner verifies checkpoint provenance against manifest content, study
dates, processing versions, analysis code, parameters and SHA-256 input/output
hashes. Legacy readable Parquets are recomputed rather than silently reused.
Ingestion rejects incompatible completed history before downloading new releases.
See [checkpoint provenance](#checkpoint-provenance) for migration and safe updates.

### Equivalent manual stages and customization

For stage-level reruns or customization, the commands below are the equivalent
stages executed by the runner. Keep this complete path layout together so the
runner and final audit see the same artifacts.

Build the manifest:

```bash
uv run clinvar-history build-manifest \
  --start 2020-01-01 \
  --end 2026-07-24 \
  --output data/full/manifest/releases.parquet
```

Ingest both source layers sequentially because they share a single-writer
DuckDB ledger:

```bash
uv run clinvar-history ingest-vcf \
  --manifest data/full/manifest/releases.parquet \
  --root data/full

uv run clinvar-history ingest-submissions \
  --manifest data/full/manifest/releases.parquet \
  --root data/full
```

Build conservative candidate attributions and write every analysis result to
the audited results directory:

```bash
uv run clinvar-history attribute-events \
  --root data/full \
  --output results/full-v1/attributions.parquet \
  --max-lag-days 35

uv run python analysis/01_reproduce_renovo.py \
  --root data/full \
  --start 2020-01-06 \
  --end 2026-07-24 \
  --output results/full-v1/renovo_quarterly.parquet
uv run python analysis/02_sampling_bias.py \
  --root data/full \
  --output results/full-v1/sampling_bias.parquet
uv run python analysis/03_transition_matrix.py \
  --root data/full \
  --output results/full-v1/transition_matrix.parquet
uv run python analysis/04_time_to_resolution.py \
  --root data/full \
  --study-end 2026-07-24 \
  --output results/full-v1/vus_time_to_resolution.parquet
uv run python analysis/05_event_mechanisms.py \
  --root data/full \
  --attributions results/full-v1/attributions.parquet \
  --output results/full-v1/event_mechanisms.parquet
uv run python analysis/06_scv_time_to_resolution.py \
  --root data/full \
  --study-end 2026-07-24 \
  --output results/full-v1/scv_vus_time_to_resolution.parquet
```

The analysis stages only read local Parquet files. They never contact NCBI.

After the runner or the equivalent manual stages finish, independently audit
the retained artifacts and regenerate the five figures used by the final
report:

```bash
uv run python analysis/07_build_final_report.py \
  --data-root data/full \
  --results-root results/full-v1 \
  --figures-dir docs/figures/full-backfill
```

This command is also offline. It verifies manifest/ledger agreement, available
official checksums, state and episode invariants, result schemas, and the shared
479,525-transition total before writing dependency-free SVG figures. Its JSON
stdout contains the recomputed headline metrics.

## Checkpoint provenance

A managed result is a Parquet plus an adjacent `.parquet.checkpoint.json` file.
The JSON records the exact identity and SHA-256 of the output. Missing, corrupt,
modified or incompatible checkpoints are recomputed into a temporary file before
publication. A failed calculation preserves the previous published result.
The manifest has its own date-bound identity. Raw readable files are never adopted
as result checkpoints automatically.

The dataset ledger records each release's complete source identity and processing
versions. Schema, normalization and ingestion versions are explicit; parser,
normalization, diff and storage source hashes additionally guard against forgotten
version bumps. Changing source dates, URLs, sizes or checksums, removing completed
history or inserting earlier releases requires a separate data root. Appending a
chronological suffix is supported. Bump the semantic version constants in
`checkpoints.py` when persisted meanings change.

For the retained **v0.1** dataset, audit and register it once:

```bash
uv run clinvar-history adopt-legacy \
  --root data/full --results results/full-v1 \
  --start 2020-01-01 --end 2026-07-24
```

Adoption validates retained artifacts and records the declared v0.1 provenance.
It cannot prove old source-code identity from Parquet alone, and it refuses to
relabel an already versioned dataset. Legacy analysis outputs remain unverified;
the first managed analysis recomputes them. Preserve the published study with a
new result directory when exploring changed analysis settings.

To recalculate from the completed local dataset without network access:

```bash
uv run clinvar-history run-study \
  --root data/full --results results/verified-study \
  --start 2020-01-01 --end 2026-07-24 --offline
```

For a later end date, run the same command without `--offline`, set the new
`--end`, and keep the original start date. A new date horizon triggers manifest
discovery; `--refresh-manifest` also refreshes an unchanged horizon. The existing
history must still match before ingestion. Result identities include the new
horizon and input content, preventing reuse of old censored outcomes.
Networked stages retry at most three times by default; configure `--attempts` and
`--retry-delay` or the corresponding runner environment variables. Disk capacity
is checked before source ingestion. Sources stream as documented;
raw releases are not retained by the managed runner.

Manual analysis commands above remain available and write ordinary Parquets.
Use `run-study` for managed reuse; manual outputs do not acquire verified
checkpoint identities simply by occupying the expected paths.

## Exploratory analyses

The [research protocol](docs/research-protocol.md) defines the observation units,
endpoints and limitations. The commands below support bounded exploratory runs
within the July study horizon.

Run the bounded first-episode competing-risk pilot using retained episodes:

```bash
uv run python analysis/08_competing_risks.py \
  --root data/full \
  --episodes results/full-v1/vus_time_to_resolution.parquet \
  --manifest data/full/manifest/releases.parquet --layer vcv \
  --cohort-start 2020-01-07 --cohort-end 2020-12-31 \
  --study-end 2026-07-24 --max-clusters 5000 --replicates 200 \
  --output results/contributor-pilot/vcv_competing_risks.parquet
uv run python analysis/08_competing_risks.py \
  --episodes results/full-v1/scv_vus_time_to_resolution.parquet \
  --manifest data/full/manifest/releases.parquet --layer scv \
  --cohort-start 2020-02-01 --cohort-end 2020-12-31 \
  --study-end 2026-07-24 --max-clusters 5000 --replicates 200 \
  --output results/contributor-pilot/scv_competing_risks.parquet
```

The pilot estimates Aalen–Johansen cumulative incidence of first exits to B/LB,
P/LP, conflict and other states, with variant-cluster bootstrap intervals.
It reports midpoint and removal sensitivities; `--root` additionally reconstructs
monthly and quarterly VCV snapshots. The first incident episode is selected before
cohort filtering. These are observation-time research estimates, not fitted
interval-censored biological event-time distributions. Optional `--stratum`,
`--horizons`, `--seed` and resource bounds are shown in `--help`.

Run the bounded current-XML evidence review, then write the report:

```bash
uv run clinvar-history review-xml \
  --root data/full --attributions results/full-v1/attributions.parquet \
  --output results/contributor-pilot/xml --limit 12 \
  --max-record-bytes 2000000 --max-total-bytes 20000000
uv run clinvar-history compare-xml-versions \
  --variation-id 978448 --old-version 1 --new-version 2 \
  --output results/contributor-pilot/xml-version-pair
uv run python analysis/09_build_contributor_report.py \
  --pilot-root results/contributor-pilot --output results/contributor-pilot/report.md
```

Use `review-xml --offline` with the same arguments to verify and reuse cached
responses. The review preserves extra evidence and retrieval provenance, separates
germline from other assertions and leaves historical causes unresolved when only
current XML is available. `compare-xml-versions` additionally checks two explicit
VCV versions and rejects fallback to another version. The report includes a
historical submission-composition example and the limitations of version-specific
XML as a release-date snapshot. XML caches and generated result files stay outside Git.

## Why higher temporal resolution matters

The quarterly baseline compares states at the beginning and end of an
interval. A path such as `VUS → benign → VUS` is therefore invisible. Weekly
events additionally recover:

- transient conflict episodes;
- repeated and reversed transitions;
- dwell times in each state;
- less coarsely interval-censored event dates;
- separate VUS-to-benign and VUS-to-pathogenic outcomes;
- temporal evidence about the SCV mechanism behind a VCV change.

The synthetic positive control contains five true transitions. A quarterly
start/end comparison recovers one (20%); monthly sampling recovers all five.

## RENOVO baseline

Bonetti et al. reported a mean quarterly VUS reclassification rate of 0.42%,
a mean VUS discovery rate of 12.40%, and a discovery/reclassification ratio of
29.76 for their 2020–March 2024 analysis. Those are comparison targets, not
hard-coded expected values:

<https://pmc.ncbi.nlm.nih.gov/articles/PMC11293099/>

`renovo_quarterly_metrics` reproduces the paper's consecutive-snapshot
formulas: reclassification rate is resolved previous-snapshot VUS divided by
current VUS count; discovery rate is the net increase in VUS divided by
current VUS count; their quotient is the discovery/reclassification ratio.
`cohort_vus_metrics` is kept as a separate, more conventional starting-cohort
analysis. Differences in release dates, normalization, GRCh38 VCF coverage,
and treatment of composite labels must be reported rather than tuned away.

## Real-data pilot

`config/pilot.toml` records the bounded January–March 2020 pilot: four weekly
VCFs and three monthly submission summaries.

- compressed network input: 232,622,241 bytes (221.8 MiB);
- retained data: 25.1 MiB; result tables: 0.19 MiB; retained raw files: zero;
- conservative peak-disk bound: below 38 MiB, including one atomic state
  replacement;
- peak RSS: 1.60 GiB for corrected BGZF VCF ingestion, 2.57 GiB during the
  largest two-state SCV comparison, and 1.89 GiB for attribution;
- final states: 616,550 VCV rows and 1,056,294 unique SCV accessions;
- classification events: 1,372 raw-label VCV changes, of which 1,366 changed
  canonical state, and 575 SCV reclassifications;
- attribution of the 1,366 canonical VCV transitions: 1,147 SCV additions,
  93 expert-panel additions, 58 multiple-candidate events, 17 condition
  changes, 15 exact-matching SCV reclassification candidates, and 36
  unexplained events (2.64%).

The pilot cadence estimates are superseded by the full-backfill cadence
analysis above and remain only as schema and restart validation evidence. The
three-week VCF pilot is not a valid estimate of the 2020–2024 quarterly trend.
Weekly sampling recovered all 1,366 canonical transitions. Monthly and
quarterly endpoint sampling both recovered 1,359 (99.49%), missing seven
recurrent transitions and five conflict episodes; their mean release-date
error was 11.36 days. Monthly and quarterly coincide only because the
observation period was three weeks.

The exact-matching SCV label is a specific lower-bound candidate, not a causal
claim: it requires the same canonical old and new states in the SCV and VCV
events. Candidate intervals `(previous snapshot, current snapshot]` must also
overlap between the weekly VCV and monthly SCV layers; `--max-lag-days` is an
additional bound on their observation dates. Joint effects from several
submissions may nevertheless be missed.

The pilot exposed and now regression-tests three historical-source behaviors:

- BGZF VCFs contain concatenated gzip members;
- old evaluation dates may use `Jun 29, 2010` rather than ISO format;
- a few old SCVs contain duplicated missing and non-missing rows.

For duplicate accession/version groups, a non-missing classification is
preferred over missing rows; identical classifications use the latest
evaluation date. Conflicting non-missing classifications also use the strictly
latest evaluation date. A tie is retained as a synthetic `conflict` state whose
raw value contains every exact source label, rather than choosing one
arbitrarily.

## Resource profile

The audited full backfill used this storage profile:

- cumulative compressed source input: 44,250,913,390 bytes (41.21 GiB);
- retained allocated dataset: 829.5 MiB;
- seven result Parquet files: 28.6 MiB;
- raw source-release files retained: zero.

Allowing 8–10 GiB of free disk and 32 GiB of RAM remains conservative
operational guidance for a rerun. These allowances are planning margins, not
observed full-run peak-disk or peak-RSS measurements.

## Interpretation limits

- VCF includes only variants with precise genomic coordinates and is not the
  complete ClinVar corpus.
- `observed_at` is a release date, so biological change time is
  interval-censored.
- the first VCF/SCV snapshot is bootstrap inventory. Bootstrap SCV additions
  are excluded from attribution candidates, and bootstrap VUS episodes are
  excluded from time-to-resolution analyses as left-censored prevalence.
- attribution requires overlapping VCV/SCV observation intervals and applies a
  configurable maximum date lag. A single candidate is still evidence, not
  proof of causality; competing candidate types are labeled
  `multiple_candidates`.
- raw classifications are retained. Unknown/composite labels are never forced
  into the ordinal B–P scale.

Publication-grade subgroup estimates still require explicit stratification
and survival/competing-risk modelling rather than interpreting raw censored
fractions as cumulative incidence.

## Data source and licensing

The pipeline downloads source releases directly from ClinVar; this repository
does not redistribute those releases. ClinVar requests attribution when its
data are copied or distributed. Its [data-use guidance](https://www.ncbi.nlm.nih.gov/clinvar/docs/maintenance_use/)
and the broader [NCBI policies and disclaimers](https://www.ncbi.nlm.nih.gov/home/about/policies/)
describe the medical-information disclaimer and molecular-data rights
limitations that apply.

The [MIT license](LICENSE) covers only this repository's original code and
documentation. It does not license the ClinVar source releases or grant rights
that may belong to third-party submitters. Do not use this software or its
outputs for direct diagnosis or medical decision-making.

## Citation

Citation metadata are provided in [CITATION.cff](CITATION.cff). On GitHub, use
the repository's **Cite this repository** control to export a formatted
citation.

See [the data dictionary](docs/data-dictionary.md) for field-level semantics.
