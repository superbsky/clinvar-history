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

## Full-backfill results

The completed 6 January 2020–24 July 2026 backfill produced these audited
headlines:

- 417 releases total (338 VCF + 79 SCV), all completed;
- 479,525 canonical VCV transitions;
- aggregate VUS count increased from 217,177 to 2,315,022 (10.66x);
- against the release-level event log, transition recall was 99.91% for the
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

Full ClinVar XML is deliberately excluded. Compressed source releases are
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

Checkpoint reuse is only structural/readability-based, not semantic validation:
a result checkpoint is reused when its Parquet is readable, and ingestion may
reuse an existing manifest. These checks do not prove that artifacts match a
changed date horizon or configuration. For such changes, use fresh
manifest/data/result paths or remove the existing manifest and checkpoints
before rerunning.

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
