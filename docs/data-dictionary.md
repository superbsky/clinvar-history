# ClinVar history data dictionary

All classifications retain both the source label and a normalized state:
`benign`, `likely_benign`, `vus`, `likely_pathogenic`, `pathogenic`,
`conflict`, `other`, or `missing`. `classification_ambiguous=true` marks a
source value that was not mapped onto the recognized germline states.

## Release manifest

`data/manifest/releases.parquet`

| Field | Type | Meaning |
|---|---|---|
| `release_date` | date | VCF filename date; first day of `source_month` for monthly summaries |
| `source_month` | string | Stable `YYYY-MM` key for monthly submission summaries |
| `kind` | string | `vcf` or `submission_summary` |
| `url` | string | Immutable NCBI archive URL |
| `size_bytes` | integer | Compressed size reported by the archive index |
| `md5` | string | Official sidecar checksum when available |

The FTP `Released` column is not used as the monthly time coordinate because
archived files may be republished later. For example, 2020 submission files
currently carry 2021 archive timestamps. `source_month` is the stable
longitudinal key; its day-one `release_date` is a nominal observation date.

## Latest VCV state

`data/state/vcv_latest.parquet`

| Field | Meaning |
|---|---|
| `variation_id` | ClinVar VariationID; aggregate stable key |
| `allele_id` | ClinVar AlleleID when present |
| `chromosome`, `position`, `reference`, `alternate` | GRCh38 VCF locus |
| `raw_classification` | Exact decoded `CLNSIG` |
| `classification` | Normalized germline state |
| `classification_ambiguous` | Unknown/composite normalization flag |
| `review_status` | Decoded `CLNREVSTAT` |
| `gene_info` | Source `GENEINFO` |
| `conditions` | `CLNDISDB`, falling back to `CLNDN` |
| `variant_type` | Source `CLNVC` |

## Latest SCV state

`data/state/scv_latest.parquet`

| Field | Meaning |
|---|---|
| `scv_accession` | SCV accession without version; longitudinal key |
| `scv_version` | Current integer version |
| `variation_id` | Aggregate VariationID represented by the submission |
| `raw_classification` | Exact `ClinicalSignificance` label |
| `classification` | Normalized germline state |
| `classification_ambiguous` | Unknown/composite normalization flag |
| `last_evaluated` | Submitter evaluation date, if reported |
| `review_status` | Submission review status |
| `conditions` | Reported phenotype, falling back to submitted phenotype |
| `submitter` | Submitting organization |
| `collection_method` | Reported collection method |

## VCV events

`data/events/vcv/release_date=YYYY-MM-DD/events.parquet`

Event types:

- `variant_added`;
- `variant_removed`;
- `classification_changed`;
- `review_status_changed`;
- `condition_set_changed`.

Each row contains `variation_id`, `event_type`, `observed_at`, and old/new
classification, raw classification, review status, conditions, gene
information, and variant type. Unchanged records are absent. The first
partition represents bootstrap inventory rather than known creation times.

## SCV events

`data/events/scv/release_date=YYYY-MM-DD/events.parquet`

Event types:

- `scv_added`;
- `scv_removed`;
- `scv_version_changed`;
- `scv_reclassified`;
- `review_status_changed`;
- `condition_set_changed`.

`scv_version_changed` is deliberately separate from `scv_reclassified`.
A reclassification requires the classification value itself to change.
Fields preserve accession, VariationID, observation date, versions, and old/new
classification, review, condition, submitter, and collection-method values.

## Attributed aggregate events

`results/attributions.parquet`

| Field | Meaning |
|---|---|
| `variation_id` | Changed aggregate record |
| `observed_at` | VCV release date |
| `primary_mechanism` | Unique candidate, `multiple_candidates`, or `unexplained` |
| `candidate_mechanisms` | All distinct mechanism types with compatible observation intervals |
| `supporting_event_ids` | Stable SCV event identifiers supporting candidates |

Candidate mechanisms are `exact_matching_scv_reclassification`,
`expert_panel_added`, `scv_added`, `scv_removed`,
`review_precedence_changed`, and `condition_changed`. The exact-match label
requires an SCV event with the same canonical old and new states as the VCV
event. It is a deliberately specific lower-bound candidate: it can miss real
SCV contributions when several submissions jointly change the aggregate.
First-snapshot `scv_added` rows are excluded because they describe baseline
inventory. A candidate requires overlap between the VCV and SCV
interval-censoring windows `(previous snapshot, current snapshot]`, plus the
configured maximum lag between observation dates. Attribution is evidentiary
and must not be reported as proven causality.

## Processing ledger

`data/clinvar.duckdb`, table `pipeline_runs`

| Field | Meaning |
|---|---|
| `kind`, `url` | Composite release key |
| `release_date` | Source release date |
| `status` | `processing`, `completed`, or `failed` |
| `started_at`, `completed_at` | Pipeline timestamps |
| `error` | Last failure message |
| `source_md5` | Checksum calculated while streaming the compressed source |
| `official_md5` | Sidecar checksum recorded when processing began |

State and event files are written through a temporary file followed by an
atomic rename. Each completed release also has an immutable state checkpoint;
the latest-state path is advanced only after its event partition is durable.
A completed URL is skipped; a failed URL is retried. Completed releases must
form a prefix of the manifest, preventing an accidental historical backfill
from rewinding state.

## Analysis outputs

Default versioned location: `results/v1/`.

- `renovo_quarterly.parquet`: consecutive-snapshot VUS counts,
  reclassification/discovery rates, and their ratio;
- `sampling_bias.parquet`: true/observed transitions, event recall, missed
  recurrent transitions, missed conflict episodes, and mean event-date error
  by cadence;
- `transition_matrix.parquet`: directed state-pair counts;
- `vus_time_to_resolution.parquet`: episode dates, duration, benign,
  pathogenic, conflict, other-exit, or censored outcome, plus review status,
  gene, condition, and variant type for stratification;
- `scv_vus_time_to_resolution.parquet`: submitter-level VUS episodes and
  outcomes, retaining submitter, collection method, review status, and
  condition for stratification;
- `event_mechanisms.parquet`: counts of conservative primary mechanisms.

State episodes use half-open intervals `[started_at, ended_at)`. Active
episodes are right-censored at the configured study end; removal is censored
at its observation date. Conflict and other non-VUS exits are explicit
observed outcomes rather than censoring. VUS inventory in the first snapshot
is excluded from resolution outputs because its true onset predates the study.
