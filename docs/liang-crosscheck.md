# Cross-check of released ClinVar reversal histories

Analysis date: 11 November 2026. This compares the release-level ClinVar history with the public code and selected reversal data accompanying Liang et al., [Temporal analysis of clinical variant reclassification (2012–2025)](https://link.springer.com/article/10.1186/s13073-026-01805-9). It evaluates source coverage, endpoint agreement, temporal alignment and estimator equivalence. It does not reproduce the entire paper.

## Study scope

Their study spans 148 monthly XML observations from 2012–2025 and already includes competing risks, schema analysis, disease-domain analyses and a two-reviewer mechanism assessment. Its own discussion identifies submonthly dynamics as a limitation. The following comparison distinguishes that breadth from our more focused temporal reconstruction.

| Dimension | Their work | Our July baseline / November extension | Assessment |
|---|---|---|---|
| Historical/format breadth | Longer monthly XML history and earlier schema eras | 2020–2026 precise-coordinate GRCh38 VCF plus monthly summaries; limited XML review | Their coverage is broader; our later endpoint does not compensate for missing earlier eras |
| Temporal resolution | Monthly | Actual VCF releases, then paired coarse reconstructions | Demonstrable complement; quantify additional observations instead of claiming every extra edge is a clinical reinterpretation |
| SCV unit | Aggregate work informed by submission metadata | Individual accession histories and first incident SCV VUS episodes | More explicit submission-level longitudinal unit in our project |
| Time-to-event question | First aggregate change after first variant appearance | First exit after observed incident VUS; recurrence supplementary | Different estimands. Neither is automatically superior |
| Competing risks | Already implemented | Also implemented; complete cohorts and SCV cluster dependence | Shared method, not a methodological novelty claim |
| Mechanism interpretation | Broader heuristic taxonomy and two-author manual validation | Conservative candidates, source-coverage qualifiers, dated/versioned XML | Complementary evidence; our automated case counts do not replace their human review |
| Cadence metrics | Monthly histories | Separate endpoint counts, actual matching edges, synthetic edges and changed cohort eligibility | Strong candidate contribution when tied to concrete recoverable examples |
| Domain/evidence analyses | Broader clinical-domain and evidence work | Mostly exploratory source-field strata | Their work is stronger in breadth; our strata are not an equivalent adjusted clinical analysis |
| Update engineering | Pinned parser's checkpoint check is file existence/size | November source/parameter/code/output identities and atomic, tested reuse | Practical engineering improvement; not proof of greater scientific depth |

The source locations are the pinned repository's `scripts/clinvar_ftp_only_v5.py` (classification mapping and `CheckpointManager`) and `revision_analyses/29_km_r3.py` (first-event definition, AJ and bootstrap), alongside our `research.py`, `research_sampling.py`, `research_statistics.py` and managed runner. Our mapping mirror retains their [MIT attribution](third-party-liang-license.txt).

## What was actually cross-checked

After reproducing the local research dataset as described in the [protocol](research-protocol.md), use a fresh reference checkout:

```bash
git clone https://github.com/DanielLiangAjj/clinvar-longitudinal.git artifacts/references/clinvar-longitudinal
git -C artifacts/references/clinvar-longitudinal checkout --detach dae628001d2b335068cef74e262e4431b288eb00
uv sync --locked --all-groups
TZ=UTC uv run --group research python analysis/19_liang_crosscheck.py
# Optional source revalidation: three bounded network streams, <=600 MB total.
TZ=UTC uv run python analysis/20_verify_transient_vcf_sources.py
```

Reuse the existing pinned checkout when it is already present; the command refuses changed reference source functions.

Their checkout is required at `dae628001d2b335068cef74e262e4431b288eb00`; the executed reference functions must have no diff from that commit. Only the inspected pure categorization/AJ function definitions are compiled, not their scripts' startup code or absolute local paths. The comparison mapping matches their original function on all 240 distinct raw VCF labels encountered. It is a comparison policy, not a proposed change to ingestion.

Inputs are their complete public 6,340-row reversal-case table (including all 38,284 encoded history entries), final-snapshot 3,594,583-ID origin map and published numerical outputs. Our reference dataset is the 418-source research copy, with its manifest/event hashes retained in the [comparison summary](research-tables/liang-crosscheck/summary.json). This is a complete comparison of the released selected-case table, not all events behind their paper. Full per-variant summary/history tables and snapshot checkpoints are not redistributed in the pinned repository.

### Coverage and selected-case classification agreement

At the nominal 1 July 2025 boundary, our last actual VCF release is 30 June. Our state contains 3,512,796 variants; **3,512,152 IDs** overlap their final XML population. That is 99.98% of our VCF population and 97.71% of their XML population. XML-only cases are not automatically parser errors: representation, genomic-coordinate coverage, removal and source vintage must be considered.

Of the 6,340 selected reversal cases, 6,212 are present in their final XML population and **5,936 are present in both final populations**. Thus 128 are absent from their final population and another 276 final XML cases lack a VCF counterpart. Among the 5,936 shared cases:

- 5,874 (98.96%) agree in their exact nine-category mapping;
- 5,879 (99.04%) agree after collapsing the benign and pathogenic families;
- 57 family disagreements remain explicit, rather than being counted as verified errors in either pipeline.

See [all selected cases](research-tables/liang-crosscheck/case_snapshot_comparison.csv) and [disagreement labels](research-tables/liang-crosscheck/snapshot_disagreement_labels.csv). These agreement percentages are conditional on the selected reversal table, not general classification accuracy over 3.5 million variants.

### Matching published reversal histories

There are **361 direct monthly B/LB↔P/LP edges** whose previous nominal month is at least January 2020 and current month at most July 2025. We matched the same direction and VariationID against the dense VCF history:

| Matching convention | Compatible published edges |
|---|---|
| Same calendar month as the nominal current XML month | 63 / 361 |
| Strict nominal start-to-start interval | 283 / 361 |
| Expanded nominal previous-month start through current-month end, capped at 1 July 2025 | 340 / 361 |

The large difference between these columns shows why nominal months must be aligned before calling a discrepancy. Expanded-window matches are temporal compatibility, not proof of identical observations or exact reproduction. Source dates can differ within a month. [Every matched/unmatched edge](research-tables/liang-crosscheck/published_reversal_edge_comparison.csv) is retained.

For **11 of the 21 edges without a direct native match**, the coarse start/end states still agree with our history. An intermediate state explains why a monthly endpoint edge need not be a direct observed edge. For example, VariationID **9785** is P→B in their September-to-October monthly entry, while our actual releases show:

| Actual VCF observation | Raw change |
|---|---|
| 12 September 2021 | Pathogenic → Conflicting interpretations of pathogenicity |
| 27 September 2021 | Conflicting interpretations of pathogenicity → Benign |

That is compatible with their monthly result and answers a different question from a direct P→B reinterpretation. The [remaining paths](research-tables/liang-crosscheck/unmatched_published_edge_paths.csv) and [endpoint comparison](research-tables/liang-crosscheck/reversal_endpoint_comparison.csv) expose unresolved cases for further source-level checking; agreement has not been forced by widening every window indefinitely.

### Additional content visible in the dense history

On the **same 3,512,152 final shared IDs**, over 6 January 2020–1 July 2025 and with their mapping:

| Quantity | Dense releases | Calendar-month endpoints |
|---|---|---|
| B/LB↔P/LP edges | 576 | 346 |
| Distinct variants with such an edge | 453 | 339 |
| Endpoint edges matching a direct native edge within the sample interval | Reference | 339 |
| Synthesized endpoint edges | Not applicable | 7 |

Count recovery is 60.07%, while matched-edge recall is 58.85% for this restricted transition group and time window. These are our internally paired cadence results, not estimated sensitivity of their whole paper.

**120 variants, carrying 237 native direct reversal edges, have no such reversal on our monthly grid and are absent from their released reversal list.** All 237 edges still qualify under our stricter clinical-family policy, so the candidate list does not arise solely from the substring mapping used for comparison. All 346 monthly reversal edges occur in variants contained in their list. Six monthly reversal variants have no direct native reversal, illustrating why variant counts also need care.

A post hoc concentration check materially qualifies this count: **117 of the 120 variants form exact raw-label round trips on 10 and 17 July 2023** (234 of the 237 edges). Their labels return within seven days. The 120 variants therefore do not represent 120 independent clinical reinterpretations, and the additional count is dominated by one batch. The 117 round trips were then independently checked against **three freshly downloaded full NCBI VCF streams** (2, 10 and 17 July): all pinned and current official MD5 checks passed, all 351 selected source rows matched the event labels, and every variant returned exactly to its pre-window label. No full raw archives were retained. See [source verification](research-tables/liang-crosscheck/source_vcf_verification.json) and [extracted source fields](research-tables/liang-crosscheck/verified_source_vcf_rows.csv). This rules out a reconstruction-only label discrepancy for the checked series, but its underlying mechanism has not been established. The [date counts](research-tables/liang-crosscheck/native_only_date_counts.csv) expose the concentration. This July 2023 series is separate from the better-investigated August 2022 review-status series; their causes must not be conflated.

The [120 candidate IDs](research-tables/liang-crosscheck/native_only_variants.csv), [237 raw transitions and both policies](research-tables/liang-crosscheck/native_only_primary_policy_check.csv), [dense edges](research-tables/liang-crosscheck/native_shared_reversals.csv) and [monthly edges](research-tables/liang-crosscheck/monthly_shared_reversals.csv) are exported. This is the strongest concrete additional content found in this comparison. It is still an aggregate-record observation: its clinical significance, the reason for every reversal and equivalence to the originally distributed XML bytes require validation.

### Estimator comparison

Their pure `aj_from_cells` routine and our estimator were applied to the identical complete VCV and SCV incident cohorts with the same durations and four competing causes. Maximum absolute point-curve differences were **0.0** in both comparisons. Their routine is generic in the number of causes; the common four-cause taxonomy was supplied for this controlled comparison.

This confirms shared AJ mathematics, not reproduction of their published all-entry cohort, six cause definitions or 15.5% five-year headline. Their clock starts at first variant appearance; ours at first observed incident VUS. Comparing those probabilities as if they measure the same outcome would be misleading. More bootstrap draws reduce Monte Carlo noise and are not a new statistical method.

## Interpretation and limits

Denser observation exposes intermediate states and concentrated reversible series that can change endpoint histories and cohort eligibility. The additional candidate set is dominated by one batch and does not establish independent clinical reinterpretations. The [full report](full-research-results.md) separately quantifies the August 2022 reversible window and its contribution to monthly net-count loss.

The comparison covers the released selected-case table, not all underlying events. The exact monthly source manifest, all 148 original XML snapshots and complete per-variant history outputs were unavailable for this cross-check. The 57 selected-case family disagreements, remaining unmatched paths and XML evidence for the July 2023 cases remain unresolved. A later-modified legacy XML archive cannot by itself establish equivalence to the authors' original source bytes.
