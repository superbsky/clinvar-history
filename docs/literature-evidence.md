# Focused primary-literature assessment

Search date: 11 November 2026. This is a targeted novelty and methods assessment, not a systematic review. Searches used the public web index, PubMed records, publisher full text and source repositories. No fixed database-level sensitivity or exhaustive screening claim is made.

## Search log

| Query | Main records retained | Decision |
|---|---|---|
| ClinVar longitudinal reclassification temporal analysis 2026 competing risks Liang | Liang 2026; source repository; Zenodo archive | Closest overlap; inspect full paper and analysis code |
| ClinVar reclassification longitudinal monthly weekly snapshot temporal sampling 2024 2025 study | Bonetti 2024; Robertson 2024; ClinVar NAR update | Compare cadence and coverage |
| site.pubmed.ncbi.nlm.nih.gov variant uncertain significance reclassification cohort 2023 2024 2025 | Kobayashi 2024; Chen 2023 | Distinguish laboratory cohorts from public aggregates |
| site.ncbi.nlm.nih.gov clinvar June 2022 conflicting classifications risk allele change | Official classification documentation | Check definitions and historical changes |

## Evidence matrix

| Primary source | Relevant design and evidence | Consequence for this project |
|---|---|---|
| [Liang et al., Genome Medicine, 23 October 2026](https://link.springer.com/article/10.1186/s13073-026-01805-9) | 148 monthly XML snapshots, 2012–2025; 3.68 million variants; era/schema analysis, competing risks and mechanism validation | General longitudinal ClinVar analysis and Aalen–Johansen estimation are already published. The complementary scope is denser observation, recovered edges, individual SCV histories and explicit endpoint definitions. |
| [Bonetti et al., 2024, RENOVO](https://pmc.ncbi.nlm.nih.gov/articles/PMC11293099/) | Quarterly snapshot comparisons and prediction evaluation; 16 snapshots through March 2024 | Retain its endpoint-comparison baseline, but distinguish count recovery from exact transition-edge recall. Our analysis does not newly validate RENOVO. |
| [Kobayashi et al., 2024](https://pubmed.ncbi.nlm.nih.gov/39504018/) | Diagnostic testing cohort; over 2 million variants, initial classifications and subsequent revisions | Clinical laboratory classifications have a different sampling frame from public VCV aggregates and SCVs. Do not compare raw percentages as if populations were interchangeable. |
| [Chen et al., 2023](https://pubmed.ncbi.nlm.nih.gov/37878314/) | Hereditary testing in approximately 1.69 million people; VUS evidence and reclassification | Patient, unique-variant and public-submission denominators differ. Mean duration among resolved variants is not the cumulative incidence among all incident VUS. |
| [Robertson et al., 2024, Vigelint](https://www.nature.com/articles/s41598-024-57976-1) | ClinVar and PanelApp histories used for epilepsy-related knowledge change | Condition-specific history tools already exist. Do not claim the first useful longitudinal viewer or equate database changes with patient diagnoses. |
| [Edwards et al., 2016](https://pmc.ncbi.nlm.nih.gov/articles/PMC5557056/) | Competing-event risk estimation and estimand interpretation | Treat first benign, pathogenic and conflict exits as competing outcomes; retain explicit censoring assumptions. |
| [NCBI classification representation](https://www.ncbi.nlm.nih.gov/clinvar/docs/clinsig/) | Aggregate classifications and conflict rules, including June 2022 updates | Preserve raw labels and separate within-family labels from between-family transitions. A date near a rule change is not proof of an artifact. |

## Closest-method comparison inspected in source

The comparison uses the Liang paper and [its public code](https://github.com/DanielLiangAjj/clinvar-longitudinal), pinned to `dae628001d2b335068cef74e262e4431b288eb00`. In particular, `revision_analyses/29_km_r3.py` implements tied-time cumulative incidence; `scripts/clinvar_ftp_only_v5.py` contains its classification mapping. Its 48-event blinded two-author validation must not be equated with our automated XML audit. Our findings do not establish superiority over that study: release density, inclusion, endpoint mapping, removal handling and estimands differ. A numerical head-to-head comparison would require harmonized populations.

The defensible question is narrower: how much apparent cadence loss is attributable to reversible batches, within-family labeling, synthesized endpoint edges and changed incident-cohort membership? This is tested against this project's retained observations. It is not a claim of methodological priority or biological truth.
