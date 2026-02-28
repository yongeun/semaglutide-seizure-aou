# Semaglutide and Epilepsy/Seizure Risk: A Pharmacoepidemiologic Study Using the All of Us Research Program

## Overview

This repository contains the analysis code for a retrospective cohort study comparing the risk of epilepsy/seizure among new users of semaglutide versus other glucose-lowering drugs (other GLP-1 receptor agonists, SGLT2 inhibitors, and other glucose-lowering drugs), using electronic health record data from the [All of Us Research Program](https://allofus.nih.gov/).

## Data Access

This study uses data from the **All of Us Research Program Researcher Workbench**. The code is designed to run within the All of Us cloud-based environment using Google BigQuery. To access the data, researchers must:

1. Register as an Authorized Data User at [workbench.researchallofus.org](https://workbench.researchallofus.org/)
2. Complete the All of Us Responsible Conduct of Research Training
3. Agree to the Data User Code of Conduct
4. Create a Workspace with an appropriate research purpose

**No participant-level data is included in this repository.** Per the All of Us Data and Statistics Dissemination Policy, all aggregate statistics with cell counts fewer than 20 have been suppressed.

## Code Structure

### Primary Analysis Scripts (run in order)

| Script | Description |
|--------|-------------|
| `code/primary/01_data_extraction.py` | Extract baseline EHR demographics and survey data from BigQuery |
| `code/primary/02_preprocessing.R` | Categorize demographic/lifestyle variables; create disease phenotype flags (4 epilepsy/seizure outcome definitions) |
| `code/primary/03_cohort_iptw_cox.R` | Build drug exposure cohorts (6 pairwise comparisons), IPTW weighting, Cox proportional hazards regression |
| `code/primary/04a_tmle_iptw_all_ages.R` | IPTW + TMLE with SuperLearner (all ages) |
| `code/primary/04b_tmle_iptw_late_onset.R` | IPTW + TMLE with SuperLearner (late-onset subgroup) |
| `code/primary/04c_tmle_bootstrap.R` | Bootstrap TMLE (1,000 iterations) for confidence intervals |
| `code/primary/05_mediation_analysis.R` | Counterfactual mediation analysis (HbA1c/BMI mediators, 12-month intervals over 48 months) |
| `code/primary/06_flowchart.R` | Study cohort flowchart generation |
| `code/primary/07_km_cumulative_incidence.R` | Kaplan-Meier cumulative incidence plots |

### Supplementary Analysis Scripts

| Script | Description |
|--------|-------------|
| `code/supplementary/S01_tmle_additional_comparisons.R` | TMLE for additional drug class comparisons |
| `code/supplementary/S02_sensitivity_analyses.R` | Sensitivity analyses: 4 outcome definitions, PSM excluding severe renal/metastatic, calendar year |
| `code/supplementary/S03_subgroup_forest_plots.R` | Subgroup analyses (age, sex, race/ethnicity, obesity) with forest plots |
| `code/supplementary/S04_lmm_longitudinal.R` | Linear mixed models for BMI/HbA1c trajectories |
| `code/supplementary/S05_tmle_histogram_figure.R` | TMLE estimated probability histograms |
| `code/supplementary/S06_bootstrap_density_plots.R` | Bootstrap risk difference density plots |
| `code/supplementary/S07_cumulative_incidence_4year.R` | 4-year cumulative incidence calculations |
| `code/supplementary/S08_proportion_mediated_figure.R` | Proportion mediated over time figures |
| `code/supplementary/S09_calendar_year_sensitivity.R` | Calendar year sensitivity analysis (COVID-period imbalance) |
| `code/supplementary/S10_t2d_concept_lookup.R` | OMOP concept ID reference for T2D phenotype |

## Drug Classes

| Class | Ingredients |
|-------|------------|
| Semaglutide | semaglutide |
| Other GLP-1RA | exenatide, liraglutide, albiglutide, dulaglutide, lixisenatide |
| SGLT2i | canagliflozin, empagliflozin, dapagliflozin, ertugliflozin |
| Other GLD | glimepiride, glipizide, glyburide, alogliptin, linagliptin, saxagliptin, sitagliptin, pioglitazone, rosiglitazone |

## Study Design

- **Design**: Retrospective new-user cohort study
- **Exposure period**: January 2018 -- October 2023
- **Washout**: 183 days (6 months) with no prior use of index or comparator drugs
- **Follow-up**: Up to 48 months
- **Primary outcome**: Epilepsy or seizure (ICD-10 G40/G41/R56)
- **Statistical methods**: IPTW, Cox regression, TMLE with SuperLearner, counterfactual mediation analysis

## Environment Requirements

### R Packages
```r
bigrquery, tidyverse, glue, survival, survminer,
MatchIt, mice, tableone, lubridate, ggplot2, survey,
tmle, SuperLearner, sandwich, lmtest, lme4, lmerTest,
cowplot, doParallel, doSNOW, foreach
```

### Python Packages
```python
pandas, numpy, subprocess
```

### Platform
- All of Us Researcher Workbench (cloud R/Python environment)
- Google BigQuery access via `WORKSPACE_CDR`, `WORKSPACE_BUCKET`, `GOOGLE_PROJECT` environment variables

## All of Us Data Use

This research was conducted using data from the All of Us Research Program. The All of Us Research Program is supported by the National Institutes of Health, Office of the Director: Regional Medical Centers: 1 OT2 OD026549; 1 OT2 OD026554; 1 OT2 OD026557; 1 OT2 OD026556; 1 OT2 OD026550; 1 OT2 OD026552; 1 OT2 OD026553; 1 OT2 OD026548; 1 OT2 OD026551; 1 OT2 OD026555; IAA #: AOD 16037; Federally Qualified Health Centers: HHSN 263201600085U; Data and Research Center: 5 U2C OD023196; Biobank: 1 U24 OD023121; The Participant Center: U24 OD023176; Participant Technology Systems Center: 1 U24 OD023163; Communications and Engagement: 3 OT2 OD023205; 3 OT2 OD023206; and Community Partners: 1 OT2 OD025277; 3 OT2 OD025315; 1 OT2 OD025337; 1 OT2 OD025276.

The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

## Compliance

This code complies with the All of Us Data User Code of Conduct (V5, 10/31/2023):
- No participant-level data is included or distributed
- All aggregate statistics with cell counts < 20 are suppressed in code output
- Code uses environment variables for platform credentials (no hardcoded secrets)
- Data access requires Authorized Data User registration

## License

This code is provided for research transparency and reproducibility. See LICENSE for details.
