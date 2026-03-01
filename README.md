# Semaglutide and Epilepsy/Seizure Risk in the All of Us Research Program

## Overview

Analysis code for a retrospective new-user cohort study comparing epilepsy/seizure risk among users of semaglutide versus other glucose-lowering drugs, using electronic health record data from the [All of Us Research Program](https://allofus.nih.gov/).

## Data Access

This study uses data from the All of Us Researcher Workbench. Researchers must register as Authorized Data Users at [workbench.researchallofus.org](https://workbench.researchallofus.org/). No participant-level data is included in this repository.

## Code

| File | Description |
|------|-------------|
| `code/01_cohort_definition.R` | Study parameters, drug class definitions, variable categorization, cohort construction |
| `code/02_primary_analysis.R` | IPTW weighting, Cox proportional hazards, TMLE with SuperLearner |
| `code/03_secondary_analyses.R` | Counterfactual mediation, sensitivity analyses, subgroup analyses, LMM trajectories |
| `code/phenotype_concept_ids.R` | OMOP concept IDs for outcome and exclusion phenotypes |

## Study Design

- **Design**: New-user retrospective cohort
- **Exposure period**: January 2018 -- October 2023
- **Washout**: 183 days
- **Follow-up**: Up to 48 months
- **Comparisons**: Semaglutide vs Other GLP-1RA, SGLT2i, and Other GLD (6 pairwise)
- **Primary outcome**: Epilepsy or seizure
- **Methods**: IPTW, Cox regression, TMLE, counterfactual mediation analysis

## Requirements

R packages: `bigrquery`, `tidyverse`, `survival`, `survey`, `mice`, `tmle`, `SuperLearner`, `lme4`, `lmerTest`

## All of Us Data Use

This research was conducted using data from the All of Us Research Program, supported by the National Institutes of Health, Office of the Director. The content is solely the responsibility of the authors and does not necessarily represent the official views of the NIH.

## License

This code is provided for research transparency and reproducibility.
