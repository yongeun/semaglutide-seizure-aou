# Semaglutide and Risk of Adult-Onset Seizure: A Target Trial Emulation

## Overview

This repository contains the analysis code for a target trial emulation examining the association between **semaglutide** (a GLP-1 receptor agonist) and the risk of **adult-onset seizure**, compared with other glucose-lowering drugs (GLDs). The study leverages electronic health record (EHR) data from the [All of Us Research Program](https://allofus.nih.gov/) (CDR version 8, Controlled Tier; 633,000+ participants, 393,596 with EHR data).

## Study Design

### Design Summary

| Parameter | Value |
|---|---|
| **Design** | Five active-comparator, new-user cohort studies using target trial emulation framework |
| **Data source** | All of Us Research Program (CDR version 8, Controlled Tier, OMOP CDM v5) |
| **Exposure window** | January 1, 2018 -- October 1, 2023 |
| **Washout period** | 183 days (6 months) before index date |
| **Maximum follow-up** | 48 months (1,460 days) |
| **Primary outcome** | Incident adult-onset seizure (excluding same-day hypoglycemic seizure events) |
| **Secondary outcome** | Late-onset seizure (age >= 60 at index date) |
| **Primary analysis** | IPTW-weighted Cox proportional hazards regression |
| **Causal inference** | TMLE (doubly robust; GLM for outcome and treatment mechanisms) |
| **Analysis approach** | Intention-to-treat |
| **Software** | R version 4.3.1 on All of Us Researcher Workbench |

### Drug Class Definitions

| Class | Drugs |
|---|---|
| **Semaglutide** | Semaglutide (oral and injectable) |
| **Other GLP-1RA** | Exenatide, liraglutide, albiglutide, dulaglutide, lixisenatide |
| **SGLT2 inhibitors** | Canagliflozin, empagliflozin, dapagliflozin, ertugliflozin |
| **Other GLD** | Glimepiride, glipizide, glyburide, alogliptin, linagliptin, sitagliptin, saxagliptin, pioglitazone, rosiglitazone |

### Five Active-Comparator Cohort Studies

**Primary comparisons:**
1. Semaglutide vs. Other GLD
2. Semaglutide vs. SGLT2 inhibitors

**Supplementary comparisons:**
3. Other GLP-1RA vs. Other GLD
4. Other GLP-1RA vs. SGLT2 inhibitors
5. SGLT2 inhibitors vs. Other GLD

### Outcome Definitions

| Outcome | Definition | Phenotype |
|---|---|---|
| **Outcome 1 (broad)** | Adult-onset seizure (broad) | ICD-10 G40 + R56 codes (sensitivity 97.7%, specificity 44.1%, PPV 96.2%) |
| **Outcome 2 (refined)** | Refined seizure | Epilepsy diagnosis (G40) OR recurrent seizure (R56 with >=2 distinct dates) |
| **Outcome 3 (specific)** | G40-only | ICD-10 G40 codes only (sensitivity 84.4%, specificity 79.4%, PPV 98.3%) |
| **Outcome 4 (primary)** | Adult-onset seizure excl. hypoglycemic | Outcome 1 minus same-day co-occurrence with hypoglycemia |

### Inclusion/Exclusion Criteria

**Inclusion:**
- Age >= 18 at index date
- EHR data available (`has_ehr_data = 1`)
- Type 2 diabetes diagnosis (excluding T1D) before index date
- New user of index drug within the exposure window

**Exclusion:**
- Any use of the index drug or active comparator in the 183-day washout window
- Pre-index epilepsy or seizure diagnosis
- Pre-index mild cognitive impairment (MCI)
- Pre-index Alzheimer's disease or related dementias (ADRD)
- Pre-index stroke

## Code Structure

| Script | Purpose | Inputs | Outputs |
|---|---|---|---|
| `01_cohort_definition.R` | Study parameters, drug class definitions, variable categorization rules, eligibility criteria | — | Study configuration |
| `01_data_extraction_preprocessing.R` | BigQuery extraction, variable categorization, disease phenotyping | OMOP CDM tables via BigQuery | Preprocessed cohort dataframe |
| `02_cohort_construction.R` | New-user design with washout, baseline covariates, outcome assembly | Preprocessed data + drug exposures | Analysis-ready cohort per comparison |
| `02_primary_analysis.R` | IPTW weighting, Cox regression, TMLE | Cohort dataframes | HR, RD, RR, OR, NNT estimates |
| `03_iptw_cox.R` | Detailed IPTW pipeline: MICE, PS estimation, balance, weighted Cox | Cohort dataframes | IPTW results, balance tables, baseline characteristics |
| `03_secondary_analyses.R` | Mediation, sensitivity, subgroup, LMM trajectory analyses | IPTW-weighted cohorts | Mediation effects, forest plots, LMM results |
| `phenotype_concept_ids.R` | All OMOP Concept IDs for outcome, exclusion, and diabetes phenotypes | — | Concept ID vectors |

### Analysis Workflow

```
┌──────────────────────────────────────────────────────────────────┐
│  Step 1: Data Extraction & Preprocessing                         │
│  01_data_extraction_preprocessing.R                              │
│  - BigQuery: demographics, PPI surveys, EHR dates, death data   │
│  - Variable categorization (race, age, smoking, AUDIT-C, etc.)  │
│  - Disease phenotype flags via OMOP concept hierarchies          │
│  Output: ehr_df → mother_df → merged_df                         │
└──────────────────────┬───────────────────────────────────────────┘
                       ▼
┌──────────────────────────────────────────────────────────────────┐
│  Step 2: Cohort Construction                                     │
│  02_cohort_construction.R                                        │
│  - New-user index dates with 183-day washout                     │
│  - Sequential user handling (assign to earlier drug)             │
│  - Baseline BMI/HbA1c, medications (1yr), comorbidities (2yr)   │
│  - Outcome assembly with censoring (ITT approach)               │
│  Output: 5 active-comparator cohorts                            │
└──────────────────────┬───────────────────────────────────────────┘
                       ▼
┌──────────────────────────────────────────────────────────────────┐
│  Step 3: Primary Analysis                                        │
│  02_primary_analysis.R + 03_iptw_cox.R                           │
│  - MICE imputation for BMI/HbA1c (m=20)                         │
│  - Propensity score estimation (46 covariates)                   │
│  - IPTW weights (ATE, trimmed at 1st/99th %ile)                 │
│  - Balance assessment (SMD < 0.1 target)                         │
│  - Cox PH regression (robust SE, clustered by person_id)        │
│  - TMLE doubly robust (GLM, 1,000 bootstrap CIs)                │
│  Output: HR, RD, RR, OR, NNT with 95% CIs                      │
└──────────────────────┬───────────────────────────────────────────┘
                       ▼
┌──────────────────────────────────────────────────────────────────┐
│  Step 4: Secondary & Sensitivity Analyses                        │
│  03_secondary_analyses.R                                         │
│  - Counterfactual mediation (HbA1c, BMI; 12-month windows)      │
│  - Sensitivity: 4 outcome definitions, calendar year, PSM,       │
│    min 2-year follow-up, leave-one-out, E-values                 │
│  - Subgroup: age (<60/≥60), sex, race/ethnicity, obesity         │
│  - LMM: HbA1c and BMI trajectories over 48 months               │
│  Output: Mediation effects, forest plots, trajectory estimates   │
└──────────────────────────────────────────────────────────────────┘
```

## Statistical Methods

### Propensity Score Variables (46 covariates)

**Demographics (11):** age, sex, race/ethnicity, income, education, insurance, smoking, alcohol category, baseline BMI category, baseline HbA1c, index year group

**Concomitant Medications (18):** Anticoagulant, Antiplatelet, Beta-blocker, Biguanide, CCB, Diuretic, Ezetimibe, MRA, Other antihypertensives, RAAS inhibitor, Statin, TZD, Insulin, Semaglutide, Other GLP-1RA, SGLT2i, Sulfonylurea, DPP-4i

**Comorbidities (17, Charlson-based):** Myocardial infarction, congestive heart failure, peripheral vascular disease, cerebrovascular disease, chronic pulmonary disease, dementia, rheumatic disease, peptic ulcer disease, hemiplegia/paraplegia, HIV, hypoglycemia, hyperglycemic emergency, renal disease severity, liver disease severity, diabetic ophthalmic complications, diabetic neurological complications, malignancy status

> **Note:** Comparison-specific PS models exclude the treatment/comparator drug variables (e.g., `SEMAGLUTIDE` and `SGLT2i` are excluded from the PS model for the Semaglutide vs. SGLT2i comparison). Zero-variance covariates within each comparison are also excluded.

### IPTW (Inverse Probability of Treatment Weighting)

- **PS model:** Logistic regression: `treatment ~ 46 covariates`
- **Weight type:** Average treatment effect (ATE): `w = T/PS + (1-T)/(1-PS)`
- **Trimming:** Weights capped at 1st and 99th percentiles
- **Standardization:** `w_std = w_trimmed * N / sum(w_trimmed)`
- **Balance criterion:** All standardized mean differences (SMD) < 0.1 after weighting

### Multiple Imputation (MICE)

- **Variables imputed:** `baseline_bmi`, `baseline_hba1c`
- **Number of imputations:** m = 20
- **Seed:** 123
- **Predictors:** Treatment, demographics, concomitant medications, comorbidities (~35 predictor variables)
- **Fallback:** Mean imputation if MICE fails to converge

### Cox Proportional Hazards Regression

- **Model:** `Surv(event_time, event) ~ treatment` with IPTW weights
- **Robust SE:** Sandwich estimator with `cluster(person_id)`
- **Output:** Hazard Ratio (HR), 95% CI, p-value
- **Scope:** 5 comparisons x 4 outcome definitions x 2 cohorts (all ages, late-onset)

### TMLE (Targeted Maximum Likelihood Estimation)

- **Method:** Doubly robust estimation; fit on unweighted cohorts
- **Outcome mechanism:** GLM (generalized linear model)
- **Treatment mechanism:** GLM (generalized linear model)
- **Bootstrap:** 1,000 resamples for confidence intervals (seed = 12345)
- **Effect measures:** Risk Difference (ATE), Risk Ratio (RR), Odds Ratio (OR), Number Needed to Treat (NNT = 1/|ATE|)
- **Time-specific RD:** Kaplan-Meier survival estimates at 6, 12, 18, and 24 months

### Counterfactual Mediation Analysis (Vansteelandt-style)

- **Mediators:** HbA1c and BMI (analyzed separately)
- **Time windows:** M1 (0-12 mo), M2 (12-24 mo), M3 (24-36 mo), M4 (36-48 mo)
- **Mediator models (chained):**
  - M1 ~ treatment + X (baseline value)
  - M2 ~ treatment + M1
  - M3 ~ treatment + M2
  - M4 ~ treatment + M3
- **Missing mediators:** Last Observation Carried Forward (LOCF)
- **Outcome model:** Piecewise Cox PH with time split at 365, 730, 1095 days: `Surv(tstart, tstop, event) ~ treatment + M_period + X + strata(period) + cluster(person_id)`
- **Counterfactual paths:**
  - S11: treatment=1, mediator under treatment=1 (natural course)
  - S10: treatment=1, mediator under treatment=0 (counterfactual)
  - S00: treatment=0, mediator under treatment=0 (natural course)
- **Monte Carlo simulation:** B = 200 draws per mediator path
- **Proportion mediated:** PM = (S11 - S10) / (S11 - S00)
- **Parametric bootstrap:** R = 100 replicates for uncertainty quantification (seed = 2025)

### Sensitivity Analyses

- **Multiple outcome definitions:** All 4 outcomes through IPTW + Cox pipeline
- **Calendar year adjustment:** Index year added to PS model; exact matching on year; pre-COVID (2018-2019) restricted analysis
- **Minimum follow-up restriction:** Restricted to participants with >= 2 years of follow-up
- **Leave-one-out:** Systematic removal of one drug ingredient at a time from comparator class to assess robustness
- **E-values:** Quantification of unmeasured confounding strength needed to explain away observed associations
- **Propensity score matching (PSM):** 1:1 and 1:5 nearest-neighbor matching (caliper = 0.2 x SD of logit PS), with optional exclusion of severe renal disease or metastatic malignancy

### Subgroup Analyses

| Subgroup | Levels |
|---|---|
| **Age** | <60 vs. >=60 |
| **Sex** | Male vs. Female |
| **Race/Ethnicity** | Non-Hispanic White, Non-Hispanic Black, Hispanic, Other |
| **Obesity** | BMI <30 vs. >=30 |

- IPTW re-estimated within each subgroup
- Interaction test: `treatment * subgroup_variable`
- Forest plots with HR and 95% CI

### Longitudinal Biomarker Trajectories (LMM)

- **Outcomes:** BMI and HbA1c trajectories over 48 months
- **Model:** `value ~ treatment * time_years + baseline_value + (1 + time_years | person_id)`
- **Key estimate:** Treatment x time interaction (differential slope)
- **Fallback:** Random intercept-only model if convergence fails

## Variable Coding Reference

| Variable | Coding |
|---|---|
| `sex_cat` | 0 = Male, 1 = Female |
| `raceethnicity_cat` | 0 = NH-White, 1 = NH-Black, 2 = Hispanic, 3 = Other |
| `income` | 0 = <$35K, 1 = $35-100K, 2 = >$100K |
| `education` | 0 = HS or below, 1 = Some college/graduate, 2 = Advanced degree |
| `insurance_category` | 0 = None, 1 = Public, 2 = Private |
| `smoking` | 0 = Never, 1 = Former, 2 = Current |
| `alcohol_category` | AUDIT-C: 0 = Low risk (0-4), 1 = Increased (5-7), 2 = High (8-10), 3 = Dependent (11-12) |
| `baseline_bmi_category` | <25, 25-30, 30-35, 35-40, >=40 |
| `index_year_grouped` | Non-COVID (2018-2019, 2023) vs. COVID (2020-2022) |
| `treatment` | 1 = index drug, 0 = comparator |
| All medication/comorbidity flags | 0 = absent, 1 = present |
| Missing values | Coded as 999 |

## Data Access

This study uses data from the **All of Us Researcher Workbench**. Researchers must register as Authorized Data Users at [workbench.researchallofus.org](https://workbench.researchallofus.org/). **No participant-level data is included in this repository.**

All queries run against the Controlled Tier CDR dataset (`WORKSPACE_CDR`) via the Google BigQuery interface in the Researcher Workbench environment.

## Requirements

### R Packages

```r
# Core
bigrquery, tidyverse, lubridate, glue

# Statistical analysis
survival, survey, mice, tableone, sandwich, lmtest, MatchIt

# Causal inference
tmle, SuperLearner

# Longitudinal models
lme4, lmerTest

# Parallel computation
future, future.apply

# Visualization
ggplot2

# E-value computation
EValue  # or manual computation
```

### Platform

- R version 4.3.1
- All of Us Researcher Workbench (Jupyter environment with R kernel)
- Google Cloud BigQuery access (provided by Workbench)
- Environment variables: `GOOGLE_PROJECT`, `WORKSPACE_CDR`, `WORKSPACE_BUCKET`

## Reproducibility Notes

- All seeds are fixed: MICE (123), TMLE bootstrap (12345), mediation bootstrap (2025)
- MICE uses m = 20 imputations; first completed dataset used for primary analysis
- IPTW weights are trimmed at 1st/99th percentiles and standardized to sample size
- TMLE uses GLM for both outcome and treatment mechanisms; 1,000 bootstrap resamples
- Mediation uses 200 Monte Carlo draws per counterfactual mediator path; 100 parametric bootstrap replicates
- Comorbidities assessed up to 2 years prior to index date; medications up to 1 year prior
- Intention-to-treat approach for primary analysis
- All scripts are designed to run sequentially (01 → 02 → 03)

## All of Us Data Use

This research was conducted using data from the All of Us Research Program, supported by the National Institutes of Health, Office of the Director. The content is solely the responsibility of the authors and does not necessarily represent the official views of the National Institutes of Health.

## License

This code is provided for research transparency and reproducibility. See the All of Us Data Use Policy for restrictions on data use.
