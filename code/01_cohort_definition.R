# ============================================================================
# 01_cohort_definition.R
# Cohort construction for semaglutide vs comparator GLD study
# Platform: All of Us Researcher Workbench (BigQuery)
# ============================================================================

library(bigrquery); library(tidyverse); library(lubridate); library(glue)

# --- Study Parameters --------------------------------------------------------

EXPOSURE_START  <- as.Date("2018-01-01")
EXPOSURE_END    <- as.Date("2023-10-01")
WASHOUT_DAYS    <- 183  # 6-month washout for new-user design
MAX_FOLLOWUP    <- 1460 # 48 months

# --- Drug Class Definitions --------------------------------------------------

drug_classes <- list(
  SEMAGLUTIDE = c("semaglutide"),
  OTHER_GLPA  = c("exenatide", "liraglutide", "albiglutide",
                   "dulaglutide", "lixisenatide"),
  SGLT2       = c("canagliflozin", "empagliflozin",
                   "dapagliflozin", "ertugliflozin"),
  OtherGLD    = c("glimepiride", "glipizide", "glyburide",
                   "alogliptin", "linagliptin", "sitagliptin",
                   "saxagliptin", "pioglitazone", "rosiglitazone")
)

# --- Pairwise Comparisons (6 total) -----------------------------------------

comparisons <- list(
  c("SEMAGLUTIDE", "OTHER_GLPA"),
  c("SEMAGLUTIDE", "SGLT2"),
  c("SEMAGLUTIDE", "OtherGLD"),
  c("OTHER_GLPA",  "SGLT2"),
  c("OTHER_GLPA",  "OtherGLD"),
  c("SGLT2",       "OtherGLD")
)

# --- 1. Data Extraction (BigQuery) -------------------------------------------
# Extract from WORKSPACE_CDR:
#   - person table: person_id, birth_datetime, race, sex, ethnicity
#   - observation table (PPI surveys): income, education, smoking, alcohol (AUDIT-C),
#     insurance type
#   - Restrict to participants with has_ehr_data = 1

# --- 2. Variable Categorization ----------------------------------------------
# Demographics:
#   age_group: 18-44, 45-64, >=65 (reference date: 2023-10-01)
#   raceethnicity_cat: 0=NH-White, 1=NH-Black, 2=Hispanic, 3=NH-Other
#   sex_cat: 0=Male, 1=Female
#
# Socioeconomic:
#   income: 0=<$35k, 1=$35-100k, 2=>$100k
#   education: 0=HS or below, 1=Some college/graduate, 2=Advanced degree
#   insurance_category: 0=None, 1=Public, 2=Private
#
# Lifestyle:
#   smoking: 0=Never, 1=Former, 2=Current (from 100-cig lifetime + current frequency)
#   alcohol_category: AUDIT-C score -> 0=Low risk (0-4), 1=Increased (5-7),
#                     2=High (8-10), 3=Dependent (11-12)

# --- 3. Disease Phenotype Flags (from cb_search_all_events) ------------------
# See phenotype_concept_ids.R for OMOP concept IDs
# Outcomes (4 definitions):
#   1. epilepsy_or_seizure: Broad (G40 + R56 codes)
#   2. epilepsy_refined: Epilepsy dx OR recurrent seizure (>=2 distinct dates)
#   3. epilepsy_g40: G40 codes only (specificity 79.4%, PPV 98.3%)
#   4. epilepsy_or_seizure_excl_hypoglycemic: Outcome 1 minus same-day hypoglycemia
#
# Exclusion phenotypes: mci, adrd, stroke
# Comorbidities (Charlson-based): See PS variable list in 02_primary_analysis.R

# --- 4. EHR Date Ranges & Death Data ----------------------------------------
# Per-table date ranges extracted: measurement, condition_occurrence,
#   drug_exposure, procedure_occurrence, observation, visit_occurrence
# Death data: aou_death table (censoring)

# --- 5. Drug Exposure Index Date (New-User Design) ---------------------------
# For each drug class:
#   index_date = first prescription date within [EXPOSURE_START, EXPOSURE_END]
#   Washout: No use of index drug OR comparator drug in prior 183 days
#   Exclusions: prior epilepsy/seizure, MCI, ADRD, stroke before index date

# --- 6. Baseline Covariates (measured at index_date) -------------------------
# BMI: closest measurement within [-365, 0] days of index
#   baseline_bmi_category: <25, 25-30, 30-35, 35-40, >=40
# HbA1c: closest measurement within [-365, 0] days of index
# Concomitant medications (binary flags, any use in [-365, 0]):
#   Anticoagulant, Antiplatelet, BB, Biguanide, CCB, Diuretic,
#   Ezetimibe, MRA, OtherHTN, RAAS, Statin, TZD, Insulin
# Comorbidity flags: binary, any occurrence before index_date

# --- 7. Outcome Ascertainment -----------------------------------------------
# Event: first occurrence of outcome phenotype AFTER index_date
# Censoring: min(death_date, EXPOSURE_END, last_EHR_date, MAX_FOLLOWUP)
# Time-to-event: event_date - index_date (days)
