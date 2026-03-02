# ============================================================================
# 02_cohort_construction.R
# New-user cohort building with washout, baseline covariates, and outcome assembly
#
# This script performs:
#   1. Drug exposure index date extraction (new-user design with 183-day washout)
#   2. Pairwise cohort assembly (exclusive + sequential user handling)
#   3. Baseline covariate attachment (BMI, HbA1c, medications, comorbidities)
#   4. Outcome event assembly with censoring rules
#
# Inputs: merged_df (from 01_data_extraction_preprocessing.R)
# Outputs: 5 active-comparator cohorts ready for analysis
# ============================================================================

library(tidyverse); library(bigrquery); library(glue); library(lubridate)

# =============================================================================
# STUDY CONFIGURATION
# =============================================================================

EXPOSURE_START  <- as.Date("2018-01-01")  # Drug initiation window start
EXPOSURE_END    <- as.Date("2023-10-01")  # Drug initiation window end / data cut
WASHOUT_DAYS    <- 183                     # 6-month washout period
MAX_FOLLOWUP    <- 1460                    # 48-month maximum follow-up (days)

# ---- Drug Class Definitions ----
# Active ingredients matched to drug_exposure table via drug_concept_id
# Ingredient-level matching includes all formulations (oral, injectable, etc.)

drug_classes <- list(
  SEMAGLUTIDE = c("semaglutide"),
  OTHER_GLPA  = c("exenatide", "liraglutide", "albiglutide",
                   "dulaglutide", "lixisenatide"),
  SGLT2       = c("canagliflozin", "empagliflozin", "dapagliflozin", "ertugliflozin"),
  OtherGLD    = c("glimepiride", "glipizide", "glyburide", "alogliptin",
                   "linagliptin", "sitagliptin", "saxagliptin",
                   "pioglitazone", "rosiglitazone")
)

# ---- Five Active-Comparator Cohort Studies ----
# Primary: Semaglutide vs Other GLD, Semaglutide vs SGLT2i
# Supplementary: Other GLP-1RA vs Other GLD, Other GLP-1RA vs SGLT2i, SGLT2i vs Other GLD
comparisons <- list(
  c("SEMAGLUTIDE", "OtherGLD"),    c("SEMAGLUTIDE", "SGLT2"),
  c("OTHER_GLPA", "OtherGLD"),     c("OTHER_GLPA", "SGLT2"),
  c("SGLT2", "OtherGLD")
)

# =============================================================================
# SECTION 1: Drug Exposure Index (New-User Design)
# =============================================================================
# For each drug class:
#   index_date = MIN(drug_exposure_start_date) within [EXPOSURE_START, EXPOSURE_END]
# Washout rule:
#   No prior use of (a) the index drug AND (b) the comparator drug
#   in the 183 days before index_date
# This ensures "new user" status with respect to both drugs in the comparison

build_exposure_idx <- function(drug_concepts, dataset, washout_concepts = NULL) {
  id_str <- paste(drug_concepts$concept_id, collapse = ", ")

  # Washout includes both the index drug AND comparator drug concept IDs
  all_washout <- if (!is.null(washout_concepts)) {
    paste(unique(c(drug_concepts$concept_id, washout_concepts$concept_id)), collapse = ", ")
  } else id_str

  sql <- glue("
    WITH first_use AS (
      SELECT person_id, MIN(drug_exposure_start_date) AS index_date
      FROM `{dataset}.drug_exposure`
      WHERE drug_concept_id IN ({id_str})
        AND drug_exposure_start_date BETWEEN '{EXPOSURE_START}' AND '{EXPOSURE_END}'
      GROUP BY person_id
    ),
    prior_use AS (
      -- Participants with ANY use of index/comparator drug in washout window
      SELECT f.person_id FROM first_use f
      JOIN `{dataset}.drug_exposure` d
        ON d.person_id = f.person_id
       AND d.drug_concept_id IN ({all_washout})
       AND d.drug_exposure_start_date
           BETWEEN DATE_SUB(f.index_date, INTERVAL {WASHOUT_DAYS} DAY)
               AND DATE_SUB(f.index_date, INTERVAL 1 DAY)
      GROUP BY f.person_id
    )
    SELECT f.person_id, f.index_date
    FROM first_use f
    LEFT JOIN prior_use p ON f.person_id = p.person_id
    WHERE p.person_id IS NULL  -- Exclude anyone with prior use (washout violation)")

  download_data(sql) %>%
    mutate(person_id = as.character(person_id), index_date = as.Date(index_date))
}

# =============================================================================
# SECTION 2: Pairwise Cohort Assembly
# =============================================================================
# For each comparison (e.g., Semaglutide vs SGLT2i):
#   - Exclusive users: person_id in only one drug class -> assigned to that drug
#   - Sequential users: person_id in both drug classes -> assigned to whichever
#     started first (treatment = 1 if index drug started first, 0 if comparator)
#     Crossover date recorded for later censoring
#   - Inclusion: T2D diagnosis before index date
#   - Exclusion: pre-index epilepsy/seizure, MCI, ADRD, stroke

process_comparison <- function(exp_idx, comp_idx, base_cohort) {
  # ---- Exclusive users ----
  exp_only  <- setdiff(exp_idx$person_id, comp_idx$person_id)
  comp_only <- setdiff(comp_idx$person_id, exp_idx$person_id)

  exp_df  <- tibble(person_id = exp_only, treatment = 1L,
                    index_date = exp_idx$index_date[match(exp_only, exp_idx$person_id)])
  comp_df <- tibble(person_id = comp_only, treatment = 0L,
                    index_date = comp_idx$index_date[match(comp_only, comp_idx$person_id)])

  # ---- Sequential users ----
  # Users with prescriptions in both classes: assign to the earlier-initiated drug
  # Record crossover_date for censoring at treatment switch
  both_ids <- intersect(exp_idx$person_id, comp_idx$person_id)
  if (length(both_ids) > 0) {
    seq_df <- tibble(person_id = both_ids) %>%
      left_join(exp_idx %>% select(person_id, exp_date = index_date), by = "person_id") %>%
      left_join(comp_idx %>% select(person_id, comp_date = index_date), by = "person_id") %>%
      mutate(treatment   = ifelse(exp_date < comp_date, 1L, 0L),
             index_date  = pmin(exp_date, comp_date),
             crossover_date = pmax(exp_date, comp_date)) %>%
      select(person_id, treatment, index_date, crossover_date)
  } else {
    seq_df <- tibble()
  }

  # ---- Combine and apply inclusion/exclusion criteria ----
  cohort <- bind_rows(exp_df, comp_df, seq_df) %>%
    left_join(base_cohort, by = "person_id") %>%

    # Inclusion: T2D diagnosis on or before index date
    filter(!is.na(t2d_excluded_t1d_start_date),
           t2d_excluded_t1d_start_date <= index_date) %>%

    # Exclusion: pre-index epilepsy/seizure (outcome must occur after index)
    filter(is.na(epilepsy_or_seizure_start_date) |
           epilepsy_or_seizure_start_date >= index_date) %>%

    # Exclusion: sequential users with crossover (assign to first drug only)
    filter(is.na(crossover_date))

  cohort
}

# =============================================================================
# SECTION 3: Baseline Covariates
# =============================================================================

# ---- 3a. Baseline BMI ----
# Closest measurement in [-365, 0] days before index_date
# Valid range: 10-100 kg/m2 (extreme values excluded)
# OMOP concept_id for BMI: 3038553

add_baseline_bmi <- function(cohort, dataset) {
  sql <- glue("
    SELECT m.person_id,
           m.value_as_number AS bmi_value,
           m.measurement_date
    FROM `{dataset}.measurement` m
    WHERE m.measurement_concept_id IN (3038553)
      AND m.value_as_number BETWEEN 10 AND 100")
  bmi_raw <- download_data(sql) %>%
    mutate(person_id = as.character(person_id))

  cohort %>%
    left_join(bmi_raw, by = "person_id") %>%
    filter(measurement_date >= index_date - 365, measurement_date <= index_date) %>%
    group_by(person_id) %>%
    slice_max(measurement_date, n = 1, with_ties = FALSE) %>%
    ungroup() %>%
    mutate(
      baseline_bmi = bmi_value,
      baseline_bmi_category = cut(baseline_bmi,
        breaks = c(0, 25, 30, 35, 40, Inf),
        labels = c("<25", "25-30", "30-35", "35-40", ">=40"),
        right = FALSE)
    )
}

# ---- 3b. Baseline HbA1c ----
# Closest measurement in [-365, 0] days before index_date
# Valid range: 2-20% (extreme values excluded)
# OMOP concept_id for HbA1c: 3004410

# ---- 3c. Concomitant Medications ----
# Binary flags for any prescription in [-365, 0] days before index_date
# Drug classes matched by ingredient name in drug_exposure table:
#   Anticoagulant, Antiplatelet, Beta-blocker (BB), Biguanide (metformin),
#   CCB, Diuretic, Ezetimibe, MRA, Other antihypertensives (OtherHTN),
#   RAAS inhibitor, Statin, TZD, Insulin,
#   + treatment-specific flags: SEMAGLUTIDE, OTHER_GLPA, SGLT2i, SU, DPP4i

# ---- 3d. Comorbidities (Charlson CCI-based) ----
# Binary flags for any diagnosis within 2 years prior to index_date ([-730, 0] days)
# Identified via OMOP condition_occurrence with concept hierarchy traversal
#   myocardial_infarction, congestive_heart_failure,
#   peripheral_vascular_disease, cerebrovascular_disease,
#   chronic_pulmonary_disease, dementia, rheumatic_disease,
#   peptic_ulcer_disease, hemiplegia_or_paraplegia, hiv_infection,
#   hypoglycemia, hyperglycemic_emergency,
#   renal_disease_severity (mild/moderate vs severe),
#   liver_disease_severity (mild vs moderate/severe),
#   diabetes_with_ophthalmic_complications,
#   diabetes_with_neurological_complications,
#   malignancy_status (any malignancy vs metastatic)

# =============================================================================
# SECTION 4: Outcome Assembly
# =============================================================================
# Time-to-event construction for survival analysis
#
# Event: first occurrence of outcome phenotype AFTER index_date,
#        within MAX_FOLLOWUP (1460 days = 48 months)
#
# Censoring hierarchy (earliest of):
#   1. Death date (from aou_death table)
#   2. EXPOSURE_END (2023-10-01, data cut-off)
#   3. Last EHR date (latest activity across all EHR domains)
#   4. MAX_FOLLOWUP (index_date + 1460 days)
#   5. Crossover date (for sequential users who switched drugs)
#
# Final variables:
#   event: 0 = censored, 1 = outcome occurred
#   event_time: days from index_date to event/censor date

assemble_outcomes <- function(cohort, outcome_var = "epilepsy_or_seizure_start_date") {
  cohort %>%
    mutate(
      outcome_date = .data[[outcome_var]],

      # Event flag: outcome after index, within follow-up window
      event = as.integer(!is.na(outcome_date) & outcome_date > index_date &
                         outcome_date <= index_date + MAX_FOLLOWUP),

      # Censoring date: earliest administrative or clinical end point
      censor_date = pmin(
        coalesce(death_date, EXPOSURE_END),
        EXPOSURE_END,
        coalesce(last_ehr_date, EXPOSURE_END),
        index_date + MAX_FOLLOWUP,
        na.rm = TRUE
      ),

      # Final event date and time
      event_date = if_else(event == 1L, outcome_date, censor_date),
      event_time = as.numeric(event_date - index_date)
    ) %>%
    # Censor at crossover for sequential users
    mutate(
      event_time = if_else(!is.na(crossover_date) & crossover_date < event_date,
                           as.numeric(crossover_date - index_date), event_time),
      event = if_else(!is.na(crossover_date) & crossover_date < event_date, 0L, event)
    ) %>%
    # Exclude participants with zero or negative follow-up
    filter(event_time > 0)
}

# =============================================================================
# SECTION 5: Late-Onset Cohort Variant
# =============================================================================
# Secondary analysis: Late-onset seizure
#   - Restricted to participants age >= 60 at index date
#   - Excludes participants with ANY pre-index seizure/epilepsy
#     regardless of age (stricter than all-ages cohort)
#
# Configuration flag: exclude_preindex_outcome_any_age_in_late = TRUE
#
# Analysis approach: Intention-to-treat (ITT)
#   Follow-up from index date to outcome, death, last observation, or Oct 1, 2023
