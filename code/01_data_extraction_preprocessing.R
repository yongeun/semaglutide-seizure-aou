# ============================================================================
# 01_data_extraction_preprocessing.R
# Data extraction from All of Us BigQuery and variable categorization
# Platform: All of Us Researcher Workbench
#
# This script performs:
#   1. EHR demographics + PPI survey extraction from BigQuery
#   2. Variable categorization (race, age, smoking, alcohol, SES)
#   3. Disease phenotype flag generation via OMOP concept hierarchies
#   4. EHR date range extraction and death data integration
#   5. BMI and HbA1c panel data extraction
#
# Outputs: ehr_df, mother_df, merged_df, bmi_panel, a1c_panel
# ============================================================================

library(bigrquery); library(tidyverse); library(lubridate); library(glue)

project <- Sys.getenv("GOOGLE_PROJECT")
cdr_ds  <- Sys.getenv("WORKSPACE_CDR")
bucket  <- Sys.getenv("WORKSPACE_BUCKET")

download_data <- function(sql) {
  bq_table_download(bq_project_query(project, sql))
}

# =============================================================================
# SECTION 1: Extract EHR Demographics + PPI Survey Data
# =============================================================================
# Source: person table + observation table (PPI survey responses)
# Restricted to participants with EHR data (has_ehr_data = 1)
#
# PPI survey Concept IDs:
#   1585386  - Health Insurance (Yes/No)
#   43528428 - Insurance Type Update (Medicare, Medicaid, Employer, etc.)
#   1585375  - Annual Income
#   1585940  - Education (Highest Grade)
#   1585857  - 100 Cigarettes Lifetime (Yes/No)
#   1585860  - Current Smoking Frequency (Every Day / Some Days / Not At All)
#   1586198  - Alcohol Participant (Yes/No)
#   1586201  - Drink Frequency Past Year
#   1586207  - Average Daily Drink Count
#   1586213  - Heavy Drinking Frequency (6+ drinks)

ehr_query <- glue("
  WITH ehr AS (
    SELECT p.person_id, p.birth_datetime AS date_of_birth,
           c_race.concept_name AS race, c_sex.concept_name AS sex,
           c_ethn.concept_name AS ethnicity
    FROM `{cdr_ds}.person` p
    LEFT JOIN `{cdr_ds}.concept` c_race ON p.race_concept_id = c_race.concept_id
    LEFT JOIN `{cdr_ds}.concept` c_sex  ON p.sex_at_birth_concept_id = c_sex.concept_id
    LEFT JOIN `{cdr_ds}.concept` c_ethn ON p.ethnicity_concept_id = c_ethn.concept_id
    WHERE p.person_id IN (
      SELECT DISTINCT person_id FROM `{cdr_ds}.cb_search_person` WHERE has_ehr_data = 1
    )
  )
  SELECT ehr.*,
    ins1.aname AS ins1, ins2.aname AS ins2,
    obs1.aname AS aname1, obs2.aname AS aname2, obs3.aname AS aname3,
    obs4.aname AS aname4, obs5.aname AS aname5, obs6.aname AS aname6,
    obs7.aname AS aname7, obs8.aname AS aname8
  FROM ehr
  -- Insurance (1585386 = Yes/No, 43528428 = Type)
  LEFT JOIN (SELECT o.person_id, a.concept_name AS aname
             FROM `{cdr_ds}.observation` o
             LEFT JOIN `{cdr_ds}.concept` a ON o.value_source_concept_id = a.concept_id
             WHERE o.observation_source_concept_id = 1585386) ins1 ON ehr.person_id = ins1.person_id
  LEFT JOIN (SELECT o.person_id, a.concept_name AS aname
             FROM `{cdr_ds}.observation` o
             LEFT JOIN `{cdr_ds}.concept` a ON o.value_source_concept_id = a.concept_id
             WHERE o.observation_source_concept_id = 43528428) ins2 ON ehr.person_id = ins2.person_id
  -- Income (1585375)
  LEFT JOIN (SELECT o.person_id, a.concept_name AS aname
             FROM `{cdr_ds}.observation` o
             LEFT JOIN `{cdr_ds}.concept` a ON o.value_source_concept_id = a.concept_id
             WHERE o.observation_source_concept_id = 1585375) obs1 ON ehr.person_id = obs1.person_id
  -- Education (1585940)
  LEFT JOIN (SELECT o.person_id, a.concept_name AS aname
             FROM `{cdr_ds}.observation` o
             LEFT JOIN `{cdr_ds}.concept` a ON o.value_source_concept_id = a.concept_id
             WHERE o.observation_source_concept_id = 1585940) obs2 ON ehr.person_id = obs2.person_id
  -- 100 Cigarettes Lifetime (1585857)
  LEFT JOIN (SELECT o.person_id, a.concept_name AS aname
             FROM `{cdr_ds}.observation` o
             LEFT JOIN `{cdr_ds}.concept` a ON o.value_source_concept_id = a.concept_id
             WHERE o.observation_source_concept_id = 1585857) obs3 ON ehr.person_id = obs3.person_id
  -- Smoking Frequency (1585860)
  LEFT JOIN (SELECT o.person_id, a.concept_name AS aname
             FROM `{cdr_ds}.observation` o
             LEFT JOIN `{cdr_ds}.concept` a ON o.value_source_concept_id = a.concept_id
             WHERE o.observation_source_concept_id = 1585860) obs4 ON ehr.person_id = obs4.person_id
  -- Alcohol Participant (1586198)
  LEFT JOIN (SELECT o.person_id, a.concept_name AS aname
             FROM `{cdr_ds}.observation` o
             LEFT JOIN `{cdr_ds}.concept` a ON o.value_source_concept_id = a.concept_id
             WHERE o.observation_source_concept_id = 1586198) obs5 ON ehr.person_id = obs5.person_id
  -- Drink Frequency Past Year (1586201)
  LEFT JOIN (SELECT o.person_id, a.concept_name AS aname
             FROM `{cdr_ds}.observation` o
             LEFT JOIN `{cdr_ds}.concept` a ON o.value_source_concept_id = a.concept_id
             WHERE o.observation_source_concept_id = 1586201) obs6 ON ehr.person_id = obs6.person_id
  -- Average Daily Drink Count (1586207)
  LEFT JOIN (SELECT o.person_id, a.concept_name AS aname
             FROM `{cdr_ds}.observation` o
             LEFT JOIN `{cdr_ds}.concept` a ON o.value_source_concept_id = a.concept_id
             WHERE o.observation_source_concept_id = 1586207) obs7 ON ehr.person_id = obs7.person_id
  -- Heavy Drinking 6+ (1586213)
  LEFT JOIN (SELECT o.person_id, a.concept_name AS aname
             FROM `{cdr_ds}.observation` o
             LEFT JOIN `{cdr_ds}.concept` a ON o.value_source_concept_id = a.concept_id
             WHERE o.observation_source_concept_id = 1586213) obs8 ON ehr.person_id = obs8.person_id
")

ehr_df <- download_data(ehr_query) %>% distinct(person_id, .keep_all = TRUE)

# =============================================================================
# SECTION 2: Variable Categorization
# =============================================================================
# All missing/unknown values coded as 999

mother_df <- ehr_df %>%
  mutate(
    # ---- Race/Ethnicity ----
    # Hispanic ethnicity takes priority over race category
    raceethnicity = case_when(
      ethnicity == "Hispanic or Latino" ~ "Hispanic",
      race == "White" ~ "Non-Hispanic White",
      race == "Black or African American" ~ "Non-Hispanic Black",
      race == "Asian" ~ "Non-Hispanic Asian",
      TRUE ~ "Other Race or Ethnicity"
    ),

    # ---- Age ----
    # Reference date: end of exposure period (2023-10-01)
    age = as.integer(floor(difftime(as.Date("2023-10-01"), date_of_birth,
                                    units = "days") / 365.25)),
    age_group = cut(age, breaks = c(18, 45, 65, Inf),
                    labels = c("18-44", "45-64", ">=65"), right = FALSE),

    # ---- Sex (binary) ----
    # 0 = Male, 1 = Female, 999 = Missing/Other
    sex_cat = ifelse(sex == "Male", 0, ifelse(sex == "Female", 1, 999)),

    # ---- Race/Ethnicity (numeric) ----
    # 0 = NH-White, 1 = NH-Black, 2 = Hispanic, 3 = Other, 999 = Missing
    raceethnicity_cat = case_when(
      ethnicity == "Hispanic or Latino" ~ 2,
      race == "White" ~ 0,
      race == "Black or African American" ~ 1,
      TRUE ~ 3
    ),

    # ---- Income ----
    # 0 = Low (<$35K), 1 = Middle ($35-100K), 2 = High (>$100K)
    income = case_when(
      aname1 %in% c("Annual Income: less 10k", "Annual Income: 10k 25k",
                     "Annual Income: 25k 35k") ~ 0,
      aname1 %in% c("Annual Income: 35k 50k", "Annual Income: 50k 75k",
                     "Annual Income: 75k 100k") ~ 1,
      aname1 %in% c("Annual Income: 100k 150k", "Annual Income: 150k 200k",
                     "Annual Income: more 200k") ~ 2,
      TRUE ~ 999
    ),

    # ---- Education ----
    # 0 = High school or below, 1 = Some college/graduate, 2 = Advanced degree
    education = case_when(
      aname2 %in% c("Highest Grade: Never Attended", "Highest Grade: One Through Four",
                     "Highest Grade: Five Through Eight", "Highest Grade: Nine Through Eleven",
                     "Highest Grade: Twelve Or GED") ~ 0,
      aname2 %in% c("Highest Grade: College One to Three",
                     "Highest Grade: College Graduate") ~ 1,
      aname2 == "Highest Grade: Advanced Degree" ~ 2,
      TRUE ~ 999
    ),

    # ---- Insurance ----
    # 0 = Uninsured, 1 = Public (Medicare/Medicaid/Military/VA/Indian),
    # 2 = Private (Employer/Purchased/Other)
    insurance_category = case_when(
      ins1 == "Health Insurance: No" | ins2 == "Insurance Type Update: None" ~ 0,
      ins2 %in% c("Insurance Type Update: Medicare", "Insurance Type Update: Medicaid",
                   "Insurance Type Update: Military", "Insurance Type Update: VA",
                   "Insurance Type Update: Indian") ~ 1,
      ins2 %in% c("Insurance Type Update: Employer Or Union",
                   "Insurance Type Update: Purchased",
                   "Insurance Type Update: Other Health Plan") ~ 2,
      TRUE ~ 999
    ),

    # ---- Smoking ----
    # Two-step derivation: (1) 100-cigarette lifetime, (2) current frequency
    # 0 = Never (< 100 cigs lifetime)
    # 1 = Former (>= 100 cigs but currently not smoking)
    # 2 = Current (>= 100 cigs and currently smoking some/every day)
    cigs = ifelse(aname3 == "100 Cigs Lifetime: Yes", 1,
                  ifelse(aname3 == "100 Cigs Lifetime: No", 0, 999)),
    cigs_frequency = case_when(
      aname4 %in% c("Smoke Frequency: Some Days", "Smoke Frequency: Every Day") ~ 1,
      aname4 == "Smoke Frequency: Not At All" ~ 0,
      TRUE ~ 999
    ),
    smoking = case_when(
      cigs_frequency == 1 ~ 2,                    # Current smoker
      cigs == 0 ~ 0,                              # Never smoker
      cigs == 1 & cigs_frequency == 0 ~ 1,        # Former smoker
      TRUE ~ 999
    ),

    # ---- Alcohol (AUDIT-C) ----
    # Component 1: alcohol_freq (0-4 from "Never" to "4+ per week")
    # Component 2: avg_daily_drink (0-4 from "1-2" to "10+")
    # Component 3: heavy_drink_freq (0-4 from "Never" to "Daily")
    # AUDIT-C total = sum of 3 components (range 0-12)
    alcohol = ifelse(aname5 == "Alcohol Participant: Yes", 1,
                     ifelse(aname5 == "Alcohol Participant: No", 0, 999)),
    alcohol_freq = case_when(
      aname6 == "Drink Frequency Past Year: Never" ~ 0,
      aname6 == "Drink Frequency Past Year: Monthly Or Less" ~ 1,
      aname6 == "Drink Frequency Past Year: 2 to 4 Per Month" ~ 2,
      aname6 == "Drink Frequency Past Year: 2 to 3 Per Week" ~ 3,
      aname6 == "Drink Frequency Past Year: 4 or More Per Week" ~ 4,
      TRUE ~ 999
    ),
    avg_daily_drink = case_when(
      aname7 == "Average Daily Drink Count: 1 or 2" ~ 0,
      aname7 == "Average Daily Drink Count: 3 or 4" ~ 1,
      aname7 == "Average Daily Drink Count: 5 or 6" ~ 2,
      aname7 == "Average Daily Drink Count: 7 to 9" ~ 3,
      aname7 == "Average Daily Drink Count: 10 or More" ~ 4,
      TRUE ~ 999
    ),
    heavy_drink_freq = case_when(
      aname8 == "6 or More Drinks Occurrence: Never In Last Year" ~ 0,
      aname8 == "6 or More Drinks Occurrence: Less Than Monthly" ~ 1,
      aname8 == "6 or More Drinks Occurrence: Monthly" ~ 2,
      aname8 == "6 or More Drinks Occurrence: Weekly" ~ 3,
      aname8 == "6 or More Drinks Occurrence: Daily" ~ 4,
      TRUE ~ 999
    )
  ) %>%
  mutate(
    # ---- AUDIT-C Score and Classification ----
    # Non-drinker = 0; if any component missing = 999
    audit_c_score = case_when(
      alcohol == 0 ~ 0,
      alcohol == 999 | alcohol_freq == 999 |
        avg_daily_drink == 999 | heavy_drink_freq == 999 ~ 999,
      TRUE ~ alcohol_freq + avg_daily_drink + heavy_drink_freq
    ),
    # 0 = Low Risk (0-4), 1 = Increased Risk (5-7),
    # 2 = High Risk (8-10), 3 = Dependent (11-12)
    alcohol_category = case_when(
      audit_c_score == 999 ~ 999,
      audit_c_score <= 4 ~ 0,
      audit_c_score <= 7 ~ 1,
      audit_c_score <= 10 ~ 2,
      TRUE ~ 3
    )
  )

# =============================================================================
# SECTION 3: Disease Phenotype Flags (BigQuery)
# =============================================================================
# Uses OMOP Concept IDs from phenotype_concept_ids.R
# Method: Descendant concept traversal via cb_criteria hierarchy table
#   1. Anchor concepts identified by concept_id + rank1 filter
#   2. All descendant concepts found via path matching
#   3. Events queried from cb_search_all_events
# Returns: person_id, {disease}_start_date (earliest event), {disease} flag (0/1)

build_cte <- function(vec, cte_name, std_flag) {
  if (!length(vec)) return("")
  ids <- paste(vec, collapse = ", ")
  glue("{cte_name} AS (
    SELECT DISTINCT c.concept_id FROM `{cdr_ds}.cb_criteria` c
    JOIN (SELECT CAST(id AS STRING) AS id FROM `{cdr_ds}.cb_criteria`
          WHERE concept_id IN ({ids}) AND full_text LIKE '%_rank1]%') anchor
      ON c.path LIKE CONCAT('%.', anchor.id, '.%')
      OR c.path LIKE CONCAT('%.', anchor.id)
      OR c.path LIKE CONCAT(anchor.id, '.%')
      OR c.path = anchor.id
    WHERE c.is_standard = {std_flag} AND c.is_selectable = 1)")
}

get_person_ids <- function(dspec, d_name) {
  ctes <- c(build_cte(dspec$include_std, "inc_std_desc", 1),
            build_cte(dspec$include_src, "inc_src_desc", 0),
            build_cte(dspec$exclude_std, "exc_std_desc", 1),
            build_cte(dspec$exclude_src, "exc_src_desc", 0))
  ctes <- ctes[ctes != ""]

  has_std <- length(dspec$include_std) > 0
  has_src <- length(dspec$include_src) > 0
  inc_filter <- paste(c(
    if (has_std) "(is_standard = 1 AND concept_id IN (SELECT concept_id FROM inc_std_desc))",
    if (has_src) "(is_standard = 0 AND concept_id IN (SELECT concept_id FROM inc_src_desc))"
  ), collapse = " OR ")

  exc_filter <- ""
  exc_parts <- c(
    if (length(dspec$exclude_std)) "(is_standard = 1 AND concept_id IN (SELECT concept_id FROM exc_std_desc))",
    if (length(dspec$exclude_src)) "(is_standard = 0 AND concept_id IN (SELECT concept_id FROM exc_src_desc))"
  )
  if (length(exc_parts)) {
    exc_filter <- glue("WHERE person_id NOT IN (
      SELECT DISTINCT person_id FROM `{cdr_ds}.cb_search_all_events`
      WHERE {paste(exc_parts, collapse = ' OR ')})")
  }

  sql <- glue("WITH {paste(ctes, collapse = ',\n')},
    qualified AS (SELECT person_id, entry_date FROM `{cdr_ds}.cb_search_all_events`
                  WHERE {inc_filter})
    SELECT person_id, MIN(entry_date) AS {d_name}_start_date
    FROM qualified {exc_filter} GROUP BY person_id")

  download_data(sql) %>% mutate(!!d_name := 1L)
}

# Variant for recurrent seizure phenotype: requires >= min_events distinct dates
get_person_ids_recurrent <- function(dspec, d_name, min_events = 2) {
  # Same CTE construction as get_person_ids, plus:
  # HAVING COUNT(DISTINCT entry_date) >= min_events
  # Returns: person_id, {d_name}_start_date, {d_name}_event_count
}

# =============================================================================
# SECTION 4: Merge Flags and Build Derived Outcomes
# =============================================================================
# Disease sets defined using concept IDs from phenotype_concept_ids.R:
#   epilepsy_or_seizure, epilepsy_dx, seizure_only (recurrent),
#   epilepsy_g40, hypoglycemia, mci, adrd, stroke, t2d_excluded_t1d,
#   + all Charlson CCI comorbidities
#
# flag_df <- imap(disease_sets, get_person_ids) |> reduce(full_join, by = "person_id")
#
# Derived outcome variables:
#   epilepsy_refined = epilepsy_dx | seizure_recurrent (>= 2 distinct dates)
#   hypoglycemic_seizure = same-day co-occurrence of epilepsy_or_seizure + hypoglycemia
#   epilepsy_or_seizure_without_hypoglycemic_seizure = outcome 1 minus hypoglycemic seizures

# =============================================================================
# SECTION 5: EHR Date Ranges and Death Data
# =============================================================================
# Extract earliest and latest dates from each EHR domain table:
#   measurement, condition_occurrence, drug_exposure,
#   procedure_occurrence, observation, visit_occurrence
#
# Survey dates: first_survey_date, last_survey_date
# Death data from aou_death table: death_date, cause_of_death
# Merged into: merged_df (mother_df + EHR dates + death info)

# =============================================================================
# SECTION 6: BMI and HbA1c Panel Data (Longitudinal)
# =============================================================================
# BMI panel: all measurements with value 10-100 (concept_id = 3038553)
# HbA1c panel: all measurements with value 2-20 (concept_id = 3004410)
# Used for: baseline covariate extraction and longitudinal trajectory analysis (LMM)
# Output: bmi_panel (person_id, measurement_date, bmi_value)
#         a1c_panel (person_id, measurement_date, a1c_value)
