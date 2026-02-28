# =============================================================================
# rev00005_TMLE_Bootstrap_outcome4.R
# TMLE Bootstrap Analysis for Outcome 4 (Epilepsy/Seizure excl Hypoglycemic)
# =============================================================================
#
# Modified from rev00005_TMLE_Bootstrap.R with the following changes:
#   1. analysis_types prefix → outcome4_epilepsy_excl_hypo_*
#   2. use_superlearner = TRUE (main TMLE + bootstrap)
#   3. n_bootstrap = 1000
#
# INPUT:
#   - outcome4_epilepsy_excl_hypo_all_ages_ipwt_full_results.rds
#   - outcome4_epilepsy_excl_hypo_late_onset_ipwt_full_results.rds
#
# OUTPUT:
#   - outcome4_epilepsy_excl_hypo_*_tmle_full_results.rds
#   - outcome4_epilepsy_excl_hypo_*_tmle_comprehensive_results_*.csv
#   - Bootstrap RDS/CSV/PNG files
# =============================================================================

# >>> CELL 01: Install & Load Libraries <<<
# install.packages("MatchIt")
# install.packages("survminer")
# install.packages("mice")
# install.packages("tableone")

# Load required packages silently
suppressPackageStartupMessages({
  library(tidyverse)
  library(bigrquery)
  library(dplyr)
  library(tidyr)
  library(glue)
  library(MatchIt)
  library(survival)
  # library(survminer)  # Not needed for TMLE bootstrap
  library(lubridate)
  library(ggplot2)
  library(survey)
  library(mice)
  library(tableone)
  library(readr)
})

# >>> CELL 02: Configuration <<<
# ---------- CONFIG ----------
get_standard_config <- function() {
  list(
    exposure_window = list(start = "2018-01-01", end = "2023-10-01"), 
    data_cut_date = as.Date("2023-10-01"),
    drug_classes = list(
      SEMAGLUTIDE = c("semaglutide"),
      OTHER_GLPA = c("exenatide", "liraglutide", "albiglutide",
                     "dulaglutide", "lixisenatide"), # "tirzepatide" excluded
      SGLT2 = c("canagliflozin", "empagliflozin", "dapagliflozin", "ertugliflozin"),
      OtherGLD = c("glimepiride", "glipizide", "glyburide",
                   "alogliptin", "linagliptin", "sitagliptin", "saxagliptin", "pioglitazone", "rosiglitazone")
    ),
    comparisons = list(
      list(name = "SEMAGLUTIDE vs OTHER_GLPA", exposure = "SEMAGLUTIDE", comparator = "OTHER_GLPA"),
      list(name = "SEMAGLUTIDE vs SGLT2", exposure = "SEMAGLUTIDE", comparator = "SGLT2"),
      list(name = "SEMAGLUTIDE vs OtherGLD", exposure = "SEMAGLUTIDE", comparator = "OtherGLD"),
      list(name = "OTHER_GLPA vs SGLT2", exposure = "OTHER_GLPA", comparator = "SGLT2"),
      list(name = "OTHER_GLPA vs OtherGLD", exposure = "OTHER_GLPA", comparator = "OtherGLD"),
      list(name = "SGLT2 vs OtherGLD", exposure = "SGLT2", comparator = "OtherGLD")
    ),
    outcomes = list(
      list(var = "epilepsy_or_seizure_start_date", label = "Epilepsy/Seizure",
           late_onset = FALSE, early_onset = FALSE),
      list(var = "epilepsy_refined_start_date", label = "Epilepsy Refined",
           late_onset = FALSE, early_onset = FALSE),
      list(var = "epilepsy_g40_start_date", label = "Epilepsy G40",
           late_onset = FALSE, early_onset = FALSE),
      list(var = "epilepsy_or_seizure_without_hypoglycemic_seizure_start_date", label = "Epilepsy/Seizure excl Hypoglycemic",
           late_onset = FALSE, early_onset = FALSE)
    ),
    all_ps_vars = c(
      # Demographics and Vitals
      "age", "sex_cat", "raceethnicity_cat", "income", "education", "insurance_category",
      "smoking", "alcohol_category",
      "baseline_bmi_category", "baseline_hba1c", "index_year_grouped",
      # Medications
      "Anticoagulant", "Antiplatelet", "BB", "Biguanide", "CCB", "Diuretic",
      "Ezetimibe", "MRA", "OtherHTN", "RAAS", "Statin", "TZD", "Insulin",
      "SEMAGLUTIDE", "OTHER_GLPA", "SGLT2i", "SU", "DPP4i",
      # Comorbidities
      "myocardial_infarction", "congestive_heart_failure", "peripheral_vascular_disease",
      "cerebrovascular_disease", "chronic_pulmonary_disease", "dementia",
      "rheumatic_disease", "peptic_ulcer_disease",
      "hemiplegia_or_paraplegia", "hiv_infection",
      "hypoglycemia", "hyperglycemic_emergency",
      # New consolidated variables
      "renal_disease_severity",
      "liver_disease_severity",
      "diabetes_with_ophthalmic_complications", "diabetes_with_neurological_complications", "malignancy_status"
    ),
    categorical_vars = c(
      "sex_cat", "raceethnicity_cat", "income", "education", "insurance_category",
      "smoking", "alcohol_category", "index_year", "index_year_grouped", "baseline_bmi_category",
      "Anticoagulant", "Antiplatelet", "BB", "Biguanide", "CCB", "Diuretic",
      "Ezetimibe", "MRA", "OtherHTN", "RAAS", "Statin", "TZD", "Insulin",
      "SEMAGLUTIDE", "OTHER_GLPA", "SGLT2i", "SU", "DPP4i",
      "myocardial_infarction", "congestive_heart_failure", "peripheral_vascular_disease",
      "cerebrovascular_disease", "dementia", "chronic_pulmonary_disease",
      "rheumatic_disease", "peptic_ulcer_disease",
      "hemiplegia_or_paraplegia", "hiv_infection",
      "hypoglycemia", "hyperglycemic_emergency",
      "liver_disease_severity","renal_disease_severity",
      "diabetes_with_ophthalmic_complications", "diabetes_with_neurological_complications", "malignancy_status"),
    # Exclusion toggles (set to TRUE/FALSE to control behavior)
    exclude_renal_severe = FALSE,
    exclude_malignancy_metastatic = FALSE,
    log_exclusions = FALSE
  )}

# >>> CELL 03: MICE Imputation Function <<<
# ============================================================================
# MICE IMPUTATION FUNCTION
# ============================================================================
perform_mice_imputation <- function(df, m = 20, seed = 123) {
  set.seed(seed)
  
  # Check if imputation is needed
  missing_vars <- c("baseline_bmi", "baseline_hba1c")
  available_vars <- intersect(missing_vars, names(df))
  
  if (length(available_vars) == 0) return(df)
  
  missing_count <- sapply(df[available_vars], function(x) sum(is.na(x)))
  if (all(missing_count == 0)) return(df)
  
  cat("MICE imputation for:", paste(names(missing_count)[missing_count > 0], collapse = ", "), "\n")
  
  # Perform imputation using only relevant variables to avoid convergence issues
  # Multiple imputation - use relevant predictors including DM complications
  imp_vars <- intersect(c(
    "baseline_bmi", "baseline_hba1c",  # Variables to impute
    # Key predictors for BMI/HbA1c
    "age", "sex_cat", "treatment", "raceethnicity_cat", "smoking", "alcohol_category",
    "income", "education", "insurance_category", 
    # Diabetes-related medications
    "Biguanide", "Insulin", "SEMAGLUTIDE", "OTHER_GLPA", "TZD", "SU", "DPP4i", "SGLT2i",
    # HTN medications - relevant for BMI/HbA1c
    "RAAS", "Diuretic", "MRA", "BB", "CCB", "OtherHTN",
    # Key comorbidities
    "myocardial_infarction", "congestive_heart_failure", "chronic_pulmonary_disease",
    "liver_disease_severity", "renal_disease_severity",
    # Diabetes complications - very relevant for BMI/HbA1c
    "diabetes_with_renal_complications",  
    "diabetes_with_ophthalmic_complications",
    "diabetes_with_neurological_complications",
    # Vascular complications  
    "peripheral_vascular_disease",
    "cerebrovascular_disease",
    # Diabetes-related events
    "hypoglycemia", "hyperglycemic_emergency"
  ), names(df))
  
  tryCatch({
    imp <- mice(df[imp_vars], m = m, seed = seed, printFlag = FALSE)
    completed <- complete(imp, 1)
    # Log imputation counts for key variables
    pre_na_bmi    <- if ("baseline_bmi"    %in% names(df)) sum(is.na(df$baseline_bmi)) else NA_integer_
    pre_na_hba1c  <- if ("baseline_hba1c"  %in% names(df)) sum(is.na(df$baseline_hba1c)) else NA_integer_
    imputed_bmi   <- if ("baseline_bmi"    %in% names(df) && "baseline_bmi"   %in% names(completed)) sum(is.na(df$baseline_bmi)   & !is.na(completed$baseline_bmi)) else NA_integer_
    imputed_hba1c <- if ("baseline_hba1c"  %in% names(df) && "baseline_hba1c" %in% names(completed)) sum(is.na(df$baseline_hba1c) & !is.na(completed$baseline_hba1c)) else NA_integer_
    post_na_bmi   <- if ("baseline_bmi"    %in% names(completed)) sum(is.na(completed$baseline_bmi)) else NA_integer_
    post_na_hba1c <- if ("baseline_hba1c"  %in% names(completed)) sum(is.na(completed$baseline_hba1c)) else NA_integer_
    log_output("PROGRESS", "MICE imputation counts",
               sprintf("BMI pre_na=%s, imputed=%s, post_na=%s; HbA1c pre_na=%s, imputed=%s, post_na=%s",
                       as.character(pre_na_bmi), as.character(imputed_bmi), as.character(post_na_bmi),
                       as.character(pre_na_hba1c), as.character(imputed_hba1c), as.character(post_na_hba1c)),
               "mice")
    
    # Replace imputed values in original dataframe
    for (var in available_vars) {
      if (var %in% names(completed)) {
        df[[var]] <- completed[[var]]
      }
    }
    return(df)
  }, error = function(e) {
    cat("MICE failed, using mean imputation as fallback\n")
    if ("baseline_bmi" %in% names(df)) {
      df$baseline_bmi[is.na(df$baseline_bmi)] <- mean(df$baseline_bmi, na.rm = TRUE)
    }
    if ("baseline_hba1c" %in% names(df)) {
      df$baseline_hba1c[is.na(df$baseline_hba1c)] <- mean(df$baseline_hba1c, na.rm = TRUE)
    }
    return(df)
  })
}

# >>> CELL 04: Logging Functions <<<
# ---------- SIMPLE LOGGING ----------
init_output_tracker <- function() {
  options(semaglutide_output_tracker = data.frame(
    timestamp = character(),
    script = character(),
    type = character(),
    message = character(),
    value = character(),
    stringsAsFactors = FALSE
  ))
}
log_output <- function(type, message, value = "", script = "") {
  if (is.null(getOption("semaglutide_output_tracker"))) init_output_tracker()
  tracker <- getOption("semaglutide_output_tracker")
  new_row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    script = script, type = type, message = message, value = as.character(value),
    stringsAsFactors = FALSE
  )
  options(semaglutide_output_tracker = rbind(tracker, new_row))
  if (type %in% c("PROGRESS", "ERROR", "WARNING")) {
    cat(sprintf("[%s] %s: %s\n", type, message, value))
    # Force flush for Jupyter/IRkernel
    try(flush(stdout()), silent = TRUE)
    try(flush(stderr()), silent = TRUE)
    if (interactive()) flush.console()
  }
}
save_output_log <- function(filename = NULL) {
  if (is.null(filename)) filename <- sprintf("output_log_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S"))
  tracker <- getOption("semaglutide_output_tracker")
  if (!is.null(tracker) && nrow(tracker) > 0) {
    write.csv(tracker, filename, row.names = FALSE)
    log_output("PROGRESS", "Output log saved", filename)
    return(filename)
  }
  NULL
}
# >>> CELL 05: Data Load & Exposure Index Builder <<<
# ---------- DATA LOAD ----------
load_base_data <- function(workspace_bucket = NULL) {
  if (is.null(workspace_bucket)) workspace_bucket <- Sys.getenv('WORKSPACE_BUCKET')
  log_output("PROGRESS", "Loading base data", "", "load_base_data")
  # merged_df, bmi_panel, a1c_panel expected in bucket/data/
  system(paste0("gsutil cp ", workspace_bucket, "/data/merged_df.csv ."), intern = TRUE)
  merged_df <- read_csv("merged_df.csv", show_col_types = FALSE)
  # normalize person_id and any date columns from CSVs to consistent types
  if ("person_id" %in% names(merged_df)) merged_df <- merged_df %>% mutate(person_id = as.character(person_id))
  date_cols_csv <- names(merged_df)[grepl("date|Date", names(merged_df), ignore.case = TRUE) | names(merged_df) %in% c("index_date")]
  for (dc in date_cols_csv) {
    merged_df[[dc]] <- tryCatch(as.Date(merged_df[[dc]]), error = function(e) merged_df[[dc]])
  }
  system(paste0("gsutil cp ", workspace_bucket, "/data/bmi_panel.csv ."), intern = TRUE)
  bmi_panel <- read_csv("bmi_panel.csv", show_col_types = FALSE) %>%
    mutate(weight_date = as.Date(weight_date)) %>%
    select(person_id, weight_date, bmi) %>%
    distinct(person_id, weight_date, .keep_all = TRUE)
  if ("person_id" %in% names(bmi_panel)) bmi_panel <- bmi_panel %>% mutate(person_id = as.character(person_id))
  system(paste0("gsutil cp ", workspace_bucket, "/data/a1c_panel.csv ."), intern = TRUE)
  a1c_panel <- read_csv("a1c_panel.csv", show_col_types = FALSE)
  if ("person_id" %in% names(a1c_panel)) a1c_panel <- a1c_panel %>% mutate(person_id = as.character(person_id))
  # convert some flags to factor if present
  cols_to_factor <- intersect(c("epilepsy_or_seizure", "adrd", "mci", "stroke", "t2d_excluded_t1d"), names(merged_df))
  merged_df[cols_to_factor] <- lapply(merged_df[cols_to_factor], function(x) factor(x, levels = c(0,1), labels = c("0","1")))
  log_output("INFO", "Data loaded", sprintf("%d rows, %d columns", nrow(merged_df), ncol(merged_df)))
  list(merged_df = merged_df, bmi_panel = bmi_panel, hba1c_panel = a1c_panel)
}
# ---------- BIGQUERY DOWNLOAD ----------
download_data <- function(sql, google_project = NULL, verbose = FALSE) {
  if (is.null(google_project)) google_project <- Sys.getenv("GOOGLE_PROJECT")
  if (verbose) log_output("DEBUG", "Executing BigQuery", substr(sql, 1, 120))
  res <- bq_table_download(bq_project_query(google_project, sql))
  # normalize person_id to character if present to avoid integer64/double joins
  if ("person_id" %in% names(res)) {
    res <- res %>% mutate(person_id = as.character(person_id))
  }
  # convert date-like columns (common names) to Date
  date_like_cols <- names(res)[grepl("date$", names(res), ignore.case = TRUE) | names(res) %in% c("index_date","exposure_date","comparator_date")]
  for (dc in date_like_cols) {
    res[[dc]] <- tryCatch(as.Date(res[[dc]]), error = function(e) res[[dc]])
  }
  res
}
# ---------- EXPOSURE INDEX BUILDER (MODIFIED) ----------
build_exposure_idx <- function(concepts_df, idx_name, exposure_window = NULL, dataset = NULL, washout_concepts_df = NULL) {
  if (is.null(dataset)) dataset <- Sys.getenv("WORKSPACE_CDR")
  if (nrow(concepts_df) == 0) return(data.frame(person_id = character(0), index_date = as.Date(character(0))))
  
  # Concept IDs for the primary drug class of interest
  id_string <- paste(concepts_df$concept_id, collapse = ", ")
  
  # --- CHANGE HERE ---
  # Combine primary and washout (comparator) concepts for the prior use check
  all_washout_concepts <- concepts_df
  if (!is.null(washout_concepts_df) && nrow(washout_concepts_df) > 0) {
    all_washout_concepts <- rbind(all_washout_concepts, washout_concepts_df)
  }
  washout_id_string <- paste(unique(all_washout_concepts$concept_id), collapse = ", ")
  # --- END CHANGE ---
  sql <- glue("
    WITH first_use AS (
      SELECT person_id, MIN(drug_exposure_start_date) AS index_date
      FROM `{dataset}.drug_exposure`
      WHERE drug_concept_id IN ({id_string}) -- This remains unchanged to define the index date correctly
        AND drug_exposure_start_date BETWEEN '{exposure_window$start}' AND '{exposure_window$end}'
      GROUP BY person_id
    ),
    prior_use AS (
      SELECT f.person_id
      FROM first_use f
      JOIN `{dataset}.drug_exposure` d
        ON d.person_id = f.person_id
       -- --- CHANGE HERE ---
       -- Now checks for prior use of BOTH the index drug class AND the comparator classes
       AND d.drug_concept_id IN ({washout_id_string})
       -- --- END CHANGE ---
       AND d.drug_exposure_start_date BETWEEN DATE_SUB(f.index_date, INTERVAL 183 DAY)
                                          AND DATE_SUB(f.index_date, INTERVAL 1 DAY)
      GROUP BY f.person_id
    )
    SELECT f.person_id, f.index_date
    FROM first_use f
    LEFT JOIN prior_use p ON f.person_id = p.person_id
    WHERE p.person_id IS NULL
  ")
  
  res <- download_data(sql)
  # keep a standard column name 'index_date' and normalize types
  if ("person_id" %in% names(res)) res <- res %>% mutate(person_id = as.character(person_id))
  if ("index_date" %in% names(res)) res <- res %>% mutate(index_date = as.Date(index_date))
  res
}
# >>> CELL 06: Medication & Comorbidity Flagging <<<
# ---------- SMALL HELPERS FOR CTEs ----------
build_cte <- function(vec, cte_name, std_flag, dataset) {
  if (!length(vec)) return("")
  ids <- glue::glue_collapse(vec, sep = ", ")
  glue::glue("
    {cte_name} AS (
      SELECT DISTINCT c.concept_id
      FROM `{dataset}.cb_criteria` c
      JOIN (SELECT CAST(id AS STRING) AS id
            FROM `{dataset}.cb_criteria`
            WHERE concept_id IN ({ids})
              AND full_text LIKE '%_rank1]%') anchor
        ON c.path LIKE CONCAT('%.', anchor.id, '.%')
        OR c.path LIKE CONCAT('%.', anchor.id)
        OR c.path LIKE CONCAT(anchor.id, '.%')
        OR c.path = anchor.id
      WHERE c.is_standard = {std_flag}
        AND c.is_selectable = 1 )
  ")
}

# ---------- MEDICATION FLAGGING (kept) ----------
add_medications <- function(cohort_df, dataset = Sys.getenv("WORKSPACE_CDR")) {
  if (nrow(cohort_df) == 0) return(cohort_df)
  # ensure cohort ids and dates are consistent (character and Date) before building SQL
  if ("person_id" %in% names(cohort_df)) cohort_df <- cohort_df %>% mutate(person_id = as.character(person_id))
  if ("index_date" %in% names(cohort_df)) cohort_df <- cohort_df %>% mutate(index_date = as.Date(index_date))
  medication_list <- list(
    SEMAGLUTIDE = c("semaglutide"),
    OTHER_GLPA = c("exenatide", "liraglutide", "albiglutide",
                   "dulaglutide", "lixisenatide"),
    Biguanide = c("metformin"),
    TZD = c("pioglitazone", "rosiglitazone"),
    Insulin = c("insulin glargine", "insulin aspart", "insulin lispro",
                "insulin detemir", "insulin degludec", "insulin glulisine",
                "insulin regular", "insulin isophane", "insulin human", "insulin"),
    SGLT2i = c("canagliflozin", "empagliflozin", "dapagliflozin", "ertugliflozin"),
    DPP4i = c("alogliptin","linagliptin", "sitagliptin", "saxagliptin"),
    SU = c("glimepiride", "glipizide", "glyburide"),
    Anticoagulant = c("warfarin", "dabigatran", "rivaroxaban", "apixaban", "edoxaban", "betrixaban"),
    Antiplatelet = c("anagrelide", "cilostazol", "clopidogrel", "dipyridamole", "prasugrel",
                     "ticagrelor", "ticlopidine", "vorapaxar", "cangrelor"),
    Statin = c("atorvastatin", "rosuvastatin", "simvastatin", "pravastatin",
               "lovastatin", "fluvastatin", "pitavastatin"),
    Ezetimibe = c("ezetimibe"),
    PCSK9 = c("alirocumab", "evolocumab", "inclisiran"),
    RAAS = c("benazepril", "captopril", "enalapril", "fosinopril", "lisinopril",
             "moexipril", "perindopril", "quinapril", "ramipril", "trandolapril",
             "azilsartan", "candesartan", "eprosartan", "irbesartan", "losartan",
             "olmesartan", "telmisartan", "valsartan", "aliskiren"),
    SacubV = c("sacubitril/valsartan"),
    Diuretic = c("bendroflumethiazide", "chlorothiazide", "chlorthalidone",
                 "hydrochlorothiazide", "indapamide", "methyclothiazide",
                 "metolazone", "bumetanide", "ethacrynate sodium", "ethacrynic acid",
                 "furosemide", "torsemide", "amiloride", "triamterene"),
    MRA = c("eplerenone", "spironolactone"),
    BB = c("atenolol", "betaxolol", "bisoprolol", "metoprolol", "nebivolol",
           "nadolol", "propranolol", "acebutolol", "pindolol", "timolol",
           "carvedilol", "labetalol", "esmolol", "sotalol"),
    CCB = c("amlodipine", "felodipine", "isradipine", "nicardipine", "nifedipine",
            "nisoldipine", "clevidipine", "nimodipine", "diltiazem", "verapamil"),
    OtherHTN = c("doxazosin", "prazosin", "terazosin", "alfuzosin", "clonidine",
                 "methyldopa", "guanfacine", "hydralazine", "minoxidil",
                 "guanethidine", "tolazoline", "sodium nitroprusside",
                 "phenoxybenzamine", "phentolamine", "fenoldopam")
  )
  # remove rows with missing person_id or index_date; format dates for safe SQL injection
  valid_rows <- cohort_df %>% filter(!is.na(person_id) & !is.na(index_date))
  if (nrow(valid_rows) == 0) {
    # nothing to query — return empty flags frame (to be joined upstream)
    cohort_cte_sql <- ""
  } else {
    # ensure index_date is Date and format as YYYY-MM-DD
    valid_rows <- valid_rows %>% mutate(index_date = as.Date(index_date))
    cohort_rows <- purrr::map2_chr(valid_rows$person_id, valid_rows$index_date,
                                   ~ glue::glue("STRUCT({.x} AS person_id, DATE('{format(.y, '%Y-%m-%d')}') AS index_date)"))
    cohort_cte_sql <- paste(cohort_rows, collapse = ", ")
  }
  case_when_clauses <- purrr::map2_chr(names(medication_list), medication_list, ~{
    like_sql <- paste0("LOWER(c.concept_name) LIKE '%", .y, "%'", collapse = " OR ")
    glue::glue("WHEN {like_sql} THEN '{.x}'")
  })
  case_when_sql <- paste(case_when_clauses, collapse = "\n ")
  pivot_clauses <- paste0("MAX(IF(me.med_class = '", names(medication_list), "', 1, 0)) AS ", names(medication_list), collapse = ",\n ")
  sql_medications <- glue::glue("
    WITH
    cohort AS (
      SELECT * FROM UNNEST(ARRAY<STRUCT<person_id INT64, index_date DATE>>[{cohort_cte_sql}])
    ),
    MedicationExposures AS (
      SELECT DISTINCT
        t.person_id,
        CASE
          {case_when_sql}
          ELSE NULL
        END AS med_class
      FROM `{dataset}.drug_exposure` de
      JOIN `{dataset}.concept` c ON de.drug_concept_id = c.concept_id
      JOIN cohort t ON de.person_id = t.person_id
      WHERE de.drug_exposure_start_date BETWEEN DATE_SUB(t.index_date, INTERVAL 365 DAY) AND t.index_date
    )
    SELECT
      t.person_id,
      t.index_date,
      {pivot_clauses}
    FROM cohort t
    LEFT JOIN MedicationExposures me ON t.person_id = me.person_id
    GROUP BY t.person_id, t.index_date
  ")
  medication_flags <- download_data(sql_medications)
  # ensure returned flags use character person_id and Date index_date
  if ("person_id" %in% names(medication_flags)) medication_flags <- medication_flags %>% mutate(person_id = as.character(person_id))
  if ("index_date" %in% names(medication_flags)) medication_flags <- medication_flags %>% mutate(index_date = as.Date(index_date))
  cohort_df <- cohort_df %>% left_join(medication_flags, by = c("person_id", "index_date"))
  for (med_name in names(medication_list)) if (med_name %in% names(cohort_df)) cohort_df[[med_name]][is.na(cohort_df[[med_name]])] <- 0
  cohort_df
}

# ---------- COMORBIDITY FLAGGING (kept) ----------
add_comorbidities <- function(cohort_df, dataset = Sys.getenv("WORKSPACE_CDR")) {
  if (nrow(cohort_df) == 0) return(cohort_df)
  # ensure cohort ids and dates are consistent before building SQL
  if ("person_id" %in% names(cohort_df)) cohort_df <- cohort_df %>% mutate(person_id = as.character(person_id))
  if ("index_date" %in% names(cohort_df)) cohort_df <- cohort_df %>% mutate(index_date = as.Date(index_date))
  comorbidity_sets <- list(
    myocardial_infarction = list(
      include_std = integer(0),
      include_src = c(44832372, 1569126, 1569130, 35207702, 44820864),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    congestive_heart_failure = list(
      include_std = integer(0),
      include_src = c(44833557, 44819695, 44823110, 44824235, 35207762, 44819693, 35207755, 44824250, 35207673, 44820856, 44819710, 44825438, 35207763, 1569178, 44821955, 35207765, 35207674, 44831230, 35207760, 44819696, 44823116, 35207761, 44819692, 35207704, 44835939, 35207764, 44823108, 35207669),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    peripheral_vascular_disease = list(
      include_std = integer(0),
      include_src = c(44835640, 35208277, 35207882, 44826657, 44834745, 35207869, 35208278, 44834813, 35225419, 44820429, 44830159, 44825446, 44819723, 35207881, 35207860, 44822297, 1569271, 35207883, 1569326, 44821962, 35208276, 1569321, 44826646, 44832392, 44835959, 44834746, 35207859, 1576287),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    cerebrovascular_disease = list(
      include_std = integer(0),
      include_src = c(1568361, 1569184, 1569190, 1569193, 1568726, 1569191, 44831252, 1568360, 44835952, 1569225, 44820872, 44820873, 44820875, 1568725, 1569227, 44835946, 1568727, 1569218, 44824253, 44830001, 44832388, 44835947, 1569221),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    dementia = list(
      include_std = integer(0),
      include_src = c(44831122, 44831078, 1568088, 44820749, 1568087, 44821814, 35207116, 35207328, 1568293, 35207360, 45553736, 44831079, 35211390, 45890911, 35207114, 44820073, 35207121, 44827645, 35207361, 44820708, 44835772, 45586320, 1568295, 35207511, 44833397, 44826538, 35207115, 35207118, 44821813, 44826537),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    chronic_pulmonary_disease = list(
      include_std = integer(0),
      include_src = c(44829014, 35208065, 44820888, 44827823, 1569496, 35208036, 1569495, 44821987, 44834771, 44829011, 1569492, 1569485, 44820887, 44827824, 1569493, 35208056, 35208063, 1569486, 44826682, 44823145, 1569487, 35208037, 35208013, 44835986, 1569488, 35208026, 44829013, 35208017, 44835981, 1569494, 44820892, 44827825, 44826681, 35208027),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    rheumatic_disease = list(
      include_std = integer(0),
      include_src = c(44831518, 35208841, 44831487, 35208832, 44836170, 1570619, 44832613, 44831474, 1570614, 35208834, 44837334, 1570612, 44828021, 44833584, 44834963, 1569966, 1570039, 44821123, 44819941, 35208820),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    peptic_ulcer_disease = list(
      include_std = integer(0),
      include_src = c(44830146, 1569564, 1569565, 44825505, 1569563, 44824313, 44833634, 1569562),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    mild_liver_disease = list(
      include_std = integer(0),
      include_src = c(35208332, 44835626, 44821690, 44833242, 35208363, 35208359, 35208368, 35208338, 35208335, 44829751, 44829754, 44832480, 44824336, 35208367, 44831320, 44826726, 1569680, 1569670, 1569671, 44830943, 44829749, 35208331, 44819379, 1567376, 44821549, 35208330, 1569675, 35208362, 44825528, 35208361, 44819803, 1569681, 35225408),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    renal_mild_or_moderate = list(
      include_std = integer(0),
      include_src = c(35209277, 44836035, 44819695, 1571474, 44832368, 1571472, 44821546, 35225404, 44837191, 35209275, 44830173, 44835923, 44832369, 35207673, 44827888, 35209279, 35207672, 44835922, 44819694, 44825538, 44819696, 44823191, 44820970, 35209276, 45543164, 44819692, 35209274, 44835924),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    hemiplegia_or_paraplegia = list(
      include_std = integer(0),
      include_src = c(44823020, 35207306, 1568408, 1568415, 44832266, 35207479, 35207480, 1568412, 35207481, 44832275, 44828872, 35207319),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    any_malignancy = list(
      include_std = integer(0),
      include_src = c(44825213, 1567501, 1567466, 1567641, 1567463, 1567492, 1567651, 44827551, 1567471, 44833301, 44822887, 1567478, 44826420, 1567462, 44828729, 44827586, 1567699, 1567461, 44824032, 44832117, 44834484, 1567470, 1567569, 35206056, 1567534, 1567674, 44830971, 44819425, 44829811, 1567472, 44834492, 44819434, 44833294, 1567537, 44829849, 44832129, 44822885, 1567476, 1567689, 44820621, 1567530, 44834486, 44829790, 44828731, 35206153, 1567493, 35206140, 44829795, 44831017, 1567573, 44835670, 1567482, 44835663, 44827582, 44832131, 44827557, 44825199, 44819422, 1567469, 35206253, 1567474, 1567633, 1567567, 44829823, 1567659, 1567668, 35206101, 44836833, 1567479, 44828781, 1567468, 35206266, 1567485, 1567502, 1567473, 44829799, 44834480, 1567705, 1567578, 44821758, 44829803, 1567715, 44825214, 44821735, 1567574, 44826411, 44826385, 1567480, 1567484, 1567568, 44824020, 44827567, 44829813, 35206185, 44830966, 35206080, 1567464, 1567465, 1567483, 1567529, 45561836, 44829815, 44825232, 1567680, 35206141, 1567528, 1567615, 1567711, 44822869, 44835671, 44833300, 44824034, 1567494, 44833283, 44820606, 1567475, 44833324, 1567481, 44827548, 44835661, 35206242, 1567566, 44826409, 35206183, 1567486, 1567565, 1567533, 1567477, 35206260, 44819441, 1567675),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    moderate_or_severe_liver_disease = list(
      include_std = integer(0),
      include_src = c(44834823, 44822040, 35207901, 1569401, 1569679, 1569678, 44837179, 1569674, 44828998, 44825455, 44825529, 35208366, 44832406, 1569672, 35208364, 35208365),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    renal_severe = list(
      include_std = integer(0),
      include_src = c(44829649, 44835497, 44824342, 35209280, 44824235, 44831231, 44833130, 45548653, 44828971, 44831232, 44819693, 44825426, 44835496, 44820856, 35209291, 44821578, 44829650, 44833560, 45596188, 35207674, 35209278, 35225436, 44829062, 1576113, 44834280, 44837192, 44821950, 44831947, 35207671, 44827889),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    hiv_infection = list(
      include_std = integer(0),
      include_src = c(4829737, 35205776),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    metastatic_solid_tumor = list(
      include_std = integer(0),
      include_src = c(1567618, 44836847, 44836850, 1567623, 45557018, 44819442, 44819443, 1567619, 35206334),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    diabetes_with_renal_complications = list(
      include_std = integer(0),
      include_src = c(44824074, 1567958, 1567975, 1567909, 1567942, 1567926),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    diabetes_with_ophthalmic_complications = list(
      include_std = integer(0),
      include_src = c(1567976, 1567910, 1567959, 1567927, 1567943, 44828794),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    diabetes_with_neurological_complications = list(
      include_std = integer(0),
      include_src = c(1567933, 1567916, 1567949, 1567982, 44827615, 1567965),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    hypoglycemia = list(
      include_std = integer(0),
      include_src = c(1567955, 1567939, 1567922, 35206889, 44831048, 44820686, 1567988, 35206888, 1567971, 35206887),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    hyperglycemic_emergency = list(
      include_std = integer(0),
      include_src = c(44831044, 1567973, 1567924, 44826459, 1567907, 1567957),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    t1d = list(
      include_std = integer(0),
      include_src = c(44820682, 44819504, 44832192, 44822934, 44829881, 44833368, 44821787, 44824071, 44832190, 44819501, 44819502, 44825264, 44836918, 44822936, 44822935, 1567940, 44832191, 44820683, 44834549, 44820684, 44831046),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    t2d = list(
      include_std = integer(0),
      include_src = c(44832194, 44828795, 44832193, 44836915, 44826461, 1567956, 44833366, 44829879, 44829882, 44824073, 44824072, 44831045, 44831047, 44836914, 44833367, 44827616, 44827617, 44836916, 44829878, 44826460, 44819500),
      exclude_std = integer(0),
      exclude_src = integer(0)
    ),
    otherdm = list(
      include_std = integer(0),
      include_src = c(1567906, 1567923, 1567972, 44832187),
      exclude_std = integer(0),
      exclude_src = integer(0)
    )
  )
  # remove rows with missing person_id or index_date; format dates for safe SQL injection
  valid_rows <- cohort_df %>% filter(!is.na(person_id) & !is.na(index_date))
  if (nrow(valid_rows) == 0) {
    # nothing to query — return empty flags frame (to be joined upstream)
    cohort_cte_sql <- ""
  } else {
    # ensure index_date is Date and format as YYYY-MM-DD
    valid_rows <- valid_rows %>% mutate(index_date = as.Date(index_date))
    cohort_rows <- purrr::map2_chr(valid_rows$person_id, valid_rows$index_date,
                                   ~ glue::glue("STRUCT({.x} AS person_id, DATE('{format(.y, '%Y-%m-%d')}') AS index_date)"))
    cohort_cte_sql <- paste(cohort_rows, collapse = ", ")
  }
  ctes <- c()
  for (name in names(comorbidity_sets)) {
    dspec <- comorbidity_sets[[name]]
    if (length(dspec$include_src) > 0) ctes <- c(ctes, build_cte(dspec$include_src, glue::glue("inc_src_desc_{name}"), 0, dataset))
  }
  with_clause <- if (length(ctes) > 0) glue::glue_collapse(ctes, sep = ",\n") else ""
  joins <- c()
  for (name in names(comorbidity_sets)) {
    dspec <- comorbidity_sets[[name]]
    if (length(dspec$include_src) == 0) next
    joins <- c(joins, glue::glue("
      LEFT JOIN (
        SELECT DISTINCT t.person_id
        FROM `{dataset}.cb_search_all_events` se
        JOIN cohort t ON se.person_id = t.person_id
        WHERE se.entry_date BETWEEN DATE_SUB(t.index_date, INTERVAL 730 DAY) AND t.index_date
          AND (se.is_standard = 0 AND se.concept_id IN (SELECT concept_id FROM inc_src_desc_{name}))
      ) com_{name} ON t.person_id = com_{name}.person_id
    "))
  }
  joins_sql <- paste(joins, collapse = "\n")
  flags <- purrr::map_chr(names(comorbidity_sets), ~ glue::glue("CASE WHEN com_{.x}.person_id IS NOT NULL THEN 1 ELSE 0 END AS {.x}"))
  flags_sql <- paste(flags, collapse = ",\n ")
  sql_comorbidities <- glue::glue("
    WITH
    cohort AS (
      SELECT * FROM UNNEST(ARRAY<STRUCT<person_id INT64, index_date DATE>>[{cohort_cte_sql}])
    ){if (with_clause != '') glue::glue(',\n{with_clause}') else ''}
    SELECT
      t.person_id,
      t.index_date,
      {flags_sql}
    FROM cohort t
    {joins_sql}
  ")
  comorbidity_flags <- download_data(sql_comorbidities)
  if ("person_id" %in% names(comorbidity_flags)) comorbidity_flags <- comorbidity_flags %>% mutate(person_id = as.character(person_id))
  if ("index_date" %in% names(comorbidity_flags)) comorbidity_flags <- comorbidity_flags %>% mutate(index_date = as.Date(index_date))
  cohort_df <- cohort_df %>% left_join(comorbidity_flags, by = c("person_id", "index_date")) %>%
    mutate(across(names(comorbidity_sets), ~ tidyr::replace_na(.x, 0)))
  cohort_df
}
# >>> CELL 07: Baseline Measurements & Comparison Processor <<<
# ---------- BASELINE MEASUREMENTS ----------
add_baseline_measurements <- function(df, bmi_panel = NULL, hba1c_panel = NULL) {
  result <- df
  if (!is.null(bmi_panel)) {
    baseline_bmi_data <- df %>%
      select(person_id, index_date) %>% distinct() %>%
      left_join(bmi_panel, by = "person_id", relationship = "many-to-many") %>%
      filter(!is.na(weight_date) & !is.na(index_date) & weight_date <= index_date) %>%
      mutate(days_from_index = as.numeric(difftime(index_date, weight_date, units = "days"))) %>%
      group_by(person_id) %>% arrange(days_from_index, desc(weight_date)) %>% slice(1) %>% ungroup() %>%
      select(person_id, baseline_bmi = bmi)
    result <- result %>% left_join(baseline_bmi_data, by = "person_id")
  }
  if (!is.null(hba1c_panel)) {
    baseline_hba1c_data <- df %>%
      select(person_id, index_date) %>% distinct() %>%
      left_join(hba1c_panel, by = "person_id", relationship = "many-to-many") %>%
      filter(!is.na(date_of_measurement) & !is.na(index_date) & date_of_measurement <= index_date) %>%
      mutate(days_from_index = as.numeric(difftime(index_date, date_of_measurement, units = "days"))) %>%
      group_by(person_id) %>% arrange(days_from_index, desc(date_of_measurement)) %>% slice(1) %>% ungroup() %>%
      select(person_id, baseline_hba1c = A1c)
    result <- result %>% left_join(baseline_hba1c_data, by = "person_id")
  }
  
  result <- result %>%
    mutate(
      index_year = factor(lubridate::year(index_date)),
      # Ensure index_year is numeric for comparison
      index_year_numeric = as.numeric(as.character(index_year)),
      
      # Create the new grouped variable
      index_year_grouped = factor(case_when(
        index_year_numeric %in% c(2018, 2019, 2023) ~ "Non-COVID",
        index_year_numeric %in% 2020:2022 ~ "COVID",
        TRUE ~ "Other" # A fallback for any other years
      )),
      
      baseline_bmi_category = factor(
        case_when(
          is.na(baseline_bmi) ~ NA_character_,
          baseline_bmi < 18.5 ~ "Underweight",
          baseline_bmi >= 18.5 & baseline_bmi < 25 ~ "Normal",
          baseline_bmi >= 25 & baseline_bmi < 30 ~ "Overweight",
          baseline_bmi >= 30 & baseline_bmi < 35 ~ "Obesity Class I",
          baseline_bmi >= 35 & baseline_bmi < 40 ~ "Obesity Class II",
          baseline_bmi >= 40 ~ "Obesity Class III"
        ),
        # This sets a logical order for the categories
        levels = c(
          "Underweight",
          "Normal",
          "Overweight",
          "Obesity Class I",
          "Obesity Class II",
          "Obesity Class III")))
  result
}
# ---------- COMPARISON GROUP PROCESSOR ----------
process_comparison_group <- function(exposure_idx, comparator_idx, base_cohort, exposure_name, comparator_name) {
  # ensure the index column is named 'index_date' (some callers previously renamed it)
  if (!"index_date" %in% names(exposure_idx) && ncol(exposure_idx) >= 2) exposure_idx <- exposure_idx %>% rename(index_date = 2)
  if (!"index_date" %in% names(comparator_idx) && ncol(comparator_idx) >= 2) comparator_idx <- comparator_idx %>% rename(index_date = 2)
  exposure_idx <- exposure_idx %>% mutate(person_id = as.character(person_id), index_date = as.Date(index_date))
  comparator_idx <- comparator_idx %>% mutate(person_id = as.character(person_id), index_date = as.Date(index_date))
  combined_ids <- unique(c(exposure_idx$person_id, comparator_idx$person_id))
  combined_cohort <- base_cohort %>% filter(person_id %in% combined_ids)
  exposure_only_ids <- setdiff(exposure_idx$person_id, comparator_idx$person_id)
  comparator_only_ids <- setdiff(comparator_idx$person_id, exposure_idx$person_id)
  exposure_only_df <- data.frame(person_id = exposure_only_ids, treatment = 1,
                                 index_date = exposure_idx$index_date[match(exposure_only_ids, exposure_idx$person_id)],
                                 crossover_date = as.Date(NA),
                                 stringsAsFactors = FALSE)
  comparator_only_df <- data.frame(person_id = comparator_only_ids, treatment = 0,
                                   index_date = comparator_idx$index_date[match(comparator_only_ids, comparator_idx$person_id)],
                                   crossover_date = as.Date(NA),
                                   stringsAsFactors = FALSE)
  both_drugs_ids <- intersect(exposure_idx$person_id, comparator_idx$person_id)
  if (length(both_drugs_ids) > 0) {
    exposure_dates <- exposure_idx %>% filter(person_id %in% both_drugs_ids) %>% select(person_id, exposure_date = index_date)
    comparator_dates <- comparator_idx %>% filter(person_id %in% both_drugs_ids) %>% select(person_id, comparator_date = index_date)
    temporal_df <- merge(exposure_dates, comparator_dates, by = "person_id")
    sequential_df <- temporal_df %>%
      mutate(treatment = if_else(exposure_date < comparator_date, 1, 0),
             index_date = as.Date(if_else(exposure_date < comparator_date, exposure_date, comparator_date)),
             crossover_date = as.Date(if_else(exposure_date < comparator_date, comparator_date, exposure_date))) %>%
      select(person_id, treatment, index_date, crossover_date)
    simultaneous_ids <- temporal_df$person_id[temporal_df$exposure_date == temporal_df$comparator_date]
    sequential_df <- sequential_df %>% filter(!(person_id %in% simultaneous_ids))
  } else {
    sequential_df <- data.frame(person_id = character(0), treatment = integer(0),
                                index_date = as.Date(character(0)), crossover_date = as.Date(character(0)),
                                stringsAsFactors = FALSE)
  }
  cohort_pair <- bind_rows(exposure_only_df, comparator_only_df, sequential_df)
  # avoid losing cohort_pair$index_date via suffixing on join: drop index_date from base_cohort if present
  if ("index_date" %in% names(base_cohort)) base_cohort <- base_cohort %>% select(-index_date)
  cohort_df <- cohort_pair %>% left_join(base_cohort, by = "person_id")
  cohort_df <- cohort_df %>%
    filter(!is.na(t2d_excluded_t1d_start_date) & t2d_excluded_t1d_start_date <= index_date) %>%
    mutate(index_date = as.Date(index_date))
  
  # Diagnostics: pre-index epilepsy/seizure exclusions (before filtering)
  pre_mask <- !is.na(cohort_df$epilepsy_or_seizure_start_date) &
    (as.Date(cohort_df$epilepsy_or_seizure_start_date) < as.Date(cohort_df$index_date))
  n_pre_total <- sum(pre_mask, na.rm = TRUE)
  if (!is.na(n_pre_total) && n_pre_total > 0) {
    comp_label <- paste0(exposure_name, " vs ", comparator_name)
    by_trt <- tryCatch({
      tbl <- cohort_df[pre_mask, c("treatment")]
      if (nrow(tbl) == 0 || !"treatment" %in% names(tbl)) data.frame(treatment = integer(0), n = integer(0)) else {
        as.data.frame(stats::aggregate(rep(1, nrow(tbl)) ~ tbl$treatment, FUN = sum))
      }
    }, error = function(e) data.frame())
    log_output("PROGRESS", "Pre-index epilepsy/seizure detected (will be excluded)",
               sprintf("%s: n_total=%d%s",
                       comp_label, n_pre_total,
                       if (nrow(by_trt) > 0) paste0("; by_treatment=", paste0(by_trt[[1]], ":", by_trt[[2]], collapse=",")) else ""),
               "process_comparison_group")
    # Save IDs for audit trail
    out_name <- sprintf("preindex_epilepsy_exclusions_%s_%s.csv",
                        gsub("[^A-Za-z0-9]+", "_", comp_label),
                        format(Sys.time(), "%Y%m%d_%H%M%S"))
    safe_cols <- intersect(c("person_id","index_date","epilepsy_or_seizure_start_date","treatment"), names(cohort_df))
    if (length(safe_cols) > 0) {
      try(utils::write.csv(cohort_df[pre_mask, safe_cols, drop = FALSE], file = out_name, row.names = FALSE), silent = TRUE)
    }
  }
  
  cohort_df <- cohort_df %>%
    mutate(preindex_epilepsy_flag = ifelse(pre_mask, 1L, 0L)) %>%
    filter(preindex_epilepsy_flag == 0L) 
  list(temporal_stats = list(final_cohort = nrow(cohort_df)), cohort_df = cohort_df)
}
create_cohort_summary_df <- function(cohort_stats, comparison_results) {
  # This part, handling the overall stats, is correct and remains unchanged.
  main_stats <- data.frame(metric = names(cohort_stats), value = unlist(cohort_stats), stringsAsFactors = FALSE)
  
  comparison_dfs <- list()
  
  # This loop is now corrected to handle potentially empty cohorts.
  for (comp_name in names(comparison_results)) {
    # Directly access the final, cleaned cohort dataframe.
    final_df <- comparison_results[[comp_name]]$cohort_df
    
    # Calculate the final cohort size. If the dataframe is NULL or empty, the count is 0.
    final_count <- if (is.null(final_df) || nrow(final_df) == 0) 0 else nrow(final_df)
    
    # Create the stats list on-the-fly using the accurate final count.
    comp_data <- list(final_cohort = final_count)
    
    # This data.frame creation will now always work because comp_data is never empty.
    comp_df <- data.frame(
      comparison = comp_name,
      metric = names(comp_data),
      value = unlist(comp_data),
      stringsAsFactors = FALSE
    )
    comparison_dfs[[comp_name]] <- comp_df
  }
  
  # The rest of the function remains the same.
  all_comparisons <- do.call(rbind, comparison_dfs)
  main_stats$comparison <- "overall"
  
  # Ensure consistent column order for binding.
  all_comparisons <- all_comparisons[, c("comparison", "metric", "value")]
  
  rbind(main_stats[, c("comparison", "metric", "value")], all_comparisons)
}

get_drug_concepts <- function(names_vec, dataset = NULL) {
  if (is.null(dataset)) {
    dataset <- Sys.getenv("WORKSPACE_CDR")
  }
  
  pattern <- paste0(
    'LOWER(c.concept_name) LIKE "%', tolower(names_vec), '%"',
    collapse = " OR "
  )
  sql <- glue("
    SELECT DISTINCT c2.concept_id
    FROM `{dataset}.concept` c
    JOIN `{dataset}.concept_ancestor` ca ON c.concept_id = ca.ancestor_concept_id
    JOIN `{dataset}.concept` c2 ON c2.concept_id = ca.descendant_concept_id
    WHERE c.concept_class_id = 'Ingredient' AND ({pattern})
  ")
  download_data(sql)
}

# ============================================================================
# BOOTSTRAP HELPERS (nonparametric; resample rows; returns RD samples + summary)
# ============================================================================
config <- get_standard_config()
perform_tmle_bootstrap <- function(df,
                                   outcome_var = "event",
                                   exclude_vars = NULL,
                                   all_ps_vars = NULL,
                                   sl_lib = c("SL.glm"),
                                   n_bootstrap = 1000,
                                   seed = 12345,
                                   min_events_per_arm = 1) {
  # Parallel bootstrap using futures
  if (!requireNamespace("future", quietly = TRUE)) install.packages("future")
  if (!requireNamespace("future.apply", quietly = TRUE)) install.packages("future.apply")
  # Prevent nested over-threading inside workers
  Sys.setenv(OMP_NUM_THREADS = "1", MKL_NUM_THREADS = "1", OPENBLAS_NUM_THREADS = "1")
  if (requireNamespace("data.table", quietly = TRUE)) data.table::setDTthreads(1)
  
  n <- nrow(df)
  set.seed(seed)
  idx_list <- lapply(seq_len(n_bootstrap), function(b) sample.int(n, size = n, replace = TRUE))
  
  future::plan(future::multisession, workers = max(1, parallel::detectCores() - 1))
  
  res <- future.apply::future_lapply(seq_along(idx_list), function(b) {
    boot_df <- df[idx_list[[b]], , drop = FALSE]
    t1 <- boot_df$treatment == 1
    t0 <- boot_df$treatment == 0
    if (sum(t1) < 2 || sum(t0) < 2 ||
        sum(boot_df[[outcome_var]][t1], na.rm = TRUE) < min_events_per_arm ||
        sum(boot_df[[outcome_var]][t0], na.rm = TRUE) < min_events_per_arm) {
      return(NA_real_)
    }
    
    fit_b <- try(
      run_tmle_binary(
        df = boot_df,
        exclude_vars = exclude_vars,
        all_ps_vars = all_ps_vars,
        sl_lib = sl_lib,
        verbose = FALSE,
        use_imputation = FALSE,
        use_superlearner = TRUE
      ),
      silent = TRUE
    )
    if (inherits(fit_b, "try-error") || is.null(fit_b) ||
        is.null(fit_b$fit$estimates$ATE$psi) || !is.finite(fit_b$fit$estimates$ATE$psi)) {
      return(NA_real_)
    }
    as.numeric(fit_b$fit$estimates$ATE$psi)
  }, future.seed = TRUE)
  
  rds <- unlist(res)
  samples <- rds[is.finite(rds)]
  fails <- sum(!is.finite(rds))
  if (length(samples) == 0L) {
    return(list(
      samples = numeric(0),
      mean = NA_real_, sd = NA_real_,
      ci_lower_2.5 = NA_real_, ci_upper_97.5 = NA_real_,
      n_success = 0L, n_fail = n_bootstrap
    ))
  }
  list(
    samples = samples,
    mean = mean(samples),
    sd = stats::sd(samples),
    ci_lower_2.5 = as.numeric(stats::quantile(samples, probs = 0.025, names = FALSE)),
    ci_upper_97.5 = as.numeric(stats::quantile(samples, probs = 0.975, names = FALSE)),
    n_success = length(samples),
    n_fail = fails
  )
}

# ============================================================================
# Save bootstrap outputs (CSV summary + RDS with full samples; optional PNG)
# ============================================================================

save_bootstrap_results <- function(bootstrap_stats,
                                   comparison,
                                   outcome,
                                   analysis_prefix,
                                   tmle_point_estimate = NA_real_,
                                   write_hist_png = TRUE) {
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  safe_comp   <- gsub("[^A-Za-z0-9]+", "_", comparison)
  safe_outc   <- gsub("[^A-Za-z0-9]+", "_", outcome)
  base_prefix <- sprintf("%s_bootstrap_%s_%s_%s", analysis_prefix, safe_comp, safe_outc, ts)
  
  # 1) RDS with full results and samples
  rds_path <- sprintf("%s.rds", base_prefix)
  saveRDS(bootstrap_stats, file = rds_path)
  
  # 2) CSV summary
  summ_df <- data.frame(
    comparison = comparison,
    outcome = outcome,
    n_success = bootstrap_stats$n_success,
    n_fail = bootstrap_stats$n_fail,
    mean_rd = bootstrap_stats$mean,
    sd_rd = bootstrap_stats$sd,
    ci_lower_2.5 = bootstrap_stats$ci_lower_2.5,
    ci_upper_97.5 = bootstrap_stats$ci_upper_97.5,
    tmle_point_estimate = tmle_point_estimate,
    stringsAsFactors = FALSE
  )
  csv_path <- sprintf("%s_summary.csv", base_prefix)
  utils::write.csv(summ_df, file = csv_path, row.names = FALSE)
  
  # 3) Optional histogram PNG
  png_path <- NA_character_
  if (write_hist_png && length(bootstrap_stats$samples) > 0) {
    png_path <- sprintf("%s_hist.png", base_prefix)
    grDevices::png(png_path, width = 1000, height = 700, res = 120)
    hist(bootstrap_stats$samples,
         main = sprintf("Bootstrap RD: %s – %s", comparison, outcome),
         xlab = "Risk Difference (ATE)", breaks = 30)
    abline(v = tmle_point_estimate, lty = 2, lwd = 2)
    abline(v = 0, lty = 3)
    legend("topright",
           legend = c("TMLE point estimate", "Null (0)"),
           lty = c(2, 3), lwd = c(2, 1))
    grDevices::dev.off()
  }
  
  list(rds = rds_path, csv = csv_path, png = png_path)
}

# ============================================================================
# Combined Visualization for Multiple Bootstrap Comparisons
# ============================================================================

visualize_combined_bootstrap <- function(comparisons,
                                         outcome = "Epilepsy/Seizure excl Hypoglycemic",
                                         analysis_prefix = "epilepsy_seizure_late_onset_tmle",
                                         colors = NULL,
                                         use_jama_style = TRUE,
                                         plot_type = "multi") {  # "single", "multi", or "both"
  jama_colors <- c("#00274C", "#C4622D", "#4B7AA7", "#969696", "#2F7F7E")
  if (is.null(colors)) {
    colors <- if (use_jama_style) jama_colors else c("darkblue", "darkred", "darkgreen", "purple", "orange")
  }
  
  if (use_jama_style) {
    old_par <- par(no.readonly = TRUE)
    par(
      family = "sans",
      cex.main = 1.1,
      cex.lab = 1.0,
      cex.axis = 0.9,
      font.main = 1,
      las = 1,
      mgp = c(2.5, 0.7, 0),
      tcl = -0.3
    )
  }
  
  all_bootstrap_data <- list()
  
  for (i in seq_along(comparisons)) {
    comparison <- comparisons[i]
    safe_comp <- gsub("[^A-Za-z0-9]+", "_", comparison)
    safe_outc <- gsub("[^A-Za-z0-9]+", "_", outcome)
    pattern <- sprintf("^%s_bootstrap_%s_%s_\\d{8}_\\d{6}\\.rds$", analysis_prefix, safe_comp, safe_outc)
    files <- list.files(pattern = pattern)
    if (length(files) == 0) {
      warning(sprintf("No RDS bootstrap files found for %s (%s) with prefix %s", comparison, outcome, analysis_prefix))
      next
    }
    mt <- file.info(files)$mtime
    rds_filename <- files[order(mt, decreasing = TRUE)][1]
    bootstrap_stats <- readRDS(rds_filename)
    all_bootstrap_data[[comparison]] <- bootstrap_stats
  }
  
  if (length(all_bootstrap_data) == 0) stop("No bootstrap data could be loaded")
  
  all_samples <- unlist(lapply(all_bootstrap_data, function(x) x$samples))
  x_range <- range(c(all_samples, 0)) * 1.1
  densities <- lapply(all_bootstrap_data, function(x) density(x$samples))
  y_max <- max(sapply(densities, function(d) max(d$y))) * 1.1
  
  if (plot_type == "single" || plot_type == "both") {
    plot(NULL,
         xlim = x_range,
         ylim = c(0, y_max),
         main = if(use_jama_style) paste("Bootstrap Distributions:", outcome) else paste("Combined Bootstrap Densities:", outcome),
         xlab = "Risk Difference",
         ylab = "Density",
         type = "n",
         axes = FALSE,
         frame.plot = FALSE)
    axis(1, col = "black", col.axis = "black", lwd = 0.5)
    axis(2, col = "black", col.axis = "black", lwd = 0.5)
    if (use_jama_style) {
      abline(h = axTicks(2), col = "gray90", lty = 1, lwd = 0.5)
      abline(v = axTicks(1), col = "gray90", lty = 1, lwd = 0.5)
    }
    box(lwd = 0.5)
    abline(v = 0, col = "black", lwd = 1.5, lty = 1)
    for (i in seq_along(all_bootstrap_data)) {
      dens <- densities[[i]]
      col_i <- colors[(i - 1) %% length(colors) + 1]
      lines(dens, col = col_i, lwd = if(use_jama_style) 2 else 3)
      polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0, length(dens$y))),
              col = adjustcolor(col_i, alpha.f = if(use_jama_style) 0.15 else 0.2), border = NA)
      abline(v = all_bootstrap_data[[i]]$mean, col = col_i, lwd = if(use_jama_style) 1 else 2, lty = 2)
    }
  }
  
  if (plot_type == "multi" || plot_type == "both") {
    par(mfrow = c(2, 2))
    plot(NULL, xlim = x_range, ylim = c(0, y_max),
         main = if(use_jama_style) "A. Density Distributions" else "Combined Bootstrap Densities",
         xlab = "Risk Difference", ylab = "Density",
         type = "n", axes = FALSE, frame.plot = FALSE)
    axis(1, col = "black", col.axis = "black", lwd = 0.5)
    axis(2, col = "black", col.axis = "black", lwd = 0.5)
    if (use_jama_style) {
      abline(h = axTicks(2), col = "gray90", lty = 1, lwd = 0.5)
      abline(v = axTicks(1), col = "gray90", lty = 1, lwd = 0.5)
    }
    box(lwd = 0.5)
    abline(v = 0, col = "black", lwd = 1.5, lty = 1)
    for (i in seq_along(all_bootstrap_data)) {
      dens <- densities[[i]]
      col_i <- colors[(i - 1) %% length(colors) + 1]
      lines(dens, col = col_i, lwd = if(use_jama_style) 2 else 3)
      polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0, length(dens$y))),
              col = adjustcolor(col_i, alpha.f = if(use_jama_style) 0.15 else 0.2), border = NA)
    }
    
    boxplot_data <- lapply(all_bootstrap_data, function(x) x$samples)
    boxplot(boxplot_data,
            names = names(all_bootstrap_data),
            col = adjustcolor(colors[seq_along(all_bootstrap_data)], alpha.f = if(use_jama_style) 0.3 else 0.5),
            border = colors[seq_along(all_bootstrap_data)],
            main = if(use_jama_style) "B. Distribution Comparison" else "Bootstrap Distributions Comparison",
            ylab = "Risk Difference", las = 2, axes = FALSE, frame.plot = FALSE,
            lwd = if(use_jama_style) 1 else 1.5, boxwex = if(use_jama_style) 0.6 else 0.8,
            staplewex = 0.5, outwex = 0.5)
    axis(1, at = seq_along(all_bootstrap_data), labels = names(all_bootstrap_data), las = 2, lwd = 0.5)
    axis(2, lwd = 0.5)
    if (use_jama_style) abline(h = axTicks(2), col = "gray90", lty = 1, lwd = 0.5)
    box(lwd = 0.5)
    abline(h = 0, col = "black", lwd = 1.5, lty = 1)
    
    plot(NULL, xlim = c(0.5, length(all_bootstrap_data) + 0.5), ylim = x_range,
         main = if(use_jama_style) "C. 95% Confidence Intervals" else "95% Confidence Intervals",
         xlab = "", ylab = "Risk Difference", xaxt = "n", axes = FALSE, frame.plot = FALSE)
    axis(1, at = seq_along(all_bootstrap_data), labels = names(all_bootstrap_data), las = 2, lwd = 0.5)
    axis(2, lwd = 0.5)
    if (use_jama_style) abline(h = axTicks(2), col = "gray90", lty = 1, lwd = 0.5)
    box(lwd = 0.5)
    abline(h = 0, col = "black", lwd = 1.5, lty = 1)
    for (i in seq_along(all_bootstrap_data)) {
      d <- all_bootstrap_data[[i]]
      col_i <- colors[(i - 1) %% length(colors) + 1]
      segments(i, d$ci_lower_2.5, i, d$ci_upper_97.5, col = col_i, lwd = if(use_jama_style) 3 else 4)
      points(i, d$mean, pch = 19, col = col_i, cex = if(use_jama_style) 1.5 else 2)
    }
    
    plot.new(); plot.window(xlim = c(0, 1), ylim = c(0, 1))
    title(if(use_jama_style) "D. Statistical Summary" else "Statistical Summary")
    y_pos <- 0.9
    for (i in seq_along(all_bootstrap_data)) {
      nm <- names(all_bootstrap_data)[i]
      d <- all_bootstrap_data[[i]]
      col_i <- colors[(i - 1) %% length(colors) + 1]
      text(0.05, y_pos, nm, col = col_i, font = if(use_jama_style) 1 else 2, adj = 0, cex = if(use_jama_style) 0.85 else 0.9)
      y_pos <- y_pos - 0.08
      text(0.1, y_pos, paste("Mean:", sprintf("%.5f", d$mean)), adj = 0, cex = if(use_jama_style) 0.75 else 0.8)
      y_pos <- y_pos - 0.06
      text(0.1, y_pos, paste("95% CI: (", sprintf("%.5f", d$ci_lower_2.5), ", ", sprintf("%.5f", d$ci_upper_97.5), ")"), adj = 0, cex = if(use_jama_style) 0.75 else 0.8)
      y_pos <- y_pos - 0.06
      text(0.1, y_pos, paste("P(RD<0):", sprintf("%.3f", mean(d$samples < 0))), adj = 0, cex = if(use_jama_style) 0.75 else 0.8)
      y_pos <- y_pos - 0.1
    }
    par(mfrow = c(1, 1))
  }
  
  if (use_jama_style && exists("old_par")) par(old_par)
  invisible(all_bootstrap_data)
}

# ============================================================================
# MAIN TMLE ANALYSIS SCRIPT (tmle package; reuses IPTW cohorts)
# ============================================================================

if (!requireNamespace("tmle", quietly = TRUE)) install.packages("tmle")
if (!requireNamespace("SuperLearner", quietly = TRUE)) install.packages("SuperLearner")
suppressPackageStartupMessages(library(tmle))

# SuperLearner library preference
tmle_sl_lib <- c("SL.glm")
if (requireNamespace("glmnet", quietly = TRUE)) tmle_sl_lib <- c("SL.glmnet", tmle_sl_lib)
if (requireNamespace("xgboost", quietly = TRUE)) tmle_sl_lib <- c("SL.xgboost", tmle_sl_lib)

# Helper: run TMLE on binary event outcome using same covariates and complete-case set
run_tmle_binary <- function(df, exclude_vars = NULL, all_ps_vars = NULL,
                            sl_lib = c("SL.glm"), verbose = FALSE,
                            use_imputation = FALSE,
                            use_superlearner = TRUE) {
  
  all_covariates <- if (is.null(all_ps_vars)) get_standard_config()$all_ps_vars else all_ps_vars
  rhs_vars <- setdiff(all_covariates, exclude_vars)
  keep_vars <- rhs_vars[sapply(df[rhs_vars], function(x) length(unique(na.omit(x))) > 1)]
  model_vars <- c("treatment", keep_vars)
  
  if (use_imputation) {
    df_impute <- perform_mice_imputation(df, m = 20, seed = 123)
  } else {
    df_impute <- df
  }
  
  complete_rows <- complete.cases(df_impute[, model_vars, drop = FALSE])
  df_complete <- df_impute[complete_rows, , drop = FALSE]
  
  W_vars <- df_complete[, keep_vars, drop = FALSE]
  if (ncol(W_vars) == 0) W_vars <- data.frame(intercept = rep(1, nrow(df_complete)))
  
  use_sl <- isTRUE(use_superlearner) && requireNamespace("SuperLearner", quietly = TRUE)
  fit <- NULL
  
  if (use_sl) {
    sl_present <- unique(sl_lib)
    fit <- try(tmle::tmle(
      Y = df_complete$event,
      A = df_complete$treatment,
      W = W_vars,
      family = "binomial",
      Q.SL.library = sl_present,
      g.SL.library = sl_present
    ), silent = TRUE)
  }
  
  if (is.null(fit) || inherits(fit, "try-error")) {
    data_for_models <- cbind(df_complete[, c("event", "treatment"), drop = FALSE], W_vars)
    
    g_form <- as.formula(paste("treatment ~", paste(colnames(W_vars), collapse = " + ")))
    g_fit <- stats::glm(g_form, data = data_for_models, family = stats::binomial())
    g1W <- suppressWarnings(pmin(pmax(as.numeric(stats::predict(g_fit, type = "response")), 1e-5), 1 - 1e-5))
    
    q_form <- as.formula(paste("event ~ treatment +", paste(colnames(W_vars), collapse = " + ")))
    q_fit <- stats::glm(q_form, data = data_for_models, family = stats::binomial())
    
    newdata1 <- data_for_models; newdata1$treatment <- 1
    newdata0 <- data_for_models; newdata0$treatment <- 0
    Q1W <- suppressWarnings(pmin(pmax(as.numeric(stats::predict(q_fit, newdata = newdata1, type = "response")), 1e-5), 1 - 1e-5))
    Q0W <- suppressWarnings(pmin(pmax(as.numeric(stats::predict(q_fit, newdata = newdata0, type = "response")), 1e-5), 1 - 1e-5))
    
    fit <- tmle::tmle(
      Y = df_complete$event,
      A = df_complete$treatment,
      W = W_vars,
      family = "binomial",
      Q = cbind(Q1W = Q1W, Q0W = Q0W),
      g1W = g1W
    )
  }
  
  list(
    fit = fit,
    cohort = df_complete,
    keep_vars = keep_vars,
    n_dropped = nrow(df) - nrow(df_complete)
  )
}

# Helper: time-specific RD (Kaplan-Meier based)
calculate_tmle_time_specific_rd <- function(df, time_points, time_labels) {
  safe_len <- min(length(time_points), length(time_labels))
  time_points <- time_points[seq_len(safe_len)]
  time_labels <- time_labels[seq_len(safe_len)]
  
  calc_surv_at <- function(data_subset, t_point) {
    fit <- survival::survfit(survival::Surv(event_time, event) ~ 1, data = data_subset)
    s <- suppressWarnings(summary(fit, times = t_point))
    se <- if (!is.null(s$std.err) && length(s$std.err) > 0) s$std.err[1] else NA_real_
    list(surv = s$surv[1], se = se)
  }
  
  out <- vector("list", length(time_points))
  for (i in seq_along(time_points)) {
    t_point <- time_points[i]
    lbl <- time_labels[i]
    
    grp0 <- df[df$treatment == 0, , drop = FALSE]
    grp1 <- df[df$treatment == 1, , drop = FALSE]
    
    s0 <- calc_surv_at(grp0, t_point)
    s1 <- calc_surv_at(grp1, t_point)
    
    risk0 <- 1 - s0$surv
    risk1 <- 1 - s1$surv
    rd <- risk1 - risk0
    
    var_rd <- s0$se^2 + s1$se^2
    se_rd <- sqrt(var_rd)
    lo <- rd - 1.96 * se_rd
    hi <- rd + 1.96 * se_rd
    
    out[[i]] <- data.frame(
      time_label = lbl,
      time_days = t_point,
      risk_treated = risk1,
      risk_control = risk0,
      rd = rd,
      rd_lo_95 = lo,
      rd_hi_95 = hi,
      stringsAsFactors = FALSE
    )
  }
  do.call(rbind, out)
}

# ---- Debug/Run Flags (flip during development for speed) ----
RUN_TMLE <- TRUE
WRITE_OUTPUTS <- TRUE
LIMIT_COMPARISONS <- NULL
COVARIATE_CAP <- NULL
IPTW_RDS_PATH_OVERRIDE <- NULL   # if set, must be a *single* RDS (use only for one analysis)

analysis_types <- list(
  list(name = "all_ages",   late_onset_flag = FALSE, prefix = "outcome4_epilepsy_excl_hypo_all_ages_ipwt"),
  list(name = "late_onset", late_onset_flag = TRUE,  prefix = "outcome4_epilepsy_excl_hypo_late_onset_ipwt")
)

# Comparisons/outcome you want
allowed_comparisons <- c("SEMAGLUTIDE vs OtherGLD", "SEMAGLUTIDE vs SGLT2")
tmle_outcome_label  <- "Epilepsy/Seizure excl Hypoglycemic"

# Run TMLE for all analysis types (all_ages + late_onset)
for (analysis in analysis_types) {
  
  if (!isTRUE(RUN_TMLE)) next
  
  log_output("PROGRESS",
             sprintf("STARTING TMLE ANALYSIS TYPE: %s", analysis$name),
             "", "tmle_main")
  
  analysis_t0 <- Sys.time()
  
  # Load the exact IPTW cohorts for THIS analysis type
  iptw_rds_path <- if (!is.null(IPTW_RDS_PATH_OVERRIDE)) {
    IPTW_RDS_PATH_OVERRIDE
  } else {
    sprintf("%s_full_results.rds", analysis$prefix)
  }
  
  if (!file.exists(iptw_rds_path)) {
    log_output("ERROR", "IPTW RDS not found", iptw_rds_path, "tmle_main")
    next
  }
  
  iptw_obj <- try(readRDS(iptw_rds_path), silent = TRUE)
  if (inherits(iptw_obj, "try-error") || is.null(iptw_obj)) {
    log_output("ERROR", "Failed to read IPTW RDS", iptw_rds_path, "tmle_main")
    next
  }
  
  # Output prefix: keep separate for all_ages vs late_onset
  tmle_prefix <- sub("_ipwt$", "_tmle", analysis$prefix)
  
  tmle_results <- list()
  tmle_summary_rows <- list()
  
  comps <- allowed_comparisons
  if (!is.null(LIMIT_COMPARISONS)) comps <- comps[seq_len(min(length(comps), LIMIT_COMPARISONS))]
  
  for (comp_name in comps) {
    
    log_output("PROGRESS",
               sprintf("TMLE: Processing comparison: %s", comp_name),
               "", "tmle_main")
    
    exclude_vars <- switch(comp_name,
                           "SEMAGLUTIDE vs OTHER_GLPA" = c("SEMAGLUTIDE", "OTHER_GLPA"),
                           "SEMAGLUTIDE vs SGLT2"      = c("SEMAGLUTIDE", "SGLT2i"),
                           "SEMAGLUTIDE vs OtherGLD"   = c("SEMAGLUTIDE", "TZD", "SU", "DPP4i"),
                           "OTHER_GLPA vs SGLT2"       = c("OTHER_GLPA", "SGLT2i"),
                           "OTHER_GLPA vs OtherGLD"    = c("OTHER_GLPA", "TZD", "SU", "DPP4i"),
                           "SGLT2 vs OtherGLD"         = c("SGLT2i", "TZD", "SU", "DPP4i"),
                           NULL
    )
    
    analysis_key <- paste(comp_name, "-", tmle_outcome_label)
    
    if (is.null(iptw_obj[[analysis_key]]) || is.null(iptw_obj[[analysis_key]]$cohort)) {
      log_output("WARNING", "Missing IPTW cohort for key", analysis_key, "tmle_main")
      next
    }
    
    survival_df <- iptw_obj[[analysis_key]]$cohort
    
    # Optional covariate cap (kept for parity with your flags)
    if (!is.null(COVARIATE_CAP) && is.numeric(COVARIATE_CAP) && COVARIATE_CAP > 0) {
      # cap is applied inside run_tmle_binary only if you implemented it there;
      # otherwise ignore here.
    }
    
    tmle_fit <- tryCatch({
      run_tmle_binary(
        df = survival_df,
        exclude_vars = exclude_vars,
        all_ps_vars = config$all_ps_vars,
        sl_lib = tmle_sl_lib,
        verbose = TRUE,
        use_imputation = FALSE,
        use_superlearner = TRUE
      )
    }, error = function(e) {
      log_output("ERROR", "run_tmle_binary failed", e$message, "tmle_main")
      NULL
    })
    if (is.null(tmle_fit)) next
    
    tmle_results[[analysis_key]] <- tmle_fit
    
    # Summary row (same structure as your original)
    ate <- tmle_fit$fit$estimates$ATE
    or  <- tmle_fit$fit$estimates$OR
    rr_obj <- try(tmle_fit$fit$estimates$RR, silent = TRUE)
    rr_est <- if (!inherits(rr_obj, "try-error") && !is.null(rr_obj$psi)) rr_obj$psi else NA_real_
    rr_lo  <- if (!inherits(rr_obj, "try-error") && !is.null(rr_obj$CI))  rr_obj$CI[1] else NA_real_
    rr_hi  <- if (!inherits(rr_obj, "try-error") && !is.null(rr_obj$CI))  rr_obj$CI[2] else NA_real_
    rr_p   <- if (!inherits(rr_obj, "try-error") && !is.null(rr_obj$pvalue)) rr_obj$pvalue else NA_real_
    nnt_val <- if (!is.null(ate$psi) && is.finite(ate$psi) && ate$psi != 0) 1 / abs(ate$psi) else Inf
    
    n_treat  <- sum(tmle_fit$cohort$treatment == 1, na.rm = TRUE)
    n_ctrl   <- sum(tmle_fit$cohort$treatment == 0, na.rm = TRUE)
    ev_treat <- sum(tmle_fit$cohort$event == 1 & tmle_fit$cohort$treatment == 1, na.rm = TRUE)
    ev_ctrl  <- sum(tmle_fit$cohort$event == 1 & tmle_fit$cohort$treatment == 0, na.rm = TRUE)
    
    out_df <- data.frame(
      method = "TMLE",
      comparison = comp_name,
      outcome = tmle_outcome_label,
      analysis_type = analysis$name,
      effect_measure = c("Risk Difference", "Odds Ratio", "Risk Ratio", "Number Needed to Treat"),
      estimate = c(ate$psi, or$psi, rr_est, nnt_val),
      lower_ci = c(ate$CI[1], or$CI[1], rr_lo, NA_real_),
      upper_ci = c(ate$CI[2], or$CI[2], rr_hi, NA_real_),
      p_value  = c(ate$pvalue, or$pvalue, rr_p, NA_real_),
      n_total  = rep(nrow(tmle_fit$cohort), 4),
      n_events = rep(sum(tmle_fit$cohort$event, na.rm = TRUE), 4),
      n_at_risk_treatment = rep(n_treat, 4),
      n_at_risk_control   = rep(n_ctrl, 4),
      events_treatment    = rep(ev_treat, 4),
      events_control      = rep(ev_ctrl, 4),
      stringsAsFactors = FALSE
    )
    tmle_summary_rows[[length(tmle_summary_rows) + 1]] <- out_df
    
    # Write per-comparison summary CSV
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    safe_comp    <- gsub("[^A-Za-z0-9]+", "_", comp_name)
    safe_outcome <- gsub("[^A-Za-z0-9]+", "_", tmle_outcome_label)
    
    if (isTRUE(WRITE_OUTPUTS)) {
      fname <- sprintf("%s_summary_%s_%s_%s.csv", tmle_prefix, safe_comp, safe_outcome, ts)
      utils::write.csv(out_df, file = fname, row.names = FALSE)
    }
    
    # Plot TMLE estimated event probabilities (pre-bootstrap)
    ey0_psi <- try(tmle_fit$fit$estimates$EY0$psi, silent = TRUE)
    ey1_psi <- try(tmle_fit$fit$estimates$EY1$psi, silent = TRUE)
    ey0_lo  <- try(tmle_fit$fit$estimates$EY0$CI[1], silent = TRUE)
    ey0_hi  <- try(tmle_fit$fit$estimates$EY0$CI[2], silent = TRUE)
    ey1_lo  <- try(tmle_fit$fit$estimates$EY1$CI[1], silent = TRUE)
    ey1_hi  <- try(tmle_fit$fit$estimates$EY1$CI[2], silent = TRUE)
    ate_obj <- try(tmle_fit$fit$estimates$ATE, silent = TRUE)
    ate_psi <- if (!inherits(ate_obj, "try-error") && !is.null(ate_obj$psi)) ate_obj$psi else NA_real_
    ate_p   <- if (!inherits(ate_obj, "try-error") && !is.null(ate_obj$pvalue)) ate_obj$pvalue else NA_real_
    
    to_num <- function(x) if (!inherits(x, "try-error")) as.numeric(x) else NA_real_
    
    prob_df <- data.frame(
      Group = c("Control", "Semaglutide"),
      Probability = c(to_num(ey0_psi), to_num(ey1_psi)),
      Lower = c(to_num(ey0_lo), to_num(ey1_lo)),
      Upper = c(to_num(ey0_hi), to_num(ey1_hi))
    )
    
    
    try({
      tmle_plot <- ggplot2::ggplot(prob_df, ggplot2::aes(x = Group, y = Probability, fill = Group)) +
        ggplot2::geom_bar(stat = "identity", alpha = 0.7) +
        ggplot2::geom_errorbar(ggplot2::aes(ymin = Lower, ymax = Upper), width = 0.2) +
        ggplot2::scale_fill_manual(values = c("Control" = "#1f77b4", "Semaglutide" = "#ff7f0e")) +
        ggplot2::labs(
          title    = paste("TMLE Estimated Event Probabilities for", tmle_outcome_label),
          subtitle = paste(comp_name, "- Risk Difference =", sprintf("%.4f", ate_psi),
                           "(p =", ifelse(is.na(ate_p), "NA", sprintf("%.3g", ate_p)), ")"),
          y        = "Probability of Event",
          x        = ""
        ) +
        ggplot2::theme_minimal() +
        ggplot2::theme(legend.position = "none")
      
      print(tmle_plot)
      
      safe_comp <- gsub("[^A-Za-z0-9]+", "_", comp_name)
      safe_outcome <- gsub("[^A-Za-z0-9]+", "_", tmle_outcome_label)
      safe_analysis <- gsub("[^A-Za-z0-9]+", "_", analysis$name)
      ggplot2::ggsave(
        filename = paste0("tmle_probs_", safe_analysis, "_", safe_comp, "_", safe_outcome, ".png"),
        plot = tmle_plot, width = 8, height = 6
      )
    }, silent = TRUE)
    
    # Bootstrap (kept prints)
    total_events <- sum(tmle_fit$cohort$event, na.rm = TRUE)
    if (total_events >= 10 &&
        exists("perform_tmle_bootstrap") && exists("save_bootstrap_results")) {
      
      cat("\n=== Starting TMLE Bootstrap (RD/ATE) ===\n")
      
      bs <- perform_tmle_bootstrap(
        df = tmle_fit$cohort,
        outcome_var = "event",
        exclude_vars = exclude_vars,
        all_ps_vars  = config$all_ps_vars,
        sl_lib       = tmle_sl_lib,
        n_bootstrap  = 1000,
        seed         = 12345
      )
      
      tmle_results[[analysis_key]]$bootstrap <- bs
      
      tmle_point <- as.numeric(tmle_fit$fit$estimates$ATE$psi)
      
      saved <- save_bootstrap_results(
        bootstrap_stats    = bs,
        comparison         = comp_name,
        outcome            = tmle_outcome_label,
        analysis_prefix    = tmle_prefix,
        tmle_point_estimate= tmle_point,
        write_hist_png     = TRUE
      )
      
      # Auto-visualize most recent bootstrap for this analysis type
      try({
        visualize_combined_bootstrap(
          comparisons     = comp_name,
          outcome         = tmle_outcome_label,
          analysis_prefix = tmle_prefix,
          use_jama_style  = TRUE,
          plot_type       = "single"
        )
      }, silent = TRUE)
      
      cat("Bootstrap saved:\n",
          " - RDS: ", saved$rds, "\n",
          " - CSV: ", saved$csv, "\n",
          if (!is.na(saved$png)) paste0(" - PNG: ", saved$png, "\n") else "",
          sep = "")
      
      cat(sprintf("Bootstrap RD 95%% CI: [%.4f, %.4f]; mean=%.4f; successes=%d; failures=%d\n",
                  bs$ci_lower_2.5, bs$ci_upper_97.5, bs$mean, bs$n_success, bs$n_fail))
      
    } else {
      tmle_results[[analysis_key]]$bootstrap <- NULL
      log_output("WARNING", "Bootstrap skipped (too few events)", total_events, "tmle_main")
    }
  }
  
  # Save full TMLE results per analysis type
  if (isTRUE(WRITE_OUTPUTS)) {
    rds_filename <- sprintf("%s_full_results.rds", tmle_prefix)
    saveRDS(tmle_results, file = rds_filename)
    log_output("PROGRESS", "TMLE results object saved", rds_filename, "tmle_main")
  }
  
  # Save aggregated summary across comparisons
  if (length(tmle_summary_rows) > 0 && isTRUE(WRITE_OUTPUTS)) {
    final_tmle_df <- do.call(rbind, tmle_summary_rows)
    ts2 <- format(Sys.time(), "%Y%m%d_%H%M%S")
    comp_results_file <- sprintf("%s_comprehensive_results_%s.csv", tmle_prefix, ts2)
    utils::write.csv(final_tmle_df, comp_results_file, row.names = FALSE)
    log_output("PROGRESS", "Aggregated TMLE comprehensive results saved", comp_results_file, "tmle_main")
  }
  
  log_output("PROGRESS",
             sprintf("TMLE ANALYSIS TYPE '%s' ELAPSED", analysis$name),
             sprintf("%.1fs", as.numeric(difftime(Sys.time(), analysis_t0, units = "secs"))),
             "tmle_main")
  
  log_output("PROGRESS",
             sprintf("TMLE ANALYSIS TYPE '%s' COMPLETE", analysis$name),
             "", "tmle_main")
}

# - For late_onset:
visualize_combined_bootstrap(
  comparisons     = c("SEMAGLUTIDE vs OtherGLD", "SEMAGLUTIDE vs SGLT2"),
  outcome         = "Epilepsy/Seizure excl Hypoglycemic",
  analysis_prefix = "outcome4_epilepsy_excl_hypo_late_onset_tmle",
  colors          = c("#00274C", "#C4622D"),
  use_jama_style  = TRUE,
  plot_type       = "single"
)