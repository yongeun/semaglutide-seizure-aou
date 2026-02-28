# >>> CELL 01: Install & Load Libraries <<<
install.packages("MatchIt")
install.packages("survminer")
install.packages("mice")
install.packages("tableone")

# Load required packages silently
suppressPackageStartupMessages({
  library(tidyverse)
  library(bigrquery)
  library(dplyr)
  library(tidyr)
  library(glue)
  library(MatchIt)
  library(survival)
  library(survminer)
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
    mutate(preindex_epilepsy_flag = ifelse(pre_mask, 1L, 0L))

  # Save pre-exclusion cohort (includes patients WITH pre-index epilepsy)
  pre_exclusion_df <- cohort_df

  cohort_df <- cohort_df %>%
    filter(preindex_epilepsy_flag == 0L)
  list(temporal_stats = list(final_cohort = nrow(cohort_df),
                             pre_exclusion_total = nrow(pre_exclusion_df)),
       cohort_df = cohort_df,
       pre_exclusion_cohort_df = pre_exclusion_df)
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

# >>> CELL 08: Main Cohort Building Workflow <<<
# --- SECTION 4: MAIN WORKFLOW ---
perform_stepwise_cohort_analysis_optimized <- function() {
  log_output("PROGRESS", "Starting Epilepsy/Seizure cohort analysis", "", "")
  config <- get_standard_config()
  data_list <- load_base_data()
  merged_df <- data_list$merged_df
  bmi_panel <- data_list$bmi_panel
  hba1c_panel <- data_list$hba1c_panel
  if ("person_id" %in% names(merged_df)) merged_df <- merged_df %>% mutate(person_id = as.character(person_id))
 
  cohort_stats <- list(total_participants = nrow(merged_df))
  diabetes_patients <- merged_df %>% filter(t2d_excluded_t1d == 1, age >= 18)
  cohort_stats$diabetes_patients <- nrow(diabetes_patients)
 
  sema_concepts <- get_drug_concepts(config$drug_classes$SEMAGLUTIDE)
  other_glpa_concepts <- get_drug_concepts(config$drug_classes$OTHER_GLPA)
  sglt2_concepts <- get_drug_concepts(config$drug_classes$SGLT2)
  othergld_concepts <- get_drug_concepts(config$drug_classes$OtherGLD)
 
  # Step 1: Build BROAD cohorts to define the overall population
  sema_idx_broad <- build_exposure_idx(sema_concepts, "sema_idx_broad", exposure_window = config$exposure_window)
  other_glpa_idx_broad <- build_exposure_idx(other_glpa_concepts, "other_glpa_idx_broad", exposure_window = config$exposure_window)
  sglt2_idx_broad <- build_exposure_idx(sglt2_concepts, "sglt2_idx_broad", exposure_window = config$exposure_window)
  othergld_idx_broad <- build_exposure_idx(othergld_concepts, "othergld_idx_broad", exposure_window = config$exposure_window)
 
  sema_idx_broad <- sema_idx_broad %>% filter(person_id %in% diabetes_patients$person_id)
  other_glpa_idx_broad <- other_glpa_idx_broad %>% filter(person_id %in% diabetes_patients$person_id)
  sglt2_idx_broad <- sglt2_idx_broad %>% filter(person_id %in% diabetes_patients$person_id)
  othergld_idx_broad <- othergld_idx_broad %>% filter(person_id %in% diabetes_patients$person_id)
 
  cohort_stats$sema_new_users_broad <- nrow(sema_idx_broad)
  cohort_stats$other_glpa_new_users_broad <- nrow(other_glpa_idx_broad)
  cohort_stats$sglt2_new_users_broad <- nrow(sglt2_idx_broad)
  cohort_stats$othergld_new_users_broad <- nrow(othergld_idx_broad)
 
  # Step 2: Create the demographics table based on this broad population
  all_drug_user_ids <- unique(c(sema_idx_broad$person_id, other_glpa_idx_broad$person_id, sglt2_idx_broad$person_id, othergld_idx_broad$person_id))
  dm_with_any_drug <- diabetes_patients %>% filter(person_id %in% all_drug_user_ids)
  demographics_to_join <- dm_with_any_drug %>%
    select(person_id, t2d_excluded_t1d_start_date, death_date, age, sex_cat, raceethnicity_cat, income, education, insurance_category, smoking, alcohol_category, EHRmaxDT_min, EHRmaxDT_max, epilepsy_or_seizure_start_date, epilepsy_refined_start_date, epilepsy_g40_start_date, epilepsy_or_seizure_without_hypoglycemic_seizure_start_date, adrd_start_date, stroke_start_date)
    
  # Step 3: Run the precise, comparison-specific cohort building
  comparison_results <- list()
  log_output("PROGRESS", "Starting comparison-specific cohort building", "", "")
 
  # --- Comparison 1: SEMAGLUTIDE vs OTHER_GLPA ---
  sema_idx_for_other_glpa <- build_exposure_idx(sema_concepts, "sema_for_other_glpa", exposure_window=config$exposure_window, washout_concepts_df=other_glpa_concepts)
  other_glpa_idx_for_sema <- build_exposure_idx(other_glpa_concepts, "other_glpa_for_sema", exposure_window=config$exposure_window, washout_concepts_df=sema_concepts)
  comparison_results$`SEMAGLUTIDE vs OTHER_GLPA` <- process_comparison_group(sema_idx_for_other_glpa, other_glpa_idx_for_sema, demographics_to_join, "SEMAGLUTIDE", "OTHER_GLPA")
 
  # --- Comparison 2: SEMAGLUTIDE vs SGLT2 ---
  sema_idx_for_sglt2 <- build_exposure_idx(sema_concepts, "sema_for_sglt2", exposure_window=config$exposure_window, washout_concepts_df=sglt2_concepts)
  sglt2_idx_for_sema <- build_exposure_idx(sglt2_concepts, "sglt2_for_sema", exposure_window=config$exposure_window, washout_concepts_df=sema_concepts)
  comparison_results$`SEMAGLUTIDE vs SGLT2` <- process_comparison_group(sema_idx_for_sglt2, sglt2_idx_for_sema, demographics_to_join, "SEMAGLUTIDE", "SGLT2")
 
  # --- Comparison 3: SEMAGLUTIDE vs OtherGLD ---
  sema_idx_for_othergld <- build_exposure_idx(sema_concepts, "sema_for_othergld", exposure_window=config$exposure_window, washout_concepts_df=othergld_concepts)
  othergld_idx_for_sema <- build_exposure_idx(othergld_concepts, "othergld_for_sema", exposure_window=config$exposure_window, washout_concepts_df=sema_concepts)
  comparison_results$`SEMAGLUTIDE vs OtherGLD` <- process_comparison_group(sema_idx_for_othergld, othergld_idx_for_sema, demographics_to_join, "SEMAGLUTIDE", "OtherGLD")
 
  # --- Comparison 4: OTHER_GLPA vs SGLT2 ---
  other_glpa_idx_for_sglt2 <- build_exposure_idx(other_glpa_concepts, "other_glpa_for_sglt2", exposure_window=config$exposure_window, washout_concepts_df=sglt2_concepts)
  sglt2_idx_for_other_glpa <- build_exposure_idx(sglt2_concepts, "sglt2_for_other_glpa", exposure_window=config$exposure_window, washout_concepts_df=other_glpa_concepts)
  comparison_results$`OTHER_GLPA vs SGLT2` <- process_comparison_group(other_glpa_idx_for_sglt2, sglt2_idx_for_other_glpa, demographics_to_join, "OTHER_GLPA", "SGLT2")
 
  # --- Comparison 5: OTHER_GLPA vs OtherGLD ---
  other_glpa_idx_for_othergld <- build_exposure_idx(other_glpa_concepts, "other_glpa_for_othergld", exposure_window=config$exposure_window, washout_concepts_df=othergld_concepts)
  othergld_idx_for_other_glpa <- build_exposure_idx(othergld_concepts, "othergld_for_other_glpa", exposure_window=config$exposure_window, washout_concepts_df=other_glpa_concepts)
  comparison_results$`OTHER_GLPA vs OtherGLD` <- process_comparison_group(other_glpa_idx_for_othergld, othergld_idx_for_other_glpa, demographics_to_join, "OTHER_GLPA", "OtherGLD")
 
  # --- Comparison 6: SGLT2 vs OtherGLD ---
  sglt2_idx_for_othergld <- build_exposure_idx(sglt2_concepts, "sglt2_for_othergld", exposure_window=config$exposure_window, washout_concepts_df=othergld_concepts)
  othergld_idx_for_sglt2 <- build_exposure_idx(othergld_concepts, "othergld_for_sglt2", exposure_window=config$exposure_window, washout_concepts_df=sglt2_concepts)
  comparison_results$`SGLT2 vs OtherGLD` <- process_comparison_group(sglt2_idx_for_othergld, othergld_idx_for_sglt2, demographics_to_join, "SGLT2", "OtherGLD")
 
  # Step 4: Batch Data Enrichment Loop
  for (comp_name in names(comparison_results)) {
    full_cohort_df <- comparison_results[[comp_name]]$cohort_df
    if (nrow(full_cohort_df) == 0) next
    batch_size <- 1000
    batches <- split(full_cohort_df, ceiling(seq_along(full_cohort_df$person_id) / batch_size))
    meds_list <- lapply(batches, add_medications, dataset = Sys.getenv("WORKSPACE_CDR"))
    combined_meds <- do.call(rbind, meds_list)
    batches2 <- split(combined_meds, ceiling(seq_along(combined_meds$person_id) / batch_size))
    com_list <- lapply(batches2, add_comorbidities, dataset = Sys.getenv("WORKSPACE_CDR"))
    combined_com <- do.call(rbind, com_list)
    combined_final <- add_baseline_measurements(combined_com, bmi_panel = bmi_panel, hba1c_panel = hba1c_panel)
    comparison_results[[comp_name]]$cohort_df <- combined_final
  }
   
  # Step 5: Final Cleaning to Remove Same-Day Initiators
  clean_comparison_results <- list()
  for (comp_name in names(comparison_results)) {
    final_df <- comparison_results[[comp_name]]$cohort_df
    if (is.null(final_df) || nrow(final_df) == 0) { clean_comparison_results[[comp_name]] <- list(cohort_df=final_df); next }
    initial_rows <- nrow(final_df)
    if (comp_name == "SEMAGLUTIDE vs OTHER_GLPA") {
      final_df <- final_df %>% filter(!((treatment == 1 & OTHER_GLPA == 1) | (treatment == 0 & SEMAGLUTIDE == 1)))
    } else if (comp_name == "SEMAGLUTIDE vs SGLT2") {
      final_df <- final_df %>% filter(!((treatment == 1 & SGLT2i == 1) | (treatment == 0 & SEMAGLUTIDE == 1)))
    } else if (comp_name == "SEMAGLUTIDE vs OtherGLD") {
      final_df <- final_df %>% filter(!((treatment == 1 & (TZD == 1 | DPP4i == 1 | SU == 1)) | (treatment == 0 & SEMAGLUTIDE == 1)))
    } else if (comp_name == "OTHER_GLPA vs SGLT2") {
      final_df <- final_df %>% filter(!((treatment == 1 & SGLT2i == 1) | (treatment == 0 & OTHER_GLPA == 1)))
    } else if (comp_name == "OTHER_GLPA vs OtherGLD") {
      final_df <- final_df %>% filter(!((treatment == 1 & (TZD == 1 | DPP4i == 1 | SU == 1)) | (treatment == 0 & OTHER_GLPA == 1)))
    } else if (comp_name == "SGLT2 vs OtherGLD") {
      final_df <- final_df %>% filter(!((treatment == 1 & (TZD == 1 | DPP4i == 1 | SU == 1)) | (treatment == 0 & SGLT2i == 1)))
    }
    final_rows <- nrow(final_df)
    log_output("INFO", paste("Cleaned cohort:", comp_name), paste("Removed", initial_rows - final_rows, "same-day initiators.", final_rows, "remain."))
    clean_comparison_results[[comp_name]] <- list(
      cohort_df = final_df,
      pre_exclusion_cohort_df = comparison_results[[comp_name]]$pre_exclusion_cohort_df
    )
  }
  # Step 7: Final Summary and Output
  results_df <- create_cohort_summary_df(cohort_stats, clean_comparison_results)
  output_file <- sprintf("semaglutide_epilepsy_seizure_cohort_analysis_results_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S"))
  write.csv(results_df, output_file, row.names = FALSE)
  save_output_log("semaglutide_epilepsy_seizure_cohort_analysis_log.csv")
  log_output("PROGRESS", "Cohort analysis complete", output_file, "")
 
  list(stats = cohort_stats, comparisons = clean_comparison_results, output_file = output_file)
}

# >>> CELL 09: Save & Run Cohort Building <<<
# ---------- SAVE FULL COHORTS TO BUCKET ----------
save_full_cohorts_to_bucket <- function(comparison_results, bucket = Sys.getenv("WORKSPACE_BUCKET")) {
 if (is.null(bucket) || bucket == "") stop("WORKSPACE_BUCKET not set")
 ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
 for (comp in names(comparison_results)) {
 cohort_df <- comparison_results[[comp]]$cohort_df
 if (is.null(cohort_df) || nrow(cohort_df) == 0) {
 log_output("WARNING", sprintf("Cohort %s empty; skipping save", comp), "", "save_full_cohorts_to_bucket")
 next
 }
 safe_comp <- gsub("[^A-Za-z0-9_\\-]", "_", comp)
 fname <- sprintf("%s_semaglutide_epilepsy_seizure_full_cohort_%s.csv", safe_comp, ts)
 readr::write_excel_csv(cohort_df, fname)
 system(paste0("gsutil cp ./", fname, " ", bucket, "/data/"), intern = TRUE)
 ls_out <- system(paste0("gsutil ls ", bucket, "/data/", fname), intern = TRUE)
 log_output("PROGRESS", sprintf("Saved cohort %s to bucket", comp), paste(ls_out, collapse = "; "), "save_full_cohorts_to_bucket")

 # Also save pre-exclusion cohort (includes patients WITH pre-index epilepsy)
 pre_excl_df <- comparison_results[[comp]]$pre_exclusion_cohort_df
 if (!is.null(pre_excl_df) && nrow(pre_excl_df) > 0) {
   fname_pre <- sprintf("%s_PRE_EXCLUSION_full_cohort_%s.csv", safe_comp, ts)
   readr::write_excel_csv(pre_excl_df, fname_pre)
   system(paste0("gsutil cp ./", fname_pre, " ", bucket, "/data/"), intern = TRUE)
   ls_out_pre <- system(paste0("gsutil ls ", bucket, "/data/", fname_pre), intern = TRUE)
   log_output("PROGRESS", sprintf("Saved PRE-EXCLUSION cohort %s to bucket", comp),
              paste(ls_out_pre, collapse = "; "), "save_full_cohorts_to_bucket")
 }
 }
 invisible(TRUE)
}
# ---------- RUN and SAVE ----------
# Run cohort building and save full cohorts for external IPWT
analysis_results <- perform_stepwise_cohort_analysis_optimized()
comparison_results <- analysis_results$comparisons
tryCatch({
 save_full_cohorts_to_bucket(comparison_results)
}, error = function(e) {
 log_output("ERROR", "Failed to save full cohorts to bucket", e$message, "main_analysis")
})


# >>> CELL 10: Load Cohort CSVs <<<
# --- PART 1: DATA LOADING ---
# Function to copy file from Google Bucket and load into a dataframe
load_from_bucket <- function(file_name) {
  if (!file.exists(file_name)) {
    stop(paste("File not found:", file_name,
               "Please ensure it's in the working directory."))
  }
  read_csv(file_name, show_col_types = FALSE)
}

# MODIFICATION: Use a fixed list of the 6 cohort files instead of a dynamic pattern.
files <- c(
  "SEMAGLUTIDE_vs_OTHER_GLPA_semaglutide_epilepsy_seizure_full_cohort_20260213_053619.csv",
  "SEMAGLUTIDE_vs_SGLT2_semaglutide_epilepsy_seizure_full_cohort_20260213_053619.csv",
  "SEMAGLUTIDE_vs_OtherGLD_semaglutide_epilepsy_seizure_full_cohort_20260213_053619.csv",
  "OTHER_GLPA_vs_SGLT2_semaglutide_epilepsy_seizure_full_cohort_20260213_053619.csv",
  "OTHER_GLPA_vs_OtherGLD_semaglutide_epilepsy_seizure_full_cohort_20260213_053619.csv",
  "SGLT2_vs_OtherGLD_semaglutide_epilepsy_seizure_full_cohort_20260213_053619.csv"
)

# Load each file into a named dataframe list
cohorts <- lapply(files, load_from_bucket)
names(cohorts) <- c(
  "semaglutide_vs_other_glpa_full_cohort",
  "semaglutide_vs_sglt2_full_cohort",
  "semaglutide_vs_othergld_full_cohort",
  "other_glpa_vs_sglt2_full_cohort",
  "other_glpa_vs_othergld_full_cohort",
  "sglt2_vs_othergld_full_cohort"
)

# >>> CELL 11: Baseline Table Generation <<<
# --- PART 2: BASELINE TABLE GENERATION ---
# Define the list of variables for the baseline table
my_vars <- c(
  "age", "sex_cat", "raceethnicity_cat", "baseline_bmi", "baseline_bmi_category", "baseline_hba1c",
  "index_year", "index_year_grouped", "income", "education", "smoking", "insurance_category",
  "alcohol_category", "SEMAGLUTIDE", "OTHER_GLPA", "Biguanide", "TZD", "Insulin", "SGLT2i",
  "DPP4i", "SU", "Anticoagulant", "Antiplatelet", "Statin", "Ezetimibe",
  "RAAS", "Diuretic", "MRA", "BB", "CCB", "OtherHTN",
  "myocardial_infarction", "congestive_heart_failure",
  "peripheral_vascular_disease", "cerebrovascular_disease", "dementia",
  "chronic_pulmonary_disease", "rheumatic_disease", "peptic_ulcer_disease", "hemiplegia_or_paraplegia",
  "hiv_infection", "hypoglycemia", "hyperglycemic_emergency",
  "renal_disease_severity", "liver_disease_severity", 
  "malignancy_status","diabetes_with_ophthalmic_complications", "diabetes_with_neurological_complications", "t1d", "t2d", "otherdm"
)

# Identify which of the above variables are categorical
categorical_vars <- my_vars[!my_vars %in% c("age", "baseline_bmi", "baseline_hba1c")]

# Function to apply labels to categorical variables
apply_variable_labels <- function(data) {
  data %>%
    mutate(
      # Sex category
      sex_cat = factor(case_when(
        sex_cat == 0 ~ "Male",
        sex_cat == 1 ~ "Female",
        sex_cat == 999 ~ "Missing",
        TRUE ~ "Missing"
      ), levels = c("Male", "Female", "Missing")),
      
      # Race/ethnicity category
      raceethnicity_cat = factor(case_when(
        raceethnicity_cat == 0 ~ "Non-Hispanic White",
        raceethnicity_cat == 1 ~ "Non-Hispanic Black",
        raceethnicity_cat == 2 ~ "Hispanic",
        raceethnicity_cat == 3 ~ "Other",
        raceethnicity_cat == 999 ~ "Missing",
        TRUE ~ "Missing"
      ), levels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other", "Missing")),
      
      # Income category
      income = factor(case_when(
        income == 0 ~ "Low (<50k)",
        income == 1 ~ "Middle (50k-100k)",
        income == 2 ~ "High (>100k)",
        income == 999 ~ "Missing",
        TRUE ~ "Missing"
      ), levels = c("Low (<50k)", "Middle (50k-100k)", "High (>100k)", "Missing")),
      
      # Education category
      education = factor(case_when(
        education == 0 ~ "High School or Less",
        education == 1 ~ "Some College",
        education == 2 ~ "Advanced Degree",
        education == 999 ~ "Missing",
        TRUE ~ "Missing"
      ), levels = c("High School or Less", "Some College", "Advanced Degree", "Missing")),
      
      # Smoking category
      smoking = factor(case_when(
        smoking == 0 ~ "Never",
        smoking == 1 ~ "Former",
        smoking == 2 ~ "Current",
        smoking == 999 ~ "Missing",
        TRUE ~ "Missing"
      ), levels = c("Never", "Former", "Current", "Missing")),
      
      # Insurance category
      insurance_category = factor(case_when(
        insurance_category == 0 ~ "None",
        insurance_category == 1 ~ "Public",
        insurance_category == 2 ~ "Private",
        insurance_category == 999 ~ "Missing",
        TRUE ~ "Missing"
      ), levels = c("None", "Public", "Private", "Missing")),
      
      # Alcohol category
      alcohol_category = factor(case_when(
        alcohol_category == 0 ~ "Low Risk",
        alcohol_category == 1 ~ "Increased Risk",
        alcohol_category == 2 ~ "High Risk",
        alcohol_category == 3 ~ "Dependent",
        alcohol_category == 999 ~ "Missing",
        TRUE ~ "Missing"
      ), levels = c("Low Risk", "Increased Risk", "High Risk", "Dependent", "Missing")),
 
        # Re-create baseline_bmi_category with all detailed obesity classes
        baseline_bmi_category = factor(
          case_when(
            is.na(baseline_bmi)                   ~ "Missing",
            baseline_bmi < 18.5                   ~ "Underweight",
            baseline_bmi >= 18.5 & baseline_bmi < 25 ~ "Normal",
            baseline_bmi >= 25 & baseline_bmi < 30   ~ "Overweight",
            baseline_bmi >= 30 & baseline_bmi < 35   ~ "Obesity Class I",
            baseline_bmi >= 35 & baseline_bmi < 40   ~ "Obesity Class II",
            baseline_bmi >= 40                    ~ "Obesity Class III"
          ),
          # This sets the correct logical order for the categories
          levels = c(
            "Underweight",
            "Normal",
            "Overweight",
            "Obesity Class I",
            "Obesity Class II",
            "Obesity Class III",
            "Missing"
          )
        ),

        # Convert index_year_grouped from character to factor
        index_year_grouped = factor(index_year_grouped),

      # Medication and comorbidity binary variables (0 = No, 1 = Yes)
across(c("SEMAGLUTIDE", "OTHER_GLPA", "Biguanide", "TZD", "Insulin", "SGLT2i", "DPP4i", "SU",
         "Anticoagulant", "Antiplatelet", "Statin", "Ezetimibe", "RAAS", "Diuretic",
         "MRA", "BB", "CCB", "OtherHTN", "myocardial_infarction",
         "congestive_heart_failure", "peripheral_vascular_disease",
         "cerebrovascular_disease", "dementia", "chronic_pulmonary_disease",
         "rheumatic_disease", "peptic_ulcer_disease", "hiv_infection",
                  "diabetes_with_ophthalmic_complications", "diabetes_with_neurological_complications", 
               "hypoglycemia", "hyperglycemic_emergency", "t1d", "t2d", "otherdm"),
       ~ factor(case_when(
         . == 0 ~ "No",
         . == 1 ~ "Yes",
         TRUE ~ "No"  # Treat NA or unexpected values as "No" instead of "Missing"
       ), levels = c("No", "Yes"))),  # Only include "No" and "Yes" levels
      
      # Renal Disease Severity
      renal_disease_severity = factor(case_when(
       renal_severe == 1 ~ "Severe",
        renal_mild_or_moderate == 1 ~ "Mild-Moderate",
        TRUE ~ "None"
      ), levels = c("None", "Mild-Moderate", "Severe")),
      
      # Liver Disease Severity
      liver_disease_severity = factor(case_when(
        moderate_or_severe_liver_disease == 1 ~ "Mod-Severe",
        mild_liver_disease == 1 ~ "Mild",
        TRUE ~ "None"
      ), levels = c("None", "Mild", "Mod-Severe")),
      
      # Malignancy Status
      malignancy_status = factor(case_when(
        metastatic_solid_tumor == 1 ~ "Metastatic",
        any_malignancy == 1 ~ "Non-Metastatic",
        TRUE ~ "None"
      ), levels = c("None", "Non-Metastatic", "Metastatic")),
      
# Individual diabetes complications as factors
        diabetes_with_ophthalmic_complications = as.factor(diabetes_with_ophthalmic_complications),
        diabetes_with_neurological_complications = as.factor(diabetes_with_neurological_complications))
}

# Loop through each cohort dataframe in the 'cohorts' list
for (cohort_name in names(cohorts)) {
  
  print(paste("Processing cohort:", cohort_name))
  
  # Get the current cohort's data from the list
  cohort_data <- cohorts[[cohort_name]]
  
  # Apply labels to categorical variables
  cohort_data <- apply_variable_labels(cohort_data)
  
  # Create the baseline table, stratified by the 'treatment' variable
  baseline_table_object <- CreateTableOne(
    vars = my_vars,
    strata = "treatment",
    data = cohort_data,
    factorVars = categorical_vars
  )
  
  # Convert the tableone object to a data frame for saving
  baseline_table_df <- as.data.frame(print(baseline_table_object,
                                           showAllLevels = TRUE,
                                           printToggle = FALSE,
                                           noSpaces = TRUE))
  
  # Get the current time and format it for the output filename
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  
  # Create the final output filename using the cohort's name
  output_filename <- paste0("epilepsy_seizure_baseline_characteristics_table_", cohort_name, "_", timestamp, ".csv")
  
  # Save the data frame to a new CSV file
  write.csv(baseline_table_df, file = output_filename, row.names = TRUE)
  
  # Print a confirmation message
  print(paste("Successfully saved baseline table as:", output_filename))
  cat("----------------------------------------------------\n")
}

# >>> CELL 12: Index Ingredient Tagging <<<
# ================== PART 3 (BOTH ARMS): TAG INDEX INGREDIENTS + MAKE TABLES ==================

# Normalize types after loading
cohorts <- lapply(cohorts, function(df) {
  df %>% dplyr::mutate(
    person_id = as.character(person_id),
    index_date = as.Date(index_date)
  )
})

# ---- A) Robust ingredient tagging on index_date for a SPECIFIC ARM ----
# - ingredient_names: character vector of ingredient NAMES to consider
# - arm_value: 1 for treatment, 0 for comparator
# - out_col: column to write (e.g., "index_ingredient_name_treat" or "_comp")
get_index_ingredient_for_arm <- function(cohort_df, ingredient_names, arm_value,
                                         out_col,
                                         dataset = Sys.getenv("WORKSPACE_CDR")) {
  # ensure types
  cohort_df <- cohort_df %>%
    dplyr::mutate(person_id = as.character(person_id),
                  index_date = as.Date(index_date))

  cdf <- cohort_df %>%
    dplyr::filter(.data$treatment == arm_value,
                  !is.na(person_id), !is.na(index_date))

  # initialize target column NOW so it's never NULL
  if (!out_col %in% names(cohort_df)) {
    cohort_df[[out_col]] <- NA_character_
  } else if (is.null(cohort_df[[out_col]])) {
    cohort_df[[out_col]] <- rep(NA_character_, nrow(cohort_df))
  }

  if (nrow(cdf) == 0) return(cohort_df)

  ing_names_sql <- paste0("'", tolower(ingredient_names), "'", collapse = ", ")
  cohort_rows <- purrr::map2_chr(cdf$person_id, cdf$index_date, ~
    glue::glue("STRUCT({.x} AS person_id, DATE('{format(.y, '%Y-%m-%d')}') AS index_date)")
  )
  cohort_cte_sql <- paste(cohort_rows, collapse = ", ")

  sql <- glue::glue("
    WITH cohort AS (
      SELECT * FROM UNNEST(ARRAY<STRUCT<person_id INT64, index_date DATE>>[{cohort_cte_sql}])
    ),
    idx_day_exposures AS (
      SELECT
        t.person_id,
        ca.ancestor_concept_id AS ingredient_id
      FROM cohort t
      JOIN `{dataset}.drug_exposure` de
        ON de.person_id = t.person_id
       AND de.drug_exposure_start_date = t.index_date
      JOIN `{dataset}.concept_ancestor` ca
        ON de.drug_concept_id = ca.descendant_concept_id
      JOIN `{dataset}.concept` anc
        ON anc.concept_id = ca.ancestor_concept_id
      WHERE anc.concept_class_id = 'Ingredient'
        AND LOWER(anc.concept_name) IN ({ing_names_sql})
    ),
    tie_broken AS (
      SELECT person_id, ingredient_id
      FROM (
        SELECT
          person_id,
          ingredient_id,
          ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY ingredient_id ASC) AS rn
        FROM idx_day_exposures
      )
      WHERE rn = 1
    )
    SELECT
      tb.person_id,
      c.concept_name AS index_ingredient_name
    FROM tie_broken tb
    JOIN `{dataset}.concept` c
      ON c.concept_id = tb.ingredient_id
  ")

  ix <- download_data(sql) %>%
    dplyr::mutate(person_id = as.character(person_id))

  out <- cohort_df %>%
    dplyr::left_join(ix, by = "person_id")

  # SAFER than ifelse(): assign only on rows for this arm
  rows <- which(out$treatment == arm_value)
  if ("index_ingredient_name" %in% names(out)) {
    out[[out_col]][rows] <- out$index_ingredient_name[rows]
    out$index_ingredient_name <- NULL
  } else {
    # no matches returned from SQL; leave as NA
  }

  out
}

# ---- B) Ingredient name sets (explicit, stable) ----
INGR_SEMA     <- c("semaglutide")
INGR_OTHERGLP <- c("exenatide","liraglutide","albiglutide","dulaglutide","lixisenatide")
INGR_SGLT2    <- c("canagliflozin","empagliflozin","dapagliflozin","ertugliflozin")
INGR_OTHERGLD <- c("glimepiride","glipizide","glyburide",
                   "alogliptin","linagliptin","sitagliptin","saxagliptin",
                   "pioglitazone","rosiglitazone")

# ---- C) Map each cohort to EXPOSURE (treatment arm) and COMPARATOR (comparator arm) sets ----
exposure_ing_map <- list(
  semaglutide_vs_other_glpa_full_cohort = INGR_SEMA,
  semaglutide_vs_sglt2_full_cohort      = INGR_SEMA,
  semaglutide_vs_othergld_full_cohort   = INGR_SEMA,
  other_glpa_vs_sglt2_full_cohort       = INGR_OTHERGLP,
  other_glpa_vs_othergld_full_cohort    = INGR_OTHERGLP,
  sglt2_vs_othergld_full_cohort         = INGR_SGLT2
)
comparator_ing_map <- list(
  semaglutide_vs_other_glpa_full_cohort = INGR_OTHERGLP,
  semaglutide_vs_sglt2_full_cohort      = INGR_SGLT2,
  semaglutide_vs_othergld_full_cohort   = INGR_OTHERGLD,
  other_glpa_vs_sglt2_full_cohort       = INGR_SGLT2,
  other_glpa_vs_othergld_full_cohort    = INGR_OTHERGLD,
  sglt2_vs_othergld_full_cohort         = INGR_OTHERGLD
)

# Helper: number of non-missing unique levels
n_levels <- function(x) length(unique(stats::na.omit(x)))

# ---- D) Build per-ingredient + class-combined tables for BOTH ARMS ----
for (cohort_name in names(cohorts)) {
  message("Tagging index ingredients & building tables for BOTH arms: ", cohort_name)
  cohort_df <- cohorts[[cohort_name]]

  # Tag treatment arm with its exposure class ingredients
  exp_names <- exposure_ing_map[[cohort_name]]
  if (!is.null(exp_names)) {
    cohort_df <- get_index_ingredient_for_arm(cohort_df, exp_names, arm_value = 1,
                                              out_col = "index_ingredient_name_treat")
  } else {
    cohort_df$index_ingredient_name_treat <- NA_character_
  }

  # Tag comparator arm with its comparator class ingredients
  comp_names <- comparator_ing_map[[cohort_name]]
  if (!is.null(comp_names)) {
    cohort_df <- get_index_ingredient_for_arm(cohort_df, comp_names, arm_value = 0,
                                              out_col = "index_ingredient_name_comp")
  } else {
    cohort_df$index_ingredient_name_comp <- NA_character_
  }

  # Apply your categorical labels
  cohort_df <- apply_variable_labels(cohort_df)

  # ---------- TREATMENT ARM ----------
  treat_df <- cohort_df %>% dplyr::filter(treatment == 1L)
  # Backstop for single-ingredient classes
  if (length(exp_names) == 1 && (all(is.na(treat_df$index_ingredient_name_treat)) || !"index_ingredient_name_treat" %in% names(treat_df))) {
    treat_df$index_ingredient_name_treat <- tolower(exp_names[[1]])
  }
  # Counts
  n_treat <- treat_df %>%
    dplyr::mutate(index_ingredient_name_treat = tolower(index_ingredient_name_treat)) %>%
    dplyr::count(index_ingredient_name_treat, name = "N") %>%
    dplyr::arrange(dplyr::desc(N))
  # All of Us small-cell suppression: mask ingredient counts < 20 per Data Dissemination Policy
  n_treat$N[n_treat$N < 20] <- NA_integer_
  n_treat$index_ingredient_name_treat[is.na(n_treat$N)] <- "<suppressed>"
  write.csv(n_treat,
            file = sprintf("N_by_ingredient_TREAT_%s_%s.csv", cohort_name, format(Sys.time(), "%Y%m%d_%H%M%S")),
            row.names = FALSE)

  # Per-ingredient TableOne for treatment (only if ≥2 levels)
  if (n_levels(treat_df$index_ingredient_name_treat) >= 2) {
    tbl_treat_ing <- CreateTableOne(
      vars       = my_vars,
      strata     = "index_ingredient_name_treat",
      data       = treat_df,
      factorVars = categorical_vars
    )
    df_treat_ing <- as.data.frame(print(tbl_treat_ing, showAllLevels = TRUE, printToggle = FALSE, noSpaces = TRUE))
    write.csv(df_treat_ing,
              file = sprintf("demographics_by_ingredient_TREAT_%s_%s.csv", cohort_name, format(Sys.time(), "%Y%m%d_%H%M%S")),
              row.names = TRUE)
  } else {
    message("Note: <2 ingredient levels in TREATMENT arm; skipping per-ingredient TableOne for ", cohort_name, ".")
  }

  # Class-combined (treatment)
  tbl_treat_cls <- CreateTableOne(
    vars       = my_vars,
    data       = treat_df,
    factorVars = categorical_vars
  )
  df_treat_cls <- as.data.frame(print(tbl_treat_cls, showAllLevels = TRUE, printToggle = FALSE, noSpaces = TRUE))
  write.csv(df_treat_cls,
            file = sprintf("demographics_by_class_TREAT_%s_%s.csv", cohort_name, format(Sys.time(), "%Y%m%d_%H%M%S")),
            row.names = TRUE)

  # ---------- COMPARATOR ARM ----------
  comp_df <- cohort_df %>% dplyr::filter(treatment == 0L)
  # Backstop for single-ingredient comparator classes (e.g., if comparator is semaglutide in some designs)
  if (!is.null(comp_names) && length(comp_names) == 1 &&
      (all(is.na(comp_df$index_ingredient_name_comp)) || !"index_ingredient_name_comp" %in% names(comp_df))) {
    comp_df$index_ingredient_name_comp <- tolower(comp_names[[1]])
  }
  # Counts
  n_comp <- comp_df %>%
    dplyr::mutate(index_ingredient_name_comp = tolower(index_ingredient_name_comp)) %>%
    dplyr::count(index_ingredient_name_comp, name = "N") %>%
    dplyr::arrange(dplyr::desc(N))
  # All of Us small-cell suppression: mask ingredient counts < 20 per Data Dissemination Policy
  n_comp$N[n_comp$N < 20] <- NA_integer_
  n_comp$index_ingredient_name_comp[is.na(n_comp$N)] <- "<suppressed>"
  write.csv(n_comp,
            file = sprintf("N_by_ingredient_COMP_%s_%s.csv", cohort_name, format(Sys.time(), "%Y%m%d_%H%M%S")),
            row.names = FALSE)

  # Per-ingredient TableOne for comparator (only if ≥2 levels)
  if (n_levels(comp_df$index_ingredient_name_comp) >= 2) {
    tbl_comp_ing <- CreateTableOne(
      vars       = my_vars,
      strata     = "index_ingredient_name_comp",
      data       = comp_df,
      factorVars = categorical_vars
    )
    df_comp_ing <- as.data.frame(print(tbl_comp_ing, showAllLevels = TRUE, printToggle = FALSE, noSpaces = TRUE))
    write.csv(df_comp_ing,
              file = sprintf("demographics_by_ingredient_COMP_%s_%s.csv", cohort_name, format(Sys.time(), "%Y%m%d_%H%M%S")),
              row.names = TRUE)
  } else {
    message("Note: <2 ingredient levels in COMPARATOR arm; skipping per-ingredient TableOne for ", cohort_name, ".")
  }

  # Class-combined (comparator)
  tbl_comp_cls <- CreateTableOne(
    vars       = my_vars,
    data       = comp_df,
    factorVars = categorical_vars
  )
  df_comp_cls <- as.data.frame(print(tbl_comp_cls, showAllLevels = TRUE, printToggle = FALSE, noSpaces = TRUE))
  write.csv(df_comp_cls,
            file = sprintf("demographics_by_class_COMP_%s_%s.csv", cohort_name, format(Sys.time(), "%Y%m%d_%H%M%S")),
            row.names = TRUE)

  message("✓ Saved TREATMENT and COMPARATOR tables for: ", cohort_name)
}

cohorts <- lapply(cohorts, function(df) {
  df %>%
    dplyr::mutate(
      person_id  = as.character(person_id),
      index_date = as.Date(index_date),
      # make sure treatment is 0/1 integer (not factor/char)
      treatment  = as.integer(as.character(treatment))
    )
})

# >>> CELL 13: IPTW Logging & Followup Functions <<<
# ============================================================================
# LOGGING FUNCTIONS
# ============================================================================
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
  if (is.null(getOption("semaglutide_output_tracker"))) {
    init_output_tracker()
  }
  tracker <- getOption("semaglutide_output_tracker")
  new_row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    script = script,
    type = type,
    message = message,
    value = as.character(value),
    stringsAsFactors = FALSE
  )
  options(semaglutide_output_tracker = rbind(tracker, new_row))
  if (type %in% c("PROGRESS", "ERROR", "WARNING")) {
    cat(sprintf("[%s] %s: %s\n", type, message, value))
  }
}

save_output_log <- function(filename = NULL) {
  if (is.null(filename)) {
    filename <- sprintf("output_log_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S"))
  }
  tracker <- getOption("semaglutide_output_tracker")
  if (!is.null(tracker) && nrow(tracker) > 0) {
    write.csv(tracker, filename, row.names = FALSE)
    log_output("PROGRESS", "Output log saved", filename)
    return(filename)
  }
  return(NULL)
}

# ============================================================================
# COHORT PROCESSING FUNCTIONS
# ============================================================================
followup_and_event <- function(df, outcome_var, data_cut_date,
                               late_onset = FALSE, early_onset = FALSE,
                               verbose = FALSE, comp_name = NULL,
                               exclude_preindex_any_age = FALSE) {
  result <- data.frame(df)
  result$age_at_event <- result$age + as.numeric(as.Date(result[[outcome_var]]) - as.Date(result$index_date)) / 365.25
  result$is_qualified_event <- (
    !is.na(result[[outcome_var]]) &
    if (late_onset) {
      result$age_at_event >= 60
    } else if (early_onset) {
      result$age_at_event < 60
    } else {
      TRUE
    }
  )

  # Detect and (for late_onset) exclude preindex/same-day raw outcomes
  raw_outcome_date <- try(as.Date(result[[outcome_var]]), silent = TRUE)
  idx_date_vec     <- try(as.Date(result$index_date), silent = TRUE)
  if (!inherits(raw_outcome_date, "try-error") && !inherits(idx_date_vec, "try-error")) {
    preidx_mask_detect <- !is.na(raw_outcome_date) & (raw_outcome_date <= idx_date_vec)
    log_output("PROGRESS",
               sprintf("Detected preindex/same-day outcomes (raw) [%s]",
                       ifelse(late_onset, "late_onset", "all_ages")),
               sum(preidx_mask_detect), "followup")
    if (isTRUE(late_onset)) {
      n_preidx <- sum(preidx_mask_detect)
      n_before <- nrow(result)
      if (n_preidx > 0) {
        try({
          out_name <- sprintf("preindex_outcomes_excluded_late_onset_%s_%s.csv",
                              outcome_var, format(Sys.time(), "%Y%m%d_%H%M%S"))
          safe_cols <- intersect(c("person_id","index_date", outcome_var, "age","treatment"), names(result))
          utils::write.csv(result[preidx_mask_detect, safe_cols, drop = FALSE], file = out_name, row.names = FALSE)
          log_output("PROGRESS", "Saved list of excluded raw preindex outcomes", out_name, "followup")
        }, silent = TRUE)
        result <- result[!preidx_mask_detect, , drop = FALSE]
        log_output("PROGRESS", "Late-onset preindex exclusion applied",
                   sprintf("n_preidx=%d; n_before=%d; n_after=%d", n_preidx, n_before, nrow(result)),
                   "followup")
      }
    }
  } else {
    log_output("WARNING", "Could not parse dates for preindex exclusion check", outcome_var, "followup")
  }

  result$outcome_date <- as.Date(ifelse(result$is_qualified_event, as.character(result[[outcome_var]]), NA), origin = "1970-01-01")
  result$censor_date <- pmin(as.Date(result$EHRmaxDT_min), as.Date(data_cut_date), na.rm = TRUE)
  result$end_fu_date <- pmin(result$outcome_date, result$censor_date, result$death_date, na.rm = TRUE)
  result$time <- as.numeric(result$end_fu_date - as.Date(result$index_date))
  result$event <- as.integer(!is.na(result$outcome_date) & result$outcome_date == result$end_fu_date)

  # Diagnostics for non-positive follow-up time
  n_pre_time <- nrow(result)
  n_nonpos   <- sum(!is.na(result$time) & result$time <= 0)
  if (!is.na(n_nonpos) && n_nonpos > 0) {
    diag_df <- result
    diag_df$analysis_group <- ifelse(late_onset, "late_onset", ifelse(early_onset, "early_onset", "all_ages"))
    diag_df$time_nonpos <- !is.na(diag_df$time) & diag_df$time <= 0
    diag_df$reason_nonpos <- NA_character_
    diag_df$had_outcome_raw <- !is.na(diag_df[[outcome_var]])
    diag_df$age_ok <- if (late_onset) { diag_df$age_at_event >= 60 } else if (early_onset) { diag_df$age_at_event < 60 } else { TRUE }
    diag_df$end_is_outcome <- !is.na(diag_df$outcome_date) & (diag_df$end_fu_date == diag_df$outcome_date)
    diag_df$end_is_censor  <- (!is.na(diag_df$censor_date)) & (diag_df$end_fu_date == diag_df$censor_date)
    diag_df$end_is_death   <- (!is.na(diag_df$death_date))  & (diag_df$end_fu_date == diag_df$death_date)
    diag_df$index_ge_end   <- as.Date(diag_df$index_date) >= diag_df$end_fu_date
    diag_df$index_gt_end   <- as.Date(diag_df$index_date) >  diag_df$end_fu_date

    diag_df$reason_nonpos[diag_df$time_nonpos & diag_df$end_is_outcome] <- "preindex_or_same_day_outcome"
    diag_df$reason_nonpos[diag_df$time_nonpos & diag_df$end_is_censor]  <- ifelse(is.na(diag_df$reason_nonpos[diag_df$time_nonpos & diag_df$end_is_censor]), "preindex_or_same_day_censor", diag_df$reason_nonpos[diag_df$time_nonpos & diag_df$end_is_censor])
    diag_df$reason_nonpos[diag_df$time_nonpos & diag_df$end_is_death]   <- ifelse(is.na(diag_df$reason_nonpos[diag_df$time_nonpos & diag_df$end_is_death]), "preindex_or_same_day_death", diag_df$reason_nonpos[diag_df$time_nonpos & diag_df$end_is_death])
    diag_df$reason_nonpos[diag_df$time_nonpos & diag_df$had_outcome_raw & !diag_df$age_ok] <- "event_failed_age_criterion"
    diag_df$reason_nonpos[diag_df$time_nonpos & is.na(diag_df$reason_nonpos) & !diag_df$had_outcome_raw] <- "no_event_recorded"
    diag_df$reason_nonpos[diag_df$time_nonpos & is.na(diag_df$reason_nonpos)] <- "other_preindex_issue"

    by_reason <- diag_df %>% dplyr::filter(time_nonpos) %>% dplyr::count(reason_nonpos, name = "n")
    by_treatment <- diag_df %>% dplyr::filter(time_nonpos) %>% dplyr::count(treatment, name = "n")
    msg_reason <- paste(paste0(by_reason$reason_nonpos, "=", by_reason$n), collapse = ", ")
    msg_trt    <- paste(paste0("trt", by_treatment$treatment, "=", by_treatment$n), collapse = ", ")
    log_output("PROGRESS", "Non-positive follow-up filtered",
               sprintf("Removed %d of %d; reasons: %s; by treatment: %s",
                       n_nonpos, n_pre_time, msg_reason, msg_trt),
               "followup")

    out_name <- sprintf("nonpos_followup_%s_%s_%s.csv",
                        diag_df$analysis_group[1], outcome_var, format(Sys.time(), "%Y%m%d_%H%M%S"))
    safe_cols <- intersect(c("person_id","index_date","outcome_date","censor_date","death_date","time","event","treatment","age","age_at_event","is_qualified_event","analysis_group","end_is_outcome","end_is_censor","end_is_death","index_ge_end","reason_nonpos"), names(diag_df))
    try({
      utils::write.csv(diag_df[diag_df$time_nonpos, safe_cols, drop = FALSE], file = out_name, row.names = FALSE)
      log_output("PROGRESS", "Saved non-positive follow-up details", out_name, "followup")
    }, silent = TRUE)

    if (is.null(getOption("semaglutide_nonpos_summary"))) {
      options(semaglutide_nonpos_summary = data.frame(
        analysis_group = character(),
        comparison     = character(),
        outcome        = character(),
        reason_nonpos  = character(),
        n              = integer(),
        stringsAsFactors = FALSE
      ))
    }
    s <- getOption("semaglutide_nonpos_summary")
    comp_val <- ifelse(is.null(comp_name), "(unknown)", comp_name)
    add_rows <- data.frame(
      analysis_group = diag_df$analysis_group[1],
      comparison     = comp_val,
      outcome        = outcome_var,
      reason_nonpos  = by_reason$reason_nonpos,
      n              = by_reason$n,
      stringsAsFactors = FALSE
    )
    options(semaglutide_nonpos_summary = rbind(s, add_rows))

    if (is.null(getOption("semaglutide_nonpos_people"))) {
      options(semaglutide_nonpos_people = data.frame(
        person_id      = character(),
        analysis_group = character(),
        comparison     = character(),
        outcome        = character(),
        reason_nonpos  = character(),
        stringsAsFactors = FALSE
      ))
    }
    p <- getOption("semaglutide_nonpos_people")
    idx <- which(diag_df$time_nonpos)
    if (length(idx) > 0) {
      ppl_rows <- data.frame(
        person_id      = as.character(diag_df$person_id[idx]),
        analysis_group = diag_df$analysis_group[idx],
        comparison     = rep(comp_val, length(idx)),
        outcome        = rep(outcome_var, length(idx)),
        reason_nonpos  = diag_df$reason_nonpos[idx],
        stringsAsFactors = FALSE
      )
      options(semaglutide_nonpos_people = rbind(p, ppl_rows))
    }
  }

  result <- result %>% dplyr::filter(time > 0)
  result$event_time <- result$time
  log_output("INFO", "followup_and_event complete",
             sprintf("%d events from %d patients", sum(result$event), nrow(result)),
             "followup_and_event")
  return(result)
}

# >>> CELL 14: IPTW Core Functions <<<
# ============================================================================
# IPTW ANALYSIS FUNCTIONS
# ============================================================================
run_ipwt_and_cox <- function(df, exclude_vars = NULL, trim_threshold = 0.01,
                             all_ps_vars = NULL, verbose = NULL, stabilize = FALSE) {
  all_covariates <- if (is.null(all_ps_vars)) {
    get_standard_config()$all_ps_vars
  } else {
    all_ps_vars
  }
  rhs_vars <- setdiff(all_covariates, exclude_vars)
  keep_vars <- rhs_vars[sapply(df[rhs_vars], function(x) length(unique(na.omit(x))) > 1)]
  ps_form <- reformulate(keep_vars, response = "treatment")
  model_vars <- c("treatment", keep_vars)
  # Diagnostic: NA counts per model variable
  if (isTRUE(verbose)) {
    na_counts <- sapply(df[model_vars], function(x) sum(is.na(x)))
    if (any(na_counts > 0)) {
      na_msg <- paste(
        paste0(names(na_counts)[na_counts > 0], "=", na_counts[na_counts > 0]),
        collapse = ", "
      )
      log_output("PROGRESS", "Missing PS covariates before complete-case filter", na_msg, "ipwt")
    }
  }
  complete_rows <- complete.cases(df[model_vars])
  df_complete <- df[complete_rows, ]
  if (verbose) {
    n_removed <- nrow(df) - nrow(df_complete)
    log_output("INFO", "IPWT sample size", sprintf("%d complete cases from %d total",
                                                   nrow(df_complete), nrow(df)))
    log_output("PROGRESS", "Complete-case filtering removed rows",
               sprintf("%d rows removed due to missing PS covariates", n_removed),
               "ipwt")
  }
  ps_model <- glm(ps_form, data = df_complete, family = binomial())
  df_complete$ps <- predict(ps_model, type = "response")
  if (isTRUE(stabilize)) {
    # Stabilized IPTW: multiply by marginal treatment probability
    p_treat <- mean(df_complete$treatment, na.rm = TRUE)
    if (verbose) {
      ps_min <- min(df_complete$ps, na.rm = TRUE)
      ps_max <- max(df_complete$ps, na.rm = TRUE)
      log_output("PROGRESS", "PS summary", sprintf("stabilize=TRUE, p_treat=%.4f, ps_min=%.4f, ps_max=%.4f", p_treat, ps_min, ps_max), "ipwt")
    }
    df_complete$ipw <- ifelse(df_complete$treatment == 1,
                              p_treat/df_complete$ps,
                              (1 - p_treat)/(1 - df_complete$ps))
  } else {
    # Unstabilized IPTW
    if (verbose) {
      ps_min <- min(df_complete$ps, na.rm = TRUE)
      ps_max <- max(df_complete$ps, na.rm = TRUE)
      log_output("PROGRESS", "PS summary", sprintf("stabilize=FALSE, ps_min=%.4f, ps_max=%.4f", ps_min, ps_max), "ipwt")
    }
    df_complete$ipw <- ifelse(df_complete$treatment == 1,
                              1/df_complete$ps,
                              1/(1 - df_complete$ps))
  }
  if (!is.null(trim_threshold)) {
    lower <- quantile(df_complete$ipw, trim_threshold)
    upper <- quantile(df_complete$ipw, 1 - trim_threshold)
    n_trim_low  <- sum(df_complete$ipw < lower, na.rm = TRUE)
    n_trim_high <- sum(df_complete$ipw > upper, na.rm = TRUE)
    df_complete$ipw_trimmed <- pmin(pmax(df_complete$ipw, lower), upper)
    if (verbose) {
      log_output("PROGRESS", "IPTW trimming",
                 sprintf("lower=%.4f, upper=%.4f; trimmed low=%d, high=%d of %d",
                         lower, upper, n_trim_low, n_trim_high, nrow(df_complete)),
                 "ipwt")
    }
  } else {
    df_complete$ipw_trimmed <- df_complete$ipw
  }
  df_complete$ipw_std <- df_complete$ipw_trimmed * nrow(df_complete) / sum(df_complete$ipw_trimmed)
  if (verbose) {
    w <- df_complete$ipw_std
    ess_total <- (sum(w)^2) / sum(w^2)
    w_treat <- w[df_complete$treatment == 1]
    w_ctrl  <- w[df_complete$treatment == 0]
    ess_treat <- (sum(w_treat)^2) / sum(w_treat^2)
    ess_ctrl  <- (sum(w_ctrl)^2) / sum(w_ctrl^2)
    pseudo_treat  <- sum(w_treat)
    pseudo_control<- sum(w_ctrl)
    log_output("PROGRESS", "Weights diagnostics",
               sprintf("ESS total=%.1f, treat=%.1f, control=%.1f; pseudo pop: treat=%.1f, control=%.1f, total=%.1f",
                       ess_total, ess_treat, ess_ctrl, pseudo_treat, pseudo_control, pseudo_treat + pseudo_control),
               "ipwt")
  }
  balance_before <- calculate_balance_metrics(df_complete, keep_vars, NULL, verbose = FALSE)
  balance_after <- calculate_balance_metrics(df_complete, keep_vars, df_complete$ipw_std, verbose = FALSE)
  weighted_cox <- coxph(Surv(event_time, event) ~ treatment,
                        data = df_complete, weights = ipw_std)
  return(list(
    balance_before = balance_before,
    balance_after = balance_after,
    cox = summary(weighted_cox),
    cohort = df_complete,
    ps_model = ps_model,
    dropped_covariates = setdiff(rhs_vars, keep_vars),
    n_dropped = nrow(df) - nrow(df_complete)
  ))
}

calculate_balance_metrics <- function(df, vars, weights, verbose = FALSE) {
  balance_df <- data.frame(
    variable = character(),
    std_diff = numeric(),
    stringsAsFactors = FALSE
  )
  multilevel_overall_vars <- c("sex_cat", "raceethnicity_cat",
                               "income", "education", "insurance_category", #"baseline_bmi_category",
                               "smoking", "alcohol_category", "index_year_grouped",
                               "renal_disease_severity", "liver_disease_severity",
                               "malignancy_status")
  for (var in vars) {
    if (!var %in% names(df)) next
    if (is.factor(df[[var]]) && var %in% multilevel_overall_vars) {
      all_levels <- levels(df[[var]])
      level_smds <- numeric()
      for (lvl in all_levels) {
        var_binary <- as.numeric(df[[var]] == lvl)
        if (is.null(weights)) {
          mean_t1 <- mean(var_binary[df$treatment == 1], na.rm = TRUE)
          mean_t0 <- mean(var_binary[df$treatment == 0], na.rm = TRUE)
          var_t1 <- var(var_binary[df$treatment == 1], na.rm = TRUE)
          var_t0 <- var(var_binary[df$treatment == 0], na.rm = TRUE)
        } else {
          w_t1 <- weights[df$treatment == 1]
          w_t0 <- weights[df$treatment == 0]
          mean_t1 <- weighted.mean(var_binary[df$treatment == 1], w_t1, na.rm = TRUE)
          mean_t0 <- weighted.mean(var_binary[df$treatment == 0], w_t0, na.rm = TRUE)
          var_t1 <- sum(w_t1 * (var_binary[df$treatment == 1] - mean_t1)^2, na.rm = TRUE) / sum(w_t1)
          var_t0 <- sum(w_t0 * (var_binary[df$treatment == 0] - mean_t0)^2, na.rm = TRUE) / sum(w_t0)
        }
        pooled_sd <- sqrt((var_t1 + var_t0) / 2)
        level_smd <- if(pooled_sd > 0) (mean_t1 - mean_t0) / pooled_sd else 0
        level_smds <- c(level_smds, level_smd)
      }
      overall_smd <- ifelse(length(level_smds) > 0, max(abs(level_smds), na.rm = TRUE), 0)
      if(length(level_smds) > 0) {
        max_idx <- which.max(abs(level_smds))
        overall_smd <- overall_smd * sign(level_smds[max_idx])
      }
      balance_df <- rbind(balance_df, data.frame(variable = var, std_diff = overall_smd))
    } else if (is.factor(df[[var]])) {
      for (lvl in levels(df[[var]])[-1]) {
        var_binary <- as.numeric(df[[var]] == lvl)
        if (is.null(weights)) {
          mean_t1 <- mean(var_binary[df$treatment == 1], na.rm = TRUE)
          mean_t0 <- mean(var_binary[df$treatment == 0], na.rm = TRUE)
          var_t1 <- var(var_binary[df$treatment == 1], na.rm = TRUE)
          var_t0 <- var(var_binary[df$treatment == 0], na.rm = TRUE)
        } else {
          w_t1 <- weights[df$treatment == 1]
          w_t0 <- weights[df$treatment == 0]
          mean_t1 <- weighted.mean(var_binary[df$treatment == 1], w_t1, na.rm = TRUE)
          mean_t0 <- weighted.mean(var_binary[df$treatment == 0], w_t0, na.rm = TRUE)
          var_t1 <- sum(w_t1 * (var_binary[df$treatment == 1] - mean_t1)^2, na.rm = TRUE) / sum(w_t1)
          var_t0 <- sum(w_t0 * (var_binary[df$treatment == 0] - mean_t0)^2, na.rm = TRUE) / sum(w_t0)
        }
        pooled_sd <- sqrt((var_t1 + var_t0) / 2)
        std_diff <- if(pooled_sd > 0) (mean_t1 - mean_t0) / pooled_sd else 0
        balance_df <- rbind(balance_df, data.frame(variable = paste0(var, lvl), std_diff = std_diff))
      }
    } else {
      if (is.null(weights)) {
        mean_t1 <- mean(df[[var]][df$treatment == 1], na.rm = TRUE)
        mean_t0 <- mean(df[[var]][df$treatment == 0], na.rm = TRUE)
        var_t1 <- var(df[[var]][df$treatment == 1], na.rm = TRUE)
        var_t0 <- var(df[[var]][df$treatment == 0], na.rm = TRUE)
      } else {
        w_t1 <- weights[df$treatment == 1]
        w_t0 <- weights[df$treatment == 0]
        mean_t1 <- weighted.mean(df[[var]][df$treatment == 1], w_t1, na.rm = TRUE)
        mean_t0 <- weighted.mean(df[[var]][df$treatment == 0], w_t0, na.rm = TRUE)
        var_t1 <- sum(w_t1 * (df[[var]][df$treatment == 1] - mean_t1)^2, na.rm = TRUE) / sum(w_t1)
        var_t0 <- sum(w_t0 * (df[[var]][df$treatment == 0] - mean_t0)^2, na.rm = TRUE) / sum(w_t0)
      }
      pooled_sd <- sqrt((var_t1 + var_t0) / 2)
      std_diff <- if(pooled_sd > 0) (mean_t1 - mean_t0) / pooled_sd else 0
      balance_df <- rbind(balance_df, data.frame(variable = var, std_diff = std_diff))
    }
  }
  return(balance_df)
}

# >>> CELL 15: Results Export & Visualization <<<
# ============================================================================
# RESULTS EXPORT FUNCTIONS
# ============================================================================
save_analysis_results <- function(results_list, prefix = "analysis") {
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  all_results <- list()
  for (name in names(results_list)) {
    res <- results_list[[name]]
    if (is.null(res) || is.null(res$cox) || is.null(res$crude_hr_results) || is.null(res$incident_rate_results)) {
      log_output("WARNING", "Skipping results saving due to missing data", name, "save_analysis_results")
      next
    }
    split_name <- strsplit(name, " - ")[[1]]
    comparison_name <- split_name[1]
    outcome_name <- split_name[2]
    unweighted_stats <- res$crude_hr_results
    unweighted_cohort <- res$cohort
    total_fu_years_unweighted <- unweighted_cohort$time / 365.25
    total_fu_mean_unweighted <- mean(total_fu_years_unweighted)
    total_fu_sd_unweighted <- sd(total_fu_years_unweighted)
    weighted_cohort <- res$cohort
    inc_rates <- res$incident_rate_results
    fu_stats_weighted <- calculate_followup_stats(weighted_cohort, weight_var = "ipw_std")
    pseudo_pop_treatment <- sum(weighted_cohort$ipw_std[weighted_cohort$treatment == 1])
    pseudo_pop_control <- sum(weighted_cohort$ipw_std[weighted_cohort$treatment == 0])
    pseudo_pop_total <- pseudo_pop_treatment + pseudo_pop_control
    # Effective sample sizes (ESS): overall and by treatment arm
    w_all <- weighted_cohort$ipw_std
    ess_total <- (sum(w_all)^2) / sum(w_all^2)
    w_treat <- w_all[weighted_cohort$treatment == 1]
    w_ctrl  <- w_all[weighted_cohort$treatment == 0]
    ess_treat <- (sum(w_treat)^2) / sum(w_treat^2)
    ess_ctrl  <- (sum(w_ctrl)^2) / sum(w_ctrl^2)
    weighted_events_treatment <- sum(weighted_cohort$event[weighted_cohort$treatment == 1] * weighted_cohort$ipw_std[weighted_cohort$treatment == 1])
    weighted_events_control <- sum(weighted_cohort$event[weighted_cohort$treatment == 0] * weighted_cohort$ipw_std[weighted_cohort$treatment == 0])
    weighted_events_total <- weighted_events_treatment + weighted_events_control
    total_py_weighted <- inc_rates$treatment_person_years + inc_rates$control_person_years
    total_ir_weighted <- (weighted_events_total / total_py_weighted) * 1000
    cox_model <- res$cox
    results_df_row <- data.frame(
      comparison = comparison_name,
      outcome = outcome_name,
      unweighted_patients_total = sprintf("%d", unweighted_stats$n_total),
      unweighted_events_total = sprintf("%d", unweighted_stats$n_events),
      unweighted_events_patients_treatment = sprintf("%d / %d", unweighted_stats$events_treatment, unweighted_stats$n_treatment),
      unweighted_events_patients_control = sprintf("%d / %d", unweighted_stats$events_control, unweighted_stats$n_control),
      unweighted_fu_years_total_mean_sd = sprintf("%.2f (%.2f)", total_fu_mean_unweighted, total_fu_sd_unweighted),
      unweighted_fu_years_treatment_mean_sd = sprintf("%.2f (%.2f)", unweighted_stats$followup_years_treatment_mean, unweighted_stats$followup_years_treatment_sd),
      unweighted_fu_years_control_mean_sd = sprintf("%.2f (%.2f)", unweighted_stats$followup_years_control_mean, unweighted_stats$followup_years_control_sd),
      weighted_pseudo_pop_total = sprintf("%.1f", pseudo_pop_total),
      weighted_pseudo_pop_treatment = sprintf("%.1f", pseudo_pop_treatment),
      weighted_pseudo_pop_control = sprintf("%.1f", pseudo_pop_control),
      weighted_pseudo_events_total = sprintf("%.1f", weighted_events_total),
      weighted_pseudo_events_treatment = sprintf("%.1f", weighted_events_treatment),
      weighted_pseudo_events_control = sprintf("%.1f", weighted_events_control),
      weighted_fu_years_treatment_mean_sd = sprintf("%.2f (%.2f)", fu_stats_weighted$treatment_mean, fu_stats_weighted$treatment_sd),
      weighted_fu_years_control_mean_sd = sprintf("%.2f (%.2f)", fu_stats_weighted$control_mean, fu_stats_weighted$control_sd),
      ir_1000py_total = sprintf("%.2f", total_ir_weighted),
      ir_1000py_treatment_ci = sprintf("%.2f (%.2f to %.2f)", inc_rates$treatment_rate, inc_rates$treatment_rate_ci[1], inc_rates$treatment_rate_ci[2]),
      ir_1000py_control_ci = sprintf("%.2f (%.2f to %.2f)", inc_rates$control_rate, inc_rates$control_rate_ci[1], inc_rates$control_rate_ci[2]),
      rate_diff_1000py_ci = sprintf("%.2f (%.2f to %.2f)", inc_rates$rate_difference, inc_rates$rate_difference_ci[1], inc_rates$rate_difference_ci[2]),
      hr_crude_ci = sprintf("%.2f (%.2f to %.2f)", unweighted_stats$hr, unweighted_stats$hr_lower, unweighted_stats$hr_upper),
      hr_crude_p_value = unweighted_stats$p_value,
      hr_weighted_ci = sprintf("%.2f (%.2f to %.2f)", cox_model$conf.int[1, "exp(coef)"], cox_model$conf.int[1, "lower .95"], cox_model$conf.int[1, "upper .95"]),
      hr_p_value = cox_model$coefficients[1, "Pr(>|z|)"]
    )
    all_results[[name]] <- results_df_row
  }
  if (length(all_results) > 0) {
    final_results_df <- do.call(rbind, all_results)
    results_file <- sprintf("%s_comprehensive_results_%s.csv", prefix, timestamp)
    # Log columns present to verify ESS fields are included
    log_output("PROGRESS", "Comprehensive results columns",
               paste(colnames(final_results_df), collapse = ", "),
               "save_analysis_results")
    write.csv(final_results_df, results_file, row.names = FALSE)
    log_output("PROGRESS", "Comprehensive results saved", results_file)
  } else {
    log_output("WARNING", "No results were generated to save.", "", "save_analysis_results")
    results_file <- NULL
  }
  balance_list <- lapply(names(results_list), function(name) {
    res <- results_list[[name]]
    if (!is.null(res)) {
      before <- res$balance_before; after <- res$balance_after
      before$analysis <- name; before$weighting <- "before"
      after$analysis <- name; after$weighting <- "after"
      rbind(before, after)
    }
  })
  if (length(balance_list) > 0 && !all(sapply(balance_list, is.null))) {
    balance_df <- do.call(rbind, balance_list)
    balance_file <- sprintf("%s_balance_%s.csv", prefix, timestamp)
    write.csv(balance_df, balance_file, row.names = FALSE)
    log_output("PROGRESS", "Balance metrics saved", balance_file)
  } else {
    balance_file <- NULL
  }
  return(list(results = results_file, balance = balance_file))
}

# ============================================================================
# VISUALIZATION FUNCTIONS
# ============================================================================
create_balance_plot <- function(balance_before, balance_after, title = "SMD Comparison") {
  balance_data <- rbind(
    balance_before %>% mutate(type = "Before"),
    balance_after %>% mutate(type = "After")
  )
  p <- ggplot(balance_data, aes(x = abs(std_diff), y = reorder(variable, abs(std_diff)),
                                color = type)) +
    geom_point(size = 2, alpha = 0.7) +
    geom_vline(xintercept = 0.1, linetype = "dashed", color = "red") +
    scale_color_manual(values = c("Before" = "#E74C3C", "After" = "#27AE60")) +
    labs(x = "Absolute Standardized Mean Difference",
         y = "Variable",
         title = title,
         color = "Weighting") +
    theme_minimal() +
    theme(legend.position = "bottom")
  return(p)
}

# ============================================================================
# SURVIVAL PLOT FUNCTIONS
# ============================================================================
create_survival_plot <- function(ipwt_result, comparison_name, outcome_label, analysis_name, save_plot = TRUE) {
  if (!requireNamespace("survminer", quietly = TRUE)) {
    install.packages("survminer")
  }
  library(survminer)
  library(ggplot2)
  weighted_df <- ipwt_result$cohort
  if (comparison_name == "SEMAGLUTIDE vs OTHER_GLPA") {
    treatment_labels <- c("Other GLP-1 Agonists", "Semaglutide")
  } else if (comparison_name == "SEMAGLUTIDE vs SGLT2") {
    treatment_labels <- c("SGLT-2 Inhibitors", "Semaglutide")
  } else if (comparison_name == "SEMAGLUTIDE vs OtherGLD") {
    treatment_labels <- c("Other GLDs", "Semaglutide")
  } else if (comparison_name == "OTHER_GLPA vs SGLT2") {
    treatment_labels <- c("SGLT-2 Inhibitors", "Other GLP-1 Agonists")
  } else if (comparison_name == "OTHER_GLPA vs OtherGLD") {
    treatment_labels <- c("Other GLDs", "Other GLP-1 Agonists")
  } else if (comparison_name == "SGLT2 vs OtherGLD") {
    treatment_labels <- c("Other GLDs", "SGLT-2 Inhibitors")
  } else {
    treatment_labels <- c("Comparator", "Treatment")
  }
  weighted_df$event_time_months <- weighted_df$event_time / 30.44
  cum_incidence_fit <- survfit(Surv(event_time_months, event) ~ treatment,
                               data = weighted_df,
                               weights = weighted_df$ipw_std)
  surv_plot <- ggsurvplot(
    fit = cum_incidence_fit,
    data = weighted_df,
    fun = "event",
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.fontsize = 3,
    tables.height = 0.3,
    legend.title = "Treatment Group",
    legend.labs = treatment_labels,
    xlab = "Follow-up Time (months)",
    ylab = "Cumulative Incidence",
    title = paste("Cumulative Incidence of", outcome_label, "for", analysis_name,
                  "\nStratified by Treatment Group (IPTW-weighted)"),
    break.time.by = 12,
    xlim = c(0, 48),
    lwd = 1.0,
    censor = FALSE,
    palette = "jama",
    ggtheme = theme_bw() +
      theme(panel.grid.major.x = element_blank(),
            panel.grid.minor.x = element_blank())
  )
  if (save_plot) {
    timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
    plot_filename <- paste0("ipwt_survival_",
                            analysis_name, "_",
                            gsub(" ", "_", comparison_name), "_",
                            gsub("[/ ]", "_", outcome_label), "_",
                            timestamp, ".png")
    tryCatch({
      png(filename = plot_filename, width = 12, height = 10, units = "in", res = 300)
      print(surv_plot)
      dev.off()
      log_output("INFO", "Survival plot saved", plot_filename)
    }, error = function(e) {
      log_output("ERROR", "Failed to save survival plot", e$message)
    })
  }
  return(surv_plot)
}

# >>> CELL 16: Statistics Calculators <<<
# ============================================================================
# COMPREHENSIVE STATISTICS CALCULATORS
# ============================================================================
calculate_incident_rates <- function(cohort_df, weight_var = NULL) {
  if (!is.null(weight_var) && weight_var %in% names(cohort_df)) {
    treatment_df <- cohort_df %>% filter(treatment == 1)
    control_df <- cohort_df %>% filter(treatment == 0)
    treatment_py <- sum(treatment_df$time * treatment_df[[weight_var]] / 365.25)
    control_py <- sum(control_df$time * control_df[[weight_var]] / 365.25)
    treatment_events <- sum(treatment_df$event * treatment_df[[weight_var]])
    control_events <- sum(control_df$event * control_df[[weight_var]])
  } else {
    treatment_df <- cohort_df %>% filter(treatment == 1)
    control_df <- cohort_df %>% filter(treatment == 0)
    treatment_py <- sum(treatment_df$time) / 365.25
    control_py <- sum(control_df$time) / 365.25
    treatment_events <- sum(treatment_df$event)
    control_events <- sum(control_df$event)
  }
  treatment_rate <- (treatment_events / treatment_py) * 1000
  control_rate <- (control_events / control_py) * 1000
  rate_diff <- treatment_rate - control_rate
  treatment_rate_se <- sqrt(treatment_events) / treatment_py * 1000
  control_rate_se <- sqrt(control_events) / control_py * 1000
  treatment_rate_lower <- treatment_rate - 1.96 * treatment_rate_se
  treatment_rate_upper <- treatment_rate + 1.96 * treatment_rate_se
  control_rate_lower <- control_rate - 1.96 * control_rate_se
  control_rate_upper <- control_rate + 1.96 * control_rate_se
  rate_diff_se <- sqrt((treatment_events / treatment_py^2) + (control_events / control_py^2)) * 1000
  rate_diff_lower <- rate_diff - 1.96 * rate_diff_se
  rate_diff_upper <- rate_diff + 1.96 * rate_diff_se
  return(list(
    treatment_rate = treatment_rate,
    treatment_rate_ci = c(treatment_rate_lower, treatment_rate_upper),
    control_rate = control_rate,
    control_rate_ci = c(control_rate_lower, control_rate_upper),
    rate_difference = rate_diff,
    rate_difference_ci = c(rate_diff_lower, rate_diff_upper),
    treatment_person_years = treatment_py,
    control_person_years = control_py
  ))
}

calculate_followup_stats <- function(cohort_df, weight_var = NULL) {
  cohort_df$followup_years <- cohort_df$time / 365.25
  if (!is.null(weight_var) && weight_var %in% names(cohort_df)) {
    treatment_df <- cohort_df %>% filter(treatment == 1)
    control_df <- cohort_df %>% filter(treatment == 0)
    treatment_mean <- weighted.mean(treatment_df$followup_years, treatment_df[[weight_var]])
    control_mean <- weighted.mean(control_df$followup_years, control_df[[weight_var]])
    treatment_sd <- sqrt(sum(treatment_df[[weight_var]] * (treatment_df$followup_years - treatment_mean)^2) /
                           sum(treatment_df[[weight_var]]))
    control_sd <- sqrt(sum(control_df[[weight_var]] * (control_df$followup_years - control_mean)^2) /
                         sum(control_df[[weight_var]]))
  } else {
    treatment_df <- cohort_df %>% filter(treatment == 1)
    control_df <- cohort_df %>% filter(treatment == 0)
    treatment_mean <- mean(treatment_df$followup_years)
    control_mean <- mean(control_df$followup_years)
    treatment_sd <- sd(treatment_df$followup_years)
    control_sd <- sd(control_df$followup_years)
  }
  return(list(
    treatment_mean = treatment_mean,
    treatment_sd = treatment_sd,
    control_mean = control_mean,
    control_sd = control_sd
  ))
}

# ============================================================================
# CRUDE HR CALCULATOR
# ============================================================================
calculate_crude_hr <- function(cohort_df, comparison_name, outcome_label) {
  tryCatch({
    crude_cox <- coxph(Surv(event_time, event) ~ treatment,
                       data = cohort_df)
    hr_crude <- exp(coef(crude_cox)["treatment"])
    ci_crude <- exp(confint(crude_cox)["treatment", ])
    p_crude <- summary(crude_cox)$coefficients["treatment", "Pr(>|z|)"]
    followup_stats <- calculate_followup_stats(cohort_df, weight_var = NULL)
    result <- data.frame(
      comparison = comparison_name,
      outcome = outcome_label,
      analysis_type = "Crude_Unadjusted",
      n_total = nrow(cohort_df),
      n_events = sum(cohort_df$event),
      n_treatment = sum(cohort_df$treatment == 1),
      n_control = sum(cohort_df$treatment == 0),
      events_treatment = sum(cohort_df$treatment == 1 & cohort_df$event == 1),
      events_control = sum(cohort_df$treatment == 0 & cohort_df$event == 1),
      followup_years_treatment_mean = followup_stats$treatment_mean,
      followup_years_treatment_sd = followup_stats$treatment_sd,
      followup_years_control_mean = followup_stats$control_mean,
      followup_years_control_sd = followup_stats$control_sd,
      hr = hr_crude,
      hr_lower = ci_crude[1],
      hr_upper = ci_crude[2],
      p_value = p_crude,
      stringsAsFactors = FALSE
    )
    log_output("INFO", "Crude HR calculated",
               sprintf("%s - HR=%.3f (%.3f-%.3f), p=%.4f",
                       outcome_label, hr_crude, ci_crude[1], ci_crude[2], p_crude),
               "crude_hr")
    return(result)
  }, error = function(e) {
    log_output("ERROR", "Crude HR calculation failed", e$message, "crude_hr")
    return(NULL)
  })
}

# >>> CELL 17: Variable Labels & Helpers <<<
# ============================================================================
# VARIABLE LABELING FUNCTION
# ============================================================================
apply_variable_labels <- function(data) {
  data %>%
    mutate(
      sex_cat = factor(case_when(
        sex_cat == 0 ~ "Male",
        sex_cat == 1 ~ "Female",
        sex_cat == 999 ~ "Missing",
        TRUE ~ "Missing"
      ), levels = c("Male", "Female", "Missing")),
      raceethnicity_cat = factor(case_when(
        raceethnicity_cat == 0 ~ "Non-Hispanic White",
        raceethnicity_cat == 1 ~ "Non-Hispanic Black",
        raceethnicity_cat == 2 ~ "Hispanic",
        raceethnicity_cat == 3 ~ "Other",
        raceethnicity_cat == 999 ~ "Missing",
        TRUE ~ "Missing"
      ), levels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other", "Missing")),
      income = factor(case_when(
        income == 0 ~ "Low (<50k)",
        income == 1 ~ "Middle (50k-100k)",
        income == 2 ~ "High (>100k)",
        income == 999 ~ "Missing",
        TRUE ~ "Missing"
      ), levels = c("Low (<50k)", "Middle (50k-100k)", "High (>100k)", "Missing")),
      education = factor(case_when(
        education == 0 ~ "High School or Less",
        education == 1 ~ "Some College",
        education == 2 ~ "Advanced Degree",
        education == 999 ~ "Missing",
        TRUE ~ "Missing"
      ), levels = c("High School or Less", "Some College", "Advanced Degree", "Missing")),
      smoking = factor(case_when(
        smoking == 0 ~ "Never",
        smoking == 1 ~ "Former",
        smoking == 2 ~ "Current",
        smoking == 999 ~ "Missing",
        TRUE ~ "Missing"
      ), levels = c("Never", "Former", "Current", "Missing")),
      insurance_category = factor(case_when(
        insurance_category == 0 ~ "None",
        insurance_category == 1 ~ "Public",
        insurance_category == 2 ~ "Private",
        insurance_category == 999 ~ "Missing",
        TRUE ~ "Missing"
      ), levels = c("None", "Public", "Private", "Missing")),
      alcohol_category = factor(case_when(
        alcohol_category == 0 ~ "Low Risk",
        alcohol_category == 1 ~ "Increased Risk",
        alcohol_category == 2 ~ "High Risk",
        alcohol_category == 3 ~ "Dependent",
        alcohol_category == 999 ~ "Missing",
        TRUE ~ "Missing"
      ), levels = c("Low Risk", "Increased Risk", "High Risk", "Dependent", "Missing")),
across(c("SEMAGLUTIDE", "OTHER_GLPA", "Biguanide", "TZD", "Insulin", "SGLT2i", "DPP4i", "SU",
         "Anticoagulant", "Antiplatelet", "Statin", "Ezetimibe", "RAAS", "Diuretic",
         "MRA", "BB", "CCB", "OtherHTN", "myocardial_infarction",
         "congestive_heart_failure", "peripheral_vascular_disease",
         "cerebrovascular_disease", "dementia", "chronic_pulmonary_disease",
         "rheumatic_disease", "peptic_ulcer_disease", "hiv_infection",
         "hypoglycemia", "hyperglycemic_emergency",
                "diabetes_with_ophthalmic_complications", "diabetes_with_neurological_complications", 
          "t1d", "t2d", "otherdm"),
       ~ factor(case_when(
         . == 0 ~ "No",
         . == 1 ~ "Yes",
         TRUE ~ "No"  # Treat NA or unexpected values as "No" instead of "Missing"
       ), levels = c("No", "Yes")))  # Only include "No" and "Yes" levels
    )
}

# ---- Ingredient name sets (same as earlier) ----
INGR_SEMA     <- c("semaglutide")
INGR_OTHERGLP <- c("exenatide","liraglutide","albiglutide","dulaglutide","lixisenatide")
INGR_SGLT2    <- c("canagliflozin","empagliflozin","dapagliflozin","ertugliflozin")
INGR_OTHERGLD <- c("glimepiride","glipizide","glyburide",
                   "alogliptin","linagliptin","sitagliptin","saxagliptin",
                   "pioglitazone","rosiglitazone")

# ---- Map comparison name -> exposure/comparator ingredient sets ----
get_exp_ing_by_comp <- function(comp_name) {
  switch(comp_name,
    "SEMAGLUTIDE vs OTHER_GLPA" = INGR_SEMA,
    "SEMAGLUTIDE vs SGLT2"      = INGR_SEMA,
    "SEMAGLUTIDE vs OtherGLD"   = INGR_SEMA,
    "OTHER_GLPA vs SGLT2"       = INGR_OTHERGLP,
    "OTHER_GLPA vs OtherGLD"    = INGR_OTHERGLP,
    "SGLT2 vs OtherGLD"         = INGR_SGLT2,
    NULL
  )
}
get_comp_ing_by_comp <- function(comp_name) {
  switch(comp_name,
    "SEMAGLUTIDE vs OTHER_GLPA" = INGR_OTHERGLP,
    "SEMAGLUTIDE vs SGLT2"      = INGR_SGLT2,
    "SEMAGLUTIDE vs OtherGLD"   = INGR_OTHERGLD,
    "OTHER_GLPA vs SGLT2"       = INGR_SGLT2,
    "OTHER_GLPA vs OtherGLD"    = INGR_OTHERGLD,
    "SGLT2 vs OtherGLD"         = INGR_OTHERGLD,
    NULL
  )
}

# Helper: number of non-missing unique levels
n_levels <- function(x) length(unique(stats::na.omit(x)))

# ---- Core routine used inside the IPTW script (run once per comparison/analysis) ----
tag_and_make_tables_for_both_arms <- function(cohort_df, comp_name, analysis_prefix) {
  message("Tagging index ingredients & building tables for BOTH arms (IPTW cohort prep): ", comp_name)

  # Normalize types
  cohort_df <- cohort_df %>%
    dplyr::mutate(
      person_id  = as.character(person_id),
      index_date = as.Date(index_date),
      treatment  = as.integer(as.character(treatment))
    )

  # Ingredient name sets (assumes helpers exist in your env)
  exp_names  <- get_exp_ing_by_comp(comp_name)
  comp_names <- get_comp_ing_by_comp(comp_name)

  # Tag treatment arm
  if (!is.null(exp_names)) {
    cohort_df <- get_index_ingredient_for_arm(
      cohort_df, exp_names, arm_value = 1, out_col = "index_ingredient_name_treat"
    )
  } else {
    cohort_df$index_ingredient_name_treat <- NA_character_
  }

  # Tag comparator arm
  if (!is.null(comp_names)) {
    cohort_df <- get_index_ingredient_for_arm(
      cohort_df, comp_names, arm_value = 0, out_col = "index_ingredient_name_comp"
    )
  } else {
    cohort_df$index_ingredient_name_comp <- NA_character_
  }

  # Apply labels
  cohort_df <- apply_variable_labels(cohort_df)

  # ---------- TREATMENT ARM ----------
  treat_df <- cohort_df %>% dplyr::filter(treatment == 1L)

  # Backstop for single-ingredient classes
  if (!is.null(exp_names) && length(exp_names) == 1 &&
      (all(is.na(treat_df$index_ingredient_name_treat)) ||
       !"index_ingredient_name_treat" %in% names(treat_df))) {
    treat_df$index_ingredient_name_treat <- tolower(exp_names[[1]])
  }

  # Counts by ingredient (treatment)
  n_treat <- treat_df %>%
    dplyr::mutate(index_ingredient_name_treat = tolower(index_ingredient_name_treat)) %>%
    dplyr::count(index_ingredient_name_treat, name = "N") %>%
    dplyr::arrange(dplyr::desc(N))
  # All of Us small-cell suppression: mask ingredient counts < 20 per Data Dissemination Policy
  n_treat$N[n_treat$N < 20] <- NA_integer_
  n_treat$index_ingredient_name_treat[is.na(n_treat$N)] <- "<suppressed>"

  write.csv(
    n_treat,
    file = sprintf("N_by_ingredient_TREAT_%s_%s_%s.csv",
                   gsub("[^A-Za-z0-9]+","_", comp_name),
                   analysis_prefix,
                   format(Sys.time(), "%Y%m%d_%H%M%S")),
    row.names = FALSE
  )

  # Per-ingredient TableOne for TREATMENT (kept)
  if (n_levels(treat_df$index_ingredient_name_treat) >= 2) {
    tbl_treat_ing <- CreateTableOne(
      vars       = my_vars,
      strata     = "index_ingredient_name_treat",
      data       = treat_df,
      factorVars = categorical_vars
    )
    df_treat_ing <- as.data.frame(print(tbl_treat_ing, showAllLevels = TRUE, printToggle = FALSE, noSpaces = TRUE))
    write.csv(
      df_treat_ing,
      file = sprintf("demographics_by_ingredient_TREAT_%s_%s_%s.csv",
                     gsub("[^A-Za-z0-9]+","_", comp_name),
                     analysis_prefix,
                     format(Sys.time(), "%Y%m%d_%H%M%S")),
      row.names = TRUE
    )
  } else {
    message("Note: <2 ingredient levels in TREATMENT arm; skipping per-ingredient TableOne for ", comp_name, ".")
  }

  # Class-combined (treatment) — kept
  tbl_treat_cls <- CreateTableOne(
    vars       = my_vars,
    data       = treat_df,
    factorVars = categorical_vars
  )
  df_treat_cls <- as.data.frame(print(tbl_treat_cls, showAllLevels = TRUE, printToggle = FALSE, noSpaces = TRUE))
  write.csv(
    df_treat_cls,
    file = sprintf("demographics_by_class_TREAT_%s_%s_%s.csv",
                   gsub("[^A-Za-z0-9]+","_", comp_name),
                   analysis_prefix,
                   format(Sys.time(), "%Y%m%d_%H%M%S")),
    row.names = TRUE
  )

  # ---------- COMPARATOR ARM ----------
  comp_df <- cohort_df %>% dplyr::filter(treatment == 0L)

  if (!is.null(comp_names) && length(comp_names) == 1 &&
      (all(is.na(comp_df$index_ingredient_name_comp)) ||
       !"index_ingredient_name_comp" %in% names(comp_df))) {
    comp_df$index_ingredient_name_comp <- tolower(comp_names[[1]])
  }

  # Keep ONLY the comparator **counts** file; remove comparator demographics outputs
  n_comp <- comp_df %>%
    dplyr::mutate(index_ingredient_name_comp = tolower(index_ingredient_name_comp)) %>%
    dplyr::count(index_ingredient_name_comp, name = "N") %>%
    dplyr::arrange(dplyr::desc(N))
  # All of Us small-cell suppression: mask ingredient counts < 20 per Data Dissemination Policy
  n_comp$N[n_comp$N < 20] <- NA_integer_
  n_comp$index_ingredient_name_comp[is.na(n_comp$N)] <- "<suppressed>"

  write.csv(
    n_comp,
    file = sprintf("N_by_ingredient_COMP_%s_%s_%s.csv",
                   gsub("[^A-Za-z0-9]+","_", comp_name),
                   analysis_prefix,
                   format(Sys.time(), "%Y%m%d_%H%M%S")),
    row.names = FALSE
  )

  # NOTE:
  # Removed:
  # - Per-ingredient comparator TableOne & write.csv (demographics_by_ingredient_COMP_*.csv)
  # - Class-combined comparator TableOne & write.csv (demographics_by_class_COMP_*.csv)

  # Return tagged cohort for downstream use
  cohort_df
}

# ============================================================================
# MAIN IPTW ANALYSIS SCRIPT (patched)
# - Fix 1: Tagging/tables run on post-exclusion set (align counts with analysis N)
# - Fix 2: Baseline table reattaches descriptives + normalizes binaries to No/Yes
# ============================================================================

# Load required packages
library(tidyverse)
library(survival)
# If not already loaded elsewhere:
# library(tableone)

# -------------------- Helpers (add once, near the top) -----------------------
# Detect 0/1 or logical or yes/no-like character/factor
is_binary_like <- function(x) {
  if (is.logical(x)) return(TRUE)
  if (is.numeric(x)) {
    ux <- unique(na.omit(x))
    return(length(ux) <= 2 && all(ux %in% c(0, 1)))
  }
  if (is.factor(x) || is.character(x)) {
    xc <- tolower(trimws(as.character(x)))
    ok <- xc %in% c("0","1","no","yes","n","y","false","true","t","f")
    return(all(ok | is.na(xc)))
  }
  FALSE
}

# Normalize a binary-like column to factor("No","Yes")
to_yesno <- function(x) {
  if (is.logical(x)) {
    return(factor(ifelse(x, "Yes", "No"), levels = c("No","Yes")))
  }
  if (is.numeric(x)) {
    return(factor(ifelse(x == 1, "Yes", "No"), levels = c("No","Yes")))
  }
  xc <- tolower(trimws(as.character(x)))
  yes_vals <- c("1","yes","y","true","t")
  no_vals  <- c("0","no","n","false","f")
  out <- ifelse(xc %in% yes_vals, "Yes",
                ifelse(xc %in% no_vals, "No", NA))
  factor(out, levels = c("No","Yes"))
}
# ---------------------------------------------------------------------------

# >>> CELL 18: Main IPTW Analysis Loop <<<
# Load configuration
config <- get_standard_config()

# Define the different analyses to run
analysis_types <- list(
  list(name = "all_ages",  late_onset_flag = FALSE, prefix = "epilepsy_seizure_all_ages_ipwt")
)

# Create comparison_results from the 'cohorts' list
comparison_results <- lapply(cohorts, function(df) {
  list(cohort_df = df)
})
names(comparison_results) <- c(
  "SEMAGLUTIDE vs OTHER_GLPA",
  "SEMAGLUTIDE vs SGLT2",
  "SEMAGLUTIDE vs OtherGLD",
  "OTHER_GLPA vs SGLT2",
  "OTHER_GLPA vs OtherGLD",
  "SGLT2 vs OtherGLD"
)

# Outer loop for each analysis type
for (analysis in analysis_types) {
  log_output("PROGRESS", sprintf("STARTING ANALYSIS TYPE: %s", analysis$name), "", "main_script")
  ipwt_results <- list()

  # Process each comparison
  for (comp_name in names(comparison_results)) {
    log_output("PROGRESS", sprintf("Processing comparison: %s", comp_name), "", "main_analysis")
    cohort_df <- comparison_results[[comp_name]]$cohort_df

    # Apply new summary variables and labels (on the raw comparison cohort)
    cohort_df <- cohort_df %>%
      mutate(
        renal_disease_severity = factor(case_when(
          renal_severe == 1 ~ "Severe",
          renal_mild_or_moderate == 1 ~ "Mild-Moderate",
          TRUE ~ "None"
        ), levels = c("None", "Mild-Moderate", "Severe")),
        liver_disease_severity = factor(case_when(
          moderate_or_severe_liver_disease == 1 ~ "Mod-Severe",
          mild_liver_disease == 1 ~ "Mild",
          TRUE ~ "None"
        ), levels = c("None", "Mild", "Mod-Severe")),
        malignancy_status = factor(case_when(
          metastatic_solid_tumor == 1 ~ "Metastatic",
          any_malignancy == 1 ~ "Non-Metastatic",
          TRUE ~ "None"
        ), levels = c("None", "Non-Metastatic", "Metastatic"))#,
     #   diabetic_complication_count = diabetes_with_renal_complications +
      #    diabetes_with_ophthalmic_complications +
       #   diabetes_with_neurological_complications,
        #diabetic_complication_count_cat = factor(
         # case_when(
          #  diabetic_complication_count >= 2 ~ "2+",
           # diabetic_complication_count == 1 ~ "1",
            #TRUE ~ "0"
         # ), levels = c("0", "1", "2+")
       # )
      ) %>%
      apply_variable_labels()

    # Ensure categorical variables are factors (from config)
    for (var in config$categorical_vars) {
      if (var %in% names(cohort_df)) {
        cohort_df[[var]] <- as.factor(cohort_df[[var]])
      }
    }

    # ------------------------ IMPORTANT CHANGE ------------------------------
    # We NO LONGER run tag_and_make_tables_for_both_arms() here on the raw cohort.
    # We will run it later on the POST-EXCLUSION cohort inside the outcome loop.
    # ------------------------------------------------------------------------

    # Exclude vars for PS model by comparison
    exclude_vars <- switch(
      comp_name,
      "SEMAGLUTIDE vs OTHER_GLPA" = c("SEMAGLUTIDE", "OTHER_GLPA"),
      "SEMAGLUTIDE vs SGLT2"      = c("SEMAGLUTIDE", "SGLT2i"),
      "SEMAGLUTIDE vs OtherGLD"   = c("SEMAGLUTIDE", "TZD", "SU", "DPP4i"),
      "OTHER_GLPA vs SGLT2"       = c("OTHER_GLPA", "SGLT2i"),
      "OTHER_GLPA vs OtherGLD"    = c("OTHER_GLPA", "TZD", "SU", "DPP4i"),
      "SGLT2 vs OtherGLD"         = c("SGLT2i", "TZD", "SU", "DPP4i"),
      NULL
    )

    # Do tagging/tables once per comparison per analysis after exclusions
    did_tags_for_comp <- FALSE

    # Loop through each outcome
    for (outcome in config$outcomes) {
      outcome_var   <- outcome$var
      outcome_label <- outcome$label

      log_output("INFO",
                 sprintf("Processing outcome: %s for %s (%s)", outcome_label, comp_name, analysis$name),
                 "", "main_analysis")

      survival_df <- followup_and_event(
        df            = cohort_df,
        outcome_var   = outcome_var,
        data_cut_date = config$data_cut_date,
        late_onset    = analysis$late_onset_flag,
        verbose       = TRUE,
        comp_name     = comp_name
      )
      # Exclude severe renal disease and/or metastatic malignancy from analysis cohort (configurable)
      pre_n <- nrow(survival_df)
      cond_renal <- !is.na(survival_df$renal_disease_severity) & survival_df$renal_disease_severity == "Severe"
      cond_malig <- !is.na(survival_df$malignancy_status) & survival_df$malignancy_status == "Metastatic"
      # Build exclusion mask based on config toggles
      exclude_mask <- rep(FALSE, nrow(survival_df))
      if (isTRUE(config$exclude_renal_severe)) exclude_mask <- exclude_mask | cond_renal
      if (isTRUE(config$exclude_malignancy_metastatic)) exclude_mask <- exclude_mask | cond_malig
      n_excl_renal <- sum(cond_renal)
      n_excl_malig <- sum(cond_malig)
      n_excl_enabled <- sum(exclude_mask)
      if (isTRUE(config$log_exclusions)) {
        log_output("PROGRESS", "Exclusion counts (renal/metastatic)",
                   sprintf("renal Severe=%d, malignancy Metastatic=%d, enabled_exclusions=%d",
                           n_excl_renal, n_excl_malig, n_excl_enabled), "main_analysis")
      }
      survival_df <- survival_df %>%
        dplyr::filter(!exclude_mask) %>%
        droplevels()
      post_n <- nrow(survival_df)
      if (isTRUE(config$log_exclusions)) {
        log_output("PROGRESS", "Cohort size after exclusions",
                   sprintf("%d -> %d (removed %d)", pre_n, post_n, pre_n - post_n),
                   "main_analysis")
      }

      if (nrow(survival_df) < 10) {
        log_output("WARNING", "Skipping analysis due to insufficient data", nrow(survival_df)); next
      }
      n_t1 <- sum(survival_df$treatment == 1)
      n_t0 <- sum(survival_df$treatment == 0)
      if (n_t1 < 2 || n_t0 < 2) {
        log_output("ERROR", "Insufficient observations in one treatment group",
                   sprintf("Treatment=1: %d, Treatment=0: %d", n_t1, n_t0), "main_analysis")
        next
      }

    # MICE imputation for baseline_bmi and baseline_hba1c
      survival_df <- perform_mice_imputation(survival_df, m = 20, seed = 123)
    # Recompute BMI categories after imputation so PS complete-case does not drop imputed rows
      if ("baseline_bmi" %in% names(survival_df)) {
        survival_df <- survival_df %>% mutate(
          baseline_bmi_category = factor(
            case_when(
              is.na(baseline_bmi)                      ~ "Missing",
              baseline_bmi < 18.5                      ~ "Underweight",
              baseline_bmi >= 18.5 & baseline_bmi < 25 ~ "Normal",
              baseline_bmi >= 25   & baseline_bmi < 30 ~ "Overweight",
              baseline_bmi >= 30   & baseline_bmi < 35 ~ "Obesity Class I",
              baseline_bmi >= 35   & baseline_bmi < 40 ~ "Obesity Class II",
              baseline_bmi >= 40                       ~ "Obesity Class III"
            ),
            levels = c(
              "Underweight",
              "Normal",
              "Overweight",
              "Obesity Class I",
              "Obesity Class II",
              "Obesity Class III",
              "Missing"
            )
          )
        )
      }  
        
      # ------------------------ IMPORTANT CHANGE ------------------------------
      # Align everything to the ACTUAL analysis set (post-exclusions)
      keep_ids <- unique(survival_df$person_id)
      cohort_df_post <- cohort_df %>%
        dplyr::semi_join(tibble(person_id = keep_ids), by = "person_id")

      # Run tagging/tables exactly once per comparison per analysis (post-exclusion)
      if (!did_tags_for_comp) {
        invisible(
          tag_and_make_tables_for_both_arms(
            cohort_df       = cohort_df_post,
            comp_name       = comp_name,
            analysis_prefix = analysis$prefix
          )
        )
        did_tags_for_comp <- TRUE
      }
      # ------------------------------------------------------------------------

      ipwt_result <- tryCatch({
        run_ipwt_and_cox(
          df             = survival_df,
          exclude_vars   = exclude_vars,
          all_ps_vars    = config$all_ps_vars,
          trim_threshold = 0.01,
          verbose        = TRUE,
          stabilize      = FALSE
        )
      }, error = function(e) {
        log_output("ERROR", sprintf("run_ipwt_and_cox failed for %s", comp_name), e$message, "main_analysis")
        return(NULL)
      })
      if (is.null(ipwt_result)) next

      crude_hr <- tryCatch({
        calculate_crude_hr(
          cohort_df       = ipwt_result$cohort,
          comparison_name = comp_name,
          outcome_label   = outcome_label
        )
      }, error = function(e) {
        log_output("ERROR", sprintf("calculate_crude_hr failed for %s", comp_name), e$message, "main_analysis")
        return(NULL)
      })

      incident_rates <- tryCatch({
        calculate_incident_rates(
          cohort_df = ipwt_result$cohort,
          weight_var = "ipw_std"
        )
      }, error = function(e) {
        log_output("ERROR", sprintf("calculate_incident_rates failed for %s", comp_name), e$message, "main_analysis")
        return(NULL)
      })

      ipwt_result$crude_hr_results      <- crude_hr
      ipwt_result$incident_rate_results <- incident_rates
      analysis_key <- paste(comp_name, "-", outcome_label)
      ipwt_results[[analysis_key]] <- ipwt_result

      # --- Create baseline (demographics) table for the exact IPTW analysis cohort ---
      # Use the same variables and formatting as the pre-IPTW baseline table
      try({
        # Start from the exact modeling cohort
        iptw_cohort <- ipwt_result$cohort

        # Re-attach descriptives from the post-exclusion cohort set
        # Choose the variables you actually want to describe
        desc_vars <- unique(c(my_vars, categorical_vars))
        desc_vars <- intersect(desc_vars, names(cohort_df_post))
        # Don't overwrite existing columns in iptw_cohort
        keep_desc <- setdiff(desc_vars, names(iptw_cohort))

        if (length(keep_desc) > 0) {
          iptw_cohort <- iptw_cohort %>%
            dplyr::left_join(
              cohort_df_post %>% dplyr::select(person_id, dplyr::all_of(keep_desc)),
              by = "person_id"
            )
        }

        # Normalize binary-like columns to factor("No","Yes") once, right before the table
        bin_candidates <- names(iptw_cohort)[vapply(iptw_cohort, is_binary_like, logical(1))]
        if (length(bin_candidates) > 0) {
          iptw_cohort[bin_candidates] <- lapply(iptw_cohort[bin_candidates], to_yesno)
        }

        # Build the variable lists present in this cohort
        iptw_vars     <- intersect(my_vars, names(iptw_cohort))
        iptw_cat_vars <- intersect(categorical_vars, iptw_vars)

        # Ensure non-binary categoricals are factors
        for (v in iptw_cat_vars) {
          if (!is.factor(iptw_cohort[[v]])) {
            iptw_cohort[[v]] <- factor(iptw_cohort[[v]])
          }
        }

        # Mark normalized binaries as factor vars too (those now exactly c("No","Yes"))
        bin_factor_vars <- intersect(bin_candidates, iptw_vars)

        if (length(iptw_vars) > 0) {
          iptw_table_obj <- CreateTableOne(
            vars       = iptw_vars,
            strata     = "treatment",
            data       = iptw_cohort,
            factorVars = unique(c(iptw_cat_vars, bin_factor_vars))
          )
          iptw_table_df <- as.data.frame(print(
            iptw_table_obj,
            showAllLevels = TRUE,
            printToggle   = FALSE,
            noSpaces      = TRUE
          ))
          ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
          safe_comp    <- gsub("[^A-Za-z0-9]+", "_", comp_name)
          safe_outcome <- gsub("[^A-Za-z0-9]+", "_", outcome_label)
          fname <- sprintf(
            "epilepsy_seizure_baseline_characteristics_table_IPTW_%s_%s_%s_%s.csv",
            analysis$prefix, safe_comp, safe_outcome, ts
          )
          write.csv(iptw_table_df, file = fname, row.names = TRUE)
          log_output("PROGRESS", "Saved IPTW analysis cohort baseline table", fname, "baseline_table_IPTW")
        } else {
          log_output("WARNING", "No overlapping variables for IPTW baseline table; skipping",
                     paste(comp_name, outcome_label), "baseline_table_IPTW")
        }
      }, silent = TRUE)

      balance_plot <- tryCatch({
        create_balance_plot(
          balance_before = ipwt_result$balance_before,
          balance_after  = ipwt_result$balance_after,
          title          = sprintf("Balance for %s - %s (%s)", comp_name, outcome_label, analysis$name)
        )
      }, error = function(e) {
        log_output("ERROR", sprintf("create_balance_plot failed for %s", comp_name), e$message, "main_analysis")
        return(NULL)
      })
      if (!is.null(balance_plot)) {
        timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
        balance_plot_filename <- sprintf("%s_balance_%s_%s_%s.png",
                                         analysis$prefix,
                                         gsub(" ", "_", comp_name),
                                         gsub("[/ ]", "_", outcome_label),
                                         timestamp)
        ggsave(balance_plot_filename, balance_plot)
      }

      surv_plot <- tryCatch({
        create_survival_plot(
          ipwt_result    = ipwt_result,
          comparison_name= comp_name,
          outcome_label  = outcome_label,
          analysis_name  = analysis$name,
          save_plot      = TRUE
        )
      }, error = function(e) {
        log_output("ERROR", sprintf("create_survival_plot failed for %s", comp_name), e$message, "main_analysis")
        return(NULL)
      })
    } # end outcomes loop
  }   # end comparisons loop

  saved_files <- tryCatch({
    save_analysis_results(ipwt_results, prefix = analysis$prefix)
  }, error = function(e) {
    log_output("ERROR", "save_analysis_results failed", e$message, "main_analysis")
    NULL
  })

  rds_filename <- sprintf("%s_full_results.rds", analysis$prefix)
  saveRDS(ipwt_results, file = rds_filename)
  log_output("PROGRESS", "Full results object saved to", rds_filename, "main_script")

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  save_output_log(sprintf("%s_log_%s.csv", analysis$prefix, timestamp))

  log_output("PROGRESS", sprintf("ANALYSIS TYPE '%s' COMPLETE", analysis$name), "", "main_script")
}

message("All analyses for Epilepsy/Seizure completed successfully.")

# >>> CELL 19: Combined Results Table <<<
# =============================================================================
# COMBINED IPTW RESULTS TABLE ACROSS ALL OUTCOMES & ANALYSIS TYPES
# =============================================================================
combined_rows <- list()
for (analysis in analysis_types) {
  rds_file <- sprintf("%s_full_results.rds", analysis$prefix)
  if (!file.exists(rds_file)) {
    message("RDS not found: ", rds_file)
    next
  }
  ipwt_obj <- readRDS(rds_file)
  for (key in names(ipwt_obj)) {
    res <- ipwt_obj[[key]]
    if (is.null(res) || is.null(res$cox)) next
    split_name <- strsplit(key, " - ")[[1]]
    comp  <- trimws(split_name[1])
    outc  <- trimws(split_name[2])
    cox_s <- summary(res$cox)
    hr    <- cox_s$conf.int[1, "exp(coef)"]
    lo    <- cox_s$conf.int[1, "lower .95"]
    hi    <- cox_s$conf.int[1, "upper .95"]
    pval  <- cox_s$coefficients[1, "Pr(>|z|)"]
    n_total  <- nrow(res$cohort)
    n_events <- sum(res$cohort$event, na.rm = TRUE)
    n_treat  <- sum(res$cohort$treatment == 1)
    n_ctrl   <- sum(res$cohort$treatment == 0)
    ev_treat <- sum(res$cohort$event[res$cohort$treatment == 1], na.rm = TRUE)
    ev_ctrl  <- sum(res$cohort$event[res$cohort$treatment == 0], na.rm = TRUE)
    combined_rows[[length(combined_rows) + 1]] <- data.frame(
      analysis_type = analysis$name,
      comparison    = comp,
      outcome       = outc,
      n_total       = n_total,
      n_treatment   = n_treat,
      n_control     = n_ctrl,
      events_total  = n_events,
      events_treat  = ev_treat,
      events_ctrl   = ev_ctrl,
      hr_weighted   = round(hr, 3),
      hr_lower_95   = round(lo, 3),
      hr_upper_95   = round(hi, 3),
      hr_ci         = sprintf("%.2f (%.2f-%.2f)", hr, lo, hi),
      p_value       = signif(pval, 4),
      stringsAsFactors = FALSE
    )
  }
}
if (length(combined_rows) > 0) {
  combined_df <- do.call(rbind, combined_rows)
  combined_df <- combined_df[order(combined_df$analysis_type, combined_df$comparison, combined_df$outcome), ]
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  combined_file <- sprintf("IPTW_all_outcomes_combined_results_%s.csv", ts)
  write.csv(combined_df, combined_file, row.names = FALSE)
  message("Combined IPTW results saved to: ", combined_file)
  print(combined_df[, c("analysis_type", "comparison", "outcome", "n_total",
                         "n_treatment", "n_control",
                         "events_total", "events_treat", "events_ctrl",
                         "hr_ci", "p_value")])
}

res_all  <- readRDS("epilepsy_seizure_all_ages_ipwt_full_results.rds")
res_late <- readRDS("epilepsy_seizure_late_onset_ipwt_full_results.rds")

ids_all  <- unique(unlist(lapply(res_all, function(x) x$cohort$person_id)))
ids_late <- unique(unlist(lapply(res_late, function(x) x$cohort$person_id)))

cat("All-ages N:", length(ids_all), "\n")
cat("Late-onset N:", length(ids_late), "\n")
cat("Difference:", length(setdiff(ids_late, ids_all)), "added;",
    length(setdiff(ids_all, ids_late)), "lost\n")

# >>> CELL 19: Combined Results Table <<<
# =============================================================================
# COMBINED IPTW RESULTS TABLE ACROSS ALL OUTCOMES & ANALYSIS TYPES
# =============================================================================
combined_rows <- list()
for (analysis in analysis_types) {
  rds_file <- sprintf("%s_full_results.rds", analysis$prefix)
  if (!file.exists(rds_file)) {
    message("RDS not found: ", rds_file)
    next
  }
  ipwt_obj <- readRDS(rds_file)
  for (key in names(ipwt_obj)) {
    res <- ipwt_obj[[key]]
    if (is.null(res) || is.null(res$cox)) next
    split_name <- strsplit(key, " - ")[[1]]
    comp  <- trimws(split_name[1])
    outc  <- trimws(split_name[2])
    cox_s <- res$cox  # already a summary.coxph object
    hr    <- cox_s$conf.int[1, "exp(coef)"]
    lo    <- cox_s$conf.int[1, "lower .95"]
    hi    <- cox_s$conf.int[1, "upper .95"]
    pval  <- cox_s$coefficients[1, "Pr(>|z|)"]
    n_total  <- nrow(res$cohort)
    n_events <- sum(res$cohort$event, na.rm = TRUE)
    n_treat  <- sum(res$cohort$treatment == 1)
    n_ctrl   <- sum(res$cohort$treatment == 0)
    ev_treat <- sum(res$cohort$event[res$cohort$treatment == 1], na.rm = TRUE)
    ev_ctrl  <- sum(res$cohort$event[res$cohort$treatment == 0], na.rm = TRUE)
    combined_rows[[length(combined_rows) + 1]] <- data.frame(
      analysis_type = analysis$name,
      comparison    = comp,
      outcome       = outc,
      n_total       = n_total,
      n_treatment   = n_treat,
      n_control     = n_ctrl,
      events_total  = n_events,
      events_treat  = ev_treat,
      events_ctrl   = ev_ctrl,
      hr_weighted   = round(hr, 3),
      hr_lower_95   = round(lo, 3),
      hr_upper_95   = round(hi, 3),
      hr_ci         = sprintf("%.2f (%.2f-%.2f)", hr, lo, hi),
      p_value       = signif(pval, 4),
      stringsAsFactors = FALSE
    )
  }
}
if (length(combined_rows) > 0) {
  combined_df <- do.call(rbind, combined_rows)
  combined_df <- combined_df[order(combined_df$analysis_type, combined_df$comparison, combined_df$outcome), ]
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
  combined_file <- sprintf("IPTW_all_outcomes_combined_results_%s.csv", ts)
  write.csv(combined_df, combined_file, row.names = FALSE)
  message("Combined IPTW results saved to: ", combined_file)
  print(combined_df[, c("analysis_type", "comparison", "outcome", "n_total",
                         "n_treatment", "n_control",
                         "events_total", "events_treat", "events_ctrl",
                         "hr_ci", "p_value")])
}

res_all  <- readRDS("epilepsy_seizure_all_ages_ipwt_full_results.rds")
res_late <- readRDS("epilepsy_seizure_late_onset_ipwt_full_results.rds")

ids_all  <- unique(unlist(lapply(res_all, function(x) x$cohort$person_id)))
ids_late <- unique(unlist(lapply(res_late, function(x) x$cohort$person_id)))

cat("All-ages N:", length(ids_all), "\n")
cat("Late-onset N:", length(ids_late), "\n")
cat("Difference:", length(setdiff(ids_late, ids_all)), "added;",
    length(setdiff(ids_all, ids_late)), "lost\n")



