# >>> CELL 01: INSTALL PACKAGES
install.packages("MatchIt")
install.packages("survminer")
install.packages("mice")
install.packages("tableone")

# >>> CELL 02: LOAD LIBRARIES
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

# >>> CELL 03: CONFIG
# ---------- CONFIG ----------
# Outcome 4: Epilepsy/Seizure excluding Hypoglycemic Seizure
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
      list(var = "epilepsy_or_seizure_without_hypoglycemic_seizure_start_date",
           label = "Epilepsy/Seizure excl Hypoglycemic",
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

# >>> CELL 04: MICE IMPUTATION FUNCTION
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

# >>> CELL 05: LOGGING FUNCTIONS
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

# >>> CELL 06: DATA LOAD FUNCTIONS
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

# >>> CELL 07: EXPOSURE INDEX BUILDER
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

# >>> CELL 08: CTE HELPER
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

# >>> CELL 09: MEDICATION FLAGGING
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
    # nothing to query -- return empty flags frame (to be joined upstream)
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

# >>> CELL 10: COMORBIDITY FLAGGING
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
    # nothing to query -- return empty flags frame (to be joined upstream)
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

# >>> CELL 11: BASELINE MEASUREMENTS
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

# >>> CELL 12: PROCESS COMPARISON GROUP (MODIFIED FOR OUTCOME 4)
# ---------- COMPARISON GROUP PROCESSOR ----------
# Modified for Outcome 4: uses epilepsy_or_seizure_without_hypoglycemic_seizure_start_date
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
  # OUTCOME 4: Use epilepsy_or_seizure_without_hypoglycemic_seizure_start_date
  pre_mask <- !is.na(cohort_df$epilepsy_or_seizure_without_hypoglycemic_seizure_start_date) &
              (as.Date(cohort_df$epilepsy_or_seizure_without_hypoglycemic_seizure_start_date) < as.Date(cohort_df$index_date))
  n_pre_total <- sum(pre_mask, na.rm = TRUE)
  if (!is.na(n_pre_total) && n_pre_total > 0) {
    comp_label <- paste0(exposure_name, " vs ", comparator_name)
    by_trt <- tryCatch({
      tbl <- cohort_df[pre_mask, c("treatment")]
      if (nrow(tbl) == 0 || !"treatment" %in% names(tbl)) data.frame(treatment = integer(0), n = integer(0)) else {
        as.data.frame(stats::aggregate(rep(1, nrow(tbl)) ~ tbl$treatment, FUN = sum))
      }
    }, error = function(e) data.frame())
    log_output("PROGRESS", "Pre-index epilepsy/seizure (excl hypoglycemic) detected (will be excluded)",
               sprintf("%s: n_total=%d%s",
                       comp_label, n_pre_total,
                       if (nrow(by_trt) > 0) paste0("; by_treatment=", paste0(by_trt[[1]], ":", by_trt[[2]], collapse=",")) else ""),
               "process_comparison_group")
    # Save IDs for audit trail
    out_name <- sprintf("preindex_epilepsy_excl_hypo_exclusions_%s_%s.csv",
                        gsub("[^A-Za-z0-9]+", "_", comp_label),
                        format(Sys.time(), "%Y%m%d_%H%M%S"))
    safe_cols <- intersect(c("person_id","index_date","epilepsy_or_seizure_without_hypoglycemic_seizure_start_date","treatment"), names(cohort_df))
    if (length(safe_cols) > 0) {
      try(utils::write.csv(cohort_df[pre_mask, safe_cols, drop = FALSE], file = out_name, row.names = FALSE), silent = TRUE)
    }
  }

  cohort_df <- cohort_df %>%
    mutate(preindex_epilepsy_flag = ifelse(pre_mask, 1L, 0L)) %>%
    filter(preindex_epilepsy_flag == 0L)
  list(temporal_stats = list(final_cohort = nrow(cohort_df)), cohort_df = cohort_df)
}

# >>> CELL 13: COHORT SUMMARY AND DRUG CONCEPTS
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

# >>> CELL 14: FOLLOWUP AND EVENT FUNCTION
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

# >>> CELL 15: MAIN WORKFLOW - COHORT BUILDING WITH FLOW CHART COUNTS
# --- SECTION 4: MAIN WORKFLOW ---
# Outcome 4: Epilepsy/Seizure excluding Hypoglycemic Seizure - Flow Chart Version
perform_stepwise_cohort_analysis_optimized <- function() {
  log_output("PROGRESS", "Starting Epilepsy/Seizure excl Hypoglycemic cohort analysis (Flow Chart)", "", "")
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
    select(person_id, t2d_excluded_t1d_start_date, death_date, age, sex_cat, raceethnicity_cat, income, education, insurance_category, smoking, alcohol_category, EHRmaxDT_min, EHRmaxDT_max, epilepsy_or_seizure_without_hypoglycemic_seizure_start_date, adrd_start_date, stroke_start_date)

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
    clean_comparison_results[[comp_name]] <- list(cohort_df = final_df)
  }
  # Step 7: Final Summary and Output
  results_df <- create_cohort_summary_df(cohort_stats, clean_comparison_results)
  output_file <- sprintf("outcome4_flowchart_cohort_analysis_results_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S"))
  write.csv(results_df, output_file, row.names = FALSE)
  save_output_log("outcome4_flowchart_cohort_analysis_log.csv")
  log_output("PROGRESS", "Cohort analysis complete", output_file, "")

  list(stats = cohort_stats, comparisons = clean_comparison_results, output_file = output_file)
}

# >>> CELL 16: SAVE FULL COHORTS TO BUCKET
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
 }
 invisible(TRUE)
}

# >>> CELL 17: RUN COHORT BUILDING AND SAVE
# ---------- RUN and SAVE ----------
# Run cohort building and save full cohorts for external IPWT
analysis_results <- perform_stepwise_cohort_analysis_optimized()
comparison_results <- analysis_results$comparisons
tryCatch({
 save_full_cohorts_to_bucket(comparison_results)
}, error = function(e) {
 log_output("ERROR", "Failed to save full cohorts to bucket", e$message, "main_analysis")
})

# >>> CELL 18: LOAD COHORT FILES FROM DISK
# --- PART 1: DATA LOADING ---
# Function to copy file from Google Bucket and load into a dataframe
load_from_bucket <- function(file_name) {
  if (!file.exists(file_name)) {
    stop(paste("File not found:", file_name,
               "Please ensure it's in the working directory."))
  }
  read_csv(file_name, show_col_types = FALSE)
}

# Outcome 4 cohort files
# NOTE: Update the COHORT_TIMESTAMP below to match your output file timestamps
COHORT_TIMESTAMP <- "20260213_053619"
files <- sprintf(
  "%s_semaglutide_epilepsy_seizure_full_cohort_%s.csv",
  c("SEMAGLUTIDE_vs_OTHER_GLPA", "SEMAGLUTIDE_vs_SGLT2", "SEMAGLUTIDE_vs_OtherGLD",
    "OTHER_GLPA_vs_SGLT2", "OTHER_GLPA_vs_OtherGLD", "SGLT2_vs_OtherGLD"),
  COHORT_TIMESTAMP
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

# >>> CELL 19: NORMALIZE COHORT TYPES
cohorts <- lapply(cohorts, function(df) {
  df %>% dplyr::mutate(
    person_id = as.character(person_id),
    index_date = as.Date(index_date),
    # make sure treatment is 0/1 integer (not factor/char)
    treatment  = as.integer(as.character(treatment))
  )
})

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

# >>> CELL 20: FLOW CHART COUNTING LOGIC
# ============================================================================
# FLOW CHART: Count participants at each cohort-building step
# Outcome 4: Epilepsy/Seizure excl Hypoglycemic
# ============================================================================

cat("\n")
cat("================================================================\n")
cat("  FLOW CHART COUNTS: Epilepsy/Seizure excl Hypoglycemic\n")
cat("================================================================\n\n")

config <- get_standard_config()
outcome_var <- "epilepsy_or_seizure_without_hypoglycemic_seizure_start_date"
outcome_label <- "Epilepsy/Seizure excl Hypoglycemic"

# Initialize flow chart tracking data frame
flowchart_rows <- list()

for (comp_name in names(comparison_results)) {
  cat(sprintf("\n--- %s ---\n", comp_name))

  cohort_df <- comparison_results[[comp_name]]$cohort_df

  # Step 0: Total loaded from cohort CSV (already passed washout + T2D + pre-index + crossover)
  n_total_loaded <- nrow(cohort_df)
  n_treat_loaded <- sum(cohort_df$treatment == 1)
  n_comp_loaded  <- sum(cohort_df$treatment == 0)

  cat(sprintf("  [Step 0] Loaded from cohort CSV:                  N = %d (Treat = %d, Comp = %d)\n",
              n_total_loaded, n_treat_loaded, n_comp_loaded))

  # Step 1: After same-day initiator exclusion
  # Apply the same-day initiator filter
  if (comp_name == "SEMAGLUTIDE vs OTHER_GLPA") {
    cohort_clean <- cohort_df %>% filter(!((treatment == 1 & OTHER_GLPA == 1) | (treatment == 0 & SEMAGLUTIDE == 1)))
  } else if (comp_name == "SEMAGLUTIDE vs SGLT2") {
    cohort_clean <- cohort_df %>% filter(!((treatment == 1 & SGLT2i == 1) | (treatment == 0 & SEMAGLUTIDE == 1)))
  } else if (comp_name == "SEMAGLUTIDE vs OtherGLD") {
    cohort_clean <- cohort_df %>% filter(!((treatment == 1 & (TZD == 1 | DPP4i == 1 | SU == 1)) | (treatment == 0 & SEMAGLUTIDE == 1)))
  } else if (comp_name == "OTHER_GLPA vs SGLT2") {
    cohort_clean <- cohort_df %>% filter(!((treatment == 1 & SGLT2i == 1) | (treatment == 0 & OTHER_GLPA == 1)))
  } else if (comp_name == "OTHER_GLPA vs OtherGLD") {
    cohort_clean <- cohort_df %>% filter(!((treatment == 1 & (TZD == 1 | DPP4i == 1 | SU == 1)) | (treatment == 0 & OTHER_GLPA == 1)))
  } else if (comp_name == "SGLT2 vs OtherGLD") {
    cohort_clean <- cohort_df %>% filter(!((treatment == 1 & (TZD == 1 | DPP4i == 1 | SU == 1)) | (treatment == 0 & SGLT2i == 1)))
  }

  n_after_sameday <- nrow(cohort_clean)
  n_removed_sameday <- n_total_loaded - n_after_sameday
  n_treat_after_sameday <- sum(cohort_clean$treatment == 1)
  n_comp_after_sameday  <- sum(cohort_clean$treatment == 0)

  cat(sprintf("  [Step 1] After same-day initiator removal:        N = %d (Treat = %d, Comp = %d) [removed %d]\n",
              n_after_sameday, n_treat_after_sameday, n_comp_after_sameday, n_removed_sameday))

  # Step 2: Apply follow-up and event processing (removes non-positive follow-up)
  survival_df <- followup_and_event(
    df            = cohort_clean,
    outcome_var   = outcome_var,
    data_cut_date = config$data_cut_date,
    late_onset    = FALSE,
    verbose       = TRUE,
    comp_name     = comp_name
  )

  n_after_followup <- nrow(survival_df)
  n_removed_followup <- n_after_sameday - n_after_followup
  n_treat_after_followup <- sum(survival_df$treatment == 1)
  n_comp_after_followup  <- sum(survival_df$treatment == 0)

  cat(sprintf("  [Step 2] After follow-up filtering (time > 0):    N = %d (Treat = %d, Comp = %d) [removed %d]\n",
              n_after_followup, n_treat_after_followup, n_comp_after_followup, n_removed_followup))

  # Step 3: Count events
  n_events_total <- sum(survival_df$event)
  n_events_treat <- sum(survival_df$event[survival_df$treatment == 1])
  n_events_comp  <- sum(survival_df$event[survival_df$treatment == 0])

  cat(sprintf("  [Events] Total = %d (Treat = %d, Comp = %d)\n",
              n_events_total, n_events_treat, n_events_comp))

  # Final summary
  cat(sprintf("  [FINAL]  N = %d | Treatment = %d | Comparator = %d | Events = %d\n",
              n_after_followup, n_treat_after_followup, n_comp_after_followup, n_events_total))

  # Collect row for CSV output
  flowchart_rows[[comp_name]] <- data.frame(
    comparison = comp_name,
    outcome = outcome_label,
    n_loaded_total = n_total_loaded,
    n_loaded_treatment = n_treat_loaded,
    n_loaded_comparator = n_comp_loaded,
    n_removed_sameday_initiators = n_removed_sameday,
    n_after_sameday_total = n_after_sameday,
    n_after_sameday_treatment = n_treat_after_sameday,
    n_after_sameday_comparator = n_comp_after_sameday,
    n_removed_nonpos_followup = n_removed_followup,
    n_final_total = n_after_followup,
    n_final_treatment = n_treat_after_followup,
    n_final_comparator = n_comp_after_followup,
    n_events_total = n_events_total,
    n_events_treatment = n_events_treat,
    n_events_comparator = n_events_comp,
    stringsAsFactors = FALSE
  )
}

# >>> CELL 21: SAVE FLOW CHART COUNTS TO CSV
# Combine and save flow chart counts
flowchart_df <- do.call(rbind, flowchart_rows)
rownames(flowchart_df) <- NULL

output_filename <- "outcome4_flowchart_counts.csv"
write.csv(flowchart_df, file = output_filename, row.names = FALSE)
cat(sprintf("\nFlow chart counts saved to: %s\n", output_filename))

# Print final summary table
cat("\n================================================================\n")
cat("  SUMMARY TABLE\n")
cat("================================================================\n")
print(flowchart_df %>% select(comparison, n_loaded_total, n_after_sameday_total, n_final_total,
                               n_final_treatment, n_final_comparator, n_events_total))

cat("\n")
cat("Flow chart analysis complete for Outcome 4: Epilepsy/Seizure excl Hypoglycemic.\n")

# Save output log
save_output_log(sprintf("outcome4_flowchart_log_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S")))
