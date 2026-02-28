# >>> CELL 01: Install & Load Libraries <<<
install.packages("MatchIt")
# install.packages("survminer")  # Optional: only needed for survival plots
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
  # library(survminer)  # Optional: only needed for survival plots
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

# >>> CELL 05: Data Load Functions <<<
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

# >>> CELL 06: Exposure Index Builder <<<
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

# >>> CELL 07: CTE Builder <<<
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

# >>> CELL 08: Medication Flagging <<<
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

# >>> CELL 09: Comorbidity Flagging <<<
# ---------- COMORBIDITY FLAGGING (kept) ----------
add_comorbidities <- function(cohort_df, dataset = Sys.getenv("WORKSPACE_CDR")) {
  if (nrow(cohort_df) == 0) return(cohort_df)
  if ("person_id" %in% names(cohort_df)) cohort_df <- cohort_df %>% mutate(person_id = as.character(person_id))
  if ("index_date" %in% names(cohort_df)) cohort_df <- cohort_df %>% mutate(index_date = as.Date(index_date))
  comorbidity_sets <- list(
    myocardial_infarction = list(include_std = integer(0), include_src = c(44832372, 1569126, 1569130, 35207702, 44820864), exclude_std = integer(0), exclude_src = integer(0)),
    congestive_heart_failure = list(include_std = integer(0), include_src = c(44833557, 44819695, 44823110, 44824235, 35207762, 44819693, 35207755, 44824250, 35207673, 44820856, 44819710, 44825438, 35207763, 1569178, 44821955, 35207765, 35207674, 44831230, 35207760, 44819696, 44823116, 35207761, 44819692, 35207704, 44835939, 35207764, 44823108, 35207669), exclude_std = integer(0), exclude_src = integer(0)),
    peripheral_vascular_disease = list(include_std = integer(0), include_src = c(44835640, 35208277, 35207882, 44826657, 44834745, 35207869, 35208278, 44834813, 35225419, 44820429, 44830159, 44825446, 44819723, 35207881, 35207860, 44822297, 1569271, 35207883, 1569326, 44821962, 35208276, 1569321, 44826646, 44832392, 44835959, 44834746, 35207859, 1576287), exclude_std = integer(0), exclude_src = integer(0)),
    cerebrovascular_disease = list(include_std = integer(0), include_src = c(1568361, 1569184, 1569190, 1569193, 1568726, 1569191, 44831252, 1568360, 44835952, 1569225, 44820872, 44820873, 44820875, 1568725, 1569227, 44835946, 1568727, 1569218, 44824253, 44830001, 44832388, 44835947, 1569221), exclude_std = integer(0), exclude_src = integer(0)),
    dementia = list(include_std = integer(0), include_src = c(44831122, 44831078, 1568088, 44820749, 1568087, 44821814, 35207116, 35207328, 1568293, 35207360, 45553736, 44831079, 35211390, 45890911, 35207114, 44820073, 35207121, 44827645, 35207361, 44820708, 44835772, 45586320, 1568295, 35207511, 44833397, 44826538, 35207115, 35207118, 44821813, 44826537), exclude_std = integer(0), exclude_src = integer(0)),
    chronic_pulmonary_disease = list(include_std = integer(0), include_src = c(44829014, 35208065, 44820888, 44827823, 1569496, 35208036, 1569495, 44821987, 44834771, 44829011, 1569492, 1569485, 44820887, 44827824, 1569493, 35208056, 35208063, 1569486, 44826682, 44823145, 1569487, 35208037, 35208013, 44835986, 1569488, 35208026, 44829013, 35208017, 44835981, 1569494, 44820892, 44827825, 44826681, 35208027), exclude_std = integer(0), exclude_src = integer(0)),
    rheumatic_disease = list(include_std = integer(0), include_src = c(44831518, 35208841, 44831487, 35208832, 44836170, 1570619, 44832613, 44831474, 1570614, 35208834, 44837334, 1570612, 44828021, 44833584, 44834963, 1569966, 1570039, 44821123, 44819941, 35208820), exclude_std = integer(0), exclude_src = integer(0)),
    peptic_ulcer_disease = list(include_std = integer(0), include_src = c(44830146, 1569564, 1569565, 44825505, 1569563, 44824313, 44833634, 1569562), exclude_std = integer(0), exclude_src = integer(0)),
    mild_liver_disease = list(include_std = integer(0), include_src = c(35208332, 44835626, 44821690, 44833242, 35208363, 35208359, 35208368, 35208338, 35208335, 44829751, 44829754, 44832480, 44824336, 35208367, 44831320, 44826726, 1569680, 1569670, 1569671, 44830943, 44829749, 35208331, 44819379, 1567376, 44821549, 35208330, 1569675, 35208362, 44825528, 35208361, 44819803, 1569681, 35225408), exclude_std = integer(0), exclude_src = integer(0)),
    renal_mild_or_moderate = list(include_std = integer(0), include_src = c(35209277, 44836035, 44819695, 1571474, 44832368, 1571472, 44821546, 35225404, 44837191, 35209275, 44830173, 44835923, 44832369, 35207673, 44827888, 35209279, 35207672, 44835922, 44819694, 44825538, 44819696, 44823191, 44820970, 35209276, 45543164, 44819692, 35209274, 44835924), exclude_std = integer(0), exclude_src = integer(0)),
    hemiplegia_or_paraplegia = list(include_std = integer(0), include_src = c(44823020, 35207306, 1568408, 1568415, 44832266, 35207479, 35207480, 1568412, 35207481, 44832275, 44828872, 35207319), exclude_std = integer(0), exclude_src = integer(0)),
    any_malignancy = list(include_std = integer(0), include_src = c(44825213, 1567501, 1567466, 1567641, 1567463, 1567492, 1567651, 44827551, 1567471, 44833301, 44822887, 1567478, 44826420, 1567462, 44828729, 44827586, 1567699, 1567461, 44824032, 44832117, 44834484, 1567470, 1567569, 35206056, 1567534, 1567674, 44830971, 44819425, 44829811, 1567472, 44834492, 44819434, 44833294, 1567537, 44829849, 44832129, 44822885, 1567476, 1567689, 44820621, 1567530, 44834486, 44829790, 44828731, 35206153, 1567493, 35206140, 44829795, 44831017, 1567573, 44835670, 1567482, 44835663, 44827582, 44832131, 44827557, 44825199, 44819422, 1567469, 35206253, 1567474, 1567633, 1567567, 44829823, 1567659, 1567668, 35206101, 44836833, 1567479, 44828781, 1567468, 35206266, 1567485, 1567502, 1567473, 44829799, 44834480, 1567705, 1567578, 44821758, 44829803, 1567715, 44825214, 44821735, 1567574, 44826411, 44826385, 1567480, 1567484, 1567568, 44824020, 44827567, 44829813, 35206185, 44830966, 35206080, 1567464, 1567465, 1567483, 1567529, 45561836, 44829815, 44825232, 1567680, 35206141, 1567528, 1567615, 1567711, 44822869, 44835671, 44833300, 44824034, 1567494, 44833283, 44820606, 1567475, 44833324, 1567481, 44827548, 44835661, 35206242, 1567566, 44826409, 35206183, 1567486, 1567565, 1567533, 1567477, 35206260, 44819441, 1567675), exclude_std = integer(0), exclude_src = integer(0)),
    moderate_or_severe_liver_disease = list(include_std = integer(0), include_src = c(44834823, 44822040, 35207901, 1569401, 1569679, 1569678, 44837179, 1569674, 44828998, 44825455, 44825529, 35208366, 44832406, 1569672, 35208364, 35208365), exclude_std = integer(0), exclude_src = integer(0)),
    renal_severe = list(include_std = integer(0), include_src = c(44829649, 44835497, 44824342, 35209280, 44824235, 44831231, 44833130, 45548653, 44828971, 44831232, 44819693, 44825426, 44835496, 44820856, 35209291, 44821578, 44829650, 44833560, 45596188, 35207674, 35209278, 35225436, 44829062, 1576113, 44834280, 44837192, 44821950, 44831947, 35207671, 44827889), exclude_std = integer(0), exclude_src = integer(0)),
    hiv_infection = list(include_std = integer(0), include_src = c(4829737, 35205776), exclude_std = integer(0), exclude_src = integer(0)),
    metastatic_solid_tumor = list(include_std = integer(0), include_src = c(1567618, 44836847, 44836850, 1567623, 45557018, 44819442, 44819443, 1567619, 35206334), exclude_std = integer(0), exclude_src = integer(0)),
    diabetes_with_renal_complications = list(include_std = integer(0), include_src = c(44824074, 1567958, 1567975, 1567909, 1567942, 1567926), exclude_std = integer(0), exclude_src = integer(0)),
    diabetes_with_ophthalmic_complications = list(include_std = integer(0), include_src = c(1567976, 1567910, 1567959, 1567927, 1567943, 44828794), exclude_std = integer(0), exclude_src = integer(0)),
    diabetes_with_neurological_complications = list(include_std = integer(0), include_src = c(1567933, 1567916, 1567949, 1567982, 44827615, 1567965), exclude_std = integer(0), exclude_src = integer(0)),
    hypoglycemia = list(include_std = integer(0), include_src = c(1567955, 1567939, 1567922, 35206889, 44831048, 44820686, 1567988, 35206888, 1567971, 35206887), exclude_std = integer(0), exclude_src = integer(0)),
    hyperglycemic_emergency = list(include_std = integer(0), include_src = c(44831044, 1567973, 1567924, 44826459, 1567907, 1567957), exclude_std = integer(0), exclude_src = integer(0)),
    t1d = list(include_std = integer(0), include_src = c(44820682, 44819504, 44832192, 44822934, 44829881, 44833368, 44821787, 44824071, 44832190, 44819501, 44819502, 44825264, 44836918, 44822936, 44822935, 1567940, 44832191, 44820683, 44834549, 44820684, 44831046), exclude_std = integer(0), exclude_src = integer(0)),
    t2d = list(include_std = integer(0), include_src = c(44832194, 44828795, 44832193, 44836915, 44826461, 1567956, 44833366, 44829879, 44829882, 44824073, 44824072, 44831045, 44831047, 44836914, 44833367, 44827616, 44827617, 44836916, 44829878, 44826460, 44819500), exclude_std = integer(0), exclude_src = integer(0)),
    otherdm = list(include_std = integer(0), include_src = c(1567906, 1567923, 1567972, 44832187), exclude_std = integer(0), exclude_src = integer(0))
  )
  valid_rows <- cohort_df %>% filter(!is.na(person_id) & !is.na(index_date))
  if (nrow(valid_rows) == 0) {
    cohort_cte_sql <- ""
  } else {
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

# >>> CELL 10: Baseline Measurements <<<
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
    index_year_numeric = as.numeric(as.character(index_year)),
    index_year_grouped = factor(case_when(
      index_year_numeric %in% c(2018, 2019, 2023) ~ "Non-COVID",
      index_year_numeric %in% 2020:2022 ~ "COVID",
      TRUE ~ "Other"
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
      levels = c("Underweight", "Normal", "Overweight", "Obesity Class I", "Obesity Class II", "Obesity Class III")))
      result
}

# >>> CELL 11: Comparison Group Processor <<<
process_comparison_group <- function(exposure_idx, comparator_idx, base_cohort, exposure_name, comparator_name) {
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
                                 crossover_date = as.Date(NA), stringsAsFactors = FALSE)
  comparator_only_df <- data.frame(person_id = comparator_only_ids, treatment = 0,
                                   index_date = comparator_idx$index_date[match(comparator_only_ids, comparator_idx$person_id)],
                                   crossover_date = as.Date(NA), stringsAsFactors = FALSE)
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
                                index_date = as.Date(character(0)), crossover_date = as.Date(character(0)), stringsAsFactors = FALSE)
  }
  cohort_pair <- bind_rows(exposure_only_df, comparator_only_df, sequential_df)
  if ("index_date" %in% names(base_cohort)) base_cohort <- base_cohort %>% select(-index_date)
  cohort_df <- cohort_pair %>% left_join(base_cohort, by = "person_id")
  cohort_df <- cohort_df %>%
    filter(!is.na(t2d_excluded_t1d_start_date) & t2d_excluded_t1d_start_date <= index_date) %>%
    mutate(index_date = as.Date(index_date),
           pre_index_outcome = ifelse(!is.na(epilepsy_or_seizure_without_hypoglycemic_seizure_start_date) & epilepsy_or_seizure_without_hypoglycemic_seizure_start_date < index_date, 1, 0)) %>%
    filter(pre_index_outcome == 0) %>%
    filter(is.na(crossover_date))
  list(temporal_stats = list(final_cohort = nrow(cohort_df)), cohort_df = cohort_df)
}

# >>> CELL 12: Load Cohort CSVs from Bucket <<<
my_bucket <- Sys.getenv('WORKSPACE_BUCKET')
load_from_bucket <- function(file_name) {
  system(paste0("gsutil cp ", my_bucket, "/data/", file_name, " ."), intern=T)
  read_csv(file_name, show_col_types = FALSE)
}
files <- c(
  "SEMAGLUTIDE_vs_OTHER_GLPA_semaglutide_epilepsy_seizure_full_cohort_20260216_024506.csv",
  "SEMAGLUTIDE_vs_SGLT2_semaglutide_epilepsy_seizure_full_cohort_20260216_024506.csv",
  "SEMAGLUTIDE_vs_OtherGLD_semaglutide_epilepsy_seizure_full_cohort_20260216_024506.csv",
  "OTHER_GLPA_vs_SGLT2_semaglutide_epilepsy_seizure_full_cohort_20260216_024506.csv",
  "OTHER_GLPA_vs_OtherGLD_semaglutide_epilepsy_seizure_full_cohort_20260216_024506.csv",
  "SGLT2_vs_OtherGLD_semaglutide_epilepsy_seizure_full_cohort_20260216_024506.csv"
)
cohorts <- lapply(files, load_from_bucket)
names(cohorts) <- c("semaglutide_vs_other_glpa_full_cohort","semaglutide_vs_sglt2_full_cohort","semaglutide_vs_othergld_full_cohort","other_glpa_vs_sglt2_full_cohort","other_glpa_vs_othergld_full_cohort","sglt2_vs_othergld_full_cohort")

# Reload with local file check
load_from_bucket <- function(file_name) {
  if (!file.exists(file_name)) stop(paste("File not found:", file_name))
  read_csv(file_name, show_col_types = FALSE)
}
files <- c(
  "SEMAGLUTIDE_vs_OTHER_GLPA_semaglutide_epilepsy_seizure_full_cohort_20260216_024506.csv",
  "SEMAGLUTIDE_vs_SGLT2_semaglutide_epilepsy_seizure_full_cohort_20260216_024506.csv",
  "SEMAGLUTIDE_vs_OtherGLD_semaglutide_epilepsy_seizure_full_cohort_20260216_024506.csv",
  "OTHER_GLPA_vs_SGLT2_semaglutide_epilepsy_seizure_full_cohort_20260216_024506.csv",
  "OTHER_GLPA_vs_OtherGLD_semaglutide_epilepsy_seizure_full_cohort_20260216_024506.csv",
  "SGLT2_vs_OtherGLD_semaglutide_epilepsy_seizure_full_cohort_20260216_024506.csv"
)
cohorts <- lapply(files, load_from_bucket)
names(cohorts) <- c("semaglutide_vs_other_glpa_full_cohort","semaglutide_vs_sglt2_full_cohort","semaglutide_vs_othergld_full_cohort","other_glpa_vs_sglt2_full_cohort","other_glpa_vs_othergld_full_cohort","sglt2_vs_othergld_full_cohort")

# >>> CELL 13: Variable Labels & Baseline Table Vars <<<
my_vars <- c("age","sex_cat","raceethnicity_cat","baseline_bmi","baseline_bmi_category","baseline_hba1c","index_year","index_year_grouped","income","education","smoking","insurance_category","alcohol_category","SEMAGLUTIDE","OTHER_GLPA","Biguanide","TZD","Insulin","SGLT2i","DPP4i","SU","Anticoagulant","Antiplatelet","Statin","Ezetimibe","RAAS","Diuretic","MRA","BB","CCB","OtherHTN","myocardial_infarction","congestive_heart_failure","peripheral_vascular_disease","cerebrovascular_disease","dementia","chronic_pulmonary_disease","rheumatic_disease","peptic_ulcer_disease","hemiplegia_or_paraplegia","hiv_infection","hypoglycemia","hyperglycemic_emergency","renal_disease_severity","liver_disease_severity","malignancy_status","diabetes_with_ophthalmic_complications","diabetes_with_neurological_complications","t1d","t2d","otherdm")
categorical_vars <- my_vars[!my_vars %in% c("age","baseline_bmi","baseline_hba1c")]

apply_variable_labels <- function(data) {
  data %>% mutate(
    sex_cat = factor(case_when(sex_cat==0~"Male",sex_cat==1~"Female",sex_cat==999~"Missing",TRUE~"Missing"), levels=c("Male","Female","Missing")),
    raceethnicity_cat = factor(case_when(raceethnicity_cat==0~"Non-Hispanic White",raceethnicity_cat==1~"Non-Hispanic Black",raceethnicity_cat==2~"Hispanic",raceethnicity_cat==3~"Other",raceethnicity_cat==999~"Missing",TRUE~"Missing"), levels=c("Non-Hispanic White","Non-Hispanic Black","Hispanic","Other","Missing")),
    income = factor(case_when(income==0~"Low (<50k)",income==1~"Middle (50k-100k)",income==2~"High (>100k)",income==999~"Missing",TRUE~"Missing"), levels=c("Low (<50k)","Middle (50k-100k)","High (>100k)","Missing")),
    education = factor(case_when(education==0~"High School or Less",education==1~"Some College",education==2~"Advanced Degree",education==999~"Missing",TRUE~"Missing"), levels=c("High School or Less","Some College","Advanced Degree","Missing")),
    smoking = factor(case_when(smoking==0~"Never",smoking==1~"Former",smoking==2~"Current",smoking==999~"Missing",TRUE~"Missing"), levels=c("Never","Former","Current","Missing")),
    insurance_category = factor(case_when(insurance_category==0~"None",insurance_category==1~"Public",insurance_category==2~"Private",insurance_category==999~"Missing",TRUE~"Missing"), levels=c("None","Public","Private","Missing")),
    alcohol_category = factor(case_when(alcohol_category==0~"Low Risk",alcohol_category==1~"Increased Risk",alcohol_category==2~"High Risk",alcohol_category==3~"Dependent",alcohol_category==999~"Missing",TRUE~"Missing"), levels=c("Low Risk","Increased Risk","High Risk","Dependent","Missing")),
    across(c("SEMAGLUTIDE","OTHER_GLPA","Biguanide","TZD","Insulin","SGLT2i","DPP4i","SU","Anticoagulant","Antiplatelet","Statin","Ezetimibe","RAAS","Diuretic","MRA","BB","CCB","OtherHTN","myocardial_infarction","congestive_heart_failure","peripheral_vascular_disease","cerebrovascular_disease","dementia","chronic_pulmonary_disease","rheumatic_disease","peptic_ulcer_disease","hiv_infection","hypoglycemia","hyperglycemic_emergency","diabetes_with_ophthalmic_complications","diabetes_with_neurological_complications","t1d","t2d","otherdm"), ~ factor(case_when(.==0~"No",.==1~"Yes",TRUE~"No"), levels=c("No","Yes")))
  )
}

# >>> CELL 14: Normalize Types & Ingredient Helpers <<<
cohorts <- lapply(cohorts, function(df) { df %>% dplyr::mutate(person_id=as.character(person_id), index_date=as.Date(index_date), treatment=as.integer(as.character(treatment))) })

get_index_ingredient_for_arm <- function(cohort_df, ingredient_names, arm_value, out_col, dataset=Sys.getenv("WORKSPACE_CDR")) {
  cohort_df <- cohort_df %>% dplyr::mutate(person_id=as.character(person_id), index_date=as.Date(index_date))
  cdf <- cohort_df %>% dplyr::filter(.data$treatment==arm_value, !is.na(person_id), !is.na(index_date))
  if (!out_col %in% names(cohort_df)) { cohort_df[[out_col]] <- NA_character_ } else if (is.null(cohort_df[[out_col]])) { cohort_df[[out_col]] <- rep(NA_character_, nrow(cohort_df)) }
  if (nrow(cdf)==0) return(cohort_df)
  ing_names_sql <- paste0("'", tolower(ingredient_names), "'", collapse=", ")
  cohort_rows <- purrr::map2_chr(cdf$person_id, cdf$index_date, ~ glue::glue("STRUCT({.x} AS person_id, DATE('{format(.y, '%Y-%m-%d')}') AS index_date)"))
  cohort_cte_sql <- paste(cohort_rows, collapse=", ")
  sql <- glue::glue("
    WITH cohort AS (SELECT * FROM UNNEST(ARRAY<STRUCT<person_id INT64, index_date DATE>>[{cohort_cte_sql}])),
    idx_day_exposures AS (
      SELECT t.person_id, ca.ancestor_concept_id AS ingredient_id
      FROM cohort t JOIN `{dataset}.drug_exposure` de ON de.person_id=t.person_id AND de.drug_exposure_start_date=t.index_date
      JOIN `{dataset}.concept_ancestor` ca ON de.drug_concept_id=ca.descendant_concept_id
      JOIN `{dataset}.concept` anc ON anc.concept_id=ca.ancestor_concept_id
      WHERE anc.concept_class_id='Ingredient' AND LOWER(anc.concept_name) IN ({ing_names_sql})
    ),
    tie_broken AS (SELECT person_id, ingredient_id FROM (SELECT person_id, ingredient_id, ROW_NUMBER() OVER (PARTITION BY person_id ORDER BY ingredient_id ASC) AS rn FROM idx_day_exposures) WHERE rn=1)
    SELECT tb.person_id, c.concept_name AS index_ingredient_name FROM tie_broken tb JOIN `{dataset}.concept` c ON c.concept_id=tb.ingredient_id
  ")
  ix <- download_data(sql) %>% dplyr::mutate(person_id=as.character(person_id))
  out <- cohort_df %>% dplyr::left_join(ix, by="person_id")
  rows <- which(out$treatment==arm_value)
  if ("index_ingredient_name" %in% names(out)) { out[[out_col]][rows] <- out$index_ingredient_name[rows]; out$index_ingredient_name <- NULL }
  out
}

INGR_SEMA     <- c("semaglutide")
INGR_OTHERGLP <- c("exenatide","liraglutide","albiglutide","dulaglutide","lixisenatide")
INGR_SGLT2    <- c("canagliflozin","empagliflozin","dapagliflozin","ertugliflozin")
INGR_OTHERGLD <- c("glimepiride","glipizide","glyburide","alogliptin","linagliptin","sitagliptin","saxagliptin","pioglitazone","rosiglitazone")
get_exp_ing_by_comp <- function(comp_name) { switch(comp_name, "SEMAGLUTIDE vs OTHER_GLPA"=INGR_SEMA, "SEMAGLUTIDE vs SGLT2"=INGR_SEMA, "SEMAGLUTIDE vs OtherGLD"=INGR_SEMA, "OTHER_GLPA vs SGLT2"=INGR_OTHERGLP, "OTHER_GLPA vs OtherGLD"=INGR_OTHERGLP, "SGLT2 vs OtherGLD"=INGR_SGLT2, NULL) }
get_comp_ing_by_comp <- function(comp_name) { switch(comp_name, "SEMAGLUTIDE vs OTHER_GLPA"=INGR_OTHERGLP, "SEMAGLUTIDE vs SGLT2"=INGR_SGLT2, "SEMAGLUTIDE vs OtherGLD"=INGR_OTHERGLD, "OTHER_GLPA vs SGLT2"=INGR_SGLT2, "OTHER_GLPA vs OtherGLD"=INGR_OTHERGLD, "SGLT2 vs OtherGLD"=INGR_OTHERGLD, NULL) }
n_levels <- function(x) length(unique(stats::na.omit(x)))

tag_and_make_tables_for_both_arms <- function(cohort_df, comp_name, analysis_prefix) {
  message("Tagging index ingredients & building tables for BOTH arms: ", comp_name)
  cohort_df <- cohort_df %>% dplyr::mutate(person_id=as.character(person_id), index_date=as.Date(index_date), treatment=as.integer(as.character(treatment)))
  exp_names <- get_exp_ing_by_comp(comp_name); comp_names <- get_comp_ing_by_comp(comp_name)
  if (!is.null(exp_names)) { cohort_df <- get_index_ingredient_for_arm(cohort_df, exp_names, arm_value=1, out_col="index_ingredient_name_treat") } else { cohort_df$index_ingredient_name_treat <- NA_character_ }
  if (!is.null(comp_names)) { cohort_df <- get_index_ingredient_for_arm(cohort_df, comp_names, arm_value=0, out_col="index_ingredient_name_comp") } else { cohort_df$index_ingredient_name_comp <- NA_character_ }
  cohort_df <- apply_variable_labels(cohort_df)
  treat_df <- cohort_df %>% dplyr::filter(treatment==1L)
  if (!is.null(exp_names) && length(exp_names)==1 && (all(is.na(treat_df$index_ingredient_name_treat)) || !"index_ingredient_name_treat" %in% names(treat_df))) { treat_df$index_ingredient_name_treat <- tolower(exp_names[[1]]) }
  n_treat <- treat_df %>% dplyr::mutate(index_ingredient_name_treat=tolower(index_ingredient_name_treat)) %>% dplyr::count(index_ingredient_name_treat, name="N") %>% dplyr::arrange(dplyr::desc(N))
  # All of Us small-cell suppression: mask ingredient counts < 20 per Data Dissemination Policy
  n_treat$N[n_treat$N < 20] <- NA_integer_
  n_treat$index_ingredient_name_treat[is.na(n_treat$N)] <- "<suppressed>"
  write.csv(n_treat, file=sprintf("N_by_ingredient_TREAT_%s_%s_%s.csv", gsub("[^A-Za-z0-9]+","_",comp_name), analysis_prefix, format(Sys.time(),"%Y%m%d_%H%M%S")), row.names=FALSE)
  if (n_levels(treat_df$index_ingredient_name_treat)>=2) { tbl <- CreateTableOne(vars=my_vars, strata="index_ingredient_name_treat", data=treat_df, factorVars=categorical_vars); write.csv(as.data.frame(print(tbl, showAllLevels=TRUE, printToggle=FALSE, noSpaces=TRUE)), file=sprintf("demographics_by_ingredient_TREAT_%s_%s_%s.csv", gsub("[^A-Za-z0-9]+","_",comp_name), analysis_prefix, format(Sys.time(),"%Y%m%d_%H%M%S")), row.names=TRUE) }
  tbl_cls <- CreateTableOne(vars=my_vars, data=treat_df, factorVars=categorical_vars)
  write.csv(as.data.frame(print(tbl_cls, showAllLevels=TRUE, printToggle=FALSE, noSpaces=TRUE)), file=sprintf("demographics_by_class_TREAT_%s_%s_%s.csv", gsub("[^A-Za-z0-9]+","_",comp_name), analysis_prefix, format(Sys.time(),"%Y%m%d_%H%M%S")), row.names=TRUE)
  comp_df <- cohort_df %>% dplyr::filter(treatment==0L)
  if (!is.null(comp_names) && length(comp_names)==1 && (all(is.na(comp_df$index_ingredient_name_comp)) || !"index_ingredient_name_comp" %in% names(comp_df))) { comp_df$index_ingredient_name_comp <- tolower(comp_names[[1]]) }
  n_comp <- comp_df %>% dplyr::mutate(index_ingredient_name_comp=tolower(index_ingredient_name_comp)) %>% dplyr::count(index_ingredient_name_comp, name="N") %>% dplyr::arrange(dplyr::desc(N))
  # All of Us small-cell suppression: mask ingredient counts < 20 per Data Dissemination Policy
  n_comp$N[n_comp$N < 20] <- NA_integer_
  n_comp$index_ingredient_name_comp[is.na(n_comp$N)] <- "<suppressed>"
  write.csv(n_comp, file=sprintf("N_by_ingredient_COMP_%s_%s_%s.csv", gsub("[^A-Za-z0-9]+","_",comp_name), analysis_prefix, format(Sys.time(),"%Y%m%d_%H%M%S")), row.names=FALSE)
  cohort_df
}

# >>> CELL 15: Binary Helpers & Logging (re-init) <<<
is_binary_like <- function(x) {
  if (is.logical(x)) return(TRUE)
  if (is.numeric(x)) { ux <- unique(na.omit(x)); return(length(ux)<=2 && all(ux %in% c(0,1))) }
  if (is.factor(x)||is.character(x)) { xc <- tolower(trimws(as.character(x))); ok <- xc %in% c("0","1","no","yes","n","y","false","true","t","f"); return(all(ok|is.na(xc))) }
  FALSE
}
to_yesno <- function(x) {
  if (is.logical(x)) return(factor(ifelse(x,"Yes","No"), levels=c("No","Yes")))
  if (is.numeric(x)) return(factor(ifelse(x==1,"Yes","No"), levels=c("No","Yes")))
  xc <- tolower(trimws(as.character(x)))
  out <- ifelse(xc %in% c("1","yes","y","true","t"), "Yes", ifelse(xc %in% c("0","no","n","false","f"), "No", NA))
  factor(out, levels=c("No","Yes"))
}

# Re-init logging for analysis phase
init_output_tracker <- function() { options(semaglutide_output_tracker = data.frame(timestamp=character(), script=character(), type=character(), message=character(), value=character(), stringsAsFactors=FALSE)) }
log_output <- function(type, message, value="", script="") {
  if (is.null(getOption("semaglutide_output_tracker"))) init_output_tracker()
  tracker <- getOption("semaglutide_output_tracker")
  new_row <- data.frame(timestamp=format(Sys.time(),"%Y-%m-%d %H:%M:%S"), script=script, type=type, message=message, value=as.character(value), stringsAsFactors=FALSE)
  options(semaglutide_output_tracker=rbind(tracker, new_row))
  if (type %in% c("PROGRESS","ERROR","WARNING")) cat(sprintf("[%s] %s: %s\n", type, message, value))
}
save_output_log <- function(filename=NULL) {
  if (is.null(filename)) filename <- sprintf("output_log_%s.csv", format(Sys.time(),"%Y%m%d_%H%M%S"))
  tracker <- getOption("semaglutide_output_tracker")
  if (!is.null(tracker) && nrow(tracker)>0) { write.csv(tracker, filename, row.names=FALSE); log_output("PROGRESS","Output log saved",filename); return(filename) }
  return(NULL)
}

# >>> CELL 16: followup_and_event <<<
followup_and_event <- function(df, outcome_var, data_cut_date, late_onset=FALSE, early_onset=FALSE, verbose=FALSE, comp_name=NULL, exclude_preindex_any_age=FALSE) {
  result <- data.frame(df)
  result$age_at_event <- result$age + as.numeric(as.Date(result[[outcome_var]]) - as.Date(result$index_date)) / 365.25
  result$is_qualified_event <- (!is.na(result[[outcome_var]]) & if (late_onset) { result$age_at_event >= 60 } else if (early_onset) { result$age_at_event < 60 } else { TRUE })
  raw_outcome_date <- try(as.Date(result[[outcome_var]]), silent=TRUE)
  idx_date_vec <- try(as.Date(result$index_date), silent=TRUE)
  if (!inherits(raw_outcome_date,"try-error") && !inherits(idx_date_vec,"try-error")) {
    preidx_mask_detect <- !is.na(raw_outcome_date) & (raw_outcome_date <= idx_date_vec)
    log_output("PROGRESS", sprintf("Detected preindex/same-day outcomes (raw) [%s]", ifelse(late_onset,"late_onset","all_ages")), sum(preidx_mask_detect), "followup")
    if (isTRUE(late_onset)) {
      n_preidx <- sum(preidx_mask_detect); n_before <- nrow(result)
      if (n_preidx > 0) {
        try({ out_name <- sprintf("preindex_outcomes_excluded_late_onset_%s_%s.csv", outcome_var, format(Sys.time(),"%Y%m%d_%H%M%S")); safe_cols <- intersect(c("person_id","index_date",outcome_var,"age","treatment"), names(result)); utils::write.csv(result[preidx_mask_detect, safe_cols, drop=FALSE], file=out_name, row.names=FALSE) }, silent=TRUE)
        result <- result[!preidx_mask_detect, , drop=FALSE]
        log_output("PROGRESS","Late-onset preindex exclusion applied", sprintf("n_preidx=%d; n_before=%d; n_after=%d", n_preidx, n_before, nrow(result)), "followup")
      }
    }
  }
  result$outcome_date <- as.Date(ifelse(result$is_qualified_event, as.character(result[[outcome_var]]), NA), origin="1970-01-01")
  result$censor_date <- pmin(as.Date(result$EHRmaxDT_min), as.Date(data_cut_date), na.rm=TRUE)
  result$end_fu_date <- pmin(result$outcome_date, result$censor_date, result$death_date, na.rm=TRUE)
  result$time <- as.numeric(result$end_fu_date - as.Date(result$index_date))
  result$event <- as.integer(!is.na(result$outcome_date) & result$outcome_date == result$end_fu_date)
  n_pre_time <- nrow(result); n_nonpos <- sum(!is.na(result$time) & result$time <= 0)
  if (!is.na(n_nonpos) && n_nonpos > 0) {
    diag_df <- result; diag_df$analysis_group <- ifelse(late_onset,"late_onset",ifelse(early_onset,"early_onset","all_ages"))
    diag_df$time_nonpos <- !is.na(diag_df$time) & diag_df$time <= 0
    diag_df$reason_nonpos <- NA_character_; diag_df$had_outcome_raw <- !is.na(diag_df[[outcome_var]])
    diag_df$age_ok <- if (late_onset) { diag_df$age_at_event >= 60 } else if (early_onset) { diag_df$age_at_event < 60 } else { TRUE }
    diag_df$end_is_outcome <- !is.na(diag_df$outcome_date) & (diag_df$end_fu_date == diag_df$outcome_date)
    diag_df$end_is_censor <- (!is.na(diag_df$censor_date)) & (diag_df$end_fu_date == diag_df$censor_date)
    diag_df$end_is_death <- (!is.na(diag_df$death_date)) & (diag_df$end_fu_date == diag_df$death_date)
    diag_df$reason_nonpos[diag_df$time_nonpos & diag_df$end_is_outcome] <- "preindex_or_same_day_outcome"
    diag_df$reason_nonpos[diag_df$time_nonpos & diag_df$end_is_censor & is.na(diag_df$reason_nonpos)] <- "preindex_or_same_day_censor"
    diag_df$reason_nonpos[diag_df$time_nonpos & diag_df$end_is_death & is.na(diag_df$reason_nonpos)] <- "preindex_or_same_day_death"
    diag_df$reason_nonpos[diag_df$time_nonpos & diag_df$had_outcome_raw & !diag_df$age_ok & is.na(diag_df$reason_nonpos)] <- "event_failed_age_criterion"
    diag_df$reason_nonpos[diag_df$time_nonpos & is.na(diag_df$reason_nonpos) & !diag_df$had_outcome_raw] <- "no_event_recorded"
    diag_df$reason_nonpos[diag_df$time_nonpos & is.na(diag_df$reason_nonpos)] <- "other_preindex_issue"
    by_reason <- diag_df %>% dplyr::filter(time_nonpos) %>% dplyr::count(reason_nonpos, name="n")
    by_treatment <- diag_df %>% dplyr::filter(time_nonpos) %>% dplyr::count(treatment, name="n")
    msg_reason <- paste(paste0(by_reason$reason_nonpos,"=",by_reason$n), collapse=", ")
    msg_trt <- paste(paste0("trt",by_treatment$treatment,"=",by_treatment$n), collapse=", ")
    log_output("PROGRESS","Non-positive follow-up filtered", sprintf("Removed %d of %d; reasons: %s; by treatment: %s", n_nonpos, n_pre_time, msg_reason, msg_trt), "followup")
    out_name <- sprintf("nonpos_followup_%s_%s_%s.csv", diag_df$analysis_group[1], outcome_var, format(Sys.time(),"%Y%m%d_%H%M%S"))
    safe_cols <- intersect(c("person_id","index_date","outcome_date","censor_date","death_date","time","event","treatment","age","age_at_event","is_qualified_event","analysis_group","end_is_outcome","end_is_censor","end_is_death","reason_nonpos"), names(diag_df))
    try({ utils::write.csv(diag_df[diag_df$time_nonpos, safe_cols, drop=FALSE], file=out_name, row.names=FALSE) }, silent=TRUE)
  }
  result <- result %>% dplyr::filter(time > 0)
  result$event_time <- result$time
  log_output("INFO","followup_and_event complete", sprintf("%d events from %d patients", sum(result$event), nrow(result)), "followup_and_event")
  return(result)
}

# >>> CELL 17: IPTW & Cox Functions <<<
run_ipwt_and_cox <- function(df, exclude_vars=NULL, trim_threshold=0.01, all_ps_vars=NULL, verbose=NULL, stabilize=FALSE) {
  all_covariates <- if (is.null(all_ps_vars)) get_standard_config()$all_ps_vars else all_ps_vars
  rhs_vars <- setdiff(all_covariates, exclude_vars)
  keep_vars <- rhs_vars[sapply(df[rhs_vars], function(x) length(unique(na.omit(x)))>1)]
  ps_form <- reformulate(keep_vars, response="treatment"); model_vars <- c("treatment", keep_vars)
  if (isTRUE(verbose)) { na_counts <- sapply(df[model_vars], function(x) sum(is.na(x))); if (any(na_counts>0)) { na_msg <- paste(paste0(names(na_counts)[na_counts>0],"=",na_counts[na_counts>0]), collapse=", "); log_output("PROGRESS","Missing PS covariates before complete-case filter", na_msg, "ipwt") } }
  complete_rows <- complete.cases(df[model_vars]); df_complete <- df[complete_rows, ]
  if (verbose) { n_removed <- nrow(df)-nrow(df_complete); log_output("INFO","IPWT sample size", sprintf("%d complete cases from %d total", nrow(df_complete), nrow(df))); log_output("PROGRESS","Complete-case filtering removed rows", sprintf("%d rows removed due to missing PS covariates", n_removed), "ipwt") }
  ps_model <- glm(ps_form, data=df_complete, family=binomial()); df_complete$ps <- predict(ps_model, type="response")
  if (isTRUE(stabilize)) { p_treat <- mean(df_complete$treatment, na.rm=TRUE); df_complete$ipw <- ifelse(df_complete$treatment==1, p_treat/df_complete$ps, (1-p_treat)/(1-df_complete$ps)) } else { df_complete$ipw <- ifelse(df_complete$treatment==1, 1/df_complete$ps, 1/(1-df_complete$ps)) }
  if (!is.null(trim_threshold)) { lower <- quantile(df_complete$ipw, trim_threshold); upper <- quantile(df_complete$ipw, 1-trim_threshold); df_complete$ipw_trimmed <- pmin(pmax(df_complete$ipw, lower), upper) } else { df_complete$ipw_trimmed <- df_complete$ipw }
  df_complete$ipw_std <- df_complete$ipw_trimmed * nrow(df_complete) / sum(df_complete$ipw_trimmed)
  balance_before <- calculate_balance_metrics(df_complete, keep_vars, NULL, verbose=FALSE)
  balance_after <- calculate_balance_metrics(df_complete, keep_vars, df_complete$ipw_std, verbose=FALSE)
  weighted_cox <- coxph(Surv(event_time, event) ~ treatment, data=df_complete, weights=ipw_std)
  return(list(balance_before=balance_before, balance_after=balance_after, cox=summary(weighted_cox), cohort=df_complete, ps_model=ps_model, dropped_covariates=setdiff(rhs_vars, keep_vars), n_dropped=nrow(df)-nrow(df_complete)))
}

# >>> CELL 18: Balance Metrics <<<
calculate_balance_metrics <- function(df, vars, weights, verbose=FALSE) {
  balance_df <- data.frame(variable=character(), std_diff=numeric(), stringsAsFactors=FALSE)
  multilevel_overall_vars <- c("sex_cat","raceethnicity_cat","income","education","insurance_category","smoking","alcohol_category","index_year_grouped","renal_disease_severity","liver_disease_severity","malignancy_status")
  for (var in vars) {
    if (!var %in% names(df)) next
    if (is.factor(df[[var]]) && var %in% multilevel_overall_vars) {
      all_levels <- levels(df[[var]]); level_smds <- numeric()
      for (lvl in all_levels) {
        var_binary <- as.numeric(df[[var]]==lvl)
        if (is.null(weights)) { mean_t1 <- mean(var_binary[df$treatment==1], na.rm=TRUE); mean_t0 <- mean(var_binary[df$treatment==0], na.rm=TRUE); var_t1 <- var(var_binary[df$treatment==1], na.rm=TRUE); var_t0 <- var(var_binary[df$treatment==0], na.rm=TRUE) } else { w_t1 <- weights[df$treatment==1]; w_t0 <- weights[df$treatment==0]; mean_t1 <- weighted.mean(var_binary[df$treatment==1], w_t1, na.rm=TRUE); mean_t0 <- weighted.mean(var_binary[df$treatment==0], w_t0, na.rm=TRUE); var_t1 <- sum(w_t1*(var_binary[df$treatment==1]-mean_t1)^2, na.rm=TRUE)/sum(w_t1); var_t0 <- sum(w_t0*(var_binary[df$treatment==0]-mean_t0)^2, na.rm=TRUE)/sum(w_t0) }
        pooled_sd <- sqrt((var_t1+var_t0)/2); level_smd <- if(pooled_sd>0) (mean_t1-mean_t0)/pooled_sd else 0; level_smds <- c(level_smds, level_smd)
      }
      overall_smd <- ifelse(length(level_smds)>0, max(abs(level_smds), na.rm=TRUE), 0)
      if(length(level_smds)>0) { max_idx <- which.max(abs(level_smds)); overall_smd <- overall_smd * sign(level_smds[max_idx]) }
      balance_df <- rbind(balance_df, data.frame(variable=var, std_diff=overall_smd))
    } else if (is.factor(df[[var]])) {
      for (lvl in levels(df[[var]])[-1]) {
        var_binary <- as.numeric(df[[var]]==lvl)
        if (is.null(weights)) { mean_t1 <- mean(var_binary[df$treatment==1], na.rm=TRUE); mean_t0 <- mean(var_binary[df$treatment==0], na.rm=TRUE); var_t1 <- var(var_binary[df$treatment==1], na.rm=TRUE); var_t0 <- var(var_binary[df$treatment==0], na.rm=TRUE) } else { w_t1 <- weights[df$treatment==1]; w_t0 <- weights[df$treatment==0]; mean_t1 <- weighted.mean(var_binary[df$treatment==1], w_t1, na.rm=TRUE); mean_t0 <- weighted.mean(var_binary[df$treatment==0], w_t0, na.rm=TRUE); var_t1 <- sum(w_t1*(var_binary[df$treatment==1]-mean_t1)^2, na.rm=TRUE)/sum(w_t1); var_t0 <- sum(w_t0*(var_binary[df$treatment==0]-mean_t0)^2, na.rm=TRUE)/sum(w_t0) }
        pooled_sd <- sqrt((var_t1+var_t0)/2); std_diff <- if(pooled_sd>0) (mean_t1-mean_t0)/pooled_sd else 0
        balance_df <- rbind(balance_df, data.frame(variable=paste0(var, lvl), std_diff=std_diff))
      }
    } else {
      if (is.null(weights)) { mean_t1 <- mean(df[[var]][df$treatment==1], na.rm=TRUE); mean_t0 <- mean(df[[var]][df$treatment==0], na.rm=TRUE); var_t1 <- var(df[[var]][df$treatment==1], na.rm=TRUE); var_t0 <- var(df[[var]][df$treatment==0], na.rm=TRUE) } else { w_t1 <- weights[df$treatment==1]; w_t0 <- weights[df$treatment==0]; mean_t1 <- weighted.mean(df[[var]][df$treatment==1], w_t1, na.rm=TRUE); mean_t0 <- weighted.mean(df[[var]][df$treatment==0], w_t0, na.rm=TRUE); var_t1 <- sum(w_t1*(df[[var]][df$treatment==1]-mean_t1)^2, na.rm=TRUE)/sum(w_t1); var_t0 <- sum(w_t0*(df[[var]][df$treatment==0]-mean_t0)^2, na.rm=TRUE)/sum(w_t0) }
      pooled_sd <- sqrt((var_t1+var_t0)/2); std_diff <- if(pooled_sd>0) (mean_t1-mean_t0)/pooled_sd else 0
      balance_df <- rbind(balance_df, data.frame(variable=var, std_diff=std_diff))
    }
  }
  return(balance_df)
}

# >>> CELL 19: Statistics & Visualization Functions <<<
calculate_incident_rates <- function(cohort_df, weight_var=NULL) {
  if (!is.null(weight_var) && weight_var %in% names(cohort_df)) { treatment_df <- cohort_df %>% filter(treatment==1); control_df <- cohort_df %>% filter(treatment==0); treatment_py <- sum(treatment_df$time*treatment_df[[weight_var]]/365.25); control_py <- sum(control_df$time*control_df[[weight_var]]/365.25); treatment_events <- sum(treatment_df$event*treatment_df[[weight_var]]); control_events <- sum(control_df$event*control_df[[weight_var]]) } else { treatment_df <- cohort_df %>% filter(treatment==1); control_df <- cohort_df %>% filter(treatment==0); treatment_py <- sum(treatment_df$time)/365.25; control_py <- sum(control_df$time)/365.25; treatment_events <- sum(treatment_df$event); control_events <- sum(control_df$event) }
  treatment_rate <- (treatment_events/treatment_py)*1000; control_rate <- (control_events/control_py)*1000; rate_diff <- treatment_rate-control_rate
  treatment_rate_se <- sqrt(treatment_events)/treatment_py*1000; control_rate_se <- sqrt(control_events)/control_py*1000
  rate_diff_se <- sqrt((treatment_events/treatment_py^2)+(control_events/control_py^2))*1000
  return(list(treatment_rate=treatment_rate, treatment_rate_ci=c(treatment_rate-1.96*treatment_rate_se, treatment_rate+1.96*treatment_rate_se), control_rate=control_rate, control_rate_ci=c(control_rate-1.96*control_rate_se, control_rate+1.96*control_rate_se), rate_difference=rate_diff, rate_difference_ci=c(rate_diff-1.96*rate_diff_se, rate_diff+1.96*rate_diff_se), treatment_person_years=treatment_py, control_person_years=control_py))
}
calculate_followup_stats <- function(cohort_df, weight_var=NULL) {
  cohort_df$followup_years <- cohort_df$time/365.25
  if (!is.null(weight_var) && weight_var %in% names(cohort_df)) { treatment_df <- cohort_df %>% filter(treatment==1); control_df <- cohort_df %>% filter(treatment==0); treatment_mean <- weighted.mean(treatment_df$followup_years, treatment_df[[weight_var]]); control_mean <- weighted.mean(control_df$followup_years, control_df[[weight_var]]); treatment_sd <- sqrt(sum(treatment_df[[weight_var]]*(treatment_df$followup_years-treatment_mean)^2)/sum(treatment_df[[weight_var]])); control_sd <- sqrt(sum(control_df[[weight_var]]*(control_df$followup_years-control_mean)^2)/sum(control_df[[weight_var]])) } else { treatment_df <- cohort_df %>% filter(treatment==1); control_df <- cohort_df %>% filter(treatment==0); treatment_mean <- mean(treatment_df$followup_years); control_mean <- mean(control_df$followup_years); treatment_sd <- sd(treatment_df$followup_years); control_sd <- sd(control_df$followup_years) }
  return(list(treatment_mean=treatment_mean, treatment_sd=treatment_sd, control_mean=control_mean, control_sd=control_sd))
}
calculate_crude_hr <- function(cohort_df, comparison_name, outcome_label) {
  tryCatch({ crude_cox <- coxph(Surv(event_time, event) ~ treatment, data=cohort_df); hr_crude <- exp(coef(crude_cox)["treatment"]); ci_crude <- exp(confint(crude_cox)["treatment", ]); p_crude <- summary(crude_cox)$coefficients["treatment","Pr(>|z|)"]; followup_stats <- calculate_followup_stats(cohort_df, weight_var=NULL)
    data.frame(comparison=comparison_name, outcome=outcome_label, analysis_type="Crude_Unadjusted", n_total=nrow(cohort_df), n_events=sum(cohort_df$event), n_treatment=sum(cohort_df$treatment==1), n_control=sum(cohort_df$treatment==0), events_treatment=sum(cohort_df$treatment==1 & cohort_df$event==1), events_control=sum(cohort_df$treatment==0 & cohort_df$event==1), followup_years_treatment_mean=followup_stats$treatment_mean, followup_years_treatment_sd=followup_stats$treatment_sd, followup_years_control_mean=followup_stats$control_mean, followup_years_control_sd=followup_stats$control_sd, hr=hr_crude, hr_lower=ci_crude[1], hr_upper=ci_crude[2], p_value=p_crude, stringsAsFactors=FALSE)
  }, error=function(e) { log_output("ERROR","Crude HR calculation failed", e$message, "crude_hr"); return(NULL) })
}

save_analysis_results <- function(results_list, prefix="analysis") {
  timestamp <- format(Sys.time(),"%Y%m%d_%H%M%S"); all_results <- list()
  for (name in names(results_list)) {
    res <- results_list[[name]]
    if (is.null(res)||is.null(res$cox)||is.null(res$crude_hr_results)||is.null(res$incident_rate_results)) { log_output("WARNING","Skipping results saving due to missing data", name, "save_analysis_results"); next }
    split_name <- strsplit(name," - ")[[1]]; comparison_name <- split_name[1]; outcome_name <- split_name[2]
    unweighted_stats <- res$crude_hr_results; unweighted_cohort <- res$cohort
    total_fu_years_unweighted <- unweighted_cohort$time/365.25; total_fu_mean_unweighted <- mean(total_fu_years_unweighted); total_fu_sd_unweighted <- sd(total_fu_years_unweighted)
    weighted_cohort <- res$cohort; inc_rates <- res$incident_rate_results; fu_stats_weighted <- calculate_followup_stats(weighted_cohort, weight_var="ipw_std")
    pseudo_pop_treatment <- sum(weighted_cohort$ipw_std[weighted_cohort$treatment==1]); pseudo_pop_control <- sum(weighted_cohort$ipw_std[weighted_cohort$treatment==0]); pseudo_pop_total <- pseudo_pop_treatment+pseudo_pop_control
    weighted_events_treatment <- sum(weighted_cohort$event[weighted_cohort$treatment==1]*weighted_cohort$ipw_std[weighted_cohort$treatment==1]); weighted_events_control <- sum(weighted_cohort$event[weighted_cohort$treatment==0]*weighted_cohort$ipw_std[weighted_cohort$treatment==0]); weighted_events_total <- weighted_events_treatment+weighted_events_control
    total_py_weighted <- inc_rates$treatment_person_years+inc_rates$control_person_years; total_ir_weighted <- (weighted_events_total/total_py_weighted)*1000; cox_model <- res$cox
    results_df_row <- data.frame(comparison=comparison_name, outcome=outcome_name, unweighted_patients_total=sprintf("%d", unweighted_stats$n_total), unweighted_events_total=sprintf("%d", unweighted_stats$n_events), unweighted_events_patients_treatment=sprintf("%d / %d", unweighted_stats$events_treatment, unweighted_stats$n_treatment), unweighted_events_patients_control=sprintf("%d / %d", unweighted_stats$events_control, unweighted_stats$n_control), unweighted_fu_years_total_mean_sd=sprintf("%.2f (%.2f)", total_fu_mean_unweighted, total_fu_sd_unweighted), hr_crude_ci=sprintf("%.2f (%.2f to %.2f)", unweighted_stats$hr, unweighted_stats$hr_lower, unweighted_stats$hr_upper), hr_crude_p_value=unweighted_stats$p_value, hr_weighted_ci=sprintf("%.2f (%.2f to %.2f)", cox_model$conf.int[1,"exp(coef)"], cox_model$conf.int[1,"lower .95"], cox_model$conf.int[1,"upper .95"]), hr_p_value=cox_model$coefficients[1,"Pr(>|z|)"], stringsAsFactors=FALSE)
    all_results[[name]] <- results_df_row
  }
  if (length(all_results)>0) { final_results_df <- do.call(rbind, all_results); results_file <- sprintf("%s_comprehensive_results_%s.csv", prefix, timestamp); write.csv(final_results_df, results_file, row.names=FALSE); log_output("PROGRESS","Comprehensive results saved", results_file) } else { results_file <- NULL }
  balance_list <- lapply(names(results_list), function(name) { res <- results_list[[name]]; if (!is.null(res)) { before <- res$balance_before; after <- res$balance_after; before$analysis <- name; before$weighting <- "before"; after$analysis <- name; after$weighting <- "after"; rbind(before, after) } })
  if (length(balance_list)>0 && !all(sapply(balance_list, is.null))) { balance_df <- do.call(rbind, balance_list); balance_file <- sprintf("%s_balance_%s.csv", prefix, timestamp); write.csv(balance_df, balance_file, row.names=FALSE) } else { balance_file <- NULL }
  return(list(results=results_file, balance=balance_file))
}

create_balance_plot <- function(balance_before, balance_after, title="SMD Comparison") {
  balance_data <- rbind(balance_before %>% mutate(type="Before"), balance_after %>% mutate(type="After"))
  p <- ggplot(balance_data, aes(x=abs(std_diff), y=reorder(variable, abs(std_diff)), color=type)) + geom_point(size=2, alpha=0.7) + geom_vline(xintercept=0.1, linetype="dashed", color="red") + scale_color_manual(values=c("Before"="#E74C3C","After"="#27AE60")) + labs(x="Absolute Standardized Mean Difference", y="Variable", title=title, color="Weighting") + theme_minimal() + theme(legend.position="bottom")
  return(p)
}

create_survival_plot <- function(ipwt_result, comparison_name, outcome_label, analysis_name, save_plot=TRUE) {
  # Check if survminer is available; skip if not installed
  if (!requireNamespace("survminer", quietly=TRUE)) {
    log_output("WARNING", "survminer not available, skipping survival plot", comparison_name)
    return(NULL)
  }
  library(survminer); library(ggplot2); weighted_df <- ipwt_result$cohort
  treatment_labels <- switch(comparison_name, "SEMAGLUTIDE vs OTHER_GLPA"=c("Other GLP-1 Agonists","Semaglutide"), "SEMAGLUTIDE vs SGLT2"=c("SGLT-2 Inhibitors","Semaglutide"), "SEMAGLUTIDE vs OtherGLD"=c("Other GLDs","Semaglutide"), "OTHER_GLPA vs SGLT2"=c("SGLT-2 Inhibitors","Other GLP-1 Agonists"), "OTHER_GLPA vs OtherGLD"=c("Other GLDs","Other GLP-1 Agonists"), "SGLT2 vs OtherGLD"=c("Other GLDs","SGLT-2 Inhibitors"), c("Comparator","Treatment"))
  weighted_df$event_time_months <- weighted_df$event_time/30.44
  cum_incidence_fit <- survfit(Surv(event_time_months, event) ~ treatment, data=weighted_df, weights=weighted_df$ipw_std)
  surv_plot <- ggsurvplot(fit=cum_incidence_fit, data=weighted_df, fun="event", conf.int=TRUE, risk.table=TRUE, risk.table.col="strata", risk.table.fontsize=3, tables.height=0.3, legend.title="Treatment Group", legend.labs=treatment_labels, xlab="Follow-up Time (months)", ylab="Cumulative Incidence", title=paste("Cumulative Incidence of", outcome_label, "for", analysis_name, "\nStratified by Treatment Group (IPTW-weighted)"), break.time.by=12, xlim=c(0,48), lwd=1.0, censor=FALSE, palette="jama", ggtheme=theme_bw()+theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()))
  if (save_plot) { timestamp <- format(Sys.time(),"%Y%m%d_%H%M%S"); plot_filename <- paste0("ipwt_survival_", analysis_name, "_", gsub(" ","_",comparison_name), "_", gsub("[/ ]","_",outcome_label), "_", timestamp, ".png"); tryCatch({ png(filename=plot_filename, width=12, height=10, units="in", res=300); print(surv_plot); dev.off(); log_output("INFO","Survival plot saved", plot_filename) }, error=function(e) { log_output("ERROR","Failed to save survival plot", e$message) }) }
  return(surv_plot)
}

# >>> CELL 20: Main IPTW Analysis Loop (Outcome 4 - LATE ONSET) <<<
config <- get_standard_config()

# Late-onset analysis for Outcome 4
analysis_types <- list(
  list(name = "late_onset", late_onset_flag = TRUE, prefix = "outcome4_epilepsy_excl_hypo_late_onset_ipwt")
)

comparison_results <- lapply(cohorts, function(df) { list(cohort_df = df) })
names(comparison_results) <- c("SEMAGLUTIDE vs OTHER_GLPA","SEMAGLUTIDE vs SGLT2","SEMAGLUTIDE vs OtherGLD","OTHER_GLPA vs SGLT2","OTHER_GLPA vs OtherGLD","SGLT2 vs OtherGLD")

for (analysis in analysis_types) {
  log_output("PROGRESS", sprintf("STARTING ANALYSIS TYPE: %s", analysis$name), "", "main_script")
  ipwt_results <- list()

  for (comp_name in names(comparison_results)) {
    log_output("PROGRESS", sprintf("Processing comparison: %s", comp_name), "", "main_analysis")
    cohort_df <- comparison_results[[comp_name]]$cohort_df

    cohort_df <- cohort_df %>%
      mutate(
        renal_disease_severity = factor(case_when(renal_severe==1~"Severe", renal_mild_or_moderate==1~"Mild-Moderate", TRUE~"None"), levels=c("None","Mild-Moderate","Severe")),
        liver_disease_severity = factor(case_when(moderate_or_severe_liver_disease==1~"Mod-Severe", mild_liver_disease==1~"Mild", TRUE~"None"), levels=c("None","Mild","Mod-Severe")),
        malignancy_status = factor(case_when(metastatic_solid_tumor==1~"Metastatic", any_malignancy==1~"Non-Metastatic", TRUE~"None"), levels=c("None","Non-Metastatic","Metastatic"))
      ) %>% apply_variable_labels()

    for (var in config$categorical_vars) { if (var %in% names(cohort_df)) cohort_df[[var]] <- as.factor(cohort_df[[var]]) }

    exclude_vars <- switch(comp_name,
      "SEMAGLUTIDE vs OTHER_GLPA" = c("SEMAGLUTIDE","OTHER_GLPA"),
      "SEMAGLUTIDE vs SGLT2"      = c("SEMAGLUTIDE","SGLT2i"),
      "SEMAGLUTIDE vs OtherGLD"   = c("SEMAGLUTIDE","TZD","SU","DPP4i"),
      "OTHER_GLPA vs SGLT2"       = c("OTHER_GLPA","SGLT2i"),
      "OTHER_GLPA vs OtherGLD"    = c("OTHER_GLPA","TZD","SU","DPP4i"),
      "SGLT2 vs OtherGLD"         = c("SGLT2i","TZD","SU","DPP4i"), NULL)

    did_tags_for_comp <- FALSE

    for (outcome in config$outcomes) {
      outcome_var <- outcome$var; outcome_label <- outcome$label
      log_output("INFO", sprintf("Processing outcome: %s for %s (%s)", outcome_label, comp_name, analysis$name), "", "main_analysis")

      survival_df <- followup_and_event(df=cohort_df, outcome_var=outcome_var, data_cut_date=config$data_cut_date, late_onset=analysis$late_onset_flag, verbose=TRUE, comp_name=comp_name)

      pre_n <- nrow(survival_df)
      cond_renal <- !is.na(survival_df$renal_disease_severity) & survival_df$renal_disease_severity=="Severe"
      cond_malig <- !is.na(survival_df$malignancy_status) & survival_df$malignancy_status=="Metastatic"
      exclude_mask <- rep(FALSE, nrow(survival_df))
      if (isTRUE(config$exclude_renal_severe)) exclude_mask <- exclude_mask | cond_renal
      if (isTRUE(config$exclude_malignancy_metastatic)) exclude_mask <- exclude_mask | cond_malig
      survival_df <- survival_df %>% dplyr::filter(!exclude_mask) %>% droplevels()

      if (nrow(survival_df) < 10) { log_output("WARNING","Skipping analysis due to insufficient data", nrow(survival_df)); next }
      n_t1 <- sum(survival_df$treatment==1); n_t0 <- sum(survival_df$treatment==0)
      if (n_t1 < 2 || n_t0 < 2) { log_output("ERROR","Insufficient observations in one treatment group", sprintf("Treatment=1: %d, Treatment=0: %d", n_t1, n_t0), "main_analysis"); next }

      survival_df <- perform_mice_imputation(survival_df, m=20, seed=123)
      if ("baseline_bmi" %in% names(survival_df)) {
        survival_df <- survival_df %>% mutate(baseline_bmi_category = factor(case_when(is.na(baseline_bmi)~"Missing", baseline_bmi<18.5~"Underweight", baseline_bmi>=18.5&baseline_bmi<25~"Normal", baseline_bmi>=25&baseline_bmi<30~"Overweight", baseline_bmi>=30&baseline_bmi<35~"Obesity Class I", baseline_bmi>=35&baseline_bmi<40~"Obesity Class II", baseline_bmi>=40~"Obesity Class III"), levels=c("Underweight","Normal","Overweight","Obesity Class I","Obesity Class II","Obesity Class III","Missing")))
      }

      keep_ids <- unique(survival_df$person_id)
      cohort_df_post <- cohort_df %>% dplyr::semi_join(tibble(person_id=keep_ids), by="person_id")
      if (!did_tags_for_comp) { invisible(tag_and_make_tables_for_both_arms(cohort_df=cohort_df_post, comp_name=comp_name, analysis_prefix=analysis$prefix)); did_tags_for_comp <- TRUE }

      ipwt_result <- tryCatch({ run_ipwt_and_cox(df=survival_df, exclude_vars=exclude_vars, all_ps_vars=config$all_ps_vars, trim_threshold=0.01, verbose=TRUE, stabilize=FALSE) }, error=function(e) { log_output("ERROR", sprintf("run_ipwt_and_cox failed for %s", comp_name), e$message, "main_analysis"); return(NULL) })
      if (is.null(ipwt_result)) next

      crude_hr <- tryCatch({ calculate_crude_hr(cohort_df=ipwt_result$cohort, comparison_name=comp_name, outcome_label=outcome_label) }, error=function(e) { log_output("ERROR", sprintf("calculate_crude_hr failed for %s", comp_name), e$message, "main_analysis"); return(NULL) })
      incident_rates <- tryCatch({ calculate_incident_rates(cohort_df=ipwt_result$cohort, weight_var="ipw_std") }, error=function(e) { log_output("ERROR", sprintf("calculate_incident_rates failed for %s", comp_name), e$message, "main_analysis"); return(NULL) })

      ipwt_result$crude_hr_results <- crude_hr; ipwt_result$incident_rate_results <- incident_rates
      analysis_key <- paste(comp_name, "-", outcome_label)
      ipwt_results[[analysis_key]] <- ipwt_result

      # Baseline table for IPTW cohort
      try({
        iptw_cohort <- ipwt_result$cohort
        desc_vars <- unique(c(my_vars, categorical_vars)); desc_vars <- intersect(desc_vars, names(cohort_df_post)); keep_desc <- setdiff(desc_vars, names(iptw_cohort))
        if (length(keep_desc)>0) { iptw_cohort <- iptw_cohort %>% dplyr::left_join(cohort_df_post %>% dplyr::select(person_id, dplyr::all_of(keep_desc)), by="person_id") }
        bin_candidates <- names(iptw_cohort)[vapply(iptw_cohort, is_binary_like, logical(1))]
        if (length(bin_candidates)>0) { iptw_cohort[bin_candidates] <- lapply(iptw_cohort[bin_candidates], to_yesno) }
        iptw_vars <- intersect(my_vars, names(iptw_cohort)); iptw_cat_vars <- intersect(categorical_vars, iptw_vars)
        for (v in iptw_cat_vars) { if (!is.factor(iptw_cohort[[v]])) iptw_cohort[[v]] <- factor(iptw_cohort[[v]]) }
        bin_factor_vars <- intersect(bin_candidates, iptw_vars)
        if (length(iptw_vars)>0) {
          iptw_table_obj <- CreateTableOne(vars=iptw_vars, strata="treatment", data=iptw_cohort, factorVars=unique(c(iptw_cat_vars, bin_factor_vars)))
          iptw_table_df <- as.data.frame(print(iptw_table_obj, showAllLevels=TRUE, printToggle=FALSE, noSpaces=TRUE))
          ts <- format(Sys.time(),"%Y%m%d_%H%M%S"); safe_comp <- gsub("[^A-Za-z0-9]+","_",comp_name); safe_outcome <- gsub("[^A-Za-z0-9]+","_",outcome_label)
          fname <- sprintf("outcome4_epilepsy_excl_hypo_baseline_table_IPTW_%s_%s_%s_%s.csv", analysis$prefix, safe_comp, safe_outcome, ts)
          write.csv(iptw_table_df, file=fname, row.names=TRUE); log_output("PROGRESS","Saved IPTW analysis cohort baseline table", fname, "baseline_table_IPTW")
        }
      }, silent=TRUE)

      balance_plot <- tryCatch({ create_balance_plot(balance_before=ipwt_result$balance_before, balance_after=ipwt_result$balance_after, title=sprintf("Balance for %s - %s (%s)", comp_name, outcome_label, analysis$name)) }, error=function(e) { NULL })
      if (!is.null(balance_plot)) { timestamp <- format(Sys.time(),"%Y%m%d_%H%M%S"); ggsave(sprintf("%s_balance_%s_%s_%s.png", analysis$prefix, gsub(" ","_",comp_name), gsub("[/ ]","_",outcome_label), timestamp), balance_plot) }

      surv_plot <- tryCatch({ create_survival_plot(ipwt_result=ipwt_result, comparison_name=comp_name, outcome_label=outcome_label, analysis_name=analysis$name, save_plot=TRUE) }, error=function(e) { NULL })
    } # end outcomes loop
  }   # end comparisons loop

  saved_files <- tryCatch({ save_analysis_results(ipwt_results, prefix=analysis$prefix) }, error=function(e) { log_output("ERROR","save_analysis_results failed", e$message, "main_analysis"); NULL })

  rds_filename <- sprintf("%s_full_results.rds", analysis$prefix)
  saveRDS(ipwt_results, file=rds_filename)
  log_output("PROGRESS","Full results object saved to", rds_filename, "main_script")

  timestamp <- format(Sys.time(),"%Y%m%d_%H%M%S")
  save_output_log(sprintf("%s_log_%s.csv", analysis$prefix, timestamp))
  log_output("PROGRESS", sprintf("ANALYSIS TYPE '%s' COMPLETE", analysis$name), "", "main_script")
}

message("IPTW analysis for Outcome 4 - LATE ONSET (Epilepsy/Seizure excl Hypoglycemic) completed.")

# >>> CELL 21: Bootstrap Helpers <<<
perform_tmle_bootstrap <- function(df, outcome_var="event", exclude_vars=NULL, all_ps_vars=NULL, sl_lib=c("SL.glm"), n_bootstrap=1000, seed=12345, min_events_per_arm=1) {
  if (!requireNamespace("future", quietly=TRUE)) install.packages("future")
  if (!requireNamespace("future.apply", quietly=TRUE)) install.packages("future.apply")
  Sys.setenv(OMP_NUM_THREADS="1", MKL_NUM_THREADS="1", OPENBLAS_NUM_THREADS="1")
  if (requireNamespace("data.table", quietly=TRUE)) data.table::setDTthreads(1)
  n <- nrow(df); set.seed(seed); idx_list <- lapply(seq_len(n_bootstrap), function(b) sample.int(n, size=n, replace=TRUE))
  future::plan(future::multisession, workers=max(1, parallel::detectCores()-1))
  res <- future.apply::future_lapply(seq_along(idx_list), function(b) {
    boot_df <- df[idx_list[[b]], , drop=FALSE]; t1 <- boot_df$treatment==1; t0 <- boot_df$treatment==0
    if (sum(t1)<2 || sum(t0)<2 || sum(boot_df[[outcome_var]][t1], na.rm=TRUE)<min_events_per_arm || sum(boot_df[[outcome_var]][t0], na.rm=TRUE)<min_events_per_arm) return(NA_real_)
    fit_b <- try(run_tmle_binary(df=boot_df, exclude_vars=exclude_vars, all_ps_vars=all_ps_vars, sl_lib=sl_lib, verbose=FALSE), silent=TRUE)
    if (inherits(fit_b,"try-error") || is.null(fit_b) || is.null(fit_b$fit$estimates$ATE$psi) || !is.finite(fit_b$fit$estimates$ATE$psi)) return(NA_real_)
    as.numeric(fit_b$fit$estimates$ATE$psi)
  }, future.seed=TRUE)
  rds <- unlist(res); samples <- rds[is.finite(rds)]; fails <- sum(!is.finite(rds))
  if (length(samples)==0L) return(list(samples=numeric(0), mean=NA_real_, sd=NA_real_, ci_lower_2.5=NA_real_, ci_upper_97.5=NA_real_, n_success=0L, n_fail=n_bootstrap))
  list(samples=samples, mean=mean(samples), sd=stats::sd(samples), ci_lower_2.5=as.numeric(stats::quantile(samples, probs=0.025, names=FALSE)), ci_upper_97.5=as.numeric(stats::quantile(samples, probs=0.975, names=FALSE)), n_success=length(samples), n_fail=fails)
}

save_bootstrap_results <- function(bootstrap_stats, comparison, outcome, analysis_prefix, tmle_point_estimate=NA_real_, write_hist_png=TRUE) {
  ts <- format(Sys.time(),"%Y%m%d_%H%M%S"); safe_comp <- gsub("[^A-Za-z0-9]+","_",comparison); safe_outc <- gsub("[^A-Za-z0-9]+","_",outcome)
  base_prefix <- sprintf("%s_bootstrap_%s_%s_%s", analysis_prefix, safe_comp, safe_outc, ts)
  rds_path <- sprintf("%s.rds", base_prefix); saveRDS(bootstrap_stats, file=rds_path)
  summ_df <- data.frame(comparison=comparison, outcome=outcome, n_success=bootstrap_stats$n_success, n_fail=bootstrap_stats$n_fail, mean_rd=bootstrap_stats$mean, sd_rd=bootstrap_stats$sd, ci_lower_2.5=bootstrap_stats$ci_lower_2.5, ci_upper_97.5=bootstrap_stats$ci_upper_97.5, tmle_point_estimate=tmle_point_estimate, stringsAsFactors=FALSE)
  csv_path <- sprintf("%s_summary.csv", base_prefix); utils::write.csv(summ_df, file=csv_path, row.names=FALSE)
  png_path <- NA_character_
  if (write_hist_png && length(bootstrap_stats$samples)>0) {
    png_path <- sprintf("%s_hist.png", base_prefix); grDevices::png(png_path, width=1000, height=700, res=120)
    hist(bootstrap_stats$samples, main=sprintf("Bootstrap RD: %s - %s", comparison, outcome), xlab="Risk Difference (ATE)", breaks=30)
    abline(v=tmle_point_estimate, lty=2, lwd=2); abline(v=0, lty=3); legend("topright", legend=c("TMLE point estimate","Null (0)"), lty=c(2,3), lwd=c(2,1)); grDevices::dev.off()
  }
  list(rds=rds_path, csv=csv_path, png=png_path)
}

# >>> CELL 22: Combined Bootstrap Visualization <<<
visualize_combined_bootstrap <- function(comparisons, outcome="Epilepsy/Seizure excl Hypoglycemic", analysis_prefix="outcome4_epilepsy_excl_hypo_late_onset_tmle", colors=NULL, use_jama_style=TRUE, plot_type="multi") {
  jama_colors <- c("#00274C","#C4622D","#4B7AA7","#969696","#2F7F7E")
  if (is.null(colors)) colors <- if (use_jama_style) jama_colors else c("darkblue","darkred","darkgreen","purple","orange")
  if (use_jama_style) { old_par <- par(no.readonly=TRUE); par(family="sans", cex.main=1.1, cex.lab=1.0, cex.axis=0.9, font.main=1, las=1, mgp=c(2.5,0.7,0), tcl=-0.3) }
  all_bootstrap_data <- list()
  for (i in seq_along(comparisons)) {
    comparison <- comparisons[i]; safe_comp <- gsub("[^A-Za-z0-9]+","_",comparison); safe_outc <- gsub("[^A-Za-z0-9]+","_",outcome)
    pattern <- sprintf("^%s_bootstrap_%s_%s_\\d{8}_\\d{6}\\.rds$", analysis_prefix, safe_comp, safe_outc)
    files <- list.files(pattern=pattern)
    if (length(files)==0) { warning(sprintf("No RDS bootstrap files found for %s (%s)", comparison, outcome)); next }
    mt <- file.info(files)$mtime; rds_filename <- files[order(mt, decreasing=TRUE)][1]
    all_bootstrap_data[[comparison]] <- readRDS(rds_filename)
  }
  if (length(all_bootstrap_data)==0) stop("No bootstrap data could be loaded")
  all_samples <- unlist(lapply(all_bootstrap_data, function(x) x$samples)); x_range <- range(c(all_samples, 0))*1.1
  densities <- lapply(all_bootstrap_data, function(x) density(x$samples)); y_max <- max(sapply(densities, function(d) max(d$y)))*1.1
  if (plot_type=="single" || plot_type=="both") {
    plot(NULL, xlim=x_range, ylim=c(0,y_max), main=paste("Bootstrap Distributions:", outcome), xlab="Risk Difference", ylab="Density", type="n", axes=FALSE, frame.plot=FALSE)
    axis(1, col="black", lwd=0.5); axis(2, col="black", lwd=0.5); box(lwd=0.5); abline(v=0, col="black", lwd=1.5, lty=1)
    for (i in seq_along(all_bootstrap_data)) { dens <- densities[[i]]; col_i <- colors[(i-1)%%length(colors)+1]; lines(dens, col=col_i, lwd=2); polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0, length(dens$y))), col=adjustcolor(col_i, alpha.f=0.15), border=NA); abline(v=all_bootstrap_data[[i]]$mean, col=col_i, lwd=1, lty=2) }
  }
  if (plot_type=="multi" || plot_type=="both") {
    par(mfrow=c(2,2))
    plot(NULL, xlim=x_range, ylim=c(0,y_max), main="A. Density Distributions", xlab="Risk Difference", ylab="Density", type="n", axes=FALSE, frame.plot=FALSE)
    axis(1, lwd=0.5); axis(2, lwd=0.5); box(lwd=0.5); abline(v=0, col="black", lwd=1.5, lty=1)
    for (i in seq_along(all_bootstrap_data)) { dens <- densities[[i]]; col_i <- colors[(i-1)%%length(colors)+1]; lines(dens, col=col_i, lwd=2); polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0,length(dens$y))), col=adjustcolor(col_i, alpha.f=0.15), border=NA) }
    boxplot_data <- lapply(all_bootstrap_data, function(x) x$samples)
    boxplot(boxplot_data, names=names(all_bootstrap_data), col=adjustcolor(colors[seq_along(all_bootstrap_data)], alpha.f=0.3), border=colors[seq_along(all_bootstrap_data)], main="B. Distribution Comparison", ylab="Risk Difference", las=2, axes=FALSE, frame.plot=FALSE)
    axis(1, at=seq_along(all_bootstrap_data), labels=names(all_bootstrap_data), las=2, lwd=0.5); axis(2, lwd=0.5); box(lwd=0.5); abline(h=0, col="black", lwd=1.5, lty=1)
    plot(NULL, xlim=c(0.5, length(all_bootstrap_data)+0.5), ylim=x_range, main="C. 95% Confidence Intervals", xlab="", ylab="Risk Difference", xaxt="n", axes=FALSE, frame.plot=FALSE)
    axis(1, at=seq_along(all_bootstrap_data), labels=names(all_bootstrap_data), las=2, lwd=0.5); axis(2, lwd=0.5); box(lwd=0.5); abline(h=0, col="black", lwd=1.5, lty=1)
    for (i in seq_along(all_bootstrap_data)) { d <- all_bootstrap_data[[i]]; col_i <- colors[(i-1)%%length(colors)+1]; segments(i, d$ci_lower_2.5, i, d$ci_upper_97.5, col=col_i, lwd=3); points(i, d$mean, pch=19, col=col_i, cex=1.5) }
    plot.new(); plot.window(xlim=c(0,1), ylim=c(0,1)); title("D. Statistical Summary"); y_pos <- 0.9
    for (i in seq_along(all_bootstrap_data)) { nm <- names(all_bootstrap_data)[i]; d <- all_bootstrap_data[[i]]; col_i <- colors[(i-1)%%length(colors)+1]; text(0.05, y_pos, nm, col=col_i, font=1, adj=0, cex=0.85); y_pos <- y_pos-0.08; text(0.1, y_pos, paste("Mean:", sprintf("%.5f", d$mean)), adj=0, cex=0.75); y_pos <- y_pos-0.06; text(0.1, y_pos, paste("95% CI: (", sprintf("%.5f", d$ci_lower_2.5), ",", sprintf("%.5f", d$ci_upper_97.5), ")"), adj=0, cex=0.75); y_pos <- y_pos-0.06; text(0.1, y_pos, paste("P(RD<0):", sprintf("%.3f", mean(d$samples<0))), adj=0, cex=0.75); y_pos <- y_pos-0.1 }
    par(mfrow=c(1,1))
  }
  if (use_jama_style && exists("old_par")) par(old_par)
  invisible(all_bootstrap_data)
}

# >>> CELL 23: TMLE Core Functions <<<
if (!requireNamespace("tmle", quietly=TRUE)) install.packages("tmle")
if (!requireNamespace("SuperLearner", quietly=TRUE)) install.packages("SuperLearner")
suppressPackageStartupMessages(library(tmle))

COVARIATE_CAP <- NULL

run_tmle_binary <- function(df, exclude_vars=NULL, all_ps_vars=NULL, sl_lib=c("SL.glm"), verbose=FALSE, use_imputation=TRUE, use_superlearner=FALSE) {
  t_start <- Sys.time()
  all_covariates <- if (is.null(all_ps_vars)) get_standard_config()$all_ps_vars else all_ps_vars
  rhs_vars <- setdiff(all_covariates, exclude_vars)
  keep_vars <- rhs_vars[sapply(df[rhs_vars], function(x) length(unique(na.omit(x)))>1)]
  model_vars <- c("treatment", keep_vars)
  if (verbose) log_output("PROGRESS","TMLE: begin run_tmle_binary", sprintf("n=%d; covariates_kept=%d", nrow(df), length(keep_vars)), "tmle_fit")
  if (isTRUE(use_imputation)) { df_impute <- perform_mice_imputation(df, m=20, seed=123) } else { df_impute <- df }
  if (isTRUE(use_imputation) && "baseline_bmi" %in% names(df_impute)) {
    df_impute <- df_impute %>% mutate(baseline_bmi_category = factor(case_when(is.na(baseline_bmi)~"Missing", baseline_bmi<18.5~"Underweight", baseline_bmi>=18.5&baseline_bmi<25~"Normal", baseline_bmi>=25&baseline_bmi<30~"Overweight", baseline_bmi>=30&baseline_bmi<35~"Obesity Class I", baseline_bmi>=35&baseline_bmi<40~"Obesity Class II", baseline_bmi>=40~"Obesity Class III"), levels=c("Underweight","Normal","Overweight","Obesity Class I","Obesity Class II","Obesity Class III","Missing")))
  }
  complete_rows <- complete.cases(df_impute[, model_vars, drop=FALSE]); df_complete <- df_impute[complete_rows, , drop=FALSE]
  if (verbose) log_output("INFO", sprintf("TMLE: complete cases retained = %d / %d", nrow(df_complete), nrow(df)))
  W_vars <- df_complete[, keep_vars, drop=FALSE]
  if (!is.null(COVARIATE_CAP) && ncol(W_vars)>COVARIATE_CAP) W_vars <- W_vars[, seq_len(COVARIATE_CAP), drop=FALSE]
  if (ncol(W_vars)==0) W_vars <- data.frame(intercept=rep(1, nrow(df_complete)))
  use_sl <- isTRUE(use_superlearner) && requireNamespace("SuperLearner", quietly=TRUE)
  fit <- NULL
  if (use_sl) {
    sl_present <- unique(sl_lib)
    fit <- try(tmle::tmle(Y=df_complete$event, A=df_complete$treatment, W=W_vars, family="binomial", Q.SL.library=sl_present, g.SL.library=sl_present), silent=TRUE)
  }
  if (is.null(fit) || inherits(fit,"try-error")) {
    data_for_models <- cbind(df_complete[, c("event","treatment"), drop=FALSE], W_vars)
    g_form <- as.formula(paste("treatment ~", paste(colnames(W_vars), collapse=" + ")))
    g_fit <- try(stats::glm(g_form, data=data_for_models, family=stats::binomial()), silent=TRUE)
    g1W <- suppressWarnings(pmin(pmax(as.numeric(stats::predict(g_fit, type="response")), 1e-5), 1-1e-5))
    q_form <- as.formula(paste("event ~ treatment +", paste(colnames(W_vars), collapse=" + ")))
    q_fit <- try(stats::glm(q_form, data=data_for_models, family=stats::binomial()), silent=TRUE)
    newdata1 <- data_for_models; newdata1$treatment <- 1; newdata0 <- data_for_models; newdata0$treatment <- 0
    Q1W <- suppressWarnings(pmin(pmax(as.numeric(stats::predict(q_fit, newdata=newdata1, type="response")), 1e-5), 1-1e-5))
    Q0W <- suppressWarnings(pmin(pmax(as.numeric(stats::predict(q_fit, newdata=newdata0, type="response")), 1e-5), 1-1e-5))
    Q_init <- cbind(Q1W=Q1W, Q0W=Q0W)
    fit <- tmle::tmle(Y=df_complete$event, A=df_complete$treatment, W=W_vars, family="binomial", Q=Q_init, g1W=g1W)
  }
  list(fit=fit, cohort=df_complete, keep_vars=keep_vars, n_dropped=nrow(df)-nrow(df_complete))
}

tmle_sl_lib <- c("SL.glm")
if (requireNamespace("glmnet", quietly=TRUE)) tmle_sl_lib <- c("SL.glmnet", tmle_sl_lib)
if (requireNamespace("xgboost", quietly=TRUE)) tmle_sl_lib <- c("SL.xgboost", tmle_sl_lib)

calculate_tmle_time_specific_rd <- function(df, time_points, time_labels) {
  if (!all(c("event_time","event","treatment") %in% names(df))) return(NULL)
  safe_len <- min(length(time_points), length(time_labels)); if (safe_len==0) return(NULL)
  time_points <- time_points[seq_len(safe_len)]; time_labels <- time_labels[seq_len(safe_len)]
  calc_surv_at <- function(data_subset, t_point) {
    if (nrow(data_subset)<1) return(list(surv=NA_real_, se=NA_real_))
    fit <- survival::survfit(survival::Surv(event_time, event) ~ 1, data=data_subset)
    s <- try(suppressWarnings(summary(fit, times=t_point)), silent=TRUE)
    if (inherits(s,"try-error") || length(s$surv)==0) return(list(surv=NA_real_, se=NA_real_))
    se <- if (!is.null(s$std.err) && length(s$std.err)>0) s$std.err[1] else NA_real_
    list(surv=s$surv[1], se=se)
  }
  out <- vector("list", length(time_points))
  for (i in seq_along(time_points)) {
    t_point <- time_points[i]; lbl <- time_labels[i]
    grp0 <- df[df$treatment==0, , drop=FALSE]; grp1 <- df[df$treatment==1, , drop=FALSE]
    s0 <- calc_surv_at(grp0, t_point); s1 <- calc_surv_at(grp1, t_point)
    risk0 <- if (!is.na(s0$surv)) 1-s0$surv else NA_real_; risk1 <- if (!is.na(s1$surv)) 1-s1$surv else NA_real_
    rd <- if (!is.na(risk0) && !is.na(risk1)) risk1-risk0 else NA_real_
    var_rd <- sum(c(s0$se, s1$se)^2, na.rm=TRUE); se_rd <- if (is.finite(var_rd)) sqrt(var_rd) else NA_real_
    lo <- if (!is.na(rd) && !is.na(se_rd)) rd-1.96*se_rd else NA_real_; hi <- if (!is.na(rd) && !is.na(se_rd)) rd+1.96*se_rd else NA_real_
    out[[i]] <- data.frame(time_label=lbl, time_days=t_point, risk_treated=risk1, risk_control=risk0, rd=rd, rd_lo_95=lo, rd_hi_95=hi, stringsAsFactors=FALSE)
  }
  do.call(rbind, out)
}

# >>> CELL 24: Main TMLE Analysis Loop (Outcome 4 - LATE ONSET) <<<
for (analysis in analysis_types) {
  log_output("PROGRESS", sprintf("STARTING TMLE ANALYSIS TYPE: %s", analysis$name), "", "tmle_main")
  analysis_t0 <- Sys.time()
  if (!identical(analysis$name, "late_onset")) next
  tmle_results <- list(); tmle_summary_rows <- list()

  iptw_rds_path <- if (exists("IPTW_RDS_PATH_OVERRIDE") && !is.null(IPTW_RDS_PATH_OVERRIDE)) IPTW_RDS_PATH_OVERRIDE else sprintf("%s_full_results.rds", analysis$prefix)
  iptw_obj <- try(readRDS(iptw_rds_path), silent=TRUE)

  allowed_comparisons <- c("SEMAGLUTIDE vs OtherGLD", "SEMAGLUTIDE vs SGLT2")
  comps <- intersect(names(comparison_results), allowed_comparisons)
  if (exists("LIMIT_COMPARISONS") && !is.null(LIMIT_COMPARISONS)) comps <- comps[seq_len(min(length(comps), LIMIT_COMPARISONS))]

  for (comp_name in comps) {
    if (!(comp_name %in% allowed_comparisons)) next
    comp_t0 <- Sys.time()
    log_output("PROGRESS", sprintf("TMLE: Processing comparison: %s", comp_name), "", "tmle_main")
    cohort_df <- comparison_results[[comp_name]]$cohort_df
    if (is.null(cohort_df) || nrow(cohort_df)==0) { log_output("WARNING","Empty cohort; skipping", comp_name, "tmle_main"); next }

    cohort_df <- cohort_df %>%
      mutate(
        renal_disease_severity = factor(case_when(renal_severe==1~"Severe", renal_mild_or_moderate==1~"Mild-Moderate", TRUE~"None"), levels=c("None","Mild-Moderate","Severe")),
        liver_disease_severity = factor(case_when(moderate_or_severe_liver_disease==1~"Mod-Severe", mild_liver_disease==1~"Mild", TRUE~"None"), levels=c("None","Mild","Mod-Severe")),
        malignancy_status = factor(case_when(metastatic_solid_tumor==1~"Metastatic", any_malignancy==1~"Non-Metastatic", TRUE~"None"), levels=c("None","Non-Metastatic","Metastatic"))
      ) %>% apply_variable_labels()
    for (var in config$categorical_vars) { if (var %in% names(cohort_df)) cohort_df[[var]] <- as.factor(cohort_df[[var]]) }

    exclude_vars <- switch(comp_name,
      "SEMAGLUTIDE vs OTHER_GLPA"=c("SEMAGLUTIDE","OTHER_GLPA"), "SEMAGLUTIDE vs SGLT2"=c("SEMAGLUTIDE","SGLT2i"),
      "SEMAGLUTIDE vs OtherGLD"=c("SEMAGLUTIDE","TZD","SU","DPP4i"), "OTHER_GLPA vs SGLT2"=c("OTHER_GLPA","SGLT2i"),
      "OTHER_GLPA vs OtherGLD"=c("OTHER_GLPA","TZD","SU","DPP4i"), "SGLT2 vs OtherGLD"=c("SGLT2i","TZD","SU","DPP4i"), NULL)

    outs <- config$outcomes
    if (exists("LIMIT_OUTCOMES") && !is.null(LIMIT_OUTCOMES)) outs <- outs[seq_len(min(length(outs), LIMIT_OUTCOMES))]

    for (outcome in outs) {
      out_t0 <- Sys.time(); outcome_var <- outcome$var; outcome_label <- outcome$label
      analysis_key <- paste(comp_name, "-", outcome_label)
      log_output("INFO", sprintf("TMLE: Processing outcome: %s for %s (%s)", outcome_label, comp_name, analysis$name), "", "tmle_main")

      if (!inherits(iptw_obj,"try-error") && !is.null(iptw_obj[[analysis_key]]) && !is.null(iptw_obj[[analysis_key]]$cohort)) {
        survival_df <- iptw_obj[[analysis_key]]$cohort
        log_output("PROGRESS","TMLE: reusing exact IPTW cohort", sprintf("n=%d", nrow(survival_df)), "tmle_main")
      } else {
        survival_df <- followup_and_event(df=cohort_df, outcome_var=outcome_var, data_cut_date=config$data_cut_date, late_onset=analysis$late_onset_flag, verbose=TRUE, comp_name=comp_name)
        pre_n <- nrow(survival_df)
        cond_renal <- !is.na(survival_df$renal_disease_severity) & survival_df$renal_disease_severity=="Severe"
        cond_malig <- !is.na(survival_df$malignancy_status) & survival_df$malignancy_status=="Metastatic"
        exclude_mask <- rep(FALSE, nrow(survival_df))
        if (isTRUE(config$exclude_renal_severe)) exclude_mask <- exclude_mask | cond_renal
        if (isTRUE(config$exclude_malignancy_metastatic)) exclude_mask <- exclude_mask | cond_malig
        if (any(exclude_mask)) survival_df <- survival_df %>% dplyr::filter(!exclude_mask) %>% droplevels()
      }

      if (nrow(survival_df)<10) { log_output("WARNING","Skipping TMLE due to insufficient data", nrow(survival_df), "tmle_main"); next }
      n_t1 <- sum(survival_df$treatment==1); n_t0 <- sum(survival_df$treatment==0)
      if (n_t1<2 || n_t0<2) { log_output("ERROR","Insufficient observations", sprintf("Treatment=1: %d, Treatment=0: %d", n_t1, n_t0), "tmle_main"); next }

      tmle_fit <- tryCatch({
        run_tmle_binary(df=survival_df, exclude_vars=exclude_vars, all_ps_vars=config$all_ps_vars, sl_lib=tmle_sl_lib, verbose=TRUE, use_imputation=FALSE, use_superlearner=TRUE)
      }, error=function(e) { log_output("ERROR", sprintf("run_tmle_binary failed for %s", comp_name), e$message, "tmle_main"); return(NULL) })
      if (is.null(tmle_fit)) next

      tmle_results[[analysis_key]] <- tmle_fit
      log_output("PROGRESS","TMLE: fit stored", sprintf("analysis=%s; comp=%s; outcome=%s", analysis$name, comp_name, outcome_label), "tmle_main")

      # Plot TMLE estimated event probabilities
      ey0_psi <- try(tmle_fit$fit$estimates$EY0$psi, silent=TRUE); ey1_psi <- try(tmle_fit$fit$estimates$EY1$psi, silent=TRUE)
      ey0_lo <- try(tmle_fit$fit$estimates$EY0$CI[1], silent=TRUE); ey0_hi <- try(tmle_fit$fit$estimates$EY0$CI[2], silent=TRUE)
      ey1_lo <- try(tmle_fit$fit$estimates$EY1$CI[1], silent=TRUE); ey1_hi <- try(tmle_fit$fit$estimates$EY1$CI[2], silent=TRUE)
      ate_obj <- try(tmle_fit$fit$estimates$ATE, silent=TRUE)
      ate_psi <- if (!inherits(ate_obj,"try-error") && !is.null(ate_obj$psi)) ate_obj$psi else NA_real_
      ate_p <- if (!inherits(ate_obj,"try-error") && !is.null(ate_obj$pvalue)) ate_obj$pvalue else NA_real_
      to_num <- function(x) if (!inherits(x,"try-error")) as.numeric(x) else NA_real_
      prob_df <- data.frame(Group=c("Control","Semaglutide"), Probability=c(to_num(ey0_psi), to_num(ey1_psi)), Lower=c(to_num(ey0_lo), to_num(ey1_lo)), Upper=c(to_num(ey0_hi), to_num(ey1_hi)))
      try({
        tmle_plot <- ggplot2::ggplot(prob_df, ggplot2::aes(x=Group, y=Probability, fill=Group)) + ggplot2::geom_bar(stat="identity", alpha=0.7) + ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower, ymax=Upper), width=0.2) + ggplot2::scale_fill_manual(values=c("Control"="#1f77b4","Semaglutide"="#ff7f0e")) + ggplot2::labs(title=paste("TMLE Estimated Event Probabilities for", outcome_label), subtitle=paste(comp_name,"- Risk Difference =", sprintf("%.4f", ate_psi), "(p =", ifelse(is.na(ate_p),"NA", sprintf("%.3g", ate_p)), ")"), y="Probability of Event", x="") + ggplot2::theme_minimal() + ggplot2::theme(legend.position="none")
        print(tmle_plot); safe_comp <- gsub("[^A-Za-z0-9]+","_",comp_name); safe_outcome <- gsub("[^A-Za-z0-9]+","_",outcome_label)
        ggplot2::ggsave(filename=paste0("tmle_probs_", safe_comp, "_", safe_outcome, ".png"), plot=tmle_plot, width=8, height=6)
      }, silent=TRUE)

      time_points <- c(182.5, 365, 547.5, 730); time_labels <- c("6_months","1_year","18_months","2_years")
      time_specific_results <- tryCatch({ calculate_tmle_time_specific_rd(tmle_fit$cohort, time_points, time_labels) }, error=function(e) { NULL })

      ts <- format(Sys.time(),"%Y%m%d_%H%M%S"); safe_comp <- gsub("[^A-Za-z0-9]+","_",comp_name); safe_outcome <- gsub("[^A-Za-z0-9]+","_",outcome_label)
      tmle_prefix <- sub("_ipwt$","_tmle", analysis$prefix)
      ate <- tmle_fit$fit$estimates$ATE; or <- tmle_fit$fit$estimates$OR
      EY1 <- try(tmle_fit$fit$estimates$EY1$psi, silent=TRUE); EY0 <- try(tmle_fit$fit$estimates$EY0$psi, silent=TRUE)
      rr_obj <- try(tmle_fit$fit$estimates$RR, silent=TRUE)
      rr_est <- if (!inherits(rr_obj,"try-error") && !is.null(rr_obj$psi)) rr_obj$psi else if (!inherits(EY1,"try-error") && !inherits(EY0,"try-error") && !is.null(EY0) && EY0>0) EY1/EY0 else NA_real_
      rr_lo <- if (!inherits(rr_obj,"try-error") && !is.null(rr_obj$CI)) rr_obj$CI[1] else NA_real_
      rr_hi <- if (!inherits(rr_obj,"try-error") && !is.null(rr_obj$CI)) rr_obj$CI[2] else NA_real_
      rr_p <- if (!inherits(rr_obj,"try-error") && !is.null(rr_obj$pvalue)) rr_obj$pvalue else NA_real_
      nnt_val <- if (!is.null(ate$psi) && is.finite(ate$psi) && ate$psi!=0) 1/abs(ate$psi) else Inf
      n_treat <- sum(tmle_fit$cohort$treatment==1, na.rm=TRUE); n_ctrl <- sum(tmle_fit$cohort$treatment==0, na.rm=TRUE)
      ev_treat <- sum(tmle_fit$cohort$event==1 & tmle_fit$cohort$treatment==1, na.rm=TRUE); ev_ctrl <- sum(tmle_fit$cohort$event==1 & tmle_fit$cohort$treatment==0, na.rm=TRUE)

      out_df <- data.frame(method="TMLE", comparison=comp_name, outcome=outcome_label,
        effect_measure=c("Risk Difference","Odds Ratio","Risk Ratio","Number Needed to Treat"),
        estimate=c(ate$psi, or$psi, rr_est, nnt_val), lower_ci=c(ate$CI[1], or$CI[1], rr_lo, NA_real_), upper_ci=c(ate$CI[2], or$CI[2], rr_hi, NA_real_), p_value=c(ate$pvalue, or$pvalue, rr_p, NA_real_),
        n_total=rep(nrow(tmle_fit$cohort),4), n_events=rep(sum(tmle_fit$cohort$event, na.rm=TRUE),4),
        n_at_risk_treatment=rep(n_treat,4), n_at_risk_control=rep(n_ctrl,4), events_treatment=rep(ev_treat,4), events_control=rep(ev_ctrl,4), stringsAsFactors=FALSE)
      tmle_summary_rows[[length(tmle_summary_rows)+1]] <- out_df
      fname <- sprintf("%s_summary_%s_%s_%s.csv", tmle_prefix, safe_comp, safe_outcome, ts)
      if (isTRUE(WRITE_OUTPUTS)) try(utils::write.csv(out_df, file=fname, row.names=FALSE), silent=TRUE)
      if (!is.null(time_specific_results) && isTRUE(WRITE_OUTPUTS)) { fname_time <- sprintf("%s_time_specific_RD_%s_%s_%s.csv", tmle_prefix, safe_comp, safe_outcome, ts); try(utils::write.csv(time_specific_results, file=fname_time, row.names=FALSE), silent=TRUE) }

      # Baseline table for TMLE cohort
      try({
        tmle_cohort <- tmle_fit$cohort; tmle_vars <- intersect(my_vars, names(tmle_cohort)); tmle_cat_vars <- intersect(categorical_vars, tmle_vars)
        for (v in tmle_cat_vars) { if (!is.factor(tmle_cohort[[v]])) tmle_cohort[[v]] <- as.factor(tmle_cohort[[v]]) }
        if (length(tmle_vars)>0) {
          tmle_table_obj <- CreateTableOne(vars=tmle_vars, strata="treatment", data=tmle_cohort, factorVars=tmle_cat_vars)
          tmle_table_df <- as.data.frame(print(tmle_table_obj, showAllLevels=TRUE, printToggle=FALSE, noSpaces=TRUE))
          ts2 <- format(Sys.time(),"%Y%m%d_%H%M%S")
          fname2 <- sprintf("outcome4_epilepsy_excl_hypo_baseline_table_TMLE_%s_%s_%s_%s.csv", tmle_prefix, safe_comp, safe_outcome, ts2)
          if (isTRUE(WRITE_OUTPUTS)) { write.csv(tmle_table_df, file=fname2, row.names=TRUE); log_output("PROGRESS","Saved TMLE baseline table", fname2, "baseline_table_TMLE") }
        }
      }, silent=TRUE)
    }
  }

  tmle_prefix <- sub("_ipwt$","_tmle", analysis$prefix)
  rds_filename <- sprintf("%s_full_results.rds", tmle_prefix)
  if (isTRUE(WRITE_OUTPUTS)) { saveRDS(tmle_results, file=rds_filename); log_output("PROGRESS","TMLE results object saved to", rds_filename, "tmle_main") }
  timestamp <- format(Sys.time(),"%Y%m%d_%H%M%S")
  if (isTRUE(WRITE_OUTPUTS)) save_output_log(sprintf("%s_log_%s.csv", tmle_prefix, timestamp))
  if (length(tmle_summary_rows)>0 && isTRUE(WRITE_OUTPUTS)) {
    final_tmle_df <- tryCatch({ do.call(rbind, tmle_summary_rows) }, error=function(e) NULL)
    if (!is.null(final_tmle_df)) { comp_results_file <- sprintf("%s_comprehensive_results_%s.csv", analysis$prefix, timestamp); utils::write.csv(final_tmle_df, comp_results_file, row.names=FALSE); log_output("PROGRESS","Aggregated TMLE comprehensive results saved", comp_results_file, "tmle_main") }
  }
  log_output("PROGRESS", sprintf("TMLE ANALYSIS TYPE '%s' COMPLETE", analysis$name), "", "tmle_main")
}

# >>> CELL 25: Final Bootstrap Visualization <<<
visualize_combined_bootstrap(
  comparisons = c("SEMAGLUTIDE vs OtherGLD", "SEMAGLUTIDE vs SGLT2"),
  outcome = "Epilepsy/Seizure excl Hypoglycemic",
  analysis_prefix = "outcome4_epilepsy_excl_hypo_late_onset_tmle",
  colors = c("#00274C", "#C4622D"),
  use_jama_style = TRUE,
  plot_type = "single"
)
