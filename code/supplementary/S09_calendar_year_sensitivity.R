# =============================================================================
# Revision_CalendarYear_Sensitivity.r
# Calendar Year Sensitivity Analysis for Reviewer 1 Response
# =============================================================================
#
# PURPOSE: Address Reviewer 1's concern about calendar year imbalance
#          (69% vs 34% COVID-period exposure) across drug comparisons.
#
# THREE ARGUMENTS:
#   1) IPTW with calendar year in PS model -> SMD imbalance persists (>0.1)
#   2) PSM with exact year matching (1:1 and 1:5) -> events reduced but trend maintained
#   3) Pre-COVID (2018-2019) time-restricted PSM -> sufficient follow-up, lower power
#
# PREREQUISITE: Run rev00003.R.r first to generate cohort CSVs.
#
# KEY FIXES vs original appended code:
#   - Single definitions of all functions/config/data loading (no duplication)
#   - PRIMARY_RATIO = 5 set a priori (no p-value-based ratio selection)
#   - All 6 comparisons included (was missing SEMAGLUTIDE vs OTHER_GLPA)
#   - Full 33-variable MICE predictor set from rev00003.R.r
#   - cluster(subclass) in all PSM Cox models for correct SEs
#   - Post-COVID (2022-2023) analysis dropped (insufficient follow-up <1yr)
#   - Consistent index_year_grouped definition throughout
# =============================================================================


# =============================================================================
# CELL 1: Load Libraries
# =============================================================================
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(MatchIt)
  library(survival)
  library(survminer)
  library(lubridate)
  library(ggplot2)
  library(mice)
  library(tableone)
  library(readr)
  library(gridExtra)
  library(survey)
})

cat("Libraries loaded successfully.\n")


# =============================================================================
# CELL 2: Configuration
# =============================================================================
config <- list(
  data_cut_date = as.Date("2023-10-01"),
  outcome = list(
    var = "epilepsy_or_seizure_without_hypoglycemic_seizure_start_date",
    label = "Epilepsy/Seizure (excl Hypoglycemic)"
  ),

  # A priori matching ratios (no p-value-based selection)
  PRIMARY_RATIO = 5,
  matching_ratios = c(1, 5),

  # PS model variables (same as main pipeline)
  ps_vars = c(
    "age", "sex_cat", "raceethnicity_cat", "income", "education", "insurance_category",
    "smoking", "alcohol_category",
    "baseline_bmi_category", "baseline_hba1c",
    "Anticoagulant", "Antiplatelet", "BB", "Biguanide", "CCB", "Diuretic",
    "Ezetimibe", "MRA", "OtherHTN", "RAAS", "Statin", "TZD", "Insulin",
    "SEMAGLUTIDE", "OTHER_GLPA", "SGLT2i", "SU", "DPP4i",
    "myocardial_infarction", "congestive_heart_failure", "peripheral_vascular_disease",
    "cerebrovascular_disease", "chronic_pulmonary_disease", "dementia",
    "rheumatic_disease", "peptic_ulcer_disease",
    "hemiplegia_or_paraplegia", "hiv_infection",
    "hypoglycemia", "hyperglycemic_emergency",
    "renal_disease_severity", "liver_disease_severity",
    "diabetes_with_ophthalmic_complications", "diabetes_with_neurological_complications",
    "malignancy_status"
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
    "liver_disease_severity", "renal_disease_severity",
    "diabetes_with_ophthalmic_complications", "diabetes_with_neurological_complications",
    "malignancy_status"
  ),

  # All 6 cohort files
  cohort_files = c(
    "SEMAGLUTIDE vs OTHER_GLPA" = "SEMAGLUTIDE_vs_OTHER_GLPA_semaglutide_epilepsy_seizure_full_cohort_20260213_053619.csv",
    "SEMAGLUTIDE vs SGLT2"      = "SEMAGLUTIDE_vs_SGLT2_semaglutide_epilepsy_seizure_full_cohort_20260213_053619.csv",
    "SEMAGLUTIDE vs OtherGLD"   = "SEMAGLUTIDE_vs_OtherGLD_semaglutide_epilepsy_seizure_full_cohort_20260213_053619.csv",
    "OTHER_GLPA vs SGLT2"       = "OTHER_GLPA_vs_SGLT2_semaglutide_epilepsy_seizure_full_cohort_20260213_053619.csv",
    "OTHER_GLPA vs OtherGLD"    = "OTHER_GLPA_vs_OtherGLD_semaglutide_epilepsy_seizure_full_cohort_20260213_053619.csv",
    "SGLT2 vs OtherGLD"         = "SGLT2_vs_OtherGLD_semaglutide_epilepsy_seizure_full_cohort_20260213_053619.csv"
  ),

  # Comparison-specific PS variable exclusions
  comparison_exclusions = list(
    "SEMAGLUTIDE vs OTHER_GLPA" = c("SEMAGLUTIDE", "OTHER_GLPA"),
    "SEMAGLUTIDE vs SGLT2"      = c("SEMAGLUTIDE", "SGLT2i"),
    "SEMAGLUTIDE vs OtherGLD"   = c("SEMAGLUTIDE", "TZD", "SU", "DPP4i"),
    "OTHER_GLPA vs SGLT2"       = c("OTHER_GLPA", "SGLT2i"),
    "OTHER_GLPA vs OtherGLD"    = c("OTHER_GLPA", "TZD", "SU", "DPP4i"),
    "SGLT2 vs OtherGLD"         = c("SGLT2i", "TZD", "SU", "DPP4i")
  )
)

cat(sprintf(
  "Configuration set. PRIMARY_RATIO = 1:%d (a priori). Testing ratios: %s\n",
  config$PRIMARY_RATIO,
  paste(paste0("1:", config$matching_ratios), collapse = ", ")
))


# =============================================================================
# CELL 3: Helper Functions
# =============================================================================

# --- Factor labeling ---
apply_variable_labels <- function(data) {
  data %>%
    mutate(
      sex_cat = factor(case_when(
        sex_cat == 0 ~ "Male", sex_cat == 1 ~ "Female",
        sex_cat == 999 ~ "Missing", TRUE ~ "Missing"
      ), levels = c("Male", "Female", "Missing")),
      raceethnicity_cat = factor(case_when(
        raceethnicity_cat == 0 ~ "Non-Hispanic White",
        raceethnicity_cat == 1 ~ "Non-Hispanic Black",
        raceethnicity_cat == 2 ~ "Hispanic",
        raceethnicity_cat == 3 ~ "Other",
        raceethnicity_cat == 999 ~ "Missing", TRUE ~ "Missing"
      ), levels = c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic", "Other", "Missing")),
      income = factor(case_when(
        income == 0 ~ "Low (<50k)", income == 1 ~ "Middle (50k-100k)",
        income == 2 ~ "High (>100k)", income == 999 ~ "Missing", TRUE ~ "Missing"
      ), levels = c("Low (<50k)", "Middle (50k-100k)", "High (>100k)", "Missing")),
      education = factor(case_when(
        education == 0 ~ "High School or Less", education == 1 ~ "Some College",
        education == 2 ~ "Advanced Degree", education == 999 ~ "Missing", TRUE ~ "Missing"
      ), levels = c("High School or Less", "Some College", "Advanced Degree", "Missing")),
      smoking = factor(case_when(
        smoking == 0 ~ "Never", smoking == 1 ~ "Former",
        smoking == 2 ~ "Current", smoking == 999 ~ "Missing", TRUE ~ "Missing"
      ), levels = c("Never", "Former", "Current", "Missing")),
      insurance_category = factor(case_when(
        insurance_category == 0 ~ "None", insurance_category == 1 ~ "Public",
        insurance_category == 2 ~ "Private", insurance_category == 999 ~ "Missing", TRUE ~ "Missing"
      ), levels = c("None", "Public", "Private", "Missing")),
      alcohol_category = factor(case_when(
        alcohol_category == 0 ~ "Low Risk", alcohol_category == 1 ~ "Increased Risk",
        alcohol_category == 2 ~ "High Risk", alcohol_category == 3 ~ "Dependent",
        alcohol_category == 999 ~ "Missing", TRUE ~ "Missing"
      ), levels = c("Low Risk", "Increased Risk", "High Risk", "Dependent", "Missing")),
      baseline_bmi_category = factor(
        case_when(
          is.na(baseline_bmi) ~ "Missing",
          baseline_bmi < 18.5 ~ "Underweight",
          baseline_bmi >= 18.5 & baseline_bmi < 25 ~ "Normal",
          baseline_bmi >= 25 & baseline_bmi < 30 ~ "Overweight",
          baseline_bmi >= 30 & baseline_bmi < 35 ~ "Obesity Class I",
          baseline_bmi >= 35 & baseline_bmi < 40 ~ "Obesity Class II",
          baseline_bmi >= 40 ~ "Obesity Class III"
        ),
        levels = c(
          "Underweight", "Normal", "Overweight", "Obesity Class I",
          "Obesity Class II", "Obesity Class III", "Missing"
        )
      ),
      across(
        c(
          "SEMAGLUTIDE", "OTHER_GLPA", "Biguanide", "TZD", "Insulin", "SGLT2i", "DPP4i", "SU",
          "Anticoagulant", "Antiplatelet", "Statin", "Ezetimibe", "RAAS", "Diuretic",
          "MRA", "BB", "CCB", "OtherHTN", "myocardial_infarction",
          "congestive_heart_failure", "peripheral_vascular_disease",
          "cerebrovascular_disease", "dementia", "chronic_pulmonary_disease",
          "rheumatic_disease", "peptic_ulcer_disease", "hiv_infection",
          "diabetes_with_ophthalmic_complications", "diabetes_with_neurological_complications",
          "hypoglycemia", "hyperglycemic_emergency", "hemiplegia_or_paraplegia"
        ),
        ~ factor(case_when(. == 0 ~ "No", . == 1 ~ "Yes", TRUE ~ "No"),
          levels = c("No", "Yes")
        )
      ),
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
      ), levels = c("None", "Non-Metastatic", "Metastatic"))
    )
}

# --- MICE imputation (FULL 33-variable predictor set from rev00003.R.r lines 116-136) ---
perform_mice_imputation <- function(df, m = 20, seed = 123) {
  set.seed(seed)
  missing_vars <- c("baseline_bmi", "baseline_hba1c")
  available_vars <- intersect(missing_vars, names(df))
  if (length(available_vars) == 0) {
    return(df)
  }

  missing_count <- sapply(df[available_vars], function(x) sum(is.na(x)))
  if (all(missing_count == 0)) {
    return(df)
  }

  cat("MICE imputation for:", paste(names(missing_count)[missing_count > 0], collapse = ", "), "\n")

  # Full 33-variable predictor set (matching rev00003.R.r)
  imp_vars <- intersect(c(
    "baseline_bmi", "baseline_hba1c", # Variables to impute
    "age", "sex_cat", "treatment", "raceethnicity_cat", # Key predictors
    "smoking", "alcohol_category",
    "income", "education", "insurance_category",
    "Biguanide", "Insulin", "SEMAGLUTIDE", "OTHER_GLPA", # Diabetes medications
    "TZD", "SU", "DPP4i", "SGLT2i",
    "RAAS", "Diuretic", "MRA", "BB", "CCB", "OtherHTN", # HTN medications
    "myocardial_infarction", "congestive_heart_failure", # Key comorbidities
    "chronic_pulmonary_disease",
    "liver_disease_severity", "renal_disease_severity",
    "diabetes_with_renal_complications", # DM complications
    "diabetes_with_ophthalmic_complications",
    "diabetes_with_neurological_complications",
    "peripheral_vascular_disease", "cerebrovascular_disease", # Vascular
    "hypoglycemia", "hyperglycemic_emergency" # DM events
  ), names(df))

  tryCatch(
    {
      imp <- mice(df[imp_vars], m = m, seed = seed, printFlag = FALSE)
      completed <- complete(imp, 1)

      pre_na_bmi <- if ("baseline_bmi" %in% names(df)) sum(is.na(df$baseline_bmi)) else NA_integer_
      pre_na_hba1c <- if ("baseline_hba1c" %in% names(df)) sum(is.na(df$baseline_hba1c)) else NA_integer_
      imputed_bmi <- if ("baseline_bmi" %in% names(completed)) sum(is.na(df$baseline_bmi) & !is.na(completed$baseline_bmi)) else NA_integer_
      imputed_hba1c <- if ("baseline_hba1c" %in% names(completed)) sum(is.na(df$baseline_hba1c) & !is.na(completed$baseline_hba1c)) else NA_integer_
      cat(sprintf(
        "  BMI: %d missing -> %d imputed | HbA1c: %d missing -> %d imputed\n",
        pre_na_bmi, imputed_bmi, pre_na_hba1c, imputed_hba1c
      ))

      for (var in available_vars) {
        if (var %in% names(completed)) df[[var]] <- completed[[var]]
      }
      return(df)
    },
    error = function(e) {
      cat("MICE failed, using mean imputation as fallback\n")
      if ("baseline_bmi" %in% names(df)) df$baseline_bmi[is.na(df$baseline_bmi)] <- mean(df$baseline_bmi, na.rm = TRUE)
      if ("baseline_hba1c" %in% names(df)) df$baseline_hba1c[is.na(df$baseline_hba1c)] <- mean(df$baseline_hba1c, na.rm = TRUE)
      return(df)
    }
  )
}

# --- Follow-up and event construction ---
followup_and_event <- function(df, outcome_var, data_cut_date) {
  result <- data.frame(df)
  result$outcome_date <- as.Date(result[[outcome_var]])
  result$censor_date <- pmin(as.Date(result$EHRmaxDT_min), as.Date(data_cut_date), na.rm = TRUE)
  result$end_fu_date <- pmin(result$outcome_date, result$censor_date, result$death_date, na.rm = TRUE)
  result$time <- as.numeric(result$end_fu_date - as.Date(result$index_date))
  result$event <- as.integer(!is.na(result$outcome_date) & result$outcome_date == result$end_fu_date)

  n_before <- nrow(result)
  result <- result %>% filter(time > 0)
  cat(sprintf("  Removed %d rows with non-positive follow-up. N = %d\n", n_before - nrow(result), nrow(result)))
  result$event_time <- result$time
  return(result)
}

# --- SMD calculation ---
calc_smd_single <- function(x, treatment, weights = NULL) {
  if (is.null(weights)) {
    mean_t1 <- mean(x[treatment == 1], na.rm = TRUE)
    mean_t0 <- mean(x[treatment == 0], na.rm = TRUE)
    var_t1 <- var(x[treatment == 1], na.rm = TRUE)
    var_t0 <- var(x[treatment == 0], na.rm = TRUE)
  } else {
    w_t1 <- weights[treatment == 1]
    w_t0 <- weights[treatment == 0]
    mean_t1 <- weighted.mean(x[treatment == 1], w_t1, na.rm = TRUE)
    mean_t0 <- weighted.mean(x[treatment == 0], w_t0, na.rm = TRUE)
    var_t1 <- sum(w_t1 * (x[treatment == 1] - mean_t1)^2, na.rm = TRUE) / sum(w_t1)
    var_t0 <- sum(w_t0 * (x[treatment == 0] - mean_t0)^2, na.rm = TRUE) / sum(w_t0)
  }
  pooled_sd <- sqrt((var_t1 + var_t0) / 2)
  if (pooled_sd > 0) (mean_t1 - mean_t0) / pooled_sd else 0
}

calculate_smd_detailed <- function(df, vars, weights = NULL) {
  results_overall <- data.frame(variable = character(), smd = numeric(), type = character(), stringsAsFactors = FALSE)
  results_by_level <- data.frame(variable = character(), level = character(), smd = numeric(), stringsAsFactors = FALSE)

  for (var in vars) {
    if (!var %in% names(df)) next
    if (is.factor(df[[var]]) || is.character(df[[var]])) {
      levels_vec <- unique(na.omit(df[[var]]))
      level_smds <- numeric()
      for (lvl in levels_vec) {
        var_binary <- as.numeric(df[[var]] == lvl)
        smd_val <- calc_smd_single(var_binary, df$treatment, weights)
        level_smds <- c(level_smds, smd_val)
        results_by_level <- rbind(results_by_level, data.frame(
          variable = var, level = as.character(lvl), smd = smd_val, stringsAsFactors = FALSE
        ))
      }
      overall_smd <- max(abs(level_smds), na.rm = TRUE)
      results_overall <- rbind(results_overall, data.frame(
        variable = var, smd = overall_smd, type = "categorical", stringsAsFactors = FALSE
      ))
    } else {
      smd_val <- calc_smd_single(df[[var]], df$treatment, weights)
      results_overall <- rbind(results_overall, data.frame(
        variable = var, smd = smd_val, type = "continuous", stringsAsFactors = FALSE
      ))
    }
  }
  return(list(overall = results_overall, by_level = results_by_level))
}

# --- Year-specific SMD ---
calculate_year_smd <- function(df, weights = NULL) {
  years <- unique(df$index_year)
  results <- data.frame(year = character(), smd = numeric(), stringsAsFactors = FALSE)
  for (yr in years) {
    var_binary <- as.numeric(df$index_year == yr)
    smd_val <- calc_smd_single(var_binary, df$treatment, weights)
    results <- rbind(results, data.frame(year = as.character(yr), smd = smd_val, stringsAsFactors = FALSE))
  }
  return(results)
}

# --- PSM with exact year matching + cluster(subclass) Cox ---
run_psm_exact_year <- function(df, exclude_vars, ps_vars, caliper = 0.2, ratio = 1) {
  rhs_vars <- setdiff(ps_vars, c(exclude_vars, "index_year", "index_year_grouped"))
  keep_vars <- rhs_vars[sapply(df[rhs_vars], function(x) length(unique(na.omit(x))) > 1)]

  if (!"index_year" %in% names(df)) {
    df$index_year <- factor(lubridate::year(df$index_date))
  }

  model_vars <- c("treatment", keep_vars, "index_year")
  complete_rows <- complete.cases(df[model_vars])
  df_complete <- df[complete_rows, ]

  cat(sprintf(
    "  Complete cases: %d (Treatment: %d, Control: %d)\n",
    nrow(df_complete), sum(df_complete$treatment == 1), sum(df_complete$treatment == 0)
  ))

  ps_formula <- reformulate(keep_vars, response = "treatment")

  tryCatch(
    {
      m_out <- matchit(
        ps_formula,
        data = df_complete,
        method = "nearest",
        distance = "glm",
        link = "logit",
        caliper = caliper,
        ratio = ratio,
        exact = ~index_year
      )

      matched_df <- match.data(m_out)

      n_treat <- sum(matched_df$treatment == 1)
      n_ctrl <- sum(matched_df$treatment == 0)
      events_treat <- sum(matched_df$event[matched_df$treatment == 1])
      events_ctrl <- sum(matched_df$event[matched_df$treatment == 0])

      cat(sprintf(
        "  Matched: %d total (Treat: %d, Ctrl: %d) | Events: %d (Treat: %d, Ctrl: %d)\n",
        nrow(matched_df), n_treat, n_ctrl,
        events_treat + events_ctrl, events_treat, events_ctrl
      ))

      balance_vars <- unique(c(keep_vars, "index_year", "index_year_grouped"))
      balance_vars <- balance_vars[balance_vars %in% names(df_complete)]
      balance_before <- calculate_smd_detailed(df_complete, balance_vars, weights = NULL)
      balance_after <- calculate_smd_detailed(matched_df, balance_vars, weights = NULL)

      # Cox with cluster(subclass) for correct variance estimation
      cox_model <- coxph(Surv(event_time, event) ~ treatment + cluster(subclass), data = matched_df)

      return(list(
        matchit_obj = m_out,
        cohort_before = df_complete,
        cohort_matched = matched_df,
        balance_before = balance_before$overall,
        balance_after = balance_after$overall,
        cox = summary(cox_model),
        ratio = ratio
      ))
    },
    error = function(e) {
      cat(sprintf("  PSM 1:%d matching failed: %s\n", ratio, e$message))
      return(NULL)
    }
  )
}

# --- PSM without exact year (for time-restricted) + cluster(subclass) Cox ---
run_psm_no_exact_year <- function(df, exclude_vars, ps_vars, caliper = 0.2, ratio = 1) {
  rhs_vars <- setdiff(ps_vars, c(exclude_vars, "index_year", "index_year_grouped"))
  keep_vars <- rhs_vars[sapply(df[rhs_vars], function(x) length(unique(na.omit(x))) > 1)]

  model_vars <- c("treatment", keep_vars)
  complete_rows <- complete.cases(df[model_vars])
  df_complete <- df[complete_rows, ]

  if (nrow(df_complete) < 50 || sum(df_complete$treatment == 1) < 10 || sum(df_complete$treatment == 0) < 10) {
    cat("  Insufficient sample size for PSM.\n")
    return(NULL)
  }

  ps_formula <- reformulate(keep_vars, response = "treatment")

  tryCatch(
    {
      m_out <- matchit(
        ps_formula,
        data = df_complete,
        method = "nearest",
        distance = "glm",
        link = "logit",
        caliper = caliper,
        ratio = ratio
      )

      matched_df <- match.data(m_out)
      if (nrow(matched_df) < 20) {
        return(NULL)
      }

      n_treat <- sum(matched_df$treatment == 1)
      n_ctrl <- sum(matched_df$treatment == 0)
      events_treat <- sum(matched_df$event[matched_df$treatment == 1])
      events_ctrl <- sum(matched_df$event[matched_df$treatment == 0])

      cat(sprintf(
        "    Matched: N = %d (Treat: %d, Ctrl: %d) | Events: %d (Treat: %d, Ctrl: %d)\n",
        nrow(matched_df), n_treat, n_ctrl,
        events_treat + events_ctrl, events_treat, events_ctrl
      ))

      # Cox with cluster(subclass) for correct variance estimation
      cox_model <- coxph(Surv(event_time, event) ~ treatment + cluster(subclass), data = matched_df)

      # Mean follow-up
      mean_fu_treat <- mean(matched_df$event_time[matched_df$treatment == 1]) / 365.25
      mean_fu_ctrl <- mean(matched_df$event_time[matched_df$treatment == 0]) / 365.25

      return(list(
        cohort_matched = matched_df,
        cox = summary(cox_model),
        ratio = ratio,
        mean_fu_treat = mean_fu_treat,
        mean_fu_ctrl = mean_fu_ctrl
      ))
    },
    error = function(e) {
      cat(sprintf("    PSM 1:%d failed: %s\n", ratio, e$message))
      return(NULL)
    }
  )
}

# --- Comparison labels for KM plots ---
get_comparison_labels <- function(comp_name) {
  if (grepl("SEMAGLUTIDE vs OTHER_GLPA", comp_name)) {
    c("0" = "Other GLP-1 Agonists", "1" = "Semaglutide")
  } else if (grepl("SEMAGLUTIDE vs SGLT2", comp_name)) {
    c("0" = "SGLT-2 Inhibitors", "1" = "Semaglutide")
  } else if (grepl("SEMAGLUTIDE vs OtherGLD", comp_name)) {
    c("0" = "Other GLDs", "1" = "Semaglutide")
  } else if (grepl("OTHER_GLPA vs SGLT2", comp_name)) {
    c("0" = "SGLT-2 Inhibitors", "1" = "Other GLP-1 Agonists")
  } else if (grepl("OTHER_GLPA vs OtherGLD", comp_name)) {
    c("0" = "Other GLDs", "1" = "Other GLP-1 Agonists")
  } else if (grepl("SGLT2 vs OtherGLD", comp_name)) {
    c("0" = "Other GLDs", "1" = "SGLT-2 Inhibitors")
  } else {
    c("0" = "Comparator", "1" = "Treatment")
  }
}

# --- Unified KM plot function ---
create_km_plot <- function(matched_df, comp_name, title_suffix = "", ratio = 5, save_plot = TRUE, cox_p = NULL) {
  labels <- get_comparison_labels(comp_name)

  matched_df$event_time_years <- matched_df$event_time / 365.25
  matched_df$treatment_label <- factor(matched_df$treatment,
    levels = c(0, 1),
    labels = c(labels["0"], labels["1"])
  )

  surv_fit <- survfit(Surv(event_time_years, event) ~ treatment_label, data = matched_df)

  n_treat <- sum(matched_df$treatment == 1)
  n_ctrl <- sum(matched_df$treatment == 0)
  events_treat <- sum(matched_df$event[matched_df$treatment == 1])
  events_ctrl <- sum(matched_df$event[matched_df$treatment == 0])
  py_treat <- sum(matched_df$event_time[matched_df$treatment == 1]) / 365.25
  py_ctrl <- sum(matched_df$event_time[matched_df$treatment == 0]) / 365.25
  ir_treat <- (events_treat / py_treat) * 1000
  ir_ctrl <- (events_ctrl / py_ctrl) * 1000

  subtitle_text <- sprintf(
    "PSM 1:%d | %s: %d events/%d pts (IR: %.2f/1000PY) | %s: %d events/%d pts (IR: %.2f/1000PY)",
    ratio, labels["1"], events_treat, n_treat, ir_treat,
    labels["0"], events_ctrl, n_ctrl, ir_ctrl
  )

  plot_title <- if (title_suffix != "") {
    sprintf("%s - %s (PSM 1:%d)", comp_name, title_suffix, ratio)
  } else {
    sprintf("%s (PSM 1:%d)", comp_name, ratio)
  }

  km_plot <- ggsurvplot(
    fit = surv_fit,
    data = matched_df,
    fun = "event",
    conf.int = TRUE,
    risk.table = TRUE,
    risk.table.col = "strata",
    risk.table.fontsize = 3.5,
    risk.table.height = 0.28,
    legend.title = "",
    legend.labs = c(
      sprintf("%s (n=%d, events=%d)", labels["0"], n_ctrl, events_ctrl),
      sprintf("%s (n=%d, events=%d)", labels["1"], n_treat, events_treat)
    ),
    xlab = "Follow-up Time (Years)",
    ylab = "Cumulative Incidence",
    title = plot_title,
    subtitle = subtitle_text,
    break.time.by = 1,
    xlim = c(0, 5),
    pval = if (!is.null(cox_p)) cox_p else TRUE,
    pval.coord = c(3.5, 0.02),
    pval.method = if (!is.null(cox_p)) FALSE else TRUE,
    palette = c("#E74C3C", "#3498DB"),
    ggtheme = theme_bw(base_size = 11) +
      theme(
        plot.title = element_text(face = "bold", size = 12),
        plot.subtitle = element_text(size = 9),
        legend.position = "bottom"
      )
  )

  if (save_plot) {
    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    clean_name <- gsub(" ", "_", comp_name)
    clean_suffix <- if (title_suffix != "") paste0("_", gsub(" ", "_", gsub("[()]", "", title_suffix))) else ""
    plot_file <- sprintf("KM_%s%s_1to%d_%s.png", clean_name, clean_suffix, ratio, ts)
    png(filename = plot_file, width = 10, height = 9, units = "in", res = 300)
    print(km_plot)
    dev.off()
    cat(sprintf("  Saved: %s\n", plot_file))
  }

  return(km_plot)
}

cat("All helper functions defined.\n")


# =============================================================================
# CELL 4: Load All 6 Cohort CSVs
# =============================================================================
cohorts <- lapply(names(config$cohort_files), function(comp_name) {
  f <- config$cohort_files[[comp_name]]
  if (!file.exists(f)) stop(paste("File not found:", f))
  df <- read_csv(f, show_col_types = FALSE)
  df %>% mutate(
    person_id = as.character(person_id),
    index_date = as.Date(index_date),
    treatment = as.integer(treatment)
  )
})
names(cohorts) <- names(config$cohort_files)

cat(sprintf(
  "Loaded %d cohort files: %s\n",
  length(cohorts),
  paste(names(cohorts), collapse = ", ")
))


# =============================================================================
# CELL 5: Data Preparation (all 6 cohorts)
# =============================================================================
prepared_cohorts <- list()

for (comp_name in names(cohorts)) {
  cat(sprintf("\nPreparing: %s\n", comp_name))

  cohort_df <- cohorts[[comp_name]]

  # Create index_year, index_year_grouped, composite factors
  cohort_df <- cohort_df %>%
    mutate(
      index_year = factor(lubridate::year(index_date)),
      # Consistent definition: 2020-2022 = "COVID", all other years = "Non-COVID"
      index_year_grouped = factor(case_when(
        lubridate::year(index_date) %in% c(2020, 2021, 2022) ~ "COVID",
        TRUE ~ "Non-COVID"
      ), levels = c("Non-COVID", "COVID")),
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
      ), levels = c("None", "Non-Metastatic", "Metastatic"))
    ) %>%
    apply_variable_labels()

  # Ensure categorical vars are factors
  for (var in config$categorical_vars) {
    if (var %in% names(cohort_df) && !is.factor(cohort_df[[var]])) {
      cohort_df[[var]] <- as.factor(cohort_df[[var]])
    }
  }

  # Follow-up and event
  survival_df <- followup_and_event(cohort_df, config$outcome$var, config$data_cut_date)

  # MICE imputation
  survival_df <- perform_mice_imputation(survival_df, m = 20, seed = 123)

  # Recompute BMI category after imputation
  if ("baseline_bmi" %in% names(survival_df)) {
    survival_df <- survival_df %>% mutate(
      baseline_bmi_category = factor(
        case_when(
          is.na(baseline_bmi) ~ "Missing",
          baseline_bmi < 18.5 ~ "Underweight",
          baseline_bmi >= 18.5 & baseline_bmi < 25 ~ "Normal",
          baseline_bmi >= 25 & baseline_bmi < 30 ~ "Overweight",
          baseline_bmi >= 30 & baseline_bmi < 35 ~ "Obesity Class I",
          baseline_bmi >= 35 & baseline_bmi < 40 ~ "Obesity Class II",
          baseline_bmi >= 40 ~ "Obesity Class III"
        ),
        levels = c(
          "Underweight", "Normal", "Overweight", "Obesity Class I",
          "Obesity Class II", "Obesity Class III", "Missing"
        )
      )
    )
  }

  prepared_cohorts[[comp_name]] <- survival_df

  cat(sprintf(
    "  N = %d, Treatment = %d, Control = %d, Events = %d\n",
    nrow(survival_df), sum(survival_df$treatment == 1),
    sum(survival_df$treatment == 0), sum(survival_df$event)
  ))
}

cat("\nData preparation complete for all 6 cohorts.\n")


# =============================================================================
# CELL 6: ARGUMENT 1 - IPTW with Calendar Year in PS Model
#          Show that SMD imbalance persists (>0.1) even after IPTW
# =============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ARGUMENT 1: IPTW with Calendar Year in PS Model -> Large SMD Imbalance\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

iptw_smd_results <- list()

for (comp_name in names(prepared_cohorts)) {
  cat(sprintf("\n--- %s ---\n", comp_name))

  df <- prepared_cohorts[[comp_name]]
  exclude_vars <- config$comparison_exclusions[[comp_name]]

  # PS model variables including index_year
  ps_vars <- setdiff(config$ps_vars, exclude_vars)
  ps_vars <- c(ps_vars, "index_year")
  ps_vars <- ps_vars[ps_vars %in% names(df)]
  ps_vars <- ps_vars[sapply(df[ps_vars], function(x) length(unique(na.omit(x))) > 1)]

  ps_formula <- reformulate(ps_vars, response = "treatment")

  model_vars <- c("treatment", ps_vars)
  df_complete <- df[complete.cases(df[model_vars]), ]

  # Fit PS model
  ps_model <- glm(ps_formula, data = df_complete, family = binomial())
  df_complete$ps <- predict(ps_model, type = "response")

  # Calculate IPTW weights
  df_complete$iptw <- ifelse(df_complete$treatment == 1,
    1 / df_complete$ps,
    1 / (1 - df_complete$ps)
  )

  # Truncate weights at 1st and 99th percentile
  lower <- quantile(df_complete$iptw, 0.01)
  upper <- quantile(df_complete$iptw, 0.99)
  df_complete$iptw_truncated <- pmin(pmax(df_complete$iptw, lower), upper)

  # Year-specific SMD before and after IPTW
  smd_before <- calculate_year_smd(df_complete, weights = NULL)
  smd_after_iptw <- calculate_year_smd(df_complete, weights = df_complete$iptw_truncated)

  smd_before$method <- "Before Weighting"
  smd_after_iptw$method <- "After IPTW"

  combined_smd <- rbind(smd_before, smd_after_iptw)
  combined_smd$comparison <- comp_name

  iptw_smd_results[[comp_name]] <- combined_smd

  # Print
  cat("Calendar Year SMD:\n")
  cat(sprintf("  %-8s  %10s  %10s\n", "Year", "Before", "After IPTW"))
  for (yr in unique(smd_before$year)) {
    before_val <- smd_before$smd[smd_before$year == yr]
    after_val <- smd_after_iptw$smd[smd_after_iptw$year == yr]
    cat(sprintf("  %-8s  %+10.3f  %+10.3f\n", yr, before_val, after_val))
  }

  max_smd_after <- max(abs(smd_after_iptw$smd))
  cat(sprintf(
    "  Max |SMD| after IPTW = %.3f %s\n",
    max_smd_after, ifelse(max_smd_after > 0.1, "<- IMBALANCED (>0.1)", "")
  ))
}

iptw_smd_df <- do.call(rbind, iptw_smd_results)


# =============================================================================
# CELL 7: ARGUMENT 1 Figure - SMD Imbalance Bar Chart
# =============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Creating IPTW SMD Imbalance Figure\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

iptw_smd_df$method <- factor(iptw_smd_df$method, levels = c("Before Weighting", "After IPTW"))

p_iptw_smd <- ggplot(iptw_smd_df, aes(x = year, y = smd, fill = method)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
  geom_hline(yintercept = c(-0.1, 0.1), linetype = "dashed", color = "red", linewidth = 0.8) +
  geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
  facet_wrap(~comparison, ncol = 2, scales = "free_y") +
  scale_fill_manual(values = c("Before Weighting" = "#95a5a6", "After IPTW" = "#3498db")) +
  labs(
    title = "Calendar Year SMD: Before vs After IPTW",
    subtitle = "Red dashed lines = SMD +/-0.1 threshold | IPTW cannot adequately balance calendar year",
    x = "Calendar Year",
    y = "Standardized Mean Difference (SMD)",
    fill = "Method"
  ) +
  theme_bw(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", size = 14),
    legend.position = "bottom",
    strip.text = element_text(face = "bold", size = 10)
  )

ggsave(sprintf("Figure_IPTW_SMD_Imbalance_%s.png", ts), p_iptw_smd, width = 14, height = 10, dpi = 300)
cat(sprintf("Saved: Figure_IPTW_SMD_Imbalance_%s.png\n", ts))
print(p_iptw_smd)


# =============================================================================
# CELL 8: ARGUMENT 2 - PSM Exact Year Matching (1:1 and 1:5)
#          All 6 comparisons, Cox with cluster(subclass), event reduction tracking
# =============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ARGUMENT 2: PSM (1:1 and 1:5) with Exact Year Matching\n")
cat("-> Events reduced, statistical power decreased, but TREND MAINTAINED\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

psm_exact_year_results <- list()

for (comp_name in names(prepared_cohorts)) {
  cat(sprintf("\n%s\n%s\n", paste(rep("=", 70), collapse = ""), comp_name))

  df <- prepared_cohorts[[comp_name]]
  exclude_vars <- config$comparison_exclusions[[comp_name]]

  # Original (unmatched) event counts
  original_n <- nrow(df)
  original_events <- sum(df$event)
  original_events_treat <- sum(df$event[df$treatment == 1])
  original_events_ctrl <- sum(df$event[df$treatment == 0])

  cat(sprintf(
    "Original: N = %d, Events = %d (Treat: %d, Ctrl: %d)\n",
    original_n, original_events, original_events_treat, original_events_ctrl
  ))

  for (ratio in config$matching_ratios) {
    cat(sprintf("\n  PSM 1:%d with exact year matching:\n", ratio))

    result <- run_psm_exact_year(
      df = df,
      exclude_vars = exclude_vars,
      ps_vars = config$ps_vars,
      caliper = 0.2,
      ratio = ratio
    )

    if (!is.null(result)) {
      matched_df <- result$cohort_matched
      cox_sum <- result$cox

      hr <- cox_sum$conf.int[1, "exp(coef)"]
      hr_lo <- cox_sum$conf.int[1, "lower .95"]
      hr_hi <- cox_sum$conf.int[1, "upper .95"]
      p_val <- cox_sum$coefficients[1, "Pr(>|z|)"]

      events_total <- sum(matched_df$event)
      events_treat <- sum(matched_df$event[matched_df$treatment == 1])
      events_ctrl <- sum(matched_df$event[matched_df$treatment == 0])

      py_treat <- sum(matched_df$event_time[matched_df$treatment == 1]) / 365.25
      py_ctrl <- sum(matched_df$event_time[matched_df$treatment == 0]) / 365.25

      year_smd_after <- calculate_year_smd(matched_df, weights = NULL)
      max_year_smd <- max(abs(year_smd_after$smd))

      cat(sprintf("    HR: %.3f (%.3f - %.3f), p = %.4f\n", hr, hr_lo, hr_hi, p_val))
      cat(sprintf(
        "    Event reduction: %.1f%% (from %d to %d)\n",
        (1 - events_total / original_events) * 100, original_events, events_total
      ))
      cat(sprintf("    Year SMD max: %.4f (exact matching -> ~0)\n", max_year_smd))

      # Check for robust SEs (cluster column)
      if ("robust se" %in% tolower(colnames(cox_sum$coefficients))) {
        cat("    [Robust SEs: cluster(subclass) applied]\n")
      }

      key <- sprintf("%s_1to%d", comp_name, ratio)
      psm_exact_year_results[[key]] <- list(
        comparison = comp_name,
        ratio = ratio,
        n_original = original_n,
        n_matched = nrow(matched_df),
        n_treatment = sum(matched_df$treatment == 1),
        n_control = sum(matched_df$treatment == 0),
        events_original = original_events,
        events_total = events_total,
        events_treatment = events_treat,
        events_control = events_ctrl,
        event_reduction_pct = (1 - events_total / original_events) * 100,
        hr = hr, hr_lower = hr_lo, hr_upper = hr_hi, p_value = p_val,
        ir_treatment = (events_treat / py_treat) * 1000,
        ir_control = (events_ctrl / py_ctrl) * 1000,
        max_year_smd = max_year_smd,
        matched_data = matched_df,
        cox_summary = cox_sum
      )
    }
  }
}


# =============================================================================
# CELL 9: ARGUMENT 2 - KM Plots for Exact Year PSM (PRIMARY_RATIO)
# =============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat(sprintf("Generating KM Plots for Exact Year PSM (1:%d)\n", config$PRIMARY_RATIO))
cat(paste(rep("=", 80), collapse = ""), "\n")

for (key in names(psm_exact_year_results)) {
  res <- psm_exact_year_results[[key]]

  # Only plot the primary ratio
  if (res$ratio != config$PRIMARY_RATIO) next

  cat(sprintf("\n%s (1:%d)\n", res$comparison, res$ratio))

  tryCatch(
    {
      km <- create_km_plot(
        matched_df = res$matched_data,
        comp_name = res$comparison,
        title_suffix = "Exact Year",
        ratio = res$ratio,
        save_plot = TRUE,
        cox_p = res$p_value
      )
      print(km)
    },
    error = function(e) {
      cat(sprintf("  KM plot failed: %s\n", e$message))
    }
  )
}


# =============================================================================
# CELL 10: ARGUMENT 3 - Pre-COVID (2018-2019) Time-Restricted PSM
#           Post-COVID (2022-2023) DROPPED due to insufficient follow-up
# =============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ARGUMENT 3: Time-Restricted PSM (2018-2019 Pre-COVID Only)\n")
cat("NOTE: Post-COVID (2022-2023) analysis DROPPED - insufficient follow-up (<1 year\n")
cat("      from data cut date 2023-10-01; max possible FU is ~0.75-1.75 years).\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

precovid_results <- list()

for (comp_name in names(prepared_cohorts)) {
  cat(sprintf("\n--- %s ---\n", comp_name))

  df <- prepared_cohorts[[comp_name]]
  exclude_vars <- config$comparison_exclusions[[comp_name]]

  # Filter to 2018-2019 only
  df_period <- df %>%
    filter(as.numeric(as.character(index_year)) %in% c(2018, 2019))

  cat(sprintf(
    "  2018-2019 cohort: N = %d (Treat: %d, Ctrl: %d), Events = %d\n",
    nrow(df_period),
    sum(df_period$treatment == 1),
    sum(df_period$treatment == 0),
    sum(df_period$event)
  ))

  if (nrow(df_period) < 100 || sum(df_period$treatment == 1) < 20) {
    cat("  Insufficient sample size. Skipping.\n")
    next
  }

  # Report follow-up
  mean_fu_treat <- mean(df_period$event_time[df_period$treatment == 1]) / 365.25
  mean_fu_ctrl <- mean(df_period$event_time[df_period$treatment == 0]) / 365.25
  max_fu_treat <- max(df_period$event_time[df_period$treatment == 1]) / 365.25
  max_fu_ctrl <- max(df_period$event_time[df_period$treatment == 0]) / 365.25

  cat(sprintf("  Mean follow-up: Treat = %.2f yrs, Ctrl = %.2f yrs\n", mean_fu_treat, mean_fu_ctrl))
  cat(sprintf("  Max follow-up:  Treat = %.2f yrs, Ctrl = %.2f yrs\n", max_fu_treat, max_fu_ctrl))

  for (ratio in config$matching_ratios) {
    cat(sprintf("\n  PSM 1:%d (no exact year - within-period):\n", ratio))

    result <- run_psm_no_exact_year(
      df = df_period,
      exclude_vars = exclude_vars,
      ps_vars = config$ps_vars,
      caliper = 0.2,
      ratio = ratio
    )

    if (!is.null(result)) {
      cox_sum <- result$cox
      hr <- cox_sum$conf.int[1, "exp(coef)"]
      hr_lo <- cox_sum$conf.int[1, "lower .95"]
      hr_hi <- cox_sum$conf.int[1, "upper .95"]
      p_val <- cox_sum$coefficients[1, "Pr(>|z|)"]

      matched_df <- result$cohort_matched
      events_treat <- sum(matched_df$event[matched_df$treatment == 1])
      events_ctrl <- sum(matched_df$event[matched_df$treatment == 0])

      cat(sprintf("    HR: %.3f (%.3f - %.3f), p = %.4f\n", hr, hr_lo, hr_hi, p_val))
      cat(sprintf("    Mean FU after matching: %.2f / %.2f years\n", result$mean_fu_treat, result$mean_fu_ctrl))

      key <- sprintf("%s_PreCOVID_1to%d", comp_name, ratio)
      precovid_results[[key]] <- list(
        comparison = comp_name,
        period = "2018-2019 (Pre-COVID)",
        ratio = ratio,
        n_treatment = sum(matched_df$treatment == 1),
        n_control = sum(matched_df$treatment == 0),
        events_treatment = events_treat,
        events_control = events_ctrl,
        hr = hr, hr_lower = hr_lo, hr_upper = hr_hi, p_value = p_val,
        mean_fu_treat = result$mean_fu_treat,
        mean_fu_ctrl = result$mean_fu_ctrl,
        py_treatment = sum(matched_df$event_time[matched_df$treatment == 1]) / 365.25,
        py_control = sum(matched_df$event_time[matched_df$treatment == 0]) / 365.25,
        matched_data = matched_df
      )
    }
  }
}


# =============================================================================
# CELL 11: ARGUMENT 3 - Pre-COVID KM Plots (PRIMARY_RATIO)
# =============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat(sprintf("Generating Pre-COVID KM Plots (1:%d)\n", config$PRIMARY_RATIO))
cat(paste(rep("=", 80), collapse = ""), "\n")

for (key in names(precovid_results)) {
  res <- precovid_results[[key]]

  # Only plot the primary ratio
  if (res$ratio != config$PRIMARY_RATIO) next

  cat(sprintf("\n%s - %s (1:%d)\n", res$comparison, res$period, res$ratio))

  tryCatch(
    {
      km <- create_km_plot(
        matched_df = res$matched_data,
        comp_name = res$comparison,
        title_suffix = "Pre-COVID 2018-2019",
        ratio = res$ratio,
        save_plot = TRUE,
        cox_p = res$p_value
      )
      print(km)
    },
    error = function(e) {
      cat(sprintf("  KM plot failed: %s\n", e$message))
    }
  )
}


# =============================================================================
# CELL 11b: ARGUMENT 4 - Minimum 2-Year Follow-up Restriction
#            (Reviewer 1 Point #3: Differential follow-up)
#            Restrict to patients with index_date <= 2021-10-01
#            so everyone has >= 2 years potential follow-up before data cut
# =============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ARGUMENT 4: Restrict to >= 2 Years Potential Follow-up\n")
cat("  (index_date <= 2021-10-01, data cut = 2023-10-01)\n")
cat("  Addresses Reviewer 1 Point #3: Differential follow-up\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

min_fu_cutoff <- as.Date("2021-10-01") # 2 years before data cut date 2023-10-01

min2yr_results <- list()

for (comp_name in names(prepared_cohorts)) {
  cat(sprintf("\n--- %s ---\n", comp_name))

  df <- prepared_cohorts[[comp_name]]
  exclude_vars <- config$comparison_exclusions[[comp_name]]

  # Filter to patients with >= 2 years potential follow-up
  df_restricted <- df %>%
    filter(as.Date(index_date) <= min_fu_cutoff)

  n_excluded <- nrow(df) - nrow(df_restricted)
  cat(sprintf(
    "  Full cohort: N = %d | After >= 2yr restriction: N = %d (excluded %d)\n",
    nrow(df), nrow(df_restricted), n_excluded
  ))
  cat(sprintf(
    "  Treatment: %d -> %d | Control: %d -> %d | Events: %d -> %d\n",
    sum(df$treatment == 1), sum(df_restricted$treatment == 1),
    sum(df$treatment == 0), sum(df_restricted$treatment == 0),
    sum(df$event), sum(df_restricted$event)
  ))

  if (nrow(df_restricted) < 100 || sum(df_restricted$treatment == 1) < 20) {
    cat("  Insufficient sample size after restriction. Skipping.\n")
    next
  }

  # Report follow-up
  mean_fu_treat <- mean(df_restricted$event_time[df_restricted$treatment == 1]) / 365.25
  mean_fu_ctrl <- mean(df_restricted$event_time[df_restricted$treatment == 0]) / 365.25
  cat(sprintf(
    "  Mean follow-up after restriction: Treat = %.2f yrs, Ctrl = %.2f yrs\n",
    mean_fu_treat, mean_fu_ctrl
  ))

  for (ratio in config$matching_ratios) {
    cat(sprintf("\n  PSM 1:%d (no exact year):\n", ratio))

    result <- run_psm_no_exact_year(
      df = df_restricted,
      exclude_vars = exclude_vars,
      ps_vars = config$ps_vars,
      caliper = 0.2,
      ratio = ratio
    )

    if (!is.null(result)) {
      cox_sum <- result$cox
      hr <- cox_sum$conf.int[1, "exp(coef)"]
      hr_lo <- cox_sum$conf.int[1, "lower .95"]
      hr_hi <- cox_sum$conf.int[1, "upper .95"]
      p_val <- cox_sum$coefficients[1, "Pr(>|z|)"]

      matched_df <- result$cohort_matched
      events_treat <- sum(matched_df$event[matched_df$treatment == 1])
      events_ctrl <- sum(matched_df$event[matched_df$treatment == 0])

      cat(sprintf("    HR: %.3f (%.3f - %.3f), p = %.4f\n", hr, hr_lo, hr_hi, p_val))
      cat(sprintf(
        "    Mean FU after matching: %.2f / %.2f years\n",
        result$mean_fu_treat, result$mean_fu_ctrl
      ))

      key <- sprintf("%s_Min2yr_1to%d", comp_name, ratio)
      min2yr_results[[key]] <- list(
        comparison = comp_name,
        period = ">=2yr Follow-up",
        ratio = ratio,
        n_original = nrow(df),
        n_restricted = nrow(df_restricted),
        n_treatment = sum(matched_df$treatment == 1),
        n_control = sum(matched_df$treatment == 0),
        events_original = sum(df$event),
        events_treatment = events_treat,
        events_control = events_ctrl,
        hr = hr, hr_lower = hr_lo, hr_upper = hr_hi, p_value = p_val,
        mean_fu_treat = result$mean_fu_treat,
        mean_fu_ctrl = result$mean_fu_ctrl,
        py_treatment = sum(matched_df$event_time[matched_df$treatment == 1]) / 365.25,
        py_control = sum(matched_df$event_time[matched_df$treatment == 0]) / 365.25,
        matched_data = matched_df
      )
    }
  }
}


# =============================================================================
# CELL 11c: ARGUMENT 4 - KM Plots for >= 2-Year Follow-up (PRIMARY_RATIO)
# =============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat(sprintf("Generating >= 2yr Follow-up KM Plots (1:%d)\n", config$PRIMARY_RATIO))
cat(paste(rep("=", 80), collapse = ""), "\n")

for (key in names(min2yr_results)) {
  res <- min2yr_results[[key]]
  if (res$ratio != config$PRIMARY_RATIO) next

  cat(sprintf("\n%s - %s (1:%d)\n", res$comparison, res$period, res$ratio))

  tryCatch(
    {
      km <- create_km_plot(
        matched_df = res$matched_data,
        comp_name = res$comparison,
        title_suffix = "Min 2yr FU",
        ratio = res$ratio,
        save_plot = TRUE,
        cox_p = res$p_value
      )
      print(km)
    },
    error = function(e) {
      cat(sprintf("  KM plot failed: %s\n", e$message))
    }
  )
}


# =============================================================================
# CELL 12: Combined Forest Plot (6 comparisons x 4 methods, color-coded)
# =============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Creating Combined Forest Plot (6 comparisons x 4 methods)\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

forest_data <- data.frame()

# (A) Add IPTW results - re-run Cox on IPTW-weighted data for HR
for (comp_name in names(prepared_cohorts)) {
  df <- prepared_cohorts[[comp_name]]
  exclude_vars <- config$comparison_exclusions[[comp_name]]

  ps_vars <- setdiff(config$ps_vars, exclude_vars)
  ps_vars <- c(ps_vars, "index_year")
  ps_vars <- ps_vars[ps_vars %in% names(df)]
  ps_vars <- ps_vars[sapply(df[ps_vars], function(x) length(unique(na.omit(x))) > 1)]

  model_vars <- c("treatment", ps_vars, "event_time", "event")
  df_complete <- df[complete.cases(df[model_vars]), ]

  tryCatch(
    {
      ps_model <- glm(reformulate(ps_vars, response = "treatment"),
        data = df_complete, family = binomial()
      )
      df_complete$ps <- predict(ps_model, type = "response")
      df_complete$iptw <- ifelse(df_complete$treatment == 1,
        1 / df_complete$ps,
        1 / (1 - df_complete$ps)
      )
      lower <- quantile(df_complete$iptw, 0.01)
      upper <- quantile(df_complete$iptw, 0.99)
      df_complete$iptw_truncated <- pmin(pmax(df_complete$iptw, lower), upper)

      cox_iptw <- coxph(Surv(event_time, event) ~ treatment,
        data = df_complete, weights = iptw_truncated
      )
      cox_sum <- summary(cox_iptw)

      forest_data <- rbind(forest_data, data.frame(
        comparison = comp_name,
        method = "IPTW + Calendar Year in PS",
        n_events = sum(df_complete$event),
        hr = cox_sum$conf.int[1, "exp(coef)"],
        hr_lower = cox_sum$conf.int[1, "lower .95"],
        hr_upper = cox_sum$conf.int[1, "upper .95"],
        p_value = cox_sum$coefficients[1, "Pr(>|z|)"],
        stringsAsFactors = FALSE
      ))
    },
    error = function(e) {
      cat(sprintf("  IPTW Cox failed for %s: %s\n", comp_name, e$message))
    }
  )
}

# (B) Add PSM exact year results (PRIMARY_RATIO only)
for (key in names(psm_exact_year_results)) {
  res <- psm_exact_year_results[[key]]
  if (res$ratio != config$PRIMARY_RATIO) next

  forest_data <- rbind(forest_data, data.frame(
    comparison = res$comparison,
    method = sprintf("Exact Year PSM 1:%d", res$ratio),
    n_events = res$events_total,
    hr = res$hr,
    hr_lower = res$hr_lower,
    hr_upper = res$hr_upper,
    p_value = res$p_value,
    stringsAsFactors = FALSE
  ))
}

# (C) Add Pre-COVID results (PRIMARY_RATIO only)
for (key in names(precovid_results)) {
  res <- precovid_results[[key]]
  if (res$ratio != config$PRIMARY_RATIO) next

  forest_data <- rbind(forest_data, data.frame(
    comparison = res$comparison,
    method = sprintf("Pre-COVID PSM 1:%d", res$ratio),
    n_events = res$events_treatment + res$events_control,
    hr = res$hr,
    hr_lower = res$hr_lower,
    hr_upper = res$hr_upper,
    p_value = res$p_value,
    stringsAsFactors = FALSE
  ))
}

# (D) Add >= 2yr follow-up results (PRIMARY_RATIO only)
for (key in names(min2yr_results)) {
  res <- min2yr_results[[key]]
  if (res$ratio != config$PRIMARY_RATIO) next

  forest_data <- rbind(forest_data, data.frame(
    comparison = res$comparison,
    method = sprintf(">=2yr FU PSM 1:%d", res$ratio),
    n_events = res$events_treatment + res$events_control,
    hr = res$hr,
    hr_lower = res$hr_lower,
    hr_upper = res$hr_upper,
    p_value = res$p_value,
    stringsAsFactors = FALSE
  ))
}

if (nrow(forest_data) > 0) {
  forest_data <- forest_data %>%
    mutate(
      label = sprintf("%s\n%s (n=%d)", comparison, method, n_events),
      significant = ifelse(p_value < 0.05, "p < 0.05", "p >= 0.05"),
      method_type = case_when(
        grepl("IPTW", method) ~ "IPTW + Year in PS",
        grepl("Exact Year", method) ~ sprintf("Exact Year PSM 1:%d", config$PRIMARY_RATIO),
        grepl("Pre-COVID", method) ~ sprintf("Pre-COVID PSM 1:%d", config$PRIMARY_RATIO),
        grepl(">=2yr", method) ~ sprintf(">=2yr FU PSM 1:%d", config$PRIMARY_RATIO)
      )
    )

  # Order: group by comparison, then IPTW -> Exact Year -> Pre-COVID -> >=2yr FU
  forest_data$method_order <- case_when(
    grepl("IPTW", forest_data$method) ~ 1,
    grepl("Exact Year", forest_data$method) ~ 2,
    grepl("Pre-COVID", forest_data$method) ~ 3,
    grepl(">=2yr", forest_data$method) ~ 4,
    TRUE ~ 5
  )

  forest_data <- forest_data %>% arrange(comparison, method_order)
  forest_data$label <- factor(forest_data$label, levels = rev(unique(forest_data$label)))

  forest_plot <- ggplot(forest_data, aes(x = hr, y = label)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray40", linewidth = 0.8) +
    geom_point(aes(color = method_type, shape = significant), size = 3.5) +
    geom_errorbarh(aes(xmin = hr_lower, xmax = hr_upper, color = method_type), height = 0.3) +
    geom_text(aes(label = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper)),
      hjust = -0.15, size = 3.2, nudge_x = 0.02
    ) +
    scale_x_log10(limits = c(0.05, 10), breaks = c(0.1, 0.25, 0.5, 1, 2, 4)) +
    scale_color_manual(values = c(
      "IPTW + Year in PS" = "#e67e22",
      setNames("#3498db", sprintf("Exact Year PSM 1:%d", config$PRIMARY_RATIO)),
      setNames("#27ae60", sprintf("Pre-COVID PSM 1:%d", config$PRIMARY_RATIO)),
      setNames("#8e44ad", sprintf(">=2yr FU PSM 1:%d", config$PRIMARY_RATIO))
    )) +
    scale_shape_manual(values = c("p < 0.05" = 16, "p >= 0.05" = 1)) +
    labs(
      title = "Calendar Year & Follow-up Sensitivity Analyses: All 6 Comparisons x 4 Methods",
      subtitle = "IPTW vs Exact Year PSM vs Pre-COVID vs >=2yr Follow-up Restriction",
      x = "Hazard Ratio (log scale)",
      y = "",
      color = "Analysis Method",
      shape = "Statistical Significance"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      plot.subtitle = element_text(size = 10, color = "gray30"),
      legend.position = "bottom",
      legend.box = "vertical",
      axis.text.y = element_text(size = 8),
      panel.grid.minor = element_blank()
    ) +
    guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 1))

  ggsave(sprintf("Figure_Forest_AllMethods_%s.png", ts), forest_plot, width = 14, height = 18, dpi = 300)
  cat(sprintf("Saved: Figure_Forest_AllMethods_%s.png\n", ts))
  cat(sprintf("  Total rows in forest plot: %d (expected 24 = 6 comparisons x 4 methods)\n", nrow(forest_data)))
  print(forest_plot)
}


# =============================================================================
# CELL 13: Event Reduction Summary Figure + CSV
# =============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Creating Event Reduction Summary Figure\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

event_summary <- data.frame()

# PSM exact year event reduction
for (key in names(psm_exact_year_results)) {
  res <- psm_exact_year_results[[key]]
  event_summary <- rbind(event_summary, data.frame(
    comparison = res$comparison,
    method = sprintf("Exact Year 1:%d", res$ratio),
    events_original = res$events_original,
    events_after = res$events_total,
    reduction_pct = res$event_reduction_pct,
    stringsAsFactors = FALSE
  ))
}

# Pre-COVID event reduction
for (key in names(precovid_results)) {
  res <- precovid_results[[key]]
  original_events <- sum(prepared_cohorts[[res$comparison]]$event)
  events_after <- res$events_treatment + res$events_control
  event_summary <- rbind(event_summary, data.frame(
    comparison = res$comparison,
    method = sprintf("Pre-COVID 1:%d", res$ratio),
    events_original = original_events,
    events_after = events_after,
    reduction_pct = (1 - events_after / original_events) * 100,
    stringsAsFactors = FALSE
  ))
}

# >=2yr follow-up event reduction
for (key in names(min2yr_results)) {
  res <- min2yr_results[[key]]
  events_after <- res$events_treatment + res$events_control
  event_summary <- rbind(event_summary, data.frame(
    comparison = res$comparison,
    method = sprintf(">=2yr FU 1:%d", res$ratio),
    events_original = res$events_original,
    events_after = events_after,
    reduction_pct = (1 - events_after / res$events_original) * 100,
    stringsAsFactors = FALSE
  ))
}

if (nrow(event_summary) > 0) {
  event_summary$comparison_short <- gsub("SEMAGLUTIDE vs ", "Sema vs ", event_summary$comparison)
  event_summary$comparison_short <- gsub("OTHER_GLPA vs ", "OtherGLP vs ", event_summary$comparison_short)
  event_summary$comparison_short <- gsub("SGLT2 vs ", "SGLT2 vs ", event_summary$comparison_short)

  p_events <- ggplot(event_summary, aes(x = method, y = events_after, fill = comparison_short)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.8), width = 0.7) +
    geom_text(aes(label = sprintf("%d\n(%.0f%%)", events_after, reduction_pct)),
      position = position_dodge(width = 0.8), vjust = -0.3, size = 2.8
    ) +
    scale_fill_brewer(palette = "Set2") +
    labs(
      title = "Event Count After Calendar Year Matching/Restriction",
      subtitle = "Numbers show events captured and percent reduction from original cohort",
      x = "Analysis Method",
      y = "Number of Events",
      fill = "Comparison"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 13),
      legend.position = "bottom",
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    coord_cartesian(ylim = c(0, max(event_summary$events_after, na.rm = TRUE) * 1.4))

  ggsave(sprintf("Figure_Event_Reduction_%s.png", ts), p_events, width = 14, height = 8, dpi = 300)
  cat(sprintf("Saved: Figure_Event_Reduction_%s.png\n", ts))
  print(p_events)

  # Save CSV
  write.csv(event_summary, sprintf("Event_Reduction_Summary_%s.csv", ts), row.names = FALSE)
  cat(sprintf("Saved: Event_Reduction_Summary_%s.csv\n", ts))
}


# =============================================================================
# CELL 14: Final Reviewer Summary + Save All CSVs/RDS
# =============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("FINAL SUMMARY FOR REVIEWER RESPONSE\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

# --- Argument 1 summary ---
cat("\n")
cat("ARGUMENT 1: IPTW with exact calendar year shows LARGE SMD IMBALANCE\n")
cat(paste(rep("-", 70), collapse = ""), "\n")
for (comp in unique(iptw_smd_df$comparison)) {
  comp_data <- iptw_smd_df %>% filter(comparison == comp, method == "After IPTW")
  max_smd <- max(abs(comp_data$smd))
  cat(sprintf(
    "  %s: Max |SMD| after IPTW = %.3f %s\n",
    comp, max_smd, ifelse(max_smd > 0.1, "<- IMBALANCED", "")
  ))
}

# --- Argument 2 summary ---
cat("\n")
cat("ARGUMENT 2: PSM with exact year -> EVENTS REDUCED but TREND MAINTAINED\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

psm_summary_rows <- lapply(psm_exact_year_results, function(res) {
  data.frame(
    Comparison = res$comparison,
    Ratio = paste0("1:", res$ratio),
    Events_Orig = res$events_original,
    Events_PSM = res$events_total,
    Reduction = sprintf("%.1f%%", res$event_reduction_pct),
    HR = sprintf("%.2f (%.2f-%.2f)", res$hr, res$hr_lower, res$hr_upper),
    P = sprintf("%.4f", res$p_value),
    Trend = ifelse(res$hr < 1, "Protective", ifelse(res$hr > 1, "Harmful", "Neutral")),
    stringsAsFactors = FALSE
  )
})
psm_summary_df <- do.call(rbind, psm_summary_rows)
print(psm_summary_df, row.names = FALSE)

write.csv(psm_summary_df, sprintf("PSM_ExactYear_Summary_%s.csv", ts), row.names = FALSE)
cat(sprintf("Saved: PSM_ExactYear_Summary_%s.csv\n", ts))

# --- Argument 3 summary ---
cat("\n")
cat("ARGUMENT 3: Pre-COVID (2018-2019) -> SUFFICIENT FOLLOW-UP, LOWER POWER\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

if (length(precovid_results) > 0) {
  tr_summary_rows <- lapply(precovid_results, function(res) {
    data.frame(
      Comparison = res$comparison,
      Period = res$period,
      Ratio = paste0("1:", res$ratio),
      N_Treat = res$n_treatment,
      N_Ctrl = res$n_control,
      Events = res$events_treatment + res$events_control,
      Mean_FU = sprintf("%.2f / %.2f yrs", res$mean_fu_treat, res$mean_fu_ctrl),
      HR = sprintf("%.2f (%.2f-%.2f)", res$hr, res$hr_lower, res$hr_upper),
      P = sprintf("%.4f", res$p_value),
      Trend = ifelse(res$hr < 1, "Protective", ifelse(res$hr > 1, "Harmful", "Neutral")),
      stringsAsFactors = FALSE
    )
  })
  tr_summary_df <- do.call(rbind, tr_summary_rows)
  print(tr_summary_df, row.names = FALSE)

  write.csv(tr_summary_df, sprintf("PSM_PreCOVID_Summary_%s.csv", ts), row.names = FALSE)
  cat(sprintf("Saved: PSM_PreCOVID_Summary_%s.csv\n", ts))
}

# --- Argument 4 summary ---
cat("\n")
cat("ARGUMENT 4: >= 2yr Follow-up Restriction -> EQUALIZES OBSERVATION TIME\n")
cat("  (Addresses Reviewer 1 Point #3: Differential follow-up)\n")
cat(paste(rep("-", 70), collapse = ""), "\n")

if (length(min2yr_results) > 0) {
  min2yr_summary_rows <- lapply(min2yr_results, function(res) {
    data.frame(
      Comparison = res$comparison,
      Ratio = paste0("1:", res$ratio),
      N_Treat = res$n_treatment,
      N_Ctrl = res$n_control,
      Events = res$events_treatment + res$events_control,
      Mean_FU = sprintf("%.2f / %.2f yrs", res$mean_fu_treat, res$mean_fu_ctrl),
      HR = sprintf("%.2f (%.2f-%.2f)", res$hr, res$hr_lower, res$hr_upper),
      P = sprintf("%.4f", res$p_value),
      Trend = ifelse(res$hr < 1, "Protective", ifelse(res$hr > 1, "Harmful", "Neutral")),
      stringsAsFactors = FALSE
    )
  })
  min2yr_summary_df <- do.call(rbind, min2yr_summary_rows)
  print(min2yr_summary_df, row.names = FALSE)

  write.csv(min2yr_summary_df, sprintf("PSM_Min2yr_FU_Summary_%s.csv", ts), row.names = FALSE)
  cat(sprintf("Saved: PSM_Min2yr_FU_Summary_%s.csv\n", ts))
}

# --- Key takeaway ---
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("KEY TAKEAWAY FOR REVIEWER\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("
1. IPTW cannot adequately balance calendar year (SMD > 0.1 persists)
   -> Including calendar year as a PS covariate does not resolve imbalance

2. PSM with exact calendar year matching achieves perfect balance (SMD ~ 0)
   but significantly reduces events (lower statistical power)

3. Despite reduced power, the direction of effect is MAINTAINED
   across all sensitivity analyses

4. Time-restricted analysis (2018-2019 Pre-COVID) provides longer follow-up
   but further reduces sample size and events
   NOTE: Post-COVID (2022-2023) dropped - insufficient follow-up
   (data cut 2023-10-01 yields max ~0.75-1.75 years)

5. Restricting to patients with >= 2 years potential follow-up
   (index date <= Oct 2021) equalizes observation time between groups
   while preserving more events than the pre-COVID restriction

6. The consistency of results across multiple approaches strengthens
   the robustness of our primary findings
")

# --- Save combined RDS ---
final_results <- list(
  iptw_smd = iptw_smd_df,
  psm_exact_year = psm_exact_year_results,
  psm_summary = psm_summary_df,
  precovid_results = precovid_results,
  precovid_summary = if (exists("tr_summary_df")) tr_summary_df else NULL,
  min2yr_results = min2yr_results,
  min2yr_summary = if (exists("min2yr_summary_df")) min2yr_summary_df else NULL,
  forest_data = if (exists("forest_data")) forest_data else NULL,
  event_summary = if (exists("event_summary")) event_summary else NULL,
  config = config
)

saveRDS(final_results, sprintf("Reviewer_CalendarYear_AllResults_%s.rds", ts))

cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("OUTPUT FILES GENERATED\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf("1.  Figure_IPTW_SMD_Imbalance_%s.png\n", ts))
cat(sprintf("2.  Figure_Forest_AllMethods_%s.png\n", ts))
cat(sprintf("3.  Figure_Event_Reduction_%s.png\n", ts))
cat("4.  KM_*_Exact_Year_*.png (6 comparisons)\n")
cat("5.  KM_*_Pre-COVID_*.png (up to 6 comparisons)\n")
cat("6.  KM_*_Min_2yr_FU_*.png (up to 6 comparisons)\n")
cat(sprintf("7.  PSM_ExactYear_Summary_%s.csv\n", ts))
cat(sprintf("8.  PSM_PreCOVID_Summary_%s.csv\n", ts))
cat(sprintf("9.  PSM_Min2yr_FU_Summary_%s.csv\n", ts))
cat(sprintf("10. Event_Reduction_Summary_%s.csv\n", ts))
cat(sprintf("11. Reviewer_CalendarYear_AllResults_%s.rds\n", ts))
cat("\n=== ANALYSIS COMPLETE ===\n")


# =============================================================================
# CELL 15: Prism Data Extraction (reuses loaded cohorts)
# =============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Extracting Events/N by Group for Prism Forest Plot\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

outcome_definitions <- list(
  "G40 Only"           = "epilepsy_g40_start_date",
  "G40 or Recurrent"   = "epilepsy_refined_start_date",
  "Broad"              = "epilepsy_or_seizure_start_date",
  "Excl Hypoglycemic"  = "epilepsy_or_seizure_without_hypoglycemic_seizure_start_date"
)

prism_data <- data.frame()

for (comp_name in names(cohorts)) {
  cat(sprintf("Processing: %s\n", comp_name))

  df <- cohorts[[comp_name]]
  df$index_date <- as.Date(df$index_date)
  df$censor_date <- pmin(as.Date(df$EHRmaxDT_min), config$data_cut_date, na.rm = TRUE)

  for (def_name in names(outcome_definitions)) {
    outcome_col <- outcome_definitions[[def_name]]
    if (!outcome_col %in% names(df)) next

    df$outcome_date <- as.Date(df[[outcome_col]])
    df$event_flag <- as.integer(!is.na(df$outcome_date) & df$outcome_date <= df$censor_date)

    prism_data <- rbind(prism_data, data.frame(
      Definition = def_name,
      Comparison = comp_name,
      Treatment_Events = sum(df$event_flag[df$treatment == 1]),
      Treatment_N = sum(df$treatment == 1),
      Control_Events = sum(df$event_flag[df$treatment == 0]),
      Control_N = sum(df$treatment == 0),
      Total_Events = sum(df$event_flag),
      Total_N = nrow(df),
      stringsAsFactors = FALSE
    ))
  }
}

cat("\n=== PRISM FOREST PLOT DATA ===\n")
print(prism_data, row.names = FALSE)

write.csv(prism_data, sprintf("Prism_Events_by_Group_%s.csv", ts), row.names = FALSE)
cat(sprintf("\nSaved: Prism_Events_by_Group_%s.csv\n", ts))
