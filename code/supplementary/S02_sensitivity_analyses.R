# =============================================================================
# rev00007_sensitivity.R.r
# Comprehensive Sensitivity Analyses with Forest Plots
# =============================================================================
#
# THREE ANALYSES, EACH PRODUCING A SUPPLEMENTAL FOREST PLOT:
#
#   1) 4 Outcome Definitions x All Ages / Late Onset (IPTW from RDS)
#   2) 4-Panel Sensitivity: PSM + IPTW excl severe renal & metastatic
#   3) Calendar Year Sensitivity: 4 methods x 6 comparisons
#
# DATA SOURCE: Pre-computed IPTW RDS files containing $cohort (fully enriched/
#              imputed dataframes with event, event_time, treatment, all PS
#              covariates) and $cox (summary(coxph(...))).
#
# All functions inlined (no source() calls).
# =============================================================================


# >>> CELL 01: Load Libraries <<<
suppressPackageStartupMessages({
  library(dplyr)
  library(ggplot2)
  library(survival)
  library(MatchIt)
  library(gridExtra)
  library(cowplot)
  library(lubridate)
})

cat("Libraries loaded successfully.\n")


# >>> CELL 02: Configuration <<<

# ---------- RDS file paths ----------
rds_paths <- list(
  all_ages_4outcomes = "epilepsy_seizure_all_ages_ipwt_full_results.rds",
  late_onset_4outcomes = "epilepsy_seizure_late_onset_ipwt_full_results.rds",
  outcome4_all_ages = "outcome4_epilepsy_excl_hypo_all_ages_ipwt_full_results.rds",
  outcome4_late_onset = "outcome4_epilepsy_excl_hypo_late_onset_ipwt_full_results.rds"
)

# ---------- Outcome definitions ----------
outcome_defs <- list(
  list(var = "epilepsy_or_seizure_start_date",
       label = "Broad (Epilepsy/Seizure)"),
  list(var = "epilepsy_refined_start_date",
       label = "Refined Epilepsy"),
  list(var = "epilepsy_g40_start_date",
       label = "G40 Epilepsy"),
  list(var = "epilepsy_or_seizure_without_hypoglycemic_seizure_start_date",
       label = "Excl Hypoglycemic")
)

# ---------- 6 comparisons ----------
comparisons_6 <- c(
  "SEMAGLUTIDE vs OtherGLD",
  "SEMAGLUTIDE vs SGLT2",
  "SGLT2 vs OtherGLD",
  "OTHER_GLPA vs OtherGLD",
  "OTHER_GLPA vs SGLT2",
  "SEMAGLUTIDE vs OTHER_GLPA"
)

# ---------- 5 comparisons for Analysis 2 (matching image) ----------
comparisons_5 <- c(
  "SEMAGLUTIDE vs OtherGLD",
  "SEMAGLUTIDE vs SGLT2",
  "SGLT2 vs OtherGLD",
  "OTHER_GLPA vs OtherGLD",
  "OTHER_GLPA vs SGLT2"
)

# ---------- Display labels for comparisons ----------
comparison_display <- c(
  "SEMAGLUTIDE vs OTHER_GLPA" = "Semaglutide vs Other GLP-1RAs",
  "SEMAGLUTIDE vs SGLT2"      = "Semaglutide vs SGLT2i",
  "SEMAGLUTIDE vs OtherGLD"   = "Semaglutide vs Other GLDs",
  "OTHER_GLPA vs SGLT2"       = "Other GLP-1RAs vs SGLT2i",
  "OTHER_GLPA vs OtherGLD"    = "Other GLP-1RAs vs Other GLDs",
  "SGLT2 vs OtherGLD"         = "SGLT2i vs Other GLDs"
)

# ---------- PS variables (33 vars) ----------
ps_vars <- c(
  "age", "sex_cat", "raceethnicity_cat", "income", "education", "insurance_category",
  "smoking", "alcohol_category",
  "baseline_bmi_category", "baseline_hba1c", "index_year_grouped",
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
)

# ---------- Comparison-specific PS variable exclusions ----------
comparison_exclusions <- list(
  "SEMAGLUTIDE vs OTHER_GLPA" = c("SEMAGLUTIDE", "OTHER_GLPA"),
  "SEMAGLUTIDE vs SGLT2"      = c("SEMAGLUTIDE", "SGLT2i"),
  "SEMAGLUTIDE vs OtherGLD"   = c("SEMAGLUTIDE", "TZD", "SU", "DPP4i"),
  "OTHER_GLPA vs SGLT2"       = c("OTHER_GLPA", "SGLT2i"),
  "OTHER_GLPA vs OtherGLD"    = c("OTHER_GLPA", "TZD", "SU", "DPP4i"),
  "SGLT2 vs OtherGLD"         = c("SGLT2i", "TZD", "SU", "DPP4i")
)

# ---------- Categorical vars ----------
categorical_vars <- c(
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
)

PRIMARY_RATIO <- 5

cat("Configuration set.\n")


# >>> CELL 03: Helper Functions <<<

# --- SMD calculation (single variable) ---
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

# --- SMD detailed (all variables) ---
calculate_smd_detailed <- function(df, vars, weights = NULL) {
  results_overall <- data.frame(variable = character(), smd = numeric(),
                                 type = character(), stringsAsFactors = FALSE)
  for (var in vars) {
    if (!var %in% names(df)) next
    if (is.factor(df[[var]]) || is.character(df[[var]])) {
      levels_vec <- unique(na.omit(df[[var]]))
      level_smds <- numeric()
      for (lvl in levels_vec) {
        var_binary <- as.numeric(df[[var]] == lvl)
        smd_val <- calc_smd_single(var_binary, df$treatment, weights)
        level_smds <- c(level_smds, smd_val)
      }
      overall_smd <- max(abs(level_smds), na.rm = TRUE)
      results_overall <- rbind(results_overall, data.frame(
        variable = var, smd = overall_smd, type = "categorical",
        stringsAsFactors = FALSE
      ))
    } else {
      smd_val <- calc_smd_single(df[[var]], df$treatment, weights)
      results_overall <- rbind(results_overall, data.frame(
        variable = var, smd = smd_val, type = "continuous",
        stringsAsFactors = FALSE
      ))
    }
  }
  return(list(overall = results_overall))
}

# --- Balance metrics (for IPTW) ---
calculate_balance_metrics <- function(df, vars, weights, verbose = FALSE) {
  balance_df <- data.frame(variable = character(), std_diff = numeric(),
                           stringsAsFactors = FALSE)
  multilevel_overall_vars <- c("sex_cat", "raceethnicity_cat", "income",
                                "education", "insurance_category", "smoking",
                                "alcohol_category", "index_year_grouped",
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
        level_smd <- if (pooled_sd > 0) (mean_t1 - mean_t0) / pooled_sd else 0
        level_smds <- c(level_smds, level_smd)
      }
      overall_smd <- ifelse(length(level_smds) > 0,
                             max(abs(level_smds), na.rm = TRUE), 0)
      if (length(level_smds) > 0) {
        max_idx <- which.max(abs(level_smds))
        overall_smd <- overall_smd * sign(level_smds[max_idx])
      }
      balance_df <- rbind(balance_df, data.frame(variable = var,
                                                  std_diff = overall_smd))
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
        std_diff <- if (pooled_sd > 0) (mean_t1 - mean_t0) / pooled_sd else 0
        balance_df <- rbind(balance_df, data.frame(variable = paste0(var, lvl),
                                                    std_diff = std_diff))
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
      std_diff <- if (pooled_sd > 0) (mean_t1 - mean_t0) / pooled_sd else 0
      balance_df <- rbind(balance_df, data.frame(variable = var,
                                                  std_diff = std_diff))
    }
  }
  return(balance_df)
}

# --- IPTW and Cox ---
run_iptw_and_cox <- function(df, exclude_vars = NULL, trim_threshold = 0.01,
                             all_ps_vars = NULL, verbose = FALSE,
                             stabilize = FALSE) {
  all_covariates <- if (is.null(all_ps_vars)) ps_vars else all_ps_vars
  rhs_vars <- setdiff(all_covariates, exclude_vars)
  keep_vars <- rhs_vars[sapply(df[rhs_vars], function(x)
    length(unique(na.omit(x))) > 1)]
  ps_form <- reformulate(keep_vars, response = "treatment")
  model_vars <- c("treatment", keep_vars)

  complete_rows <- complete.cases(df[model_vars])
  df_complete <- df[complete_rows, ]

  if (isTRUE(verbose)) {
    cat(sprintf("  IPTW: %d complete cases from %d total (%d removed)\n",
                nrow(df_complete), nrow(df), nrow(df) - nrow(df_complete)))
  }

  ps_model <- glm(ps_form, data = df_complete, family = binomial())
  df_complete$ps <- predict(ps_model, type = "response")

  if (isTRUE(stabilize)) {
    p_treat <- mean(df_complete$treatment, na.rm = TRUE)
    df_complete$ipw <- ifelse(df_complete$treatment == 1,
                              p_treat / df_complete$ps,
                              (1 - p_treat) / (1 - df_complete$ps))
  } else {
    df_complete$ipw <- ifelse(df_complete$treatment == 1,
                              1 / df_complete$ps,
                              1 / (1 - df_complete$ps))
  }

  if (!is.null(trim_threshold)) {
    lower <- quantile(df_complete$ipw, trim_threshold)
    upper <- quantile(df_complete$ipw, 1 - trim_threshold)
    df_complete$ipw_trimmed <- pmin(pmax(df_complete$ipw, lower), upper)
  } else {
    df_complete$ipw_trimmed <- df_complete$ipw
  }

  df_complete$ipw_std <- df_complete$ipw_trimmed * nrow(df_complete) /
    sum(df_complete$ipw_trimmed)

  balance_before <- calculate_balance_metrics(df_complete, keep_vars, NULL,
                                              verbose = FALSE)
  balance_after <- calculate_balance_metrics(df_complete, keep_vars,
                                             df_complete$ipw_std,
                                             verbose = FALSE)

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

# --- PSM with exact year matching + cluster(subclass) Cox ---
run_psm_exact_year <- function(df, exclude_vars, ps_vars_arg,
                               caliper = 0.2, ratio = 1) {
  rhs_vars <- setdiff(ps_vars_arg,
                       c(exclude_vars, "index_year", "index_year_grouped"))
  keep_vars <- rhs_vars[sapply(df[rhs_vars], function(x)
    length(unique(na.omit(x))) > 1)]

  if (!"index_year" %in% names(df)) {
    df$index_year <- factor(lubridate::year(df$index_date))
  }

  model_vars <- c("treatment", keep_vars, "index_year")
  complete_rows <- complete.cases(df[model_vars])
  df_complete <- df[complete_rows, ]

  cat(sprintf("  Complete cases: %d (Treatment: %d, Control: %d)\n",
              nrow(df_complete), sum(df_complete$treatment == 1),
              sum(df_complete$treatment == 0)))

  ps_formula <- reformulate(keep_vars, response = "treatment")

  tryCatch({
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

    cat(sprintf("  Matched: %d total (Treat: %d, Ctrl: %d) | Events: %d (Treat: %d, Ctrl: %d)\n",
                nrow(matched_df), n_treat, n_ctrl,
                events_treat + events_ctrl, events_treat, events_ctrl))

    cox_model <- coxph(Surv(event_time, event) ~ treatment + cluster(subclass),
                       data = matched_df)

    return(list(
      cohort_matched = matched_df,
      cox = summary(cox_model),
      ratio = ratio
    ))
  }, error = function(e) {
    cat(sprintf("  PSM 1:%d matching failed: %s\n", ratio, e$message))
    return(NULL)
  })
}

# --- PSM matching original approach (no exact year, logit-PS distance) ---
# Replicates run_psm_and_cox() from 111225_ITT-sensitivity_analysis.r:
#   - index_year_grouped INCLUDED in PS covariates (not excluded)
#   - Pre-computed logit-PS passed as distance
#   - std.caliper = TRUE, replace = FALSE
#   - NO exact matching on index_year
run_psm_original <- function(df, exclude_vars, ps_vars_arg,
                             caliper = 0.2, ratio = 1) {
  rhs_vars <- setdiff(ps_vars_arg, exclude_vars)
  keep_vars <- rhs_vars[sapply(df[rhs_vars], function(x)
    length(unique(na.omit(x))) > 1)]

  model_vars <- c("treatment", keep_vars)
  complete_rows <- complete.cases(df[model_vars])
  df_complete <- df[complete_rows, ]

  cat(sprintf("  Complete cases: %d (Treatment: %d, Control: %d)\n",
              nrow(df_complete), sum(df_complete$treatment == 1),
              sum(df_complete$treatment == 0)))

  if (nrow(df_complete) < 10 ||
      sum(df_complete$treatment == 1) < 2 ||
      sum(df_complete$treatment == 0) < 2) {
    cat("  Insufficient data for matching.\n")
    return(NULL)
  }

  ps_formula <- reformulate(keep_vars, response = "treatment")

  tryCatch({
    # Fit PS model and compute logit-PS (matching original approach)
    ps_fit <- glm(ps_formula, data = df_complete, family = binomial())
    df_complete$ps <- predict(ps_fit, type = "response")
    df_complete$ps_logit <- predict(ps_fit, type = "link")

    m_out <- matchit(
      ps_formula,
      data = df_complete,
      method = "nearest",
      ratio = ratio,
      distance = df_complete$ps_logit,
      caliper = caliper,
      std.caliper = TRUE,
      replace = FALSE
    )

    matched_df <- match.data(m_out)

    n_treat <- sum(matched_df$treatment == 1)
    n_ctrl <- sum(matched_df$treatment == 0)
    events_treat <- sum(matched_df$event[matched_df$treatment == 1])
    events_ctrl <- sum(matched_df$event[matched_df$treatment == 0])

    cat(sprintf("  Matched: %d total (Treat: %d, Ctrl: %d) | Events: %d (Treat: %d, Ctrl: %d)\n",
                nrow(matched_df), n_treat, n_ctrl,
                events_treat + events_ctrl, events_treat, events_ctrl))

    cox_model <- coxph(Surv(event_time, event) ~ treatment + cluster(subclass),
                       data = matched_df)

    return(list(
      cohort_matched = matched_df,
      cox = summary(cox_model),
      ratio = ratio
    ))
  }, error = function(e) {
    cat(sprintf("  PSM 1:%d matching failed: %s\n", ratio, e$message))
    return(NULL)
  })
}

# --- PSM without exact year (for time-restricted) + cluster(subclass) Cox ---
run_psm_no_exact_year <- function(df, exclude_vars, ps_vars_arg,
                                  caliper = 0.2, ratio = 1) {
  rhs_vars <- setdiff(ps_vars_arg,
                       c(exclude_vars, "index_year", "index_year_grouped"))
  keep_vars <- rhs_vars[sapply(df[rhs_vars], function(x)
    length(unique(na.omit(x))) > 1)]

  model_vars <- c("treatment", keep_vars)
  complete_rows <- complete.cases(df[model_vars])
  df_complete <- df[complete_rows, ]

  if (nrow(df_complete) < 50 ||
      sum(df_complete$treatment == 1) < 10 ||
      sum(df_complete$treatment == 0) < 10) {
    cat("  Insufficient sample size for PSM.\n")
    return(NULL)
  }

  ps_formula <- reformulate(keep_vars, response = "treatment")

  tryCatch({
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
    if (nrow(matched_df) < 20) return(NULL)

    n_treat <- sum(matched_df$treatment == 1)
    n_ctrl <- sum(matched_df$treatment == 0)
    events_treat <- sum(matched_df$event[matched_df$treatment == 1])
    events_ctrl <- sum(matched_df$event[matched_df$treatment == 0])

    cat(sprintf("    Matched: N = %d (Treat: %d, Ctrl: %d) | Events: %d (Treat: %d, Ctrl: %d)\n",
                nrow(matched_df), n_treat, n_ctrl,
                events_treat + events_ctrl, events_treat, events_ctrl))

    cox_model <- coxph(Surv(event_time, event) ~ treatment + cluster(subclass),
                       data = matched_df)

    return(list(
      cohort_matched = matched_df,
      cox = summary(cox_model),
      ratio = ratio
    ))
  }, error = function(e) {
    cat(sprintf("    PSM 1:%d failed: %s\n", ratio, e$message))
    return(NULL)
  })
}

# --- Ensure composite variables exist in cohort ---
ensure_composite_vars <- function(df) {
  # renal_disease_severity
  if (!"renal_disease_severity" %in% names(df)) {
    if (all(c("renal_severe", "renal_mild_or_moderate") %in% names(df))) {
      df$renal_disease_severity <- factor(case_when(
        df$renal_severe == 1 | df$renal_severe == "Yes" ~ "Severe",
        df$renal_mild_or_moderate == 1 |
          df$renal_mild_or_moderate == "Yes" ~ "Mild-Moderate",
        TRUE ~ "None"
      ), levels = c("None", "Mild-Moderate", "Severe"))
    }
  }

  # liver_disease_severity
  if (!"liver_disease_severity" %in% names(df)) {
    if (all(c("moderate_or_severe_liver_disease", "mild_liver_disease") %in%
            names(df))) {
      df$liver_disease_severity <- factor(case_when(
        df$moderate_or_severe_liver_disease == 1 |
          df$moderate_or_severe_liver_disease == "Yes" ~ "Mod-Severe",
        df$mild_liver_disease == 1 |
          df$mild_liver_disease == "Yes" ~ "Mild",
        TRUE ~ "None"
      ), levels = c("None", "Mild", "Mod-Severe"))
    }
  }

  # malignancy_status
  if (!"malignancy_status" %in% names(df)) {
    if (all(c("metastatic_solid_tumor", "any_malignancy") %in% names(df))) {
      df$malignancy_status <- factor(case_when(
        df$metastatic_solid_tumor == 1 |
          df$metastatic_solid_tumor == "Yes" ~ "Metastatic",
        df$any_malignancy == 1 |
          df$any_malignancy == "Yes" ~ "Non-Metastatic",
        TRUE ~ "None"
      ), levels = c("None", "Non-Metastatic", "Metastatic"))
    }
  }

  # index_year
  if (!"index_year" %in% names(df) && "index_date" %in% names(df)) {
    df$index_year <- factor(lubridate::year(as.Date(df$index_date)))
  }

  # index_year_grouped
  if (!"index_year_grouped" %in% names(df) && "index_year" %in% names(df)) {
    yr_num <- as.numeric(as.character(df$index_year))
    df$index_year_grouped <- factor(case_when(
      yr_num %in% c(2020, 2021, 2022) ~ "COVID",
      TRUE ~ "Non-COVID"
    ), levels = c("Non-COVID", "COVID"))
  }

  # Ensure categorical vars are factors
  for (var in categorical_vars) {
    if (var %in% names(df) && !is.factor(df[[var]])) {
      df[[var]] <- as.factor(df[[var]])
    }
  }

  return(df)
}

# --- Extract HR from stored cox summary ---
extract_hr <- function(cox_summary) {
  hr <- cox_summary$conf.int[1, "exp(coef)"]
  hr_lo <- cox_summary$conf.int[1, "lower .95"]
  hr_hi <- cox_summary$conf.int[1, "upper .95"]
  p_val <- cox_summary$coefficients[1, "Pr(>|z|)"]
  list(hr = hr, hr_lower = hr_lo, hr_upper = hr_hi, p_value = p_val)
}

# --- Parse RDS key into comparison + outcome ---
parse_rds_key <- function(key) {
  # Keys formatted as "COMPARISON - OUTCOME_LABEL"
  parts <- strsplit(key, " - ", fixed = TRUE)[[1]]
  if (length(parts) >= 2) {
    comp <- parts[1]
    outcome <- paste(parts[-1], collapse = " - ")
  } else {
    comp <- key
    outcome <- "Unknown"
  }
  list(comparison = comp, outcome = outcome)
}

cat("All helper functions defined.\n")


# =============================================================================
# >>> CELL 04: Analysis 1 - SKIPPED (only running Analysis 2) <<<
# =============================================================================
if (FALSE) {
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ANALYSIS 1: IPTW Results Across 4 Outcome Definitions\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

forest_data_1 <- data.frame()

for (age_group in c("All Ages", "Late Onset")) {
  rds_file <- if (age_group == "All Ages") {
    rds_paths$all_ages_4outcomes
  } else {
    rds_paths$late_onset_4outcomes
  }

  cat(sprintf("\nLoading %s: %s\n", age_group, rds_file))

  rds_obj <- tryCatch({
    readRDS(rds_file)
  }, error = function(e) {
    cat(sprintf("  Failed to load: %s\n", e$message))
    NULL
  })

  if (is.null(rds_obj)) next

  cat(sprintf("  Found %d analyses\n", length(rds_obj)))

  # Reorder keys: SEMA vs OtherGLD first, SEMA vs SGLT2 second, then rest
  priority_comps <- c("SEMAGLUTIDE vs OtherGLD", "SEMAGLUTIDE vs SGLT2",
                      "SGLT2 vs OtherGLD", "OTHER_GLPA vs OtherGLD",
                      "OTHER_GLPA vs SGLT2", "SEMAGLUTIDE vs OTHER_GLPA")
  all_keys <- names(rds_obj)
  ordered_keys <- character(0)
  for (pc in priority_comps) {
    ordered_keys <- c(ordered_keys, all_keys[grepl(pc, all_keys, fixed = TRUE)])
  }
  # Add any remaining keys not matched
  ordered_keys <- c(ordered_keys, setdiff(all_keys, ordered_keys))

  for (key in ordered_keys) {
    res <- rds_obj[[key]]
    parsed <- parse_rds_key(key)
    comp_name <- parsed$comparison
    outcome_label <- parsed$outcome

    # Extract HR from stored cox summary
    hr_info <- tryCatch({
      extract_hr(res$cox)
    }, error = function(e) {
      cat(sprintf("  Could not extract HR from '%s': %s\n", key, e$message))
      NULL
    })

    if (is.null(hr_info)) next

    # Extract event counts from cohort
    cohort <- res$cohort
    n_total <- nrow(cohort)
    n_events <- sum(cohort$event)
    events_treat <- sum(cohort$event[cohort$treatment == 1])
    events_ctrl <- sum(cohort$event[cohort$treatment == 0])
    n_treat <- sum(cohort$treatment == 1)
    n_ctrl <- sum(cohort$treatment == 0)

    display_comp <- ifelse(comp_name %in% names(comparison_display),
                           comparison_display[comp_name], comp_name)

    forest_data_1 <- rbind(forest_data_1, data.frame(
      comparison = comp_name,
      comparison_display = display_comp,
      outcome = outcome_label,
      age_group = age_group,
      hr = hr_info$hr,
      hr_lower = hr_info$hr_lower,
      hr_upper = hr_info$hr_upper,
      p_value = hr_info$p_value,
      n_total = n_total,
      n_events = n_events,
      events_treat = events_treat,
      events_ctrl = events_ctrl,
      n_treat = n_treat,
      n_ctrl = n_ctrl,
      stringsAsFactors = FALSE
    ))

    cat(sprintf("  %s | HR: %.3f (%.3f-%.3f), p=%.4f, treat=%d/%d, ctrl=%d/%d\n",
                key, hr_info$hr, hr_info$hr_lower, hr_info$hr_upper,
                hr_info$p_value, events_treat, n_treat, events_ctrl, n_ctrl))
  }
}

cat(sprintf("\nForest data 1: %d rows\n", nrow(forest_data_1)))

# --- Forest Plot 1: Stacked Panels (Supplementary Figure Style) ---
if (nrow(forest_data_1) > 0) {
  forest_data_1 <- forest_data_1 %>%
    mutate(
      treat_name = sub(" vs .*", "", comparison_display),
      ctrl_name = sub(".* vs ", "", comparison_display),
      arm_label = sprintf("%s (%d/%d) vs %s (%d/%d)",
                          treat_name, events_treat, n_treat,
                          ctrl_name, events_ctrl, n_ctrl),
      hr_text = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper),
      p_text = ifelse(p_value < 0.001, "< 0.001", sprintf("%.3f", p_value))
    )

  # Shape map per comparison (6 comparisons)
  shape_map_1 <- c(
    "SEMAGLUTIDE vs OTHER_GLPA" = 16,  # filled circle
    "SEMAGLUTIDE vs SGLT2"      = 15,  # filled square
    "SEMAGLUTIDE vs OtherGLD"   = 17,  # filled triangle
    "OTHER_GLPA vs SGLT2"       = 18,  # filled diamond
    "OTHER_GLPA vs OtherGLD"    = 8,   # asterisk
    "SGLT2 vs OtherGLD"         = 4    # cross
  )

  # Comparison ordering
  comp_order_1 <- c(
    "SEMAGLUTIDE vs OtherGLD", "SEMAGLUTIDE vs SGLT2",
    "SGLT2 vs OtherGLD", "OTHER_GLPA vs OtherGLD",
    "OTHER_GLPA vs SGLT2", "SEMAGLUTIDE vs OTHER_GLPA"
  )

  # Panels: age_group x outcome = up to 8 panels
  # Group by outcome definition, each panel showing all comparisons
  outcome_labels <- c(
    "Broad (Epilepsy/Seizure)", "Refined Epilepsy",
    "G40 Epilepsy", "Excl Hypoglycemic"
  )

  panel_counter <- 0
  panel_letters <- LETTERS
  fp1_panels <- list()

  for (ag in c("All Ages", "Late Onset")) {
    for (oc in outcome_labels) {
      panel_counter <- panel_counter + 1
      pl <- panel_letters[panel_counter]

      pdf <- forest_data_1 %>%
        filter(age_group == ag, outcome == oc) %>%
        mutate(comparison = factor(comparison, levels = comp_order_1)) %>%
        arrange(comparison)

      if (nrow(pdf) == 0) next

      pdf$arm_label <- factor(pdf$arm_label, levels = rev(pdf$arm_label))

      x_max <- 2.0
      is_bottom <- (ag == "Late Onset" && oc == tail(outcome_labels, 1))

      p <- ggplot(pdf, aes(x = hr, y = arm_label)) +
        geom_vline(xintercept = 1, linetype = "dashed", color = "black",
                   linewidth = 0.5) +
        geom_errorbarh(aes(xmin = hr_lower, xmax = pmin(hr_upper, x_max)),
                       height = 0.3, linewidth = 0.5) +
        geom_point(aes(shape = comparison), size = 3, fill = "black") +
        scale_shape_manual(values = shape_map_1, guide = "none") +
        geom_text(aes(x = x_max + 0.15, label = hr_text),
                  hjust = 0, size = 3, fontface = "plain") +
        geom_text(aes(x = x_max + 0.85, label = p_text),
                  hjust = 0, size = 3, fontface = "plain") +
        scale_x_continuous(
          limits = c(0, x_max), breaks = seq(0, 2.0, 0.5),
          expand = expansion(mult = c(0.02, 0))
        ) +
        coord_cartesian(clip = "off") +
        labs(
          title = sprintf("(%s) IPTW — %s [%s]", pl, oc, ag),
          x = if (is_bottom) "Hazard Ratio" else NULL,
          y = NULL
        ) +
        theme_bw(base_size = 10) +
        theme(
          plot.title = element_text(face = "bold", size = 9, hjust = 0,
                                    margin = margin(b = 4)),
          plot.title.position = "plot",
          axis.text.y = element_text(size = 8, color = "black"),
          axis.text.x = element_text(size = 8),
          axis.title.x = if (is_bottom) element_text(size = 9) else element_blank(),
          panel.grid.major.y = element_blank(),
          panel.grid.minor = element_blank(),
          plot.margin = margin(t = 8, r = 120, b = 4, l = 5)
        )

      fp1_panels[[pl]] <- p
    }
  }

  if (length(fp1_panels) > 0) {
    header1 <- ggdraw() +
      draw_label("Comparison (No. of patients with event/total)",
                 x = 0.01, hjust = 0, fontface = "bold", size = 9) +
      draw_label("HR (95% CI)",
                 x = 0.72, hjust = 0, fontface = "bold", size = 9) +
      draw_label("P Value",
                 x = 0.88, hjust = 0, fontface = "bold", size = 9)

    fp1 <- plot_grid(
      header1,
      plotlist = fp1_panels,
      ncol = 1,
      rel_heights = c(0.03, rep(1, length(fp1_panels))),
      align = "v", axis = "lr"
    )

    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    ggsave(sprintf("sensitivity_ForestPlot1_Outcomes_%s.png", ts), fp1,
           width = 12, height = 3.5 * length(fp1_panels), dpi = 300)
    cat(sprintf("Saved: sensitivity_ForestPlot1_Outcomes_%s.png\n", ts))
    print(fp1)
  }
}
} # end if (FALSE) — skip Analysis 1


# =============================================================================
# >>> CELL 05: Analysis 2 - 4-Panel Sensitivity (PSM + IPTW excl) <<<
# =============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ANALYSIS 2: PSM + IPTW Excluding Severe Renal & Metastatic (4 Panels)\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

forest_data_2 <- data.frame()

# Load both RDS files
rds_all <- tryCatch(readRDS(rds_paths$outcome4_all_ages), error = function(e) {
  cat(sprintf("Failed to load all_ages RDS: %s\n", e$message)); NULL
})
rds_late <- tryCatch(readRDS(rds_paths$outcome4_late_onset), error = function(e) {
  cat(sprintf("Failed to load late_onset RDS: %s\n", e$message)); NULL
})

# Function to find cohort in RDS for a given comparison
find_cohort <- function(rds_obj, comp_name) {
  if (is.null(rds_obj)) return(NULL)
  # Try exact key match patterns
  for (key in names(rds_obj)) {
    if (grepl(comp_name, key, fixed = TRUE)) {
      cohort <- rds_obj[[key]]$cohort
      if (!is.null(cohort) && nrow(cohort) > 0) {
        return(ensure_composite_vars(cohort))
      }
    }
  }
  return(NULL)
}

# --- Panel definitions ---
# Panel A/C: PSM 1:1 matching (matching figure with equal N per arm)
# Panel B/D: IPTW on cohort excluding severe renal + metastatic
PSM_RATIO_ANALYSIS2 <- 1

panels <- list(
  list(
    label = "A",
    title = "PSM analysis with the original cohort (Adult-Onset)",
    rds_source = "all_ages",
    method = "PSM",
    exclusion = FALSE
  ),
  list(
    label = "B",
    title = "IPTW excl severe renal & metastatic (Adult-Onset)",
    rds_source = "all_ages",
    method = "IPTW",
    exclusion = TRUE
  ),
  list(
    label = "C",
    title = "PSM analysis with the original cohort (Late-Onset)",
    rds_source = "late_onset",
    method = "PSM",
    exclusion = FALSE
  ),
  list(
    label = "D",
    title = "IPTW excl severe renal & metastatic (Late-Onset)",
    rds_source = "late_onset",
    method = "IPTW",
    exclusion = TRUE
  )
)

for (panel in panels) {
  cat(sprintf("\n--- Panel %s: %s ---\n", panel$label, panel$title))

  rds_source <- if (panel$rds_source == "all_ages") rds_all else rds_late

  for (comp_name in comparisons_5) {
    cat(sprintf("\n  %s\n", comp_name))

    cohort <- find_cohort(rds_source, comp_name)
    if (is.null(cohort)) {
      cat("    Cohort not found in RDS. Skipping.\n")
      next
    }

    exclude_vars <- comparison_exclusions[[comp_name]]
    display_comp <- comparison_display[comp_name]

    if (panel$exclusion) {
      # Exclude severe renal AND metastatic
      n_before <- nrow(cohort)
      cohort <- cohort %>%
        filter(renal_disease_severity != "Severe" &
               malignancy_status != "Metastatic")
      n_excluded <- n_before - nrow(cohort)
      cat(sprintf("    Excluded %d patients (severe renal + metastatic). N: %d -> %d\n",
                  n_excluded, n_before, nrow(cohort)))

      if (nrow(cohort) < 50 ||
          sum(cohort$treatment == 1) < 10 ||
          sum(cohort$treatment == 0) < 10) {
        cat("    Insufficient sample size. Skipping.\n")
        next
      }
    }

    if (panel$method == "PSM") {
      # PSM 1:1 using original approach (no exact year, logit-PS distance)
      result <- tryCatch({
        run_psm_original(
          df = cohort,
          exclude_vars = exclude_vars,
          ps_vars_arg = ps_vars,
          caliper = 0.2,
          ratio = PSM_RATIO_ANALYSIS2
        )
      }, error = function(e) {
        cat(sprintf("    PSM failed: %s\n", e$message))
        NULL
      })

      if (!is.null(result)) {
        hr_info <- extract_hr(result$cox)
        matched_df <- result$cohort_matched
        events_treat <- sum(matched_df$event[matched_df$treatment == 1])
        events_ctrl <- sum(matched_df$event[matched_df$treatment == 0])
        n_treat <- sum(matched_df$treatment == 1)
        n_ctrl <- sum(matched_df$treatment == 0)

        forest_data_2 <- rbind(forest_data_2, data.frame(
          panel = panel$label,
          panel_title = panel$title,
          comparison = comp_name,
          comparison_display = display_comp,
          hr = hr_info$hr,
          hr_lower = hr_info$hr_lower,
          hr_upper = hr_info$hr_upper,
          p_value = hr_info$p_value,
          events_treat = events_treat,
          events_ctrl = events_ctrl,
          n_treat = n_treat,
          n_ctrl = n_ctrl,
          stringsAsFactors = FALSE
        ))

        cat(sprintf("    HR: %.3f (%.3f-%.3f), p=%.4f | Events: treat=%d/%d, ctrl=%d/%d\n",
                    hr_info$hr, hr_info$hr_lower, hr_info$hr_upper,
                    hr_info$p_value, events_treat, n_treat, events_ctrl, n_ctrl))
      }
    } else {
      # IPTW
      result <- tryCatch({
        run_iptw_and_cox(
          df = cohort,
          exclude_vars = exclude_vars,
          trim_threshold = 0.01,
          all_ps_vars = ps_vars,
          verbose = TRUE,
          stabilize = FALSE
        )
      }, error = function(e) {
        cat(sprintf("    IPTW failed: %s\n", e$message))
        NULL
      })

      if (!is.null(result)) {
        hr_info <- extract_hr(result$cox)
        iptw_cohort <- result$cohort
        events_treat <- sum(iptw_cohort$event[iptw_cohort$treatment == 1])
        events_ctrl <- sum(iptw_cohort$event[iptw_cohort$treatment == 0])
        n_treat <- sum(iptw_cohort$treatment == 1)
        n_ctrl <- sum(iptw_cohort$treatment == 0)

        forest_data_2 <- rbind(forest_data_2, data.frame(
          panel = panel$label,
          panel_title = panel$title,
          comparison = comp_name,
          comparison_display = display_comp,
          hr = hr_info$hr,
          hr_lower = hr_info$hr_lower,
          hr_upper = hr_info$hr_upper,
          p_value = hr_info$p_value,
          events_treat = events_treat,
          events_ctrl = events_ctrl,
          n_treat = n_treat,
          n_ctrl = n_ctrl,
          stringsAsFactors = FALSE
        ))

        cat(sprintf("    HR: %.3f (%.3f-%.3f), p=%.4f | Events: treat=%d/%d, ctrl=%d/%d\n",
                    hr_info$hr, hr_info$hr_lower, hr_info$hr_upper,
                    hr_info$p_value, events_treat, n_treat, events_ctrl, n_ctrl))
      }
    }
  }
}

cat(sprintf("\nForest data 2: %d rows\n", nrow(forest_data_2)))

# --- Forest Plot 2: 4-Panel Stacked (Supplementary Figure 4 Style) ---
if (nrow(forest_data_2) > 0) {

  # Split display labels into treatment vs comparator names
  # e.g., "Semaglutide vs Other GLDs" -> treat_name="Semaglutide", ctrl_name="Other GLDs"
  forest_data_2 <- forest_data_2 %>%
    mutate(
      treat_name = sub(" vs .*", "", comparison_display),
      ctrl_name = sub(".* vs ", "", comparison_display),
      # Label: "Semaglutide (15/2369) vs Other GLDs (52/2369)"
      arm_label = sprintf("%s (%d/%d) vs %s (%d/%d)",
                          treat_name, events_treat, n_treat,
                          ctrl_name, events_ctrl, n_ctrl),
      hr_text = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper),
      p_text = ifelse(p_value < 0.001, "< 0.001",
                      sprintf("%.3f", p_value))
    )

  # Assign distinct shapes per comparison (matching figure style)
  shape_map <- c(
    "SEMAGLUTIDE vs OtherGLD"  = 16,  # filled circle
    "SEMAGLUTIDE vs SGLT2"     = 15,  # filled square
    "SGLT2 vs OtherGLD"        = 17,  # filled triangle
    "OTHER_GLPA vs OtherGLD"   = 18,  # filled diamond
    "OTHER_GLPA vs SGLT2"      = 8    # asterisk/star
  )

  # Panel titles matching figure exactly
  panel_titles <- c(
    "A" = "(A) PSM analysis with the original cohort for adult-onset seizure disorder",
    "B" = "(B) IPTW analysis with a cohort excluding severe renal disease and metastatic cancer for adult-onset seizure disorder",
    "C" = "(C) PSM analysis with the original cohort for late-onset seizure disorder",
    "D" = "(D) IPTW analysis with a cohort excluding severe renal disease and metastatic cancer for late-onset seizure disorder"
  )

  # Comparison ordering (same as figure, top to bottom)
  comp_order <- comparisons_5

  make_panel_plot <- function(panel_df, panel_lbl, is_bottom = FALSE) {
    # Order comparisons
    panel_df <- panel_df %>%
      mutate(comparison = factor(comparison, levels = comp_order)) %>%
      arrange(comparison)
    panel_df$arm_label <- factor(panel_df$arm_label,
                                  levels = rev(panel_df$arm_label))

    # Map shapes
    panel_df$pt_shape <- shape_map[as.character(panel_df$comparison)]

    # Fixed x-axis range
    x_max <- 2.0

    p <- ggplot(panel_df, aes(x = hr, y = arm_label)) +
      # Reference line at HR=1
      geom_vline(xintercept = 1, linetype = "dashed", color = "black",
                 linewidth = 0.5) +
      # CI error bars
      geom_errorbarh(aes(xmin = hr_lower, xmax = pmin(hr_upper, x_max)),
                     height = 0.3, linewidth = 0.5) +
      # Points with distinct shapes
      geom_point(aes(shape = comparison), size = 3, fill = "black") +
      scale_shape_manual(values = shape_map, guide = "none") +
      # HR text annotation (right of plot area)
      geom_text(aes(x = x_max + 0.15, label = hr_text),
                hjust = 0, size = 3, fontface = "plain") +
      # P value annotation (further right)
      geom_text(aes(x = x_max + 0.85, label = p_text),
                hjust = 0, size = 3, fontface = "plain") +
      # X-axis
      scale_x_continuous(
        limits = c(0, x_max),
        breaks = seq(0, 2.0, 0.5),
        expand = expansion(mult = c(0.02, 0))
      ) +
      coord_cartesian(clip = "off") +
      labs(
        title = panel_lbl,
        x = if (is_bottom) "Hazard Ratio" else NULL,
        y = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 9, hjust = 0,
                                  margin = margin(b = 4)),
        plot.title.position = "plot",
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8),
        axis.title.x = if (is_bottom) element_text(size = 9) else element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 8, r = 120, b = 4, l = 5)
      )

    # Add column headers for HR and P Value (only on first panel)
    return(p)
  }

  panel_plots <- list()
  panel_order <- c("A", "B", "C", "D")
  for (i in seq_along(panel_order)) {
    pl <- panel_order[i]
    pdf <- forest_data_2 %>% filter(panel == pl)
    if (nrow(pdf) > 0) {
      is_bottom <- (i == length(panel_order))
      panel_plots[[pl]] <- make_panel_plot(pdf, panel_titles[pl], is_bottom)
    }
  }

  if (length(panel_plots) > 0) {
    # Add column headers as a top annotation
    header <- ggdraw() +
      draw_label("Comparison (No. of patients with event/total)",
                 x = 0.01, hjust = 0, fontface = "bold", size = 9) +
      draw_label("HR (95% CI)",
                 x = 0.72, hjust = 0, fontface = "bold", size = 9) +
      draw_label("P Value",
                 x = 0.88, hjust = 0, fontface = "bold", size = 9)

    # Stack panels vertically (ncol = 1)
    fp2 <- plot_grid(
      header,
      plotlist = panel_plots,
      ncol = 1,
      rel_heights = c(0.04, rep(1, length(panel_plots))),
      align = "v",
      axis = "lr"
    )

    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    ggsave(sprintf("sensitivity_ForestPlot2_4Panel_%s.png", ts), fp2,
           width = 12, height = 14, dpi = 300)
    cat(sprintf("Saved: sensitivity_ForestPlot2_4Panel_%s.png\n", ts))
    print(fp2)
  }
}


# =============================================================================
# >>> CELL 06: Analysis 3 - SKIPPED (only running Analysis 2) <<<
# =============================================================================
if (FALSE) {
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("ANALYSIS 3: Calendar Year Sensitivity (4 Methods x 6 Comparisons)\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

forest_data_3 <- data.frame()

# Load cohorts from outcome4 all_ages RDS
if (is.null(rds_all)) {
  rds_all <- tryCatch(readRDS(rds_paths$outcome4_all_ages), error = function(e) {
    cat(sprintf("Failed to load: %s\n", e$message)); NULL
  })
}

for (comp_name in comparisons_6) {
  cat(sprintf("\n=== %s ===\n", comp_name))

  cohort <- find_cohort(rds_all, comp_name)
  if (is.null(cohort)) {
    cat("  Cohort not found. Skipping.\n")
    next
  }

  exclude_vars <- comparison_exclusions[[comp_name]]
  display_comp <- comparison_display[comp_name]

  # --- Method 1: PSM Exact Year 1:1 ---
  cat("  Method 1: PSM Exact Year 1:1\n")
  result_1 <- tryCatch({
    run_psm_exact_year(cohort, exclude_vars, ps_vars, caliper = 0.2, ratio = 1)
  }, error = function(e) {
    cat(sprintf("    Failed: %s\n", e$message)); NULL
  })

  if (!is.null(result_1)) {
    hr_info <- extract_hr(result_1$cox)
    md <- result_1$cohort_matched
    forest_data_3 <- rbind(forest_data_3, data.frame(
      comparison = comp_name,
      comparison_display = display_comp,
      method = "PSM Exact Year 1:1",
      hr = hr_info$hr, hr_lower = hr_info$hr_lower,
      hr_upper = hr_info$hr_upper, p_value = hr_info$p_value,
      n_matched = nrow(md), n_events = sum(md$event),
      events_treat = sum(md$event[md$treatment == 1]),
      events_ctrl = sum(md$event[md$treatment == 0]),
      n_treat = sum(md$treatment == 1),
      n_ctrl = sum(md$treatment == 0),
      stringsAsFactors = FALSE
    ))
    cat(sprintf("    HR: %.3f (%.3f-%.3f) | Events: treat=%d/%d, ctrl=%d/%d\n",
                hr_info$hr, hr_info$hr_lower, hr_info$hr_upper,
                sum(md$event[md$treatment == 1]), sum(md$treatment == 1),
                sum(md$event[md$treatment == 0]), sum(md$treatment == 0)))
  }

  # --- Method 2: PSM Exact Year 1:5 ---
  cat("  Method 2: PSM Exact Year 1:5\n")
  result_2 <- tryCatch({
    run_psm_exact_year(cohort, exclude_vars, ps_vars, caliper = 0.2, ratio = 5)
  }, error = function(e) {
    cat(sprintf("    Failed: %s\n", e$message)); NULL
  })

  if (!is.null(result_2)) {
    hr_info <- extract_hr(result_2$cox)
    md <- result_2$cohort_matched
    forest_data_3 <- rbind(forest_data_3, data.frame(
      comparison = comp_name,
      comparison_display = display_comp,
      method = "PSM Exact Year 1:5",
      hr = hr_info$hr, hr_lower = hr_info$hr_lower,
      hr_upper = hr_info$hr_upper, p_value = hr_info$p_value,
      n_matched = nrow(md), n_events = sum(md$event),
      events_treat = sum(md$event[md$treatment == 1]),
      events_ctrl = sum(md$event[md$treatment == 0]),
      n_treat = sum(md$treatment == 1),
      n_ctrl = sum(md$treatment == 0),
      stringsAsFactors = FALSE
    ))
    cat(sprintf("    HR: %.3f (%.3f-%.3f) | Events: treat=%d/%d, ctrl=%d/%d\n",
                hr_info$hr, hr_info$hr_lower, hr_info$hr_upper,
                sum(md$event[md$treatment == 1]), sum(md$treatment == 1),
                sum(md$event[md$treatment == 0]), sum(md$treatment == 0)))
  }

  # --- Method 3: Pre-COVID PSM 1:5 (2018-2019 only) ---
  cat("  Method 3: Pre-COVID PSM 1:5 (2018-2019)\n")
  cohort_precovid <- cohort %>%
    filter(as.numeric(as.character(index_year)) %in% c(2018, 2019))

  cat(sprintf("    Pre-COVID cohort: N=%d (Treat=%d, Ctrl=%d)\n",
              nrow(cohort_precovid),
              sum(cohort_precovid$treatment == 1),
              sum(cohort_precovid$treatment == 0)))

  result_3 <- tryCatch({
    run_psm_no_exact_year(cohort_precovid, exclude_vars, ps_vars,
                          caliper = 0.2, ratio = 5)
  }, error = function(e) {
    cat(sprintf("    Failed: %s\n", e$message)); NULL
  })

  if (!is.null(result_3)) {
    hr_info <- extract_hr(result_3$cox)
    md <- result_3$cohort_matched
    forest_data_3 <- rbind(forest_data_3, data.frame(
      comparison = comp_name,
      comparison_display = display_comp,
      method = "Pre-COVID PSM 1:5",
      hr = hr_info$hr, hr_lower = hr_info$hr_lower,
      hr_upper = hr_info$hr_upper, p_value = hr_info$p_value,
      n_matched = nrow(md), n_events = sum(md$event),
      events_treat = sum(md$event[md$treatment == 1]),
      events_ctrl = sum(md$event[md$treatment == 0]),
      n_treat = sum(md$treatment == 1),
      n_ctrl = sum(md$treatment == 0),
      stringsAsFactors = FALSE
    ))
    cat(sprintf("    HR: %.3f (%.3f-%.3f) | Events: treat=%d/%d, ctrl=%d/%d\n",
                hr_info$hr, hr_info$hr_lower, hr_info$hr_upper,
                sum(md$event[md$treatment == 1]), sum(md$treatment == 1),
                sum(md$event[md$treatment == 0]), sum(md$treatment == 0)))
  }

  # --- Method 4: >=2yr Follow-up PSM 1:5 (index_date <= 2021-10-01) ---
  cat("  Method 4: >=2yr Follow-up PSM 1:5\n")
  cohort_2yr <- cohort %>%
    filter(as.Date(index_date) <= as.Date("2021-10-01"))

  cat(sprintf("    >=2yr FU cohort: N=%d (Treat=%d, Ctrl=%d)\n",
              nrow(cohort_2yr),
              sum(cohort_2yr$treatment == 1),
              sum(cohort_2yr$treatment == 0)))

  result_4 <- tryCatch({
    run_psm_no_exact_year(cohort_2yr, exclude_vars, ps_vars,
                          caliper = 0.2, ratio = 5)
  }, error = function(e) {
    cat(sprintf("    Failed: %s\n", e$message)); NULL
  })

  if (!is.null(result_4)) {
    hr_info <- extract_hr(result_4$cox)
    md <- result_4$cohort_matched
    forest_data_3 <- rbind(forest_data_3, data.frame(
      comparison = comp_name,
      comparison_display = display_comp,
      method = ">=2yr FU PSM 1:5",
      hr = hr_info$hr, hr_lower = hr_info$hr_lower,
      hr_upper = hr_info$hr_upper, p_value = hr_info$p_value,
      n_matched = nrow(md), n_events = sum(md$event),
      events_treat = sum(md$event[md$treatment == 1]),
      events_ctrl = sum(md$event[md$treatment == 0]),
      n_treat = sum(md$treatment == 1),
      n_ctrl = sum(md$treatment == 0),
      stringsAsFactors = FALSE
    ))
    cat(sprintf("    HR: %.3f (%.3f-%.3f) | Events: treat=%d/%d, ctrl=%d/%d\n",
                hr_info$hr, hr_info$hr_lower, hr_info$hr_upper,
                sum(md$event[md$treatment == 1]), sum(md$treatment == 1),
                sum(md$event[md$treatment == 0]), sum(md$treatment == 0)))
  }
}

cat(sprintf("\nForest data 3: %d rows\n", nrow(forest_data_3)))

# --- Forest Plot 3: Stacked Panels (Supplementary Figure Style) ---
if (nrow(forest_data_3) > 0) {
  forest_data_3 <- forest_data_3 %>%
    mutate(
      treat_name = sub(" vs .*", "", comparison_display),
      ctrl_name = sub(".* vs ", "", comparison_display),
      arm_label = sprintf("%s (%d/%d) vs %s (%d/%d)",
                          treat_name, events_treat, n_treat,
                          ctrl_name, events_ctrl, n_ctrl),
      hr_text = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper),
      p_text = ifelse(p_value < 0.001, "< 0.001", sprintf("%.3f", p_value))
    )

  # Shape map per comparison (6 comparisons)
  shape_map_3 <- c(
    "SEMAGLUTIDE vs OTHER_GLPA" = 16,
    "SEMAGLUTIDE vs SGLT2"      = 15,
    "SEMAGLUTIDE vs OtherGLD"   = 17,
    "OTHER_GLPA vs SGLT2"       = 18,
    "OTHER_GLPA vs OtherGLD"    = 8,
    "SGLT2 vs OtherGLD"         = 4
  )

  comp_order_3 <- comparisons_6

  method_levels <- c(
    "PSM Exact Year 1:1", "PSM Exact Year 1:5",
    "Pre-COVID PSM 1:5", ">=2yr FU PSM 1:5"
  )

  method_titles <- c(
    "PSM Exact Year 1:1" = "(A) PSM with exact calendar year matching (1:1)",
    "PSM Exact Year 1:5" = "(B) PSM with exact calendar year matching (1:5)",
    "Pre-COVID PSM 1:5"  = "(C) Pre-COVID era PSM (2018-2019, 1:5)",
    ">=2yr FU PSM 1:5"   = "(D) Minimum 2-year follow-up PSM (1:5)"
  )

  fp3_panels <- list()
  for (i in seq_along(method_levels)) {
    m <- method_levels[i]

    pdf <- forest_data_3 %>%
      filter(method == m) %>%
      mutate(comparison = factor(comparison, levels = comp_order_3)) %>%
      arrange(comparison)

    if (nrow(pdf) == 0) next

    pdf$arm_label <- factor(pdf$arm_label, levels = rev(pdf$arm_label))

    x_max <- 2.0
    is_bottom <- (i == length(method_levels))

    p <- ggplot(pdf, aes(x = hr, y = arm_label)) +
      geom_vline(xintercept = 1, linetype = "dashed", color = "black",
                 linewidth = 0.5) +
      geom_errorbarh(aes(xmin = hr_lower, xmax = pmin(hr_upper, x_max)),
                     height = 0.3, linewidth = 0.5) +
      geom_point(aes(shape = comparison), size = 3, fill = "black") +
      scale_shape_manual(values = shape_map_3, guide = "none") +
      geom_text(aes(x = x_max + 0.15, label = hr_text),
                hjust = 0, size = 3, fontface = "plain") +
      geom_text(aes(x = x_max + 0.85, label = p_text),
                hjust = 0, size = 3, fontface = "plain") +
      scale_x_continuous(
        limits = c(0, x_max), breaks = seq(0, 2.0, 0.5),
        expand = expansion(mult = c(0.02, 0))
      ) +
      coord_cartesian(clip = "off") +
      labs(
        title = method_titles[m],
        x = if (is_bottom) "Hazard Ratio" else NULL,
        y = NULL
      ) +
      theme_bw(base_size = 10) +
      theme(
        plot.title = element_text(face = "bold", size = 9, hjust = 0,
                                  margin = margin(b = 4)),
        plot.title.position = "plot",
        axis.text.y = element_text(size = 8, color = "black"),
        axis.text.x = element_text(size = 8),
        axis.title.x = if (is_bottom) element_text(size = 9) else element_blank(),
        panel.grid.major.y = element_blank(),
        panel.grid.minor = element_blank(),
        plot.margin = margin(t = 8, r = 120, b = 4, l = 5)
      )

    fp3_panels[[m]] <- p
  }

  if (length(fp3_panels) > 0) {
    header3 <- ggdraw() +
      draw_label("Comparison (No. of patients with event/total)",
                 x = 0.01, hjust = 0, fontface = "bold", size = 9) +
      draw_label("HR (95% CI)",
                 x = 0.72, hjust = 0, fontface = "bold", size = 9) +
      draw_label("P Value",
                 x = 0.88, hjust = 0, fontface = "bold", size = 9)

    fp3 <- plot_grid(
      header3,
      plotlist = fp3_panels,
      ncol = 1,
      rel_heights = c(0.03, rep(1, length(fp3_panels))),
      align = "v", axis = "lr"
    )

    ts <- format(Sys.time(), "%Y%m%d_%H%M%S")
    ggsave(sprintf("sensitivity_ForestPlot3_CalendarYear_%s.png", ts), fp3,
           width = 12, height = 3.5 * length(fp3_panels), dpi = 300)
    cat(sprintf("Saved: sensitivity_ForestPlot3_CalendarYear_%s.png\n", ts))
    print(fp3)
  }
}
} # end if (FALSE) — skip Analysis 3


# =============================================================================
# >>> CELL 07: Save Analysis 2 Results <<<
# =============================================================================
cat("\n", paste(rep("=", 80), collapse = ""), "\n")
cat("Saving Analysis 2 Results\n")
cat(paste(rep("=", 80), collapse = ""), "\n")

ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

# --- Combined RDS ---
all_results <- list(
  analysis_2_panels = if (exists("forest_data_2")) forest_data_2 else NULL
)

saveRDS(all_results, sprintf("outcome4_sensitivity_Analysis2_%s.rds", ts))
cat(sprintf("Saved: outcome4_sensitivity_Analysis2_%s.rds\n", ts))

# --- Forest data CSV ---
if (exists("forest_data_2") && nrow(forest_data_2) > 0) {
  write.csv(forest_data_2,
            sprintf("sensitivity_forest_data_2_panels_%s.csv", ts),
            row.names = FALSE)
  cat(sprintf("Saved: sensitivity_forest_data_2_panels_%s.csv\n", ts))
}

# --- Summary ---
cat("\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat("OUTPUT FILES GENERATED\n")
cat(paste(rep("=", 70), collapse = ""), "\n")
cat(sprintf(" 1. outcome4_sensitivity_Analysis2_%s.rds\n", ts))
cat(sprintf(" 2. sensitivity_ForestPlot2_4Panel_*.png\n"))
cat(sprintf(" 3. sensitivity_forest_data_2_panels_%s.csv\n", ts))

cat("\n=== ANALYSIS 2 COMPLETE ===\n")
