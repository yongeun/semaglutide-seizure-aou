# ============================================================================
# 03_iptw_cox.R
# IPTW weighting, balance assessment, and Cox proportional hazards regression
#
# This script provides the full implementation of the IPTW-Cox pipeline:
#   1. MICE imputation for missing baseline BMI and HbA1c
#   2. Standardized mean difference (SMD) calculation
#   3. IPTW weight computation with trimming and standardization
#   4. Weighted Cox regression with robust standard errors
#   5. End-to-end pipeline function per comparison
#
# Inputs: Cohort dataframes from 02_cohort_construction.R
# Outputs: HR estimates, balance tables, baseline characteristics
# ============================================================================

library(survival); library(survey); library(mice)

# =============================================================================
# PROPENSITY SCORE VARIABLES (46 covariates)
# =============================================================================

ps_vars <- c(
  "age", "sex_cat", "raceethnicity_cat", "income", "education",
  "insurance_category", "smoking", "alcohol_category",
  "baseline_bmi_category", "baseline_hba1c", "index_year_grouped",
  "Anticoagulant", "Antiplatelet", "BB", "Biguanide", "CCB",
  "Diuretic", "Ezetimibe", "MRA", "OtherHTN", "RAAS", "Statin",
  "TZD", "Insulin", "SEMAGLUTIDE", "OTHER_GLPA", "SGLT2i", "SU", "DPP4i",
  "myocardial_infarction", "congestive_heart_failure",
  "peripheral_vascular_disease", "cerebrovascular_disease",
  "chronic_pulmonary_disease", "dementia", "rheumatic_disease",
  "peptic_ulcer_disease", "hemiplegia_or_paraplegia", "hiv_infection",
  "hypoglycemia", "hyperglycemic_emergency",
  "renal_disease_severity", "liver_disease_severity",
  "diabetes_with_ophthalmic_complications",
  "diabetes_with_neurological_complications", "malignancy_status"
)

# =============================================================================
# SECTION 1: MICE Imputation for Missing Baseline Values
# =============================================================================
# Missing data pattern:
#   - baseline_bmi: missing if no BMI measurement in [-365, 0] days
#   - baseline_hba1c: missing if no HbA1c measurement in [-365, 0] days
#
# MICE parameters:
#   m = 20 imputations (first completed dataset used)
#   seed = 123
#   printFlag = FALSE
#
# Predictor variables for imputation (~35 variables):
#   baseline_bmi, baseline_hba1c (to impute),
#   treatment, age, sex_cat, raceethnicity_cat, smoking, alcohol_category,
#   income, education, insurance_category,
#   Biguanide, Insulin, SEMAGLUTIDE, OTHER_GLPA, TZD, SU, DPP4i, SGLT2i,
#   RAAS, Diuretic, MRA, BB, CCB, OtherHTN,
#   myocardial_infarction, congestive_heart_failure, chronic_pulmonary_disease,
#   liver_disease_severity, renal_disease_severity,
#   diabetes_with_ophthalmic_complications,
#   diabetes_with_neurological_complications,
#   peripheral_vascular_disease, cerebrovascular_disease,
#   hypoglycemia, hyperglycemic_emergency

perform_mice_imputation <- function(df, m = 20, seed = 123) {
  set.seed(seed)
  imp_vars <- intersect(c(
    "baseline_bmi", "baseline_hba1c",
    "age", "sex_cat", "treatment", "raceethnicity_cat", "smoking",
    "alcohol_category", "income", "education", "insurance_category",
    "Biguanide", "Insulin", "SEMAGLUTIDE", "OTHER_GLPA", "TZD", "SU",
    "DPP4i", "SGLT2i", "RAAS", "Diuretic", "MRA", "BB", "CCB", "OtherHTN",
    "myocardial_infarction", "congestive_heart_failure",
    "chronic_pulmonary_disease", "liver_disease_severity",
    "renal_disease_severity", "diabetes_with_ophthalmic_complications",
    "diabetes_with_neurological_complications",
    "peripheral_vascular_disease", "cerebrovascular_disease",
    "hypoglycemia", "hyperglycemic_emergency"
  ), names(df))

  imp <- mice(df[imp_vars], m = m, seed = seed, printFlag = FALSE)
  completed <- complete(imp, 1)

  # Replace missing values with imputed values
  for (var in c("baseline_bmi", "baseline_hba1c")) {
    if (var %in% names(completed)) df[[var]] <- completed[[var]]
  }
  df
}

# =============================================================================
# SECTION 2: Standardized Mean Difference (SMD)
# =============================================================================
# SMD = (mean_treated - mean_control) / sqrt((var_treated + var_control) / 2)
#
# Computed twice:
#   (a) Before weighting: raw data (unweighted)
#   (b) After weighting: using IPTW weights via weighted.mean()
# Target: all |SMD| < 0.1 after weighting indicates acceptable balance

calculate_smd <- function(df, vars, weights = NULL) {
  smd_results <- map_dfr(vars, function(var) {
    if (!var %in% names(df)) return(NULL)
    x1 <- df[[var]][df$treatment == 1]
    x0 <- df[[var]][df$treatment == 0]

    if (is.null(weights)) {
      m1 <- mean(x1, na.rm = TRUE); m0 <- mean(x0, na.rm = TRUE)
      v1 <- var(x1, na.rm = TRUE);  v0 <- var(x0, na.rm = TRUE)
    } else {
      w1 <- weights[df$treatment == 1]; w0 <- weights[df$treatment == 0]
      m1 <- weighted.mean(x1, w1, na.rm = TRUE)
      m0 <- weighted.mean(x0, w0, na.rm = TRUE)
      v1 <- sum(w1 * (x1 - m1)^2, na.rm = TRUE) / sum(w1)
      v0 <- sum(w0 * (x0 - m0)^2, na.rm = TRUE) / sum(w0)
    }

    pooled_sd <- sqrt((v1 + v0) / 2)
    tibble(variable = var, smd = ifelse(pooled_sd > 0, (m1 - m0) / pooled_sd, 0))
  })
  smd_results
}

# =============================================================================
# SECTION 3: IPTW Weight Computation
# =============================================================================
# Steps:
#   1. Identify comparison-specific PS variables (exclude treatment drug flags)
#   2. Drop zero-variance covariates (only one unique value within comparison)
#   3. Fit logistic regression: treatment ~ covariates
#   4. Compute ATE weights: treated -> 1/PS, control -> 1/(1-PS)
#   5. Trim weights at 1st and 99th percentiles
#   6. Standardize: ipw_std = ipw_trimmed * N / sum(ipw_trimmed)

compute_iptw <- function(df, exclude_vars = NULL, trim = 0.01) {
  rhs_vars <- setdiff(ps_vars, exclude_vars)

  # Drop zero-variance covariates
  keep_vars <- rhs_vars[sapply(df[rhs_vars], function(x) length(unique(na.omit(x))) > 1)]

  # Propensity score model (logistic regression)
  ps_form <- as.formula(paste("treatment ~", paste(keep_vars, collapse = " + ")))
  ps_model <- glm(ps_form, data = df, family = binomial())
  df$ps <- predict(ps_model, type = "response")

  # ATE weights
  df$ipw <- ifelse(df$treatment == 1, 1 / df$ps, 1 / (1 - df$ps))

  # Trim at specified percentiles to reduce extreme weight influence
  lower <- quantile(df$ipw, trim)
  upper <- quantile(df$ipw, 1 - trim)
  df$ipw_trimmed <- pmin(pmax(df$ipw, lower), upper)

  # Standardize to preserve effective sample size
  df$ipw_std <- df$ipw_trimmed * nrow(df) / sum(df$ipw_trimmed)

  list(cohort = df, ps_model = ps_model, keep_vars = keep_vars)
}

# =============================================================================
# SECTION 4: Weighted Cox Proportional Hazards Regression
# =============================================================================
# Model: Surv(event_time, event) ~ treatment
# Weights: IPTW standardized weights (ipw_std)
# Robust SE: sandwich estimator (robust = TRUE)
#
# Output:
#   HR = exp(coefficient for treatment)
#   95% CI = exp(coef +/- 1.96 * robust SE)
#   p-value from robust Wald test

run_weighted_cox <- function(df) {
  fit <- coxph(Surv(event_time, event) ~ treatment,
               data = df, weights = ipw_std, robust = TRUE)
  summary(fit)
}

# =============================================================================
# SECTION 5: End-to-End Analysis Pipeline (Per Comparison)
# =============================================================================
# Runs the full IPTW-Cox pipeline for one comparison:
#   1. MICE imputation -> 2. IPTW weights -> 3. Balance check -> 4. Cox model
#
# Arguments:
#   cohort: analysis-ready dataframe with treatment, covariates, event_time, event
#   comparison_name: string label (e.g., "SEMAGLUTIDE_vs_OtherGLD")
#   outcome_var: outcome column name (e.g., "epilepsy_or_seizure_start_date")
#   exclude_vars: PS variables to exclude for this comparison
#
# Returns: list with cox results, SMD tables, weighted cohort, sample size

run_iptw_cox_pipeline <- function(cohort, comparison_name, outcome_var,
                                  exclude_vars = NULL) {
  # Step 1: MICE imputation (m=20, seed=123)
  cohort <- perform_mice_imputation(cohort)

  # Step 2: IPTW weighting (trim at 1st/99th percentiles)
  iptw_result <- compute_iptw(cohort, exclude_vars = exclude_vars)
  df_weighted <- iptw_result$cohort

  # Step 3: Balance check (SMD before and after weighting)
  smd_before <- calculate_smd(df_weighted, iptw_result$keep_vars)
  smd_after  <- calculate_smd(df_weighted, iptw_result$keep_vars,
                               weights = df_weighted$ipw_std)

  # Step 4: Cox regression with IPTW weights
  cox_result <- run_weighted_cox(df_weighted)

  list(
    cox = cox_result,
    smd_before = smd_before,
    smd_after = smd_after,
    cohort = df_weighted,
    n = nrow(df_weighted)
  )
}
