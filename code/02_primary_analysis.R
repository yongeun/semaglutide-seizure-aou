# ============================================================================
# 02_primary_analysis.R
# Primary statistical analysis: IPTW, Cox regression, TMLE
#
# This script performs:
#   1. Multiple imputation (MICE) for missing BMI and HbA1c
#   2. IPTW (Inverse Probability of Treatment Weighting)
#   3. Weighted Cox proportional hazards regression
#   4. TMLE (Targeted Maximum Likelihood Estimation) with SuperLearner
#   5. Kaplan-Meier cumulative incidence estimation
#
# Inputs: Cohort dataframes from 02_cohort_construction.R
# Outputs: HR, RD, RR, OR, NNT with 95% confidence intervals
# ============================================================================

library(survival); library(survey); library(mice)
library(tmle); library(SuperLearner)

# =============================================================================
# PROPENSITY SCORE VARIABLES (46 covariates)
# =============================================================================
# These variables are included in the propensity score logistic regression model
# to estimate the probability of receiving the index drug vs comparator.
#
# Note: For each comparison, the treatment/comparator drug flags are excluded
# from the PS model (e.g., SEMAGLUTIDE and SGLT2i excluded for the
# Semaglutide vs SGLT2i comparison). Zero-variance variables are also
# automatically removed.

ps_vars <- c(
  # ---- Demographics (11 variables) ----
  "age",                          # Continuous, years
  "sex_cat",                      # 0=Male, 1=Female
  "raceethnicity_cat",            # 0=NH-White, 1=NH-Black, 2=Hispanic, 3=Other
  "income",                       # 0=<$35K, 1=$35-100K, 2=>$100K
  "education",                    # 0=HS or below, 1=Some college, 2=Advanced
  "insurance_category",           # 0=None, 1=Public, 2=Private
  "smoking",                      # 0=Never, 1=Former, 2=Current
  "alcohol_category",             # AUDIT-C: 0=Low, 1=Increased, 2=High, 3=Dependent
  "baseline_bmi_category",        # <25, 25-30, 30-35, 35-40, >=40
  "baseline_hba1c",               # Continuous, %
  "index_year_grouped",           # Non-COVID (2018-2019, 2023) vs COVID (2020-2022)

  # ---- Concomitant Medications (18 binary flags) ----
  # Any prescription in [-365, 0] days before index_date
  "Anticoagulant", "Antiplatelet", "BB", "Biguanide", "CCB",
  "Diuretic", "Ezetimibe", "MRA", "OtherHTN", "RAAS", "Statin",
  "TZD", "Insulin",
  # Drug class flags (excluded from PS when they are treatment/comparator)
  "SEMAGLUTIDE", "OTHER_GLPA", "SGLT2i", "SU", "DPP4i",

  # ---- Comorbidities (17 variables, Charlson CCI-based) ----
  "myocardial_infarction",        # Binary
  "congestive_heart_failure",     # Binary
  "peripheral_vascular_disease",  # Binary
  "cerebrovascular_disease",      # Binary
  "chronic_pulmonary_disease",    # Binary
  "dementia",                     # Binary
  "rheumatic_disease",            # Binary
  "peptic_ulcer_disease",         # Binary
  "hemiplegia_or_paraplegia",     # Binary
  "hiv_infection",                # Binary
  "hypoglycemia",                 # Binary
  "hyperglycemic_emergency",      # Binary
  # Consolidated severity variables (ordinal: 0=none, 1=mild, 2=severe)
  "renal_disease_severity",
  "liver_disease_severity",
  "diabetes_with_ophthalmic_complications",    # Binary
  "diabetes_with_neurological_complications",  # Binary
  "malignancy_status"                          # 0=none, 1=any, 2=metastatic
)

# ---- Comparison-specific PS variable exclusions ----
# Drug class variables that overlap with treatment/comparator are excluded
# Five comparisons per manuscript (target trial emulation framework)
ps_exclude_map <- list(
  "SEMAGLUTIDE_vs_OtherGLD"   = c("SEMAGLUTIDE", "TZD", "SU", "DPP4i"),
  "SEMAGLUTIDE_vs_SGLT2"      = c("SEMAGLUTIDE", "SGLT2i"),
  "OTHER_GLPA_vs_OtherGLD"    = c("OTHER_GLPA", "TZD", "SU", "DPP4i"),
  "OTHER_GLPA_vs_SGLT2"       = c("OTHER_GLPA", "SGLT2i"),
  "SGLT2_vs_OtherGLD"         = c("SGLT2i", "TZD", "SU", "DPP4i")
)

# =============================================================================
# SECTION 1: Multiple Imputation (MICE)
# =============================================================================
# Variables imputed: baseline_bmi, baseline_hba1c
# Method: Multiple Imputation by Chained Equations
# Parameters:
#   m = 20 imputations
#   seed = 123 (reproducibility)
#   printFlag = FALSE
# Predictor matrix: treatment + demographics + medications + comorbidities (~35 vars)
# Completed dataset: first imputation used for primary analysis
# Fallback: if MICE fails to converge, mean imputation applied

# =============================================================================
# SECTION 2: IPTW (Inverse Probability of Treatment Weighting)
# =============================================================================
# PS model: logistic regression
#   glm(treatment ~ ps_vars, family = binomial())
#
# Weight calculation (ATE weights):
#   Treated:  w = 1 / PS
#   Control:  w = 1 / (1 - PS)
#
# Weight trimming:
#   Lower bound: 1st percentile of raw weights
#   Upper bound: 99th percentile of raw weights
#   ipw_trimmed = pmin(pmax(ipw, lower), upper)
#
# Weight standardization (preserves sample size):
#   ipw_std = ipw_trimmed * N / sum(ipw_trimmed)
#
# Balance assessment:
#   Standardized Mean Difference (SMD) for each covariate:
#     SMD = (mean_treated - mean_control) / sqrt((var_treated + var_control) / 2)
#   Target: all |SMD| < 0.1 after weighting

# =============================================================================
# SECTION 3: Cox Proportional Hazards Regression
# =============================================================================
# For each of 5 comparisons x 4 outcome definitions x 2 cohorts:
#
# Model specification:
#   coxph(Surv(event_time, event) ~ treatment,
#         data = cohort, weights = ipw_std, robust = TRUE)
#
# Robust standard errors via sandwich estimator with cluster(person_id)
# to account for potential within-person correlation
#
# Output per comparison:
#   - Hazard Ratio (HR) = exp(coef)
#   - 95% CI = exp(coef +/- 1.96 * robust SE)
#   - p-value from robust Wald test

# =============================================================================
# SECTION 4: TMLE (Targeted Maximum Likelihood Estimation)
# =============================================================================
# Doubly robust causal inference method that combines outcome regression
# and treatment model via a targeting step for efficient estimation.
# Fit on UNWEIGHTED cohorts (TMLE handles confounding internally).

# ---- TMLE Mechanism Specification ----
# Per manuscript: GLM for both outcome and treatment mechanisms
SL_library <- c("SL.glm")  # GLM for both Q and g mechanisms

# ---- TMLE Specification ----
# For each comparison:
#   tmle(Y = binary_outcome_at_timepoint,
#        A = treatment,
#        W = ps_vars (covariate matrix),
#        Q.SL.library = SL_library,  # Outcome model (GLM)
#        g.SL.library = SL_library,  # Treatment model (GLM)
#        family = "binomial")
#
# The tmle() function internally:
#   1. Fits outcome model Q(A,W) = E[Y|A,W] via GLM
#   2. Fits treatment model g(W) = P(A=1|W) via GLM
#   3. Computes clever covariate H(A,W) = A/g(W) - (1-A)/(1-g(W))
#   4. Fits targeting submodel: logit(Q*) = logit(Q) + epsilon * H
#   5. Derives plug-in estimates of ATE from updated Q*

# ---- Bootstrap Confidence Intervals ----
# n_bootstrap = 1000
# seed = 12345
# Method: Nonparametric bootstrap (resample with replacement)
#   Each bootstrap iteration re-runs full TMLE pipeline
#   Parallel execution: future_lapply() with multisession workers
#
# From bootstrap distribution, derive:
#   - Risk Difference (ATE): mean of bootstrap ATE estimates, 2.5th/97.5th percentile CI
#   - Risk Ratio (RR): E[Y(1)] / E[Y(0)]
#   - Odds Ratio (OR): [E[Y(1)]/(1-E[Y(1)])] / [E[Y(0)]/(1-E[Y(0)])]
#   - Number Needed to Treat (NNT): 1 / |ATE|

# ---- Time-Specific Risk Differences ----
# IPTW-weighted Kaplan-Meier survival at discrete time points
# Time points (days): 182.5, 365, 547.5, 730
# Labels: "6_months", "1_year", "18_months", "2_years"
# RD at time t = Risk_treated(t) - Risk_control(t)
#             = [1 - S_treated(t)] - [1 - S_control(t)]

# =============================================================================
# SECTION 5: Kaplan-Meier & Cumulative Incidence
# =============================================================================
# IPTW-weighted Kaplan-Meier survival curves
# survfit(Surv(event_time, event) ~ treatment,
#         data = cohort, weights = ipw_std)
#
# 95% CI: Greenwood formula with log-log transformation
# 4-year cumulative incidence reported: 1 - S(1460)
