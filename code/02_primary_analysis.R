# ============================================================================
# 02_primary_analysis.R
# Primary statistical analysis: IPTW, Cox regression, TMLE
# ============================================================================

library(survival); library(survey); library(mice)
library(tmle); library(SuperLearner)

# --- Propensity Score Variables (45 covariates) ------------------------------

ps_vars <- c(
  # Demographics
  "age", "sex_cat", "raceethnicity_cat", "income", "education",
  "insurance_category", "smoking", "alcohol_category",
  "baseline_bmi_category", "baseline_hba1c", "index_year_grouped",
  # Concomitant medications
  "Anticoagulant", "Antiplatelet", "BB", "Biguanide", "CCB",
  "Diuretic", "Ezetimibe", "MRA", "OtherHTN", "RAAS", "Statin",
  "TZD", "Insulin", "SEMAGLUTIDE", "OTHER_GLPA", "SGLT2i", "SU", "DPP4i",
  # Comorbidities (Charlson-based)
  "myocardial_infarction", "congestive_heart_failure",
  "peripheral_vascular_disease", "cerebrovascular_disease",
  "chronic_pulmonary_disease", "dementia", "rheumatic_disease",
  "peptic_ulcer_disease", "hemiplegia_or_paraplegia", "hiv_infection",
  "hypoglycemia", "hyperglycemic_emergency",
  # Consolidated severity variables
  "renal_disease_severity", "liver_disease_severity",
  "diabetes_with_ophthalmic_complications",
  "diabetes_with_neurological_complications", "malignancy_status"
)

# --- 1. Multiple Imputation (MICE) ------------------------------------------
# Variables imputed: baseline_bmi, baseline_hba1c
# Method: MICE (m = 20 imputations, seed = 123)
# Predictors: treatment, demographics, medications, comorbidities
# Fallback: mean imputation if MICE fails to converge

# --- 2. IPTW (Inverse Probability of Treatment Weighting) -------------------
# PS model: logistic regression, treatment ~ ps_vars (all 45 covariates)
# Weight: ATE weights = treated/ps + (1-treated)/(1-ps)
# Trimming: weights capped at 1st and 99th percentiles
# Balance: standardized mean differences (SMD) < 0.1 target

# --- 3. Cox Proportional Hazards Regression ----------------------------------
# For each of 6 pairwise comparisons × 4 outcome definitions:
#   Model: coxph(Surv(time, event) ~ treatment,
#                weights = iptw, data = comparison_cohort)
#   Robust SE: cluster(person_id)
#   Output: HR, 95% CI, p-value

# --- 4. TMLE (Targeted Maximum Likelihood Estimation) -----------------------
# SuperLearner library:
SL_library <- c("SL.glm", "SL.glm.interaction", "SL.step",
                "SL.step.interaction", "SL.mean")

# For each comparison:
#   tmle(Y = binary_outcome,
#        A = treatment,
#        W = ps_vars,
#        Q.SL.library = SL_library,
#        g.SL.library = SL_library,
#        family = "binomial")
#
# Bootstrap: 1000 resamples for confidence intervals
# Effect measures: Risk Difference, Risk Ratio, Odds Ratio, NNT

# --- 5. Kaplan-Meier & Cumulative Incidence ----------------------------------
# IPTW-weighted survfit() for cumulative incidence curves
# 95% CI via Greenwood formula with log-log transformation
# 4-year cumulative incidence reported for all comparisons
