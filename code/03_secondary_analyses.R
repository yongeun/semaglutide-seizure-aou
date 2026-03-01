# ============================================================================
# 03_secondary_analyses.R
# Secondary, sensitivity, and supplementary analyses
# ============================================================================

library(survival); library(survey); library(lme4); library(lmerTest)

# =============================================================================
# A. COUNTERFACTUAL MEDIATION ANALYSIS
# =============================================================================
# Design: Vansteelandt-style counterfactual mediation
# Mediators: HbA1c, BMI (analyzed separately)
# Time windows: 4 periods of 12 months each (0-12, 12-24, 24-36, 36-48 months)
# Baseline covariate (X): pre-index HbA1c or BMI value

# Mediator models (chained across periods):
#   M1 ~ treatment + X
#   M2 ~ treatment + M1
#   M3 ~ treatment + M2
#   M4 ~ treatment + M3
# Missing mediator values: LOCF (last observation carried forward)

# Outcome model (piecewise Cox):
#   Surv(tstart, tstop, event) ~ treatment + M_period + X + strata(period)
#   + cluster(person_id), ties = "breslow"
# Time split at: 365, 730, 1095 days

# Counterfactual survival curves:
#   S11: treatment=1, mediator path under treatment=1 (natural course)
#   S10: treatment=1, mediator path under treatment=0 (counterfactual)
#   S00: treatment=0, mediator path under treatment=0 (natural course)
# Monte Carlo simulation: B = 200 draws per mediator path
# Proportion mediated: (S11 - S10) / (S11 - S00)

# =============================================================================
# B. SENSITIVITY ANALYSES
# =============================================================================

# B1. Four outcome definitions (see 01_cohort_definition.R)
#     All run through IPTW + Cox pipeline for all 6 comparisons

# B2. PSM excluding severe comorbidities
#     Exclusion criteria: severe renal disease OR metastatic solid tumor
#     Matching: nearest-neighbor, caliper = 0.2 × SD(logit PS)
#     Ratios: 1:1 and 1:5

# B3. Calendar year sensitivity (COVID-period imbalance)
#     - IPTW with index_year added to PS model
#     - PSM with exact matching on index_year
#     - Time-restricted analysis: pre-COVID (2018-2019) cohort only

# =============================================================================
# C. SUBGROUP ANALYSES
# =============================================================================
# Subgroups:
#   Age: <60 vs >=60
#   Sex: Male vs Female
#   Race/ethnicity: NH-White, NH-Black, Hispanic, Other
#   Obesity: BMI <30 vs >=30
#
# Method: IPTW within each subgroup (weights re-estimated)
#   svycoxph(Surv(time, event) ~ treatment, design = svydesign(...))
# Interaction test: treatment * subgroup_variable
# Visualization: forest plots with HR and 95% CI

# =============================================================================
# D. LONGITUDINAL BIOMARKER TRAJECTORIES (LMM)
# =============================================================================
# Outcomes: BMI and HbA1c trajectories over 48 months
# Model: lmer(value ~ treatment * time_years + baseline_value
#              + (time_years | person_id), data = panel)
# Key estimate: treatment × time interaction (differential slope)
# If convergence fails: random intercept only model as fallback
