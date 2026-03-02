# ============================================================================
# 03_secondary_analyses.R
# Secondary, sensitivity, and supplementary analyses
#
# This script performs:
#   A. Counterfactual mediation analysis (Vansteelandt-style)
#   B. Sensitivity analyses (outcome definitions, calendar year, PSM,
#      min follow-up, leave-one-out, E-values)
#   C. Subgroup analyses with forest plots
#   D. Longitudinal biomarker trajectory analysis (LMM)
#
# Inputs: IPTW-weighted cohorts from 03_iptw_cox.R
# Outputs: Mediation effects, sensitivity HR tables, forest plots, LMM results
# ============================================================================

library(survival); library(survey); library(lme4); library(lmerTest)

# =============================================================================
# A. COUNTERFACTUAL MEDIATION ANALYSIS
# =============================================================================
# Design: Vansteelandt-style counterfactual mediation using piecewise Cox models
# Goal: Decompose total treatment effect into direct and indirect (mediated) effects
#
# Mediators (analyzed separately):
#   - HbA1c: glycemic control pathway
#   - BMI: weight-related pathway
#
# Temporal framework:
#   48-month follow-up divided into 4 periods of 12 months each
#   M1: 0-365 days    (0-12 months)
#   M2: 365-730 days   (12-24 months)
#   M3: 730-1095 days  (24-36 months)
#   M4: 1095-1460 days (36-48 months)
#
# ---- Mediator Models (sequential / chained) ----
# Each period's mediator depends on the previous period:
#   M1 ~ treatment + X                (X = baseline HbA1c or BMI)
#   M2 ~ treatment + M1
#   M3 ~ treatment + M2
#   M4 ~ treatment + M3
# Fitted via linear regression (lm)
# Missing mediator values handled by LOCF (Last Observation Carried Forward)
# Method for mediator aggregation: mean value within each 12-month window
#
# ---- Outcome Model (piecewise Cox) ----
# Time axis split at: 365, 730, 1095 days (creating 4 intervals)
# Model:
#   Surv(tstart, tstop, event) ~ treatment + M_period + X
#                                + strata(period)
#                                + cluster(person_id)
# where:
#   tstart, tstop = interval start/end (from survSplit)
#   M_period = mediator value for the corresponding time period
#   X = baseline mediator value (pre-index)
#   strata(period) = separate baseline hazard per interval
#   cluster(person_id) = robust SE for within-person correlation
#   ties = "breslow"
#
# ---- Counterfactual Survival Curves ----
# Three counterfactual scenarios:
#   S11: treatment=1, mediator path from treatment=1 (natural course, treated)
#   S10: treatment=1, mediator path from treatment=0 (controlled direct effect)
#   S00: treatment=0, mediator path from treatment=0 (natural course, control)
#
# To construct S10 (the counterfactual):
#   1. Predict mediator values under A=0 using the chained mediator models
#   2. Plug these "counterfactual mediator values" into outcome model with A=1
#   3. Compute survival probability
#
# Monte Carlo simulation: B = 200 draws per mediator path
#   For each draw:
#     - Sample residuals from mediator model fits
#     - Generate mediator path: M_cf = predicted + residual * rnorm(1)
#     - Compute baseline hazard from Cox model
#     - Calculate cumulative survival
#
# Seed: 123 for A=0 paths, 124 for A=1 paths
#
# ---- Effect Decomposition ----
# Total Effect (TE)    = S11 - S00 (difference in survival)
# Natural Direct Effect (NDE) = S10 - S00 (effect not through mediator)
# Natural Indirect Effect (NIE) = S11 - S10 (effect through mediator)
# Proportion Mediated (PM) = NIE / TE = (S11 - S10) / (S11 - S00)
#
# ---- Parametric Bootstrap for Uncertainty Quantification ----
# R = 100 bootstrap replicates (resampling cohort with replacement)
# B = 200 mediator simulations per replicate (Monte Carlo draws)
# Base seed = 2025 (incremented per replicate)
# Output: HR at 12, 24, 36, 48 months with bootstrap 95% CIs

# =============================================================================
# B. SENSITIVITY ANALYSES
# =============================================================================

# ---- B1. Multiple Outcome Definitions ----
# All 4 outcome definitions run through full IPTW + Cox pipeline:
#   Outcome 1: Broad (epilepsy_or_seizure)
#   Outcome 2: Refined (epilepsy_refined: epilepsy dx OR recurrent seizure >=2 dates)
#   Outcome 3: G40-only (epilepsy_g40: ICD-10 G40 codes only)
#   Outcome 4: Excl. hypoglycemic (primary: Outcome 1 minus same-day hypoglycemia)
# Each applied to all 5 comparisons x 2 cohorts (all_ages, late_onset)

# ---- B2. PSM Excluding Severe Comorbidities ----
# Additional analysis restricting to participants without:
#   - Severe renal disease (renal_disease_severity == 2)
#   - Metastatic solid tumor (malignancy_status == 2)
#
# Propensity score matching (PSM):
#   Method: Nearest-neighbor matching
#   Caliper: 0.2 x SD(logit(PS))
#   Ratios: 1:1 and 1:5 matching
#   Package: MatchIt

# ---- B3. Calendar Year Sensitivity (COVID-period Imbalance) ----
# Rationale: COVID-19 pandemic may have altered healthcare utilization patterns
# and drug prescribing. The index_year_grouped variable captures this:
#   Non-COVID: 2018, 2019, 2023
#   COVID: 2020, 2021, 2022
#
# Three approaches:
#   (a) IPTW with index_year added to PS model
#       (already included in primary analysis as index_year_grouped)
#   (b) PSM with exact matching on index_year
#       Ensures treated/control have identical calendar year distribution
#   (c) Time-restricted analysis: pre-COVID cohort only (index_date < 2020-01-01)
#       Eliminates COVID-era confounding entirely

# ---- B4. Minimum Follow-Up Restriction ----
# Restricted to participants with >= 2 years (730 days) of follow-up
# Rationale: ensures sufficient observation time for outcome ascertainment
#   filter(event_time >= 730 | event == 1)  # keep if >= 2yr FU or had event

# ---- B5. Leave-One-Out Analysis ----
# Systematic removal of one drug ingredient at a time from comparator class
# Rationale: assess whether a single drug disproportionately drives the result
# For each ingredient in comparator class:
#   Rebuild cohort excluding that ingredient, re-run IPTW + Cox
#   Report HR range across all leave-one-out iterations
# Example for Semaglutide vs Other GLD:
#   Drop glimepiride -> re-run; Drop glipizide -> re-run; etc.

# ---- B6. E-values ----
# Quantification of unmeasured confounding strength needed to explain away
# the observed association
# E-value = HR + sqrt(HR * (HR - 1)) for HR < 1 (protective association)
# Lower confidence limit E-value computed similarly
# Interpretation: an unmeasured confounder would need to be associated with
# both treatment and outcome by at least E-value to nullify the result

# =============================================================================
# C. SUBGROUP ANALYSES
# =============================================================================
# Goal: Assess treatment effect heterogeneity across pre-specified subgroups
#
# Subgroup definitions:
#   Age:            <60 vs >=60 years (at index date)
#   Sex:            Male (sex_cat=0) vs Female (sex_cat=1)
#   Race/Ethnicity: Non-Hispanic White (0), Non-Hispanic Black (1),
#                   Hispanic (2), Other (3)
#   Obesity:        BMI <30 vs >=30 (at baseline)
#
# Method:
#   IPTW re-estimated WITHIN each subgroup:
#     - PS model re-fit using only participants in the subgroup
#     - Variables with single level are automatically excluded
#       (e.g., sex_cat in sex-specific subgroups)
#     - Weights trimmed and standardized within subgroup
#
# Cox regression within each subgroup:
#   svycoxph(Surv(event_time, event) ~ treatment,
#            design = svydesign(ids = ~1, weights = ~ipw_std, data = subgroup))
#
# Interaction test:
#   Full model: Surv ~ treatment * subgroup_variable (with IPTW)
#   p-value for interaction term indicates effect modification
#
# Forest plot visualization:
#   - Point estimates (HR) with 95% CI per subgroup level
#   - Overall estimate for reference
#   - p-value for heterogeneity (interaction test)
#   - Styling: JAMA-style, blue color scheme, sans-serif font

# =============================================================================
# D. LONGITUDINAL BIOMARKER TRAJECTORIES (Linear Mixed Models)
# =============================================================================
# Goal: Assess differential trajectories of HbA1c and BMI between treatment groups
#
# Data preparation:
#   - Baseline measurement: pre-index value (time = 0)
#   - Follow-up measurements: post-index to end of follow-up
#   - Time converted from months to years for interpretable slope estimates
#   - Panel data: bmi_panel.csv, a1c_panel.csv (from 01_data_extraction_preprocessing.R)
#
# Model specification:
#   lmer(value ~ treatment * time_years + baseline_value
#        + (1 + time_years | person_id),
#        data = panel_data)
#
# Fixed effects:
#   treatment: group difference at baseline (after adjustment)
#   time_years: overall temporal trend
#   treatment * time_years: DIFFERENTIAL SLOPE (key estimate)
#     -> Positive interaction = treatment group increases more over time
#     -> Negative interaction = treatment group decreases more over time
#   baseline_value: adjustment for baseline level
#
# Random effects:
#   (1 + time_years | person_id)
#     Random intercept: individual-level baseline variation
#     Random slope: individual-level rate of change
#     Correlation between intercept and slope estimated
#
# Fallback (if convergence fails):
#   Random intercept-only model:
#   lmer(value ~ treatment * time_years + baseline_value
#        + (1 | person_id), data = panel_data)
#
# Summary metrics:
#   - Estimated slope per treatment group
#   - SE and 95% CI for treatment x time interaction
#   - p-value for differential trajectory
