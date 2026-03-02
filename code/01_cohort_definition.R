# ============================================================================
# 01_cohort_definition.R
# Study parameters, drug class definitions, and cohort eligibility criteria
# for target trial emulation of semaglutide and adult-onset seizure risk
# Platform: All of Us Researcher Workbench (CDR version 8, OMOP CDM v5)
#
# This script defines the complete study configuration including:
#   1. Study parameters (exposure window, washout, follow-up)
#   2. Drug class definitions (4 classes, 20 active ingredients)
#   3. Five active-comparator cohort study structure
#   4. Variable categorization rules (demographics, SES, lifestyle)
#   5. Outcome phenotype definitions (4 outcome variants)
#   6. Inclusion/exclusion criteria
#   7. Baseline covariate specifications
#   8. Censoring rules
# ============================================================================

library(bigrquery); library(tidyverse); library(lubridate); library(glue)

# =============================================================================
# SECTION 1: Study Parameters
# =============================================================================

EXPOSURE_START  <- as.Date("2018-01-01")  # Semaglutide FDA approval: Dec 2017
EXPOSURE_END    <- as.Date("2023-10-01")  # Data cut-off / censoring date
WASHOUT_DAYS    <- 183                     # 6-month washout for new-user design
MAX_FOLLOWUP    <- 1460                    # 48-month maximum follow-up (days)
REFERENCE_DATE  <- as.Date("2023-10-01")  # Age calculation reference date

# =============================================================================
# SECTION 2: Drug Class Definitions
# =============================================================================
# All drugs matched by active ingredient name in the OMOP drug_exposure table
# Ingredient-level matching captures all formulations, doses, and routes

drug_classes <- list(
  SEMAGLUTIDE = c("semaglutide"),
  # Includes: Ozempic (injectable), Rybelsus (oral), Wegovy (injectable)

  OTHER_GLPA  = c("exenatide",       # Byetta, Bydureon
                   "liraglutide",     # Victoza, Saxenda
                   "albiglutide",     # Tanzeum (discontinued)
                   "dulaglutide",     # Trulicity
                   "lixisenatide"),   # Adlyxin

  SGLT2       = c("canagliflozin",   # Invokana
                   "empagliflozin",   # Jardiance
                   "dapagliflozin",   # Farxiga
                   "ertugliflozin"),  # Steglatro

  OtherGLD    = c("glimepiride",     # Sulfonylurea (SU)
                   "glipizide",      # Sulfonylurea (SU)
                   "glyburide",      # Sulfonylurea (SU)
                   "alogliptin",     # DPP-4 inhibitor (DPP4i)
                   "linagliptin",    # DPP-4 inhibitor (DPP4i)
                   "sitagliptin",    # DPP-4 inhibitor (DPP4i)
                   "saxagliptin",    # DPP-4 inhibitor (DPP4i)
                   "pioglitazone",   # Thiazolidinedione (TZD)
                   "rosiglitazone")  # Thiazolidinedione (TZD)
)

# =============================================================================
# SECTION 3: Five Active-Comparator Cohort Studies
# =============================================================================
# Target trial emulation framework with active-comparator new-user design
# treatment = 1 for first drug, treatment = 0 for second drug
#
# Primary comparisons:
#   1. Semaglutide vs Other GLDs
#   2. Semaglutide vs SGLT2i
# Supplementary comparisons:
#   3. Other GLP-1RAs vs Other GLDs
#   4. Other GLP-1RAs vs SGLT2i
#   5. SGLT2i vs Other GLDs

comparisons <- list(
  c("SEMAGLUTIDE", "OtherGLD"),    # Primary: Semaglutide vs older agents
  c("SEMAGLUTIDE", "SGLT2"),       # Primary: Semaglutide vs SGLT2i
  c("OTHER_GLPA",  "OtherGLD"),    # Supplementary: Other GLP-1RA vs older agents
  c("OTHER_GLPA",  "SGLT2"),       # Supplementary: Other GLP-1RA vs SGLT2i
  c("SGLT2",       "OtherGLD")     # Supplementary: SGLT2i vs older agents
)

# PS variable exclusions per comparison (avoid conditioning on treatment)
ps_exclude_map <- list(
  "SEMAGLUTIDE_vs_OtherGLD"   = c("SEMAGLUTIDE", "TZD", "SU", "DPP4i"),
  "SEMAGLUTIDE_vs_SGLT2"      = c("SEMAGLUTIDE", "SGLT2i"),
  "OTHER_GLPA_vs_OtherGLD"    = c("OTHER_GLPA", "TZD", "SU", "DPP4i"),
  "OTHER_GLPA_vs_SGLT2"       = c("OTHER_GLPA", "SGLT2i"),
  "SGLT2_vs_OtherGLD"         = c("SGLT2i", "TZD", "SU", "DPP4i")
)

# =============================================================================
# SECTION 4: Variable Categorization Rules
# =============================================================================
# All missing/unknown values coded as 999

# ---- Demographics ----
# age: continuous, calculated as floor(difftime(REFERENCE_DATE, DOB) / 365.25)
# age_group: "18-44", "45-64", ">=65" (cut points: 18, 45, 65)
# sex_cat: 0=Male, 1=Female, 999=Missing/Other
# raceethnicity_cat: 0=NH-White, 1=NH-Black, 2=Hispanic, 3=Other
#   Hispanic ethnicity takes precedence over race category

# ---- Socioeconomic Status ----
# income:
#   0 = Low (<$35K): "less 10k", "10k 25k", "25k 35k"
#   1 = Middle ($35-100K): "35k 50k", "50k 75k", "75k 100k"
#   2 = High (>$100K): "100k 150k", "150k 200k", "more 200k"
# education:
#   0 = HS or below: "Never Attended" through "Twelve Or GED"
#   1 = Some college/graduate: "College One to Three", "College Graduate"
#   2 = Advanced degree: "Advanced Degree"
# insurance_category:
#   0 = Uninsured: "Health Insurance: No" or "Insurance Type Update: None"
#   1 = Public: Medicare, Medicaid, Military, VA, Indian Health Service
#   2 = Private: Employer/Union, Purchased, Other Health Plan

# ---- Lifestyle ----
# smoking (two-step derivation):
#   Step 1: 100 cigarettes lifetime? (PPI concept 1585857)
#   Step 2: Current frequency? (PPI concept 1585860)
#   0 = Never (<100 cigs lifetime)
#   1 = Former (>=100 cigs but currently "Not At All")
#   2 = Current (>=100 cigs and currently "Some Days" or "Every Day")
#
# alcohol_category (AUDIT-C composite):
#   Components (3 PPI questions, each scored 0-4):
#     - Drink frequency past year (concept 1586201)
#     - Average daily drinks (concept 1586207)
#     - Heavy drinking 6+ drinks (concept 1586213)
#   AUDIT-C total = sum of 3 components (range 0-12)
#   0 = Low Risk (score 0-4)
#   1 = Increased Risk (score 5-7)
#   2 = High Risk (score 8-10)
#   3 = Dependent (score 11-12)
#   Non-drinkers (concept 1586198 = "No") assigned score 0

# =============================================================================
# SECTION 5: Outcome Phenotype Definitions (4 variants)
# =============================================================================
# See phenotype_concept_ids.R for complete OMOP Concept ID lists
#
# Outcome 1 - Broad (epilepsy_or_seizure):
#   Any epilepsy (ICD-10 G40.x) OR seizure (ICD-10 R56.x) diagnosis
#   Sensitivity 97.7%, specificity 44.1%, PPV 96.2%
#   Most sensitive; may include single provoked seizures
#
# Outcome 2 - Refined (epilepsy_refined):
#   Epilepsy diagnosis (G40.x only)
#   OR Recurrent seizure (R56.x with >= 2 distinct dates)
#   Balances sensitivity and specificity
#
# Outcome 3 - G40-only (epilepsy_g40):
#   ICD-10-CM G40.0-G40.9, G40.A, G40.B only
#   Sensitivity 84.4%, specificity 79.4%, PPV 98.3%
#   Most specific; excludes all seizure-NOS codes
#
# Outcome 4 - Adult-onset seizure excluding hypoglycemic seizures (PRIMARY):
#   Outcome 1 MINUS events where seizure occurs on the
#   SAME DATE as a hypoglycemia diagnosis
#   Rationale: hypoglycemic seizures are metabolically provoked,
#   not true epileptic events; their inclusion would bias toward
#   drugs that cause more hypoglycemia (e.g., sulfonylureas)

# =============================================================================
# SECTION 6: Inclusion/Exclusion Criteria
# =============================================================================
# INCLUSION:
#   1. Age >= 18 at index date
#   2. EHR data available (has_ehr_data = 1 in cb_search_person)
#   3. Type 2 diabetes diagnosis before/on index date
#      (T2D phenotype: T2D codes present AND T1D codes absent)
#   4. New user of index drug within exposure window [2018-01-01, 2023-10-01]
#   5. Clean washout: no use of index drug OR comparator drug in prior 183 days
#
# EXCLUSION:
#   1. Pre-index epilepsy or seizure (outcome start_date < index_date)
#   2. Pre-index mild cognitive impairment (MCI)
#   3. Pre-index Alzheimer's disease or related dementias (ADRD)
#   4. Pre-index stroke
#   5. [Late-onset cohort only] Pre-index seizure at any age
#
# SECONDARY ANALYSIS:
#   Late-onset seizure: restricted to participants age >= 60 at index date
#
# ANALYSIS APPROACH: Intention-to-treat (ITT)

# =============================================================================
# SECTION 7: Baseline Covariate Specifications
# =============================================================================
# All covariates measured relative to the index date

# BMI: closest measurement in [-365, 0] days before index_date
#   Source: measurement table, concept_id = 3038553
#   Valid range: 10-100 kg/m2 (values outside range excluded)
#   Categories: <25, 25-30, 30-35, 35-40, >=40

# HbA1c: closest measurement in [-365, 0] days before index_date
#   Source: measurement table, concept_id = 3004410
#   Valid range: 2-20% (values outside range excluded)
#   Used as continuous variable in PS model

# Concomitant medications: binary flags for ANY prescription in [-365, 0] days
#   (up to 1 year prior to index date)
#   Source: drug_exposure table, matched by ingredient name
#   14 drug classes: Anticoagulant, Antiplatelet, BB, Biguanide, CCB,
#     Diuretic, Ezetimibe, MRA, OtherHTN, RAAS, Statin, TZD, Insulin
#   + 5 GLD class flags: SEMAGLUTIDE, OTHER_GLPA, SGLT2i, SU, DPP4i

# Comorbidities: binary flags for ANY diagnosis in [-730, 0] days
#   (up to 2 years prior to index date)
#   Source: condition_occurrence table via cb_search_all_events
#   17 variables (see PS variable list in 02_primary_analysis.R)

# Index year: grouped as Non-COVID (2018, 2019, 2023) vs COVID (2020-2022)
#   Captures potential confounding from pandemic-era prescribing changes

# =============================================================================
# SECTION 8: Censoring Rules
# =============================================================================
# Event time = event_date - index_date (days)
# Censoring at earliest of:
#   1. Death date (from aou_death table)
#   2. EXPOSURE_END (2023-10-01, administrative censoring)
#   3. Last EHR date (latest record across all EHR domain tables)
#   4. MAX_FOLLOWUP (index_date + 1460 days = 48 months)
#   5. Crossover date (for sequential users: date of switching to comparator)
# Participants with event_time <= 0 excluded (no follow-up)
