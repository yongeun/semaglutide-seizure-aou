# ============================================================================
# phenotype_concept_ids.R
# OMOP Concept IDs used for phenotype definitions
# Platform: All of Us Researcher Workbench (CDR version 8, OMOP CDM v5)
#
# Concept hierarchy:
#   - Standard concepts (is_standard = 1): SNOMED, LOINC, RxNorm
#   - Source concepts (is_standard = 0): ICD-10-CM, ICD-9-CM, CPT4
#   - Descendant concept traversal via cb_criteria table (path matching)
#
# Method:
#   1. Anchor concepts identified by concept_id AND full_text LIKE '%_rank1]%'
#   2. Descendants found via path matching in cb_criteria hierarchy
#   3. Events queried from cb_search_all_events using discovered concept_ids
#   4. Include and exclude sets combined to define final phenotype
# ============================================================================

# =============================================================================
# OUTCOME PHENOTYPES
# =============================================================================

# --- Outcome 1: Adult-onset seizure (broad) -----------------------------------
# Includes all epilepsy diagnoses (G40) and seizure codes (R56)
# Validation: sensitivity 97.7%, specificity 44.1%, PPV 96.2%
# Most sensitive definition; may include single provoked seizures
epilepsy_seizure_std <- c(
  37208117, 4261957, 4047903, 4150299, 4026922, 43530665, 3654658,
  380378, 374023, 45757050, 4029498, 4196708, 377091, 4326435
)
epilepsy_seizure_src <- c(1572257, 1568300)
# 1572257 = ICD-10 R56.x (Convulsions, not elsewhere classified)
# 1568300 = ICD-10 G40.x (Epilepsy and recurrent seizures)

# --- Outcome 2 Components: Epilepsy diagnoses only ---------------------------
# Epilepsy-specific ICD-10 G40 codes (no R56 seizure codes)
epilepsy_dx_std <- c(
  37208117, 4261957, 4047903, 4150299, 380378, 374023, 45757050, 4326435
)
epilepsy_dx_src <- c(1568300)

# --- Outcome 2 Components: Seizure-only codes (requires >=2 dates) -----------
# R56 codes (convulsions/seizures) that are NOT epilepsy diagnoses
# Requires >= 2 distinct event dates to count as "recurrent seizure"
# Outcome 2 = epilepsy_dx UNION seizure_recurrent
seizure_only_std <- c(43530665, 3654658, 4029498, 4196708, 377091, 4026922)
seizure_only_src <- c(1572257)

# --- Outcome 3: G40-only epilepsy (high specificity) -------------------------
# Restricted to ICD-10-CM G40.0-G40.9, G40.A, G40.B only
# Validation metrics:
#   Sensitivity: 84.4%
#   Specificity: 79.4% (95% CI: 62.1-91.3%)
#   PPV: 98.3% (95% CI: 96.6-99.3%)
# Most specific definition; excludes non-epilepsy seizure codes
epilepsy_g40_src <- c(1568300)

# --- Outcome 4 (primary): Hypoglycemia exclusion (same-day) -------------------
# Outcome 4 = Outcome 1 MINUS events where seizure occurs on the
# SAME DATE as a hypoglycemia diagnosis
# Rationale: hypoglycemic seizures are metabolic, not epileptic
# This is the PRIMARY outcome definition used in the manuscript
hypoglycemia_src <- c(
  1567955,  # E16.0 - Drug-induced hypoglycemia without coma
  1567939,  # E16.1 - Other hypoglycemia
  1567922,  # E16.2 - Hypoglycemia, unspecified
  35206889, # E11.64x - T2D with hypoglycemia
  44831048, # E11.641 - T2D with hypoglycemia with coma
  44820686, # E11.649 - T2D with hypoglycemia without coma
  1567988,  # E13.64x - Other specified DM with hypoglycemia
  35206888, # E13.641 - Other DM with hypoglycemia with coma
  1567971,  # E13.649 - Other DM with hypoglycemia without coma
  35206887  # E10.64x - T1D with hypoglycemia (for co-occurring T1D patients)
)

# =============================================================================
# EXCLUSION PHENOTYPES
# =============================================================================
# Participants with these conditions BEFORE index date are excluded from cohort

# --- MCI (Mild Cognitive Impairment) -----------------------------------------
# Exclusion rationale: MCI is associated with seizure risk and may confound
mci_src <- c(44823010, 45553737, 45595932, 37402458)
# ICD-10: G31.84 (Mild cognitive impairment)

# --- ADRD (Alzheimer's Disease & Related Dementias) ---------------------------
# Exclusion rationale: Dementia is a strong independent risk factor for seizures
adrd_std <- c(
  380701,   # Alzheimer's disease
  43021816, # Frontotemporal dementia
  4046090,  # Dementia with Lewy bodies
  37109056, # Vascular dementia
  4043378,  # Dementia
  4047747,  # Pick's disease
  378419,   # Senile dementia
  37018688  # Dementia due to other causes
)
adrd_src <- c(
  44820709, 44831078, 44824152, 44829914, 44836959, 44835772, 44835825,
  44819534, 44820749, 44821811, 44829917, 44827641, 44824106, 44827644,
  44835773, 44826538, 44836954, 44821810, 44820073, 44821813, 44831083,
  44826537, 44834581, 44827645, 45591073, 35207361, 45591076, 35207358,
  35207359, 35207116, 35207328, 45538103, 35207511, 35207360, 45595842,
  45553736, 35211390, 35207356, 45605533, 45600684, 35207115, 35207357,
  45595843, 35207118, 45538000, 35207121, 44831122, 45600683
)
# ICD-10: G30.x (Alzheimer), F01.x (Vascular dementia), F02.x (Dementia in
# other diseases), F03.x (Unspecified dementia), G31.0x (Frontotemporal)

# --- Stroke -------------------------------------------------------------------
# Exclusion rationale: Post-stroke seizures are a distinct clinical entity
stroke_std <- c(
  4045734,  # Ischemic stroke
  437847,   # Cerebral infarction
  4048785,  # Intracerebral hemorrhage
  3656941,  # Non-traumatic intracerebral hemorrhage
  3656857,  # Non-traumatic subarachnoid hemorrhage
  4043734,  # Cerebral hemorrhage
  437544,   # Cerebrovascular disease
  437540,   # Cerebral ischemia
  4110192,  # Stroke
  372924,   # Cerebral embolism
  46270031, # Ischemic stroke (SNOMED)
  4110190,  # Cerebral thrombosis
  4269209,  # Cerebrovascular accident
  45772786, # Acute stroke
  4332369,  # Transient cerebral ischemia
  43531607, # Acute ischemic stroke
  373503,   # Occlusion of cerebral arteries
  4111714,  # Cerebral artery occlusion with infarction
  4045736,  # Lacunar infarction
  762933,   # Stroke sequelae
  4334245,  # Old cerebral infarction
  443454,   # Cerebral infarction unspecified
  4108356   # Cerebral artery stenosis
)
stroke_src <- c(
  1569193,  # I63.x - Cerebral infarction
  44824253, # I61.x - Intracerebral hemorrhage
  44830001, # I60.x - Subarachnoid hemorrhage
  1568724,  # I67.x - Other cerebrovascular diseases
  45571838, # I65.x - Occlusion/stenosis of precerebral arteries
  1568360,  # I66.x - Occlusion/stenosis of cerebral arteries
  44834648, # G45.x - Transient cerebral ischemic attacks
  44820873  # I69.x - Sequelae of cerebrovascular disease
)

# =============================================================================
# DIABETES PHENOTYPE
# =============================================================================
# Used for inclusion criterion: T2D diagnosis before index date
# T1D codes used as exclusion set (T2D phenotype requires absence of T1D)

# --- T2D (Type 2 Diabetes, excluding Type 1) ----------------------------------
t2d_include_std <- c(
  43531578, 443729, 4228443, 35626070, 36714116, 45770928, 43530689,
  43530690, 45770880, 45773064, 40485020, 4193704, 43531010, 201530,
  443732, 43531651, 4222415, 43531564, 4200875, 40482801, 201826, 45757277,
  4129519, 43530685, 37018728, 45769905, 4221495, 43531563
)
t2d_include_src <- c(
  44828795, 44832194, 44832193, 44836915, 44826461, 1567956, 44833366,
  44829882, 44829879, 44824073, 44831045, 44824072, 44831047, 44836914,
  44833367, 44827617, 44827616, 44836916, 44829878, 44819500, 44826460
)
# ICD-10: E11.x (T2D with various complications)

# --- T1D Exclusion Set --------------------------------------------------------
t2d_exclude_std <- c(
  37018566, 42535540, 45769873, 201254, 4224254, 45773688, 4225656, 435216,
  4225055, 45757266, 45769830, 43531009, 318712, 40484648, 45770902, 35626069,
  45763584, 201531, 40484649, 4228112, 43531008, 45757507, 443412, 45757432
)
t2d_exclude_src <- c(
  44820682, 44819504, 44832192, 44829881, 44822934, 44833368, 44824071,
  44821787, 44832190, 44819501, 44819502, 44825264, 44836918, 44822936,
  44822935, 1567940, 44832191, 44831046, 44820684, 44834549, 44820683
)
# ICD-10: E10.x (T1D with various complications)
# Method: Person must have T2D code AND NOT have T1D code

# =============================================================================
# COMORBIDITY PHENOTYPES (Charlson Comorbidity Index)
# =============================================================================
# Binary flags created for each comorbidity, assessed up to 2 years prior to index date
# Used as covariates in propensity score model (46 total covariates)
#
# Comorbidity categories:
#   myocardial_infarction         - AMI, old MI
#   congestive_heart_failure      - CHF (systolic, diastolic, combined)
#   peripheral_vascular_disease   - PAD, aortic aneurysm
#   cerebrovascular_disease       - Stroke, TIA, cerebral ischemia
#   chronic_pulmonary_disease     - COPD, asthma, emphysema
#   dementia                      - All-cause dementia
#   rheumatic_disease             - RA, SLE, other connective tissue disease
#   peptic_ulcer_disease          - Gastric/duodenal ulcer
#   hemiplegia_or_paraplegia      - Hemiplegia, paraplegia, quadriplegia
#   hiv_infection                 - HIV/AIDS
#   hypoglycemia                  - Hypoglycemic episodes
#   hyperglycemic_emergency       - DKA, HHS
#   renal_disease_severity        - 0=none, 1=mild/moderate, 2=severe/ESRD
#   liver_disease_severity        - 0=none, 1=mild, 2=moderate/severe
#   diabetes_with_ophthalmic_complications     - Diabetic retinopathy
#   diabetes_with_neurological_complications   - Diabetic neuropathy
#   malignancy_status             - 0=none, 1=any malignancy, 2=metastatic

# =============================================================================
# PHENOTYPE ASCERTAINMENT METHOD
# =============================================================================
# All phenotypes ascertained via cb_search_all_events table in OMOP CDM.
# Concept hierarchies traversed using cb_criteria table (descendant lookup):
#   1. Identify anchor concepts matching the provided IDs with rank1 flag
#   2. Traverse path hierarchy to find all descendant concepts
#   3. Query events separately for standard (is_standard=1) and
#      source (is_standard=0) concepts
#   4. Combine via UNION
# Exclusion phenotypes applied via NOT IN subquery on person_id.
#
# Event date: entry_date from cb_search_all_events
# Start date: MIN(entry_date) per person for each phenotype
# Flag: 1 if person has any qualifying event, 0 otherwise
