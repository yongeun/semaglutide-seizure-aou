library("bigrquery")
library("tidyverse")
library("lubridate")
library("dplyr")

# Copy the file from the Google Bucket to the local workspace
name_of_file_in_bucket <- "ehr_df.csv"
my_bucket <- Sys.getenv("WORKSPACE_BUCKET")
system(paste0("gsutil cp ", my_bucket, "/data/", name_of_file_in_bucket, " ."), intern = TRUE)

# Load the EHR CSV file (with lowercase variable names)
ehr_df <- read_csv(name_of_file_in_bucket)
ehr_df <- ehr_df[!duplicated(ehr_df$person_id), ]

# Create a new column for combined race and ethnicity
ehr_df$raceethnicity <- ifelse(ehr_df$ethnicity == "Hispanic or Latino", "Hispanic",
  ifelse(ehr_df$race == "White", "Non-Hispanic White",
    ifelse(ehr_df$race == "Black or African American", "Non-Hispanic Black",
      ifelse(ehr_df$race == "Asian", "Non-Hispanic Asian",
        ifelse(ehr_df$race == "Middle Eastern or North African", "Non-Hispanic Middle Eastern or North African",
          ifelse(ehr_df$race == "Native Hawaiian or Other Pacific Islander", "Non-Hispanic Native Hawaiian or Other Pacific Islander",
            ifelse(ehr_df$race == "American Indian or Alaska Native", "Non-Hispanic American Indian or Alaska Native",
              "Other Race or Ethnicity"
            )
          )
        )
      )
    )
  )
)

# Create mother_df and add age and age_group variables
mother_df <- ehr_df %>%
  mutate(
    age = as.integer(floor(difftime(as.Date("2023-10-01"), date_of_birth, units = "days") / 365.25)),
    age_group = cut(age, breaks = c(18, 45, 65, Inf), labels = c("18-44", "45-64", ">=65"), right = FALSE, include.lowest = TRUE),
    age_group_code = as.numeric(as.factor(age_group))
  )

mother_df <- mother_df %>%
  mutate(
    sex_cat = case_when(
      sex == "Male" ~ 0,
      sex == "Female" ~ 1,
      is.na(sex) ~ 999,
      TRUE ~ 999
    ),
    race_cat = case_when(
      race == "White" ~ 0,
      race == "Black or African American" ~ 1,
      race == "Asian" ~ 2,
      race == "Middle Eastern or North African" ~ 3,
      race == "Native Hawaiian or Other Pacific Islander" ~ 4,
      race == "American Indian or Alaska Native" ~ 5,
      race == "More than one population" ~ 6,
      is.na(race) ~ 999,
      TRUE ~ 999
    ),
    ethnicity_cat = case_when(
      ethnicity == "Not Hispanic or Latino" ~ 0,
      ethnicity == "Hispanic or Latino" ~ 1,
      is.na(ethnicity) ~ 999,
      TRUE ~ 999
    ),
    raceethnicity_cat = case_when(
      ethnicity_cat == 1 ~ 2,
      ethnicity_cat == 0 & race_cat == 0 ~ 0,
      ethnicity_cat == 0 & race_cat == 1 ~ 1,
      ethnicity_cat == 0 & race_cat %in% c(2, 3, 4, 5, 6) ~ 3,
      TRUE ~ 999
    )
  )

mother_df <- mother_df %>%
  mutate(
    income = case_when(
      aname1 == "Annual Income: less 10k" ~ 0,
      aname1 == "Annual Income: 10k 25k" ~ 0,
      aname1 == "Annual Income: 25k 35k" ~ 0,
      aname1 == "Annual Income: 35k 50k" ~ 1,
      aname1 == "Annual Income: 50k 75k" ~ 1,
      aname1 == "Annual Income: 75k 100k" ~ 1,
      aname1 == "Annual Income: 100k 150k" ~ 2,
      aname1 == "Annual Income: 150k 200k" ~ 2,
      aname1 == "Annual Income: more 200k" ~ 2,
      is.na(aname1) ~ 999,
      TRUE ~ 999
    ),
    education = case_when(
      aname2 %in% c("Highest Grade: Never Attended", "Highest Grade: One Through Four", "Highest Grade: Five Through Eight", "Highest Grade: Nine Through Eleven", "Highest Grade: Twelve Or GED") ~ 0,
      aname2 == "Highest Grade: College One to Three" ~ 1,
      aname2 == "Highest Grade: College Graduate" ~ 1,
      aname2 == "Highest Grade: Advanced Degree" ~ 2,
      is.na(aname2) ~ 999,
      TRUE ~ 999
    ),
    cigs = case_when(
      aname3 == "100 Cigs Lifetime: No" ~ 0,
      aname3 == "100 Cigs Lifetime: Yes" ~ 1,
      is.na(aname3) ~ 999,
      TRUE ~ 999
    ),
    cigs_frequency = case_when(
      aname4 == "Smoke Frequency: Some Days" ~ 1,
      aname4 == "Smoke Frequency: Every Day" ~ 1,
      aname4 == "Smoke Frequency: Not At All" ~ 0,
      is.na(aname4) ~ 999,
      TRUE ~ 999
    ),
    smoking = case_when(
      cigs_frequency == 1 ~ 2, # Override: If cigs_frequency == 1, they are an active smoker (2)
      cigs == 0 ~ 0, # If cigs == 0, they are never smoker (0)
      cigs == 1 & cigs_frequency == 0 ~ 1, # If cigs == 1 and cigs_frequency == 0, they are a former smoker (1)
      TRUE ~ 999 # For all other cases, assign 999
    )
  )

mother_df <- mother_df %>%
  mutate(
    insurance_category = case_when(
      # None (0)
      ins1 == "Health Insurance: No" ~ 0,
      ins2 == "Insurance Type Update: None" ~ 0,

      # Public (1)
      ins2 %in% c(
        "Insurance Type Update: Medicare",
        "Insurance Type Update: Medicaid",
        "Insurance Type Update: Military",
        "Insurance Type Update: VA",
        "Insurance Type Update: Indian"
      ) ~ 1,

      # Private (2)
      ins2 %in% c(
        "Insurance Type Update: Employer Or Union",
        "Insurance Type Update: Purchased",
        "Insurance Type Update: Other Health Plan"
      ) ~ 2,

      # Missing in ins2 (998)
      is.na(ins2) ~ 999,

      # Invalid or unrecognized in ins2 (999)
      ins2 %in% c("PMI: Skip", "Invalid") ~ 999,

      # Unknown Type in ins1 (999)
      ins1 == "Health Insurance: Yes" &
        (is.na(ins2) | !ins2 %in% c(
          "Insurance Type Update: None",
          "Insurance Type Update: Medicare",
          "Insurance Type Update: Medicaid",
          "Insurance Type Update: Military",
          "Insurance Type Update: VA",
          "Insurance Type Update: Indian",
          "Insurance Type Update: Employer Or Union",
          "Insurance Type Update: Purchased",
          "Insurance Type Update: Other Health Plan"
        )) ~ 999,

      # Missing in ins1 (999)
      is.na(ins1) ~ 999,

      # Default (999 for any other edge case)
      TRUE ~ 999
    )
  )

# Define alcohol and drinking frequency categories
mother_df <- mother_df %>%
  mutate(
    alcohol = case_when(
      aname5 == "Alcohol Participant: No" ~ 0,
      aname5 == "Alcohol Participant: Yes" ~ 1,
      is.na(aname5) ~ 999,
      TRUE ~ 999
    ),
    alcohol_freq = case_when(
      aname6 == "Drink Frequency Past Year: Never" ~ 0,
      aname6 == "Drink Frequency Past Year: Monthly Or Less" ~ 1,
      aname6 == "Drink Frequency Past Year: 2 to 4 Per Month" ~ 2,
      aname6 == "Drink Frequency Past Year: 2 to 3 Per Week" ~ 3,
      aname6 == "Drink Frequency Past Year: 4 or More Per Week" ~ 4,
      is.na(aname6) ~ 999,
      TRUE ~ 999
    ),
    avg_daily_drink = case_when(
      aname7 == "Average Daily Drink Count: 1 or 2" ~ 0,
      aname7 == "Average Daily Drink Count: 3 or 4" ~ 1,
      aname7 == "Average Daily Drink Count: 5 or 6" ~ 2,
      aname7 == "Average Daily Drink Count: 7 to 9" ~ 3,
      aname7 == "Average Daily Drink Count: 10 or More" ~ 4,
      is.na(aname7) ~ 999,
      TRUE ~ 999
    ),
    heavy_drink_freq = case_when(
      aname8 == "6 or More Drinks Occurrence: Never In Last Year" ~ 0,
      aname8 == "6 or More Drinks Occurrence: Less Than Monthly" ~ 1,
      aname8 == "6 or More Drinks Occurrence: Monthly" ~ 2,
      aname8 == "6 or More Drinks Occurrence: Weekly" ~ 3,
      aname8 == "6 or More Drinks Occurrence: Daily" ~ 4,
      is.na(aname8) ~ 999,
      TRUE ~ 999
    ),
    missing_alcohol_data = case_when(
      alcohol == 999 ~ 1, # Missing or non-informative alcohol use
      alcohol == 1 & (avg_daily_drink == 999 | alcohol_freq == 999) ~ 1, # Insufficient data for weekly_alcohol_grams
      TRUE ~ 0 # Not missing
    )
  )

mother_df <- mother_df %>%
  mutate(
    missing_alcohol_data = case_when(
      alcohol == 999 ~ 1, # Missing alcohol use indicator
      alcohol == 1 & (alcohol_freq == 999 | avg_daily_drink == 999 | heavy_drink_freq == 999) ~ 1, # Insufficient data for drinkers
      TRUE ~ 0 # Not missing (includes non-drinkers)
    ),
    audit_c_score = case_when(
      alcohol == 0 ~ 0, # Non-drinkers get score 0
      missing_alcohol_data == 1 ~ 999, # Missing/incomplete data
      TRUE ~ alcohol_freq + avg_daily_drink + heavy_drink_freq # Valid drinkers get sum (0-12)
    ),
    alcohol_category = case_when(
      audit_c_score == 999 ~ 999, # Missing
      audit_c_score <= 4 ~ 0, # Low Risk
      audit_c_score <= 7 ~ 1, # Increased Risk
      audit_c_score <= 10 ~ 2, # High Risk
      TRUE ~ 3 # Dependent
    )
  )

generate_freq_table <- function(variable, variable_name) {
  freq_table <- table(variable, useNA = "ifany")
  freq_df <- as.data.frame(freq_table)
  # Ensure that the frequency table has exactly two columns
  if (ncol(freq_df) == 2) {
    colnames(freq_df) <- c(variable_name, "frequency")
  } else {
    warning("Frequency table does not have two columns; please check your variable.")
  }
  print(freq_df)
}

cat("### Frequency Tables: Original and Categorized Variables ###\n")

cat("\n--- SEX ---\n")
generate_freq_table(mother_df$sex, "original sex")
generate_freq_table(mother_df$sex_cat, "sex (categorized)")

cat("\n--- RACE ---\n")
generate_freq_table(mother_df$race, "original race")
generate_freq_table(mother_df$race_cat, "race (categorized)")

cat("\n--- ETHNICITY ---\n")
generate_freq_table(mother_df$ethnicity, "original ethnicity")
generate_freq_table(mother_df$ethnicity_cat, "ethnicity (categorized)")

cat("\n--- RACE & ETHNICITY ---\n")
generate_freq_table(mother_df$raceethnicity, "original race & ethnicity")
generate_freq_table(mother_df$raceethnicity_cat, "race & ethnicity (categorized)")

cat("\n--- AGE ---\n")
generate_freq_table(mother_df$age_group, "age group (categorized)")
generate_freq_table(mother_df$age_group_code, "age group code")

cat("\n--- INSURANCE ---\n")
generate_freq_table(mother_df$ins1, "original insurance")
generate_freq_table(mother_df$ins2, "original insurance 2")
generate_freq_table(mother_df$insurance_category, "insurance category")

cat("\n--- INCOME ---\n")
generate_freq_table(mother_df$aname1, "original income")
generate_freq_table(mother_df$income, "income (categorized)")

cat("\n--- EDUCATION ---\n")
generate_freq_table(mother_df$aname2, "original education")
generate_freq_table(mother_df$education, "education (categorized)")

cat("\n--- CIGARETTE USE ---\n")
generate_freq_table(mother_df$aname3, "original cigarette use")
generate_freq_table(mother_df$cigs, "cigarette use (categorized)")

generate_freq_table(mother_df$aname4, "original cigarette use frequency")
generate_freq_table(mother_df$cigs_frequency, "cigarette use frequency (categorized)")
generate_freq_table(mother_df$smoking, "smoking status (categorized)")

cat("\n--- ALCOHOL USE ---\n")
generate_freq_table(mother_df$aname5, "original alcohol use")
generate_freq_table(mother_df$alcohol, "alcohol use (categorized)")

cat("\n--- ALCOHOL FREQUENCY ---\n")
generate_freq_table(mother_df$aname6, "original alcohol frequency")
generate_freq_table(mother_df$alcohol_freq, "alcohol frequency (categorized)")

cat("\n--- AVERAGE DAILY DRINK ---\n")
generate_freq_table(mother_df$aname7, "original average daily drink")
generate_freq_table(mother_df$avg_daily_drink, "average daily drink (categorized)")

cat("\n--- HEAVY DRINK FREQUENCY ---\n")
generate_freq_table(mother_df$aname8, "original heavy drink frequency")
generate_freq_table(mother_df$heavy_drink_freq, "heavy drink frequency (categorized)")

cat("\n--- ALCOHOL EXCLUSION ---\n")
generate_freq_table(mother_df$audit_c_score, "Audit C Score")
generate_freq_table(mother_df$alcohol_category, "Alcohol Category")

cat("Rows in mother_df:", nrow(mother_df), "\n")

backup_df <- mother_df

library(bigrquery)
library(dplyr)
library(tidyr)
library(lubridate)
library(glue)
library(purrr)
library(tidyverse)

project <- Sys.getenv("GOOGLE_PROJECT")
cdr_ds <- Sys.getenv("WORKSPACE_CDR") # e.g. "fc-aou-control-ctr-2025"

###############################################################################
## 1b. G40 epilepsy: single source concept ID 1568300                       ##
##     G40-only definition — specificity 79.4%, PPV 98.3%                    ##
##     (Kee et al. / reviewer Table 3)                                       ##
###############################################################################
g40_concept_ids <- c(1568300)

disease_sets <- list(
  # ————————————————————————————————————————————————————————————————
  # 1. Epilepsy / Seizure
  # ————————————————————————————————————————————————————————————————
  epilepsy_or_seizure = list(
    include_std = c(
      37208117, 4261957, 4047903, 4150299, 4026922, 43530665, 3654658,
      380378, 374023, 45757050, 4029498, 4196708, 377091, 4326435
    ),
    include_src = c(1572257, 1568300),
    exclude_std = integer(0),
    exclude_src = integer(0)
  ),

  # 1b. Epilepsy-specific diagnoses (any single occurrence qualifies)
  epilepsy_dx = list(
    include_std = c(37208117, 4261957, 4047903, 4150299, 380378, 374023, 45757050, 4326435),
    include_src = c(1568300),
    exclude_std = integer(0),
    exclude_src = integer(0)
  ),

  # 1c. Seizure-only codes (requires >=2 distinct event dates — handled by get_person_ids_recurrent)
  seizure_only = list(
    include_std = c(43530665, 3654658, 4029498, 4196708, 377091, 4026922),
    include_src = c(1572257),
    exclude_std = integer(0),
    exclude_src = integer(0)
  ),

  # 1d. Validated G40-only epilepsy (ICD-10-CM G40.0–G40.9, G40.A, G40.B)
  #     Specificity 79.4% (62.1–91.3), PPV 98.3% (96.6–99.3)
  #     concept IDs looked up dynamically above
  epilepsy_g40 = list(
    include_std = integer(0),
    include_src = g40_concept_ids,
    exclude_std = integer(0),
    exclude_src = integer(0)
  ),

  # ————————————————————————————————————————————————————————————————
  # 2. Mild Cognitive Impairment (MCI)
  # ————————————————————————————————————————————————————————————————
  mci = list(
    include_std = integer(0),
    include_src = c(44823010, 45553737, 45595932, 37402458),
    exclude_std = integer(0),
    exclude_src = integer(0)
  ),

  # ————————————————————————————————————————————————————————————————
  # 3. Alzheimer's disease & related dementias (ADRD)
  # ————————————————————————————————————————————————————————————————
  adrd = list(
    include_std = c(
      380701, 43021816, 4046090, 37109056, 4043378, 4047747, 378419, 37018688
    ),
    include_src = c(
      44820709, 44831078, 44824152, 44829914, 44836959, 44835772, 44835825, 44819534, 44820749, 44821811, 44829917, 44827641, 44824106, 44827644, 44835773, 44826538, 44836954, 44821810, 44820073, 44821813, 44831083, 44826537, 44834581, 44827645,
      45591073, 35207361, 45591076, 35207358, 35207359, 35207116, 35207328, 45538103, 35207511, 35207360, 45595842, 45553736, 35211390, 35207356, 45605533, 45600684, 35207115, 35207357, 45595843, 35207118, 45538000, 35207121,
      44831122, 45600683
    ),
    exclude_std = integer(0),
    exclude_src = integer(0)
  ),

  # ————————————————————————————————————————————————————————————————
  # 4. Stroke
  # ————————————————————————————————————————————————————————————————
  stroke = list(
    include_std = c(
      4045734, 437847, 4048785, 3656941, 3656857, 4043734, 437544, 437540,
      4110192, 372924, 46270031, 4110190, 4269209, 45772786, 4332369, 43531607,
      373503, 4111714, 4045736, 762933, 4334245, 443454, 4108356
    ),
    include_src = c(
      1569193, 44824253, 44830001, 1568724, 45571838, 1568360, 44834648, 44820873
    ),
    exclude_std = integer(0),
    exclude_src = integer(0)
  ),

  # ————————————————————————————————————————————————————————————————
  # 6. Diabetes
  # ————————————————————————————————————————————————————————————————
  diabetes = list(
    include_std = integer(0),
    include_src = c(
      44819503, 44820685, 44832189, 1567906, 44831044, 1567956, 44828791,
      1567923, 1567940, 1567972, 44836911, 44835747
    ),
    exclude_std = integer(0),
    exclude_src = integer(0)
  ),

  # ————————————————————————————————————————————————————————————————
  # 7. T2D excluded T1D
  # ————————————————————————————————————————————————————————————————
  t2d_excluded_t1d = list(
    include_std = c(
      43531578, 443729, 4228443, 35626070, 36714116, 45770928, 43530689,
      43530690, 45770880, 45773064, 40485020, 4193704, 43531010, 201530,
      443732, 43531651, 4222415, 43531564, 4200875, 40482801, 201826, 45757277,
      4129519, 43530685, 37018728, 45769905, 4221495, 43531563
    ),
    include_src = c(
      44828795, 44832194, 44832193, 44836915, 44826461, 1567956, 44833366,
      44829882, 44829879, 44824073, 44831045, 44824072, 44831047, 44836914,
      44833367, 44827617, 44827616, 44836916, 44829878, 44819500, 44826460
    ),
    exclude_std = c(
      37018566, 42535540, 45769873, 201254, 4224254, 45773688, 4225656, 435216,
      4225055, 45757266, 45769830, 43531009, 318712, 40484648, 45770902, 35626069,
      45763584, 201531, 40484649, 4228112, 43531008, 45757507, 443412, 45757432
    ),
    exclude_src = c(
      44820682, 44819504, 44832192, 44829881, 44822934, 44833368, 44824071,
      44821787, 44832190, 44819501, 44819502, 44825264, 44836918, 44822936,
      44822935, 1567940, 44832191, 44831046, 44820684, 44834549, 44820683
    )
  ),

  # ————————————————————————————————————————————————————————————————
  # 8. Hypoglycemia (for co-occurrence with seizure)
  # ————————————————————————————————————————————————————————————————
  hypoglycemia = list(
    include_std = integer(0),
    include_src = c(
      1567955, 1567939, 1567922, 35206889, 44831048, 44820686,
      1567988, 35206888, 1567971, 35206887
    ),
    exclude_std = integer(0),
    exclude_src = integer(0)
  )
)

###############################################################################
## 2.  Helper: build a "descendants" CTE for any concept vector ###############
###############################################################################
build_cte <- function(vec, cte_name, std_flag) {
  if (!length(vec)) {
    return("")
  } # nothing to do if vector is empty
  ids <- glue_collapse(vec, sep = ", ")
  glue("
    {cte_name} AS (
      SELECT DISTINCT c.concept_id
      FROM `{cdr_ds}.cb_criteria`  c
      JOIN (SELECT CAST(id AS STRING) AS id
            FROM `{cdr_ds}.cb_criteria`
            WHERE concept_id IN ({ids})
              AND full_text LIKE '%_rank1]%') anchor
        ON  c.path LIKE CONCAT('%.', anchor.id, '.%')
        OR  c.path LIKE CONCAT('%.', anchor.id)
        OR  c.path LIKE CONCAT(anchor.id, '.%')
        OR  c.path =  anchor.id
      WHERE c.is_standard  = {std_flag}
        AND c.is_selectable = 1 )
  ")
}

###############################################################################
## 3.  Core: pull person-IDs for one disease ##################################
###############################################################################
get_person_ids <- function(dspec, d_name) { # .x = dspec, .y = name

  ## ---- 3a.  Build up to four CTEs -----------------------------------------
  ctes <- c(
    build_cte(dspec$include_std, "inc_std_desc", 1),
    build_cte(dspec$include_src, "inc_src_desc", 0),
    build_cte(dspec$exclude_std, "exc_std_desc", 1),
    build_cte(dspec$exclude_src, "exc_src_desc", 0)
  )
  ctes <- ctes[ctes != ""] # drop empties
  with_clause <- glue_collapse(ctes, sep = ",\n")

  ## ---- 3b.  Inclusion and exclusion filters -------------------------------
  has_std <- length(dspec$include_std) > 0
  has_src <- length(dspec$include_src) > 0

  inc_filter <- glue_collapse(c(
    if (has_std) {
      " (is_standard = 1 AND concept_id IN (SELECT concept_id FROM inc_std_desc)) "
    },
    if (has_src) {
      " (is_standard = 0 AND concept_id IN (SELECT concept_id FROM inc_src_desc)) "
    }
  ), sep = " OR ")

  exc_parts <- c()
  if (length(dspec$exclude_std)) {
    exc_parts <- c(
      exc_parts,
      " (is_standard = 1 AND concept_id IN (SELECT concept_id FROM exc_std_desc)) "
    )
  }
  if (length(dspec$exclude_src)) {
    exc_parts <- c(
      exc_parts,
      " (is_standard = 0 AND concept_id IN (SELECT concept_id FROM exc_src_desc)) "
    )
  }

  ## ---- KEY CHANGE: exc_filter starts with WHERE, not AND ------------------
  exc_filter <- if (length(exc_parts)) glue("
      WHERE person_id NOT IN (
        SELECT DISTINCT person_id
        FROM `{cdr_ds}.cb_search_all_events`
        WHERE {glue_collapse(exc_parts, sep = ' OR ')}
      )") else ""

  ## ---- 3c.  Full SQL -------------------------------------------------------
  sql <- glue("
  WITH
  {with_clause},

  qualified AS (
    SELECT person_id, entry_date
    FROM `{cdr_ds}.cb_search_all_events`
    WHERE {inc_filter}
  )

  SELECT
    person_id,
    MIN(entry_date) AS {d_name}_start_date
  FROM qualified
  {exc_filter}            -- now begins with WHERE (or is empty)
  GROUP BY person_id
  ")

  bq_table_download(
    bq_project_query(project, as.character(sql))
  ) %>%
    mutate(!!d_name := 1L)
}

###############################################################################
## 3b. Variant: pull person-IDs with recurrent events (>=min_events dates) ####
###############################################################################
get_person_ids_recurrent <- function(dspec, d_name, min_events = 2) {
  ## ---- Build CTEs (identical to get_person_ids) ----------------------------
  ctes <- c(
    build_cte(dspec$include_std, "inc_std_desc", 1),
    build_cte(dspec$include_src, "inc_src_desc", 0),
    build_cte(dspec$exclude_std, "exc_std_desc", 1),
    build_cte(dspec$exclude_src, "exc_src_desc", 0)
  )
  ctes <- ctes[ctes != ""]
  with_clause <- glue_collapse(ctes, sep = ",\n")

  ## ---- Inclusion filter ----------------------------------------------------
  has_std <- length(dspec$include_std) > 0
  has_src <- length(dspec$include_src) > 0

  inc_filter <- glue_collapse(c(
    if (has_std) {
      " (is_standard = 1 AND concept_id IN (SELECT concept_id FROM inc_std_desc)) "
    },
    if (has_src) {
      " (is_standard = 0 AND concept_id IN (SELECT concept_id FROM inc_src_desc)) "
    }
  ), sep = " OR ")

  ## ---- Exclusion filter ----------------------------------------------------
  exc_parts <- c()
  if (length(dspec$exclude_std)) {
    exc_parts <- c(
      exc_parts,
      " (is_standard = 1 AND concept_id IN (SELECT concept_id FROM exc_std_desc)) "
    )
  }
  if (length(dspec$exclude_src)) {
    exc_parts <- c(
      exc_parts,
      " (is_standard = 0 AND concept_id IN (SELECT concept_id FROM exc_src_desc)) "
    )
  }

  exc_filter <- if (length(exc_parts)) glue("
      WHERE person_id NOT IN (
        SELECT DISTINCT person_id
        FROM `{cdr_ds}.cb_search_all_events`
        WHERE {glue_collapse(exc_parts, sep = ' OR ')}
      )") else ""

  ## ---- SQL with event count and HAVING clause ------------------------------
  sql <- glue("
  WITH
  {with_clause},

  qualified AS (
    SELECT person_id, entry_date
    FROM `{cdr_ds}.cb_search_all_events`
    WHERE {inc_filter}
  )

  SELECT
    person_id,
    MIN(entry_date)            AS {d_name}_start_date,
    COUNT(DISTINCT entry_date) AS {d_name}_event_count
  FROM qualified
  {exc_filter}
  GROUP BY person_id
  HAVING COUNT(DISTINCT entry_date) >= {min_events}
  ")

  bq_table_download(
    bq_project_query(project, as.character(sql))
  ) %>%
    mutate(!!d_name := 1L)
}

###############################################################################
## 4.  Run for every disease & combine ########################################
###############################################################################

# 4a. Run get_person_ids on all disease_sets EXCEPT seizure_only
standard_sets <- disease_sets[setdiff(names(disease_sets), "seizure_only")]
flag_df <- purrr::imap(standard_sets, get_person_ids) |>
  reduce(full_join, by = "person_id") |>
  mutate(across(where(is.numeric), ~ replace_na(.x, 0L)))

# 4b. Run get_person_ids_recurrent on seizure_only (>=2 distinct dates)
seizure_recurrent_df <- get_person_ids_recurrent(
  dspec = disease_sets[["seizure_only"]],
  d_name = "seizure_recurrent",
  min_events = 2
)

# 4c. Combine standard flags with seizure_recurrent
flag_df <- flag_df %>%
  full_join(seizure_recurrent_df, by = "person_id") %>%
  mutate(across(where(is.numeric), ~ replace_na(.x, 0L)))

# 4d. Compute derived columns
flag_df <- flag_df %>%
  mutate(
    epilepsy_refined = as.integer((epilepsy_dx == 1) | (seizure_recurrent == 1)),
    epilepsy_refined_start_date = pmin(epilepsy_dx_start_date,
      seizure_recurrent_start_date,
      na.rm = TRUE
    ),
    hypoglycemic_seizure = as.integer(
      (epilepsy_or_seizure == 1) &
        (hypoglycemia == 1) &
        (!is.na(hypoglycemia_start_date)) &
        (!is.na(epilepsy_or_seizure_start_date)) &
        (hypoglycemia_start_date == epilepsy_or_seizure_start_date)
    ),
    epilepsy_or_seizure_without_hypoglycemic_seizure = as.integer(
      (epilepsy_or_seizure == 1) & (hypoglycemic_seizure == 0L)
    ),
    epilepsy_or_seizure_without_hypoglycemic_seizure_start_date = if_else(
      epilepsy_or_seizure_without_hypoglycemic_seizure == 1L,
      epilepsy_or_seizure_start_date,
      as.Date(NA)
    )
  )

###############################################################################
## 5.  Merge into mother_df ###################################################
###############################################################################
# 'backup_df' must already be in your workspace (it was in your earlier snippet)
disease_flag_cols <- setdiff(names(disease_sets), "seizure_only")

mother_df <- backup_df %>%
  left_join(flag_df, by = "person_id") %>%
  mutate(across(
    c(
      all_of(disease_flag_cols), "seizure_recurrent",
      "epilepsy_refined", "hypoglycemic_seizure",
      "epilepsy_or_seizure_without_hypoglycemic_seizure"
    ),
    ~ replace_na(.x, 0L)
  ))

## Quick sanity check: total counts for each disease flag
mother_df %>% summarise(across(all_of(disease_flag_cols), sum))

###############################################################################
## 5b. Diagnostic output: refined epilepsy vs original ########################
###############################################################################
cat("\n=== Cross-tabulation: epilepsy_or_seizure vs epilepsy_refined ===\n")
print(table(
  original = mother_df$epilepsy_or_seizure,
  refined = mother_df$epilepsy_refined
))

cat("\n=== Component breakdown ===\n")
cat(
  "  epilepsy_dx only (no recurrent seizure):",
  sum(mother_df$epilepsy_dx == 1 & mother_df$seizure_recurrent == 0), "\n"
)
cat(
  "  seizure_recurrent only (no epilepsy dx):",
  sum(mother_df$epilepsy_dx == 0 & mother_df$seizure_recurrent == 1), "\n"
)
cat(
  "  Both epilepsy_dx AND seizure_recurrent :",
  sum(mother_df$epilepsy_dx == 1 & mother_df$seizure_recurrent == 1), "\n"
)
cat(
  "  Original but NOT refined              :",
  sum(mother_df$epilepsy_or_seizure == 1 & mother_df$epilepsy_refined == 0), "\n"
)

cat("\n=== Validated G40-only definition (high specificity) ===\n")
cat("  epilepsy_g40 count                    :", sum(mother_df$epilepsy_g40 == 1), "\n")
print(table(
  original = mother_df$epilepsy_or_seizure,
  g40_only = mother_df$epilepsy_g40
))
print(table(
  refined = mother_df$epilepsy_refined,
  g40_only = mother_df$epilepsy_g40
))
cat(
  "  G40 but NOT in epilepsy_refined       :",
  sum(mother_df$epilepsy_g40 == 1 & mother_df$epilepsy_refined == 0), "\n"
)
cat(
  "  epilepsy_refined but NOT G40          :",
  sum(mother_df$epilepsy_refined == 1 & mother_df$epilepsy_g40 == 0), "\n"
)

cat("\n=== Hypoglycemic seizure ===\n")
cat("  Hypoglycemia count          :", sum(mother_df$hypoglycemia == 1), "\n")
cat("  Hypoglycemic seizure count  :", sum(mother_df$hypoglycemic_seizure == 1), "\n")
cat(
  "  Proportion of seizure cohort:",
  round(sum(mother_df$hypoglycemic_seizure == 1) /
    max(sum(mother_df$epilepsy_or_seizure == 1), 1), 4), "\n"
)

# Temporal ordering: hypoglycemia before vs after seizure
hypo_seizure_temporal <- mother_df %>%
  filter(
    hypoglycemic_seizure == 1,
    !is.na(hypoglycemia_start_date),
    !is.na(epilepsy_or_seizure_start_date)
  ) %>%
  mutate(
    hypo_first = hypoglycemia_start_date < epilepsy_or_seizure_start_date,
    same_day = hypoglycemia_start_date == epilepsy_or_seizure_start_date
  )
cat("  Hypoglycemia BEFORE seizure :", sum(hypo_seizure_temporal$hypo_first, na.rm = TRUE), "\n")
cat("  Same day                    :", sum(hypo_seizure_temporal$same_day, na.rm = TRUE), "\n")
cat(
  "  Hypoglycemia AFTER seizure  :",
  sum(!hypo_seizure_temporal$hypo_first & !hypo_seizure_temporal$same_day, na.rm = TRUE), "\n"
)

###############################################################################
## 6. Save into mother_df.csv file  ###########################################
###############################################################################
# This snippet assumes that you run setup first

# This code saves your dataframe into a csv file in a "data" folder in Google Bucket

# Replace df with THE NAME OF YOUR DATAFRAME
my_dataframe <- mother_df

# Replace 'test.csv' with THE NAME of the file you're going to store in the bucket (don't delete the quotation marks)
destination_filename <- "mother_df.csv"

########################################################################
##
################# DON'T CHANGE FROM HERE ###############################
##
########################################################################

# store the dataframe in current workspace
write_excel_csv(my_dataframe, destination_filename)

# Get the bucket name
my_bucket <- Sys.getenv("WORKSPACE_BUCKET")

# Copy the file from current workspace to the bucket
system(paste0("gsutil cp ./", destination_filename, " ", my_bucket, "/data/"), intern = T)

# Check if file is in the bucket
system(paste0("gsutil ls ", my_bucket, "/data/*.csv"), intern = T)
