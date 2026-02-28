# ================================================================
# Outcome 4: Epilepsy/Seizure excl Hypoglycemic
# Design: LMM slope comparison HbA1c and BMI with baseline values
# Modified to load from IPTW results RDS file
# ================================================================

library(lme4)
library(lmerTest)
library(ggplot2)
library(dplyr)
library(readr)
library(tidyr)

# Prevent function conflicts
select <- dplyr::select
mutate <- dplyr::mutate
filter <- dplyr::filter
arrange <- dplyr::arrange
group_by <- dplyr::group_by
summarise <- dplyr::summarise
ungroup <- dplyr::ungroup

# ===============================================================================
# HELPER FUNCTION: followup_and_event (if not already defined)
# ===============================================================================
followup_and_event <- function(df, outcome_date_col, censor_date, use_death = FALSE) {
  df %>%
    mutate(
      outcome_date = as.Date(.data[[outcome_date_col]]),
      censor_date = as.Date(censor_date),
      end_fu_date = case_when(
        !is.na(outcome_date) & outcome_date <= censor_date ~ outcome_date,
        TRUE ~ censor_date
      ),
      event = as.integer(!is.na(outcome_date) & outcome_date <= censor_date),
      event_time = as.numeric(end_fu_date - index_date)
    ) %>%
    filter(event_time > 0)
}

# ===============================================================================
# LONGITUDINAL LMM ANALYSIS FUNCTION
# ===============================================================================
analyze_longitudinal_biomarker <- function(cohort_list,
                                           panel_data,
                                           biomarker_name = "BMI",
                                           date_col = "weight_date",
                                           value_col = "bmi",
                                           baseline_col = "baseline_bmi") {

  results_list <- list()

  for (comp_name in names(cohort_list)) {
    cat(sprintf("\n=== %s Longitudinal Analysis: %s ===\n", biomarker_name, comp_name))

    # 1. 코호트 데이터 준비
    cohort_item <- cohort_list[[comp_name]]

    # Handle different data structures
    if (is.data.frame(cohort_item)) {
      df <- cohort_item
    } else if (is.list(cohort_item) && "cohort" %in% names(cohort_item)) {
      df <- cohort_item$cohort
    } else {
      cat("Unknown cohort structure, skipping...\n")
      next
    }

    cat(sprintf("Cohort dimensions: %d x %d\n", nrow(df), ncol(df)))

    survival_df <- followup_and_event(df, "epilepsy_or_seizure_start_date",
                                      as.Date("2023-10-01"), FALSE)

    if (nrow(survival_df) < 20) {
      cat("Insufficient cohort data\n")
      next
    }

    # 2. Longitudinal 데이터 생성 (baseline + follow-up values)
    baseline_data <- survival_df %>%
      select(person_id, treatment, index_date, end_fu_date, !!sym(baseline_col)) %>%
      mutate(person_id = as.character(person_id)) %>%
      filter(!is.na(!!sym(baseline_col))) %>%
      mutate(
        time_months = 0,
        measurement_date = index_date,
        value = !!sym(baseline_col),
        measurement_type = "baseline",
        baseline_covariate = !!sym(baseline_col)
      ) %>%
      select(person_id, treatment, index_date, end_fu_date,
             time_months, measurement_date, value, measurement_type, baseline_covariate)

    followup_data <- survival_df %>%
      select(person_id, treatment, index_date, end_fu_date, !!sym(baseline_col)) %>%
      mutate(person_id = as.character(person_id)) %>%
      left_join(
        panel_data %>% mutate(person_id = as.character(person_id)),
        by = "person_id"
      ) %>%
      filter(
        !is.na(!!sym(value_col)),
        !is.na(!!sym(date_col)),
        !!sym(date_col) > index_date,
        !!sym(date_col) <= end_fu_date
      ) %>%
      mutate(
        time_months = as.numeric(!!sym(date_col) - index_date) / 30.44,
        measurement_date = !!sym(date_col),
        value = !!sym(value_col),
        measurement_type = "follow_up",
        baseline_covariate = !!sym(baseline_col)
      ) %>%
      select(person_id, treatment, index_date, end_fu_date,
             time_months, measurement_date, value, measurement_type, baseline_covariate)

    longitudinal_data <- bind_rows(baseline_data, followup_data) %>%
      filter(!is.na(value), !is.na(baseline_covariate)) %>%
      mutate(
        treatment_label = factor(treatment, levels = c(0, 1), labels = c("Control", "Treatment")),
        time_years = time_months / 12
      ) %>%
      arrange(person_id, time_months)

    if (nrow(longitudinal_data) < 50) {
      cat(sprintf("Insufficient longitudinal measurements (%d)\n", nrow(longitudinal_data)))
      next
    }

    # 3. 데이터 요약
    patient_summary <- longitudinal_data %>%
      group_by(person_id, treatment_label) %>%
      summarise(
        n_measurements = n(),
        has_baseline = any(measurement_type == "baseline"),
        has_followup = any(measurement_type == "follow_up"),
        min_time = min(time_months),
        max_time = max(time_months),
        baseline_value = value[time_months == 0][1],
        baseline_covariate_value = first(baseline_covariate),
        .groups = "drop"
      )

    max_followup_years <- max(longitudinal_data$time_months) / 12

    cat(sprintf("Patients with data: %d\n", nrow(patient_summary)))
    cat(sprintf("Total measurements: %d\n", nrow(longitudinal_data)))
    cat(sprintf("Max follow-up: %.2f years\n", max_followup_years))

    # 4. LMM 분석
    cat(sprintf("\n=== LMM Analysis for %s ===\n", biomarker_name))

    lmm_model <- tryCatch({
      lmer(value ~ baseline_covariate + time_years * treatment_label + (1 + time_years | person_id),
           data = longitudinal_data,
           control = lmerControl(optimizer = "bobyqa"))
    }, error = function(e1) {
      cat("Random slope model failed, trying intercept only...\n")
      tryCatch({
        lmer(value ~ baseline_covariate + time_years * treatment_label + (1 | person_id),
             data = longitudinal_data)
      }, error = function(e2) {
        cat("LMM failed, using linear model...\n")
        lm(value ~ baseline_covariate + time_years * treatment_label, data = longitudinal_data)
      })
    })

    if (!is.null(lmm_model)) {
      print(summary(lmm_model))

      if (inherits(lmm_model, "lmerMod")) {
        coeffs <- fixef(lmm_model)
        coeffs_summary <- summary(lmm_model)$coefficients
      } else {
        coeffs <- coef(lmm_model)
        coeffs_summary <- summary(lmm_model)$coefficients
      }

      slope_control <- coeffs["time_years"]
      interaction_term <- "time_years:treatment_labelTreatment"

      if (interaction_term %in% names(coeffs)) {
        slope_treatment <- slope_control + coeffs[interaction_term]
        if (interaction_term %in% rownames(coeffs_summary)) {
          interaction_p <- coeffs_summary[interaction_term, ncol(coeffs_summary)]
        } else {
          interaction_p <- NA
        }
      } else {
        slope_treatment <- slope_control
        interaction_p <- NA
      }

      baseline_diff <- if ("treatment_labelTreatment" %in% names(coeffs)) {
        coeffs["treatment_labelTreatment"]
      } else { NA }

      baseline_p <- if ("treatment_labelTreatment" %in% rownames(coeffs_summary)) {
        coeffs_summary["treatment_labelTreatment", ncol(coeffs_summary)]
      } else { NA }

      baseline_coef <- if ("baseline_covariate" %in% names(coeffs)) {
        coeffs["baseline_covariate"]
      } else { NA }

      baseline_coef_p <- if ("baseline_covariate" %in% rownames(coeffs_summary)) {
        coeffs_summary["baseline_covariate", ncol(coeffs_summary)]
      } else { NA }

      cat(sprintf("\n=== %s RESULTS ===\n", biomarker_name))
      cat(sprintf("Control slope: %.4f units/year\n", slope_control))
      cat(sprintf("Treatment slope: %.4f units/year\n", slope_treatment))
      cat(sprintf("Slope difference: %.4f units/year (p=%.4f)\n",
                  slope_treatment - slope_control, interaction_p))

      results_list[[comp_name]] <- list(
        model = lmm_model,
        biomarker = biomarker_name,
        baseline_col = baseline_col,
        baseline_coef = baseline_coef,
        baseline_coef_p = baseline_coef_p,
        baseline_diff = baseline_diff,
        baseline_p = baseline_p,
        slope_control = slope_control,
        slope_treatment = slope_treatment,
        slope_difference = slope_treatment - slope_control,
        interaction_p = interaction_p,
        n_patients = nrow(patient_summary),
        n_measurements = nrow(longitudinal_data),
        max_followup_years = max_followup_years
      )
    }

    cat(paste(rep("=", 80), collapse = ""), "\n")
  }

  return(results_list)
}

# ===============================================================================
# LOAD DATA - Outcome 4 IPTW Results
# ===============================================================================
cat("\n========== LOADING OUTCOME 4 IPTW RESULTS ==========\n")

# Load IPTW results (contains cohort data)
ipwt_results <- readRDS("outcome4_epilepsy_excl_hypo_all_ages_ipwt_full_results.rds")

# Extract cohorts for analysis (Semaglutide comparisons only)
cohorts <- list(
  "SEMAGLUTIDE vs SGLT2" = ipwt_results[["SEMAGLUTIDE vs SGLT2 - Epilepsy/Seizure excl Hypoglycemic"]],
  "SEMAGLUTIDE vs OtherGLD" = ipwt_results[["SEMAGLUTIDE vs OtherGLD - Epilepsy/Seizure excl Hypoglycemic"]]
)

cat(sprintf("Loaded %d cohorts for analysis\n", length(cohorts)))

# ===============================================================================
# LOAD PANEL DATA
# ===============================================================================
cat("\n========== LOADING PANEL DATA ==========\n")

# BMI panel
cat("Loading BMI panel data...\n")
bmi_panel <- read_csv("bmi_panel.csv", show_col_types = FALSE) %>%
  mutate(
    person_id = as.character(person_id),
    weight_date = as.Date(weight_date)
  ) %>%
  select(person_id, weight_date, bmi) %>%
  distinct(person_id, weight_date, .keep_all = TRUE) %>%
  filter(!is.na(bmi), !is.na(weight_date))

cat(sprintf("BMI panel: %d measurements from %d patients\n",
            nrow(bmi_panel), length(unique(bmi_panel$person_id))))

# HbA1c panel
cat("Loading HbA1c panel data...\n")
hba1c_panel <- read_csv("a1c_panel.csv", show_col_types = FALSE) %>%
  mutate(
    person_id = as.character(person_id),
    date_of_measurement = as.Date(date_of_measurement)
  ) %>%
  select(person_id, date_of_measurement, A1c) %>%
  distinct(person_id, date_of_measurement, .keep_all = TRUE) %>%
  filter(!is.na(A1c), !is.na(date_of_measurement))

cat(sprintf("HbA1c panel: %d measurements from %d patients\n",
            nrow(hba1c_panel), length(unique(hba1c_panel$person_id))))

# ===============================================================================
# RUN LMM ANALYSES
# ===============================================================================
cat("\n========== BMI LMM ANALYSIS ==========\n")
bmi_lmm_results <- analyze_longitudinal_biomarker(
  cohort_list = cohorts,
  panel_data = bmi_panel,
  biomarker_name = "BMI",
  date_col = "weight_date",
  value_col = "bmi",
  baseline_col = "baseline_bmi"
)

cat("\n========== HbA1c LMM ANALYSIS ==========\n")
hba1c_lmm_results <- analyze_longitudinal_biomarker(
  cohort_list = cohorts,
  panel_data = hba1c_panel,
  biomarker_name = "HbA1c",
  date_col = "date_of_measurement",
  value_col = "A1c",
  baseline_col = "baseline_hba1c"
)

# ===============================================================================
# CREATE SUPPLEMENTARY TABLE FORMAT (like eTable 10, 11)
# Format: Term | Estimate (β) | 95% CI | Standard error | P-value
# ===============================================================================
create_supplementary_table <- function(lmm_result, biomarker_name, comparison_name) {

  model <- lmm_result$model

  if (inherits(model, "lmerMod")) {
    # Get coefficients with confidence intervals
    coef_summary <- summary(model)$coefficients

    # Get confidence intervals
    ci <- tryCatch({
      confint(model, method = "Wald", parm = "beta_")
    }, error = function(e) {
      # Fallback: compute CI manually
      se <- coef_summary[, "Std. Error"]
      est <- coef_summary[, "Estimate"]
      cbind("2.5 %" = est - 1.96 * se, "97.5 %" = est + 1.96 * se)
    })

    # Extract values
    terms <- rownames(coef_summary)
    estimates <- coef_summary[, "Estimate"]
    se <- coef_summary[, "Std. Error"]
    p_values <- coef_summary[, ncol(coef_summary)]  # Last column is p-value

  } else {
    # For regular lm model
    coef_summary <- summary(model)$coefficients
    ci <- confint(model)

    terms <- rownames(coef_summary)
    estimates <- coef_summary[, "Estimate"]
    se <- coef_summary[, "Std. Error"]
    p_values <- coef_summary[, "Pr(>|t|)"]
  }

  # Format p-values with significance stars
  format_pvalue <- function(p) {
    if (is.na(p)) return("NA")
    stars <- ifelse(p < 0.001, "***", ifelse(p < 0.01, "**", ifelse(p < 0.05, "*", "")))
    if (p < 0.001) {
      return(paste0("<0.001", stars))
    } else {
      return(paste0(sprintf("%.3f", p), stars))
    }
  }

  # Rename terms for clarity
  rename_term <- function(term, biomarker) {
    term <- gsub("\\(Intercept\\)", "(Intercept)", term)
    term <- gsub("time_years", "Time (years)", term)
    term <- gsub("baseline_covariate", paste0("Baseline ", biomarker), term)
    term <- gsub("treatment_labelTreatment", "Semaglutide", term)
    term <- gsub("Time \\(years\\):Semaglutide", "Time × Semaglutide", term)
    term <- gsub("Semaglutide:Time \\(years\\)", "Time × Semaglutide", term)
    return(term)
  }

  # Create table
  table_df <- data.frame(
    Term = sapply(terms, function(t) rename_term(t, biomarker_name)),
    `Estimate (β)` = sprintf("%.3f", estimates),
    `95% CI` = sprintf("%.3f  %.3f", ci[, 1], ci[, 2]),
    `Standard error` = sprintf("%.3f", se),
    `P-value` = sapply(p_values, format_pvalue),
    stringsAsFactors = FALSE,
    check.names = FALSE
  )

  return(table_df)
}

# ===============================================================================
# CREATE SUMMARY AND SAVE RESULTS
# ===============================================================================
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- "lmm_final_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- Create Supplementary Tables in the exact format ----

# For each comparison, create a combined table with HbA1c and BMI sections
for (comp_name in names(cohorts)) {
  cat(sprintf("\n=== Creating Supplementary Table for %s ===\n", comp_name))

  tables_list <- list()

  # HbA1c section
  if (comp_name %in% names(hba1c_lmm_results)) {
    hba1c_table <- create_supplementary_table(
      hba1c_lmm_results[[comp_name]],
      "HbA1c",
      comp_name
    )
    tables_list[["HbA1c"]] <- hba1c_table
    cat("HbA1c table created\n")
  }

  # BMI section
  if (comp_name %in% names(bmi_lmm_results)) {
    bmi_table <- create_supplementary_table(
      bmi_lmm_results[[comp_name]],
      "BMI",
      comp_name
    )
    tables_list[["BMI"]] <- bmi_table
    cat("BMI table created\n")
  }

  # Combine tables with section headers
  if (length(tables_list) > 0) {
    combined_rows <- list()

    for (biomarker in names(tables_list)) {
      # Add section header row
      header_row <- data.frame(
        Term = paste0("--- ", ifelse(biomarker == "HbA1c", "Glycated Hemoglobin (HbA1c, %)", "Body Mass Index (BMI)"), " ---"),
        `Estimate (β)` = "",
        `95% CI` = "",
        `Standard error` = "",
        `P-value` = "",
        stringsAsFactors = FALSE,
        check.names = FALSE
      )
      combined_rows[[length(combined_rows) + 1]] <- header_row
      combined_rows[[length(combined_rows) + 1]] <- tables_list[[biomarker]]
    }

    combined_table <- do.call(rbind, combined_rows)

    # Save as CSV
    safe_name <- gsub(" ", "_", comp_name)
    filename <- file.path(output_dir,
                          paste0("Supplementary_Table_LMM_", safe_name, "_", timestamp, ".csv"))
    write.csv(combined_table, filename, row.names = FALSE)
    cat(sprintf("Saved: %s\n", basename(filename)))
  }
}

# ---- Also create a simple summary table for quick reference ----
all_results <- list()
for (comp_name in names(bmi_lmm_results)) {
  all_results[[paste(comp_name, "BMI", sep = " - ")]] <- bmi_lmm_results[[comp_name]]
}
for (comp_name in names(hba1c_lmm_results)) {
  all_results[[paste(comp_name, "HbA1c", sep = " - ")]] <- hba1c_lmm_results[[comp_name]]
}

if (length(all_results) > 0) {
  summary_df <- data.frame(
    Comparison = names(all_results),
    Biomarker = sapply(all_results, function(x) x$biomarker),
    N_Patients = sapply(all_results, function(x) x$n_patients),
    N_Measurements = sapply(all_results, function(x) x$n_measurements),
    Slope_Control = sapply(all_results, function(x) sprintf("%.3f", x$slope_control)),
    Slope_Treatment = sapply(all_results, function(x) sprintf("%.3f", x$slope_treatment)),
    Slope_Difference = sapply(all_results, function(x) sprintf("%.3f", x$slope_difference)),
    Interaction_P = sapply(all_results, function(x) {
      p <- x$interaction_p
      if (p < 0.001) "<0.001" else sprintf("%.3f", p)
    }),
    stringsAsFactors = FALSE
  )

  write.csv(summary_df,
            file.path(output_dir, paste0("Table_LMM_Summary_", timestamp, ".csv")),
            row.names = FALSE)
}

# ---- Save RDS for future use ----
lmm_results <- list(
  bmi_results = bmi_lmm_results,
  hba1c_results = hba1c_lmm_results,
  summary = summary_df
)
saveRDS(lmm_results,
        file.path(output_dir, paste0("lmm_analysis_results_", timestamp, ".rds")))

# ---- Print final summary ----
cat("\n")
cat("============================================================\n")
cat("LMM ANALYSIS COMPLETE - SUPPLEMENTARY TABLE FORMAT\n")
cat("============================================================\n")
cat("Output directory:", output_dir, "\n\n")

cat("FILES SAVED:\n")
cat("  Supplementary Tables (eTable 10/11 format):\n")
for (comp_name in names(cohorts)) {
  safe_name <- gsub(" ", "_", comp_name)
  cat("    - Supplementary_Table_LMM_", safe_name, "_", timestamp, ".csv\n", sep="")
}
cat("\n  Summary:\n")
cat("    - Table_LMM_Summary_", timestamp, ".csv\n", sep="")
cat("    - lmm_analysis_results_", timestamp, ".rds\n", sep="")

cat("\n============================================================\n")
cat("KEY RESULTS (Time × Semaglutide interaction):\n")
cat("============================================================\n")
if (exists("summary_df")) {
  for (i in 1:nrow(summary_df)) {
    row <- summary_df[i, ]
    cat(sprintf("%s:\n", row$Comparison))
    cat(sprintf("  Slope diff: %s/year (p=%s)\n\n", row$Slope_Difference, row$Interaction_P))
  }
}

cat("\n=== LMM ANALYSIS COMPLETE ===\n")
