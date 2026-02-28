# =============================================================================
# rev00009_subgroup_forest_outcome4.R.r
# Outcome 4 Subgroup Analysis with Forest Plots
# Modified to load from IPTW results RDS file
# =============================================================================
#
# INPUT: outcome4_epilepsy_excl_hypo_all_ages_ipwt_full_results.rds
#
# OUTPUT:
#   - outcome4_forest_semaglutide_vs_othergld.pdf
#   - outcome4_forest_semaglutide_vs_sglt2i.pdf
#   - outcome4_subgroup_results_{timestamp}.csv
#   - subgroup_forest_results_{timestamp}.rds
# =============================================================================

# >>> CELL 01: Load Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(survival)
  library(ggplot2)
  library(readr)
  library(survey)
})

# Prevent function conflicts
select <- dplyr::select
mutate <- dplyr::mutate
filter <- dplyr::filter
arrange <- dplyr::arrange
group_by <- dplyr::group_by
summarise <- dplyr::summarise
ungroup <- dplyr::ungroup

cat("Libraries loaded successfully.\n")

# >>> CELL 02: Configuration
config <- list(
  data_cut_date = as.Date("2023-10-01"),
  outcome_var = "epilepsy_or_seizure_start_date",
  outcome_label = "Epilepsy/Seizure excl Hypoglycemic",

  # PS variables for IPTW (excluding comparison-specific drugs)
  base_ps_vars = c(
    "age", "sex_cat", "raceethnicity_cat", "income", "education", "insurance_category",
    "smoking", "alcohol_category", "baseline_bmi_category", "baseline_hba1c", "index_year_grouped",
    "Anticoagulant", "Antiplatelet", "BB", "Biguanide", "CCB", "Diuretic",
    "Ezetimibe", "MRA", "OtherHTN", "RAAS", "Statin", "TZD", "Insulin",
    "myocardial_infarction", "congestive_heart_failure", "peripheral_vascular_disease",
    "cerebrovascular_disease", "chronic_pulmonary_disease", "dementia",
    "rheumatic_disease", "peptic_ulcer_disease", "hemiplegia_or_paraplegia", "hiv_infection",
    "hypoglycemia", "hyperglycemic_emergency",
    "renal_disease_severity", "liver_disease_severity",
    "diabetes_with_ophthalmic_complications", "diabetes_with_neurological_complications", "malignancy_status"
  ),

  # Comparison-specific exclusions
  comparison_exclusions = list(
    "SEMAGLUTIDE vs OtherGLD" = c("SEMAGLUTIDE", "SU", "DPP4i", "TZD"),
    "SEMAGLUTIDE vs SGLT2"    = c("SEMAGLUTIDE", "SGLT2i")
  )
)

# >>> CELL 03: Helper Functions
followup_and_event <- function(df, outcome_date_col, censor_date) {
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

# IPTW function
compute_iptw <- function(df, ps_vars, trim_pct = 0.01) {
  available_vars <- intersect(ps_vars, names(df))

  # Remove variables with only one level (e.g., sex_cat in sex-specific subgroups)
  valid_vars <- available_vars[sapply(available_vars, function(v) {
    if (is.factor(df[[v]]) || is.character(df[[v]])) {
      length(unique(na.omit(df[[v]]))) > 1
    } else {
      TRUE
    }
  })]

  ps_formula <- as.formula(paste("treatment ~", paste(valid_vars, collapse = " + ")))

  ps_model <- tryCatch({
    glm(ps_formula, data = df, family = binomial())
  }, error = function(e) {
    cat("  PS model failed, using simplified model (age only)\n")
    glm(treatment ~ age, data = df, family = binomial())
  })

  df$ps <- predict(ps_model, type = "response")
  df <- df %>%
    mutate(
      iptw = ifelse(treatment == 1, 1 / ps, 1 / (1 - ps))
    )

  # Trim extreme weights
  lower <- quantile(df$iptw, trim_pct, na.rm = TRUE)
  upper <- quantile(df$iptw, 1 - trim_pct, na.rm = TRUE)
  df$iptw <- pmin(pmax(df$iptw, lower), upper)

  # Normalize weights
  df <- df %>%
    group_by(treatment) %>%
    mutate(iptw = iptw / mean(iptw, na.rm = TRUE)) %>%
    ungroup()

  return(df)
}

# Weighted Cox model
fit_weighted_cox <- function(df) {
  design <- svydesign(ids = ~1, weights = ~iptw, data = df)
  model <- tryCatch({
    svycoxph(Surv(event_time, event) ~ treatment, design = design)
  }, error = function(e) {
    cat("  Weighted Cox failed, using unweighted\n")
    coxph(Surv(event_time, event) ~ treatment, data = df)
  })
  return(model)
}

# >>> CELL 03b: Interaction P-value Function
compute_interaction_pvalues <- function(cohort_df, ps_vars, config) {
  # Prepare the full cohort with follow-up/event and IPTW
  full_df <- cohort_df %>%
    followup_and_event(config$outcome_var, config$data_cut_date) %>%
    compute_iptw(ps_vars)

  # Create subgroup variables for interaction testing
  full_df <- full_df %>%
    mutate(
      age_group = factor(ifelse(age >= 60, ">=60", "<60")),
      sex_group = factor(as.character(sex_cat)),
      race_group = factor(case_when(
        raceethnicity_cat == "Non-Hispanic White" ~ "Non-Hispanic White",
        raceethnicity_cat == "Non-Hispanic Black" ~ "Non-Hispanic Black",
        raceethnicity_cat == "Hispanic"           ~ "Hispanic",
        TRUE                                      ~ "Others"
      ))
    )

  # Define interaction tests: group_name -> variable in full_df
  interaction_vars <- list(
    "Age"            = "age_group",
    "Sex"            = "sex_group",
    "Race/Ethnicity" = "race_group"
  )

  interaction_pvals <- data.frame(
    group = character(), p_interaction = numeric(),
    stringsAsFactors = FALSE
  )

  for (grp_name in names(interaction_vars)) {
    var_name <- interaction_vars[[grp_name]]
    cat(sprintf("  Interaction test: treatment x %s (%s)...\n", var_name, grp_name))

    p_int <- tryCatch({
      design_int <- survey::svydesign(ids = ~1, weights = ~iptw, data = full_df)
      # Fit model with interaction
      fml_int <- as.formula(sprintf(
        "Surv(event_time, event) ~ treatment * %s", var_name
      ))
      model_int <- survey::svycoxph(fml_int, design = design_int)

      # Fit model without interaction
      fml_main <- as.formula(sprintf(
        "Surv(event_time, event) ~ treatment + %s", var_name
      ))
      model_main <- survey::svycoxph(fml_main, design = design_int)

      # Wald test for interaction terms via anova (LRT-like comparison)
      # Use regTermTest for the interaction term(s)
      rt <- survey::regTermTest(model_int,
                                sprintf("treatment:%s", var_name),
                                method = "Wald")
      rt$p
    }, error = function(e) {
      cat(sprintf("    Interaction test failed: %s\n", e$message))
      NA_real_
    })

    interaction_pvals <- rbind(interaction_pvals, data.frame(
      group = grp_name, p_interaction = p_int,
      stringsAsFactors = FALSE
    ))
    cat(sprintf("    P for interaction = %s\n",
                ifelse(is.na(p_int), "NA", sprintf("%.4f", p_int))))
  }

  return(interaction_pvals)
}

# >>> CELL 04: Subgroup Analysis Function
run_subgroup_analysis <- function(cohort_df, comp_name, exclude_vars, config) {

  # Define subgroups (Age, Sex, Race/Ethnicity - no BMI/Obesity)
  # Use key-based filtering instead of filter_expr for sex to handle conversion properly
  subgroup_definitions <- list(
    list(key = "age_ge_60",     group = "Age",            label = ">=60"),
    list(key = "age_lt_60",     group = "Age",            label = "<60"),
    list(key = "sex_male",      group = "Sex",            label = "Male"),
    list(key = "sex_female",    group = "Sex",            label = "Female"),
    list(key = "race_nh_white", group = "Race/Ethnicity", label = "Non-Hispanic White"),
    list(key = "race_nh_black", group = "Race/Ethnicity", label = "Non-Hispanic Black"),
    list(key = "race_hispanic", group = "Race/Ethnicity", label = "Hispanic"),
    list(key = "race_others",   group = "Race/Ethnicity", label = "Others")
  )

  # PS vars excluding comparison-specific drugs
  ps_vars <- setdiff(config$base_ps_vars, exclude_vars)

  # Compute interaction P-values on the full cohort
  cat("  Computing interaction P-values...\n")
  interaction_pvals <- compute_interaction_pvalues(cohort_df, ps_vars, config)

  # Results collector
  results <- data.frame(
    comparison = character(), group = character(), label = character(),
    n_total = integer(), n_treat = integer(), n_ctrl = integer(),
    n_events = integer(), n_events_treat = integer(), n_events_ctrl = integer(),
    hr = numeric(), hr_lower = numeric(), hr_upper = numeric(), p_value = numeric(),
    p_interaction = numeric(),
    stringsAsFactors = FALSE
  )

  # Overall analysis first
  cat(sprintf("  Overall analysis...\n"))
  overall_df <- cohort_df %>%
    followup_and_event(config$outcome_var, config$data_cut_date) %>%
    compute_iptw(ps_vars)

  overall_model <- fit_weighted_cox(overall_df)
  overall_coef <- summary(overall_model)$coefficients

  results <- rbind(results, data.frame(
    comparison = comp_name,
    group = "Overall",
    label = "Overall",
    n_total = nrow(overall_df),
    n_treat = sum(overall_df$treatment == 1),
    n_ctrl = sum(overall_df$treatment == 0),
    n_events = sum(overall_df$event),
    n_events_treat = sum(overall_df$event[overall_df$treatment == 1]),
    n_events_ctrl = sum(overall_df$event[overall_df$treatment == 0]),
    hr = exp(overall_coef[1, 1]),
    hr_lower = exp(overall_coef[1, 1] - 1.96 * overall_coef[1, 3]),
    hr_upper = exp(overall_coef[1, 1] + 1.96 * overall_coef[1, 3]),
    p_value = overall_coef[1, ncol(overall_coef)],
    p_interaction = NA_real_,
    stringsAsFactors = FALSE
  ))

  # Subgroup analyses
  for (sg_def in subgroup_definitions) {
    sg <- sg_def$key
    cat(sprintf("  %s: %s...\n", sg_def$group, sg_def$label))

    sg_df <- tryCatch({
      df_sub <- cohort_df

      # Convert sex_cat to character for proper filtering (as in reference code)
      # Add as column for dplyr filter access
      df_sub$sex_chr <- if ("sex_cat" %in% names(df_sub)) as.character(df_sub$sex_cat) else NA_character_

      # Apply subgroup-specific filter
      if (sg == "age_ge_60") {
        df_sub <- df_sub %>% dplyr::filter(age >= 60)
      } else if (sg == "age_lt_60") {
        df_sub <- df_sub %>% dplyr::filter(age < 60)
      } else if (sg == "sex_female") {
        df_sub <- df_sub %>% dplyr::filter(sex_chr == "Female")
      } else if (sg == "sex_male") {
        df_sub <- df_sub %>% dplyr::filter(sex_chr == "Male")
      } else if (sg == "race_nh_white") {
        df_sub <- df_sub %>% dplyr::filter(raceethnicity_cat == "Non-Hispanic White")
      } else if (sg == "race_nh_black") {
        df_sub <- df_sub %>% dplyr::filter(raceethnicity_cat == "Non-Hispanic Black")
      } else if (sg == "race_hispanic") {
        df_sub <- df_sub %>% dplyr::filter(raceethnicity_cat == "Hispanic")
      } else if (sg == "race_others") {
        df_sub <- df_sub %>% dplyr::filter(!raceethnicity_cat %in% c("Non-Hispanic White", "Non-Hispanic Black", "Hispanic"))
      }

      df_sub %>% followup_and_event(config$outcome_var, config$data_cut_date)
    }, error = function(e) {
      cat(sprintf("    Filter failed: %s\n", e$message))
      return(NULL)
    })

    if (is.null(sg_df) || nrow(sg_df) < 20 || sum(sg_df$event) < 3) {
      cat(sprintf("    Skipped (insufficient data: N=%d, events=%d)\n",
                  ifelse(is.null(sg_df), 0, nrow(sg_df)),
                  ifelse(is.null(sg_df), 0, sum(sg_df$event))))
      next
    }

    sg_df <- compute_iptw(sg_df, ps_vars)
    sg_model <- tryCatch({
      fit_weighted_cox(sg_df)
    }, error = function(e) {
      cat(sprintf("    Model failed: %s\n", e$message))
      return(NULL)
    })

    if (is.null(sg_model)) next

    sg_coef <- summary(sg_model)$coefficients

    # Look up interaction P-value for this group
    grp_pint <- interaction_pvals$p_interaction[interaction_pvals$group == sg_def$group]
    if (length(grp_pint) == 0) grp_pint <- NA_real_

    results <- rbind(results, data.frame(
      comparison = comp_name,
      group = sg_def$group,
      label = sg_def$label,
      n_total = nrow(sg_df),
      n_treat = sum(sg_df$treatment == 1),
      n_ctrl = sum(sg_df$treatment == 0),
      n_events = sum(sg_df$event),
      n_events_treat = sum(sg_df$event[sg_df$treatment == 1]),
      n_events_ctrl = sum(sg_df$event[sg_df$treatment == 0]),
      hr = exp(sg_coef[1, 1]),
      hr_lower = exp(sg_coef[1, 1] - 1.96 * sg_coef[1, 3]),
      hr_upper = exp(sg_coef[1, 1] + 1.96 * sg_coef[1, 3]),
      p_value = sg_coef[1, ncol(sg_coef)],
      p_interaction = grp_pint,
      stringsAsFactors = FALSE
    ))
  }

  return(results)
}

# >>> CELL 05: Forest Plot Function
create_forest_plot <- function(results_df, title, save_path = NULL) {

  # Order subgroups (no Obesity)
  group_order <- c("Overall", "Age", "Sex", "Race/Ethnicity")
  results_df$group <- factor(results_df$group, levels = group_order)
  results_df <- results_df %>% arrange(group, label)

  # Build P-interaction label per group (show once per group header)
  pint_lookup <- results_df %>%
    filter(!is.na(p_interaction)) %>%
    distinct(group, .keep_all = TRUE) %>%
    mutate(
      pint_text = ifelse(
        p_interaction < 0.001, "P interaction <0.001",
        sprintf("P interaction = %.3f", p_interaction)
      )
    ) %>%
    select(group, pint_text)

  # Create row data: interleave group headers (with P-interaction) and subgroup rows
  # Assign numeric y-position for precise ordering
  plot_data <- data.frame(
    row_label = character(), hr = numeric(), hr_lower = numeric(), hr_upper = numeric(),
    hr_text = character(), events_text = character(), pint_text = character(),
    is_header = logical(), y_pos = numeric(),
    stringsAsFactors = FALSE
  )

  y <- 0

  for (grp in rev(group_order)) {
    grp_rows <- results_df %>% filter(group == grp)
    if (nrow(grp_rows) == 0) next

    if (grp == "Overall") {
      # Overall row (no header)
      y <- y + 1
      row <- grp_rows[1, ]
      plot_data <- rbind(plot_data, data.frame(
        row_label = "Overall",
        hr = row$hr, hr_lower = row$hr_lower, hr_upper = row$hr_upper,
        hr_text = sprintf("%.2f (%.2f-%.2f)", row$hr, row$hr_lower, row$hr_upper),
        events_text = sprintf("%d/%d", row$n_events, row$n_total),
        pint_text = "",
        is_header = FALSE, y_pos = y,
        stringsAsFactors = FALSE
      ))
    } else {
      # Group header row
      y <- y + 1
      pint_label <- pint_lookup$pint_text[pint_lookup$group == grp]
      if (length(pint_label) == 0) pint_label <- ""
      plot_data <- rbind(plot_data, data.frame(
        row_label = as.character(grp),
        hr = NA, hr_lower = NA, hr_upper = NA,
        hr_text = "", events_text = "",
        pint_text = pint_label,
        is_header = TRUE, y_pos = y,
        stringsAsFactors = FALSE
      ))

      # Subgroup rows
      for (i in seq_len(nrow(grp_rows))) {
        y <- y + 1
        row <- grp_rows[i, ]
        plot_data <- rbind(plot_data, data.frame(
          row_label = paste0("    ", row$label),
          hr = row$hr, hr_lower = row$hr_lower, hr_upper = row$hr_upper,
          hr_text = sprintf("%.2f (%.2f-%.2f)", row$hr, row$hr_lower, row$hr_upper),
          events_text = sprintf("%d/%d", row$n_events, row$n_total),
          pint_text = "",
          is_header = FALSE, y_pos = y,
          stringsAsFactors = FALSE
        ))
      }
    }

    # Add spacing between groups
    y <- y + 0.4
  }

  # Determine x position for right-side annotations
  x_hr_text   <- 3.2
  x_pint_text <- 3.2

  # Create plot
  p <- ggplot(plot_data, aes(y = y_pos)) +
    geom_vline(xintercept = 1, linetype = "dashed", color = "gray50") +
    # Points and error bars (non-header rows only)
    geom_point(
      data = plot_data %>% filter(!is_header),
      aes(x = hr), size = 3, color = "#2E86AB"
    ) +
    geom_errorbarh(
      data = plot_data %>% filter(!is_header),
      aes(xmin = hr_lower, xmax = hr_upper),
      height = 0.2, color = "#2E86AB"
    ) +
    # Row labels on the left
    geom_text(
      aes(x = 0.03, label = row_label),
      hjust = 0, size = 3.5,
      fontface = ifelse(plot_data$is_header, "bold", "plain")
    ) +
    # HR text on the right
    geom_text(
      data = plot_data %>% filter(!is_header),
      aes(x = x_hr_text, label = hr_text),
      hjust = 0, size = 3, color = "gray30"
    ) +
    # P-interaction text on header rows
    geom_text(
      data = plot_data %>% filter(is_header & pint_text != ""),
      aes(x = x_pint_text, label = pint_text),
      hjust = 0, size = 3, fontface = "italic", color = "gray40"
    ) +
    scale_x_log10(
      breaks = c(0.1, 0.25, 0.5, 1, 2),
      labels = c("0.1", "0.25", "0.5", "1", "2"),
      limits = c(0.02, 6)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.05))) +
    labs(
      title = title,
      x = "Hazard Ratio (95% CI)",
      y = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank()
    )

  if (!is.null(save_path)) {
    ggsave(save_path, plot = p, width = 12, height = 8, dpi = 300)
    cat(sprintf("Forest plot saved: %s\n", save_path))
  }

  return(p)
}

# >>> CELL 06: Load Data from IPTW Results
cat("\n========== LOADING IPTW RESULTS ==========\n")

ipwt_results <- readRDS("outcome4_epilepsy_excl_hypo_all_ages_ipwt_full_results.rds")

# Extract cohorts
cohorts <- list(
  "SEMAGLUTIDE vs OtherGLD" = ipwt_results[["SEMAGLUTIDE vs OtherGLD - Epilepsy/Seizure excl Hypoglycemic"]]$cohort,
  "SEMAGLUTIDE vs SGLT2" = ipwt_results[["SEMAGLUTIDE vs SGLT2 - Epilepsy/Seizure excl Hypoglycemic"]]$cohort
)

cat(sprintf("Loaded %d cohorts:\n", length(cohorts)))
for (comp_name in names(cohorts)) {
  df <- cohorts[[comp_name]]
  cat(sprintf("  %s: N = %d (Treat: %d, Control: %d)\n",
              comp_name, nrow(df),
              sum(df$treatment == 1),
              sum(df$treatment == 0)))

  # Check sex_cat values for debugging
  cat(sprintf("    sex_cat values: %s\n",
              paste(unique(df$sex_cat), collapse = ", ")))
  cat(sprintf("    sex_cat class: %s\n", class(df$sex_cat)))
}

# >>> CELL 07: Run Subgroup Analyses
cat("\n========== RUNNING SUBGROUP ANALYSES ==========\n")

all_results <- list()

for (comp_name in names(cohorts)) {
  cat(sprintf("\n=== %s ===\n", comp_name))

  exclude_vars <- config$comparison_exclusions[[comp_name]]
  results <- run_subgroup_analysis(cohorts[[comp_name]], comp_name, exclude_vars, config)
  all_results[[comp_name]] <- results
}

# >>> CELL 08: Create Forest Plots and Save Results
cat("\n========== CREATING FOREST PLOTS ==========\n")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- "subgroup_forest_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# Create forest plots
for (comp_name in names(all_results)) {
  clean_name <- tolower(gsub(" ", "_", gsub(" vs ", "_vs_", comp_name)))
  save_path <- file.path(output_dir, sprintf("outcome4_forest_%s_%s.pdf", clean_name, timestamp))
  title <- sprintf("Subgroup Analysis: %s\nOutcome: Adult-Onset Seizure Disorder", comp_name)
  create_forest_plot(all_results[[comp_name]], title, save_path)
}

# Combine all results
combined_results <- do.call(rbind, all_results)
rownames(combined_results) <- NULL

# Format for paper
combined_results <- combined_results %>%
  mutate(
    HR_formatted = sprintf("%.2f (%.2f-%.2f)", hr, hr_lower, hr_upper),
    P_formatted = ifelse(p_value < 0.001, "<0.001",
                         ifelse(p_value < 0.01, sprintf("%.3f", p_value),
                                sprintf("%.2f", p_value))),
    P_interaction_formatted = ifelse(
      is.na(p_interaction), "",
      ifelse(p_interaction < 0.001, "<0.001",
             sprintf("%.3f", p_interaction))
    )
  )

# Save CSV
write.csv(combined_results,
          file.path(output_dir, sprintf("outcome4_subgroup_results_%s.csv", timestamp)),
          row.names = FALSE)

# Save RDS
saveRDS(list(
  results = combined_results,
  by_comparison = all_results
), file.path(output_dir, sprintf("subgroup_forest_results_%s.rds", timestamp)))

# >>> CELL 09: Summary
cat("\n")
cat("============================================================\n")
cat("SUBGROUP FOREST ANALYSIS COMPLETE\n")
cat("============================================================\n")
cat("Output directory:", output_dir, "\n\n")

cat("FILES SAVED:\n")
for (comp_name in names(all_results)) {
  clean_name <- tolower(gsub(" ", "_", gsub(" vs ", "_vs_", comp_name)))
  cat(sprintf("  - outcome4_forest_%s_%s.pdf\n", clean_name, timestamp))
}
cat(sprintf("  - outcome4_subgroup_results_%s.csv\n", timestamp))
cat(sprintf("  - subgroup_forest_results_%s.rds\n", timestamp))

cat("\n============================================================\n")
cat("RESULTS SUMMARY:\n")
cat("============================================================\n")
print(combined_results %>%
        select(comparison, group, label, n_total, n_events, HR_formatted, P_formatted, P_interaction_formatted))

cat("\n=== SUBGROUP FOREST ANALYSIS COMPLETE ===\n")
