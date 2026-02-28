# =============================================================================
# rev00014_4year_cumulative_incidence.R
# Calculate 4-Year Cumulative Incidence for All Group Comparisons
# =============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
  library(tidyr)
})

# =============================================================================
# Configuration
# =============================================================================

# RDS file paths - modify these paths as needed
rds_files <- list(
  all_ages = "epilepsy_seizure_all_ages_ipwt_full_results.rds",
  late_onset = "epilepsy_seizure_late_onset_ipwt_full_results.rds"
)

# Target time point (48 months = 4 years)
TARGET_MONTHS <- 48

# =============================================================================
# Function: Calculate Cumulative Incidence at Specific Time Point
# =============================================================================
calculate_cumulative_incidence <- function(cohort_df, time_months = 48, weight_var = "ipw_std") {

  # Convert event_time to months if in days
  if (!"event_time_months" %in% names(cohort_df)) {
    cohort_df$event_time_months <- cohort_df$event_time / 30.44
  }

  # Check if weight variable exists
  if (!weight_var %in% names(cohort_df)) {
    warning(sprintf("Weight variable '%s' not found, using unweighted analysis", weight_var))
    cohort_df[[weight_var]] <- 1
  }

  results <- list()

  for (trt in c(0, 1)) {
    subset_df <- cohort_df %>% filter(treatment == trt)

    if (nrow(subset_df) < 2) {
      results[[as.character(trt)]] <- list(
        n = nrow(subset_df),
        events = sum(subset_df$event),
        cum_inc = NA,
        ci_lower = NA,
        ci_upper = NA,
        se = NA
      )
      next
    }

    # Fit weighted Kaplan-Meier
    fit <- survfit(
      Surv(event_time_months, event) ~ 1,
      data = subset_df,
      weights = subset_df[[weight_var]]
    )

    # Get summary at target time
    summ <- summary(fit, times = time_months, extend = TRUE)

    # Cumulative incidence = 1 - Survival
    if (length(summ$surv) > 0) {
      surv <- summ$surv[1]
      se_surv <- summ$std.err[1]

      cum_inc <- 1 - surv
      # SE for cumulative incidence is same as SE for survival
      se_cum_inc <- se_surv

      # 95% CI (using log-log transformation for stability)
      ci_lower <- 1 - min(1, exp(log(surv) + 1.96 * se_surv / surv))
      ci_upper <- 1 - max(0, exp(log(surv) - 1.96 * se_surv / surv))

      # Ensure bounds are valid
      ci_lower <- max(0, ci_lower)
      ci_upper <- min(1, ci_upper)

    } else {
      cum_inc <- NA
      se_cum_inc <- NA
      ci_lower <- NA
      ci_upper <- NA
    }

    results[[as.character(trt)]] <- list(
      n = nrow(subset_df),
      events = sum(subset_df$event),
      cum_inc = cum_inc,
      ci_lower = ci_lower,
      ci_upper = ci_upper,
      se = se_cum_inc
    )
  }

  return(results)
}

# =============================================================================
# Function: Format Cumulative Incidence as Percentage String
# =============================================================================
format_cum_inc <- function(cum_inc_result, digits = 2) {
  if (is.na(cum_inc_result$cum_inc)) {
    return("NA")
  }

  pct <- cum_inc_result$cum_inc * 100
  ci_low <- cum_inc_result$ci_lower * 100
  ci_high <- cum_inc_result$ci_upper * 100

  sprintf("%.*f%% (%.*f-%.*f)",
          digits, pct,
          digits, ci_low,
          digits, ci_high)
}

# =============================================================================
# Function: Calculate Risk Difference
# =============================================================================
calculate_risk_difference <- function(cum_inc_results) {
  trt1 <- cum_inc_results[["1"]]
  trt0 <- cum_inc_results[["0"]]

  if (is.na(trt1$cum_inc) || is.na(trt0$cum_inc)) {
    return(list(rd = NA, rd_ci_lower = NA, rd_ci_upper = NA, rd_se = NA))
  }

  rd <- trt1$cum_inc - trt0$cum_inc
  rd_se <- sqrt(trt1$se^2 + trt0$se^2)
  rd_ci_lower <- rd - 1.96 * rd_se
  rd_ci_upper <- rd + 1.96 * rd_se

  return(list(
    rd = rd,
    rd_ci_lower = rd_ci_lower,
    rd_ci_upper = rd_ci_upper,
    rd_se = rd_se
  ))
}

# =============================================================================
# Main Analysis
# =============================================================================

cat("\n=============================================================================\n")
cat("4-Year Cumulative Incidence Analysis\n")
cat("=============================================================================\n\n")

# Storage for all results
all_results <- list()

for (analysis_name in names(rds_files)) {
  rds_path <- rds_files[[analysis_name]]

  cat(sprintf("\n--- Analysis: %s ---\n", toupper(analysis_name)))
  cat(sprintf("Loading: %s\n", rds_path))

  # Check if file exists
  if (!file.exists(rds_path)) {
    cat(sprintf("  [WARNING] File not found: %s\n", rds_path))
    next
  }

  # Load IPTW results
  ipwt_results <- readRDS(rds_path)

  cat(sprintf("  Found %d comparisons\n", length(ipwt_results)))

  analysis_results <- data.frame(
    Analysis = character(),
    Comparison = character(),
    Outcome = character(),
    N_Treatment = integer(),
    N_Comparator = integer(),
    Events_Treatment = integer(),
    Events_Comparator = integer(),
    CumInc_Treatment = character(),
    CumInc_Comparator = character(),
    Risk_Difference = character(),
    stringsAsFactors = FALSE
  )

  for (key in names(ipwt_results)) {
    cat(sprintf("\n  Processing: %s\n", key))

    result <- ipwt_results[[key]]

    # Check if cohort data exists
    if (is.null(result$cohort)) {
      cat("    [WARNING] No cohort data found\n")
      next
    }

    cohort_df <- result$cohort

    # Parse comparison name and outcome from key
    # Key format: "SEMAGLUTIDE vs OTHER_GLPA - Epilepsy/Seizure excl Hypoglycemic"
    parts <- strsplit(key, " - ")[[1]]
    comparison_name <- parts[1]
    outcome_name <- if (length(parts) > 1) parts[2] else "Unknown"

    # Get treatment labels
    treatment_labels <- switch(comparison_name,
      "SEMAGLUTIDE vs OTHER_GLPA" = c("Other GLP-1RA", "Semaglutide"),
      "SEMAGLUTIDE vs SGLT2" = c("SGLT2i", "Semaglutide"),
      "SEMAGLUTIDE vs OtherGLD" = c("Other GLDs", "Semaglutide"),
      "OTHER_GLPA vs SGLT2" = c("SGLT2i", "Other GLP-1RA"),
      "OTHER_GLPA vs OtherGLD" = c("Other GLDs", "Other GLP-1RA"),
      "SGLT2 vs OtherGLD" = c("Other GLDs", "SGLT2i"),
      c("Comparator", "Treatment")
    )

    # Calculate cumulative incidence
    cum_inc <- calculate_cumulative_incidence(cohort_df, time_months = TARGET_MONTHS)

    # Calculate risk difference
    rd <- calculate_risk_difference(cum_inc)

    # Format results
    trt_label <- treatment_labels[2]  # Treatment (coded as 1)
    comp_label <- treatment_labels[1]  # Comparator (coded as 0)

    # Format risk difference
    if (!is.na(rd$rd)) {
      rd_pct <- rd$rd * 100
      rd_ci_low <- rd$rd_ci_lower * 100
      rd_ci_high <- rd$rd_ci_upper * 100
      rd_str <- sprintf("%.2f%% (%.2f to %.2f)", rd_pct, rd_ci_low, rd_ci_high)
    } else {
      rd_str <- "NA"
    }

    # Add to results table
    analysis_results <- rbind(analysis_results, data.frame(
      Analysis = analysis_name,
      Comparison = comparison_name,
      Outcome = outcome_name,
      N_Treatment = cum_inc[["1"]]$n,
      N_Comparator = cum_inc[["0"]]$n,
      Events_Treatment = cum_inc[["1"]]$events,
      Events_Comparator = cum_inc[["0"]]$events,
      CumInc_Treatment = format_cum_inc(cum_inc[["1"]]),
      CumInc_Comparator = format_cum_inc(cum_inc[["0"]]),
      Risk_Difference = rd_str,
      stringsAsFactors = FALSE
    ))

    # Print summary
    cat(sprintf("    %s: N=%d, Events=%d, 4-yr CumInc=%s\n",
                trt_label, cum_inc[["1"]]$n, cum_inc[["1"]]$events,
                format_cum_inc(cum_inc[["1"]])))
    cat(sprintf("    %s: N=%d, Events=%d, 4-yr CumInc=%s\n",
                comp_label, cum_inc[["0"]]$n, cum_inc[["0"]]$events,
                format_cum_inc(cum_inc[["0"]])))
    cat(sprintf("    Risk Difference: %s\n", rd_str))
  }

  all_results[[analysis_name]] <- analysis_results
}

# =============================================================================
# Combine and Export Results
# =============================================================================

cat("\n\n=============================================================================\n")
cat("SUMMARY TABLE: 4-Year Cumulative Incidence\n")
cat("=============================================================================\n\n")

# Combine all results
final_results <- do.call(rbind, all_results)

if (nrow(final_results) > 0) {
  # Print formatted table
  print(final_results, row.names = FALSE)

  # Save to CSV
  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  csv_filename <- sprintf("4year_cumulative_incidence_%s.csv", timestamp)
  write.csv(final_results, file = csv_filename, row.names = FALSE)
  cat(sprintf("\n\nResults saved to: %s\n", csv_filename))

  # Also create a formatted version for publication
  pub_table <- final_results %>%
    mutate(
      `Treatment (N)` = sprintf("%s (N=%d)",
                                 sapply(Comparison, function(x) {
                                   switch(x,
                                     "SEMAGLUTIDE vs OTHER_GLPA" = "Semaglutide",
                                     "SEMAGLUTIDE vs SGLT2" = "Semaglutide",
                                     "SEMAGLUTIDE vs OtherGLD" = "Semaglutide",
                                     "OTHER_GLPA vs SGLT2" = "Other GLP-1RA",
                                     "OTHER_GLPA vs OtherGLD" = "Other GLP-1RA",
                                     "SGLT2 vs OtherGLD" = "SGLT2i",
                                     "Treatment"
                                   )
                                 }), N_Treatment),
      `Comparator (N)` = sprintf("%s (N=%d)",
                                  sapply(Comparison, function(x) {
                                    switch(x,
                                      "SEMAGLUTIDE vs OTHER_GLPA" = "Other GLP-1RA",
                                      "SEMAGLUTIDE vs SGLT2" = "SGLT2i",
                                      "SEMAGLUTIDE vs OtherGLD" = "Other GLDs",
                                      "OTHER_GLPA vs SGLT2" = "SGLT2i",
                                      "OTHER_GLPA vs OtherGLD" = "Other GLDs",
                                      "SGLT2 vs OtherGLD" = "Other GLDs",
                                      "Comparator"
                                    )
                                  }), N_Comparator),
      `4-Year Cumulative Incidence (Treatment)` = CumInc_Treatment,
      `4-Year Cumulative Incidence (Comparator)` = CumInc_Comparator,
      `Risk Difference (95% CI)` = Risk_Difference
    ) %>%
    select(Analysis, Comparison,
           `Treatment (N)`, `Comparator (N)`,
           `4-Year Cumulative Incidence (Treatment)`,
           `4-Year Cumulative Incidence (Comparator)`,
           `Risk Difference (95% CI)`)

  pub_csv_filename <- sprintf("4year_cumulative_incidence_publication_%s.csv", timestamp)
  write.csv(pub_table, file = pub_csv_filename, row.names = FALSE)
  cat(sprintf("Publication-ready table saved to: %s\n", pub_csv_filename))

} else {
  cat("No results to display. Please check RDS file paths.\n")
}

cat("\n=== Done ===\n")
