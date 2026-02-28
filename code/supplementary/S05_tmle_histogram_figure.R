# =============================================================================
# rev00011_tmle_histogram_figure.R.r
# TMLE Late-Onset Results - Histogram Comparison Figure
# =============================================================================
#
# INPUT: outcome4_epilepsy_excl_hypo_late_onset_tmle_full_results.rds
#
# OUTPUT:
#   - tmle_figures/tmle_late_onset_histogram_comparison.pdf
#   - tmle_figures/tmle_late_onset_histogram_comparison.png
#   - tmle_figures/tmle_late_onset_semaglutide_vs_othergld.pdf/png
#   - tmle_figures/tmle_late_onset_semaglutide_vs_sglt2.pdf/png
#   - tmle_figures/tmle_late_onset_summary.csv
# =============================================================================

# >>> CELL 01: Load Libraries
suppressPackageStartupMessages({
  library(tidyverse)
  library(ggplot2)
  library(cowplot)
})

cat("Libraries loaded successfully.\n")

# >>> CELL 02: Configuration
config <- list(
  # Use main TMLE results RDS (contains SEMAGLUTIDE comparisons with EY0/EY1)
  input_rds = "outcome4_epilepsy_excl_hypo_late_onset_tmle_full_results.rds",
  # Fallback to CSV if RDS not available
  input_csv = "outcome4_epilepsy_excl_hypo_late_onset_tmle_comprehensive_results_20260220_222048.csv",
  output_dir = "tmle_figures",

  # Color scheme
  colors = list(
    control = "#1f77b4",    # Blue
    treatment = "#ff7f0e"   # Orange
  ),

  # Only include these comparisons
  target_comparisons = c("SEMAGLUTIDE vs OtherGLD", "SEMAGLUTIDE vs SGLT2"),

  # Comparison labels for display
  comparison_labels = list(
    "SEMAGLUTIDE vs OtherGLD" = "Semaglutide vs Other GLD",
    "SEMAGLUTIDE vs SGLT2" = "Semaglutide vs SGLT2i"
  )
)

# Create output directory
dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)

# >>> CELL 03: Load TMLE Results from RDS
cat("\n========== LOADING TMLE RESULTS ==========\n")

# Try RDS first, fallback to CSV
use_rds <- file.exists(config$input_rds)

if (use_rds) {
  cat(sprintf("Loading RDS file: %s\n", config$input_rds))
  tmle_results <- readRDS(config$input_rds)

  cat(sprintf("Loaded TMLE results with %d analyses\n", length(tmle_results)))
  cat("Available analyses:\n")
  for (name in names(tmle_results)) {
    cat(sprintf("  - %s\n", name))
  }
} else {
  cat("RDS file not found, using CSV fallback\n")
}

# >>> CELL 04: Extract Probability Estimates from RDS
cat("\n========== EXTRACTING PROBABILITY ESTIMATES ==========\n")

extract_probabilities <- function(tmle_fit, comp_name) {
  tryCatch({
    # Extract EY0 (control) and EY1 (treatment) estimates
    ey0_psi <- tmle_fit$fit$estimates$EY0$psi
    ey1_psi <- tmle_fit$fit$estimates$EY1$psi
    ey0_ci <- tmle_fit$fit$estimates$EY0$CI
    ey1_ci <- tmle_fit$fit$estimates$EY1$CI

    # Risk difference and p-value
    ate <- tmle_fit$fit$estimates$ATE
    rd <- ate$psi
    rd_ci <- ate$CI
    p_value <- ate$pvalue

    # Get sample sizes from cohort
    cohort <- tmle_fit$cohort
    n_treat <- sum(cohort$treatment == 1)
    n_ctrl <- sum(cohort$treatment == 0)
    events_treat <- sum(cohort$event[cohort$treatment == 1])
    events_ctrl <- sum(cohort$event[cohort$treatment == 0])

    data.frame(
      comparison = comp_name,
      ey0 = ey0_psi,
      ey0_lower = ey0_ci[1],
      ey0_upper = ey0_ci[2],
      ey1 = ey1_psi,
      ey1_lower = ey1_ci[1],
      ey1_upper = ey1_ci[2],
      rd = rd,
      rd_lower = rd_ci[1],
      rd_upper = rd_ci[2],
      p_val = p_value,
      n_treat = n_treat,
      n_ctrl = n_ctrl,
      events_treat = events_treat,
      events_ctrl = events_ctrl,
      stringsAsFactors = FALSE
    )
  }, error = function(e) {
    cat(sprintf("  Error extracting %s: %s\n", comp_name, e$message))
    return(NULL)
  })
}

# Extract from target comparisons
prob_data_list <- list()

for (analysis_name in names(tmle_results)) {
  # Parse comparison name from analysis name
  comp_name <- gsub(" - .*$", "", analysis_name)

  # Check if this comparison is in target list
  if (!comp_name %in% config$target_comparisons) {
    cat(sprintf("Skipping: %s\n", comp_name))
    next
  }

  cat(sprintf("Processing: %s\n", comp_name))
  tmle_fit <- tmle_results[[analysis_name]]

  probs <- extract_probabilities(tmle_fit, comp_name)
  if (!is.null(probs)) {
    prob_data_list[[comp_name]] <- probs
    cat(sprintf("  EY0 (Control): %.4f (%.4f - %.4f)\n",
                probs$ey0, probs$ey0_lower, probs$ey0_upper))
    cat(sprintf("  EY1 (Treatment): %.4f (%.4f - %.4f)\n",
                probs$ey1, probs$ey1_lower, probs$ey1_upper))
    cat(sprintf("  Risk Difference: %.4f (p = %.2e)\n", probs$rd, probs$p_val))
  }
}

# Combine all data
prob_data <- do.call(rbind, prob_data_list)
rownames(prob_data) <- NULL

cat(sprintf("\nExtracted probabilities for %d comparisons\n", nrow(prob_data)))

# >>> CELL 05: Create Individual Comparison Plots
cat("\n\n========== CREATING HISTOGRAM PLOTS ==========\n")

create_comparison_plot <- function(prob_row, colors, labels) {
  comp_name <- prob_row$comparison
  display_label <- ifelse(comp_name %in% names(labels), labels[[comp_name]], comp_name)

  # Parse treatment and control names
  treatment_name <- "Semaglutide"
  control_name <- ifelse(grepl("OtherGLD", comp_name), "Other GLD", "SGLT2i")

  # Create data for plotting (use TMLE EY0/EY1 estimates)
  plot_df <- data.frame(
    Group = factor(c(control_name, treatment_name), levels = c(control_name, treatment_name)),
    Probability = c(prob_row$ey0, prob_row$ey1),
    Lower = c(prob_row$ey0_lower, prob_row$ey1_lower),
    Upper = c(prob_row$ey0_upper, prob_row$ey1_upper),
    Events = c(prob_row$events_ctrl, prob_row$events_treat),
    N = c(prob_row$n_ctrl, prob_row$n_treat),
    Type = c("Control", "Treatment")
  )

  # Format p-value
  p_text <- ifelse(prob_row$p_val < 0.001, "p < 0.001",
                   sprintf("p = %.3f", prob_row$p_val))

  # Create plot
  p <- ggplot(plot_df, aes(x = Group, y = Probability * 100, fill = Type)) +
    geom_bar(stat = "identity", alpha = 0.85, width = 0.6) +
    geom_errorbar(aes(ymin = Lower * 100, ymax = Upper * 100),
                  width = 0.15, linewidth = 0.8, color = "gray30") +
    geom_text(aes(label = sprintf("%d/%d", Events, N)),
              vjust = -0.8, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = c("Control" = colors$control, "Treatment" = colors$treatment)) +
    scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
    labs(
      title = display_label,
      subtitle = sprintf("Risk Difference: %.2f%% (%.2f%% to %.2f%%), %s",
                         prob_row$rd * 100, prob_row$rd_lower * 100, prob_row$rd_upper * 100, p_text),
      y = "TMLE Estimated Probability (%)",
      x = ""
    ) +
    theme_minimal(base_size = 12) +
    theme(
      plot.title = element_text(face = "bold", size = 14, hjust = 0.5),
      plot.subtitle = element_text(size = 10, hjust = 0.5, color = "gray30"),
      legend.position = "none",
      panel.grid.major.x = element_blank(),
      axis.text.x = element_text(size = 12, face = "bold"),
      axis.title.y = element_text(size = 11)
    )

  return(p)
}

# Create individual plots
plot_list <- list()
for (i in 1:nrow(prob_data)) {
  plot_list[[i]] <- create_comparison_plot(
    prob_data[i, ],
    config$colors,
    config$comparison_labels
  )
}

# >>> CELL 06: Create Combined Figure
cat("\n========== CREATING COMBINED FIGURE ==========\n")

# Arrange plots side by side (2 plots)
n_plots <- length(plot_list)

# Create combined plot with title
combined_title <- ggdraw() +
  draw_label(
    "TMLE Analysis: Late-Onset Adult Seizure Disorder",
    fontface = "bold",
    size = 16,
    x = 0.5,
    hjust = 0.5
  )

# Combine all plots
plot_grid_obj <- plot_grid(
  plotlist = plot_list,
  ncol = 2,
  nrow = 1,
  labels = c("A", "B"),
  label_size = 14
)

# Final combined figure
final_plot <- plot_grid(
  combined_title,
  plot_grid_obj,
  ncol = 1,
  rel_heights = c(0.08, 0.92)
)

# >>> CELL 07: Save Figures
cat("\n========== SAVING FIGURES ==========\n")

# Save combined figure
pdf_path <- file.path(config$output_dir, "tmle_late_onset_histogram_comparison.pdf")
png_path <- file.path(config$output_dir, "tmle_late_onset_histogram_comparison.png")

ggsave(pdf_path, plot = final_plot, width = 12, height = 6, dpi = 300)
cat(sprintf("Saved: %s\n", pdf_path))

ggsave(png_path, plot = final_plot, width = 12, height = 6, dpi = 300)
cat(sprintf("Saved: %s\n", png_path))

# Save individual plots
for (i in 1:n_plots) {
  comp_name <- prob_data$comparison[i]
  safe_name <- tolower(gsub("[^A-Za-z0-9]+", "_", comp_name))

  ind_pdf <- file.path(config$output_dir, sprintf("tmle_late_onset_%s.pdf", safe_name))
  ind_png <- file.path(config$output_dir, sprintf("tmle_late_onset_%s.png", safe_name))

  ggsave(ind_pdf, plot = plot_list[[i]], width = 6, height = 5, dpi = 300)
  ggsave(ind_png, plot = plot_list[[i]], width = 6, height = 5, dpi = 300)
  cat(sprintf("Saved: %s\n", ind_pdf))
}

# >>> CELL 08: Create Summary Table
cat("\n========== SUMMARY TABLE ==========\n")

summary_table <- prob_data %>%
  mutate(
    Comparison = comparison,
    `Control (EY0)` = sprintf("%.2f%% (%.2f-%.2f)", ey0*100, ey0_lower*100, ey0_upper*100),
    `Treatment (EY1)` = sprintf("%.2f%% (%.2f-%.2f)", ey1*100, ey1_lower*100, ey1_upper*100),
    `Events (Treat/Ctrl)` = sprintf("%d/%d vs %d/%d", events_treat, n_treat, events_ctrl, n_ctrl),
    `Risk Difference (%)` = sprintf("%.2f (%.2f to %.2f)", rd*100, rd_lower*100, rd_upper*100),
    `P-value` = ifelse(p_val < 0.001, "<0.001", sprintf("%.2e", p_val))
  ) %>%
  select(Comparison, `Control (EY0)`, `Treatment (EY1)`, `Events (Treat/Ctrl)`, `Risk Difference (%)`, `P-value`)

print(summary_table)

# Save summary table
summary_path <- file.path(config$output_dir, "tmle_late_onset_summary.csv")
write.csv(summary_table, summary_path, row.names = FALSE)
cat(sprintf("\nSaved summary: %s\n", summary_path))

# >>> CELL 09: Print Completion Message
cat("\n")
cat("============================================================\n")
cat("TMLE HISTOGRAM FIGURE GENERATION COMPLETE\n")
cat("============================================================\n")
cat("Output directory:", config$output_dir, "\n\n")
cat("FILES SAVED:\n")
cat("  - tmle_late_onset_histogram_comparison.pdf (combined figure)\n")
cat("  - tmle_late_onset_histogram_comparison.png (combined figure)\n")
for (i in 1:n_plots) {
  comp_name <- prob_data$comparison[i]
  safe_name <- tolower(gsub("[^A-Za-z0-9]+", "_", comp_name))
  cat(sprintf("  - tmle_late_onset_%s.pdf/png (individual)\n", safe_name))
}
cat("  - tmle_late_onset_summary.csv\n")
cat("============================================================\n")
