# =============================================================================
# rev00012_suppl_figure6_from_rds.R.r
# Supplementary Figure 6: Proportion Mediated Over Time
# Uses existing mediation analysis RDS file
# =============================================================================
#
# INPUT: mediation_analysis_results_*.rds (from rev00006_mediation_final.R.r)
#
# OUTPUT:
#   - Supplementary_Figure_6_Proportion_Mediated.pdf
#   - Supplementary_Figure_6_Proportion_Mediated.png
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
  input_rds = "mediation_analysis_results_20260223_051609.rds",
  output_dir = "Mediation results"
)

# Create output directory
dir.create(config$output_dir, showWarnings = FALSE, recursive = TRUE)

# >>> CELL 03: Helper Function - Extract Proportion Mediated
.getS_at <- function(curves, lbl, t) {
  x <- curves %>%
    mutate(lbl = case_when(
      grepl("^S_11", type) ~ "S11",
      grepl("^S_10", type) ~ "S10",
      grepl("^S_00", type) ~ "S00",
      TRUE ~ NA_character_
    )) %>%
    filter(lbl == !!lbl, t_months <= !!t) %>%
    arrange(t_months)
  if (nrow(x) == 0) return(NA_real_)
  x$S[nrow(x)]
}

pm_from_S <- function(curves, clip0 = TRUE) {
  grid <- sort(unique(curves$t_months))

  S11 <- vapply(grid, function(tt) .getS_at(curves, "S11", tt), numeric(1))
  S10 <- vapply(grid, function(tt) .getS_at(curves, "S10", tt), numeric(1))
  S00 <- vapply(grid, function(tt) .getS_at(curves, "S00", tt), numeric(1))

  denom <- S11 - S00
  pm_raw <- (S11 - S10) / denom
  pm <- ifelse(is.finite(denom) & abs(denom) > 1e-12, pm_raw, NA_real_)
  if (clip0) pm <- pmax(0, pm)

  tibble(t_months = grid, PM_S = pm)
}

# >>> CELL 04: Load RDS
cat("\n========== LOADING MEDIATION RESULTS ==========\n")

mediation_results <- readRDS(config$input_rds)

cat("Loaded mediation results.\n")
cat("Available components:\n")
for (name in names(mediation_results)) {
  cat(sprintf("  - %s\n", name))
}

# >>> CELL 05: Extract Proportion Mediated Over Time
cat("\n========== EXTRACTING PROPORTION MEDIATED ==========\n")

pm_time_data <- bind_rows(
  pm_from_S(mediation_results$res_othergld_a1c$curves) %>%
    mutate(Comparison = "Semaglutide vs Other GLDs", Mediator = "HbA1c", Panel = "A"),
  pm_from_S(mediation_results$res_sglt2_a1c$curves) %>%
    mutate(Comparison = "Semaglutide vs SGLT2i", Mediator = "HbA1c", Panel = "B"),
  pm_from_S(mediation_results$res_othergld_bmi$curves) %>%
    mutate(Comparison = "Semaglutide vs Other GLDs", Mediator = "BMI", Panel = "C"),
  pm_from_S(mediation_results$res_sglt2_bmi$curves) %>%
    mutate(Comparison = "Semaglutide vs SGLT2i", Mediator = "BMI", Panel = "D")
) %>%
  mutate(
    PM_pct = PM_S * 100,
    Panel_label = paste0("(", Panel, ") ", Mediator, ": ", Comparison)
  )

# Print summary
cat("\nProportion Mediated Summary:\n")
pm_time_data %>%
  group_by(Panel, Comparison, Mediator) %>%
  summarise(
    Mean_PM = mean(PM_pct, na.rm = TRUE),
    Max_PM = max(PM_pct, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  print()

# >>> CELL 06: Create Figure
cat("\n========== CREATING SUPPLEMENTARY FIGURE 6 ==========\n")

create_pm_panel <- function(data, panel_letter) {
  panel_data <- data %>% filter(Panel == panel_letter)
  title_text <- unique(panel_data$Panel_label)

  # Calculate y-axis limit
  y_max <- max(15, max(panel_data$PM_pct, na.rm = TRUE) + 2)

  ggplot(panel_data, aes(x = t_months, y = PM_pct)) +
    geom_line(linewidth = 1.2, color = "#2E86AB") +
    geom_point(size = 2, color = "#2E86AB") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_x_continuous(
      breaks = c(0, 12, 24, 36, 48),
      limits = c(0, 48),
      expand = c(0.02, 0)
    ) +
    scale_y_continuous(
      limits = c(-2, y_max),
      expand = c(0, 0)
    ) +
    labs(
      title = title_text,
      x = "Months since index date",
      y = "Proportion mediated (%)"
    ) +
    theme_bw(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 10),
      panel.grid.minor = element_blank(),
      axis.title = element_text(size = 10),
      axis.text = element_text(size = 9)
    )
}

# Generate 4 panels
panel_A <- create_pm_panel(pm_time_data, "A")
panel_B <- create_pm_panel(pm_time_data, "B")
panel_C <- create_pm_panel(pm_time_data, "C")
panel_D <- create_pm_panel(pm_time_data, "D")

# Combine into 2x2 grid
suppl_fig6 <- plot_grid(
  panel_A, panel_B, panel_C, panel_D,
  ncol = 2, nrow = 2,
  labels = c("A", "B", "C", "D"),
  label_size = 14,
  label_fontface = "bold"
)

# Add main title
title_grob <- ggdraw() +
  draw_label(
    "Supplementary Figure 6. Proportion Mediated by HbA1c and BMI",
    fontface = "bold",
    size = 12,
    x = 0.5,
    hjust = 0.5
  )

final_fig <- plot_grid(
  title_grob,
  suppl_fig6,
  ncol = 1,
  rel_heights = c(0.05, 0.95)
)

# >>> CELL 07: Save Figure
cat("\n========== SAVING FIGURE ==========\n")

timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

pdf_path <- file.path(config$output_dir, paste0("Supplementary_Figure_6_Proportion_Mediated_", timestamp, ".pdf"))
png_path <- file.path(config$output_dir, paste0("Supplementary_Figure_6_Proportion_Mediated_", timestamp, ".png"))

ggsave(pdf_path, plot = final_fig, width = 10, height = 8, dpi = 300)
cat(sprintf("Saved: %s\n", pdf_path))

ggsave(png_path, plot = final_fig, width = 10, height = 8, dpi = 300)
cat(sprintf("Saved: %s\n", png_path))

# >>> CELL 08: Summary
cat("\n")
cat("============================================================\n")
cat("SUPPLEMENTARY FIGURE 6 GENERATION COMPLETE\n")
cat("============================================================\n")
cat("Output directory:", config$output_dir, "\n\n")
cat("FILES SAVED:\n")
cat(sprintf("  - Supplementary_Figure_6_Proportion_Mediated_%s.pdf\n", timestamp))
cat(sprintf("  - Supplementary_Figure_6_Proportion_Mediated_%s.png\n", timestamp))
cat("\nFigure description:\n")
cat("  (A) HbA1c mediation: Semaglutide vs Other GLDs\n")
cat("  (B) HbA1c mediation: Semaglutide vs SGLT2i\n")
cat("  (C) BMI mediation: Semaglutide vs Other GLDs\n")
cat("  (D) BMI mediation: Semaglutide vs SGLT2i\n")
cat("============================================================\n")
