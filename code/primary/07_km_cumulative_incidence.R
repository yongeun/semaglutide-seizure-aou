# =============================================================================
# rev00015_KM_plot_semaglutide.R
# Kaplan-Meier (Cumulative Incidence) Plots for Semaglutide Comparisons
# Base R version (no survminer required)
# =============================================================================

# Load required packages
suppressPackageStartupMessages({
  library(survival)
  library(dplyr)
})

# =============================================================================
# Configuration
# =============================================================================

# RDS file paths
rds_files <- list(
  all_ages = "epilepsy_seizure_all_ages_ipwt_full_results.rds",
  late_onset = "epilepsy_seizure_late_onset_ipwt_full_results.rds"
)

# Comparisons to plot
target_comparisons <- c("SEMAGLUTIDE vs OtherGLD", "SEMAGLUTIDE vs SGLT2")

# Output directory
output_dir <- "."

# JAMA color palette
jama_colors <- c("#00274C", "#C4622D")  # Navy (Comparator), Orange (Semaglutide)

# =============================================================================
# Function: Create IPTW-weighted Cumulative Incidence Plot (Base R - Single)
# =============================================================================
plot_cumulative_incidence <- function(cohort_df,
                                       comparison_name,
                                       analysis_name,
                                       weight_var = "ipw_std",
                                       colors = jama_colors,
                                       show_title = FALSE) {

  # Treatment labels
  treatment_labels <- switch(comparison_name,
    "SEMAGLUTIDE vs OTHER_GLPA" = c("Other GLP-1RAs", "Semaglutide"),
    "SEMAGLUTIDE vs SGLT2" = c("SGLT2i", "Semaglutide"),
    "SEMAGLUTIDE vs OtherGLD" = c("Other GLDs", "Semaglutide"),
    c("Comparator", "Treatment")
  )

  # Convert event_time to months
  if (!"event_time_months" %in% names(cohort_df)) {
    cohort_df$event_time_months <- cohort_df$event_time / 30.44
  }

  # Check weight variable
  if (!weight_var %in% names(cohort_df)) {
    cohort_df[[weight_var]] <- 1
  }

  # Fit weighted survival curves by treatment
  fit <- survfit(
    Surv(event_time_months, event) ~ treatment,
    data = cohort_df,
    weights = cohort_df[[weight_var]]
  )

  # Get number at risk at specific times
  time_points <- c(0, 12, 24, 36, 48)
  summ <- summary(fit, times = time_points, extend = TRUE)

  # Extract n.risk for each strata
  strata_names <- names(fit$strata)
  n_risk_list <- list()

  for (s in strata_names) {
    idx <- which(summ$strata == s)
    n_risk_list[[s]] <- summ$n.risk[idx]
    # Pad if needed
    if (length(n_risk_list[[s]]) < length(time_points)) {
      n_risk_list[[s]] <- c(n_risk_list[[s]],
                             rep(NA, length(time_points) - length(n_risk_list[[s]])))
    }
  }

  # Extract survival curves for each group
  get_surv_data <- function(fit, strata_name) {
    idx <- which(names(fit$strata) == strata_name)
    if (idx == 1) {
      start_idx <- 1
    } else {
      start_idx <- sum(fit$strata[1:(idx-1)]) + 1
    }
    end_idx <- sum(fit$strata[1:idx])

    time <- c(0, fit$time[start_idx:end_idx])
    surv <- c(1, fit$surv[start_idx:end_idx])
    upper <- c(1, fit$upper[start_idx:end_idx])
    lower <- c(1, fit$lower[start_idx:end_idx])

    # Cumulative incidence
    list(
      time = time,
      cum_inc = 1 - surv,
      ci_lower = 1 - upper,
      ci_upper = 1 - lower
    )
  }

  surv_data_0 <- get_surv_data(fit, "treatment=0")
  surv_data_1 <- get_surv_data(fit, "treatment=1")

  # Determine y-axis max
  y_max <- max(c(surv_data_0$ci_upper, surv_data_1$ci_upper), na.rm = TRUE) * 1.2
  y_max <- max(y_max, 0.05)  # At least 5%
  y_max <- ceiling(y_max * 100) / 100  # Round up to nearest percent

  # Set up plot layout: main plot + risk table
  layout(matrix(c(1, 2), nrow = 2), heights = c(3.5, 1))

  # === Main Plot ===
  par(
    mar = c(1, 5, 2, 2),
    family = "sans",
    cex.axis = 0.95,
    cex.lab = 1.1,
    las = 1,
    mgp = c(3, 0.7, 0)
  )

  # Empty plot
  plot(NULL,
       xlim = c(0, 48),
       ylim = c(0, y_max),
       xlab = "",
       ylab = "Cumulative Incidence (%)",
       xaxt = "n",
       yaxt = "n",
       frame.plot = FALSE)

  # Y-axis with percentage labels
  y_ticks <- pretty(c(0, y_max), n = 5)
  y_ticks <- y_ticks[y_ticks <= y_max]
  axis(2, at = y_ticks, labels = paste0(y_ticks * 100, "%"), lwd = 0.8)

  # Grid lines
  abline(h = y_ticks, col = "gray85", lty = 1, lwd = 0.5)
  abline(v = seq(0, 48, 12), col = "gray85", lty = 1, lwd = 0.5)

  # Confidence interval polygons (step function style)
  # Group 0 (Comparator)
  n0 <- length(surv_data_0$time)
  x0_poly <- c(surv_data_0$time, rev(surv_data_0$time))
  y0_lower <- surv_data_0$ci_lower
  y0_upper <- surv_data_0$ci_upper

  # Step function for CI
  x0_step <- rep(surv_data_0$time, each = 2)[-1]
  x0_step <- c(x0_step, tail(surv_data_0$time, 1))
  y0_lower_step <- rep(y0_lower, each = 2)
  y0_lower_step <- y0_lower_step[-length(y0_lower_step)]
  y0_upper_step <- rep(y0_upper, each = 2)
  y0_upper_step <- y0_upper_step[-length(y0_upper_step)]

  polygon(c(x0_step, rev(x0_step)),
          c(y0_lower_step, rev(y0_upper_step)),
          col = adjustcolor(colors[1], alpha.f = 0.2),
          border = NA)

  # Group 1 (Semaglutide)
  x1_step <- rep(surv_data_1$time, each = 2)[-1]
  x1_step <- c(x1_step, tail(surv_data_1$time, 1))
  y1_lower <- surv_data_1$ci_lower
  y1_upper <- surv_data_1$ci_upper
  y1_lower_step <- rep(y1_lower, each = 2)
  y1_lower_step <- y1_lower_step[-length(y1_lower_step)]
  y1_upper_step <- rep(y1_upper, each = 2)
  y1_upper_step <- y1_upper_step[-length(y1_upper_step)]

  polygon(c(x1_step, rev(x1_step)),
          c(y1_lower_step, rev(y1_upper_step)),
          col = adjustcolor(colors[2], alpha.f = 0.2),
          border = NA)

  # Step function lines
  # Group 0
  x0_line <- rep(surv_data_0$time, each = 2)[-1]
  x0_line <- c(x0_line, tail(surv_data_0$time, 1))
  y0_line <- rep(surv_data_0$cum_inc, each = 2)
  y0_line <- y0_line[-length(y0_line)]
  lines(x0_line, y0_line, col = colors[1], lwd = 2.5)

  # Group 1
  x1_line <- rep(surv_data_1$time, each = 2)[-1]
  x1_line <- c(x1_line, tail(surv_data_1$time, 1))
  y1_line <- rep(surv_data_1$cum_inc, each = 2)
  y1_line <- y1_line[-length(y1_line)]
  lines(x1_line, y1_line, col = colors[2], lwd = 2.5)

  # Box
  box(lwd = 0.8)

  # Legend
  legend("topleft",
         legend = rev(treatment_labels),  # Semaglutide first
         col = rev(colors),
         lwd = 2.5,
         bty = "n",
         cex = 1.0,
         seg.len = 2)

  # Title if requested
  if (show_title) {
    title(main = sprintf("%s (%s)", comparison_name, analysis_name), cex.main = 1.1)
  }

  # === Risk Table ===
  par(mar = c(4, 5, 0.5, 2))

  plot(NULL,
       xlim = c(0, 48),
       ylim = c(0.3, 2.7),
       xlab = "Follow-up Time (months)",
       ylab = "",
       xaxt = "n",
       yaxt = "n",
       frame.plot = FALSE)

  axis(1, at = seq(0, 48, 12), lwd = 0.8)

  # "No. at risk" label
  mtext("No. at risk", side = 2, line = 3.5, cex = 0.9, las = 0)

  # Group labels and numbers
  text(-4, 2, treatment_labels[2], adj = 1, cex = 0.95, col = colors[2], font = 2, xpd = TRUE)
  text(-4, 1, treatment_labels[1], adj = 1, cex = 0.95, col = colors[1], font = 2, xpd = TRUE)

  n_risk_1 <- n_risk_list[["treatment=1"]]
  n_risk_0 <- n_risk_list[["treatment=0"]]

  for (i in seq_along(time_points)) {
    # Semaglutide (treatment=1)
    val1 <- if (i <= length(n_risk_1) && !is.na(n_risk_1[i])) format(round(n_risk_1[i]), big.mark = ",") else ""
    text(time_points[i], 2, val1, cex = 0.9, col = colors[2])

    # Comparator (treatment=0)
    val0 <- if (i <= length(n_risk_0) && !is.na(n_risk_0[i])) format(round(n_risk_0[i]), big.mark = ",") else ""
    text(time_points[i], 1, val0, cex = 0.9, col = colors[1])
  }
}

# =============================================================================
# Function: Save Single Plot
# =============================================================================
save_single_plot <- function(cohort_df, comparison_name, analysis_name,
                              weight_var = "ipw_std", colors = jama_colors) {

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
  clean_comp <- gsub(" ", "_", comparison_name)

  # PNG
  png_file <- sprintf("%s/KM_%s_%s_%s.png", output_dir, analysis_name, clean_comp, timestamp)
  png(png_file, width = 7, height = 6, units = "in", res = 300)
  plot_cumulative_incidence(cohort_df, comparison_name, analysis_name, weight_var, colors)
  dev.off()
  cat(sprintf("Saved: %s\n", png_file))

  # PDF
  pdf_file <- sprintf("%s/KM_%s_%s_%s.pdf", output_dir, analysis_name, clean_comp, timestamp)
  pdf(pdf_file, width = 7, height = 6)
  plot_cumulative_incidence(cohort_df, comparison_name, analysis_name, weight_var, colors)
  dev.off()
  cat(sprintf("Saved: %s\n", pdf_file))
}

# =============================================================================
# Function: Save Combined Panel (2 plots side by side)
# =============================================================================
save_combined_panel <- function(cohort_list, comparison_names, analysis_name,
                                 weight_var = "ipw_std", colors = jama_colors) {

  timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # Treatment labels for both comparisons
  labels_list <- list(
    "SEMAGLUTIDE vs OtherGLD" = c("Other GLDs", "Semaglutide"),
    "SEMAGLUTIDE vs SGLT2" = c("SGLT2i", "Semaglutide")
  )

  # PNG
  png_file <- sprintf("%s/KM_%s_combined_%s.png", output_dir, analysis_name, timestamp)
  png(png_file, width = 12, height = 6, units = "in", res = 300)

  # Layout: 2 main plots on top, 2 risk tables on bottom
  layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE), heights = c(3.5, 1))

  for (i in seq_along(cohort_list)) {
    comp_name <- comparison_names[i]
    cohort_df <- cohort_list[[i]]
    treatment_labels <- labels_list[[comp_name]]

    if (!"event_time_months" %in% names(cohort_df)) {
      cohort_df$event_time_months <- cohort_df$event_time / 30.44
    }
    if (!weight_var %in% names(cohort_df)) {
      cohort_df[[weight_var]] <- 1
    }

    fit <- survfit(Surv(event_time_months, event) ~ treatment, data = cohort_df,
                   weights = cohort_df[[weight_var]])

    # Get survival data
    get_surv_data <- function(fit, strata_name) {
      idx <- which(names(fit$strata) == strata_name)
      start_idx <- if (idx == 1) 1 else sum(fit$strata[1:(idx-1)]) + 1
      end_idx <- sum(fit$strata[1:idx])
      time <- c(0, fit$time[start_idx:end_idx])
      surv <- c(1, fit$surv[start_idx:end_idx])
      upper <- c(1, fit$upper[start_idx:end_idx])
      lower <- c(1, fit$lower[start_idx:end_idx])
      list(time = time, cum_inc = 1 - surv, ci_lower = 1 - upper, ci_upper = 1 - lower)
    }

    surv_data_0 <- get_surv_data(fit, "treatment=0")
    surv_data_1 <- get_surv_data(fit, "treatment=1")

    y_max <- max(c(surv_data_0$ci_upper, surv_data_1$ci_upper), na.rm = TRUE) * 1.2
    y_max <- max(y_max, 0.05)
    y_max <- ceiling(y_max * 100) / 100

    # Main plot
    par(mar = c(1, 5, 3, 1), cex.axis = 0.9, cex.lab = 1.0, las = 1)

    plot(NULL, xlim = c(0, 48), ylim = c(0, y_max),
         xlab = "", ylab = "Cumulative Incidence (%)", xaxt = "n", yaxt = "n", frame.plot = FALSE)

    y_ticks <- pretty(c(0, y_max), n = 5)
    y_ticks <- y_ticks[y_ticks <= y_max]
    axis(2, at = y_ticks, labels = paste0(y_ticks * 100, "%"), lwd = 0.8)

    abline(h = y_ticks, col = "gray85", lwd = 0.5)
    abline(v = seq(0, 48, 12), col = "gray85", lwd = 0.5)

    # CI polygons and lines
    # Group 0
    x0_step <- rep(surv_data_0$time, each = 2)[-1]
    x0_step <- c(x0_step, tail(surv_data_0$time, 1))
    y0_lower_step <- rep(surv_data_0$ci_lower, each = 2); y0_lower_step <- y0_lower_step[-length(y0_lower_step)]
    y0_upper_step <- rep(surv_data_0$ci_upper, each = 2); y0_upper_step <- y0_upper_step[-length(y0_upper_step)]
    polygon(c(x0_step, rev(x0_step)), c(y0_lower_step, rev(y0_upper_step)),
            col = adjustcolor(colors[1], alpha.f = 0.2), border = NA)
    y0_line <- rep(surv_data_0$cum_inc, each = 2); y0_line <- y0_line[-length(y0_line)]
    lines(x0_step, y0_line, col = colors[1], lwd = 2.5)

    # Group 1
    x1_step <- rep(surv_data_1$time, each = 2)[-1]
    x1_step <- c(x1_step, tail(surv_data_1$time, 1))
    y1_lower_step <- rep(surv_data_1$ci_lower, each = 2); y1_lower_step <- y1_lower_step[-length(y1_lower_step)]
    y1_upper_step <- rep(surv_data_1$ci_upper, each = 2); y1_upper_step <- y1_upper_step[-length(y1_upper_step)]
    polygon(c(x1_step, rev(x1_step)), c(y1_lower_step, rev(y1_upper_step)),
            col = adjustcolor(colors[2], alpha.f = 0.2), border = NA)
    y1_line <- rep(surv_data_1$cum_inc, each = 2); y1_line <- y1_line[-length(y1_line)]
    lines(x1_step, y1_line, col = colors[2], lwd = 2.5)

    box(lwd = 0.8)

    legend("topleft", legend = rev(treatment_labels), col = rev(colors), lwd = 2.5, bty = "n", cex = 0.9)

    # Panel label
    mtext(LETTERS[i], side = 3, line = 1, adj = 0, font = 2, cex = 1.3)
  }

  # Risk tables
  for (i in seq_along(cohort_list)) {
    comp_name <- comparison_names[i]
    cohort_df <- cohort_list[[i]]
    treatment_labels <- labels_list[[comp_name]]

    if (!"event_time_months" %in% names(cohort_df)) {
      cohort_df$event_time_months <- cohort_df$event_time / 30.44
    }
    if (!weight_var %in% names(cohort_df)) cohort_df[[weight_var]] <- 1

    fit <- survfit(Surv(event_time_months, event) ~ treatment, data = cohort_df,
                   weights = cohort_df[[weight_var]])

    time_points <- c(0, 12, 24, 36, 48)
    summ <- summary(fit, times = time_points, extend = TRUE)

    n_risk_0 <- summ$n.risk[summ$strata == "treatment=0"]
    n_risk_1 <- summ$n.risk[summ$strata == "treatment=1"]
    if (length(n_risk_0) < 5) n_risk_0 <- c(n_risk_0, rep(NA, 5 - length(n_risk_0)))
    if (length(n_risk_1) < 5) n_risk_1 <- c(n_risk_1, rep(NA, 5 - length(n_risk_1)))

    par(mar = c(4, 5, 0.5, 1))
    plot(NULL, xlim = c(0, 48), ylim = c(0.3, 2.7), xlab = "Follow-up (months)",
         ylab = "", xaxt = "n", yaxt = "n", frame.plot = FALSE)
    axis(1, at = seq(0, 48, 12), lwd = 0.8)
    mtext("No. at risk", side = 2, line = 3.5, cex = 0.8, las = 0)

    text(-4, 2, treatment_labels[2], adj = 1, cex = 0.85, col = colors[2], font = 2, xpd = TRUE)
    text(-4, 1, treatment_labels[1], adj = 1, cex = 0.85, col = colors[1], font = 2, xpd = TRUE)

    for (j in seq_along(time_points)) {
      val1 <- if (!is.na(n_risk_1[j])) format(round(n_risk_1[j]), big.mark = ",") else ""
      val0 <- if (!is.na(n_risk_0[j])) format(round(n_risk_0[j]), big.mark = ",") else ""
      text(time_points[j], 2, val1, cex = 0.8, col = colors[2])
      text(time_points[j], 1, val0, cex = 0.8, col = colors[1])
    }
  }

  dev.off()
  cat(sprintf("Saved: %s\n", png_file))

  # PDF
  pdf_file <- sprintf("%s/KM_%s_combined_%s.pdf", output_dir, analysis_name, timestamp)
  pdf(pdf_file, width = 12, height = 6)
  # Repeat same plotting code for PDF...
  layout(matrix(c(1, 2, 3, 4), nrow = 2, byrow = TRUE), heights = c(3.5, 1))

  for (i in seq_along(cohort_list)) {
    comp_name <- comparison_names[i]
    cohort_df <- cohort_list[[i]]
    treatment_labels <- labels_list[[comp_name]]

    if (!"event_time_months" %in% names(cohort_df)) cohort_df$event_time_months <- cohort_df$event_time / 30.44
    if (!weight_var %in% names(cohort_df)) cohort_df[[weight_var]] <- 1

    fit <- survfit(Surv(event_time_months, event) ~ treatment, data = cohort_df, weights = cohort_df[[weight_var]])

    get_surv_data <- function(fit, strata_name) {
      idx <- which(names(fit$strata) == strata_name)
      start_idx <- if (idx == 1) 1 else sum(fit$strata[1:(idx-1)]) + 1
      end_idx <- sum(fit$strata[1:idx])
      time <- c(0, fit$time[start_idx:end_idx])
      surv <- c(1, fit$surv[start_idx:end_idx])
      upper <- c(1, fit$upper[start_idx:end_idx])
      lower <- c(1, fit$lower[start_idx:end_idx])
      list(time = time, cum_inc = 1 - surv, ci_lower = 1 - upper, ci_upper = 1 - lower)
    }

    surv_data_0 <- get_surv_data(fit, "treatment=0")
    surv_data_1 <- get_surv_data(fit, "treatment=1")

    y_max <- max(c(surv_data_0$ci_upper, surv_data_1$ci_upper), na.rm = TRUE) * 1.2
    y_max <- max(y_max, 0.05)
    y_max <- ceiling(y_max * 100) / 100

    par(mar = c(1, 5, 3, 1), cex.axis = 0.9, cex.lab = 1.0, las = 1)
    plot(NULL, xlim = c(0, 48), ylim = c(0, y_max), xlab = "", ylab = "Cumulative Incidence (%)", xaxt = "n", yaxt = "n", frame.plot = FALSE)
    y_ticks <- pretty(c(0, y_max), n = 5); y_ticks <- y_ticks[y_ticks <= y_max]
    axis(2, at = y_ticks, labels = paste0(y_ticks * 100, "%"), lwd = 0.8)
    abline(h = y_ticks, col = "gray85", lwd = 0.5); abline(v = seq(0, 48, 12), col = "gray85", lwd = 0.5)

    x0_step <- rep(surv_data_0$time, each = 2)[-1]; x0_step <- c(x0_step, tail(surv_data_0$time, 1))
    y0_lower_step <- rep(surv_data_0$ci_lower, each = 2); y0_lower_step <- y0_lower_step[-length(y0_lower_step)]
    y0_upper_step <- rep(surv_data_0$ci_upper, each = 2); y0_upper_step <- y0_upper_step[-length(y0_upper_step)]
    polygon(c(x0_step, rev(x0_step)), c(y0_lower_step, rev(y0_upper_step)), col = adjustcolor(colors[1], alpha.f = 0.2), border = NA)
    y0_line <- rep(surv_data_0$cum_inc, each = 2); y0_line <- y0_line[-length(y0_line)]
    lines(x0_step, y0_line, col = colors[1], lwd = 2.5)

    x1_step <- rep(surv_data_1$time, each = 2)[-1]; x1_step <- c(x1_step, tail(surv_data_1$time, 1))
    y1_lower_step <- rep(surv_data_1$ci_lower, each = 2); y1_lower_step <- y1_lower_step[-length(y1_lower_step)]
    y1_upper_step <- rep(surv_data_1$ci_upper, each = 2); y1_upper_step <- y1_upper_step[-length(y1_upper_step)]
    polygon(c(x1_step, rev(x1_step)), c(y1_lower_step, rev(y1_upper_step)), col = adjustcolor(colors[2], alpha.f = 0.2), border = NA)
    y1_line <- rep(surv_data_1$cum_inc, each = 2); y1_line <- y1_line[-length(y1_line)]
    lines(x1_step, y1_line, col = colors[2], lwd = 2.5)

    box(lwd = 0.8)
    legend("topleft", legend = rev(treatment_labels), col = rev(colors), lwd = 2.5, bty = "n", cex = 0.9)
    mtext(LETTERS[i], side = 3, line = 1, adj = 0, font = 2, cex = 1.3)
  }

  for (i in seq_along(cohort_list)) {
    comp_name <- comparison_names[i]
    cohort_df <- cohort_list[[i]]
    treatment_labels <- labels_list[[comp_name]]

    if (!"event_time_months" %in% names(cohort_df)) cohort_df$event_time_months <- cohort_df$event_time / 30.44
    if (!weight_var %in% names(cohort_df)) cohort_df[[weight_var]] <- 1

    fit <- survfit(Surv(event_time_months, event) ~ treatment, data = cohort_df, weights = cohort_df[[weight_var]])
    time_points <- c(0, 12, 24, 36, 48)
    summ <- summary(fit, times = time_points, extend = TRUE)
    n_risk_0 <- summ$n.risk[summ$strata == "treatment=0"]
    n_risk_1 <- summ$n.risk[summ$strata == "treatment=1"]
    if (length(n_risk_0) < 5) n_risk_0 <- c(n_risk_0, rep(NA, 5 - length(n_risk_0)))
    if (length(n_risk_1) < 5) n_risk_1 <- c(n_risk_1, rep(NA, 5 - length(n_risk_1)))

    par(mar = c(4, 5, 0.5, 1))
    plot(NULL, xlim = c(0, 48), ylim = c(0.3, 2.7), xlab = "Follow-up (months)", ylab = "", xaxt = "n", yaxt = "n", frame.plot = FALSE)
    axis(1, at = seq(0, 48, 12), lwd = 0.8)
    mtext("No. at risk", side = 2, line = 3.5, cex = 0.8, las = 0)

    text(-4, 2, treatment_labels[2], adj = 1, cex = 0.85, col = colors[2], font = 2, xpd = TRUE)
    text(-4, 1, treatment_labels[1], adj = 1, cex = 0.85, col = colors[1], font = 2, xpd = TRUE)

    for (j in seq_along(time_points)) {
      val1 <- if (!is.na(n_risk_1[j])) format(round(n_risk_1[j]), big.mark = ",") else ""
      val0 <- if (!is.na(n_risk_0[j])) format(round(n_risk_0[j]), big.mark = ",") else ""
      text(time_points[j], 2, val1, cex = 0.8, col = colors[2])
      text(time_points[j], 1, val0, cex = 0.8, col = colors[1])
    }
  }

  dev.off()
  cat(sprintf("Saved: %s\n", pdf_file))
}

# =============================================================================
# Main Analysis
# =============================================================================

cat("\n=============================================================================\n")
cat("Kaplan-Meier Plots: Semaglutide Comparisons (Base R)\n")
cat("=============================================================================\n\n")

for (analysis_name in names(rds_files)) {
  rds_path <- rds_files[[analysis_name]]

  cat(sprintf("\n=== %s ===\n", toupper(analysis_name)))

  if (!file.exists(rds_path)) {
    cat(sprintf("  [WARNING] File not found: %s\n", rds_path))
    next
  }

  ipwt_results <- readRDS(rds_path)
  cat(sprintf("  Loaded: %s\n", rds_path))

  cohort_list <- list()
  comparison_found <- c()

  for (target_comp in target_comparisons) {
    matching_keys <- grep(target_comp, names(ipwt_results), value = TRUE)

    if (length(matching_keys) == 0) {
      cat(sprintf("  [WARNING] Comparison not found: %s\n", target_comp))
      next
    }

    key <- matching_keys[1]
    cat(sprintf("\n  Processing: %s\n", key))

    result <- ipwt_results[[key]]
    if (is.null(result$cohort)) { cat("    [WARNING] No cohort data\n"); next }

    cohort_df <- result$cohort

    n_trt <- sum(cohort_df$treatment == 1)
    n_comp <- sum(cohort_df$treatment == 0)
    events_trt <- sum(cohort_df$event[cohort_df$treatment == 1])
    events_comp <- sum(cohort_df$event[cohort_df$treatment == 0])

    cat(sprintf("    Semaglutide: N=%d, Events=%d\n", n_trt, events_trt))
    cat(sprintf("    Comparator: N=%d, Events=%d\n", n_comp, events_comp))

    save_single_plot(cohort_df, target_comp, analysis_name, "ipw_std", jama_colors)

    cohort_list[[target_comp]] <- cohort_df
    comparison_found <- c(comparison_found, target_comp)
  }

  if (length(cohort_list) == 2) {
    cat("\n  Creating combined panel...\n")
    save_combined_panel(cohort_list, comparison_found, analysis_name, "ipw_std", jama_colors)
  }
}

cat("\n=============================================================================\n")
cat("Done!\n")
cat("=============================================================================\n")
