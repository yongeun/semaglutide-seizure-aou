# =============================================================================
# rev00013_redraw_density_no_legend.R
# Redraw Bootstrap Density Plots WITHOUT Legend (Original Style)
# =============================================================================

# =============================================================================
# Configuration
# =============================================================================
output_dir <- "."  # Change if needed

# JAMA colors
jama_colors <- c("#00274C", "#C4622D")

# =============================================================================
# Function: Draw Density Plot Without Legend (Original Base R Style)
# =============================================================================
draw_density_no_legend <- function(rds_path,
                                   output_prefix,
                                   outcome = "Epilepsy/Seizure excl Hypoglycemic",
                                   colors = jama_colors,
                                   use_jama_style = TRUE) {

  # Load RDS
  if (!file.exists(rds_path)) {
    stop(sprintf("File not found: %s", rds_path))
  }

  tmle_results <- readRDS(rds_path)
  cat(sprintf("Loaded: %s\n", rds_path))
  cat(sprintf("Keys: %s\n", paste(names(tmle_results), collapse = ", ")))

  # Extract bootstrap data
  all_bootstrap_data <- list()

  for (key in names(tmle_results)) {
    if (!is.null(tmle_results[[key]]$bootstrap) &&
        !is.null(tmle_results[[key]]$bootstrap$samples) &&
        length(tmle_results[[key]]$bootstrap$samples) > 0) {

      comp_name <- gsub(" - Epilepsy.*", "", key)
      all_bootstrap_data[[comp_name]] <- tmle_results[[key]]$bootstrap
      cat(sprintf("  %s: %d samples\n", comp_name, length(tmle_results[[key]]$bootstrap$samples)))
    }
  }

  if (length(all_bootstrap_data) == 0) {
    warning("No bootstrap samples found in RDS")
    return(NULL)
  }

  # Calculate ranges
  all_samples <- unlist(lapply(all_bootstrap_data, function(x) x$samples))
  x_range <- range(c(all_samples, 0)) * 1.1
  densities <- lapply(all_bootstrap_data, function(x) density(x$samples))
  y_max <- max(sapply(densities, function(d) max(d$y))) * 1.1

  # Output files
  ts <- format(Sys.time(), "%Y%m%d_%H%M%S")

  # --- PNG Output ---
  png_file <- sprintf("%s/%s_density_no_legend_%s.png", output_dir, output_prefix, ts)
  grDevices::png(png_file, width = 1200, height = 800, res = 150)

  if (use_jama_style) {
    par(
      family = "sans",
      cex.main = 1.1,
      cex.lab = 1.0,
      cex.axis = 0.9,
      font.main = 1,
      las = 1,
      mgp = c(2.5, 0.7, 0),
      tcl = -0.3,
      mar = c(5, 5, 4, 2)
    )
  }

  plot(NULL,
       xlim = x_range,
       ylim = c(0, y_max),
       main = paste("Bootstrap Distributions:", outcome),
       xlab = "Risk Difference",
       ylab = "Density",
       type = "n",
       axes = FALSE,
       frame.plot = FALSE)

  axis(1, col = "black", col.axis = "black", lwd = 0.5)
  axis(2, col = "black", col.axis = "black", lwd = 0.5)

  if (use_jama_style) {
    abline(h = axTicks(2), col = "gray90", lty = 1, lwd = 0.5)
    abline(v = axTicks(1), col = "gray90", lty = 1, lwd = 0.5)
  }

  box(lwd = 0.5)
  abline(v = 0, col = "black", lwd = 1.5, lty = 1)

  # Draw density curves
  for (i in seq_along(all_bootstrap_data)) {
    dens <- densities[[i]]
    col_i <- colors[(i - 1) %% length(colors) + 1]
    lines(dens, col = col_i, lwd = 2)
    polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0, length(dens$y))),
            col = adjustcolor(col_i, alpha.f = 0.15), border = NA)
    abline(v = all_bootstrap_data[[i]]$mean, col = col_i, lwd = 1, lty = 2)
  }

  # NO LEGEND - removed

  grDevices::dev.off()
  cat(sprintf("Saved PNG: %s\n", png_file))

  # --- PDF Output ---
  pdf_file <- sprintf("%s/%s_density_no_legend_%s.pdf", output_dir, output_prefix, ts)
  grDevices::pdf(pdf_file, width = 10, height = 7)

  if (use_jama_style) {
    par(
      family = "sans",
      cex.main = 1.1,
      cex.lab = 1.0,
      cex.axis = 0.9,
      font.main = 1,
      las = 1,
      mgp = c(2.5, 0.7, 0),
      tcl = -0.3,
      mar = c(5, 5, 4, 2)
    )
  }

  plot(NULL,
       xlim = x_range,
       ylim = c(0, y_max),
       main = paste("Bootstrap Distributions:", outcome),
       xlab = "Risk Difference",
       ylab = "Density",
       type = "n",
       axes = FALSE,
       frame.plot = FALSE)

  axis(1, col = "black", col.axis = "black", lwd = 0.5)
  axis(2, col = "black", col.axis = "black", lwd = 0.5)

  if (use_jama_style) {
    abline(h = axTicks(2), col = "gray90", lty = 1, lwd = 0.5)
    abline(v = axTicks(1), col = "gray90", lty = 1, lwd = 0.5)
  }

  box(lwd = 0.5)
  abline(v = 0, col = "black", lwd = 1.5, lty = 1)

  for (i in seq_along(all_bootstrap_data)) {
    dens <- densities[[i]]
    col_i <- colors[(i - 1) %% length(colors) + 1]
    lines(dens, col = col_i, lwd = 2)
    polygon(c(dens$x, rev(dens$x)), c(dens$y, rep(0, length(dens$y))),
            col = adjustcolor(col_i, alpha.f = 0.15), border = NA)
    abline(v = all_bootstrap_data[[i]]$mean, col = col_i, lwd = 1, lty = 2)
  }

  # NO LEGEND - removed

  grDevices::dev.off()
  cat(sprintf("Saved PDF: %s\n", pdf_file))

  return(list(png = png_file, pdf = pdf_file, data = all_bootstrap_data))
}

# =============================================================================
# Generate Plots
# =============================================================================

cat("\n=== Generating Density Plots WITHOUT Legend (Original Style) ===\n\n")

# All Ages
cat("--- All Ages ---\n")
result_all_ages <- draw_density_no_legend(
  rds_path = "outcome4_epilepsy_excl_hypo_all_ages_tmle_full_results.rds",
  output_prefix = "all_ages",
  outcome = "Epilepsy/Seizure excl Hypoglycemic (All Ages)",
  colors = c("#00274C", "#C4622D"),
  use_jama_style = TRUE
)

# Late Onset
cat("\n--- Late Onset ---\n")
result_late_onset <- draw_density_no_legend(
  rds_path = "outcome4_epilepsy_excl_hypo_late_onset_tmle_full_results.rds",
  output_prefix = "late_onset",
  outcome = "Epilepsy/Seizure excl Hypoglycemic (Late Onset)",
  colors = c("#00274C", "#C4622D"),
  use_jama_style = TRUE
)

cat("\n=== Done ===\n")
