# >>> CELL 01: ================================================================
# Outcome 4: Epilepsy/Seizure excl Hypoglycemic - FINAL PAPER VERSION
# Design: 12-month windows up to 48 months (0-12-24-36-48m)
# Counterfactual mediation (Vansteelandt-style), single mediator.
#
# OUTPUT: 논문에 필요한 Figure와 Table만 저장
#   - Figure: Counterfactual cumulative incidence plots (PDF)
#   - Table: Proportion mediated summary & HR table (CSV)
# ================================================================
suppressPackageStartupMessages({
  library(dplyr); library(readr); library(survival); library(ggplot2); library(tidyr)
})

# >>> CELL 02: ----- helpers ---------------------------------------------------
to_chr <- function(x) if (is.character(x)) x else as.character(x)

H0_at <- function(bh_df, p, t) {
  x <- bh_df %>% dplyr::filter(period == p)
  if (nrow(x) == 0) return(0)
  id <- which(x$time <= t)
  if (length(id) == 0) return(0)
  x$hazard[max(id)]
}

make_time_grid_48 <- function(bh) {
  g1 <- bh %>% filter(period==1, time>0,    time<=365)  %>% pull(time)
  g2 <- bh %>% filter(period==2, time>365,  time<=730)  %>% pull(time)
  g3 <- bh %>% filter(period==3, time>730,  time<=1095) %>% pull(time)
  g4 <- bh %>% filter(period==4, time>1095, time<=1460) %>% pull(time)
  sort(unique(c(0, g1, 365, g2, 730, g3, 1095, g4, 1460)))
}

# >>> CELL 03: ----- 1) Build baseline X and window M1..M4 (mean or last) -----
build_X_M_48 <- function(df0, panel_path,
                         date_col  = "date_of_measurement",
                         value_col = "A1c",
                         baseline_col = NULL,
                         locf = TRUE,
                         window_method = c("mean","last")) {

  method <- match.arg(window_method)

  df <- df0 %>%
    mutate(person_id = to_chr(person_id),
           index_date = as.Date(index_date),
           event_time = as.numeric(event_time),
           event      = as.integer(event),
           treatment  = as.integer(treatment))

  pan_raw <- readr::read_csv(panel_path, show_col_types = FALSE) %>%
    mutate(person_id = to_chr(person_id))
  stopifnot(all(c("person_id", date_col, value_col) %in% names(pan_raw)))

  pan <- pan_raw %>%
    transmute(person_id,
              .date = as.Date(.data[[date_col]]),
              .val  = as.numeric(.data[[value_col]])) %>%
    filter(is.finite(.val))

  long <- df %>%
    select(person_id, index_date, event_time) %>%
    inner_join(pan, by="person_id") %>%
    mutate(days_since = as.numeric(.date - index_date)) %>%
    filter(is.finite(days_since))

  if (!is.null(baseline_col) && baseline_col %in% names(df)) {
    X_df <- df %>% transmute(person_id, X = as.numeric(.data[[baseline_col]]))
  } else {
    X_df <- long %>%
      filter(days_since <= 0) %>%
      group_by(person_id) %>%
      slice_max(order_by = days_since, n = 1, with_ties = FALSE) %>%
      ungroup() %>%
      transmute(person_id, X = .val)
  }

  win_agg <- function(start, end, method = c("mean","last")) {
    method <- match.arg(method)
    tmp <- long %>%
      mutate(bound = pmin(event_time, end)) %>%
      filter(days_since > start, days_since <= bound)

    if (method == "mean") {
      tmp %>%
        group_by(person_id) %>%
        summarise(m = mean(.val, na.rm = TRUE), .groups = "drop")
    } else {
      tmp %>%
        group_by(person_id) %>%
        slice_max(order_by = days_since, n = 1, with_ties = FALSE) %>%
        ungroup() %>%
        transmute(person_id, m = .val)
    }
  }

  M1 <- win_agg(   0,  365, method) %>% rename(M1 = m)
  M2 <- win_agg( 365,  730, method) %>% rename(M2 = m)
  M3 <- win_agg( 730, 1095, method) %>% rename(M3 = m)
  M4 <- win_agg(1095, 1460, method) %>% rename(M4 = m)

  out <- df %>%
    left_join(X_df, by="person_id") %>%
    left_join(M1, by="person_id") %>%
    left_join(M2, by="person_id") %>%
    left_join(M3, by="person_id") %>%
    left_join(M4, by="person_id")

  if (locf) {
    out <- out %>%
      mutate(
        M1 = ifelse(is.na(M1), X,  M1),
        M2 = ifelse(is.na(M2), M1, M2),
        M3 = ifelse(is.na(M3), M2, M3),
        M4 = ifelse(is.na(M4), M3, M4)
      )
  }

  out %>%
    mutate(.etime = event_time, .e = event,
           event_time = pmin(event_time, 1460),
           event      = ifelse(.e == 1 & .etime <= 1460, 1L, 0L))
}

# >>> CELL 04: ----- 2) Start-stop & attach piecewise mediator ------------------
build_startstop_48 <- function(dfM) {
  survSplit(
    data   = dfM,
    cut    = c(365, 730, 1095),
    end    = "event_time",
    start  = "tstart",
    event  = "event",
    episode= "period"
  ) %>%
    filter(tstart < 1460) %>%
    mutate(period   = as.integer(period),
           M_period = case_when(
             period==1L ~ M1,
             period==2L ~ M2,
             period==3L ~ M3,
             period==4L ~ M4,
             TRUE ~ NA_real_
           ))
}

# >>> CELL 05: ----- 3) Cox: Surv ~ A + M_period + X + strata(period) ----------
fit_cox_48 <- function(ss) {
  coxph(Surv(tstart, event_time, event) ~ treatment + M_period + X +
          strata(period) + cluster(person_id),
        data = ss, ties = "breslow", robust = TRUE)
}

get_bh_by_period <- function(fit) {
  bh <- basehaz(fit, centered = FALSE)
  bh %>%
    mutate(period = as.integer(sub("period=", "", strata, fixed = TRUE))) %>%
    select(period, time, hazard) %>%
    arrange(period, time)
}

# >>> CELL 06: ----- 4) Mediator models & chained simulation --------------------
fit_med_models_chain <- function(dfM) {
  list(
    m1 = lm(M1 ~ treatment + X,  data = dfM),
    m2 = lm(M2 ~ treatment + M1, data = dfM),
    m3 = lm(M3 ~ treatment + M2, data = dfM),
    m4 = lm(M4 ~ treatment + M3, data = dfM)
  )
}

sim_paths_chain <- function(dfM, med_models, a_value, B = 200, seed = 123) {
  set.seed(seed)
  N <- nrow(dfM)
  s1 <- summary(med_models$m1)$sigma; if (!is.finite(s1)) s1 <- 0
  s2 <- summary(med_models$m2)$sigma; if (!is.finite(s2)) s2 <- 0
  s3 <- summary(med_models$m3)$sigma; if (!is.finite(s3)) s3 <- 0
  s4 <- summary(med_models$m4)$sigma; if (!is.finite(s4)) s4 <- 0

  arr <- array(NA_real_, dim = c(N, 4, B),
               dimnames = list(NULL, c("M1","M2","M3","M4"), NULL))

  for (b in 1:B) {
    mu1 <- as.numeric(predict(med_models$m1,
                              newdata = data.frame(treatment = a_value, X = dfM$X)))
    M1  <- mu1 + rnorm(N, 0, s1)

    mu2 <- as.numeric(predict(med_models$m2,
                              newdata = data.frame(treatment = a_value, M1 = M1)))
    M2  <- mu2 + rnorm(N, 0, s2)

    mu3 <- as.numeric(predict(med_models$m3,
                              newdata = data.frame(treatment = a_value, M2 = M2)))
    M3  <- mu3 + rnorm(N, 0, s3)

    mu4 <- as.numeric(predict(med_models$m4,
                              newdata = data.frame(treatment = a_value, M3 = M3)))
    M4  <- mu4 + rnorm(N, 0, s4)

    arr[, , b] <- cbind(M1, M2, M3, M4)
  }
  arr
}

# >>> CELL 07: ----- 5) Stepwise survival prediction (0-48m, 4 windows) --------
predict_survival_step_48 <- function(dfM, fit, bh, A_fixed, M_paths, times_days) {
  beta <- coef(fit)
  bA <- unname(beta["treatment"])
  bM <- unname(beta["M_period"])
  bX <- unname(beta["X"])

  N <- nrow(dfM); B <- dim(M_paths)[3]
  S_by_b <- matrix(NA_real_, nrow = B, ncol = length(times_days))

  for (b in 1:B) {
    lp1 <- bA*A_fixed + bM*M_paths[,1,b] + bX*dfM$X
    lp2 <- bA*A_fixed + bM*M_paths[,2,b] + bX*dfM$X
    lp3 <- bA*A_fixed + bM*M_paths[,3,b] + bX*dfM$X
    lp4 <- bA*A_fixed + bM*M_paths[,4,b] + bX*dfM$X

    S_it <- matrix(NA_real_, nrow=N, ncol=length(times_days))
    for (k in seq_along(times_days)) {
      t <- times_days[k]
      if (t <= 365) {
        H1 <- H0_at(bh, 1, t) - H0_at(bh, 1, 0)
        S_it[,k] <- exp(- H1*exp(lp1))
      } else if (t <= 730) {
        H1 <- H0_at(bh,1,365)-H0_at(bh,1,0)
        H2 <- H0_at(bh,2,t)  -H0_at(bh,2,365)
        S_it[,k] <- exp(- H1*exp(lp1) - H2*exp(lp2))
      } else if (t <= 1095) {
        H1 <- H0_at(bh,1,365)-H0_at(bh,1,0)
        H2 <- H0_at(bh,2,730)-H0_at(bh,2,365)
        H3 <- H0_at(bh,3,t)  -H0_at(bh,3,730)
        S_it[,k] <- exp(- H1*exp(lp1) - H2*exp(lp2) - H3*exp(lp3))
      } else {
        H1 <- H0_at(bh,1,365) -H0_at(bh,1,0)
        H2 <- H0_at(bh,2,730) -H0_at(bh,2,365)
        H3 <- H0_at(bh,3,1095)-H0_at(bh,3,730)
        H4 <- H0_at(bh,4,t)   -H0_at(bh,4,1095)
        S_it[,k] <- exp(- H1*exp(lp1) - H2*exp(lp2) - H3*exp(lp3) - H4*exp(lp4))
      }
    }
    S_by_b[b,] <- colMeans(S_it, na.rm=TRUE)
  }

  S_hat <- colMeans(S_by_b, na.rm=TRUE)
  tibble::tibble(
    t_days   = times_days,
    t_months = round(times_days/30.44),
    S        = pmin(pmax(S_hat,0),1)
  )
}

# >>> CELL 08: ----- 6) End-to-end runner (SIMPLIFIED - no file saving) ----------
run_cf_mediation_48m <- function(df0,
                                 panel_path,
                                 date_col  = "date_of_measurement",
                                 value_col = "A1c",
                                 baseline_col = NULL,
                                 B = 200, seed = 123, locf = TRUE,
                                 window_method = c("mean","last")) {

  dfM <- build_X_M_48(df0, panel_path, date_col, value_col,
                      baseline_col, locf, window_method = match.arg(window_method))

  ss <- build_startstop_48(dfM)
  cat("Events by period & treatment:\n")
  ss %>% dplyr::group_by(period, treatment) %>%
    dplyr::summarise(n=dplyr::n(), events=sum(event), .groups="drop") %>% print(n=Inf)

  fit <- fit_cox_48(ss)
  bh  <- get_bh_by_period(fit)
  grid <- make_time_grid_48(bh)
  if (length(grid) < 5) grid <- c(0,365,730,1095,1460)

  mods <- fit_med_models_chain(dfM)
  M_A0 <- sim_paths_chain(dfM, mods, a_value = 0, B = B, seed = seed)
  M_A1 <- sim_paths_chain(dfM, mods, a_value = 1, B = B, seed = seed + 1)

  S11 <- predict_survival_step_48(dfM, fit, bh, A_fixed = 1, M_paths = M_A1, times_days = grid) %>%
    dplyr::mutate(type = "S_11 (A=1, M from A=1)")
  S10 <- predict_survival_step_48(dfM, fit, bh, A_fixed = 1, M_paths = M_A0, times_days = grid) %>%
    dplyr::mutate(type = "S_10 (A=1, M from A=0)")
  S00 <- predict_survival_step_48(dfM, fit, bh, A_fixed = 0, M_paths = M_A0, times_days = grid) %>%
    dplyr::mutate(type = "S_00 (A=0, M from A=0)")

  curves <- dplyr::bind_rows(S11,S10,S00) %>%
    dplyr::mutate(cum_incidence = 1 - S)

  list(curves = curves, cox_fit = fit, basehaz = bh, dfM = dfM, ss = ss)
}

# >>> CELL 09: ----- Proportion mediated helpers --------------------------------
.getS_at <- function(curves, lbl, t){
  x <- curves |>
    dplyr::mutate(lbl = dplyr::case_when(
      grepl("^S_11", type) ~ "S11",
      grepl("^S_10", type) ~ "S10",
      grepl("^S_00", type) ~ "S00",
      TRUE ~ NA_character_
    )) |>
    dplyr::filter(lbl == !!lbl, t_months <= !!t) |>
    dplyr::arrange(t_months)
  if (nrow(x) == 0) return(NA_real_)
  x$S[nrow(x)]
}

pm_from_S <- function(curves, clip0 = TRUE) {
  grid <- sort(unique(curves$t_months))

  S11 <- vapply(grid, function(tt) .getS_at(curves, "S11", tt), numeric(1))
  S10 <- vapply(grid, function(tt) .getS_at(curves, "S10", tt), numeric(1))
  S00 <- vapply(grid, function(tt) .getS_at(curves, "S00", tt), numeric(1))

  denom <- S11 - S00
  pm <- (S11 - S10) / denom
  pm_raw <- (S11 - S10) / denom
  pm <- ifelse(is.finite(denom) & abs(denom) > 1e-12, pm_raw, NA_real_)
  if (clip0) pm <- pmax(0, pm)

  mean_pm <- mean(pm, na.rm = TRUE)
  max_pm  <- max(pm, na.rm = TRUE)
  idx_max <- which.max(pm)
  t_at_max <- if (length(idx_max) && is.finite(max_pm)) grid[idx_max] else NA_real_

  list(
    pm_by_time = tibble::tibble(t_months = grid, PM_S = pm),
    summary = tibble::tibble(
      mean_PM_S = mean_pm,
      max_PM_S  = max_pm,
      t_months_at_max = t_at_max
    )
  )
}

# >>> CELL 10: ---- Parametric bootstrap HR table ----
hr_table_parametric <- function(res,
                                t_months = c(12,24,36,48),
                                R = 300, B = 200, seed = 2025) {
  stopifnot(all(c("dfM","cox_fit","basehaz") %in% names(res)))
  dfM <- res$dfM; fit <- res$cox_fit
  set.seed(seed)
  if (!requireNamespace("MASS", quietly = TRUE))
    stop("Please install.packages('MASS') for mvrnorm().")

  mods <- fit_med_models_chain(dfM)
  per_of <- function(tm) if (tm<=12) 1L else if (tm<=24) 2L else if (tm<=36) 3L else 4L
  periods <- vapply(t_months, per_of, integer(1))

  beta <- coef(fit); V <- vcov(fit)
  need <- c("treatment","M_period","X")
  idx  <- match(need, names(beta))
  if (any(is.na(idx))) stop("Cox model must contain ", paste(need, collapse=", "))

  K <- length(t_months)
  HR_TE  <- HR_NDE <- HR_NIE <- matrix(NA_real_, nrow = R, ncol = K)

  mean_exp_lp <- function(bA, bM, bX, A_fixed, paths) {
    N <- nrow(dfM)
    out <- numeric(4)
    for (p in 1:4) {
      tmp <- sapply(1:B, function(bi) {
        lp <- bA*A_fixed + bM*paths[,p,bi] + bX*dfM$X
        mean(exp(lp), na.rm = TRUE)
      })
      out[p] <- mean(tmp, na.rm = TRUE)
    }
    out
  }

  for (r in 1:R) {
    b_draw <- MASS::mvrnorm(1, mu = beta[idx], Sigma = V[idx, idx, drop=FALSE])
    bA <- b_draw[1]; bM <- b_draw[2]; bX <- b_draw[3]

    M_A0 <- sim_paths_chain(dfM, mods, a_value = 0, B = B, seed = seed + 1000*r)
    M_A1 <- sim_paths_chain(dfM, mods, a_value = 1, B = B, seed = seed + 1000*r + 1)

    mA1M1 <- mean_exp_lp(bA,bM,bX, A_fixed=1, paths=M_A1)
    mA1M0 <- mean_exp_lp(bA,bM,bX, A_fixed=1, paths=M_A0)
    mA0M0 <- mean_exp_lp(bA,bM,bX, A_fixed=0, paths=M_A0)

    for (k in seq_len(K)) {
      p <- periods[k]
      HR_TE [r,k] <- mA1M1[p] / mA0M0[p]
      HR_NDE[r,k] <- mA1M0[p] / mA0M0[p]
      HR_NIE[r,k] <- HR_TE[r,k] / HR_NDE[r,k]
    }
  }

  summarize_cols <- function(M)
    apply(M, 2, function(v) c(est=mean(v, na.rm=TRUE),
                              lwr=unname(quantile(v, .025, na.rm=TRUE)),
                              upr=unname(quantile(v, .975, na.rm=TRUE))))
  s_TE  <- summarize_cols(HR_TE)
  s_NDE <- summarize_cols(HR_NDE)
  s_NIE <- summarize_cols(HR_NIE)

  tibble::tibble(
    t_months   = t_months,
    HR_TE      = s_TE ["est",], HR_TE_lwr  = s_TE ["lwr",], HR_TE_upr  = s_TE ["upr",],
    HR_NDE     = s_NDE["est",], HR_NDE_lwr = s_NDE["lwr",], HR_NDE_upr = s_NDE["upr",],
    HR_NIE     = s_NIE["est",], HR_NIE_lwr = s_NIE["lwr",], HR_NIE_upr = s_NIE["upr",]
  )
}

# >>> CELL 11: ============================================================
# RUN ALL ANALYSES
# ============================================================
cat("\n========== LOADING IPTW RESULTS ==========\n")
ipwt_results <- readRDS("outcome4_epilepsy_excl_hypo_all_ages_ipwt_full_results.rds")

cat("\n========== ANALYSIS 1: SEMAGLUTIDE vs SGLT2 x A1c ==========\n")
res_sglt2_a1c <- run_cf_mediation_48m(
  df0 = ipwt_results[[ "SEMAGLUTIDE vs SGLT2 - Epilepsy/Seizure excl Hypoglycemic" ]]$cohort,
  panel_path = "a1c_panel.csv",
  date_col   = "date_of_measurement",
  value_col  = "A1c",
  B = 100, seed = 2025,
  baseline_col="baseline_hba1c",
  window_method = "mean"
)

cat("\n========== ANALYSIS 2: SEMAGLUTIDE vs OtherGLD x A1c ==========\n")
res_othergld_a1c <- run_cf_mediation_48m(
  df0 = ipwt_results[[ "SEMAGLUTIDE vs OtherGLD - Epilepsy/Seizure excl Hypoglycemic" ]]$cohort,
  panel_path = "a1c_panel.csv",
  date_col   = "date_of_measurement",
  value_col  = "A1c",
  B = 100, seed = 2025,
  baseline_col="baseline_hba1c",
  window_method = "mean"
)

cat("\n========== ANALYSIS 3: SEMAGLUTIDE vs SGLT2 x BMI ==========\n")
res_sglt2_bmi <- run_cf_mediation_48m(
  df0 = ipwt_results[[ "SEMAGLUTIDE vs SGLT2 - Epilepsy/Seizure excl Hypoglycemic" ]]$cohort,
  panel_path = "bmi_panel.csv",
  date_col   = "weight_date",
  value_col  = "bmi",
  B = 100, seed = 2025,
  baseline_col="baseline_bmi",
  window_method = "mean"
)

cat("\n========== ANALYSIS 4: SEMAGLUTIDE vs OtherGLD x BMI ==========\n")
res_othergld_bmi <- run_cf_mediation_48m(
  df0 = ipwt_results[[ "SEMAGLUTIDE vs OtherGLD - Epilepsy/Seizure excl Hypoglycemic" ]]$cohort,
  panel_path = "bmi_panel.csv",
  date_col   = "weight_date",
  value_col  = "bmi",
  B = 100, seed = 2025,
  baseline_col="baseline_bmi",
  window_method = "mean"
)

# >>> CELL 12: ============================================================
# GENERATE FINAL OUTPUTS FOR PAPER
# ============================================================
timestamp <- format(Sys.time(), "%Y%m%d_%H%M%S")
output_dir <- "mediation_final_results"
dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)

# ---- 1. FIGURE: Counterfactual Cumulative Incidence Plots ----
create_mediation_figure <- function(res, title_text, filename) {
  curves <- res$curves %>%
    mutate(
      type_label = case_when(
        grepl("S_11", type) ~ "Semaglutide (actual)",
        grepl("S_10", type) ~ "Semaglutide (counterfactual)",
        grepl("S_00", type) ~ "Comparator",
        TRUE ~ type
      )
    )

  p <- ggplot(curves, aes(x = t_months, y = cum_incidence, color = type_label)) +
    geom_step(linewidth = 1) +
    scale_x_continuous(breaks = c(0, 12, 24, 36, 48), limits = c(0, 48)) +
    scale_y_continuous(labels = scales::percent_format(accuracy = 0.1)) +
    labs(
      title = title_text,
      x = "Months since index date",
      y = "Cumulative incidence",
      color = NULL
    ) +
    theme_bw(base_size = 12) +
    theme(
      legend.position = "bottom",
      panel.grid.minor = element_blank(),
      plot.title = element_text(face = "bold", size = 11)
    )

  if (requireNamespace("ggsci", quietly = TRUE)) {
    p <- p + ggsci::scale_color_jama()
  }

  ggsave(filename, plot = p, width = 7, height = 5, dpi = 300)
  return(p)
}

# Generate all 4 figures
fig1 <- create_mediation_figure(res_sglt2_a1c,
  "Semaglutide vs SGLT2i: HbA1c Mediation",
  file.path(output_dir, paste0("Figure_Mediation_SGLT2_A1c_", timestamp, ".pdf")))

fig2 <- create_mediation_figure(res_othergld_a1c,
  "Semaglutide vs Other GLDs: HbA1c Mediation",
  file.path(output_dir, paste0("Figure_Mediation_OtherGLD_A1c_", timestamp, ".pdf")))

fig3 <- create_mediation_figure(res_sglt2_bmi,
  "Semaglutide vs SGLT2i: BMI Mediation",
  file.path(output_dir, paste0("Figure_Mediation_SGLT2_BMI_", timestamp, ".pdf")))

fig4 <- create_mediation_figure(res_othergld_bmi,
  "Semaglutide vs Other GLDs: BMI Mediation",
  file.path(output_dir, paste0("Figure_Mediation_OtherGLD_BMI_", timestamp, ".pdf")))

# ---- 1b. SUPPLEMENTARY FIGURE 6: Proportion Mediated Over Time ----
cat("\n========== GENERATING SUPPLEMENTARY FIGURE 6 ==========\n")

# Extract pm_by_time for all 4 analyses
pm_time_data <- bind_rows(
  pm_from_S(res_othergld_a1c$curves)$pm_by_time %>%
    mutate(Comparison = "Semaglutide vs Other GLDs", Mediator = "HbA1c", Panel = "A"),
  pm_from_S(res_sglt2_a1c$curves)$pm_by_time %>%
    mutate(Comparison = "Semaglutide vs SGLT2i", Mediator = "HbA1c", Panel = "B"),
  pm_from_S(res_othergld_bmi$curves)$pm_by_time %>%
    mutate(Comparison = "Semaglutide vs Other GLDs", Mediator = "BMI", Panel = "C"),
  pm_from_S(res_sglt2_bmi$curves)$pm_by_time %>%
    mutate(Comparison = "Semaglutide vs SGLT2i", Mediator = "BMI", Panel = "D")
) %>%
  mutate(
    PM_pct = PM_S * 100,
    Panel_label = paste0("(", Panel, ") ", Mediator, ": ", Comparison)
  )

# Create individual panel plots
create_pm_panel <- function(data, panel_letter) {
  panel_data <- data %>% filter(Panel == panel_letter)
  title_text <- unique(panel_data$Panel_label)

  ggplot(panel_data, aes(x = t_months, y = PM_pct)) +
    geom_line(linewidth = 1, color = "#2E86AB") +
    geom_point(size = 1.5, color = "#2E86AB") +
    geom_hline(yintercept = 0, linetype = "dashed", color = "gray50") +
    scale_x_continuous(breaks = c(0, 12, 24, 36, 48), limits = c(0, 48)) +
    scale_y_continuous(limits = c(-5, max(15, max(panel_data$PM_pct, na.rm = TRUE) + 2))) +
    labs(
      title = title_text,
      x = "Months since index date",
      y = "Proportion mediated (%)"
    ) +
    theme_bw(base_size = 10) +
    theme(
      plot.title = element_text(face = "bold", size = 9),
      panel.grid.minor = element_blank()
    )
}

# Generate 4 panels
panel_A <- create_pm_panel(pm_time_data, "A")
panel_B <- create_pm_panel(pm_time_data, "B")
panel_C <- create_pm_panel(pm_time_data, "C")
panel_D <- create_pm_panel(pm_time_data, "D")

# Combine into 2x2 grid
if (requireNamespace("cowplot", quietly = TRUE)) {
  library(cowplot)
  suppl_fig6 <- plot_grid(
    panel_A, panel_B, panel_C, panel_D,
    ncol = 2, nrow = 2,
    labels = c("A", "B", "C", "D"),
    label_size = 12
  )
} else {
  if (requireNamespace("gridExtra", quietly = TRUE)) {
    library(gridExtra)
    suppl_fig6 <- grid.arrange(panel_A, panel_B, panel_C, panel_D, ncol = 2)
  }
}

# Save Supplementary Figure 6
ggsave(
  file.path(output_dir, paste0("Supplementary_Figure_6_Proportion_Mediated_", timestamp, ".pdf")),
  plot = suppl_fig6, width = 10, height = 8, dpi = 300
)
ggsave(
  file.path(output_dir, paste0("Supplementary_Figure_6_Proportion_Mediated_", timestamp, ".png")),
  plot = suppl_fig6, width = 10, height = 8, dpi = 300
)
cat("Supplementary Figure 6 saved.\n")

# ---- 2. TABLE: Proportion Mediated Summary ----
pm_summary <- bind_rows(
  pm_from_S(res_sglt2_a1c$curves)$summary %>% mutate(Comparison = "Semaglutide vs SGLT2i", Mediator = "HbA1c"),
  pm_from_S(res_othergld_a1c$curves)$summary %>% mutate(Comparison = "Semaglutide vs Other GLDs", Mediator = "HbA1c"),
  pm_from_S(res_sglt2_bmi$curves)$summary %>% mutate(Comparison = "Semaglutide vs SGLT2i", Mediator = "BMI"),
  pm_from_S(res_othergld_bmi$curves)$summary %>% mutate(Comparison = "Semaglutide vs Other GLDs", Mediator = "BMI")
) %>%
  select(Comparison, Mediator, mean_PM_S, max_PM_S, t_months_at_max) %>%
  mutate(
    mean_PM_S = round(mean_PM_S * 100, 1),
    max_PM_S = round(max_PM_S * 100, 1)
  ) %>%
  rename(
    `Mean Proportion Mediated (%)` = mean_PM_S,
    `Max Proportion Mediated (%)` = max_PM_S,
    `Time at Max (months)` = t_months_at_max
  )

write.csv(pm_summary,
  file.path(output_dir, paste0("Table_Proportion_Mediated_Summary_", timestamp, ".csv")),
  row.names = FALSE)

# ---- 3. TABLE: HR with 95% CI at key timepoints ----
cat("\n========== COMPUTING HR TABLES (this may take a few minutes) ==========\n")

hr_all <- bind_rows(
  hr_table_parametric(res_sglt2_a1c, t_months = c(12, 24, 36, 48), R = 300, B = 100, seed = 2025) %>%
    mutate(Comparison = "Semaglutide vs SGLT2i", Mediator = "HbA1c"),
  hr_table_parametric(res_othergld_a1c, t_months = c(12, 24, 36, 48), R = 300, B = 100, seed = 2025) %>%
    mutate(Comparison = "Semaglutide vs Other GLDs", Mediator = "HbA1c"),
  hr_table_parametric(res_sglt2_bmi, t_months = c(12, 24, 36, 48), R = 300, B = 100, seed = 2025) %>%
    mutate(Comparison = "Semaglutide vs SGLT2i", Mediator = "BMI"),
  hr_table_parametric(res_othergld_bmi, t_months = c(12, 24, 36, 48), R = 300, B = 100, seed = 2025) %>%
    mutate(Comparison = "Semaglutide vs Other GLDs", Mediator = "BMI")
) %>%
  select(Comparison, Mediator, t_months, everything()) %>%
  mutate(
    HR_TE_formatted = sprintf("%.2f (%.2f-%.2f)", HR_TE, HR_TE_lwr, HR_TE_upr),
    HR_NDE_formatted = sprintf("%.2f (%.2f-%.2f)", HR_NDE, HR_NDE_lwr, HR_NDE_upr),
    HR_NIE_formatted = sprintf("%.2f (%.2f-%.2f)", HR_NIE, HR_NIE_lwr, HR_NIE_upr)
  )

# Save full HR table
write.csv(hr_all,
  file.path(output_dir, paste0("Table_Mediation_HR_Full_", timestamp, ".csv")),
  row.names = FALSE)

# Save formatted HR table for paper
hr_formatted <- hr_all %>%
  select(Comparison, Mediator, t_months, HR_TE_formatted, HR_NDE_formatted, HR_NIE_formatted) %>%
  rename(
    `Time (months)` = t_months,
    `Total Effect HR (95% CI)` = HR_TE_formatted,
    `Natural Direct Effect HR (95% CI)` = HR_NDE_formatted,
    `Natural Indirect Effect HR (95% CI)` = HR_NIE_formatted
  )

write.csv(hr_formatted,
  file.path(output_dir, paste0("Table_Mediation_HR_Formatted_", timestamp, ".csv")),
  row.names = FALSE)

# >>> CELL 13: ============================================================
# SAVE RDS FOR FUTURE USE
# ============================================================
mediation_results <- list(
  res_sglt2_a1c = res_sglt2_a1c,
  res_othergld_a1c = res_othergld_a1c,
  res_sglt2_bmi = res_sglt2_bmi,
  res_othergld_bmi = res_othergld_bmi,
  pm_summary = pm_summary,
  hr_all = hr_all
)

saveRDS(mediation_results,
  file.path(output_dir, paste0("mediation_analysis_results_", timestamp, ".rds")))

cat("\nRDS file saved for future use.\n")

# >>> CELL 14: ============================================================
# SUMMARY OUTPUT
# ============================================================
cat("\n")
cat("============================================================\n")
cat("MEDIATION ANALYSIS COMPLETE - FINAL PAPER OUTPUTS\n
")
cat("============================================================\n")
cat("Output directory:", output_dir, "\n\n")

cat("FIGURES (PDF):\n")
cat("  - Figure_Mediation_SGLT2_A1c_", timestamp, ".pdf\n", sep="")
cat("  - Figure_Mediation_OtherGLD_A1c_", timestamp, ".pdf\n", sep="")
cat("  - Figure_Mediation_SGLT2_BMI_", timestamp, ".pdf\n", sep="")
cat("  - Figure_Mediation_OtherGLD_BMI_", timestamp, ".pdf\n", sep="")
cat("  - Supplementary_Figure_6_Proportion_Mediated_", timestamp, ".pdf\n", sep="")
cat("  - Supplementary_Figure_6_Proportion_Mediated_", timestamp, ".png\n", sep="")

cat("\nTABLES (CSV):\n")
cat("  - Table_Proportion_Mediated_Summary_", timestamp, ".csv\n", sep="")
cat("  - Table_Mediation_HR_Formatted_", timestamp, ".csv\n", sep="")
cat("  - Table_Mediation_HR_Full_", timestamp, ".csv\n", sep="")

cat("\nRDS (for future use):\n")
cat("  - mediation_analysis_results_", timestamp, ".rds\n", sep="")

cat("\n============================================================\n")
cat("PROPORTION MEDIATED SUMMARY:\n")
cat("============================================================\n")
print(pm_summary)

cat("\n============================================================\n")
cat("HR TABLE (48 months):\n")
cat("============================================================\n")
hr_48m <- hr_formatted %>% filter(`Time (months)` == 48)
print(hr_48m)
