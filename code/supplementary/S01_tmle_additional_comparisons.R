# >>> CELL 01: Install & Load Libraries (TMLE Only - Combined All Ages & Late Onset) <<<
# This script runs TMLE analysis for additional comparisons:
# - SGLT2 vs OtherGLD
# - OTHER_GLPA vs OtherGLD
# - OTHER_GLPA vs SGLT2
# For both ALL_AGES and LATE_ONSET

# Load required packages silently
suppressPackageStartupMessages({
  library(tidyverse)
  library(dplyr)
  library(tidyr)
  library(survival)
  library(lubridate)
  library(ggplot2)
  library(tableone)
  library(readr)
})

# Install TMLE packages if needed
if (!requireNamespace("tmle", quietly=TRUE)) install.packages("tmle")
if (!requireNamespace("SuperLearner", quietly=TRUE)) install.packages("SuperLearner")
suppressPackageStartupMessages(library(tmle))

# >>> CELL 02: Configuration <<<
get_standard_config <- function() {
  list(
    exposure_window = list(start = "2018-01-01", end = "2023-10-01"),
    data_cut_date = as.Date("2023-10-01"),
    drug_classes = list(
      SEMAGLUTIDE = c("semaglutide"),
      OTHER_GLPA = c("exenatide", "liraglutide", "albiglutide", "dulaglutide", "lixisenatide"),
      SGLT2 = c("canagliflozin", "empagliflozin", "dapagliflozin", "ertugliflozin"),
      OtherGLD = c("glimepiride", "glipizide", "glyburide", "alogliptin", "linagliptin", "sitagliptin", "saxagliptin", "pioglitazone", "rosiglitazone")
    ),
    comparisons = list(
      list(name = "SEMAGLUTIDE vs OTHER_GLPA", exposure = "SEMAGLUTIDE", comparator = "OTHER_GLPA"),
      list(name = "SEMAGLUTIDE vs SGLT2", exposure = "SEMAGLUTIDE", comparator = "SGLT2"),
      list(name = "SEMAGLUTIDE vs OtherGLD", exposure = "SEMAGLUTIDE", comparator = "OtherGLD"),
      list(name = "OTHER_GLPA vs SGLT2", exposure = "OTHER_GLPA", comparator = "SGLT2"),
      list(name = "OTHER_GLPA vs OtherGLD", exposure = "OTHER_GLPA", comparator = "OtherGLD"),
      list(name = "SGLT2 vs OtherGLD", exposure = "SGLT2", comparator = "OtherGLD")
    ),
    outcomes = list(
      list(var = "epilepsy_or_seizure_without_hypoglycemic_seizure_start_date",
           label = "Epilepsy/Seizure excl Hypoglycemic",
           late_onset = FALSE, early_onset = FALSE)
    ),
    all_ps_vars = c(
      "age", "sex_cat", "raceethnicity_cat", "income", "education", "insurance_category",
      "smoking", "alcohol_category", "baseline_bmi_category", "baseline_hba1c", "index_year_grouped",
      "Anticoagulant", "Antiplatelet", "BB", "Biguanide", "CCB", "Diuretic",
      "Ezetimibe", "MRA", "OtherHTN", "RAAS", "Statin", "TZD", "Insulin",
      "SEMAGLUTIDE", "OTHER_GLPA", "SGLT2i", "SU", "DPP4i",
      "myocardial_infarction", "congestive_heart_failure", "peripheral_vascular_disease",
      "cerebrovascular_disease", "chronic_pulmonary_disease", "dementia",
      "rheumatic_disease", "peptic_ulcer_disease", "hemiplegia_or_paraplegia", "hiv_infection",
      "hypoglycemia", "hyperglycemic_emergency",
      "renal_disease_severity", "liver_disease_severity",
      "diabetes_with_ophthalmic_complications", "diabetes_with_neurological_complications", "malignancy_status"
    ),
    categorical_vars = c(
      "sex_cat", "raceethnicity_cat", "income", "education", "insurance_category",
      "smoking", "alcohol_category", "index_year", "index_year_grouped", "baseline_bmi_category",
      "Anticoagulant", "Antiplatelet", "BB", "Biguanide", "CCB", "Diuretic",
      "Ezetimibe", "MRA", "OtherHTN", "RAAS", "Statin", "TZD", "Insulin",
      "SEMAGLUTIDE", "OTHER_GLPA", "SGLT2i", "SU", "DPP4i",
      "myocardial_infarction", "congestive_heart_failure", "peripheral_vascular_disease",
      "cerebrovascular_disease", "dementia", "chronic_pulmonary_disease",
      "rheumatic_disease", "peptic_ulcer_disease", "hemiplegia_or_paraplegia", "hiv_infection",
      "hypoglycemia", "hyperglycemic_emergency",
      "liver_disease_severity", "renal_disease_severity",
      "diabetes_with_ophthalmic_complications", "diabetes_with_neurological_complications", "malignancy_status"
    ),
    exclude_renal_severe = FALSE,
    exclude_malignancy_metastatic = FALSE
  )
}

config <- get_standard_config()

# >>> CELL 03: Logging Functions <<<
init_output_tracker <- function() {
  options(semaglutide_output_tracker = data.frame(
    timestamp=character(), script=character(), type=character(),
    message=character(), value=character(), stringsAsFactors=FALSE
  ))
}

log_output <- function(type, message, value="", script="") {
  if (is.null(getOption("semaglutide_output_tracker"))) init_output_tracker()
  tracker <- getOption("semaglutide_output_tracker")
  new_row <- data.frame(
    timestamp = format(Sys.time(), "%Y-%m-%d %H:%M:%S"),
    script = script, type = type, message = message, value = as.character(value),
    stringsAsFactors = FALSE
  )
  options(semaglutide_output_tracker = rbind(tracker, new_row))
  if (type %in% c("PROGRESS", "ERROR", "WARNING")) {
    cat(sprintf("[%s] %s: %s\n", type, message, value))
    try(flush(stdout()), silent = TRUE)
    if (interactive()) flush.console()
  }
}

save_output_log <- function(filename=NULL) {
  tracker <- getOption("semaglutide_output_tracker")
  if (is.null(tracker) || nrow(tracker)==0) return(invisible(NULL))
  if (is.null(filename)) filename <- sprintf("output_log_%s.csv", format(Sys.time(), "%Y%m%d_%H%M%S"))
  write.csv(tracker, file=filename, row.names=FALSE)
  log_output("PROGRESS", "Output log saved", filename)
}

# >>> CELL 04: TableOne Variables <<<
my_vars <- c("age", "sex_cat", "raceethnicity_cat", "income", "education", "insurance_category",
             "smoking", "alcohol_category", "baseline_bmi_category", "baseline_hba1c", "index_year_grouped",
             "Anticoagulant", "Antiplatelet", "BB", "Biguanide", "CCB", "Diuretic",
             "Ezetimibe", "MRA", "OtherHTN", "RAAS", "Statin", "TZD", "Insulin",
             "SEMAGLUTIDE", "OTHER_GLPA", "SGLT2i", "SU", "DPP4i",
             "myocardial_infarction", "congestive_heart_failure", "peripheral_vascular_disease",
             "cerebrovascular_disease", "chronic_pulmonary_disease", "dementia",
             "rheumatic_disease", "peptic_ulcer_disease", "hemiplegia_or_paraplegia", "hiv_infection",
             "hypoglycemia", "hyperglycemic_emergency",
             "renal_disease_severity", "liver_disease_severity",
             "diabetes_with_ophthalmic_complications", "diabetes_with_neurological_complications", "malignancy_status")
categorical_vars <- my_vars[!my_vars %in% c("age", "baseline_bmi", "baseline_hba1c")]

WRITE_OUTPUTS <- TRUE

# >>> CELL 05: TMLE Core Functions <<<
COVARIATE_CAP <- NULL

run_tmle_binary <- function(df, exclude_vars=NULL, all_ps_vars=NULL, sl_lib=c("SL.glm"), verbose=FALSE, use_imputation=TRUE, use_superlearner=FALSE) {
  t_start <- Sys.time()
  all_covariates <- if (is.null(all_ps_vars)) get_standard_config()$all_ps_vars else all_ps_vars
  rhs_vars <- setdiff(all_covariates, exclude_vars)
  keep_vars <- rhs_vars[sapply(df[rhs_vars], function(x) length(unique(na.omit(x)))>1)]
  model_vars <- c("treatment", keep_vars)
  if (verbose) log_output("PROGRESS","TMLE: begin run_tmle_binary", sprintf("n=%d; covariates_kept=%d", nrow(df), length(keep_vars)), "tmle_fit")

  df_impute <- df  # Skip imputation since IPTW already did it

  complete_rows <- complete.cases(df_impute[, model_vars, drop=FALSE])
  df_complete <- df_impute[complete_rows, , drop=FALSE]
  if (verbose) log_output("INFO", sprintf("TMLE: complete cases retained = %d / %d", nrow(df_complete), nrow(df)))

  W_vars <- df_complete[, keep_vars, drop=FALSE]
  if (!is.null(COVARIATE_CAP) && ncol(W_vars)>COVARIATE_CAP) W_vars <- W_vars[, seq_len(COVARIATE_CAP), drop=FALSE]
  if (ncol(W_vars)==0) W_vars <- data.frame(intercept=rep(1, nrow(df_complete)))

  use_sl <- isTRUE(use_superlearner) && requireNamespace("SuperLearner", quietly=TRUE)
  fit <- NULL

  if (use_sl) {
    sl_present <- unique(sl_lib)
    fit <- try(tmle::tmle(Y=df_complete$event, A=df_complete$treatment, W=W_vars, family="binomial", Q.SL.library=sl_present, g.SL.library=sl_present), silent=TRUE)
  }

  if (is.null(fit) || inherits(fit,"try-error")) {
    data_for_models <- cbind(df_complete[, c("event","treatment"), drop=FALSE], W_vars)
    g_form <- as.formula(paste("treatment ~", paste(colnames(W_vars), collapse=" + ")))
    g_fit <- try(stats::glm(g_form, data=data_for_models, family=stats::binomial()), silent=TRUE)
    g1W <- suppressWarnings(pmin(pmax(as.numeric(stats::predict(g_fit, type="response")), 1e-5), 1-1e-5))
    q_form <- as.formula(paste("event ~ treatment +", paste(colnames(W_vars), collapse=" + ")))
    q_fit <- try(stats::glm(q_form, data=data_for_models, family=stats::binomial()), silent=TRUE)
    newdata1 <- data_for_models; newdata1$treatment <- 1
    newdata0 <- data_for_models; newdata0$treatment <- 0
    Q1W <- suppressWarnings(pmin(pmax(as.numeric(stats::predict(q_fit, newdata=newdata1, type="response")), 1e-5), 1-1e-5))
    Q0W <- suppressWarnings(pmin(pmax(as.numeric(stats::predict(q_fit, newdata=newdata0, type="response")), 1e-5), 1-1e-5))
    Q_init <- cbind(Q1W=Q1W, Q0W=Q0W)
    fit <- tmle::tmle(Y=df_complete$event, A=df_complete$treatment, W=W_vars, family="binomial", Q=Q_init, g1W=g1W)
  }

  list(fit=fit, cohort=df_complete, keep_vars=keep_vars, n_dropped=nrow(df)-nrow(df_complete))
}

tmle_sl_lib <- c("SL.glm")
if (requireNamespace("glmnet", quietly=TRUE)) tmle_sl_lib <- c("SL.glmnet", tmle_sl_lib)
if (requireNamespace("xgboost", quietly=TRUE)) tmle_sl_lib <- c("SL.xgboost", tmle_sl_lib)

calculate_tmle_time_specific_rd <- function(df, time_points, time_labels) {
  if (!all(c("event_time","event","treatment") %in% names(df))) return(NULL)
  safe_len <- min(length(time_points), length(time_labels))
  if (safe_len==0) return(NULL)
  time_points <- time_points[seq_len(safe_len)]
  time_labels <- time_labels[seq_len(safe_len)]

  calc_surv_at <- function(data_subset, t_point) {
    if (nrow(data_subset)<1) return(list(surv=NA_real_, se=NA_real_))
    fit <- survival::survfit(survival::Surv(event_time, event) ~ 1, data=data_subset)
    s <- try(suppressWarnings(summary(fit, times=t_point)), silent=TRUE)
    if (inherits(s,"try-error") || length(s$surv)==0) return(list(surv=NA_real_, se=NA_real_))
    se <- if (!is.null(s$std.err) && length(s$std.err)>0) s$std.err[1] else NA_real_
    list(surv=s$surv[1], se=se)
  }

  out <- vector("list", length(time_points))
  for (i in seq_along(time_points)) {
    t_point <- time_points[i]; lbl <- time_labels[i]
    grp0 <- df[df$treatment==0, , drop=FALSE]
    grp1 <- df[df$treatment==1, , drop=FALSE]
    s0 <- calc_surv_at(grp0, t_point)
    s1 <- calc_surv_at(grp1, t_point)
    risk0 <- if (!is.na(s0$surv)) 1-s0$surv else NA_real_
    risk1 <- if (!is.na(s1$surv)) 1-s1$surv else NA_real_
    rd <- if (!is.na(risk0) && !is.na(risk1)) risk1-risk0 else NA_real_
    var_rd <- sum(c(s0$se, s1$se)^2, na.rm=TRUE)
    se_rd <- if (is.finite(var_rd)) sqrt(var_rd) else NA_real_
    lo <- if (!is.na(rd) && !is.na(se_rd)) rd-1.96*se_rd else NA_real_
    hi <- if (!is.na(rd) && !is.na(se_rd)) rd+1.96*se_rd else NA_real_
    out[[i]] <- data.frame(time_label=lbl, time_days=t_point, risk_treated=risk1, risk_control=risk0, rd=rd, rd_lo_95=lo, rd_hi_95=hi, stringsAsFactors=FALSE)
  }
  do.call(rbind, out)
}

# >>> CELL 06: TMLE Analysis Function <<<
run_tmle_for_analysis <- function(analysis_name, iptw_rds_path, output_prefix, allowed_comparisons) {
  message(sprintf("\n========================================"))
  message(sprintf("=== TMLE ANALYSIS: %s ===", toupper(analysis_name)))
  message(sprintf("========================================"))
  message(sprintf("Comparisons: %s", paste(allowed_comparisons, collapse=", ")))
  message(sprintf("Loading IPTW results from: %s", iptw_rds_path))

  # Load IPTW results
  iptw_obj <- tryCatch({
    readRDS(iptw_rds_path)
  }, error = function(e) {
    log_output("ERROR", sprintf("Failed to load IPTW RDS file: %s", iptw_rds_path), e$message, "tmle_main")
    return(NULL)
  })

  if (is.null(iptw_obj)) {
    message("Skipping this analysis due to missing IPTW file.")
    return(NULL)
  }

  message(sprintf("Loaded IPTW results: %d analyses found", length(iptw_obj)))

  tmle_results <- list()
  tmle_summary_rows <- list()

  for (comp_name in allowed_comparisons) {
    log_output("PROGRESS", sprintf("TMLE: Processing comparison: %s", comp_name), "", "tmle_main")

    exclude_vars <- switch(comp_name,
      "SEMAGLUTIDE vs OTHER_GLPA"=c("SEMAGLUTIDE","OTHER_GLPA"),
      "SEMAGLUTIDE vs SGLT2"=c("SEMAGLUTIDE","SGLT2i"),
      "SEMAGLUTIDE vs OtherGLD"=c("SEMAGLUTIDE","TZD","SU","DPP4i"),
      "OTHER_GLPA vs SGLT2"=c("OTHER_GLPA","SGLT2i"),
      "OTHER_GLPA vs OtherGLD"=c("OTHER_GLPA","TZD","SU","DPP4i"),
      "SGLT2 vs OtherGLD"=c("SGLT2i","TZD","SU","DPP4i"), NULL)

    for (outcome in config$outcomes) {
      outcome_var <- outcome$var
      outcome_label <- outcome$label
      analysis_key <- paste(comp_name, "-", outcome_label)
      log_output("INFO", sprintf("TMLE: Processing outcome: %s for %s (%s)", outcome_label, comp_name, analysis_name), "", "tmle_main")

      # Load cohort from IPTW results
      if (!is.null(iptw_obj[[analysis_key]]) && !is.null(iptw_obj[[analysis_key]]$cohort)) {
        survival_df <- iptw_obj[[analysis_key]]$cohort
        log_output("PROGRESS","TMLE: reusing exact IPTW cohort", sprintf("n=%d", nrow(survival_df)), "tmle_main")
      } else {
        log_output("WARNING", sprintf("No IPTW cohort found for: %s", analysis_key), "", "tmle_main")
        next
      }

      if (nrow(survival_df)<10) {
        log_output("WARNING","Skipping TMLE due to insufficient data", nrow(survival_df), "tmle_main")
        next
      }

      n_t1 <- sum(survival_df$treatment==1)
      n_t0 <- sum(survival_df$treatment==0)
      if (n_t1<2 || n_t0<2) {
        log_output("ERROR","Insufficient observations", sprintf("Treatment=1: %d, Treatment=0: %d", n_t1, n_t0), "tmle_main")
        next
      }

      # Run TMLE
      tmle_fit <- tryCatch({
        run_tmle_binary(df=survival_df, exclude_vars=exclude_vars, all_ps_vars=config$all_ps_vars, sl_lib=tmle_sl_lib, verbose=TRUE, use_imputation=FALSE, use_superlearner=TRUE)
      }, error=function(e) {
        log_output("ERROR", sprintf("run_tmle_binary failed for %s", comp_name), e$message, "tmle_main")
        return(NULL)
      })

      if (is.null(tmle_fit)) next

      tmle_results[[analysis_key]] <- tmle_fit
      log_output("PROGRESS","TMLE: fit stored", sprintf("analysis=%s; comp=%s; outcome=%s", analysis_name, comp_name, outcome_label), "tmle_main")

      # Plot TMLE estimated event probabilities
      ey0_psi <- try(tmle_fit$fit$estimates$EY0$psi, silent=TRUE)
      ey1_psi <- try(tmle_fit$fit$estimates$EY1$psi, silent=TRUE)
      ey0_lo <- try(tmle_fit$fit$estimates$EY0$CI[1], silent=TRUE)
      ey0_hi <- try(tmle_fit$fit$estimates$EY0$CI[2], silent=TRUE)
      ey1_lo <- try(tmle_fit$fit$estimates$EY1$CI[1], silent=TRUE)
      ey1_hi <- try(tmle_fit$fit$estimates$EY1$CI[2], silent=TRUE)
      ate_obj <- try(tmle_fit$fit$estimates$ATE, silent=TRUE)
      ate_psi <- if (!inherits(ate_obj,"try-error") && !is.null(ate_obj$psi)) ate_obj$psi else NA_real_
      ate_p <- if (!inherits(ate_obj,"try-error") && !is.null(ate_obj$pvalue)) ate_obj$pvalue else NA_real_

      to_num <- function(x) if (!inherits(x,"try-error")) as.numeric(x) else NA_real_
      prob_df <- data.frame(
        Group=c("Control","Treatment"),
        Probability=c(to_num(ey0_psi), to_num(ey1_psi)),
        Lower=c(to_num(ey0_lo), to_num(ey1_lo)),
        Upper=c(to_num(ey0_hi), to_num(ey1_hi))
      )

      try({
        tmle_plot <- ggplot2::ggplot(prob_df, ggplot2::aes(x=Group, y=Probability, fill=Group)) +
          ggplot2::geom_bar(stat="identity", alpha=0.7) +
          ggplot2::geom_errorbar(ggplot2::aes(ymin=Lower, ymax=Upper), width=0.2) +
          ggplot2::scale_fill_manual(values=c("Control"="#1f77b4","Treatment"="#ff7f0e")) +
          ggplot2::labs(
            title=paste("TMLE Estimated Event Probabilities for", outcome_label),
            subtitle=paste(comp_name,"- Risk Difference =", sprintf("%.4f", ate_psi), "(p =", ifelse(is.na(ate_p),"NA", sprintf("%.3g", ate_p)), ")"),
            y="Probability of Event", x=""
          ) +
          ggplot2::theme_minimal() +
          ggplot2::theme(legend.position="none")
        print(tmle_plot)
        safe_comp <- gsub("[^A-Za-z0-9]+","_",comp_name)
        safe_outcome <- gsub("[^A-Za-z0-9]+","_",outcome_label)
        ggplot2::ggsave(filename=paste0("tmle_probs_", analysis_name, "_", safe_comp, "_", safe_outcome, ".png"), plot=tmle_plot, width=8, height=6)
      }, silent=TRUE)

      # Time-specific results
      time_points <- c(182.5, 365, 547.5, 730)
      time_labels <- c("6_months","1_year","18_months","2_years")
      time_specific_results <- tryCatch({
        calculate_tmle_time_specific_rd(tmle_fit$cohort, time_points, time_labels)
      }, error=function(e) { NULL })

      # Save results
      ts <- format(Sys.time(),"%Y%m%d_%H%M%S")
      safe_comp <- gsub("[^A-Za-z0-9]+","_",comp_name)
      safe_outcome <- gsub("[^A-Za-z0-9]+","_",outcome_label)
      tmle_prefix <- output_prefix

      ate <- tmle_fit$fit$estimates$ATE
      or <- tmle_fit$fit$estimates$OR
      EY1 <- try(tmle_fit$fit$estimates$EY1$psi, silent=TRUE)
      EY0 <- try(tmle_fit$fit$estimates$EY0$psi, silent=TRUE)
      rr_obj <- try(tmle_fit$fit$estimates$RR, silent=TRUE)
      rr_est <- if (!inherits(rr_obj,"try-error") && !is.null(rr_obj$psi)) rr_obj$psi else if (!inherits(EY1,"try-error") && !inherits(EY0,"try-error") && !is.null(EY0) && EY0>0) EY1/EY0 else NA_real_
      rr_lo <- if (!inherits(rr_obj,"try-error") && !is.null(rr_obj$CI)) rr_obj$CI[1] else NA_real_
      rr_hi <- if (!inherits(rr_obj,"try-error") && !is.null(rr_obj$CI)) rr_obj$CI[2] else NA_real_
      rr_p <- if (!inherits(rr_obj,"try-error") && !is.null(rr_obj$pvalue)) rr_obj$pvalue else NA_real_
      nnt_val <- if (!is.null(ate$psi) && is.finite(ate$psi) && ate$psi!=0) 1/abs(ate$psi) else Inf
      n_treat <- sum(tmle_fit$cohort$treatment==1, na.rm=TRUE)
      n_ctrl <- sum(tmle_fit$cohort$treatment==0, na.rm=TRUE)
      ev_treat <- sum(tmle_fit$cohort$event==1 & tmle_fit$cohort$treatment==1, na.rm=TRUE)
      ev_ctrl <- sum(tmle_fit$cohort$event==1 & tmle_fit$cohort$treatment==0, na.rm=TRUE)

      out_df <- data.frame(
        method="TMLE", analysis_type=analysis_name, comparison=comp_name, outcome=outcome_label,
        effect_measure=c("Risk Difference","Odds Ratio","Risk Ratio","Number Needed to Treat"),
        estimate=c(ate$psi, or$psi, rr_est, nnt_val),
        lower_ci=c(ate$CI[1], or$CI[1], rr_lo, NA_real_),
        upper_ci=c(ate$CI[2], or$CI[2], rr_hi, NA_real_),
        p_value=c(ate$pvalue, or$pvalue, rr_p, NA_real_),
        n_total=rep(nrow(tmle_fit$cohort),4),
        n_events=rep(sum(tmle_fit$cohort$event, na.rm=TRUE),4),
        n_at_risk_treatment=rep(n_treat,4),
        n_at_risk_control=rep(n_ctrl,4),
        events_treatment=rep(ev_treat,4),
        events_control=rep(ev_ctrl,4),
        stringsAsFactors=FALSE
      )
      tmle_summary_rows[[length(tmle_summary_rows)+1]] <- out_df

      fname <- sprintf("%s_summary_%s_%s_%s.csv", tmle_prefix, safe_comp, safe_outcome, ts)
      if (isTRUE(WRITE_OUTPUTS)) try(utils::write.csv(out_df, file=fname, row.names=FALSE), silent=TRUE)

      if (!is.null(time_specific_results) && isTRUE(WRITE_OUTPUTS)) {
        fname_time <- sprintf("%s_time_specific_RD_%s_%s_%s.csv", tmle_prefix, safe_comp, safe_outcome, ts)
        try(utils::write.csv(time_specific_results, file=fname_time, row.names=FALSE), silent=TRUE)
      }

      # Baseline table for TMLE cohort
      try({
        tmle_cohort <- tmle_fit$cohort
        tmle_vars <- intersect(my_vars, names(tmle_cohort))
        tmle_cat_vars <- intersect(categorical_vars, tmle_vars)
        for (v in tmle_cat_vars) {
          if (!is.factor(tmle_cohort[[v]])) tmle_cohort[[v]] <- as.factor(tmle_cohort[[v]])
        }
        if (length(tmle_vars)>0) {
          tmle_table_obj <- CreateTableOne(vars=tmle_vars, strata="treatment", data=tmle_cohort, factorVars=tmle_cat_vars)
          tmle_table_df <- as.data.frame(print(tmle_table_obj, showAllLevels=TRUE, printToggle=FALSE, noSpaces=TRUE))
          ts2 <- format(Sys.time(),"%Y%m%d_%H%M%S")
          fname2 <- sprintf("%s_baseline_table_%s_%s_%s.csv", tmle_prefix, safe_comp, safe_outcome, ts2)
          if (isTRUE(WRITE_OUTPUTS)) {
            write.csv(tmle_table_df, file=fname2, row.names=TRUE)
            log_output("PROGRESS","Saved TMLE baseline table", fname2, "baseline_table_TMLE")
          }
        }
      }, silent=TRUE)
    }
  }

  # Save all results
  rds_filename <- sprintf("%s_full_results.rds", output_prefix)
  if (isTRUE(WRITE_OUTPUTS)) {
    saveRDS(tmle_results, file=rds_filename)
    log_output("PROGRESS","TMLE results object saved to", rds_filename, "tmle_main")
  }

  timestamp <- format(Sys.time(),"%Y%m%d_%H%M%S")
  if (isTRUE(WRITE_OUTPUTS)) save_output_log(sprintf("%s_log_%s.csv", output_prefix, timestamp))

  if (length(tmle_summary_rows)>0 && isTRUE(WRITE_OUTPUTS)) {
    final_tmle_df <- tryCatch({ do.call(rbind, tmle_summary_rows) }, error=function(e) NULL)
    if (!is.null(final_tmle_df)) {
      comp_results_file <- sprintf("%s_comprehensive_results_%s.csv", output_prefix, timestamp)
      utils::write.csv(final_tmle_df, comp_results_file, row.names=FALSE)
      log_output("PROGRESS","Aggregated TMLE comprehensive results saved", comp_results_file, "tmle_main")
    }
  }

  log_output("PROGRESS", sprintf("TMLE ANALYSIS '%s' COMPLETE", analysis_name), "", "tmle_main")
  return(tmle_results)
}

# >>> CELL 07: Run TMLE for Both Analysis Types <<<
message("\n")
message("##############################################################")
message("# TMLE ANALYSIS - ADDITIONAL COMPARISONS                     #")
message("# Comparisons: SGLT2 vs OtherGLD, OTHER_GLPA vs OtherGLD,    #")
message("#              OTHER_GLPA vs SGLT2                           #")
message("# Analysis Types: ALL_AGES and LATE_ONSET                    #")
message("##############################################################")

# Define comparisons
additional_comparisons <- c("SGLT2 vs OtherGLD", "OTHER_GLPA vs OtherGLD", "OTHER_GLPA vs SGLT2")

# Run for ALL_AGES
all_ages_results <- run_tmle_for_analysis(
  analysis_name = "all_ages",
  iptw_rds_path = "outcome4_epilepsy_excl_hypo_all_ages_ipwt_full_results.rds",
  output_prefix = "outcome4_epilepsy_excl_hypo_all_ages_tmle_additional",
  allowed_comparisons = additional_comparisons
)

# Run for LATE_ONSET
late_onset_results <- run_tmle_for_analysis(
  analysis_name = "late_onset",
  iptw_rds_path = "outcome4_epilepsy_excl_hypo_late_onset_ipwt_full_results.rds",
  output_prefix = "outcome4_epilepsy_excl_hypo_late_onset_tmle_additional",
  allowed_comparisons = additional_comparisons
)

message("\n")
message("##############################################################")
message("# ALL TMLE ANALYSES COMPLETE                                 #")
message("##############################################################")
