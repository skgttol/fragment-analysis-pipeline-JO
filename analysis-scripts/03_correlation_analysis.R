#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                             #
#         Advanced Fragment Analysis Workflow - Script 3 of 5                 #
#                   (Correlation & Integration)                               #
#                                                                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# DESCRIPTION:
# This script determines the "Instability Threshold" using AIC model selection.
#
# FEATURES:
# 1. 7-Way Model Tournament (Linear, Quad, 2-Ph, 3-Ph, Hockey, Exp, Power).
# 2. Dynamic Genotype Labeling (Uses labels from Script 01).
# 3. Robust "Brute Force" Scanning (Prevents crashes on small data).
# 4. Full Diagnostic Suite (Actual vs Predicted, Alternative Plots).
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#=============================================================================#
# PART 0: SETUP AND LOAD PROCESSED DATA                                    ####
#=============================================================================#

# 0A. Load Libraries ----
packages <- c("tidyverse", "readxl", "here", "openxlsx", "broom.mixed",
              "ggpubr", "janitor", "yaml", "logr", "lme4", "lmerTest",
              "emmeans", "ggeffects", "patchwork", "gridExtra", "grid",
              "RColorBrewer", "gtable", "ggnewscale", "GGally", "reshape2", "MASS",
              "ggrepel", "dendextend", "boot", "segmented") 

installed_packages <- packages %in% rownames(utils::installed.packages())
if (any(installed_packages == FALSE)) {
  utils::install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

# 0B. Load Configuration & Find Data ----
if (!file.exists(here::here("config.yml"))) stop("CRITICAL ERROR: 'config.yml' not found.")
config <- yaml::read_yaml(here::here("config.yml"))

# --- 1. Bulletproof Latest Directory Search ---
analysis_base_dir <- here::here(config$paths$output_dir_base)

# Get all directories, but explicitly remove the parent base directory from the list
all_dirs <- list.dirs(path = analysis_base_dir, full.names = TRUE, recursive = FALSE)
all_dirs <- setdiff(all_dirs, analysis_base_dir)

if(length(all_dirs) == 0) stop("CRITICAL ERROR: No analysis output folders found in base directory.")

# Filter strictly for folders that contain your specific naming convention string
target_dirs <- all_dirs[grepl("_Analysis_v", basename(all_dirs))]

if(length(target_dirs) == 0) stop("CRITICAL ERROR: No folders matching the '_Analysis_v' convention were found.")

# Because your format is YYYY-MM-DD...HHMM, alphabetical max() safely finds the newest run
latest_analysis_dir <- max(target_dirs)

logr::log_print(paste("Successfully identified latest analysis directory:", basename(latest_analysis_dir)), console = TRUE)

rdata_path <- file.path(latest_analysis_dir, "processing_complete.RData")
if (!file.exists(rdata_path)) stop("ERROR: 'processing_complete.RData' not found.")

load(rdata_path)
# Ensure functions are loaded

## 0B. Load Functions ----
#-------------------------#
if (!file.exists(here::here("functions.R"))) {
  stop("CRITICAL ERROR: 'functions.R' not found. Please ensure it is in the main project directory.")
}
source(here::here("functions.R"))

# func_path <- "C:/Users/skgttol/OneDrive - University College London/PhD/PhD_Thesis/02_Data/Template_RScripts/CAGsizing_Flexible/functions.R"
# if (!file.exists(func_path)) {
#   if(file.exists("functions.R")) func_path <- "functions.R" else stop("functions.R not found.")
# }
# source(func_path)

if (!exists("modeling_data")) {
  stop(paste(
    "ERROR: 'modeling_data' not found in loaded environment.",
    "Script 03 requires Script 02 to have been run first.",
    "Run 02_generate_plots.R before this script."
  ))
}

# 3. Setup Logging
try(logr::log_close(), silent = TRUE)
if (!exists("output_dir")) {
  output_dir <- file.path(latest_analysis_dir, "results") 
  dir.create(output_dir, showWarnings = FALSE)
}
log_path <- file.path(output_dir, "analysis_log.log") 
lf <- logr::log_open(log_path, show_notes = FALSE, logdir = FALSE)

logr::log_print("\n\n--- SCRIPT 03: CORRELATION & INTEGRATION INITIALIZED ---")

# Setup Directories
corr_plot_dir <- file.path(output_dir, "04_correlation")
cag_exp_plot_dir <- file.path(corr_plot_dir, "cag_vs_mode")
alternatives_dir <- file.path(cag_exp_plot_dir, "alternatives") 
matrix_plot_dir <- file.path(corr_plot_dir, "corr_matrix")
pca_plot_dir <- file.path(corr_plot_dir, "pca")
diag_plot_dir <- file.path(corr_plot_dir, "diagnostics")

dir.create(corr_plot_dir, showWarnings = FALSE)
dir.create(cag_exp_plot_dir, showWarnings = FALSE)
dir.create(alternatives_dir, showWarnings = FALSE)
dir.create(matrix_plot_dir, showWarnings = FALSE)
dir.create(pca_plot_dir, showWarnings = FALSE)
dir.create(diag_plot_dir, showWarnings = FALSE)

# 4. Define key variables
primary_var   <- config$key_variables$primary_group_var
secondary_var <- config$key_variables$secondary_group_var
time_var      <- config$key_variables$time_variable
diff_var      <- config$key_variables$optional_grouping_var 
if(is.null(secondary_var) || secondary_var == 'null') grouping_var <- primary_var else grouping_var <- secondary_var

# --- get_lbl() WRAPPER ---
# The canonical implementation lives in functions.R (Section 6).
# This local wrapper binds the script-level `config` and `modeling_data`
# so call sites don't need to pass them explicitly each time.
# NOTE: modeling_data may not exist yet at this point in the script;
# get_lbl() handles that gracefully by falling back to primary_group_var.
get_lbl <- function(style = "auto") {
  data_to_check <- if (exists("modeling_data")) modeling_data else NULL
  get_lbl_canonical(style = style, data = data_to_check, config = config)
}

plot_label_col <- get_lbl("auto")
logr::log_print(paste("Using Label Column:", plot_label_col))

# Create Lookup Map
if(exists("modeling_data")) {
  label_lookup <- modeling_data %>% 
    dplyr::distinct(!!sym(primary_var), !!sym(plot_label_col))
}

# FIX: Establish the active palette globally from the loaded RData
active_palette <- if(exists("primary_palette")) primary_palette else NULL

#=============================================================================#
# PART 1: DEFINING ROBUST STATISTICAL FUNCTIONS                           #####
#=============================================================================#

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 1A. Helper: Calculate R-Squared
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
calc_population_r2 <- function(model, data, y_col = "slope_val") {
  tryCatch({
    if(inherits(model, "lmerMod")) {
      preds <- predict(model, newdata = data, re.form = NA)
    } else {
      preds <- predict(model, newdata = data)
    }
    r2 <- cor(data[[y_col]], preds, use = "complete.obs")^2
    return(r2)
  }, error = function(e) return(NA))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 1B. Helper: Quadratic Slope Crossing Calculation (NEW)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
get_quad_crossing_point <- function(model_lin, model_quad, data, divergence_percent = 0.20) {
  # 1. Get Linear Slope (m)
  # We use the fixed effect for 'baseline_cag'
  m <- fixef(model_lin)["baseline_cag"]
  
  # 2. Get Quadratic Parameters (ax^2 + bx + c)
  cf_q <- fixef(model_quad)
  a <- cf_q[grep("I\\(", names(cf_q))] # Coefficient for x^2
  b <- cf_q["baseline_cag"]            # Coefficient for x
  
  # Safety: If curve is concave (upside down U) or flat, no acceleration threshold exists.
  if(length(a) == 0 || a <= 0) return(NA)
  
  # 3. Solve for X where Quad Slope > Linear Slope
  # 2ax + b = m  =>  2ax = m - b  =>  x = (m - b) / 2a
  crossing_x <- (m - b) / (2 * a)
  
  # 4. Return
  return(as.numeric(crossing_x))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 1C. Helper: AIC Profile Scan Plotter (ROBUST VERSION)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
generate_aic_profile_plot <- function(hockey_scan_df, seg1_scan_df, 
                                      aic_lin, aic_quad, aic_3ph, aic_exp, aic_pow, aic_loglin) {
  
  # 1. Clean the scan dataframes safely
  seg1_scan_df <- if (!is.null(seg1_scan_df)) na.omit(seg1_scan_df) else data.frame()
  hockey_scan_df <- if (!is.null(hockey_scan_df)) na.omit(hockey_scan_df) else data.frame()
  
  # If ALL models failed (very rare), abort to prevent a blank square
  if(nrow(seg1_scan_df) == 0 && nrow(hockey_scan_df) == 0 && 
     is.infinite(aic_lin) && is.infinite(aic_quad)) return(NULL)
  
  # 2. Define Reference Lines Dataframe for Legend
  ref_lines <- data.frame(
    Model = c("Linear", "Quadratic", "3-Phase (Free)", "Exponential", "Variable Power", "Log-Linear (Y-Transformed)"),
    AIC = c(aic_lin, aic_quad, aic_3ph, aic_exp, aic_pow, aic_loglin),
    Color = c("black", "forestgreen", "orange", "magenta", "darkcyan", "brown4"),
    Type = c("dotted", "longdash", "dashed", "dotdash", "twodash", "solid")
  )
  
  # Remove Infinite AICs (failed models)
  ref_lines <- ref_lines[is.finite(ref_lines$AIC), ]
  
  # 3. Base plot with horizontal reference lines
  p_scan <- ggplot() +
    geom_hline(data = ref_lines, aes(yintercept = AIC, color = Model, linetype = Model), linewidth = 0.8, alpha = 0.8)
  
  # 4. Add 2-Phase scan ONLY if valid data exists
  if (nrow(seg1_scan_df) > 0) {
    p_scan <- p_scan +
      geom_line(data = seg1_scan_df, aes(x = Threshold, y = AIC, color = "2-Phase (Free)", linetype = "2-Phase (Free)"), linewidth = 0.8) +
      geom_point(data = seg1_scan_df, aes(x = Threshold, y = AIC, color = "2-Phase (Free)"), size = 1.5, alpha = 0.6) +
      geom_point(data = seg1_scan_df[which.min(seg1_scan_df$AIC),], aes(x = Threshold, y = AIC), color = "blue", size = 4, shape = 19)
  }
  
  # 5. Add Hockey scan ONLY if valid data exists
  if (nrow(hockey_scan_df) > 0) {
    p_scan <- p_scan +
      geom_line(data = hockey_scan_df, aes(x = Threshold, y = AIC, color = "Hockey Stick (Fixed)", linetype = "Hockey Stick (Fixed)"), linewidth = 0.8) +
      geom_point(data = hockey_scan_df, aes(x = Threshold, y = AIC, color = "Hockey Stick (Fixed)"), size = 1.5, alpha = 0.6) +
      geom_point(data = hockey_scan_df[which.min(hockey_scan_df$AIC),], aes(x = Threshold, y = AIC), color = "purple", size = 4, shape = 19)
  }
  
  # 6. Scales and Labels
  p_scan <- p_scan + 
    scale_color_manual(name = "Model", values = c(
      "Linear" = "black", 
      "Quadratic" = "forestgreen", 
      "3-Phase (Free)" = "orange", 
      "Exponential" = "magenta", 
      "Variable Power" = "darkcyan", 
      "Log-Linear (Y-Transformed)" = "brown4", 
      "2-Phase (Free)" = "blue", 
      "Hockey Stick (Fixed)" = "purple"
    )) +
    scale_linetype_manual(name = "Model", values = c(
      "Linear" = "dotted", 
      "Quadratic" = "longdash", 
      "3-Phase (Free)" = "dashed", 
      "Exponential" = "dotdash",
      "Variable Power" = "twodash", 
      "Log-Linear (Y-Transformed)" = "solid", 
      "2-Phase (Free)" = "solid", 
      "Hockey Stick (Fixed)" = "solid"
    )) +
    labs(title = "AIC Profile Scan", subtitle = "Lower AIC = Better Fit", x = "Tested Threshold (CAG)", y = "AIC") +
    theme_publication() + 
    theme(legend.position = "right", legend.direction = "vertical")
  
  return(p_scan)
}
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 1D. MAIN FUNCTION: BRUTE FORCE MIXED EFFECTS TOURNAMENT (7-WAY)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
get_robust_mixed_analysis <- function(df, x_col, y_col, group_var, n_boot = 500, divergence_percent = 0.20, boundary_buffer = 10) { 
  work_df <- data.frame(baseline_cag = df[[x_col]], slope_val = df[[y_col]], group_id = df[[group_var]]) 
  work_df <- na.omit(work_df)
  min_valid <- min(work_df$baseline_cag) + boundary_buffer
  max_valid <- max(work_df$baseline_cag) - boundary_buffer
  
  # Grids for Scanning
  scan_seq <- unique(round(seq(min_valid, max_valid, length.out = 50))) 
  scan_exp_k <- seq(0.01, 0.15, length.out = 20)  
  scan_pow_p <- seq(1.5, 5.0, length.out = 20)    
  
  # Global Control to Silence Boundary/Singular Warnings
  safe_ctrl <- lmerControl(
    check.conv.singular = .makeCC(action = "ignore", tol = 1e-4),
    check.scaleX = "ignore",
    check.nobs.vs.nlev = "ignore",
    check.nobs.vs.rankZ = "ignore",
    check.nobs.vs.nRE = "ignore",
    optCtrl = list(print_level = 0) # Tells the optimizer to shut up
  )
  
  # Helper: Safe LMER Fit (SILENCED)
  fit_lmer_safe <- function(formula, data) {
    tryCatch({
      m <- suppressWarnings(suppressMessages({
        lmer(formula, data = data, REML = FALSE, control = safe_ctrl)
      }))
      
      if("(Intercept)" %in% names(fixef(m)) && fixef(m)["(Intercept)"] < 0) {
        f_char <- deparse(formula)
        f_new <- as.formula(paste(f_char, "- 1"))
        m <- suppressWarnings(suppressMessages({
          lmer(f_new, data = data, REML = FALSE, control = safe_ctrl)
        }))
      }
      return(m)
    }, error = function(e) return(NULL))
  }
  
  # --- 1. Standard Models ---
  fit_lin <- fit_lmer_safe(slope_val ~ baseline_cag + (1 | group_id), work_df)
  fit_quad <- fit_lmer_safe(slope_val ~ baseline_cag + I(baseline_cag^2) + (1 | group_id), work_df)
  
  # --- 2. 2-Phase Free ---
  best_seg1_aic <- Inf; best_seg1_thresh <- NA; fit_seg1_best <- NULL
  seg1_scan_df <- data.frame(Threshold = scan_seq, AIC = NA)
  
  for(i in 1:length(scan_seq)) {
    t <- scan_seq[i]
    work_df$term_diff <- pmax(work_df$baseline_cag - t, 0)
    fit <- fit_lmer_safe(slope_val ~ baseline_cag + term_diff + (1 | group_id), work_df)
    if(!is.null(fit)) {
      cf <- fixef(fit)
      has_int <- "(Intercept)" %in% names(cf); s1 <- if(has_int) cf["baseline_cag"] else cf[1]; s_diff <- if(has_int) cf["term_diff"] else cf[2]
      # FIX: Added is.na() checks
      if(!is.na(s1) && !is.na(s_diff) && s1 >= 0 && (s1 + s_diff) > s1) {
        this_aic <- AIC(fit)
        seg1_scan_df$AIC[i] <- this_aic
        if(this_aic < best_seg1_aic) { best_seg1_aic <- this_aic; best_seg1_thresh <- t; fit_seg1_best <- fit }
      }
    }
  }
  
  # --- 3. Hockey Stick ---
  best_hockey_aic <- Inf; best_hockey_thresh <- NA; fit_hockey_best <- NULL
  hockey_scan_df <- data.frame(Threshold = scan_seq, AIC = NA)
  
  for(i in 1:length(scan_seq)) {
    t <- scan_seq[i]
    work_df$hs_term <- pmax(work_df$baseline_cag - t, 0)
    fit <- fit_lmer_safe(slope_val ~ hs_term + (1 | group_id), work_df)
    if(!is.null(fit)) {
      cf <- fixef(fit); slope_val <- cf[length(cf)]
      # FIX: Added is.na() check
      if(!is.na(slope_val) && slope_val > 0) {
        this_aic <- AIC(fit)
        hockey_scan_df$AIC[i] <- this_aic
        if(this_aic < best_hockey_aic) { best_hockey_aic <- this_aic; best_hockey_thresh <- t; fit_hockey_best <- fit }
      }
    }
  }
  
  # --- 4. 3-Phase Free ---
  best_seg2_aic <- Inf; fit_seg2_best <- NULL; thresh_seg2_vec <- c(NA, NA)
  if(nrow(work_df) >= 10) { 
    for(t1 in scan_seq) {
      for(t2 in scan_seq) {
        if(t2 > (t1 + 4)) {
          work_df$d1 <- pmax(work_df$baseline_cag - t1, 0)
          work_df$d2 <- pmax(work_df$baseline_cag - t2, 0)
          
          fit <- tryCatch({ 
            suppressWarnings(suppressMessages(
              lmer(slope_val ~ baseline_cag + d1 + d2 + (1 | group_id), data = work_df, REML = FALSE, control = safe_ctrl)
            ))
          }, error = function(e) NULL)          
          if(is.null(fit)) fit <- tryCatch({ lm(slope_val ~ baseline_cag + d1 + d2, data = work_df) }, error = function(e) NULL)
          if(is.null(fit)) next
          
          if(inherits(fit, "lmerMod")) { if(fixef(fit)[1] < 0) fit <- lmer(slope_val ~ baseline_cag + d1 + d2 - 1 + (1 | group_id), data = work_df, REML = FALSE, control = safe_ctrl) } 
          else { if(coef(fit)[1] < 0) fit <- lm(slope_val ~ baseline_cag + d1 + d2 - 1, data = work_df) }
          
          cf <- if(inherits(fit, "lmerMod")) fixef(fit) else coef(fit)
          has_int <- "(Intercept)" %in% names(cf); idx <- if(has_int) 2 else 1
          s1 <- if("baseline_cag" %in% names(cf)) cf["baseline_cag"] else cf[idx]
          s_d1 <- if("d1" %in% names(cf)) cf["d1"] else cf[idx+1]
          s_d2 <- if("d2" %in% names(cf)) cf["d2"] else cf[idx+2]
          
          # FIX: Added is.na() checks
          if(!is.na(s1) && !is.na(s_d1) && !is.na(s_d2) && s1 >= 0 && (s1+s_d1) > s1 && (s1+s_d1+s_d2) > (s1+s_d1)) {
            if(AIC(fit) < best_seg2_aic) { best_seg2_aic <- AIC(fit); fit_seg2_best <- fit; thresh_seg2_vec <- c(t1, t2) }
          }
        }
      }
    }
  }
  
  # --- 5. Exponential Scan ---
  best_exp_aic <- Inf; best_exp_k <- NA; fit_exp_best <- NULL
  cag_min <- min(work_df$baseline_cag)
  
  for(k in scan_exp_k) {
    work_df$exp_term <- exp(k * (work_df$baseline_cag - cag_min))
    fit <- fit_lmer_safe(slope_val ~ exp_term + (1 | group_id), work_df)
    if(!is.null(fit)) {
      cf <- fixef(fit); b_val <- cf[length(cf)]
      # FIX: Added is.na() check
      if(!is.na(b_val) && b_val > 0) {
        if(AIC(fit) < best_exp_aic) { best_exp_aic <- AIC(fit); best_exp_k <- k; fit_exp_best <- fit }
      }
    }
  }
  
  # --- 6. Variable Power Scan ---
  best_pow_aic <- Inf; best_pow_p <- NA; fit_pow_best <- NULL
  for(p in scan_pow_p) {
    work_df$pow_term <- work_df$baseline_cag^p
    fit <- fit_lmer_safe(slope_val ~ pow_term + (1 | group_id), work_df)
    if(!is.null(fit)) {
      cf <- fixef(fit); b_val <- cf[length(cf)]
      # FIX: Added is.na() check
      if(!is.na(b_val) && b_val > 0) {
        if(AIC(fit) < best_pow_aic) { best_pow_aic <- AIC(fit); best_pow_p <- p; fit_pow_best <- fit }
      }
    }
  }
  
  # --- 9. Log-Linear Scan (Y-Transformed) ---
  min_val <- min(work_df$slope_val, na.rm=TRUE)
  y_shift <- if(min_val <= 0) abs(min_val) + 0.001 else 0
  
  work_df$log_y <- log(work_df$slope_val + y_shift)
  fit_loglin <- tryCatch({
    suppressWarnings(suppressMessages({
      lmer(log_y ~ baseline_cag + (1 | group_id), data = work_df, REML = FALSE, control = safe_ctrl)
    }))
  }, error = function(e) NULL)
  
  aic_loglin <- Inf
  r2_loglin <- NA
  
  if(!is.null(fit_loglin)) {
    log_preds <- predict(fit_loglin, newdata = work_df, re.form = NA)
    orig_preds <- exp(log_preds) - y_shift
    
    rss <- sum((work_df$slope_val - orig_preds)^2, na.rm = TRUE)
    n <- nrow(work_df)
    k <- length(fixef(fit_loglin)) + 1 
    
    aic_loglin <- n * log(rss/n) + 2 * k + n + n * log(2 * pi)
    r2_loglin <- cor(work_df$slope_val, orig_preds, use="complete.obs")^2
  }
  
  # --- 7. AIC Comparison ---
  aic_lin <- if(!is.null(fit_lin)) AIC(fit_lin) else Inf
  aic_quad <- if(!is.null(fit_quad)) AIC(fit_quad) else Inf
  aic_seg1 <- best_seg1_aic
  aic_hockey <- best_hockey_aic
  aic_seg2 <- best_seg2_aic 
  aic_exp <- best_exp_aic
  aic_pow <- best_pow_aic
  
  # --- 8. R-Squared ---
  r2_lin <- if(!is.null(fit_lin)) calc_population_r2(fit_lin, work_df) else NA
  r2_quad <- if(!is.null(fit_quad)) calc_population_r2(fit_quad, work_df) else NA
  r2_seg1 <- if(!is.null(fit_seg1_best)) calc_population_r2(fit_seg1_best, work_df) else NA
  r2_hockey <- if(!is.null(fit_hockey_best)) calc_population_r2(fit_hockey_best, work_df) else NA
  r2_seg2 <- if(!is.null(fit_seg2_best)) calc_population_r2(fit_seg2_best, work_df) else NA
  r2_exp <- if(!is.null(fit_exp_best)) calc_population_r2(fit_exp_best, work_df) else NA
  r2_pow <- if(!is.null(fit_pow_best)) calc_population_r2(fit_pow_best, work_df) else NA
  
  # --- 9. Determine Winner ---
  thresh_quad <- if(!is.null(fit_lin) && !is.null(fit_quad)) get_quad_crossing_point(fit_lin, fit_quad, work_df, divergence_percent) else NA
  
  best_model <- "Linear"; best_aic <- aic_lin; raw_thresh <- NA
  
  if (aic_quad < (best_aic - 2)) { best_model <- "Quadratic"; best_aic <- aic_quad; raw_thresh <- thresh_quad }
  if (aic_seg1 < (best_aic - 2)) { best_model <- "2-Phase (Free)"; best_aic <- aic_seg1; raw_thresh <- best_seg1_thresh }
  if (aic_seg2 < (best_aic - 2)) { best_model <- "3-Phase (Free)"; best_aic <- aic_seg2; raw_thresh <- thresh_seg2_vec[1] }
  #if (aic_hockey < (best_aic - 2)) { best_model <- "Hockey Stick"; best_aic <- aic_hockey; raw_thresh <- best_hockey_thresh }
  if (aic_exp < (best_aic - 2))  { best_model <- "Exponential"; best_aic <- aic_exp; raw_thresh <- NA }
  if (aic_pow < (best_aic - 2))  { best_model <- "Variable Power"; best_aic <- aic_pow; raw_thresh <- NA }
  if (aic_loglin < (best_aic - 2)) { best_model <- "Log-Linear (Y-Transformed)"; best_aic <- aic_loglin; raw_thresh <- NA }
  
  threshold_est <- NA
  if (!is.na(raw_thresh)) {
    if (raw_thresh >= min_valid && raw_thresh <= max_valid) threshold_est <- raw_thresh else best_model <- paste0(best_model, " (Edge Artifact Ignored)")
  }
  
  # --- 10. Bootstrapping ---
  threshold_ci <- c(NA, NA, NA)
  if (!is.na(threshold_est)) {
    message(paste0("...Best fit: ", best_model, ". Bootstrapping..."))
    boot_fn <- function(data, indices) {
      d_b <- data[indices, ]
      b_seq <- seq(min(d_b$baseline_cag, na.rm=T)+5, max(d_b$baseline_cag, na.rm=T)-5, length.out=15)
      tryCatch({
        val <- NA
        if (grepl("Quadratic", best_model)) {
          ml <- lmer(slope_val ~ baseline_cag + (1|group_id), data=d_b, REML=F, control=safe_ctrl)
          mq <- lmer(slope_val ~ baseline_cag + I(baseline_cag^2) + (1|group_id), data=d_b, REML=F, control=safe_ctrl)
          if(fixef(mq)[grep("I\\(", names(fixef(mq)))] > 0) val <- get_quad_crossing_point(ml, mq, d_b, divergence_percent)
        } else if (grepl("Hockey", best_model)) {
          best_a <- Inf; best_t <- NA
          for(tb in b_seq) {
            d_b$ht <- pmax(d_b$baseline_cag - tb, 0)
            fb <- tryCatch(lmer(slope_val ~ ht + (1|group_id), data=d_b, REML=F, control=safe_ctrl), error=function(e) NULL)
            if(!is.null(fb) && !is.na(fixef(fb)["ht"]) && fixef(fb)["ht"] > 0 && AIC(fb) < best_a) { best_a <- AIC(fb); best_t <- tb }
          }
          val <- best_t
        } else if (grepl("2-Phase", best_model)) {
          best_a <- Inf; best_t <- NA
          for(tb in b_seq) {
            d_b$td <- pmax(d_b$baseline_cag - tb, 0)
            fb <- tryCatch(lmer(slope_val ~ baseline_cag + td + (1|group_id), data=d_b, REML=F, control=safe_ctrl), error=function(e) NULL)
            if(!is.null(fb)) {
              cf <- fixef(fb); if(!is.na(cf[2]) && !is.na(cf[3]) && cf[2] >= 0 && (cf[2]+cf[3]) > cf[2] && AIC(fb) < best_a) { best_a <- AIC(fb); best_t <- tb }
            }
          }
          val <- best_t
        } 
        if (grepl("3-Phase", best_model)) val <- threshold_est 
        return(val)
      }, error = function(e) return(NA))
    }
    set.seed(123)
    boot_obj <- boot(data = work_df, statistic = boot_fn, R = n_boot)
    boot_vals <- na.omit(boot_obj$t)
    if(length(boot_vals) > 50) {
      threshold_ci <- quantile(boot_vals, probs = c(0.025, 0.5, 0.975))
      threshold_est <- median(boot_vals, na.rm = TRUE)
    }
  }
  
  # --- 11. Generate Lines ---
  x_seq <- seq(min(work_df$baseline_cag), max(work_df$baseline_cag), length.out = 200)
  
  # Initialize with Base Models
  lines_df <- rbind(
    data.frame(x = x_seq, y = predict(fit_lin, newdata = data.frame(baseline_cag=x_seq), re.form = NA), Model = "Linear"),
    data.frame(x = x_seq, y = predict(fit_quad, newdata = data.frame(baseline_cag=x_seq), re.form = NA), Model = "Quadratic")
  )
  
  if(!is.null(fit_seg1_best)) {
    cf <- fixef(fit_seg1_best); int <- if("(Intercept)" %in% names(cf)) cf["(Intercept)"] else 0
    s1 <- cf["baseline_cag"]; s_diff <- cf["term_diff"]
    lines_df <- rbind(lines_df, data.frame(x=x_seq, y=int + s1*x_seq + s_diff*pmax(x_seq - best_seg1_thresh, 0), Model="2-Phase (Free)"))
  }
  if(!is.null(fit_hockey_best)) {
    cf <- fixef(fit_hockey_best); int <- if("(Intercept)" %in% names(cf)) cf["(Intercept)"] else 0; s <- cf["hs_term"]
    lines_df <- rbind(lines_df, data.frame(x=x_seq, y=int + s*pmax(x_seq - best_hockey_thresh, 0), Model="Hockey Stick"))
  }
  if(!is.null(fit_seg2_best)) {
    cf <- fixef(fit_seg2_best); int <- if("(Intercept)" %in% names(cf)) cf["(Intercept)"] else 0
    lines_df <- rbind(lines_df, data.frame(x=x_seq, y=int + cf["baseline_cag"]*x_seq + cf["d1"]*pmax(x_seq - thresh_seg2_vec[1], 0) + cf["d2"]*pmax(x_seq - thresh_seg2_vec[2], 0), Model="3-Phase (Free)"))
  }
  if(!is.null(fit_exp_best)) {
    cf <- fixef(fit_exp_best); int <- if("(Intercept)" %in% names(cf)) cf["(Intercept)"] else 0; s <- cf["exp_term"]
    lines_df <- rbind(lines_df, data.frame(x=x_seq, y=int + s*exp(best_exp_k * (x_seq - cag_min)), Model="Exponential"))
  }
  if(!is.null(fit_pow_best)) {
    cf <- fixef(fit_pow_best); int <- if("(Intercept)" %in% names(cf)) cf["(Intercept)"] else 0; s <- cf["pow_term"]
    lines_df <- rbind(lines_df, data.frame(x=x_seq, y=int + s*(x_seq^best_pow_p), Model="Variable Power"))
  }
  
  if(!is.null(fit_loglin)) {
    cf <- fixef(fit_loglin)
    if(!exists("y_shift")) {
      min_v <- min(work_df$slope_val, na.rm=TRUE)
      y_shift <- if(min_v <= 0) abs(min_v) + 0.001 else 0
    }
    pred_vals <- exp(cf["(Intercept)"] + cf["baseline_cag"] * x_seq) - y_shift
    pred_vals[is.infinite(pred_vals)] <- NA
    lines_df <- rbind(lines_df, data.frame(x = x_seq, y = pred_vals, Model = "Log-Linear (Y-Transformed)"))
  }
  
  lines_df$Model <- as.character(lines_df$Model)
  
  # --- 12. Get Equation of Best Model ---
  winner_obj <- NULL
  if (grepl("Linear", best_model)) winner_obj <- fit_lin
  if (grepl("Quadratic", best_model)) winner_obj <- fit_quad
  if (grepl("2-Phase", best_model)) winner_obj <- fit_seg1_best
  if (grepl("3-Phase", best_model)) winner_obj <- fit_seg2_best
  if (grepl("Hockey", best_model)) winner_obj <- fit_hockey_best
  if (grepl("Exponential", best_model)) winner_obj <- fit_exp_best
  if (grepl("Power", best_model)) winner_obj <- fit_pow_best
  if (grepl("Log-Linear", best_model)) winner_obj <- fit_loglin
  
  best_equation <- get_model_equation(
    model_obj = winner_obj, 
    model_type = gsub(" \\(Edge Artifact Ignored\\)", "", best_model),
    params = list(exp_k = best_exp_k, pow_p = best_pow_p, y_shift = y_shift),
    thresholds = list(hockey = best_hockey_thresh, seg1 = best_seg1_thresh)
  )
  
  return(list(
    work_df = work_df, lines = lines_df, 
    best_model = best_model, 
    best_equation = best_equation, 
    threshold_est = threshold_est, threshold_ci = threshold_ci, 
    stats = list(AIC_lin = aic_lin, AIC_quad = aic_quad, AIC_seg1 = aic_seg1, AIC_seg2 = aic_seg2, AIC_hockey = aic_hockey, AIC_exp = aic_exp, AIC_pow = aic_pow, AIC_loglin = aic_loglin),
    r2 = list(lin = r2_lin, quad = r2_quad, seg1 = r2_seg1, seg2 = r2_seg2, hockey = r2_hockey, exp = r2_exp, pow = r2_pow,loglin = r2_loglin),
    thresholds = list(quad = thresh_quad, seg1 = best_seg1_thresh, seg2 = thresh_seg2_vec, hockey = best_hockey_thresh),
    params = list(exp_k = best_exp_k, pow_p = best_pow_p, y_shift = y_shift), 
    scans = list(hockey = hockey_scan_df, seg1 = seg1_scan_df),
    models = list(Linear=fit_lin, Quadratic=fit_quad, `2-Phase (Free)`=fit_seg1_best, `3-Phase (Free)`=fit_seg2_best, `Hockey Stick`=fit_hockey_best, Exponential=fit_exp_best, `Variable Power`=fit_pow_best, `Log-Linear (Y-Transformed)`=fit_loglin)
  ))
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# 1E. Helper: Extract Model Equation Strings
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
get_model_equation <- function(model_obj, model_type, params = NULL, thresholds = NULL) {
  if(is.null(model_obj)) return(NA)
  
  # Extract fixed effects (round for readability)
  cf <- fixef(model_obj)
  fmt <- function(x) formatC(x, format = "g", digits = 4)
  
  eq_str <- ""
  
  if (model_type == "Linear") {
    # y = Int + Slope*x
    eq_str <- paste0("Rate = ", fmt(cf["(Intercept)"]), " + ", fmt(cf["baseline_cag"]), " * CAG")
    
  } else if (model_type == "Quadratic") {
    # y = Int + b*x + a*x^2
    a_term <- cf[grep("I\\(", names(cf))]
    eq_str <- paste0("Rate = ", fmt(cf["(Intercept)"]), " + ", fmt(cf["baseline_cag"]), " * CAG + ", fmt(a_term), " * CAG^2")
    
  } else if (model_type == "Hockey Stick") {
    # y = Int + Slope * (x - Threshold)+
    t_val <- thresholds$hockey
    eq_str <- paste0("Rate = ", fmt(cf["(Intercept)"]), " + ", fmt(cf["hs_term"]), " * max(0, CAG - ", round(t_val, 1), ")")
    
  } else if (model_type == "2-Phase (Free)") {
    # y = Int + Slope1*x + SlopeDiff * (x - Threshold)+
    t_val <- thresholds$seg1
    eq_str <- paste0("Rate = ", fmt(cf["(Intercept)"]), " + ", fmt(cf["baseline_cag"]), " * CAG + ", 
                     fmt(cf["term_diff"]), " * max(0, CAG - ", round(t_val, 1), ")")
    
  } else if (model_type == "Exponential") {
    # y = Int + Coeff * exp(k * (x - min))
    k_val <- params$exp_k
    eq_str <- paste0("Rate = ", fmt(cf["(Intercept)"]), " + ", fmt(cf["exp_term"]), " * exp(", round(k_val, 4), " * (CAG - Min))")
    
  } else if (model_type == "Variable Power") {
    # y = Int + Coeff * x^p
    p_val <- params$pow_p
    eq_str <- paste0("Rate = ", fmt(cf["(Intercept)"]), " + ", fmt(cf["pow_term"]), " * CAG^", round(p_val, 2))
    
  } else if (model_type == "Log-Linear (Y-Transformed)") {
    # Model: ln(y + shift) = Int + Slope*x
    # Display as: y = exp(Int + Slope*x) - shift
    shift <- params$y_shift
    eq_str <- paste0("Rate = exp(", fmt(cf["(Intercept)"]), " + ", fmt(cf["baseline_cag"]), " * CAG) - ", round(shift, 4))
    
  } else {
    eq_str <- "Complex/Other"
  }
  
  return(eq_str)
}

#=============================================================================#
# PART 2 & 3: PLOTTING LOOPS                                              #####
#=============================================================================#

logr::log_print("\n--- Starting PART 2: Descriptive Correlation Plots ---")
cag_col <- if("mode" %in% names(data_per_pcr)) "mode" else response_vars[1]
baseline_map_agg <- data_per_pcr %>%
  group_by(!!sym(primary_var), !!sym(secondary_var)) %>%
  dplyr::filter(!!sym(time_var) == min(!!sym(time_var), na.rm = TRUE)) %>%
  summarise(baseline_cag = mean(!!sym(cag_col), na.rm = TRUE), .groups = "drop")

#change this to model the desired variables
cag_exp_plot_vars <- c("mode_change", "instability_index_change", "expansion_index_change")

variable_specific_exclusions <- list(
  # "response_variable_name" = c("Genotype_A", "Genotype_B"),
  #"mode_change" = c("190")
)
# --- Initialize Storage Lists ---
summary_results_list <- list() # For Table 1 (Best Model Only)
detailed_results_list <- list() # For Table 2 (All Models)


for (current_resp_var in cag_exp_plot_vars) {
  slope_col <- paste0("slope_", current_resp_var)
  if (!slope_col %in% colnames(all_group_slopes)) next
  
  y_label <- response_labels[[current_resp_var]] %||% current_resp_var
  short_name <- response_shortnames[[current_resp_var]] %||% current_resp_var
  slopes_agg <- all_group_slopes %>% 
    group_by(!!sym(primary_var), !!sym(secondary_var)) %>% 
    summarise(slope_val = mean(!!sym(slope_col), na.rm=TRUE), .groups="drop")
  
  # 1. Determine Exclusions for THIS specific variable
  # Start with global exclusions
  current_drop_list <- c("Control", "ExcludeMe") 
  
  # Add variable-specific exclusions if they exist in your list
  if (current_resp_var %in% names(variable_specific_exclusions)) {
    specific_drops <- variable_specific_exclusions[[current_resp_var]]
    current_drop_list <- c(current_drop_list, specific_drops)
    logr::log_print(paste("-> EXCLUDING specific genotypes for", current_resp_var, ":", paste(specific_drops, collapse = ", ")), console=TRUE)
  }
  
  plot_data_agg <- left_join(baseline_map_agg, slopes_agg, by = c(primary_var, secondary_var)) %>% 
    dplyr::filter(!is.na(slope_val)) %>%
    dplyr::filter(!.data[[primary_var]] %in% current_drop_list) %>%
    left_join(label_lookup, by = primary_var) # DYNAMIC LABELS
  
 
  
  if (nrow(plot_data_agg) > 3) {
    
  # 2. GENERATE GLOBAL CORRELATION PLOT
    p_global <- ggplot(plot_data_agg, aes(x = baseline_cag, y = slope_val)) +
      geom_point(aes(color = !!sym(plot_label_col)), size = 3, alpha = 0.8) +
      geom_smooth(method = "lm", se = FALSE, color="grey50", linetype="dashed") +
      labs(title = paste("Baseline vs Rate (Aggregated):", y_label), 
           x = "Baseline Modal CAG", y = "Expansion Rate", color = "Genotype") + 
      theme_publication()
    
    p_global <- apply_smart_palette(p_global, plot_label_col)
    ggsave(file.path(cag_exp_plot_dir, paste0("corr_global_", short_name, ".tiff")), 
           p_global, width = 10, height = 7, dpi = 300, compression = "lzw")
    
    # 3. GENERATE BY-GENOTYPE CORRELATION PLOT
    p_geno <- ggplot(plot_data_agg, aes(x = baseline_cag, y = slope_val, color = !!sym(plot_label_col))) +
      geom_point(size = 3) + 
      geom_smooth(method = "lm", se = FALSE) +
      labs(title = paste("Baseline vs Rate (By Genotype):", y_label), 
           x = "Baseline Modal CAG", y = "Expansion Rate", color = "Genotype") + 
      theme_publication()
    p_geno <- apply_smart_palette(p_geno, plot_label_col)
    
    ggsave(file.path(cag_exp_plot_dir, paste0("corr_geno_", short_name, ".tiff")), 
           p_geno, width = 10, height = 7, dpi = 300, compression = "lzw")
    
    # Print to viewer (per your preference)
    print(p_global)
    print(p_geno)
  }
}
#=========================================================================##
# Part 3: Robust Anaylsis ####
#=========================================================================##

logr::log_print("\n--- Starting PART 3: Robust Analysis Loop ---")
# --- FIX: Dynamic Key Selection ---
# We only want to join by Genotype and Clone. Because we use True BLUPs, 
# 'rep' and 'batch' do not exist in the slope table and should be ignored.
intended_keys_bio <- c(primary_var, secondary_var)
join_keys_bio <- intended_keys_bio[intended_keys_bio %in% colnames(all_group_slopes)]

baseline_map_bio <- summary_per_rep %>%
  group_by(across(all_of(join_keys_bio))) %>%
  dplyr::filter(!!sym(time_var) == min(!!sym(time_var), na.rm = TRUE)) %>%
  summarise(baseline_cag = mean(!!sym(cag_col), na.rm = TRUE), .groups = "drop")

stats_results_list <- list()

for (current_resp_var in cag_exp_plot_vars) {
  slope_col <- paste0("slope_", current_resp_var)
  if (!slope_col %in% colnames(all_group_slopes)) next
  y_label <- response_labels[[current_resp_var]] %||% current_resp_var
  short_name <- response_shortnames[[current_resp_var]] %||% current_resp_var
  
  if (length(join_keys_bio) > 0) {
    slopes_bio <- all_group_slopes %>% dplyr::select(all_of(join_keys_bio), slope_val = !!sym(slope_col))
    
    # 1. Determine Exclusions for THIS specific variable
    current_drop_list <- c("Control", "ExcludeMe") 
    
    if (current_resp_var %in% names(variable_specific_exclusions)) {
      specific_drops <- variable_specific_exclusions[[current_resp_var]]
      current_drop_list <- c(current_drop_list, specific_drops)
      logr::log_print(paste("-> EXCLUDING specific genotypes for", current_resp_var, ":", paste(specific_drops, collapse = ", ")))
    }
    
    # 2. Filter the Data
    plot_data_bio <- left_join(baseline_map_bio, slopes_bio, by = join_keys_bio) %>% 
      dplyr::filter(!is.na(slope_val)) %>% 
      dplyr::filter(!.data[[primary_var]] %in% current_drop_list) %>% 
      left_join(label_lookup, by = primary_var) %>%
      restore_all_factors()
    
    if (nrow(plot_data_bio) > 10) { 
      logr::log_print(paste("Analyzing:", current_resp_var))
      res_robust <- get_robust_mixed_analysis(df = plot_data_bio, x_col = "baseline_cag", y_col = "slope_val", group_var = primary_var, n_boot = 500, boundary_buffer = 10)
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
      # DATA EXPORT: Create Two Separate Tables
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
      
      winner_r2 <- NA
      if(grepl("Linear", res_robust$best_model)) winner_r2 <- res_robust$r2$lin
      if(grepl("Quadratic", res_robust$best_model)) winner_r2 <- res_robust$r2$quad
      if(grepl("2-Phase", res_robust$best_model)) winner_r2 <- res_robust$r2$seg1
      if(grepl("3-Phase", res_robust$best_model)) winner_r2 <- res_robust$r2$seg2
      if(grepl("Hockey", res_robust$best_model)) winner_r2 <- res_robust$r2$hockey
      if(grepl("Exponential", res_robust$best_model)) winner_r2 <- res_robust$r2$exp
      if(grepl("Power", res_robust$best_model)) winner_r2 <- res_robust$r2$pow
      if(grepl("Log-Linear", res_robust$best_model)) winner_r2 <- res_robust$r2$loglin
      
      summary_entry <- data.frame(
        Variable = current_resp_var,
        Winning_Model = res_robust$best_model,
        Equation = res_robust$best_equation,
        Onset_Threshold = res_robust$threshold_est, 
        CI_Lower_95 = res_robust$threshold_ci[1],
        CI_Upper_95 = res_robust$threshold_ci[3],
        R_Squared = winner_r2,
        AIC_Score = min(unlist(res_robust$stats), na.rm=TRUE) 
      )
      summary_results_list[[current_resp_var]] <- summary_entry
      
      detailed_entry <- data.frame(
        Variable = current_resp_var,
        AIC_Linear = res_robust$stats$AIC_lin,
        AIC_Quadratic = res_robust$stats$AIC_quad,
        AIC_2Phase = res_robust$stats$AIC_seg1,
        AIC_3Phase = res_robust$stats$AIC_seg2,
        AIC_Hockey = res_robust$stats$AIC_hockey,
        AIC_Exp = res_robust$stats$AIC_exp,
        AIC_Power = res_robust$stats$AIC_pow,
        AIC_LogLin = res_robust$stats$AIC_loglin,
        R2_Linear = res_robust$r2$lin,
        R2_Quadratic = res_robust$r2$quad,
        R2_2Phase = res_robust$r2$seg1,
        R2_3Phase = res_robust$r2$seg2,
        R2_Hockey = res_robust$r2$hockey,
        R2_Exp = res_robust$r2$exp,
        R2_Power = res_robust$r2$pow,
        R2_LogLin = res_robust$r2$loglin,
        Thresh_Quad_Candidate = res_robust$thresholds$quad,
        Thresh_2Ph_Candidate = res_robust$thresholds$seg1,
        Thresh_Hockey_Candidate = res_robust$thresholds$hockey,
        Param_Exp_k = res_robust$params$exp_k,
        Param_Power_p = res_robust$params$pow_p
      )
      detailed_results_list[[current_resp_var]] <- detailed_entry
      
      subtitle_stats <- paste0("Best Fit: ", res_robust$best_model)
      if (!is.na(res_robust$threshold_est)) {
        subtitle_stats <- paste0(subtitle_stats, "\nOnset Threshold: ", round(res_robust$threshold_est, 1), " CAGs")
        if (!is.na(res_robust$threshold_ci[1])) subtitle_stats <- paste0(subtitle_stats, " [95% CI: ", round(res_robust$threshold_ci[1], 1), "-", round(res_robust$threshold_ci[3], 1), "]")
      }
      
      fmt_aic <- function(val) if(is.infinite(val)) "NA" else round(val, 1)
      fmt_r2  <- function(val) if(is.na(val)) "NA" else round(val, 2)
      
      caption_text <- paste0(
        "AIC (R2): Lin=", fmt_aic(res_robust$stats$AIC_lin), "(", fmt_r2(res_robust$r2$lin), ")",
        " | 2-Ph=", fmt_aic(res_robust$stats$AIC_seg1), "(", fmt_r2(res_robust$r2$seg1), ")",
        "\n3-Ph=", fmt_aic(res_robust$stats$AIC_seg2), "(", fmt_r2(res_robust$r2$seg2), ")",
        " | Hockey=", fmt_aic(res_robust$stats$AIC_hockey), "(", fmt_r2(res_robust$r2$hockey), ")",
        "\nExp=", fmt_aic(res_robust$stats$AIC_exp), "(", fmt_r2(res_robust$r2$exp), ")",
        " | Quad=", fmt_aic(res_robust$stats$AIC_quad), "(", fmt_r2(res_robust$r2$quad), ")",
        " | Power=", fmt_aic(res_robust$stats$AIC_pow), "(", fmt_r2(res_robust$r2$pow), ")",
        " | LogLin=", fmt_aic(res_robust$stats$AIC_loglin), "(", fmt_r2(res_robust$r2$loglin), ")"
      )
      clean_best_model <- gsub(" (Edge Artifact Ignored)", "", res_robust$best_model, fixed = TRUE)
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
      # THE PALETTE FIX (Calculated ONCE for all 4 plots)
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
      n_groups_bio <- length(unique(plot_data_bio[[plot_label_col]]))
      safe_pal_bio <- tryCatch({
        create_custom_palette(plot_data_bio, config, primary_var)
      }, error = function(e) { 
        n_groups <- length(unique(plot_data_bio[[primary_var]]))
        RColorBrewer::brewer.pal(min(n_groups, 9), "Set1")[1:n_groups]
        })

      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
      # PLOT 2: PUBLICATION VIEW (WINNER)
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
      p_clean <- ggplot(plot_data_bio, aes(x = baseline_cag, y = slope_val)) +
        {if(!is.na(res_robust$threshold_est)) list(
          annotate("rect", xmin = res_robust$threshold_ci[1], xmax = res_robust$threshold_ci[3], 
                   ymin = -Inf, ymax = Inf, fill = "red", alpha = 0.1),
          geom_vline(xintercept = res_robust$threshold_est, linetype = "dashed", color = "red", linewidth = 0.5)
        )} +
        geom_point(aes(color = !!sym(plot_label_col)), size = 3, alpha = 0.6) + 
        geom_line(data = res_robust$lines %>% dplyr::filter(Model == clean_best_model), aes(x = x, y = y), color = "black", linewidth = 1.2) +
        labs(title = paste("Instability Dynamics:", y_label), subtitle = subtitle_stats, x = "Baseline Modal CAG", y = "Rate", caption = caption_text, color = str_to_title(primary_var)) +
        theme_publication()
      
      p_clean <- apply_smart_palette(p_clean, plot_label_col)
      
      ggsave(file.path(cag_exp_plot_dir, paste0("final_model_", short_name, ".tiff")), p_clean, width = 10, height = 8, compression = "lzw")
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
      # PLOT 1: COMPARISON VIEW
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
      # 1. Base plot and Points (Colored by Genotype)
      p_main_base <- ggplot(plot_data_bio, aes(x = baseline_cag, y = slope_val)) +
        geom_point(aes(color = !!sym(plot_label_col)), size = 3, alpha = 0.6) + 
        labs(
          title = paste("Instability Dynamics:", y_label), 
          subtitle = subtitle_stats, 
          x = "Baseline Modal CAG", 
          y = "Rate", 
          caption = caption_text
        ) +
        theme_publication() + 
        theme(legend.position = "right")
      
      # 2. Apply your smart palette for the points
      p_main_points <- apply_smart_palette(p_main_base, plot_label_col)
      
      # 3. Add a NEW color scale, then add the Lines (Colored by Model)
      p_main <- p_main_points + 
        ggnewscale::new_scale_color() + # <--- THE MAGIC LINE
        geom_line(data = res_robust$lines, 
                  aes(x = x, y = y, group = Model, color = Model, 
                      linewidth = Model == res_robust$best_model, 
                      linetype = Model == res_robust$best_model)) +
        
        # 4. Define the colors specifically for the new Model scale
        scale_color_manual(
          values = c("Linear" = "black", "Quadratic" = "forestgreen", 
                     "2-Phase (Free)" = "blue", "3-Phase (Free)" = "orange", 
                     "Hockey Stick" = "purple", "Exponential" = "magenta", 
                     "Variable Power" = "darkcyan", "Log-Linear (Y-Transformed)" = "brown4"), 
          name = "Model"
        ) +
        scale_linewidth_manual(values = c("TRUE" = 1.2, "FALSE" = 0.8), guide="none") +
        scale_linetype_manual(values = c("TRUE" = "solid", "FALSE" = "dotted"), guide="none")
      
      print(p_main)
      
      p_scan <- generate_aic_profile_plot(
        res_robust$scans$hockey, 
        res_robust$scans$seg1, 
        res_robust$stats$AIC_lin, 
        res_robust$stats$AIC_quad,    
        res_robust$stats$AIC_seg2, 
        res_robust$stats$AIC_exp,     
        res_robust$stats$AIC_pow,     
        res_robust$stats$AIC_loglin
      )
      
      p_combined <- if (!is.null(p_scan)) p_main + p_scan + plot_layout(widths = c(1.5, 1)) else p_main
      ggsave(file.path(cag_exp_plot_dir, paste0("model-comp_AIC_", short_name, ".tiff")), p_combined, width = 16, height = 8, compression = "lzw")
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
      # PLOT 3: SUMMARY VIEW (Genotype Means + SEM + Trendline)
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
      summary_plot_data <- plot_data_bio %>%
        dplyr::group_by(!!sym(primary_var), !!sym(plot_label_col)) %>%
        dplyr::summarise(
          mean_cag = mean(baseline_cag, na.rm = TRUE),
          mean_rate = mean(slope_val, na.rm = TRUE),
          se_rate = sd(slope_val, na.rm = TRUE) / sqrt(n()),
          n = n(),
          .groups = "drop"
        )
      
      p_summary <- ggplot(summary_plot_data, aes(x = mean_cag, y = mean_rate)) +
        geom_line(data = res_robust$lines %>% dplyr::filter(Model == clean_best_model), 
                  aes(x = x, y = y), color = "black", linewidth = 1.0, inherit.aes = FALSE) +
        geom_errorbar(aes(ymin = mean_rate - se_rate, ymax = mean_rate + se_rate, color = !!sym(plot_label_col)), 
                      width = 2, linewidth = 0.8, alpha = 0.8) +
        geom_point(aes(color = !!sym(plot_label_col)), size = 4) +
        {if(!is.na(res_robust$threshold_est)) list(
          geom_vline(xintercept = res_robust$threshold_est, linetype = "dashed", color = "red", linewidth = 0.5)
        )} +
        labs(
          title = paste("Instability Dynamics (Group Means):", y_label), 
          subtitle = paste0("Points represent Genotype Average +/- SEM.\n", subtitle_stats), 
          x = "Mean Baseline Modal CAG", 
          y = "Mean Rate", 
          caption = caption_text,
          color = str_to_title(primary_var)
        ) +
        theme_publication() + theme(legend.position = "right")
      
      p_summary <- apply_smart_palette(p_summary, plot_label_col)
      print(p_summary)
      ggsave(file.path(cag_exp_plot_dir, paste0("summary_means_", short_name, ".tiff")), p_summary, width = 10, height = 7, compression = "lzw")
      
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
      # LOOP: ACTUAL vs PREDICTED + ALTERNATIVE CLEAN PLOTS
      # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
      model_names <- c("Linear", "Quadratic", "2-Phase (Free)", "3-Phase (Free)", "Hockey Stick", "Exponential", "Variable Power", "Log-Linear (Y-Transformed)")
      
      for(m_name in model_names) {
        m_obj <- res_robust$models[[m_name]]
        if(!is.null(m_obj)) {
          
          temp_df <- res_robust$work_df
          
          if(m_name == "2-Phase (Free)") {
            t_val <- res_robust$thresholds$seg1
            temp_df$term_diff <- pmax(temp_df$baseline_cag - t_val, 0)
          } else if(m_name == "Hockey Stick") {
            t_val <- res_robust$thresholds$hockey
            temp_df$hs_term <- pmax(temp_df$baseline_cag - t_val, 0)
          } else if(m_name == "3-Phase (Free)") {
            t1 <- res_robust$thresholds$seg2[1]; t2 <- res_robust$thresholds$seg2[2]
            temp_df$d1 <- pmax(temp_df$baseline_cag - t1, 0)
            temp_df$d2 <- pmax(temp_df$baseline_cag - t2, 0)
          } else if(m_name == "Exponential") {
            k_val <- res_robust$params$exp_k; c_min <- min(temp_df$baseline_cag)
            temp_df$exp_term <- exp(k_val * (temp_df$baseline_cag - c_min))
          } else if(m_name == "Variable Power") {
            p_val <- res_robust$params$pow_p
            temp_df$pow_term <- temp_df$baseline_cag^p_val
          } 
          
          if (m_name == "Log-Linear (Y-Transformed)") {
            log_pred <- tryCatch(predict(m_obj, newdata=temp_df), error=function(e) NA)
            shift <- res_robust$params$y_shift
            temp_df$pred_val <- exp(log_pred) - shift
          } else {
            temp_df$pred_val <- tryCatch(predict(m_obj, newdata=temp_df), error=function(e) NA)
          } 
          
          if(!all(is.na(temp_df$pred_val))) {
            r2_val <- cor(temp_df$pred_val, temp_df$slope_val, use="complete.obs")^2
            p_diag <- ggplot(temp_df, aes(x = pred_val, y = slope_val)) +
              geom_point(aes(color = group_id), size = 3, alpha = 0.7) +
              geom_abline(intercept = 0, slope = 1, linetype = "dashed", color = "black") +
              labs(title = paste0("Actual vs. Predicted: ", m_name), 
                   subtitle = paste0("R2 = ", round(r2_val, 2)),
                   x = "Predicted Rate", y = "Actual Rate", color = str_to_title(primary_var)) +
              theme_publication()
            safe_name <- gsub(" ", "", gsub("[()]", "", m_name))
            ggsave(file.path(diag_plot_dir, paste0("diag_", safe_name, "_", short_name, ".tiff")), p_diag, width = 6, height = 6, compression = "lzw")
          }
          
          p_alt <- ggplot(plot_data_bio, aes(x = baseline_cag, y = slope_val)) +
            geom_point(aes(color = !!sym(plot_label_col)), size = 3, alpha = 0.6) +
            geom_line(data = res_robust$lines %>% dplyr::filter(Model == m_name), aes(x = x, y = y), color = "black", linewidth = 1.2) +
            labs(title = paste("Instability Dynamics:", y_label), subtitle = paste("Model:", m_name), x = "Baseline Modal CAG", y = "Rate", 
                 color = str_to_title(primary_var)) +
            theme_publication()
          p_alt <- apply_smart_palette(p_alt, primary_var)
          
          safe_name <- gsub(" ", "", gsub("[()]", "", m_name))
          ggsave(file.path(alternatives_dir, paste0("final_alt_", safe_name, "_", short_name, ".tiff")), p_alt, width = 10, height = 8, compression = "lzw")
        }
      }
    }
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# FINAL EXPORT: Save the Two Tables
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
if (length(summary_results_list) > 0) {
  
  # 1. Save Best Model Summary
  df_summary <- dplyr::bind_rows(summary_results_list)
  write.csv(df_summary, file.path(corr_plot_dir, "Summary_Best_Models.csv"), row.names = FALSE)
  logr::log_print("Saved Summary_Best_Models.csv")
  
  # 2. Save Detailed Diagnostics
  df_details <- dplyr::bind_rows(detailed_results_list)
  write.csv(df_details, file.path(corr_plot_dir, "Detailed_Model_Stats.csv"), row.names = FALSE)
  logr::log_print("Saved Detailed_Model_Stats.csv")
}

#=============================================================================#
# PART 1B: (NEW) FIT EXPLORATORY NON-LINEAR MODELS
#=============================================================================#

logr::log_print("\n--- Starting PART 1B: Fitting Exploratory Non-Linear Models ---", console = TRUE)

model_comparison_list <- list()
plot_list_poly <- list()
exploratory_models <- list()

# Loop through each response variable (e.g., mode_change, instability_index)
for (resp_var in response_vars) {
  logr::log_print(paste("...starting model comparison for:", resp_var))
  
  # 1. Get the original model object and data from Script 2
  model_obj_linear_reml <- all_model_outputs[[resp_var]]$model_object
  
  if (is.null(model_obj_linear_reml)) {
    logr::log_print(paste("...skipping", resp_var, "(no model object found)."))
    next
  }
  
  # Get the *exact* data frame used for the original model
  is_lmer <- inherits(model_obj_linear_reml, "lmerMod")
  
  clean_data <- if (is_lmer) {
    model_obj_linear_reml@frame # Access S4 slot for lmerMod
  } else {
    model_obj_linear_reml$model # Access S3 element for lm
  }
  
  if (is.null(clean_data)) {
    logr::log_print(paste("...skipping", resp_var, "(could not extract model data frame)."))
    next
  }
  
  orig_formula <- formula(model_obj_linear_reml)
  
  # Define the fit function and arguments for ML (for comparison)
  fit_func <- if (is_lmer) lmerTest::lmer else stats::lm
  fit_args_ml <- if (is_lmer) list(data = clean_data, REML = FALSE, na.action = na.omit) else list(data = clean_data, na.action = na.omit)
  
  # Define arguments for REML (for plotting predictions)
  fit_args_reml <- if (is_lmer) list(data = clean_data, REML = TRUE, na.action = na.omit) else list(data = clean_data, na.action = na.omit)
  
  tryCatch({
    # 2. Define the *three* formulas
    formula_str_linear <- deparse(orig_formula)
    formula_str_poly <- gsub(time_var, paste0("poly(", time_var, ", 2)"), formula_str_linear, fixed = TRUE)
    # Use log(time + 1) to handle Day 0
    formula_str_log_time <- gsub(time_var, paste0("log(", time_var, " + 1)"), formula_str_linear, fixed = TRUE)
    
    # 3. Fit all three models using Maximum Likelihood (REML = FALSE)
    logr::log_print("......fitting Linear (ML)")
    model_linear_ml <- do.call(fit_func, c(list(formula = as.formula(formula_str_linear)), fit_args_ml))
    
    logr::log_print("......fitting Polynomial (ML)")
    model_poly_ml <- do.call(fit_func, c(list(formula = as.formula(formula_str_poly)), fit_args_ml))
    
    logr::log_print("......fitting Log-Time (ML)")
    model_log_time_ml <- do.call(fit_func, c(list(formula = as.formula(formula_str_log_time)), fit_args_ml))
    
    # Store fitted models for Part 1D
    exploratory_models[[resp_var]] <- list(
      Linear = model_linear_ml,
      Polynomial = model_poly_ml,
      `Log-Time` = model_log_time_ml
    )
    
    # 4. Compare models
    logr::log_print("......extracting AIC/BIC for model comparison.")
    
    models_to_compare <- list(
      Linear = model_linear_ml, 
      Polynomial = model_poly_ml, 
      `Log-Time` = model_log_time_ml
    )
    
    model_comp_df <- purrr::map_df(models_to_compare, ~{
      tibble::tibble(
        AIC = stats::AIC(.x),
        BIC = stats::BIC(.x),
        logLik = as.numeric(stats::logLik(.x))
      )
    }, .id = "model_name") %>%
      dplyr::mutate(response_variable = resp_var) %>%
      dplyr::select(response_variable, model_name, AIC, BIC, logLik)
    
    model_comparison_list[[resp_var]] <- model_comp_df
    
    # 5. (FOR PLOTTING) Re-fit Polynomial model with REML = TRUE
    logr::log_print("......re-fitting Polynomial (REML) for plotting.")
    model_poly_reml <- do.call(fit_func, c(list(formula = as.formula(formula_str_poly)), fit_args_reml))
    
    # 6. Generate prediction plots for the Polynomial model
    pred_terms <- c(paste0(time_var, "[all]"), fixed_effects_to_plot) # fixed_effects_to_plot is from .RData
    model_preds_poly_raw <- ggeffects::ggpredict(model_poly_reml, terms = pred_terms)
    
    model_preds_poly <- as.data.frame(model_preds_poly_raw)
    if ("group" %in% colnames(model_preds_poly_raw)) {
      model_preds_poly <- model_preds_poly %>% 
        dplyr::rename(!!rlang::sym(fixed_effects_to_plot[[1]]) := group)
    }
    if ("facet" %in% colnames(model_preds_poly_raw)) {
      model_preds_poly <- model_preds_poly %>%
        dplyr::rename(!!rlang::sym(fixed_effects_to_plot[[2]]) := facet)
    }
    
    y_axis_label <- response_labels[[resp_var]] %||% stringr::str_to_title(resp_var)
    short_resp_var <- response_shortnames[[resp_var]] %||% resp_var
    
    p_poly_preds <- ggplot() +
      geom_point(
        data = modeling_data, # 'modeling_data' from .RData
        aes(
          x = !!sym(time_var),
          y = !!sym(resp_var),
          color = if(has_pseudo_clones) clone_rank else !!sym(color_var)
        ),
        alpha = 0.6, size = 2.5
      ) +
      { if (!is.null(custom_palette)) {
        scale_color_manual(values = custom_palette, name = if(has_pseudo_clones) "Clone Rank" else stringr::str_to_title(color_var))
      }} +
      ggnewscale::new_scale_color() +
      ggnewscale::new_scale_fill() +
      geom_line(
        data = model_preds_poly,
        aes(x = x, y = predicted, color = !!sym(primary_var)),
        linewidth = 1.2
      ) +
      geom_ribbon(
        data = model_preds_poly,
        aes(x = x, ymin = conf.low, ymax = conf.high, fill = !!sym(primary_var)),
        alpha = 0.2
      ) +
      { if (!is.null(x_var_palette)) {
        c(scale_color_manual(values = x_var_palette, name = "Genotype"),
          scale_fill_manual(values = x_var_palette, name = "Genotype"))
      }} +
      facet_wrap(vars(!!sym(primary_var)), scales = "free_y", axes = "all") +
      labs(
        title = paste("Exploratory Polynomial Model Fit:", y_axis_label),
        subtitle = "Lines show the non-linear (polynomial) model trend.",
        x = stringr::str_to_title(time_var),
        y = y_axis_label
      ) +
      theme_publication(base_size = 14) +
      theme(legend.position = "none")
    
    # Add legend
    if(is_clone_style_analysis && !is.null(p_legend)) {
      p_poly_preds <- p_poly_preds + p_legend + patchwork::plot_layout(widths = c(3, 1))
      plot_width_preds <- 16
    } else {
      p_poly_preds <- p_poly_preds + theme(legend.position = "right")
      plot_width_preds <- 12
    }
    
    # Per user preference: Print plot to R viewer
    print(p_poly_preds)
    
    # Per user preference: Save as TIFF to nested folder
    ggsave(
      file.path(corr_plot_dir, paste0("p_poly_fit_", short_resp_var, ".tiff")),
      p_poly_preds, width = plot_width_preds, height = 8, dpi = 300, device = 'tiff',
      compression = "lzw"
    )
    
  }, error = function(e) {
    logr::log_print(paste("...ERROR fitting exploratory models for", resp_var, ":", e$message))
  })
} # End of response variable loop


#=============================================================================#
# PART 1C: (NEW) COMPARE MODEL FIT (LINEAR VS. POLY VS. LOG-TIME)
#=============================================================================#
logr::log_print("\n--- Starting PART 1C: Comparing Model Fit Statistics ---", console = TRUE)

if (length(model_comparison_list) > 0) {
  # 1. Combine all comparison tables and save to Excel
  all_model_comps <- dplyr::bind_rows(model_comparison_list)
  
  model_export_list_fit <- list(Model_Fit_Comparison = all_model_comps)
  
  # --- (MODIFIED) Shortened file name ---
  openxlsx::write.xlsx(
    model_export_list_fit, 
    file.path(latest_analysis_dir, "Model_Fit_Comparison.xlsx"),
    rowNames = FALSE,
    overwrite = TRUE 
  )
  logr::log_print("...model fit comparison table saved to 'Model_Fit_Comparison.xlsx'.")
  
  # 2. Create a plot comparing AIC/BIC
  all_model_comps_long <- all_model_comps %>%
    dplyr::select(response_variable, model_name, AIC, BIC) %>%
    tidyr::pivot_longer(
      cols = c(AIC, BIC),
      names_to = "metric",
      values_to = "value"
    ) %>%
    dplyr::group_by(response_variable, metric) %>%
    dplyr::mutate(
      is_best = (value == min(value, na.rm = TRUE)),
      delta = value - min(value, na.rm = TRUE)
    ) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(
      # --- (MODIFIED) Removed "Factor" ---
      model_label = factor(model_name, levels = c("Linear", "Polynomial", "Log-Time"))
    )
  
  p_model_fit_comp <- ggplot(all_model_comps_long, 
                             aes(x = model_label, 
                                 y = delta, 
                                 fill = model_label)) +
    geom_col(color = "black", alpha = 0.8) +
    geom_text(
      data = . %>% dplyr::filter(is_best == TRUE),
      aes(label = "Best Fit"),
      vjust = -0.5, size = 3, fontface = "bold"
    ) +
    facet_grid(
      rows = vars(response_variable), 
      cols = vars(metric), 
      scales = "free_y"
    ) +
    labs(
      title = "Exploratory Model Fit Comparison (Lower is Better)",
      subtitle = "Shows the Delta AIC and Delta BIC (difference from the best-fitting model).",
      x = "Model Type",
      y = "Delta (Value - Best Value)",
      fill = "Model Type"
    ) +
    theme_publication() +
    theme(
      legend.position = "none",
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  # Per user preference: Print plot to R viewer
  print(p_model_fit_comp)
  
  # Per user preference: Save as TIFF to nested folder
  ggsave(
    file.path(corr_plot_dir, "p_exploratory_model_fit.tiff"),
    p_model_fit_comp, width = 8, height = 10, dpi = 300, compression = "lzw"
  )
  
} else {
  logr::log_print("...no model comparison data to plot.")
}


#=============================================================================#
# PART 1D: (NEW) VISUALIZE ALL MODEL FITS
#=============================================================================#
logr::log_print("\n--- Starting PART 1D: Visualizing All 3 Model Fits ---", console = TRUE)

# This uses the 'exploratory_models' list fitted in Part 1B

for (resp_var in names(exploratory_models)) {
  logr::log_print(paste("...generating 3-model comparison plot for:", resp_var))
  
  # 1. Get the model objects from the list
  models_to_plot <- exploratory_models[[resp_var]]
  
  # --- (MODIFIED) Removed "Factor" ---
  if (is.null(models_to_plot$Linear) || is.null(models_to_plot$Polynomial) || 
      is.null(models_to_plot$`Log-Time`)) {
    logr::log_print(paste("...skipping plot for", resp_var, "(one or more models failed in Part 1B)."))
    next
  }
  
  # 2. Define prediction terms
  pred_terms <- c(paste0(time_var, "[all]"), fixed_effects_to_plot)
  
  # 3. Generate predictions for all three models
  preds_linear <- ggeffects::ggpredict(models_to_plot$Linear, terms = pred_terms) %>%
    as.data.frame() %>%
    mutate(model_type = "Linear")
  
  preds_poly <- ggeffects::ggpredict(models_to_plot$Polynomial, terms = pred_terms) %>%
    as.data.frame() %>%
    mutate(model_type = "Polynomial")
  
  preds_log_time <- ggeffects::ggpredict(models_to_plot$`Log-Time`, terms = pred_terms) %>%
    as.data.frame() %>%
    mutate(model_type = "Log-Time")
  
  # 4. Combine prediction data
  # --- (MODIFIED) Removed "Factor" ---
  all_preds <- dplyr::bind_rows(preds_linear, preds_poly, preds_log_time) %>%
    # Rename 'group' and 'facet' to the real variable names
    { if ("group" %in% colnames(.)) rename(., !!rlang::sym(fixed_effects_to_plot[[1]]) := group) else . } %>%
    { if ("facet" %in% colnames(.)) rename(., !!rlang::sym(fixed_effects_to_plot[[2]]) := facet) else . } %>%
    mutate(
      model_type = factor(model_type, levels = c("Linear", "Polynomial", "Log-Time"))
    )
  
  # 5. Get labels
  y_axis_label <- response_labels[[resp_var]] %||% stringr::str_to_title(resp_var)
  short_resp_var <- response_shortnames[[resp_var]] %||% resp_var
  
  # 6. Create the plot
  p_all_fits <- ggplot(all_preds, aes(x = x, y = predicted, color = model_type)) +
    # Add the raw data points in the background
    geom_point(
      data = modeling_data,
      aes(x = !!sym(time_var), y = !!sym(resp_var)),
      color = "grey70", alpha = 0.5, inherit.aes = FALSE
    ) +
    # Add prediction lines
    geom_line(aes(linetype = model_type), linewidth = 1.1) +
    # Add ribbons (confidence intervals)
    geom_ribbon(aes(ymin = conf.low, ymax = conf.high, fill = model_type), 
                alpha = 0.2, linetype = 0) +
    # Facet by the main experimental group
    facet_wrap(vars(!!sym(primary_var)), scales = "free_y") +
    # --- (MODIFIED) Removed "Factor" ---
    scale_linetype_manual(values = c("Linear" = "solid", "Polynomial" = "dashed", 
                                     "Log-Time" = "dotdash")) +
    scale_color_brewer(palette = "Dark2") +
    scale_fill_brewer(palette = "Dark2") +
    labs(
      title = paste("Model Fit Comparison:", y_axis_label),
      subtitle = "Compares all 3 model fits against raw data (grey points).",
      x = stringr::str_to_title(time_var),
      y = y_axis_label,
      color = "Model Type",
      fill = "Model Type",
      linetype = "Model Type"
    ) +
    theme_publication() +
    theme(legend.position = "bottom")
  
  # Per user preference: Print plot to R viewer
  print(p_all_fits)
  
  # Per user preference: Save as TIFF to nested folder
  ggsave(
    file.path(corr_plot_dir, paste0("p_all_model_fits_", short_resp_var, ".tiff")),
    p_all_fits, width = 12, height = 9, dpi = 300, compression = "lzw"
  )
} # End of response variable loop

#=============================================================================#
# PART 2: CREATE BASE 'BY-GROUP' DATA
#=============================================================================#
logr::log_print("\n--- Starting PART 2: Creating Base 'By-Group' Data ---", console = TRUE)

# 1. 'all_group_slopes' is already loaded from the .RData file
logr::log_print("...using 'all_group_slopes' table loaded from .RData file.")

# 2. Join baseline and all slopes to create the "base" master data (BY GROUP)
# FIX: Use 'baseline_map_agg' which was created in Part 1
if(exists("baseline_map_agg") && nrow(baseline_map_agg) > 0) {
  
  # Ensure we are joining on the correct columns
  # baseline_map_agg usually has [primary_var, secondary_var, baseline_cag]
  # all_group_slopes has [primary_var, secondary_var, slope_mode, slope_ii...]
  
  join_cols <- c(primary_var, secondary_var)
  # Filter join columns to only those present in both
  valid_join <- join_cols[join_cols %in% colnames(baseline_map_agg) & join_cols %in% colnames(all_group_slopes)]
  
  master_data_by_group <- baseline_map_agg %>%
    left_join(all_group_slopes, by = valid_join)
  
  logr::log_print(paste("...base 'by_group' correlation data created (joined on:", paste(valid_join, collapse=", "), ")"))
  
} else {
  logr::log_print("...no baseline data found from Part 1, using 'all_group_slopes' as base.")
  master_data_by_group <- all_group_slopes
}
master_data_by_group <- apply_renaming_and_factors(master_data_by_group, config)
#=============================================================================#
# PART 2A: (NEW) LOAD AND JOIN GROUP-LEVEL METADATA
#=============================================================================#
logr::log_print("\n--- Starting PART 2A: Loading Group-Level Metadata ---", console = TRUE)

meta_cfg <- config$external_metadata

if (!is.null(meta_cfg$group_level_path) && meta_cfg$group_level_path != 'null') {
  group_meta_path <- here::here(meta_cfg$group_level_path)
  if (file.exists(group_meta_path)) {
    logr::log_print(paste("...Found optional group-level metadata. Loading:", basename(group_meta_path)))
    external_metadata_group <- readxl::read_excel(group_meta_path)
    
    if (all(c(primary_var, grouping_var) %in% colnames(external_metadata_group))) {
      master_data_by_group <- master_data_by_group %>%
        left_join(external_metadata_group, by = c(primary_var, grouping_var))
      logr::log_print("...Group-level metadata successfully joined to 'by_group' data.")
    } else {
      logr::log_print(paste("...WARNING: Group-level metadata skipped. File MUST contain columns:", primary_var, "and", grouping_var))
    }
  } else {
    logr::log_print(paste("...WARNING: Group-level metadata file not found at:", group_meta_path))
  }
} else {
  logr::log_print("...No group-level metadata path specified in config. Skipping.")
}


#=============================================================================#
# PART 2B: (NEW) CREATE GENOTYPE-LEVEL SUMMARY (WITH N > 3 FILTERING)
#=============================================================================#
logr::log_print("\n--- Starting PART 2B: Creating Genotype-Level Summary ---", console = TRUE)

# Define ID columns that should NOT be averaged
known_id_cols <- c(grouping_var, "rep", "pcr", "original_clone_id", "pseudo_clone_id", "clone_rank")

# Find all genotype-level metadata columns (non-numeric, non-group_id)
genotype_level_vars <- master_data_by_group %>%
  dplyr::select(
    -any_of(known_id_cols), 
    -where(is.numeric)
  ) %>%
  colnames()

if(length(genotype_level_vars) == 0) {
  logr::log_print("...no genotype-level metadata found. Summarizing by 'primary_var' only.")
  genotype_level_vars <- c(primary_var)
} else {
  logr::log_print(paste("...summarizing data by genotype, grouped by:", paste(genotype_level_vars, collapse=", ")))
}

master_data_by_genotype <- master_data_by_group %>%
  group_by(across(all_of(genotype_level_vars))) %>%
  summarise(
    across(where(is.numeric), ~mean(.x, na.rm = TRUE)),
    n_groups = n(),
    .groups = "drop"
  )

logr::log_print("...'by_genotype' summary data frame created.")
logr::log_print("...Applying filter: Genotype-level plots will only be generated for genotypes with > 3 groups (clones).")


#=============================================================================#
# PART 2C: (NEW) LOAD AND JOIN GENOTYPE-LEVEL METADATA
#=============================================================================#
logr::log_print("\n--- Starting PART 2C: Loading Genotype-Level Metadata ---", console = TRUE)

if (!is.null(meta_cfg$genotype_level_path) && meta_cfg$genotype_level_path != 'null') {
  genotype_meta_path <- here::here(meta_cfg$genotype_level_path)
  if (file.exists(genotype_meta_path)) {
    logr::log_print(paste("...Found optional genotype-level metadata. Loading:", basename(genotype_meta_path)))
    external_metadata_genotype <- readxl::read_excel(genotype_meta_path)
    
    if (primary_var %in% colnames(external_metadata_genotype)) {
      
      # --- This is the key join ---
      master_data_by_genotype <- master_data_by_genotype %>%
        left_join(external_metadata_genotype, by = primary_var)
      
      logr::log_print("...Genotype-level metadata successfully joined to 'by_genotype' data.")
    } else {
      logr::log_print(paste("...WARNING: Genotype-level metadata skipped. File MUST contain column:", primary_var))
    }
  } else {
    logr::log_print(paste("...WARNING: Genotype-level metadata file not found at:", genotype_meta_path))
  }
} else {
  logr::log_print("...No genotype-level metadata path specified in config. Skipping.")
}


#=============================================================================#
# --- (NEW) STATISTICAL ANALYSIS N-CHECK ---
#=============================================================================#
# Define the minimum number of data points (rows) needed to run
# correlations, PCA, clustering, and regression.
min_n_for_correlation <- 3 # Set to 5 for moderately robust stats

# Check 1: By-Group (Clone/Treatment)
run_by_group_stats <- nrow(master_data_by_group) >= min_n_for_correlation

# Check 2: By-Genotype (uses the *filtered* table from Part 2B)
run_by_genotype_stats <- nrow(master_data_by_genotype) >= min_n_for_correlation

if(!run_by_group_stats) {
  logr::log_print(paste("\n--- SKIPPING PARTS 3, 4, 5, 6, 7, 8 (By-Group) ---"))
  logr::log_print(paste("Reason: Not enough data points (N =", nrow(master_data_by_group), ") for correlation analysis."))
  logr::log_print(paste("Minimum required is", min_n_for_correlation, "."))
}


#=============================================================================#
# PART 3: (SMART) BUILD AND PLOT CORRELATION MATRICES (BY GROUP) ####
#=============================================================================#

if (run_by_group_stats) {
  logr::log_print("\n--- Starting PART 3: Building 'By-Group' Correlation Matrices ---", console = TRUE)
  
  # --- 1. Load Pretty Labels from Config ---
  logr::log_print("...loading 'variable_pretty_labels' from config.yml.")
  
  pretty_labels <- config$variable_pretty_labels %||% list()
  if (!primary_var %in% names(pretty_labels)) pretty_labels[[primary_var]] <- stringr::str_to_title(primary_var)
  allowed_vars <- names(pretty_labels)
  # --- 2. Prepare data for plotting (FIX: Sync Labels & Drop X) ---
  vars_for_ggpairs <- master_data_by_group %>% 
    # Bridge the raw IDs to the Publication Labels so the color palette matches!
    dplyr::left_join(label_lookup %>% dplyr::distinct(), by = primary_var) %>%
    dplyr::mutate(!!sym(primary_var) := !!sym(plot_label_col)) %>%
    dplyr::select(-!!sym(plot_label_col)) %>%
    
    # MAGIC LINE: Keep only numeric columns that are explicitly named in your config labels (plus Genotype)
    dplyr::select(where(is.numeric) & any_of(allowed_vars), !!sym(primary_var)) %>%
    
    janitor::remove_empty("cols") %>%
    dplyr::mutate(!!sym(primary_var) := as.character(!!sym(primary_var))) %>%
    dplyr::select(!!sym(primary_var), everything())
  current_col_names <- colnames(vars_for_ggpairs)
  
  # --- 3. Create Wrapped Labels array specifically for GGally ---
  final_labels_ordered <- sapply(current_col_names, function(x) {
    if (x %in% names(pretty_labels)) pretty_labels[[x]] else x
  })
  final_labels_wrapped <- unname(sapply(final_labels_ordered, function(x) stringr::str_wrap(x, width = 15)))
  
  # --- PART 3A: Numeric-Only Correlation Heatmap (BY GROUP) ---
  logr::log_print("...generating numeric correlation heatmap (by-group).")
  
  numeric_vars_to_correlate <- vars_for_ggpairs %>% 
    dplyr::select(where(is.numeric))
  
  numeric_col_names <- colnames(numeric_vars_to_correlate)
  numeric_labels <- sapply(numeric_col_names, function(x) {
    if(x %in% names(pretty_labels)) stringr::str_wrap(pretty_labels[[x]], 15) else x
  })
  names(numeric_labels) <- numeric_col_names
  
  if (ncol(numeric_vars_to_correlate) > 1) {
    cor_matrix <- cor(numeric_vars_to_correlate, use = "pairwise.complete.obs")
    
    # --- THE FIX: Mask the upper triangle ---
    # This sets all values above the diagonal to NA so they don't plot
    cor_matrix[upper.tri(cor_matrix)] <- NA 
    
    # Add na.rm = TRUE to drop those blank NAs from the plot dataframe
    melted_cor_matrix <- reshape2::melt(cor_matrix, na.rm = TRUE) 
    # ----------------------------------------
    
    p_corr_heatmap <- ggplot(melted_cor_matrix, aes(x = Var1, y = Var2, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(
        low = "#0072B2", high = "#D55E00", mid = "white",
        midpoint = 0, limit = c(-1, 1), name = "Correlation"
      ) +
      geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
      theme_minimal(base_size = 12) +
      labs(title = "Overall Correlation Matrix of Numeric Variables (By-Group)", x = "", y = "") +
      scale_x_discrete(labels = numeric_labels) +
      scale_y_discrete(labels = numeric_labels) +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 10), 
        axis.text.y = element_text(size = 10), 
        panel.grid = element_blank(),
        panel.background = element_blank() # Ensures the empty half is pure white
      ) +
      coord_fixed()
    
    print(p_corr_heatmap)
    ggsave(file.path(matrix_plot_dir, "corr_heatmap_by_group.tiff"), p_corr_heatmap, width = 9, height = 8, dpi = 300, compression = "lzw")
  } else {
    logr::log_print("...Skipping correlation heatmap (not enough numeric variables).")
  }
  
  # --- CUSTOM GEOM FUNCTIONS FOR GGPAIRS ---
  custom_barDiag <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_bar(color = "black", alpha = 0.7) +
      scale_x_discrete(drop = FALSE) 
  }
  custom_box_jitter <- function(data, mapping, ...) {
    ggplot(data = data, mapping = mapping) +
      geom_boxplot(color = "black", outlier.shape = NA, alpha = 0.6) +
      geom_jitter(width = 0.1, alpha = 0.7, height = 0, shape = 19)
  }
  
  # --- PART 3B: Mixed-Data ggpairs Plot (BY GROUP, Colored by Genotype) ---
  logr::log_print("...generating ggpairs matrix (By-Group, Colored by Genotype).")
  
  if (ncol(vars_for_ggpairs) > 1) {
    p_corr_ggpairs_color <- GGally::ggpairs(
      vars_for_ggpairs,
      columnLabels = final_labels_wrapped, 
      upper = list(continuous = GGally::wrap("cor", size = 3, use = "pairwise.complete.obs"), # FIX: Safe correlations
                   combo = custom_box_jitter),
      lower = list(continuous = GGally::wrap("points", alpha = 0.3, na.rm = TRUE), # FIX: Drop NAs
                   combo = "facetdensity"),
      diag = list(continuous = GGally::wrap("densityDiag", color = "black", alpha = 0.7, na.rm = TRUE), 
                  discrete = custom_barDiag),
      mapping = ggplot2::aes(color = .data[[primary_var]], fill = .data[[primary_var]]) # FIX: Safe pronoun mapping
    ) +
      labs(title = "By-Group Correlation Matrix (Colored by Genotype)") +
      theme_publication(base_size = 8) +
      theme(strip.text.x = element_text(size = 8), strip.text.y = element_text(size = 8), axis.text = element_text(size = 6))
    
    p_corr_ggpairs_color <- apply_smart_palette(p_corr_ggpairs_color, primary_var)
    
    print(p_corr_ggpairs_color)
    ggsave(file.path(matrix_plot_dir, "corr_matrix_by_genotype.tiff"), p_corr_ggpairs_color, width = 14, height = 14, dpi = 300, compression = "lzw")
  } else {
    logr::log_print("...Skipping by-group ggpairs plot (not enough variables).")
  }
  
  # --- PART 3C: Global ggpairs Plot (BY GROUP, All Data Combined) ---
  logr::log_print("...generating global ggpairs matrix (by-group, all data combined).")
  
  if (ncol(vars_for_ggpairs) > 1) {
    p_corr_ggpairs_global <- GGally::ggpairs(
      vars_for_ggpairs,
      columnLabels = final_labels_wrapped,
      upper = list(continuous = GGally::wrap("cor", size = 3, use = "pairwise.complete.obs"), combo = "box"),
      lower = list(continuous = GGally::wrap("points", alpha = 0.3, na.rm = TRUE), combo = "facetdensity"),
      diag = list(continuous = GGally::wrap("densityDiag", color = "black", alpha = 0.7, na.rm = TRUE), discrete = custom_barDiag)
    ) +
      labs(title = "Overall By-Group Correlation Matrix (All Genotypes Combined)") +
      theme_publication(base_size = 8) +
      theme(strip.text.x = element_text(size = 8), strip.text.y = element_text(size = 8), axis.text = element_text(size = 6))

    print(p_corr_ggpairs_global)
    ggsave(file.path(matrix_plot_dir, "corr_matrix_by_group.tiff"), p_corr_ggpairs_global, width = 14, height = 14, dpi = 300, compression = "lzw")
  } else {
    logr::log_print("...Skipping global ggpairs plot (not enough variables).")
  }
  
  # --- PART 3D: Faceted Correlation Heatmap (BY GROUP, by Genotype) ---
  logr::log_print("...generating faceted correlation heatmaps (by-group, per genotype).")
  
  faceted_cor_data <- NULL
  if (ncol(numeric_vars_to_correlate) > 1) {
    group_counts <- vars_for_ggpairs %>% count(!!sym(primary_var), name = "n_groups_per_genotype")
    if(any(group_counts$n_groups_per_genotype >= min_n_for_correlation)) {
      
      faceted_cor_data <- vars_for_ggpairs %>%
        left_join(group_counts, by = primary_var) %>%
        dplyr::filter(n_groups_per_genotype >= min_n_for_correlation) %>% 
        dplyr::select(-n_groups_per_genotype) %>%
        dplyr::group_by(!!sym(primary_var)) %>%
        dplyr::summarise(cor_matrix = list(cor(across(where(is.numeric)), use = "pairwise.complete.obs")), .groups = "drop") %>%
        rowwise() %>%
        mutate(melted_cor = list(reshape2::melt(cor_matrix))) %>%
        unnest(melted_cor) %>%
        dplyr::select(-cor_matrix) %>%
        left_join(group_counts, by = primary_var) %>% 
        mutate(
          facet_label = paste0(!!sym(primary_var), " (n=", n_groups_per_genotype, ")"),
          Var1 = factor(Var1, levels = numeric_col_names, labels = numeric_labels),
          Var2 = factor(Var2, levels = numeric_col_names, labels = numeric_labels)
        )
      
      p_corr_heatmap_faceted <- ggplot(faceted_cor_data, aes(x = Var1, y = Var2, fill = value)) +
        geom_tile(color = "white") +
        scale_fill_gradient2(low = "#0072B2", high = "#D55E00", mid = "white", midpoint = 0, limit = c(-1, 1), name = "Correlation") +
        geom_text(aes(label = round(value, 2)), color = "black", size = 2) + 
        facet_wrap(vars(facet_label)) +
        theme_minimal(base_size = 10) +
        labs(
          title = "By-Group Correlation Matrix (Faceted by Genotype)",
          subtitle = paste("Each panel shows correlations for", grouping_var, "within that genotype. (Only shows genotypes with N >=", min_n_for_correlation, ")"),
          x = "", y = ""
        ) +
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8), axis.text.y = element_text(size = 8), panel.grid = element_blank()) +
        coord_fixed()
      
      print(p_corr_heatmap_faceted)
      ggsave(file.path(matrix_plot_dir, "corr_heatmap_by_genotype.tiff"), p_corr_heatmap_faceted, width = 12, height = 10, dpi = 300, compression = "lzw")
    } else {
      logr::log_print(paste("...Skipping faceted heatmap (no single genotype has N >=", min_n_for_correlation, ")."))
    }
  } else {
    logr::log_print("...Skipping faceted heatmap (not enough numeric variables).")
  }
  
  # --- PART 3E: Faceted Bar Chart (BY GROUP, by Genotype) ---
  logr::log_print("...generating faceted bar chart (by-group, by Genotype).")
  
  if (!is.null(faceted_cor_data) && nrow(faceted_cor_data) > 0) {
    cor_barchart_data_A <- faceted_cor_data %>% dplyr::filter(as.character(Var1) != as.character(Var2))
    
    p_corr_barchart_by_geno <- ggplot(cor_barchart_data_A, aes(x = Var2, y = value, fill = !!sym(primary_var))) +
      geom_col(position = "dodge", color = "black") +
      facet_grid(rows = vars(Var1), cols = vars(facet_label), labeller = label_wrap_gen(width = 15)) +
      geom_hline(yintercept = 0, color = "grey20", linetype = "dashed") +
      scale_y_continuous(limits = c(-1, 1), breaks = c(-1, -0.5, 0, 0.5, 1)) +
      theme_publication(base_size = 10) +
      labs(
        title = "By-Group Correlation Coefficients by Genotype",
        subtitle = "Shows the Pearson's (R) correlation between all numeric variables.",
        x = "Variable 2", y = "Correlation Coefficient (R)", fill = "Genotype"
      ) +
      theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 8), legend.position = "none")
    
    p_corr_barchart_by_geno <- apply_smart_palette(p_corr_barchart_by_geno, primary_var)
    
    print(p_corr_barchart_by_geno)
    ggsave(file.path(matrix_plot_dir, "corr_bars_by_genotype.tiff"), p_corr_barchart_by_geno, width = 12, height = 12, dpi = 300, compression = "lzw")
  } else {
    logr::log_print("...Skipping faceted bar chart (no correlation data was generated in 3D).")
  }
  
  # --- PART 3F: Faceted Bar Chart (BY GROUP, by Correlation Pair) ---
  logr::log_print("...generating faceted bar chart of correlation values (by-group, by pair).")
  
  if (!is.null(faceted_cor_data) && nrow(faceted_cor_data) > 0) {
    cor_barchart_data_B <- faceted_cor_data %>% dplyr::filter(as.integer(Var1) < as.integer(Var2))
    
    p_corr_barchart_by_pair <- ggplot(cor_barchart_data_B, aes(x = !!sym(primary_var), y = value, fill = !!sym(primary_var))) +
      geom_col(position = "dodge", color = "black") +
      facet_grid(rows = vars(Var1), cols = vars(Var2), labeller = label_wrap_gen(width = 15)) +
      geom_hline(yintercept = 0, color = "grey20", linetype = "dashed") +
      scale_y_continuous(limits = c(-1, 1)) +
      theme_publication(base_size = 10) +
      labs(
        title = "By-Group Correlation Coefficients by Pair",
        subtitle = "Shows the Pearson's (R) correlation for each pair of variables, grouped by genotype.",
        x = "Genotype", y = "Correlation Coefficient (R)", fill = "Genotype"
      ) +
      theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1, size = 8), legend.position = "none")
    
    p_corr_barchart_by_pair <- apply_smart_palette(p_corr_barchart_by_pair, primary_var)
    
    print(p_corr_barchart_by_pair)
    ggsave(file.path(matrix_plot_dir, "corr_bars_by_pair.tiff"), p_corr_barchart_by_pair, width = 12, height = 12, dpi = 300, compression = "lzw")
  } else {
    logr::log_print("...Skipping faceted pair bar chart (no correlation data was generated in 3D).")
  }  
} # End of 'run_by_group_stats' block

#=============================================================================#
# PART 3G: (SMART) BUILD AND PLOT CORRELATION MATRIX (BY GENOTYPE) ####
#=============================================================================#
if (run_by_genotype_stats) {
  logr::log_print("\n--- Starting PART 3G: Building 'By-Genotype' Correlation Matrix ---", console = TRUE)
  
  master_data_by_genotype_filtered <- master_data_by_genotype %>%
    dplyr::filter(!is.na(!!sym(primary_var)) & !!sym(primary_var) != "NA")
  
  logr::log_print(paste("...Removed NA genotypes. Plotting with", nrow(master_data_by_genotype_filtered), "rows."))
  
  # FIX: Apply label mapping and drop X exactly as above
  vars_for_ggpairs_genotype <- master_data_by_genotype_filtered %>%
    dplyr::left_join(label_lookup %>% dplyr::distinct(), by = primary_var) %>%
    dplyr::mutate(!!sym(primary_var) := !!sym(plot_label_col)) %>%
    dplyr::select(-!!sym(plot_label_col)) %>%
    dplyr::select(where(is.numeric), !!sym(primary_var), -any_of(c("n_groups", "X", "X1", "x"))) %>%
    janitor::remove_empty("cols") %>%
    mutate(!!sym(primary_var) := as.character(!!sym(primary_var))) %>%
    dplyr::select(!!sym(primary_var), everything()) 
  
  current_col_names_geno <- colnames(vars_for_ggpairs_genotype)
  
  final_labels_ordered_geno <- sapply(current_col_names_geno, function(x) {
    if (x %in% names(pretty_labels)) pretty_labels[[x]] else x
  })
  final_labels_wrapped_geno <- unname(sapply(final_labels_ordered_geno, function(x) stringr::str_wrap(x, width = 15)))
  
  if (ncol(vars_for_ggpairs_genotype) > 1) {
    
    custom_cor <- function(data, mapping, ...) {
      x_var <- rlang::as_name(mapping$x)
      y_var <- rlang::as_name(mapping$y)
      cor_val <- cor(data[[x_var]], data[[y_var]], use = "pairwise.complete.obs")
      cor_label <- sprintf("Overall R = %.3f", cor_val)
      GGally::ggally_text(label = cor_label, color = "black", size = 4, ...) + theme(panel.background = element_rect(fill = "white"))
    }
    
    custom_scatter_with_trend <- function(data, mapping, ...) {
      ggplot(data = data, mapping = mapping) +
        geom_point(alpha = 0.7, size = 3) +
        geom_smooth(mapping = aes(group = 1), method = "lm", se = FALSE, color = "black", linewidth = 1, linetype = "dashed")
    }
    
    p_corr_ggpairs_genotype <- GGally::ggpairs(
      vars_for_ggpairs_genotype,
      columnLabels = final_labels_wrapped_geno, 
      upper = list(continuous = custom_cor, combo = custom_box_jitter),
      lower = list(continuous = custom_scatter_with_trend, combo = "facetdensity"),
      diag = list(continuous = GGally::wrap("densityDiag", alpha = 0.7), discrete = custom_barDiag),
      mapping = ggplot2::aes(color = .data[[primary_var]], fill = .data[[primary_var]])
    ) +
      labs(title = "Genotype-Level Correlation Matrix (All Data Averaged)") +
      theme_publication(base_size = 8) +
      theme(strip.text.x = element_text(size = 8), strip.text.y = element_text(size = 8), axis.text = element_text(size = 6), legend.position = "right")
    
   p_corr_ggpairs_genotype <- apply_smart_palette(p_corr_ggpairs_genotype, primary_var)
    
    print(p_corr_ggpairs_genotype)
    ggsave(file.path(matrix_plot_dir, "corr_matrix_by_geno.tiff"), p_corr_ggpairs_genotype, width = 14, height = 14, dpi = 300, compression = "lzw")
  } else {
    logr::log_print("...Skipping by-genotype ggpairs plot (not enough variables).")
  }
} else {
  logr::log_print("\n--- Skipping PART 3G (By-Genotype) ---")
  logr::log_print(paste("Reason: Not enough data points (N =", nrow(master_data_by_genotype), ") for correlation analysis."))
}

#=============================================================================#
# PART 4 & 5: (SMART) MULTI-VARIABLE STEPWISE REGRESSION ####
#=============================================================================#

if (run_by_group_stats) {
  logr::log_print("\n--- Starting PART 4 & 5: Stepwise Regression (Clone-Level Data) ---", console = TRUE)
  
  # 1. Identify all outcome variables (anything starting with "slope_")
  outcome_vars <- grep("^slope_", colnames(master_data_by_group), value = TRUE)
  
  # 2. Identify all potential predictors 
  # (All numeric columns EXCEPT IDs and the slope outcome columns)
  all_numeric_cols <- master_data_by_group %>%
    dplyr::select(where(is.numeric), -any_of(known_id_cols), -starts_with("slope_")) %>%
    colnames()
  
  if (length(all_numeric_cols) < 2) {
    logr::log_print("...Skipping Stepwise Regression (Not enough numeric predictor variables).")
    run_stepwise <- FALSE
  } else {
    logr::log_print(paste("...Pool of predictors:", paste(all_numeric_cols, collapse=", ")))
    run_stepwise <- TRUE
    
    # Initialize Excel workbook to store all models
    excel_file_path <- file.path(latest_analysis_dir, "Predictive_Model_Outputs.xlsx")
    wb_models <- openxlsx::createWorkbook()
    
    # List to store winning models for plotting in Part 6
    winning_models_list <- list()
  }
} else {
  run_stepwise <- FALSE
}

if (run_stepwise) {
  for (outcome_var in outcome_vars) {
    logr::log_print(paste("\n...Modeling Outcome:", outcome_var))
    
    # Drop rows with NA in the outcome variable to ensure clean modeling
    model_df <- master_data_by_group %>% 
      dplyr::filter(!is.na(!!sym(outcome_var))) %>%
      dplyr::select(all_of(outcome_var), all_of(all_numeric_cols)) %>%
      na.omit() # stepAIC requires a complete dataset without NAs
    
    if (nrow(model_df) < (length(all_numeric_cols) + 2)) {
      logr::log_print(paste("......Skipped (Not enough complete rows to model", length(all_numeric_cols), "predictors)."))
      next
    }
    
    # --- Build the Full and Null Models ---
    full_formula <- as.formula(paste(outcome_var, "~", paste(all_numeric_cols, collapse = " + ")))
    null_formula <- as.formula(paste(outcome_var, "~ 1"))
    
    full_model <- stats::lm(full_formula, data = model_df)
    null_model <- stats::lm(null_formula, data = model_df)
    
    # --- Run Stepwise Regression (Both Directions) ---
    logr::log_print("......running bidirectional stepAIC...")
    step_model <- tryCatch({
      MASS::stepAIC(full_model, 
                    scope = list(lower = null_model, upper = full_model), 
                    direction = "both", 
                    trace = 0) # trace=0 keeps the console clean
    }, error = function(e) { NULL })
    
    if (!is.null(step_model)) {
      # Save to list for Part 6
      winning_models_list[[outcome_var]] <- step_model
      
      # Tidy the output for Excel
      step_model_tidy <- broom::tidy(step_model, conf.int = TRUE) %>%
        dplyr::mutate(
          term_pretty = sapply(term, function(x) pretty_labels[[x]] %||% x),
          outcome = pretty_labels[[outcome_var]] %||% outcome_var
        ) %>%
        dplyr::select(outcome, term_pretty, estimate, std.error, p.value, conf.low, conf.high, term)
      
      # Add to Excel Workbook (Sheet name max 31 chars)
      sheet_name <- substr(paste0("Model_", gsub("slope_", "", outcome_var)), 1, 31)
      openxlsx::addWorksheet(wb_models, sheet_name)
      openxlsx::writeData(wb_models, sheet_name, step_model_tidy)
      
      logr::log_print(paste("......Optimal model found with", length(coef(step_model))-1, "predictors."))
    }
  }
  
  # Save the combined Excel file
  if (length(names(wb_models)) > 0) {
    openxlsx::saveWorkbook(wb_models, excel_file_path, overwrite = TRUE)
    logr::log_print(paste("\n...All stepwise models saved to:", basename(excel_file_path)))
  }
}

#=============================================================================#
# PART 6: (SMART) VISUALIZE STEPWISE MODEL RESULTS ####
#=============================================================================#

if (run_stepwise && length(winning_models_list) > 0) {
  logr::log_print("\n--- Starting PART 6: Visualizing Stepwise Models ---", console = TRUE)
  
  for (outcome_var in names(winning_models_list)) {
    step_model <- winning_models_list[[outcome_var]]
    pretty_outcome <- pretty_labels[[outcome_var]] %||% outcome_var
    safe_outcome_name <- gsub("slope_", "", outcome_var)
    
    # --- Plot 1: Coefficient Forest Plot ---
    step_model_tidy <- broom::tidy(step_model, conf.int = TRUE) %>%
      dplyr::filter(term != "(Intercept)")
    
    if (nrow(step_model_tidy) > 0) {
      logr::log_print(paste("...generating coefficient plot for", safe_outcome_name))
      
      step_model_plot_data <- step_model_tidy %>%
        mutate(
          label = sapply(term, function(x) stringr::str_wrap(pretty_labels[[x]] %||% x, 25)),
          label = reorder(label, estimate),
          significance = ifelse(p.value < 0.05, "Significant (p < 0.05)", "Not Significant")
        )
      
      p_coef_plot <- ggplot(step_model_plot_data, aes(x = estimate, y = label, color = significance)) +
        geom_vline(xintercept = 0, linetype = "dashed", color = "grey50", linewidth = 1) +
        geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width = 0.2, linewidth = 1) +
        geom_point(size = 4) +
        scale_color_manual(values = c("Significant (p < 0.05)" = "#D55E00", "Not Significant" = "grey60")) +
        labs(
          title = paste("Key Drivers of", pretty_outcome),
          subtitle = "Displays the estimated effect size (with 95% CI) of variables selected by the stepwise model.",
          x = "Estimated Effect on Rate",
          y = "Predictor Variable",
          color = "P-Value"
        ) +
        theme_publication() + theme(legend.position = "bottom")
      
      print(p_coef_plot)
      ggsave(file.path(corr_plot_dir, paste0("predictive_drivers_", safe_outcome_name, ".tiff")), 
             p_coef_plot, width = 10, height = 7, dpi = 300, compression = "lzw")
      
      # --- Plot 2: Partial Residual Plots for Significant Predictors ---
      # Find predictors that are actually significant
      sig_predictors <- step_model_tidy %>% dplyr::filter(p.value < 0.05) %>% pull(term)
      
      if (length(sig_predictors) > 0) {
        logr::log_print(paste("......generating partial effect plots for", length(sig_predictors), "significant predictors."))
        
        for (pred_var in sig_predictors) {
          pretty_pred <- pretty_labels[[pred_var]] %||% pred_var
          
          # Use ggeffects to calculate marginal effects
          pred_data <- ggeffects::ggpredict(step_model, terms = pred_var)
          
          p_effect <- ggplot(pred_data, aes(x = x, y = predicted)) +
            geom_line(color = "#0072B2", linewidth = 1.2) +
            geom_ribbon(aes(ymin = conf.low, ymax = conf.high), fill = "#0072B2", alpha = 0.2) +
            labs(
              title = paste("Isolated Effect:", pretty_pred, "on", pretty_outcome),
              subtitle = "Shows the predicted rate when controlling for all other model variables.",
              x = pretty_pred,
              y = paste("Predicted", pretty_outcome)
            ) +
            theme_publication()
          
          ggsave(file.path(corr_plot_dir, paste0("effect_", safe_outcome_name, "_vs_", pred_var, ".tiff")), 
                 p_effect, width = 7, height = 6, dpi = 300, compression = "lzw")
        }
      }
    } else {
      logr::log_print(paste("...Stepwise model for", safe_outcome_name, "found NO significant predictors (Null model won)."))
    }
  }
}

#=============================================================================#
# PART 7: (SMART) PRINCIPAL COMPONENT ANALYSIS (PCA)
#=============================================================================#

if (run_by_group_stats) {
  logr::log_print("\n--- Starting PART 7: Principal Component Analysis (PCA) ---", console = TRUE)
  
  # 1. Prepare data (use the 'by_group' data)
  pca_data_numeric <- master_data_by_group %>%
    dplyr::select(where(is.numeric), -any_of(known_id_cols)) %>%
    na.omit() # PCA cannot handle missing values
  
  # 2. Run PCA
  if(nrow(pca_data_numeric) > 1 && ncol(pca_data_numeric) > 1) {
    pca_fit <- prcomp(pca_data_numeric, scale. = TRUE, center = TRUE)
    
    # 3. Get metadata for plotting
    pca_data_with_meta <- pca_fit$x %>%
      as.data.frame() %>%
      dplyr::bind_cols(master_data_by_group[rownames(pca_data_numeric), ])
    
    # 4. Plot Scree Plot (Variance Explained)
    pca_var_explained <- broom::tidy(pca_fit, matrix = "pcs") %>%
      dplyr::filter(PC <= 10)
    
    p_scree <- ggplot(pca_var_explained, aes(x = PC, y = percent)) +
      geom_col(fill = "#0072B2", alpha = 0.8) +
      geom_line(aes(y = cumulative), color = "red", group = 1) +
      geom_point(aes(y = cumulative), color = "red") +
      scale_y_continuous(labels = scales::percent_format(accuracy = 1)) +
      labs(
        title = "PCA: Variance Explained by Component",
        subtitle = "Red line shows cumulative variance.",
        x = "Principal Component (PC)",
        y = "Variance Explained"
      ) +
      theme_publication()
    
    print(p_scree)
    ggsave(file.path(pca_plot_dir, "pca_scree_plot.tiff"), p_scree, 
           width = 8, height = 6, dpi = 300, compression = "lzw")
    
    # 5. Plot PCA Biplot (PC1 vs PC2)
    pc1_var <- scales::percent(pca_var_explained$percent[1], accuracy = 1)
    pc2_var <- scales::percent(pca_var_explained$percent[2], accuracy = 1)
    
    pca_loadings <- as.data.frame(pca_fit$rotation) %>%
      tibble::rownames_to_column("variable") %>%
      dplyr::select(variable, PC1, PC2)
    
    arrow_labels <- unlist(pretty_labels[names(pretty_labels) %in% pca_loadings$variable])
    
    pca_loadings <- pca_loadings %>%
      dplyr::mutate(label = dplyr::coalesce(arrow_labels[variable], variable))
    
    arrow_scale <- max(abs(pca_data_with_meta$PC1), abs(pca_data_with_meta$PC2)) / max(abs(pca_loadings$PC1), abs(pca_loadings$PC2)) * 0.8
    
    p_biplot <- ggplot(pca_data_with_meta, aes(x = PC1, y = PC2)) +
      geom_point(aes(color = !!sym(primary_var)), size = 3, alpha = 0.7) +
      ggrepel::geom_text_repel(aes(label = !!sym(grouping_var)), size = 2) +
      geom_segment(
        data = pca_loadings, 
        aes(x = 0, y = 0, xend = PC1 * arrow_scale, yend = PC2 * arrow_scale),
        arrow = arrow(length = unit(0.2, "cm")), color = "black", alpha = 0.8
      ) +
      ggrepel::geom_text_repel(
        data = pca_loadings,
        aes(x = PC1 * arrow_scale, y = PC2 * arrow_scale, label = label),
        color = "black", fontface = "bold",
        box.padding = 0.5
      ) +
      labs(
        title = "PCA Biplot (PC1 vs. PC2)",
        subtitle = paste("Points are groups (", grouping_var, "). Arrows are variables."),
        x = paste0("PC1 (", pc1_var, ")"),
        y = paste0("PC2 (", pc2_var, ")"),
        color = stringr::str_to_title(primary_var)
      ) +
      theme_publication() +
      geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey50")
    

    p_biplot <- apply_smart_palette(p_biplot, primary_var)
    
    print(p_biplot)
    ggsave(file.path(pca_plot_dir, "pca_biplot.tiff"), p_biplot, 
           width = 11, height = 9, dpi = 300, compression = "lzw")
    
    # 6. (NEW) Plot PCA Loadings Heatmap
    logr::log_print("...generating PCA loadings heatmap.")
    
    # Get the rotation matrix (loadings)
    pca_loadings_matrix <- pca_fit$rotation
    
    # Use reshape2::melt to convert matrix to long format for ggplot
    melted_pca_loadings <- reshape2::melt(pca_loadings_matrix, 
                                          varnames = c("variable", "PC"))
    
    # Get pretty labels
    heatmap_labels <- unlist(pretty_labels[names(pretty_labels) %in% melted_pca_loadings$variable])
    
    p_loadings_heatmap <- ggplot(melted_pca_loadings, 
                                 aes(x = PC, y = variable, fill = value)) +
      geom_tile(color = "white") +
      scale_fill_gradient2(
        low = "#0072B2", high = "#D55E00", mid = "white",
        midpoint = 0, limit = c(-1, 1), name = "Loading"
      ) +
      geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
      scale_y_discrete(labels = heatmap_labels) +
      labs(
        title = "PCA Component Loadings",
        subtitle = "Shows which variables contribute to each component.",
        x = "Principal Component",
        y = "Variable"
      ) +
      theme_minimal(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1)
      )
    
    # Per your preference: Print plot to R viewer
    print(p_loadings_heatmap)
    
    # Per your preference: Save as TIFF to nested folder
    ggsave(
      file.path(pca_plot_dir, "pca_loadings_heatmap.tiff"),
      p_loadings_heatmap, width = 8, height = 8, dpi = 300, compression = "lzw"
    )
    
    # 7. (NEW) Plot PCA Loadings Bar Charts
    logr::log_print("...generating PCA loadings bar charts.")
    
    # 1. Ensure pretty_labels is a standard named vector (not a list)
    pretty_labels_vec <- unlist(pretty_labels)
    
    # 2. Run the PCA tidy + label step
    pca_loadings_tidy <- broom::tidy(pca_fit, matrix = "rotation") %>%
      dplyr::mutate(
        # Look up the name in the vector; if NA (no match), keep the original 'column' ID
        label = dplyr::coalesce(pretty_labels_vec[column], column)
      )
    
    p_loadings_bars <- ggplot(pca_loadings_tidy, 
                              aes(x = reorder(label, value), y = value, fill = PC)) +
      geom_col(color = "black", alpha = 0.8) +
      coord_flip() + # Flip to horizontal bar chart
      facet_wrap(vars(PC), scales = "free_x") +
      labs(
        title = "PCA Component Loadings",
        subtitle = "Shows the contribution of each variable to each component.",
        x = "Variable",
        y = "Loading (Contribution Strength)"
      ) +
      theme_publication() +
      theme(
        legend.position = "none",
        strip.text = element_text(face = "bold")
      )
    
    # Per your preference: Print plot to R viewer
    print(p_loadings_bars)
    
    # Per your preference: Save as TIFF to nested folder
    ggsave(
      file.path(pca_plot_dir, "pca_loadings_bar_charts.tiff"),
      p_loadings_bars, width = 12, height = 8, dpi = 300, compression = "lzw"
    )
    
  } else {
    logr::log_print("...Skipping PCA (not enough data after NA removal or only 1 variable).")
  }
} # End of 'run_by_group_stats' block

#=============================================================================#
# PART 8: (SMART) HIERARCHICAL CLUSTERING
#=============================================================================#

if (run_by_group_stats) {
  # --- Re-open log file in case it timed out ---
  try(logr::log_close(), silent = TRUE)
  lf <- logr::log_open(log_path, show_notes = FALSE, logdir = FALSE)
  
  logr::log_print("\n--- Starting PART 8: Hierarchical Clustering ---", console = TRUE)
  
  # --- (FIX) This is the new, robust data preparation step ---
  
  # 1. Prepare data (use the 'by_group' data again)
  
  # Select ONLY the numeric columns for clustering,
  # EXCLUDING any numeric ID columns.
  cluster_data_numeric_with_na <- master_data_by_group %>%
    dplyr::select(where(is.numeric), -any_of(known_id_cols))
  
  # Select the metadata columns we need for labels
  cluster_data_meta_with_na <- master_data_by_group %>%
    dplyr::select(!!sym(primary_var), !!sym(grouping_var))
  
  # Create a complete table *first*, then remove rows with NAs
  # This ensures numeric data and metadata stay perfectly aligned.
  complete_data_for_cluster <- cluster_data_numeric_with_na %>%
    dplyr::bind_cols(cluster_data_meta_with_na) %>%
    na.omit() # This removes NAs from numeric columns
  
  # Now, split back into numeric and meta
  cluster_data_numeric <- complete_data_for_cluster %>%
    dplyr::select(where(is.numeric))
  
  cluster_data_meta <- complete_data_for_cluster %>%
    dplyr::select(!!sym(primary_var), !!sym(grouping_var))
  
  # --- END FIX ---
  
  # 2. Scale data and calculate distance
  if(nrow(cluster_data_numeric) > 2) {
    dist_matrix <- dist(scale(cluster_data_numeric), method = "euclidean")
    
    # 3. Perform clustering
    hclust_fit <- hclust(dist_matrix, method = "ward.D2")
    
    # 4. Create color mapping for labels
    logr::log_print("...creating colored dendrogram.")
    
    # Safely define the active palette
    active_palette <- if(exists("primary_palette")) primary_palette else if(exists("x_var_palette")) x_var_palette else NULL
    
    # Re-order our metadata to match the plot order (hclust_fit$order)
    genotypes_in_order <- cluster_data_meta[[primary_var]][hclust_fit$order]
    
    # Map genotypes to their colors from the palette (with fallback to black if palette is missing)
    if (!is.null(active_palette)) {
      colors_for_labels <- unname(active_palette[as.character(genotypes_in_order)])
    } else {
      colors_for_labels <- rep("black", length(genotypes_in_order))
    }
    
    # Get the labels (clone IDs) in the correct plot order
    labels_in_order <- cluster_data_meta[[grouping_var]][hclust_fit$order]
    
    # 5. Plot with dendextend
    dend <- as.dendrogram(hclust_fit)
    
    # Assign the correctly-ordered labels AND colors
    dend <- dend %>%
      dendextend::set("labels", labels_in_order) %>%
      dendextend::set("labels_col", colors_for_labels) %>% 
      dendextend::hang.dendrogram(hang_height = .6)
    
    # Save the plot
    tiff(file.path(corr_plot_dir, "cluster_dendrogram.tiff"), 
         width = 12, height = 8, units = "in", res = 300, compression = "lzw")
    
    plot(dend, 
         main = "Hierarchical Clustering of Groups (Colored by Genotype)",
         sub = "Based on all scaled numeric variables",
         xlab = paste("Groups (", grouping_var, ")"),
         ylab = "Height (Ward's distance)")
    
    # Add a legend safely
    if (!is.null(active_palette)) {
      legend("topright", 
             title = stringr::str_to_title(primary_var), 
             legend = names(active_palette), 
             fill = active_palette,
             bty = "n")
    }
    
    dev.off()
    
    # --- Also print to RStudio viewer ---
    plot(dend, 
         main = "Hierarchical Clustering of Groups (Colored by Genotype)",
         sub = "Based on all scaled numeric variables",
         xlab = paste("Groups (", grouping_var, ")"),
         ylab = "Height (Ward's distance)")
    
    # Add a legend safely
    if (!is.null(active_palette)) {
      legend("topright", 
             title = stringr::str_to_title(primary_var), 
             legend = names(active_palette), 
             fill = active_palette,
             bty = "n")
    } else {
    logr::log_print("...Skipping hierarchical clustering (not enough numeric variables or rows).")
  }
} # End of 'run_by_group_stats' block
}
# --- Finalize ---
logr::log_print("\n--- SCRIPT 03: CORRELATION & INTEGRATION FINISHED SUCCESSFULLY ---", console = TRUE)
logr::log_close()

