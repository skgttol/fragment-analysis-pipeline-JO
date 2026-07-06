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
func_path <- "C:/Users/skgttol/OneDrive - University College London/PhD/PhD_Thesis/02_Data/Template_RScripts/CAGsizing_Flexible/functions.R"
if (!file.exists(func_path)) {
  if(file.exists("functions.R")) func_path <- "functions.R" else stop("functions.R not found.")
}
source(func_path)

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

# --- Finalize ---
logr::log_print("\n--- SCRIPT 03: CORRELATION & INTEGRATION FINISHED SUCCESSFULLY ---", console = TRUE)
logr::log_close()

