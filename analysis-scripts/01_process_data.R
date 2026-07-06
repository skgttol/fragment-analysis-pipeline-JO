#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                             #
#            Advanced Fragment Analysis Workflow - Script 1 of 5              #
#                        (Data Processing & QC)                               #
#                                                                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# DESCRIPTION:
# This script loads the PROCESSED data (Excel output) from the external
# fragment analysis program. It performs:
#  1. Hybrid Data Loading (patches names using Platemap if needed).
#  2. Data Cleaning & Factorizing.
#  3. Normalization (Calculation of Baseline & Mode Change).
#  4. Data Summarization (Replicates, Groups).
#  5. Outlier Detection (Biological & Technical).
#  6. Export of clean data for Script 02 (Plotting).
#
# PREREQUISITES:
# - The external analysis (e.g., Romeo) must be finished.
# - The output Excel file must be saved in the path defined in 'config.yml'.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#=============================================================================#
# PART 0: SETUP
#=============================================================================#

# 0A. Load Libraries ----
packages <- c("tidyverse", "readxl", "here", "openxlsx", "janitor", "yaml", "logr", "ggrepel", "ggpubr")
invisible(lapply(packages, library, character.only = TRUE))

# 0B. Load Functions ----
#-------------------------#
if (!file.exists(here::here("functions.R"))) {
  stop("CRITICAL ERROR: 'functions.R' not found. Please ensure it is in the main project directory.")
}
source(here::here("functions.R"))

# 0C. Load Config & Logging ----
if (!file.exists(here::here("config.yml"))) {
  stop("CRITICAL ERROR: 'config.yml' not found.")
}
config <- yaml::read_yaml(here::here("config.yml"))

# Create output directory
output_dir <- here::here(
  config$paths$output_dir_base,
  paste0(format(Sys.Date(), "%Y-%m-%d"), "_Analysis_v", format(Sys.time(), "%H%M"))
)
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)

# Open Log
try(logr::log_close(), silent = TRUE)
log_path <- file.path(output_dir, "analysis_log.log")
lf <- logr::log_open(log_path, show_notes = FALSE, logdir = FALSE)

logr::log_print("--- SCRIPT 01: DATA PROCESSING INITIALIZED ---", console = TRUE)
logr::log_print(paste("Output will be saved to:", output_dir))


#=============================================================================#
# PART 1: VALIDATION & DATA LOADING ####
#=============================================================================#

## --- 1A. Validate Config & Check for External File ---####
# We set check_external_files = TRUE because this script REQUIRES the data.
validate_config(config, check_external_files = TRUE)

## --- 1B. Load Platemap (Optional) ---####
# We try to load the platemap to assist with patching missing sample names.
# If FSA files aren't found locally, we skip the metadata link but still load the platemap file.
platemap_data <- NULL

if (dir.exists(here::here(config$paths$fsa_folder)) &&
    length(list.files(here::here(config$paths$fsa_folder), pattern = "\\.fsa$", recursive = TRUE)) > 0) {
  
  logr::log_print("FSA folder found. Loading metadata for robust linking...")
  platemap_raw <- load_platemap(config)
  fsa_metadata <- get_fsa_metadata(config$paths$fsa_folder, config)
  
  # Join
  fsa_df_for_join <- fsa_metadata %>% dplyr::select(plate = plate_from_filename, well = well_from_filename, fsa_filename)
  platemap_raw_lower <- platemap_raw %>% dplyr::mutate(plate = tolower(plate), well = tolower(well))
  fsa_df_for_join_lower <- fsa_df_for_join %>% dplyr::mutate(plate = tolower(plate), well = tolower(well))
  platemap_data <- dplyr::left_join(platemap_raw_lower, fsa_df_for_join_lower, by = c("plate", "well"))
  
} else {
  logr::log_print("NOTE: No local FSA files found. Attempting to load Platemap file directly for patching...")
  # Fallback: Just load the platemap file so we can map well -> sample_name
  tryCatch({
    platemap_data <- load_platemap(config)
    # If no FSA files, we can't fully join, but we pass the raw platemap to the loader
    # The loader will handle matching based on available columns.
  }, error = function(e) {
    logr::log_print("WARNING: Could not load platemap. Data loading will rely entirely on sample names in the Excel file.")
  })
}

## --- 1C. Load & Process External Data ---####
logr::log_print("\n--- Loading and Normalizing Sizing Data ---", console = TRUE)

data_load_list <- load_and_process_data(config, platemap_data)
processed_data <- data_load_list$processed_data
excluded_data <- data_load_list$excluded_data
data_with_reasons <- data_load_list$pre_filter_data_with_reasons

## --- 1D. Apply Variable Factor Ordering ---####
all_grouping_vars <- get_grouping_vars(config)

processed_data <- apply_factor_levels(
  data = processed_data,
  config = config,
  all_grouping_vars = all_grouping_vars
)

#=============================================================================#
# PART 2: NORMALIZATION TO BASELINE ####
#=============================================================================#
logr::log_print("\n--- Normalizing Data to Baseline ---", console = TRUE)

normalization_list <- normalize_data(
  processed_data = processed_data,
  config = config,
  all_grouping_vars = all_grouping_vars
)
normalized_data <- normalization_list$normalized_data
baselines_table <- normalization_list$baselines_table
print(baselines_table)

#=============================================================================#
# PART 2B: BUILD DISPLAY LABELS & ENFORCE FACTOR ORDER ####
#=============================================================================#
logr::log_print("Building display label columns and enforcing factor order...")

label_result <- build_label_columns(
  data      = normalized_data,
  baselines = baselines_table,
  config    = config,
  extra_dfs = list()   # summary data frames not yet created; applied again in Part 3
)

normalized_data <- label_result$primary
baselines_table <- label_result$baselines
genotype_map    <- label_result$genotype_map

logr::log_print("Label columns created. Genotype map (in plot order):")
print(genotype_map)

#=============================================================================#
# PART 3: SUMMARIZATION & OUTLIER DETECTION ---- ####
#=============================================================================#
logr::log_print("\n--- Summarizing Data & Checking Outliers ---", console = TRUE)

# --- 3A. Summarize ---
summary_list <- summarize_data(
  normalized_data = normalized_data,
  config = config,
  all_grouping_vars = all_grouping_vars
)

data_per_pcr <- summary_list$data_per_pcr
summary_per_rep <- summary_list$summary_per_rep
summary_by_group <- summary_list$summary_by_group

# Propagate label columns and factor order to all summary data frames.
# genotype_map was built in Part 2B and is reused here so the order is
# guaranteed identical across every data frame in the pipeline.
label_result2 <- build_label_columns(
  data      = data_per_pcr,
  baselines = baselines_table,  # baselines_table already has labels; used only for the map
  config    = config,
  extra_dfs = list(
    summary_per_rep  = summary_per_rep,
    summary_by_group = summary_by_group
  )
)
data_per_pcr     <- label_result2$primary
summary_per_rep  <- label_result2$summary_per_rep
summary_by_group <- label_result2$summary_by_group
# baselines_table already correct from Part 2B; no need to reassign
logr::log_print("Label columns propagated to all summary data frames.")

# --- PART 3B. OUTLIER DETECTION ---####
diagnostics_dir <- file.path(output_dir, "01_QC_plots")
dir.create(diagnostics_dir, showWarnings = FALSE)

response_vars <- names(config$response_variables)
response_labels <- config$response_variable_labels %||% list()
response_shortnames <- config$response_variable_shortnames %||% list()

cfg_vars <- config$key_variables
rep_var <- cfg_vars$optional_grouping_var %||% cfg_vars$repeated_measure_var
if (is.null(rep_var) || rep_var == 'null') rep_var <- "rep"
time_var <- cfg_vars$time_variable
primary_var <- cfg_vars$primary_group_var
secondary_var <- cfg_vars$secondary_group_var

bio_outliers_list <- list()
tech_outliers_list <- list()

logr::log_print("Generating outlier visualization plots (IQR Method)...")

for (current_var in response_vars) {
  
  y_axis_label <- response_labels[[current_var]] %||% stringr::str_to_title(current_var)
  short_resp_var <- response_shortnames[[current_var]] %||% current_var
  
  ## --- Plot 1: Biological Replicate Outliers ---####
  logr::log_print(paste("...generating biological outlier plot for", current_var))
  
  # Calculate IQR and Identify Outliers (Grouped by Genotype, Clone, Time)
  # The variation here is between Biological Replicates (Wells)
  summary_with_outliers <- summary_per_rep %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(c(all_grouping_vars, time_var)))) %>%
    dplyr::mutate(
      Q1 = stats::quantile(!!sym(current_var), 0.25, na.rm = TRUE),
      Q3 = stats::quantile(!!sym(current_var), 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      is_outlier = !!sym(current_var) < (Q1 - 1.5 * IQR) | !!sym(current_var) > (Q3 + 1.5 * IQR)
    ) %>%
    dplyr::ungroup()
  
  outliers_found <- dplyr::filter(summary_with_outliers, is_outlier == TRUE)
  if (nrow(outliers_found) > 0) {
    logr::log_print(paste("Potential biological outliers for", current_var, "identified:"))
    # logr::log_print(outliers_found) # Uncomment to see in log
    bio_outliers_list[[current_var]] <- outliers_found
  }
  
  # Define facets dynamically
  facet_vars <- if (!is.null(secondary_var) && secondary_var != 'null') {
    vars(!!sym(primary_var), !!sym(secondary_var))
  } else {
    vars(!!sym(primary_var))
  }
  
  p_bio_diag <- summary_with_outliers %>%
    dplyr::mutate(!!time_var := as.numeric(as.character(!!sym(time_var)))) %>%
    ggplot(aes(x = !!sym(time_var), y = !!sym(current_var))) +
    geom_boxplot(
      aes(fill = !!sym(primary_var),
          group = interaction(!!sym(time_var), !!sym(primary_var))
      ),
      outlier.shape = NA,
      alpha = 0.4,
      width = 0.8
    ) +
    geom_point(
      aes(shape = is_outlier, color = is_outlier),
      size = 2,
      position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.1)
    ) +
    ggrepel::geom_text_repel(
      data = . %>% dplyr::filter(is_outlier == TRUE),
      aes(label = !!sym(rep_var)),
      box.padding = 0.5,
      max.overlaps = Inf
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    scale_shape_manual(values = c("TRUE" = 8, "FALSE" = 19)) +
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    facet_wrap(facet_vars, scales = "free_y", axes = "all") +
    labs(
      title = paste("Biological Outlier Detection:", y_axis_label),
      subtitle = "Outliers (1.5xIQR) are highlighted in red and labeled by biological replicate ID.",
      x = stringr::str_to_title(time_var),
      y = y_axis_label
    ) +
    theme_publication() +
    theme(legend.position = "none")
  
  print(p_bio_diag)
  ggsave(
    file.path(diagnostics_dir, paste0("p_outlier_bio_", short_resp_var, ".tiff")),
    p_bio_diag, width = 14, height = 10, dpi = 150, device = 'tiff', compression = "lzw"
  )
  
  # --- AUTO-REMOVE BIOLOGICAL OUTLIERS ---
  if (isTRUE(config$parameters$remove_bio_outliers)) {
    logr::log_print(paste("...Auto-removing biological outliers for", current_var, "(Dropping entire row)"))
    
    # Filter out the outlier rows completely, then drop the temporary calculation columns
    summary_per_rep <- summary_with_outliers %>%
      dplyr::filter(is_outlier == FALSE) %>%
      dplyr::select(-Q1, -Q3, -IQR, -is_outlier)
    
  } else {
    # Even if we don't remove them, we must clean up the temp columns for the next loop iteration
    summary_per_rep <- summary_with_outliers %>%
      dplyr::select(-Q1, -Q3, -IQR, -is_outlier)
  }
  
  # --- Plot 2: Technical Replicate (PCR) Outliers ---####
  logr::log_print(paste("...generating technical outlier plot for", current_var))
  
  # Group by Bio-Rep (Genotype + Clone + Time + Rep) to find outlier PCRs
  tech_grouping <- unique(c(all_grouping_vars, time_var, rep_var))
  tech_grouping <- tech_grouping[tech_grouping %in% colnames(data_per_pcr)]
  
  data_with_tech_outliers <- data_per_pcr %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(tech_grouping))) %>%
    dplyr::mutate(
      Q1 = stats::quantile(!!sym(current_var), 0.25, na.rm = TRUE),
      Q3 = stats::quantile(!!sym(current_var), 0.75, na.rm = TRUE),
      IQR = Q3 - Q1,
      # Only calc outlier if we have >2 points, otherwise everything is FALSE
      is_outlier = if(n() > 2) (!!sym(current_var) < (Q1 - 1.5 * IQR) | !!sym(current_var) > (Q3 + 1.5 * IQR)) else FALSE
    ) %>%
    dplyr::ungroup()
  
  tech_outliers_found <- dplyr::filter(data_with_tech_outliers, is_outlier == TRUE)
  if (nrow(tech_outliers_found) > 0) {
    logr::log_print(paste("Potential TECHNICAL (PCR) outliers for", current_var, "identified:"))
    tech_outliers_list[[current_var]] <- tech_outliers_found
  }
  
  label_col <- if ("pcr" %in% colnames(data_per_pcr)) "pcr" else rep_var
  
  # UPDATED PLOT LAYOUT
  p_tech_diag <- data_with_tech_outliers %>%
    dplyr::mutate(!!time_var := as.numeric(as.character(!!sym(time_var))), !!sym(rep_var) := as.factor(!!sym(rep_var))) %>%
    ggplot(aes(x = !!sym(time_var), y = !!sym(current_var))) +
    
    # Removed geom_boxplot entirely. Relying purely on points for n=3 visibility.
    geom_point(
      aes(fill = !!sym(rep_var), shape = is_outlier, color = is_outlier, group = !!sym(rep_var)),
      size = 2,
      position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)
    ) +
    ggrepel::geom_text_repel(
      data = function(df) dplyr::filter(df, is_outlier == TRUE),
      aes(group = !!sym(rep_var), label = !!sym(label_col)),
      box.padding = 0.5,
      max.overlaps = Inf,
      position = position_jitterdodge(dodge.width = 0.8, jitter.width = 0.2)
    ) +
    scale_x_continuous(breaks = scales::pretty_breaks(n = 10)) +
    # Use fillable shapes so the biological replicate fill shows through
    scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 21)) + 
    scale_color_manual(values = c("TRUE" = "red", "FALSE" = "black")) +
    
    # CRITICAL FIX: Removed axes = "all"
    facet_wrap(facet_vars, scales = "free_y", axes = "all") + 
    
    labs(
      title = paste("Technical (PCR) Outlier Detection:", y_axis_label),
      subtitle = "Points represent individual PCRs. Biological replicates distinguished by color.",
      x = stringr::str_to_title(time_var),
      y = y_axis_label,
      fill = "Bio. Rep",
      shape = "Outlier Detection"
    ) +
    theme_publication() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "right")
  
  ggsave(
    file.path(diagnostics_dir, paste0("p_outlier_tech_", short_resp_var, ".tiff")),
    p_tech_diag, width = 16, height = 10, dpi = 150, device = 'tiff', compression = "lzw"
  )  
  
  # --- AUTO-REMOVE TECHNICAL OUTLIERS ---
  if (isTRUE(config$parameters$remove_tech_outliers)) {
    logr::log_print(paste("...Auto-removing technical (PCR) outliers for", current_var, "(Dropping entire row)"))
    
    # Filter out the outlier rows completely, then drop the temporary calculation columns
    data_per_pcr <- data_with_tech_outliers %>%
      dplyr::filter(is_outlier == FALSE) %>%
      dplyr::select(-Q1, -Q3, -IQR, -is_outlier)
    
  } else {
    # Clean up temp columns for the next loop iteration
    data_per_pcr <- data_with_tech_outliers %>%
      dplyr::select(-Q1, -Q3, -IQR, -is_outlier)
  }
  
} # End of response variable loop

#=============================================================================#
# PART 4: PSEUDO-CLONES & PRE-EXCLUSION QC ---####
#=============================================================================#

## --- 4A. Pseudo-Clone Mapping ---####

# 1. Check Config Flag
# Default to FALSE if missing to be safe
use_pseudo_clones_flag <- config$parameters$use_pseudo_clones %||% FALSE

# 2. Identify the variable to target (Secondary or Repeated)
target_var <- cfg_vars$secondary_group_var %||% cfg_vars$repeated_measure_var

# 3. Determine if we should run the ranking logic
# Must have the flag TRUE AND a valid variable to rank
should_run_ranking <- (use_pseudo_clones_flag && !is.null(target_var) && target_var != 'null')

# Global boolean for Script 02 to see
has_pseudo_clones <- should_run_ranking
is_clone_style_analysis <- should_run_ranking # Alias

if (should_run_ranking) {
  log_print(paste("Generating ranked pseudo-IDs for variable:", target_var))
  
  # --- Execute Pseudo-Clone Ranking Logic ---
  pseudo_clone_map <- summary_per_rep %>%
    mutate(!!sym(target_var) := as.character(!!sym(target_var))) %>%
    distinct(!!sym(primary_var), !!sym(target_var)) %>%
    group_by(!!sym(primary_var)) %>%
    arrange(!!sym(target_var)) %>%
    mutate(
      clone_rank = row_number(),
      pseudo_clone_id = paste(!!sym(primary_var), clone_rank, sep = "-"),
      original_clone_id = !!sym(target_var)
    ) %>% ungroup() %>%
    mutate(clone_rank = factor(clone_rank), pseudo_clone_id = factor(pseudo_clone_id, levels = unique(pseudo_clone_id)))
  
  join_map <- dplyr::select(pseudo_clone_map, !!sym(primary_var), original_clone_id, pseudo_clone_id, clone_rank)
  join_key <- c(primary_var, setNames("original_clone_id", target_var)) 
  
  # Apply map to all dataframes
  data_per_pcr <- data_per_pcr %>% mutate(!!sym(target_var) := as.character(!!sym(target_var))) %>% left_join(join_map, by = join_key)
  summary_per_rep <- summary_per_rep %>% mutate(!!sym(target_var) := as.character(!!sym(target_var))) %>% left_join(join_map, by = join_key)
  summary_by_group <- summary_by_group %>% mutate(!!sym(target_var) := as.character(!!sym(target_var))) %>% left_join(join_map, by = join_key)
  
  # Set coloring to the new Rank ID
  color_var <- "pseudo_clone_id" 
  export_to_excel(list('Clone_ID_Rank_Map' = pseudo_clone_map), file.path(output_dir, "Clone_ID_Rank_Mapping.xlsx"))
  
} else {
  # --- Treatment / No-Rank Logic ---
  log_print("Pseudo-clone mapping disabled (flag is FALSE or no variable found).")
  
  # If we have a secondary variable (e.g. Treatment), use it for color
  if (!is.null(cfg_vars$secondary_group_var) && cfg_vars$secondary_group_var != 'null') {
    color_var <- cfg_vars$secondary_group_var
    log_print(paste("Plots will be colored by secondary variable:", color_var))
  } else {
    # Fallback to Primary (Genotype)
    color_var <- primary_var
    log_print(paste("Plots will be colored by primary variable:", color_var))
  }
}

## --- 4B. Pre-Exclusion QC Plots ---####
log_print("\n--- Starting PART 4B: Generating Pre-Exclusion QC Plots ---", console = TRUE)

# --- FIX: Explicitly define variables from config to prevent missing object errors ---
time_var <- config$key_variables$time_variable
primary_var <- config$key_variables$primary_group_var
re_cross <- config$key_variables$optional_crossed_effect

# Define key variables for plotting
time_var_sym <- sym(time_var)
primary_var_sym <- sym(primary_var)
peak_height_cutoff <- config$parameters$peak_height_threshold %||% 0

# Check if we have a crossed effect (Batch) to plot
# Ensure re_cross is NULL if it's the string "null"
if (!is.null(re_cross) && re_cross == 'null') re_cross <- NULL

qc_plot_data <- data_with_reasons %>%
  dplyr::filter(  # Force R to use the dplyr version
    !is.na(as.character(.data[[primary_var]])) &
      !grepl("^NA$|^null$", as.character(.data[[primary_var]]), ignore.case = TRUE)
  ) %>%
  dplyr::mutate(  # Force R to use the dplyr version
    is_excluded = !is.na(exclusion_reason),
    "{time_var}" := as.numeric(.data[[time_var]])
  )

# --- Create custom Y-axis log scale ---
y_breaks_config <- config$plot_settings$y_axis_log_breaks
if (!is.null(y_breaks_config) && all(y_breaks_config != 'null')) {
  log_print("...using custom Y-axis log breaks from config.")
  y_vals <- as.numeric(unlist(y_breaks_config))
  scale_y_log_custom <- scale_y_continuous(breaks = log10(y_vals), labels = y_vals)
} else {
  log_print("...using default Y-axis log breaks.")
  scale_y_log_custom <- scale_y_continuous(labels = scales::math_format(10^.x))
}


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot 1: Peak Height vs. Time
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
log_print("...plotting peak height vs. time (pre-exclusion)...")

# 1. Define Base Aesthetics
p1_aes <- aes(
  x = !!time_var_sym, 
  y = log10(target_peak_height), 
  color = exclusion_reason
)

# 2. Define Jitter Aesthetics (Add shape if batch exists)
jitter_aes <- aes()
if (!is.null(re_cross)) {
  jitter_aes <- aes(shape = !!sym(re_cross))
}

p_height_vs_time <- ggplot(qc_plot_data, p1_aes) +
  # Use dynamic jitter aesthetics
  geom_jitter(mapping = jitter_aes, width = 0.1, alpha = 0.7) +
  scale_y_log_custom + 
  facet_wrap(vars(!!primary_var_sym)) +
  labs(
    title = "Pre-Exclusion: Peak Height vs. Time",
    subtitle = "Points are colored by exclusion reason.",
    x = stringr::str_to_title(time_var),
    y = "Peak Height (Log10 Scale)",
    color = "Exclusion Reason",
    shape = if(!is.null(re_cross)) "Batch" else NULL
  ) +
  theme_publication() +
  theme(legend.position = "bottom", legend.text = element_text(size = 8))

print(p_height_vs_time)
ggsave(
  file.path(diagnostics_dir, "p_qc_peak_height_vs_time.tiff"),
  p_height_vs_time, width = 12, height = 9, dpi = 300, compression = "lzw"
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot 2: Peak Height vs. Mode
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
log_print("...plotting peak height vs. mode (good vs. excluded)...")

# 1. Define Base Point Aesthetics
# We always color by primary_var and shape by exclusion status
p2_point_aes <- aes(
  color = !!primary_var_sym, 
  shape = is_excluded
)

# 2. Add 'group' aesthetic for proper dodging if batch exists
if (!is.null(re_cross)) {
  # If batch exists, grouping by it helps jitter stay consistent
  p2_point_aes <- utils::modifyList(p2_point_aes, aes(group = !!sym(re_cross)))
}

p_height_vs_mode <- ggplot(qc_plot_data, aes(x = mode, y = log10(target_peak_height))) +
  geom_point(
    mapping = p2_point_aes,
    alpha = 0.6,
    position = position_jitter(width=0.1)
  ) +
  geom_smooth(data = function(df) dplyr::filter(df, is_excluded == FALSE), 
              method = "lm", se = FALSE, color = "black", linetype = "solid") +
  geom_smooth(data = function(df) dplyr::filter(df, is_excluded == TRUE), 
              method = "lm", se = FALSE, color = "red", linetype = "dashed") +
  # Stat correlations
  ggpubr::stat_cor(data = function(df) dplyr::filter(df, is_excluded == FALSE), 
                   aes(label = paste("R(good) =", ..r..)),
                   method = "pearson", label.x.npc = 0.05, label.y.npc = 0.95, color = "black", size = 3) +
  ggpubr::stat_cor(aes(label = paste("R(overall) =", ..r..)),
                   method = "pearson", label.x.npc = 0.05, label.y.npc = 0.85, color = "red", size = 3) +
  scale_y_log_custom + 
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19), name = "Is Excluded") + 
  facet_wrap(vars(!!primary_var_sym), scales = "free_x") +
  labs(
    title = "Pre-Exclusion: Peak Height vs. CAG Length",
    x = "Modal CAG Size",
    y = "Peak Height (Log10 Scale)",
    color = stringr::str_to_title(primary_var)
  ) +
  theme_publication()

print(p_height_vs_mode)
ggsave(
  file.path(diagnostics_dir, "p_qc_peak_height_vs_mode.tiff"),
  p_height_vs_mode, width = 12, height = 9, dpi = 300, compression = "lzw"
)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Plot 3: Global Peak Height (with Cutoff)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
log_print("...plotting global peak height vs. mode (with cutoff)...")

# Use the same flexible aesthetics as Plot 2
# Note: shape is mapped to is_excluded, but we might want batch grouping if present
p3_point_aes <- aes(shape = is_excluded)
if (!is.null(re_cross)) {
  p3_point_aes <- utils::modifyList(p3_point_aes, aes(group = !!sym(re_cross)))
}

p_height_global <- ggplot(qc_plot_data, 
                          aes(x = mode, y = log10(target_peak_height), color = !!primary_var_sym)) +
  geom_hline(yintercept = log10(peak_height_cutoff), color = "red", linetype = "dashed", linewidth = 1) +
  annotate("text", x = max(qc_plot_data$mode, na.rm = TRUE), y = log10(peak_height_cutoff), 
           label = paste("Cutoff =", peak_height_cutoff), hjust = 1.1, vjust = -0.5, color = "red", size = 3.5) +
  
  geom_point(
    mapping = p3_point_aes, 
    alpha = 0.6,
    position = position_jitter(width=0.1)
  ) +
  
  scale_y_log_custom + 
  scale_shape_manual(values = c("TRUE" = 4, "FALSE" = 19), name = "Is Excluded") +
  labs(
    title = "Global Peak Height vs. CAG Length",
    subtitle = paste("Dashed red line shows the peak height cutoff of", peak_height_cutoff),
    x = "Modal CAG Size",
    y = "Peak Height (Log10 Scale)",
    color = stringr::str_to_title(primary_var)
  ) +
  theme_publication()

print(p_height_global)
ggsave(
  file.path(diagnostics_dir, "p_qc_peak_NET_global.tiff"),
  p_height_global, width = 11, height = 8, dpi = 300, compression = "lzw"
)


#=============================================================================#
# PART 5: EXPORT DATA AND SAVE ENVIRONMENT ####
#=============================================================================#
log_print("\n--- Starting PART 5: Exporting Data Tables ---", console = TRUE)

all_bio_outliers <- bind_rows(bio_outliers_list, .id = "response_variable")
all_tech_outliers <- bind_rows(tech_outliers_list, .id = "response_variable")

# --- UPDATE EXCLUDED_DATA WITH STATISTICAL OUTLIERS ---
# If the auto-remove flags were set to TRUE, append the dropped rows to the 
# master excluded_data log so downstream scripts see them.

if (isTRUE(config$parameters$remove_bio_outliers) && nrow(all_bio_outliers) > 0) {
  bio_excl_formatted <- all_bio_outliers %>%
    dplyr::mutate(exclusion_reason = paste("Bio Outlier (IQR):", response_variable)) %>%
    # CRITICAL FIX: Force all columns to character to prevent bind_rows type crashes
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character)) 
  
  excluded_data <- dplyr::bind_rows(
    excluded_data %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.character)), 
    bio_excl_formatted
  )
}

if (isTRUE(config$parameters$remove_tech_outliers) && nrow(all_tech_outliers) > 0) {
  tech_excl_formatted <- all_tech_outliers %>%
    dplyr::mutate(exclusion_reason = paste("Tech Outlier (IQR):", response_variable)) %>%
    # CRITICAL FIX: Force all columns to character to prevent bind_rows type crashes
    dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
  
  excluded_data <- dplyr::bind_rows(
    excluded_data %>% dplyr::mutate(dplyr::across(dplyr::everything(), as.character)), 
    tech_excl_formatted
  )
}

export_list_long <- list(
  'Data_Per_PCR' = data_per_pcr,
  'Summary_Per_Rep' = summary_per_rep,
  'Summary_Baseline_Only' = baselines_table,
  'Summary_By_Group' = summary_by_group,
  'Data_Excluded' = excluded_data,
  'All_Bio_Outliers' = all_bio_outliers,
  'All_Tech_Outliers' = all_tech_outliers
)
export_to_excel(export_list_long, file.path(output_dir, "Complete_Analysis_Data_Long.xlsx"))


# --- Export Wide Formats (GraphPad Prism) ---
log_print("Creating wide-format files for GraphPad Prism.")

# 1. Bio-Replicate Level (Well Reps)
graphpad_export_list_well <- list()

# Define robust ID columns (check what actually exists)
possible_bio_ids <- c("Genotype_Pub", all_grouping_vars, "rep")
valid_bio_ids <- intersect(colnames(summary_per_rep), possible_bio_ids)

# 1. Grab the batch variable from your config
batch_var <- config$key_variables$optional_crossed_effect

# 2. Add it to your valid IDs (only if it actually exists in your config)
if (!is.null(batch_var) && batch_var != 'null') {
  valid_bio_ids <- unique(c(valid_bio_ids, batch_var))
}

for (current_resp_var in response_vars) {
  
  if(!current_resp_var %in% colnames(summary_per_rep)) next
  
  temp_wide_df <- summary_per_rep %>%
    dplyr::select(dplyr::any_of(valid_bio_ids), !!sym(time_var), value = !!sym(current_resp_var)) %>%
    dplyr::group_by(dplyr::across(dplyr::any_of(intersect(valid_bio_ids, colnames(.))))) %>%
    dplyr::arrange(!!sym(time_var), .by_group = TRUE) %>%
    ungroup() %>%
    tidyr::pivot_wider(
      id_cols = dplyr::all_of(valid_bio_ids),
      names_from = !!sym(time_var), 
      values_from = value
    ) %>%
    # Reorder columns: IDs first, then Timepoints ordered numerically
    dplyr::select(
      dplyr::all_of(valid_bio_ids),
      dplyr::all_of({
        # Find all columns that aren't IDs (these are your timepoints)
        time_cols <- setdiff(colnames(.), valid_bio_ids)
        # Sort them numerically so Day 10 comes after Day 2, not Day 1
        time_cols[order(as.numeric(time_cols))]
      })
    )
  graphpad_export_list_well[[current_resp_var]] <- temp_wide_df
}
export_to_excel(graphpad_export_list_well, file.path(output_dir, "GraphPad_Data_Well_Reps.xlsx"))


# 2. Tech-Replicate Level (PCR Reps)
graphpad_export_list_pcr <- list()

# Define robust ID columns for PCR level
possible_pcr_ids <- c("Genotype_Pub", all_grouping_vars, "rep", "pcr")
valid_pcr_ids <- intersect(colnames(data_per_pcr), possible_pcr_ids)

for (current_resp_var in response_vars) {
  
  if(!current_resp_var %in% colnames(data_per_pcr)) next
  
  temp_wide_df_pcr <- data_per_pcr %>%
    # Create temp ID for pivoting (ensures uniqueness)
    #dplyr::mutate(pcr_pivot_id = paste(!!!syms(valid_pcr_ids), sep="_")) %>% 
    dplyr::select(dplyr::all_of(valid_pcr_ids), !!sym(time_var), value = !!sym(current_resp_var)) %>%
    dplyr::group_by(Genotype_Pub) %>%
    dplyr::arrange(!!sym(time_var), .by_group = TRUE) %>%
    ungroup() %>%
    tidyr::pivot_wider(
      id_cols = c(dplyr::all_of(valid_pcr_ids)), # Use dynamic list here
      names_from = !!sym(time_var), 
      values_from = value,
      values_fn = mean,
      values_fill = NA 
    ) %>%
    #dplyr::select(-c(fsa_filename)) %>% # Remove temp ID 
    dplyr::select(
      dplyr::all_of(intersect(valid_pcr_ids, colnames(.))),
      dplyr::all_of({
        cn <- setdiff(colnames(.), valid_pcr_ids)
        cn[order(suppressWarnings(as.numeric(cn)))]
      })
    )
  
  graphpad_export_list_pcr[[current_resp_var]] <- temp_wide_df_pcr
}
export_to_excel(graphpad_export_list_pcr, file.path(output_dir, "GraphPad_Data_PCR_Reps.xlsx"))

cfg_vars <- config$key_variables
primary_var <- cfg_vars$primary_group_var

# Create label_lookup for downstream scripts (ensures Script 05 can run independently of 02)
label_lookup <- genotype_map %>%
  dplyr::select(
    !!sym(config$key_variables$primary_group_var),
    dplyr::any_of(c("Label_Pub", "Label_Exp", "Genotype_Pub", "Genotype_Exp"))
  ) %>%
  dplyr::distinct()

# --- Save Final Environment ---
rdata_path <- file.path(output_dir, "processing_complete.RData")
save.image(file = rdata_path) 
log_print(paste("R environment saved to:", rdata_path))
log_print("This file is required for 02_generate_plots.R")

# --- Finalize ---
log_print("\n--- SCRIPT 01: DATA PROCESSING FINISHED SUCCESSFULLY ---", console = TRUE)
log_close()
