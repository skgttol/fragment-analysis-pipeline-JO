#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                             #
#         Advanced Fragment Analysis Workflow - Script 2 of 5                 #
#                   (Statistical Modeling & Plotting)                         #
#                                                                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# DESCRIPTION:
# This script is OPTIONAL. It loads the "processing_complete.RData"
# environment saved by '01_process_data.R' and performs all statistical
# modeling and plot generation.
#
# WORKFLOW:
#  - PART 0: Setup and Load Processed Data
#  - PART 1: Generate Baseline Data, Palettes, & Slope Table
#  - PART 2: Generate Descriptive Plots
#    - PLOT 1: Group Averages
#    - PLOT 2: Detailed View
#    - PLOT 3: Individual Clone Profiles (Faceted by Genotype)
#    - PLOT 4: (NEW) Individual Rep Profiles (Faceted by Clone)
#    - PLOT 2B: Plate QC Heatmaps
#    - PLOT 2C: Replication Counts
#  - PART 3: Perform Statistical Modeling
#  - PART 4: Generate Model-Based Plots
#  - PART 5: Export Statistical Model Results & Save RData
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#


#=============================================================================#
# PART 0: SETUP AND LOAD PROCESSED DATA ####
#=============================================================================#

## 0A. Load Libraries ----
#------------------------#
packages <- c("tidyverse", "readxl", "here", "openxlsx", "broom.mixed",
              "ggpubr", "janitor", "yaml", "logr", "lme4", "lmerTest",
              "emmeans", "ggeffects", "patchwork", "gridExtra", "grid",
              "RColorBrewer", "gtable", "ggnewscale", "GGally", "reshape2", "MASS", "cowplot")
# Install missing
installed_packages <- packages %in% rownames(utils::installed.packages())
if (any(installed_packages == FALSE)) {
  utils::install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

## 0B. Load Functions ----
#-------------------------#
if (!file.exists(here::here("functions.R"))) {
  stop("CRITICAL ERROR: 'functions.R' not found. Please ensure it is in the main project directory.")
}
source(here::here("functions.R"))

# if (!file.exists("C:/Users/skgttol/OneDrive - University College London/PhD/PhD_Thesis/02_Data/Template_RScripts/CAGsizing_Flexible/functions.R")) {
#   stop("CRITICAL ERROR: 'functions.R' not found. Please ensure it is in the main project directory.")
# }
# source("C:/Users/skgttol/OneDrive - University College London/PhD/PhD_Thesis/02_Data/Template_RScripts/CAGsizing_Flexible/functions.R")

## 0C. Load Configuration and Find Latest Data ----
#--------------------------------------------------#
if (!file.exists(here::here("config.yml"))) {
  stop("CRITICAL ERROR: 'config.yml' not found. Please create it before running.")
}
  config <- yaml::read_yaml(here::here("config.yml"))

# --- DETERMINE EXPERIMENT STYLE (THE SWITCH) ---
exp_style <- config$key_variables$experiment_style
if (is.null(exp_style) || exp_style == 'null') {
  exp_style <- "clone" # Fallback so older non-treatment configs never break
}
is_treatment_exp <- (exp_style == "treatment")
logr::log_print(paste("Experiment Plotting Style set to:", stringr::str_to_title(exp_style)))

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
if (!file.exists(rdata_path)) {
  stop("ERROR: 'processing_complete.RData' not found in:", latest_analysis_dir)
}

# 1. Load the saved environment from Script 01 FIRST.
load(rdata_path)
print(paste("Data loaded successfully from:", rdata_path))

# This loads all data, config, output_dir, palettes, etc.
source("C:/Users/skgttol/OneDrive - University College London/PhD/PhD_Thesis/02_Data/Template_RScripts/CAGsizing_Flexible/functions.R")
config <- yaml::read_yaml(here::here("config.yml"))
# 2. NOW, (re)open the log file.
try(logr::log_close(), silent = TRUE)
log_path <- file.path(output_dir, "analysis_log.log") # output_dir is loaded from .RData
lf <- logr::log_open(log_path, show_notes = FALSE, logdir = FALSE) # Append

logr::log_print("\n\n--- SCRIPT 02: MODELING & PLOTTING INITIALIZED ---", console = TRUE)
logr::log_print("Log re-opened in append mode.")

# --- (NEW) INITIALIZE PLOT DATABASE ---
# If plot_database doesn't exist from Script 01, create it.
if (!exists("plot_database")) {
  plot_database <- list()
}

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

logr::log_print(paste("Default plot label style:", config$parameters$active_label_style %||% "exploratory"))

#=============================================================================#
# PART 1: GENERATE PLOTTING DATA, PALETTES, & SLOPE TABLE                  ####
#=============================================================================#
logr::log_print("\n--- Starting PART 1: Creating Plotting Data & Palettes ---", console = TRUE)

# --- CRITICAL FIX: Ensure 'clone_rank' exists in all source data ---
# We do this BEFORE creating baseline data so the rank propagates to modeling_data
if (exists("pseudo_clone_map") && has_pseudo_clones) {
  
  logr::log_print("...Patching missing 'clone_rank' into source dataframes.")
  
  # Define join keys based on config
  cfg_vars <- config$key_variables
  primary_var <- cfg_vars$primary_group_var
  clone_var <- cfg_vars$repeated_measure_var
  
  # The map uses 'original_clone_id', the data uses whatever 'clone_var' is (e.g. "clone")
  join_key <- c(primary_var, setNames("original_clone_id", clone_var)) 
  
  # 1. Fix summary_per_rep
  if (!"clone_rank" %in% colnames(summary_per_rep)) {
    summary_per_rep <- summary_per_rep %>%
      mutate(!!sym(clone_var) := as.character(!!sym(clone_var))) %>%
      left_join(pseudo_clone_map, by = join_key)
  }
  
  # 2. Fix data_per_pcr
  if (!"clone_rank" %in% colnames(data_per_pcr)) {
    data_per_pcr <- data_per_pcr %>%
      mutate(!!sym(clone_var) := as.character(!!sym(clone_var))) %>%
      left_join(pseudo_clone_map, by = join_key)
  }
}

## --- 1A. Create specialized data for plotting ----
plot_data_list <- create_baseline_data(
  data_per_pcr = data_per_pcr,
  summary_per_rep = summary_per_rep,
  config = config,
  all_grouping_vars = all_grouping_vars
)
modeling_data <- plot_data_list$modeling_data
shared_baseline_pcr_points <- plot_data_list$shared_baseline_pcr_points

if (exists("pseudo_clone_map") && has_pseudo_clones) {
  if (!"clone_rank" %in% colnames(modeling_data)) {
    modeling_data <- modeling_data %>%
      mutate(!!sym(clone_var) := as.character(!!sym(clone_var))) %>%
      left_join(pseudo_clone_map, by = join_key)
  }
  if (!"clone_rank" %in% colnames(shared_baseline_pcr_points)) {
    shared_baseline_pcr_points <- shared_baseline_pcr_points %>%
      mutate(!!sym(clone_var) := as.character(!!sym(clone_var))) %>%
      left_join(pseudo_clone_map, by = join_key)
  }
}


# --- NEW: GLOBAL VALUE RENAMING (e.g., FKO -> FAN1 KO) ---
if ("value_renaming" %in% names(config) && length(config$value_renaming) > 0) {
  logr::log_print("...Applying global value renaming from config.yml")
  
  
  # Convert the config list to a named vector: c("Old" = "New")
  rename_map <- unlist(config$value_renaming)
  # Helper function to safely find and replace values while preserving factor order
  apply_renaming <- function(df) {
    if (is.null(df) || nrow(df) == 0) return(df)
    
    df %>%
      dplyr::mutate(dplyr::across(
        dplyr::where(~ is.character(.) || is.factor(.)),
        ~ {
          if (is.factor(.)) {
            # 1. Capture the original ordered levels
            old_levels <- levels(.)
            
            # 2. Translate the levels themselves using the map
            new_levels <- dplyr::coalesce(unname(rename_map[old_levels]), old_levels)
            
            # 3. Translate the actual data values
            vec <- as.character(.)
            new_vec <- dplyr::coalesce(unname(rename_map[vec]), vec)
            
            # 4. Rebuild the factor using the translated, ordered levels
            # unique() prevents errors if two different old levels map to the same new name
            factor(new_vec, levels = unique(new_levels))
            
          } else {
            # For standard character columns, just replace the text
            vec <- as.character(.)
            dplyr::coalesce(unname(rename_map[vec]), vec)
          }
        }
      ))
  } 
  
  # Apply to all master dataframes
  modeling_data <- apply_renaming(modeling_data)
  summary_per_rep <- apply_renaming(summary_per_rep)
  data_per_pcr <- apply_renaming(data_per_pcr)
  if (exists("shared_baseline_pcr_points")) shared_baseline_pcr_points <- apply_renaming(shared_baseline_pcr_points)
}

# NEW: APPLY GLOBAL FACTOR LEVELS FROM CONFIG TO ALL DATAFRAMES 
if ("factor_levels" %in% names(config)) {
  logr::log_print("...Applying global factor levels from config.yml")
  
  # Grab the rename map so we can translate the config levels to match the data
  rename_map <- if ("value_renaming" %in% names(config)) unlist(config$value_renaming) else NULL
  
  for (col_name in names(config$factor_levels)) {
    # Get the explicitly ordered vector from the config
    raw_levels <- as.character(config$factor_levels[[col_name]])
    
    # Translate the config levels using the rename map so they match!
    ordered_levels <- if (!is.null(rename_map)) {
      dplyr::coalesce(unname(rename_map[raw_levels]), raw_levels)
    } else {
      raw_levels
    }
    
    # Force unique to prevent crashes if two old levels map to one new level
    ordered_levels <- unique(ordered_levels)
    
    # Safely apply it to every master dataframe if the column exists
    if (col_name %in% names(modeling_data)) {
      modeling_data[[col_name]] <- factor(as.character(modeling_data[[col_name]]), levels = ordered_levels)
    }
    if (col_name %in% names(summary_per_rep)) {
      summary_per_rep[[col_name]] <- factor(as.character(summary_per_rep[[col_name]]), levels = ordered_levels)
    }
    if (col_name %in% names(data_per_pcr)) {
      data_per_pcr[[col_name]] <- factor(as.character(data_per_pcr[[col_name]]), levels = ordered_levels)
    }
    if (exists("shared_baseline_pcr_points") && col_name %in% names(shared_baseline_pcr_points)) {
      shared_baseline_pcr_points[[col_name]] <- factor(as.character(shared_baseline_pcr_points[[col_name]]), levels = ordered_levels)
    }
  }
}
# Apply this to all master dataframes
primary_var_name <- config$key_variables$primary_group_var
modeling_data <- sync_label_factors(modeling_data, primary_var_name)
summary_per_rep <- sync_label_factors(summary_per_rep, primary_var_name)
data_per_pcr <- sync_label_factors(data_per_pcr, primary_var_name)
if (exists("shared_baseline_pcr_points")) {
  shared_baseline_pcr_points <- sync_label_factors(shared_baseline_pcr_points, primary_var_name)
}
logr::log_print("Derived label columns synced to primary variable factor order.")

## --- 1B. Define color variables and distinct palettes ----
cfg_vars <- config$key_variables

# (NEW) Define batch variable
re_cross <- cfg_vars$optional_crossed_effect
is_crossed_model <- (!is.null(re_cross) && re_cross != 'null')
if(is_crossed_model) {
  logr::log_print(paste("...Found crossed effect (e.g., batch):", re_cross, ". Will add as 'shape' to plots."))
}

# (NEW) Define rep variable
rep_var <- cfg_vars$optional_grouping_var %||% cfg_vars$repeated_measure_var
if (is.null(rep_var) || rep_var == 'null') { rep_var <- "rep" }

# --- 1. Primary Palette (e.g., Genotype) ---
primary_palette <- create_custom_palette(modeling_data, config, cfg_vars$primary_group_var)
logr::log_print(paste("Primary palette generated for:", cfg_vars$primary_group_var))

# --- 2. Secondary Palette (e.g., Treatment, Clone) ---
secondary_palette <- NULL
secondary_var <- cfg_vars$secondary_group_var
if (!is.null(secondary_var) && secondary_var != "null") {
  secondary_palette <- create_custom_palette(modeling_data, config, secondary_var)
  
  # Fallback: If no color mapping in config, assign Set2 to secondary groups
  if (is.null(config$color_mapping) && length(secondary_palette) > 0) {
    n_sec_levels <- length(secondary_palette)
    secondary_palette[] <- RColorBrewer::brewer.pal(max(3, n_sec_levels), "Set2")[1:n_sec_levels]
  }
  logr::log_print(paste("Secondary palette generated for:", secondary_var))
}

# --- 3. Rank-Based Palette & Legend (For Pseudo-Clones) ---
custom_palette <- NULL
p_legend <- NULL 

if (has_pseudo_clones) {
  logr::log_print("Creating rank-based color palette for clones.")
  max_rank <- max(as.numeric(levels(droplevels(pseudo_clone_map$clone_rank))))
  n_colors_to_request <- max(3, min(9, max_rank))
  base_palette <- RColorBrewer::brewer.pal(n_colors_to_request, "Set1")
  
  if (max_rank > length(base_palette)) {
    final_rank_palette <- rep(base_palette, length.out = max_rank)
  } else {
    final_rank_palette <- base_palette[1:max_rank]
  }
  names(final_rank_palette) <- as.character(1:max_rank)
  custom_palette <- final_rank_palette # Retained for legacy plot logic
  logr::log_print(paste("Rank palette created with", max_rank, "colors."))
  
  # --- Create the Custom Hierarchical Legend Plot ---
  logr::log_print("Building custom hierarchical legend...")
  legend_data <- pseudo_clone_map %>%
    dplyr::arrange(!!rlang::sym(cfg_vars$primary_group_var), clone_rank) %>%
    dplyr::mutate(original_clone_id = factor(original_clone_id, levels = rev(unique(.$original_clone_id))))
  
  p_legend <- ggplot(legend_data,
                     aes(x = 0, y = original_clone_id, color = clone_rank)) +
    geom_point(size = 4) +
    geom_text(aes(label = original_clone_id), x = 0.1, hjust = 0, color = "black", size = 3.5) +
    scale_color_manual(values = custom_palette) +
    facet_grid(rows = vars(!!rlang::sym(cfg_vars$primary_group_var)),
               scales = "free_y", space = "free_y", switch = "y") +
    xlim(-0.1, 1) +
    theme_minimal(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid = element_blank(),
      axis.text.x = element_blank(),
      axis.title = element_blank(),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      strip.text.y.left = element_text(
        angle = 0, face = "bold.italic", size = 11, hjust = 0.5
      ),
      plot.margin = margin(t = 20)
    )
  if (!is.null(p_legend)) {
    plot_database[["p_hierarchical_legend"]] <- p_legend
  }
} else {
  # Legacy fallback for older plot logic if pseudo_clones are off
  custom_palette <- if(!is.null(secondary_palette)) secondary_palette else primary_palette
}

# Get plot settings
plot_settings <- config$plot_settings
x_breaks <- plot_settings$x_axis_breaks

#=============================================================================#
# PART 2: GENERATE DESCRIPTIVE PLOTS                                      #####
#=============================================================================#
logr::log_print("\n--- Starting PART 2: Generating Descriptive Plots ---", console = TRUE)
plots_dir <- file.path(output_dir, "02_descriptive_plots")
dir.create(plots_dir, showWarnings = FALSE)


for (resp_var in response_vars) {
  y_axis_label <- response_labels[[resp_var]] %||% stringr::str_to_title(resp_var)
  short_resp_var <- response_shortnames[[resp_var]] %||% resp_var
  y_scale_setting <- if (resp_var == "mode") "free_y" else "fixed"
  
  logr::log_print(paste("...generating descriptive plots for", resp_var))
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## PLOT 1: Group / Treatment Averages & Trendline
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  if (is_treatment_exp || is_clone_style_analysis) {
    logr::log_print(paste("...generating averages plot for", resp_var))
    label_exp <- get_lbl("exploratory")
    
    if (is_treatment_exp) {
      # --- TREATMENT STYLE PLOT ---
      color_col <- config$key_variables$secondary_group_var # e.g., "treatment"
      facet_vars <- list(rlang::sym(label_exp))
      if (is_crossed_model) facet_vars <- append(facet_vars, list(rlang::sym(re_cross)))
      
      plot_title <- paste("Treatment Averages:", y_axis_label)
      plot_sub <- "Lines are linear model fits. Points are treatment means +/- SEM."
      
    } else {
      # ---  CLONE STYLE PLOT ---
      color_col <- label_exp
      facet_vars <- list(rlang::sym(label_exp))
      plot_title <- paste("Group Averages & Trendline:", y_axis_label)
      plot_sub <- "Lines are linear model fits. Points are group means +/- SEM."
    }
    
    p_avg <- ggplot(
      data = modeling_data,
      aes(x = !!sym(time_var), y = !!sym(resp_var), color = !!sym(color_col), group = !!sym(color_col))
    ) +
      # Zero reference line for change-from-baseline variables
      { if (grepl("change|_change", resp_var, ignore.case = TRUE)) geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) } +
      stat_summary(fun.data = "mean_se", geom = "errorbar", width = 1, linewidth = 0.7) +
      stat_summary(fun = mean, geom = "point", size = 3) +
      geom_smooth(method = "lm", se = TRUE, linewidth = 1, aes(fill = !!sym(color_col)), alpha = 0.15) +
      facet_wrap(facet_vars, scales = y_scale_setting, axes = "all") +
      labs(title = plot_title, subtitle = plot_sub,
           x = stringr::str_to_title(time_var), y = y_axis_label,
           color = stringr::str_to_title(color_col),
           fill  = stringr::str_to_title(color_col)) +
      theme_publication() +
      theme(legend.position = if(is_treatment_exp) "right" else "none") 
    
    # Pass the plot and the color variable to our universal helper
    p_avg <- apply_smart_palette(p_avg, color_col)
    
    if (!is.null(x_breaks) && all(x_breaks != 'null')) p_avg <- p_avg + scale_x_continuous(breaks = as.numeric(x_breaks))
    
    print(p_avg)
    ggsave(file.path(plots_dir, paste0("01_avg_", short_resp_var, ".tiff")), p_avg, width = if(is_treatment_exp) 12 else 10, height = 7, dpi = 300, device = 'tiff', compression = "lzw")
    plot_database[[paste0("p_avg_", short_resp_var)]] <- p_avg
    assign(paste0("p_avg_", short_resp_var), p_avg)
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  ## PLOT 2: Detailed View with All Individual Replicates
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  logr::log_print(paste("...generating detailed view plot for", resp_var))
  label_exp <- get_lbl("exploratory")
  
  if (is_treatment_exp) {
    color_col_string <- config$key_variables$secondary_group_var # "treatment"
    facet_vars <- list(rlang::sym(label_exp))
    if (is_crossed_model) facet_vars <- append(facet_vars, list(rlang::sym(re_cross)))
  } else {
    color_col_string <- if(has_pseudo_clones) "clone_rank" else label_exp
    facet_vars <- list(rlang::sym(label_exp))
  }
  
  p_detail <- ggplot(
    data = modeling_data,
    aes(x = !!sym(time_var), y = !!sym(resp_var))
  ) 
  
  # 1. Add Jitter Points (Raw Data)
  base_jitter_aes <- aes(y = !!sym(resp_var), color = !!sym(color_col_string))
  final_jitter_aes <- create_plot_aes(base_jitter_aes, re_cross)
  
  p_detail <- p_detail + geom_jitter(
    data = dplyr::bind_rows(data_per_pcr, shared_baseline_pcr_points),
    mapping = final_jitter_aes,
    width = 0.15, alpha = 0.18, size = 1
    )
  
  # 2. Add Mean Points (Bio-Reps)
  base_point_aes <- aes(y = !!sym(resp_var), color = !!sym(color_col_string))
  final_point_aes <- create_plot_aes(base_point_aes, re_cross)
  
  p_detail <- p_detail + geom_point(
    data = summary_per_rep,
    mapping = final_point_aes,
    alpha = 0.85, size = 2.5
    )
  
  # 3. Add Trendlines & Faceting

    p_detail <- p_detail + 
    geom_smooth(
      method = "lm", se = TRUE, linewidth = 1.1,
      aes(color = !!sym(color_col_string), fill = !!sym(color_col_string)), alpha = 0.15) +
    facet_wrap(facet_vars, scales = y_scale_setting, axes = "all") +
    labs(
      title = paste("Detailed View with All Replicates:", y_axis_label),
      subtitle = "Faint points = individual PCRs. Bold points = biological replicates. Band = 95% CI.",
      x = stringr::str_to_title(time_var),
      y = y_axis_label, 
      color = stringr::str_to_title(color_col_string),
      fill  = stringr::str_to_title(color_col_string),
      shape = if (!is.null(re_cross) && re_cross != 'null') stringr::str_to_title(re_cross) else NULL
    ) +
    theme_publication() +
    theme(legend.position = "right")
  
    p_detail <- apply_smart_palette(p_detail, color_col_string)
  
  if (!is.null(x_breaks) && all(x_breaks != 'null')) {
    p_detail <- p_detail + scale_x_continuous(breaks = as.numeric(x_breaks))
  }
  
  assign(paste0("p_detail_", short_resp_var), p_detail)
  plot_database[[paste0("p_detail_", short_resp_var)]] <- p_detail
  
  print(p_detail)
  ggsave(
    file.path(plots_dir, paste0("02_detail_", short_resp_var, ".tiff")),
    p_detail, width = if(is_treatment_exp) 12 else 10, height = 8, dpi = 300, device = 'tiff',
    compression = "lzw"
  )  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 3: (Conditional) Individual Profiles (Faceted by Genotype)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  secondary_var <- cfg_vars$secondary_group_var
  if (!is.null(secondary_var) && secondary_var != 'null') {
    logr::log_print(paste("...generating individual profiles plot for", secondary_var))
    
    p_traj <- ggplot(
      data = summary_per_rep,
      aes(
        x = !!sym(time_var),
        y = !!sym(resp_var),
        color = if(has_pseudo_clones) clone_rank else !!sym(color_var),
        group = !!sym(color_var)
      )
    ) +
      { if (grepl("change|_change", resp_var, ignore.case = TRUE)) geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) } +
      geom_line(alpha = 0.7, linewidth = 0.7) +
      geom_point(
        mapping = create_plot_aes(aes(color = if(has_pseudo_clones) clone_rank else !!sym(color_var)), re_cross),
        alpha = 0.85, size = 2
      )
    
    p_traj <- p_traj + 
      facet_wrap(vars(!!sym(label_exp)), scales = y_scale_setting, axes = "all") +
      labs(
        title = paste("Individual Profiles:", y_axis_label),
        subtitle = paste("Each colored line represents a unique", stringr::str_to_title(secondary_var)),
        x = stringr::str_to_title(time_var),
        y = y_axis_label
      ) +
      theme_publication(base_size = 14) +
      theme(legend.position = "none")
    
    p_traj <- apply_smart_palette(p_traj, if(has_pseudo_clones) "clone_rank" else color_var)
    
    if (!is.null(x_breaks) && all(x_breaks != 'null')) {
      p_traj <- p_traj + scale_x_continuous(breaks = as.numeric(x_breaks))
    }
    
    if(is_clone_style_analysis && !is.null(p_legend)) {
      combined_plot_traj <- p_traj + p_legend + patchwork::plot_layout(widths = c(3, 1))
      plot_width_traj <- 16
    } else {
      combined_plot_traj <- p_traj
      plot_width_traj <- 12
    }
    
    # Store the plot object
    assign(paste0("p_profiles_", short_resp_var), combined_plot_traj)
    # --- Add to plot database ---
    plot_database[[paste0("p_profiles_", short_resp_var)]] <- combined_plot_traj
    
    print(combined_plot_traj)
    ggsave(
      file.path(plots_dir, paste0("03_profiles_", short_resp_var, ".tiff")),
      combined_plot_traj, width = plot_width_traj, height = 8, dpi = 300, device = 'tiff',
      compression = "lzw"
    )
  }
}

# --- 2B. (MODIFIED) Plate QC Heatmaps ---
logr::log_print("...generating plate QC heatmaps (one plot per metric).")

# (NEW) Loop by RESPONSE VARIABLE first
for (current_resp_var in response_vars) {
  
  y_axis_label <- response_labels[[current_resp_var]] %||% stringr::str_to_title(current_resp_var)
  short_resp_var <- response_shortnames[[current_resp_var]] %||% current_resp_var
  logr::log_print(paste("...generating heatmap for metric:", y_axis_label))
  
  # (NEW) Calculate the global min/max for THIS variable
  var_limits_data <- data_per_pcr %>%
    dplyr::summarise(
      min_val = min(!!sym(current_resp_var), na.rm = TRUE),
      max_val = max(!!sym(current_resp_var), na.rm = TRUE)
    )
  val_limits <- c(var_limits_data$min_val, var_limits_data$max_val)
  
  # (NEW) Check for valid data
  if (is.infinite(val_limits[1]) || is.infinite(val_limits[2])) {
    logr::log_print(paste("...WARNING: Skipping heatmap for '", y_axis_label, "' because it contains no valid (non-NA) data."))
    logr::log_print("...Please check your 'config.yml' response_variables and your external data file.")
    next # Skip to the next variable
  }
  
  # (NEW) Call the heatmap function ONCE for this metric
  # This now uses 'plate_raw' and 'well_raw' from functions.R
  p_heatmap_metric <- generate_plate_heatmap(
    data_to_plot = data_per_pcr, # Pass ALL data
    config = config,
    current_resp_var = current_resp_var, # Tell it which variable to plot
    current_label = y_axis_label,
    val_limits = val_limits # Pass the scale limits
  )
  
  # (NEW) DYNAMIC LAYOUT CALCULATION
  n_plates <- length(unique(data_per_pcr$plate_raw)) # Use plate_raw for count
  plot_ncol_est <- ceiling(sqrt(n_plates))
  plot_nrow_est <- ceiling(n_plates / plot_ncol_est)
  base_panel_width <- 6 
  base_panel_height <- 4.5 
  plot_width_dynamic <- (base_panel_width * plot_ncol_est) + 1.5
  plot_height_dynamic <- (base_panel_height * plot_nrow_est) + 0.5
  
  # Store the plot object
  assign(paste0("p_heatmap_", short_resp_var), p_heatmap_metric)
  # --- (NEW) Add to plot database ---
  plot_database[[paste0("p_heatmap_", short_resp_var)]] <- p_heatmap_metric
  
  print(p_heatmap_metric)
  ggsave(
    filename = file.path(diagnostics_dir, paste0("hm_faceted_", short_resp_var, ".tiff")),
    plot = p_heatmap_metric, 
    width = plot_width_dynamic, height = plot_height_dynamic, 
    dpi = 300, device = 'tiff',
    compression = "lzw"
  )
} # --- End of new loop ---

#=============================================================================#
# PART 2C: REPLICATION COUNT SUMMARY (Universal + Bubble)
#=============================================================================#
logr::log_print("\n--- Starting PART 2C: Generating Replication Count Summary ---", console = TRUE)

# --- Explicitly define variables for this block ---
time_var <- config$key_variables$time_variable
clone_var <- config$key_variables$secondary_group_var 
batch_var <- config$key_variables$optional_crossed_effect 
rep_var <- config$key_variables$optional_grouping_var %||% config$key_variables$repeated_measure_var

if (is.null(clone_var) || clone_var == 'null') clone_var <- config$key_variables$primary_group_var
if (is.null(rep_var) || rep_var == 'null') rep_var <- "rep"

diagnostics_dir <- file.path(output_dir, "01_QC_plots")
dir.create(diagnostics_dir, showWarnings = FALSE)

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 1. SIDE-BY-SIDE BAR CHART (Always Run - Aggregates all batches)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
logr::log_print("...Generating standard side-by-side replication plot (Aggregated).")

# Group ONLY by Clone and Time (ignoring Batch to aggregate)
rep_count_summary <- data_per_pcr %>%
  dplyr::group_by(!!sym(clone_var), !!sym(time_var)) %>%
  dplyr::summarise(
    n_pcr_reps = n(), 
    n_well_reps = dplyr::n_distinct(!!sym(rep_var)), 
    .groups = "drop"
  ) %>%
  tidyr::pivot_longer(
    cols = c(n_pcr_reps, n_well_reps), 
    names_to = "Rep_Type_Raw", 
    values_to = "Count"
  ) %>%
  dplyr::mutate(
    Rep_Type = factor(Rep_Type_Raw, 
                      levels = c("n_pcr_reps", "n_well_reps"), 
                      labels = c("Technical (PCR)", "Biological (Well)"))
  )

p_rep_counts <- ggplot(rep_count_summary, aes(x = as.factor(!!sym(time_var)), y = Count, fill = Rep_Type)) +
  geom_col(position = position_dodge(width = 0.8), color = "black") +
  geom_text(
    aes(label = Count, group = Rep_Type), 
    position = position_dodge(width = 0.8), 
    vjust = -0.5, size = 3
  ) +
  facet_wrap(vars(!!sym(clone_var)), scales = "free_x", ncol = 6, axes = "all") +
  scale_fill_brewer(palette = "Pastel1") +
  labs(
    title = "Replication Structure per Clone", 
    subtitle = "Total counts (aggregating all batches).",
    x = stringr::str_to_title(time_var), 
    y = "Count (N)",
    fill = "Unit"
  ) +
  theme_publication() + 
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

print(p_rep_counts)
ggsave(file.path(diagnostics_dir, "p_qc_rep_counts.tiff"), p_rep_counts, width = 18, height = 12, dpi = 300, compression = "lzw")
plot_database[["p_qc_rep_counts"]] <- p_rep_counts


# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# 2. BUBBLE PLOT (Conditional - Only if Batch exists)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
if (is_crossed_model && !is.null(batch_var)) {
  logr::log_print("...Batch variable found. Generating BUBBLE plot for breakdown.")
  
  pcr_counts_by_batch <- data_per_pcr %>%
    dplyr::group_by(!!sym(clone_var), !!sym(time_var), !!sym(batch_var)) %>%
    dplyr::summarise(n_pcr_reps = n(), .groups = "drop")
  
  well_counts_by_batch <- data_per_pcr %>%
    dplyr::distinct(!!sym(clone_var), !!sym(time_var), !!sym(batch_var), !!sym(rep_var)) %>%
    dplyr::group_by(!!sym(clone_var), !!sym(time_var), !!sym(batch_var)) %>%
    dplyr::summarise(n_well_reps = n(), .groups = "drop")
  
  plot_rep_count_data <- pcr_counts_by_batch %>%
    dplyr::full_join(well_counts_by_batch, by = c(clone_var, time_var, batch_var)) %>%
    tidyr::pivot_longer(cols = c(n_pcr_reps, n_well_reps), names_to = "Rep_Type_Raw", values_to = "Count") %>%
    dplyr::filter(!is.na(Count) & Count > 0) %>%
    dplyr::mutate(
      Rep_Type = factor(Rep_Type_Raw, levels = c("n_pcr_reps", "n_well_reps"), labels = c("Technical (PCR)", "Biological (Well)")),
      !!sym(batch_var) := as.factor(!!sym(batch_var))
    )
  
  p_rep_counts_bubble <- ggplot(plot_rep_count_data, aes(x = as.factor(!!sym(time_var)), y = Rep_Type, size = Count)) +
    geom_point(aes(color = !!sym(batch_var)), position = position_dodge(width = 0.3), alpha = 0.7) +
    geom_text(aes(label = Count, group = !!sym(batch_var)), color = "black", position = position_dodge(width = 0.3), size = 3) +
    facet_wrap(vars(!!sym(clone_var)), ncol = 6, axes = "all") +
    scale_size_continuous(range = c(5, 15)) + scale_color_brewer(palette = "Set2") +
    labs(title = "Replication Structure: Contribution by Batch (Bubble)", x = stringr::str_to_title(time_var), y = "Unit") +
    theme_publication() + theme(axis.text.x = element_text(angle = 90))
  
  print(p_rep_counts_bubble)
  ggsave(file.path(diagnostics_dir, "p_qc_rep_counts_bubble.tiff"), p_rep_counts_bubble, width = 18, height = 12, dpi = 300, compression = "lzw")
  plot_database[["p_qc_rep_counts_bubble"]] <- p_rep_counts_bubble
}
#=============================================================================#
# PART 3: PERFORM STATISTICAL MODELING (Main + QC)                        #####
#=============================================================================#
logr::log_print("\n--- Starting PART 3: Statistical Modeling (Main + QC) ---", console = TRUE)

preds_plot_dir <- file.path(output_dir, "03_model_plots")
dir.create(preds_plot_dir, showWarnings = FALSE)

all_model_outputs <- list()

for (current_resp_var in response_vars) {
  message(paste("\n--- Modeling:", current_resp_var, "---"))
  
## 3A. Main Statistical Model (Group Effects) ----
# ------------------------------------------ #
  model_outputs <- run_statistical_model(
    data_summary = modeling_data,
    config = config,
    response_variable = current_resp_var
  )
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # --- GHOST ASSASSIN & NA SCRUBBER (CLEAN EXPORT TABLES) ---
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  fixed_vars <- config$key_variables$model_fixed_effects
  if (!is.null(fixed_vars)) {
    
    # 1. Find the true combinations that physically exist in your raw data
    real_combos <- modeling_data %>%
      dplyr::select(dplyr::any_of(fixed_vars)) %>%
      dplyr::distinct() %>%
      dplyr::mutate(dplyr::across(dplyr::everything(), as.character))
    
    for (tbl_name in names(model_outputs$results_tables)) {
      tbl <- model_outputs$results_tables[[tbl_name]]
      
      # 2. Kill structural ghosts using semi_join
      shared_cols <- intersect(names(tbl), names(real_combos))
      if (length(shared_cols) > 0) {
        tbl <- tbl %>%
          dplyr::mutate(dplyr::across(dplyr::all_of(shared_cols), as.character)) %>%
          dplyr::semi_join(real_combos, by = shared_cols)
      }
      
      # 3. THE NA SCRUBBER: Drop any rows where the math failed (NA estimates)
      # This cleans up the Rate_Comparisons and any missing data combinations
      if ("Expansion_Rate" %in% names(tbl))  tbl <- tbl %>% dplyr::filter(!is.na(Expansion_Rate))
      if ("Rate_Difference" %in% names(tbl)) tbl <- tbl %>% dplyr::filter(!is.na(Rate_Difference))
      if ("Estimate" %in% names(tbl))        tbl <- tbl %>% dplyr::filter(!is.na(Estimate))
      if ("P_Value" %in% names(tbl))         tbl <- tbl %>% dplyr::filter(!is.na(P_Value))
      
      # Save the cleaned table back to the list
      model_outputs$results_tables[[tbl_name]] <- tbl
    }
    logr::log_print("   -> Cleaned model tables of ghosts and NA values.")
  }
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #

  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # --- THE BASELINE REVEALER: Rename (Intercept) ---
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  if ("Model_Coefficients" %in% names(model_outputs$results_tables)) {
    logr::log_print("   -> Renaming (Intercept) to expose the hidden baseline...")
    
    # 1. Start with Time = 0
    baseline_info <- c(paste0(config$key_variables$time_variable, "=0"))
    
    # 2. Dynamically pull the exact reference levels for your specific run
    fixed_vars <- config$key_variables$model_fixed_effects
    for (var in fixed_vars) {
      if (var %in% names(modeling_data)) {
        if (is.factor(modeling_data[[var]])) {
          ref_level <- levels(modeling_data[[var]])[1]
        } else {
          # Fallback if it's just a character column
          ref_level <- sort(unique(as.character(modeling_data[[var]])))[1] 
        }
        baseline_info <- c(baseline_info, paste0(stringr::str_to_title(var), " [", ref_level, "]"))
      }
    }
    
    # 3. Stitch it together into a beautiful string
    baseline_string <- paste(baseline_info, collapse = " | ")
    new_intercept_name <- paste0("Baseline Starting Value (", baseline_string, ")")
    
    # 4. Overwrite (Intercept) in the actual table
    model_outputs$results_tables$Model_Coefficients <- model_outputs$results_tables$Model_Coefficients %>%
      dplyr::mutate(
        Parameter = ifelse(Parameter == "(Intercept)", new_intercept_name, Parameter)
      )
  }
  
  # ---------------------------------------------------------------------------#
  # 3A-BIS. QC: CLONE HETEROGENEITY (Do Clones within a Genotype differ?) ----
  # ---------------------------------------------------------------------------#
  if (is_clone_style_analysis) {
    logr::log_print(paste("   Running Clone Heterogeneity Test for", current_resp_var))
    geno_groups <- summary_per_rep %>% distinct(!!sym(primary_var))
    clone_qc_list <- list()
    
    time_var_string <- config$key_variables$time_variable
    
    # Define variable names for robust matching
    clone_col <- cfg_vars$secondary_group_var
    
    if (is.null(clone_col) || clone_col == 'null' || !clone_col %in% colnames(summary_per_rep)) {
      logr::log_print("   Skipping clone heterogeneity: secondary_group_var not found in data.")
      clone_col <- NULL
    }
    
    if (!is.null(clone_col)) {
    for(i in 1:nrow(geno_groups)) {
      cur_geno <- geno_groups[[primary_var]][i]
      sub_data <- summary_per_rep %>% dplyr::filter(!!sym(primary_var) == cur_geno)
      
      if (length(unique(sub_data[[clone_col]])) > 1) {
        try({
          # Fit Model: Response ~ Time * Clone
          form <- as.formula(paste(current_resp_var, "~", time_var_string, "*", clone_col))
          mod <- stats::lm(form, data = sub_data, na.action = na.omit)
          
          anova_res <- stats::anova(mod)
          
          # --- ROBUST P-VALUE EXTRACTION ---
          term1 <- paste0(time_var_string, ":", clone_col)
          term2 <- paste0(clone_col, ":", time_var_string)
          
          p_val <- NA
          
          if (term1 %in% rownames(anova_res)) {
            p_val <- anova_res[term1, "Pr(>F)"]
          } else if (term2 %in% rownames(anova_res)) {
            p_val <- anova_res[term2, "Pr(>F)"]
          } else {
            p_val <- tail(na.omit(anova_res$`Pr(>F)`), 1)
          }
          
          # FIX: Create dataframe with DYNAMIC column name for Genotype
          res_row <- data.frame(
            P_Value_Clones_Diff = p_val, 
            Significance = ifelse(!is.na(p_val) & p_val < 0.05, "YES", "No")
          )
          res_row[[primary_var]] <- cur_geno # Assign 'genotype' dynamically
          
          clone_qc_list[[as.character(cur_geno)]] <- res_row
          
        }, silent = TRUE)
      }
    }
    }
    if(length(clone_qc_list) > 0) model_outputs$results_tables$QC_Clone_Differences <- bind_rows(clone_qc_list)
  }
  
  # ---------------------------------------------------------------------------#
  # 3B. QC: SLOPE HETEROGENEITY (Replicates within Clones) & OUTLIER DETECTION ----
  # ---------------------------------------------------------------------------#
  qc_results_df <- NULL
  
  if (is_clone_style_analysis) {
    logr::log_print(paste("   Running Slope Heterogeneity & Outlier Test for", current_resp_var))
    
    loop_data <- summary_per_rep %>%
      dplyr::mutate(facet_label = paste(!!sym(primary_var), !!sym(cfg_vars$secondary_group_var), sep = " | "))
    
    unique_facets <- loop_data %>% distinct(facet_label, !!sym(primary_var), !!sym(cfg_vars$secondary_group_var))
    qc_list <- list()
    
    for(i in 1:nrow(unique_facets)) {
      f_label <- unique_facets$facet_label[i]
      c_geno <- unique_facets[[primary_var]][i]
      c_clone <- unique_facets[[cfg_vars$secondary_group_var]][i]
      
      # FIX: Safe UID Creation matching Plot Labels
      subset_data <- loop_data %>%
        dplyr::filter(!!sym(primary_var) == c_geno, !!sym(cfg_vars$secondary_group_var) == c_clone)
      
      if(is_crossed_model) {
        subset_data <- subset_data %>% mutate(uid = as.factor(paste(!!sym(re_cross), !!sym(rep_var), sep="-")))
      } else {
        subset_data <- subset_data %>% mutate(uid = as.factor(!!sym(rep_var)))
      }
      
      p_val <- NA
      outlier_rep <- NA
      model_type_used <- "None"
      
      n_lines <- length(unique(subset_data$uid))
      n_rows <- nrow(subset_data)
      n_params_full <- 2 * n_lines
      
      # -- FIT MODEL --
      if (n_lines > 1 && n_rows >= n_params_full + 1) {
        
        # ATTEMPT 1: Full Interaction (Response ~ Time * UID)
        qc_mod <- tryCatch({
          stats::lm(as.formula(paste(current_resp_var, "~", time_var_string, "* uid")), 
                    data = subset_data, na.action = na.omit)
        }, error = function(e) NULL)
        
        if (!is.null(qc_mod)) {
          # Success: Extract Interaction P-value
          aov_res <- stats::anova(qc_mod)
          interaction_term_name <- paste(time_var_string, ":uid", sep="")
          
          if (interaction_term_name %in% rownames(aov_res)) {
            p_val <- aov_res[interaction_term_name, "Pr(>F)"]
          } else {
            p_val <- tail(na.omit(aov_res$`Pr(>F)`), 1)
          }
          model_type_used <- "Full Interaction"
          
          # -- OUTLIER DETECTION (Method: Most Divergent Slope) --
          if (!is.na(p_val) && p_val < 0.05) {
            # Calculate slope for EACH UID (SAFE VERSION)
            slopes <- subset_data %>%
              # 1. Remove NAs for the specific variables
              dplyr::filter(!is.na(!!sym(current_resp_var)) & !is.na(!!sym(time_var_string))) %>%
              group_by(uid) %>%
              # 2. Ensure at least 2 points exist to draw a line
              dplyr::filter(n() >= 2) %>%
              # 3. Use tryCatch to prevent isolated model crashes
              do({
                fit <- tryCatch({
                  stats::lm(as.formula(paste(current_resp_var, "~", time_var_string)), data = .)
                }, error = function(e) NULL)
                if(!is.null(fit)) broom::tidy(fit) else data.frame()
              }) %>%
              ungroup() 
            
            # 4. Safely extract driver
            if(nrow(slopes) > 0 && time_var_string %in% slopes$term) {
              slopes <- slopes %>%
                dplyr::filter(term == time_var_string) %>%
                dplyr::select(uid, slope = estimate)
              
              if(nrow(slopes) > 0) {
                mean_slope <- mean(slopes$slope, na.rm = TRUE)
                slopes$diff <- abs(slopes$slope - mean_slope)
                
                driver <- slopes %>% arrange(desc(diff)) %>% slice(1)
                outlier_rep <- as.character(driver$uid)
              }
            }
          }
        } else {
          # ATTEMPT 2: Parallel Slopes (Response ~ Time + UID)
          n_params_red <- n_lines + 1
          if (n_rows >= n_params_red + 1) {
            qc_mod_red <- tryCatch({
              stats::lm(as.formula(paste(current_resp_var, "~", time_var_string, "+ uid")), 
                        data = subset_data, na.action = na.omit)
            }, error = function(e) NULL)
            
            if (!is.null(qc_mod_red)) {
              aov_res <- stats::anova(qc_mod_red)
              uid_row <- which(grepl("uid", rownames(aov_res)))[1]
              if(!is.na(uid_row)) {
                p_val <- aov_res[uid_row, "Pr(>F)"]
                model_type_used <- "Parallel Fallback"
              }
            }
          }
        }
      } # End DF Check
      
      qc_list[[as.character(f_label)]] <- data.frame(
        Genotype = c_geno, Clone = c_clone, facet_label = f_label,
        P_Value_Het = p_val, Driver_Rep = outlier_rep, Model_Used = model_type_used
      )
    } # End Loop (Clones)
    
    qc_results_df <- dplyr::bind_rows(qc_list)
  } # End QC Block Condition
  
  # 3C. Store Results
  model_outputs$results_tables$QC_Slope_Heterogeneity <- qc_results_df
  all_model_outputs[[current_resp_var]] <- model_outputs
  
} # End Main Loop (Response Vars)

# ================================================================= #
## 3C: INFLUENCE DIAGNOSTICS (Standardized Pearson Residuals) ----
# ================================================================= #


logr::log_print("\n--- Starting PART 3C: Running Influence Diagnostics ---", console = TRUE)

# Define label columns safely
lbl_geno  <- config$key_variables$primary_group_var
lbl_clone <- config$key_variables$secondary_group_var
lbl_time  <- config$key_variables$time_variable
lbl_rep   <- config$key_variables$optional_grouping_var %||% config$key_variables$repeated_measure_var
if (is.null(lbl_rep) || lbl_rep == 'null') lbl_rep <- "rep"

all_drivers_list <- list()

for (var_name in names(all_model_outputs)) {
  
  # 1. Map the var_name to its short version
  short_resp_var <- config$response_variable_shortnames[[var_name]] %||% var_name
  
  # 2. Build the nested directory paths using the short name
  var_plot_dir <- file.path(preds_plot_dir, short_resp_var)
  model_qc_dir <- file.path(var_plot_dir, "model_qc")
  
  # 3. Create the nested directories safely
  dir.create(var_plot_dir, showWarnings = FALSE, recursive = TRUE)
  dir.create(model_qc_dir, showWarnings = FALSE, recursive = TRUE)
  
  model_outputs <- all_model_outputs[[var_name]]
  if (is.null(model_outputs) || is.null(model_outputs$model_object)) next
  
  model_obj <- model_outputs$model_object
  logr::log_print(paste("...assessing influence for:", var_name))
  
  tryCatch({
    # 1. Re-align Data to match the model's exact internal frame
    cols_in_model <- unique(c(var_name, lbl_time, lbl_geno, lbl_clone, lbl_rep)) 
    data_aligned <- modeling_data %>% tidyr::drop_na(any_of(cols_in_model))
    
    # 2. Extract Standardized Pearson Residuals 
    # (This safely accounts for unequal variance margins set by gls/lme!)
    std_resids <- as.numeric(stats::residuals(model_obj, type = "pearson"))
    
    if (length(std_resids) != nrow(data_aligned)) {
      logr::log_print(paste("   -> Data length mismatch. Skipping", var_name))
      next
    }
    
    # 3. Build diagnostic dataframe
    influence_data <- data_aligned %>%
      dplyr::mutate(
        obs_id = dplyr::row_number(),
        std_resid = std_resids,
        abs_resid = abs(std_resids),
        Meta_Geno = !!sym(lbl_geno),
        Meta_Clone = if(!is.null(lbl_clone) && lbl_clone %in% names(.)) !!sym(lbl_clone) else "N/A",
        Meta_Time = !!sym(lbl_time),
        Meta_Rep = if (!is.null(lbl_rep) && lbl_rep %in% names(.)) !!sym(lbl_rep) else NA_character_,
        ) %>%
      dplyr::mutate(
        plot_label = paste(Meta_Geno, Meta_Clone, paste0("D", Meta_Time), Meta_Rep, sep = " | ")
      )
    
    # 4. Identify Drivers (Absolute Standardized Residual > 3 is the statistical gold standard)
    residual_threshold <- 3.0
    
    drivers <- influence_data %>%
      dplyr::filter(abs_resid > residual_threshold) %>%
      dplyr::arrange(desc(abs_resid))
    
    if (nrow(drivers) > 0) {
      drivers$Driven_Metric <- var_name 
      all_drivers_list[[var_name]] <- drivers
    }
    
    # 5. PLOT
    pretty_name <- response_labels[[var_name]] %||% stringr::str_to_title(var_name)
    
    p_influence <- ggplot(influence_data, aes(x = obs_id, y = abs_resid)) +
      geom_bar(stat = "identity", width = 0.2, color = "grey60") +
      geom_point(aes(color = abs_resid > residual_threshold), size = 2) +
      geom_hline(yintercept = residual_threshold, linetype = "dashed", color = "firebrick") +
      ggrepel::geom_text_repel(
        data = head(drivers, 10), aes(label = plot_label), 
        size = 3, color = "firebrick", fontface = "bold", 
        box.padding = 0.5, max.overlaps = Inf, nudge_y = 0.2
      ) +
      scale_color_manual(values = c("FALSE" = "black", "TRUE" = "firebrick")) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.1))) +
      labs(
        title = paste("Influence Diagnostics:", pretty_name), 
        subtitle = "Points above the red line (Standardized Residual > 3) are statistical outliers driving the model.",
        x = "Observation Index", 
        y = "Absolute Standardized Residual"
      ) +
      theme_publication(base_size = 12) + theme(legend.position = "none") 
    
    ggsave(file.path(model_qc_dir, paste0("inf_resids_", var_name, ".tiff")), p_influence, width = 10, height = 6, dpi=300, compression="lzw")
    
    # 6. SAVE Single Metric CSV
    if (nrow(drivers) > 0) {
      drivers_export <- drivers %>%
        dplyr::select(any_of(c("plot_label", "plate_raw", "well_raw", "Meta_Geno", "Meta_Clone", "Meta_Time", "Meta_Rep", var_name, "std_resid", "abs_resid"))) %>%
        dplyr::rename(Pearson_Residual = std_resid, Severity = abs_resid)
      write.csv(drivers_export, file.path(model_qc_dir, paste0("drivers_", var_name, ".csv")), row.names = FALSE)
    }
    
  }, error = function(e) {
    logr::log_print(paste("   WARNING: Influence check failed for", var_name, ":", e$message))
  })
}

# --- AGGREGATE MULTI-METRIC DRIVERS ---
if (length(all_drivers_list) > 0) {
  logr::log_print("...Calculating Multi-Metric Drivers table...")
  
  all_drivers_combined <- dplyr::bind_rows(all_drivers_list)
  
  multi_driver_summary <- all_drivers_combined %>%
    dplyr::group_by(plot_label, Meta_Geno, Meta_Clone, Meta_Time, Meta_Rep, across(any_of(c("plate_raw", "well_raw")))) %>%
    dplyr::summarise(
      N_Metrics_Affected = dplyr::n_distinct(Driven_Metric),
      List_of_Metrics = paste(unique(Driven_Metric), collapse = ", "),
      Max_Severity = max(abs_resid),
      .groups = "drop"
    ) %>%
    dplyr::filter(N_Metrics_Affected > 1) %>% 
    dplyr::arrange(desc(N_Metrics_Affected), desc(Max_Severity))
  
  if (nrow(multi_driver_summary) > 0) {
    save_path <- file.path(output_dir, "MULTI_METRIC_DRIVERS_SUMMARY.csv")
    write.csv(multi_driver_summary, save_path, row.names = FALSE)
    logr::log_print(paste("   -> Found", nrow(multi_driver_summary), "samples driving multiple metrics. Saved summary CSV."))
  } else {
    logr::log_print("   -> No samples found that drive multiple metrics simultaneously.")
  }
}

#=============================================================================#
# PART 4: GENERATE MODEL-BASED PLOTS                                       ####
#=============================================================================#
logr::log_print("\n--- Starting PART 4: Generating Model-Based Plots ---", console = TRUE)
preds_plot_dir <- file.path(output_dir, "03_model_plots")
dir.create(preds_plot_dir, showWarnings = FALSE)

# --- MAP LOOKUP & PALETTE SYNC ---
label_pub <- get_lbl("publication")
label_lookup <- modeling_data %>% 
  dplyr::distinct(!!sym(config$key_variables$primary_group_var), !!sym(label_pub))

# 1. Capture the official factor levels
pub_levels <- if(is.factor(modeling_data[[label_pub]])) levels(modeling_data[[label_pub]]) else unique(as.character(modeling_data[[label_pub]]))

for (current_resp_var in response_vars) {
  
  model_outputs <- all_model_outputs[[current_resp_var]]
  if (is.null(model_outputs) || is.null(model_outputs$model_object)) {
    logr::log_print(paste("Skipping plots for", current_resp_var, "as model fitting failed."))
    next
  }
  
  logr::log_print(paste("...generating model plots for", current_resp_var))
  
  model_fit_obj <- model_outputs$model_object
  y_axis_label <- response_labels[[current_resp_var]] %||% stringr::str_to_title(current_resp_var)
  fixed_effects_to_plot <- cfg_vars$model_fixed_effects %||% cfg_vars$primary_group_var
  short_resp_var <- response_shortnames[[current_resp_var]] %||% current_resp_var
  y_scale_setting <- if (current_resp_var == "mode") "free_y" else "fixed"
  
  var_plot_dir <- file.path(preds_plot_dir, short_resp_var)
  model_qc_dir <- file.path(preds_plot_dir, short_resp_var, "model_qc")
  
  dir.create(var_plot_dir, showWarnings = FALSE)
  dir.create(model_qc_dir, showWarnings = FALSE)
  
  time_bin_width <- as.numeric(config$plot_settings$time_bin_width)  # Aggregation window for smoothing timepoints
  
  # ============================================================================= #
  # --- MODEL DIAGNOSTICS & QC PLOTS (Homoscedasticity & Normality) ---
  # ============================================================================= #
  logr::log_print(paste("Generating Model QC Plots for", current_resp_var, "..."), console = TRUE)
  
  # Create a nested folder structure specifically for QC diagnostics
  diagnostics_dir <- file.path(output_dir, "01_QC_plots")
  dir.create(diagnostics_dir, showWarnings = FALSE)
  
  # 1. Extract the raw model object from your results list for the CURRENT variable
  model_obj <- all_model_outputs[[current_resp_var]]$model_object 
  
  if (is.null(model_obj)) {
    logr::log_print(paste("Skipping QC plots for", current_resp_var, "- No model object found."))
  } else {
    
    # 2. Align the metadata with the model residuals
    used_indices <- names(resid(model_obj))
    model_data_subset <- modeling_data[used_indices, ]
    
    # 3. Build the synchronized QC dataframe
    qc_data <- data.frame(
      Fitted    = as.numeric(fitted(model_obj)),
      Residuals = as.numeric(resid(model_obj, type = "pearson")),
      Genotype  = model_data_subset[[primary_var]] # Guaranteed to match lengths!
    )
    
    # --- 4. Determine Model Type Label ---
    model_class <- class(model_obj)[1]
    is_robust <- model_class %in% c("gls", "lme", "lme4") 
    
    model_type_text <- if(is_robust) {
      "Model: Robust (Adjusted for Unequal Variance)"
    } else {
      "Model: Standard (Assumes Equal Variance)"
    }
    
    # --- 5. Plot A: Homoscedasticity (Residuals vs Fitted) ---
    p_resid <- ggplot(qc_data, aes(x = Fitted, y = Residuals, color = Genotype)) +
      geom_hline(yintercept = 0, linetype = "dashed", color = "black") +
      geom_point(alpha = 0.6, size = 2.5) +
      scale_color_manual(values = primary_palette) + 
      geom_smooth(method = "loess", se = FALSE, color = if(is_robust) "darkgreen" else "red", linewidth = 1) +
      labs(
        title = paste("Residuals vs Fitted:", current_resp_var),
        subtitle = paste("Homoscedasticity |", model_type_text),
        x = "Fitted Values", 
        y = "Pearson Residuals",
        color = stringr::str_to_title(primary_var)
      ) +
      theme_cowplot(font_size = 14) +
      theme(panel.grid.major = element_line(color = "grey95"), legend.position = "none") 
    
    # --- 6. Plot B: Normal Q-Q Plot ---
    p_qq <- ggplot(qc_data, aes(sample = Residuals, color = Genotype)) +
      stat_qq(alpha = 0.5, size = 2.5) +
      stat_qq_line(aes(sample = Residuals), color = "red", linewidth = 1, inherit.aes = FALSE) + 
      scale_color_manual(values = primary_palette) +
      labs(
        title = paste("Normal Q-Q Plot:", current_resp_var),
        subtitle = "Normality of Residuals",
        x = "Theoretical Quantiles", 
        y = "Standardized Residuals",
        color = stringr::str_to_title(primary_var)
      ) +
      theme_cowplot(font_size = 14) +
      theme(
        panel.grid.major = element_line(color = "grey95"),
        legend.position = "right" 
      )
    
    # --- 7. Combine, Print, and Save ---
    p_combined_qc <- cowplot::plot_grid(p_resid, p_qq, ncol = 2, labels = c("A", "B"), rel_widths = c(1, 1.2), align = "hv", axis = "tblr")
    
    base::print(p_combined_qc)
    
    file_name <- paste0("QC_Diagnostics_", current_resp_var, ".tiff")
    ggsave(
      filename = file.path(model_qc_dir, file_name),
      plot = p_combined_qc, 
      width = 14, 
      height = 6, 
      dpi = 300, 
      compression = "lzw"
    )
    
    logr::log_print(paste("   -> QC plots saved for", current_resp_var))
  }
  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # BASE MODEL PREDICTIONS (For Ribbons/Lines)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  combo_vars_preds <- c(primary_var)
  if (is_treatment_exp) combo_vars_preds <- c(combo_vars_preds, cfg_vars$secondary_group_var)
  
  real_combos_preds <- modeling_data %>% 
    dplyr::distinct(dplyr::across(dplyr::all_of(combo_vars_preds))) %>% 
    mutate(across(everything(), as.character))
  
  if (is_treatment_exp) {
    pred_terms <- c(paste0(time_var, "[all]"), cfg_vars$secondary_group_var, primary_var)
    if (is_crossed_model) pred_terms <- c(pred_terms, re_cross)
  } else {
    pred_terms <- c(paste0(time_var, "[all]"), primary_var)
  }
  
  # Ensure we only ask for terms that actually exist in the model
  model_vars <- all.vars(stats::formula(model_fit_obj))
  valid_pred_terms <- pred_terms[gsub("\\[.*?\\]", "", pred_terms) %in% model_vars]
  
  # NEW CODE (Calculates Marginal Means)
  model_preds_raw <- tryCatch({
    ggeffects::ggemmeans(model_fit_obj, terms = valid_pred_terms, weights = "proportional")
  }, error = function(e) {
    logr::log_print(paste("GGEFFECTS ERROR for", current_resp_var, ":", e$message))
    return(NULL)
  })
  
  if (is.null(model_preds_raw)) {
    logr::log_print(paste("...WARNING: Predictions failed for", current_resp_var, "- Skipping to next plot."))
    next
  }
  
  model_preds <- as.data.frame(model_preds_raw)
  
  # --- THE BULLETPROOF PREDICTED COLUMN FIX 1 ---
  if (!"conf.low" %in% names(model_preds) && "lower" %in% names(model_preds)) {
    model_preds$conf.low <- model_preds$lower
    model_preds$conf.high <- model_preds$upper
  }
  if (!"predicted" %in% names(model_preds)) {
    if (!is.null(model_preds_raw$predicted)) {
      model_preds$predicted <- model_preds_raw$predicted
    } else if ("fit" %in% names(model_preds)) {
      model_preds$predicted <- model_preds$fit
    } else if ("estimate" %in% names(model_preds)) {
      model_preds$predicted <- model_preds$estimate
    } else if ("panel" %in% names(model_preds)) {
      model_preds$predicted <- model_preds$panel
    } else if (current_resp_var %in% names(model_preds)) {
      model_preds$predicted <- model_preds[[current_resp_var]]
    } else if ("conf.low" %in% names(model_preds) && "conf.high" %in% names(model_preds)) {
      model_preds$predicted <- (model_preds$conf.low + model_preds$conf.high) / 2
    }
  }
  
  # Force columns to be numeric (using as.character first is a safety 
  # measure in case R secretly turned them into factors!)
  model_preds$predicted <- as.numeric(as.character(model_preds$predicted))
  model_preds$conf.low  <- as.numeric(as.character(model_preds$conf.low))
  model_preds$conf.high <- as.numeric(as.character(model_preds$conf.high))
  
  if (is_treatment_exp) {
    if("group" %in% names(model_preds)) model_preds <- model_preds %>% dplyr::rename(!!sym(cfg_vars$secondary_group_var) := group)
    if("facet" %in% names(model_preds)) model_preds <- model_preds %>% dplyr::rename(!!sym(primary_var) := facet)
  } else {
    if("group" %in% names(model_preds)) model_preds <- model_preds %>% dplyr::rename(!!sym(primary_var) := group)
  }
  
  model_preds_filtered <- model_preds %>%
    mutate(across(any_of(combo_vars_preds), as.character)) %>%
    dplyr::semi_join(real_combos_preds, by = combo_vars_preds) %>%
    left_join(label_lookup %>% mutate(across(everything(), as.character)), by = primary_var) %>%
    restore_all_factors()
  
  if (is_treatment_exp) {
    color_col_preds <- cfg_vars$secondary_group_var 
    facet_vars_preds <- list(rlang::sym(label_pub))
  } else {
    color_col_preds <- label_pub
    facet_vars_preds <- list(rlang::sym(label_pub))
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 1: Model Predictions (Raw Data + Lines) -
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # 1. Safely determine variables and create clean labels for the legends
  shape_var_preds <- if (!is.null(rep_var) && rep_var != "null") rep_var else cfg_vars$secondary_group_var
  clean_color_preds <- ifelse(color_col_preds %in% c("Genotype_Pub", "Genotype_Exp", "genotype"), "Genotype", stringr::str_to_title(color_col_preds))
  clean_shape_preds <- ifelse(shape_var_preds == "rep", "Replicate", stringr::str_to_title(shape_var_preds))
  
  plot_model_preds <- ggplot() +
    geom_point(data = modeling_data, 
               aes(x = !!sym(time_var), y = !!sym(current_resp_var), 
                   color = !!sym(if(is_treatment_exp) color_col_preds else if(has_pseudo_clones) "clone_rank" else color_var), 
                   shape = !!sym(shape_var_preds)), 
               alpha = 0.8, size = 2.5)
  
  if (!is_treatment_exp) {
    # If the custom table exists, hide this redundant dots legend
    show_clone_leg <- if(is_clone_style_analysis && !is.null(p_legend)) "none" else "legend"
    if (!is.null(custom_palette)) plot_model_preds <- plot_model_preds + scale_color_manual(values = unname(custom_palette), name = "Clone", guide = show_clone_leg)
    plot_model_preds <- plot_model_preds + ggnewscale::new_scale_color() + ggnewscale::new_scale_fill()
  }
  
  plot_model_preds <- plot_model_preds +
    geom_line(data = model_preds_filtered, aes(x = x, y = predicted, color = !!sym(color_col_preds)), linewidth = 1.2) +
    geom_ribbon(data = model_preds_filtered, aes(x = x, ymin = conf.low, ymax = conf.high, fill = !!sym(color_col_preds)), alpha = 0.2) +
    facet_wrap(facet_vars_preds, scales = y_scale_setting, axes = "all") +
    # Use our clean names here!
    labs(title = paste("Modeled Trend of", y_axis_label), 
         subtitle = "Lines show overall trend. Points are individual bio-reps.", 
         x = stringr::str_to_title(time_var), y = y_axis_label, 
         color = clean_color_preds, fill = clean_color_preds, shape = clean_shape_preds) +
    theme_publication(base_size = 14) + theme(legend.position = "right")
  
  plot_model_preds <- apply_smart_palette(plot_model_preds, color_col_preds)
  
  if (!is.null(x_breaks) && all(x_breaks != 'null')) plot_model_preds <- plot_model_preds + scale_x_continuous(breaks = as.numeric(x_breaks))
  
  # --- NEW: Attach the Clone Rank Table ---
  combined_plot_1 <- if(is_clone_style_analysis && !is.null(p_legend)) plot_model_preds + p_legend + patchwork::plot_layout(widths = c(4, 1)) else plot_model_preds
  plot_width_1 <- if(is_clone_style_analysis && !is.null(p_legend)) 14 else 10
  
  assign(paste0("p_preds_", short_resp_var), combined_plot_1)
  plot_database[[paste0("p_preds_", short_resp_var)]] <- combined_plot_1
  base::print(combined_plot_1)
  ggsave(file.path(var_plot_dir, paste0("01_preds_", short_resp_var, ".tiff")), combined_plot_1, width = if(is_treatment_exp) 12 else plot_width_1, height = 8, dpi = 300, compression = "lzw")  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 1A: Model Predictions
  # (Clonal Means - One Point per Clone per Timepoint)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  clone_grp_vars <- c(primary_var, cfg_vars$secondary_group_var, "time_binned")
  if ("clone_rank" %in% names(modeling_data)) clone_grp_vars <- c(clone_grp_vars, "clone_rank")
  
  point_data_clones <- modeling_data %>%
    dplyr::mutate(time_binned = round(!!sym(time_var) / time_bin_width) * time_bin_width) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(clone_grp_vars))) %>%
    dplyr::summarise(!!rlang::sym(current_resp_var) := mean(!!rlang::sym(current_resp_var), na.rm = TRUE), .groups = "drop") %>%
    dplyr::left_join(label_lookup %>% dplyr::mutate(across(everything(), as.character)), by = primary_var) %>% restore_all_factors()
  
  plot_model_1a <- ggplot() +
    geom_point(data = point_data_clones, 
               aes(x = time_binned, y = !!sym(current_resp_var), color = !!sym(if(has_pseudo_clones) "clone_rank" else cfg_vars$secondary_group_var)), alpha = 0.7, size = 2.5)
  
  if (!is_treatment_exp && !is.null(custom_palette)) {
    # If the custom table exists, hide this redundant dots legend
    show_clone_leg <- if(is_clone_style_analysis && !is.null(p_legend)) "none" else "legend"
    plot_model_1a <- plot_model_1a + scale_color_manual(values = unname(custom_palette), name = "Clone", guide = show_clone_leg)
  }
  
  plot_model_1a <- plot_model_1a + ggnewscale::new_scale_color() + ggnewscale::new_scale_fill() +
    geom_line(data = model_preds_filtered, aes(x = x, y = predicted, color = !!sym(color_col_preds)), linewidth = 1.2) +
    geom_ribbon(data = model_preds_filtered, aes(x = x, ymin = conf.low, ymax = conf.high, fill = !!sym(color_col_preds)), alpha = 0.2) +
    facet_wrap(facet_vars_preds, scales = y_scale_setting, axes = "all") +
    # Use our clean names here!
    labs(title = paste("Modeled Trend:", y_axis_label), 
         subtitle = "Lines: Modeled trend. Points: Individual clones (rep & batch averaged).",
         x = stringr::str_to_title(time_var), y = y_axis_label, 
         color = clean_color_preds, fill = clean_color_preds) +
    theme_publication(base_size = 14) + theme(legend.position = "right")
  
  plot_model_1a <- apply_smart_palette(plot_model_1a, color_col_preds)
  
  if (!is.null(x_breaks) && all(x_breaks != 'null')) plot_model_1a <- plot_model_1a + scale_x_continuous(breaks = as.numeric(x_breaks))
  
  # --- NEW: Attach the Clone Rank Table ---
  combined_plot_1a <- if(is_clone_style_analysis && !is.null(p_legend)) plot_model_1a + p_legend + patchwork::plot_layout(widths = c(4, 1)) else plot_model_1a
  plot_width_1a <- if(is_clone_style_analysis && !is.null(p_legend)) 14 else 10
  
  assign(paste0("p_preds_avg_", short_resp_var), combined_plot_1a)
  plot_database[[paste0("p_preds_avg_", short_resp_var)]] <- combined_plot_1a
  base::print(combined_plot_1a)
  ggsave(file.path(var_plot_dir, paste0("01a_preds_avg_", short_resp_var, ".tiff")), combined_plot_1a, width = if(is_treatment_exp) 12 else plot_width_1a, height = 8, dpi = 300, compression = "lzw")
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 1B: Model Predictions (Genotype Means - One Point per Genotype per Time)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  geno_grp_vars <- c(primary_var, "time_binned")
  if (is_treatment_exp) geno_grp_vars <- c(geno_grp_vars, cfg_vars$secondary_group_var)
  
  point_data_geno <- modeling_data %>%
    dplyr::mutate(time_binned = round(!!sym(time_var) / time_bin_width) * time_bin_width) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(geno_grp_vars))) %>%
    dplyr::summarise(!!rlang::sym(current_resp_var) := mean(!!rlang::sym(current_resp_var), na.rm = TRUE), .groups = "drop") %>%
    dplyr::left_join(label_lookup %>% dplyr::mutate(across(everything(), as.character)), by = primary_var) %>% restore_all_factors()
  
  plot_model_1b <- ggplot() +
    # FIX: Map x to time_binned
    geom_point(data = point_data_geno, aes(x = time_binned, y = !!sym(current_resp_var), color = !!sym(color_col_preds)), size = 3.5, alpha = 0.9) +
    geom_line(data = model_preds_filtered, aes(x = x, y = predicted, color = !!sym(color_col_preds)), linewidth = 1.2) +
    geom_ribbon(data = model_preds_filtered, aes(x = x, ymin = conf.low, ymax = conf.high, fill = !!sym(color_col_preds)), alpha = 0.2) +
    facet_wrap(facet_vars_preds, scales = y_scale_setting, axes = "all") +
    labs(title = paste("Modeled Trend:", y_axis_label), subtitle = "Lines: Modeled genotype trend. Points: Genotype means.", 
         x = stringr::str_to_title(time_var), 
         y = y_axis_label,
         color = stringr::str_to_title(primary_var), 
         fill = stringr::str_to_title(primary_var)) +
    theme_publication(base_size = 14) + theme(legend.position = "right")
  
  plot_model_1b <- apply_smart_palette(plot_model_1b, color_col_preds)
  if (!is.null(x_breaks) && all(x_breaks != 'null')) plot_model_1b <- plot_model_1b + scale_x_continuous(breaks = as.numeric(x_breaks))
  
  assign(paste0("p_preds_1b_geno_", short_resp_var), plot_model_1b)
  plot_database[[paste0("p_preds_1b_geno_", short_resp_var)]] <- plot_model_1b
  base::print(plot_model_1b)
  ggsave(file.path(var_plot_dir, paste0("01b_preds_geno_", short_resp_var, ".tiff")), plot_model_1b, width = if(is_treatment_exp) 12 else 10, height = 8, dpi = 300, compression = "lzw")
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 2 Setup: Overall Model Predictions Summary (No Facets)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  if (is_treatment_exp) {
    pred_terms_overall <- c(paste0(time_var, "[all]"), cfg_vars$secondary_group_var, primary_var)
    color_col_overall <- cfg_vars$secondary_group_var
  } else {
    pred_terms_overall <- c(paste0(time_var, "[all]"), primary_var)
    color_col_overall <- label_pub
  }
  
  valid_pred_terms_overall <- pred_terms_overall[gsub("\\[.*?\\]", "", pred_terms_overall) %in% model_vars]
  
  # NEW CODE (Calculates Marginal Means)
  model_preds_overall_raw <- tryCatch({
    ggeffects::ggemmeans(model_fit_obj, terms = valid_pred_terms_overall, weights = "proportional")
  }, error = function(e) NULL)
  
  if(!is.null(model_preds_overall_raw)) {
    model_preds_overall <- as.data.frame(model_preds_overall_raw)
    
    # --- THE BULLETPROOF PREDICTED COLUMN FIX 2 ---
    if (!"conf.low" %in% names(model_preds_overall) && "lower" %in% names(model_preds_overall)) {
      model_preds_overall$conf.low <- model_preds_overall$lower
      model_preds_overall$conf.high <- model_preds_overall$upper
    }
    if (!"predicted" %in% names(model_preds_overall)) {
      if (!is.null(model_preds_overall_raw$predicted)) {
        model_preds_overall$predicted <- model_preds_overall_raw$predicted
      } else if ("fit" %in% names(model_preds_overall)) {
        model_preds_overall$predicted <- model_preds_overall$fit
      } else if ("estimate" %in% names(model_preds_overall)) {
        model_preds_overall$predicted <- model_preds_overall$estimate
      } else if ("panel" %in% names(model_preds_overall)) {
        model_preds_overall$predicted <- model_preds_overall$panel
      } else if (current_resp_var %in% names(model_preds_overall)) {
        model_preds_overall$predicted <- model_preds_overall[[current_resp_var]]
      } else if ("conf.low" %in% names(model_preds_overall) && "conf.high" %in% names(model_preds_overall)) {
        model_preds_overall$predicted <- (model_preds_overall$conf.low + model_preds_overall$conf.high) / 2
      }
    }

    if (is_treatment_exp) {
      if ("group" %in% names(model_preds_overall)) model_preds_overall <- model_preds_overall %>% dplyr::rename(!!sym(cfg_vars$secondary_group_var) := group)
      if ("facet" %in% names(model_preds_overall)) model_preds_overall <- model_preds_overall %>% dplyr::rename(!!sym(primary_var) := facet)
    } else {
      if ("group" %in% names(model_preds_overall)) model_preds_overall <- model_preds_overall %>% dplyr::rename(!!sym(primary_var) := group)
    }
    
    combo_vars_overall <- c(primary_var, if(is_treatment_exp) cfg_vars$secondary_group_var else NULL)
    real_combos_overall <- modeling_data %>% 
      dplyr::distinct(dplyr::across(dplyr::all_of(combo_vars_overall))) %>% 
      dplyr::mutate(across(everything(), as.character))
    
    model_preds_overall <- model_preds_overall %>%
      dplyr::mutate(across(any_of(combo_vars_overall), as.character)) %>%
      dplyr::semi_join(real_combos_overall, by = combo_vars_overall) %>%
      dplyr::left_join(label_lookup %>% dplyr::mutate(across(everything(), as.character)), by = primary_var) %>% restore_all_factors()
    
    # 1. Start a list with your mandatory base grouping variable
    group_vars <- c(color_col_overall)
    
    # 2. Conditionally add re_cross if it exists and is valid
    if (!is.null(re_cross) && re_cross != "null" && re_cross %in% names(model_preds_overall)) {
      group_vars <- c(group_vars, re_cross)
    }
    
    # 3. Conditionally add genotype if it is present in the data
    if ("genotype" %in% names(model_preds_overall)) {
      group_vars <- c(group_vars, "genotype")
    }
    
    # 4. Build the expression dynamically
    if (length(group_vars) > 1) {
      # rlang::syms converts the text vector into a list of symbol objects
      # !!! (triple-bang) splices those symbols directly into the interaction() function
      group_expr <- rlang::expr(interaction(!!!rlang::syms(group_vars)))
    } else {
      # If it's just one variable, no interaction() needed
      group_expr <- rlang::expr(!!sym(group_vars[1]))
    }
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # PLOT 2A: Overall Summary (No Facets, Clonal Means)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    var_name <- rlang::as_name(rlang::ensym(primary_var))
    n_groups <- length(unique(point_data_clones[[var_name]]))
    
    # 1. Conditionally build the aesthetics
    if (is_treatment_exp) {
      pt_aes <- aes(x = time_binned, y = !!sym(current_resp_var), color = !!sym(color_col_overall), shape = !!sym(primary_var))
      ln_aes <- aes(x = x, y = predicted, color = !!sym(color_col_overall), linetype = !!sym(primary_var), group = !!group_expr)
    } else {
      pt_aes <- aes(x = time_binned, y = !!sym(current_resp_var), color = !!sym(color_col_overall))
      ln_aes <- aes(x = x, y = predicted, color = !!sym(color_col_overall), group = !!group_expr)
    }
    
    # 2. Build the base plot
    p_overall_2a <- ggplot() +
      geom_point(data = point_data_clones, 
                 mapping = pt_aes, 
                 alpha = 0.6, size = 2.5, position = position_jitter(width = 0.5, height = 0, seed = 123),
                 key_glyph = draw_key_shifted_point) + 
      
      geom_line(data = model_preds_overall, 
                mapping = ln_aes, 
                linewidth = 1.5,
                key_glyph = draw_key_short_line) +      
      
      geom_ribbon(data = model_preds_overall, aes(x = x, ymin = conf.low, ymax = conf.high, fill = !!sym(color_col_overall), group = !!group_expr), alpha = 0.2) +
      
      labs(title = paste("Overall Modeled Trend for", y_axis_label), 
           subtitle = "Lines: Modeled overall trend. Points: Individual clones (unfaceted).", 
           x = stringr::str_to_title(cfg_vars$time_variable), 
           y = y_axis_label, 
           color = stringr::str_to_title(primary_var), 
           fill = stringr::str_to_title(primary_var)) +
      theme_publication(base_size = 14) + 
      theme(
        legend.position = "right",
        legend.key.width = unit(1.5, "cm") 
      )
    
    # 3. Conditionally add shape/linetype labels and guide overrides
    if (is_treatment_exp) {
      p_overall_2a <- p_overall_2a + 
        labs(linetype = stringr::str_to_title(primary_var),
             shape = stringr::str_to_title(primary_var)) +
        guides(
          shape = guide_legend(override.aes = list(size = 3.5, linewidth = 1.5)),
          linetype = "legend" 
        )
    }
    
    p_overall_2a <- apply_smart_palette(p_overall_2a, color_col_overall)
    
    if (!is.null(x_breaks) && all(x_breaks != 'null')) p_overall_2a <- p_overall_2a + scale_x_continuous(breaks = as.numeric(x_breaks))
    
    assign(paste0("p_overall_2a_", short_resp_var), p_overall_2a)
    plot_database[[paste0("p_overall_2a_", short_resp_var)]] <- p_overall_2a
    base::print(p_overall_2a)
    ggsave(file.path(var_plot_dir, paste0("02a_overall_", short_resp_var, ".tiff")), p_overall_2a, width = if(is_treatment_exp) 12 else 10, height = 7, dpi = 300, compression = "lzw")
    
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # PLOT 2B: Overall Summary (No Facets, Genotype Means)
    # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
    # 1. Conditionally build the aesthetics
    if (is_treatment_exp) {
      ln_aes <- aes(x = x, y = predicted, color = !!sym(color_col_overall), linetype = !!sym(primary_var), group = !!group_expr)
    } else {
      ln_aes <- aes(x = x, y = predicted, color = !!sym(color_col_overall), group = !!group_expr)
    }
    
    p_overall_2b <- ggplot() +
      # FIX: Map x to time_binned
      geom_point(data = point_data_geno, aes(x = time_binned, y = !!sym(current_resp_var), color = !!sym(color_col_overall)), alpha = 0.9, size = 3) +
      geom_line(data = model_preds_overall, mapping = ln_aes, linewidth = 1.5) +
      geom_ribbon(data = model_preds_overall, aes(x = x, ymin = conf.low, ymax = conf.high, fill = !!sym(color_col_overall), group = !!group_expr), alpha = 0.2) +
      labs(title = paste("Overall Modeled Trend for", y_axis_label), subtitle = "Lines: Modeled overall trend. Points: Genotype means (unfaceted).", 
           x = stringr::str_to_title(time_var), 
           y = y_axis_label, 
           color = stringr::str_to_title(primary_var), 
           fill = stringr::str_to_title(primary_var)) +
      theme_publication(base_size = 14) + theme(legend.position = "right")
    
    # 3. Conditionally add shape/linetype labels and guide overrides
    if (is_treatment_exp) {
      p_overall_2a <- p_overall_2a + 
        labs(linetype = stringr::str_to_title(primary_var)) +
        guides(
          shape = guide_legend(override.aes = list(size = 3.5, linewidth = 1.5))
        )
    }
    
    p_overall_2b <- apply_smart_palette(p_overall_2b, color_col_overall)
    
    if (!is.null(x_breaks) && all(x_breaks != 'null')) p_overall_2b <- p_overall_2b + scale_x_continuous(breaks = as.numeric(x_breaks))
    
    assign(paste0("p_overall_2b_", short_resp_var), p_overall_2b)
    plot_database[[paste0("p_overall_2b_", short_resp_var)]] <- p_overall_2b
    base::print(p_overall_2b)
    ggsave(file.path(var_plot_dir, paste0("02b_overall_", short_resp_var, ".tiff")), p_overall_2b, width = if(is_treatment_exp) 12 else 10, height = 7, dpi = 300, compression = "lzw")
  }  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 3: Individual Clone Trends
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  if (is_clone_style_analysis) {
    p_clone_trends <- ggplot(modeling_data, aes(x = !!sym(time_var), y = !!sym(current_resp_var), color = clone_rank, group = !!sym(color_var))) +
      geom_jitter(data = data_per_pcr, 
                  mapping = create_plot_aes(aes(color = clone_rank, group = !!sym(color_var)), re_cross), 
                  width = 0.1, alpha = 0.2, size = 1) +
      geom_point(mapping = create_plot_aes(aes(), re_cross), alpha = 0.7, size = 2.5) +
      geom_smooth(method = "lm", se = FALSE, linewidth = 0.8) +
      facet_wrap(vars(!!sym(label_pub)), scales = y_scale_setting, axes = "all") +
      labs(title = paste("Individual Clone Trends for", y_axis_label), subtitle = "Lines are linear trends. Points are bio-reps (shaped by batch).", x = stringr::str_to_title(time_var), y = y_axis_label) +
      theme_publication(base_size = 12) + theme(legend.position = "none")
    
    p_clone_trends <- apply_smart_palette(p_clone_trends, "clone_rank")
    if (!is.null(x_breaks) && all(x_breaks != 'null')) p_clone_trends <- p_clone_trends + scale_x_continuous(breaks = as.numeric(x_breaks))
    
    combined_plot <- if(is_clone_style_analysis && !is.null(p_legend)) p_clone_trends + p_legend + patchwork::plot_layout(widths = c(3, 1)) else p_clone_trends
    plot_width <- if(is_clone_style_analysis && !is.null(p_legend)) 16 else 12
    
    assign(paste0("p_clone_trends_", short_resp_var), combined_plot)
    plot_database[[paste0("p_clone_trends_", short_resp_var)]] <- combined_plot
    base::print(combined_plot)
    ggsave(file.path(var_plot_dir, paste0("03_clone_trends_", short_resp_var, ".tiff")), combined_plot, width = plot_width, height = 9, dpi = 300, compression = "lzw")
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 3B: Clone Heterogeneity Profiles 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  secondary_var <- cfg_vars$secondary_group_var
  
  if (is_clone_style_analysis && !is.null(secondary_var) && secondary_var != 'null') {
    plot_data_clones <- summary_per_rep %>% arrange(!!sym(primary_var), clone_rank) %>% 
      dplyr::select(-dplyr::any_of(label_pub)) %>%
      dplyr::left_join(label_lookup, by = primary_var) %>% 
      restore_all_factors()
    
    qc_diff_df <- model_outputs$results_tables$QC_Clone_Differences
    driver_labels <- list(); unique_genos <- unique(plot_data_clones[[primary_var]])
    
    for(g in unique_genos) {
      sub_g <- plot_data_clones %>% dplyr::filter(!!sym(primary_var) == g)
      tryCatch({
        if(length(unique(sub_g[[secondary_var]])) > 1) {
          slopes <- sub_g %>% dplyr::filter(!is.na(!!sym(current_resp_var)) & !is.na(!!sym(time_var))) %>% group_by(!!sym(secondary_var)) %>% dplyr::filter(n() >= 2) %>% do({ fit <- tryCatch({ stats::lm(as.formula(paste(current_resp_var, "~", time_var)), data = .) }, error = function(e) NULL); if(!is.null(fit)) broom::tidy(fit) else data.frame() }) %>% ungroup()
          if(nrow(slopes) > 0 && time_var %in% slopes$term) {
            slopes <- slopes %>% dplyr::filter(term == time_var) %>% dplyr::select(Clone = !!sym(secondary_var), slope = estimate)
            if(nrow(slopes) > 0) {
              avg_slope <- mean(slopes$slope, na.rm=TRUE); driver_row <- slopes %>% mutate(diff = abs(slope - avg_slope)) %>% arrange(desc(diff)) %>% slice(1)
              driver_labels[[as.character(g)]] <- as.character(driver_row$Clone)
            } else driver_labels[[as.character(g)]] <- NA
          } else driver_labels[[as.character(g)]] <- NA
        } else driver_labels[[as.character(g)]] <- NA
      }, error = function(e) { driver_labels[[as.character(g)]] <- NA })
    }
    
    stats_annot <- data.frame(Genotype_Temp = unique_genos) %>% mutate(driver = sapply(Genotype_Temp, function(x) driver_labels[[as.character(x)]]))
    names(stats_annot)[names(stats_annot) == "Genotype_Temp"] <- primary_var
    if(!is.null(qc_diff_df)) {
      if (!primary_var %in% names(qc_diff_df) && "Genotype" %in% names(qc_diff_df)) qc_diff_df <- qc_diff_df %>% dplyr::rename(!!sym(primary_var) := Genotype)
      stats_annot <- stats_annot %>% left_join(qc_diff_df, by = primary_var) %>% mutate(label_text = case_when(is.na(P_Value_Clones_Diff) ~ "", P_Value_Clones_Diff < 0.05 ~ sprintf("p = %.3f*", P_Value_Clones_Diff), TRUE ~ sprintf("p = %.2f (ns)", P_Value_Clones_Diff)), final_label = label_text, text_color = ifelse(!is.na(P_Value_Clones_Diff) & P_Value_Clones_Diff < 0.05, "red", "black"))
    } else { stats_annot$final_label <- ""; stats_annot$text_color <- "black" }
    
    label_map <- plot_data_clones %>% distinct(!!sym(primary_var), !!sym(label_pub))
    stats_annot <- stats_annot %>% left_join(label_map, by = primary_var) %>% dplyr::select(-dplyr::any_of(label_pub)) %>% restore_all_factors()
    
    label_data_clones <- plot_data_clones %>% group_by(!!sym(primary_var), !!sym(label_pub), !!sym(color_var), !!sym(secondary_var)) %>% do({ sub_d <- .; mod <- stats::lm(as.formula(paste(current_resp_var, "~", time_var)), data = sub_d); max_t <- max(sub_d[[time_var]], na.rm = TRUE); pred_y <- predict(mod, newdata = setNames(data.frame(max_t), time_var)); data.frame(x_final = max_t, y_final = as.numeric(pred_y), clone_rank = first(sub_d$clone_rank)) }) %>% ungroup() %>% dplyr::rename(!!sym(time_var) := x_final, !!sym(current_resp_var) := y_final)
    driver_df <- data.frame(G = names(driver_labels), D = unlist(driver_labels)) %>% setNames(c(primary_var, "Driver_Name"))
    label_data_clones <- label_data_clones %>% left_join(driver_df, by = primary_var) %>% mutate(is_driver = (as.character(!!sym(secondary_var)) == as.character(Driver_Name)), plot_label = ifelse(is_driver, paste0(!!sym(secondary_var), "*"), as.character(!!sym(secondary_var))), font_face = ifelse(is_driver, "bold", "plain")) %>% restore_all_factors()
    
    p_clone_het <- ggplot(plot_data_clones, aes(x = !!sym(time_var), y = !!sym(current_resp_var))) + 
      geom_smooth(method = "lm", alpha = 0.15, se = TRUE, linewidth = 0.8, aes(group = !!sym(color_var), color = clone_rank, fill = clone_rank), fullrange = TRUE) + 
      geom_point(alpha = 0.4, size = 1, aes(color = clone_rank)) + 
      ggrepel::geom_text_repel(data = label_data_clones, aes(label = plot_label, color = clone_rank, fontface = font_face), size = 2.5, nudge_x = 10, direction = "y", segment.size = 0.2, max.overlaps = Inf, min.segment.length = 0, show.legend = FALSE) + 
      geom_text(data = stats_annot, aes(label = final_label, color = I(text_color)), x = -Inf, y = Inf, hjust = -0.1, vjust = 1.2, inherit.aes = FALSE, size = 3, fontface = "bold") + 
      facet_wrap(vars(!!sym(primary_var)), scales = y_scale_setting, axes = "all") + 
      scale_x_continuous(breaks = as.numeric(x_breaks), expand = expansion(mult = c(0.05, 0.1))) + 
      labs(title = paste("Clone Heterogeneity within Genotypes:", y_axis_label), 
           subtitle = "p: Tests if clones have significantly different slopes. (*) denotes the 'Driver' clone.", 
           x = stringr::str_to_title(time_var), y = y_axis_label) + 
      theme_publication(base_size = 12) + 
      theme(legend.position = "none", strip.text = element_text(size = rel(1.1), face = "bold"))
    
    p_clone_het <- apply_smart_palette(p_clone_het, "clone_rank")
    
    assign(paste0("p_clone_het_", short_resp_var), p_clone_het)
    plot_database[[paste0("p_clone_het_", short_resp_var)]] <- p_clone_het
    base::print(p_clone_het)
    ggsave(file.path(var_plot_dir, paste0("03b_clone_het_", short_resp_var, ".tiff")), p_clone_het, width = 14, height = 10, dpi = 300, compression = "lzw")
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 4: Individual Replicate Profiles 
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  if (is_clone_style_analysis && !is.null(secondary_var) && secondary_var != 'null') {
    plot_data <- summary_per_rep %>% mutate(facet_label = paste(!!sym(primary_var), !!sym(secondary_var), sep = " | ")) %>% 
      arrange(!!sym(primary_var), clone_rank) %>% mutate(facet_label = factor(facet_label, levels = unique(.$facet_label)))
    if(is_crossed_model) plot_data <- plot_data %>% mutate(uid = as.factor(paste(!!sym(re_cross), !!sym(rep_var), sep="-"))) else plot_data <- plot_data %>% 
        mutate(uid = as.factor(!!sym(rep_var)))
    label_data <- plot_data %>% group_by(facet_label, uid) %>% dplyr::filter(!!sym(time_var) == max(!!sym(time_var), na.rm=TRUE)) %>% ungroup()
    
    qc_df <- model_outputs$results_tables$QC_Slope_Heterogeneity
    if (!is.null(qc_df)) {
      stat_labels <- qc_df %>% mutate(base_p = case_when(is.na(P_Value_Het) ~ "", P_Value_Het < 0.001 ~ "p < 0.001***", P_Value_Het < 0.05 ~ sprintf("p = %.3f*", P_Value_Het), TRUE ~ sprintf("p = %.2f (ns)", P_Value_Het)), final_label = ifelse(!is.na(P_Value_Het) & P_Value_Het < 0.05 & !is.na(Driver_Rep), paste0(base_p, "\n(Driver: ", Driver_Rep, ")"), base_p), text_color = ifelse(!is.na(P_Value_Het) & P_Value_Het < 0.05, "red", "black"))
    } else { stat_labels <- data.frame(facet_label = unique(plot_data$facet_label), final_label="", text_color="black") }
    
    p_rep_traj <- ggplot(plot_data, aes(x = !!sym(time_var), y = !!sym(current_resp_var), group = uid, color = !!sym(rep_var), fill = !!sym(rep_var))) + 
      { if (grepl("change|_change", current_resp_var, ignore.case = TRUE))
        geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) } +
      geom_smooth(method = "lm", alpha = 0.12, se = TRUE, linewidth = 0.8) + 
      geom_point(mapping = create_plot_aes(aes(), re_cross), alpha = 0.7, size = 2) + 
      ggrepel::geom_text_repel(data=label_data, aes(label=uid), size=2.5, nudge_x=5, direction="y", segment.size=0.2, show.legend=FALSE) + 
      geom_text(data=stat_labels, aes(label=final_label, color=I(text_color), group=facet_label), x=-Inf, y=Inf, hjust=-0.1, vjust=1.2, inherit.aes=FALSE, size=2.8, fontface="bold") + 
      facet_wrap(vars(facet_label), scales="fixed", axes="all") + 
      scale_x_continuous(breaks = as.numeric(x_breaks), expand = expansion(mult = c(0.05, 0.25))) + 
      labs(title = paste("Individual Replicate Trends:", y_axis_label), 
           subtitle = "Each line is one biological replicate. Red p-value = significant slope divergence.",
           x = stringr::str_to_title(time_var),
           y = y_axis_label,
           color = "Replicate",
           fill  = "Replicate",
           shape = if (is_crossed_model) "Batch" else NULL
      ) +
      { if(is_crossed_model) labs(shape = "Batch") } + 
      theme_publication(base_size = 12) + 
      theme(legend.position = "none", strip.text = element_text(size = rel(0.7), face="bold"))
    
    assign(paste0("p_rep_profiles_", short_resp_var), p_rep_traj)
    plot_database[[paste0("p_rep_profiles_", short_resp_var)]] <- p_rep_traj
    base::print(p_rep_traj)
    ggsave(file.path(var_plot_dir, paste0("04_rep_profiles_", short_resp_var, ".tiff")), p_rep_traj, width=16, height=12, dpi=300, compression="lzw")
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 5: Slope Coefficients (Faceted & Averaged Versions)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  logr::log_print(paste("...generating slope plots for", current_resp_var), console=TRUE)
  
  slope_col_name <- "Expansion_Rate"
  
  # --- NEW: Safe Batch Variable Logic ---
  # 1. Start with the crossed effect (batch)
  safe_batch_var <- re_cross 
  
  # 2. If null, fall back to the optional grouping var (as requested)
  if (is.null(safe_batch_var) || safe_batch_var == "null") {
    safe_batch_var <- cfg_vars$optional_grouping_var
  }
  
  # 3. Final safety net (if both are null in config)
  if (is.null(safe_batch_var) || safe_batch_var == "null") {
    safe_batch_var <- NULL
  }
  
  # --- 1. Define Variables based on Experiment Style ---
  if (is_treatment_exp) {
    plot_x_var     <- cfg_vars$secondary_group_var # Treatment
    plot_fill_var  <- label_pub
    plot_facet_var <- label_pub                    # Facet by Genotype
    batch_var      <- safe_batch_var               # Safe Fallback Variable
    color_label    <- cfg_vars$secondary_group_var
  } else {
    plot_x_var     <- label_pub                    # Genotype
    plot_fill_var  <- label_pub
    plot_facet_var <- NULL
    batch_var      <- safe_batch_var
    color_label    <- primary_var
  }
  
  # --- 2. Prepare Data Layers ---
  
  # A. Individual Slopes (The Dots)
  dot_data <- data.frame()
  if (!is.null(model_outputs$results_tables$Individual_Slopes_BLUPs)) {
    dot_data <- model_outputs$results_tables$Individual_Slopes_BLUPs %>%
      dplyr::mutate(!!sym(primary_var) := as.character(!!sym(primary_var))) %>%
      dplyr::left_join(label_lookup %>% mutate(across(everything(), as.character)), by = primary_var) %>%
      restore_all_factors()
  }
  
  # B. Overall Data (The Bars)
  # any_of() ensures it doesn't crash if batch_var was marginalized out of the model
  bar_data_master <- model_outputs$results_tables$Group_Expansion_Rates %>%
    mutate(across(any_of(c(primary_var, batch_var, plot_x_var)), as.character)) %>%
    dplyr::left_join(label_lookup %>% mutate(across(everything(), as.character)), by = primary_var) %>%
    restore_all_factors()
  
  # Apply factor levels from config, translating through value_renaming so
  # renamed values (e.g. "FAN1 KO") match the level strings.
  bar_data_master <- apply_renaming_and_factors(bar_data_master, config)
  
  # --- 3. PLOT 5A: FACETED (Batch-Specific) ---
  logr::log_print("......creating faceted batch-specific plot")
  
  # CRITICAL FIX: Ensure batch_var actually exists in the data before attempting to facet by it!
  valid_batch_facet <- if (!is.null(batch_var) && batch_var %in% names(bar_data_master)) batch_var else NULL
  
  # Build the facet formula dynamically based on what survives the checks
  if (!is.null(plot_facet_var) && !is.null(valid_batch_facet)) {
    facet_formula <- as.formula(paste("~", plot_facet_var, "+", valid_batch_facet))
  } else if (!is.null(plot_facet_var)) {
    facet_formula <- as.formula(paste("~", plot_facet_var))
  } else if (!is.null(valid_batch_facet)) {
    facet_formula <- as.formula(paste("~", valid_batch_facet))
  } else {
    facet_formula <- NULL
  }
  
  p_slopes_faceted <- ggplot(bar_data_master, aes(x = !!sym(plot_x_var), y = !!sym(slope_col_name), fill = !!sym(plot_fill_var))) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_col(color = "black", alpha = 0.8, position = position_dodge(width = 0.9)) +
    geom_errorbar(aes(ymin = CI_Lower_95, ymax = CI_Upper_95), width = 0.3, position = position_dodge(width = 0.9)) +
    labs(title = paste("Faceted Expansion Rates:", y_axis_label),
         subtitle = "Showing modeled individual rates.",
         y = paste0("Rate (Slope per ", stringr::str_to_title(config$key_variables$time_variable), ")"),
         x = stringr::str_to_title(color_var)) +
    theme_publication(base_size = 14) + theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")
  
  # Apply facets safely
  if (!is.null(facet_formula)) {
    p_slopes_faceted <- p_slopes_faceted + facet_wrap(facet_formula, axes = "all")
  }
  
  # --- 4. PLOT 5B: AVERAGED (Aggregated across Batches) ---
  logr::log_print("......creating averaged (aggregated) plot")
  
  bar_data_avg <- bar_data_master %>%
    dplyr::group_by(dplyr::across(dplyr::any_of(c(plot_facet_var, plot_x_var, plot_fill_var)))) %>%
    dplyr::summarise(
      Mean_Rate = mean(!!sym(slope_col_name), na.rm = TRUE),
      Lower = mean(CI_Lower_95, na.rm = TRUE),
      Upper = mean(CI_Upper_95, na.rm = TRUE),
      .groups = "drop"
    ) %>%
    # FIX 1: Hard-convert X to a factor in the data so ggplot doesn't get confused
    dplyr::mutate(!!sym(plot_x_var) := as.factor(!!sym(plot_x_var)))
  
  # Also hard-convert the dot data if it exists
  if(nrow(dot_data) > 0) {
    dot_data <- dot_data %>% dplyr::mutate(!!sym(plot_x_var) := as.factor(!!sym(plot_x_var)))
  }
  
  # --- NEW: Extract ALL LMER Significance Brackets ---
  bracket_data <- data.frame()
  
  if (!is.null(model_outputs$results_tables$Rate_Comparisons)) {
    comps <- model_outputs$results_tables$Rate_Comparisons
    
    # UPDATE 1: Check for the new column name
    if ("p_value_adj" %in% names(comps)) {
      
      bracket_data <- comps %>%
        # 1. Split the 'contrast' column (e.g., "WT - S126D") into two targets
        dplyr::mutate(
          group1 = trimws(sapply(strsplit(as.character(contrast), " - "), `[`, 1)),
          group2 = trimws(sapply(strsplit(as.character(contrast), " - "), `[`, 2))
        ) %>%
        # UPDATE 2: Rename from the new column name
        dplyr::rename(p.adj = p_value_adj) %>%
        # 2. Map the adjusted p-value to stars
        dplyr::mutate(
          stars = dplyr::case_when(
            p.adj <= 0.001 ~ "***",
            p.adj <= 0.01  ~ "**",
            p.adj <= 0.05  ~ "*",
            TRUE ~ "ns"
          )
        ) %>%
        # 3. Filter out "ns" to prevent drawing a messy web of non-significant brackets
        dplyr::filter(stars != "ns")
      
      # 4. Safety Check: Translate original names to the renamed labels on the plot X-axis
      if (exists("label_lookup") && nrow(bracket_data) > 0) {
        translate_dict <- setNames(as.character(label_lookup[[label_pub]]), as.character(label_lookup[[primary_var]]))
        
        bracket_data <- bracket_data %>%
          dplyr::mutate(
            group1 = ifelse(!is.na(translate_dict[group1]), translate_dict[group1], group1),
            group2 = ifelse(!is.na(translate_dict[group2]), translate_dict[group2], group2)
          )
      }
      
      # 5. Calculate staggering Y positions so the brackets stack cleanly and don't overlap
      if(nrow(bracket_data) > 0) {
        y_max_limit <- max(bar_data_avg$Upper, na.rm = TRUE)
        step_increase <- y_max_limit * 0.15 # Spacing between stacked brackets
        
        bracket_data <- bracket_data %>%
          dplyr::mutate(y.position = y_max_limit + (dplyr::row_number() * step_increase))
      }
    }
  }
  
  # --- PLOT 5B: AVERAGED (Aggregated across Batches) ---
  
  p_slopes_avg <- ggplot(bar_data_avg, aes(
    x = !!sym(plot_x_var), 
    y = Mean_Rate, 
    fill = !!sym(plot_fill_var),
    group = !!sym(plot_fill_var) 
  )) +
    geom_hline(yintercept = 0, color = "black", linewidth = 0.5) +
    geom_col(width = 0.75, color = "black", alpha = 0.8, position = position_dodge(width = 0.85))
  
  if(nrow(dot_data) > 0) {
    p_slopes_avg <- p_slopes_avg + 
      geom_point(data = dot_data, aes(y = Individual_Slope, x = !!sym(plot_x_var), fill = !!sym(plot_fill_var)), 
                 inherit.aes = FALSE, shape = 21, color = "black", alpha = 0.7, size = 2,
                 position = position_jitterdodge(jitter.width = 0.2, dodge.width = 0.85))
  }
  
  p_slopes_avg <- p_slopes_avg +
    geom_errorbar(aes(ymin = Lower, ymax = Upper), width = 0.25, position = position_dodge(width = 0.85))
  
  # --- NEW: ADD SIGNIFICANCE BRACKETS ---
  if (nrow(bracket_data) > 0) {
    p_slopes_avg <- p_slopes_avg +
      ggpubr::stat_pvalue_manual(
        data = bracket_data, 
        label = "stars",
        y.position = "y.position",
        bracket.size = 0.6,
        label.size = 6,
        tip.length = 0.02
      )
  }
  # --------------------------------------
  
  p_slopes_avg <- p_slopes_avg +
    labs(title = paste("Averaged Expansion Rates:", y_axis_label),
         subtitle = "Aggregated group rates. Error bars show 95% CIs. Brackets show Tukey-adjusted LMER significance.",
         y = paste0("Mean Rate (Slope per ", stringr::str_to_title(config$key_variables$time_variable), ")"),
         x = stringr::str_to_title(color_label)) +
    theme_publication(base_size = 14) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), axis.title.x = element_blank(), legend.position = "none")  # --- 5. Apply Palettes and Save ---
  plot_list <- list(faceted = p_slopes_faceted, averaged = p_slopes_avg)
  
  for(type in names(plot_list)) {
    p <- plot_list[[type]]
    p <- apply_smart_palette(p, plot_fill_var)
    
    file_name <- paste0("p_slopes_", type, "_", short_resp_var)
    plot_database[[file_name]] <- p
    assign(paste0(file_name), p_rep_traj)
    base::print(p)
    ggsave(file.path(var_plot_dir, paste0("05_", file_name, ".tiff")), p, width = 12, height = 8, dpi = 300, compression = "lzw")
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 6: (QC) Slope Calculation Comparison
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  if (is_treatment_exp) {
    qc_x_var <- cfg_vars$secondary_group_var 
    qc_facet_vars <- c(primary_var)
    if (is_crossed_model) qc_facet_vars <- unique(c(qc_facet_vars, re_cross))
  } else {
    qc_x_var <- tail(fixed_effects_to_plot, 1); qc_facet_vars <- head(fixed_effects_to_plot, -1) 
  }
  fixed_effects_to_compare <- unique(c(qc_x_var, qc_facet_vars))
  
  slopes_from_model <- all_model_outputs[[current_resp_var]]$results_tables$Group_Expansion_Rates %>% dplyr::select(dplyr::all_of(fixed_effects_to_compare), slope = !!rlang::sym(slope_col_name)) %>% dplyr::mutate(method = "Formal Model (fits all data)")
  avg_modeling_data <- modeling_data %>% dplyr::group_by(dplyr::across(dplyr::all_of(c(fixed_effects_to_compare, time_var)))) %>% dplyr::summarise(!!rlang::sym(current_resp_var) := mean(!!rlang::sym(current_resp_var), na.rm = TRUE), .groups = "drop")
  slopes_from_trendline <- avg_modeling_data %>% dplyr::group_by(dplyr::across(dplyr::all_of(fixed_effects_to_compare))) %>% dplyr::filter(n() > 1) %>% dplyr::do(broom::tidy(stats::lm(as.formula(paste(current_resp_var, "~", time_var)), data = .))) %>% dplyr::ungroup() %>% dplyr::filter(term == time_var) %>% dplyr::select(dplyr::all_of(fixed_effects_to_compare), slope = estimate) %>% dplyr::mutate(method = "Simple Trendline (fits averages)")
  comparison_data <- dplyr::bind_rows(slopes_from_model, slopes_from_trendline)
  combo_vars_qc <- fixed_effects_to_compare
  real_combos_qc <- modeling_data %>% dplyr::distinct(dplyr::across(dplyr::all_of(combo_vars_qc))) %>% dplyr::mutate(across(everything(), as.character))
  
  comparison_data <- comparison_data %>% 
    dplyr::mutate(across(any_of(combo_vars_qc), as.character)) %>% 
    dplyr::semi_join(real_combos_qc, by = combo_vars_qc) %>% 
    dplyr::left_join(label_lookup %>% mutate(across(everything(), as.character)), by = primary_var) %>% 
    restore_all_factors()
  
  # --- EXTRACT FORMULAS FOR PRINTING ---
  
  # 1. Simple Trendline Formula (Grouping happens before this is applied)
  simple_formula_str <- paste(current_resp_var, "~", time_var)
  
  # 2. Formal Model Formula (Extract from the saved model object)
  # Replace 'model_fit' with the actual name of the saved model object in your list
  formal_model_obj <- all_model_outputs[[current_resp_var]]$model_object
  
  if (!is.null(formal_model_obj)) {
    # Extract the formula and collapse it into a single clean string
    formal_formula_str <- paste(trimws(deparse(formula(formal_model_obj))), collapse = " ")
  } else {
    formal_formula_str <- "[Formal formula not found in model_outputs]"
  }
  
  plot_qc_x_var <- qc_x_var
  plot_qc_facet_vars <- qc_facet_vars
  if (!is_treatment_exp && qc_x_var == primary_var) plot_qc_x_var <- label_pub
  plot_qc_facet_vars[plot_qc_facet_vars == primary_var] <- label_pub
  
  p_compare <- ggplot(comparison_data, aes(x = !!sym(plot_qc_x_var), y = slope, fill = method)) + 
    geom_col(position = position_dodge(width = 0.9), color = "black", alpha = 0.8) + scale_fill_brewer(palette = "Set1") + 
    labs(title = paste("Calculation QC for", y_axis_label), 
         subtitle = "Compares the complex 'Formal Model' (red) to a 'Simple Trendline' of the group averages (blue).", 
         x = stringr::str_to_title(qc_x_var), 
         y = paste0("Rate of Change (Slope per ", stringr::str_to_title(stringr::str_remove(time_var, "s$")), ")"), 
         fill = "Calculation Method",
         caption = paste0("Formal Model (Red): ", formal_formula_str, 
                          "\nSimple Trendline (Blue): ", simple_formula_str, " (applied to group averages)")) + 
    theme_publication(base_size = 14) + 
    theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "bottom", 
          plot.caption = element_text(hjust = 0, size = 10, color = "grey30", margin = margin(t = 15)))
  
  if (length(plot_qc_facet_vars) > 0) p_compare <- p_compare + facet_wrap(vars(!!!lapply(plot_qc_facet_vars, rlang::sym)), axes = "all")
  base::print(p_compare)
  ggsave(file.path(var_plot_dir, paste0("p_qc_slope_compare_", short_resp_var, ".tiff")), p_compare, width = if(is_treatment_exp) 12 else 11, height = 8, dpi = 300, compression = "lzw")
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 6B: (QC) Visualizing the Shift: Naive Lines vs. Formal Model Line
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  tryCatch({
    if(is_crossed_model) vis_data <- modeling_data %>% dplyr::mutate(uid = as.factor(paste(!!sym(primary_var), !!sym(re_cross), !!sym(rep_var), sep="_"))) else vis_data <- modeling_data %>% dplyr::mutate(uid = as.factor(paste(!!sym(primary_var), !!sym(rep_var), sep="_")))
    
    pred_terms_qc <- c(paste0(time_var, "[all]"), primary_var)
    valid_pred_terms_qc <- pred_terms_qc[gsub("\\[.*?\\]", "", pred_terms_qc) %in% model_vars]
    
    # NEW CODE (Calculates Marginal Means)
    formal_preds_raw <- tryCatch({
      ggeffects::ggemmeans(model_fit_obj, terms = valid_pred_terms_qc, weights = "proportional")
    }, error = function(e) NULL)
    
    if(!is.null(formal_preds_raw)) {
      formal_preds <- as.data.frame(formal_preds_raw)
      
      # --- THE BULLETPROOF PREDICTED COLUMN FIX 3 ---
      if (!"conf.low" %in% names(formal_preds) && "lower" %in% names(formal_preds)) {
        formal_preds$conf.low <- formal_preds$lower; formal_preds$conf.high <- formal_preds$upper
      }
      if (!"predicted" %in% names(formal_preds)) {
        if (!is.null(formal_preds_raw$predicted)) { formal_preds$predicted <- formal_preds_raw$predicted
        } else if ("fit" %in% names(formal_preds)) { formal_preds$predicted <- formal_preds$fit
        } else if ("estimate" %in% names(formal_preds)) { formal_preds$predicted <- formal_preds$estimate
        } else if ("panel" %in% names(formal_preds)) { formal_preds$predicted <- formal_preds$panel
        } else if (current_resp_var %in% names(formal_preds)) { formal_preds$predicted <- formal_preds[[current_resp_var]]
        } else if ("conf.low" %in% names(formal_preds) && "conf.high" %in% names(formal_preds)) {
          formal_preds$predicted <- (formal_preds$conf.low + formal_preds$conf.high) / 2
        }
      }

      formal_preds <- formal_preds %>% dplyr::rename(!!sym(primary_var) := group) %>% dplyr::left_join(label_lookup %>% mutate(across(everything(), as.character)), by = primary_var) %>% restore_all_factors()
      p_vis_shift <- ggplot() + 
        geom_smooth( data = vis_data, aes(
          x = !!sym(time_var),
          y = !!sym(current_resp_var),
          group = uid
        ), method = "lm", se = FALSE, color = "grey60", linewidth = 0.4, alpha = 0.6) + 
        geom_point(data = vis_data, aes(
          x = !!sym(time_var),
          y = !!sym(current_resp_var)
        ), color = "grey80", size = 1.5, alpha = 0.5) + 
        geom_line(data = formal_preds,
                    aes(x = x, y = predicted, color = !!sym(label_pub)), linewidth = 2) + 
        facet_wrap(vars(!!sym(label_pub)), scales = y_scale_setting) + 
        labs(
          title = paste("Understanding the Model Adjustment:", y_axis_label),
          subtitle = "Thin lines = Naive Fits. Thick colored line = Adjusted Model Fit.",
          x = stringr::str_to_title(time_var),
          y = y_axis_label,
          color = "Formal Model Estimate"
        ) + 
        theme_publication(base_size = 14) + theme(legend.position = "bottom")
      
      p_vis_shift <- apply_smart_palette(p_vis_shift, label_pub)
      base::print(p_vis_shift)
      ggsave(file.path(var_plot_dir, paste0("p_qc_vis_shift_", short_resp_var, ".tiff")), p_vis_shift, width = 10, height = 6, dpi = 300, compression = "lzw")
    }
  }, error = function(e) logr::log_print(paste("...WARNING: Could not calculate visualization shift plot. Error:", e$message)))
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 7: (QC) Variance Components / Effect Size (LMER or LM)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  p_var_comp <- NULL
  if (inherits(model_fit_obj, "lmerMod")) {
    if ("Random_Effects_Variance" %in% names(model_outputs$results_tables)) {
      var_comp_data <- model_outputs$results_tables$Random_Effects_Variance 
      
      p_var_comp <- ggplot(var_comp_data, aes(x = "Variance Explained", y = ICC_Percent, fill = Source)) + 
        geom_col(color = "white", linewidth = 0.4, width = 0.55) +
        geom_text(aes(label = paste0(round(ICC_Percent, 1), "%")), position = position_stack(vjust = 0.5), size = 3.5, fontface = "bold") + 
        scale_y_continuous(labels = scales::percent_format(scale = 1), expand = expansion(mult = c(0, 0.03))) + 
        scale_fill_brewer(palette = "Set3") +
        coord_flip() + 
        labs(title = paste("Sources of Variance:", y_axis_label), 
             subtitle = "From LMER model. Shows % of total variance explained by each random effect.", 
             x = NULL, 
             y = "Percent of Total Variance (%)", 
             fill = "Variance Source") + 
        theme_publication(base_size = 12) + 
        theme(axis.text.y = element_blank(), axis.ticks.y = element_blank(), panel.grid.major.y = element_blank(), legend.position = "right")
    }
  } else if (inherits(model_fit_obj, "lm")) {
    aov_res <- stats::anova(model_fit_obj)
    ss_total <- sum(aov_res$`Sum Sq`, na.rm = TRUE)
    var_data <- data.frame(Source = rownames(aov_res), SS = aov_res$`Sum Sq`) %>% 
      dplyr::mutate(Pct_Explained = (SS / ss_total) * 100, Source = stringr::str_replace_all(Source, ":", " × ")) %>% 
      dplyr::arrange(Pct_Explained) %>% 
      dplyr::mutate(Source = factor(Source, levels = Source),         
                    # Colour by whether this is an interaction term
                    Term_Type = ifelse(grepl("×", Source), "Interaction", "Main Effect")) 
    
    p_var_comp <- ggplot(var_data,
                         aes(x = Source, y = Pct_Explained, fill = Term_Type)) +
      geom_col(color = "black", alpha = 0.85) +
      geom_text(aes(label = sprintf("%.1f%%", Pct_Explained)),
                hjust = -0.15, size = 3.5) +
      scale_y_continuous(expand = expansion(mult = c(0, 0.2))) +
      scale_fill_manual(values = c("Main Effect" = "#0072B2", "Interaction" = "#CC79A7")) +
      coord_flip() +
      labs(
        title    = paste("Variance Explained by Fixed Effects:", y_axis_label),
        subtitle = "ANOVA eta-squared. Percentage of total variance explained by each model term.",
        x = "Model Term", y = "% of Total Variance Explained",
        fill = "Term Type"
      ) +
      theme_publication(base_size = 12) +
      theme(legend.position = "right")
  }
  if (!is.null(p_var_comp)) { base::print(p_var_comp); ggsave(file.path(model_qc_dir, paste0("p_qc_varcomp_", short_resp_var, ".tiff")), p_var_comp, width = 9, height = 6, dpi = 300, compression = "lzw") }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 8: Automated Slope Fold Change (Ratio) vs. Control
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  logr::log_print(paste("...generating Fold Change plot(s) for", current_resp_var), console=TRUE)
  
  slopes_data <- model_outputs$results_tables$Group_Expansion_Rates
  
  # Safely pull control variables from the config (with fallbacks just in case)
  fc_cfg <- config$parameters$fold_change_controls
  ctrl_geno <- ifelse(!is.null(fc_cfg$baseline_genotype), fc_cfg$baseline_genotype, "WT")
  ctrl_dose <- ifelse(!is.null(fc_cfg$baseline_dose), as.character(fc_cfg$baseline_dose), "0")

  if (!is.null(slopes_data) && nrow(slopes_data) > 0) {
    
    # --- 1. Determine Plotting Strategy Based on Experiment Type ---
    if (is_treatment_exp) {
      plot_x_var    <- cfg_vars$secondary_group_var 
      plot_fill_var <- label_pub                    
      grouping_cols <- c(primary_var, plot_x_var)
      methods_to_run <- c("genotype_matched", "global_absolute")
      color_label <- cfg_vars$secondary_group_var
    } else {
      plot_x_var    <- label_pub
      plot_fill_var <- label_pub
      grouping_cols <- c(primary_var)
      methods_to_run <- c("clone_baseline")
      color_label <- primary_var
    }
    
    # Clean the incoming data
    combo_vars_fc <- c(primary_var, if(is_treatment_exp) cfg_vars$secondary_group_var else NULL)
    overall_slopes <- slopes_data %>% 
      dplyr::mutate(dplyr::across(dplyr::any_of(combo_vars_fc), as.character)) %>% 
      dplyr::semi_join(
        modeling_data %>% 
          dplyr::distinct(across(all_of(combo_vars_fc))) %>% 
          dplyr::mutate(dplyr::across(dplyr::everything(), as.character)), 
        by = combo_vars_fc
      )
    
    # --- 2. Loop Through and Generate the Required Plots ---
    for (method in methods_to_run) {
      tryCatch({
        
        # Calculate Math based on the specific method
        if (method == "clone_baseline") {
          logr::log_print(paste("......calculating clone baseline vs", ctrl_geno))
          
          baseline_val <- overall_slopes %>% 
            dplyr::filter(!!sym(primary_var) == ctrl_geno) %>% 
            dplyr::summarise(Base = mean(Expansion_Rate, na.rm = TRUE)) %>% dplyr::pull(Base)
          
          fc_math <- overall_slopes %>% 
            dplyr::group_by(dplyr::across(dplyr::any_of(grouping_cols))) %>% 
            dplyr::summarise(Mean = mean(Expansion_Rate, na.rm=TRUE), L = mean(CI_Lower_95, na.rm=TRUE), U = mean(CI_Upper_95, na.rm=TRUE), .groups="drop") %>%
            dplyr::mutate(ratio = Mean/baseline_val[1], CI_Lower_95 = L/baseline_val[1], CI_Upper_95 = U/baseline_val[1])
          
          sub_text <- paste("Ratio relative to the", ctrl_geno, "average.")
          file_suffix <- "vs_control_clone"
          
        } else if (method == "global_absolute") {
          logr::log_print(paste("......calculating global baseline vs", ctrl_geno, "at dose", ctrl_dose))
          
          baseline_val <- overall_slopes %>% 
            dplyr::filter(!!sym(primary_var) == ctrl_geno & !!sym(plot_x_var) == ctrl_dose) %>% 
            dplyr::summarise(Base = mean(Expansion_Rate, na.rm = TRUE)) %>% dplyr::pull(Base)
          
          fc_math <- overall_slopes %>% 
            dplyr::group_by(dplyr::across(dplyr::any_of(grouping_cols))) %>% 
            dplyr::summarise(Mean = mean(Expansion_Rate, na.rm=TRUE), L = mean(CI_Lower_95, na.rm=TRUE), U = mean(CI_Upper_95, na.rm=TRUE), .groups="drop") %>%
            dplyr::mutate(ratio = Mean/baseline_val[1], CI_Lower_95 = L/baseline_val[1], CI_Upper_95 = U/baseline_val[1])
          
          sub_text <- paste("Ratio relative to the Untreated", ctrl_geno, "baseline.")
          file_suffix <- "global_baseline"
          
        } else if (method == "genotype_matched") {
          logr::log_print(paste("......calculating genotype-matched baseline at dose", ctrl_dose))
          
          baseline_lookup <- overall_slopes %>% 
            dplyr::filter(!!sym(plot_x_var) == ctrl_dose) %>% 
            dplyr::group_by(!!sym(primary_var)) %>% 
            dplyr::summarise(Base = mean(Expansion_Rate, na.rm = TRUE), .groups = "drop")
          
          fc_math <- overall_slopes %>% 
            dplyr::group_by(dplyr::across(dplyr::any_of(grouping_cols))) %>% 
            dplyr::summarise(Mean = mean(Expansion_Rate, na.rm=TRUE), L = mean(CI_Lower_95, na.rm=TRUE), U = mean(CI_Upper_95, na.rm=TRUE), .groups="drop") %>%
            dplyr::left_join(baseline_lookup, by = primary_var) %>%
            dplyr::mutate(ratio = Mean/Base, CI_Lower_95 = L/Base, CI_Upper_95 = U/Base)
          
          sub_text <- "Ratio relative to the Untreated control within each genotype."
          file_suffix <- "genotype_matched"
        }
        
        # Apply labels and ordered factors
        fold_change_data <- fc_math %>%
          dplyr::left_join(label_lookup %>% mutate(across(everything(), as.character)), by = primary_var) %>% 
          restore_all_factors() %>%
          # FIX 1: Hard-convert the X-axis to a factor inside the dataset
          dplyr::mutate(!!sym(plot_x_var) := as.factor(!!sym(plot_x_var)))
        
        # Apply factor levels, translating through value_renaming so renamed
        # values match the level strings and bars are not silently dropped.
        fold_change_data <- apply_renaming_and_factors(fold_change_data, config)
        
        # --- 3. Build and Save the Plot ---
        # FIX 2: Added explicit 'group' aesthetic 
        p_fold_change <- ggplot(fold_change_data, aes(
          x = !!sym(plot_x_var), 
          y = ratio, 
          fill = !!sym(plot_fill_var),
          group = !!sym(plot_fill_var)
        )) + 
          geom_hline(yintercept = 1.0, linetype = "dashed", color = "blue", linewidth = 1) + 
          
          # FIX 3: Explicitly bound the bar width (0.7) to match the dodge width (0.7)
          geom_col(width = 0.8, color = "black", alpha = 0.8, position = position_dodge(width = 0.85)) + 
          geom_errorbar(aes(ymin = CI_Lower_95, ymax = CI_Upper_95), width = 0.25, position = position_dodge(width = 0.85)) + 
          
          scale_y_continuous(expand = expansion(mult = c(0, 0.05))) + 
          labs(
            title = paste("Relative Rate of Change in", y_axis_label), 
            subtitle = sub_text, 
            y = "Rate Fold Change (Ratio)", 
            x = stringr::str_to_title(color_label)
          ) + 
          theme_publication(base_size = 14) + 
          theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
        
        p_fold_change <- apply_smart_palette(p_fold_change, plot_fill_var)
        
        # Save dynamically named file
        base::print(p_fold_change)
        file_name <- paste0("p_slope_FC_", file_suffix, "_", short_resp_var, ".tiff")
        ggsave(file.path(var_plot_dir, file_name), p_fold_change, width = 10, height = 7, dpi = 300, compression = "lzw")    
        
      }, error = function(e) {
        logr::log_print(paste("Error generating fold change plot (", method, "):", e$message))
      })
    }
  }  
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 9: Coefficients Plot (Forest Plot of Effects)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  coef_data <- model_outputs$results_tables$Model_Coefficients
  if (!is.null(coef_data) && nrow(coef_data) > 0) {
    plot_coefs <- coef_data %>%
      dplyr::filter(!grepl("Baseline", Parameter)) %>%
      dplyr::mutate(
        P_Num = suppressWarnings(as.numeric(gsub("< ", "", P_Value))),
        Significance = dplyr::case_when(
          P_Num < 0.001 ~ "p < 0.001",
          P_Num < 0.05  ~ "p < 0.05",
          TRUE          ~ "n.s."
        ),
        # Human-readable labels using the time_var dynamically
        Clean_Label = Parameter %>%
          stringr::str_replace_all(paste0(time_var, " × "), "Speed: ") %>%
          stringr::str_replace_all(paste0("^", time_var, "$"), "Baseline Expansion Speed") %>%
          stringr::str_trim(),
        Clean_Label = stats::reorder(Clean_Label, Estimate),
        # CIs
        CI_Lower_95 = suppressWarnings(as.numeric(CI_Lower_95)),
        CI_Upper_95 = suppressWarnings(as.numeric(CI_Upper_95))
      ) %>%
      dplyr::filter(!is.na(Estimate), !is.na(CI_Lower_95))
    
    p_coef <- ggplot(plot_coefs,
                     aes(x = Estimate, y = Clean_Label, color = Significance)) +
      geom_vline(xintercept = 0, linetype = "dashed", color = "grey40", linewidth = 0.6) +
      annotate("text", x = 0, y = -Inf, label = "No effect", vjust = -0.4,
               hjust = 0.5, size = 3, color = "grey40") +
      geom_errorbarh(
        aes(xmin = CI_Lower_95, xmax = CI_Upper_95),
        height = 0.25, linewidth = 0.9
      ) +
      geom_point(size = 4) +
      scale_color_manual(
        values = c("p < 0.001" = "firebrick", "p < 0.05" = "#E69F00", "n.s." = "grey60")
      ) +
      labs(
        title    = paste("Estimated Effects on", y_axis_label),
        subtitle = "Points are model coefficients (effect size). Bars show 95% confidence intervals.",
        caption  = "Terms not crossing 0 are statistically significant.",
        x = "Coefficient Estimate (Effect Size)",
        y = "Predictor",
        color = "Significance"
      ) +
      theme_publication(base_size = 12) +
      theme(
        legend.position      = "right",
        panel.grid.major.y   = element_line(color = "grey93"),
        panel.grid.major.x   = element_line(color = "grey93")
      )
    base::print(p_coef)
    ggsave(file.path(var_plot_dir, paste0("p_coef_", short_resp_var, ".tiff")),
           p_coef, width = 9, height = 6, dpi = 300, compression = "lzw")
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 10: Clonal Waterfall Plot (Individual Responder Profiles)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  if (!is.null(model_outputs$results_tables$Individual_Slopes_BLUPs)) {
    logr::log_print(paste("...generating Waterfall plot for", current_resp_var))
    
    # Identify the column holding your clone names
    clone_col <- cfg_vars$secondary_group_var 
    
    wf_data <- model_outputs$results_tables$Individual_Slopes_BLUPs %>%
      dplyr::mutate(!!sym(primary_var) := as.character(!!sym(primary_var))) %>%
      dplyr::left_join(label_lookup %>% mutate(across(everything(), as.character)), by = primary_var) %>%
      restore_all_factors() %>%
      # Sort by slope to create the waterfall effect
      dplyr::arrange(desc(Individual_Slope)) %>%
      # Make Rank a character string so ggplot treats it as a discrete axis
      dplyr::mutate(Rank = as.character(row_number()), 
                    Direction = ifelse(Individual_Slope > 0, "Expansion", "Contraction/Stable"))
    
    color_var_wf <- label_pub
    
    # Create a "dictionary" that safely matches the Rank (1, 2, 3) to the Clone Name
    clone_labels <- setNames(as.character(wf_data[[clone_col]]), wf_data$Rank)
    
    p_waterfall <- ggplot(wf_data,
                          aes(x = reorder(Rank, -Individual_Slope),
                              y = Individual_Slope,
                              fill = !!sym(color_var_wf))) +
      geom_col(color = "black", linewidth = 0.2, alpha = 0.85) +
      geom_hline(yintercept = 0, color = "black", linewidth = 0.9) +
      # Annotate the zero line
      annotate("text", x = Inf, y = 0, label = "No expansion  ",
               hjust = 1, vjust = -0.5, size = 3, color = "grey40") +
      labs(
        title    = paste("Individual Responder Waterfall:", y_axis_label),
        subtitle = "Each bar is one clone/replicate ranked by its rate of change.",
        x        = NULL,
        y        = paste0("Rate of Change per ", stringr::str_to_title(config$key_variables$time_variable)),
        fill     = stringr::str_to_title(color_var_wf)
      ) +
      scale_x_discrete(labels = clone_labels) +
      theme_publication(base_size = 14) +
      theme(
        axis.text.x        = element_text(angle = 45, hjust = 1, size = 9),
        axis.ticks.x       = element_line(color = "black"),
        panel.grid.major.x = element_blank()
      )
    
    p_waterfall <- apply_smart_palette(p_waterfall, color_var_wf)
    
    assign(paste0("p_waterfall_", short_resp_var), p_waterfall)
    plot_database[[paste0("p_waterfall_", short_resp_var)]] <- p_waterfall
    base::print(p_waterfall)
    ggsave(file.path(var_plot_dir, paste0("p_waterfall_", short_resp_var, ".tiff")), p_waterfall, width = 12, height = 6, dpi = 300, compression = "lzw")
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 11: Time-Point "Pseudo-Raincloud" (Distribution Shifts)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  logr::log_print(paste("...generating Distribution Shift plot for", current_resp_var))
  
  rc_data <- modeling_data %>%
    # 1. Bin the time variable (e.g., group Day 13, 14, 15 into "14")
    dplyr::mutate(time_binned = round(!!sym(time_var) / time_bin_width) * time_bin_width) %>%
    # 2. Convert the BINNED time to a factor so we get clean, grouped buckets on the x-axis
    dplyr::mutate(Time_Factor = factor(time_binned)) %>%
    # 3. Prevent duplicate column errors
    dplyr::select(-dplyr::any_of(label_pub)) %>%
    dplyr::left_join(label_lookup %>% mutate(across(everything(), as.character)), by = primary_var) %>%
    restore_all_factors()
  
  color_var_rc <- if(is_treatment_exp) cfg_vars$secondary_group_var else label_pub
  
  p_raincloud <- ggplot(rc_data,
                        aes(x = Time_Factor, y = !!sym(current_resp_var),
                            fill = !!sym(color_var_rc), color = !!sym(color_var_rc))) +
    { if (grepl("change|_change", current_resp_var, ignore.case = TRUE))
      geom_hline(yintercept = 0, linetype = "dashed", colour = "grey50", linewidth = 0.5) } +
    # Violin showing full distribution shape
    geom_violin(alpha = 0.25, color = NA, trim = TRUE,
                position = position_dodge(width = 0.85)) +
    # Boxplot for median + IQR — narrow so it sits inside the violin
    geom_boxplot(width = 0.15, alpha = 0.7, outlier.shape = NA,
                 color = "black",
                 position = position_dodge(width = 0.85)) +
    # Individual points jittered within dodge
    geom_point(alpha = 0.4, size = 1.2,
               position = position_jitterdodge(jitter.width = 0.08, dodge.width = 0.85)) +
    labs(
      title    = paste("Population Distribution Shift over Time:", y_axis_label),
      subtitle = paste("Data binned into", time_bin_width, time_var,
                       "intervals. Violin = distribution; box = IQR; points = bio-reps."),
      x     = paste(stringr::str_to_title(time_var), "(Binned)"),
      y     = y_axis_label,
      fill  = stringr::str_to_title(color_var_rc),
      color = stringr::str_to_title(color_var_rc)
    ) +
    theme_publication(base_size = 14) +
    theme(legend.position = "right")
  
  # 1. Apply the Facet Structure (Only if it's a treatment experiment)
  if (is_treatment_exp) {
    p_raincloud <- p_raincloud + facet_wrap(vars(!!sym(label_pub)))
  }
  
  # 2. Apply the Smart Palette (Handles all the color logic automatically)
  p_raincloud <- apply_smart_palette(p_raincloud, color_var_rc)
  
  assign(paste0("p_raincloud_", short_resp_var), p_raincloud)
  plot_database[[paste0("p_raincloud_", short_resp_var)]] <- p_raincloud
  base::print(p_raincloud)
  ggsave(file.path(var_plot_dir, paste0("p_raincloud_", short_resp_var, ".tiff")), p_raincloud, width = 12, height = 7, dpi = 300, compression = "lzw")
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  # PLOT 12: Volcano Plot (Effect Size vs. Significance)
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
  if (!is.null(model_outputs$results_tables$Rate_Comparisons)) {
    logr::log_print(paste("...generating Volcano plot for", current_resp_var))
    
    volc_data <- model_outputs$results_tables$Rate_Comparisons %>%
      dplyr::filter(!is.na(p_value_adj) & !is.na(Rate_Difference)) %>%
      dplyr::mutate(
        p_value_adj         = suppressWarnings(as.numeric(p_value_adj)),
        Rate_Difference = suppressWarnings(as.numeric(Rate_Difference)),
        logP            = -log10(p_value_adj),
        effect_thresh   = 0.01,
        Sig_Label       = dplyr::case_when(
          p_value_adj < 0.05 & Rate_Difference >  effect_thresh ~ "Significant Increase",
          p_value_adj < 0.05 & Rate_Difference < -effect_thresh ~ "Significant Decrease",
          TRUE ~ "Not Significant"
        ),
        # FIX: Directly target the raw 'contrast' column and swap the spaced-hyphen
        Contrast_Label  = stringr::str_replace(contrast, " - ", " vs ")
      )
    
    if (nrow(volc_data) > 0) {
      # Symmetric x-axis range
      max_diff <- max(abs(volc_data$Rate_Difference), na.rm = TRUE) * 1.1
      
      p_volcano <- ggplot(volc_data,
                          aes(x = Rate_Difference, y = logP, color = Sig_Label)) +
        # Threshold lines
        geom_vline(xintercept = c(-0.01, 0.01),
                   linetype = "dashed", color = "grey55", linewidth = 0.5) +
        geom_hline(yintercept = -log10(0.05),
                   linetype = "dashed", color = "grey55", linewidth = 0.5) +
        # Threshold annotations
        annotate("text", x =  0.01, y = 0, label = "Effect threshold",
                 hjust = -0.1, vjust = -0.5, size = 2.8, color = "grey45") +
        annotate("text", x = -max_diff, y = -log10(0.05),
                 label = "p = 0.05", hjust = -0.1, vjust = -0.5,
                 size = 2.8, color = "grey45") +
        geom_point(size = 3.5, alpha = 0.85) +
        ggrepel::geom_text_repel(
          data = dplyr::filter(volc_data, Sig_Label != "Not Significant"),
          aes(label = Contrast_Label),
          size = 3, show.legend = FALSE, max.overlaps = 15,
          box.padding = 0.4, min.segment.length = 0.2
        ) +
        scale_color_manual(values = c(
          "Significant Increase" = "firebrick",
          "Significant Decrease" = "dodgerblue3",
          "Not Significant"      = "grey80"
        )) +
        scale_x_continuous(limits = c(-max_diff, max_diff)) +
        labs(
          title    = paste("Rate Comparison Volcano Plot:", y_axis_label),
          subtitle = "Effect size (x) vs. statistical significance (y). Dashed lines = thresholds.",
          x        = "Difference in Expansion Rate",
          y        =  expression(-log[10]("adj. p-value")),
          color    = "Direction"
        ) +
        theme_publication(base_size = 14) +
        theme(legend.position = "bottom")
      p_volcano <- p_volcano + 
        
        # 1. Expand the Axes Padding (Fixes the sliced dots)
        # 'mult' adds percentage-based padding. 
        # c(0.05, 0.15) means 5% padding at the bottom, and 15% at the top.
        scale_y_continuous(expand = expansion(mult = c(0.05, 0.15))) +

        # 2. Turn off Clipping (A safety net for edge-cases)
        # This tells ggplot: "Even if a dot or label slightly crosses the axis line, don't chop it off."
        coord_cartesian(clip = "off") +
        
        # 3. Add Physical Plot Margins (Fixes the cut-off text on the right)
        # margin(top, right, bottom, left)
        theme(plot.margin = margin(t = 10, r = 30, b = 10, l = 10))
      assign(paste0("p_volcano_", short_resp_var), p_volcano)
      plot_database[[paste0("p_volcano_", short_resp_var)]] <- p_volcano
      base::print(p_volcano)
      ggsave(file.path(var_plot_dir, paste0("p_volcano_", short_resp_var, ".tiff")), p_volcano, width = 10, height = 8, dpi = 300, compression = "lzw")
    }
  }
  
} # End Main Response Variable Loop

#=============================================================================#
# PART 4B: CROSS-METRIC CORRELATION (Mode vs. Instability) ####
#=============================================================================#
logr::log_print("\n--- Generating Cross-Metric Correlation Plots ---", console = TRUE)

var_x <- "mode_change"
var_y <- "instability_index_change"

# Check if both variables were modeled and produced BLUPs
if (!is.null(all_model_outputs[[var_x]]$results_tables$Individual_Slopes_BLUPs) && 
    !is.null(all_model_outputs[[var_y]]$results_tables$Individual_Slopes_BLUPs)) {
  
  blup_x <- all_model_outputs[[var_x]]$results_tables$Individual_Slopes_BLUPs %>%
    dplyr::rename(Slope_X = Individual_Slope) %>%
    dplyr::select(dplyr::any_of(c(primary_var, cfg_vars$secondary_group_var, re_cross, "Method")), Slope_X)
  
  blup_y <- all_model_outputs[[var_y]]$results_tables$Individual_Slopes_BLUPs %>%
    dplyr::rename(Slope_Y = Individual_Slope) %>%
    dplyr::select(dplyr::any_of(c(primary_var, cfg_vars$secondary_group_var, re_cross, "Method")), Slope_Y)
  
  # Join them together based on their shared grouping variables
  join_cols <- intersect(names(blup_x)[names(blup_x) != "Slope_X"], names(blup_y)[names(blup_y) != "Slope_Y"])
  corr_data <- dplyr::inner_join(blup_x, blup_y, by = join_cols) %>%
    dplyr::mutate(!!sym(primary_var) := as.character(!!sym(primary_var))) %>%
    dplyr::left_join(label_lookup %>% mutate(across(everything(), as.character)), by = primary_var) %>%
    restore_all_factors()
  
  if (nrow(corr_data) > 3) {
    color_var_corr <- if(is_treatment_exp) cfg_vars$secondary_group_var else label_pub
    color_var_label <- if(is_treatment_exp) cfg_vars$secondary_group_var else cfg_vars$primary_group_var
    
    # Use pretty labels for axis titles if available
    label_x <- response_labels[["mode_change"]] %||% "Mode Change Slope"
    label_y <- response_labels[["instability_index"]] %||% "Instability Index Slope"
    
    p_cross_corr <- ggplot(corr_data, aes(x = Slope_X, y = Slope_Y)) +
      # Overall regression (dashed, behind points)
      geom_smooth(method = "lm", color = "grey40", linetype = "dashed",
                  fill = "grey85", alpha = 0.3, linewidth = 0.8) +
      # Per-group regression lines (solid, coloured)
      geom_smooth(aes(color = !!sym(color_var_corr)),
                  method = "lm", se = FALSE, linewidth = 0.8, alpha = 0.7) +
      geom_point(aes(fill = !!sym(color_var_corr)),
                 shape = 21, color = "black", size = 4, alpha = 0.85) +
      # Overall Pearson r
      ggpubr::stat_cor(method = "pearson", label.x.npc = "left",
                       label.y.npc = "top", size = 4.5, color = "grey30") +
      labs(
        title    = paste("Cross-Metric Correlation:",
                         gsub(" Slope", "", label_x), "vs.",
                         gsub(" Slope", "", label_y)),
        subtitle = "Do clones that expand fastest also become the most unstable? Dashed = overall trend.",
        x    = paste("Rate of Change in", label_x),
        y    = paste("Rate of Change in", label_y),
        fill  = stringr::str_to_title(color_var_label),
        color = stringr::str_to_title(color_var_label)
      ) +
      theme_publication(base_size = 14) +
      theme(legend.position = "right")
    
    p_cross_corr <- apply_smart_palette(p_cross_corr, color_var_corr)
    
    if (is_treatment_exp) p_cross_corr <- p_cross_corr + facet_wrap(vars(!!sym(label_pub))) 
    
    plot_database[["p_cross_corr_mode_vs_iic"]] <- p_cross_corr
    base::print(p_cross_corr)
    
    corr_dir <- file.path(output_dir, "04_correlation")
    dir.create(corr_dir, showWarnings = FALSE)
    ggsave(file.path(corr_dir, "p_cross_corr_mode_vs_iic.tiff"), p_cross_corr, width = 10, height = 7, dpi = 300, compression = "lzw")
  }
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# --- THE ULTIMATE EXCEL POLISH: Clean, Format, and Humanize Tables --- ####
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
logr::log_print("...Applying final polish to statistical tables before Excel export.", console = TRUE)

format_p_val <- function(p) {
  dplyr::case_when(
    is.na(p) ~ NA_character_,
    p < 0.001 ~ "< 0.001",
    TRUE ~ sprintf("%.4f", p)
  )
}

# Provide a safety net of every possible name statistical packages use for p-values
pval_cols <- c("P_Value", "p.value", "p-value", "Pr(>F)", "Pr(>|t|)", "Pr(>Chisq)", "p_value")
pval_adj_cols <- c(
  # Standard base R & broom naming
  "p.adj", "p_adj", "padj", 
  "p.adjusted", "p_adjusted", 
  
  # Emmeans & Custom explicit naming
  "p.value.adj", "p_value_adj", "P_Value_Adj",
  "p_value_tukey", "p_value_bonferroni",
  
  # Bioinformatics packages (limma, DESeq2, edgeR)
  "adj.P.Val", "adj.p.value", "adj_p_value",
  
  # False Discovery Rate / Q-value equivalents
  "FDR", "fdr", "q.value", "qvalue", "q_value"
)


for (current_resp_var in response_vars) {
  model_outputs <- all_model_outputs[[current_resp_var]]
  if (is.null(model_outputs)) next
  
  if ("Model_Coefficients" %in% names(model_outputs$results_tables)) {
    mc <- model_outputs$results_tables$Model_Coefficients
    baseline_info <- c(paste0(config$key_variables$time_variable, "=0"))
    for (var in config$key_variables$model_fixed_effects) {
      if (var %in% names(modeling_data)) {
        ref <- if(is.factor(modeling_data[[var]])) levels(modeling_data[[var]])[1] else sort(unique(as.character(modeling_data[[var]])))[1]
        baseline_info <- c(baseline_info, paste0(stringr::str_to_title(var), " [", ref, "]"))
      }
    }
    baseline_string <- paste0("Baseline: ", paste(baseline_info, collapse = " | "))
    
    mc <- mc %>%
      dplyr::filter(!is.na(Estimate)) %>%
      dplyr::mutate(Parameter = ifelse(Parameter == "(Intercept)", baseline_string, Parameter)) %>%
      dplyr::mutate(Parameter = stringr::str_replace_all(Parameter, ":", " × ")) %>%
      # FIX: Crash-proof p-value formatting
      dplyr::mutate(dplyr::across(dplyr::any_of(c(pval_cols, pval_adj_cols)), ~ format_p_val(as.numeric(.))))
    
    # Standardize the column name back to P_Value for Excel
    names(mc)[names(mc) %in% pval_cols] <- "P_Value"
    model_outputs$results_tables$Model_Coefficients <- mc
  }
  
  if ("ANOVA_Summary" %in% names(model_outputs$results_tables)) {
    aov_tbl <- model_outputs$results_tables$ANOVA_Summary
    if ("Sum Sq" %in% names(aov_tbl)) {
      total_ss <- sum(aov_tbl$`Sum Sq`, na.rm = TRUE)
      aov_tbl <- aov_tbl %>%
        dplyr::mutate(Pct_Variance_Explained = sprintf("%.1f%%", (`Sum Sq` / total_ss) * 100)) %>%
        dplyr::relocate(Pct_Variance_Explained, .after = `Sum Sq`)
    }
    
    aov_tbl <- aov_tbl %>%
      dplyr::mutate(Factor = stringr::str_replace_all(Factor, ":", " × ")) %>%
      # FIX: Crash-proof p-value formatting
      dplyr::mutate(dplyr::across(dplyr::any_of(pval_cols), ~ format_p_val(as.numeric(.))))
    
    names(aov_tbl)[names(aov_tbl) %in% pval_cols] <- "P_Value"
    model_outputs$results_tables$ANOVA_Summary <- aov_tbl
  }
  
  if ("Rate_Comparisons" %in% names(model_outputs$results_tables)) {
    rc <- model_outputs$results_tables$Rate_Comparisons
    rc <- rc %>%
      dplyr::filter(!is.na(Rate_Difference)) %>%
      dplyr::mutate(
        # FIX: Look for " - " (space-hyphen-space) so it ignores the hyphen in Wild-Type
        Group_A = stringr::str_trim(stringr::str_split_fixed(contrast, " - ", n = 2)[, 1]),
        Group_B = stringr::str_trim(stringr::str_split_fixed(contrast, " - ", n = 2)[, 2])
      ) %>%
      dplyr::relocate(Group_A, Group_B, .after = contrast) %>%
      dplyr::select(-contrast) %>%
      # FIX: Crash-proof p-value formatting
      dplyr::mutate(dplyr::across(dplyr::any_of(pval_cols), ~ format_p_val(as.numeric(.))))
    
    names(rc)[names(rc) %in% pval_cols] <- "P_Value"
    names(rc)[names(rc) %in% pval_adj_cols] <- "Adj_P_Value"
    model_outputs$results_tables$Rate_Comparisons <- rc
  }
  
  for (tbl_name in names(model_outputs$results_tables)) {
    tbl <- model_outputs$results_tables[[tbl_name]]
    if ("Expansion_Rate" %in% names(tbl)) tbl <- tbl %>% dplyr::filter(!is.na(Expansion_Rate))
    model_outputs$results_tables[[tbl_name]] <- tbl
  }
  
  if ("Random_Effects_Variance" %in% names(model_outputs$results_tables)) {
    rev_tbl <- model_outputs$results_tables$Random_Effects_Variance
    rev_tbl <- rev_tbl %>%
      dplyr::mutate(
        ICC_Percent = ifelse(!is.na(ICC_Percent), sprintf("%.1f%%", ICC_Percent), NA_character_),
        Source = stringr::str_replace(Source, "\\.\\(Intercept\\)", " (Baseline Variance)"),
        Source = stringr::str_replace(Source, "Residual", "Residual (Unexplained Noise)")
      ) %>%
      dplyr::mutate(dplyr::across(dplyr::where(is.numeric), ~ round(.x, 4)))
    model_outputs$results_tables$Random_Effects_Variance <- rev_tbl
  }
  
  if ("Individual_Slopes_BLUPs" %in% names(model_outputs$results_tables)) {
    blup_tbl <- model_outputs$results_tables$Individual_Slopes_BLUPs
    blup_tbl <- blup_tbl %>%
      # FIX: Use any_of() so it doesn't crash if the Naive method was used
      dplyr::mutate(dplyr::across(dplyr::any_of(c("Individual_Slope", "Intercept_Deviation", "Slope_Deviation")), ~ round(.x, 5)))
    model_outputs$results_tables$Individual_Slopes_BLUPs <- blup_tbl
  }
  
  all_model_outputs[[current_resp_var]] <- model_outputs
}

#=============================================================================#
# PART 5: EXPORT STATISTICAL MODEL RESULTS & SAVE RDATA ####
#=============================================================================#
logr::log_print("\n--- Starting PART 5: Exporting Statistical Model Results ---", console = TRUE)

library(openxlsx)
wb <- createWorkbook()
title_style <- createStyle(textDecoration = "bold", fontSize = 14)

for (var_name in names(all_model_outputs)) {
  model_obj <- all_model_outputs[[var_name]]
  if (is.null(model_obj)) next
  
  addWorksheet(wb, var_name)
  row_pos <- 1
  results_for_var <- model_obj$results_tables
  
  for (table_name in names(results_for_var)) {
    tbl <- results_for_var[[table_name]]
    if (is.null(tbl) || nrow(tbl) == 0) next
    
    writeData(wb, sheet = var_name, x = table_name, startRow = row_pos, colNames = FALSE)
    addStyle(wb, sheet = var_name, style = title_style, rows = row_pos, cols = 1, gridExpand = TRUE)
    row_pos <- row_pos + 1
    writeData(wb, sheet = var_name, x = tbl, startRow = row_pos, withFilter = FALSE)
    row_pos <- row_pos + nrow(tbl) + 3
  }
}

saveWorkbook(wb, file.path(output_dir, "Statistical_Model_Output.xlsx"), overwrite = TRUE)

# ================================================================= #
# 5C. CREATE "FINAL_SLOPES_WIDE" (The Replacement for all_group_slopes)
# ================================================================= #
tryCatch({ logr::log_print("...Creating consolidated wide slope table (BLUPs) for next script.") }, error = function(e) message("...Creating consolidated wide slope table (BLUPs) for next script."))

all_group_slopes <- NULL

for (var_name in names(all_model_outputs)) {
  if (!is.null(all_model_outputs[[var_name]]) && !is.null(all_model_outputs[[var_name]]$results_tables$Individual_Slopes_BLUPs)) {
    
    blup_table <- all_model_outputs[[var_name]]$results_tables$Individual_Slopes_BLUPs
    
    # --- THE FIX: Exclude ALL model-specific statistical columns ---
    cols_to_exclude <- c(
      "Individual_Slope", 
      "Intercept_Deviation", 
      "Slope_Deviation",     # <-- Was missing
      "Group_Base_Slope",    # <-- Was missing
      "Method"
    )
    current_group_cols <- setdiff(names(blup_table), cols_to_exclude)
    # ---------------------------------------------------------------
    
    if ("value_renaming" %in% names(config)) {
      rename_map <- unlist(config$value_renaming)
      blup_table <- blup_table %>%
        dplyr::mutate(dplyr::across(dplyr::any_of(current_group_cols), ~ dplyr::coalesce(unname(rename_map[as.character(.)]), as.character(.))))
    }
    
    to_merge <- blup_table %>%
      dplyr::mutate(across(dplyr::any_of(current_group_cols), as.character)) %>%
      dplyr::select(dplyr::any_of(current_group_cols), Individual_Slope) %>%
      dplyr::rename(!!paste0("slope_", var_name) := Individual_Slope)
    
    if (is.null(all_group_slopes)) {
      all_group_slopes <- to_merge
    } else {
      join_keys <- intersect(names(all_group_slopes), names(to_merge))
      all_group_slopes <- dplyr::left_join(all_group_slopes, to_merge, by = join_keys)
    }
  }
}
# --- Save RData for Script 03 ---
rdata_path <- file.path(output_dir, "processing_complete.RData")
save.image(file = rdata_path)
logr::log_print(paste("R environment with all plot objects saved to:", rdata_path))

# --- Finalize ---
logr::log_print("\n--- SCRIPT 02: MODELING & PLOTTING FINISHED SUCCESSFULLY ---", console = TRUE)
logr::log_close()
