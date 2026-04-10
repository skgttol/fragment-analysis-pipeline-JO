#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                             #
#         Advanced Fragment Analysis Workflow - Script 4 of 4                 #
#                   (Combined Figure Generation)                              #
#                                                                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# DESCRIPTION:
# This script LOADS the final "processing_complete.RData" environment.
# It automatically "harvests" all plots from Scripts 02 & 03 into a 
# single list ('all_plots') to make combining figures easy.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#=============================================================================#
# PART 0: SETUP AND LOAD DATA
#=============================================================================#

# 0A. Load Libraries ----
packages <- c("tidyverse", "here", "yaml", "logr", "patchwork", 
              "ggdendro", "dendextend", "ggplot2")
installed_packages <- packages %in% rownames(utils::installed.packages())
if (any(installed_packages == FALSE)) utils::install.packages(packages[!installed_packages])
invisible(lapply(packages, library, character.only = TRUE))

# 0B. Load Configuration & Data ----
if (!file.exists(here::here("config.yml"))) stop("config.yml not found.")
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
print(paste("Data loaded successfully from:", rdata_path))

# Ensure functions are loaded
func_path <- "C:/Users/skgttol/OneDrive - University College London/PhD/PhD_Thesis/02_Data/Template_RScripts/CAGsizing_Flexible/functions.R"
if (!file.exists(func_path)) {
  if(file.exists("functions.R")) func_path <- "functions.R" else stop("functions.R not found.")
}
source(func_path)

# (Re)open log
try(logr::log_close(), silent = TRUE)
log_path <- file.path(output_dir, "analysis_log.log") 
lf <- logr::log_open(log_path, show_notes = FALSE, logdir = FALSE) 

logr::log_print("\n\n--- SCRIPT 04: FIGURE ASSEMBLY INITIALIZED ---", console = TRUE)

figure_dir <- here::here(latest_analysis_dir, "combined_figures")
dir.create(figure_dir, showWarnings = FALSE)


#=============================================================================#
# PART 1: THE "PLOT HARVESTER" (Builds the Database)
#=============================================================================#
logr::log_print("\n--- Building Plot Registry (Database) ---", console = TRUE)

# 1. Initialize Master List
all_plots <- list(
  descriptive = list(), # From Script 02
  modeling    = list(), # From Script 02
  qc          = list(), # From Script 02
  exploratory = list()  # From Script 03 (PCA, Corr)
)

# 2. Harvest Script 02 Plots 
if (exists("plot_database")) {
  names_db <- names(plot_database)
  
  # Descriptive (+ Rainclouds)
  desc_names <- names_db[grep("^p_avg|^p_detail|^p_profiles|^p_heatmap|^p_raincloud", names_db)]
  all_plots$descriptive <- plot_database[desc_names]
  
  # Modeling (+ Waterfalls)
  model_names <- names_db[grep("^p_preds|^p_overall|^p_slopes|^p_clone_trends|^p_rep_profiles|^p_waterfall", names_db)]
  all_plots$modeling <- plot_database[model_names]
  
  # QC (+ Volcano)
  qc_names <- names_db[grep("^p_qc|^p_clone_het|^p_volcano", names_db)]
  all_plots$qc <- plot_database[qc_names]
  
  # Cross-Correlation
  cross_names <- names_db[grep("^p_cross_corr", names_db)]
  all_plots$exploratory <- plot_database[cross_names]
}

# 3. Harvest Script 03 Plots (Loose objects in environment)
# Define known Script 03 plot names
script3_patterns <- c("^p_biplot", "^p_scree", "^p_corr", "^p_model_fit", "^p_loadings", "^dend")

# Find them in the loaded environment
env_objects <- ls()
for (pat in script3_patterns) {
  matches <- env_objects[grep(pat, env_objects)]
  for (m in matches) {
    obj <- get(m)
    # Check if it's a plot or dendrogram
    if (inherits(obj, "ggplot") || inherits(obj, "dendrogram") || inherits(obj, "ggmatrix")) {
      all_plots$exploratory[[m]] <- obj
    }
  }
}

# 4. PRINT THE MENU
logr::log_print("-----------------------------------------------------------", console=TRUE)
logr::log_print("               AVAILABLE PLOT DATABASE                     ", console=TRUE)
logr::log_print("-----------------------------------------------------------", console=TRUE)
logr::log_print("Use 'all_plots$category$plot_name' to select graphs.", console=TRUE)
logr::log_print("", console=TRUE)

#print(str(all_plots, max.level = 2, list.len = 100))

logr::log_print("-----------------------------------------------------------", console=TRUE)


#=============================================================================#
# PART 2: AUTOMATIC SINGLE-VARIABLE FIGURES (Trend + Slope)
#=============================================================================#
logr::log_print("\n--- Generating Single-Variable Combined Figures ---", console = TRUE)

# We will look for the "Averaged" versions of the plots we made in Script 2
for (target_var in response_vars) {
  short_target <- response_shortnames[[target_var]] %||% target_var
  
  # Try to find the plots in the database (using the names from the new Script 2)
  plot_trend_name <- paste0("p_overall_2a_", short_target)      # Averaged Trend
  plot_slope_name <- paste0("p_slopes_averaged_", short_target) # Averaged Slopes
  
  # # Fallbacks just in case older plot names exist
  # if (is.null(all_plots$modeling[[plot_trend_name]])) plot_trend_name <- paste0("p_preds_", short_target)
  # if (is.null(all_plots$modeling[[plot_slope_name]])) plot_slope_name <- paste0("p_slopes_", short_target)
  
  p_trend <- all_plots$modeling[[plot_trend_name]]
  p_slope <- all_plots$modeling[[plot_slope_name]]
  
  if (!is.null(p_trend) && !is.null(p_slope)) {
    logr::log_print(paste("...Building combined figure for:", target_var))
    
    # Wrap titles to prevent runoff
    p_trend <- p_trend + labs(
      title = stringr::str_wrap(p_trend$labels$title, width = 50),
      subtitle = stringr::str_wrap(p_trend$labels$subtitle, width = 60)
    ) + theme_publication(base_size = 18)
    
    p_slope <- p_slope + labs(
      title = stringr::str_wrap(p_slope$labels$title, width = 50),
      subtitle = stringr::str_wrap(p_slope$labels$subtitle, width = 60)
    ) + theme_publication(base_size = 18) +
      theme(legend.position = "none")
    
    # Combine using Patchwork
    fig_single <- (p_trend + p_slope) + 
      plot_layout(widths = c(2, 1), guides = "collect")
      
    
    ggsave(file.path(figure_dir, paste0("Fig_Combined_", short_target, ".tiff")), 
           fig_single, width = 20, height = 8, dpi = 300, compression = "lzw")
  } else {
    logr::log_print(paste("...Skipping", target_var, "- Required plots not found in database."))
  }
}

#=============================================================================#
# PART 3: AUTOMATIC PAIRED FIGURES (4-Panel Layouts)
#=============================================================================#
logr::log_print("\n--- Generating Paired-Variable (4-Panel) Figures ---", console = TRUE)

# Define which variables you want to plot together as pairs. 
# You can add as many pairs to this list as you want!
paired_analyses <- list(
  c("mode_change", "instability_index_change")
  # c("another_var", "yet_another_var") # Add more pairs here if needed
)

for (pair in paired_analyses) {
  var_A <- pair[1]
  var_B <- pair[2]
  
  short_A <- response_shortnames[[var_A]] %||% var_A
  short_B <- response_shortnames[[var_B]] %||% var_B
  
  # Grab the 4 plots
  p1 <- all_plots$modeling[[paste0("p_overall_2a_", short_A)]] %||% all_plots$modeling[[paste0("p_preds_", short_A)]]
  p2 <- all_plots$modeling[[paste0("p_slopes_averaged_", short_A)]] %||% all_plots$modeling[[paste0("p_slopes_", short_A)]]
  p3 <- all_plots$modeling[[paste0("p_overall_2a_", short_B)]] %||% all_plots$modeling[[paste0("p_preds_", short_B)]]
  p4 <- all_plots$modeling[[paste0("p_slopes_averaged_", short_B)]] %||% all_plots$modeling[[paste0("p_slopes_", short_B)]]
  
  if (!is.null(p1) && !is.null(p2) && !is.null(p3) && !is.null(p4)) {
    logr::log_print(paste("...Building paired 4-panel figure for:", var_A, "&", var_B))
    
    safe_wrap <- function(p) {
      p + labs(
        title = if(!is.null(p$labels$title)) stringr::str_wrap(p$labels$title, width = 50) else NULL,
        subtitle = if(!is.null(p$labels$subtitle)) stringr::str_wrap(p$labels$subtitle, width = 60) else NULL
      )
    }
    
    p1 <- safe_wrap(p1) + theme_publication(base_size = 18)
    p2 <- safe_wrap(p2) + theme_publication(base_size = 18) + theme(legend.position = "none")
    p3 <- safe_wrap(p3) + theme_publication(base_size = 18)
    p4 <- safe_wrap(p4) + theme_publication(base_size = 18) + theme(legend.position = "none")
    
    # Note: Unifying Y-axes between completely different metrics (like Mode and IIC) 
    # is statistically dangerous because they are on different scales. 
    # This layout aligns them visually without forcing the Y-axes to match.
    
    fig_paired <- (p1 + p2) / (p3 + p4) + 
      plot_layout(widths = c(1.8, 1.2), guides = "collect") +
      theme(
        legend.position = "right", 
        legend.direction = "vertical",
        legend.justification = "top",
        legend.box.margin = margin(l = 10, t = 50) 
      )
      
    
    save_path <- file.path(figure_dir, paste0("Fig_Paired_", short_A, "_vs_", short_B, ".tiff"))
    ggsave(save_path, fig_paired, width = 22, height = 14, dpi = 300, compression = "lzw")
  } else {
    logr::log_print(paste("...Skipping paired figure for", var_A, "&", var_B, "- missing plots."))
  }
}

logr::log_print("\n--- SCRIPT 04 COMPLETE ---", console = TRUE)
logr::log_close()
