#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                             #
#         Advanced Fragment Analysis Workflow - Script 6                      #
#                           (HTML Report Generation)                          #
#                                                                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#=============================================================================#
# PART 0: SETUP AND LOAD DATA
#=============================================================================#

# 0A. Load Libraries ----
packages <- c("tidyverse", "here", "yaml", "logr", "rmarkdown", "plotly", "kableExtra", "knitr", "htmltools", "ragg")
installed_packages <- packages %in% rownames(utils::installed.packages())
if (any(installed_packages == FALSE)) {
  utils::install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

# 0B. Load Configuration & Find Data ----
if (!file.exists(here::here("config.yml"))) stop("CRITICAL ERROR: 'config.yml' not found.")
config <- yaml::read_yaml(here::here("config.yml"))

# Bulletproof Latest Directory Search (Matches Script 05)
analysis_base_dir <- here::here(config$paths$output_dir_base)
all_dirs <- list.dirs(path = analysis_base_dir, full.names = TRUE, recursive = FALSE)
all_dirs <- setdiff(all_dirs, analysis_base_dir)

if(length(all_dirs) == 0) stop("CRITICAL ERROR: No analysis output folders found.")
target_dirs <- all_dirs[grepl("_Analysis_v", basename(all_dirs))]
if(length(target_dirs) == 0) stop("CRITICAL ERROR: No '_Analysis_v' folders found.")

latest_analysis_dir <- max(target_dirs)

# Load the Environment
rdata_path <- file.path(latest_analysis_dir, "processing_complete.RData")
if (!file.exists(rdata_path)) stop("ERROR: 'processing_complete.RData' not found.")
load(rdata_path)

# (Re)open log to append
try(logr::log_close(), silent = TRUE)
log_path <- file.path(output_dir, "analysis_log.log")
lf <- logr::log_open(log_path, show_notes = FALSE, logdir = FALSE)

logr::log_print("\n\n--- SCRIPT 06: HTML REPORT GENERATION INITIALIZED ---", console = TRUE)

#=============================================================================#
# PART 1: RENDER HTML REPORT
#=============================================================================#
logr::log_print("\n--- Rendering HTML Interactive Report ---", console = TRUE)

# Define the path to your R Markdown template
report_path <- here::here("Fragment_Analysis_Report.Rmd")


if (file.exists(report_path)) {
  
  # Give it a clean, timestamped name
  output_filename <- paste0("Interactive_Report_", format(Sys.time(), "%Y%m%d_%H%M"), ".html")
  
  tryCatch({
    rmarkdown::render(
      input = report_path, 
      output_dir = output_dir,        # Saves it right into today's analysis folder
      output_file = output_filename,
      envir = globalenv(),            # CRITICAL: Passes the loaded RData environment to the markdown!
      quiet = TRUE
    )
    logr::log_print(paste("✅ Successfully generated HTML report:", output_filename), console = TRUE)
    logr::log_print(paste("Saved in:", output_dir), console = TRUE)
    
  }, error = function(e) {
    logr::log_print(paste("❌ ERROR rendering HTML report:", e$message), console = TRUE)
  })
  
} else {
  logr::log_print(paste("❌ WARNING: RMarkdown template not found at:", report_path), console = TRUE)
  logr::log_print("Please ensure 'Fragment_Analysis_Report.Rmd' is in your main project directory.", console = TRUE)
}

logr::log_print("\n--- SCRIPT 06 COMPLETE ---", console = TRUE)
logr::log_close()
