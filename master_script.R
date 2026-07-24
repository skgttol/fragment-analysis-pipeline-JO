#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                 MODULAR PIPELINE MASTER CONTROL SCRIPT                      #
#                 (Phase 1: Collect -> Phase 2: Execute)                      #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

# =============================================================================
# GLOBAL SETUP
# =============================================================================
if (!requireNamespace("here", quietly = TRUE)) install.packages("here")

message("Checking local environment...")
if (!file.exists(here::here("config.yml"))) {
  stop("❌ CRITICAL ERROR: 'config.yml' not found. Please ensure your working directory is set correctly.")
}

FORCE_UPDATE_SCRIPTS <- FALSE
# --- GITHUB REPOSITORY SETTINGS ---
# Replace 'YOUR_USERNAME' and 'YOUR_REPO' with your actual GitHub details.
# Ensure the branch name ('main' or 'master') is correct.
github_base_url <- "https://raw.githubusercontent.com/skgttol/fragment-analysis-pipeline-JO/main/"

# Master list of all required scripts
pipeline_scripts <- c(
  "functions.R",
  "analysis-scripts/00_generate_settings.R",
  "analysis-scripts/01_process_data.R",
  "analysis-scripts/02_generate_plots.R",
  "analysis-scripts/03_correlation_analysis.R",
  "analysis-scripts/04_figure_generation.R",
  "analysis-scripts/05_trace-from-txt.R",
  "analysis-scripts/06_render_html.R",
  "Fragment_Analysis_Report.Rmd"
)

# =============================================================================
# PHASE 1: COLLECT & VERIFY ALL SCRIPTS
# =============================================================================
script_dir <- here::here("analysis-scripts")
dir.create(script_dir, recursive = TRUE, showWarnings = FALSE)
github_script_dir <- here::here()


# Check if any scripts are physically missing locally
missing_locally <- any(!file.exists(file.path(github_script_dir, pipeline_scripts)))

if (missing_locally || FORCE_UPDATE_SCRIPTS) {
  message("\n--- PHASE 1: COLLECTING/UPDATING SCRIPTS FROM GITHUB ---")
  
  for (script in pipeline_scripts) {
    url <- paste0(github_base_url, script)
    dest <- file.path(github_script_dir, script)
    
    message(sprintf("Downloading %s...", script))
    
    dl_status <- try(download.file(url, destfile = dest, mode = "wb", # Changed to 'wb' to prevent Windows line break bugs
                                   quiet = TRUE), silent = TRUE)
    
    if (inherits(dl_status, "try-error") || dl_status != 0) {
      stop(sprintf("\n❌ ERROR: Failed to download '%s'.\nCheck your internet, URL, or ensure your token hasn't expired.", script))
    }
  }
  message("✅ Scripts updated and verified.\n")
} else {
  message("\n--- PHASE 1: SKIPPED DOWNLOAD (Using local script cache) ---")
}

# =============================================================================
# PHASE 2: EXECUTION ENGINE
# =============================================================================
message("\n--- PHASE 2: EXECUTING PIPELINE ---")

# Helper function to run a local script safely in an isolated environment
# Helper function to run a local script safely in an isolated environment
run_local_script <- function(script_name) {
  script_path <- file.path(script_dir, script_name)
  message("=======================================================")
  message(sprintf("▶ RUNNING: %s", script_name))
  message("=======================================================")
  
  # 1. Take a snapshot of the environment BEFORE the script runs
  pre_existing_vars <- ls(envir = .GlobalEnv)
  
  # 2. Run the script normally (local = FALSE)
  # This ensures save.image(), saveRDS(), and your helper functions work perfectly.
  source(script_path, local = FALSE)
  
  # 3. Identify all NEW variables created by this specific script
  current_vars <- ls(envir = .GlobalEnv)
  vars_to_delete <- setdiff(current_vars, pre_existing_vars)
  
  # 4. Aggressive Cleanup
  # We delete ONLY the variables this script just made, leaving your base setup intact
  rm(list = vars_to_delete, envir = .GlobalEnv) 
  
  gc(verbose = FALSE)                                 # Force Garbage Collection (Free RAM)
  
  # Suppress warnings for connections in case none are open
  suppressWarnings(closeAllConnections())             # Close dangling file read/writes
  graphics.off()                                      # Close all open plot devices
  
  message(sprintf("✅ SUCCESS: %s completed cleanly and memory was freed.\n", script_name))
}

# =============================================================================
# PIPELINE CONTROL PANEL
# Toggle these switches to TRUE or FALSE to control which scripts execute.
# =============================================================================

# [Setup] Usually run only once per project to generate the initial settings
RUN_00_SETTINGS    <- FALSE  

# [Core Pipeline]
RUN_01_PROCESS     <- TRUE   # Data Processing & QC
RUN_02_PLOTS       <- TRUE   # Statistical Modeling & Plotting
RUN_03_CORRELATION <- TRUE   # Correlation & Integration
RUN_04_FIGURES     <- TRUE   # Combined Figure Generation
RUN_05_TRACES      <- TRUE   # Trace Reconstruction & PDF Reports
RUN_06_HTML        <- TRUE   # Integrated HTML report


# =============================================================================
# EXECUTION BLOCKS
# =============================================================================
message("\n=======================================================")
message("🚀 STARTING PIPELINE EXECUTION")
message("=======================================================\n")


# --- 00: Settings Generation ---
if (RUN_00_SETTINGS) {
  tryCatch({
    run_local_script("00_generate_settings.R")
  }, error = function(e) stop(sprintf("\n❌ PIPELINE HALTED AT 00_settings_generation.R:\n%s", e$message), call. = FALSE))
} else {
  message("⏭️ SKIPPED: 00_generate_settings.R (Run-once setup)")
}

# --- 01: Process Data ---
if (RUN_01_PROCESS) {
  tryCatch({
    run_local_script("01_process_data.R")
  }, error = function(e) stop(sprintf("\n❌ PIPELINE HALTED AT 01_process_data.R:\n%s", e$message), call. = FALSE))
} else {
  message("⏭️ SKIPPED: 01_process_data.R")
}

# --- 02: Generate Plots ---
if (RUN_02_PLOTS) {
  tryCatch({
    run_local_script("02_generate_plots.R")
  }, error = function(e) stop(sprintf("\n❌ PIPELINE HALTED AT 02_generate_plots.R:\n%s", e$message), call. = FALSE))
} else {
  message("⏭️ SKIPPED: 02_generate_plots.R")
}

# --- 03: Correlation Analysis ---
if (RUN_03_CORRELATION) {
  tryCatch({
    run_local_script("03_correlation_analysis.R")
  }, error = function(e) stop(sprintf("\n❌ PIPELINE HALTED AT 03_correlation_analysis.R:\n%s", e$message), call. = FALSE))
} else {
  message("⏭️ SKIPPED: 03_correlation_analysis.R")
}

# --- 04: Figure Generation ---
if (RUN_04_FIGURES) {
  tryCatch({
    run_local_script("04_figure_generation.R")
  }, error = function(e) stop(sprintf("\n❌ PIPELINE HALTED AT 04_figure_generation.R:\n%s", e$message), call. = FALSE))
} else {
  message("⏭️ SKIPPED: 04_figure_generation.R")
}

# --- 05: Trace Generation ---
if (RUN_05_TRACES) {
  tryCatch({
    run_local_script("05_trace-from-txt.R")
  }, error = function(e) stop(sprintf("\n❌ PIPELINE HALTED AT 05_trace-from-txt.R:\n%s", e$message), call. = FALSE))
} else {
  message("⏭️ SKIPPED: 05_trace-from-txt.R")
}

# --- 06: HTML Report Generation ---
if (RUN_06_HTML) {
  tryCatch({
    run_local_script("06_render_html.R")
  }, error = function(e) stop(sprintf("\n❌ PIPELINE HALTED AT 06_render_html.R:\n%s", e$message), call. = FALSE))
} else {
  message("⏭️ SKIPPED: 06_render_html.R")
}

message("\n🎉 PIPELINE EXECUTION COMPLETELY FINISHED 🎉\n")
