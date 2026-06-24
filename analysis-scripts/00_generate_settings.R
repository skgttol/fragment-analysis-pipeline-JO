#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                             #
#            Advanced Fragment Analysis Workflow - Script 0 of 2              #
#                        (Settings File Generation)                           #
#                                                                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# DESCRIPTION:
# This script is OPTIONAL. It is used to generate the "Settings File" required
# by the external fragment analysis program (e.g., Romeo).
#
# USAGE:
# Run this script ONLY if you need to generate a new settings file based on
# your FSA files and Platemap.
#
# OUTPUT:
# An Excel file containing the sample list and analysis parameters, saved
# to the output directory.
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#=============================================================================#
# PART 0: SETUP
#=============================================================================#

# 0A. Load Libraries ----
packages <- c("tidyverse", "readxl", "here", "openxlsx", "janitor", "yaml", "logr")
installed_packages <- packages %in% rownames(utils::installed.packages())
if (any(installed_packages == FALSE)) {
  utils::install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

# 0B. Load Functions ----
#-------------------------#
if (!file.exists(here::here("functions.R"))) {
  stop("CRITICAL ERROR: 'functions.R' not found. Please ensure it is in the main project directory.")
}
source(here::here("functions.R"))

# if (!file.exists("C:/Users/skgttol/OneDrive - University College London/PhD/PhD_Thesis/02_Data/Template_RScripts/CAGsizing_Flexible/functions.R")) {
#   stop("CRITICAL ERROR: 'functions.R' not found. Please ensure it is in the main project directory.")
# }
# source("C:/Users/skgttol/OneDrive - University College London/PhD/PhD_Thesis/02_Data/Template_RScripts/CAGsizing_Flexible/functions.R")

# 0C. Load Config & Logging ----
if (!file.exists(here::here("config.yml"))) {
  stop("CRITICAL ERROR: 'config.yml' not found.")
}
config <- yaml::read_yaml(here::here("config.yml"))

# Create output directory
settings_output_dir <- here::here(paste0("2a_settings_generation"))
dir.create(settings_output_dir, recursive = TRUE, showWarnings = FALSE)

# Open Log
try(logr::log_close(), silent = TRUE)
log_path <- file.path(settings_output_dir, "settings_log.log")
lf <- logr::log_open(log_path, show_notes = FALSE, logdir = FALSE)

logr::log_print("--- SCRIPT 00: SETTINGS GENERATION INITIALIZED ---", console = TRUE)

# --- VALIDATION ---
# Check inputs, but FALSE for check_external_files (since we haven't made it yet)
validate_config(config, check_external_files = FALSE)


#=============================================================================#
# PART 1: LOAD PLATEMAP & FSA METADATA
#=============================================================================#
logr::log_print("\n--- Loading Platemap and FSA Metadata ---", console = TRUE)

# 1. Load Platemap
platemap_raw <- load_platemap(config)

# 2. Load FSA Metadata
fsa_metadata <- get_fsa_metadata(config$paths$fsa_folder, config)

# 3. Join Data
# Prepare for join by standardizing columns
fsa_df_for_join <- fsa_metadata %>%
  dplyr::select(plate = plate_from_filename, well = well_from_filename, fsa_filename)

platemap_raw_lower <- platemap_raw %>%
  dplyr::mutate(plate = tolower(plate), well = tolower(well))

fsa_df_for_join_lower <- fsa_df_for_join %>%
  dplyr::mutate(plate = tolower(plate), well = tolower(well))

# Join
platemap_data <- dplyr::left_join(platemap_raw_lower, fsa_df_for_join_lower, by = c("plate", "well"))

# Check for unmatched samples
unmatched_samples <- dplyr::filter(platemap_data, is.na(fsa_filename))
if (nrow(unmatched_samples) > 0) {
  logr::log_print("WARNING: The following samples from the platemap did not have a matching .fsa file:")
  logr::log_print(unmatched_samples)
} else {
  logr::log_print("All platemap samples successfully matched to FSA files.")
}


#=============================================================================#
# PART 2: GENERATE SETTINGS FILE
#=============================================================================#
logr::log_print("\n--- Generating Settings File ---", console = TRUE)

settings_file <- generate_custom_settings_file(
  platemap_data = platemap_data,
  config = config,
  output_dir = settings_output_dir
)

if (!is.null(settings_file)) {
  logr::log_print(
    c("*****************************************************************",
      "** SETTINGS FILE GENERATED SUCCESSFULLY **",
      "** **",
      "** Location: **",
      paste("** ", settings_file),
      "** **",
      "** NEXT STEPS: **",
      "** 1. Load this file into your fragment analysis program (e.g. Romeo). **",
      "** 2. Run the analysis. **",
      "** 3. Save the output Excel file to the location specified in config.yml. **",
      "** 4. Run '01_process_data.R'. **",
      "*****************************************************************"),
    console = TRUE, hide_notes = TRUE
  )
} else {
  logr::log_print("ERROR: Settings file generation failed. Check log for details.", console = TRUE)
}

logr::log_close()
