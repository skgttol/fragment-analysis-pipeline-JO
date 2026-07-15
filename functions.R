#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#                                                                             #
#             Advanced Fragment Analysis Workflow - Functions                 #
#                                                                             #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#
# DESCRIPTION:
# This script contains all reusable functions for the fragment analysis
# pipeline (01_process_data.R, 02_generate_plots.R, 03_correlation_analysis.R).
#
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#=============================================================================#
# SECTION 0: HELPER & UTILITY FUNCTIONS
#=============================================================================#
# --- 1. Environment Setup ---####
setup_environment <- function() {
  # Bioconductor packages notoriously overwrite basic dplyr functions. 
  # This locks them in so your data wrangling never crashes.
  if (!requireNamespace("conflicted", quietly = TRUE)) install.packages("conflicted")
  library(conflicted)
# --- 1. Core Data Wrangling (Always prefer Tidyverse) ---
conflicts_prefer(
  dplyr::filter,
  dplyr::select,
  dplyr::lag,
  dplyr::collapse,
  dplyr::combine,
  tidyr::expand,
  tidyr::pack,
  tidyr::unpack
)

# --- 2. Statistical Modeling (Always prefer lmerTest & lme4) ---
conflicts_prefer(
  lmerTest::lmer,
  lmerTest::step,
  lme4::lmList
)

# --- 3. Visualization & Clustering (Tailored for Script 4 & 5) ---
conflicts_prefer(
  dendextend::rotate,     # Protects your Ward clustering dendrograms
  dendextend::labels,
  cowplot::get_legend,    # Protects your publication figure grids
  cowplot::get_title,
  cowplot::align_plots,
  cowplot::stamp,
  patchwork::area
)

# --- 4. Quietly silence benign dataset/utility overlaps ---
conflicts_prefer(
  janitor::chisq.test,
  janitor::fisher.test,
  base::intersect,
  base::setdiff,
  .quiet = TRUE
)
}


`%||%` <- function(lhs, rhs) if (!is.null(lhs)) lhs else rhs

#' Get All Grouping Variables
#'
#' @description Safely parses the config to get a vector of all unique grouping
#' variables, explicitly ignoring any that are set to 'null' or are NULL.
#'
#' @param config The loaded config.yml list.
#' @return A character vector of all valid, non-null grouping variable names.
get_grouping_vars <- function(config) {
  cfg_vars <- config$key_variables
  
  if (is.null(cfg_vars$primary_group_var)) {
    stop("Config Error: 'primary_group_var' must be defined in 'key_variables'.")
  }
  
  all_vars <- list(
    cfg_vars$primary_group_var,
    cfg_vars$secondary_group_var,
    cfg_vars$optional_grouping_var
  )
  
  valid_vars <- unique(unlist(all_vars))
  valid_vars <- valid_vars[!is.null(valid_vars) & valid_vars != 'null']
  
  return(valid_vars)
}

#' Apply Custom Factor Levels
#'
#' @description (MODIFIED) This function now performs two steps:
#' 1. (NEW) Filters the data. If 'factor_levels' are defined in the config
#'    for a column, it is treated as an INCLUSION LIST. Any rows with
#'    values *not* in that list are removed.
#' 2. (EXISTING) Converts character columns to factors, applying
#'    custom levels from the config file if they are provided.
#'
#' @param data The data frame to modify (e.g., processed_data).
#' @param config The loaded config.yml list.
#' @param all_grouping_vars A character vector of all grouping columns.
#' @return The filtered and factored data frame.
apply_factor_levels <- function(data, config, all_grouping_vars) {
  cfg_vars <- config$key_variables
  
  # (This part is unchanged)
  cols_to_factor_list <- list(
    all_grouping_vars,
    cfg_vars$repeated_measure_var,
    cfg_vars$optional_random_effect,
    cfg_vars$optional_crossed_effect # <-- (NEW) Add crossed effect here
  )
  cols_to_factor <- unique(unlist(cols_to_factor_list))
  cols_to_factor <- cols_to_factor[!is.null(cols_to_factor) & cols_to_factor != 'null']
  cols_to_factor <- cols_to_factor[cols_to_factor %in% colnames(data)]
  
  factor_level_definitions <- config$factor_levels
  data_filtered <- data # Start with the original data
  
  # --- (NEW) PART 1: FILTERING BASED ON factor_levels ---
  logr::log_print("Applying factor_levels inclusion filter (Script 01 - Part 2C)...")
  
  if (!is.null(factor_level_definitions)) {
    for (col_name in names(factor_level_definitions)) {
      
      custom_levels <- factor_level_definitions[[col_name]]
      
      # Check if the column exists in the data and the levels are defined
      if (col_name %in% colnames(data_filtered) && !is.null(custom_levels) && all(custom_levels != "null")) {
        
        logr::log_print(paste0("...Filtering '", col_name, "' to ONLY include levels specified in config.yml."))
        
        rows_before <- nrow(data_filtered)
        
        # Ensure column is character for comparison
        data_filtered <- data_filtered %>%
          dplyr::mutate(!!rlang::sym(col_name) := as.character(!!rlang::sym(col_name))) 
        
        # Apply the filter
        data_filtered <- data_filtered %>%
          dplyr::filter(!!rlang::sym(col_name) %in% custom_levels)
        
        rows_after <- nrow(data_filtered)
        logr::log_print(paste0("...Removed ", (rows_before - rows_after), " rows not matching levels for '", col_name, "'."))
        
        # (NEW) Add a check for samples that were *in* the list but not in the data
        data_levels_found <- unique(data_filtered[[col_name]])
        levels_not_found <- setdiff(custom_levels, data_levels_found)
        if (length(levels_not_found) > 0) {
          logr::log_print(paste0("...WARNING: The following levels for '", col_name, "' were in the config but not found in the data: ", 
                                 paste(levels_not_found, collapse = ", ")))
        }
        
      } else if (col_name %in% colnames(data_filtered)) {
        logr::log_print(paste0("...No levels defined for '", col_name, "'. All existing values will be kept."))
      }
    }
  }
  # --- END NEW PART 1 ---
  
  
  # --- (EXISTING) PART 2: CONVERTING TO FACTORS ---
  logr::log_print("Converting key grouping variables to factors...")
  
  data_factored <- data_filtered %>% # (MODIFIED: Use data_filtered)
    dplyr::mutate(dplyr::across(dplyr::all_of(cols_to_factor), ~{
      col_name <- dplyr::cur_column()
      
      if (!is.null(factor_level_definitions) && col_name %in% names(factor_level_definitions)) {
        custom_levels <- factor_level_definitions[[col_name]]
        
        if (!is.null(custom_levels) && all(custom_levels != "null")) {
          logr::log_print(paste0("...applying custom factor order for '", col_name, "'."))
          # This is now safe because all non-listed values have been removed
          return(factor(., levels = custom_levels))
        }
      }
      
      logr::log_print(paste0("...applying default alphabetical order for '", col_name, "'."))
      return(factor(.))
    }))
  
  return(data_factored)
}


#' Export All Results to a Multi-Sheet Excel File
#'
#' @description Writes a named list of data frames to a single .xlsx file.
#'
#' @param export_list A named list of data frames.
#' @param file_path The full file path for the output Excel file.
export_to_excel <- function(export_list, file_path) {
  valid_list <- export_list[!sapply(export_list, is.null)]
  valid_list <- valid_list[sapply(valid_list, nrow) > 0]
  
  if (length(valid_list) == 0) {
    message(paste("No data to export to:", file_path))
    return(invisible(NULL))
  }
  
  openxlsx::write.xlsx(valid_list, file = file_path, rowNames = FALSE)
  message(paste("Data successfully saved to:", file_path))
}

#' Get a list of files from a config path
#'
#' @description Handles a config path that can be a single file, 
#' a folder, or a list of files.
#' @param path_config The value from config$...$path
#' @param file_pattern A regex pattern to filter by if path_config is a folder.
#' @return A character vector of file paths.
get_file_paths <- function(path_config, file_pattern) {
  if (is.null(path_config) || all(path_config == 'null')) {
    stop("Path is not defined in config.")
  }
  
  file_list <- c()
  
  # Case 1: It's a list from YAML. Unlist it to a character vector.
  if (is.list(path_config)) {
    file_list <- unlist(path_config)
    
    # Case 2: It's a character vector (could be length 1 or more)
  } else if (is.character(path_config)) {
    
    # Sub-case 2a: It's a single string, which could be a file OR a folder
    if (length(path_config) == 1) {
      if (dir.exists(path_config)) {
        # It's a folder
        file_list <- list.files(
          path = path_config,
          pattern = file_pattern,
          full.names = TRUE,
          recursive = FALSE
        )
        if (length(file_list) == 0) {
          warning(paste("No files matching pattern", file_pattern, "found in folder:", path_config))
        }
      } else if (file.exists(path_config)) {
        # It's a single file
        file_list <- c(path_config)
      } else {
        stop(paste("Path not found:", path_config))
      }
    } else {
      # Sub-case 2b: It's already a character vector of multiple paths
      # (This is the case that was failing)
      file_list <- path_config
    }
  } else {
    stop("Path config is not in a recognized format (must be a single path, a folder, a list, or a vector of paths).")
  }
  
  # Final check: do all files in the list exist?
  if (length(file_list) > 0) {
    files_exist <- file.exists(file_list)
    if (!all(files_exist)) {
      stop(paste("The following files listed in the config path do not exist:\n", 
                 paste(file_list[!files_exist], collapse = "\n")))
    }
  }
  return(file_list)
}

#' Apply a Simple or Composite Parsing Rule
#'
#' @description Applies one or more parsing rules and combines them.
#' Provides backward compatibility for single-rule format.
#'
#' @param strings A character vector of filenames.
#' @param rule_config The config list for this rule (e.g., config$fsa...$plate_number)
#' @return A character vector of the extracted & combined values.
apply_composite_parsing_rule <- function(strings, rule_config) {
  # Check for old format (has 'method' at the top level)
  if (!is.null(rule_config$method)) {
    logr::log_print("...using single parsing rule (legacy format).")
    return(apply_parsing_rule(strings, rule_config))
  }
  
  # New composite format
  separator <- rule_config$separator %||% "_"
  rules_list <- rule_config[!names(rule_config) %in% "separator"]
  
  if (length(rules_list) == 0) {
    stop("Composite parsing rule defined, but no sub-rules (e.g., part_1, part_2) found.")
  }
  
  logr::log_print(paste("...using composite parsing rule with", length(rules_list), "parts."))
  
  # Run each sub-rule
  extracted_parts <- purrr::map(rules_list, ~ apply_parsing_rule(strings, .x))
  
  # Combine the results
  combined_strings <- do.call(paste, c(extracted_parts, sep = separator))
  
  return(combined_strings)
}


#=============================================================================#
# SECTION 1: DATA LOADING AND PRE-PROCESSING (SCRIPT 01)
#=============================================================================#

#' Apply a Config-Based Parsing Rule
#'
#' @param strings A character vector of filenames.
#' @param rule A list containing method and parameters.
#' @return A character vector of the extracted values.
apply_parsing_rule <- function(strings, rule) {
  if (is.null(rule) || is.null(rule$method)) {
    logr::log_print("WARNING: A parsing rule is null or missing a 'method'. Returning NA.", console = TRUE)
    return(NA_character_)
  }
  
  if (rule$method == "split") {
    if (is.null(rule$delimiter) || is.null(rule$position)) {
      stop("For 'split' method, 'delimiter' and 'position' must be defined in config.yml.")
    }
    return(stringr::str_split_i(strings, rule$delimiter, rule$position))
    
  } else if (rule$method == "regex") {
    if (is.null(rule$pattern)) {
      stop("For 'regex' method, 'pattern' must be defined in config.yml.")
    }
    return(stringr::str_extract(strings, rule$pattern))
    
  } else {
    stop(paste("Invalid parsing method specified:", rule$method, ". Must be 'split' or 'regex'."))
  }
}

#' Find and Parse Metadata from all FSA Files
#'
#' @param fsa_base_folder The path to the top-level folder.
#' @param config The loaded configuration list.
#' @return A tidy data frame with metadata for each .fsa file.
get_fsa_metadata <- function(fsa_base_folder, config) {
  if (is.null(fsa_base_folder) || fsa_base_folder == 'null') {
    stop("Config Error: 'paths:fsa_folder' is not defined.")
  }
  
  all_fsa_paths <- list.files(
    path = fsa_base_folder,
    pattern = "\\.fsa$",
    recursive = TRUE,
    full.names = TRUE
  )
  
  if (length(all_fsa_paths) == 0) {
    stop("No .fsa files were found in: ", fsa_base_folder)
  }
  
  parsing_rules <- config$fsa_filename_parsing
  if (is.null(parsing_rules)) {
    stop("Config Error: `fsa_filename_parsing` section is missing from config.yml.")
  }
  
  fsa_metadata_df <- tibble::tibble(full_path = all_fsa_paths) %>%
    dplyr::mutate(
      fsa_filename = basename(full_path),
      
      # --- MODIFIED: Use composite parser ---
      plate_extracted = apply_composite_parsing_rule(fsa_filename, parsing_rules$plate_number),
      well_extracted = apply_composite_parsing_rule(fsa_filename, parsing_rules$well_name),
      
      # --- (NEW) Apply the consistent cleaning rule ---
      plate_from_filename = clean_plate_name(plate_extracted),
      well_from_filename = tolower(well_extracted) # Wells just need to be lowercase
    ) %>%
    dplyr::select(-plate_extracted, -well_extracted) # Drop the intermediate columns
  
  logr::log_print(paste("Found and parsed metadata for", nrow(fsa_metadata_df), ".fsa files."))
  return(fsa_metadata_df)
}

#' Load a Tidy Platemap from CSV or Excel
#'
#' @description (MODIFIED) Preserves original plate/well names.
#' @param path The file path to the CSV or Excel file.
#' @return A single tidy data frame with platemap info.
load_tidy_platemap <- function(path) {
  if (!file.exists(path)) stop("Tidy platemap file not found: ", path)
  
  if (tools::file_ext(path) == "csv") {
    platemap_data <- readr::read_csv(path, show_col_types = FALSE)
  } else {
    platemap_data <- readxl::read_excel(path)
  }
  
  required_cols <- c("plate", "well", "sample_name")
  if (!all(required_cols %in% colnames(platemap_data))) {
    stop("Tidy platemap must contain columns: ", paste(required_cols, collapse = ", "))
  }
  
  platemap_tidy <- platemap_data %>%
    # --- (MODIFIED) Preserve raw names before cleaning ---
    dplyr::mutate(
      plate_raw = as.character(plate), # Preserve original plate name
      well_raw = as.character(well),   # Preserve original well name
      plate = clean_plate_name(plate),
      well = tolower(as.character(well))
    ) %>%
    # Select all original columns plus the new raw ones
    dplyr::select(plate, well, sample_name, plate_raw, well_raw, dplyr::everything())
  
  logr::log_print(paste("Successfully loaded", nrow(platemap_tidy), "samples from tidy file:", basename(path)))
  return(platemap_tidy)
}

#' Load and Tidy a Stacked Excel Platemap
#'
#' @description (MODIFIED) Preserves original plate/well names.
#' @param path The file path to the single Excel file.
#' @param config The loaded config.yml list.
#' @return A single tidy data frame with platemap info.
load_stacked_platemap_excel <- function(path, config) {
  if (!file.exists(path)) stop("Stacked platemap not found: ", path)
  
  identifier <- config$platemap$stacked_excel_identifier %||% "^Plate"
  logr::log_print(paste("...using stacked identifier:", identifier))
  
  raw_sheet <- readxl::read_excel(path, col_names = FALSE)
  plate_starts <- which(stringr::str_detect(raw_sheet[[1]], identifier))
  
  if (length(plate_starts) == 0) {
    stop(paste("No rows matching '", identifier, "' found in Excel file:", path))
  }
  
  platemap_tidy <- purrr::map_df(plate_starts, function(i) {
    name_raw <- raw_sheet[[1]][i] # <-- CAPTURE RAW NAME HERE
    grid <- raw_sheet[(i + 1):(i + 8), 1:12]
    name_clean <- clean_plate_name(name_raw) 
    
    grid %>%
      tibble::as_tibble(.name_repair = "minimal") %>%
      dplyr::rename_with(~ as.character(1:12)) %>%
      dplyr::mutate(Row = LETTERS[1:8]) %>%
      tidyr::pivot_longer(cols = -Row, names_to = "Column", values_to = "sample_name") %>%
      dplyr::mutate(
        Column = sprintf("%02d", as.numeric(Column)),
        well_raw = paste0(Row, Column), # <-- CAPTURE RAW WELL HERE
        well = tolower(well_raw),
        plate_raw = name_raw,           # <-- ADD RAW PLATE NAME
        plate = name_clean              # <-- ADD CLEAN PLATE NAME
      ) %>%
      dplyr::select(plate, well, sample_name, plate_raw, well_raw) %>% # <-- SELECT NEW COLUMNS
      dplyr::filter(!is.na(sample_name) & sample_name != "") # Filter for non-empty
  })
  
  logr::log_print(paste("Successfully loaded and tidied", length(plate_starts), "plates from stacked Excel:", basename(path)))
  return(platemap_tidy)
}

#' Load Platemap (Wrapper Function)
#'
#' @description (MODIFIED) Manages loading of raw plate/well names.
#' @param config The loaded config.yml list.
#' @return A single tidy data frame with platemap info.
load_platemap <- function(config) {
  platemap_cfg <- config$platemap
  if (is.null(platemap_cfg) || is.null(platemap_cfg$path) || is.null(platemap_cfg$format)) {
    stop("Config Error: 'platemap' section with 'path' and 'format' must be defined.")
  }
  
  file_paths <- get_file_paths(platemap_cfg$path, "\\.(csv|xlsx|xls)$")
  if (length(file_paths) == 0) {
    stop("No platemap files found at specified path.")
  }
  logr::log_print(paste("Found", length(file_paths), "platemap file(s)."))
  
  # --- (MODIFIED) Loader functions now need to accept config ---
  if (platemap_cfg$format == "stacked_excel") {
    platemap_raw <- purrr::map_df(file_paths, ~ load_stacked_platemap_excel(.x, config))
  } else if (platemap_cfg$format == "tidy_csv") {
    platemap_raw <- purrr::map_df(file_paths, ~ load_tidy_platemap(.x))
  } else if (platemap_cfg$format == "tidy_excel") {
    platemap_raw <- purrr::map_df(file_paths, ~ load_tidy_platemap(.x))
  } else {
    stop("Invalid platemap format in config.yml. Must be 'stacked_excel', 'tidy_csv', or 'tidy_excel'.")
  }
  # --- END MODIFICATION ---
  
  platemap_raw_cleaned <- platemap_raw %>%
    dplyr::distinct() 
  
  logr::log_print(paste("Total platemap samples loaded and combined:", nrow(platemap_raw_cleaned)))
  return(platemap_raw_cleaned)
}

#' Generate Settings File for Custom Analysis (e.g. Romeo/Juliet)
#'
#' @description Creates an Excel file with specific columns required for external tools.
#' Includes logic to split sample names and supports "custom" baseline grouping.
#'
#' @param platemap_data The processed platemap dataframe.
#' @param config The config list.
#' @param output_dir Directory to save the file.
generate_custom_settings_file <- function(platemap_data, config, output_dir) {
  
  # 1. Filter valid data
  settings_data <- platemap_data %>%
    dplyr::filter(!is.na(fsa_filename))
  
  if (nrow(settings_data) == 0) {
    stop("No matched FSA files found in the platemap data. Cannot generate settings file.")
  }
  
  # 2. Parse Sample Names Safely (prevents duplicate column name collisions)
  sample_cols <- config$sample_name_columns
  if (is.null(sample_cols)) stop("Config Error: 'sample_name_columns' must be defined.")
  
  settings_data_parsed <- settings_data %>%
    tidyr::separate(
      col = sample_name, 
      into = sample_cols, 
      sep = "_", 
      remove = FALSE, 
      extra = "warn", 
      fill = "right"
    )
  
  # 3. Define variables
  # Note: using rlang::`%||%` or base R fallback if %||% is not loaded
  `%||%` <- rlang::`%||%`
  
  time_var <- config$key_variables$time_variable
  primary_var <- config$key_variables$primary_group_var
  secondary_var <- config$key_variables$secondary_group_var
  batch_var <- config$key_variables$optional_crossed_effect %||% 
    config$key_variables$optional_grouping_var
  
  if (!time_var %in% colnames(settings_data_parsed)) {
    logr::log_print(paste("WARNING: Time variable", time_var, "not found. Skipping settings file."))
    return(NULL)
  }
  
  # 4. Logic Flags
  has_secondary <- !is.null(secondary_var) && secondary_var != "null" && secondary_var %in% colnames(settings_data_parsed)
  has_batch <- !is.null(batch_var) && batch_var != "null" && batch_var %in% colnames(settings_data_parsed)
  
  baseline_method <- config$parameters$baseline_grouping_level %||% "group"
  
  # 5. Create Group Column Logic (Replaced rlang/case_when with safe conditional logic)
  tryCatch({
    
    # Pre-extract the time column
    settings_prepared <- settings_data_parsed %>%
      dplyr::mutate(
        time_col = as.numeric(stringr::str_extract(as.character(.data[[time_var]]), "\\d+\\.?\\d*"))
      )
    
    # Apply Grouping Logic
    if (baseline_method == "global_control") {
      settings_prepared$group_col <- config$parameters$global_control_group$level
      
    } else if (baseline_method == "genotype") {
      settings_prepared$group_col <- as.character(settings_prepared[[primary_var]])
      
    } else if (baseline_method == "custom") {
      custom_vars <- config$parameters$baseline_group_vars
      if(is.null(custom_vars)) stop("baseline_group_vars must be defined when baseline_grouping_level is 'custom'")
      
      missing <- setdiff(custom_vars, colnames(settings_prepared))
      if(length(missing) > 0) stop(paste("Missing custom baseline vars:", paste(missing, collapse=", ")))
      
      # tidyr::unite cleanly pastes dynamic column lists together
      settings_prepared <- settings_prepared %>%
        tidyr::unite("group_col", dplyr::all_of(custom_vars), sep = "_", remove = FALSE)
      
    } else {
      # Standard Grouping (Primary + Optional Secondary + Optional Batch)
      group_cols <- primary_var
      if (has_secondary) group_cols <- c(group_cols, secondary_var)
      if (has_batch) group_cols <- c(group_cols, batch_var)
      
      settings_prepared <- settings_prepared %>%
        tidyr::unite("group_col", dplyr::all_of(group_cols), sep = "_", remove = FALSE)
    }
    
    # 6. Map Final Columns (Safe Extraction prevents NULLs from dropping columns)
    cs <- config$custom_settings %||% list()
    
    safe_val <- function(val, default = NA) {
      if (is.null(val)) default else val
    }
    
    final_settings_file <- settings_prepared %>%
      dplyr::transmute(
        sample_name = fsa_filename,
        flank = safe_val(cs$flank),
        repeat_unit = safe_val(cs$repeat_unit),
        correction = safe_val(cs$correction),
        pathogenic_threshold = safe_val(cs$pathogenic_threshold),
        ii_threshold_method = safe_val(cs$ii_threshold_method, "Height"),
        ii_threshold_height = safe_val(cs$ii_threshold_height),
        ii_threshold_area = safe_val(cs$ii_threshold_area),
        target_channel = safe_val(cs$target_channel),
        target_start = safe_val(cs$target_start),
        target_end = safe_val(cs$target_end),
        plot_start = safe_val(cs$plot_start),
        plot_end = safe_val(cs$plot_end),
        control_rpt = safe_val(cs$control_rpt, ""),
        group = group_col,
        group_control_sample = safe_val(cs$group_control_sample_prefix, ""),
        time = time_col,
        exclude = ""
      )
    
    # 7. Add defaults if they are completely missing from the transmute mapping
    defaults <- list(shared_sample_donor = NA, shared_sample_recipient = NA, unique_id = NA)
    for (col in names(defaults)) {
      if (!col %in% colnames(final_settings_file)) {
        final_settings_file[[col]] <- safe_val(cs[[col]], defaults[[col]])
      }
    }
    
    # 8. Save output
    if (!dir.exists(output_dir)) dir.create(output_dir, recursive = TRUE)
    file_name <- paste0(format(Sys.Date(), "%Y-%m-%d"), "_", format(Sys.time(), "%H%M"), "_romeo_settings.xlsx")
    full_path <- file.path(output_dir, file_name)
    
    if (requireNamespace("openxlsx", quietly = TRUE)) {
      openxlsx::write.xlsx(final_settings_file, full_path)
      logr::log_print(paste("Saved custom settings file to:", full_path))
    } else {
      logr::log_print("WARNING: 'openxlsx' package not installed. Skipping Excel export.")
    }
    
    return(full_path)
    
  }, error = function(e) {
    logr::log_print(paste("Error generating custom settings file:", e$message))
    return(NULL)
  })
}

#' Standardize a Plate Name for Joining
#'
#' @description Cleans plate names to ensure robust matching.
#' CHANGE: Now normalizes spaces/hyphens to UNDERSCORES instead of removing them.
#' Example: "Plate 1" -> "plate_1", "1.1_220925" -> "1.1_220925"
#'
#' @param plate_string A character vector of plate names.
#' @return A cleaned, standardized character vector.
clean_plate_name <- function(plate_string) {
  plate_string %>%
    as.character() %>%
    tolower() %>%
    # Remove the word "plate" if present
    stringr::str_remove_all("plate") %>%
    stringr::str_trim() %>%
    # Replace spaces and hyphens with underscores to maintain structure
    stringr::str_replace_all("[ -]", "_") %>%
    # Remove duplicate underscores
    stringr::str_replace_all("_+", "_") %>%
    # Remove leading/trailing underscores
    stringr::str_remove("^_+|_+$")
}

#' Load and Process Data from External Program (Hybrid + Robust Exclusions)
#'
#' @description Loads external Excel files with Hybrid Strategy.
#' 1. FIXES Romeo column names.
#' 2. ALWAYS parses fsa_filename to generate plate/well/plate_raw/well_raw.
#' 3. PATCHES missing sample names using the platemap (if needed) on a file-by-file basis.
#'    - Uses dynamic joining (Filename vs Plate/Well) depending on what's available.
#' 4. APPLIES robust exclusions.
#'
#' @param config The loaded config.yml list.
#' @param platemap_data The tidy platemap data frame.
#' @return A list containing `processed_data`, `excluded_data`, and `pre_filter_data_with_reasons`.
load_and_process_data <- function(config, platemap_data = NULL) {
  
  # --- Get list of files ---
  ext_path_config <- config$paths$external_data_path
  file_paths <- get_file_paths(ext_path_config, "\\.(xlsx|xls)$")
  if (length(file_paths) == 0) stop("No external data files found.")
  
  logr::log_print(paste("Found", length(file_paths), "external data file(s)."))
  
  # --- 1. Load Each File & Process Individually ---
  resp_vars_cfg <- config$response_variables
  external_numeric_cols <- c(as.character(unlist(resp_vars_cfg)), "target_peak_height", "time")
  
  # Extract parsing rules once to use inside the loop
  parsing_rules <- config$fsa_filename_parsing
  
  data_list <- purrr::map(file_paths, ~ {
    current_path <- .x
    logr::log_print(paste("...loading:", basename(current_path)))
    
    # Read header
    headers <- readxl::read_excel(current_path, sheet = "height", n_max = 0) %>% colnames()
    col_types <- rep("guess", length(headers))
    if ("exclude" %in% tolower(headers)) col_types[which(tolower(headers) == "exclude")] <- "text"
    
    data <- readxl::read_excel(current_path, sheet = "height", col_types = col_types) %>%
      janitor::clean_names()
    
    # --- A. COLUMN SWAP FIX (Romeo format) ---
    if ("sample_name" %in% colnames(data) && "name" %in% colnames(data)) {
      first_val <- as.character(data$sample_name[1])
      if (!is.na(first_val) && grepl("\\.fsa$", first_val, ignore.case = TRUE)) {
        data <- data %>%
          dplyr::rename(fsa_filename = sample_name) %>%
          dplyr::rename(sample_name = name)
      }
    }
    
    # --- B. Standardize FSA Column ---
    if (!"fsa_filename" %in% colnames(data)) {
      if (grepl("\\.fsa$", as.character(data[[1]][1]), ignore.case = TRUE)) {
        data <- dplyr::rename(data, fsa_filename = 1)
      }
    }
    
    # Ensure fsa_filename is character
    if("fsa_filename" %in% colnames(data)) {
      data$fsa_filename <- as.character(data$fsa_filename)
    }
    
    # --- C. (NEW LOGIC) ALWAYS PARSE FILENAME FOR PLATE/WELL ---
    if ("fsa_filename" %in% colnames(data) && !is.null(parsing_rules)) {
      data <- data %>%
        dplyr::mutate(
          plate_extracted = apply_composite_parsing_rule(fsa_filename, parsing_rules$plate_number),
          well_extracted  = apply_composite_parsing_rule(fsa_filename, parsing_rules$well_name),
          # Generate clean versions (Using updated clean_plate_name which preserves underscores)
          plate = clean_plate_name(plate_extracted),
          well  = tolower(well_extracted),
          # Generate the RAW versions
          plate_raw = plate_extracted,
          well_raw  = well_extracted
        ) %>%
        dplyr::select(-plate_extracted, -well_extracted)
    }
    
    # --- D. (UPDATED) HYBRID PATCHING (Per File) ---
    needs_patching <- FALSE
    if (!"sample_name" %in% colnames(data)) {
      needs_patching <- TRUE
      data$sample_name <- NA_character_
    } else if (all(is.na(data$sample_name) | trimws(data$sample_name) == "")) {
      needs_patching <- TRUE
    }
    
    if (needs_patching) {
      if (!is.null(platemap_data)) {
        logr::log_print(paste("...file missing sample names. Patching from platemap."))
        
        # --- DYNAMIC JOIN STRATEGY ---
        if ("fsa_filename" %in% colnames(platemap_data)) {
          # Strategy A: Robust Link (Platemap has filenames)
          logr::log_print("......joining by 'fsa_filename'.")
          
          lookup <- platemap_data %>% 
            dplyr::select(fsa_filename, plate_sample_name = sample_name) %>%
            dplyr::mutate(fsa_filename = tolower(basename(fsa_filename)))
          
          data <- data %>%
            dplyr::mutate(fsa_lower = tolower(basename(fsa_filename))) %>%
            dplyr::left_join(lookup, by = c("fsa_lower" = "fsa_filename")) %>%
            dplyr::select(-fsa_lower)
          
        } else {
          # Strategy B: Position Link (Platemap only has Plate/Well)
          logr::log_print("......joining by 'plate' and 'well'.")
          
          if(!all(c("plate", "well") %in% colnames(data))) {
            logr::log_print("WARNING: Cannot patch by Plate/Well because filename parsing failed (columns missing). Check config rules.")
          } else {
            # Prepare Lookup: Force character types and use updated cleaner
            lookup <- platemap_data %>%
              dplyr::mutate(
                plate = clean_plate_name(plate), # Re-run clean to ensure consistency
                well = tolower(as.character(well))
              ) %>%
              dplyr::select(plate, well, plate_sample_name = sample_name)
            
            # Perform Join
            data <- data %>%
              dplyr::left_join(lookup, by = c("plate", "well"))
          }
        }
        
        # Apply the patch
        if("plate_sample_name" %in% colnames(data)) {
          n_patched <- sum(!is.na(data$plate_sample_name))
          logr::log_print(paste0("......patched ", n_patched, " sample names."))
          
          data <- data %>%
            dplyr::mutate(sample_name = dplyr::coalesce(sample_name, plate_sample_name)) %>%
            dplyr::select(-plate_sample_name)
        }
        
      } else {
        logr::log_print("WARNING: Sample names missing and NO platemap provided for patching.")
      }
    }
    
    # Force numeric types safely
    data %>%
      dplyr::mutate(dplyr::across(any_of(external_numeric_cols), as.numeric))
  })
  
  # Combine all data
  pre_filter_data <- dplyr::bind_rows(data_list)
  
  # --- 3. Rename Response Variables (Config Mapping) ---
  response_vars_rename <- setNames(as.character(resp_vars_cfg), names(resp_vars_cfg))
  existing_cols <- colnames(pre_filter_data)
  
  for (new_name in names(response_vars_rename)) {
    old_name <- response_vars_rename[[new_name]]
    clean_old <- janitor::make_clean_names(old_name)
    if (clean_old %in% existing_cols) {
      pre_filter_data <- dplyr::rename(pre_filter_data, !!new_name := !!clean_old)
    }
  }
  
  # --- 4. Split Sample Name ---
  sample_cols <- config$sample_name_columns
  if (is.null(sample_cols)) stop("Config Error: 'sample_name_columns' must be defined.")
  
  pre_filter_data <- pre_filter_data %>%
    tidyr::separate(sample_name, into = sample_cols, sep = "_", remove = FALSE, fill = "right")
  
  # --- 5. Tag Rows for Exclusion ---
  peak_thresh <- config$parameters$peak_height_threshold %||% 0
  height_col <- if("target_peak_height" %in% colnames(pre_filter_data)) "target_peak_height" else NULL
  blank_part_regex <- "^_|_$|__|_\\s*_"
  
  data_with_reasons <- pre_filter_data %>%
    dplyr::mutate(
      exclusion_reason = dplyr::case_when(
        is.na(sample_name) ~ "Missing sample_name (Patch failed)",
        trimws(tolower(sample_name)) %in% c("", "blank") ~ "Blank or 'Blank' name",
        stringr::str_detect(sample_name, blank_part_regex) ~ "Malformed name (e.g., __)",
        !is.null(height_col) & .[[height_col]] < peak_thresh ~ paste("Peak height <", peak_thresh),
        "exclude" %in% names(.) & !is.na(exclude) & tolower(exclude) != "n" ~ "Marked 'exclude' in file",
        TRUE ~ NA_character_
      )
    )
  
  processed_data_temp <- data_with_reasons %>% dplyr::filter(is.na(exclusion_reason))
  excluded_data_1 <- data_with_reasons %>% dplyr::filter(!is.na(exclusion_reason))
  
  # --- 6. Apply Config Exclusions (Round 2) ---
  config_exclusions <- config$exclusions
  processed_data_final <- processed_data_temp
  excluded_data_2 <- data.frame()
  
  if (!is.null(config_exclusions) && all(config_exclusions != 'null')) {
    logr::log_print("Applying exclusion rules from config...")
    
    processed_data_temp <- processed_data_temp %>% dplyr::mutate(temp_row_id = dplyr::row_number())
    ids_to_exclude <- c()
    
    for (rule in config_exclusions) {
      clean_rule <- Filter(Negate(is.null), rule)
      if (length(clean_rule) == 0) next
      
      tryCatch({
        rule_df <- as.data.frame(clean_rule, stringsAsFactors = FALSE) %>%
          dplyr::mutate(dplyr::across(everything(), as.character))
        
        matched_ids <- processed_data_temp %>%
          dplyr::mutate(dplyr::across(names(rule_df), as.character)) %>%
          dplyr::inner_join(rule_df, by = names(rule_df)) %>%
          dplyr::pull(temp_row_id)
        
        ids_to_exclude <- c(ids_to_exclude, matched_ids)
      }, error = function(e) {
        logr::log_print(paste("WARNING: Skipped malformed exclusion rule:", paste(names(rule), collapse=", ")))
      })
    }
    
    ids_to_exclude <- unique(ids_to_exclude)
    
    if (length(ids_to_exclude) > 0) {
      config_excluded_rows <- processed_data_temp %>%
        dplyr::filter(temp_row_id %in% ids_to_exclude) %>%
        dplyr::select(-temp_row_id) %>%
        dplyr::mutate(exclusion_reason = "Excluded by config.yml")
      
      processed_data_final <- processed_data_temp %>%
        dplyr::filter(!temp_row_id %in% ids_to_exclude) %>%
        dplyr::select(-temp_row_id)
      
      excluded_data_2 <- config_excluded_rows
    } else {
      processed_data_final <- processed_data_temp %>% dplyr::select(-temp_row_id)
    }
  }
  
  # --- 7. Finalize ---
  excluded_data <- dplyr::bind_rows(excluded_data_1, excluded_data_2)
  
  time_var <- config$key_variables$time_variable
  if(time_var != "null" && time_var %in% colnames(processed_data_final)) {
    processed_data_final[[time_var]] <- as.numeric(gsub("[^0-9\\.]", "", processed_data_final[[time_var]]))
  }
  
  logr::log_print(paste("External data loaded. Valid rows:", nrow(processed_data_final), "| Excluded:", nrow(excluded_data)))
  
  return(list(
    processed_data = processed_data_final,
    excluded_data = excluded_data,
    pre_filter_data_with_reasons = data_with_reasons
  ))
}


#' Normalize Data to Baseline (Time-Aligned)
#'
#' @description Calculates baselines and creates 'mode_change' and 'instability_index_change'.
#' UPDATED: Automatically aligns timepoints so every group starts at Day 0.
#' Preserves the original time values in a new column '{time_var}_raw'.
#' NOW INCLUDES: Change in Instability Index.
#'
#' @param processed_data The clean data frame.
#' @param config The loaded config.yml list.
#' @param all_grouping_vars A character vector of default grouping columns.
#' @return A list containing `normalized_data` and `baselines_table`.
normalize_data <- function(processed_data, config, all_grouping_vars) {
  
  cfg_vars <- config$key_variables
  time_var <- cfg_vars$time_variable
  baseline_method <- config$parameters$baseline_grouping_level %||% "group"
  
  baselines_table <- data.frame() 
  
  # Check if instability_index is present in the data
  has_ii <- "instability_index" %in% colnames(processed_data)
  has_expi <- "expansion_index" %in% colnames(processed_data)
  has_contri <- "contraction_index" %in% colnames(processed_data)
  
  # --- 1. Determine Grouping Variables ---
  if (baseline_method == "global_control") {
    logr::log_print("Using 'global_control' baseline method.")
    
    control_group_cfg <- config$parameters$global_control_group
    if (is.null(control_group_cfg) || is.null(control_group_cfg$variable) || is.null(control_group_cfg$level)) {
      stop("For 'global_control' method, 'global_control_group' with 'variable' and 'level' must be defined.")
    }
    # Use standard grouping vars for time alignment logic
    baseline_grouping_vars <- all_grouping_vars
    
  } else if (baseline_method == "genotype") {
    baseline_grouping_vars <- cfg_vars$primary_group_var
    logr::log_print(paste("Using 'genotype' method: Baselines/Time averaged per", cfg_vars$primary_group_var))
    
  } else if (baseline_method == "custom") {
    custom_vars <- config$parameters$baseline_group_vars
    if (is.null(custom_vars)) stop("Baseline method 'custom' requires 'baseline_group_vars' in config.")
    baseline_grouping_vars <- unlist(custom_vars)
    logr::log_print(paste("Using 'custom' method. Grouped by:", paste(baseline_grouping_vars, collapse = ", ")))
    
  } else {
    # Default: Group
    logr::log_print("Using 'group' method.")
    baseline_grouping_vars <- all_grouping_vars 
    re_cross <- cfg_vars$optional_crossed_effect
    if (!is.null(re_cross) && re_cross != 'null' && re_cross %in% colnames(processed_data)) {
      baseline_grouping_vars <- unique(c(baseline_grouping_vars, re_cross))
    }
  }
  
  # --- 2. Calculate Per-Group Start Time ---
  # Find the minimum time recorded for EACH specific group
  group_stats <- processed_data %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(baseline_grouping_vars))) %>%
    dplyr::mutate(
      group_start_time = min(!!rlang::sym(time_var), na.rm = TRUE)
    ) %>%
    dplyr::ungroup() 
  
  # --- 3. Align Time & Retain Raw ---
  logr::log_print("Aligning timelines: Subtracting group-specific start time from time variable.")
  
  normalized_data_aligned <- group_stats %>%
    dplyr::mutate(
      # Save the raw time first (e.g. "days_raw")
      !!paste0(time_var, "_raw") := !!rlang::sym(time_var),
      # Overwrite main time with aligned version (starts at 0)
      !!rlang::sym(time_var) := !!rlang::sym(time_var) - group_start_time
    )
  
  # --- 4. Calculate Baseline Values ---
  
  if (baseline_method == "global_control") {
    # Special calculation for global control
    control_var <- config$parameters$global_control_group$variable
    control_level <- config$parameters$global_control_group$level
    
    # Get subset for control group at Time 0
    baseline_subset <- normalized_data_aligned %>%
      dplyr::filter(!!rlang::sym(control_var) == control_level) %>%
      dplyr::filter(!!rlang::sym(time_var) == 0)
    
    global_baseline_val <- mean(baseline_subset$mode, na.rm = TRUE)
    global_baseline_ii <- if(has_ii) mean(baseline_subset$instability_index, na.rm = TRUE) else NA_real_
    
    logr::log_print(paste0("Global Control Baseline (Mode): ", round(global_baseline_val, 2)))
    
    # Apply to everyone
    normalized_data <- normalized_data_aligned %>%
      dplyr::mutate(
        final_baseline = global_baseline_val,
        mode_change = mode - final_baseline
      )
    
    if(has_ii) {
      normalized_data <- normalized_data %>%
        dplyr::mutate(
          final_baseline_ii = global_baseline_ii,
          instability_index_change = instability_index - final_baseline_ii
        )
    }
    
    if(has_expi) {
      normalized_data <- normalized_data %>%
        dplyr::mutate(
          final_baseline_expi = global_baseline_expi,
          expansion_index_change = expansion_index - final_baseline_expi
        )
    }
    
    if(has_contri) {
      normalized_data <- normalized_data %>%
        dplyr::mutate(
          final_baseline_contri = global_baseline_contri,
          contraction_index_change = contraction_index - final_baseline_contri
        )
    }
    
    
    baselines_table <- data.frame(
      Type = "Global Control",
      Group = control_level,
      Baseline_Value = global_baseline_val,
      Baseline_II = if(has_ii) global_baseline_ii else NA,
      Baseline_ExpI = if(has_expi) global_baseline_expi else NA,
      Baseline_ContrI = if(has_contri) global_baseline_contri else NA
    )
    
  } else {
    # Standard Dynamic Calculation (Group/Genotype/Custom)
    baselines <- normalized_data_aligned %>%
      dplyr::group_by(dplyr::across(dplyr::all_of(baseline_grouping_vars))) %>%
      dplyr::filter(!!rlang::sym(time_var) == 0) %>%
      dplyr::summarise(
        final_baseline = mean(mode, na.rm = TRUE),
        # Calculate II baseline if column exists
        final_baseline_ii = if(has_ii) mean(instability_index, na.rm = TRUE) else NA_real_,
        final_baseline_expi = if(has_expi) mean(expansion_index, na.rm = TRUE) else NA_real_,
        final_baseline_contri = if(has_contri) mean(contraction_index, na.rm = TRUE) else NA_real_,
        original_start_day = dplyr::first(group_start_time),
        n_baseline_samples = n(),
        .groups = "drop"
      )
    
    baselines_table <- baselines
    
    # Join back
    normalized_data <- normalized_data_aligned %>%
      dplyr::left_join(baselines, by = baseline_grouping_vars) %>%
      dplyr::mutate(mode_change = mode - final_baseline)
    
    # Calculate II Change if available
    if(has_ii) {
      normalized_data <- normalized_data %>%
        dplyr::mutate(instability_index_change = instability_index - final_baseline_ii)
    }
    
    if(has_expi) {
      normalized_data <- normalized_data %>%
        dplyr::mutate(expansion_index_change = expansion_index - final_baseline_expi)
    }
    
    if(has_contri) {
      normalized_data <- normalized_data %>%
        dplyr::mutate(contraction_index_change = contraction_index - final_baseline_contri)
    }
  }
  
  logr::log_print("Data successfully normalized (Mode & Instability Index) and time-aligned.")
  
  return(list(
    normalized_data = normalized_data, 
    baselines_table = baselines_table
  ))
}

#' Summarize Processed Data
#'
#' @description (MODIFIED) Creates `data_per_pcr`, `summary_per_rep`, and `summary_by_group`.
#' FIXED: Robustly handles 'plate_raw' and 'well_raw' using any_of() to prevent crashes if missing.
#'
#' @param normalized_data The data from `normalize_data()`.
#' @param config The loaded configuration list.
#' @param all_grouping_vars A character vector of all grouping columns.
#' @return A named list containing the three summary data frames.
summarize_data <- function(normalized_data, config, all_grouping_vars) {
  
  response_vars <- names(config$response_variables)
  time_var <- config$key_variables$time_variable
  re_cross <- config$key_variables$optional_crossed_effect
  
  # Define metadata columns we want to keep IF they exist
  meta_cols <- c("plate", "well", "plate_raw", "well_raw", "fsa_filename", "run_date")
  
  # --- NEW: Preserve Label Columns ---
  # Check if these columns exist (created in Script 01 Part 2B)
  possible_labels <- c("Genotype_Pub", "Genotype_Exp")
  found_labels <- possible_labels[possible_labels %in% colnames(normalized_data)]
  data_per_pcr <- normalized_data
  
  # --- 1. Bio-Rep Summary ---
  bio_rep_grouping <- c(all_grouping_vars, time_var, "rep", re_cross, found_labels)
  bio_rep_grouping <- bio_rep_grouping[bio_rep_grouping %in% colnames(data_per_pcr)]
  
  summary_per_rep <- data_per_pcr %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(bio_rep_grouping))) %>%
    dplyr::summarise(
      dplyr::across(dplyr::all_of(response_vars), ~mean(.x, na.rm = TRUE)),
      # FIX: Use any_of() instead of calling plate_raw directly
      dplyr::across(dplyr::any_of(meta_cols), dplyr::first),
      .groups = "drop"
    )
  
  logr::log_print(paste("Created summary_per_rep with", nrow(summary_per_rep), "rows."))
  
  # --- 2. Experimental Group Summary ---
  exp_grouping <- c(all_grouping_vars, time_var, re_cross, found_labels)
  # Remove "rep" or "pcr" if they accidentally ended up in the grouping vars
  exp_grouping <- exp_grouping[!exp_grouping %in% c("rep", "pcr")]
  exp_grouping <- exp_grouping[exp_grouping %in% colnames(summary_per_rep)] # Safety check
  
  summary_by_group <- summary_per_rep %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(exp_grouping))) %>%
    dplyr::summarise(
      dplyr::across(
        dplyr::all_of(response_vars),
        list(avg = ~mean(.x, na.rm = TRUE),
             sd = ~sd(.x, na.rm = TRUE),
             sem = ~sd(.x, na.rm = TRUE) / sqrt(dplyr::n())),
        .names = "{.col}_{.fn}"
      ),
      # FIX: Use any_of() here as well
      dplyr::across(dplyr::any_of(meta_cols), dplyr::first),
      .groups = "drop"
    )
  
  logr::log_print(paste("Created summary_by_group with", nrow(summary_by_group), "rows."))
  
  return(list(
    data_per_pcr = data_per_pcr,
    summary_per_rep = summary_per_rep,
    summary_by_group = summary_by_group,
    exp_grouping = exp_grouping
  ))
}
#=============================================================================#
# SECTION 2: PLOTTING & MODELING (SCRIPT 02)
#=============================================================================#

#' Create Baseline Data for Plotting and Modeling
#'
#' @description Universal function to handle shared vs. individual baselines.
#' Adapts based on `baseline_grouping_level` and `share_baselines` in config.yml.
#'
#' @param data_per_pcr The per-PCR replicate data frame.
#' @param summary_per_rep The per-biological replicate data frame.
#' @param config The loaded configuration list.
#' @param all_grouping_vars A character vector of all grouping columns.
#' @return A named list containing `modeling_data` and `shared_baseline_pcr_points`.
#' Create Baseline Data for Plotting and Modeling
#'
#' @description Universal function to handle shared vs. individual baselines.
#' Adapts based on `baseline_grouping_level` and `share_baselines` in config.yml.
create_baseline_data <- function(data_per_pcr, summary_per_rep, config, all_grouping_vars) {
  
  time_var <- config$key_variables$time_variable
  baseline_method <- config$parameters$baseline_grouping_level %||% "group"
  
  share_baselines <- config$parameters$share_baselines %||% FALSE
  if (is.character(share_baselines)) share_baselines <- tolower(share_baselines) %in% c("true", "yes", "t", "1")
  
  # ── CRITICAL FIX: Identify label columns to strip before broadcasting ──
  label_cols <- c("Genotype_Pub", "Genotype_Exp", "label_text")
  found_labels <- label_cols[label_cols %in% colnames(data_per_pcr)]
  
  # Ensure the master combinations table HOLDS the correct labels
  all_group_combinations <- data_per_pcr %>%
    dplyr::distinct(dplyr::across(dplyr::all_of(c(all_grouping_vars, found_labels))))
  
  # =====================================================================#
  # PATH A: Global Control 
  # =====================================================================#
  if (baseline_method == "global_control") {
    logr::log_print("Creating shared baseline data for 'global_control' (Dropping all grouping vars at Day 0).")
    
    day0_pcr <- data_per_pcr %>%
      dplyr::filter(!!rlang::sym(time_var) == 0) %>%
      dplyr::select(-dplyr::any_of(c(all_grouping_vars, found_labels))) %>%
      dplyr::distinct()
    
    shared_baseline_pcr_points <- tidyr::crossing(all_group_combinations, day0_pcr)
    
    day0_bio <- summary_per_rep %>%
      dplyr::filter(!!rlang::sym(time_var) == 0) %>%
      dplyr::select(-dplyr::any_of(c(all_grouping_vars, found_labels))) %>%
      dplyr::distinct()
    
    expanded_baseline_bio_reps <- tidyr::crossing(all_group_combinations, day0_bio)
    non_baseline_bio_reps <- summary_per_rep %>% dplyr::filter(!!rlang::sym(time_var) > 0)
    
    modeling_data <- dplyr::bind_rows(expanded_baseline_bio_reps, non_baseline_bio_reps)
    
    # =====================================================================#
    # PATH B: Custom Method
    # =====================================================================#
  } else if (baseline_method == "custom") {
    
    if (isTRUE(share_baselines)) {
      # ── Backwards Compatibility Check ──
      drop_vars_raw <- config$parameters$pool_baselines_across
      
      # If the new variable is missing, look for the legacy variable name
      if (is.null(drop_vars_raw)) {
        drop_vars_raw <- config$parameters$baseline_drop_vars
        
        # If found under the old name, print a gentle deprecation notice
        if (!is.null(drop_vars_raw) && length(drop_vars_raw) > 0 && drop_vars_raw[1] != 'null') {
          logr::log_print(
            "NOTE: 'baseline_drop_vars' is deprecated in config.yml. Please update your config to use 'pool_baselines_across'. Continuing execution with legacy settings.", 
            console = TRUE
          )
        }
      }
      
      if (is.null(drop_vars_raw) || length(drop_vars_raw) == 0 || drop_vars_raw[1] == 'null') {
        logr::log_print("WARNING: 'custom' baseline requested with 'share_baselines: true', but 'baseline_drop_vars' is missing. Defaulting to individual baselines.")
        modeling_data <- summary_per_rep
        shared_baseline_pcr_points <- data_per_pcr %>% dplyr::filter(!!rlang::sym(time_var) == 0)
        
      } else {
        drop_vars <- unlist(drop_vars_raw)
        drop_vars <- drop_vars[drop_vars %in% colnames(data_per_pcr)]
        join_keys <- setdiff(all_grouping_vars, drop_vars)
        
        logr::log_print(paste("Creating custom shared baselines. Dropping:", paste(drop_vars, collapse = ", "), "| Segregating strictly by:", paste(join_keys, collapse = ", ")))
        
        # ── CRITICAL FIX: Strip old labels from the raw Day 0 points ──
        day0_pcr <- data_per_pcr %>%
          dplyr::filter(!!rlang::sym(time_var) == 0) %>%
          dplyr::select(-dplyr::any_of(c(drop_vars, found_labels))) %>% 
          dplyr::distinct() 
        
        shared_baseline_pcr_points <- all_group_combinations %>%
          dplyr::left_join(day0_pcr, by = join_keys, relationship = "many-to-many") %>%
          dplyr::filter(!is.na(!!sym(time_var))) 
        
        # ── CRITICAL FIX: Strip old labels from the raw Day 0 Bio-Reps ──
        day0_bio <- summary_per_rep %>%
          dplyr::filter(!!rlang::sym(time_var) == 0) %>%
          dplyr::select(-dplyr::any_of(c(drop_vars, found_labels))) %>%
          dplyr::distinct()
        
        expanded_baseline_bio_reps <- all_group_combinations %>%
          dplyr::left_join(day0_bio, by = join_keys, relationship = "many-to-many") %>%
          dplyr::filter(!is.na(!!sym(time_var)))
        
        non_baseline_bio_reps <- summary_per_rep %>% dplyr::filter(!!rlang::sym(time_var) > 0)
        modeling_data <- dplyr::bind_rows(expanded_baseline_bio_reps, non_baseline_bio_reps)
      }
      
    } else {
      logr::log_print("Method is 'custom', but 'share_baselines' is FALSE. Using exact individual Day 0 points.")
      modeling_data <- summary_per_rep
      shared_baseline_pcr_points <- data_per_pcr %>% dplyr::filter(!!rlang::sym(time_var) == 0)
    }
    
    # =====================================================================#
    # PATH C: Standard Individual Baselines
    # =====================================================================#
  } else {
    logr::log_print(paste("Using individual baselines. Method specified:", baseline_method))
    modeling_data <- summary_per_rep
    shared_baseline_pcr_points <- data_per_pcr %>% dplyr::filter(!!rlang::sym(time_var) == 0)
  }
  
  return(list(
    modeling_data = modeling_data,
    shared_baseline_pcr_points = shared_baseline_pcr_points
  ))
}

#' Perform Adaptive Statistical Modeling (Tournament Version)
#'
#' @description Fits models in a tournament style:
#' 1. Structure Test: Nested vs Simple (if applicable)
#' 2. Slope Test: Random Slope vs Random Intercept
#' 2.5 Batch Slope Test: Tests if batches grew at different speeds
#' 3. Final Test: Winner vs Simple LM
#' 4. Variance Test (AIC): Equal Variance vs Unequal Variance (gls/lme)
#'
#' @param data_summary The per-replicate summary data frame.
#' @param config The loaded configuration list.
#' @param response_variable The name of the column for the outcome variable.
#' @return A list containing `results_tables` and `model_object`.
run_statistical_model <- function(data_summary, config, response_variable) {
  
  # --- Helper Function: Extract Variance ---
  get_variance_summary <- function(model_obj) {
    if (inherits(model_obj, "lmerMod")) {
      var_corr_df <- as.data.frame(lme4::VarCorr(model_obj))
      total_var <- sum(var_corr_df$vcov, na.rm = TRUE)
      variance_summary_df <- var_corr_df %>%
        dplyr::select(Source = grp, Variance = vcov, Std.Dev = sdcor) %>%
        dplyr::mutate(
          ICC_Percent = (Variance / total_var) * 100,
          Note = ifelse(Variance < 1e-8, "Boundary fit: Variance is effectively zero.", NA)
        )
      return(variance_summary_df)
    }
    return(NULL) 
  }
  
  message(paste("Running Statistical Model (Tournament) for:", response_variable))
  interpretation_log <- c()
  interpretation_log <- c(interpretation_log, paste("--- Analysis Log for:", response_variable, "---"))
  
  # --- 1. Define variables ---
  cfg_vars <- config$key_variables
  resp_var <- response_variable
  time_var <- cfg_vars$time_variable
  primary_group <- cfg_vars$primary_group_var
  fixed_effects_to_test <- cfg_vars$model_fixed_effects %||% primary_group
  
  clean_fixed_effects <- unique(trimws(unlist(strsplit(fixed_effects_to_test, "[\\*:]"))))
  clean_fixed_effects <- clean_fixed_effects[clean_fixed_effects %in% colnames(data_summary)]
  
  base_formula_string <- paste(resp_var, "~", paste(c(time_var, fixed_effects_to_test), collapse = " * "))
  
  re_cross <- cfg_vars$optional_crossed_effect
  fixed_effects_formula <- base_formula_string 
  if (!is.null(re_cross) && re_cross != 'null') {
    message(paste("...adding '", re_cross, "' as a fixed effect."))
    interpretation_log <- c(interpretation_log, paste("Fixed Effect Added:", re_cross, "(Crossed Effect)"))
    fixed_effects_formula <- paste(base_formula_string, "+", re_cross)
    clean_fixed_effects <- unique(c(clean_fixed_effects, re_cross))
  }
  
  # --- SMART RANDOM EFFECT PARSING ---
  re_config_str <- cfg_vars$nesting_structure
  
  if (!is.null(re_config_str) && re_config_str != "null" && re_config_str != "") {
    random_effect_vars <- trimws(unlist(strsplit(re_config_str, "[/:]")))
    # Safely drop any variables that don't actually exist in the data
    random_effect_vars <- random_effect_vars[random_effect_vars %in% colnames(data_summary)]
    
    is_mixed_model <- length(random_effect_vars) > 0
    is_nested_model <- length(random_effect_vars) > 1
    
    primary_re <- if(length(random_effect_vars) > 0) random_effect_vars[1] else NULL
    
    # Rebuild the max structure safely based only on columns that exist
    structure_max <- paste(random_effect_vars, collapse = "/")
  } else {
    re1 <- cfg_vars$repeated_measure_var
    re2 <- cfg_vars$optional_random_effect
    
    if(is.null(re1) || re1 == 'null') re1 <- NULL
    if(is.null(re2) || re2 == 'null') re2 <- NULL
    
    is_mixed_model <- !is.null(re1) || !is.null(re2)
    is_nested_model <- !is.null(re1) && !is.null(re2)
    
    random_effect_vars <- c()
    if(!is.null(re1)) random_effect_vars <- c(random_effect_vars, re1)
    if(!is.null(re2)) random_effect_vars <- c(random_effect_vars, re2)
    
    primary_re <- re1 %||% re2
    structure_max <- if(is_nested_model) paste0(re1, "/", re2) else primary_re 
  }
  
  # --- 2. Data Validation & UID Creation ---
  required_vars <- unique(c(resp_var, time_var, clean_fixed_effects, random_effect_vars))
  
  if (!is.null(re_cross) && re_cross != 'null' && !is.null(re2) && re2 %in% names(data_summary)) {
    data_summary <- data_summary %>% dplyr::mutate(!!sym(re2) := as.factor(paste0("B", !!sym(re_cross), "_R", !!sym(re2))))
  }
  
  model_data_base <- data_summary %>%
    dplyr::ungroup() %>%
    dplyr::select(dplyr::all_of(required_vars)) %>%
    dplyr::filter(dplyr::if_all(dplyr::where(is.numeric), ~ !is.infinite(.))) %>%
    stats::na.omit() %>%
    as.data.frame()
  
  if (is_mixed_model) {
    model_data_base <- model_data_base %>% dplyr::mutate(dplyr::across(dplyr::all_of(random_effect_vars), as.factor))
  }
  
  # --- THE SAFETY CHECK (Avoid 1-Level Crash) ---
  if (is_mixed_model) {
    num_re_levels <- length(unique(model_data_base[[primary_re]]))
    if (num_re_levels < 2) {
      message("SAFETY TRIGGER: Random effect dropped to < 2 levels after NA omission. Reverting to LM.")
      interpretation_log <- c(interpretation_log, "SAFETY TRIGGER: Insufficient random effect levels. Forced simple LM.")
      is_mixed_model <- FALSE
    }
  }
  
  # --- 3. The Tournament ---
  model_fit <- NULL
  model_results <- list()
  model_formula_lm <- stats::as.formula(fixed_effects_formula)
  bobyqa_ctrl <- lme4::lmerControl(optimizer = "bobyqa", optCtrl = list(maxfun = 1e5))
  
  if (is_mixed_model) {
    
    formula_max_slope <- stats::as.formula(paste(fixed_effects_formula, "+ (1 +", time_var, "|", structure_max, ")"))
    formula_max_int   <- stats::as.formula(paste(fixed_effects_formula, "+ (1 |", structure_max, ")"))
    
    # RESTORED: Initialization Block
    message("Initializing: Fitting maximal model to identify clean data...")
    
    temp_fit <- tryCatch({
      lmerTest::lmer(formula_max_slope, data = model_data_base, REML = FALSE, na.action = na.omit, control = bobyqa_ctrl)
    }, error = function(e) { NULL })
    
    if(is.null(temp_fit)) {
      temp_fit <- tryCatch({
        lmerTest::lmer(formula_max_int, data = model_data_base, REML = FALSE, na.action = na.omit, control = bobyqa_ctrl)
      }, error = function(e) { NULL })
    }
    
    if(is.null(temp_fit)){
      message("NOTE: All LMER initialization failed. Falling back to simple LM.")
      interpretation_log <- c(interpretation_log, "ERROR: Mixed Model (LMER) initialization failed on full dataset. Reverting to Simple Linear Model (LM).")
      model_data_clean <- model_data_base 
    } else {
      model_data_clean <- temp_fit@frame
      message(paste("...Clean data identified:", nrow(model_data_clean), "rows."))
      interpretation_log <- c(interpretation_log, paste("Data Cleaning: Final dataset has", nrow(model_data_clean), "valid rows (missing values removed)."))
    }
    
    fit_LM <- stats::lm(model_formula_lm, data = model_data_clean)
    
    if(!is.null(temp_fit)) {
      
      # --- ROUND 1: Structure ---
      # CRITICAL FIX 1: Use primary_re instead of re1 as the baseline
      current_best_structure <- primary_re 
      
      if(is_nested_model) {
        message("--- Round 1: Nested vs Simple Structure ---")
        interpretation_log <- c(interpretation_log, "--- Tournament Round 1: Random Effect Structure ---")
        
        # CRITICAL FIX 2: Use primary_re here as well
        f_simple <- stats::as.formula(paste(fixed_effects_formula, "+ (1 |", primary_re, ")"))
        m_simple <- try(lmerTest::lmer(f_simple, data = model_data_clean, REML = FALSE, control = bobyqa_ctrl), silent=TRUE)
        
        f_nested <- stats::as.formula(paste(fixed_effects_formula, "+ (1 |", structure_max, ")"))
        m_nested <- try(lmerTest::lmer(f_nested, data = model_data_clean, REML = FALSE, control = bobyqa_ctrl), silent=TRUE)
        
        if(!inherits(m_simple, "try-error") && !inherits(m_nested, "try-error")) {
          lrt_struct <- anova(m_simple, m_nested)
          p_struct <- lrt_struct$`Pr(>Chisq)`[2]
          model_results[["Test_Round1_Structure"]] <- as.data.frame(lrt_struct)
          
          if(!is.na(p_struct) && p_struct < 0.05) {
            message(paste("...Result: Nested Structure (", structure_max, ") wins."))
            current_best_structure <- structure_max
            interpretation_log <- c(interpretation_log, paste("Result: Nested Structure (", structure_max, ") is significantly better (p =", round(p_struct, 4), ")."))
          } else {
            message(paste("...Result: Simple Structure (", primary_re, ") is sufficient."))
            current_best_structure <- primary_re
            interpretation_log <- c(interpretation_log, paste("Result: Simple Structure (", primary_re, ") is sufficient. Nested structure was not significant (p =", round(p_struct, 4), ")."))
          }      
        }
      } else {
        message("--- Round 1 Skipped: No nested structure detected ---")
        # CRITICAL FIX 3: Update the log message so it doesn't print a blank space if re1 is empty
        interpretation_log <- c(interpretation_log, paste("Round 1 Skipped: Only simple random effect requested (", primary_re, ")."))
      }      
      
      # --- ROUND 2: Slope vs Intercept ---
      message("--- Round 2: Random Slope vs Intercept ---")
      interpretation_log <- c(interpretation_log, "--- Tournament Round 2: Random Slope vs. Random Intercept ---")
      
      f_int <- stats::as.formula(paste(fixed_effects_formula, "+ (1 |", current_best_structure, ")"))
      m_int <- lmerTest::lmer(f_int, data = model_data_clean, REML = FALSE, control = bobyqa_ctrl)
      
      f_slope <- stats::as.formula(paste(fixed_effects_formula, "+ (1 +", time_var, "|", current_best_structure, ")"))
      
      m_slope <- tryCatch({
        lmerTest::lmer(f_slope, data = model_data_clean, REML = FALSE, control = bobyqa_ctrl)
      }, error = function(e) {
        message("...Slope model errored out. Using Intercept.")
        return(NULL)
      })
      
      best_lmer_model <- m_int
      best_lmer_formula <- f_int
      lmer_desc <- "Random Intercept"
      random_effect_string <- paste("~ 1 |", current_best_structure)
      
      if(!is.null(m_slope)) {
        lrt_slope <- anova(m_int, m_slope)
        p_slope <- lrt_slope$`Pr(>Chisq)`[2]
        model_results[["Test_Round2_Slope"]] <- as.data.frame(lrt_slope)
        
        if(!is.na(p_slope) && p_slope < 0.05) {
          message("...Result: Random Slope wins.")
          best_lmer_model <- m_slope
          best_lmer_formula <- f_slope
          lmer_desc <- "Random Slope + Intercept"
          random_effect_string <- paste("~", time_var, "|", current_best_structure)
          interpretation_log <- c(interpretation_log, paste("Result: Random Slope model is significantly better (p =", round(p_slope, 4), ")."))
        } else {
          message("...Result: Random Slope not significant. Keeping Intercept only.")
          interpretation_log <- c(interpretation_log, paste("Result: Random Intercept model is sufficient. Random Slope was not significant (p =", round(p_slope, 4), ")."))
        }
      }
      
      # --- ROUND 2.5: Batch Random Slope vs Batch Random Intercept ---
      if (!is.null(re_cross) && re_cross != "null") {
        message("--- Round 2.5: Testing Batch Random Slope ---")
        interpretation_log <- c(interpretation_log, "--- Tournament Round 2.5: Testing Batch Random Slope ---")
        
        current_formula_str <- paste(deparse(best_lmer_formula, width.cutoff = 500), collapse = " ")        
        complex_batch_form <- stats::as.formula(paste(current_formula_str, "+ (0 +", time_var, "|", re_cross, ")"))
        
        m_batch_slope <- tryCatch({
          lmerTest::lmer(complex_batch_form, data = model_data_clean, REML = FALSE, control = bobyqa_ctrl)
        }, error = function(e) { NULL })
        
        if (!is.null(m_batch_slope) && !lme4::isSingular(m_batch_slope)) {
          batch_anova <- anova(best_lmer_model, m_batch_slope)
          batch_p_val <- batch_anova$`Pr(>Chisq)`[2]
          
          model_results[["Test_Round2.5_Batch_Slope"]] <- as.data.frame(batch_anova)
          
          if (!is.na(batch_p_val) && batch_p_val < 0.05) {
            message("...Result: Batch Random Slope wins (Batches grew at different speeds).")
            best_lmer_formula <- complex_batch_form
            best_lmer_model <- m_batch_slope
            lmer_desc <- paste(lmer_desc, "+ Batch Slope")
            random_effect_string <- paste("list(", current_best_structure, " = ", random_effect_string, ", ", re_cross, " = ~", time_var, " - 1)")
            interpretation_log <- c(interpretation_log, paste("Result: Adding a random slope for", re_cross, "significantly improved the model (p =", round(batch_p_val, 4), ")."))
          } else {
            message("...Result: Batch Random Slope not significant. Keeping fixed batch intercept.")
            interpretation_log <- c(interpretation_log, paste("Result: Batch random slope was not significant (p =", round(batch_p_val, 4), ")."))
          }
        } else {
          message("...Result: Batch Random Slope model failed or was singular. Skipping.")
          interpretation_log <- c(interpretation_log, paste("Result: Batch Random Slope model was singular/unstable. Maintained previous model."))
        }
      }
      
      # --- ROUND 3: LMER vs LM ---
      message("--- Round 3: Final LMER vs LM ---")
      interpretation_log <- c(interpretation_log, "--- Tournament Round 3: Mixed Model vs. Fixed Effects Only ---")
      
      lrt_final <- anova(best_lmer_model, fit_LM)
      p_final <- lrt_final$`Pr(>Chisq)`[2]
      model_results[["Test_Round3_LM_vs_LMER"]] <- as.data.frame(lrt_final)
      
      if(!is.na(p_final) && p_final < 0.05) {
        message(paste("...Result: LMER wins. Final Model:", lmer_desc))
        interpretation_log <- c(interpretation_log, paste("Result: Mixed Model (", lmer_desc, ") explains significantly more variance than LM (p =", round(p_final, 4), ")."))
        
        message("...Refitting winning LMER with REML=TRUE")
        winning_equal_var_model <- lmerTest::lmer(best_lmer_formula, data = model_data_clean, REML = TRUE, control = bobyqa_ctrl)
        winning_type <- "lmer"
        interpretation_log <- c(interpretation_log, paste("Final Formula (REML):", deparse(best_lmer_formula)))
        
      } else {
        message("...Result: LMER not better than LM. Using LM.")
        interpretation_log <- c(interpretation_log, paste("Result: Mixed Model was NOT significantly better than simple LM (p =", round(p_final, 4), ")."))
        winning_equal_var_model <- fit_LM
        winning_type <- "lm"
      }
    } else {
      winning_equal_var_model <- fit_LM
      winning_type <- "lm"
    }
  } else {
    model_data_clean <- model_data_base
    winning_equal_var_model <- stats::lm(model_formula_lm, data = model_data_clean)
    winning_type <- "lm"
    
    # --- NEW: LM FALLBACK LOGGING ---
    interpretation_log <- c(interpretation_log, 
                            "--- Mixed Model Tournament Skipped ---",
                            "Reason: No random effects were configured, or the random effect dropped to fewer than 2 levels after missing data was removed.",
                            "Action: Fit a standard Fixed-Effects Linear Model (LM).",
                            paste("Formula:", paste(deparse(model_formula_lm, width.cutoff = 500), collapse = " "))
    )

    message("=== QC: Checking factor levels in model_data_base ===")
    
    for (col in names(model_data_base)) {
      if (is.factor(model_data_base[[col]])) {
        lvls <- levels(model_data_base[[col]])
        actual <- unique(as.character(model_data_base[[col]]))
        message(paste0("  FACTOR '", col, "': levels=", paste(lvls, collapse=","), 
                       " | actual values=", paste(actual, collapse=",")))
      }
    }
    message("=== END QC ===")
  }
  
  # --- ROUND 4: THE VARIANCE TOURNAMENT (AIC Check) ---
  message("--- Round 4: Equal vs Unequal Variance (AIC Tournament) ---")
  final_model_fit <- winning_equal_var_model
  
  # FORCE-FEED the exact formula to nlme using do.call
  unequal_var_model <- tryCatch({
    if (winning_type == "lm") {
      do.call(nlme::gls, list(
        model = model_formula_lm, 
        data = model_data_clean, 
        weights = nlme::varIdent(form = stats::as.formula(paste("~ 1 |", primary_group)))
      ))
    } else {
      do.call(nlme::lme, list(
        fixed = model_formula_lm, 
        random = stats::as.formula(random_effect_string), 
        data = model_data_clean, 
        weights = nlme::varIdent(form = stats::as.formula(paste("~ 1 |", primary_group)))
      ))
    }
  }, error = function(e) { NULL })
  
  if (!is.null(unequal_var_model)) {
    aic_equal <- stats::AIC(winning_equal_var_model)
    aic_unequal <- stats::AIC(unequal_var_model)
    
    # --- NEW: Save results to the model_results list ---
    variance_test_df <- data.frame(
      Variance_Model = c("Equal Variance (Standard)", "Unequal Variance (gls/lme)"),
      AIC_Value = c(aic_equal, aic_unequal),
      AIC_Difference = c(0, aic_equal - aic_unequal), # Positive means unequal is better
      Decision = if (aic_unequal < (aic_equal - 2)) c("Runner-up", "Winner (Selected)") else c("Winner (Selected)", "Runner-up")
    )
    model_results[["Test_Round4_Variance_AIC"]] <- variance_test_df

    if (aic_unequal < (aic_equal - 2)) {
      message(paste("...Result: Unequal Variance model wins! (AIC drop:", round(aic_equal - aic_unequal, 1), ")"))
      interpretation_log <- c(interpretation_log, "Result: Unequal Variance model (gls/lme) selected via AIC.")
      final_model_fit <- unequal_var_model
    } else {
      message("...Result: Equal Variance sufficient. Keeping standard model.")
      interpretation_log <- c(interpretation_log, "Result: Equal Variance model maintained via AIC.")
    }
  } else {
    message("...Result: Unequal Variance model failed to converge. Keeping standard model.")
    interpretation_log <- c(interpretation_log, "Result: Unequal Variance model failed to converge. Keeping standard model.")
    
    # --- NEW: Log the failure cleanly in the exported tables ---
    variance_test_df <- data.frame(
      Variance_Model = c("Equal Variance (Standard)", "Unequal Variance (gls/lme)"),
      AIC_Value = c(stats::AIC(winning_equal_var_model), NA),
      AIC_Difference = c(0, NA),
      Decision = c("Winner (Selected)", "Failed to Converge")
    )
    model_results[["Test_Round4_Variance_AIC"]] <- variance_test_df
  }
  
  # --- LOG THE UNIVERSAL FINAL FORMULA ---
  try({
    # 1. Extract the base formula safely (collapsing multi-line formulas into one string)
    final_form_str <- paste(trimws(deparse(stats::formula(final_model_fit), width.cutoff = 500)), collapse = " ")
    
    # 2. nlme::lme objects hide their random effects from the standard formula() call, so we append them manually
    if (inherits(final_model_fit, "lme")) {
      final_form_str <- paste(final_form_str, "| Random:", random_effect_string)
    }
    
    # 3. If the Unequal Variance model won Round 4, append the variance structure to the log
    if (!is.null(unequal_var_model) && identical(final_model_fit, unequal_var_model)) {
      final_form_str <- paste(final_form_str, "| Variance: varIdent(form = ~ 1 |", primary_group, ")")
    }
    
    # 4. Push to the log
    interpretation_log <- c(interpretation_log, "--------------------------------------------------")
    interpretation_log <- c(interpretation_log, paste(">>> FINAL DEPLOYED MODEL FORMULA:", final_form_str))
    interpretation_log <- c(interpretation_log, "--------------------------------------------------")
  }, silent = TRUE)
  
  # --- 4. Results Extraction & Formatting ---
  format_table <- function(df) {
    if(is.null(df)) return(NULL)
    df %>%
      dplyr::mutate(across(where(is.numeric), ~ round(., 5))) %>%
      dplyr::rename_with(~ dplyr::case_when(
        . == "term" ~ "Parameter", . == "estimate" ~ "Estimate", . == "std.error" ~ "SE",
        . == "statistic" ~ "t_Stat", . == "p.value" ~ "P_Value", . == "df" ~ "DF",
        . == "lower.CL" ~ "CI_Lower_95", . == "upper.CL" ~ "CI_Upper_95",
        . == "asymp.LCL" ~ "CI_Lower_95", . == "asymp.UCL" ~ "CI_Upper_95",
        . == "conf.low" ~ "CI_Lower_95", . == "conf.high" ~ "CI_Upper_95", TRUE ~ .
      ))
  }
  
  # A. Coefficients
  raw_coefs <- broom.mixed::tidy(final_model_fit, effects = "fixed", conf.int = TRUE)
  model_results[["Model_Coefficients"]] <- format_table(raw_coefs)
  
  # B. ANOVA
  try({
    raw_anova <- as.data.frame(anova(final_model_fit)) %>% tibble::rownames_to_column("Factor")
    model_results[["ANOVA_Summary"]] <- format_table(raw_anova) %>% dplyr::rename_with(~ gsub("Pr(>F)", "P_Value", ., fixed=TRUE))
    
    # --- NEW: DYNAMIC SIGNIFICANCE LOGGING ---
    interpretation_log <- c(interpretation_log, "--- Final Fixed Effects Significance (ANOVA) ---")
    
    # Smart-search for the p-value column (handles lm, lmer, and gls outputs)
    pval_col <- names(raw_anova)[grepl("Pr|p.value|p-value", names(raw_anova), ignore.case = TRUE)][1]
    
    if (!is.na(pval_col)) {
      for (i in 1:nrow(raw_anova)) {
        factor_name <- raw_anova$Factor[i]
        p_val <- as.numeric(raw_anova[[pval_col]][i])
        
        # Ignore intercept/residuals and only log actual biological factors
        if (!is.na(p_val) && !grepl("Residual|Intercept", factor_name, ignore.case = TRUE)) {
          if (p_val <= 0.001) {
            interpretation_log <- c(interpretation_log, paste("[***] HIGHLY SIGNIFICANT:", factor_name, "(p =", round(p_val, 4), ")"))
          } else if (p_val <= 0.05) {
            interpretation_log <- c(interpretation_log, paste("[ * ] SIGNIFICANT:", factor_name, "(p =", round(p_val, 4), ")"))
          } else {
            interpretation_log <- c(interpretation_log, paste("[ - ] Not Significant:", factor_name, "(p =", round(p_val, 4), ")"))
          }
        }
      }
    } else {
      interpretation_log <- c(interpretation_log, "Note: Could not automatically parse ANOVA p-values for the log.")
    }
    # -----------------------------------------
  }, silent = TRUE)
  
  # C. Group Slopes (For overall bars)
  try({
    trend_specs <- stats::as.formula(paste("~", paste(setdiff(clean_fixed_effects, re_cross), collapse = " + ")))
    
    em_trends <- emmeans::emtrends(final_model_fit, specs = trend_specs, var = time_var, 
                                   data = model_data_clean, adjust = "tukey", lmer.df = "asymptotic")
    emm_df <- as.data.frame(em_trends)
    
    trend_col_name <- paste0(time_var, ".trend")
    if(trend_col_name %in% names(emm_df)) names(emm_df)[names(emm_df) == trend_col_name] <- "Expansion_Rate"
    
    # 1. Format the table first so we have standardized column names (CI_Lower_95, etc.)
    formatted_emm_df <- format_table(emm_df)
    
    # --- NEW: CALCULATE ADVANCED METRICS ---
    # A. Find the total time duration of the experiment
    total_time_duration <- max(model_data_clean[[time_var]], na.rm = TRUE) - min(model_data_clean[[time_var]], na.rm = TRUE)
    
    # B. Identify the grouping column (usually the first column, e.g., 'Genotype')
    group_col <- names(formatted_emm_df)[1]
    
    # C. Safely extract the WT baseline rate
    wt_rate <- formatted_emm_df %>%
      dplyr::filter(grepl("WT|Wildtype", !!sym(group_col), ignore.case = TRUE)) %>%
      dplyr::pull(Expansion_Rate) %>%
      mean(na.rm = TRUE) # Uses mean just in case there are multiple WT-like rows
    
    # If no WT is found, set to NA so the script doesn't crash
    if (is.nan(wt_rate) || length(wt_rate) == 0) wt_rate <- NA_real_
    
    # D. Calculate all metrics
    formatted_emm_df <- formatted_emm_df %>%
      dplyr::mutate(
        # 1. Total Expansion Over Time
        Total_Expansion_Est = Expansion_Rate * total_time_duration,
        Total_Expansion_CI_Lower = CI_Lower_95 * total_time_duration,
        Total_Expansion_CI_Upper = CI_Upper_95 * total_time_duration,
        
        # 2. Time to expand exactly 1 CAG (only if expanding)
        Rate_per_CAG = dplyr::case_when(
          Expansion_Rate > 0 ~ 1 / Expansion_Rate,
          TRUE ~ NA_real_  # Returns NA if the rate is 0 or negative (contracting)
        ),
        
        # 3. Rate Change compared to WT
        Rate_Difference_vs_WT = Expansion_Rate - wt_rate,
        Rate_FoldChange_vs_WT = Expansion_Rate / wt_rate
      ) %>%
      # Round the new numeric columns to 4 decimal places for clean viewing
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(.x, 4)))
    # ---------------------------------------
    
    # Save it back to the results list
    model_results[["Group_Expansion_Rates"]] <- formatted_emm_df
    
    # --- NEW: GET BOTH RAW AND ADJUSTED P-VALUES ---
    # 1. Calculate raw p-values (No adjustment)
    pairs_raw <- as.data.frame(summary(graphics::pairs(em_trends), adjust = "none")) %>%
      dplyr::rename(p_value = p.value)
    
    # 2. Calculate adjusted p-values (Tukey adjustment for multiple comparisons)
    pairs_adj <- as.data.frame(summary(graphics::pairs(em_trends), adjust = "tukey")) %>%
      dplyr::rename(p_value_adj = p.value)
    
    # 3. Merge them perfectly by their contrast name
    emm_pairs_combined <- pairs_raw %>%
      dplyr::left_join(pairs_adj %>% dplyr::select(contrast, p_value_adj), by = "contrast")
    
    # 4. Format and save
    model_results[["Rate_Comparisons"]] <- format_table(emm_pairs_combined) %>% 
      dplyr::rename(Rate_Difference = Estimate)
    # -----------------------------------------------
    
    # --- MARGINAL EFFECTS FOR BLUP CENTERING ---
    marginal_effects <- setdiff(clean_fixed_effects, re_cross)
    if(length(marginal_effects) > 0) {
      specs_marg <- stats::as.formula(paste("~", paste(rev(marginal_effects), collapse = " | ")))
      emm_marg_raw <- emmeans::emtrends(final_model_fit, specs = specs_marg, var = time_var, 
                                        data = model_data_clean, lmer.df = "asymptotic")
      emm_df_marg <- as.data.frame(emm_marg_raw)
      if(trend_col_name %in% names(emm_df_marg)) names(emm_df_marg)[names(emm_df_marg) == trend_col_name] <- "Expansion_Rate"
    } else {
      emm_df_marg <- emm_df
    }
    
    # -------------------------------------------------------------------------#
    # D. TRUE BLUP EXTRACTION
    # -------------------------------------------------------------------------#
    individual_slopes_final <- NULL
    grouping_var_blup <- if(exists("primary_re")) primary_re else primary_group   
    naive_grouping <- unique(c(primary_group, cfg_vars$secondary_group_var, grouping_var_blup))
    naive_grouping <- naive_grouping[!is.null(naive_grouping) & naive_grouping != 'null' & naive_grouping %in% colnames(data_summary)]
    
    try({
      if (length(naive_grouping) > 0) {
        individual_slopes_final <- data_summary %>%
          dplyr::ungroup() %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(naive_grouping))) %>%
          dplyr::summarise(
            Individual_Slope = {
              x_vals <- .data[[time_var]]
              y_vals <- .data[[resp_var]]
              valid_idx <- !is.na(x_vals) & !is.na(y_vals)
              x_clean <- x_vals[valid_idx]
              y_clean <- y_vals[valid_idx]
              if (length(x_clean) > 1 && stats::var(x_clean) > 0) { stats::cov(x_clean, y_clean) / stats::var(x_clean) } else { NA_real_ }
            }, .groups = "drop"
          ) %>%
          dplyr::filter(!is.na(Individual_Slope)) %>%
          dplyr::mutate(Method = "Naive (Pooled by Clone)")
      }
    }, silent = TRUE)
    
    if (winning_type == "lmer") {
      try({
        ranef_list <- lme4::ranef(winning_equal_var_model)
        exact_match <- names(ranef_list)[names(ranef_list) == grouping_var_blup]
        target_coef_name <- if(length(exact_match) > 0) exact_match[1] else names(ranef_list)[grepl(grouping_var_blup, names(ranef_list))][1]
        
        if (!is.na(target_coef_name)) {
          blup_df <- ranef_list[[target_coef_name]]
          
          if (time_var %in% colnames(blup_df)) {
            blup_df <- blup_df %>%
              tibble::rownames_to_column("compound_id") %>%
              dplyr::select(compound_id, Intercept_Deviation = `(Intercept)`, Slope_Deviation = !!sym(time_var))            
            
            if (grepl(":", target_coef_name)) {
              split_names <- unlist(strsplit(target_coef_name, ":"))
              blup_df <- blup_df %>% tidyr::separate(compound_id, into = split_names, sep = ":")
            } else {
              blup_df <- blup_df %>% dplyr::rename(!!sym(target_coef_name) := compound_id)
            }
            
            join_cols <- intersect(names(blup_df), names(data_summary))
            group_map <- data_summary %>%
              dplyr::ungroup() %>% dplyr::select(dplyr::all_of(naive_grouping)) %>%
              dplyr::distinct() %>% dplyr::mutate(across(all_of(join_cols), as.character))
            
            blup_df <- blup_df %>% dplyr::mutate(across(all_of(join_cols), as.character))
            mapped_deviations <- blup_df %>% dplyr::inner_join(group_map, by = join_cols)
            
            emm_join_cols <- intersect(names(mapped_deviations), names(emm_df_marg))
            emm_clean <- emm_df_marg %>% dplyr::mutate(across(all_of(emm_join_cols), as.character))
            mapped_deviations <- mapped_deviations %>% dplyr::mutate(across(all_of(emm_join_cols), as.character))
            
            individual_slopes_final <- mapped_deviations %>%
              dplyr::left_join(emm_clean %>% dplyr::select(all_of(emm_join_cols), Group_Base_Slope = Expansion_Rate), by = emm_join_cols) %>%
              dplyr::mutate(Individual_Slope = Group_Base_Slope + Slope_Deviation, Method = "True BLUP (Model-Adjusted)") %>%
              dplyr::filter(!is.na(Individual_Slope))
            
          } else { message("   -> NOTE: Model dropped Random Slopes. Keeping Naive Slopes.") }
        }
      }, silent = TRUE)
    }
    
    if(!is.null(individual_slopes_final) && nrow(individual_slopes_final) > 0) {
      model_results[["Individual_Slopes_BLUPs"]] <- individual_slopes_final
      message(paste0("   -> Individual dots calculated using: ", individual_slopes_final$Method[1]))
    }
    
  }, silent = FALSE)
  
  # E. Variance Summary 
  variance_summary_df <- get_variance_summary(winning_equal_var_model)
  if (!is.null(variance_summary_df)) {
    model_results[["Random_Effects_Variance"]] <- variance_summary_df %>% dplyr::mutate(across(where(is.numeric), ~ round(., 4)))
  }
  
  model_results[["Analysis_Log"]] <- data.frame(Log_Entry = interpretation_log)
  return(list(results_tables = model_results, model_object = final_model_fit))
}
#=============================================================================#
# SECTION 3: PLOTTING THEMES & HELPERS (SCRIPT 01 & 02)
#=============================================================================#

#' Custom ggplot Theme for Publication-Ready Plots
#'
#' @param base_size The base font size.
#' @return A ggplot theme object.
theme_publication <- function(base_size = 12) {
  ggplot2::theme_classic(base_size = base_size) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1.2), hjust = 0, margin = ggplot2::margin(b = 5)),
      plot.subtitle = ggplot2::element_text(size = ggplot2::rel(1), hjust = 0, margin = ggplot2::margin(b = 10)),
      plot.caption = ggplot2::element_text(hjust = 1, size = ggplot2::rel(0.8), color = "grey30"),
      plot.background = ggplot2::element_rect(fill = "white", colour = NA),
      plot.margin = ggplot2::margin(10, 10, 10, 10),
      panel.background = ggplot2::element_rect(fill = "white", colour = NA),
      panel.grid.major = ggplot2::element_line(colour = "grey92", linewidth = 0.5),
      panel.grid.minor = ggplot2::element_blank(),
      axis.line = ggplot2::element_line(colour = "black", linewidth = 0.7, lineend = "square"),
      axis.ticks = ggplot2::element_line(colour = "black", linewidth = 0.5),
      axis.title = ggplot2::element_text(face = "bold", size = ggplot2::rel(1)),
      axis.text = ggplot2::element_text(size = ggplot2::rel(0.9)),
      legend.title = ggplot2::element_text(face = "plain"),
      legend.background = ggplot2::element_rect(fill = "white", colour = NA),
      legend.key = ggplot2::element_rect(fill = "white", colour = NA),
      legend.position = "right",
      strip.background = ggplot2::element_blank(),
      strip.text = ggplot2::element_text(face = "bold", size = ggplot2::rel(1), colour = "black", hjust = 0, margin = ggplot2::margin(2, 2, 2, 2))
    )
}

#' Define a Custom Color Palette (Smart Version)
#'
#' @description Priorities: 
#' 1. Config 'color_mapping'
#' 2. High-contrast set (if n <= 4)
#' 3. Standard extended palette
#'
#' @param data The dataframe containing the grouping variable.
#' @param config The loaded config object.
#' @param group_var The name of the column that contains the group labels.
#' @return A named vector of hexadecimal color codes.
create_custom_palette <- function(data, config, group_var) {
  if (is.null(group_var) || !group_var %in% colnames(data)) return(NULL)
  
  # Ensure factor
  if (!is.factor(data[[group_var]])) data[[group_var]] <- factor(data[[group_var]])
  group_levels <- levels(data[[group_var]])
  n_levels <- length(group_levels)
  
  if (n_levels == 0) return(NULL)
  
  # Initialize palette vector
  palette <- setNames(rep(NA, n_levels), group_levels)
  
  # 1. Check Config for Manual Mappings
  manual_colors <- config$color_mapping
  if (!is.null(manual_colors)) {
    for (lvl in group_levels) {
      if (lvl %in% names(manual_colors)) {
        palette[lvl] <- manual_colors[[lvl]]
      }
    }
  }
  
  # 2. Identify unassigned levels
  unassigned_levels <- names(palette)[is.na(palette)]
  n_unassigned <- length(unassigned_levels)
  
  if (n_unassigned > 0) {
    # 3. Choose Auto-Palette based on remaining count
    auto_colors <- c("#0072B2", "#CC79A7", "#F0E442", "#009E73", "#D55E00", "#000000", "#7570B3", "#A6761D")
    std_colors <- c(
      "#0072B2", "#CC79A7", "#009E73", "#F0E442", "#56B4E9",
      "#E69F00", "#000000", "#7570B3", "#E7298A", "#66A61E", "#E6AB02",
      "#A6761D", "#666666", "#1B9E77", "#D95F02", "#D55E00"
    )
    
    base_set <- if (n_unassigned <= 4) auto_colors else std_colors
    
    # Filter out colors already used manually
    if (!is.null(manual_colors)) {
      base_set <- setdiff(base_set, unlist(manual_colors))
    }
    
    # Safety Net: If manual colors consumed the entire base set
    if (length(base_set) == 0) {
      base_set <- c("#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd")
    }
    
    # INTERPOLATION FIX: Dynamically stretch the palette if we need more colors than exist
    if (n_unassigned > length(base_set)) {
      selected_colors <- grDevices::colorRampPalette(base_set)(n_unassigned)
    } else {
      selected_colors <- base_set[1:n_unassigned]
    }
    
    # Assign them
    for (i in 1:n_unassigned) {
      palette[unassigned_levels[i]] <- selected_colors[i]
    }
  }
  
  return(palette)
}

#' Generate a 96-well Plate Heatmap Visualization (Faceted by Plate)
#'
#' @description (MODIFIED) Generates a single figure for ONE metric,
#' faceted to show all plates. All panels share the same color scale.
#' Now uses 'plate_raw' and 'well_raw' for labels and facets.
#'
#' @param data_to_plot The *entire* data_per_pcr data frame.
#' @param config The loaded configuration list.
#' @param current_resp_var The name of the *one* response variable to plot (e.g., "mode").
#' @param current_label The "pretty label" for that variable (e.g., "Modal CAG").
#' @param val_limits A numeric vector c(min, max) for the color scale.
#' @return A ggplot object representing the faceted plate heatmaps.
generate_plate_heatmap <- function(data_to_plot, config, 
                                   current_resp_var, current_label, 
                                   val_limits) { 
  
  primary_var <- config$key_variables$primary_group_var
  
  # --- 1. Prepare data for this metric (Using raw IDs) ---
  plot_data_metric <- data_to_plot %>%
    dplyr::select(
      plate_raw, well_raw, !!sym(primary_var), # <-- (MODIFIED)
      value = !!sym(current_resp_var)
    ) %>%
    dplyr::mutate(
      # Use well_raw to get Row/Column
      Row_char = toupper(stringr::str_extract(well_raw, stringr::regex("[A-H]", ignore_case = TRUE))),
      Column = as.integer(stringr::str_extract(well_raw, "\\d+")),
      Row = factor(Row_char, levels = rev(LETTERS[1:8]))
    ) %>%
    dplyr::select(-Row_char) 
  
  # --- 2. Create the labels (Genotype or Well) ---
  label_data <- plot_data_metric %>%
    dplyr::select(plate_raw, well_raw, !!sym(primary_var), Row, Column) %>% # <-- (MODIFIED)
    dplyr::distinct() 
  
  if (!is.null(primary_var) && primary_var != 'null' && primary_var %in% colnames(label_data)) {
    label_data$label_text <- as.character(label_data[[primary_var]])
  } else {
    label_data$label_text <- as.character(label_data$well_raw) # <-- (MODIFIED)
  }
  
  # --- 3. Create the plot ---
  heatmap_plot <- ggplot2::ggplot(
    plot_data_metric, 
    ggplot2::aes(x = Column, y = Row)
  ) +
    ggplot2::geom_tile(ggplot2::aes(fill = value), color = "white", linewidth = 0.5) +
    ggplot2::geom_text(
      data = label_data, 
      ggplot2::aes(label = label_text), 
      size = 2.0, lineheight = .8, vjust = 0.5, color = "black"
    ) +
    ggplot2::scale_fill_viridis_c(
      option = "plasma", 
      name = "Value",
      limits = val_limits
    ) +
    ggplot2::scale_x_continuous(breaks = 1:12, expand = c(0, 0)) +
    ggplot2::coord_fixed() +
    
    # (MODIFIED) Facet by the raw plate ID
    facet_wrap(vars(plate_raw)) + 
    
    ggplot2::labs(
      title = paste("Plate Layout Heatmap (Original IDs):", current_label), # <-- (MODIFIED)
      subtitle = "Color indicates value. All plates shown share the same color scale.",
      x = "Plate Column",
      y = "Plate Row"
    ) +
    ggplot2::theme_minimal(base_size = 12) +
    ggplot2::theme(
      plot.title = ggplot2::element_text(face = "bold"),
      strip.text = ggplot2::element_text(face = "bold"),
      panel.grid = ggplot2::element_blank()
    )
  
  return(heatmap_plot)
}

#' Create Plot Aesthetics with Conditional Shape (Robust Version)
#'
#' @description Merges a base aesthetic mapping with a conditional shape aesthetic.
#' This simplifies adding 'shape' for a crossed effect (e.g., batch).
#' This version uses base R utils::modifyList for broad compatibility.
#'
#' @param base_aes The base ggplot aesthetic mapping (e.g., aes(x = col1, y = col2)).
#' @param shape_var A character string (e.g., "batch") for the shape variable.
#'                  If NULL or "null", no shape is added.
#' @return A merged aesthetic mapping.
#'
create_plot_aes <- function(base_aes, shape_var = NULL) {
  if (is.null(shape_var) || shape_var == 'null') {
    return(base_aes) # Return the original if no shape var
  }
  
  # Create the new aesthetic to add
  shape_aes <- ggplot2::aes(shape = !!rlang::sym(shape_var))
  
  # Use modifyList to merge the two lists.
  # This will safely add 'shape' or overwrite it if it already exists.
  final_aes_list <- utils::modifyList(base_aes, shape_aes)
  
  # Re-class the merged list so ggplot recognizes it as an aesthetic
  class(final_aes_list) <- "uneval"
  return(final_aes_list)
}


#=============================================================================#
# SECTION 4: CONFIGURATION VALIDATION
#=============================================================================#

#' Validate Configuration File
#'
#' @description Performs a "pre-flight check" on the loaded config list to
#' catch common errors before the pipeline runs.
#' FIXED: Correctly handles interaction terms (e.g. "genotype * treatment")
#'
#' @param config The config list read from yaml::read_yaml().
#' @param check_external_files A logical flag. If TRUE, it will check that
#' the `external_data_path` file(s) exist. This is set to TRUE for
#' Run 2 of Script 01 and for Script 02.
#'
#' @return Invisibly returns TRUE if successful, or stops with an error.
validate_config <- function(config, check_external_files = FALSE) {
  
  # Helper function for consistent error messages
  val_stop <- function(...) {
    msg <- paste(..., collapse = "")
    stop(paste("\n\n-- CONFIG VALIDATION FAILED ----------------\n\n", msg,
               "\n\n--------------------------------------------\n"),
         call. = FALSE)
  }
  
  message("--- Validating config.yml ---")
  
  # --- 1. Validate `paths` and `platemap` ---
  if (!dir.exists(here::here(config$paths$fsa_folder))) {
    val_stop("paths:fsa_folder: Directory not found.\n",
             "Checked: ", here::here(config$paths$fsa_folder))
  }
  
  if (is.null(config$paths$output_dir_base)) {
    val_stop("paths:output_dir_base: This field cannot be null.")
  }
  
  if (is.null(config$platemap$format) || !config$platemap$format %in% c("stacked_excel", "tidy_excel", "tidy_csv")) {
    val_stop("platemap:format: Must be one of 'stacked_excel', 'tidy_excel', or 'tidy_csv'.")
  }
  
  # Try to find platemap files
  platemap_files <- tryCatch({
    get_file_paths(config$platemap$path, "\\.(csv|xlsx|xls)$")
  }, error = function(e) { val_stop("platemap:path: ", e$message) })
  
  if (length(platemap_files) == 0) {
    val_stop("platemap:path: No platemap files were found at the specified path.")
  }
  
  # --- 2. Validate external data path (if required) ---
  if (check_external_files) {
    ext_files <- tryCatch({
      get_file_paths(config$paths$external_data_path, "\\.(xlsx|xls)$")
    }, error = function(e) { val_stop("paths:external_data_path: ", e$message) })
    
    if (length(ext_files) == 0) {
      val_stop("paths:external_data_path: No external data files were found.\n",
               "This is required for Run 2 of Script 01 and for Script 02.")
    }
    message(paste("...found", length(ext_files), "external data file(s)."))
  } else {
    message("...skipping external data file check (for Script 01, Run 1).")
  }
  
  # --- 3. Validate platemap columns ---
  message("...checking platemap columns.")
  
  loader_func <- switch(config$platemap$format,
                        "stacked_excel" = function(path) load_stacked_platemap_excel(path, config),
                        "tidy_csv" = load_tidy_platemap,
                        "tidy_excel" = load_tidy_platemap)
  
  sample_platemap <- tryCatch({
    loader_func(platemap_files[1]) # Just load first file to check structure
  }, error = function(e) { val_stop("platemap: Could not load sample platemap.\n", e$message) })
  
  required_platemap_cols <- c("plate", "well", "sample_name")
  if (!all(required_platemap_cols %in% colnames(sample_platemap))) {
    val_stop("platemap: Platemap file is missing required columns.\n",
             "Must contain: 'plate', 'well', and 'sample_name'.\n",
             "Found: ", paste(colnames(sample_platemap), collapse = ", "))
  }
  
  # --- 4. Validate `sample_name_columns` and `key_variables` ---
  message("...checking key variable definitions.")
  sample_cols <- config$sample_name_columns
  if (is.null(sample_cols)) {
    val_stop("sample_name_columns: This section must be defined.")
  }
  
  key_vars <- config$key_variables
  check_var <- function(var_name, var_value) {
    # It's okay for optional vars to be null
    if (is.null(var_value) || var_value == 'null') return(TRUE)
    
    # If it's not null, it must be in sample_cols
    if (!var_value %in% sample_cols) {
      val_stop("key_variables:", var_name, " ('", var_value, "') is not listed in 'sample_name_columns'.")
    }
  }
  
  # Check required vars
  check_var("time_variable", key_vars$time_variable)
  check_var("primary_group_var", key_vars$primary_group_var)
  
  # Check optional vars
  check_var("secondary_group_var", key_vars$secondary_group_var)
  check_var("optional_grouping_var", key_vars$optional_grouping_var)
  check_var("repeated_measure_var", key_vars$repeated_measure_var)
  check_var("optional_random_effect", key_vars$optional_random_effect)
  check_var("optional_crossed_effect", key_vars$optional_crossed_effect)
  
  # --- CHECK MODEL FIXED EFFECTS (FIXED) ---
  fx_effects_raw <- unlist(key_vars$model_fixed_effects)
  
  if (length(fx_effects_raw) > 0) {
    # Split interaction terms (e.g. "genotype * treatment" -> "genotype", "treatment")
    fx_effects_split <- unlist(strsplit(fx_effects_raw, "[\\*:]"))
    fx_effects_clean <- trimws(fx_effects_split)
    
    # Check if these individual parts exist in columns
    missing_fx <- setdiff(fx_effects_clean, sample_cols)
    
    if (length(missing_fx) > 0) {
      val_stop("key_variables:model_fixed_effects: The following variables are not listed in 'sample_name_columns':\n",
               paste(missing_fx, collapse = ", "), "\n",
               "Note: Ensure you list base columns (e.g. 'genotype'), not just the interaction term.")
    }
  }
  
  # --- CHECK THAT MODEL VARS ARE IN GROUPING VARS ---
  message("...verifying model effects survive data summarization.")
  
  # 1. Define the pool of variables that survive the Script 01 summarization
  grouping_pool <- unique(c(
    key_vars$time_variable,
    key_vars$primary_group_var,
    key_vars$secondary_group_var,
    key_vars$optional_grouping_var,
    key_vars$repeated_measure_var,
    key_vars$optional_crossed_effect
  ))
  grouping_pool <- grouping_pool[!is.null(grouping_pool) & grouping_pool != 'null']
  
  # 2. Extract all random effects required by the model
  re_effects_raw <- c()
  if (!is.null(key_vars$nesting_structure) && key_vars$nesting_structure != 'null' && key_vars$nesting_structure != "") {
    re_effects_raw <- trimws(unlist(strsplit(key_vars$nesting_structure, "[/:]")))
  } else {
    re_effects_raw <- c(key_vars$repeated_measure_var, key_vars$optional_random_effect)
  }
  re_effects_clean <- unique(re_effects_raw[!is.null(re_effects_raw) & re_effects_raw != 'null'])
  
  # 3. Combine fixed and random effects (fx_effects_clean is calculated just above this in your script)
  fx_effects_clean_safe <- if (exists("fx_effects_clean")) fx_effects_clean else c()
  all_model_vars <- unique(c(fx_effects_clean_safe, re_effects_clean))
  
  # 4. Check if any model variable is missing from the grouping pool
  missing_from_pool <- setdiff(all_model_vars, grouping_pool)
  
  if (length(missing_from_pool) > 0) {
    val_stop("key_variables: The following modeling variables are NOT assigned to a grouping role:\n",
             paste(missing_from_pool, collapse = ", "), "\n\n",
             "Any variable used in 'model_fixed_effects', 'optional_random_effect', or 'nesting_structure' ",
             "MUST also be assigned to a structural role (e.g., primary_group_var, secondary_group_var, ",
             "optional_grouping_var, optional_crossed_effect, or repeated_measure_var) ",
             "so it is not dropped when the data is summarized.")
  }
  
  # --- 5. Validate `parameters` ---
  bl_level <- config$parameters$baseline_grouping_level
  if (!bl_level %in% c("group", "genotype", "global_control", "custom")) {
    val_stop("parameters:baseline_grouping_level: Must be one of 'group', 'genotype', 'custom' or 'global_control'.")
  }
  
  if (bl_level == "global_control") {
    if (is.null(config$parameters$global_control_group$variable) || is.null(config$parameters$global_control_group$level)) {
      val_stop("parameters:global_control_group: 'variable' and 'level' must be set when using 'global_control' baseline.")
    }
  }
  
  # --- 6. Validate `response_variables` ---
  if (is.null(config$response_variables) || is.null(names(config$response_variables))) {
    val_stop("response_variables: This section is missing or malformed.")
  }
  
  if (!"mode" %in% names(config$response_variables)) {
    val_stop("response_variables: Must include an entry for 'mode'.\n",
             "Example: mode: 'target_peak_rpt'\n",
             "(This is required by the normalization function.)")
  }
  
  if (is.null(config$response_variable_shortnames)) {
    message("...Warning: 'response_variable_shortnames' is not defined. File names for plots may be very long.")
  }
  
  message("--- Validation Passed Successfully ---")
  return(invisible(TRUE))
}

#=============================================================================#
# SECTION 5: TRACE RECONSTRUCTION (NEW)
#=============================================================================#

#' Reconstruct Trace from Peak Table
#'
#' @description Takes a GeneMapper peak table (Size, Height, Area) and 
#' reconstructs a synthetic Gaussian trace for plotting.
#'
#' @param peak_data A dataframe containing peaks for ONE sample.
#' @param min_x Minimum x-axis value (CAG units).
#' @param max_x Maximum x-axis value (CAG units).
#' @param resolution Points per unit (default 10).
#' @param width_factor Controls the width of the peaks (visual smoothing).
#' @return A dataframe with columns: x (CAG), y (Height).
reconstruct_trace <- function(peak_data, min_x, max_x, resolution = 10, width_factor = 1.5) {
  
  # Create X sequence
  x_seq <- seq(min_x, max_x, by = 1/resolution)
  y_seq <- numeric(length(x_seq))
  
  # For each peak, add a Gaussian curve to the Y sequence
  # We assume 'Size' in peak_data is already converted to CAG units if needed,
  # OR we convert it here. The prompt implies input is bp, x is size/3.
  
  for (i in 1:nrow(peak_data)) {
    peak_center <- peak_data$Size[i] / 3 # Convert bp to CAG
    peak_height <- peak_data$Height[i]
    
    # Estimate Sigma from Area/Height ratio or use fixed visual width
    # A standard Gaussian: Height * exp(-0.5 * ((x - mu) / sigma)^2)
    # Visual width factor allows user to make peaks "chunkier" for visibility
    sigma <- width_factor * 0.3 # 0.3 is a base scaling factor
    
    # Calculate Gaussian
    curve <- peak_height * exp(-0.5 * ((x_seq - peak_center) / sigma)^2)
    
    # Add to total (superposition)
    y_seq <- y_seq + curve
  }
  
  return(data.frame(x = x_seq, y = y_seq))
}

#' Load and Process GeneMapper Peaks
#'
#' @description Loads the raw text file, cleans names, and nests peaks by sample.
#'
#' @param file_path Path to the exported text file.
#' @param sample_name_parser Function or logic to parse 'Sample File Name'.
#' @return A nested tibble with columns: sample_id, peaks (nested df).
load_genemapper_peaks <- function(file_path, config) {
  
  raw_df <- read.table(file_path, header = TRUE, sep = "\t", stringsAsFactors = FALSE, quote = "\"", fill = TRUE)
  
  # Basic Cleanup
  clean_df <- raw_df %>%
    dplyr::select(Sample_File = `Sample.File.Name`, Size, Height, Area) %>%
    dplyr::filter(!is.na(Size), !is.na(Height)) %>%
    # Parse the filename to get metadata matches (This relies on Script 01 logic)
    dplyr::mutate(fsa_filename = basename(Sample_File))
  
  return(clean_df)
}

# --- CUSTOM LEGEND GLYPHS (Side-by-Side Line & Point) ---
draw_key_short_line <- function(data, params, size) {
  lty <- if (!is.null(data$linetype)) data$linetype else 1
  col <- if (!is.null(data$colour)) data$colour else "black"
  lwd <- if (!is.null(data$linewidth)) data$linewidth else 1
  
  # Draws the line from 5% to 60% of the box width
  grid::linesGrob(
    x = c(0.05, 0.6), y = c(0.5, 0.5),
    gp = grid::gpar(col = col, lwd = lwd * ggplot2::.pt, lty = lty)
  )
}

draw_key_shifted_point <- function(data, params, size) {
  shape <- if (!is.null(data$shape)) data$shape else 19
  col <- if (!is.null(data$colour)) data$colour else "black"
  sz <- if (!is.null(data$size)) data$size else 3
  
  # Draws the point at 85% of the box width (next to the line)
  grid::pointsGrob(
    x = 0.85, y = 0.5, pch = shape,
    gp = grid::gpar(col = col, fill = col, fontsize = sz * ggplot2::.pt)
  )
}

# --- UNIVERSAL FACTOR RESTORER ---
# This function looks at the "master" modeling_data and forces any incoming 
# dataframe to adopt the exact same factor levels for all overlapping columns.
restore_all_factors <- function(df, reference_df = modeling_data) {
  if (is.null(df) || nrow(df) == 0) return(df)
  
  shared_cols <- intersect(names(df), names(reference_df))
  
  for (col in shared_cols) {
    if (is.factor(reference_df[[col]])) {
      # Re-apply the exact levels from the reference dataset
      df[[col]] <- factor(as.character(df[[col]]), levels = levels(reference_df[[col]]))
    }
  }
  return(df)
}

# =============================================================================#
# UNIVERSAL COLOR APPLIER
# =============================================================================#
# Smartly routes the correct palette based on the variable being mapped.
# Uses unname() to map colors strictly by factor order, completely bypassing
# the need to translate raw names to publication labels!
apply_smart_palette <- function(p, color_var_name) {
  
  # 1. Determine which palette to use
  active_pal <- primary_palette
  if (color_var_name == "clone_rank") {
    active_pal <- custom_palette
  } else if (!is.null(config$key_variables$secondary_group_var) && color_var_name == config$key_variables$secondary_group_var) {
    active_pal <- secondary_palette
  }
  
  # 2. Apply it to the plot
  if (!is.null(active_pal)) {
    
    n_groups <- length(unique(p$data[[color_var_name]]))
    
    if (length(active_pal) < n_groups) {
      active_pal <- grDevices::colorRampPalette(active_pal)(n_groups)
    }
    
    p <- p +
      scale_color_manual(values = unname(active_pal)) +
      scale_fill_manual(values = unname(active_pal))
  }
  
  return(p)
}



# --- SYNC DERIVED LABELS TO PRIMARY VARIABLE ORDER ---
# Because bind_rows() strips factor attributes, label columns (like Genotype_Exp) 
# revert to alphabetical text. This forces them to perfectly match the 
# explicitly ordered primary variable.

sync_label_factors <- function(df, primary_col) {
  if (is.null(df) || nrow(df) == 0 || !primary_col %in% names(df)) return(df)
  
  label_cols <- intersect(c("Genotype_Pub", "Genotype_Exp"), names(df))
  
  for (lbl_col in label_cols) {
    correct_order <- df %>%
      dplyr::arrange(!!rlang::sym(primary_col)) %>%
      dplyr::pull(!!rlang::sym(lbl_col)) %>%
      unique() %>%
      as.character()
    
    df[[lbl_col]] <- factor(as.character(df[[lbl_col]]), levels = correct_order)
  }
  return(df)
}

#=============================================================================#
# SECTION 6: LABEL & FACTOR MANAGEMENT  (single source of truth)
#=============================================================================#

#' Resolve the display label column name for a given style
#'
#' @description Canonical, single implementation used by ALL scripts.
#' Reads the active style from config$label_parameters$active_label_style,
#' validates the column exists in `data`, and falls back to primary_group_var
#' if the label columns have not yet been created (e.g. early in Script 01).
#'
#' Column name map:
#'   "publication" -> "Genotype_Pub"   (clean display name only, e.g. "Wild-Type")
#'   "exploratory" -> "Genotype_Exp"   (name + CAG size, e.g. "Wild-Type (140Q)")
#'
#' @param style One of "auto", "publication", or "exploratory".
#' @param data  The data frame to check for column existence (required).
#' @param config The config list (required).
#' @return A single character string: the column name to use.
get_lbl_canonical <- function(style = "auto", data, config) {
  # 1. Resolve style
  if (style == "auto") {
    style <- config$label_parameters$active_label_style %||% "exploratory"
  }
  # 2. Map style -> column name
  col_name <- switch(style,
                     "publication" = "Genotype_Pub",
                     "exploratory" = "Genotype_Exp",
                     "Genotype_Exp"   # safe fallback for any unrecognised string
  )
  # 3. Validate: fall back to raw primary var if label columns don't exist yet
  primary_var <- config$key_variables$primary_group_var
  if (is.null(data) || !col_name %in% names(data)) {
    if (!is.null(data)) {
      message(paste0(
        "get_lbl(): column '", col_name, "' not found in data. ",
        "Falling back to primary variable '", primary_var, "'. ",
        "Ensure build_label_columns() has been called first."
      ))
    }
    return(primary_var)
  }
  return(col_name)
}

#' Build the Genotype_Pub / Genotype_Exp label columns on a data frame
#'
#' @description Creates two display-label factor columns on `data` and on every
#' additional data frame passed in `extra_dfs`.  This is the ONLY place
#' Genotype_Pub and Genotype_Exp are constructed.
#'
#' Order of operations (matters!):
#'   1. value_renaming is applied first so labels contain the final display name.
#'   2. config$label_formats templates are filled with the renamed geno name and
#'      the calculated average CAG baseline.
#'   3. Factor levels are set in the order dictated by:
#'      a. config$factor_levels (if defined for primary_group_var), or
#'      b. ascending average CAG baseline (auto-sort).
#'   4. The primary grouping column itself is also re-levelled to match.
#'
#' @param data          The primary data frame (typically normalized_data).
#' @param baselines     The baselines_table data frame (must have final_baseline col).
#' @param config        The loaded config list.
#' @param extra_dfs     Named list of additional data frames to apply labels to.
#' @return A named list with the same names as extra_dfs plus "primary" for `data`,
#'         each having Genotype_Pub and Genotype_Exp columns added/updated, and
#'         an additional element "genotype_map" holding the lookup table.
#'         
#'  # build_label_columns() (defined in functions.R) is the SINGLE place that:
#   1. Applies value_renaming to the primary grouping column FIRST
#   2. Calculates average CAG baseline per group (for label templates)
#   3. Determines display order (config factor_levels OR ascending CAG)
#   4. Builds Genotype_Pub and Genotype_Exp factor columns in that order
#   5. Re-levels the primary column itself to match
#
# By doing renaming BEFORE label construction, the display labels will always
# contain the final human-readable name (e.g. "Wild-Type (140Q)") rather than
# the internal code (e.g. "WT (140Q)").
build_label_columns <- function(data, baselines, config, extra_dfs = list()) {
  
  primary_var   <- config$key_variables$primary_group_var
  secondary_var <- config$key_variables$secondary_group_var
  sym_primary   <- rlang::sym(primary_var)
  
  # Define grouping variables (include secondary if it exists in data)
  map_group_vars <- primary_var
  if (!is.null(secondary_var) && secondary_var != "null" && secondary_var %in% names(data)) {
    map_group_vars <- c(map_group_vars, secondary_var)
  }
  
  # ── Step 1: Apply value_renaming BEFORE building labels ──
  rename_map <- NULL
  if (!is.null(config$value_renaming) && length(config$value_renaming) > 0) {
    rename_map <- unlist(config$value_renaming)
  }
  
  rename_col <- function(vec) {
    if (is.null(rename_map)) return(vec)
    was_factor <- is.factor(vec)
    old_levels <- if (was_factor) levels(vec) else NULL
    char_vec   <- as.character(vec)
    new_vec    <- dplyr::coalesce(unname(rename_map[char_vec]), char_vec)
    if (was_factor) {
      new_levels <- dplyr::coalesce(unname(rename_map[old_levels]), old_levels)
      return(factor(new_vec, levels = new_levels))
    }
    new_vec
  }
  
  apply_rename_all <- function(df) {
    if (is.null(df)) return(df)
    for (col in map_group_vars) {
      if (col %in% names(df)) df[[col]] <- rename_col(df[[col]])
    }
    df
  }
  
  data      <- apply_rename_all(data)
  baselines <- apply_rename_all(baselines)
  extra_dfs <- lapply(extra_dfs, apply_rename_all)
  
  # ── Step 2: Calculate baseline CAG and build alias dictionary ──
  genotype_map <- data %>%
    dplyr::filter(!is.na(!!sym_primary)) %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(map_group_vars))) %>%
    dplyr::summarise(
      avg_baseline = round(mean(final_baseline, na.rm = TRUE), 0),
      .groups = "drop"
    ) %>%
    dplyr::mutate(
      # Universal dictionary variables available inside config glue templates
      base      = avg_baseline,
      primary   = as.character(.data[[primary_var]]),
      geno      = as.character(.data[[primary_var]]),
      secondary = if (!is.null(secondary_var) && secondary_var %in% names(.)) as.character(.data[[secondary_var]]) else "",
      treat     = if (!is.null(secondary_var) && secondary_var %in% names(.)) as.character(.data[[secondary_var]]) else ""
    )
  
  # ── Step 3: Determine display hierarchy order ──
  genotype_map <- genotype_map %>%
    dplyr::arrange(
      if (!is.null(config$factor_levels[[primary_var]])) match(as.character(.data[[primary_var]]), config$factor_levels[[primary_var]]) else avg_baseline,
      if (!is.null(secondary_var) && secondary_var %in% names(.) && !is.null(config$factor_levels[[secondary_var]])) match(as.character(.data[[secondary_var]]), config$factor_levels[[secondary_var]]) else 1
    )
  
  # ── Step 4: Build label strings using config templates ──
  fmt_pub <- config$label_formats$publication %||% "{geno}"
  fmt_exp <- config$label_formats$exploratory %||% "{geno} ({base}Q)"
  
  genotype_map <- genotype_map %>%
    dplyr::mutate(
      Label_Pub = as.character(glue::glue_data(., fmt_pub)),
      Label_Exp = as.character(glue::glue_data(., fmt_exp))
    )
  
  logr::log_print("build_label_columns(): label lookup map created:")
  logr::log_print(as.data.frame(genotype_map[, c(map_group_vars, "Label_Pub", "Label_Exp")]))
  
  # ── Step 5: Attach label columns and set factor levels ──
  ordered_primary <- unique(as.character(genotype_map[[primary_var]]))
  ordered_pub     <- unique(genotype_map$Label_Pub)
  ordered_exp     <- unique(genotype_map$Label_Exp)
  
  attach_labels <- function(df) {
    if (is.null(df) || !primary_var %in% names(df)) return(df)
    
    # Determine which grouping keys actually exist in this specific dataframe
    join_keys <- intersect(map_group_vars, names(df))
    
    df %>%
      dplyr::select(-dplyr::any_of(c("Genotype_Pub", "Genotype_Exp", "avg_baseline", "Label_Pub", "Label_Exp"))) %>%
      dplyr::left_join(
        genotype_map %>% dplyr::select(dplyr::all_of(join_keys), avg_baseline, Label_Pub, Label_Exp) %>% dplyr::distinct(),
        by = join_keys
      ) %>%
      dplyr::mutate(
        !!sym_primary := factor(as.character(!!sym_primary), levels = ordered_primary),
        Genotype_Pub   = factor(Label_Pub, levels = ordered_pub),
        Genotype_Exp   = factor(Label_Exp, levels = ordered_exp)
      ) %>%
      dplyr::select(-avg_baseline, -Label_Pub, -Label_Exp)
  }
  
  result <- c(
    list(primary = attach_labels(data),
         baselines = attach_labels(baselines),
         genotype_map = genotype_map),
    lapply(extra_dfs, attach_labels)
  )
  
  return(result)
}

#' Apply value_renaming AND factor levels to all columns of a data frame
#'
#' @description General-purpose helper used in Script 02 to refresh all
#' character/factor columns after loading from .RData.  Does NOT touch
#' Genotype_Pub / Genotype_Exp — those are managed by build_label_columns().
#'
#' @param df        The data frame to modify.
#' @param config    The loaded config list.
#' @return The modified data frame.
apply_renaming_and_factors <- function(df, config) {
  if (is.null(df) || nrow(df) == 0) return(df)
  
  rename_map <- NULL
  if (!is.null(config$value_renaming) && length(config$value_renaming) > 0) {
    rename_map <- unlist(config$value_renaming)
  }
  
  # ── Step A: rename values while preserving factor level order ──
  if (!is.null(rename_map)) {
    df <- df %>%
      dplyr::mutate(dplyr::across(
        # Skip the pre-built label columns — they are already correct
        dplyr::where(~ is.character(.) || is.factor(.)) &
          !dplyr::matches("^Genotype_(Pub|Exp)$"),
        ~ {
          was_factor <- is.factor(.)
          old_levels <- if (was_factor) levels(.) else NULL
          vec        <- as.character(.)
          new_vec    <- dplyr::coalesce(unname(rename_map[vec]), vec)
          if (was_factor) {
            new_levels <- dplyr::coalesce(unname(rename_map[old_levels]), old_levels)
            factor(new_vec, levels = new_levels)
          } else {
            new_vec
          }
        }
      ))
  }
  
  # ── Step B: re-apply config factor_levels (translated through rename_map) ──
  if (!is.null(config$factor_levels)) {
    for (col_name in names(config$factor_levels)) {
      if (!col_name %in% names(df)) next
      raw_levels <- as.character(config$factor_levels[[col_name]])
      ordered_levels <- if (!is.null(rename_map)) {
        dplyr::coalesce(unname(rename_map[raw_levels]), raw_levels)
      } else {
        raw_levels
      }
      df[[col_name]] <- factor(as.character(df[[col_name]]), levels = ordered_levels)
    }
  }
  
  df
}

generate_rounded_peaks <- function(df, width_cag = 0.45, resolution = 51) {
  if(nrow(df) == 0) return(df)
  
  # 1. High-resolution grid for perfectly smooth drawing
  bell_curve <- data.frame(
    offset = seq(-width_cag, width_cag, length.out = resolution)
  )
  
  # 2. Gaussian curve math (Authentic fragment peak shape)
  # Setting sigma to width/3 ensures the tails smoothly hit 0 at the edges
  sigma <- width_cag / 3
  bell_curve$mult <- exp(-0.5 * (bell_curve$offset / sigma)^2)
  
  # 3. Apply this shape to every single row in the data
  df %>%
    merge(bell_curve, by = NULL) %>% 
    dplyr::mutate(
      CAG = CAG + offset,
      Height = Height * mult
    ) %>%
    dplyr::arrange(CAG) %>%
    dplyr::select(-offset, -mult)
}

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


