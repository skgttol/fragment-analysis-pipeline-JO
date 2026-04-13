#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
#           Advanced Fragment Analysis Workflow - Script 5                    #
#    (PDF Reports + Publication Figures | Zoomed Summaries & Smoothing)       #
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#=============================================================================#
# PART 0: SETUP & PATHS
#=============================================================================#
# Packages
packages <- c("tidyverse", "here", "yaml", "ggridges", "cowplot", "logr", "magick", "grid", "data.table", "writexl", "ggh4x")
installed_packages <- packages %in% rownames(utils::installed.packages())
if (any(installed_packages == FALSE)) {
  utils::install.packages(packages[!installed_packages])
}
invisible(lapply(packages, library, character.only = TRUE))

func_path <- "C:/Users/skgttol/OneDrive - University College London/PhD/PhD_Thesis/02_Data/Template_RScripts/CAGsizing_Flexible/functions.R"
if (!file.exists(func_path)) stop("CRITICAL ERROR: 'functions.R' not found.")
source(func_path)

config <- yaml::read_yaml(here::here("config.yml"))

# --- Check Multiple Input Files ---
raw_paths_config <- config$trace_settings$genemapper_file_path
if (is.list(raw_paths_config)) raw_paths_config <- unlist(raw_paths_config)
combined_file_paths <- here::here(raw_paths_config)

if (any(!file.exists(combined_file_paths))) {
  missing_files <- combined_file_paths[!file.exists(combined_file_paths)]
  stop(paste("CRITICAL ERROR: The following data files were not found:\n", 
             paste(missing_files, collapse = "\n")))
}

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

load(file.path(latest_analysis_dir, "processing_complete.RData"))

config <- yaml::read_yaml(here::here("config.yml"))

trace_output_dir <- file.path(latest_analysis_dir, "trace_plots")
dir.create(trace_output_dir, showWarnings = FALSE)

# --- Organized Sub-folders ---
pdf_base_dir   <- file.path(trace_output_dir, "pdf_reports") 
tiff_dir       <- file.path(trace_output_dir, "publication_figures") 
excel_dir      <- file.path(trace_output_dir, "data_tables") # New folder for Excel

# Specific PDF Folders
pdf_dir_exp    <- file.path(pdf_base_dir, "Expanded")
pdf_dir_wt     <- file.path(pdf_base_dir, "WT")
pdf_dir_target <- file.path(pdf_base_dir, "Target")
pdf_dir_smooth <- file.path(pdf_base_dir, "Smoothed_Inspection")

# Create all directories recursively
dir.create(tiff_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(excel_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(pdf_dir_exp, showWarnings = FALSE, recursive = TRUE)
dir.create(pdf_dir_wt, showWarnings = FALSE, recursive = TRUE)
dir.create(pdf_dir_target, showWarnings = FALSE, recursive = TRUE)
dir.create(pdf_dir_smooth, showWarnings = FALSE, recursive = TRUE)


try(logr::log_close(), silent = TRUE)
lf <- logr::log_open(file.path(trace_output_dir, "trace_log.log"), show_notes = FALSE, logdir = FALSE)
logr::log_print("\n--- SCRIPT 05: FINAL VISUALIZATION ---", console = TRUE)

get_lbl <- function(style = "auto") {
  data_to_check <- if (exists("modeling_data")) modeling_data else NULL
  get_lbl_canonical(style = style, data = data_to_check, config = config)
}

#=============================================================================#
# PART 1: LOAD & OPTIMIZE RAW DATA
#=============================================================================#
logr::log_print("Loading and Optimizing Data...", console=TRUE)

# --- 1. CONFIG VARIABLES ---
flank  <- as.numeric(config$custom_settings$flank)
corr   <- as.numeric(config$custom_settings$correction)
limits_exp <- config$trace_settings$x_limits
limits_wt  <- config$trace_settings$wt_x_limits

rfu_threshold <- as.numeric(config$parameters$peak_height_threshold)
if(is.null(rfu_threshold)) rfu_threshold <- 50 

col_mode <- "mode"
col_ii   <- "instability_index" 

if(is.null(flank) | is.null(corr) | is.null(limits_exp)) stop("CRITICAL ERROR: Primary Config variables missing.")

# 2. Fast Load
raw_file_paths <- config$trace_settings$genemapper_file_path
if (is.list(raw_file_paths)) raw_file_paths <- unlist(raw_file_paths)

raw_trace_data <- purrr::map_df(raw_file_paths, function(fp) {
  dt <- data.table::fread(here::here(fp), fill = TRUE, data.table = FALSE)
  colnames(dt) <- make.names(colnames(dt))
  return(dt)
})

# 3. Identify Sample Column
possible_cols <- c("Sample.File.Name", "Sample.Name", "File.Name", "fsa_filename")
join_col <- intersect(names(raw_trace_data), possible_cols)[1]
if (is.na(join_col)) stop("ERROR: Could not find Sample Name column in the combined data.")

logr::log_print(paste("Successfully loaded", nrow(raw_trace_data), "rows combined."), console = TRUE)

# 4. PRE-PROCESSING
min_keep <- limits_exp[1]
max_keep <- limits_exp[2]
if(!is.null(limits_wt)) {
  min_keep <- min(min_keep, limits_wt[1])
  max_keep <- max(max_keep, limits_wt[2])
}
buffer <- 10

optimized_data <- raw_trace_data %>%
  dplyr::select(fsa_filename = all_of(join_col), Size, Height) %>%
  dplyr::mutate(across(c(Size, Height), as.numeric)) %>%
  dplyr::filter(!is.na(Size) & !is.na(Height)) %>%
  dplyr::mutate(CAG = ((Size - flank) / 3) + corr) %>%
  dplyr::filter(CAG >= (min_keep - buffer) & CAG <= (max_keep + buffer))

# 5. CREATE LOOKUP LIST
trace_lookup <- split(optimized_data, optimized_data$fsa_filename)
logr::log_print("Data optimized.", console=TRUE)


#=============================================================================#
# PART 2: METADATA MERGE & FACTOR ORDERING
#=============================================================================#
cols_to_keep <- c("fsa_filename", 
                  "pcr", "rep",
                  col_mode, col_ii, "ii_threshold_abs", 
                  "target_start", "target_end",        
                  config$key_variables$primary_group_var, 
                  config$key_variables$time_variable, 
                  config$key_variables$secondary_group_var,
                  config$key_variables$optional_grouping_var,
                  config$key_variables$optional_crossed_effect)

cols_to_keep <- cols_to_keep[!sapply(cols_to_keep, is.null)]
cols_to_keep <- cols_to_keep[cols_to_keep != "null"]
cols_to_keep <- unique(cols_to_keep)

meta_clean <- processed_data %>% 
  dplyr::select(any_of(cols_to_keep)) %>%
  mutate(is_excluded = FALSE, exclusion_reason = NA_character_) %>%
  mutate(!!sym(config$key_variables$time_variable) := as.numeric(!!sym(config$key_variables$time_variable)))

if (exists("excluded_data") && nrow(excluded_data) > 0) {
  meta_excl <- excluded_data %>%
    dplyr::select(any_of(c(cols_to_keep, "exclusion_reason"))) %>%
    mutate(is_excluded = TRUE) %>%
    mutate(!!sym(config$key_variables$time_variable) := as.numeric(!!sym(config$key_variables$time_variable)))
  
  common_cols <- intersect(names(meta_clean), names(meta_excl))
  all_metadata <- bind_rows(
    meta_clean %>% dplyr::select(any_of(common_cols)),
    meta_excl %>% dplyr::select(any_of(common_cols))
  )
} else {
  all_metadata <- meta_clean
}

all_metadata <- apply_renaming_and_factors(all_metadata, config)

# --- ROBUST BIO_REP_ID CREATION ---
potential_id_cols <- c(
  config$key_variables$primary_group_var, 
  config$key_variables$secondary_group_var, 
  config$key_variables$optional_grouping_var,
  config$key_variables$optional_crossed_effect
)

potential_id_cols <- unique(potential_id_cols)

valid_id_cols <- potential_id_cols[
  !is.null(potential_id_cols) & potential_id_cols != "null" & potential_id_cols %in% colnames(all_metadata)
]
peaks_annotated <- all_metadata %>%
  tidyr::unite("bio_id", dplyr::all_of(valid_id_cols), sep = "_", remove = FALSE)

if(!"pcr" %in% names(peaks_annotated)) {
  peaks_annotated <- peaks_annotated %>% mutate(pcr = as.numeric(str_extract(fsa_filename, "(?<=_)\\d+(?=.fsa)")))
}

plot_label_col <- get_lbl("publication") # Get the column name (e.g., "Genotype_Pub")
logr::log_print(paste("Using Fancy Label Column:", plot_label_col))

# B. Join Fancy Labels to peaks_annotated
# We check if 'modeling_data' (loaded from RData) has the map.
if(exists("modeling_data") && plot_label_col %in% names(modeling_data)) {
  
  # 1. Create unique map
  label_map <- modeling_data %>% 
    dplyr::distinct(!!sym(config$key_variables$primary_group_var), !!sym(plot_label_col))
  
  # 2. Join to peaks_annotated
  peaks_annotated <- peaks_annotated %>%
    left_join(label_map, by = config$key_variables$primary_group_var)
  
  # 3. Apply Factor Ordering (Crucial for Plots)
  # This ensures "HD 100" doesn't appear before "HD 20" alphabetically
  if(is.factor(peaks_annotated[[config$key_variables$primary_group_var]])) {
    # Get the order of the raw levels
    raw_levels <- levels(peaks_annotated[[config$key_variables$primary_group_var]])
    
    # Arrange map by that order to get the matching fancy order
    ordered_map <- label_map %>%
      mutate(raw_fct = factor(!!sym(config$key_variables$primary_group_var), levels = raw_levels)) %>%
      arrange(raw_fct)
    
    fancy_levels <- unique(ordered_map[[plot_label_col]])
    
    # Apply levels
    peaks_annotated[[plot_label_col]] <- factor(peaks_annotated[[plot_label_col]], levels = fancy_levels)
  }
} else {
  # Fallback if map missing: Create a dummy column matching primary var
  peaks_annotated[[plot_label_col]] <- peaks_annotated[[config$key_variables$primary_group_var]]
}

# ---------------------------------------------------------

baseline_table <- peaks_annotated %>%
  dplyr::filter(!is_excluded) %>% 
  dplyr::filter(!!sym(config$key_variables$time_variable) == min(!!sym(config$key_variables$time_variable), na.rm=TRUE)) %>%
  group_by(bio_id) %>%
  summarise(baseline_cag = mean(as.numeric(!!sym(col_mode)), na.rm = TRUE), .groups="drop")


#=============================================================================#
# PART 3: PDF REPORTS (Organized & Vectorized with facet_grid)
#=============================================================================#
logr::log_print("Generating PDF Reports...", console=TRUE)

# --- 1. Join Target Metadata ---
cols_to_join <- c("fsa_filename", "target_start", "target_end")
cols_to_join <- intersect(cols_to_join, colnames(processed_data))
if(length(cols_to_join) > 1 && !"target_start" %in% colnames(peaks_annotated)) {
  # Line 244 — change to:
  meta_lookup <- processed_data %>% 
    dplyr::select(dplyr::all_of(cols_to_join)) %>% 
    dplyr::distinct(fsa_filename, .keep_all = TRUE)  # ← force one row per file
  
  peaks_annotated <- peaks_annotated %>% dplyr::left_join(meta_lookup, by="fsa_filename")
}

# --- 2. Update bio_rep_id for Reports ---
possible_rep_cols <- c("rep", "bio_rep", "replicate", "well_id", "clone_rep")
found_rep_col <- intersect(possible_rep_cols, colnames(peaks_annotated))
rep_col <- if(length(found_rep_col) > 0) found_rep_col[1] else NULL

if (!is.null(rep_col)) {
  # If a replicate column is found, combine bio_id and the replicate number
  peaks_annotated <- peaks_annotated %>% 
    tidyr::unite("bio_rep_id", dplyr::all_of(c("bio_id", rep_col)), sep = "_", remove = FALSE)
} else {
  # Fallback: If no replicate column exists in the metadata, just mirror bio_id
  peaks_annotated <- peaks_annotated %>% 
    dplyr::mutate(bio_rep_id = bio_id)
}
# Run after line 255 (second bio_rep_id creation)
dup_check <- peaks_annotated %>% 
  group_by(fsa_filename) %>% 
  dplyr::filter(n() > 1)
cat("Duplicate fsa_filename rows:", nrow(dup_check), "\n")
if(nrow(dup_check) > 0) print(head(dup_check %>% select(fsa_filename, bio_rep_id, pcr, everything())))

# --- 3. Helper Functions ---
triangulate_data <- function(df, width_cag = 0.4) {
  if(nrow(df) == 0) return(df)
  bind_rows(df %>% mutate(CAG=CAG-width_cag, Height=0), df, df %>% mutate(CAG=CAG+width_cag, Height=0)) %>% arrange(CAG)
}

plot_modes <- list(list(name="Expanded", limits=limits_exp))
if(!is.null(limits_wt)) plot_modes[[length(plot_modes)+1]] <- list(name="WT", limits=limits_wt)
plot_modes[[length(plot_modes)+1]] <- list(name="Target", limits=NULL)

unique_genos <- unique(peaks_annotated[[config$key_variables$primary_group_var]])
unique_genos <- unique_genos[!is.na(unique_genos)]
time_var <- config$key_variables$time_variable

for (geno in unique_genos) {
  geno_data <- peaks_annotated %>% dplyr::filter(!!sym(config$key_variables$primary_group_var) == geno)
  unique_bioreps <- unique(geno_data$bio_rep_id)
  
  if(length(unique_bioreps) == 0) next
  logr::log_print(paste0("Genotype: ", geno, " | Pages: ", length(unique_bioreps)), console=TRUE)
  
  for(mode in plot_modes) {
    mode_name <- mode$name
    static_limits <- mode$limits
    
    curr_save_dir <- switch(mode_name,
                            "Expanded" = pdf_dir_exp,
                            "WT"       = pdf_dir_wt,
                            "Target"   = pdf_dir_target,
                            pdf_base_dir)
    
    safe_geno <- gsub("[^A-Za-z0-9_]", "_", geno)
    fname <- file.path(curr_save_dir, paste0("TraceReport_", mode_name, "_", safe_geno, ".pdf"))
    
    # OPEN PDF ONCE PER FILE
    pdf(file = fname, width = 8.5, height = 11)
    pb <- txtProgressBar(min = 0, max = length(unique_bioreps), style = 3, char = "=")
    
    for (i in seq_along(unique_bioreps)) {
      biorep <- unique_bioreps[i]
      rep_data <- geno_data %>% dplyr::filter(bio_rep_id == biorep)
      
      # --- PART 3 UPDATE: Dynamic Header for Batch ---
      curr_clone <- unique(rep_data[[config$key_variables$secondary_group_var]])[1]
      
      # Determine if a Batch exists
      batch_var_name <- config$key_variables$optional_crossed_effect
      curr_batch <- if(!is.null(batch_var_name) && batch_var_name %in% names(rep_data)) {
        unique(rep_data[[batch_var_name]])[1]
      } else ""
      
      curr_opt <- if(!is.null(config$key_variables$optional_grouping_var) && 
                     config$key_variables$optional_grouping_var %in% names(rep_data)) {
        unique(rep_data[[config$key_variables$optional_grouping_var]])[1]
      } else unique(rep_data[[found_rep_col]])[1]
      
      label_parts <- c()
      if(curr_batch != "") label_parts <- c(label_parts, paste("Batch:", curr_batch))
      if(curr_opt != "")   label_parts <- c(label_parts, paste("Rep:", curr_opt)) 
      
      label_final <- if(length(label_parts) > 0) paste(label_parts, collapse=" | ") else "Rep: 1"
      
      # 1. Extract the parent bio_id for the current replicate
      curr_bio_id <- unique(rep_data$bio_id)[1]
      
      # 2. Look up the baseline using bio_id instead of bio_rep_id
      base_val <- if(mode_name != "WT") { 
        val <- baseline_table %>% 
          dplyr::filter(bio_id == curr_bio_id) %>% 
          pull(baseline_cag)
        
        if(length(val) > 0) val[1] else NA 
      } else NA      
      # Assemble Page Data
      page_traces <- list()
      page_meta <- list()
      
      safe_time_vec <- as.numeric(as.character(rep_data[[time_var]]))
      unique_days <- sort(unique(na.omit(safe_time_vec)))
      
      # Dynamically find the maximum number of PCRs (fallback to 3 if missing)
      max_pcr <- max(as.numeric(as.character(rep_data$pcr)), na.rm = TRUE)
      if (is.infinite(max_pcr) || is.na(max_pcr)) max_pcr <- 3
      
      grid_skeleton <- expand.grid(Day = unique_days, PCR = 1:max_pcr, stringsAsFactors = FALSE)
      
      for (row_idx in 1:nrow(grid_skeleton)) {
        # FIX: Use unique variable names to prevent dplyr masking
        curr_day <- grid_skeleton$Day[row_idx]
        curr_pcr <- grid_skeleton$PCR[row_idx]
        
        # 2. BULLETPROOF FIX: Use Base R to find the row (Immune to dplyr masking)
        safe_pcr_vec <- as.numeric(as.character(rep_data$pcr))
        match_idx <- which(safe_time_vec == curr_day & safe_pcr_vec == curr_pcr)
        
        if(length(match_idx) > 0) {
          match_row <- rep_data[match_idx[1], , drop = FALSE]
        } else {
          match_row <- rep_data[0, , drop = FALSE] # Empty dataframe
        }
        
        meta <- data.frame(Day = curr_day, PCR = curr_pcr, mode_val = NA, ii_val = NA, ii_thresh = NA, w_min = NA, w_max = NA, 
                           is_excluded = FALSE, has_data = FALSE, center_x = NA, excl_reason = "")
        
        if (nrow(match_row) > 0) {
          meta$has_data <- TRUE
          meta$is_excluded <- match_row$is_excluded
          if ("exclusion_reason" %in% names(match_row) && !is.na(match_row$exclusion_reason)) {
            meta$excl_reason <- match_row$exclusion_reason
          } else if (match_row$is_excluded) {
            meta$excl_reason <- "Manual/Unknown Exclusion"
          }
          
          meta$mode_val <- if(col_mode %in% names(match_row)) as.numeric(match_row[[col_mode]]) else NA
          meta$ii_val <- if(col_ii %in% names(match_row)) as.numeric(match_row[[col_ii]]) else NA
          meta$ii_thresh <- if("ii_threshold_abs" %in% names(match_row)) as.numeric(match_row[["ii_threshold_abs"]]) else NA
          w_min_s <- if("target_start" %in% names(match_row)) as.numeric(match_row[["target_start"]]) else NA
          w_max_s <- if("target_end" %in% names(match_row)) as.numeric(match_row[["target_end"]]) else NA
          
          if (mode_name == "Target") {
            meta$w_min <- w_min_s; meta$w_max <- w_max_s
            use_min <- if(!is.na(w_min_s)) w_min_s - 5 else limits_exp[1]
            use_max <- if(!is.na(w_max_s)) w_max_s + 5 else limits_exp[2]
          } else {
            use_min <- static_limits[1]; use_max <- static_limits[2]
            if (mode_name == "Expanded") { meta$w_min <- w_min_s; meta$w_max <- w_max_s }
          }
          
          meta$center_x <- mean(c(use_min, use_max), na.rm = TRUE)
          
          clean_fname <- as.character(match_row$fsa_filename)  # ← explicit fix
          trace <- trace_lookup[[clean_fname]]
          if(is.null(trace)) trace <- trace_lookup[[sub("\\.fsa$", "", match_row$fsa_filename)]]
          
          if (!is.null(trace) && nrow(trace) > 0) {
            trace <- trace %>% dplyr::filter(CAG >= use_min & CAG <= use_max)
            if(nrow(trace) > 0) {
              trace_tri <- triangulate_data(trace)
              # FIX: Use curr_day and curr_pcr here too
              trace_tri$Day <- curr_day; trace_tri$PCR <- curr_pcr
              page_traces[[length(page_traces) + 1]] <- trace_tri
            } else { meta$has_data <- FALSE }
          } else { meta$has_data <- FALSE }
        } else {
          # Provide dummy center for missing data
          use_min <- if(mode_name == "Target") limits_exp[1] else static_limits[1]
          use_max <- if(mode_name == "Target") limits_exp[2] else static_limits[2]
          meta$center_x <- mean(c(use_min, use_max), na.rm = TRUE)
        }
        
        lbl <- paste0("Size: ", ifelse(is.na(meta$mode_val), "NA", round(meta$mode_val, 2)))
        if(!is.na(meta$ii_val)) lbl <- paste0(lbl, " | II: ", round(meta$ii_val, 2))
        meta$label <- lbl
        page_meta[[length(page_meta) + 1]] <- meta
      }
      
      # 3. FIX: Lock the Day_Label into a strictly ordered numerical factor
      ordered_day_labels <- paste(stringr::str_to_title(time_var), unique_days)
      
      df_meta <- dplyr::bind_rows(page_meta) %>%
        dplyr::mutate(
          PCR_Label = paste("PCR", PCR), 
          Day_Label = factor(paste(stringr::str_to_title(time_var), Day), levels = ordered_day_labels)
        )
      
      # Construct Plot
      p_page <- ggplot()
      
      # 2. Safely check for excluded data BEFORE adding the layer
      excluded_meta <- df_meta %>% dplyr::filter(is_excluded == TRUE)
      
      if(nrow(excluded_meta) > 0) {
        p_page <- p_page + 
          geom_rect(data = excluded_meta, 
                    aes(xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf), 
                    fill = "pink", alpha = 0.3)
      }
      
      p_page <- p_page +
        geom_text(data = df_meta %>% dplyr::filter(is_excluded == TRUE | has_data == FALSE),
                  aes(x = center_x, y = Inf, label = ifelse(has_data == FALSE, "No Data", paste0("EXCLUDED\n(", excl_reason, ")"))), 
                  color = "red", size = 2, fontface = "bold", vjust = 3, na.rm = TRUE) +
        geom_vline(data = df_meta, aes(xintercept = w_min), linetype = "dashed", color = "darkgreen", alpha = 0.5, na.rm = TRUE) +
        geom_vline(data = df_meta, aes(xintercept = w_max), linetype = "dashed", color = "darkgreen", alpha = 0.5, na.rm = TRUE) +
        geom_vline(data = df_meta, aes(xintercept = mode_val), linetype = "dotdash", color = "blue", alpha = 0.6, na.rm = TRUE) +
        geom_hline(data = df_meta, aes(yintercept = ii_thresh), linetype = "dashed", color = "purple", alpha = 0.7, na.rm = TRUE) +
        geom_hline(yintercept = rfu_threshold, linetype = "dotted", color = "red") +
        {if(!is.na(base_val)) geom_vline(xintercept = base_val, linetype = "dashed", color = "grey50")} +
        geom_text(data = df_meta, aes(x = -Inf, y = Inf, label = label), hjust = -0.05, vjust = 1.5, size = 2, fontface = "bold") +
        ggh4x::facet_grid2(
          Day_Label ~ PCR_Label, 
          scales = "free_y", 
          independent = "y"  # <--- This is the magic ggh4x argument for grids
        ) + 
        labs(title = paste0(geno, " | ", curr_clone, " | ", label_final, " | ", mode_name)) +
        theme_minimal() +
        theme(
          axis.title = element_blank(), 
          panel.grid = element_blank(), 
          panel.border = element_rect(colour = "grey80", fill = NA),
          strip.text = element_text(size = 8, face = "bold"), # The "Days" and "PCR" labels
          plot.title = element_text(size = 14, face = "bold", hjust = 0.5),
          
          # --- NEW FONT CONTROLS ---
          axis.text.x = element_text(size = 5), # Shrinks the CAG numbers (120, 140, 160)
          axis.text.y = element_text(size = 5), # Shrinks the RFU numbers (0, 500, 1000)
          axis.title.x = element_text(size = 7) # Shrinks the "CAG Length" title if you add one
        )
      
      if(length(page_traces) > 0) {
        df_traces <- dplyr::bind_rows(page_traces) %>% 
          dplyr::mutate(
            PCR_Label = paste("PCR", PCR), 
            Day_Label = factor(paste(stringr::str_to_title(time_var), Day), 
                               levels = levels(df_meta$Day_Label))  # ← use df_meta's levels explicitly
          )
        
        p_page <- p_page + geom_area(data = df_traces, aes(x = CAG, y = Height), fill = "blue", alpha = 0.2) +
          geom_line(data = df_traces, aes(x = CAG, y = Height), color = "blue", linewidth = 0.2)
      }
      
      if(mode_name != "Target") p_page <- p_page + coord_cartesian(xlim = c(static_limits[1], static_limits[2]))
      
      # Get one filename per Day/PCR combo from the traces
      trace_labels <- df_traces %>%
        dplyr::distinct(Day_Label, PCR_Label, fsa_filename)
      
      p_page <- p_page + 
        geom_text(data = trace_labels, 
                  aes(x = -Inf, y = Inf, label = fsa_filename), 
                  hjust = -0.01, vjust = 4, size = 1.5, color = "darkred", inherit.aes = FALSE)
      
      print(p_page)
      setTxtProgressBar(pb, i)
    }
    close(pb)
    dev.off()
  }
}


#=============================================================================#
# PART 4: PUBLICATION FIGURES
#=============================================================================#
try({ while(!is.null(dev.list())) dev.off() }, silent=TRUE)
logr::log_print("Generating Publication Figures (TIFFs)...", console=TRUE)
#use label_lookup as this has genotype_pub

plot_data_full <- optimized_data %>%
  left_join(peaks_annotated, by = "fsa_filename") %>% 
  dplyr::filter(!is_excluded) %>%
  dplyr::filter(!is.na(!!sym(config$key_variables$primary_group_var))) %>%
  dplyr::filter(!is.na(target_start) & !is.na(target_end)) %>% 
  dplyr::filter(CAG >= target_start & CAG <= target_end) %>%
  group_by(fsa_filename) %>%
  mutate(Rel_Height = Height / max(Height, na.rm = TRUE)) %>% 
  ungroup() %>% restore_all_factors()

plot_data_mean <- plot_data_full %>%
  group_by(!!sym(config$key_variables$primary_group_var), !!sym(config$key_variables$time_variable), CAG) %>%
  summarise(Mean_Rel_Height = mean(Rel_Height, na.rm = TRUE), .groups="drop") %>%
  left_join(label_lookup, by = config$key_variables$primary_group_var)

baseline_visual_peaks <- plot_data_mean %>%
  dplyr::filter(!!sym(config$key_variables$time_variable) == min(!!sym(config$key_variables$time_variable), na.rm=TRUE)) %>%
  group_by(!!sym(config$key_variables$primary_group_var)) %>%
  
  # 1. Identify the max height for this group
  mutate(Max_H = max(Mean_Rel_Height)) %>%
  
  # 2. dplyr::filter to keep ALL rows that match this max height (handling ties)
  dplyr::filter(Mean_Rel_Height == Max_H) %>%
  
  # 3. Summarise to get the Median CAG of those tied rows
  summarise(Center_CAG = median(CAG), .groups="drop") %>%
  
  dplyr::select(!!sym(config$key_variables$primary_group_var), Center_CAG) %>%
  left_join(label_lookup, by = config$key_variables$primary_group_var)


window_limits <- baseline_visual_peaks %>%
  mutate(
    Win_Min = case_when(Center_CAG < 36 ~ Center_CAG - 8, TRUE ~ Center_CAG - 15),
    Win_Max = case_when(Center_CAG < 36 ~ Center_CAG + 8, TRUE ~ Center_CAG + 30)
  )

start_end_raw <- plot_data_mean %>%
  group_by(!!sym(config$key_variables$primary_group_var)) %>%
  dplyr::filter(!!sym(config$key_variables$time_variable) == min(!!sym(config$key_variables$time_variable)) | 
           !!sym(config$key_variables$time_variable) == max(!!sym(config$key_variables$time_variable))) %>%
  mutate(Timepoint_Label = ifelse(
    !!sym(config$key_variables$time_variable) == min(!!sym(config$key_variables$time_variable)), 
    "Start", "End"
  )) %>%
  mutate(Timepoint_Label = factor(Timepoint_Label, levels = c("Start", "End"))) %>%
  ungroup()

aligned_data <- start_end_raw %>%
  inner_join(window_limits, by = config$key_variables$primary_group_var) %>%
  dplyr::filter(CAG >= Win_Min & CAG <= Win_Max) %>%
  group_by(!!sym(config$key_variables$primary_group_var), Timepoint_Label, CAG) %>%
  summarise(Mean_Rel_Height = mean(Mean_Rel_Height, na.rm=TRUE), .groups="drop") %>%
  group_by(!!sym(config$key_variables$primary_group_var)) %>%
  mutate(Facet_Norm_Height = Mean_Rel_Height / max(Mean_Rel_Height, na.rm=TRUE)) %>%
  ungroup() %>%
  left_join(label_lookup, by = config$key_variables$primary_group_var) %>% restore_all_factors()

force_width_df <- window_limits %>% tidyr::pivot_longer(cols = c(Win_Min, Win_Max), values_to = "CAG") %>% dplyr::select(!!sym(config$key_variables$primary_group_var), CAG)

pub_label_col <- get_lbl("publication") # Get the column name (e.g., "Genotype_Pub")
logr::log_print(paste("Using Fancy Label Column:", pub_label_col))

if(nrow(aligned_data) > 0) {
  p_aligned <- ggplot(aligned_data, aes(x = CAG)) +
    geom_area(aes(y = Facet_Norm_Height, fill = Timepoint_Label, group = Timepoint_Label), alpha = 0.5, position = "identity") +
    geom_vline(data = baseline_visual_peaks, aes(xintercept = Center_CAG), linetype = "dotted", color = "black", alpha = 0.8) +
    geom_blank(data = force_width_df, aes(x = CAG)) + 
    facet_wrap(vars(!!sym(pub_label_col)), ncol=1, scales="free") + 
    scale_fill_manual(values = c("Start" = "grey60", "End" = "red")) + 
    labs(title = "Total Shift: Start vs End (Aligned)", y = "Relative Height", fill = "Timepoint") +
    theme_cowplot()
  ggsave(file.path(tiff_dir, "Publication_Start_vs_End_ALIGNED.tiff"), p_aligned, width = 8, height = 12, compression = "lzw")
}

calc_smooth_curve <- function(cag, height, spar_val=0.4) {
  if(length(cag) < 8) return(data.frame(CAG=cag, Height=height))
  tryCatch({
    fit <- smooth.spline(x = cag, y = height, spar = spar_val)
    pred_x <- seq(min(cag), max(cag), length.out = 200)
    pred_y <- predict(fit, pred_x)$y; pred_y[pred_y < 0] <- 0 
    return(data.frame(CAG = pred_x, Height = pred_y))
  }, error = function(e) return(data.frame(CAG=cag, Height=height)))
}
smoothed_data <- aligned_data %>%
  group_by(!!sym(config$key_variables$primary_group_var), Timepoint_Label) %>%
  reframe(calc_smooth_curve(CAG, Facet_Norm_Height)) %>%
  group_by(!!sym(config$key_variables$primary_group_var)) %>%
  mutate(Norm_Height_Smooth = Height / max(Height, na.rm=TRUE)) %>%
  ungroup()

if(nrow(smoothed_data) > 0) {
  p_smooth_aligned <- ggplot(smoothed_data, aes(x = CAG)) +
    geom_area(aes(y = Norm_Height_Smooth, fill = Timepoint_Label), alpha = 0.4, position = "identity") +
    geom_line(aes(y = Norm_Height_Smooth, color = Timepoint_Label), linewidth=0.5) +
    geom_vline(data = baseline_visual_peaks, aes(xintercept = Center_CAG), linetype = "dotted") +
    geom_blank(data = force_width_df, aes(x = CAG)) + 
    facet_wrap(vars(!!sym(pub_label_col)), ncol=1, scales="free") + 
    scale_fill_manual(values = c("Start" = "grey60", "End" = "red")) + scale_color_manual(values = c("Start" = "grey40", "End" = "darkred")) + 
    labs(title = "Total Shift: Start vs End (Smoothed)", y = "Relative Height") + theme_cowplot()
  ggsave(file.path(tiff_dir, "Publication_Start_vs_End_SMOOTHED.tiff"), p_smooth_aligned, width = 8, height = 12, compression = "lzw")
}

ridge_data <- plot_data_mean %>%
  inner_join(window_limits, by = config$key_variables$primary_group_var) %>%
  dplyr::filter(CAG >= Win_Min & CAG <= Win_Max) %>%
  group_by(!!sym(config$key_variables$primary_group_var), !!sym(config$key_variables$time_variable)) %>%
  mutate(Norm_Height = Mean_Rel_Height / max(Mean_Rel_Height, na.rm = TRUE)) %>%
  ungroup()

if(nrow(ridge_data) > 0) {
  ridge_data$Time_Factor <- factor(ridge_data[[config$key_variables$time_variable]], levels = sort(unique(ridge_data[[config$key_variables$time_variable]]), decreasing = TRUE))
  p_ridge <- ggplot(ridge_data, aes(x = CAG, y = Time_Factor, height = Norm_Height, fill = ..x..)) +
    geom_density_ridges_gradient(stat = "identity", scale = 1.5, rel_min_height = 0.01) +
    scale_fill_viridis_c(name = "CAG Length", option = "C") +
    geom_vline(data = baseline_visual_peaks, aes(xintercept = Center_CAG), linetype="dotted", alpha=0.5) +
    facet_wrap(vars(!!sym(pub_label_col)), scales = "free", ncol = 4) +
    labs(title = "Expansion Dynamics Over Time (Ridge Plots)", x = "CAG Length", y = "Timepoint") + theme_cowplot() + theme(legend.position = "none", axis.text.y = element_text(size = 8))
  ggsave(file.path(tiff_dir, "Publication_Ridge_Plots.tiff"), p_ridge, width = 16, height = 10, compression = "lzw")
}

# --- PLOT 4D: OVERVIEW GRID ---
logr::log_print("  -> Generating Faceted Trace Overview...", console=TRUE)
for (geno in unique_genos) {
  sub_data <- plot_data_full %>% 
    dplyr::filter(!!sym(config$key_variables$primary_group_var) == geno) %>%
    inner_join(window_limits %>% dplyr::filter(!!sym(config$key_variables$primary_group_var) == geno), by = config$key_variables$primary_group_var) %>%
    dplyr::filter(CAG >= Win_Min & CAG <= Win_Max)
  if(nrow(sub_data) > 0) {
    p_facet <- ggplot(sub_data, aes(x = CAG, y = Rel_Height)) + 
      geom_area(fill="grey80") + 
      geom_line(size=0.3) +
      facet_grid(rows = vars(!!sym(config$key_variables$secondary_group_var)), cols = vars(!!sym(config$key_variables$time_variable))) +
      labs(title = paste("Trace Overview:", geno), subtitle = "Rows = Clone/BioRep, Cols = Timepoint", y = "Intensity") + 
      theme_cowplot() +
      theme(strip.text = element_text(size=8), axis.text = element_blank(), axis.ticks = element_blank())
    safe_geno <- gsub("[^A-Za-z0-9_]", "_", geno)
    ggsave(file.path(tiff_dir, paste0("Overview_Grid_", safe_geno, ".tiff")), p_facet, width = 16, height = 12, compression = "lzw")
  }
}

# --- PLOT 4E: UNALIGNED START vs END (FULL RANGE, STACKED) ---
logr::log_print("  -> Generating Stacked Unaligned Start vs End Plot...", console=TRUE)

# 1. Filter Data: Start/End Only, Full Range
unaligned_data <- plot_data_mean %>%
  # Filter to configured expanded limits
  dplyr::filter(CAG >= limits_exp[1] & CAG <= limits_exp[2]) %>%
  # Keep only Start and End timepoints
  group_by(!!sym(config$key_variables$primary_group_var)) %>%
  dplyr::filter(!!sym(config$key_variables$time_variable) == min(!!sym(config$key_variables$time_variable)) | 
           !!sym(config$key_variables$time_variable) == max(!!sym(config$key_variables$time_variable))) %>%
  mutate(Timepoint_Label = ifelse(
    !!sym(config$key_variables$time_variable) == min(!!sym(config$key_variables$time_variable)), 
    "Start", "End"
  )) %>%
  mutate(Timepoint_Label = factor(Timepoint_Label, levels = c("Start", "End"))) %>%
  ungroup() %>% restore_all_factors()

# 3. Generate Plot
if(nrow(unaligned_data) > 0) {
  p_unaligned_stack <- ggplot(unaligned_data, aes(x = CAG, y = Mean_Rel_Height)) +
    # Fill areas for Start/End
    geom_area(aes(fill = Timepoint_Label), alpha = 0.5, position = "identity") +
    # Add line for definition
    geom_line(aes(color = Timepoint_Label), size = 0.3, alpha = 0.8) +
    
    # Facet: Stacked vertically (ncol=1)
    facet_wrap(vars(!!sym(pub_label_col)), ncol = 1, strip.position = "right") +
    
    # Aesthetics matching Aligned plots
    scale_fill_manual(values = c("Start" = "grey60", "End" = "red")) +
    scale_color_manual(values = c("Start" = "grey40", "End" = "darkred")) +
    
    scale_y_continuous(breaks = NULL, name = "Relative Height") +
    labs(title = "Unaligned Shift: Start vs End (Full Range)",
         subtitle = "Averaged Traces per Genotype",
         x = "CAG Length", 
         fill = "Timepoint", color = "Timepoint") +
    
    theme_cowplot() +
    theme(
      strip.text.y = element_text(angle = 0, face = "bold", size = 10), # Horizontal text on right
      panel.spacing.y = unit(0.2, "lines"), # Reduce gap between rows
      axis.line.y = element_blank() # Remove Y axis line for cleanliness
    )
  print(p_unaligned_stack)
  # Save
  ggsave(file.path(tiff_dir, "Publication_Start_vs_End_UNALIGNED_STACKED.tiff"), 
         p_unaligned_stack, width = 10, height = 14, compression = "lzw")
}

# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
# PLOT 4F: GENOTYPE BASELINE OVERVIEW (Replicates & Averages)
# ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
logr::log_print("  -> Generating Genotype Baseline Overview Plots...", console=TRUE)

# 1. Setup Directory
overview_dir <- file.path(trace_output_dir, "genotype_overview")
dir.create(overview_dir, showWarnings = FALSE, recursive = TRUE)

# 2. Define X-Axis Limits (Auto-detect to include WT + Expanded)
# This ensures we see the short WT peak and the long Expanded peak 
wt_limits_cfg <- config$trace_settings$wt_x_limits
exp_limits_cfg <- config$trace_settings$x_limits

x_min_ov <- if(!is.null(wt_limits_cfg)) min(wt_limits_cfg[1], exp_limits_cfg[1]) else exp_limits_cfg[1]
x_max_ov <- if(!is.null(wt_limits_cfg)) max(wt_limits_cfg[2], exp_limits_cfg[2]) else exp_limits_cfg[2]

# 3. Filter Data for Baseline Only
# We use 'plot_data_full' which already has the fancy labels joined from Part 4 Setup
baseline_full_data <- plot_data_full %>%
  dplyr::filter(!!sym(config$key_variables$time_variable) == min(!!sym(config$key_variables$time_variable), na.rm = TRUE))

if(nrow(baseline_full_data) > 0) {
  
  # --- A. AVERAGE PLOT (Mean Trace per Genotype) ---
  # Bin data to 0.1 CAG units to calculate a smooth average trace
  avg_trace_data <- baseline_full_data %>%
    mutate(CAG_Bin = round(CAG, 1)) %>%
    group_by(!!sym(plot_label_col), CAG_Bin) %>%
    summarise(Height = mean(Height, na.rm=TRUE), .groups="drop") %>%
    rename(CAG = CAG_Bin)
  
  p_avg_overview <- ggplot(avg_trace_data, aes(x = CAG, y = Height)) +
    geom_area(fill = "black", alpha = 0.8) + 
    geom_line(color = "black", linewidth = 0.3) +
    # Use free_y to account for WT peaks being much higher than mutant peaks
    facet_grid(rows = vars(!!sym(plot_label_col)), scales = "free_y", switch = "y") +
    scale_x_continuous(limits = c(x_min_ov, x_max_ov)) +
    labs(
      title = "Genotype Baseline Overview: Group Average",
      subtitle = "Composite average of all biological replicates at baseline timepoint.",
      x = "CAG Repeat Length",
      y = "Intensity (RFU)"
    ) +
    theme_cowplot() +
    theme(
      strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1, size = 10),
      axis.text.y = element_blank(), # Clean look
      axis.ticks.y = element_blank(),
      axis.line.y = element_blank(),
      panel.spacing = unit(0.2, "lines"),
      strip.background = element_blank()
    )
  
  ggsave(file.path(overview_dir, "Genotype_Overview_Average.tiff"), p_avg_overview, width = 8, height = 10, compression = "lzw")
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # --- B. SINGLE REPLICATE SELECTION (Auto-Fallback or External Override) ---
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  target_rep_id <- "1" 
  rep_var <- config$key_variables$optional_grouping_var %||% config$key_variables$repeated_measure_var
  if (is.null(rep_var) || rep_var == 'null') rep_var <- "rep"
  
  clone_var <- config$key_variables$secondary_group_var
  primary_var <- config$key_variables$primary_group_var
  time_var <- config$key_variables$time_variable
  
  # Check for external override file
  override_file_path <- file.path(analysis_base_dir, "representative_traces.csv")
  
  if (file.exists(override_file_path)) {
    logr::log_print("Using external 'representative_traces.csv' for plot selection...", console=TRUE)
    manual_selection <- read.csv(override_file_path, stringsAsFactors = FALSE)
    
    # Filter metadata to strictly match the rows provided in the CSV
    selected_meta <- peaks_annotated %>%
      dplyr::filter(!is_excluded) %>%
      semi_join(manual_selection, by = intersect(colnames(manual_selection), colnames(peaks_annotated)))
    
  } else {
    logr::log_print("Auto-selecting 1st clone and target replicate for overview plots...", console=TRUE)
    
    selected_meta <- peaks_annotated %>%
      dplyr::filter(!is_excluded) %>%
      group_by(!!sym(primary_var)) %>%
      # 1. Sort clones alphabetically/factor-order, then reps numerically
      arrange(!!sym(clone_var), as.numeric(as.character(!!sym(rep_var)))) %>%
      # 2. Isolate just the FIRST clone for this genotype
      dplyr::filter(!!sym(clone_var) == first(!!sym(clone_var))) %>%
      # 3. Prefer target_rep_id; if missing, fall back to the first available rep (min)
      dplyr::filter(if(any(!!sym(rep_var) == target_rep_id)) !!sym(rep_var) == target_rep_id else !!sym(rep_var) == first(!!sym(rep_var))) %>%
      slice(1) %>% # Ensure absolutely only 1 row per genotype (unless external file used)
      ungroup()
  }
  
  # Extract the bio_rep_ids we selected (so we can track them across timepoints)
  target_bioreps <- unique(selected_meta$bio_rep_id)
  
  # Create a Display Label (Genotype + Clone) just in case the external file specifies multiple clones per genotype
  plot_data_reps <- optimized_data %>%
    left_join(peaks_annotated, by = "fsa_filename") %>%
    dplyr::filter(bio_rep_id %in% target_bioreps) %>%
    mutate(Display_Label = paste(!!sym(pub_label_col), "|", !!sym(clone_var)))
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # --- C. BASELINE SPLIT PLOT (WT vs Expanded) ---
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # Filter to baseline timepoint only
  baseline_rep_data <- plot_data_reps %>% 
    dplyr::filter(!!sym(time_var) == min(!!sym(time_var), na.rm=TRUE))
  
  wt_lims <- c(10, 25)
  exp_lims <- c(120, 160)
  
  data_wt <- baseline_rep_data %>%
    dplyr::filter(CAG >= wt_lims[1] & CAG <= wt_lims[2]) %>%
    mutate(Panel_Region = "WT Allele")
  
  data_exp <- baseline_rep_data %>%
    dplyr::filter(CAG >= exp_lims[1] & CAG <= exp_lims[2]) %>%
    mutate(Panel_Region = "Expanded Allele")
  
  single_rep_split <- bind_rows(data_wt, data_exp) %>%
    mutate(Panel_Region = factor(Panel_Region, levels = c("WT Allele", "Expanded Allele")))
  
  if(nrow(single_rep_split) > 0) {
    # Base theme
    common_theme <- theme_cowplot() +
      theme(
        strip.text.y.left = element_text(angle = 0, face = "bold", hjust = 1, size = 10),
        strip.text.x = element_text(face = "bold", size = 11),
        strip.background = element_blank(),
        axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank()
      )
    
    # Plot WT
    p_wt_smooth <- ggplot(single_rep_split %>% dplyr::filter(Panel_Region == "WT Allele"), aes(x = CAG, y = Height)) +
      geom_area(aes(fill = !!sym(pub_label_col)), alpha = 0.4) + 
      geom_line(color = "black", linewidth = 0.3) +
      facet_grid(rows = vars(Display_Label), switch = "y") +
      scale_x_continuous(limits = wt_lims, breaks = seq(0, 100, by = 5)) +
      labs(title = NULL, x = NULL, y = "Intensity (RFU)") +
      scale_fill_manual(values = unname(create_custom_palette(single_rep_split, config, pub_label_col))) +
      common_theme + guides(fill = "none") + theme(plot.margin = margin(r = 5))
    
    # Plot Expanded
    p_exp_smooth <- ggplot(single_rep_split %>% dplyr::filter(Panel_Region == "Expanded Allele"), aes(x = CAG, y = Height)) +
      geom_area(aes(fill = !!sym(pub_label_col)), alpha = 0.4) + 
      geom_line(color = "black", linewidth = 0.3) +
      facet_grid(rows = vars(Display_Label), scales = "free_y") +
      scale_x_continuous(limits = exp_lims, breaks = seq(0, 1000, by = 20)) +
      labs(title = NULL, x = NULL, y = NULL) +
      scale_fill_manual(values = unname(create_custom_palette(single_rep_split, config, pub_label_col))) +
      common_theme + guides(fill = "none") + theme(strip.text.y = element_blank(), plot.margin = margin(l = 5))
    
    # Stitch
    final_grid_smooth <- plot_grid(p_wt_smooth, p_exp_smooth, align = "h", axis = "bt", nrow = 1, rel_widths = c(1, 2.5))
    final_plot_smooth <- ggdraw() +
      draw_plot(final_grid_smooth, x = 0, y = 0.05, width = 1, height = 0.9) +
      draw_label("Genotype Baseline Overview: Representative Trace", x = 0.5, y = 0.97, fontface = "bold", size = 14) +
      draw_label("CAG Repeat Length", x = 0.5, y = 0.02, fontface = "bold")
    
    ggsave(file.path(overview_dir, "Genotype_SingleRep_Smoothed.tiff"), final_plot_smooth, width = 8, height = 10, compression = "lzw")
  }
  
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  # --- D. NEW PLOT: SINGLE REP START VS END OVERLAY ---
  # ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~ #
  logr::log_print("  -> Generating Start vs End Representative Overlay...", console=TRUE)
  
  start_end_overlay_data <- plot_data_reps %>%
    dplyr::filter(CAG >= exp_lims[1] & CAG <= exp_lims[2]) %>%
    # Filter for Min and Max timepoints
    dplyr::filter(!!sym(time_var) == min(!!sym(time_var), na.rm=TRUE) | 
                    !!sym(time_var) == max(!!sym(time_var), na.rm=TRUE)) %>%
    mutate(Timepoint_Label = ifelse(
      !!sym(time_var) == min(!!sym(time_var), na.rm=TRUE), "Start", "End"
    )) %>%
    mutate(Timepoint_Label = factor(Timepoint_Label, levels = c("Start", "End"))) %>%
    
    # THE FIX: Isolate exactly one PCR per biological replicate per timepoint
    group_by(bio_rep_id, Timepoint_Label) %>%
    dplyr::filter(pcr == min(as.numeric(as.character(pcr)), na.rm = TRUE)) %>%
    ungroup() %>%
    
    # Normalize height per file
    group_by(fsa_filename) %>%
    mutate(Norm_Height = Height / max(Height, na.rm = TRUE)) %>%
    ungroup()
  
  if(nrow(start_end_overlay_data) > 0) {
    p_rep_overlay <- ggplot(start_end_overlay_data, aes(x = CAG, y = Norm_Height)) +
      geom_area(aes(fill = Timepoint_Label), alpha = 0.4, position = "identity") +
      geom_line(aes(color = Timepoint_Label), linewidth = 0.5) +
      # Facet by Display Label (Genotype + Clone) just in case the CSV loads multiple clones
      facet_wrap(vars(Display_Label), ncol = 1, scales = "free_y", strip.position = "right") +
      scale_fill_manual(values = c("Start" = "grey60", "End" = "red")) +
      scale_color_manual(values = c("Start" = "black", "End" = "darkred")) +
      labs(title = "Representative Single Trace: Start vs End Expansion",
           subtitle = "Showing Expanded Allele Region Only",
           x = "CAG Length", y = "Normalized Intensity", fill = "Timepoint", color = "Timepoint") +
      theme_cowplot() +
      theme(
        strip.text.y = element_text(angle = 0, face = "bold"),
        axis.text.y = element_blank(), axis.ticks.y = element_blank()
      )
    
    ggsave(file.path(overview_dir, "Representative_Start_vs_End_Overlay.tiff"), p_rep_overlay, width = 12, height = 8, compression = "lzw")
  }
  
} # <-- End of the if(nrow(baseline_full_data) > 0) block

#=============================================================================#
# PART 5: EXCEL DATA EXPORT & SUMMARY PLOT (With Raw Comparisons)
#=============================================================================#
logr::log_print("Generating Summary Excel Export & Peak Plot...", console=TRUE)

if(exists("aligned_data") && nrow(aligned_data) > 0) {
  
  # 1. Calculate Raw Modes (Average of GeneMapper Modes) for Comparison
  # -------------------------------------------------------------------------
  # This gets the actual numeric average of the replicates, not the peak of the trace
  raw_modes_wide <- peaks_annotated %>%
    dplyr::filter(!is_excluded) %>%
    group_by(!!sym(config$key_variables$primary_group_var)) %>%
    # Filter for Start and End timepoints only
    dplyr::filter(!!sym(config$key_variables$time_variable) == min(!!sym(config$key_variables$time_variable), na.rm=TRUE) | 
             !!sym(config$key_variables$time_variable) == max(!!sym(config$key_variables$time_variable), na.rm=TRUE)) %>%
    mutate(Timepoint_Label = ifelse(
      !!sym(config$key_variables$time_variable) == min(!!sym(config$key_variables$time_variable), na.rm=TRUE), 
      "Start", "End"
    )) %>%
    group_by(!!sym(config$key_variables$primary_group_var), Timepoint_Label) %>%
    summarise(
      Raw_Mean_Mode = mean(as.numeric(!!sym(col_mode)), na.rm = TRUE), 
      .groups = "drop"
    ) %>%
    tidyr::pivot_wider(names_from = Timepoint_Label, values_from = Raw_Mean_Mode, names_prefix = "Raw_Mode_")
  
  # 2. Generate Summary Export Table (Merged)
  # -------------------------------------------------------------------------
  summary_export_df <- aligned_data %>%
    group_by(!!sym(config$key_variables$primary_group_var), Timepoint_Label) %>%
    summarise(
      Plotted_Modal_Peak = CAG[which.max(Facet_Norm_Height)],
      .groups = "drop"
    ) %>%
    # Pivot Plotted Peaks
    tidyr::pivot_wider(names_from = Timepoint_Label, values_from = Plotted_Modal_Peak, names_prefix = "Plotted_Peak_") %>%
    
    # Join Window Limits & Baseline Info
    left_join(window_limits, by = config$key_variables$primary_group_var) %>%
    rename(Baseline_Peak_Used = Center_CAG) %>%
    
    # Join the Raw Modes calculated above
    left_join(raw_modes_wide, by = config$key_variables$primary_group_var) %>%
    
    mutate(Plot_Type = "Aligned_Start_vs_End") %>%
    
    # Select and Order Columns cleanly
    dplyr::select(
      !!sym(config$key_variables$primary_group_var),
      Plotted_Peak_Start,
      Plotted_Peak_End,
      Raw_Mode_Start,   # <--- New Column
      Raw_Mode_End,     # <--- New Column
      Win_Min, 
      Win_Max, 
      Baseline_Peak_Used
    )
  
  # 3. Export to Excel
  # -------------------------------------------------------------------------
  xlsx_path <- file.path(excel_dir, "Summary_Data_Aligned_Peaks.xlsx")
  tryCatch({
    writexl::write_xlsx(summary_export_df, xlsx_path)
    logr::log_print(paste("Saved Summary Excel to:", xlsx_path))
  }, error = function(e) {
    write.csv(summary_export_df, file.path(excel_dir, "Summary_Data_Aligned_Peaks.csv"), row.names = FALSE)
  })
  
  # 4. Generate Summary Peak Plot (Visualizing the Shift)
  # -------------------------------------------------------------------------
  # We need a long format version just for plotting the "Plotted" peaks
  summary_long_for_plot <- aligned_data %>%
    group_by(!!sym(config$key_variables$primary_group_var), Timepoint_Label) %>%
    summarise(Plotted_Modal_Peak = CAG[which.max(Facet_Norm_Height)], .groups = "drop") %>%
    left_join(window_limits, by = config$key_variables$primary_group_var) %>%
    rename(Baseline_Peak_Used = Center_CAG)
  
  p_peak_summary <- ggplot(summary_long_for_plot, aes(x = Plotted_Modal_Peak, y = Genotype_Pub)) +
    geom_errorbarh(aes(xmin = Win_Min, xmax = Win_Max), height = 0.2, color = "grey80", size = 2) +
    geom_point(aes(x = Baseline_Peak_Used), shape = 3, size = 3, color = "black") + 
    geom_line(aes(group = !!sym(config$key_variables$primary_group_var)), color = "grey50", size = 0.5) +
    geom_point(aes(color = Timepoint_Label), size = 4) +
    scale_color_manual(values = c("Start" = "grey60", "End" = "red")) +
    labs(title = "Summary of Aligned Peak Shifts",
         subtitle = "Grey Bar = Window Used | Cross (+) = Baseline Peak | Dots = Plotted Peaks",
         x = "CAG Length", y = "Genotype", color = "Timepoint") +
    theme_cowplot() +
    theme(panel.grid.major.x = element_line(color = "grey95"))
  
  ggsave(file.path(tiff_dir, "Summary_Aligned_Peaks_Plot.tiff"), p_peak_summary, width = 10, height = 8, compression = "lzw")
}

#=============================================================================#
# PART 6: STATISTICS PLOTS
#=============================================================================#
logr::log_print("Generating Statistical Summary Plots...", console=TRUE)

summary_data <- peaks_annotated %>%
  dplyr::filter(!is_excluded) %>%
  mutate(
    Genotype = as.factor(!!sym(plot_label_col)),
    Time = as.factor(!!sym(config$key_variables$time_variable))
  )

if (!is.null(config$factor_levels) && primary_var %in% names(config$factor_levels)) {
  summary_data$Genotype <- factor(summary_data$Genotype, levels = config$factor_levels[[primary_var]])
}

if(col_ii %in% names(summary_data)) {
  p_ii <- ggplot(summary_data, aes(x = Time, y = !!sym(col_ii), fill = Genotype)) +
    geom_boxplot(alpha = 0.7, outlier.shape = NA) + 
    geom_point(position = position_jitterdodge(jitter.width = 0.2), size = 2, alpha = 0.8) +
    labs(title = "Somatic Instability Quantification", subtitle = "Instability Index", x = "Timepoint", y = "Instability Index") +
    theme_cowplot() + scale_fill_brewer(palette = "Set1") + theme(panel.grid.major.y = element_line(color="grey90"))
  ggsave(file.path(tiff_dir, "Summary_Instability_Index.tiff"), p_ii, width = 10, height = 6, compression = "lzw")
}

if(col_mode %in% names(summary_data)) {
  baseline_modes <- summary_data %>%
    dplyr::filter(!!sym(config$key_variables$time_variable) == min(as.numeric(as.character(Time)), na.rm=TRUE)) %>%
    dplyr::select(bio_rep_id, Baseline_Mode = !!sym(col_mode))
  shift_data <- summary_data %>% left_join(baseline_modes, by = "bio_rep_id") %>% mutate(Mode_Shift = !!sym(col_mode) - Baseline_Mode)
  
  p_shift <- ggplot(shift_data, aes(x = Time, y = Mode_Shift, color = Genotype, group = Genotype)) +
    stat_summary(fun = mean, geom = "line", size = 1, position = position_dodge(width = 0.5)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5), alpha = 0.5) +
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey50") +
    labs(title = "Net Expansion Over Time", x = "Timepoint", y = "Change in CAG Length") +
    theme_cowplot() + scale_color_brewer(palette = "Set1") + theme(panel.grid.major.y = element_line(color="grey90"))
  ggsave(file.path(tiff_dir, "Summary_Mode_Shift.tiff"), p_shift, width = 10, height = 6, compression = "lzw")
}

if(nrow(optimized_data) > 0) {
  calc_smooth_mode <- function(cag_vec, height_vec) {
    if(length(cag_vec) < 10) return(NA)
    tryCatch({
      fit <- smooth.spline(x = cag_vec, y = height_vec, spar = 0.25)
      pred_x <- seq(min(cag_vec), max(cag_vec), length.out = 500)
      return(pred_x[which.max(predict(fit, pred_x)$y)])
    }, error = function(e) return(NA))
  }
  window_lookup <- peaks_annotated %>% dplyr::select(fsa_filename, target_start, target_end) %>% distinct(fsa_filename, .keep_all = TRUE)
  smoothed_modes <- optimized_data %>% left_join(window_lookup, by = "fsa_filename") %>% dplyr::filter(CAG >= target_start & CAG <= target_end) %>% group_by(fsa_filename) %>% summarise(Smoothed_Mode = calc_smooth_mode(CAG, Height), .groups="drop")
  smoothed_summary <- summary_data %>% left_join(smoothed_modes, by = "fsa_filename") %>% dplyr::filter(!is.na(Smoothed_Mode)) 
  baseline_smooth <- smoothed_summary %>% dplyr::filter(as.numeric(as.character(Time)) == min(as.numeric(as.character(Time)), na.rm=TRUE)) %>% group_by(bio_rep_id) %>% summarise(Baseline_Smooth = mean(Smoothed_Mode, na.rm=TRUE), .groups="drop")
  shift_data_smooth <- smoothed_summary %>% left_join(baseline_smooth, by = "bio_rep_id") %>% dplyr::filter(!is.na(Baseline_Smooth)) %>% mutate(Mode_Shift = Smoothed_Mode - Baseline_Smooth)
  
  p_smooth <- ggplot(shift_data_smooth, aes(x = Time, y = Mode_Shift, color = Genotype, group = Genotype)) +
    stat_summary(fun = mean, geom = "line", linewidth = 1, position = position_dodge(width = 0.5)) +
    stat_summary(fun.data = mean_se, geom = "errorbar", width = 0.2, position = position_dodge(width = 0.5)) +
    geom_point(position = position_dodge(width = 0.5), alpha = 0.3) +
    labs(title = "Net Expansion (Smoothed Mode)", x = "Timepoint", y = "Change in CAG Length") + theme_cowplot() + scale_color_brewer(palette = "Set1")
  ggsave(file.path(tiff_dir, "Summary_Mode_Shift_Smoothed.tiff"), p_smooth, width = 10, height = 6, compression = "lzw")
}


#=============================================================================#
# PART 7: SMOOTHED TRACE INSPECTION REPORTS (PDF)
#=============================================================================#
logr::log_print("Generating Smoothed Trace Inspection Reports...", console=TRUE)

plot_smooth_check <- function(df, x_min, x_max, win_min=NA, win_max=NA, 
                              spar_val=0.25, baseline_val=NA, ii_thresh_val=NA, orig_mode=NA) {
  
  calc_min <- if(!is.na(win_min)) win_min else x_min
  calc_max <- if(!is.na(win_max)) win_max else x_max
  df_calc <- df %>% dplyr::filter(CAG >= calc_min & CAG <= calc_max) 
  df_vis  <- df %>% dplyr::filter(CAG >= x_min & CAG <= x_max)        
  
  img <- image_graph(width = 400, height = 200, res = 96)
  
  if(nrow(df_calc) < 10) {
    print(ggplot() + annotate("text", x=0.5, y=0.5, label="Insufficient Data") + theme_void())
  } else {
    plot_obj <- tryCatch({
      fit <- smooth.spline(x = df_calc$CAG, y = df_calc$Height, spar = spar_val)
      pred_grid <- seq(min(df_calc$CAG), max(df_calc$CAG), length.out = 500)
      pred_data <- data.frame(CAG = pred_grid, Height = predict(fit, pred_grid)$y)
      list(pred = pred_data, peak = pred_data$CAG[which.max(pred_data$Height)], error = FALSE)
    }, error = function(e) list(error = TRUE, msg = e$message))
    
    label_txt <- paste0("Smooth: ", round(plot_obj$peak, 2))
    if(!is.na(orig_mode)) label_txt <- paste0(label_txt, " | Raw: ", round(orig_mode, 2))
    
    p <- ggplot() +
      geom_area(data = df_vis, aes(x = CAG, y = Height), fill = "grey80", alpha = 0.5) +
      {if(!is.na(win_min)) geom_vline(xintercept = win_min, linetype = "dotted", color = "darkgreen")} +
      {if(!is.na(win_max)) geom_vline(xintercept = win_max, linetype = "dotted", color = "darkgreen")} +
      {if(!is.na(baseline_val)) geom_vline(xintercept = baseline_val, linetype = "dashed", color = "grey50")} +
      {if(!is.na(orig_mode)) geom_vline(xintercept = orig_mode, linetype = "dotdash", color = "blue", alpha = 0.6)} +
      {if(!is.na(ii_thresh_val)) geom_hline(yintercept = ii_thresh_val, linetype = "dashed", color = "purple", alpha = 0.7)} +
      scale_x_continuous(limits = c(x_min, x_max)) +
      theme_minimal() +
      theme(axis.title = element_blank(), panel.grid = element_blank(), panel.border = element_rect(colour = "grey80", fill = NA))
    
    if(!plot_obj$error) {
      p <- p +
        geom_line(data = plot_obj$pred, aes(x = CAG, y = Height), color = "red", linewidth = 0.5) +
        geom_vline(xintercept = plot_obj$peak, linetype = "dashed", color = "blue") +
        annotate("text", x = -Inf, y = Inf, label = label_txt, hjust = -0.1, vjust = 1.5, size = 3, fontface = "bold", color = "blue")
    } else {
      p <- p + annotate("text", x=0.5, y=0.5, label="Fit Error", color="red")
    }
    print(p)
  }
  dev.off()
  return(img)
}

for (geno in unique_genos) {
  geno_data <- peaks_annotated %>% dplyr::filter(!!sym(config$key_variables$primary_group_var) == geno)
  unique_bioreps <- unique(geno_data$bio_rep_id)
  if(length(unique_bioreps) == 0) next
  logr::log_print(paste0("Generating Smoothed PDF: ", geno, "..."), console=TRUE)
  pb <- txtProgressBar(min = 0, max = length(unique_bioreps), style = 3, char = "=")
  pdf_pages_list <- list()
  
  for (i in seq_along(unique_bioreps)) {
    biorep <- unique_bioreps[i]
    rep_data <- geno_data %>% dplyr::filter(bio_rep_id == biorep)
    curr_clone <- unique(rep_data[[config$key_variables$secondary_group_var]])
    # 1. Extract the parent bio_id for the current replicate
    curr_bio_id <- unique(rep_data$bio_id)[1]
    
    # 2. Look up the baseline using bio_id instead of bio_rep_id
    base_val <- if(mode_name != "WT") { 
      val <- baseline_table %>% 
        dplyr::filter(bio_id == curr_bio_id) %>% 
        pull(baseline_cag)
      
      if(length(val) > 0) val[1] else NA 
    } else NA          
    rows_list <- list()
    
    # --- DYNAMIC PCR COUNT ---
    # Find the max PCRs for this specific replicate to adapt the layout
    max_pcr <- max(as.numeric(as.character(rep_data$pcr)), na.rm = TRUE)
    if (is.infinite(max_pcr) || is.na(max_pcr)) max_pcr <- 3
    
    for (day in sort(unique(rep_data[[config$key_variables$time_variable]]))) {
      day_data <- rep_data %>% dplyr::filter(!!sym(config$key_variables$time_variable) == day)
      pcr_imgs <- list()
      
      for (p in 1:max_pcr) {
        match_row <- day_data %>% dplyr::filter(pcr == p) %>% slice(1)
        if (nrow(match_row) > 0) {
          
          # --- THE FACTOR INDEXING FIX ---
          # Force the filename to character so R searches by actual file name
          clean_fname <- as.character(match_row$fsa_filename)
          trace <- trace_lookup[[clean_fname]]
          if(is.null(trace)) trace <- trace_lookup[[sub("\\.fsa$", "", clean_fname)]]
          
          w_min_sample <- if("target_start" %in% names(match_row)) as.numeric(match_row[["target_start"]]) else NA
          w_max_sample <- if("target_end" %in% names(match_row)) as.numeric(match_row[["target_end"]]) else NA
          curr_thresh_ii <- if("ii_threshold_abs" %in% names(match_row)) as.numeric(match_row[["ii_threshold_abs"]]) else NA
          curr_mode <- if(col_mode %in% names(match_row)) as.numeric(match_row[[col_mode]]) else NA
          
          if (!is.null(trace) && nrow(trace) > 0) {
            img <- plot_smooth_check(trace, limits_exp[1], limits_exp[2], w_min_sample, w_max_sample, spar_val = 0.25, baseline_val = base_val, ii_thresh_val = curr_thresh_ii, orig_mode = curr_mode)
            if(match_row$is_excluded) img <- image_colorize(img, opacity=25, color="pink") %>% image_annotate("EXCLUDED", size=18, color="red", gravity="center")
            pcr_imgs[[p]] <- img
          } else { 
            pcr_imgs[[p]] <- image_blank(width=400, height=200, color="white") %>% image_annotate("No Data", color="grey80", gravity="center", size=20) 
          }
        } else { 
          pcr_imgs[[p]] <- image_blank(width=400, height=200, color="white") %>% image_annotate("Missing", color="grey90", gravity="center", size=20) 
        }
      }
      
      # FIX: Dynamically capitalize the time variable name from config
      rows_list <- c(rows_list, list(image_append(c(image_blank(width=60, height=200, color="white") %>% 
                                                      image_annotate(paste0(stringr::str_to_title(config$key_variables$time_variable), " ", day), size=20, gravity="center", degrees=-90), image_append(image_join(pcr_imgs), stack = FALSE)), stack = FALSE)))
    }
    full_grid <- image_append(image_join(rows_list), stack = TRUE)
    
    # Look up fancy label for the current raw 'geno'
    curr_fancy_label <- peaks_annotated %>% 
      dplyr::filter(!!sym(config$key_variables$primary_group_var) == geno) %>% 
      pull(!!sym(plot_label_col)) %>% 
      unique() %>% .[1]
    
    if(is.na(curr_fancy_label)) curr_fancy_label <- geno
    
    title_img <- image_blank(width=image_info(full_grid)$width, height=50, color="white") %>% 
      image_annotate(paste0(curr_fancy_label, " | ", curr_clone, " | SMOOTHED CHECK"), size=24, gravity="center", weight=700)
    pdf_pages_list <- c(pdf_pages_list, list(image_append(c(title_img, full_grid), stack = TRUE)))
    gc(verbose = FALSE); setTxtProgressBar(pb, i); flush.console()
  }
  close(pb)
  if(length(pdf_pages_list) > 0) {
    safe_geno <- gsub("[^A-Za-z0-9_]", "_", geno)
    fname <- paste0("TraceReport_Smoothed_", safe_geno, ".pdf")
    image_write(image_join(pdf_pages_list), path = file.path(pdf_dir_smooth, fname), format = "pdf")
  }
}

logr::log_print("\n--- SCRIPT 05 COMPLETE ---", console = TRUE)
logr::log_close()