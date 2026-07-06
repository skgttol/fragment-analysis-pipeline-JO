#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#
# Advanced Fragment Analysis Workflow - Script 4 of 4
# (Extrinsic Multi-Omics Integration & Dimensionality Reduction)
#~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~#

#=============================================================================#
# PART 0: SETUP & ANCHOR RE-HYDRATION                                      ####
#=============================================================================#

packages <- c("tidyverse", "readxl", "here", "openxlsx", "broom",
              "ggpubr", "janitor", "yaml", "logr", "emmeans", 
              "ggeffects", "patchwork", "gridExtra", "grid",
              "RColorBrewer", "gtable", "ggnewscale", "GGally", "MASS",
              "ggrepel", "dendextend") 

installed_packages <- packages %in% rownames(utils::installed.packages())
if (any(installed_packages == FALSE)) utils::install.packages(packages[!installed_packages])
invisible(lapply(packages, library, character.only = TRUE))

if (!file.exists(here::here("config.yml"))) stop("CRITICAL ERROR: 'config.yml' not found.")
config <- yaml::read_yaml(here::here("config.yml"))

# --- Bulletproof Newest Directory Search ---
analysis_base_dir <- here::here(config$paths$output_dir_base)
all_dirs <- setdiff(list.dirs(path = analysis_base_dir, full.names = TRUE, recursive = FALSE), analysis_base_dir)
target_dirs <- all_dirs[grepl("_Analysis_v", basename(all_dirs))]
if(length(target_dirs) == 0) stop("CRITICAL ERROR: No analysis output folders found.")

latest_analysis_dir <- max(target_dirs)
rdata_path <- file.path(latest_analysis_dir, "processing_complete.RData")
if (!file.exists(rdata_path)) stop("CRITICAL ERROR: 'processing_complete.RData' not found.")

load(rdata_path)

# Ensure helper functions are attached
if(file.exists("functions.R")) source("functions.R")

# Setup Logging & Output Architecture
output_dir      <- file.path(latest_analysis_dir)
corr_plot_dir   <- file.path(output_dir, "04_external_correlations")
matrix_plot_dir <- file.path(corr_plot_dir, "correlation_matrices")
bivariate_dir   <- file.path(corr_plot_dir, "bivariate_scatterplots")
step_model_dir  <- file.path(corr_plot_dir, "stepwise_regressions")
pca_plot_dir    <- file.path(corr_plot_dir, "pca_projections")
clust_plot_dir  <- file.path(corr_plot_dir, "hierarchical_clustering")

invisible(lapply(c(matrix_plot_dir, bivariate_dir, step_model_dir, pca_plot_dir, clust_plot_dir), dir.create, recursive = TRUE, showWarnings = FALSE))

try(logr::log_close(), silent = TRUE)
lf <- logr::log_open(file.path(output_dir, "extrinsic_integration.log"), show_notes = FALSE, logdir = FALSE)
logr::log_print("\n\n--- SCRIPT 04: EXTRINSIC INTEGRATION INITIALIZED ---")

# Define primary mapping anchors
primary_var   <- config$key_variables$primary_group_var
secondary_var <- config$key_variables$secondary_group_var
time_var      <- config$key_variables$time_variable
grouping_var  <- if(is.null(secondary_var) || secondary_var == 'null') primary_var else secondary_var

plot_label_col <- if(exists("data_per_pcr") && "Genotype_Exp" %in% colnames(data_per_pcr)) "Genotype_Exp" else primary_var

# NEW (Bulletproof, Tidyverse-compliant)
if(exists("modeling_data")) {
  label_lookup <- modeling_data %>% 
    dplyr::select(!!sym(primary_var), dplyr::any_of(c("Genotype_Pub", "Genotype_Exp", plot_label_col))) %>% 
    dplyr::distinct()
  
  logr::log_print(paste("Label lookup matrix re-hydrated. Columns mapped:", paste(colnames(label_lookup), collapse = ", ")))
} else { 
  label_lookup <- NULL 
}

cag_col <- if("mode" %in% names(data_per_pcr)) "mode" else response_vars[1]

baseline_map_agg <- data_per_pcr %>%
  group_by(!!sym(primary_var), !!sym(secondary_var)) %>%
  dplyr::filter(!!sym(time_var) == min(!!sym(time_var), na.rm = TRUE)) %>%
  summarise(baseline_cag = mean(!!sym(cag_col), na.rm = TRUE), .groups = "drop")

best_models_path <- file.path(output_dir, "03_instability_models/Summary_Best_Models.csv")
if(file.exists(best_models_path)) cag_threshold_table <- read.csv(best_models_path)

#=============================================================================#
# PART 1: EXTRINSIC METADATA INGESTION & MASTER TRUTH TABLE                ####
#=============================================================================#
logr::log_print("\n--- Starting PART 1: Ingesting Extrinsic Metadata ---")

join_cols <- c(primary_var, secondary_var)
valid_join <- join_cols[join_cols %in% colnames(baseline_map_agg) & join_cols %in% colnames(all_group_slopes)]

master_data_by_group <- baseline_map_agg %>% left_join(all_group_slopes, by = valid_join)
if(exists("apply_renaming_and_factors")) master_data_by_group <- apply_renaming_and_factors(master_data_by_group, config)

meta_cfg <- config$external_metadata
imported_metadata_cols <- c() # The Memory Tracker

# 1A. Build Master 'By-Group' Data & Ingest Group Metadata (List-Safe) ----
if (!is.null(meta_cfg$group_level_path) && meta_cfg$group_level_path != 'null') {
  
  # 1. Use your custom get_file_paths function to unlist the YAML configuration safely
  group_meta_paths <- tryCatch({
    get_file_paths(meta_cfg$group_level_path, "\\.(xlsx|xls)$")
  }, error = function(e) { NULL })
  
  if (!is.null(group_meta_paths) && length(group_meta_paths) > 0) {
    logr::log_print(paste("...Found", length(group_meta_paths), "group-level metadata files to combine."))
    
    # 2. Loop through all files, fix scientific notation, and bind rows
    external_metadata_group <- purrr::map_df(group_meta_paths, function(current_path) {
      logr::log_print(paste("......Reading metadata file:", basename(current_path)))
      
      # Read as text to preserve formatting digits
      raw_file <- readxl::read_excel(current_path, col_types = "text") %>%
        readr::type_convert(guess_integer = FALSE) # Prevents conversion of long IDs to scientific double types
      
      # --- YOUR FIX: Reverse Excel's Scientific Notation Corruption ---
      if (grouping_var %in% colnames(raw_file)) {
        raw_file <- raw_file %>%
          dplyr::mutate(
            !!sym(grouping_var) := as.character(!!sym(grouping_var)),
            # Strip artificial decimal zeros (1.00E+10 -> 1E+10)
            !!sym(grouping_var) := stringr::str_replace_all(!!sym(grouping_var), "\\.0+E", "E"),
            # Strip optional '+' and leading zeros in exponent (E+06 -> E6)
            !!sym(grouping_var) := stringr::str_replace_all(!!sym(grouping_var), "E\\+?0*([0-9]+)", "E\\1")
          )
      }
      return(raw_file)
    })
    external_metadata_group <- apply_renaming_and_factors(external_metadata_group, config)
    
    # 3. Proceed with the Smart Merge on the combined dataset
    core_keys <- c(primary_var, grouping_var)
    batch_var <- config$key_variables$optional_crossed_effect %||% "differentiation" 
    
    if (all(core_keys %in% colnames(external_metadata_group))) {
      
      can_exact_match <- batch_var %in% colnames(master_data_by_group) && batch_var %in% colnames(external_metadata_group)
      
      if (can_exact_match) {
        logr::log_print(paste("...High-Resolution Match: Joining on", primary_var, ",", grouping_var, ", and", batch_var))
        master_data_by_group <- master_data_by_group %>% 
          dplyr::left_join(external_metadata_group, by = c(core_keys, batch_var))
        
      } else {
        logr::log_print(paste("...Low-Resolution Match: Averaging combined metadata by", primary_var, "and", grouping_var))        
        external_metadata_avg <- external_metadata_group %>%
          dplyr::group_by(dplyr::across(dplyr::all_of(core_keys))) %>%
          dplyr::summarise(
            dplyr::across(dplyr::where(is.numeric), ~mean(.x, na.rm = TRUE)),
            n_metadata_batches = dplyr::n(),
            .groups = "drop"
          )
        
        master_data_by_group <- master_data_by_group %>% 
          dplyr::left_join(external_metadata_avg, by = core_keys)
      }
      
      # --- THE CLEANUP FIX: Targeted removal ---
      # Explicitly remove the columns you don't want using the 'any_of' helper
      # This handles both the physical_wells column and the n_metadata_batches 
      # column, even if they were created dynamically.
      master_data_by_group <- master_data_by_group %>% 
        dplyr::select(-dplyr::any_of(c("physical_wells_sampled", "n_metadata_batches")))
      
      logr::log_print("Successfully attached external variables and cleaned unwanted metadata columns.")
    }
  } else {
    logr::log_print("WARNING: Config paths for group_level_path could not be resolved or files do not exist.")
  }
}        


# Keep strictly the imported columns that are numeric for the Univariate sweep
numeric_imported_cols <- master_data_by_group %>% 
  dplyr::select(dplyr::any_of(imported_metadata_cols)) %>% 
  dplyr::select(dplyr::where(is.numeric), -dplyr::any_of(c("n_metadata_batches"))) %>% 
  colnames()

#=============================================================================#
# PART 2: UNIVARIATE METRIC SWEEPS (1-ON-1)                                ####
#=============================================================================#
logr::log_print("\n--- Starting PART 3H: Automated 1-on-1 Extrinsic Metric Sweep ---")

# --- STATISTICAL N-CHECK ---
min_n_for_correlation <- 3 
run_by_group_stats <- nrow(master_data_by_group) >= min_n_for_correlation
run_by_genotype_stats <- nrow(master_data_by_genotype) >= min_n_for_correlation

if(!run_by_group_stats) stop("CRITICAL ABORT: Insufficient data points (N < 3) to execute correlation matrix.")


# 1. Automatically grab the primary outcome rate (e.g. slope_mode_change)
target_rate_col <- grep("^slope_", colnames(master_data_by_group), value = TRUE)[1]

# 2. Identify all numeric columns that were imported from outside
# (We explicitly subtract the slope outcome columns and base CAG variables)
extrinsic_metrics <- master_data_by_group %>%
  dplyr::select(where(is.numeric), -any_of(known_id_cols), -starts_with("slope_"), -starts_with("baseline_")) %>%
  colnames()

if(length(extrinsic_metrics) == 0) {
  logr::log_print("NOTE: No standalone numeric extrinsic columns found in metadata to run 1-on-1 sweep.")
} else {
  logr::log_print(paste("Sweeping Target [", target_rate_col, "] against", length(extrinsic_metrics), "extrinsic metrics..."))
  
  bivariate_summary_list <- list()
  y_rate_label <- pretty_labels[[target_rate_col]] %||% target_rate_col
  
  for(metric in extrinsic_metrics) {
    x_metric_label <- pretty_labels[[metric]] %||% metric
    
    # Isolate valid complete pairs for this specific metric
    clean_pair_df <- master_data_by_group %>%
      dplyr::select(!!sym(primary_var), x_val = !!sym(metric), y_val = !!sym(target_rate_col)) %>%
      tidyr::drop_na(x_val, y_val)
    
    if(!is.null(label_lookup)) {
      clean_pair_df <- clean_pair_df %>% left_join(label_lookup %>% distinct(), by = primary_var)
      plot_group_var <- if(plot_label_col %in% colnames(clean_pair_df)) plot_label_col else primary_var
    } else { plot_group_var <- primary_var }
    
    if(nrow(clean_pair_df) >= 3) {
      stats_test <- tryCatch(cor.test(clean_pair_df$x_val, clean_pair_df$y_val, method = "pearson"), error = function(e) NULL)
      
      r_est <- if(!is.null(stats_test)) stats_test$estimate else NA
      p_est <- if(!is.null(stats_test)) stats_test$p.value else NA
      
      # Log entry for the master ranking table
      bivariate_summary_list[[metric]] <- data.frame(
        Extrinsic_Metric = metric,
        Target_Outcome   = target_rate_col,
        Valid_Samples_N  = nrow(clean_pair_df),
        Pearson_R        = r_est,
        R_Squared        = r_est^2,
        P_Value          = p_est,
        Significance     = case_when(
          is.na(p_est)  ~ "Failed",
          p_est < 0.001 ~ "***",
          p_est < 0.01  ~ "**",
          p_est < 0.05  ~ "*",
          TRUE          ~ "ns"
        )
      )
      
      # Generate Individual Bivariate Plot
      p_ind <- ggplot(clean_pair_df, aes(x = x_val, y = y_val)) +
        geom_point(aes(color = !!sym(plot_group_var)), size = 3.5, alpha = 0.75) +
        geom_smooth(method = "lm", color = "black", linetype = "dashed", alpha = 0.15) +
        labs(
          title = paste("Bivariate Fit:", x_metric_label),
          subtitle = sprintf("Target: %s\nPearson R = %.3f (R² = %.3f) | p-value = %.4g", y_rate_label, r_est, r_est^2, p_est),
          x = x_metric_label,
          y = y_rate_label,
          color = "Genotype"
        ) +
        theme_publication()
      
      if(exists("apply_smart_palette")) p_ind <- apply_smart_palette(p_ind, plot_group_var)
      print(p_ind)
      
      safe_filename <- compress_filename(metric)
      ggsave(file.path(bivariate_dir, paste0(safe_filename, ".tiff")), p_ind, device = "tiff", width = 8, height = 6, dpi = 300, compression = "lzw")
      
    } else { logr::log_print(paste("Skipped bivariate sweep for", metric, "- insufficient pairs.")) }
  }
  
  # Save Master Sweep CSV
  if(length(bivariate_summary_list) > 0) {
    master_bivariate_table <- dplyr::bind_rows(bivariate_summary_list) %>% arrange(P_Value)
    write.csv(master_bivariate_table, file.path(bivariate_dir, "Compiled_Bivariate_Metric_Sweep.csv"), row.names = FALSE)
    logr::log_print("Successfully exported 'Compiled_Bivariate_Metric_Sweep.csv'.")
  }
}

#=============================================================================#
# PART 3: MULTI-OMICS CORRELATION MATRICES                                 ####
#=============================================================================#

logr::log_print("\n--- Starting PART 3: Global Correlation Matrices ---")

pretty_labels <- config$variable_pretty_labels %||% list()
if (!primary_var %in% names(pretty_labels)) pretty_labels[[primary_var]] <- stringr::str_to_title(primary_var)

# Bridge IDs safely
if(!is.null(label_lookup)) {
  vars_for_ggpairs <- master_data_by_group %>% 
    dplyr::left_join(label_lookup %>% dplyr::distinct(), by = primary_var) %>%
    dplyr::mutate(!!sym(primary_var) := !!sym(plot_label_col)) %>%
    dplyr::select(-!!sym(plot_label_col))
} else {
  vars_for_ggpairs <- master_data_by_group
}

vars_for_ggpairs <- vars_for_ggpairs %>% 
  dplyr::select(where(is.numeric), !!sym(primary_var)) %>%
  janitor::remove_empty("cols")

numeric_vars_to_correlate <- vars_for_ggpairs %>% dplyr::select(where(is.numeric))
numeric_col_names <- colnames(numeric_vars_to_correlate)

numeric_labels <- sapply(numeric_col_names, function(x) {
  if(x %in% names(pretty_labels)) stringr::str_wrap(pretty_labels[[x]], 15) else x
})

# 3A. Pure Tidyverse Heatmap (No reshape2)
if (ncol(numeric_vars_to_correlate) > 1) {
  cor_matrix <- cor(numeric_vars_to_correlate, use = "pairwise.complete.obs")
  cor_matrix[upper.tri(cor_matrix)] <- NA 
  
  melted_cor_matrix <- cor_matrix %>%
    as.data.frame() %>%
    tibble::rownames_to_column(var = "Var1") %>%
    tidyr::pivot_longer(cols = -Var1, names_to = "Var2", values_to = "value", values_drop_na = TRUE) %>%
    mutate(
      Var1 = factor(Var1, levels = numeric_col_names, labels = numeric_labels),
      Var2 = factor(Var2, levels = numeric_col_names, labels = numeric_labels)
    )
  
  p_corr_heatmap <- ggplot(melted_cor_matrix, aes(x = Var1, y = Var2, fill = value)) +
    geom_tile(color = "white") +
    scale_fill_gradient2(low = "#0072B2", high = "#D55E00", mid = "white", midpoint = 0, limit = c(-1, 1), name = "Pearson R") +
    geom_text(aes(label = round(value, 2)), color = "black", size = 3) +
    labs(title = "Global Bivariate Correlation Matrix", x = "", y = "") +
    theme_minimal(base_size = 12) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1), panel.grid = element_blank()) +
    coord_fixed()
  
  print(p_corr_heatmap)
  ggsave(file.path(matrix_plot_dir, "corr_heatmap_global.tiff"), p_corr_heatmap, device = "tiff", width = 10, height = 9, dpi = 300, compression = "lzw")
}


#=============================================================================#
# PART 4: MULTIVARIATE STEPWISE REGRESSION                                 ####
#=============================================================================#

logr::log_print("\n--- Starting PART 4, 5 & 6: Stepwise AIC Regression ---")

outcome_vars     <- grep("^slope_", colnames(master_data_by_group), value = TRUE)
all_numeric_cols <- master_data_by_group %>% dplyr::select(where(is.numeric), -any_of(known_id_cols), -starts_with("slope_")) %>% colnames()

if (length(all_numeric_cols) >= 2) {
  wb_models <- openxlsx::createWorkbook()
  
  for (outcome_var in outcome_vars) {
    model_df <- master_data_by_group %>% dplyr::select(all_of(outcome_var), all_of(all_numeric_cols)) %>% na.omit()
    if (nrow(model_df) < (length(all_numeric_cols) + 2)) next
    
    full_model <- stats::lm(as.formula(paste(outcome_var, "~", paste(all_numeric_cols, collapse = " + "))), data = model_df)
    null_model <- stats::lm(as.formula(paste(outcome_var, "~ 1")), data = model_df)
    
    step_model <- tryCatch(MASS::stepAIC(full_model, scope=list(lower=null_model, upper=full_model), direction="both", trace=0), error = function(e) NULL)
    
    if (!is.null(step_model)) {
      step_tidy <- broom::tidy(step_model, conf.int = TRUE) %>%
        dplyr::mutate(term_pretty = sapply(term, function(x) pretty_labels[[x]] %||% x))
      
      openxlsx::addWorksheet(wb_models, substr(outcome_var, 1, 31))
      openxlsx::writeData(wb_models, substr(outcome_var, 1, 31), step_tidy)
      
      # Plot Drivers
      plot_df <- step_tidy %>% dplyr::filter(term != "(Intercept)") %>% mutate(label = reorder(term_pretty, estimate))
      if(nrow(plot_df) > 0) {
        p_coef <- ggplot(plot_df, aes(x = estimate, y = label, color = p.value < 0.05)) +
          geom_vline(xintercept = 0, linetype="dashed", color="grey50") +
          geom_errorbar(aes(xmin = conf.low, xmax = conf.high), width=0.2, linewidth=1) +
          geom_point(size = 4) +
          scale_color_manual(values = c("TRUE" = "#D55E00", "FALSE" = "grey60"), name="p < 0.05") +
          labs(title = paste("Stepwise Predictors of", outcome_var), x = "Effect on Rate", y = "") +
          theme_publication()
        
        print(p_coef)
        ggsave(file.path(step_model_dir, paste0("stepwise_coefs_", outcome_var, ".tiff")), p_coef, device="tiff", width=8, height=6, compression="lzw")
      }
    }
  }
  openxlsx::saveWorkbook(wb_models, file.path(step_model_dir, "Stepwise_Model_Results.xlsx"), overwrite = TRUE)
}

#=============================================================================#
# PART 5: DIMENSIONALITY REDUCTION (PCA)                                   ####
#=============================================================================#
logr::log_print("\n--- Starting PART 7: Principal Component Analysis ---")

# CRITICAL FIX: Attach a permanent numeric Row ID before passing to na.omit()
pca_source_data <- master_data_by_group %>% 
  tibble::rowid_to_column("safe_row_id")

pca_data_clean <- pca_source_data %>%
  dplyr::select(safe_row_id, where(is.numeric), -any_of(known_id_cols)) %>%
  na.omit() 

if(nrow(pca_data_clean) > 3 && ncol(pca_data_clean) > 2) {
  
  # Run PCA strictly on the math columns (excluding the row ID)
  pca_fit <- prcomp(pca_data_clean %>% dplyr::select(-safe_row_id), scale. = TRUE, center = TRUE)
  
  # Join metadata back safely using the exact tracked ID
  pca_coords_with_meta <- pca_fit$x %>%
    as.data.frame() %>%
    dplyr::mutate(safe_row_id = pca_data_clean$safe_row_id) %>%
    dplyr::inner_join(pca_source_data, by = "safe_row_id")
  
  pca_var <- broom::tidy(pca_fit, matrix = "pcs")
  
  p_biplot <- ggplot(pca_coords_with_meta, aes(x = PC1, y = PC2)) +
    geom_hline(yintercept = 0, linetype="dashed", color="grey60") +
    geom_vline(xintercept = 0, linetype="dashed", color="grey60") +
    geom_point(aes(color = !!sym(primary_var)), size = 3.5, alpha = 0.8) +
    labs(
      title = "Integrated PCA Biplot", 
      x = sprintf("PC1 (%.1f%%)", pca_var$percent[1]*100), 
      y = sprintf("PC2 (%.1f%%)", pca_var$percent[2]*100),
      color = "Genotype"
    ) +
    theme_publication()
  
  if(exists("apply_smart_palette")) p_biplot <- apply_smart_palette(p_biplot, primary_var)
  print(p_biplot)
  ggsave(file.path(pca_plot_dir, "pca_biplot.tiff"), p_biplot, device="tiff", width=10, height=8, compression="lzw")
}

#=============================================================================#
# PART 6: PHENOTYPIC CLUSTERING                                            ####
#=============================================================================#

logr::log_print("\n--- Starting PART 8: Hierarchical Clustering ---")

if(exists("pca_data_clean") && nrow(pca_data_clean) > 3) {
  
  # Re-use our safe tracked IDs for the distance matrix
  dist_matrix <- dist(scale(pca_data_clean %>% dplyr::select(-safe_row_id)), method = "euclidean")
  hclust_fit  <- hclust(dist_matrix, method = "ward.D2")
  
  # Map metadata precisely to dendrogram tip order
  tip_order_ids <- pca_data_clean$safe_row_id[hclust_fit$order]
  meta_ordered  <- tibble(safe_row_id = tip_order_ids) %>% left_join(pca_source_data, by = "safe_row_id")
  
  active_palette <- if(exists("primary_palette")) primary_palette else NULL
  tip_colors <- if(!is.null(active_palette)) unname(active_palette[as.character(meta_ordered[[primary_var]])]) else rep("black", nrow(meta_ordered))
  
  dend <- as.dendrogram(hclust_fit) %>%
    dendextend::set("labels", as.character(meta_ordered[[grouping_var]])) %>%
    dendextend::set("labels_col", tip_colors) %>% 
    dendextend::hang.dendrogram(hang_height = 0.5)
  
  # Export to viewer
  plot(dend, main="Hierarchical Ward Clustering", ylab="Ward's Distance")
  
  # Export to File safely
  tiff(file.path(clust_plot_dir, "hierarchical_ward_tree.tiff"), width = 12, height = 8, units = "in", res = 300, compression = "lzw")
  plot(dend, main="Hierarchical Ward Clustering", ylab="Ward's Distance")
  dev.off()
}

logr::log_print("\n=== SCRIPT 04: MULTI-OMICS INTEGRATION COMPLETE ===")
logr::log_close()
