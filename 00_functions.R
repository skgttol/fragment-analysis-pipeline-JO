# ==============================================================================#
# 00_functions.R
# Core Functions for Hybrid DEP + Limma Proteomics Pipeline
# ==============================================================================#

# --- 1. Environment Setup ---####
setup_environment <- function() {
  if (isTRUE(getOption("pipeline_env_ready"))) {
    message("   -> Environment already initialised, skipping.")
    return(invisible(NULL))
  }
  # --- 1. Define Required Packages ---
  # Standard R packages for data wrangling and visualization
  cran_pkgs <- c(
    "tidyverse", "data.table", "pheatmap", "ggrepel", "cowplot", 
    "yaml", "openxlsx", "conflicted", "ggpubr", "ComplexUpset", "UpSetR", "magick",
    "patchwork", "RColorBrewer", "rmarkdown", "knitr", "kableExtra", "here", "base64enc"
  )
  
  # Bioconductor packages for proteomics, normalization, and functional GSEA
  bioc_pkgs <- c(
    "DEP", "limma", "SummarizedExperiment", "Biostrings", "ComplexHeatmap", 
    "vsn", "clusterProfiler", "org.Hs.eg.db", "msigdbr", "enrichplot", "STRINGdb"
  )
  
  # --- 2. Install BiocManager if missing ---
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  
  # --- 3. Loop, Check, and Install ---
  for (pkg in c(cran_pkgs, bioc_pkgs)) {
    if (!require(pkg, character.only = TRUE, quietly = TRUE)) {
      message(paste("Installing missing package:", pkg))
      if (pkg %in% bioc_pkgs) {
        BiocManager::install(pkg, update = FALSE, ask = FALSE)
      } else {
        install.packages(pkg, dependencies = TRUE)
      }
      # Load immediately after installation
      library(pkg, character.only = TRUE)
    }
  }
  
  # --- 4. Enforce Namespace Rules ---
  # Bioconductor packages notoriously overwrite basic dplyr functions. 
  # This locks them in so your data wrangling never crashes.
  library(conflicted)
  conflicts_prefer(dplyr::filter)
  conflicts_prefer(dplyr::select)
  conflicts_prefer(dplyr::rename)
  conflicts_prefer(dplyr::mutate)
  conflicts_prefer(base::as.factor)    # Add this
  conflicts_prefer(base::as.numeric)   # Add this too — same issue likely lurks here
  conflicts_prefer(base::intersect)    # Already used explicitly but worth locking in
  conflicts_prefer(base::setdiff)
  
  options(pipeline_env_ready = TRUE)
  message("--- Environment Ready ---")
}

# --- 2. Helper: Save and Print Plots --- ####
save_and_print_plot <- function(plot_obj, file_path, width = 10, height = 8) {
  
  # Sanitise the basename to remove special characters from contrast names
  # e.g. "WT - (IgG + BeadOnly)/2" -> "WT_vs_IgGBeadOnly2"
  safe_base      <- make_safe_contrast_name(basename(file_path))
  safe_file_path <- file.path(dirname(file_path), safe_base)
  
  # Determine how to render the object
  # pheatmap and recordedplot objects need grid rendering
  # ggplot and standard R plots use print()
  is_pheatmap <- inherits(plot_obj, "pheatmap")
  is_gtable   <- inherits(plot_obj, "gtable")
  is_recorded <- inherits(plot_obj, "recordedplot")
  
  render_plot <- function() {
    if (is_pheatmap) {
      grid::grid.newpage()
      grid::grid.draw(plot_obj$gtable)
    } else if (is_gtable) {
      grid::grid.newpage()
      grid::grid.draw(plot_obj)
    } else if (is_recorded) {
      replayPlot(plot_obj)
    } else {
      # ggplot, base R plots, and anything else
      print(plot_obj)
    }
  }
  
  # Print to screen
  render_plot()
  
  # Save TIFF at 300 dpi for publication
  tiff(paste0(safe_file_path, ".tiff"),
       width = width, height = height, units = "in", res = 300, compression = "lzw", type = "cairo")
  render_plot()
  dev.off()
  }

# 3. Contaminatn FIltering ####

parse_crap_fasta <- function(fasta_path, ms_data_list = NULL) {
  lines   <- readr::read_lines(fasta_path)
  headers <- lines[startsWith(lines, ">")]
  
  message(paste0("   -> Parsing cRAP FASTA: ", length(headers), " protein entries"))
  
  # Split header by pipe to capture all components
  # >sp|cRAP001|P00330|ADH1_YEAST Alcohol dehydrogenase...
  # x[1] = >sp, x[2] = cRAP001, x[3] = P00330, x[4] = ADH1_YEAST
  clean_headers <- sub("^>", "", headers)
  
  # Extract components into a list of vectors
  parts <- strsplit(clean_headers, "|", fixed = TRUE)
  
  # Collect all IDs from positions 2, 3, and 4
  all_ids <- unique(unlist(lapply(parts, function(x) {
    # Component 2 (cRAP ID), 3 (Accession), 4 (Entry Name + Description)
    ids <- x[2:4] 
    # Clean up the 4th component to remove description after the space
    ids[3] <- sub(" .*$", "", ids[3])
    return(ids)
  })))
  
  all_ids <- all_ids[!is.na(all_ids) & all_ids != ""]
  
  # Add gene symbol prefixes (e.g. ADH1 from ADH1_YEAST)
  gene_prefixes <- unique(stringr::str_match(all_ids, "^([A-Z0-9]+)_[A-Z]+$")[, 2])
  all_ids <- unique(c(all_ids, gene_prefixes[!is.na(gene_prefixes)]))
  
  # Cross-reference with MS data
  if (!is.null(ms_data_list)) {
    accessions_found <- c()
    for (exp_name in names(ms_data_list)) {
      df <- ms_data_list[[exp_name]]
      header_col <- intersect(c("T: Protein", "T: Majority protein IDs"), colnames(df))[1]
      if (!is.na(header_col)) {
        headers_in_data <- df[[header_col]]
        for (id in all_ids) {
          # Use fixed matching to avoid regex issues with special chars
          matching_rows <- headers_in_data[!is.na(headers_in_data) & stringr::str_detect(headers_in_data, stringr::fixed(id))]
          if (length(matching_rows) > 0) {
            # Extract accessions (e.g., P00330)
            acc <- stringr::str_extract(matching_rows, "[A-Z][0-9][A-Z0-9]{3}[0-9]|[A-Z][0-9][A-Z0-9]{3}[0-9][A-Z0-9]{2}")
            accessions_found <- c(accessions_found, acc[!is.na(acc)])
          }
        }
      }
    }
    all_ids <- unique(c(all_ids, accessions_found))
  }
  
  message(paste0("   -> cRAP ID set: ", length(all_ids), " unique identifiers"))
  return(all_ids)
}

filter_contaminants <- function(df, crap_ids, crapome_list, bioid_list, is_bioid = FALSE) {
  
  # Identify exactly WHY proteins are being flagged
  removed_crap <- df %>% dplyr::filter(ID %in% crap_ids) %>% 
    dplyr::mutate(Reason = "cRAP Contaminant")
  
  removed_crapome <- df %>% dplyr::filter(!ID %in% crap_ids & name %in% crapome_list) %>% 
    dplyr::mutate(Reason = "CRAPome (High Frequency)")
  
  message(paste("Loaded CRAPome. Found", length(crapome_list), "contaminants appearing in >", config$analysis_params$crapome_threshold, "experiments. These will be removed."))
  
  
  removed_bioid <- tibble::tibble()
  if (is_bioid) {
    removed_bioid <- df %>% dplyr::filter(!ID %in% crap_ids & !name %in% crapome_list & name %in% bioid_list) %>% 
      dplyr::mutate(Reason = "BioID Contaminant")
  }
  
  # Perform the actual filtering
  df_filt <- df %>% dplyr::filter(!ID %in% crap_ids) %>% dplyr::filter(!name %in% crapome_list)
  if (is_bioid) df_filt <- df_filt %>% dplyr::filter(!name %in% bioid_list)
  
  # Return both the clean data and the detailed log
  return(list(
    filtered = df_filt,
    log = dplyr::bind_rows(removed_crap, removed_crapome, removed_bioid) %>% dplyr::select(ID, name, Reason)
  ))
}

# --- 4. Helper: Per-Sample Imputation Plot --- ####
plot_imputation_per_sample <- function(se_norm, se_imputed, title) {
  df_norm <- as.data.frame(SummarizedExperiment::assay(se_norm)) %>%
    tibble::rownames_to_column("Protein") %>%
    tidyr::pivot_longer(-Protein, names_to = "Sample", values_to = "Intensity") %>%
    dplyr::filter(!is.na(Intensity)) %>%
    dplyr::mutate(State = "1_Before_Imputation")
  
  df_imp <- as.data.frame(SummarizedExperiment::assay(se_imputed)) %>%
    tibble::rownames_to_column("Protein") %>%
    tidyr::pivot_longer(-Protein, names_to = "Sample", values_to = "Intensity") %>%
    dplyr::mutate(State = "2_After_Imputation")
  
  df_combined <- dplyr::bind_rows(df_norm, df_imp)
  
  p <- ggplot(df_combined, aes(x = Intensity, fill = State)) +
    geom_density(alpha = 0.5, color = NA) +
    facet_wrap(~ Sample) + 
    scale_fill_manual(values = c("1_Before_Imputation" = "steelblue", "2_After_Imputation" = "firebrick")) +
    theme_bw() +
    labs(title = title, x = "VSN Normalized Intensity", y = "Density") +
    theme(strip.text = element_text(size = 8), legend.position = "bottom")
  
  return(p)
}

# --- 5. Custom Mixed Imputation (Clone-Aware MNAR + KNN) ---
impute_mixed_mar_mnar <- function(se, min_missing_pct = 0.75) {
  # Explicitly scope SummarizedExperiment calls
  mat <- SummarizedExperiment::assay(se)
  
  if ("clone_id" %in% colnames(SummarizedExperiment::colData(se))) {
    grouping_var <- SummarizedExperiment::colData(se)$clone_id
    message("   -> Applying mixed imputation at the CLONE level")
  } else {
    grouping_var <- SummarizedExperiment::colData(se)$condition
    message("   -> Applying mixed imputation at the CONDITION level")
  }
  
  # --- FIX 7: Compute noise parameters per sample, not globally ---
  # Each sample has its own detection floor and noise spread.
  # Using a global quantile would over-impute high-abundance samples
  # and under-impute low-abundance ones.
  noise_floor_vec <- apply(mat, 2, quantile, 0.01, na.rm = TRUE)
  sd_noise_vec    <- apply(mat, 2, sd, na.rm = TRUE) * 0.1
  
  # Safety: replace any NA/zero SDs (e.g. from near-empty samples) with a small fallback
  sd_noise_vec[is.na(sd_noise_vec) | sd_noise_vec == 0] <- 0.01
  
  for (grp in unique(grouping_var)) {
    col_idx      <- which(grouping_var == grp)
    missing_frac <- rowSums(is.na(mat[, col_idx, drop = FALSE])) / length(col_idx)
    mnar_rows    <- which(missing_frac >= min_missing_pct)
    
    for (r in mnar_rows) {
      na_cols <- col_idx[is.na(mat[r, col_idx])]
      # Each missing value is drawn from that specific sample's noise distribution
      mat[r, na_cols] <- rnorm(
        n    = length(na_cols),
        mean = noise_floor_vec[na_cols],
        sd   = sd_noise_vec[na_cols]
      )
    }
  }
  
  # Explicitly scope the reassignment and the DEP::impute call
  SummarizedExperiment::assay(se) <- mat
  se_final <- DEP::impute(se, fun = "knn", rowmax = 0.90)
  
  return(se_final)
}

# --- 6. Downstream Wrappers (PCA, Limma & Volcano) --- ####
generate_pca_plot <- function(mat, design, title) {
  pca      <- prcomp(t(na.omit(mat)), scale. = TRUE)
  pct_var  <- round(100 * summary(pca)$importance[2, 1:2], 1)
  
  pca_df <- data.frame(pca$x, label = design$label, condition = design$condition)
  
  p <- ggplot2::ggplot(pca_df, ggplot2::aes(x = PC1, y = PC2, colour = condition, label = label)) +
    ggplot2::geom_point(size = 3) +
    ggrepel::geom_text_repel(size = 3) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = title,
      x     = paste0("PC1 (", pct_var[1], "% variance)"),
      y     = paste0("PC2 (", pct_var[2], "% variance)")
    )
  return(p)
}

generate_pca_full <- function(mat, design, title, n_loadings = 20, output_dir) {
  
  # Remove rows with any NA (PCA requires complete matrix)
  mat_complete <- na.omit(mat)
  pca          <- prcomp(t(mat_complete), scale. = TRUE)
  pct_var      <- round(100 * summary(pca)$importance[2, ], 1)
  
  # ── 1. Scree plot (variance explained per component) ─────────────────────────
  scree_df <- data.frame(
    PC      = paste0("PC", 1:min(10, length(pct_var))),
    Var_Pct = pct_var[1:min(10, length(pct_var))]
  ) %>% dplyr::mutate(PC = factor(PC, levels = PC))
  
  p_scree <- ggplot2::ggplot(scree_df, ggplot2::aes(x = PC, y = Var_Pct)) +
    ggplot2::geom_col(fill = "steelblue", alpha = 0.8, colour = "white") +
    ggplot2::geom_line(ggplot2::aes(group = 1), colour = "black", linewidth = 0.5) +
    ggplot2::geom_point(size = 2) +
    ggplot2::geom_text(ggplot2::aes(label = paste0(Var_Pct, "%")), 
                       vjust = -0.5, size = 3) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title    = paste(title, "- Scree Plot"),
      subtitle = "Variance explained by each principal component",
      x        = "Principal Component",
      y        = "% Variance Explained"
    )
  
  save_and_print_plot(p_scree, file.path(output_dir, "PCA_1_Scree_Plot"))
  
  # ── 2. Scores plot (PC1 vs PC2, the standard PCA plot) ───────────────────────
  scores_df <- data.frame(
    pca$x[, 1:min(4, ncol(pca$x))],
    label     = design$label,
    condition = design$condition
  )
  
  p_scores <- ggplot2::ggplot(scores_df, 
                              ggplot2::aes(x = PC1, y = PC2, colour = condition)) +
    ggplot2::geom_point(size = 3.5, alpha = 0.9) +
    ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf) +
    {
      # Only draw ellipses for conditions with enough points
      ellipse_groups <- scores_df %>%
        dplyr::group_by(condition) %>%
        dplyr::summarise(n = dplyr::n(), .groups = "drop") %>%
        dplyr::filter(n >= 4) %>%
        dplyr::pull(condition)
      
      if (length(ellipse_groups) > 0) {
        ggplot2::stat_ellipse(
          data = dplyr::filter(scores_df, condition %in% ellipse_groups),
          ggplot2::aes(group = condition),
          type = "norm", linetype = "dashed", alpha = 0.4
        )
      }
    }    
  ggplot2::theme_bw() +
    ggplot2::labs(
      title    = paste(title, "- PCA Scores"),
      x        = paste0("PC1 (", pct_var[1], "%)"),
      y        = paste0("PC2 (", pct_var[2], "%)"),
      colour   = "Condition"
    )
  
  save_and_print_plot(suppressWarnings(p_scores), file.path(output_dir, "PCA_2_Scores_PC1vsPC2"))
  
  # Also plot PC2 vs PC3 — sometimes the biological effect is not in PC1
  if (ncol(pca$x) >= 3) {
    p_scores_23 <- ggplot2::ggplot(scores_df,
                                   ggplot2::aes(x = PC2, y = PC3, colour = condition)) +
      ggplot2::geom_point(size = 3.5, alpha = 0.9) +
      ggrepel::geom_text_repel(aes(label = label), size = 3, max.overlaps = Inf) +
      ggplot2::stat_ellipse(ggplot2::aes(group = condition),
                            type = "norm", linetype = "dashed", alpha = 0.4) +
      ggplot2::theme_bw() +
      ggplot2::labs(
        title  = paste(title, "- PCA Scores"),
        x      = paste0("PC2 (", pct_var[2], "%)"),
        y      = paste0("PC3 (", pct_var[3], "%)"),
        colour = "Condition"
      )
    save_and_print_plot(suppressWarnings(p_scores_23), file.path(output_dir, "PCA_3_Scores_PC2vsPC3"))
  }
  
  # ── 3. Loadings plot (which proteins drive each PC) ───────────────────────────
  # Loadings are in pca$rotation — rows are proteins, columns are PCs
  loadings_df <- as.data.frame(pca$rotation[, 1:min(4, ncol(pca$rotation))]) %>%
    tibble::rownames_to_column("Protein")
  
  for (pc in paste0("PC", 1:min(4, ncol(pca$rotation)))) {
    
    top_loadings <- loadings_df %>%
      dplyr::select(Protein, Loading = !!rlang::sym(pc)) %>%
      dplyr::mutate(abs_loading = abs(Loading),
                    Direction   = ifelse(Loading > 0, "Positive", "Negative")) %>%
      dplyr::slice_max(order_by = abs_loading, n = n_loadings)
    
    p_load <- ggplot2::ggplot(
      top_loadings,
      ggplot2::aes(x = Loading, 
                   y = reorder(Protein, Loading), 
                   fill = Direction)
    ) +
      ggplot2::geom_col(alpha = 0.85) +
      ggplot2::scale_fill_manual(values = c("Positive" = "#C0392B", "Negative" = "#2980B9")) +
      ggplot2::geom_vline(xintercept = 0, linewidth = 0.4) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none") +
      ggplot2::labs(
        title    = paste(title, "- Top", n_loadings, "Loadings for", pc),
        subtitle = paste0(pc, " explains ", pct_var[which(paste0("PC", 1:length(pct_var)) == pc)], 
                          "% of variance. Proteins furthest from zero drive this component most."),
        x        = paste("Loading on", pc),
        y        = "Protein"
      )
    
    save_and_print_plot(p_load, 
                        file.path(output_dir, paste0("PCA_4_Loadings_", pc)),
                        width = 8, height = 7)
  }
  
  # ── 4. Correlation of PCs with experimental metadata ─────────────────────────
  # This directly answers "is PC1 capturing genotype or something else?"
  scores_meta <- data.frame(pca$x[, 1:min(4, ncol(pca$x))]) %>%
    tibble::rownames_to_column("label") %>%
    dplyr::left_join(design, by = "label")
  
  # For each PC, compute correlation with condition as a numeric factor
  meta_cor_df <- lapply(paste0("PC", 1:min(4, ncol(pca$x))), function(pc) {
    lapply(c("condition"), function(meta_var) {
      vals <- as.numeric(as.factor(scores_meta[[meta_var]]))
      if (length(unique(vals)) < 2) return(NULL)
      cor_val <- cor(scores_meta[[pc]], vals, method = "spearman")
      data.frame(PC = pc, Metadata = meta_var, Spearman_r = round(cor_val, 3))
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()
  
  message("   -> PCA metadata correlations:")
  message(paste(capture.output(print(meta_cor_df)), collapse = "\n"))
  readr::write_csv(meta_cor_df, file.path(output_dir, "PCA_5_Metadata_Correlations.csv"))  
  # Return the PCA object for any downstream use
  return(invisible(pca))
}

# --- Helper: MDS Plot (PCA Alternative for Missing Values) --- ####
generate_mds_plot <- function(mat, design, title, output_dir) {
  
  # Calculate pairwise Pearson correlation, ignoring NAs
  cor_mat <- cor(mat, method = "pearson", use = "pairwise.complete.obs")
  
  # Convert correlation to a distance matrix (1 - correlation)
  # High correlation (1) = Distance of 0. Low correlation = larger distance.
  dist_mat <- as.dist(1 - cor_mat)
  
  # Perform Multi-Dimensional Scaling
  mds_fit <- cmdscale(dist_mat, eig = TRUE, k = 2)
  
  # Calculate percentage of variance explained by the dimensions
  eig <- mds_fit$eig
  pct_var <- round(100 * eig / sum(eig[eig > 0]), 1)
  
  mds_df <- data.frame(
    Dim1      = mds_fit$points[, 1],
    Dim2      = mds_fit$points[, 2],
    label     = design$label,
    condition = design$condition
  )
  
  p_mds <- ggplot2::ggplot(mds_df, ggplot2::aes(x = Dim1, y = Dim2, colour = condition)) +
    ggplot2::geom_point(size = 3.5, alpha = 0.9) +
    ggrepel::geom_text_repel(ggplot2::aes(label = label), size = 3, max.overlaps = Inf) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title    = paste(title, "- MDS Plot (PCA Alternative)"),
      subtitle = "Distance based on pairwise Pearson correlation (handles missing values)",
      x        = paste0("Dimension 1 (", pct_var[1], "%)"),
      y        = paste0("Dimension 2 (", pct_var[2], "%)")
    )
  
  save_and_print_plot(p_mds, file.path(output_dir, "QC7_MDS_Plot"))
}

generate_ma_plot <- function(de_result, title, fdr_thr = 0.05, lfc_thr = 1.0) {
  plot_data <- de_result %>%
    dplyr::mutate(
      Significant = adj.P.Val < fdr_thr & abs(logFC) > lfc_thr,
      # AveExpr is produced by limma's topTable automatically
    )
  
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = AveExpr, y = logFC)) +
    ggplot2::geom_point(ggplot2::aes(colour = Significant), alpha = 0.5, size = 1.2) +
    ggplot2::geom_hline(yintercept = 0,           colour = "black",  linewidth = 0.5) +
    ggplot2::geom_hline(yintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", colour = "grey40") +
    ggplot2::geom_smooth(method = "loess", se = FALSE, colour = "red", linewidth = 0.8) +
    ggplot2::scale_colour_manual(values = c("TRUE" = "#C0392B", "FALSE" = "grey75")) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title    = paste("MA Plot:", title),
      subtitle = "A loess trend near y=0 confirms normalisation is unbiased",
      x        = "Average Expression (AveExpr)",
      y        = "Log\u2082 Fold Change",
      colour   = "Significant"
    )
  return(suppressWarnings(p))
  }

perform_limma_analysis <- function(data_matrix, sample_metadata, contrast_definitions, blocking_variable, grouping_column) {
  # Force valid R names for the design matrix to prevent Limma parsing errors
  group_names <- make.names(as.character(sample_metadata[[grouping_column]]))
  group <- factor(group_names)  
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  
  block_ids <- make.names(blocking_variable)
  
  # Suppress duplicateCorrelation warning if blocking isn't necessary
  if(length(unique(block_ids)) == length(block_ids)) {
    fit <- limma::lmFit(data_matrix, design)
  } else {
    corfit <- limma::duplicateCorrelation(data_matrix, design, block = block_ids)
    fit <- limma::lmFit(data_matrix, design, block = block_ids, correlation = corfit$consensus)
  }
  
  # Safely translate the YAML contrasts to match the make.names() format
  # This ensures "WT - (IgG + Bead)/2" resolves correctly
  valid_contrasts <- sapply(contrast_definitions, function(x) {
    parsed_x <- x
    for(lvl in unique(as.character(sample_metadata[[grouping_column]]))) {
      # Replace the raw name with the safe name, ensuring word boundaries
      parsed_x <- gsub(paste0("\\b", lvl, "\\b"), make.names(lvl), parsed_x)
    }
    return(parsed_x)
  })
  
  cm <- limma::makeContrasts(contrasts = valid_contrasts, levels = design)
  fit2 <- limma::contrasts.fit(fit, cm)
  fit2 <- limma::eBayes(fit2)
  
  res_list <- list()
  for (i in seq_along(valid_contrasts)) {
    orig_name <- names(contrast_definitions)[i]
    if (is.null(orig_name) || orig_name == "") orig_name <- contrast_definitions[i]
    
    res_list[[orig_name]] <- limma::topTable(fit2, coef = i, number = Inf, sort.by = "none") %>% 
      tibble::rownames_to_column("Protein")
  }
  return(res_list)
}

generate_universal_volcano_plot <- function(res_df, 
                                            contrast_str, 
                                            exp_name, 
                                            fdr_thr = 0.05, 
                                            lfc_thr = 1.0, 
                                            poi_list = c(), 
                                            n_label = 15, 
                                            blacklist_proteins = c(),
                                            use_nominal_p = FALSE) {
  
  message("      -> Rendering universal volcano plot...")
  
  # 1. Standardize column names across Limma ('name', 'adj.P.Val') and Raw ('Protein', 'P_Value')
  plot_data <- res_df
  if (!"Protein" %in% colnames(plot_data) && "name" %in% colnames(plot_data)) {
    plot_data$Protein <- plot_data$name
  }
  if (!"Status" %in% colnames(plot_data)) {
    plot_data$Status <- "Quantitative"
  }
  
  # Determine which P-value column to use for the Y-axis
  if (use_nominal_p || !("adj.P.Val" %in% colnames(plot_data))) {
    p_col <- if ("P_Value" %in% colnames(plot_data)) "P_Value" else "P.Value"
    y_label_text <- "-Log10(Unadjusted P-Value)"
  } else {
    p_col <- "adj.P.Val"
    y_label_text <- "-Log10(FDR Adjusted P-Value)"
  }
  
  # 2. Extract condition names from contrast string (e.g., "WT - FAN1KO" -> Target: WT, Control: FAN1KO)
  cond_parts   <- stringr::str_split(contrast_str, "\\s*-\\s*")[[1]]
  target_cond  <- if (length(cond_parts) >= 1) cond_parts[1] else "Target"
  control_cond <- if (length(cond_parts) >= 2) cond_parts[2] else "Control"
  
  # 3. Calculate finite Y-axis ceiling
  finite_pvals <- plot_data[[p_col]][plot_data[[p_col]] > 0 & !is.na(plot_data[[p_col]])]
  max_quant_y  <- if (length(finite_pvals) > 0) max(-log10(finite_pvals), na.rm = TRUE) else 4
  max_quant_y  <- ifelse(is.infinite(max_quant_y) | is.nan(max_quant_y), 4, max_quant_y)
  
  # Synthetic ceiling for binary dropouts (15% above the highest quantitative hit)
  binary_ceiling <- max_quant_y * 1.15
  
  # 4. Build Plotting Variables & Apply Jitter
  set.seed(42) # Lock seed for reproducible jittering
  plot_data <- plot_data %>% 
    dplyr::filter(Status != "Not Enough Data", !is.na(logFC)) %>%
    dplyr::mutate(
      Raw_MinusLog10P = -log10(.data[[p_col]]),
      
      # Apply Y-jitter ONLY if the hit is a binary dropout
      Plot_Y = dplyr::case_when(
        grepl("Exclusive to", Status) ~ binary_ceiling + runif(dplyr::n(), -0.3, 0.3),
        is.infinite(Raw_MinusLog10P)  ~ binary_ceiling,
        TRUE                          ~ Raw_MinusLog10P
      ),
      Plot_X = logFC,
      
      # Classify significance categories dynamically (Upgraded for Two-Tier Matching)
      Significance = dplyr::case_when(
        Protein %in% poi_list ~ "Protein of Interest (POI)",
        
        # Use grepl to catch both "(High Confidence)" and "(Exploratory Candidate)" suffixes
        grepl(paste("Exclusive to", target_cond), Status, fixed = TRUE) ~ paste("Binary On:", target_cond),
        grepl(paste("Exclusive to", control_cond), Status, fixed = TRUE) ~ paste("Binary Off:", control_cond),
        
        # FIX: Use !!sym() instead of .data[[]] for maximum backward compatibility
        Status == "Quantitative" & !!sym(p_col) < fdr_thr & logFC > lfc_thr  ~ paste("Enriched in", target_cond),
        Status == "Quantitative" & !!sym(p_col) < fdr_thr & logFC < -lfc_thr ~ paste("Enriched in", control_cond),
        TRUE ~ "Not Significant / Background"
      ),
      Is_POI = Protein %in% poi_list,
      
      # NEW: Map precise transparency based on confidence tier
      Plot_Alpha = dplyr::case_when(
        Is_POI ~ 1.0,                                                   # POIs are always solid
        grepl("(High Confidence)", Status, fixed = TRUE) ~ 0.9,         # Tier 1 Binary: Solid
        grepl("(Exploratory Candidate)", Status, fixed = TRUE) ~ 0.4,   # Tier 2 Binary: Faded
        Significance != "Not Significant / Background" ~ 0.8,           # Quant Hits: Mostly solid
        TRUE ~ 0.3                                                      # Background: Very faded
      )
    )
  
  # Apply statistical blacklist masking if requested
  if (length(blacklist_proteins) > 0) {
    plot_data <- plot_data %>% dplyr::filter(!Protein %in% blacklist_proteins)
  }
  
  # 5. Build Smart Labeling Subset
  label_data <- plot_data %>%
    dplyr::filter(
      Is_POI | 
        grepl("Binary", Significance) | 
        (Significance != "Not Significant / Background" & !!sym(p_col) < 0.01) # FIX
    ) %>%
    dplyr::group_by(Significance) %>%
    dplyr::slice_max(order_by = abs(Plot_X), n = n_label, with_ties = FALSE) %>%
    dplyr::ungroup() %>%
    dplyr::bind_rows(dplyr::filter(plot_data, Is_POI)) %>%
    dplyr::distinct(Protein, .keep_all = TRUE)
  
  # 6. Axis boundaries and palette
  x_max <- max(abs(plot_data$Plot_X), na.rm = TRUE) * 1.1
  x_max <- ifelse(x_max < 5, 5, x_max)
  y_max <- max(plot_data$Plot_Y, na.rm = TRUE) * 1.1
  
  custom_colors <- c(
    "Protein of Interest (POI)"    = "#984EA3", # Purple
    setNames("#E41A1C", paste("Binary On:", target_cond)),          # Bright Red
    setNames("#377EB8", paste("Binary Off:", control_cond)),        # Blue
    setNames("#FB8072", paste("Enriched in", target_cond)),         # Soft Red
    setNames("#80B1D3", paste("Enriched in", control_cond)),        # Soft Blue
    "Not Significant / Background" = "#CCCCCC"                      # Grey
  )
  
  # Check if binary hits actually exist in this dataset to decide whether to draw the top band
  has_binary <- any(grepl("Binary", plot_data$Significance))
  
  # 7. Construct ggplot
  p_volc <- ggplot2::ggplot(plot_data, ggplot2::aes(x = Plot_X, y = Plot_Y)) +
    ggplot2::geom_hline(yintercept = -log10(fdr_thr), linetype = "dashed", color = "grey40", alpha = 0.7) +
    ggplot2::geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", color = "grey40", alpha = 0.7) +
    
    # Optional binary band guideline
    (if (has_binary) ggplot2::geom_hline(yintercept = binary_ceiling - 0.5, linetype = "dotted", color = "grey60", alpha = 0.5) else NULL) +
    (if (has_binary) ggplot2::annotate("text", x = 0, y = binary_ceiling, label = "Binary On/Off Dropouts (Jittered)", fontface = "italic", color = "grey30", size = 3) else NULL) +
    
    # Points layers (Now mapping alpha dynamically!)
    ggplot2::geom_point(data = dplyr::filter(plot_data, Significance == "Not Significant / Background"),
                        color = "#CCCCCC", ggplot2::aes(alpha = Plot_Alpha), size = 1.5) +
    ggplot2::geom_point(data = dplyr::filter(plot_data, Significance != "Not Significant / Background" & !Is_POI),
                        ggplot2::aes(color = Significance, alpha = Plot_Alpha), size = 2.5) +
    ggplot2::geom_point(data = dplyr::filter(plot_data, Is_POI),
                        ggplot2::aes(color = Significance, alpha = Plot_Alpha), size = 3.5, shape = 17) +
    
    # Directional Side Banners
    ggplot2::annotate("label", x = x_max * 0.85, y = y_max, 
                      label = paste0("Enriched in ", target_cond, " \u2192"), 
                      fontface = "bold", color = "#B2182B", fill = "#FDDBC7", size = 3.5, label.size = 0.5) +
    ggplot2::annotate("label", x = -x_max * 0.85, y = y_max, 
                      label = paste0("\u2190 Enriched in ", control_cond), 
                      fontface = "bold", color = "#2166AC", fill = "#D1E5F0", size = 3.5, label.size = 0.5) +
    
    # Smart Labeling (Fades exploratory text too!)
    ggrepel::geom_text_repel(
      data = label_data,
      ggplot2::aes(label = Protein, alpha = Plot_Alpha, fontface = ifelse(Is_POI, "bold", "plain")),
      size          = ifelse(label_data$Is_POI, 3.8, 3.0),
      color         = ifelse(label_data$Is_POI, "#000000", "#333333"),
      box.padding   = 0.4, point.padding = 0.3, max.overlaps = 25,
      min.segment.length = 0.1, segment.color = "grey50", segment.size = 0.4,
      show.legend   = FALSE
    ) +
    
    # Styling & Scales
    ggplot2::scale_color_manual(values = custom_colors) +
    ggplot2::scale_alpha_identity() +   # <--- CRITICAL FIX: Applies the specific numeric alphas
    ggplot2::scale_x_continuous(limits = c(-x_max, x_max), expand = ggplot2::expansion(mult = 0.05)) +
    ggplot2::scale_y_continuous(limits = c(0, y_max * 1.05), expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      legend.position = "bottom", legend.title = ggplot2::element_blank(),
      legend.text = ggplot2::element_text(size = 10),
      plot.title = ggplot2::element_text(face = "bold", size = 14),
      plot.subtitle = ggplot2::element_text(color = "grey30", size = 10),
      panel.grid.minor = ggplot2::element_blank(),
      panel.border = ggplot2::element_rect(color = "black", linewidth = 1)
    ) +
    ggplot2::labs(
      title    = paste("Differential Abundance:", exp_name, "|", contrast_str),
      subtitle = paste("Thresholds: P/FDR <", fdr_thr, "| |Log2FC| >", lfc_thr, "| Triangles = POIs"),
      x        = paste0("Log2 Fold Change (", target_cond, " / ", control_cond, ")"),
      y        = ifelse(has_binary, paste(y_label_text, "[Binary Hits Jittered at Top]"), y_label_text)
    )
  
  return(p_volc)
}

generate_uncompressed_sidecar_plot <- function(res_df, 
                                               contrast_str, 
                                               exp_name,
                                               fdr_thr        = 0.05, 
                                               lfc_thr        = 1.0, 
                                               poi_list       = c(), 
                                               use_nominal_p  = FALSE,
                                               noise_floor    = 16.5,
                                               n_label_nonsig = 15,
                                               blacklist_proteins = c()) {
  
  # ── 1. Extract & Clean Quantitative Data for the Universal Volcano ──────────
  quant_df <- res_df %>% 
    dplyr::filter(Status == "Quantitative" & !is.na(P_Value) & !is.na(logFC))
  
  # ── 2. Build Center Panel using your NEW Universal Volcano Function ─────────
  p_center <- generate_universal_volcano_plot(
    res_df             = quant_df,
    contrast_str       = contrast_str,
    exp_name           = exp_name,
    fdr_thr            = fdr_thr,
    lfc_thr            = lfc_thr,
    poi_list           = poi_list,
    n_label            = n_label_nonsig,
    blacklist_proteins = blacklist_proteins,
    use_nominal_p      = use_nominal_p
  )
  
  # Strip the local title/subtitle so we can add a global master title
  p_center <- p_center + ggplot2::ggtitle(NULL, subtitle = NULL)
  
  # ── 3. Filter & Rank Binary Dropouts by Physical Abundance ──────────────────
  binary_wt <- res_df %>% 
    dplyr::filter(grepl("Exclusive to", Status) & logFC > 0) %>%
    dplyr::mutate(
      Is_True_Dropout = Mean_Target >= (noise_floor + 1.5),
      Is_POI          = Protein %in% poi_list,
      Hit_Tier        = ifelse(grepl("(High Confidence)", Status, fixed = TRUE), "Tier 1", "Tier 2"),
    ) %>%
    dplyr::filter(Is_True_Dropout) %>%
    dplyr::arrange(dplyr::desc(Mean_Target))
  
  binary_mut <- res_df %>% 
    dplyr::filter(grepl("Exclusive to", Status) & logFC < 0) %>%
    dplyr::mutate(
      Is_True_Dropout = Mean_Control >= (noise_floor + 1.5),
      Is_POI          = Protein %in% poi_list,
      Hit_Tier        = ifelse(grepl("(High Confidence)", Status, fixed = TRUE), "Tier 1", "Tier 2"),

    ) %>%
    dplyr::filter(Is_True_Dropout) %>%
    dplyr::arrange(dplyr::desc(Mean_Control))
  
  # ── 4. Build Right Sidecar (Binary On in Target) ────────────────────────────
  p_right <- ggplot2::ggplot(binary_wt, ggplot2::aes(x = 1, y = Mean_Target)) +
    ggplot2::geom_point(ggplot2::aes(shape = Is_POI, fill = Is_POI, alpha = Hit_Tier), colour = "#E41A1C", size = 3, stroke = 0.8) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = Protein, alpha = Hit_Tier, fontface = ifelse(Hit_Tier == "Tier 1" | Is_POI, "bold", "plain")), 
      size = 2.8, direction = "y", nudge_x = 0.2, hjust = 0, segment.size = 0.3, max.overlaps = Inf
    ) +
    ggplot2::scale_shape_manual(values = c("TRUE" = 23, "FALSE" = 21), guide = "none") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "gold", "FALSE" = "#E41A1C"), guide = "none") +
    ggplot2::scale_alpha_manual(values = c("Tier 1" = 1.0, "Tier 2" = 0.4), guide = "none") +  # Fades the exploratory hits
    ggplot2::scale_x_continuous(limits = c(0.8, 2.2), expand = c(0,0)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(x = "Binary On", y = "Target Abundance (VSN)", subtitle = "Exclusive to Target") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), plot.subtitle = ggplot2::element_text(size = 9, face = "bold", colour = "#E41A1C", hjust = 0.5))
  
  # ── 5. Build Left Sidecar (Binary Off in Target) ────────────────────────────
  p_left <- ggplot2::ggplot(binary_mut, ggplot2::aes(x = 1, y = Mean_Control)) +
    ggplot2::geom_point(ggplot2::aes(shape = Is_POI, fill = Is_POI, alpha = Hit_Tier), colour = "#377EB8", size = 3, stroke = 0.8) +
    ggrepel::geom_text_repel(
      ggplot2::aes(label = Protein, alpha = Hit_Tier, fontface = ifelse(Hit_Tier == "Tier 1" | Is_POI, "bold", "plain")), 
      size = 2.8, direction = "y", nudge_x = -0.2, hjust = 1, segment.size = 0.3, max.overlaps = Inf
    ) +
    ggplot2::scale_shape_manual(values = c("TRUE" = 23, "FALSE" = 21), guide = "none") +
    ggplot2::scale_fill_manual(values = c("TRUE" = "gold", "FALSE" = "#377EB8"), guide = "none") +
    ggplot2::scale_alpha_manual(values = c("Tier 1" = 1.0, "Tier 2" = 0.4), guide = "none") +  # Fades the exploratory hits
    ggplot2::scale_x_continuous(limits = c(-0.2, 1.2), expand = c(0,0)) +
    ggplot2::theme_bw(base_size = 11) +
    ggplot2::labs(x = "Binary Off", y = "Control Abundance (VSN)", subtitle = "Exclusive to Control") +
    ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank(), panel.grid.major.x = ggplot2::element_blank(), plot.subtitle = ggplot2::element_text(size = 9, face = "bold", colour = "#377EB8", hjust = 0.5))
  
  # ── 6. Stitch Together & Add Global Master Title ────────────────────────────
  # FIX: Construct the title string from the passed variables before plotting
  title_str    <- paste("Differential Abundance:", exp_name, "|", contrast_str)
  subtitle_str <- paste("Thresholds: P/FDR <", fdr_thr, "| |Log2FC| >", lfc_thr, "| Triangles = POIs")
  
  final_layout <- (p_left + p_center + p_right) + 
    patchwork::plot_layout(widths = c(1, 4.5, 1)) +
    patchwork::plot_annotation(
      title = title_str,
      subtitle = subtitle_str,
      theme = ggplot2::theme(plot.title = ggplot2::element_text(hjust = 0.5, face = "bold", size = 14),
                             plot.subtitle = ggplot2::element_text(hjust = 0.5, color = "grey30", size = 10))
    )
  
  return(final_layout)
}

# Keeps your advanced, highly-customized Volcano plotting logic
generate_volcano_plot <- function(plot_data, 
                                  title, 
                                  fdr_thr      = 0.05, 
                                  lfc_thr      = 1.0, 
                                  poi_list     = c(),
                                  show_nominal_p = FALSE,   # Set TRUE if n < 200 proteins
                                  n_label_nonsig = 15) {    # How many top non-POI hits to label
  
  # --- 1. Choose which p-value column to use for the y-axis ---
  if (show_nominal_p) {
    plot_data <- plot_data %>% dplyr::mutate(y_val = -log10(P.Value))
    y_label   <- "-log\u2081\u2080 (Nominal P-value)"
    y_thr     <- -log10(0.05)
    thr_label <- "p = 0.05"
  } else {
    plot_data <- plot_data %>% dplyr::mutate(y_val = -log10(adj.P.Val))
    y_label   <- "-log\u2081\u2080 (FDR-adjusted P-value)"
    y_thr     <- -log10(fdr_thr)
    thr_label <- paste0("FDR = ", fdr_thr)
  }
  
  # --- 2. Classify significance using FDR regardless of y-axis choice ---
  #        (You always want the colour to reflect the real threshold)
  plot_data <- plot_data %>%
    dplyr::mutate(
      Significance = dplyr::case_when(
        adj.P.Val < fdr_thr & logFC >  lfc_thr ~ "Enriched",
        adj.P.Val < fdr_thr & logFC < -lfc_thr ~ "Depleted",
        TRUE                                    ~ "NS"
      ),
      Is_POI = Protein %in% poi_list
    )
  
  # --- 3. Count significant hits for the subtitle ---
  n_total   <- nrow(plot_data)
  n_up      <- sum(plot_data$Significance == "Enriched")
  n_down    <- sum(plot_data$Significance == "Depleted")
  
  subtitle_text <- paste0(
    "n = ", n_total, " proteins  |  ",
    "Enriched: ", n_up, "  |  Depleted: ", n_down,
    "  |  LFC > \u00B1", lfc_thr, ", ", thr_label
  )
  
  # --- 4. Identify proteins to label ---
  # Top hits by significance that are NOT POIs (avoid double-labelling)
  top_hits <- plot_data %>%
    dplyr::filter(Significance != "NS" & !Is_POI) %>%
    dplyr::slice_max(order_by = y_val, n = n_label_nonsig)
  
  poi_data     <- dplyr::filter(plot_data, Is_POI)
  poi_sig_data <- dplyr::filter(plot_data, Is_POI & Significance != "NS")
  
  # --- 5. Colour palette ---
  sig_colors <- c(
    "Enriched" = "#C0392B",   # Clean red
    "Depleted" = "#2980B9",   # Clean blue
    "NS"       = "grey80"
  )
  
  # --- 6. Build the plot ---
  p <- ggplot2::ggplot(plot_data, ggplot2::aes(x = logFC, y = y_val)) +
    
    # Threshold lines — draw first so they sit behind points
    ggplot2::geom_hline(yintercept = y_thr,          linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    ggplot2::geom_vline(xintercept = c(-lfc_thr, lfc_thr), linetype = "dashed", colour = "grey40", linewidth = 0.4) +
    
    # NS points (draw first so significant sit on top)
    ggplot2::geom_point(
      data  = dplyr::filter(plot_data, Significance == "NS" & !Is_POI),
      colour = "grey80", alpha = 0.45, size = 1.2
    ) +
    
    # Significant points
    ggplot2::geom_point(
      data = dplyr::filter(plot_data, Significance != "NS" & !Is_POI),
      ggplot2::aes(colour = Significance), alpha = 0.75, size = 1.8
    ) +
    ggplot2::scale_colour_manual(
      values = sig_colors,
      labels = c("Enriched" = paste0("Enriched (n=", n_up, ")"),
                 "Depleted" = paste0("Depleted (n=", n_down, ")"),
                 "NS"       = "Not significant"),
      Protein   = NULL
    ) +
    
    # Top-hit labels (non-POI significant proteins)
    ggrepel::geom_text_repel(
      data          = top_hits,
      ggplot2::aes(label = Protein, colour = Significance),
      size          = 2.8,
      max.overlaps  = 20,
      segment.size  = 0.3,
      segment.alpha = 0.5,
      show.legend   = FALSE
    ) +
    
    # POI points — drawn on top of everything else as gold stars
    ggplot2::geom_point(
      data   = poi_data,
      shape  = 23,          # Diamond
      size   = 3.5,
      fill   = "gold",
      colour = "black",
      stroke = 0.6
    ) +
    
    # POI labels — only label POIs that are actually significant to avoid clutter
    # Change to `poi_data` if you want all POIs labelled regardless
    ggrepel::geom_text_repel(
      data          = poi_sig_data,
      ggplot2::aes(label = Protein),
      colour        = "black",
      fontface      = "bold",
      size          = 3.2,
      max.overlaps  = Inf,
      box.padding   = 0.5,
      segment.size  = 0.4,
      segment.colour = "black"
    ) +
    
    # Annotate the FDR/p threshold line directly on the plot
    ggplot2::annotate(
      "text",
      x     = Inf,
      y     = y_thr + 0.15,
      label = thr_label,
      hjust = 1.05,
      size  = 2.8,
      colour = "grey40",
      fontface = "italic"
    ) +
    
    # Axes and labels
    ggplot2::labs(
      title    = title,
      subtitle = subtitle_text,
      x        = "Log\u2082 Fold Change",
      y        = y_label
    ) +
    
    ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title       = ggplot2::element_text(hjust = 0.5, face = "bold", size = 13),
      plot.subtitle    = ggplot2::element_text(hjust = 0.5, size = 9, colour = "grey40"),
      legend.position  = "bottom",
      legend.text      = ggplot2::element_text(size = 9),
      panel.grid.minor = ggplot2::element_blank(),
      panel.grid.major = ggplot2::element_line(colour = "grey92"),
      axis.title       = ggplot2::element_text(size = 11)
    )
  
  return(suppressWarnings(p))
  }

# --- Helper: Generate Color-Coded Imputation Excel ---
export_imputation_map_excel <- function(se_norm, imputed_matrix, exp_design, output_file, 
                                        min_missing_pct = 0.75) {
  
  message("   -> Generating Color-Coded Imputation Map...")
  
  # 1. Extract the pre-imputation data (with NAs)
  mat_norm  <- SummarizedExperiment::assay(se_norm)
  is_na_mat <- is.na(mat_norm)
  
  # 2. Build a blank "Status" matrix
  status_mat <- matrix(
    "Observed",
    nrow     = nrow(mat_norm),
    ncol     = ncol(mat_norm),
    dimnames = dimnames(mat_norm)
  )
  
  # 3. Determine Grouping Variable (Clone vs Condition)
  grouping_var <- if ("clone_id" %in% colnames(exp_design)) exp_design$clone_id else exp_design$condition
  
  # 4. Re-run the MNAR logic using the same per-sample noise parameters as impute_mixed_mar_mnar()
  #    This guarantees the colour map exactly reflects what imputation actually did.
  #    Previously this used a global noise floor — now it mirrors the per-sample approach.
  noise_floor_vec <- apply(mat_norm, 2, quantile, 0.01, na.rm = TRUE)
  sd_noise_vec    <- apply(mat_norm, 2, sd, na.rm = TRUE) * 0.1
  sd_noise_vec[is.na(sd_noise_vec) | sd_noise_vec == 0] <- 0.01  # Safety fallback
  
  for (grp in unique(grouping_var)) {
    col_idx      <- which(grouping_var == grp)
    missing_frac <- rowSums(is_na_mat[, col_idx, drop = FALSE]) / length(col_idx)
    mnar_rows    <- which(missing_frac >= min_missing_pct)
    
    for (r in mnar_rows) {
      na_cols <- col_idx[is_na_mat[r, col_idx]]
      status_mat[r, na_cols] <- "MNAR_Imputed"
    }
  }
  
  # 5. Any NAs not captured by the MNAR logic were filled by KNN (MAR)
  status_mat[is_na_mat & status_mat == "Observed"] <- "KNN_Imputed"
  
  # 6. Initialize Excel Workbook
  wb         <- openxlsx::createWorkbook()
  sheet_name <- "Imputed_Data_Map"
  openxlsx::addWorksheet(wb, sheet_name)
  
  # Write the numeric imputed values
  out_df <- as.data.frame(imputed_matrix) %>% tibble::rownames_to_column("Protein")
  openxlsx::writeData(wb, sheet_name, out_df)
  
  # 7. Define cell styles
  style_mnar   <- openxlsx::createStyle(fgFill = "#ff9999")  # Light red  — MNAR (noise floor)
  style_knn    <- openxlsx::createStyle(fgFill = "#99ccff")  # Light blue — KNN  (MAR)
  style_legend <- openxlsx::createStyle(textDecoration = "bold")
  
  # 8. Apply styles column by column
  #    +1 to rows for the Excel header row; +1 to cols for the "Protein" column
  for (c in 1:ncol(status_mat)) {
    mnar_rows_xl <- which(status_mat[, c] == "MNAR_Imputed") + 1
    knn_rows_xl  <- which(status_mat[, c] == "KNN_Imputed")  + 1
    
    if (length(mnar_rows_xl) > 0) openxlsx::addStyle(wb, sheet_name, style_mnar, rows = mnar_rows_xl, cols = c + 1)
    if (length(knn_rows_xl)  > 0) openxlsx::addStyle(wb, sheet_name, style_knn,  rows = knn_rows_xl,  cols = c + 1)
  }
  
  # 9. Add a legend at the top right of the sheet
  legend_col <- ncol(out_df) + 2
  openxlsx::writeData(wb, sheet_name, "LEGEND:",              startCol = legend_col, startRow = 1)
  openxlsx::addStyle(wb, sheet_name, style_legend, rows = 1, cols = legend_col)
  
  openxlsx::writeData(wb, sheet_name, "MNAR (Noise Floor)",   startCol = legend_col, startRow = 2)
  openxlsx::addStyle(wb, sheet_name, style_mnar,   rows = 2,  cols = legend_col)
  
  openxlsx::writeData(wb, sheet_name, "KNN (MAR Borrowing)",  startCol = legend_col, startRow = 3)
  openxlsx::addStyle(wb, sheet_name, style_knn,    rows = 3,  cols = legend_col)
  
  openxlsx::writeData(wb, sheet_name, paste0("MNAR threshold: ", min_missing_pct * 100, "% missing per group"),
                      startCol = legend_col, startRow = 5)
  
  # Save
  openxlsx::saveWorkbook(wb, output_file, overwrite = TRUE)
  message(paste("   -> Saved imputation map:", output_file))
}

# --- Helper: Safely convert a human-readable contrast string to a valid R name ---
# e.g. "WT - (IgG + BeadOnly)/2"  -->  "WT_vs_IgGBeadOnly2"
make_safe_contrast_name <- function(x) {
  gsub("[^[:alnum:]_]", "", gsub(" - ", "_vs_", x))
}

# --- Helper: Load a checkpoint and guarantee de_results are named ---
load_checkpoint <- function(ckpt_file) {
  if (!file.exists(ckpt_file)) {
    message(paste("   -> Skipping: No checkpoint found at", ckpt_file))
    return(NULL)
  }
  ckpt <- readRDS(ckpt_file)
  
  # Use the stored safe names if the live names were dropped
  if (is.null(names(ckpt$de_results)) || any(names(ckpt$de_results) == "")) {
    message("   -> Warning: de_results names were missing. Restoring from stored safe names.")
    if (!is.null(ckpt$de_results_safe_names)) {
      names(ckpt$de_results) <- ckpt$de_results_safe_names
    } else {
      stop(paste("Cannot restore names: checkpoint is missing de_results_safe_names.", ckpt_file))
    }
  }
  return(ckpt)
}

add_rank_plot <- function(imputed_matrix, exp_design, bait, poi_list, output_path, exp_name) {
  
  # Use only non-background samples
  sample_cols <- exp_design$label
  
  rank_df <- data.frame(
    Protein    = rownames(imputed_matrix),
    MeanIntensity = rowMeans(imputed_matrix[, sample_cols, drop = FALSE], na.rm = TRUE)
  ) %>%
    dplyr::arrange(dplyr::desc(MeanIntensity)) %>%
    dplyr::mutate(
      Rank  = dplyr::row_number(),
      Class = dplyr::case_when(
        Protein == bait         ~ "Bait",
        Protein %in% poi_list   ~ "POI",
        TRUE                    ~ "Other"
      )
    )
  
  label_df <- dplyr::filter(rank_df, Class != "Other")
  
  p <- ggplot2::ggplot(rank_df, ggplot2::aes(x = Rank, y = MeanIntensity)) +
    ggplot2::geom_point(ggplot2::aes(colour = Class), alpha = 0.5, size = 1.2) +
    ggrepel::geom_text_repel(
      data         = label_df,
      ggplot2::aes(label = Protein, colour = Class),
      size         = 3, fontface = "bold", max.overlaps = Inf
    ) +
    ggplot2::scale_colour_manual(values = c("Bait" = "red", "POI" = "navy", "Other" = "grey70")) +
    ggplot2::theme_bw() +
    ggplot2::labs(
      title = paste(exp_name, "- Protein Abundance Rank"),
      x     = "Protein Rank (by Mean Intensity)",
      y     = "Mean VSN Intensity"
    )
  
  save_and_print_plot(p, output_path)
}



generate_qc_report <- function(exp_name, output_base_dir, config) {
  
  if (!requireNamespace("htmltools", quietly = TRUE)) {
    message("   -> htmltools not installed. Skipping HTML report.")
    return(invisible(NULL))
  }
  
  safe_path <- function(x) gsub("\\\\", "/", x)
  
  exp_dir  <- file.path(output_base_dir, exp_name)
  qc_dir   <- file.path(exp_dir, "1_QC_and_Normalization")
  de_dir   <- file.path(exp_dir, "2_Differential_Expression")
  logs_dir <- file.path(exp_dir, "0_Logs_and_Exclusions")
  
  ckpt_path <- file.path(exp_dir, paste0(exp_name, "_checkpoint.rds"))
  if (!file.exists(ckpt_path)) {
    message(paste("   -> Skipping QC report: no checkpoint found for", exp_name))
    return(invisible(NULL))
  }
  ckpt <- readRDS(ckpt_path)
  
  # ── Helper: encode TIFF as base64 PNG for self-contained embedding ─────────
  if (!requireNamespace("magick", quietly = TRUE)) {
    install.packages("magick")
  }
  
  encode_image <- function(path) {
    if (is.null(path) || !file.exists(path)) return(NULL)
    
    # Read the TIFF from disk
    img <- tryCatch(magick::image_read(path), error = function(e) NULL)
    if (is.null(img)) return(NULL)
    
    # Convert to PNG in memory and write to a raw vector
    img_png <- magick::image_convert(img, "png")
    raw_png <- magick::image_write(img_png, format = "png")
    
    paste0("data:image/png;base64,", base64enc::base64encode(raw_png))
  }
  
  make_image <- function(tiff_path, width = "90%") {
    encoded <- encode_image(tiff_path)
    if (is.null(encoded)) return("<p><em>Image not found or conversion failed.</em></p>")
    paste0('<img src="', encoded, '" style="width:', width, '; display:block; margin:auto;">')
  }
  
  # ── Collect tiff plots ─────────────────────────────────────────────────────────
  qc_tiffs  <- sort(list.files(qc_dir, pattern = "\\.tiff$",           full.names = TRUE))
  de_tiffs  <- sort(list.files(de_dir, pattern = "^Volcano.*\\.tiff$", full.names = TRUE))
  ma_tiffs  <- sort(list.files(de_dir, pattern = "^MA_.*\\.tiff$",     full.names = TRUE))

  # ── Precompute all tables ─────────────────────────────────────────────────────
  n_samples   <- nrow(ckpt$exp_design_clean)
  n_proteins  <- nrow(ckpt$imputed_matrix)
  n_contrasts <- length(ckpt$de_results)
  
  excl_files   <- list.files(logs_dir, pattern = "Excluded_.*\\.csv$", full.names = TRUE)
  excl_summary <- lapply(excl_files, function(f) {
    df <- tryCatch(readr::read_csv(f, show_col_types = FALSE), error = function(e) NULL)
    if (is.null(df)) return(NULL)
    data.frame(Filter_Step = gsub("Excluded_|\\.csv", "", basename(f)),
               Proteins_Removed = nrow(df))
  }) %>% dplyr::bind_rows()
  
  threshold_log_path <- file.path(logs_dir, "Filter_Thresholds_MissingValues.csv")
  threshold_log <- if (file.exists(threshold_log_path))
    readr::read_csv(threshold_log_path, show_col_types = FALSE) else NULL
  
  # --- Replace the sig_summary block around line ~660 ---
  if (length(ckpt$de_results) > 0 && !is.null(names(ckpt$de_results))) {
    sig_summary <- lapply(names(ckpt$de_results), function(cn) {
      df <- ckpt$de_results[[cn]]
      
      data.frame(
        Contrast = cn,
        Quant_Hits = sum(df$Status == "Quantitative" & df$adj.P.Val < config$thresholds$fdr & abs(df$logFC) > config$thresholds$lfc, na.rm=TRUE),
        Tier1_HighConf = sum(grepl("(High Confidence)", df$Status, fixed = TRUE), na.rm=TRUE),
        Tier2_Exploratory = sum(grepl("(Exploratory Candidate)", df$Status, fixed = TRUE), na.rm=TRUE)
      )
    }) %>% dplyr::bind_rows() %>% 
      # CRITICAL FIX: Exclude Tier 2 from formal significance summation
      dplyr::mutate(Total_Significant = Quant_Hits + Tier1_HighConf)
  } else {
    sig_summary <- NULL
  }
  
  cv_summary <- tryCatch({
    lapply(unique(ckpt$exp_design_clean$condition), function(cond) {
      cols <- ckpt$exp_design_clean$label[ckpt$exp_design_clean$condition == cond]
      if (length(cols) < 2) return(NULL)
      mat <- ckpt$norm_matrix[, cols, drop = FALSE]
      cvs <- (apply(mat, 1, sd, na.rm=TRUE) / rowMeans(mat, na.rm=TRUE)) * 100
      data.frame(Condition = cond,
                 Median_CV = round(median(cvs, na.rm=TRUE), 1),
                 Pct_Below_20 = round(mean(cvs < 20, na.rm=TRUE) * 100, 1))
    }) %>% dplyr::bind_rows()
  }, error = function(e) NULL)
  
  mv_summary <- lapply(unique(ckpt$exp_design_clean$condition), function(cond) {
    cols <- ckpt$exp_design_clean$label[ckpt$exp_design_clean$condition == cond]
    mat  <- ckpt$norm_matrix[, cols, drop = FALSE]
    data.frame(Condition = cond, N_Samples = length(cols),
               Median_Missing_Pct = round(median(rowMeans(is.na(mat)) * 100), 1))
  }) %>% dplyr::bind_rows()
  
  poi_list     <- config$poi_list
  poi_detected <- base::intersect(poi_list, rownames(ckpt$imputed_matrix))
  poi_missing  <- base::setdiff(poi_list, rownames(ckpt$imputed_matrix))
  
  poi_stats <- if (length(poi_detected) > 0) {
    lapply(poi_detected, function(prot) {
      lapply(names(ckpt$de_results), function(cn) {
        df  <- ckpt$de_results[[cn]]
        hit <- df[df$Protein == prot, ]
        if (nrow(hit) == 0) return(data.frame(Protein=prot, Contrast=cn, logFC=NA, FDR=NA, Significant=FALSE))
        data.frame(Protein=prot, Contrast=cn,
                   logFC=round(hit$logFC[1],2), FDR=signif(hit$adj.P.Val[1],3),
                   Significant=hit$adj.P.Val[1] < config$thresholds$fdr & abs(hit$logFC[1]) > config$thresholds$lfc)
      }) %>% dplyr::bind_rows()
    }) %>% dplyr::bind_rows()
  } else NULL
  
  # ── HTML helper functions ─────────────────────────────────────────────────────
  
  # Convert a data.frame to an HTML table string
  df_to_html_table <- function(df, highlight_col = NULL) {
    if (is.null(df) || nrow(df) == 0) return("<p><em>No data available.</em></p>")
    
    header <- paste0("<thead><tr>",
                     paste0("<th>", colnames(df), "</th>", collapse=""),
                     "</tr></thead>")
    
    rows <- apply(df, 1, function(row) {
      highlight <- !is.null(highlight_col) && 
        highlight_col %in% names(row) && 
        isTRUE(as.logical(row[[highlight_col]]))
      style <- if (highlight) ' style="background:#FFF2CC"' else ""
      paste0("<tr", style, ">",
             paste0("<td>", row, "</td>", collapse=""),
             "</tr>")
    })
    
    paste0('<div class="table-wrapper"><table class="summary-table"><', 
           header, "<tbody>", paste(rows, collapse=""), 
           "</tbody></table></div>")
  }
  
  # Wrap a section in a collapsible <details> block
  # This directly answers Q1 — code/content folding per section
  make_section <- function(title, content, guidance = NULL, open = FALSE) {
    open_attr <- if (open) " open" else ""
    guidance_html <- if (!is.null(guidance)) {
      paste0('<div class="guidance"><p>', guidance, '</p></div>')
    } else ""
    paste0(
      '<details', open_attr, '>',
      '<summary class="section-header">', title, '</summary>',
      '<div class="section-content">',
      guidance_html,
      content,
      '</div>',
      '</details>'
    )
  }
  
  # Embed a PNG as a base64 image — completely self-contained, no path resolution needed
  make_image <- function(png_path, width = "90%") {
    encoded <- encode_image(png_path)
    if (is.null(encoded)) return("<p><em>Image not found.</em></p>")
    paste0('<img src="', encoded, '" style="width:', width, '; display:block; margin:auto;">')
  }
  
  # ── QC guidance text ──────────────────────────────────────────────────────────
  qc_guidance <- list(
    "QC0"  = "Raw intensity distributions before normalisation. Boxes should be roughly aligned in median and spread. A sample sitting more than ~1 unit below the median of others is a candidate for removal.",
    "QC1"  = "Protein detection numbers per sample. Bars should be consistent within conditions. A sample detecting less than 60% of the proteins seen in its replicates likely has a preparation problem.",
    "QC2"  = "Missing value pattern. Grey tiles are missing values. A well-behaved experiment shows missing values concentrated in low-abundance proteins. Column-wide stripes indicate a failed sample.",
    "QC3"  = "VSN normalisation density check. Density curves should overlap tightly after normalisation. Persistent offsets suggest batch effects requiring additional correction.",
    "QC4"  = "VSN mean-variance fit. The running SD should be approximately flat across the intensity range. An upward trend at low intensities indicates incomplete variance stabilisation.",
    "QC5"  = "Per-sample imputation shift. The red curve should show a small left-shifted shoulder. A very large red peak means a high fraction of values were imputed.",
    "QC6"  = "Top 100 variable proteins heatmap (z-scored). Replicates of the same condition should cluster together.",
    "PCA"  = "PCA scores plot. Replicates should cluster tightly; conditions should separate. Check the loadings plots to understand which proteins drive the separation.",
    "QC8"  = "Spearman correlation matrices. Replicate correlations should be above 0.95 for a clean IP-MS experiment. The delta heatmap shows imputation-induced inflation.",
    "QC9"  = "Replicate CV. The red dashed line marks 20% CV (industry benchmark). Median CVs below 20% indicate good quantitative reproducibility.",
    "QC10" = "Protein abundance rank plot. The bait protein (red) should appear near the top. POIs (navy) are labelled."
  )
  
  # ── Build HTML body ───────────────────────────────────────────────────────────
  
  # Overview section — open by default
  overview_content <- paste0(
    '<table class="summary-table"><tbody>',
    '<tr><td><strong>Experiment</strong></td><td>', exp_name, '</td></tr>',
    '<tr><td><strong>Date</strong></td><td>', format(Sys.Date(), "%d %B %Y"), '</td></tr>',
    '<tr><td><strong>Imputation</strong></td><td>', config$experiments[[exp_name]]$imputation, '</td></tr>',
    '<tr><td><strong>Bait protein</strong></td><td>', config$experimental_design$bait_protein, '</td></tr>',
    '<tr><td><strong>FDR threshold</strong></td><td>', config$thresholds$fdr, '</td></tr>',
    '<tr><td><strong>LFC threshold</strong></td><td>&plusmn;', config$thresholds$lfc, '</td></tr>',
    '<tr><td><strong>Samples passing QC</strong></td><td>', n_samples, '</td></tr>',
    '<tr><td><strong>Proteins quantified</strong></td><td>', n_proteins, '</td></tr>',
    '<tr><td><strong>Contrasts tested</strong></td><td>', n_contrasts, '</td></tr>',
    '</tbody></table>'
  )
  
  sections <- list()
  sections[["overview"]] <- make_section("Experiment Overview", overview_content, open = TRUE)
  
  # Filtering funnel
  sections[["filtering"]] <- make_section(
    "Filtering Funnel",
    df_to_html_table(excl_summary),
    guidance = "Proteins removed at each QC step, from the initial detected set to the final quantified matrix."
  )
  
  # Missing value thresholds
  sections[["mv_thresholds"]] <- make_section(
    "Missing Value Filter Thresholds",
    df_to_html_table(threshold_log),
    guidance = "Exact proportional thresholds applied per clone. A protein was retained if it met the valid-value requirement in enough clones within at least one genotype."
  )
  
  # Missing value summary
  sections[["mv_summary"]] <- make_section(
    "Missing Value Summary Per Condition",
    df_to_html_table(mv_summary),
    guidance = "Median percentage missing per condition before imputation. Values above 50% suggest the condition may be underpowered."
  )
  
  # Replicate quality
  sections[["cv"]] <- make_section(
    "Replicate Quality (CV)",
    df_to_html_table(cv_summary),
    guidance = "Median CV < 20% is the standard mass spectrometry benchmark. Pct_Below_20 shows what fraction of proteins meet this threshold."
  )
  
  # QC plots — one collapsible section per plot
  qc_plot_sections <- lapply(qc_tiffs, function(png_path) {
    plot_base  <- tools::file_path_sans_ext(basename(png_path))
    plot_title <- gsub("_", " ", plot_base)
    qc_key     <- regmatches(plot_base, regexpr("QC[0-9]+|PCA", plot_base))
    guidance   <- if (length(qc_key) > 0 && qc_key %in% names(qc_guidance))
      qc_guidance[[qc_key]] else NULL
    make_section(plot_title, make_image(png_path), guidance = guidance)
  })
  sections[["qc_plots"]] <- paste0(
    '<details><summary class="section-header top-level">QC Plots</summary>',
    '<div class="section-content">',
    paste(qc_plot_sections, collapse = "\n"),
    '</div></details>'
  )
  
  # DE summary table
  sections[["de_summary"]] <- make_section(
    "Differential Expression Summary",
    df_to_html_table(sig_summary),
    guidance = paste0("Significant hits per contrast (FDR < ", config$thresholds$fdr,
                      ", |LFC| > ", config$thresholds$lfc, ").")
  )
  
  # Volcano plots — one collapsible per plot
  volcano_sections <- lapply(de_tiffs, function(png_path) {
    plot_title <- gsub("_", " ", tools::file_path_sans_ext(basename(png_path)))
    make_section(plot_title, make_image(png_path))
  })
  sections[["volcanos"]] <- paste0(
    '<details><summary class="section-header top-level">Volcano Plots</summary>',
    '<div class="section-content">',
    '<p>Gold diamonds = proteins of interest. Dashed lines = FDR and LFC thresholds.</p>',
    paste(volcano_sections, collapse = "\n"),
    '</div></details>'
  )
  
  # MA plots — one collapsible per plot
  if (length(ma_tiffs) > 0) {
    ma_sections <- lapply(ma_tiffs, function(png_path) {
      plot_title <- gsub("_", " ", tools::file_path_sans_ext(basename(png_path)))
      make_section(plot_title, make_image(png_path),
                   guidance = "The loess line should sit near y=0. An upward curve at low intensities indicates intensity-dependent bias.")
    })
    sections[["ma_plots"]] <- paste0(
      '<details><summary class="section-header top-level">MA Plots</summary>',
      '<div class="section-content">',
      paste(ma_sections, collapse = "\n"),
      '</div></details>'
    )
  }
  
  # POI summary
  poi_content <- paste0(
    "<p><strong>Detected (", length(poi_detected), "/", length(poi_list), "):</strong> ",
    if (length(poi_detected) > 0) paste(poi_detected, collapse = ", ") else "none", "</p>"
  )
  if (length(poi_missing) > 0) {
    poi_content <- paste0(poi_content,
                          '<p class="warning"><strong>Not detected:</strong> ', 
                          paste(poi_missing, collapse = ", "), 
                          " — absent after all filtering steps.</p>")
  }
  if (!is.null(poi_stats) && nrow(poi_stats) > 0) {
    poi_content <- paste0(poi_content, df_to_html_table(poi_stats, highlight_col = "Significant"))
  }
  sections[["poi"]] <- make_section(
    "Proteins of Interest",
    poi_content,
    guidance = paste0("Rows highlighted in yellow are significant (FDR < ", 
                      config$thresholds$fdr, ", |LFC| > ", config$thresholds$lfc, ").")
  )
  
  # Session info
  si_text <- paste(capture.output(sessionInfo()), collapse = "\n")
  sections[["session"]] <- make_section(
    "Reproducibility",
    paste0('<pre class="session-info">', si_text, '</pre>')
  )
  
  # ── CSS ───────────────────────────────────────────────────────────────────────
  css <- '
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
           max-width: 1100px; margin: 0 auto; padding: 20px; color: #333; }
    h1   { color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 8px; }
    
    details { margin: 6px 0; border: 1px solid #ddd; border-radius: 6px; overflow: hidden; }
    details > details { margin: 4px 8px; border-color: #e8e8e8; }
    
    summary.section-header {
      padding: 10px 16px; cursor: pointer; font-weight: 600;
      background: #f0f4f8; color: #2c3e50;
      list-style: none; display: flex; align-items: center; gap: 8px;
    }
    summary.section-header.top-level {
      background: #2c3e50; color: white; font-size: 1.05em; padding: 12px 16px;
    }
    summary.section-header::before { content: "▶"; font-size: 0.8em; transition: transform 0.2s; }
    details[open] > summary.section-header::before { transform: rotate(90deg); }
    
    .section-content { padding: 12px 16px; }
    .guidance { background: #f8f9fa; border-left: 4px solid #4a90d9;
                padding: 8px 12px; margin-bottom: 12px; border-radius: 0 4px 4px 0;
                font-size: 0.92em; color: #555; }
    .warning  { background: #fff3cd; border-left: 4px solid #ffc107;
                padding: 8px 12px; border-radius: 0 4px 4px 0; }
    
    .table-wrapper { overflow-x: auto; margin: 8px 0; }
    table.summary-table { border-collapse: collapse; width: 100%; font-size: 0.9em; }
    table.summary-table th { background: #2c3e50; color: white; padding: 8px 12px; text-align: left; }
    table.summary-table td { padding: 6px 12px; border-bottom: 1px solid #eee; }
    table.summary-table tr:hover td { background: #f5f5f5; }
    
    img { border-radius: 4px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); margin: 8px 0; }
    pre.session-info { background: #f8f8f8; padding: 12px; border-radius: 4px;
                       font-size: 0.8em; overflow-x: auto; white-space: pre-wrap; }
  '
  
  # ── Assemble final HTML ───────────────────────────────────────────────────────
  html <- paste0(
    '<!DOCTYPE html><html lang="en"><head>',
    '<meta charset="UTF-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1.0">',
    '<title>QC Report: ', exp_name, '</title>',
    '<style>', css, '</style>',
    '</head><body>',
    '<h1>Proteomics QC Report: ', exp_name, '</h1>',
    '<p style="color:#888; font-size:0.9em;">Generated ', 
    format(Sys.Date(), "%d %B %Y"), ' &bull; Hybrid DEP/Limma Pipeline</p>',
    paste(sections, collapse = "\n"),
    '</body></html>'
  )
  
  # ── Write ─────────────────────────────────────────────────────────────────────
  out_path <- file.path(exp_dir, paste0(exp_name, "_QC_Report.html"))
  writeLines(html, out_path)
  message(paste("   -> QC Report generated:", safe_path(out_path)))
}

generate_correlation_matrices <- function(se_norm, se_imputed, exp_design, 
                                          output_dir, exp_name) {
  
  mat_norm    <- SummarizedExperiment::assay(se_norm)
  mat_imputed <- SummarizedExperiment::assay(se_imputed)
  
  # ── Helper: build annotated heatmap ──────────────────────────────────────────
  make_cor_heatmap <- function(mat, title, method = "spearman", 
                               use_pairwise = TRUE) {
    
    use_arg <- if (use_pairwise) "pairwise.complete.obs" else "complete.obs"
    cor_mat <- cor(mat, method = method, use = use_arg)
    
    # Annotation bar showing condition per sample
    anno_df  <- data.frame(
      Condition = exp_design$condition,
      row.names = colnames(mat)
    )
    
    # Colour palette: white at ~0.9 so differences in the high-correlation
    # range are visible — don't use the full 0-1 range or everything looks similar
    col_fun <- colorRampPalette(c("#2166AC", "#F7F7F7", "#D6604D"))(100)
    
    # Round values for display
    display_mat <- round(cor_mat, 2)
    
    p <- pheatmap::pheatmap(
      cor_mat,
      main            = title,
      display_numbers = display_mat,
      number_format   = "%.2f",
      fontsize_number = 7,
      color           = col_fun,
      breaks          = seq(0.7, 1.0, length.out = 101),
      annotation_col  = anno_df,
      clustering_method = "ward.D2",
      angle_col       = 45,
      border_color    = "white",
      silent          = TRUE
    )
    return(p)
  }
  
  # ── 1. Pre-imputation Spearman (pairwise — keeps real observed data only) ─────
  message("   -> Generating pre-imputation Spearman correlation matrix...")
  p_pre_spearman <- make_cor_heatmap(
    mat    = mat_norm,
    title  = paste(exp_name, "- Pre-Imputation (Spearman, pairwise)"),
    method = "spearman",
    use_pairwise = TRUE
  )
  
  save_and_print_plot(p_pre_spearman, 
                      file.path(output_dir, "QC8a_Correlation_PreImputation_Spearman"),
                      width = 14, height = 12)
  
  # ── 2. Post-imputation Spearman (complete matrix — the standard QC plot) ──────
  message("   -> Generating post-imputation Spearman correlation matrix...")
  p_post_spearman <- make_cor_heatmap(
    mat    = mat_imputed,
    title  = paste(exp_name, "- Post-Imputation (Spearman, complete)"),
    method = "spearman",
    use_pairwise = FALSE
  )
  
  save_and_print_plot(p_post_spearman,
                      file.path(output_dir, "QC8b_Correlation_PostImputation_Spearman"),
                      width = 14, height = 12)
  
  # ── 3. Pre vs Post Pearson delta — detects imputation-inflated correlations ───
  # If imputation artificially inflates replicate similarity, Pearson correlations
  # will increase substantially between pre and post matrices.
  # A delta > 0.05 for any pair warrants a note in your methods.
  message("   -> Calculating imputation correlation inflation...")
  
  cor_pre_pearson  <- cor(mat_norm,    method = "pearson", use = "pairwise.complete.obs")
  cor_post_pearson <- cor(mat_imputed, method = "pearson", use = "complete.obs")
  
  delta_mat <- cor_post_pearson - cor_pre_pearson
  
  # Only report non-diagonal pairs
  pairs_df <- as.data.frame(as.table(delta_mat)) %>%
    dplyr::rename(Sample_A = Var1, Sample_B = Var2, Delta_r = Freq) %>%
    dplyr::filter(as.character(Sample_A) < as.character(Sample_B)) %>%
    dplyr::left_join(
      exp_design %>% dplyr::select(label, condition) %>% 
        dplyr::rename(Sample_A = label, Cond_A = condition),
      by = "Sample_A"
    ) %>%
    dplyr::left_join(
      exp_design %>% dplyr::select(label, condition) %>% 
        dplyr::rename(Sample_B = label, Cond_B = condition),
      by = "Sample_B"
    ) %>%
    dplyr::mutate(
      Pair_Type = ifelse(Cond_A == Cond_B, "Within-condition replicate", "Cross-condition"),
      Flag      = Delta_r > 0.05
    ) %>%
    dplyr::arrange(dplyr::desc(Delta_r))
  
  readr::write_csv(pairs_df, 
                   file.path(output_dir, "QC8c_Imputation_Correlation_Inflation.csv"))
  
  n_flagged <- sum(pairs_df$Flag)
  within_flagged <- sum(pairs_df$Flag & pairs_df$Pair_Type == "Within-condition replicate")
  cross_flagged  <- sum(pairs_df$Flag & pairs_df$Pair_Type == "Cross-condition")
  
  if (n_flagged > 0) {
    message(paste0("   -> Warning: ", n_flagged, " sample pairs show >0.05 Pearson r inflation ",
                   "after imputation (",
                   within_flagged, " within-condition, ",
                   cross_flagged,  " cross-condition). ",
                   "Check QC8c_Imputation_Correlation_Inflation.csv"))
  }
  
  # ── 4. Delta heatmap ──────────────────────────────────────────────────────────
  p_delta <- pheatmap::pheatmap(
    delta_mat,
    main              = paste(exp_name, "- Correlation Inflation from Imputation (Pearson delta)"),
    display_numbers   = round(delta_mat, 3),
    number_format     = "%.3f",
    fontsize_number   = 7,
    color             = colorRampPalette(c("white", "#FFF2CC", "#E41A1C"))(100),
    breaks            = seq(0, 0.1, length.out = 101),
    clustering_method = "ward.D2",
    angle_col         = 45,
    border_color      = "white",
    silent            = TRUE
  )
  
  save_and_print_plot(p_delta,
                      file.path(output_dir, "QC8d_Correlation_Inflation_Delta"),
                      width = 9, height = 8)
}

#==============================================================================#
# Non-Imputed Dual-Track Statistical Testing (With POI Labeling)
#==============================================================================#

calculate_raw_abundance_stats <- function(mat_scrubbed, exp_design, contrast_str, 
                                          output_base_dir, exp_name, id_map = NULL,
                                          poi_list = c(), n_label_nonsig = 15, blacklist_proteins = c(), generate_plots = TRUE) {
  
  groups <- stringr::str_trim(stringr::str_split(contrast_str, " - ")[[1]])
  target_cond  <- groups[1]
  control_cond <- paste(groups[-1], collapse = " - ") # Safely collapses if multiple "-" exist
  
  # 2. Smart Parser to handle composite contrasts like "(IgG + BeadOnly)/2"
  parse_conditions <- function(cond_str) {
    # Extract any alphanumeric word from the string
    words <- stringr::str_extract_all(cond_str, "[a-zA-Z0-9_]+")[[1]]
    # Match against actual condition names in the design
    valid_conds <- base::intersect(words, unique(exp_design$condition))
    return(exp_design$label[exp_design$condition %in% valid_conds])
  }
  
  target_cols  <- parse_conditions(target_cond)
  control_cols <- parse_conditions(control_cond)
  
  if (length(target_cols) == 0 || length(control_cols) == 0) {
    stop(paste("Could not parse valid samples from contrast:", contrast_str))
  }
  
  safe_contrast <- make_safe_contrast_name(contrast_str)
  nested_out_dir <- file.path(output_base_dir, exp_name, "2_Differential_Abundance")
  dir.create(nested_out_dir, recursive = TRUE, showWarnings = FALSE)
  
    # --- ADD THESE BEFORE THE LOOP ---
    # 1. Dynamically estimate the empirical instrument noise floor (1st quantile of all data)
    all_valid_vals <- as.numeric(mat_scrubbed)[!is.na(as.numeric(mat_scrubbed)) & as.numeric(mat_scrubbed) > 0]
    noise_floor    <- quantile(all_valid_vals, 0.01, na.rm = TRUE)
    
    # 2. Define parameters
    min_binary_dist      <- config$pipeline_parameters$min_binary_distance %||% 2.0
    min_exploratory_reps <- config$pipeline_parameters$min_exploratory_replicates %||% 2
    
    res_list <- list()
    
    # --- START THE LOOP ---
    for (i in 1:nrow(mat_scrubbed)) {
      prot  <- rownames(mat_scrubbed)[i]
      vals_T <- as.numeric(mat_scrubbed[i, target_cols])
      vals_C <- as.numeric(mat_scrubbed[i, control_cols])
      
      clean_T <- vals_T[!is.na(vals_T)]
      clean_C <- vals_C[!is.na(vals_C)]
      
      n_T <- length(clean_T)
      n_C <- length(clean_C)
      
      mean_T <- mean(clean_T)
      mean_C <- mean(clean_C)
      logFC  <- mean_T - mean_C 
      pval   <- NA
      status <- "Not Enough Data"
      
      # Track 1: Quantitative Shift
      if (n_T >= 2 && n_C >= 2) {
        if (var(clean_T) == 0 && var(clean_C) == 0) {
          pval <- 1
        } else {
          t_res <- tryCatch(t.test(clean_T, clean_C, var.equal = FALSE), 
                            error = function(e) list(p.value = NA))
          pval <- t_res$p.value
        }
        status <- "Quantitative"
      } 
      
      # Track 2a: Tier 1 - High-Confidence Binary On/Off (100% Replicates + Noise Floor Gate)
      else if (n_T == length(target_cols) && n_C == 0 && (mean_T - noise_floor) >= min_binary_dist) {
        status <- paste("Exclusive to", target_cond, "(High Confidence)")
        logFC  <- (mean_T - noise_floor)  
        pval   <- 0   
      } 
      else if (n_T == 0 && n_C == length(control_cols) && (mean_C - noise_floor) >= min_binary_dist) {
        status <- paste("Exclusive to", control_cond, "(High Confidence)")
        logFC  <- -(mean_C - noise_floor) 
        pval   <- 0   
      }
      
      # Track 2b: Tier 2 - Exploratory Candidate Binary On/Off (Partial Replicates OR Near Noise Floor)
      else if (n_T >= min_exploratory_reps && n_C == 0) {
        status <- paste("Exclusive to", target_cond, "(Exploratory Candidate)")
        logFC  <- (mean_T - noise_floor)
        pval   <- 0.01 # Nominal p-value ensures they don't get FDR-starred as validated hits
      } 
      else if (n_T == 0 && n_C >= min_exploratory_reps) {
        status <- paste("Exclusive to", control_cond, "(Exploratory Candidate)")
        logFC  <- -(mean_C - noise_floor)   
        pval   <- 0.01 
      }
      
      res_list[[i]] <- data.frame(
        Protein = prot,
        Status  = status,
        Target_Condition = target_cond,
        Control_Condition = control_cond,
        N_Target = n_T,
        N_Control = n_C,
        Mean_Target = mean_T,
        Mean_Control = mean_C,
        logFC = logFC,
        P_Value = pval
      )
    }
  
  res_df <- dplyr::bind_rows(res_list)
  
  if (!is.null(id_map)) {
    res_df <- res_df %>% dplyr::left_join(id_map, by = "Protein")
  }
  
  res_df$adj.P.Val <- NA
  quant_idx <- res_df$Status == "Quantitative" & !is.na(res_df$P_Value)
  
  if (sum(quant_idx) > 0) {
    res_df$adj.P.Val[quant_idx] <- p.adjust(res_df$P_Value[quant_idx], method = "BH")
  }
  
  # 1. High-Confidence dropouts act as validated hits (FDR = 0)
  res_df$adj.P.Val[grepl("(High Confidence)", res_df$Status, fixed = TRUE)] <- 0
  
  # 2. Exploratory dropouts get FDR = 1 so they fail strict statistical thresholds
  #    but remain available for visual inspection in Sidecar plots
  res_df$adj.P.Val[grepl("(Exploratory Candidate)", res_df$Status, fixed = TRUE)] <- 1
  
  readr::write_csv(res_df, file.path(nested_out_dir, paste0(safe_contrast, "_Raw_Stats.csv")))
  
  # ============================================================================
  # --- Visualization ---
  # ============================================================================
  if (generate_plots) {
    
    # Render the uncompressed 3-panel Sidecar Plot instead of the squashed 1-panel plot
    p_sidecar <- generate_uncompressed_sidecar_plot(
      res_df             = res_df,
      contrast_str       = contrast_str,
      exp_name           = exp_name,
      fdr_thr            = config$thresholds$fdr,
      lfc_thr            = config$thresholds$lfc,
      poi_list           = poi_list,
      use_nominal_p      = FALSE,
      noise_floor        = noise_floor,       # Automatically passes the matrix noise floor!
      n_label_nonsig     = n_label_nonsig,
      blacklist_proteins = blacklist_proteins
    )
    
    p_volc <- generate_universal_volcano_plot(
      res_df             = res_df,
      contrast_str       = contrast_str,
      exp_name           = exp_name,
      fdr_thr            = config$thresholds$fdr,
      lfc_thr            = config$thresholds$lfc,
      poi_list           = poi_list,
      use_nominal_p      = TRUE,
      n_label            = n_label_nonsig,
      blacklist_proteins = blacklist_proteins
    )
    
    # Print to R view and save to nested TIFF folder
    safe_name <- make_safe_contrast_name(contrast_str)
    
    # Increase the width from 10 to 14 to give the sidecars breathing room
    save_and_print_plot(p_volc, file.path(nested_out_dir, paste0("Volcano_", safe_name)), width = 14, height = 8)
    save_and_print_plot(p_sidecar, file.path(nested_out_dir, paste0("Sidecar_", safe_name)), width = 14, height = 8)
    
  }
  
  return(res_df)
}

#==============================================================================#
# Presence-threshold scrubbing + non-imputed abundance reporting
# for interactome (IP-MS) data
#
# Rationale: imputation assumes missingness is technical (MAR/MNAR noise to be
# filled in). For an interactome experiment, absence of a bait interactor in a
# given genotype/clone can be biologically real. Instead of imputing, this
# enforces a minimum-evidence rule at two nested levels and reports raw
# (VSN-normalized, log2-scale) abundances -- no fold-change, no limma contrast.
#==============================================================================#

#' Apply clone- and genotype-level presence thresholds to a normalized matrix
#'
#' Level 1 (technical replicates within a clone):
#'   If fewer than `min_valid_prop` of a clone's replicates have a valid
#'   (non-NA) value for a protein, ALL values for that clone/protein are set
#'   to NA. A single surviving replicate is not treated as representative.
#'
#' Level 2 (clones within a genotype):
#'   After Level 1, if fewer than `min_clone_prop` of the genotype's clones
#'   still have any valid value for that protein, the ENTIRE genotype is set
#'   to NA for that protein.
#'
#' @param mat Numeric matrix, proteins (rows) x samples (columns). Must be the
#'   VSN-normalized but NOT imputed matrix (i.e. real missingness preserved).
#'   vsn output is already on an approximately log2 (glog) scale, so no
#'   further transform is needed.
#' @param exp_design Data frame with `label` (matching colnames(mat)),
#'   `condition` (genotype), and ideally `clone_id`. If `clone_id` is absent,
#'   each replicate is treated as its own clone (Level 1 becomes a no-op;
#'   only the genotype-level threshold has effect).
#' @param min_valid_prop Minimum proportion of replicates required within a
#'   clone. Default 0.5 (i.e. >=2 of 3 must be present).
#' @param min_clone_prop Minimum proportion of clones required within a
#'   genotype. Default 0.5 (i.e. >=2 of 3 clones must have survived Level 1).
#'
#' @return list(mat_scrubbed = matrix with extra NAs applied,
#'              na_log = long data frame documenting every scrub event)
scrub_by_clone_genotype_presence <- function(mat, exp_design,
                                             min_valid_prop = 0.5,
                                             min_clone_prop = 0.5) {
  
  stopifnot(all(exp_design$label %in% colnames(mat)))
  
  has_clones <- "clone_id" %in% colnames(exp_design)
  if (!has_clones) {
    message("   -> No clone_id column found; treating each replicate as its own clone ",
            "(only the genotype-level threshold will have any effect).")
    exp_design <- exp_design
    exp_design$clone_id <- exp_design$label
  }
  
  mat_out <- mat
  na_log  <- vector("list", 0)
  genotypes <- unique(exp_design$condition)
  
  for (geno in genotypes) {
    
    geno_design <- exp_design[exp_design$condition == geno, ]
    clones      <- unique(geno_design$clone_id)
    
    # Per-protein, per-clone survival flag after Level 1
    clone_survives <- matrix(TRUE, nrow = nrow(mat), ncol = length(clones),
                             dimnames = list(rownames(mat), clones))
    
    # ---- Level 1: replicate presence within each clone ----
    for (cl in clones) {
      cl_cols <- geno_design$label[geno_design$clone_id == cl]
      cl_cols <- cl_cols[cl_cols %in% colnames(mat_out)]
      if (length(cl_cols) == 0) next
      
      present  <- rowSums(!is.na(mat_out[, cl_cols, drop = FALSE]))
      required <- ceiling(length(cl_cols) * min_valid_prop)
      fails    <- present < required
      
      if (any(fails)) {
        mat_out[fails, cl_cols]   <- NA
        clone_survives[fails, cl] <- FALSE
        
        na_log[[length(na_log) + 1]] <- data.frame(
          Protein     = rownames(mat)[fails],
          Genotype    = geno,
          Clone       = cl,
          Level       = "Replicate (within clone)",
          N_Present   = present[fails],
          N_Required  = required,
          N_Total     = length(cl_cols)
        )
      }
    }
    
    # ---- Level 2: clone presence within the genotype ----
    n_clones_total  <- length(clones)
    required_clones <- ceiling(n_clones_total * min_clone_prop)
    n_clones_valid  <- rowSums(clone_survives)
    geno_fails      <- n_clones_valid < required_clones
    
    if (any(geno_fails)) {
      geno_cols <- geno_design$label[geno_design$label %in% colnames(mat_out)]
      mat_out[geno_fails, geno_cols] <- NA
      
      na_log[[length(na_log) + 1]] <- data.frame(
        Protein     = rownames(mat)[geno_fails],
        Genotype    = geno,
        Clone       = "ALL",
        Level       = "Clone (within genotype)",
        N_Present   = n_clones_valid[geno_fails],
        N_Required  = required_clones,
        N_Total     = n_clones_total
      )
    }
  }
  
  na_log_df <- if (length(na_log) > 0) {
    dplyr::bind_rows(na_log)
  } else {
    data.frame(Protein = character(), Genotype = character(), Clone = character(),
               Level = character(), N_Present = integer(), N_Required = integer(),
               N_Total = integer())
  }
  
  list(mat_scrubbed = mat_out, na_log = na_log_df)
}


#' Export scrubbed abundance values as a raw abundance report (no fold-change)
#'
#' @param mat_scrubbed Output of scrub_by_clone_genotype_presence()$mat_scrubbed
#' @param exp_design Data frame with label, condition, clone_id, replicate
#' @param na_log Output of scrub_by_clone_genotype_presence()$na_log
#' @param output_dir Directory to write CSVs to
#' @param exp_name Experiment name, used in output filenames
#' @param id_map Optional data frame with Protein/ID -> gene name mapping to
#'   join in for readability (e.g. df_unique %>% select(ID, name))
export_interactome_abundance_report <- function(mat_scrubbed, exp_design, na_log,
                                                output_dir, exp_name,
                                                id_map = NULL) {
  
  dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
  
  long_df <- as.data.frame(mat_scrubbed) %>%
    tibble::rownames_to_column("Protein") %>%
    tidyr::pivot_longer(-Protein, names_to = "label", values_to = "log2_abundance") %>%
    dplyr::left_join(
      exp_design %>% dplyr::select(label, condition, clone_id, replicate),
      by = "label"
    )
  
  if (!is.null(id_map)) {
    long_df <- long_df %>% dplyr::left_join(id_map, by = "Protein")
  }
  
  readr::write_csv(long_df,
                   file.path(output_dir, paste0(exp_name, "_Abundance_LongFormat.csv")))
  
  # Per-clone mean log2 abundance -- a summary view; still raw abundance, no contrast
  clone_means <- long_df %>%
    dplyr::group_by(Protein, condition, clone_id) %>%
    dplyr::summarise(
      mean_log2_abundance = mean(log2_abundance, na.rm = TRUE),
      n_valid_reps        = sum(!is.na(log2_abundance)),
      n_total_reps        = dplyr::n(),
      .groups = "drop"
    ) %>%
    dplyr::mutate(mean_log2_abundance = ifelse(is.nan(mean_log2_abundance),
                                               NA, mean_log2_abundance))
  
  readr::write_csv(clone_means,
                   file.path(output_dir, paste0(exp_name, "_Abundance_CloneMeans.csv")))
  
  wide_out <- as.data.frame(mat_scrubbed) %>% tibble::rownames_to_column("Protein")
  readr::write_csv(wide_out,
                   file.path(output_dir, paste0(exp_name, "_Abundance_Matrix_Wide.csv")))
  
  readr::write_csv(na_log,
                   file.path(output_dir, paste0(exp_name, "_Threshold_NA_Log.csv")))
  
  invisible(list(long = long_df, clone_means = clone_means, wide = wide_out))
}

#' Per-protein abundance strip plot across genotype, colored by clone
#'
#' Shows every surviving replicate value as a jittered point, plus a mean +/-
#' SEM crossbar per genotype. Genotypes that were fully NA'd by
#' scrub_by_clone_genotype_presence() simply show no points -- the gap itself
#' is the result, not something to be filled in.
#'
#' @param long_df Long-format data frame from export_interactome_abundance_report()
#'   (columns: Protein, label, log2_abundance, condition, clone_id, replicate,
#'   optionally name/gene symbol from id_map)
#' @param proteins Character vector of Protein IDs (or gene names, if `name`
#'   column present and name_col = "name") to plot. Keep this short (<= ~12)
#'   for a legible facet grid; use the heatmap for large-scale overviews.
#' @param name_col Column to match `proteins` against and to use as facet
#'   labels. Defaults to "Protein"; set to "name" to match by gene symbol.
#' @param exp_name Used in the plot title
#'
#' @return a ggplot object
plot_protein_abundance_strip <- function(long_df, proteins, name_col = "Protein",
                                         exp_name = "") {
  
  plot_data <- long_df %>%
    dplyr::filter(.data[[name_col]] %in% proteins)
  
  if (nrow(plot_data) == 0) {
    stop("None of the requested proteins were found in long_df's '", name_col, "' column.")
  }
  
  summary_data <- plot_data %>%
    dplyr::group_by(.data[[name_col]], condition) %>%
    dplyr::summarise(
      mean_abund = mean(log2_abundance, na.rm = TRUE),
      sem        = sd(log2_abundance, na.rm = TRUE) / sqrt(sum(!is.na(log2_abundance))),
      .groups = "drop"
    ) %>%
    dplyr::filter(!is.nan(mean_abund))
  
  ggplot2::ggplot(plot_data, ggplot2::aes(x = condition, y = log2_abundance)) +
    ggplot2::geom_jitter(ggplot2::aes(color = clone_id), width = 0.15, height = 0,
                         size = 2, alpha = 0.8, na.rm = TRUE) +
    ggplot2::geom_crossbar(
      data = summary_data,
      ggplot2::aes(y = mean_abund, ymin = mean_abund - sem, ymax = mean_abund + sem),
      width = 0.4, fill = NA, color = "black", linewidth = 0.4
    ) +
    ggplot2::facet_wrap(stats::as.formula(paste0("~", name_col)), scales = "free_y") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = paste(exp_name, "- Presence-Thresholded Abundance"),
      subtitle = "Points = replicates surviving threshold; bar = mean +/- SEM; empty genotype = failed presence threshold",
      y = "log2 abundance (VSN-normalized, non-imputed)",
      x = "", color = "Clone"
    )
}


#' Genotype x clone abundance heatmap, NA cells shown as blank
#'
#' Built on clone-mean log2 abundance (from export_interactome_abundance_report()'s
#' clone_means output), so this reflects the same presence-thresholded values as
#' the strip plot, just at proteome scale.
#'
#' @param clone_means Data frame from export_interactome_abundance_report()
#'   (columns: Protein, condition, clone_id, mean_log2_abundance, ...)
#' @param proteins Optional character vector to subset to (e.g. a POI list or
#'   a set passed the row-level missing-value filter). If NULL, uses all
#'   proteins in clone_means (can be slow/illegible for very large sets).
#' @param exp_name Used in the plot title
#' @param cluster_rows Whether to hierarchically cluster proteins (default TRUE)
#'
#' @return the pheatmap object (also drawn as a side effect)
plot_abundance_heatmap <- function(clone_means, proteins = NULL, exp_name = "",
                                   cluster_rows = TRUE) {
  
  data <- clone_means
  if (!is.null(proteins)) {
    data <- data %>% dplyr::filter(Protein %in% proteins)
  }
  
  wide_mat <- data %>%
    dplyr::mutate(sample_col = paste(condition, clone_id, sep = "_")) %>%
    dplyr::select(Protein, sample_col, mean_log2_abundance) %>%
    tidyr::pivot_wider(names_from = sample_col, values_from = mean_log2_abundance) %>%
    tibble::column_to_rownames("Protein") %>%
    as.matrix()
  
  # --- FIX: Custom Distance Matrix for NAs ---
  # hclust will crash if dist() returns NA. We temporarily substitute NAs 
  # with the matrix minimum just for the distance math.
  if (cluster_rows && nrow(wide_mat) > 1) {
    mat_for_dist <- wide_mat
    min_val <- min(mat_for_dist, na.rm = TRUE)
    mat_for_dist[is.na(mat_for_dist)] <- min_val - 1
    custom_dist <- dist(mat_for_dist)
  } else {
    custom_dist <- "euclidean"
  }
  
  pheatmap::pheatmap(
    wide_mat,
    main              = paste(exp_name, "- Presence-Thresholded Abundance (clone means)"),
    cluster_rows      = cluster_rows,
    clustering_distance_rows = custom_dist,
    cluster_cols      = FALSE,
    na_col            = "grey85",
    show_rownames     = (nrow(wide_mat) <= 100),
    fontsize_row      = 6,             # <-- Shrink font to fit 100 rows cleanly
    fontsize_col      = 8,             # <-- Optional: adjust column header size
    color             = colorRampPalette(c("navy", "white", "firebrick3"))(100),
    border_color      = NA
  )
}



# --- 7. The Master Hybrid DEP Wrapper ---
run_hybrid_dep_workflow <- function(raw_data, exp_name, col_mapping, exp_design, active_design = NULL,
                                    crap_ids, crapome_list, bioid_list, is_bioid = FALSE,
                                    contrasts, output_base_dir, imputation_method = "mixed",
                                    fdr_thr = 0.05, lfc_thr = 1.0, poi_list = c(),
                                    condition_col = "condition", background_controls = c(),
                                    config) {  # Fix 6: config now passed explicitly
  
  message(paste("\n========================================================"))
  message(paste("--- Running Hybrid Pipeline for:", exp_name, "---"))
  message(paste("========================================================"))
  
  out_dir  <- file.path(output_base_dir, exp_name)
  logs_dir <- file.path(out_dir, "0_Logs_and_Exclusions")
  qc_dir   <- file.path(out_dir, "1_QC_and_Normalization")
  de_dir   <- file.path(out_dir, "2_Differential_Expression")
  
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_dir,   recursive = TRUE, showWarnings = FALSE)
  dir.create(de_dir,   recursive = TRUE, showWarnings = FALSE)
  
  # ===========================================================================#
  # 1. INITIAL SETUP & DESIGN WRANGLING ####
  # ===========================================================================#
  colnames(exp_design) <- tolower(gsub('^\xef\xbb\xbf', '', colnames(exp_design)))
  cond_col_lower <- tolower(condition_col)
  
  exp_design <- exp_design %>% dplyr::filter(!is.na(label) & label != "")
  
  if (cond_col_lower != "condition" && cond_col_lower %in% colnames(exp_design)) {
    exp_design <- exp_design %>% dplyr::rename(condition = !!sym(cond_col_lower))
  } else if (!"condition" %in% colnames(exp_design)) {
    stop(paste("CRITICAL ERROR: Column '", condition_col, "' not found in design CSV."))
  }
  
  req_cols <- c("label", "condition", "replicate")
  if (!all(req_cols %in% colnames(exp_design))) stop("exp_design missing required columns.")
  
  exp_design <- exp_design %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(replicate = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(orig_label = label) %>% # <--- FIX: Preserve the original raw headers!
    as.data.frame()
  
  # ===========================================================================#
  # 2. RAW DATA ALIGNMENT & CLEANUP ####
  # ===========================================================================#
  df <- raw_data %>%
    dplyr::rename(ID = !!sym(col_mapping$id_col), name_raw = !!sym(col_mapping$name_col)) %>%
    dplyr::mutate(
      name = stringr::str_split_i(name_raw, ";", 1),
      name = ifelse(is.na(name) | name == "", ID, name)
    )
  
  lfq_cols <- match(exp_design$label, colnames(df))
  if (any(is.na(lfq_cols))) stop("Design CSV labels do not match raw data headers.")
  
  # BULLETPROOF FIX: Use native Base R matrix operations instead of dplyr::across.
  # This completely prevents the "Predicate must return TRUE or FALSE" crash.
  mat_clean <- as.matrix(df[, lfq_cols])
  mat_clean <- apply(mat_clean, 2, as.numeric) # Force numeric
  mat_clean[is.nan(mat_clean) | is.infinite(mat_clean) | mat_clean == 0] <- NA
  df[, lfq_cols] <- mat_clean
  
  # ===========================================================================#
  # 3. EXCLUDE FAILED SAMPLES (With Control-Shielding Safety Net) ####
  # ===========================================================================#
  has_clones_in_design <- "clone_id" %in% colnames(exp_design)
  protein_counts        <- colSums(!is.na(df[, lfq_cols, drop = FALSE]))
  sample_drop_threshold <- config$pipeline_parameters$sample_drop_threshold
  cutoff                <- max(protein_counts) * sample_drop_threshold
  
  # ── 1. Gather all configured control condition names ────────────────────────
  configured_controls <- c(
    config$experimental_design$background_controls,
    config$experimental_design$blank_controls
  )
  
  # ── 2. Build a bulletproof shield for control samples ───────────────────────
  # Protects by exact YAML match OR by keyword pattern in the sample label/condition
  protected_labels <- exp_design$label[
    exp_design$condition %in% configured_controls |
      grepl("Background|IgG|Bead|Blank|Control|Buffer", exp_design$label, ignore.case = TRUE) |
      grepl("Background|IgG|Bead|Blank|Control|Buffer", exp_design$condition, ignore.case = TRUE)
  ]
  
  # ── 3. Drop only true biological samples that fail the quality cutoff ───────
  low_quality_samples <- names(protein_counts[
    protein_counts < cutoff & 
      !names(protein_counts) %in% protected_labels
  ])
  
  if (length(low_quality_samples) > 0) {
    message(paste(">>> Dropping poor biological samples (<", round(cutoff), "proteins):",
                  paste(low_quality_samples, collapse = ", ")))
    
    excluded_samples_df <- data.frame(
      Sample             = low_quality_samples,
      Proteins_Detected  = protein_counts[low_quality_samples],
      Required_Threshold = round(cutoff),
      Reason             = paste("Failed", sample_drop_threshold * 100, "% quality threshold")
    )
    readr::write_csv(excluded_samples_df, file.path(logs_dir, "Excluded_Samples.csv"))
    
    exp_design <- exp_design %>% dplyr::filter(!label %in% low_quality_samples)
    lfq_cols   <- match(exp_design$label, colnames(df))
  } else {
    message(paste(">>> All biological samples passed quality threshold (>= ", round(cutoff), " proteins). Control samples shielded."))
  }
  
  # Add after the low_quality_samples removal block in section 3:
  if (has_clones_in_design) {
    # Check for any clones that have lost all replicates after sample dropping
    surviving_clones <- exp_design %>%
      dplyr::group_by(condition, clone_id) %>%
      dplyr::summarise(n_remaining = dplyr::n(), .groups = "drop")
    
    lost_clones <- surviving_clones %>% dplyr::filter(n_remaining == 0)
    sole_rep_clones <- surviving_clones %>% dplyr::filter(n_remaining == 1)
    
    if (nrow(sole_rep_clones) > 0) {
      message("   -> Warning: The following clones now have only 1 replicate after sample QC:")
      print(sole_rep_clones)
    }
  }
  
  # Snapshot the pristine data before any proteins are deleted
  df_full <- df

  # ===========================================================================#
  # 4. EXCLUDE CONTAMINANTS, BLANKS & MATRIX BACKGROUND ####
  # ===========================================================================#
  filt_res    <- filter_contaminants(df, crap_ids, crapome_list, bioid_list, is_bioid)
  df_filtered <- filt_res$filtered
  contam_log  <- filt_res$log
  
  lfq_col_names <- colnames(df)[lfq_cols]
  
  if (nrow(contam_log) > 0) {
    dropped_contaminants <- df %>%
      dplyr::inner_join(contam_log %>% dplyr::select(ID, Reason), by = "ID") %>%
      dplyr::select(ID, name, Full_Description = name_raw, Reason,
                    dplyr::all_of(lfq_col_names)) %>%
      dplyr::mutate(Samples_Detected = rowSums(!is.na(
        dplyr::select(., dplyr::all_of(lfq_col_names))
      ))) %>%
      dplyr::arrange(dplyr::desc(Samples_Detected))
    
    readr::write_csv(dropped_contaminants,
                     file.path(logs_dir, "Excluded_Contaminants.csv"))
  }
  
  df_unique <- DEP::make_unique(df_filtered, "name", "ID", delim = ";")
  
  # ── TIER 1: Absolute Blank Exclusion (100% Presence in System Blanks) ──────
  blank_conds  <- config$experimental_design$blank_controls
  blank_design <- exp_design %>% dplyr::filter(condition %in% blank_conds)
  
  if (nrow(blank_design) > 0) {
    message(paste(">>> Tier 1: Screening system blanks for 100% ambient carryover:",
                  paste(unique(blank_design$condition), collapse = ", ")))
    
    blank_cols <- match(blank_design$label, colnames(df_unique))
    
    # Calculate how many blank samples each protein appears in
    blank_detections <- rowSums(!is.na(df_unique[, blank_cols, drop = FALSE]))
    
    # DROP RULE: Present in 100% of blank tubes (e.g., 3 out of 3)
    is_system_contam <- blank_detections == length(blank_cols)
    system_contams   <- df_unique$ID[is_system_contam]
    
    if (length(system_contams) > 0) {
      message(paste("   -> Dropping", length(system_contams),
                    "proteins detected in 100% of system blank controls."))
      
      blank_col_names <- colnames(df_unique)[blank_cols]
      
      # FIX: Hardcode length(blank_cols) instead of passing the 450-element vector
      dropped_blanks <- df_unique %>%
        dplyr::filter(ID %in% system_contams) %>%
        dplyr::select(ID, name, Full_Description = name_raw, dplyr::all_of(blank_col_names)) %>%
        dplyr::mutate(
          Reason = paste0("Tier 1 Blank Exclusion (Present in ", length(blank_cols), "/", length(blank_cols), " blank tubes)")
        )
      
      readr::write_csv(dropped_blanks, file.path(logs_dir, "Excluded_Tier1_SystemBlanks.csv"))
      df_unique <- df_unique %>% dplyr::filter(!ID %in% system_contams)
    }
    
    # Remove system blanks from exp_design so they do not enter downstream normalization/Limma
    exp_design <- exp_design %>% dplyr::filter(!condition %in% blank_conds)
  }
  
  # ── TIER 2: Max-of-Means Quantitative Gate (vs Matrix Controls) ────────────
  bg_conds  <- config$experimental_design$background_controls
  bg_design <- exp_design %>% dplyr::filter(condition %in% bg_conds)
  
  if (nrow(bg_design) > 0) {
    message(paste(">>> Tier 2: Running Max-of-Means quantitative gate against:",
                  paste(unique(bg_design$condition), collapse = ", ")))
    
    bg_cols  <- match(bg_design$label, colnames(df_unique))
    exp_cols <- match(exp_design$label[!exp_design$condition %in% bg_conds], colnames(df_unique))
    
    # 1. Calculate Linear Mean Abundance in Matrix Controls
    bg_matrix <- as.matrix(df_unique[, bg_cols, drop = FALSE])
    bg_means  <- rowMeans(bg_matrix, na.rm = TRUE)
    bg_means[is.na(bg_means) | is.nan(bg_means)] <- 0 
    
    # 2. Calculate MAXIMUM Linear Mean Abundance across any experimental condition
    exp_matrix <- as.matrix(df_unique[, exp_cols, drop = FALSE])
    exp_conds  <- exp_design$condition[!exp_design$condition %in% bg_conds]
    
    max_exp_means <- apply(exp_matrix, 1, function(row) {
      cond_means <- tapply(as.numeric(row), exp_conds, mean, na.rm = TRUE)
      max_val    <- max(cond_means, na.rm = TRUE)
      return(ifelse(is.na(max_val) | is.infinite(max_val), 0, max_val))
    })
    
    # 3. CRITICAL FIX: Compute true Log2 Fold Enrichment over background safely
    enrichment_fc <- dplyr::case_when(
      bg_means == 0 & max_exp_means > 0 ~ Inf,  # Clean specific binder (no background signal)
      bg_means == 0 & max_exp_means == 0 ~ 0,    # Absent in both experimental and background tubes
      max_exp_means == 0                 ~ -Inf, # Present only in background tubes
      TRUE ~ log2(max_exp_means / bg_means)      # True Log2 Fold Change
    )
    
    # 4. Apply Quantitative Threshold
    min_enrichment_thr <- config$pipeline_parameters$min_background_enrichment %||% 1.0
    quantitative_fails <- enrichment_fc < min_enrichment_thr
    
    matrix_contams <- df_unique$ID[quantitative_fails]
    
    if (length(matrix_contams) > 0) {
      message(paste("   -> Dropping", length(matrix_contams),
                    "proteins failing Tier 2 Max-of-Means enrichment (<", 
                    min_enrichment_thr, "log2 FC over matrix controls)."))
      
      bg_col_names  <- colnames(df_unique)[bg_cols]
      exp_col_names <- colnames(df_unique)[exp_cols]
      
      dropped_matrix <- df_unique %>%
        dplyr::filter(ID %in% matrix_contams) %>%
        dplyr::select(ID, name, Full_Description = name_raw,
                      dplyr::all_of(bg_col_names), dplyr::all_of(exp_col_names)) %>%
        dplyr::mutate(
          Contaminant_Type = ifelse(bg_means[quantitative_fails] == 0, 
                                    "Sparse Ghost Dropout (Absent in Exp & Bg)", 
                                    "True Matrix Binder (< 1.0 log2FC over Bg)"),
          Reason        = paste0("Tier 2 Matrix Gate (< ", min_enrichment_thr, " log2FC over IgG/Beads)"),
          Bg_Mean_Log2  = round(ifelse(bg_means[quantitative_fails] == 0, 0, log2(bg_means[quantitative_fails])), 2),
          Max_Exp_Log2  = round(ifelse(max_exp_means[quantitative_fails] == 0, 0, log2(max_exp_means[quantitative_fails])), 2),
          Enrichment_FC = round(ifelse(is.infinite(enrichment_fc[quantitative_fails]), 99.99, enrichment_fc[quantitative_fails]), 2)
        ) %>%
        dplyr::arrange(Enrichment_FC)
      
      readr::write_csv(dropped_matrix, file.path(logs_dir, "Excluded_Tier2_BackgroundEnrichment.csv"))
    }
    
    df_unique <- df_unique %>% dplyr::filter(!ID %in% matrix_contams)
  }
  
  # Re-align column indices to match the cleaned design matrix
  lfq_cols <- match(exp_design$label, colnames(df_unique))
  
  # ── Create a visualisation-only design that excludes background controls ──────
  # Background controls are retained in the full exp_design for limma contrasts
  # but should not appear in QC plots, PCA, heatmaps, or correlation matrices
  vis_design <- exp_design %>%
    dplyr::filter(!condition %in% config$experimental_design$background_controls)
  
  vis_labels <- vis_design$label
  
  # ===========================================================================#
  # 5. DEP log2 OBJECT CREATION & CLONE-AWARE MISSING VALUE FILTER ####
  # ===========================================================================#
  safe_labels         <- paste0("sample_", seq_along(lfq_cols))
  colnames(df_unique)[lfq_cols] <- safe_labels
  exp_design$label    <- safe_labels
  
  #Log2 transformed values
  se <- DEP::make_se(df_unique, lfq_cols, exp_design)
  
  pretty_names <- if ("clone_id" %in% colnames(exp_design)) {
    paste(exp_design$clone_id, exp_design$replicate, sep = "_")
  } else {
    paste(exp_design$condition, exp_design$replicate, sep = "_")
  }
  
  colnames(se)                           <- pretty_names
  SummarizedExperiment::colData(se)$name <- pretty_names
  SummarizedExperiment::colData(se)$ID   <- pretty_names
  exp_design$label                       <- pretty_names
  
  save_and_print_plot(DEP::plot_numbers(se),
                      file.path(qc_dir, "QC1_Protein_Numbers"))
  
  # Fix 4: Clone-aware proportional missing value filter
  min_valid_prop <- config$pipeline_parameters$min_valid_proportion %||% 0.50
  min_clone_prop <- config$pipeline_parameters$min_clone_proportion %||% 0.50
  
  is_na_mat  <- is.na(SummarizedExperiment::assay(se))
  col_data   <- as.data.frame(SummarizedExperiment::colData(se))
  has_clones <- "clone_id" %in% colnames(col_data)
  
  if (has_clones) {
    message("   -> Applying clone-aware proportional missing value filter...")
    message(paste0("      Requires >= ", min_valid_prop * 100,
                   "% valid values per clone in >= ",
                   min_clone_prop * 100,
                   "% of clones, in at least one genotype"))
    
    # Pre-compute design structure once for efficiency
    design_structure <- lapply(
      split(seq_len(nrow(col_data)), col_data$condition),
      function(geno_idx) split(geno_idx, col_data$clone_id[geno_idx])
    )
    
    # Log the detected structure
    message("   -> Detected design structure:")
    for (geno in names(design_structure)) {
      clone_sizes <- sapply(design_structure[[geno]], length)
      message(paste0("      ", geno, ": ", length(clone_sizes), " clone(s) — ",
                     paste(names(clone_sizes), "(n=", clone_sizes, ")", collapse = ", ")))
    }
    
    keep_rows <- apply(!is_na_mat, 1, function(present) {
      any(sapply(design_structure, function(clones_in_geno) {
        n_clones_total   <- length(clones_in_geno)
        n_clones_passing <- sum(sapply(clones_in_geno, function(clone_idx) {
          n_present <- sum(present[clone_idx])
          required  <- ceiling(length(clone_idx) * min_valid_prop)
          n_present >= required
        }))
        required_clones <- ceiling(n_clones_total * min_clone_prop)
        n_clones_passing >= required_clones
      }))
    })
    
    # Write threshold log for methods reporting and QC report
    threshold_log <- lapply(names(design_structure), function(geno) {
      lapply(names(design_structure[[geno]]), function(clone) {
        n_reps   <- length(design_structure[[geno]][[clone]])
        required <- ceiling(n_reps * min_valid_prop)
        data.frame(
          Genotype           = geno,
          Clone              = clone,
          N_Replicates       = n_reps,
          Required_Valid     = required,
          Proportion_Applied = min_valid_prop
        )
      }) %>% dplyr::bind_rows()
    }) %>%
      dplyr::bind_rows() %>%
      dplyr::mutate(
        N_Clones_In_Genotype = sapply(Genotype, function(g) length(design_structure[[g]])),
        Required_Clones      = ceiling(N_Clones_In_Genotype * min_clone_prop)
      )
    
    message("   -> Resolved thresholds per clone:")
    print(threshold_log)
    readr::write_csv(threshold_log,
                     file.path(logs_dir, "Filter_Thresholds_MissingValues.csv"))
    
  } else {
    message("   -> No clone_id found. Applying condition-level proportional filter...")
    message(paste0("      Requires >= ", min_valid_prop * 100,
                   "% valid values in at least one condition"))
    
    keep_rows <- apply(!is_na_mat, 1, function(present) {
      any(tapply(present, col_data$condition, function(x) {
        required <- ceiling(length(x) * min_valid_prop)
        sum(x) >= required
      }))
    })
  }
  
  message(paste0("   -> Retained ", sum(keep_rows),  " / ", length(keep_rows),
                 " proteins after missing value filter"))
  message(paste0("   -> Removed  ", sum(!keep_rows), " proteins"))
  
  # Log dropped proteins
  dropped_by_missingness <- SummarizedExperiment::rowData(se)$name[!keep_rows]
  if (length(dropped_by_missingness) > 0) {
    dropped_mv_df <- df_unique %>%
      dplyr::filter(name %in% dropped_by_missingness) %>%
      dplyr::select(ID, name, Full_Description = name_raw,
                    dplyr::all_of(safe_labels)) %>%
      dplyr::mutate(
        Reason = paste0("Failed clone-aware missing value filter (",
                        min_valid_prop * 100, "% per clone, ",
                        min_clone_prop * 100, "% of clones)"),
        Total_Detections = rowSums(!is.na(dplyr::select(., dplyr::all_of(safe_labels))))
      ) %>%
      dplyr::arrange(dplyr::desc(Total_Detections))
    
    readr::write_csv(dropped_mv_df,
                     file.path(logs_dir, "Excluded_Proteins_MissingValues.csv"))
  }
  
  se_filtered <- se[keep_rows, ]
  save_and_print_plot(DEP::plot_missval(se_filtered),
                      file.path(qc_dir, "QC2_Missing_Values"))
  
  # ===========================================================================#
  # 6. NORMALISATION & IMPUTATION ####
  # ===========================================================================#
  
  # Fix 1: Raw intensity boxplot appears once only, here in section 6
  raw_long <- as.data.frame(SummarizedExperiment::assay(se_filtered)) %>%
    tibble::rownames_to_column("Protein") %>%
    tidyr::pivot_longer(-Protein, names_to = "Sample", values_to = "Intensity") %>%
    dplyr::filter(!is.na(Intensity)) %>%
    dplyr::left_join(
      as.data.frame(SummarizedExperiment::colData(se_filtered)) %>%
        tibble::rownames_to_column("Sample") %>%
        dplyr::select(Sample, condition),
      by = "Sample"
    )
  
  raw_long <- raw_long %>% 
    dplyr::filter(!condition %in% config$experimental_design$background_controls)
  
  p_raw_box <- ggplot2::ggplot(
    raw_long,
    ggplot2::aes(x = Sample, y = Intensity, fill = condition)
  ) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(
      title = paste(exp_name, "- Raw Intensity Distribution (Pre-Normalisation)"),
      y     = "Raw Intensity",
      x     = ""
    )
  save_and_print_plot(p_raw_box,
                      file.path(qc_dir, "QC0_Raw_Intensity_Distribution"))
  
  se_norm  <- DEP::normalize_vsn(se_filtered)
  save_and_print_plot(DEP::plot_normalization(se_filtered, se_norm),
                      file.path(qc_dir, "QC3_VSN_Density"))
  
  message("   -> Generating VSN Mean-Variance Plot...")
  norm_mat     <- SummarizedExperiment::assay(se_norm)
  vsn_plot_obj <- vsn::meanSdPlot(norm_mat, plot = FALSE)
  p_vsn        <- vsn_plot_obj$gg +
    ggplot2::theme_bw() +
    ggplot2::ggtitle(paste(exp_name, "- VSN Mean-Variance Fit"))
  
  save_and_print_plot(p_vsn, file.path(qc_dir, "QC4_VSN_MeanSdPlot"), width = 7, height = 7)
  
  analysis_mode <- config$experiments[[exp_name]]$analysis_mode
  
  if (!is.null(analysis_mode) && analysis_mode == "abundance_only") {
    
    message(">>> Interactome mode: skipping imputation and differential expression <<<")
    message(">>> Applying clone/genotype presence thresholds to normalized data <<<")
    
    norm_mat <- SummarizedExperiment::assay(se_norm)
    
    scrub_res <- scrub_by_clone_genotype_presence(
      mat            = norm_mat,
      exp_design     = exp_design,
      min_valid_prop = config$pipeline_parameters$min_valid_proportion %||% 0.5,
      min_clone_prop = config$pipeline_parameters$min_clone_proportion %||% 0.5
    )
    
    # Rownames in the matrix are 'name', so we map that to 'Protein' for the join
    id_map <- df_unique %>% dplyr::select(Protein = name, ID) %>% dplyr::distinct()
    
    report <- export_interactome_abundance_report(
      mat_scrubbed = scrub_res$mat_scrubbed,
      exp_design   = exp_design,
      na_log       = scrub_res$na_log,
      output_dir   = de_dir,
      exp_name     = exp_name,
      id_map       = id_map
    )
    
    # --- NEW: Reorder contrasts so the blacklist is processed first ---
    bl_contrast <- config$comparative_analysis$blacklist_contrast
    ordered_contrasts <- contrasts
    if (!is.null(bl_contrast) && bl_contrast %in% contrasts) {
      ordered_contrasts <- c(bl_contrast, setdiff(contrasts, bl_contrast))
    }
    
    message("   -> Calculating non-imputed raw statistics across contrasts...")
    raw_stats_list <- list()
    ko_binders     <- c()
    
    for (contrast in ordered_contrasts) {
      tryCatch({
        # If this is the blacklist contrast, ko_binders is empty (plots everything).
        # If this is a biological contrast, ko_binders masks the background from the plot.
        res_df <- calculate_raw_abundance_stats(
          mat_scrubbed       = scrub_res$mat_scrubbed,
          exp_design         = exp_design,
          contrast_str       = contrast,
          output_base_dir    = output_base_dir,
          exp_name           = exp_name,
          id_map             = id_map,
          poi_list           = poi_list,
          blacklist_proteins = ko_binders 
        )
        raw_stats_list[[contrast]] <- res_df
        
        # If we just processed the blacklist contrast, extract the hits to mask subsequent plots
        if (!is.null(bl_contrast) && contrast == bl_contrast) {
          target_cond <- stringr::str_trim(stringr::str_split(contrast, "-")[[1]][1])
          ko_binders <- res_df %>%
            dplyr::filter(
              (Status == "Quantitative" & P_Value < config$thresholds$fdr & logFC > config$thresholds$lfc) |
                Status == paste("Exclusive to", target_cond)
            ) %>%
            dplyr::pull(Protein)
          message(paste("      -> Extracted", length(ko_binders), "background binders from", bl_contrast, "to mask in true comparisons."))
        }
      }, error = function(e) {
        message("      -> Skipping stats for ", contrast, ": ", e$message)
      })
    }
    
    # Visualisations (Leveraging save_and_print_plot for R viewing and TIFF generation)
    if (length(poi_list) > 0) {
      p_strip <- plot_protein_abundance_strip(
        long_df  = report$long, 
        proteins = poi_list, 
        name_col = "Protein", 
        exp_name = exp_name
      )
      save_and_print_plot(p_strip, file.path(de_dir, paste0(exp_name, "_Abundance_StripPlot_POI")))
    }
    
    p_heat <- plot_abundance_heatmap(report$clone_means, exp_name = exp_name)
    save_and_print_plot(p_heat, file.path(de_dir, paste0(exp_name, "_Abundance_Heatmap_AllProteins")), width = 10, height = 14)
    
    # --- RESCUED QC PLOTS FOR ABUNDANCE MODE ---
    message("   -> Generating QC Plots (MDS, CVs, Rank, Correlation) for Raw Data...")
    
    # Establish visualisation bounds (excluding background controls like IgG)
    vis_design <- exp_design %>% dplyr::filter(!condition %in% config$experimental_design$background_controls)
    vis_labels <- vis_design$label
    vis_matrix <- norm_mat[, vis_labels, drop = FALSE]
    
    # 1. MDS Plot (Handles NAs)
    generate_mds_plot(
      mat        = vis_matrix,
      design     = vis_design,
      title      = exp_name,
      output_dir = qc_dir
    )
    
    # 2. Replicate CVs (Handles NAs natively via na.rm = TRUE)
    cv_data <- data.frame(Protein = rownames(vis_matrix))
    for (cond in unique(vis_design$condition)) {
      cols <- vis_design$label[vis_design$condition == cond]
      if (length(cols) > 1) {
        means <- rowMeans(vis_matrix[, cols, drop = FALSE], na.rm = TRUE)
        sds   <- apply(vis_matrix[, cols, drop = FALSE], 1, sd, na.rm = TRUE)
        cv_data[[cond]] <- (sds / means) * 100
      }
    }
    cv_melt <- tidyr::pivot_longer(cv_data, cols = -Protein, names_to = "Condition", values_to = "CV") %>%
      dplyr::filter(!is.na(CV) & !is.infinite(CV))
    
    p_cv <- ggplot2::ggplot(cv_melt, ggplot2::aes(x = Condition, y = CV, fill = Condition)) +
      ggplot2::geom_violin(alpha = 0.7, trim = FALSE) +
      ggplot2::geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
      ggplot2::geom_hline(yintercept = 20, linetype = "dashed", color = "red") +
      ggplot2::theme_bw() +
      ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
      ggplot2::labs(title = paste(exp_name, "- Replicate CV (%)"), y = "Coefficient of Variation (%)", x = "")
    save_and_print_plot(p_cv, file.path(qc_dir, "QC9_Replicate_CVs"), width = 8, height = 6)
    
    # 3. Protein Rank Plot (Split by Genotype)
    add_rank_plot(
      imputed_matrix = vis_matrix,
      exp_design     = vis_design,
      bait           = config$experimental_design$bait_protein,
      poi_list       = poi_list,
      output_path    = file.path(qc_dir, "QC10_Protein_Rank_Plot"),
      exp_name       = exp_name
    )
    
    # 4. Pre-Imputation Correlation Matrix
    cor_mat <- cor(vis_matrix, method = "spearman", use = "pairwise.complete.obs")
    anno_df <- data.frame(Condition = vis_design$condition, row.names = colnames(vis_matrix))
    p_cor <- pheatmap::pheatmap(
      cor_mat,
      main            = paste(exp_name, "- Pre-Imputation Spearman Correlation"),
      display_numbers = round(cor_mat, 2),
      number_format   = "%.2f",
      fontsize_number = 7,
      color           = colorRampPalette(c("#2166AC", "#F7F7F7", "#D6604D"))(100),
      breaks          = seq(0.7, 1.0, length.out = 101),
      annotation_col  = anno_df,
      clustering_method = "ward.D2",
      angle_col       = 45,
      border_color    = "white",
      silent          = TRUE
    )
    save_and_print_plot(p_cor, file.path(qc_dir, "QC8_Correlation_Matrix"), width = 9, height = 8)

    # Bypass Limma DE and Imputation matrices so the checkpoint doesn't crash
    se_imputed <- se_norm 
    limma_res  <- NULL
    
  } else {
    
  # Read MNAR threshold from config
  mnar_threshold <- config$pipeline_parameters$imputation$mnar_threshold
  
  if (is.null(imputation_method) || length(imputation_method) == 0) {
    imputation_method <- "mixed"
    message("   -> Warning: No imputation method specified. Defaulting to 'mixed'.")
  }
  
  if (imputation_method == "mixed") {
    message(">>> Applying Mixed Imputation (Condition/Clone-aware MNAR + KNN) <<<")
    se_imputed <- impute_mixed_mar_mnar(se_norm, min_missing_pct = mnar_threshold)
  } else if (imputation_method == "knn") {
    message(">>> Applying Standard KNN Imputation <<<")
    se_imputed <- DEP::impute(se_norm, fun = "knn", rowmax = 0.9)
  } else {
    message(">>> Applying MinProb (Left-Shift) Imputation <<<")
    se_imputed <- DEP::impute(se_norm, fun = "MinProb",
                              q = config$pipeline_parameters$imputation$minprob_q)
  }
  
  message("   -> Generating Per-Sample Imputation Plot...")
  p_imp <- plot_imputation_per_sample(
    se_norm, se_imputed,
    paste(exp_name, "- Imputation Shift per Sample")
  )
  save_and_print_plot(p_imp,
                      file.path(qc_dir, "QC5_Imputation_Shift_PerSample"),
                      width = 12, height = 10)
  
  imputed_matrix <- SummarizedExperiment::assay(se_imputed)
  
  vis_design <- exp_design %>%
    dplyr::filter(!condition %in% config$experimental_design$background_controls)
  vis_labels <- vis_design$label
  
  vis_matrix <- imputed_matrix[, vis_labels, drop = FALSE]  # used after imputation
  
  export_imputation_map_excel(
    se_norm         = se_norm,
    imputed_matrix  = imputed_matrix,
    exp_design      = exp_design,
    output_file     = file.path(qc_dir, paste0("QC5b_Imputation_ColorMap_",
                                               exp_name, ".xlsx")),
    min_missing_pct = mnar_threshold
  )
  
  # ===========================================================================#
  # 6b. SCOPE INTERCEPT (Late-Filtering for PhD Subsets)
  # ===========================================================================#
  # VSN and Imputation used the full batch degrees of freedom. Now, slice the 
  # metadata and matrices down to ONLY the active scope before DE and diagnostics!
  if (!is.null(active_design)) {
    target_conditions <- unique(active_design$condition)
    
    message(paste("   -> Applying Scope Intercept: Slicing from", 
                  length(unique(exp_design$condition)), "total batch conditions down to", 
                  length(target_conditions), "active scope conditions:", 
                  paste(target_conditions, collapse = ", ")))
    
    # Filter metadata (which already has clean replicate numbers and pretty_names!)
    exp_design <- exp_design %>% dplyr::filter(condition %in% target_conditions)
    
    # Slice the SummarizedExperiment and Matrix objects cleanly by sample label
    se_norm    <- se_norm[, exp_design$label]
    se_imputed <- se_imputed[, exp_design$label]
    imputed_matrix <- SummarizedExperiment::assay(se_imputed)
    norm_mat       <- SummarizedExperiment::assay(se_norm)
  }
  
  # ============================================================================#
  # 7. HEATMAP, PCA, RANK PLOT, CORRELATION & CVs ####
  # ============================================================================#

  message("   -> Generating Heatmap for Top Variable Proteins...")
  vars            <- apply(vis_matrix, 1, var)
  top_var_indices <- order(vars, decreasing = TRUE)[1:min(100, nrow(vis_matrix))]
  heatmap_mat_top <- t(scale(t(vis_matrix[top_var_indices, ])))
  
  p_heat <- pheatmap::pheatmap(
    heatmap_mat_top,
    main              = paste(exp_name, "- Top 100 Variable Proteins (Z-score)"),
    clustering_method = "ward.D2",
    show_colnames     = TRUE,
    show_rownames     = (nrow(heatmap_mat_top) <= 100),
    fontsize_row      = 6,             # <-- Shrink font to fit 100 rows cleanly
    fontsize_col      = 8,             # <-- Optional: adjust column header size
    color             = colorRampPalette(c("navy", "white", "firebrick3"))(100)
  )
  
  save_and_print_plot(p_heat, file.path(qc_dir, "QC6_Top_Variable_Proteins"), width = 10, height = 16)
  
  # Save the exact names of the top 100 variable proteins
  readr::write_csv(
    data.frame(Rank = 1:length(top_var_indices), Protein = rownames(heatmap_mat_top)),
    file.path(qc_dir, "QC6_Top_Variable_Proteins_List.csv")
  )
  
  generate_pca_full(
    mat        = vis_matrix,
    design     = vis_design,
    title      = paste(exp_name, "PCA"),
    n_loadings = 20,
    output_dir = qc_dir
  )
  
  add_rank_plot(
    imputed_matrix = vis_matrix,
    exp_design     = vis_design,
    bait           = config$experimental_design$bait_protein,
    poi_list       = poi_list,
    output_path    = file.path(qc_dir, "QC10_Protein_Rank_Plot"),
    exp_name       = exp_name
  )
  
  generate_correlation_matrices(
    se_norm    = se_norm[,    vis_labels],   # subset to non-background samples
    se_imputed = se_imputed[, vis_labels],
    exp_design = vis_design,
    output_dir = qc_dir,
    exp_name   = exp_name
  )
  
  # Data matrix exports
  message("   -> Exporting formatted data matrices...")
  readr::write_csv(
    as.data.frame(SummarizedExperiment::assay(se_filtered)) %>%
      tibble::rownames_to_column("Protein"),
    file.path(qc_dir, "Matrix_1_Filtered_Raw.csv")
  )
  readr::write_csv(
    as.data.frame(norm_mat) %>% tibble::rownames_to_column("Protein"),
    file.path(qc_dir, "Matrix_2_VSN_Normalized.csv")
  )
  readr::write_csv(
    as.data.frame(imputed_matrix) %>% tibble::rownames_to_column("Protein"),
    file.path(qc_dir, "Matrix_3_Fully_Imputed.csv")
  )
  
  # Replicate CVs
  message("   -> Calculating Replicate CVs...")
  cv_data <- data.frame(Protein = rownames(norm_mat))
  for (cond in unique(vis_design$condition)) {
    cols <- vis_design$label[vis_design$condition == cond]
    if (length(cols) > 1) {
      means              <- rowMeans(norm_mat[, cols, drop = FALSE], na.rm = TRUE)
      sds                <- apply(norm_mat[, cols, drop = FALSE], 1, sd, na.rm = TRUE)
      cv_data[[cond]]    <- (sds / means) * 100
    }
  }
  
  cv_melt <- tidyr::pivot_longer(cv_data, cols = -Protein,
                                 names_to = "Condition", values_to = "CV") %>%
    dplyr::filter(!is.na(CV) & !is.infinite(CV))
  
  p_cv <- ggplot2::ggplot(cv_melt,
                          ggplot2::aes(x = Condition, y = CV, fill = Condition)) +
    ggplot2::geom_violin(alpha = 0.7, trim = FALSE) +
    ggplot2::geom_boxplot(width = 0.2, fill = "white", outlier.shape = NA) +
    ggplot2::geom_hline(yintercept = 20, linetype = "dashed", color = "red") +
    ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = paste(exp_name, "- Replicate CV (%)"),
                  y = "Coefficient of Variation (%)", x = "")
  save_and_print_plot(p_cv, file.path(qc_dir, "QC9_Replicate_CVs"),
                      width = 8, height = 6)
  
  # ============================================================================#
  # 8. DIFFERENTIAL EXPRESSION ####
  # ============================================================================#
  limma_res <- perform_limma_analysis(
    imputed_matrix, exp_design, contrasts, exp_design$clone_id, "condition"
  )
  
  bl_contrast <- config$comparative_analysis$blacklist_contrast
  black_safe  <- if (!is.null(bl_contrast)) make_safe_contrast_name(bl_contrast) else NULL
  ko_binders  <- c()
  
  # --- 1. Identify background binders FIRST ---
  if (!is.null(black_safe) && black_safe %in% names(limma_res)) {
    ko_binders_df <- limma_res[[black_safe]] %>% dplyr::filter(adj.P.Val < fdr_thr & logFC > 0)
    ko_binders    <- ko_binders_df$Protein  # FIX
    
    if (length(ko_binders) > 0) {
      message(paste("   -> Extracted", length(ko_binders), "background binders from", bl_contrast, "to mask in true comparisons."))
      
      # Write the exclusion log directly here
      id_map <- df_unique %>% dplyr::select(ID, Protein = name) %>% dplyr::distinct() # Map cleanly
      blacklist_log <- ko_binders_df %>%
        dplyr::left_join(id_map, by = "Protein") %>%
        dplyr::select(ID, Protein, logFC, adj.P.Val) %>%
        dplyr::mutate(Reason = "Statistical Blacklist (IP-MS Background)") %>%
        dplyr::arrange(dplyr::desc(logFC))
      readr::write_csv(blacklist_log, file.path(logs_dir, "Excluded_Statistical_Blacklist.csv"))
    }
  }
  
  # --- 2. Loop through and Plot ---
  for (contrast_name in names(limma_res)) {
    safe_name <- make_safe_contrast_name(contrast_name)
    plot_df   <- limma_res[[contrast_name]]
    
    # Save the FULL raw stats to CSV (unmasked)
    readr::write_csv(plot_df, file.path(de_dir, paste0(safe_name, "_Results.csv")))
    
    # ---> NEW: Mask the KO binders from the visualization dataframe <---
    if (!is.null(black_safe) && safe_name != black_safe && length(ko_binders) > 0) {
      plot_df <- plot_df %>% dplyr::filter(!Protein %in% ko_binders)
    }
    
    n_proteins    <- nrow(plot_df)
    use_nominal_p <- n_proteins < 200
    if (use_nominal_p) message(paste0("   -> Warning: Only ", n_proteins, " proteins in '", contrast_name, "'. Switching volcano to nominal p-value."))
    
    p_volc <- generate_universal_volcano_plot(
      res_df             = plot_df,
      contrast_str       = contrast_name,
      exp_name           = exp_name,
      fdr_thr            = fdr_thr,
      lfc_thr            = lfc_thr,
      poi_list           = poi_list,
      blacklist_proteins = ko_binders,
      use_nominal_p      = use_nominal_p
    )
    
    # Print to R view and save to nested TIFF folder
    save_and_print_plot(p_volc, file.path(de_dir, paste0("Volcano_", safe_name)), width = 10, height = 8)
    
    p_ma <- generate_ma_plot(
      de_result = plot_df,
      title     = paste(exp_name, "-", contrast_name),
      fdr_thr   = fdr_thr,
      lfc_thr   = lfc_thr
    )
    save_and_print_plot(p_ma, file.path(de_dir, paste0("MA_", safe_name)))
  }
  
  # Fix 3: Blacklist contrast name uses make_safe_contrast_name()
  if (!is.null(config$comparative_analysis$blacklist_contrast)) {
    black_safe <- make_safe_contrast_name(config$comparative_analysis$blacklist_contrast)
    
    if (black_safe %in% names(limma_res)) {
      ko_binders_df <- limma_res[[black_safe]] %>%
        dplyr::filter(adj.P.Val < fdr_thr & logFC > 0)
      
      if (nrow(ko_binders_df) > 0) {
        id_map <- df_unique %>%
          dplyr::select(ID, Protein) %>%
          dplyr::distinct()
        
        blacklist_log <- ko_binders_df %>%
          dplyr::left_join(id_map, by = "name") %>%
          dplyr::select(ID, name, logFC, adj.P.Val) %>%
          dplyr::mutate(Reason = "Statistical Blacklist (IP-MS Background)") %>%
          dplyr::arrange(dplyr::desc(logFC))
        
        readr::write_csv(blacklist_log,
                         file.path(logs_dir, "Excluded_Statistical_Blacklist.csv"))
      }
    }
  }
  
  # Combined DE Excel workbook
  wb_de <- openxlsx::createWorkbook()
  for (contrast_name in names(limma_res)) {
    safe_name  <- make_safe_contrast_name(contrast_name)
    sheet_name <- substr(safe_name, 1, 31)
    
    openxlsx::addWorksheet(wb_de, sheet_name)
    
    de_sheet <- limma_res[[contrast_name]] %>%
      dplyr::arrange(adj.P.Val) %>%
      dplyr::mutate(dplyr::across(where(is.numeric), ~ round(., 4)))
    
    openxlsx::writeDataTable(wb_de, sheet_name, de_sheet,
                             tableStyle = "TableStyleMedium9")
    
    sig_rows <- which(
      de_sheet$adj.P.Val < fdr_thr & abs(de_sheet$logFC) > lfc_thr
    ) + 1
    
    if (length(sig_rows) > 0) {
      openxlsx::addStyle(
        wb_de, sheet_name,
        openxlsx::createStyle(fgFill = "#FFF2CC"),
        rows = sig_rows, cols = 1:ncol(de_sheet), gridExpand = TRUE
      )
    }
  }
  
  openxlsx::saveWorkbook(
    wb_de,
    file.path(de_dir, paste0(exp_name, "_All_DE_Results.xlsx")),
    overwrite = TRUE
  )
  
  }
  
  # ============================================================================
  # 9. SAVE PIPELINE CHECKPOINT
  # ============================================================================
  message("   -> Saving environment checkpoint for downstream scripts...")
  
  if (!is.null(analysis_mode) && analysis_mode == "abundance_only") {
    # Abundance Mode: Pass our custom dual-track frames as de_results
    checkpoint_de_res <- raw_stats_list
    names(checkpoint_de_res) <- make_safe_contrast_name(names(raw_stats_list))
    safe_names <- names(checkpoint_de_res)
  } else {
    # Standard Mode: Run the Limma safeguard
    if (is.null(names(limma_res)) || any(names(limma_res) == "")) {
      stop(paste(
        "CRITICAL: limma results have missing names for experiment:", exp_name,
        "\nThis would corrupt all downstream scripts. Check perform_limma_analysis()."
      ))
    }
    checkpoint_de_res <- limma_res
    names(checkpoint_de_res) <- make_safe_contrast_name(names(limma_res))
    safe_names <- names(checkpoint_de_res)
  }
  
  final_imputed_mat <- SummarizedExperiment::assay(se_imputed)
  final_norm_mat    <- SummarizedExperiment::assay(se_norm)
  
  checkpoint <- list(
    exp_design_clean      = exp_design,
    raw_df_aligned        = df_full,
    imputed_matrix        = final_imputed_mat, 
    norm_matrix           = final_norm_mat,
    de_results            = checkpoint_de_res, # Universally formatted list
    de_results_safe_names = safe_names
  )
  
  saveRDS(checkpoint,
          file.path(out_dir, paste0(exp_name, "_checkpoint.rds")))
  
  message(paste(">>> Pipeline Complete for:", exp_name, "<<<"))
  
  # Generate the HTML QC report
  generate_qc_report(exp_name, output_base_dir, config)
  
  return(list(se_imputed = se_imputed, de_results = checkpoint_de_res))
}

preprocess_and_normalize <- function(raw_data, base_exp_name, col_mapping, exp_design,
                                     crap_ids, crapome_list, bioid_list, is_bioid = FALSE,
                                     condition_col = "condition", config, output_base_dir) {
  
  shared_logs_dir <- file.path(output_base_dir, base_exp_name, "0_Shared_Preprocessing_Logs")
  dir.create(shared_logs_dir, recursive = TRUE, showWarnings = FALSE)
  
  # ===========================================================================#
  # 1. INITIAL SETUP & DESIGN WRANGLING
  # ===========================================================================#
  colnames(exp_design) <- tolower(gsub('^\xef\xbb\xbf', '', colnames(exp_design)))
  cond_col_lower <- tolower(condition_col)
  exp_design <- exp_design %>% dplyr::filter(!is.na(label) & label != "")
  
  if (cond_col_lower != "condition" && cond_col_lower %in% colnames(exp_design)) {
    exp_design <- exp_design %>% dplyr::rename(condition = !!sym(cond_col_lower))
  } else if (!"condition" %in% colnames(exp_design)) {
    stop(paste("CRITICAL ERROR: Column '", condition_col, "' not found in design CSV."))
  }
  
  req_cols <- c("label", "condition", "replicate")
  if (!all(req_cols %in% colnames(exp_design))) stop("exp_design missing required columns.")
  
  exp_design <- exp_design %>%
    dplyr::group_by(condition) %>%
    dplyr::mutate(replicate = dplyr::row_number()) %>%
    dplyr::ungroup() %>%
    dplyr::mutate(orig_label = label) %>%
    as.data.frame()
  
  # ===========================================================================#
  # 2. RAW DATA ALIGNMENT & CLEANUP
  # ===========================================================================#
  df <- raw_data %>%
    dplyr::rename(ID = !!sym(col_mapping$id_col), name_raw = !!sym(col_mapping$name_col)) %>%
    dplyr::mutate(
      name = stringr::str_split_i(name_raw, ";", 1),
      name = ifelse(is.na(name) | name == "", ID, name)
    )
  
  lfq_cols <- match(exp_design$label, colnames(df))
  if (any(is.na(lfq_cols))) stop("Design CSV labels do not match raw data headers.")
  
  mat_clean <- as.matrix(df[, lfq_cols])
  mat_clean <- apply(mat_clean, 2, as.numeric)
  mat_clean[is.nan(mat_clean) | is.infinite(mat_clean) | mat_clean == 0] <- NA
  df[, lfq_cols] <- mat_clean
  
  # ===========================================================================#
  # 3. EXCLUDE FAILED SAMPLES
  # ===========================================================================#
  protein_counts        <- colSums(!is.na(df[, lfq_cols, drop = FALSE]))
  sample_drop_threshold <- config$pipeline_parameters$sample_drop_threshold
  cutoff                <- max(protein_counts) * sample_drop_threshold
  
  configured_controls <- c(config$experimental_design$background_controls, config$experimental_design$blank_controls)
  protected_labels <- exp_design$label[
    exp_design$condition %in% configured_controls |
      grepl("Background|IgG|Bead|Blank|Control|Buffer", exp_design$label, ignore.case = TRUE) |
      grepl("Background|IgG|Bead|Blank|Control|Buffer", exp_design$condition, ignore.case = TRUE)
  ]
  
  low_quality_samples <- names(protein_counts[protein_counts < cutoff & !names(protein_counts) %in% protected_labels])
  
  if (length(low_quality_samples) > 0) {
    message(paste(">>> Dropping poor biological samples (<", round(cutoff), "proteins):", paste(low_quality_samples, collapse = ", ")))
    excluded_samples_df <- data.frame(
      Sample = low_quality_samples, Proteins_Detected = protein_counts[low_quality_samples],
      Required_Threshold = round(cutoff), Reason = paste("Failed", sample_drop_threshold * 100, "% quality threshold")
    )
    readr::write_csv(excluded_samples_df, file.path(shared_logs_dir, "Excluded_Samples.csv"))
    
    exp_design <- exp_design %>% dplyr::filter(!label %in% low_quality_samples)
    lfq_cols   <- match(exp_design$label, colnames(df))
  }
  
  df_full <- df
  
  # ===========================================================================#
  # 4. EXCLUDE CONTAMINANTS, BLANKS & MATRIX BACKGROUND
  # ===========================================================================#
  filt_res    <- filter_contaminants(df, crap_ids, crapome_list, bioid_list, is_bioid)
  df_filtered <- filt_res$filtered
  contam_log  <- filt_res$log
  lfq_col_names <- colnames(df)[lfq_cols]
  
  if (nrow(contam_log) > 0) {
    dropped_contaminants <- df %>%
      dplyr::inner_join(contam_log %>% dplyr::select(ID, Reason), by = "ID") %>%
      dplyr::select(ID, name, Full_Description = name_raw, Reason, dplyr::all_of(lfq_col_names)) %>%
      dplyr::mutate(Samples_Detected = rowSums(!is.na(dplyr::select(., dplyr::all_of(lfq_col_names))))) %>%
      dplyr::arrange(dplyr::desc(Samples_Detected))
    readr::write_csv(dropped_contaminants, file.path(shared_logs_dir, "Excluded_Contaminants.csv"))
  }
  
  df_unique <- DEP::make_unique(df_filtered, "name", "ID", delim = ";")
  
  # ── TIER 1: System Blanks ──
  blank_conds  <- config$experimental_design$blank_controls
  blank_design <- exp_design %>% dplyr::filter(condition %in% blank_conds)
  if (nrow(blank_design) > 0) {
    blank_cols <- match(blank_design$label, colnames(df_unique))
    blank_detections <- rowSums(!is.na(df_unique[, blank_cols, drop = FALSE]))
    system_contams   <- df_unique$ID[blank_detections == length(blank_cols)]
    
    if (length(system_contams) > 0) {
      blank_col_names <- colnames(df_unique)[blank_cols]
      dropped_blanks <- df_unique %>%
        dplyr::filter(ID %in% system_contams) %>%
        dplyr::select(ID, name, Full_Description = name_raw, dplyr::all_of(blank_col_names)) %>%
        dplyr::mutate(Reason = paste0("Tier 1 Blank Exclusion (Present in ", length(blank_cols), "/", length(blank_cols), " blank tubes)"))
      readr::write_csv(dropped_blanks, file.path(shared_logs_dir, "Excluded_Tier1_SystemBlanks.csv"))
      df_unique <- df_unique %>% dplyr::filter(!ID %in% system_contams)
    }
    exp_design <- exp_design %>% dplyr::filter(!condition %in% blank_conds)
  }
  
  # ── TIER 2: Matrix Background ──
  bg_conds  <- config$experimental_design$background_controls
  bg_design <- exp_design %>% dplyr::filter(condition %in% bg_conds)
  if (nrow(bg_design) > 0) {
    bg_cols  <- match(bg_design$label, colnames(df_unique))
    exp_cols <- match(exp_design$label[!exp_design$condition %in% bg_conds], colnames(df_unique))
    
    bg_matrix <- as.matrix(df_unique[, bg_cols, drop = FALSE])
    bg_means  <- rowMeans(bg_matrix, na.rm = TRUE)
    bg_means[is.na(bg_means) | is.nan(bg_means)] <- 0 
    
    exp_matrix <- as.matrix(df_unique[, exp_cols, drop = FALSE])
    exp_conds  <- exp_design$condition[!exp_design$condition %in% bg_conds]
    
    max_exp_means <- apply(exp_matrix, 1, function(row) {
      cond_means <- tapply(as.numeric(row), exp_conds, mean, na.rm = TRUE)
      max_val    <- max(cond_means, na.rm = TRUE)
      return(ifelse(is.na(max_val) | is.infinite(max_val), 0, max_val))
    })
    
    enrichment_fc <- dplyr::case_when(
      bg_means == 0 & max_exp_means > 0 ~ Inf, 
      bg_means == 0 & max_exp_means == 0 ~ 0,   
      max_exp_means == 0                 ~ -Inf, 
      TRUE ~ log2(max_exp_means / bg_means)      
    )
    
    min_enrichment_thr <- config$pipeline_parameters$min_background_enrichment %||% 1.0
    quantitative_fails <- enrichment_fc < min_enrichment_thr
    matrix_contams <- df_unique$ID[quantitative_fails]
    
    if (length(matrix_contams) > 0) {
      bg_col_names  <- colnames(df_unique)[bg_cols]
      exp_col_names <- colnames(df_unique)[exp_cols]
      dropped_matrix <- df_unique %>%
        dplyr::filter(ID %in% matrix_contams) %>%
        dplyr::select(ID, name, Full_Description = name_raw, dplyr::all_of(bg_col_names), dplyr::all_of(exp_col_names)) %>%
        dplyr::mutate(
          Contaminant_Type = ifelse(bg_means[quantitative_fails] == 0, "Sparse Ghost Dropout", "True Matrix Binder"),
          Reason        = paste0("Tier 2 Matrix Gate (< ", min_enrichment_thr, " log2FC over IgG/Beads)"),
          Bg_Mean_Log2  = round(ifelse(bg_means[quantitative_fails] == 0, 0, log2(bg_means[quantitative_fails])), 2),
          Max_Exp_Log2  = round(ifelse(max_exp_means[quantitative_fails] == 0, 0, log2(max_exp_means[quantitative_fails])), 2),
          Enrichment_FC = round(ifelse(is.infinite(enrichment_fc[quantitative_fails]), 99.99, enrichment_fc[quantitative_fails]), 2)
        )
      readr::write_csv(dropped_matrix, file.path(shared_logs_dir, "Excluded_Tier2_BackgroundEnrichment.csv"))
    }
    df_unique <- df_unique %>% dplyr::filter(!ID %in% matrix_contams)
  }
  
  lfq_cols <- match(exp_design$label, colnames(df_unique))
  
  # ===========================================================================#
  # 5. DEP log2 OBJECT CREATION & VSN NORMALIZATION
  # ===========================================================================#
  safe_labels <- paste0("sample_", seq_along(lfq_cols))
  colnames(df_unique)[lfq_cols] <- safe_labels
  exp_design$label <- safe_labels
  
  se <- DEP::make_se(df_unique, lfq_cols, exp_design)
  pretty_names <- if ("clone_id" %in% colnames(exp_design)) paste(exp_design$clone_id, exp_design$replicate, sep = "_") else paste(exp_design$condition, exp_design$replicate, sep = "_")
  colnames(se) <- pretty_names
  SummarizedExperiment::colData(se)$name <- pretty_names
  SummarizedExperiment::colData(se)$ID   <- pretty_names
  exp_design$label <- pretty_names
  
  min_valid_prop <- config$pipeline_parameters$min_valid_proportion %||% 0.50
  min_clone_prop <- config$pipeline_parameters$min_clone_proportion %||% 0.50
  
  is_na_mat <- is.na(SummarizedExperiment::assay(se))
  col_data  <- as.data.frame(SummarizedExperiment::colData(se))
  
  if ("clone_id" %in% colnames(col_data)) {
    design_structure <- lapply(split(seq_len(nrow(col_data)), col_data$condition), function(geno_idx) split(geno_idx, col_data$clone_id[geno_idx]))
    keep_rows <- apply(!is_na_mat, 1, function(present) {
      any(sapply(design_structure, function(clones_in_geno) {
        n_clones_total   <- length(clones_in_geno)
        n_clones_passing <- sum(sapply(clones_in_geno, function(clone_idx) {
          sum(present[clone_idx]) >= ceiling(length(clone_idx) * min_valid_prop)
        }))
        n_clones_passing >= ceiling(n_clones_total * min_clone_prop)
      }))
    })
    
    threshold_log <- lapply(names(design_structure), function(geno) {
      lapply(names(design_structure[[geno]]), function(clone) {
        data.frame(Genotype = geno, Clone = clone, N_Replicates = length(design_structure[[geno]][[clone]]),
                   Required_Valid = ceiling(length(design_structure[[geno]][[clone]]) * min_valid_prop), Proportion_Applied = min_valid_prop)
      }) %>% dplyr::bind_rows()
    }) %>% dplyr::bind_rows() %>% dplyr::mutate(N_Clones_In_Genotype = sapply(Genotype, function(g) length(design_structure[[g]])), Required_Clones = ceiling(N_Clones_In_Genotype * min_clone_prop))
    readr::write_csv(threshold_log, file.path(shared_logs_dir, "Filter_Thresholds_MissingValues.csv"))
  } else {
    keep_rows <- apply(!is_na_mat, 1, function(present) {
      any(tapply(present, col_data$condition, function(x) sum(x) >= ceiling(length(x) * min_valid_prop)))
    })
  }
  
  dropped_by_missingness <- SummarizedExperiment::rowData(se)$name[!keep_rows]
  if (length(dropped_by_missingness) > 0) {
    dropped_mv_df <- df_unique %>% dplyr::filter(name %in% dropped_by_missingness) %>%
      dplyr::select(ID, name, Full_Description = name_raw, dplyr::all_of(safe_labels)) %>%
      dplyr::mutate(Reason = "Failed clone-aware missing value filter", Total_Detections = rowSums(!is.na(dplyr::select(., dplyr::all_of(safe_labels)))))
    readr::write_csv(dropped_mv_df, file.path(shared_logs_dir, "Excluded_Proteins_MissingValues.csv"))
  }
  
  se_filtered <- se[keep_rows, ]
  se_norm     <- DEP::normalize_vsn(se_filtered)
  
  return(list(df_full = df_full, df_unique = df_unique, exp_design_clean = exp_design, se_filtered = se_filtered, se_norm = se_norm))
}

execute_analysis_mode <- function(prepped_data, base_exp_name, mode, scope_folder, active_design,
                                  contrasts, output_base_dir, imputation_method = "mixed",
                                  fdr_thr = 0.05, lfc_thr = 1.0, poi_list = c(), config) {
  
  message(paste("\n---> Executing Downstream Branch:", base_exp_name, "| Mode:", mode, "| Scope:", scope_folder))
  
  # ── Setup Nested Hierarchy ──────────────────────────────────────────────────
  mode_folder <- ifelse(mode == "abundance_only", "Abundance_Only", "Standard_Imputed")
  out_dir  <- file.path(output_base_dir, base_exp_name, mode_folder, scope_folder, base_exp_name)
  
  logs_dir <- file.path(out_dir, "0_Logs_and_Exclusions")
  qc_dir   <- file.path(out_dir, "1_QC_and_Normalization")
  de_dir   <- file.path(out_dir, ifelse(mode == "abundance_only", "2_Differential_Abundance", "2_Differential_Expression"))
  
  dir.create(logs_dir, recursive = TRUE, showWarnings = FALSE)
  dir.create(qc_dir,   recursive = TRUE, showWarnings = FALSE)
  dir.create(de_dir,   recursive = TRUE, showWarnings = FALSE)
  
  # Copy shared pre-processing logs into this specific mode's folder so generate_qc_report works seamlessly
  shared_logs_dir <- file.path(output_base_dir, base_exp_name, "0_Shared_Preprocessing_Logs")
  if (dir.exists(shared_logs_dir)) {
    file.copy(list.files(shared_logs_dir, full.names = TRUE), logs_dir, overwrite = TRUE)
  }
  
  # ── The Scope Intercept (Late-Filtering) ────────────────────────────────────
  target_conditions <- unique(active_design$condition)
  exp_design <- prepped_data$exp_design_clean %>% dplyr::filter(condition %in% target_conditions)
  
  se_filtered <- prepped_data$se_filtered[, exp_design$label]
  se_norm     <- prepped_data$se_norm[, exp_design$label]
  df_unique   <- prepped_data$df_unique
  
  # ── Localized QC0-QC4 Plots ─────────────────────────────────────────────────
  save_and_print_plot(DEP::plot_numbers(se_filtered), file.path(qc_dir, "QC1_Protein_Numbers"))
  save_and_print_plot(DEP::plot_missval(se_filtered), file.path(qc_dir, "QC2_Missing_Values"))
  save_and_print_plot(DEP::plot_normalization(se_filtered, se_norm), file.path(qc_dir, "QC3_VSN_Density"))
  
  norm_mat     <- SummarizedExperiment::assay(se_norm)
  vsn_plot_obj <- vsn::meanSdPlot(norm_mat, plot = FALSE)
  save_and_print_plot(vsn_plot_obj$gg + ggplot2::theme_bw() + ggplot2::ggtitle(paste(base_exp_name, "- VSN Mean-Variance Fit")), 
                      file.path(qc_dir, "QC4_VSN_MeanSdPlot"), width = 7, height = 7)
  
  raw_long <- as.data.frame(SummarizedExperiment::assay(se_filtered)) %>%
    tibble::rownames_to_column("Protein") %>% tidyr::pivot_longer(-Protein, names_to = "Sample", values_to = "Intensity") %>%
    dplyr::filter(!is.na(Intensity)) %>% dplyr::left_join(exp_design, by = c("Sample" = "label")) %>%
    dplyr::filter(!condition %in% config$experimental_design$background_controls)
  
  p_raw_box <- ggplot2::ggplot(raw_long, ggplot2::aes(x = Sample, y = Intensity, fill = condition)) +
    ggplot2::geom_boxplot(outlier.shape = NA, alpha = 0.7) + ggplot2::theme_bw() +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 45, hjust = 1)) +
    ggplot2::labs(title = paste(base_exp_name, "- Raw Intensity Distribution"), y = "Raw Intensity", x = "")
  save_and_print_plot(p_raw_box, file.path(qc_dir, "QC0_Raw_Intensity_Distribution"))
  
  # ── Branch: Abundance Only Mode ─────────────────────────────────────────────
  if (mode == "abundance_only") {
    scrub_res <- scrub_by_clone_genotype_presence(
      mat = norm_mat, exp_design = exp_design,
      min_valid_prop = config$pipeline_parameters$min_valid_proportion %||% 0.5,
      min_clone_prop = config$pipeline_parameters$min_clone_proportion %||% 0.5
    )
    
    id_map <- df_unique %>% dplyr::select(Protein = name, ID) %>% dplyr::distinct()
    report <- export_interactome_abundance_report(scrub_res$mat_scrubbed, exp_design, scrub_res$na_log, de_dir, base_exp_name, id_map)
    
    ordered_contrasts <- contrasts
    bl_contrast <- config$comparative_analysis$blacklist_contrast
    if (!is.null(bl_contrast) && bl_contrast %in% contrasts) ordered_contrasts <- c(bl_contrast, setdiff(contrasts, bl_contrast))
    
    raw_stats_list <- list()
    ko_binders <- c()
    
    for (contrast in ordered_contrasts) {
      tryCatch({
        res_df <- calculate_raw_abundance_stats(
          mat_scrubbed = scrub_res$mat_scrubbed, exp_design = exp_design, contrast_str = contrast,
          output_base_dir = out_dir, exp_name = base_exp_name, id_map = id_map, poi_list = poi_list, blacklist_proteins = ko_binders
        )
        raw_stats_list[[contrast]] <- res_df
        
        if (!is.null(bl_contrast) && contrast == bl_contrast) {
          target_cond <- stringr::str_trim(stringr::str_split(contrast, "-")[[1]][1])
          ko_binders <- res_df %>% dplyr::filter((Status == "Quantitative" & P_Value < fdr_thr & logFC > lfc_thr) | Status == paste("Exclusive to", target_cond)) %>% dplyr::pull(Protein)
        }
      }, error = function(e) message("      -> Skipping stats for ", contrast, ": ", e$message))
    }
    
    if (length(poi_list) > 0) save_and_print_plot(plot_protein_abundance_strip(report$long, poi_list, "Protein", base_exp_name), file.path(de_dir, paste0(base_exp_name, "_Abundance_StripPlot_POI")))
    save_and_print_plot(plot_abundance_heatmap(report$clone_means, exp_name = base_exp_name), file.path(de_dir, paste0(base_exp_name, "_Abundance_Heatmap_AllProteins")), width = 10, height = 14)
    
    vis_design <- exp_design %>% dplyr::filter(!condition %in% config$experimental_design$background_controls)
    vis_matrix <- norm_mat[, vis_design$label, drop = FALSE]
    
    generate_mds_plot(vis_matrix, vis_design, base_exp_name, qc_dir)
    add_rank_plot(vis_matrix, vis_design, config$experimental_design$bait_protein, poi_list, file.path(qc_dir, "QC10_Protein_Rank_Plot"), base_exp_name)
    
    cor_mat <- cor(vis_matrix, method = "spearman", use = "pairwise.complete.obs")
    p_cor <- pheatmap::pheatmap(cor_mat, main = paste(base_exp_name, "- Pre-Imputation Spearman"), display_numbers = round(cor_mat, 2), number_format = "%.2f", fontsize_number = 7, color = colorRampPalette(c("#2166AC", "#F7F7F7", "#D6604D"))(100), breaks = seq(0.7, 1.0, length.out = 101), annotation_col = data.frame(Condition = vis_design$condition, row.names = colnames(vis_matrix)), clustering_method = "ward.D2", angle_col = 45, silent = TRUE)
    save_and_print_plot(p_cor, file.path(qc_dir, "QC8_Correlation_Matrix"), width = 9, height = 8)
    
    se_imputed <- se_norm 
    checkpoint_de_res <- raw_stats_list
    names(checkpoint_de_res) <- make_safe_contrast_name(names(raw_stats_list))
    
  } else {
    # ── Branch: Standard Imputed Mode ───────────────────────────────────────────
    if (imputation_method == "mixed") {
      se_imputed <- impute_mixed_mar_mnar(se_norm, min_missing_pct = config$pipeline_parameters$imputation$mnar_threshold)
    } else if (imputation_method == "knn") {
      se_imputed <- DEP::impute(se_norm, fun = "knn", rowmax = 0.9)
    } else {
      se_imputed <- DEP::impute(se_norm, fun = "MinProb", q = config$pipeline_parameters$imputation$minprob_q)
    }
    
    save_and_print_plot(plot_imputation_per_sample(se_norm, se_imputed, paste(base_exp_name, "- Imputation Shift per Sample")), file.path(qc_dir, "QC5_Imputation_Shift_PerSample"), width = 12, height = 10)
    
    imputed_matrix <- SummarizedExperiment::assay(se_imputed)
    vis_design <- exp_design %>% dplyr::filter(!condition %in% config$experimental_design$background_controls)
    vis_matrix <- imputed_matrix[, vis_design$label, drop = FALSE]
    
    export_imputation_map_excel(se_norm, imputed_matrix, exp_design, file.path(qc_dir, paste0("QC5b_Imputation_ColorMap_", base_exp_name, ".xlsx")), config$pipeline_parameters$imputation$mnar_threshold)
    
    vars <- apply(vis_matrix, 1, var)
    heatmap_mat_top <- t(scale(t(vis_matrix[order(vars, decreasing = TRUE)[1:min(100, nrow(vis_matrix))], ])))
    save_and_print_plot(pheatmap::pheatmap(heatmap_mat_top, main = paste(base_exp_name, "- Top 100 Variable Proteins (Z-score)"), clustering_method = "ward.D2", show_rownames = nrow(heatmap_mat_top) <= 100, fontsize_row = 6, color = colorRampPalette(c("navy", "white", "firebrick3"))(100)), file.path(qc_dir, "QC6_Top_Variable_Proteins"), width = 10, height = 16)
    
    generate_pca_full(vis_matrix, vis_design, paste(base_exp_name, "PCA"), 20, qc_dir)
    add_rank_plot(vis_matrix, vis_design, config$experimental_design$bait_protein, poi_list, file.path(qc_dir, "QC10_Protein_Rank_Plot"), base_exp_name)
    generate_correlation_matrices(se_norm[, vis_design$label], se_imputed[, vis_design$label], vis_design, qc_dir, base_exp_name)
    
    readr::write_csv(as.data.frame(norm_mat) %>% tibble::rownames_to_column("Protein"), file.path(qc_dir, "Matrix_2_VSN_Normalized.csv"))
    readr::write_csv(as.data.frame(imputed_matrix) %>% tibble::rownames_to_column("Protein"), file.path(qc_dir, "Matrix_3_Fully_Imputed.csv"))
    
    limma_res <- perform_limma_analysis(imputed_matrix, exp_design, contrasts, exp_design$clone_id, "condition")
    
    bl_contrast <- config$comparative_analysis$blacklist_contrast
    black_safe  <- if (!is.null(bl_contrast)) make_safe_contrast_name(bl_contrast) else NULL
    ko_binders  <- c()
    
    if (!is.null(black_safe) && black_safe %in% names(limma_res)) {
      ko_binders_df <- limma_res[[black_safe]] %>% dplyr::filter(adj.P.Val < fdr_thr & logFC > 0)
      ko_binders    <- ko_binders_df$Protein
      if (length(ko_binders) > 0) {
        id_map <- df_unique %>% dplyr::select(ID, Protein = name) %>% dplyr::distinct()
        readr::write_csv(ko_binders_df %>% dplyr::left_join(id_map, by = "Protein") %>% dplyr::select(ID, Protein, logFC, adj.P.Val) %>% dplyr::mutate(Reason = "Statistical Blacklist (IP-MS Background)"), file.path(logs_dir, "Excluded_Statistical_Blacklist.csv"))
      }
    }
    
    wb_de <- openxlsx::createWorkbook()
    for (contrast_name in names(limma_res)) {
      safe_name <- make_safe_contrast_name(contrast_name)
      plot_df   <- limma_res[[contrast_name]]
      
      readr::write_csv(plot_df, file.path(de_dir, paste0(safe_name, "_Results.csv")))
      if (!is.null(black_safe) && safe_name != black_safe && length(ko_binders) > 0) plot_df <- plot_df %>% dplyr::filter(!Protein %in% ko_binders)
      
      save_and_print_plot(generate_universal_volcano_plot(plot_df, contrast_name, base_exp_name, fdr_thr, lfc_thr, poi_list, 15, ko_binders, nrow(plot_df) < 200), file.path(de_dir, paste0("Volcano_", safe_name)), width = 10, height = 8)
      save_and_print_plot(generate_ma_plot(plot_df, paste(base_exp_name, "-", contrast_name), fdr_thr, lfc_thr), file.path(de_dir, paste0("MA_", safe_name)))
      
      openxlsx::addWorksheet(wb_de, substr(safe_name, 1, 31))
      openxlsx::writeDataTable(wb_de, substr(safe_name, 1, 31), plot_df %>% dplyr::arrange(adj.P.Val), tableStyle = "TableStyleMedium9")
    }
    openxlsx::saveWorkbook(wb_de, file.path(de_dir, paste0(base_exp_name, "_All_DE_Results.xlsx")), overwrite = TRUE)
    
    checkpoint_de_res <- limma_res
    names(checkpoint_de_res) <- make_safe_contrast_name(names(limma_res))
  }
  
  # ── Save Standardized Checkpoint ────────────────────────────────────────────
  checkpoint <- list(
    exp_design_clean      = exp_design,
    raw_df_aligned        = prepped_data$df_full,
    imputed_matrix        = SummarizedExperiment::assay(se_imputed), 
    norm_matrix           = norm_mat,
    de_results            = checkpoint_de_res,
    de_results_safe_names = names(checkpoint_de_res)
  )
  
  saveRDS(checkpoint, file.path(out_dir, paste0(base_exp_name, "_checkpoint.rds")))
  
  # Run HTML QC Report utilizing the localized output path
  generate_qc_report(base_exp_name, file.path(output_base_dir, base_exp_name, mode_folder, scope_folder), config)
  
  return(list(se_imputed = se_imputed, de_results = checkpoint_de_res))
}

generate_integrated_report <- function(output_base_dir, config, run_log = NULL) {
  
  if (!requireNamespace("base64enc", quietly = TRUE)) {
    install.packages("base64enc")
    library(base64enc)
  }
  
  safe_path <- function(x) gsub("\\\\", "/", x)
  
  integ_dir  <- file.path(output_base_dir, "Integrated_Analysis")
  comp_dir   <- file.path(integ_dir, "Method_Comparison")
  func_dir   <- file.path(integ_dir, "Functional_Enrichment")
  net_dir    <- file.path(integ_dir, "STRING_Networks")
  quad_dir   <- file.path(integ_dir, "Integration_Global_vs_IP")
  stoich_dir <- file.path(integ_dir, "Stoichiometry")
  
  # ── Helper: encode TIFF as base64 PNG for self-contained embedding ───────────
  if (!requireNamespace("magick", quietly = TRUE)) {
    install.packages("magick")
  }
  
  encode_image <- function(path) {
    if (is.null(path) || !file.exists(path)) return(NULL)
    
    # Read the TIFF from disk
    img <- tryCatch(magick::image_read(path), error = function(e) NULL)
    if (is.null(img)) return(NULL)
    
    # Convert to PNG in memory and write to a raw vector
    img_png <- magick::image_convert(img, "png")
    raw_png <- magick::image_write(img_png, format = "png")
    
    # Encode the raw vector
    paste0("data:image/png;base64,", base64enc::base64encode(raw_png))
  }
  
  make_image <- function(tiff_path, width = "90%") {
    encoded <- encode_image(tiff_path)
    if (is.null(encoded)) return("<p><em>Image not found or conversion failed.</em></p>")
    paste0('<img src="', encoded, '" style="width:', width,
           '; display:block; margin:auto; border-radius:4px;',
           ' box-shadow:0 2px 8px rgba(0,0,0,0.1); margin-bottom:12px;">')
  }
  
  # ── Helper: data.frame to HTML table ─────────────────────────────────────────
  df_to_html_table <- function(df, highlight_col = NULL) {
    if (is.null(df) || nrow(df) == 0) return("<p><em>No data available.</em></p>")
    header <- paste0("<thead><tr>",
                     paste0("<th>", colnames(df), "</th>", collapse = ""),
                     "</tr></thead>")
    rows <- apply(df, 1, function(row) {
      highlight <- !is.null(highlight_col) &&
        highlight_col %in% names(row) &&
        isTRUE(as.logical(row[[highlight_col]]))
      style <- if (highlight) ' style="background:#FFF2CC"' else ""
      paste0("<tr", style, ">",
             paste0("<td>", row, "</td>", collapse = ""),
             "</tr>")
    })
    paste0('<div class="table-wrapper"><table class="summary-table">',
           header, "<tbody>", paste(rows, collapse = ""),
           "</tbody></table></div>")
  }
  
  # ── Helper: collapsible section ───────────────────────────────────────────────
  make_section <- function(title, content, guidance = NULL,
                           open = FALSE, top_level = FALSE) {
    open_attr    <- if (open) " open" else ""
    summary_class <- if (top_level) "section-header top-level" else "section-header"
    guidance_html <- if (!is.null(guidance)) {
      paste0('<div class="guidance"><p>', guidance, '</p></div>')
    } else ""
    paste0(
      '<details', open_attr, '>',
      '<summary class="', summary_class, '">', title, '</summary>',
      '<div class="section-content">',
      guidance_html, content,
      '</div></details>'
    )
  }
  
  # ── Collect tiff plots (absolute paths) ───────────────────────────────────────
  comp_tiffs   <- sort(list.files(comp_dir,   pattern = "\\.tiff$", full.names = TRUE))
  func_tiffs   <- sort(list.files(func_dir,   pattern = "\\.tiff$", full.names = TRUE))
  net_tiffs    <- sort(list.files(net_dir,    pattern = "\\.tiff$", full.names = TRUE, recursive = TRUE))
  quad_tiffs   <- sort(list.files(quad_dir,   pattern = "\\.tiff$", full.names = TRUE))
  stoich_tiffs <- sort(list.files(stoich_dir, pattern = "\\.tiff$", full.names = TRUE))
  sens_pngs    <- sort(list.files(file.path(output_base_dir, "Integrated_Analysis/Imputation_Sensitivity"), pattern = "\\.tiff?$", full.names = TRUE))  
  # ── Precompute tables ─────────────────────────────────────────────────────────
  
  # Run log
  run_log_df <- if (!is.null(run_log)) {
    dplyr::bind_rows(lapply(names(run_log), function(s) {
      data.frame(Stage   = s,
                 Status  = run_log[[s]]$status,
                 Time    = run_log[[s]]$time,
                 Message = run_log[[s]]$message)
    }))
  } else NULL
  
  # DE summary across all experiments
  de_summary <- lapply(names(config$experiments), function(exp_name) {
    ckpt_file <- file.path(output_base_dir, exp_name,
                           paste0(exp_name, "_checkpoint.rds"))
    if (!file.exists(ckpt_file)) return(NULL)
    ckpt <- tryCatch(load_checkpoint(ckpt_file), error = function(e) NULL)
    if (is.null(ckpt)) return(NULL)
    lapply(names(ckpt$de_results), function(cn) {
      df <- ckpt$de_results[[cn]]
      data.frame(
        Experiment = exp_name, Contrast = cn,
        N_Proteins = nrow(df),
        N_Enriched = sum(df$adj.P.Val < config$thresholds$fdr &
                           df$logFC >  config$thresholds$lfc, na.rm = TRUE),
        N_Depleted = sum(df$adj.P.Val < config$thresholds$fdr &
                           df$logFC < -config$thresholds$lfc, na.rm = TRUE)
      )
    }) %>% dplyr::bind_rows()
  }) %>% dplyr::bind_rows()
  
  # Cleaned interactome summary
  interactome_files <- list.files(comp_dir,
                                  pattern = "Final_Cleaned_Interactome.*\\.csv$",
                                  full.names = TRUE)
  interactome_summary <- if (length(interactome_files) > 0) {
    lapply(interactome_files, function(f) {
      df  <- tryCatch(readr::read_csv(f, show_col_types = FALSE),
                      error = function(e) NULL)
      if (is.null(df)) return(NULL)
      exp <- gsub("Final_Cleaned_Interactome_|\\.csv", "", basename(f))
      data.frame(
        Experiment    = exp,
        N_Interactors = nrow(df),
        Top_Hit       = if (nrow(df) > 0) df$Protein[1]              else NA,
        Top_Hit_logFC = if (nrow(df) > 0) round(df$logFC[1], 2)   else NA,
        Top_Hit_FDR   = if (nrow(df) > 0) signif(df$adj.P.Val[1], 3) else NA
      )
    }) %>% dplyr::bind_rows()
  } else NULL
  
  # Stoichiometry summary
  stoich_files <- list.files(stoich_dir,
                             pattern = "Stoichiometry_Stats.*\\.csv$",
                             full.names = TRUE)
  stoich_summary <- if (length(stoich_files) > 0) {
    lapply(stoich_files, function(f) {
      df <- tryCatch(readr::read_csv(f, show_col_types = FALSE),
                     error = function(e) NULL)
      if (is.null(df)) return(NULL)
      data.frame(
        Contrast   = gsub("Stoichiometry_Stats_|\\.csv", "", basename(f)),
        N_Proteins = nrow(df),
        N_Sig_Up   = sum(df$adj.P.Val < config$thresholds$fdr &
                           df$logFC >  config$thresholds$lfc, na.rm = TRUE),
        N_Sig_Down = sum(df$adj.P.Val < config$thresholds$fdr &
                           df$logFC < -config$thresholds$lfc, na.rm = TRUE)
      )
    }) %>% dplyr::bind_rows()
  } else NULL
  
  # ── Build sections ────────────────────────────────────────────────────────────
  sections <- list()
  
  # 1. Run summary
  run_log_html <- if (!is.null(run_log_df)) {
    df_to_html_table(run_log_df)
  } else "<p><em>Run log not available.</em></p>"
  
  sections[["run_summary"]] <- make_section(
    "Pipeline Run Summary",
    paste0(
      '<p>Bait: <strong>', config$experimental_design$bait_protein,
      '</strong> | FDR: <strong>', config$thresholds$fdr,
      '</strong> | LFC: <strong>&plusmn;', config$thresholds$lfc, '</strong></p>',
      run_log_html
    ),
    open = TRUE, top_level = TRUE
  )
  
  # 2. DE overview
  sections[["de_overview"]] <- make_section(
    "Differential Expression Overview",
    df_to_html_table(de_summary),
    guidance = paste0(
      "Significant hits per contrast across all experiments ",
      "(FDR < ", config$thresholds$fdr,
      ", |LFC| > ", config$thresholds$lfc, ")."
    ),
    top_level = TRUE
  )
  
  # 3. Cleaned interactome
  interactome_html <- if (!is.null(interactome_summary)) {
    df_to_html_table(interactome_summary)
  } else "<p><em>No cleaned interactome files found. Run 02_comparative_analysis.R first.</em></p>"
  
  sections[["interactome"]] <- make_section(
    "Cleaned Interactome Summary",
    interactome_html,
    guidance = paste0(
      "Proteins identified as true interactors after background subtraction. ",
      "Blacklist contrast: <strong>", config$comparative_analysis$blacklist_contrast,
      "</strong>. Target contrast: <strong>",
      config$comparative_analysis$target_contrast, "</strong>."
    ),
    top_level = TRUE
  )
  
  # 4. Method comparison plots — one collapsible per plot
  if (length(comp_tiffs) > 0) {
    comp_inner <- paste(sapply(comp_tiffs, function(p) {
      title <- gsub("_", " ", tools::file_path_sans_ext(basename(p)))
      make_section(title, make_image(p))
    }), collapse = "\n")
    sections[["method_comparison"]] <- make_section(
      "Method Comparison",
      comp_inner,
      guidance = "UpSet plots and known interactor recovery rates after background subtraction.",
      top_level = TRUE
    )
  }
  
  # 5. Functional enrichment — grouped by experiment then one per plot
  if (length(func_tiffs) > 0) {
    
    # Group plots by experiment name for nested collapsibles
    exp_names_in_plots <- unique(sapply(func_tiffs, function(p) {
      # Extract experiment name — filename pattern: GSEA_X_ExpName_Contrast.png
      parts <- stringr::str_match(basename(p), "GSEA_[^_]+_([^_]+(?:_[^_]+)?)_")
      if (!is.na(parts[1,2])) parts[1,2] else "Other"
    }))
    
    func_inner <- paste(sapply(func_tiffs, function(p) {
      title <- gsub("_", " ", tools::file_path_sans_ext(basename(p)))
      make_section(title, make_image(p))
    }), collapse = "\n")
    
    sections[["gsea"]] <- make_section(
      "Functional Enrichment (GSEA)",
      func_inner,
      guidance = paste0(
        "GSEA results for DNA damage response and cell cycle gene sets. ",
        "NES = Normalised Enrichment Score. Positive NES = enrichment in the ",
        "first group of the contrast. Only contrasts with at least one pathway ",
        "at FDR < 0.25 produce cnetplot and heatplot outputs."
      ),
      top_level = TRUE
    )
  }
  
  # 6. STRING networks — one collapsible per plot
  if (length(net_tiffs) > 0) {
    net_inner <- paste(sapply(net_tiffs, function(p) {
      title <- gsub("_", " ", tools::file_path_sans_ext(basename(p)))
      make_section(title, make_image(p))
    }), collapse = "\n")
    sections[["string_networks"]] <- make_section(
      "Protein Interaction Networks (STRING)",
      net_inner,
      guidance = paste0(
        "PPI networks for the top 20 significant hits per contrast, ",
        "coloured by log2 fold change. Edges reflect STRING combined score > 400."
      ),
      top_level = TRUE
    )
  }
  
  # 7. Quadrant plots
  if (length(quad_tiffs) > 0) {
    quad_inner <- paste(sapply(quad_tiffs, function(p) {
      title <- gsub("_", " ", tools::file_path_sans_ext(basename(p)))
      make_section(title, make_image(p))
    }), collapse = "\n")
    sections[["quadrant"]] <- make_section(
      "Global Proteome vs IP-MS Integration",
      quad_inner,
      guidance = paste0(
        "<strong>Affinity Driven</strong>: significant in IP but not global ",
        "(enriched by pulldown, not abundance change). ",
        "<strong>Abundance Driven</strong>: significant in the same direction in both."
      ),
      top_level = TRUE
    )
  }
  
  # 8. Stoichiometry
  stoich_content <- if (!is.null(stoich_summary)) {
    df_to_html_table(stoich_summary)
  } else ""
  if (length(stoich_tiffs) > 0) {
    stoich_plots_html <- paste(sapply(stoich_tiffs, function(p) {
      title <- gsub("_", " ", tools::file_path_sans_ext(basename(p)))
      make_section(title, make_image(p))
    }), collapse = "\n")
    stoich_content <- paste0(stoich_content, stoich_plots_html)
  }
  if (nchar(stoich_content) > 0) {
    sections[["stoichiometry"]] <- make_section(
      "Stoichiometry Analysis",
      stoich_content,
      guidance = paste0(
        "Log2(Prey/Bait) ratios relative to <strong>",
        config$experimental_design$bait_protein,
        "</strong>. A ratio of 0 means the prey tracks perfectly with the bait. ",
        "Significant changes indicate altered stoichiometry independent of bait abundance."
      ),
      top_level = TRUE
    )
  }
  
  # 7b. Imputation Sensitivity Analysis
  if (length(sens_pngs) > 0) {
    sens_inner <- paste(sapply(sens_pngs, function(p) {
      title <- gsub("_", " ", tools::file_path_sans_ext(basename(p)))
      make_section(title, make_image(p))
    }), collapse = "\n")
    sections[["imputation_sensitivity"]] <- make_section(
      "Global Imputation Sensitivity Audit",
      sens_inner,
      guidance = "<strong>Robust Hits</strong> remain invariant across pipelines. <strong>Ghost Hits</strong> are false-positive artifacts generated by imputation noise filling missing rows. <strong>Binary Dropouts</strong> are hidden biological switches showing full presence/absence dropouts.",
      top_level = TRUE
    )
  }
  
  # 9. Links to per-experiment QC reports
  qc_links <- paste(sapply(names(config$experiments), function(exp_name) {
    html_path <- file.path(output_base_dir, exp_name,
                           paste0(exp_name, "_QC_Report.html"))
    if (file.exists(html_path)) {
      abs_path <- safe_path(normalizePath(html_path, mustWork = FALSE))
      paste0('<li><a href="', abs_path, '" target="_blank">',
             exp_name, ' QC Report</a></li>')
    } else {
      paste0('<li>', exp_name, ' \u2014 QC report not found</li>')
    }
  }), collapse = "\n")
  
  sections[["qc_links"]] <- make_section(
    "Per-Experiment QC Reports",
    paste0('<ul>', qc_links, '</ul>'),
    top_level = TRUE
  )
  
  # 10. Session info
  si_text <- paste(utils::capture.output(sessionInfo()), collapse = "\n")
  sections[["session"]] <- make_section(
    "Reproducibility",
    paste0('<pre class="session-info">', si_text, '</pre>'),
    top_level = TRUE
  )
  
  # ── CSS ───────────────────────────────────────────────────────────────────────
  css <- '
    body { font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", sans-serif;
           max-width: 1200px; margin: 0 auto; padding: 20px; color: #333; }
    h1   { color: #2c3e50; border-bottom: 2px solid #2c3e50; padding-bottom: 8px; }

    details { margin: 6px 0; border: 1px solid #ddd; border-radius: 6px; overflow: hidden; }
    details > details { margin: 4px 8px; border-color: #e8e8e8; }

    summary.section-header {
      padding: 10px 16px; cursor: pointer; font-weight: 600;
      background: #f0f4f8; color: #2c3e50;
      list-style: none; display: flex; align-items: center; gap: 8px;
      user-select: none;
    }
    summary.section-header.top-level {
      background: #2c3e50; color: white; font-size: 1.05em; padding: 12px 16px;
    }
    summary.section-header::before {
      content: "\25B6"; font-size: 0.8em; transition: transform 0.2s;
    }
    details[open] > summary.section-header::before { transform: rotate(90deg); }

    .section-content { padding: 12px 16px; }
    .guidance {
      background: #f8f9fa; border-left: 4px solid #4a90d9;
      padding: 8px 12px; margin-bottom: 12px; border-radius: 0 4px 4px 0;
      font-size: 0.92em; color: #555;
    }

    .table-wrapper { overflow-x: auto; margin: 8px 0; }
    table.summary-table {
      border-collapse: collapse; width: 100%; font-size: 0.88em;
    }
    table.summary-table th {
      background: #2c3e50; color: white; padding: 8px 12px; text-align: left;
    }
    table.summary-table td { padding: 6px 12px; border-bottom: 1px solid #eee; }
    table.summary-table tr:hover td { background: #f5f5f5; }

    img {
      border-radius: 4px; box-shadow: 0 2px 8px rgba(0,0,0,0.1); margin: 8px 0;
    }
    a { color: #4a90d9; }
    ul { padding-left: 20px; line-height: 2; }
    pre.session-info {
      background: #f8f8f8; padding: 12px; border-radius: 4px;
      font-size: 0.78em; overflow-x: auto; white-space: pre-wrap;
    }
  '
  
  # ── Assemble and write ────────────────────────────────────────────────────────
  html <- paste0(
    '<!DOCTYPE html><html lang="en"><head>',
    '<meta charset="UTF-8">',
    '<meta name="viewport" content="width=device-width, initial-scale=1.0">',
    '<title>Integrated Analysis Report: FAN1 Proteomics</title>',
    '<style>', css, '</style>',
    '</head><body>',
    '<h1>Integrated Analysis Report: FAN1 Proteomics</h1>',
    '<p style="color:#888; font-size:0.9em;">Generated ',
    format(Sys.Date(), "%d %B %Y"),
    ' &bull; Hybrid DEP/Limma Pipeline</p>',
    paste(sections, collapse = "\n"),
    '</body></html>'
  )
  
  out_path <- file.path(output_base_dir, "Integrated_Analysis_Report.html")
  writeLines(html, out_path)
  message(paste("   -> Integrated report generated:", safe_path(out_path)))
}

