# =============================================================================
 Data Retrieval, Pre-processing, Batch Effect Correction & PCA
# =============================================================================

# --- 1. Install and load required packages ---
pkgs <- c("PharmacoGx", "dplyr", "ggplot2", "VIM", "sva", "biomaRt", "Biobase", "SummarizedExperiment")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) {
    if (p %in% c("PharmacoGx", "sva", "biomaRt")) {
      BiocManager::install(p, update = FALSE, ask = FALSE)
    } else {
      install.packages(p)
    }
  }
  library(p, character.only = TRUE)
}))

# Logging function
log_message <- function(message, file = "logs/preprocess.log") {
  dir.create("logs", showWarnings = FALSE)
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, message), file = file, append = TRUE)
}

log_message("Starting data preprocessing and batch correction")

set.seed(123)

# --- 2. Load raw datasets (place in data/raw/) ---
log_message("Loading raw PharmacoSets")
CCLE <- readRDS("CCLE.rds")
GDSC <- readRDS("GDSC2.rds")

CCLE <- updateObject(CCLE)
GDSC <- updateObject(GDSC)

log_message("Datasets loaded and updated successfully")

# --- 3. Helper functions ---
validate_genes <- function(genes, dataset_name) {
  log_message(sprintf("Validating gene IDs for %s", dataset_name))
  tryCatch({
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    gene_mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                          filters = "ensembl_gene_id",
                          values = head(genes, 10), mart = mart)
    valid_count <- sum(gene_mapping$hgnc_symbol != "")
    log_message(sprintf("Validated %d/10 sampled gene IDs for %s", valid_count, dataset_name))
  }, error = function(e) {
    log_message(sprintf("Gene validation skipped for %s: %s", dataset_name, conditionMessage(e)))
  })
}

normalize_data <- function(expr) {
  expr <- expr + 1
  expr <- log2(expr)
  t(scale(t(expr)))
}

filter_features <- function(data) {
  log_message("Applying MAD-based feature filtering on combined data")
  initial_genes <- nrow(data)
  
  feature_mad <- apply(data, 1, function(x) {
    if (sum(!is.na(x)) < 2) return(NA)
    med <- median(x, na.rm = TRUE)
    if (is.na(med)) return(NA)
    mad_val <- median(abs(x - med), na.rm = TRUE)
    if (is.na(mad_val) || mad_val == 0) return(NA)
    mad_val
  })
  
  valid_rows <- !is.na(feature_mad)
  data <- data[valid_rows, , drop = FALSE]
  feature_mad <- feature_mad[valid_rows]
  
  mad_threshold <- quantile(feature_mad, 0.05, na.rm = TRUE)
  keep_rows <- feature_mad > mad_threshold
  data <- data[keep_rows, , drop = FALSE]
  
  log_message(sprintf("MAD filtering: %d → %d genes", initial_genes, nrow(data)))
  data
}

# --- 4. Preprocess individual dataset (RNA expression) ---
preprocess_dataset <- function(pset, dataset_name, mDataType = "rna") {
  log_message(sprintf("Preprocessing RNA expression for %s", dataset_name))
  
  se <- summarizeMolecularProfiles(pset, mDataType = mDataType, fill.missing = FALSE)
  expr <- assay(se)
  
  log_message(sprintf("%s: %d genes × %d cell lines", dataset_name, nrow(expr), ncol(expr)))
  validate_genes(rownames(expr), dataset_name)
  
  # Simple KNN imputation if needed
  if (any(is.na(expr))) {
    log_message(sprintf("Imputing missing values in %s using kNN", dataset_name))
    expr_t <- t(expr)
    expr_imp <- VIM::kNN(as.data.frame(expr_t), k = 5, imp_var = FALSE)
    expr <- t(as.matrix(expr_imp))
  }
  
  list(data = as.data.frame(t(expr)), cells = colnames(expr), genes = rownames(expr))
}

# --- 5. Batch correction + PCA ---
process_and_plot_pca <- function(gdsc_data, ccle_data) {
  log_message("Starting batch effect correction using ComBat")
  
  common_genes <- intersect(colnames(gdsc_data$data), colnames(ccle_data$data))
  log_message(sprintf("Common genes: %d", length(common_genes)))
  
  gdsc_mat <- gdsc_data$data[, common_genes, drop = FALSE]
  ccle_mat <- ccle_data$data[, common_genes, drop = FALSE]
  
  # Combine for ComBat
  combined_data <- rbind(t(gdsc_mat), t(ccle_mat))
  
  combined_data <- filter_features(combined_data)
  
  if (nrow(combined_data) < 100) {
    log_message("Too few genes after filtering. Skipping ComBat.")
    return(NULL)
  }
  
  # ComBat batch correction
  batch_labels <- c(rep("GDSC", nrow(gdsc_mat)), rep("CCLE", nrow(ccle_mat)))
  corrected_data <- ComBat(dat = t(combined_data), batch = batch_labels, mod = NULL)
  
  # Normalize
  corrected_data <- normalize_data(corrected_data)
  
  # Plot PCA before and after
  plot_pca(t(combined_data), t(corrected_data), nrow(gdsc_mat), nrow(ccle_mat))
  
  # Save processed matrices (key output files)
  write.csv(t(combined_data), "before_combat_GDSC_CCLE.csv", row.names = TRUE)
  write.csv(t(corrected_data), "after_combat_GDSC_CCLE.csv", row.names = TRUE)
  
  log_message("Batch correction and PCA completed. Files saved.")
}

plot_pca <- function(before_data, after_data, gdsc_n, ccle_n) {
  log_message("Generating PCA plots")
  
  pca_before <- prcomp(before_data, scale. = TRUE)
  pca_after  <- prcomp(after_data,  scale. = TRUE)
  
  var_before <- summary(pca_before)$importance[2, 1:2] * 100
  var_after  <- summary(pca_after)$importance[2, 1:2] * 100
  
  df_before <- data.frame(PC1 = pca_before$x[,1], PC2 = pca_before$x[,2],
                          Batch = c(rep("GDSC", gdsc_n), rep("CCLE", ccle_n)))
  df_after  <- data.frame(PC1 = pca_after$x[,1],  PC2 = pca_after$x[,2],
                          Batch = c(rep("GDSC", gdsc_n), rep("CCLE", ccle_n)))
  
  p_before <- ggplot(df_before, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 1.8) + theme_minimal() +
    labs(title = "PCA Before ComBat", 
         x = sprintf("PC1 (%.1f%%)", var_before[1]),
         y = sprintf("PC2 (%.1f%%)", var_before[2]))
  
  p_after <- ggplot(df_after, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 1.8) + theme_minimal() +
    labs(title = "PCA After ComBat", 
         x = sprintf("PC1 (%.1f%%)", var_after[1]),
         y = sprintf("PC2 (%.1f%%)", var_after[2]))
  
  ggsave("PCA_Before_ComBat_GDSC_CCLE.png", p_before, width = 10, height = 8, dpi = 600)
  ggsave("PCA_After_ComBat_GDSC_CCLE.png",  p_after,  width = 10, height = 8, dpi = 600)
  
  print(p_before)
  print(p_after)
}

# --- 6. Run preprocessing ---
gdsc_processed <- preprocess_dataset(GDSC, "GDSC")
ccle_processed <- preprocess_dataset(CCLE, "CCLE")

if (!is.null(gdsc_processed) && !is.null(ccle_processed)) {
  process_and_plot_pca(gdsc_processed, ccle_processed)
}

log_message("Preprocessing, batch correction, and PCA plotting completed successfully")
