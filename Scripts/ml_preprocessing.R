# =============================================================================
# Data Retrieval, Pre-processing, Global Batch Correction & PCA
# First correct whole data for PCA, then per-drug split
# =============================================================================

# --- 1. Install and load packages ---
pkgs <- c("PharmacoGx", "dplyr", "ggplot2", "VIM", "sva", "biomaRt")
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

# Logging
log_message <- function(message, file = "logs/preprocess.log") {
  dir.create("logs", showWarnings = FALSE)
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, message), file = file, append = TRUE)
}

log_message("Starting global preprocessing, batch correction and PCA")

set.seed(123)

# --- 2. Load raw data ---
log_message("Loading raw PharmacoSets")
CCLE <- readRDS("CCLE.rds")
GDSC <- readRDS("GDSC2.rds")

CCLE <- updateObject(CCLE)
GDSC <- updateObject(GDSC)

log_message("Datasets loaded successfully")

# Ranked drugs
ranked_drugs <- c("Sorafenib", "Erlotinib", "Topotecan", "Paclitaxel", "Lapatinib", "Selumetinib")

# --- 3. Helper functions ---
validate_genes <- function(genes, name) {
  log_message(sprintf("Validating genes for %s", name))
  tryCatch({
    mart <- useMart("ensembl", dataset = "hsapiens_gene_ensembl")
    mapping <- getBM(attributes = c("ensembl_gene_id", "hgnc_symbol"),
                     filters = "ensembl_gene_id", values = head(genes, 10), mart = mart)
    log_message(sprintf("Validated %d/10 genes for %s", sum(mapping$hgnc_symbol != ""), name))
  }, error = function(e) log_message("Gene validation skipped"))
}

normalize_data <- function(expr) {
  expr <- expr + 1
  log2(expr) %>% t() %>% scale() %>% t()
}

filter_features <- function(data) {
  log_message("MAD-based filtering on combined data")
  initial <- nrow(data)
  mad_vals <- apply(data, 1, function(x) {
    if (sum(!is.na(x)) < 2) return(NA)
    med <- median(x, na.rm = TRUE)
    if (is.na(med)) return(NA)
    mad_val <- median(abs(x - med), na.rm = TRUE)
    if (is.na(mad_val) || mad_val == 0) return(NA)
    mad_val
  })
  valid <- !is.na(mad_vals)
  data <- data[valid, , drop = FALSE]
  mad_vals <- mad_vals[valid]
  threshold <- quantile(mad_vals, 0.05, na.rm = TRUE)
  data <- data[mad_vals > threshold, , drop = FALSE]
  log_message(sprintf("MAD filter: %d → %d genes", initial, nrow(data)))
  data
}

# Global preprocess (RNA only)
preprocess_global <- function(pset, name) {
  log_message(sprintf("Global RNA preprocessing for %s", name))
  se <- summarizeMolecularProfiles(pset, mDataType = "rna", fill.missing = FALSE)
  expr <- assay(se)
  log_message(sprintf("%s: %d genes × %d cells", name, nrow(expr), ncol(expr)))
  validate_genes(rownames(expr), name)
  
  if (any(is.na(expr))) {
    log_message(sprintf("KNN imputation for %s", name))
    expr_t <- t(expr)
    expr_imp <- VIM::kNN(as.data.frame(expr_t), k = 5, imp_var = FALSE)
    expr <- t(as.matrix(expr_imp))
  }
  list(data = as.data.frame(t(expr)), cells = colnames(expr), genes = rownames(expr))
}

# PCA plot function
plot_pca <- function(before, after, gdsc_n, ccle_n) {
  log_message("Generating PCA plots (global)")
  pca_b <- prcomp(before, scale. = TRUE)
  pca_a <- prcomp(after,  scale. = TRUE)
  
  var_b <- summary(pca_b)$importance[2, 1:2] * 100
  var_a <- summary(pca_a)$importance[2, 1:2] * 100
  
  df_b <- data.frame(PC1 = pca_b$x[,1], PC2 = pca_b$x[,2],
                     Batch = c(rep("GDSC", gdsc_n), rep("CCLE", ccle_n)))
  df_a <- data.frame(PC1 = pca_a$x[,1], PC2 = pca_a$x[,2],
                     Batch = c(rep("GDSC", gdsc_n), rep("CCLE", ccle_n)))
  
  p_b <- ggplot(df_b, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 1.8) + theme_minimal() +
    labs(title = "PCA Before ComBat (Global)", 
         x = sprintf("PC1 (%.1f%%)", var_b[1]),
         y = sprintf("PC2 (%.1f%%)", var_b[2]))
  
  p_a <- ggplot(df_a, aes(x = PC1, y = PC2, color = Batch)) +
    geom_point(size = 1.8) + theme_minimal() +
    labs(title = "PCA After ComBat (Global)", 
         x = sprintf("PC1 (%.1f%%)", var_a[1]),
         y = sprintf("PC2 (%.1f%%)", var_a[2]))
  
  ggsave("PCA_Before_ComBat_GDSC_CCLE.png", p_b, width = 10, height = 8, dpi = 600)
  ggsave("PCA_After_ComBat_GDSC_CCLE.png",  p_a, width = 10, height = 8, dpi = 600)
  print(p_b); print(p_a)
}

# --- 4. Global batch correction for PCA ---
gdsc_global <- preprocess_global(GDSC, "GDSC")
ccle_global <- preprocess_global(CCLE, "CCLE")

if (!is.null(gdsc_global) && !is.null(ccle_global)) {
  common_genes <- intersect(colnames(gdsc_global$data), colnames(ccle_global$data))
  log_message(sprintf("Global common genes: %d", length(common_genes)))
  
  gdsc_mat <- gdsc_global$data[, common_genes, drop = FALSE]
  ccle_mat <- ccle_global$data[, common_genes, drop = FALSE]
  
  combined <- rbind(t(gdsc_mat), t(ccle_mat))
  combined <- filter_features(combined)
  
  if (nrow(combined) >= 100) {
    batch_lab <- c(rep("GDSC", nrow(gdsc_mat)), rep("CCLE", nrow(ccle_mat)))
    corrected <- ComBat(dat = t(combined), batch = batch_lab, mod = NULL)
    corrected <- normalize_data(corrected)
    
    plot_pca(t(combined), t(corrected), nrow(gdsc_mat), nrow(ccle_mat))
    
    # Save global corrected matrices
    write.csv(t(combined),  "before_combat_GDSC_CCLE.csv", row.names = TRUE)
    write.csv(t(corrected), "after_combat_GDSC_CCLE.csv", row.names = TRUE)
    
    log_message("Global batch correction and PCA completed")
  }
}

# --- 5. Per-drug split after global correction (train = GDSC, valid = CCLE) ---
preprocessed_data <- list()
summary_stats <- data.frame(Drug = character(), Dataset = character(),
                            Genes = integer(), Cells = integer(),
                            FilteredGenes = integer(), CommonGenes = integer(),
                            stringsAsFactors = FALSE)

for (drug in ranked_drugs) {
  log_message(sprintf("=== Per-drug processing for %s ===", drug))
  
  train_result <- preprocess_drug(drug, GDSC, "GDSC")   # GDSC as train
  valid_result <- preprocess_drug(drug, CCLE, "CCLE")   # CCLE as valid
  
  if (is.null(train_result) || is.null(valid_result)) {
    log_message(sprintf("Skipping %s", drug))
    next
  }
  
  common_genes <- intersect(train_result$genes, valid_result$genes)
  if (length(common_genes) < 100) {
    log_message(sprintf("Skipping %s: too few common genes", drug))
    next
  }
  
  # Use global correction logic but per-drug response
  corrected_data <- correct_batch_effects(
    train_result$data, valid_result$data,
    train_result$drug_response, valid_result$drug_response,
    drug, common_genes
  )
  
  if (is.null(corrected_data)) next
  
  # Save per-drug files (exact names kept)
  saveRDS(corrected_data$train, paste0("preprocessed_", drug, "_train.rds"))
  saveRDS(corrected_data$valid, paste0("preprocessed_", drug, "_valid.rds"))
  
  write.csv(corrected_data$train, paste0("preprocessed_", drug, "_train.csv"), row.names = TRUE)
  write.csv(corrected_data$valid, paste0("preprocessed_", drug, "_valid.csv"), row.names = TRUE)
  
  summary_stats <- rbind(summary_stats,
    data.frame(Drug = drug, Dataset = "GDSC", Genes = nrow(train_result$data),
               Cells = length(train_result$cells), FilteredGenes = ncol(corrected_data$train)-1,
               CommonGenes = length(common_genes)),
    data.frame(Drug = drug, Dataset = "CCLE", Genes = nrow(valid_result$data),
               Cells = length(valid_result$cells), FilteredGenes = ncol(corrected_data$valid)-1,
               CommonGenes = length(common_genes))
  )
}

write.csv(summary_stats, "results/tables/preprocessing_summary.csv", row.names = FALSE)
saveRDS(summary_stats, "results/tables/preprocessing_summary.rds")

