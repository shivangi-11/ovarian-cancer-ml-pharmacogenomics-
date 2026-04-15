# ==============================
# Data Preprocessing & PCA Analysis
# ==============================

suppressPackageStartupMessages({
  library(PharmacoGx)
  library(sva)
  library(ggplot2)
  library(VIM)
})

set.seed(123)

# ------------------------------
# Load datasets (user-defined paths)
# ------------------------------
CCLE <- readRDS("path_to_CCLE.rds")
GDSC <- readRDS("path_to_GDSC.rds")

CCLE <- updateObject(CCLE)
GDSC <- updateObject(GDSC)

# ------------------------------
# Extract RNA expression
# ------------------------------
extract_expression <- function(dataset) {
  se <- summarizeMolecularProfiles(dataset, mDataType="rna")
  expr <- assay(se)
  
  # Log transform + scaling
  expr <- log2(expr + 1)
  expr <- t(scale(t(expr)))
  
  return(expr)
}

expr_ccle <- extract_expression(CCLE)
expr_gdsc <- extract_expression(GDSC)

# ------------------------------
# Common genes
# ------------------------------
common_genes <- intersect(rownames(expr_ccle), rownames(expr_gdsc))

expr_ccle <- expr_ccle[common_genes, ]
expr_gdsc <- expr_gdsc[common_genes, ]

# ------------------------------
# Combine datasets
# ------------------------------
combined <- cbind(expr_gdsc, expr_ccle)

# ------------------------------
# Batch correction (ComBat)
# ------------------------------
batch <- c(rep("GDSC", ncol(expr_gdsc)),
           rep("CCLE", ncol(expr_ccle)))

combat_data <- ComBat(dat=combined, batch=batch)

# ------------------------------
# PCA before and after
# ------------------------------
plot_pca <- function(data, title, filename) {
  pca <- prcomp(t(data), scale.=TRUE)
  
  df <- data.frame(PC1=pca$x[,1], PC2=pca$x[,2], Batch=batch)
  
  p <- ggplot(df, aes(PC1, PC2, color=Batch)) +
    geom_point(size=2) +
    theme_minimal() +
    ggtitle(title)
  
  ggsave(filename, p, width=8, height=6, dpi=600)
}

plot_pca(combined, "Before ComBat", "plots/pca_before.png")
plot_pca(combat_data, "After ComBat", "plots/pca_after.png")

# Save processed data
saveRDS(combat_data, "results/processed_expression.rds")
