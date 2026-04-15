# ==============================
# Biomarker Validation using TCGA-OV
# ==============================

suppressPackageStartupMessages({
  library(TCGAbiolinks)
  library(SummarizedExperiment)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
})

# ------------------------------
# STEP 1: Download TCGA-OV data
# ------------------------------
query <- GDCquery(
  project = "TCGA-OV",
  data.category = "Transcriptome Profiling",
  data.type = "Gene Expression Quantification",
  workflow.type = "STAR - Counts"
)

GDCdownload(query)
data <- GDCprepare(query)

# ------------------------------
# STEP 2: Expression matrix
# ------------------------------
expr <- assay(data)
expr <- log2(expr + 1)

# Clean ENSEMBL IDs
rownames(expr) <- gsub("\\..*", "", rownames(expr))

# ------------------------------
# STEP 3: Map genes
# ------------------------------
gene_symbols <- mapIds(
  org.Hs.eg.db,
  keys = rownames(expr),
  column = "SYMBOL",
  keytype = "ENSEMBL",
  multiVals = "first"
)

rownames(expr) <- gene_symbols
expr <- expr[!is.na(rownames(expr)), ]

# ------------------------------
# STEP 4: Extract biomarkers
# ------------------------------
genes <- c("CCR10","PLEKHH2","WNT9B","ITPRID1","CHRNG","DIRAS3")

expr_subset <- expr[rownames(expr) %in% genes, ]
expr_subset <- t(expr_subset)

# ------------------------------
# STEP 5: Boxplot
# ------------------------------
png("plots/biomarker_expression.png", width=4000, height=2500, res=600)

par(mar=c(8,5,4,2))
boxplot(expr_subset,
        main="Biomarker Expression in TCGA-OV",
        las=2,
        col="lightblue",
        outline=FALSE)

dev.off()

# ------------------------------
# STEP 6: Correlation heatmap
# ------------------------------
cor_mat <- cor(expr_subset)

png("plots/biomarker_correlation.png", width=4000, height=3000, res=600)

heatmap(cor_mat,
        col=colorRampPalette(c("navy","white","red"))(100),
        scale="none",
        margins=c(8,8),
        main="Biomarker Correlation")

dev.off()