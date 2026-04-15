# ==============================
# Pharmacogenomics Data Processing & Concordance Analysis
# ==============================

# Load required libraries
suppressPackageStartupMessages({
  library(PharmacoGx)
  library(dplyr)
  library(ggplot2)
})

# Set seed for reproducibility
set.seed(123)

# ==============================
# STEP 1: Load datasets
# ==============================
# Note: Users should provide local paths to PharmacoSet objects

ccle <- readRDS("path_to_CCLE.rds")
gdsc <- readRDS("path_to_GDSC.rds")

# Ensure compatibility
ccle <- updateObject(ccle)
gdsc <- updateObject(gdsc)

# ==============================
# STEP 2: Identify common datasets
# ==============================

common_pset <- intersectPSet(
  pSets = list(CCLE = ccle, GDSC = gdsc),
  intersectOn = c("cell.lines", "drugs"),
  strictIntersect = TRUE
)

# ==============================
# STEP 3: Subset ovarian cancer cell lines
# ==============================

cell_info <- cellInfo(common_pset$GDSC)

ovary_cells <- rownames(cell_info)[
  grepl("Ovary", cell_info$tissueid, ignore.case = TRUE)
]

GDSC_OVARY <- subsetTo(common_pset$GDSC, cells = ovary_cells)
CCLE_OVARY <- subsetTo(common_pset$CCLE, cells = ovary_cells)

# ==============================
# STEP 4: Compute drug sensitivity (AAC)
# ==============================

GDSC.aac <- summarizeSensitivityProfiles(
  GDSC_OVARY,
  sensitivity.measure = "aac_recomputed",
  summary.stat = "median"
)

CCLE.aac <- summarizeSensitivityProfiles(
  CCLE_OVARY,
  sensitivity.measure = "aac_recomputed",
  summary.stat = "median"
)

# ==============================
# STEP 5: Merge datasets for comparison
# ==============================

common_drugs <- intersect(rownames(GDSC.aac), rownames(CCLE.aac))

GDSC.aac <- GDSC.aac[common_drugs, ]
CCLE.aac <- CCLE.aac[common_drugs, ]

# ==============================
# STEP 6: Correlation analysis
# ==============================

cor_results <- data.frame()

for (drug in common_drugs) {
  
  x <- as.numeric(GDSC.aac[drug, ])
  y <- as.numeric(CCLE.aac[drug, ])
  
  valid <- complete.cases(x, y)
  
  if (sum(valid) > 2) {
    cor_results <- rbind(cor_results, data.frame(
      Drug = drug,
      Pearson = cor(x[valid], y[valid], method = "pearson"),
      Spearman = cor(x[valid], y[valid], method = "spearman")
    ))
  }
}

print(cor_results)

# ==============================
# STEP 7: Concordance Index (mCI)
# ==============================

# NOTE: Requires mCI package from bhklab

c_index <- mc_index <- NULL

for (drug in common_drugs) {
  
  c_val <- mCI::paired.concordance.index(
    GDSC.aac[drug, ], CCLE.aac[drug, ],
    delta.pred = 0, delta.obs = 0
  )$cindex
  
  mc_val <- mCI::paired.concordance.index(
    GDSC.aac[drug, ], CCLE.aac[drug, ],
    delta.pred = 0.2, delta.obs = 0.2,
    logic.operator = "or"
  )$cindex
  
  c_index <- c(c_index, c_val)
  mc_index <- c(mc_index, mc_val)
}

# Combine results
ranking_df <- data.frame(
  Drug = common_drugs,
  CI = c_index,
  mCI = mc_index
)

# Rank drugs
ranking_df <- ranking_df[order(ranking_df$mCI, decreasing = TRUE), ]

print(ranking_df)

# ==============================
# STEP 8: Save outputs
# ==============================

write.csv(ranking_df, "drug_concordance_results.csv", row.names = FALSE)

# ==============================
# END
# ==============================