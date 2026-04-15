# ============================================
# Biomarker Identification using PharmacoGx
# ============================================

suppressPackageStartupMessages({
  library(PharmacoGx)
  library(SummarizedExperiment)
  library(dplyr)
})

set.seed(123)

# ==============================
# STEP 1: Load processed dataset
# ==============================
data <- readRDS("data/common_pharmaco_datasets.rds")

# Extract datasets
GDSC <- data$GDSC
CCLE <- data$CCLE

# ==============================
# STEP 2: Define drugs of interest
# ==============================
drugs <- c("Sorafenib","Erlotinib","Topotecan",
           "Paclitaxel","Lapatinib","Selumetinib")

results_list <- list()

# ==============================
# STEP 3: Biomarker analysis
# ==============================
for (dataset_name in c("GDSC","CCLE")) {
  
  cat("Processing:", dataset_name, "\n")
  
  pset <- if (dataset_name == "GDSC") GDSC else CCLE
  
  # Run drug sensitivity association
  res <- drugSensitivitySig(
    pset,
    mDataType = "rna",
    sensitivity.measure = "aac_recomputed",
    drugs = drugs,
    verbose = FALSE
  )
  
  # Convert to dataframe
  df <- as.data.frame(res)
  
  # Filter significant biomarkers
  df_sig <- df %>%
    filter(fdr < 0.05)
  
  results_list[[dataset_name]] <- df_sig
}

# ==============================
# STEP 4: Combine results
# ==============================
biomarkers <- bind_rows(results_list, .id = "Dataset")

# ==============================
# STEP 5: Save results
# ==============================
saveRDS(biomarkers, "results/biomarkers.rds")
write.csv(biomarkers, "results/biomarkers.csv", row.names = FALSE)

cat("Biomarker identification completed.\n")