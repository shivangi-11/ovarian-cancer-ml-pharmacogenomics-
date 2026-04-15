# ============================================
# Drug Sensitivity Data Extraction (AAC & IC50)
# ============================================

suppressPackageStartupMessages({
  library(PharmacoGx)
})

set.seed(123)

# ------------------------------
# Load datasets
# ------------------------------
CCLE <- readRDS("data/CCLE.rds")
GDSC <- readRDS("data/GDSC.rds")

CCLE <- updateObject(CCLE)
GDSC <- updateObject(GDSC)

# ------------------------------
# Find common cell lines & drugs
# ------------------------------
common_pset <- intersectPSet(
  pSets = list(CCLE = CCLE, GDSC = GDSC),
  intersectOn = c("cell.lines", "drugs"),
  strictIntersect = TRUE
)

# ------------------------------
# Subset ovarian cell lines
# ------------------------------
cell_info <- cellInfo(common_pset$GDSC)

ovary_cells <- rownames(cell_info)[
  grepl("ovary", cell_info$tissueid, ignore.case = TRUE)
]

GDSC_OVARY <- subsetTo(common_pset$GDSC, cells = ovary_cells)
CCLE_OVARY <- subsetTo(common_pset$CCLE, cells = ovary_cells)

# ------------------------------
# Extract AAC & IC50
# ------------------------------
GDSC_aac <- summarizeSensitivityProfiles(
  GDSC_OVARY,
  sensitivity.measure = "aac_recomputed"
)

CCLE_aac <- summarizeSensitivityProfiles(
  CCLE_OVARY,
  sensitivity.measure = "aac_recomputed"
)

GDSC_ic50 <- summarizeSensitivityProfiles(
  GDSC_OVARY,
  sensitivity.measure = "ic50_recomputed"
)

CCLE_ic50 <- summarizeSensitivityProfiles(
  CCLE_OVARY,
  sensitivity.measure = "ic50_recomputed"
)

# ------------------------------
# Save outputs
# ------------------------------
write.csv(GDSC_aac, "GDSC_filtered_AAC_ovary.csv")
write.csv(CCLE_aac, "CCLE_filtered_AAC_ovary.csv")

write.csv(GDSC_ic50, "GDSC_filtered_ic50_ovary.csv")
write.csv(CCLE_ic50, "CCLE_filtered_ic50_ovary.csv")

# =============================================================================
# Consistency across GDSC and CCLE (CI / mCI) + AAC scatterplots with correlations
# =============================================================================

# --- 1. Install & load required packages (run once) ---
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")

# Install bhklab/mci (required for paired.concordance.index)
if (!requireNamespace("mCI", quietly = TRUE)) {
  devtools::install_github("bhklab/mci", ref = "light", quiet = TRUE)
}

# CRAN packages
pkgs <- c("reshape2", "ggplot2", "dplyr", "patchwork")
invisible(lapply(pkgs, function(p) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
  library(p, character.only = TRUE)
}))

# Set seed for reproducibility
set.seed(123)


# --- 2. Load pre-processed AAC data (from GitHub data) ---
cat("Loading pre-processed AAC matrices...\n")
GDSC.aac <- as.matrix(read.csv("data/GDSC_filtered_AAC_ovary.csv", 
                               row.names = 1, check.names = FALSE))
CCLE.aac <- as.matrix(read.csv("data/processed/CCLE_filtered_AAC_ovary.csv", 
                               row.names = 1, check.names = FALSE))

# Keep only common drugs
common_drugs <- intersect(rownames(GDSC.aac), rownames(CCLE.aac))
GDSC.aac <- GDSC.aac[common_drugs, , drop = FALSE]
CCLE.aac <- CCLE.aac[common_drugs, , drop = FALSE]
drugs <- common_drugs
cat("Using", length(drugs), "common drugs and", ncol(GDSC.aac), "ovarian cell lines.\n")

# --- 3. Consistency Analysis (2.1.2) - CI and mCI ---
cat("Computing Concordance Index (CI) and modified Concordance Index (mCI)...\n")

c_index <- mc_index <- NULL
for (drug in drugs) {
  # Standard CI (delta = 0)
  tt <- mCI::paired.concordance.index(GDSC.aac[drug, ], CCLE.aac[drug, ],
                                      delta.pred = 0, delta.obs = 0, 
                                      alternative = "greater")
  c_index <- c(c_index, tt$cindex)
  
  # mCI (delta = 0.2, logic "or")
  tt <- mCI::paired.concordance.index(GDSC.aac[drug, ], CCLE.aac[drug, ],
                                      delta.pred = 0.2, delta.obs = 0.2, 
                                      alternative = "greater", logic.operator = "or")
  mc_index <- c(mc_index, tt$cindex)
}

names(c_index) <- names(mc_index) <- drugs

# Create ranking
ranking_df <- data.frame(
  Drug = drugs,
  CI    = c_index,
  mCI   = mc_index
)
ranking_df <- ranking_df[order(ranking_df$mCI, decreasing = TRUE), ]

# Top drugs (mCI > 0.7 as per reference line in original analysis)
top_drugs <- ranking_df[ranking_df$mCI > 0.7, ]

# Save results with EXACT original file names
write.csv(ranking_df, file = "all_drugs_ranking.csv",   row.names = FALSE)


# --- 4. Modelling Sensitivity Data - AAC scatterplots + correlations ---

# Long format for merging
GDSC_long <- as.data.frame(as.matrix(GDSC.aac))
GDSC_long$Drug <- rownames(GDSC.aac)
GDSC_long <- reshape2::melt(GDSC_long, id.vars = "Drug", 
                            variable.name = "Cell_Line", value.name = "GDSC")

CCLE_long <- as.data.frame(as.matrix(CCLE.aac))
CCLE_long$Drug <- rownames(CCLE.aac)
CCLE_long <- reshape2::melt(CCLE_long, id.vars = "Drug", 
                            variable.name = "Cell_Line", value.name = "CCLE")

merged_data <- merge(GDSC_long, CCLE_long, by = c("Cell_Line", "Drug"))

write.csv(merged_data, file = "aac_data_common_merged.csv", row.names = TRUE)

# Remove NAs and save cleaned version
merged_data_cleaned <- merged_data[complete.cases(merged_data), ]
write.csv(merged_data_cleaned, file = "aac_data_common_merged_removed_na.csv", row.names = TRUE)

# Correlations per drug
correlation_results <- merged_data_cleaned %>%
  dplyr::group_by(Drug) %>%
  dplyr::summarise(
    Pearson = cor(GDSC, CCLE, method = "pearson", use = "complete.obs"),
    Spearman = cor(GDSC, CCLE, method = "spearman", use = "complete.obs"),
    n = dplyr::n()
  ) %>%
  dplyr::filter(n > 1)

# Global axis limits for consistent scaling across plots
x_limits <- range(merged_data_cleaned$GDSC, na.rm = TRUE)
y_limits <- range(merged_data_cleaned$CCLE, na.rm = TRUE)

# Per-drug scatterplots (Figure 2 style)
plot_list <- list()
for (drug in unique(correlation_results$Drug)) {
  subset_data <- merged_data_cleaned %>% dplyr::filter(Drug == drug)
  if (nrow(subset_data) < 2) next
  
  pearson_corr  <- cor(subset_data$GDSC, subset_data$CCLE, method = "pearson",  use = "complete.obs")
  spearman_corr <- cor(subset_data$GDSC, subset_data$CCLE, method = "spearman", use = "complete.obs")
  
  p <- ggplot(subset_data, aes(x = GDSC, y = CCLE)) +
    geom_point(color = "#1F77B4", size = 1.5, alpha = 0.7) +
    geom_smooth(method = "lm", color = "#D62728", se = FALSE, size = 0.8) +
    scale_x_continuous(limits = x_limits) +
    scale_y_continuous(limits = y_limits) +
    theme_minimal(base_size = 14) +
    labs(title = drug,
         subtitle = paste("Pearson:", round(pearson_corr, 2), 
                         "| Spearman:", round(spearman_corr, 2)),
         x = "GDSC (AAC)", y = "CCLE (AAC)") +
    theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 12),
          plot.subtitle = element_text(hjust = 0.5, size = 10, color = "grey30"),
          axis.title = element_text(size = 12, face = "bold"))
  
  plot_list[[drug]] <- p
}

# Combine into grid and save
if (length(plot_list) > 0) {
  plot_grid <- patchwork::wrap_plots(plot_list, ncol = 3) +
    patchwork::plot_annotation
  
  ggsave("scatterplot_correlation.jpg", plot = plot_grid, 
         width = 15, height = 12, dpi = 600, bg = "white")
  print(plot_grid)
  cat("Saved: scatterplot_correlation.jpg\n")
}
