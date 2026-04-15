# =============================================================================
#  Biomarker Identification using drugSensitivitySig
# =============================================================================

suppressPackageStartupMessages({
  library(PharmacoGx)
  library(SummarizedExperiment)
  library(dplyr)
  library(ggplot2)
  library(ggrepel)
  library(ggsci)
})

# Logging
log_message <- function(message, file = "logs/biomarker_identification.log") {
  dir.create("logs", showWarnings = FALSE)
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, message), file = file, append = TRUE)
}

log_message("Starting biomarker identification")

set.seed(123)

# --- 1. Load pre-processed common datasets (from previous scripts) ---
input_file <- "common_pharmaco_datasets.rds"
if (!file.exists(input_file)) {
  stop("Input file missing: ", input_file, ". Run previous preprocessing scripts first.")
}

log_message("Loading common PharmacoSet objects")
common <- readRDS(input_file)

# Drugs of interest
drugs <- c("Sorafenib", "Erlotinib", "Topotecan", "Paclitaxel", "Lapatinib", "Selumetinib")

# Clean publication theme
theme_pub <- theme_minimal(base_size = 14) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5),
    plot.subtitle = element_text(hjust = 0.5),
    axis.text.x = element_text(angle = 45, hjust = 1),
    strip.text = element_text(face = "bold")
  )

dir.create("results", showWarnings = FALSE, recursive = TRUE)

# --- 2. Run drugSensitivitySig for biomarker discovery ---
log_message("Running drugSensitivitySig modeling")
sensitivity_sigs <- list()

for (dataset in c("GDSC", "CCLE")) {
  log_message(paste("Processing", dataset))
  pset <- common[[dataset]]
  
  # Fix RNA profile metadata
  if ("rna" %in% names(molecularProfilesSlot(pset))) {
    mp <- molecularProfilesSlot(pset)$rna
    assayNames(mp) <- "exprs"
    metadata(mp)$annotation <- "rna"
    molecularProfilesSlot(pset)$rna <- mp
  } else {
    log_message(paste("RNA profile missing in", dataset))
    next
  }
  
  # Get valid drugs
  available_drugs <- drugNames(pset)
  valid_drugs <- character(0)
  for (d in drugs) {
    match <- available_drugs[grepl(tolower(d), tolower(available_drugs), fixed = TRUE)]
    if (length(match) > 0) valid_drugs <- c(valid_drugs, match[1])
  }
  
  if (length(valid_drugs) == 0) {
    log_message(paste("No valid drugs found in", dataset))
    next
  }
  
  # Run drugSensitivitySig
  tryCatch({
    sig <- drugSensitivitySig(
      object = pset,
      mDataType = "rna",
      drugs = valid_drugs,
      sensitivity.measure = "aac_recomputed",
      molecular.summary.stat = "median",
      sensitivity.summary.stat = "median",
      verbose = FALSE
    )
    sensitivity_sigs[[dataset]] <- sig
    log_message(paste("drugSensitivitySig completed for", dataset))
  }, error = function(e) {
    log_message(paste("drugSensitivitySig failed for", dataset, ":", conditionMessage(e)))
  })
}

# --- 3. Process and save results ---
sig_results <- data.frame()

for (dataset in names(sensitivity_sigs)) {
  if (is.null(sensitivity_sigs[[dataset]])) next
  sig <- sensitivity_sigs[[dataset]]
  
  for (drug in dimnames(sig)[[2]]) {
    for (gene in dimnames(sig)[[1]]) {
      pval <- sig[gene, drug, "pvalue"]
      fdr  <- sig[gene, drug, "fdr"]
      est  <- sig[gene, drug, "estimate"]
      
      if (!is.na(pval)) {
        sig_results <- rbind(sig_results, data.frame(
          Dataset = dataset,
          Drug    = drug,
          Gene    = gene,
          Estimate = est,
          PValue   = pval,
          FDR      = fdr
        ))
      }
    }
  }
}

# Save all and significant results
saveRDS(sig_results, "gene_drug_sensitivity_all.rds")
sig_filtered <- sig_results %>% filter(FDR < 0.05)
saveRDS(sig_filtered, "gene_drug_sensitivity_signatures.rds")

write.csv(sig_results, "gene_drug_sensitivity_all.csv", row.names = FALSE)
write.csv(sig_filtered, "/gene_drug_sensitivity_signatures.csv", row.names = FALSE)

log_message(paste("Saved", nrow(sig_results), "total associations,", nrow(sig_filtered), "significant"))

# --- 4. Visualizations ---
if (nrow(sig_results) > 0) {
  sig_results$logP <- -log10(sig_results$PValue)
  sig_results$Significance <- ifelse(sig_results$FDR < 0.05, "Significant (FDR < 0.05)", "Not Significant")
  
  # Top genes for labeling
  top_genes <- sig_results %>%
    group_by(Drug, Dataset) %>%
    arrange(FDR) %>%
    slice_head(n = 5) %>%
    ungroup() 
  
  # Bar plot of significant biomarkers
  if (nrow(sig_filtered) > 0) {
    bar_plot <- ggplot(sig_filtered, aes(x = reorder(Gene, -Estimate), y = Estimate, fill = Dataset)) +
      geom_bar(stat = "identity", position = "dodge") +
      facet_wrap(~ Drug, scales = "free_x", ncol = 3) +
      scale_fill_manual(values = pal_npg()(2)) +
      labs(title = "Significant Biomarkers (FDR < 0.05)",
           x = "Gene", y = "Effect Size", fill = "Dataset") +
      theme_pub +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))
    
    ggsave("results/plots/biomarker_bar_plot.pdf", bar_plot, width = 14, height = 9, dpi = 400)
  }
  
  log_message("Biomarker visualizations saved")
}

log_message("Biomarker identification completed successfully")
