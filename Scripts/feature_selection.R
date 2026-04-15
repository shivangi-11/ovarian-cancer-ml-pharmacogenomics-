# =============================================================================
Feature Selection using mRMR
# =============================================================================

suppressPackageStartupMessages({
  library(mRMRe)
  library(dplyr)
  library(biomaRt)
})

# Logging function
log_message <- function(message, file = "logs/feature_selection.log") {
  dir.create("logs", showWarnings = FALSE)
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, message), file = file, append = TRUE)
}

log_message("Starting feature selection with mRMR")

set.seed(123)

# Ranked drugs (consistent with preprocessing)
ranked_drugs <- c("Sorafenib", "Erlotinib", "Topotecan", "Paclitaxel", "Lapatinib", "Selumetinib")

# Storage
selected_features <- list()
top_20_features   <- list()

for (drug in ranked_drugs) {
  log_message(sprintf("Processing feature selection for: %s", drug))
  
  # Load preprocessed data from Script 01
  train_data <- tryCatch({
    readRDS(paste0("/preprocessed_", drug, "_train.rds"))
  }, error = function(e) {
    log_message(sprintf("Failed to load train data for %s: %s", drug, conditionMessage(e)))
    return(NULL)
  })
  
  valid_data <- tryCatch({
    readRDS(paste0("preprocessed_", drug, "_valid.rds"))
  }, error = function(e) {
    log_message(sprintf("Failed to load valid data for %s: %s", drug, conditionMessage(e)))
    return(NULL)
  })
  
  if (is.null(train_data) || is.null(valid_data)) {
    log_message(sprintf("Skipping %s - data loading failed", drug))
    next
  }
  
  log_message(sprintf("Train: %d samples, %d features | Valid: %d samples, %d features", 
                      nrow(train_data), ncol(train_data)-1, nrow(valid_data), ncol(valid_data)-1))
  
  # Prepare data for mRMR (features + response from train set)
  features <- train_data[, setdiff(colnames(train_data), "drug_response"), drop = FALSE]
  data_df  <- cbind(features, drug_response = train_data$drug_response)
  
  # Median imputation for missing values
  na_count <- sum(is.na(data_df))
  if (na_count > 0) {
    log_message(sprintf("Imputing %d missing values with median for %s", na_count, drug))
    for (col in colnames(data_df)) {
      if (any(is.na(data_df[[col]]))) {
        data_df[[col]][is.na(data_df[[col]])] <- median(data_df[[col]], na.rm = TRUE)
      }
    }
  }
  
  # Ensure all numeric
  if (any(!sapply(data_df, is.numeric))) {
    log_message(sprintf("Non-numeric columns detected in %s. Skipping.", drug))
    next
  }
  
  # Run mRMR
  mrmr_data <- tryCatch(mRMR.data(data = data_df), 
                        error = function(e) {
                          log_message(sprintf("mRMR data prep error for %s: %s", drug, conditionMessage(e)))
                          return(NULL)
                        })
  
  if (is.null(mrmr_data)) {
    log_message(sprintf("Skipping %s: mRMR data preparation failed", drug))
    next
  }
  
  mrmr_result <- tryCatch({
    mRMR.classic(data = mrmr_data, target_indices = ncol(data_df), feature_count = 100)
  }, error = function(e) {
    log_message(sprintf("mRMR failed for %s: %s", drug, conditionMessage(e)))
    return(NULL)
  })
  
  if (is.null(mrmr_result)) {
    log_message(sprintf("Skipping %s: mRMR execution failed", drug))
    next
  }
  
  # Extract selected features
  selected_indices <- solutions(mrmr_result)[[1]]
  selected_features[[drug]] <- colnames(features)[selected_indices]
  
  log_message(sprintf("Selected %d features for %s", length(selected_features[[drug]]), drug))
  
  # Subset train and valid to selected features
  train_selected <- train_data[, c(selected_features[[drug]], "drug_response"), drop = FALSE]
  valid_selected <- valid_data[, c(selected_features[[drug]], "drug_response"), drop = FALSE]
  
  saveRDS(train_selected, paste0("/preprocessed_", drug, "_train_selected.rds"))
  saveRDS(valid_selected, paste0("preprocessed_", drug, "_valid_selected.rds"))
  
  write.csv(train_selected, paste0("preprocessed_", drug, "_train_selected.csv"), row.names = TRUE)
  write.csv(valid_selected, paste0("preprocessed_", drug, "_valid_selected.csv"), row.names = TRUE)

}

# Save overall lists
saveRDS(selected_features, "selected_features.rds")

log_message("Feature selection completed successfully")

