# ==============================
# Feature Selection using mRMR
# ==============================

suppressPackageStartupMessages({
  library(mRMRe)
})

set.seed(123)

drugs <- c("Sorafenib","Erlotinib","Topotecan",
           "Paclitaxel","Lapatinib","Selumetinib")

selected_features <- list()

for (drug in drugs) {
  
  train_data <- readRDS(paste0("results/preprocessed_", drug, "_train.rds"))
  
  features <- train_data[, setdiff(colnames(train_data), "drug_response")]
  target <- train_data$drug_response
  
  data_df <- cbind(features, drug_response=target)
  
  # mRMR feature selection
  mrmr_data <- mRMR.data(data=data_df)
  
  mrmr_result <- mRMR.classic(
    data=mrmr_data,
    target_indices=ncol(data_df),
    feature_count=100
  )
  
  selected <- solutions(mrmr_result)[[1]]
  
  selected_features[[drug]] <- colnames(features)[selected]
}

# Save results
saveRDS(selected_features, "results/selected_features.rds")