# ==============================
# Machine Learning Models
# ==============================

suppressPackageStartupMessages({
  library(caret)
  library(glmnet)
})

set.seed(123)

drugs <- c("Sorafenib","Erlotinib","Topotecan",
           "Paclitaxel","Lapatinib","Selumetinib")

results <- data.frame()

for (drug in drugs) {
  
  train_data <- readRDS(paste0("results/preprocessed_", drug, "_train_selected.rds"))
  valid_data <- readRDS(paste0("results/preprocessed_", drug, "_valid_selected.rds"))
  
  ctrl <- trainControl(method="cv", number=5)
  
  model <- train(
    drug_response ~ .,
    data=train_data,
    method="glmnet",
    trControl=ctrl,
    tuneLength=10
  )
  
  pred <- predict(model, newdata=valid_data)
  
  pearson <- cor(pred, valid_data$drug_response)
  spearman <- cor(pred, valid_data$drug_response, method="spearman")
  
  results <- rbind(results, data.frame(
    Drug=drug,
    Pearson=pearson,
    Spearman=spearman
  ))
}

write.csv(results, "results/model_performance.csv", row.names=FALSE)