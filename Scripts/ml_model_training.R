# =============================================================================
#  Model Training
# Trains ElasticNet, Ridge, and Lasso on selected features
# Evaluates with Pearson & Spearman correlations
# Generates one combined performance plot
# =============================================================================

suppressPackageStartupMessages({
  library(caret)
  library(glmnet)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(wesanderson)
})

# Logging function
log_message <- function(message, file = "logs/train_models.log") {
  dir.create("logs", showWarnings = FALSE)
  timestamp <- format(Sys.time(), "%Y-%m-%d %H:%M:%S")
  cat(sprintf("[%s] %s\n", timestamp, message), file = file, append = TRUE)
}

log_message("Starting model training")

set.seed(123)

# Ranked drugs (consistent with previous scripts)
ranked_drugs <- c("Sorafenib", "Erlotinib", "Topotecan", "Paclitaxel", "Lapatinib", "Selumetinib")

# Storage
model_metrics <- data.frame(
  Drug = character(), Model = character(), 
  Pearson = numeric(), Spearman = numeric(), 
  stringsAsFactors = FALSE
)


# Pastel palette for plots
pastel_palette <- wes_palette("Zissou1", n = 2, type = "continuous")

# Clean plot theme
plot_theme <- theme_minimal(base_size = 12) +
  theme(
    plot.title = element_text(face = "bold", hjust = 0.5, size = 14),
    axis.text.x = element_text(angle = 45, hjust = 1),
    legend.position = "top",
    strip.text = element_text(face = "bold")
  )

dir.create("results/", showWarnings = FALSE, recursive = TRUE)


for (drug in ranked_drugs) {
  log_message(sprintf("Training models for %s", drug))
  
  # Load selected feature data
  train_data <- tryCatch({
    readRDS(paste0("preprocessed_", drug, "_train_selected.rds"))
  }, error = function(e) {
    log_message(sprintf("Failed to load train data for %s", drug))
    return(NULL)
  })
  
  valid_data <- tryCatch({
    readRDS(paste0("/preprocessed_", drug, "_valid_selected.rds"))
  }, error = function(e) {
    log_message(sprintf("Failed to load valid data for %s", drug))
    return(NULL)
  })
  
  if (is.null(train_data) || is.null(valid_data)) {
    log_message(sprintf("Skipping %s: data missing", drug))
    next
  }
  
  log_message(sprintf("Train: %d samples, %d features | Valid: %d samples", 
                      nrow(train_data), ncol(train_data)-1, nrow(valid_data)))
  
  # Median imputation
  if (sum(is.na(train_data)) > 0) {
    log_message(sprintf("Imputing missing values in train data for %s", drug))
    for (col in colnames(train_data)) {
      if (any(is.na(train_data[[col]]))) {
        train_data[[col]][is.na(train_data[[col]])] <- median(train_data[[col]], na.rm = TRUE)
      }
    }
  }
  
  # Basic checks
  if (!"drug_response" %in% colnames(train_data) || !is.numeric(train_data$drug_response)) {
    log_message(sprintf("Invalid drug_response in %s", drug))
    next
  }
  
  # Remove near-zero variance and highly correlated features (optional cleanup)
  features <- train_data[, setdiff(colnames(train_data), "drug_response")]
  nzv <- nearZeroVar(features, saveMetrics = TRUE)
  if (sum(nzv$nzv) > 0) {
    features <- features[, !nzv$nzv]
  }
  high_corr <- findCorrelation(cor(features), cutoff = 0.8)
  if (length(high_corr) > 0) {
    features <- features[, -high_corr]
  }
  
  train_data <- cbind(features, drug_response = train_data$drug_response)
  valid_data <- valid_data[, c(colnames(features), "drug_response")]
  
  # Define models
  models <- list(
    ElasticNet = list(alpha = 0.5),
    Ridge      = list(alpha = 0),
    Lasso      = list(alpha = 1)
  )
  
  for (model_name in names(models)) {
    log_message(sprintf("Training %s for %s", model_name, drug))
    
    model_fit <- tryCatch({
      train(
        drug_response ~ .,
        data = train_data,
        method = "glmnet",
        trControl = trainControl(method = "cv", number = 5, savePredictions = "final"),
        tuneGrid = expand.grid(alpha = models[[model_name]]$alpha, 
                               lambda = 10^seq(-4, 1, length = 50)),
        metric = "RMSE",
        preProcess = c("center", "scale"),
        na.action = na.exclude
      )
    }, error = function(e) {
      log_message(sprintf("Training %s failed for %s: %s", model_name, drug, conditionMessage(e)))
      return(NULL)
    })
    
    if (is.null(model_fit)) next
    
    # Predict on validation set
    predictions <- predict(model_fit, newdata = valid_data)
    actual <- valid_data$drug_response
    
    pearson  <- cor(predictions, actual, method = "pearson",  use = "complete.obs")
    spearman <- cor(predictions, actual, method = "spearman", use = "complete.obs")
    
    # Store metrics
    model_metrics <- rbind(model_metrics, 
      data.frame(Drug = drug, Model = model_name, Pearson = pearson, Spearman = spearman)
    )
    
    log_message(sprintf("%s for %s → Pearson: %.3f | Spearman: %.3f", 
                        model_name, drug, pearson, spearman))
  }
}

# Combined performance plot
metrics_df <- model_metrics %>%
  pivot_longer(cols = c(Pearson, Spearman), names_to = "Metric", values_to = "Value")

bar_plot <- ggplot(metrics_df, aes(x = Model, y = Value, fill = Metric)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.4), width = 0.35) +
  facet_wrap(~ Drug, ncol = 3, scales = "free_y") +
  scale_fill_manual(values = pastel_palette, labels = c("Pearson", "Spearman")) +
  labs(title = "model_performance_combined",
       subtitle = "Pearson and Spearman Correlations on Validation Data",
       x = "Model", y = "Correlation Value") +
  plot_theme

print(bar_plot)

dev.off()

# Save outputs
saveRDS(model_metrics, "/model_metrics.rds")
write.csv(model_metrics, "model_metrics.csv", row.names = FALSE)

log_message("Model training completed successfully")

cat("\n=== Model Training Complete ===\n")




## Similarly done taking cross- datasets (CCLE data as training and GDSC for validation and generated results similarly)##
