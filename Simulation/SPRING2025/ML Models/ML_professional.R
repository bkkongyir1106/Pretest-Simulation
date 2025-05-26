# Load required libraries
rm(list = ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest, gbm, lawstat, infotheo, ineq, caret, pROC, ROCR, randomForest, evd, discretization, nnet, ggplot2)

# Set working directory
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/MACHINE LEARNING APPROACH")

# Load external user-defined functions
source("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/MACHINE LEARNING APPROACH/fun.R")

# Feature Calculation Functions
calculate_zero_crossing_rate <- function(samples) sum(diff(samples > 0)) / length(samples)
calculate_gini_coefficient <- function(samples) ineq(samples, type = "Gini")
calculate_outliers <- function(samples) {
  qnt <- quantile(samples, probs = c(0.25, 0.75))
  H <- 1.5 * IQR(samples)
  sum(samples < (qnt[1] - H) | samples > (qnt[2] + H))
}
calculate_shapiro_wilk <- function(samples) shapiro.test(samples)$statistic
calculate_shapiro_francia <- function(samples) as.numeric(lawstat::sf.test(samples)$statistic)
calculate_lilliefors <- function(samples) as.numeric(nortest::lillie.test(samples)$statistic)
calculate_cramer_von_mises <- function(samples) as.numeric(nortest::cvm.test(samples)$statistic)
calculate_entropy <- function(samples) infotheo::entropy(discretize(samples, nbins = 10))
calculate_mad <- function(samples) mean(abs(samples - mean(samples)))
calculate_peak_to_trough <- function(samples) max(samples) / abs(min(samples))
calculate_sample_size <- function(samples) length(samples)
calculate_range <- function(samples) max(samples) - min(samples)
calculate_cv <- function(samples) sd(samples) / mean(samples)
calculate_energy <- function(samples) sum(samples^2)

# Feature Extraction Function
calculate_features <- function(samples) {
  data.frame(
    Skewness = e1071::skewness(samples),
    Kurtosis = e1071::kurtosis(samples),
    JB_Statistic = as.numeric(tseries::jarque.bera.test(samples)$statistic),
    AD_Statistic = as.numeric(nortest::ad.test(samples)$statistic),
    #Zero_Crossing_Rate = calculate_zero_crossing_rate(samples),
    Outliers = calculate_outliers(samples),
    Shapiro_Wilk = calculate_shapiro_wilk(samples),
    Liliefors = calculate_lilliefors(samples),
    Cramer_Von_Mises = calculate_cramer_von_mises(samples),
    Sample_Size = calculate_sample_size(samples),
    Range = calculate_range(samples)
    #Coefficient_of_Variation = calculate_cv(samples)
    #Energy = calculate_energy(samples)
  )
}

# Normalization Function
min_max_normalize <- function(data) (data - min(data)) / (max(data) - min(data))
normalize_data <- function(data) {
  numeric_features <- data[sapply(data, is.numeric)]
  normalized_features <- as.data.frame(lapply(numeric_features, min_max_normalize))
  cbind(normalized_features, Label = data$Label)
}

# Data Generation Function
generate_data <- function(samples_size, N, dist = "normal", label) {
  do.call(rbind, lapply(1:N, function(x) {
    samples <- generate_samples(samples_size, dist)
    features <- calculate_features(samples)
    features$Label <- label
    return(features)
  }))
}

# Model Training and Evaluation Function
train_and_evaluate <- function(model_name, train_data, test_data, method, ...) {
  model <- train(Label ~ ., data = train_data, method = method, ...)
  pred <- predict(model, newdata = test_data)
  pred <- factor(pred, levels = levels(test_data$Label))
  conf_matrix <- confusionMatrix(pred, test_data$Label)
  list(model = model, conf_matrix = conf_matrix)
}

# ROC Curve and AUC Calculation Function
plot_roc_curve <- function(models, test_data) {
  roc_curves <- lapply(models, function(model) {
    pred_prob <- predict(model, newdata = test_data, type = "prob")[, 2]
    roc(test_data$Label, pred_prob, levels = rev(levels(test_data$Label)))
  })
  auc_values <- sapply(roc_curves, auc)
  print(ggroc(roc_curves) + ggtitle("ROC Curves for Different Models") + theme_minimal())
  return(auc_values)
}

# Main Execution
runtime <- system.time({
  set.seed(123)
  normal_data <- generate_data(8, 2100, "normal", "Normal")
  non_normal_data <- do.call(rbind, lapply(c("LogNormal", "Chi-Square", "Exponential", "Weibull", "Laplace", "Gamma", "Uniform"), 
                  function(dist) generate_data(8, 300, dist, "Non_Normal")))
  
  # Combine data
  data <- rbind(normal_data, non_normal_data)
  data$Label <- as.factor(data$Label)
  
  # Split Data into Training and Test Sets
  set.seed(123)
  trainIndex <- createDataPartition(data$Label, p = .7, list = FALSE, times = 1)
  train_data <- data[trainIndex, ]
  test_data  <- data[-trainIndex, ]
  train_data <- train_data[sample(nrow(train_data)), ]
  test_data <- test_data[sample(nrow(test_data)), ]
  
  Normalized_train_data <- normalize_data(train_data)
  Normalized_test_data <- normalize_data(test_data)
  
  models <- list(
    "Logistic Regression" = train_and_evaluate("Logistic Regression", train_data, test_data, "glm", family = "binomial"),
    "Random Forest" = train_and_evaluate("Random Forest", train_data, test_data, "rf"),
    "ANN" = train_and_evaluate("ANN", train_data, test_data, "nnet", trace = FALSE),
    "GBM" = train_and_evaluate("GBM", train_data, test_data, "gbm", verbose = FALSE),
    "SVM" = train_and_evaluate("SVM", train_data, test_data, "svmRadial", trControl = trainControl(classProbs = TRUE)),
    "KNN" = train_and_evaluate("KNN", train_data, test_data, "knn")
  )
  
  auc_results <- plot_roc_curve(lapply(models, `[[`, "model"), test_data)
  print(auc_results)
  
  # Variable Importance for Random Forest
  rf_var_imp <- varImp(models[["Random Forest"]]$model)
  print(rf_var_imp)
  
  # Plot Variable Importance using ggplot2
  var_imp_data <- rf_var_imp$importance
  var_imp_data$Feature <- rownames(var_imp_data)
  var_imp_data <- var_imp_data[order(-var_imp_data$Overall), ]
  
  ggplot(var_imp_data, aes(x = reorder(Feature, Overall), y = Overall)) +
    geom_bar(stat = "identity", fill = "steelblue") +
    coord_flip() +
    labs(title = "Variable Importance - Random Forest",
         x = "Feature",
         y = "Importance") +
    theme_minimal()
  
  # Validation Data
  validation_data <- do.call(rbind, lapply(c("normal", "F", "Gumbel", "logistic"), 
                                           function(dist) generate_data(8, 300, dist, ifelse(dist == "normal", "Normal", "Non_Normal"))))
  validation_data$Label <- as.factor(validation_data$Label)
  validation_data <- validation_data[sample(nrow(validation_data)), ]
  
  # Test Evaluation
  for (model_name in names(models)) {
    pred_val <- predict(models[[model_name]]$model, newdata = test_data)
    pred_val <- factor(pred_val, levels = levels(test_data$Label))
    conf_matrix_val <- confusionMatrix(pred_val, test_data$Label)
    print(paste("Test Confusion Matrix for", model_name))
    print(conf_matrix_val)
  }
  
  # Validation Evaluation
  for (model_name in names(models)) {
    pred_val <- predict(models[[model_name]]$model, newdata = validation_data)
    pred_val <- factor(pred_val, levels = levels(validation_data$Label))
    conf_matrix_val <- confusionMatrix(pred_val, validation_data$Label)
    print(paste("Validation Confusion Matrix for", model_name))
    print(conf_matrix_val)
  }
})

print(runtime)

