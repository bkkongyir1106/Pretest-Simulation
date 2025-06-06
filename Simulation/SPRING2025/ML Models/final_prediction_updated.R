# Required Libraries
rm(list = ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest, gbm, lawstat, infotheo, ineq, caret, pROC, ROCR, randomForest, evd, discretization, nnet, ggplot2)

# Set directories 
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/MACHINE LEARNING APPROACH")

# Read-in-external user-defined functions
source("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/ML Models/fun.R")

# Zero-crossing rate
calculate_zero_crossing_rate <- function(samples) {
  signs <- samples > 0
  zero_crossings <- sum(abs(diff(signs)))
  zero_crossing_rate <- zero_crossings / (length(samples) - 1)
  return(zero_crossing_rate)
}

# GINI Coefficient 
calculate_gini_coefficient <- function(samples) {
  samples_abs <- abs(samples - min(samples))  # Shift to non-negative
  gini_value <- ineq(samples_abs, type = "Gini")
  return(gini_value)
}

# Outliers
calculate_outliers <- function(samples) {
  qnt <- quantile(samples, probs = c(0.25, 0.75))
  H <- 1.5 * IQR(samples)
  outliers <- sum(samples < (qnt[1] - H) | samples > (qnt[2] + H))
  return(outliers)
}

# SW test
calculate_shapiro_wilk <- function(samples) {
  shapiro_stat <- shapiro.test(samples)$statistic
  return(shapiro_stat)
}

# Shapiro-Francia Test
calculate_shapiro_francia <- function(samples) {
  SF.test <- lawstat::sf.test(samples)$statistic
  return(as.numeric(SF.test))
}

# Liliefors Test 
calculate_lilliefors <- function(samples) {
  LF.test <- nortest::lillie.test(samples)$statistic
  return(as.numeric(LF.test))
}

# Cramer-von Mises Test
calculate_cramer_von_mises <- function(samples) {
  CVM.test <- nortest::cvm.test(samples)$statistic
  return(as.numeric(CVM.test))
}

# Entropy
calculate_entropy <- function(samples) {
  entropy_value <- infotheo::entropy(discretize(samples, nbins = 10))
  return(entropy_value)
}

# Mean Absolute Deviation
calculate_mad <- function(samples) {
  mad_value <- mean(abs(samples - mean(samples)))
  return(mad_value)
}

# Peak-to-Trough Ratio
calculate_peak_to_trough <- function(samples) {
  peak_to_trough_ratio <- max(samples) / abs(min(samples))
  return(peak_to_trough_ratio)
}

# Sample size
calculate_sample_size <- function(samples) {
  return(length(samples))
}

# Range
calculate_range <- function(samples) {
  return(max(samples) - min(samples))
}

# CV
calculate_cv <- function(samples) {
  return(sd(samples) / mean(samples))
}

# Energy Sum of squares
calculate_energy <- function(samples) {
  return(sum(samples^2))
}

calculate_features <- function(samples) {
  skewness <- e1071::skewness(samples)
  kurtosis <- e1071::kurtosis(samples)
  jb_stat <- tseries::jarque.bera.test(samples)$statistic
  ad_stat <- nortest::ad.test(samples)$statistic
  zero_crossing_rate <- calculate_zero_crossing_rate(samples)
  gini_coefficient <- calculate_gini_coefficient(samples)
  outliers <- calculate_outliers(samples)
  shapiro_wilk <- calculate_shapiro_wilk(samples)
  lilliefors_stat <- calculate_lilliefors(samples)
  cramer_von_mises <- calculate_cramer_von_mises(samples)
  sample_size <- calculate_sample_size(samples)
  range <- calculate_range(samples)
  cv <- calculate_cv(samples)
  energy <- calculate_energy(samples)
  
  features <- data.frame(
    Skewness = skewness,
    Kurtosis = kurtosis,
    JB_Statistic = as.numeric(jb_stat),
    AD_Statistic = as.numeric(ad_stat),
    Zero_Crossing_Rate = zero_crossing_rate,
    Outliers = outliers,
    Shapiro_Wilk = as.numeric(shapiro_wilk),
    Liliefors = lilliefors_stat,
    Cramer_Von_Mises = cramer_von_mises,
    #Sample_Size = sample_size,
    Range = range,
    Coefficient_of_Variation = cv,
    Energy = energy
  )
  
  return(features)
}

# Min-Max normalization function
min_max_normalize <- function(data) {
  return((data - min(data)) / (max(data) - min(data)))
}

# Normalize all numeric columns except the Label column
normalize_data <- function(data) {
  numeric_features <- data[sapply(data, is.numeric)]  # Extract numeric features
  normalized_features <- as.data.frame(lapply(numeric_features, min_max_normalize))  # Normalize features
  normalized_data <- cbind(normalized_features, Label = data$Label)  # Re-attach the Label column
  return(normalized_data)
}

# Step 2: Prepare Data for Machine Learning Models
generate_data <- function(samples_size, N, dist = "normal", label) {
  data <- do.call(rbind, lapply(1:N, function(x) {
    samples <- generate_samples(samples_size, dist)
    features <- calculate_features(samples)
    features$Label <- label
    return(features)
  }))
  return(data)
}

  # Generate data
  set.seed(12345)
  normal_data <- generate_data(8, 2000, "normal", "Normal")
  lognormal <- generate_data(8, 200, "LogNormal", "Non_Normal")
  chisq_data <- generate_data(8, 200, "Chi-Square", "Non_Normal")
  exp_data <- generate_data(8, 200, "Exponential", "Non_Normal")
  Weibull <- generate_data(8, 200, "Weibull", "Non_Normal")
  Laplace <- generate_data(8, 200, "Laplace", "Non_Normal")
  Gamma <- generate_data(8, 200, "Gamma", "Non_Normal")
  #spike <- generate_data(8, 200, "spike", "Non_Normal")
  cauchy <- generate_data(8, 200, "cauchy", "Non_Normal")
  Uniform <- generate_data(8, 200, "Uniform", "Non_Normal")
  t_data <- generate_data(8, 200, "t", "Non_Normal")
  
  # Combine data
  non_normal_data <- rbind(lognormal, chisq_data, exp_data, Weibull, Laplace, Gamma, cauchy, Uniform)
  data <- rbind(normal_data, non_normal_data)
  data$Label <- as.factor(data$Label)
  # Split Data into Training and Test Sets
  set.seed(123)
  trainIndex <- createDataPartition(data$Label, p = .7, list = FALSE, times = 1)
  train_data <- data[trainIndex, ]
  test_data  <- data[-trainIndex, ]
  train_data <- train_data[sample(nrow(train_data)), ]
  test_data <- test_data[sample(nrow(test_data)), ]
  
  # Normalize the training and test data
  Normalized_train_data <- normalize_data(train_data)
  Normalized_test_data <- normalize_data(test_data)
  
  # Step 4: Build Machine Learning Models
  # Logistic Regression Model
  log_model <- train(Label ~ ., data = Normalized_train_data, method = "glm", family = "binomial")
  log_pred <- predict(log_model, newdata = Normalized_test_data)
  log_pred <- factor(log_pred, levels = levels(Normalized_test_data$Label))
  log_conf_matrix <- confusionMatrix(log_pred, Normalized_test_data$Label)
  cat("Test Results for:", "Logistic Regression Model" ,"\n")
  print(log_conf_matrix)
  
  # Random Forest Model
  rf_model <- train(Label ~ ., data = Normalized_train_data, method = "rf")
  rf_pred <- predict(rf_model, newdata = Normalized_test_data)
  rf_pred <- factor(rf_pred, levels = levels(Normalized_test_data$Label))
  rf_conf_matrix <- confusionMatrix(rf_pred, Normalized_test_data$Label)
  cat("Test Results for:", "Random Forest Model" ,"\n")
  print(rf_conf_matrix)
  
  # ANN Model
  ann_model <- train(Label ~ ., data = Normalized_train_data, method = "nnet", trace = FALSE)
  ann_pred <- predict(ann_model, newdata = Normalized_test_data)
  ann_pred <- factor(ann_pred, levels = levels(Normalized_test_data$Label))
  ann_conf_matrix <- confusionMatrix(ann_pred, Normalized_test_data$Label)
  print(ann_conf_matrix)
  
  # Gradient Boosting Trees (GBM)
  gbm_model <- train(Label ~ ., data = Normalized_train_data, method = "gbm", verbose = FALSE)
  gbm_pred <- predict(gbm_model, newdata = Normalized_test_data)
  gbm_pred <- factor(gbm_pred, levels = levels(Normalized_test_data$Label))
  gbm_conf_matrix <- confusionMatrix(gbm_pred, Normalized_test_data$Label)
  cat("Test Results for:", "Gradient Boosting Trees" ,"\n")
  print(gbm_conf_matrix)
  
  # Support Vector Machines (SVM)
  svm_model <- train(Label ~ ., data = Normalized_train_data, method = "svmRadial", trControl = trainControl(classProbs = TRUE))
  svm_pred <- predict(svm_model, newdata = Normalized_test_data)
  svm_pred <- factor(svm_pred, levels = levels(Normalized_test_data$Label))
  svm_conf_matrix <- confusionMatrix(svm_pred, Normalized_test_data$Label)
  cat("Test Results for:", "Support Vector Machines" ,"\n")
  print(svm_conf_matrix)
  
  # K-Nearest Neighbors (KNN)
  knn_model <- train(Label ~ ., data = Normalized_train_data, method = "knn")
  knn_pred <- predict(knn_model, newdata = Normalized_test_data)
  knn_pred <- factor(knn_pred, levels = levels(Normalized_test_data$Label))
  knn_conf_matrix <- confusionMatrix(knn_pred, Normalized_test_data$Label)
  cat("Test Results for:", "K-Nearest Neighbors" ,"\n")
  print(knn_conf_matrix)
  
  # Create a function to plot ROC curves and compute AUC values
  plot_roc_curve <- function(models, Normalized_test_data) {
    roc_curves <- list()
    auc_values <- list()
    for (model_name in names(models)) {
      model <- models[[model_name]]
      pred_prob <- predict(model, newdata = Normalized_test_data, type = "prob")[,2]
      roc_curve <- roc(Normalized_test_data$Label, pred_prob, levels = rev(levels(Normalized_test_data$Label)))
      roc_curves[[model_name]] <- roc_curve
      auc_values[[model_name]] <- auc(roc_curve)
    }
    # Plot ROC curves if available
    if (length(roc_curves) > 0) {
      print(ggroc(roc_curves) +
              ggtitle("ROC Curves for Different Models") +
              theme_minimal())
    }
    
    # Return AUC values
    return(auc_values)
  }
  
  # Store models in a list
  models <- list(
    "Logistic Regression" = log_model,
    "Random Forest" = rf_model,
    "ANN" = ann_model,
    "GBM" = gbm_model,
    "SVM" = svm_model,
    "KNN" = knn_model
  )
  
  # Plot the ROC curves 
  auc_results <- plot_roc_curve(models, Normalized_test_data)
  # Print AUC values
  print(auc_results)
  
  # Variable Importance Plot for Random Forest
  rf_var_imp <- varImp(rf_model)
  print(rf_var_imp)
  
  # Plot the Variable Importance
  plot(rf_var_imp, main = "Variable Importance - Random Forest")
  
  # Generate data from distribution different from those 
  set.seed(12345)
  normal <- generate_data(8, 200, "normal", "Normal")
  extremeskew <- generate_data(8, 200, "extremeskew", "Non_Normal")
  t_data <- generate_data(8, 100, "t", "Non_Normal")

  # Combine the validation data
  validation_data <- rbind(normal, extremeskew)
  
  # Ensure that Label is a factor
  validation_data$Label <- as.factor(validation_data$Label)
  
  # Shuffle the validation data
  validation_data <- validation_data[sample(nrow(validation_data)), ]
  
  # Normalize the training and test data
  Normalized_train_data <- normalize_data(validation_data)
  
  # Ensure Label is a factor in the validation data
  Normalized_train_data$Label <- factor(Normalized_train_data$Label, levels = levels(Normalized_train_data$Label))
  
  # logistic Regression
  log_pred_val <- predict(log_model, newdata = Normalized_train_data)
  log_pred_val <- factor(log_pred_val, levels = levels(Normalized_train_data$Label))
  log_conf_matrix_val <- confusionMatrix(log_pred_val, Normalized_train_data$Label)
  cat("Validation Results for:", "logistic Regression" ,"\n")
  print(log_conf_matrix_val)
  
  # Random Forest
  rf_pred_val <- predict(rf_model, newdata = Normalized_train_data)
  rf_pred_val <- factor(rf_pred_val, levels = levels(Normalized_train_data$Label))
  rf_conf_matrix_val <- confusionMatrix(rf_pred_val, Normalized_train_data$Label)
  cat("Validation Results for:", "Random Forest" ,"\n")
  print(rf_conf_matrix_val)
  
  # ANN 
  ann_pred_val <- predict(ann_model, newdata = Normalized_train_data)
  ann_pred_val <- factor(ann_pred_val, levels = levels(Normalized_train_data$Label))
  ann_conf_matrix_val <- confusionMatrix(ann_pred_val, Normalized_train_data$Label)
  cat("Validation Results for:", "ANN" ,"\n")
  print(ann_conf_matrix_val)
  
  # GB
  gbm_pred_val <- predict(gbm_model, newdata = Normalized_train_data)
  gbm_pred_val <- factor(gbm_pred_val, levels = levels(Normalized_train_data$Label))
  gbm_conf_matrix_val <- confusionMatrix(gbm_pred_val, Normalized_train_data$Label)
  cat("Validation Results for:", "GB" ,"\n")
  print(gbm_conf_matrix_val)
  
  # SVM
  svm_pred_val <- predict(svm_model, newdata = Normalized_train_data)
  svm_pred_val <- factor(svm_pred_val, levels = levels(Normalized_train_data$Label))
  svm_conf_matrix_val <- confusionMatrix(svm_pred_val, Normalized_train_data$Label)
  cat("Validation Results for:", "SVM" ,"\n")
  print(svm_conf_matrix_val)
  
  # KNN
  knn_pred_val <- predict(knn_model, newdata = Normalized_train_data)
  knn_pred_val <- factor(knn_pred_val, levels = levels(Normalized_train_data$Label))
  knn_conf_matrix_val <- confusionMatrix(knn_pred_val, Normalized_train_data$Label)
  cat("Validation Results for:", "SKNN" ,"\n")
  print(knn_conf_matrix_val)

