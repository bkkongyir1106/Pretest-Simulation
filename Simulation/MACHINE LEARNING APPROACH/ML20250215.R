# --------------------------------------------------------------------------------
# Project: Machine Learning Approach for Simulation Data
# Description: Train ML models on 70% of the combined dataset and test on 30%,
# then validate the models on a completely different dataset.
# --------------------------------------------------------------------------------

# ---- Setup Environment ----

# Clear workspace
rm(list = ls())

# Load required libraries using pacman
if (!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest, gbm, lawstat, infotheo, ineq,
               caret, pROC, ROCR, randomForest, evd, discretization, nnet, ggplot2)

# Set working directory (adjust path as needed)
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/MACHINE LEARNING APPROACH")

# Source external user-defined functions (if any)
source("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/ML Models/fun.R")


# ---- Feature Engineering Functions ----

calculate_zero_crossing_rate <- function(samples) {
  signs <- samples > 0
  zero_crossings <- sum(abs(diff(signs)))
  return(zero_crossings / (length(samples) - 1))
}

calculate_gini_coefficient <- function(samples) {
  samples_abs <- abs(samples - min(samples))  # Shift to non-negative
  return(ineq(samples_abs, type = "Gini"))
}

calculate_outliers <- function(samples) {
  qnt <- quantile(samples, probs = c(0.25, 0.75))
  H <- 1.5 * IQR(samples)
  return(sum(samples < (qnt[1] - H) | samples > (qnt[2] + H)))
}

calculate_shapiro_wilk <- function(samples) {
  return(as.numeric(shapiro.test(samples)$statistic))
}

calculate_shapiro_francia <- function(samples) {
  return(as.numeric(lawstat::sf.test(samples)$statistic))
}

calculate_lilliefors <- function(samples) {
  return(as.numeric(nortest::lillie.test(samples)$statistic))
}

calculate_cramer_von_mises <- function(samples) {
  return(as.numeric(nortest::cvm.test(samples)$statistic))
}

calculate_entropy <- function(samples) {
  return(infotheo::entropy(discretize(samples, nbins = 10)))
}

calculate_mad <- function(samples) {
  return(mean(abs(samples - mean(samples))))
}

calculate_peak_to_trough <- function(samples) {
  return(max(samples) / abs(min(samples)))
}

calculate_sample_size <- function(samples) {
  return(length(samples))
}

calculate_range <- function(samples) {
  return(max(samples) - min(samples))
}

calculate_cv <- function(samples) {
  return(sd(samples) / mean(samples))
}

calculate_energy <- function(samples) {
  return(sum(samples^2))
}

# Combine all feature calculations into one function
calculate_features <- function(samples) {
  features <- data.frame(
    Skewness = e1071::skewness(samples),
    Kurtosis = e1071::kurtosis(samples),
    JB_Statistic = as.numeric(tseries::jarque.bera.test(samples)$statistic),
    AD_Statistic = as.numeric(nortest::ad.test(samples)$statistic),
    Zero_Crossing_Rate = calculate_zero_crossing_rate(samples),
    #Outliers = calculate_outliers(samples),
    Shapiro_Wilk = calculate_shapiro_wilk(samples),
    #Liliefors = calculate_lilliefors(samples),
    Cramer_Von_Mises = calculate_cramer_von_mises(samples),
    Range = calculate_range(samples),
    Coefficient_of_Variation = calculate_cv(samples)
    #Energy = calculate_energy(samples)
  )
  return(features)
}


# ---- Data Normalization Functions ----

min_max_normalize <- function(data) {
  return((data - min(data)) / (max(data) - min(data)))
}

normalize_data <- function(data) {
  numeric_features <- data[sapply(data, is.numeric)]  # Extract numeric features
  normalized_features <- as.data.frame(lapply(numeric_features, min_max_normalize))
  normalized_data <- cbind(normalized_features, Label = data$Label)  # Re-attach Label column
  return(normalized_data)
}

# ---- Data Generation Function ----
generate_data <- function(samples_size, N, dist = "normal", label) {
  data <- do.call(rbind, lapply(1:N, function(x) {
    samples <- generate_samples(samples_size, dist)
    features <- calculate_features(samples)
    features$Label <- label
    return(features)
  }))
  return(data)
}


# ---- Data Preparation and Model Training ----

# Set seed for reproducibility
set.seed(12345)

# Generate primary dataset
normal_data   <- generate_data(10, 1000, "normal",      "Normal")
lognormal     <- generate_data(10, 200,  "LogNormal",   "Non_Normal")
chisq_data    <- generate_data(10, 200,  "Chi-Square",  "Non_Normal")
exp_data      <- generate_data(10, 200,  "Exponential", "Non_Normal")
Weibull       <- generate_data(10, 200,  "Weibull",     "Non_Normal")
Laplace       <- generate_data(10, 200,  "Laplace",     "Non_Normal")
Gamma         <- generate_data(10, 200,  "Gamma",       "Non_Normal")
spike         <- generate_data(10, 200,  "spike",       "Non_Normal")
cauchy        <- generate_data(10, 200,  "cauchy",      "Non_Normal")
Uniform       <- generate_data(10, 200,  "Uniform",     "Non_Normal")

# Combine datasets: Normal vs. Non_Normal (using a subset of distributions)
non_normal_data <- rbind(lognormal, chisq_data, exp_data, Laplace, Uniform)
combined_data <- rbind(normal_data, non_normal_data)
combined_data$Label <- as.factor(combined_data$Label)

# Split combined data into training (70%) and testing (30%) sets
set.seed(123)
trainIndex <- createDataPartition(combined_data$Label, p = 0.7, list = FALSE)
train_data <- combined_data[trainIndex, ]
test_data  <- combined_data[-trainIndex, ]

# Shuffle the data rows
train_data <- train_data[sample(nrow(train_data)), ]
test_data  <- test_data[sample(nrow(test_data)), ]

# Normalize the datasets
Normalized_train_data <- normalize_data(train_data)
Normalized_test_data  <- normalize_data(test_data)


# ---- Model Training & Evaluation on Test Data ----

# Logistic Regression Model
log_model <- train(Label ~ ., data = Normalized_train_data, method = "glm",
                   family = "binomial")
log_pred <- predict(log_model, newdata = Normalized_test_data)
log_conf_matrix <- confusionMatrix(factor(log_pred, levels = levels(Normalized_test_data$Label)),
                                   Normalized_test_data$Label)
cat("Test Results for: Logistic Regression Model\n")
print(log_conf_matrix)

# Random Forest Model
rf_model <- train(Label ~ ., data = Normalized_train_data, method = "rf")
rf_pred <- predict(rf_model, newdata = Normalized_test_data)
rf_conf_matrix <- confusionMatrix(factor(rf_pred, levels = levels(Normalized_test_data$Label)),
                                  Normalized_test_data$Label)
cat("Test Results for: Random Forest Model\n")
print(rf_conf_matrix)

# Artificial Neural Network (ANN) Model
ann_model <- train(Label ~ ., data = Normalized_train_data, method = "nnet", trace = FALSE)
ann_pred <- predict(ann_model, newdata = Normalized_test_data)
ann_conf_matrix <- confusionMatrix(factor(ann_pred, levels = levels(Normalized_test_data$Label)),
                                   Normalized_test_data$Label)
cat("Test Results for: ANN Model\n")
print(ann_conf_matrix)

# Gradient Boosting Machine (GBM)
gbm_model <- train(Label ~ ., data = Normalized_train_data, method = "gbm", verbose = FALSE)
gbm_pred <- predict(gbm_model, newdata = Normalized_test_data)
gbm_conf_matrix <- confusionMatrix(factor(gbm_pred, levels = levels(Normalized_test_data$Label)),
                                   Normalized_test_data$Label)
cat("Test Results for: Gradient Boosting Trees\n")
print(gbm_conf_matrix)

# Support Vector Machines (SVM)
svm_model <- train(Label ~ ., data = Normalized_train_data, method = "svmRadial",
                   trControl = trainControl(classProbs = TRUE))
svm_pred <- predict(svm_model, newdata = Normalized_test_data)
svm_conf_matrix <- confusionMatrix(factor(svm_pred, levels = levels(Normalized_test_data$Label)),
                                   Normalized_test_data$Label)
cat("Test Results for: Support Vector Machines\n")
print(svm_conf_matrix)

# K-Nearest Neighbors (KNN)
knn_model <- train(Label ~ ., data = Normalized_train_data, method = "knn")
knn_pred <- predict(knn_model, newdata = Normalized_test_data)
knn_conf_matrix <- confusionMatrix(factor(knn_pred, levels = levels(Normalized_test_data$Label)),
                                   Normalized_test_data$Label)
cat("Test Results for: K-Nearest Neighbors\n")
print(knn_conf_matrix)


# ---- ROC Curves and AUC Calculation ----

plot_roc_curve <- function(models, data_test) {
  roc_curves <- list()
  auc_values <- list()
  for (model_name in names(models)) {
    model <- models[[model_name]]
    # Obtain predicted probabilities (assuming binary classification)
    pred_prob <- predict(model, newdata = data_test, type = "prob")[,2]
    roc_curve <- roc(data_test$Label, pred_prob, levels = rev(levels(data_test$Label)))
    roc_curves[[model_name]] <- roc_curve
    auc_values[[model_name]] <- auc(roc_curve)
  }
  if (length(roc_curves) > 0) {
    print(ggroc(roc_curves) +
            ggtitle("ROC Curves for Different Models") +
            theme_minimal())
  }
  return(auc_values)
}

# Store models in a list and plot ROC curves
models_list <- list("Logistic Regression" = log_model,
                    "Random Forest" = rf_model,
                    "ANN" = ann_model,
                    "GBM" = gbm_model,
                    "SVM" = svm_model,
                    "KNN" = knn_model)
auc_results <- plot_roc_curve(models_list, Normalized_test_data)
print(auc_results)

# ---- Variable Importance (Random Forest) ----

rf_var_imp <- varImp(rf_model)
print(rf_var_imp)
plot(rf_var_imp, main = "Variable Importance - Random Forest")


# ---- Model Validation on a Different Dataset ----

# Generate validation dataset using different distributions
set.seed(12345)
normal_valid     <- generate_data(10, 100, "normal", "Normal")
extremeskew_valid<- generate_data(10, 100, "extremeskew", "Non_Normal")
t_data          <- generate_data(10, 100, "t", "Non_Normal")

# Combine and shuffle validation data
validation_data <- rbind(normal_valid, t_data)
validation_data$Label <- as.factor(validation_data$Label)
validation_data <- validation_data[sample(nrow(validation_data)), ]

# Normalize validation data
Normalized_validation_data <- normalize_data(validation_data)
Normalized_validation_data$Label <- factor(Normalized_validation_data$Label,
                                           levels = levels(Normalized_validation_data$Label))

# Validate each model on the new dataset

# Logistic Regression Validation
log_pred_val <- predict(log_model, newdata = Normalized_validation_data)
log_conf_matrix_val <- confusionMatrix(factor(log_pred_val, levels = levels(Normalized_validation_data$Label)),
                                       Normalized_validation_data$Label)
cat("Validation Results for: Logistic Regression\n")
print(log_conf_matrix_val)

# Random Forest Validation
rf_pred_val <- predict(rf_model, newdata = Normalized_validation_data)
rf_conf_matrix_val <- confusionMatrix(factor(rf_pred_val, levels = levels(Normalized_validation_data$Label)),
                                      Normalized_validation_data$Label)
cat("Validation Results for: Random Forest\n")
print(rf_conf_matrix_val)

# ANN Validation
ann_pred_val <- predict(ann_model, newdata = Normalized_validation_data)
ann_conf_matrix_val <- confusionMatrix(factor(ann_pred_val, levels = levels(Normalized_validation_data$Label)),
                                       Normalized_validation_data$Label)
cat("Validation Results for: ANN\n")
print(ann_conf_matrix_val)

# GBM Validation
gbm_pred_val <- predict(gbm_model, newdata = Normalized_validation_data)
gbm_conf_matrix_val <- confusionMatrix(factor(gbm_pred_val, levels = levels(Normalized_validation_data$Label)),
                                       Normalized_validation_data$Label)
cat("Validation Results for: GBM\n")
print(gbm_conf_matrix_val)

# SVM Validation
svm_pred_val <- predict(svm_model, newdata = Normalized_validation_data)
svm_conf_matrix_val <- confusionMatrix(factor(svm_pred_val, levels = levels(Normalized_validation_data$Label)),
                                       Normalized_validation_data$Label)
cat("Validation Results for: SVM\n")
print(svm_conf_matrix_val)

# KNN Validation
knn_pred_val <- predict(knn_model, newdata = Normalized_validation_data)
knn_conf_matrix_val <- confusionMatrix(factor(knn_pred_val, levels = levels(Normalized_validation_data$Label)),
                                       Normalized_validation_data$Label)
cat("Validation Results for: KNN\n")
print(knn_conf_matrix_val)

