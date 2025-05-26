# ---------------------------
# Clean and Improved R Code for Training & Single Dataset Prediction
# ---------------------------
rm(list = ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest, gbm, lawstat, infotheo, ineq, caret, pROC, ROCR, 
               randomForest, evd, discretization, nnet, ggplot2)

# ---------------------------
# Set Directories & Load Functions
# ---------------------------
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/MACHINE LEARNING APPROACH")
source("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/MACHINE LEARNING APPROACH/fun.R")  

# ---------------------------
# Define Feature Extraction Functions
# ---------------------------
calculate_features <- function(samples) {
  features <- data.frame(
    Skewness = e1071::skewness(samples),
    Kurtosis = e1071::kurtosis(samples),
    JB_Statistic = as.numeric(tseries::jarque.bera.test(samples)$statistic),
    AD_Statistic = as.numeric(nortest::ad.test(samples)$statistic),
    Zero_Crossing_Rate = sum(abs(diff(samples > 0))) / (length(samples) - 1),
    Outliers = sum(samples < quantile(samples, 0.25) - 1.5 * IQR(samples) |
                     samples > quantile(samples, 0.75) + 1.5 * IQR(samples)),
    Shapiro_Wilk = as.numeric(shapiro.test(samples)$statistic),
    Liliefors = as.numeric(nortest::lillie.test(samples)$statistic),
    Cramer_Von_Mises = as.numeric(nortest::cvm.test(samples)$statistic),
    Range = max(samples) - min(samples),
    Coefficient_of_Variation = sd(samples) / mean(samples),
    Energy = sum(samples^2)
  )
  return(features)
}

# ---------------------------
# Standardization and Normalization
# ---------------------------
preprocess_data <- function(train_data, other_data) {
  numeric_train <- train_data[, sapply(train_data, is.numeric)]
  preProcValues <- preProcess(numeric_train, method = c("center", "scale"))
  train_processed <- predict(preProcValues, numeric_train)
  numeric_other <- other_data[, sapply(other_data, is.numeric)]
  other_processed <- predict(preProcValues, numeric_other)
  train_processed$Label <- train_data$Label
  other_processed$Label <- other_data$Label
  return(list(train = train_processed, other = other_processed, preProcValues = preProcValues))
}

# ---------------------------
# Data Generation
# ---------------------------
generate_data <- function(sample_size, N, dist = "normal", label) {
  data <- do.call(rbind, lapply(1:N, function(x) {
    samples <- generate_samples(sample_size, dist)  
    features <- calculate_features(samples)
    features$Label <- label
    return(features)
  }))
  return(data)
}

# ---------------------------
# Generate and Prepare Data
# ---------------------------
set.seed(12345)
sample_size <- 8  
N <- 200            
normal_data <- generate_data(sample_size, N, "normal", "Normal")
non_normal_data <- rbind(
  generate_data(sample_size, N, "LogNormal", "Non_Normal"),
  generate_data(sample_size, N, "Chi-Square", "Non_Normal"),
  generate_data(sample_size, N, "Exponential", "Non_Normal"),
  generate_data(sample_size, N, "Weibull", "Non_Normal"),
  generate_data(sample_size, N, "Laplace", "Non_Normal"),
  generate_data(sample_size, N, "Gamma", "Non_Normal"),
  generate_data(sample_size, N, "beta", "Non_Normal"),
  generate_data(sample_size, N, "Uniform", "Non_Normal"),
  generate_data(sample_size, N, "t", "Non_Normal")
)

data_all <- rbind(normal_data, non_normal_data)
data_all$Label <- as.factor(data_all$Label)

# Split Data
set.seed(123)
trainIndex <- createDataPartition(data_all$Label, p = 0.7, list = FALSE)
train_data <- data_all[trainIndex, ]
test_data  <- data_all[-trainIndex, ]
train_data <- train_data[sample(nrow(train_data)), ]
test_data  <- test_data[sample(nrow(test_data)), ]

# Apply Preprocessing
norm_result <- preprocess_data(train_data, test_data)
train_norm <- norm_result$train
test_norm  <- norm_result$other
preProcValues <- norm_result$preProcValues  

# ---------------------------
# Train Machine Learning Models
# ---------------------------
ctrl <- trainControl(method = "cv", number = 10, summaryFunction = twoClassSummary, classProbs = TRUE, savePredictions = "final")
models_list <- list(
  "Logistic Regression" = train(Label ~ ., data = train_norm, method = "glm", family = "binomial", trControl = ctrl, metric = "ROC"),
  "Random Forest" = train(Label ~ ., data = train_norm, method = "rf", trControl = ctrl, metric = "ROC"),
  "ANN" = train(Label ~ ., data = train_norm, method = "nnet", trControl = ctrl, metric = "ROC", trace = FALSE),
  "GBM" = train(Label ~ ., data = train_norm, method = "gbm", trControl = ctrl, metric = "ROC", verbose = FALSE),
  "SVM" = train(Label ~ ., data = train_norm, method = "svmRadial", trControl = ctrl, metric = "ROC"),
  "KNN" = train(Label ~ ., data = train_norm, method = "knn", trControl = ctrl, metric = "ROC")
)

# ---------------------------
# Prediction Function for Single Dataset
# ---------------------------
predict_single_dataset <- function(file_path = NULL, trained_models, preProcValues) {
  if(is.null(file_path)) file_path <- file.choose()
  user_data <- read.csv(file_path, header = TRUE)
  samples <- user_data[[1]]
  user_features <- calculate_features(samples)
  user_features_norm <- predict(preProcValues, user_features)
  predictions <- sapply(trained_models, function(model) predict(model, newdata = user_features_norm))
  return(predictions)
}

# ---------------------------
# Example Usage for Single Dataset Prediction
# ---------------------------
single_predictions <- predict_single_dataset(file_path = NULL, trained_models = models_list, preProcValues = preProcValues)
print(single_predictions)
