rm(list = ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest, gbm, lawstat, infotheo, ineq, caret, pROC, ROCR, 
               randomForest, evd, discretization, nnet, ggplot2)

# ---------------------------
# Set Directories & Load Functions
# ---------------------------
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/MACHINE LEARNING APPROACH")
source("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/ML Models/fun.R") 

# ---------------------------
# Define Feature Extraction Functions
# ---------------------------
calculate_zero_crossing_rate <- function(samples) {
  signs <- samples > 0
  zero_crossings <- sum(abs(diff(signs)))
  return(zero_crossings / (length(samples) - 1))
}

calculate_gini_coefficient <- function(samples) {
  samples_abs <- abs(samples - min(samples))  
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

calculate_features <- function(samples) {
  skewness <- e1071::skewness(samples)
  kurtosis <- e1071::kurtosis(samples)
  jb_stat <- as.numeric(tseries::jarque.bera.test(samples)$statistic)
  ad_stat <- as.numeric(nortest::ad.test(samples)$statistic)
  zero_crossing_rate <- calculate_zero_crossing_rate(samples)
  gini_coefficient <- calculate_gini_coefficient(samples)
  outliers <- calculate_outliers(samples)
  shapiro_wilk <- calculate_shapiro_wilk(samples)
  lilliefors_stat <- calculate_lilliefors(samples)
  cramer_von_mises <- calculate_cramer_von_mises(samples)
  range_val <- calculate_range(samples)
  cv <- calculate_cv(samples)
  energy <- calculate_energy(samples)
  
  features <- data.frame(
    #JB_Statistics incorporate both skewness and kurtosis, so we can exclude them
    Skewness = skewness,
    Kurtosis = kurtosis,
    JB_Statistic = jb_stat,
    AD_Statistic = ad_stat,
    Zero_Crossing_Rate = zero_crossing_rate,
    Outliers = outliers,
    Shapiro_Wilk = shapiro_wilk,
    Liliefors = lilliefors_stat,
    Cramer_Von_Mises = cramer_von_mises,
    Range = range_val,
    Coefficient_of_Variation = cv,
    Energy = energy
  )
  return(features)
}

# ---------------------------
# Data Normalization using caret's preProcess
# ---------------------------
normalize_dataset <- function(train_data, other_data) {
  numeric_train <- train_data[, sapply(train_data, is.numeric)]
  preProcValues <- preProcess(numeric_train, method = c("center", "scale"))
  train_norm <- predict(preProcValues, numeric_train) 
  numeric_other <- other_data[, sapply(other_data, is.numeric)]
  other_norm <- predict(preProcValues, numeric_other)
  
  train_norm$Label <- train_data$Label
  other_norm$Label <- other_data$Label
  return(list(train = train_norm, other = other_norm, preProcValues = preProcValues))
}

# ---------------------------
# Data Generation Function (balanced sample size across groups)
# ---------------------------
generate_data <- function(sample_size, N, dist = "normal", label) {
  data <- do.call(rbind, lapply(1:N, function(x) {
    samples <- generate_samples(sample_size, dist)  # generate_samples() must be defined in fun.R
    features <- calculate_features(samples)
    features$Label <- label
    return(features)
  }))
  return(data)
}

# ---------------------------
# Set Seed for Reproducibility and Generate Training Data
# ---------------------------
set.seed(12345)
sample_size <- 8  # consistent sample size per dataset
N <- 500            # number of samples per distribution

normal_data1 <- generate_data(sample_size, 5*N, "normal", "Normal")
normal_data2 <- generate_data(sample_size, 5*N, "normal", "Normal")
# Generate non-normal data from several distributions
lognormal <- generate_data(sample_size, N, "LogNormal", "Non_Normal")
chisq_data   <- generate_data(sample_size, N, "Chi-Square", "Non_Normal")
exp_data     <- generate_data(sample_size, N, "Exponential", "Non_Normal")
Weibull      <- generate_data(sample_size, N, "Weibull", "Non_Normal")
Pareto      <- generate_data(sample_size, N, "Pareto", "Non_Normal")
Laplace      <- generate_data(sample_size, N, "Laplace", "Non_Normal")
Gamma        <- generate_data(sample_size, N, "Gamma", "Non_Normal")
#Gumbel       <- generate_data(sample_size, N, "Gumbel", "Non_Normal")
Uniform      <- generate_data(sample_size, N, "Uniform", "Non_Normal")
t_data       <- generate_data(sample_size, N, "t", "Non_Normal")
t10_data       <- generate_data(sample_size, N, "t", "Non_Normal")

non_normal_data <- rbind(lognormal, chisq_data, exp_data, Weibull, Pareto, Laplace, Gamma, t10_data, Uniform, t_data)
data_all <- rbind(normal_data1, normal_data2, non_normal_data)
data_all$Label <- as.factor(data_all$Label)

# ---------------------------
# Split Data into Training and Test Sets (Stratified)
# ---------------------------
set.seed(123)
trainIndex <- createDataPartition(data_all$Label, p = 0.7, list = FALSE)
train_data <- data_all[trainIndex, ]
test_data  <- data_all[-trainIndex, ]
train_data <- train_data[sample(nrow(train_data)), ]
test_data  <- test_data[sample(nrow(test_data)), ]

# ---------------------------
# Normalize Data (Learn normalization from training set)
# ---------------------------
norm_result <- normalize_dataset(train_data, test_data)
train_norm <- norm_result$train
test_norm  <- norm_result$other
preProcValues <- norm_result$preProcValues  # Save for future predictions

# ---------------------------
# Define Common Training Control for Cross-Validation
# ---------------------------
ctrl <- trainControl(method = "cv", number = 10, 
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE, 
                     savePredictions = "final")

# ---------------------------
# Train Machine Learning Models
# ---------------------------
log_model <- train(Label ~ .,
                   data = train_norm,
                   method = "glm", 
                   family = "binomial", 
                   trControl = ctrl, 
                   metric = "ROC")
log_pred <- predict(log_model, newdata = test_norm)
cat("Test Results for: Logistic Regression Model\n")
print(confusionMatrix(log_pred, test_norm$Label))

rf_model <- train(Label ~ ., 
                  data = train_norm, 
                  method = "rf", 
                  trControl = ctrl, 
                  metric = "ROC")
rf_pred <- predict(rf_model, newdata = test_norm)
cat("Test Results for: Random Forest Model\n")
print(confusionMatrix(rf_pred, test_norm$Label))

ann_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "nnet", 
                   trControl = ctrl, 
                   metric = "ROC", 
                   trace = FALSE)
ann_pred <- predict(ann_model, newdata = test_norm)
cat("Test Results for: ANN Model\n")
print(confusionMatrix(ann_pred, test_norm$Label))

gbm_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "gbm", 
                   trControl = ctrl, 
                   metric = "ROC", 
                   verbose = FALSE)
gbm_pred <- predict(gbm_model, newdata = test_norm)
cat("Test Results for: GBM Model\n")
print(confusionMatrix(gbm_pred, test_norm$Label))

svm_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "svmRadial", 
                   trControl = ctrl, 
                   metric = "ROC")
svm_pred <- predict(svm_model, newdata = test_norm)
cat("Test Results for: SVM Model\n")
print(confusionMatrix(svm_pred, test_norm$Label))

knn_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "knn", 
                   trControl = ctrl, 
                   metric = "ROC")
knn_pred <- predict(knn_model, newdata = test_norm)
cat("Test Results for: KNN Model\n")
print(confusionMatrix(knn_pred, test_norm$Label))

# Collect all models into a list for later use
models_list <- list(
  "Logistic Regression" = log_model,
  "Random Forest" = rf_model,
  "ANN" = ann_model,
  "GBM" = gbm_model,
  "SVM" = svm_model,
  "KNN" = knn_model
)

# ---------------------------
# Plot ROC Curves and Compute AUC for All Models
# ---------------------------
plot_roc_curve <- function(models, test_data) {
  roc_curves <- list()
  auc_values <- list()
  for (model_name in names(models)) {
    model <- models[[model_name]]
    # Assumes "Non_Normal" is the second level in Label factor
    pred_prob <- predict(model, newdata = test_data, type = "prob")[, "Normal"]
    roc_curve <- roc(test_data$Label, pred_prob, levels = rev(levels(test_data$Label)))
    roc_curves[[model_name]] <- roc_curve
    auc_values[[model_name]] <- auc(roc_curve)
  }
  roc_plot <- ggroc(roc_curves, legacy.axes = TRUE) +
    ggtitle("ROC Curves for Different Models") +
    theme_minimal()
  print(roc_plot)
  return(auc_values)
}

auc_results <- plot_roc_curve(models_list, test_norm)
cat("AUC values:\n")
print(auc_results)

# ---------------------------
# Variable Importance (Example: Random Forest)
# ---------------------------
rf_var_imp <- varImp(rf_model)
cat("Variable Importance - Random Forest:\n")
print(rf_var_imp)
plot(rf_var_imp, main = "Variable Importance - Random Forest")


# ---------------------------
# Function to Predict a Single Dataset from a CSV File
# ---------------------------
predict_single_dataset <- function(file_path = NULL, trained_models, preProcValues) {
  # If file_path is not provided, let the user choose interactively
  if(is.null(file_path)){
    cat("Please select the CSV file containing your data...\n")
    file_path <- file.choose()
  }
  # Read the CSV file (assuming header exists)
  user_data <- read.csv(file_path, header = FALSE)
  # If multiple columns, use the first column by default
  if(ncol(user_data) > 1) {
    cat("Multiple columns detected. Using the first column for prediction.\n")
  }
  samples <- user_data[[1]]
  #standardized data
  #standardized_samples <- scale(samples)
  
  # Extract features from the user's dataset
  user_features <- calculate_features(samples)
  
  # Normalize the features using the preProcValues learned from training
  user_features_norm <- predict(preProcValues, user_features)
  
  # Predict using each model and store the results
  predictions <- list()
  for(model_name in names(trained_models)) {
    pred <- predict(trained_models[[model_name]], newdata = user_features_norm)
    predictions[[model_name]] <- pred
  }
  return(predictions)
}

# ---------------------------
# Example Usage for Single Dataset Prediction
# ---------------------------
single_predictions <- predict_single_dataset(file_path = NULL,
                                              trained_models = models_list,
                                              preProcValues = preProcValues)
cat("Prediction Results for the provided dataset:\n")
print(single_predictions)
