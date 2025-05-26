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
    #Coefficient_of_Variation = cv,
    Energy = energy
  )
  return(features)
}

# ---------------------------
# Standardization & Normalization Function
# ---------------------------
preprocess_data <- function(train_data, test_data) {
  numeric_train <- train_data[, sapply(train_data, is.numeric)]
  
  # Step 1: Standardization (Z-score scaling: mean = 0, std = 1)
  preProcStandard <- preProcess(numeric_train, method = c("center", "scale"))
  train_std <- predict(preProcStandard, numeric_train)
  test_std  <- predict(preProcStandard, test_data[, sapply(test_data, is.numeric)])
  
  # Step 2: Normalization (Rescale to [0,1])
  preProcNorm <- preProcess(train_std, method = "range")
  train_norm <- predict(preProcNorm, train_std)
  test_norm  <- predict(preProcNorm, test_std)
  
  train_norm$Label <- train_data$Label
  test_norm$Label  <- test_data$Label
  
  return(list(train = train_norm, test = test_norm, preProcStandard = preProcStandard, preProcNorm = preProcNorm))
}

# ---------------------------
# Data Generation Function 
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
# Set Seed for Reproducibility and Generate Training Data
# ---------------------------
set.seed(12345)
sample_size <- 8  # consistent sample size per dataset
N <- 100            # number of samples per distribution

normal_data1 <- generate_data(sample_size, 5*N, "normal", "Normal")
normal_data2 <- generate_data(sample_size, 5*N, "normal_25", "Normal")
# Generate non-normal data from several distributions
lognormal <- generate_data(sample_size, N, "LogNormal", "Non_Normal")
chisq_data   <- generate_data(sample_size, N, "Chi-Square", "Non_Normal")
exp_data     <- generate_data(sample_size, N, "Exponential", "Non_Normal")
Weibull      <- generate_data(sample_size, N, "Weibull", "Non_Normal")
Pareto      <- generate_data(sample_size, N, "Pareto", "Non_Normal")
Laplace      <- generate_data(sample_size, N, "Laplace", "Non_Normal")
Gamma        <- generate_data(sample_size, N, "Gamma", "Non_Normal")
Uniform      <- generate_data(sample_size, N, "Uniform", "Non_Normal")
t_data       <- generate_data(sample_size, N, "t", "Non_Normal")
beta       <- generate_data(sample_size, N, "beta", "Non_Normal")
t15_data       <- generate_data(sample_size, N, "t_15", "Non_Normal")

non_normal_data <- rbind(lognormal, chisq_data, exp_data, Weibull, Pareto, Laplace, Gamma, Uniform, t_data, t15_data)
data_all <- rbind(normal_data1, normal_data2, non_normal_data)
data_all$Label <- as.factor(data_all$Label)


# ---------------------------
# Split Data (70% Train, 30% Test)
# ---------------------------
set.seed(123)
trainIndex <- createDataPartition(data_all$Label, p = 0.7, list = FALSE)
train_data <- data_all[trainIndex, ]
test_data  <- data_all[-trainIndex, ]

# ---------------------------
# Standardize & Normalize Data
# ---------------------------
norm_result <- preprocess_data(train_data, test_data)
train_norm <- norm_result$train
test_norm  <- norm_result$test

# define appropriate class reference
train_norm$Label <- relevel(train_norm$Label, ref = "Non_Normal")
test_norm$Label  <- relevel(test_norm$Label, ref = "Non_Normal")

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
# Logistic Regression
log_model <- train(Label ~ .,
                   data = train_norm,
                   method = "glm", 
                   family = "binomial", 
                   trControl = ctrl, 
                   metric = "ROC")
log_pred <- predict(log_model, newdata = test_norm)
cat("Test Results for: Logistic Regression Model\n")
print(confusionMatrix(log_pred, test_norm$Label))

# Random Forest
rf_model <- train(Label ~ ., 
                  data = train_norm, 
                  method = "rf", 
                  trControl = ctrl, 
                  metric = "ROC")
rf_pred <- predict(rf_model, newdata = test_norm)
cat("Test Results for: Random Forest Model\n")
print(confusionMatrix(rf_pred, test_norm$Label))

# Artificial Neural Network
ann_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "nnet", 
                   trControl = ctrl, 
                   metric = "ROC", 
                   trace = FALSE)
ann_pred <- predict(ann_model, newdata = test_norm)
cat("Test Results for: ANN Model\n")
print(confusionMatrix(ann_pred, test_norm$Label))

# Gradient Boosting
gbm_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "gbm", 
                   trControl = ctrl, 
                   metric = "ROC", 
                   verbose = FALSE)
gbm_pred <- predict(gbm_model, newdata = test_norm)
cat("Test Results for: GBM Model\n")
print(confusionMatrix(gbm_pred, test_norm$Label))

# Support Vector Machines
svm_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "svmRadial", 
                   trControl = ctrl, 
                   metric = "ROC")
svm_pred <- predict(svm_model, newdata = test_norm)
cat("Test Results for: SVM Model\n")
print(confusionMatrix(svm_pred, test_norm$Label))

# K-Nearest Neighbors
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

# ----------------------------------------
# Variable Importance : Random Forest
# ----------------------------------------
rf_var_imp <- varImp(rf_model)
cat("Variable Importance - Random Forest:\n")
print(rf_var_imp)
plot(rf_var_imp, main = "Variable Importance - Random Forest")


# -----------------------------------------------
# Looping Prediction over Multiple Simulations
# ----------------------------------------------
simulate_predictions <- function(n_iter = 1000, n = 10, dist = "normal",
                                 trained_models, preProcStandard, preProcNorm) {
  # Initialize accuracy counts for each model
  accuracies <- setNames(rep(0, length(trained_models)), names(trained_models))
  # Create a list to store prediction results for each model
  results_list <- lapply(names(trained_models), function(x) {
    data.frame(Iteration = integer(), True_Class = character(), Predicted_Class = character(), stringsAsFactors = FALSE)
  })
  names(results_list) <- names(trained_models)
  
  get_true_class <- function(dist) {
    if (startsWith(dist, "normal")) return("Normal") else return("Non_Normal")
  }

  for (i in 1:n_iter) {
    # generate data
    samples <- generate_samples(n, dist)
    # get true class
    true_class <- get_true_class(dist)
    
    # calculate features and preprocess them
    user_features <- calculate_features(samples)
    user_features_std <- predict(preProcStandard, user_features)
    user_features_norm <- predict(preProcNorm, user_features_std)
    
    # Loop through each model and collect its prediction
    for (model_name in names(trained_models)) {
      pred_class <- as.character(predict(trained_models[[model_name]], newdata = user_features_norm))
      
      # Append result for this model
      results_list[[model_name]] <- rbind(results_list[[model_name]], 
                                          data.frame(Iteration = i, True_Class = true_class, Predicted_Class = pred_class, stringsAsFactors = FALSE))
      
      if (pred_class == true_class) {
        accuracies[model_name] <- accuracies[model_name] + 1
      }
    }
  }
  
  # Calculate accuracies as proportions
  accuracies <- accuracies / n_iter
  
  return(list(
    accuracy = accuracies,
    prediction_table = results_list
  ))
}

# Example usage:
# Simulation for Normal distribution
distribution <- c("normal", "Exponential", "t", "Uniform")
for (j in seq_along(distribution)) {
  dist <- distribution[j]
  result[j]<- simulate_predictions(
    n_iter = 1000,
    n = 10,
    dist = dist,
    trained_models = models_list,
    preProcStandard = norm_result$preProcStandard,
    preProcNorm = norm_result$preProcNorm
  )
}

cat("Overall prediction accuracies for", dist,"distribution:\n")
print(result$accuracy)
print(result$prediction_table)

# Simulation for LogNormal distribution
dist <- "normal_5"
result <- simulate_predictions(
  n_iter = 1000,
  n = 10,
  dist = dist,
  trained_models = models_list,
  preProcStandard = norm_result$preProcStandard,
  preProcNorm = norm_result$preProcNorm
)
cat("Overall prediction accuracies for", dist, "distribution:\n")
print(result$accuracy)
print(result$prediction_table)

  
  
  
  
