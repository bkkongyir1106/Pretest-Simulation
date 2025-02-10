# Required Libraries
rm(list = ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest, lawstat, infotheo, ineq, caret, pROC, ROCR, randomForest, evd, discretization, nnet, ggplot2)

# Clear working environment
rm(list = ls())
# Set directories in local computer
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/MACHINE LEARNING APPROACH")
source("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/MACHINE LEARNING APPROACH/fun.R")

calculate_zero_crossing_rate <- function(samples) {
  zero_crossings <- sum(diff(samples > 0))
  return(zero_crossings / length(samples))
}

calculate_gini_coefficient <- function(samples) {
  gini_value <- ineq(samples, type = "Gini")
  return(gini_value)
}

calculate_outliers <- function(samples) {
  qnt <- quantile(samples, probs = c(0.25, 0.75))
  H <- 1.5 * IQR(samples)
  outliers <- sum(samples < (qnt[1] - H) | samples > (qnt[2] + H))
  return(outliers)
}

calculate_shapiro_wilk <- function(samples) {
  shapiro_stat <- shapiro.test(samples)$statistic
  return(shapiro_stat)
}

# Shapiro-Francia Test
calculate_shapiro_francia <- function(samples) {
  result <- tryCatch({
    lawstat::sf.test(samples)$statistic
  }, error = function(e) {
    warning("Shapiro-Francia test failed")
    return(NA)
  })
  return(as.numeric(result))
}

# Liliefors Test (Kolmogorov-Smirnov for estimated parameters)
calculate_lilliefors <- function(samples) {
  result <- tryCatch({
    nortest::lillie.test(samples)$statistic
  }, error = function(e) {
    warning("Lilliefors test failed")
    return(NA)
  })
  return(as.numeric(result))
}

# Cramer-von Mises Test
calculate_cramer_von_mises <- function(samples) {
  result <- tryCatch({
    nortest::cvm.test(samples)$statistic
  }, error = function(e) {
    warning("Cramer-von Mises test failed")
    return(NA)
  })
  return(as.numeric(result))
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
  
  features <- data.frame(
    Skewness = skewness,
    Kurtosis = kurtosis,
    JB_Statistic = as.numeric(jb_stat),
    AD_Statistic = as.numeric(ad_stat),
    #Zero_Crossing_Rate = zero_crossing_rate,
    Outliers = outliers,
    Shapiro_Wilk = as.numeric(shapiro_wilk),
    Liliefors = lilliefors_stat,
    Cramer_Von_Mises = cramer_von_mises
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
set.seed(123)
normal_data <- generate_data(10, 1500, "normal", "Normal")
lognormal <- generate_data(10, 300, "LogNormal", "Non_Normal")
chisq_data <- generate_data(10, 300, "Chi-Square", "Non_Normal")
exp_data <- generate_data(10, 500, "Exponential", "Non_Normal")
Weibull <- generate_data(10, 500, "Weibull", "Non_Normal")

#combine data
non_normal_data <- rbind(lognormal, chisq_data, Weibull, exp_data)

# Combine and shuffle data, excluding entropy and tail index columns
data <- rbind(normal_data, non_normal_data)
data <- data[sample(nrow(data)), ]

# Step 3: Split Data into Training and Test Sets
set.seed(123)
trainIndex <- createDataPartition(data$Label, p = .8, list = FALSE, times = 1)
train_data <- data[trainIndex, ]
test_data  <- data[-trainIndex, ]

# Ensure that Label is a factor
train_data$Label <- as.factor(train_data$Label)
test_data$Label <- as.factor(test_data$Label)

# Normalize the training and test data
Normalized_train_data <- normalize_data(train_data)
Normalized_test_data <- normalize_data(test_data)

# Step 4: Build Machine Learning Models

# Logistic Regression Model
log_model <- train(Label ~ ., data = train_data, method = "glm", family = "binomial")
log_pred <- predict(log_model, newdata = test_data)
log_pred <- factor(log_pred, levels = levels(test_data$Label))
log_conf_matrix <- confusionMatrix(log_pred, test_data$Label)
print(log_conf_matrix)

# Random Forest Model
rf_model <- train(Label ~ ., data = train_data, method = "rf")
rf_pred <- predict(rf_model, newdata = test_data)
rf_pred <- factor(rf_pred, levels = levels(test_data$Label))
rf_conf_matrix <- confusionMatrix(rf_pred, test_data$Label)
print(rf_conf_matrix)

# ANN Model
ann_model <- train(Label ~ ., data = train_data, method = "nnet", trace = FALSE)
ann_pred <- predict(ann_model, newdata = test_data)
ann_pred <- factor(ann_pred, levels = levels(test_data$Label))
ann_conf_matrix <- confusionMatrix(ann_pred, test_data$Label)
print(ann_conf_matrix)

# Step 5: Create a function to plot ROC curves with error handling
plot_roc_curve <- function(models, test_data) {
  roc_curves <- list()
  
  for (model_name in names(models)) {
    model <- models[[model_name]]
    
    # Predict probabilities with error handling
    pred_prob <- predict(model, newdata = test_data, type = "prob")[,2]
    # If prediction is successful, calculate ROC curve
    if (!is.null(pred_prob)) {
      roc_curve <- roc(test_data$Label, pred_prob, levels = rev(levels(test_data$Label)))
      roc_curves[[model_name]] <- roc_curve
    }
  }
  
  # Plot ROC curves if available
  if (length(roc_curves) > 0) {
    ggroc(roc_curves) +
      ggtitle("ROC Curves for Different Models") +
      theme_minimal()
  } else {
    print("No ROC curves to display.")
  }
}

# Store models in a list
models <- list(
  "Logistic Regression" = log_model,
  "Random Forest" = rf_model,
  "ANN" = ann_model
)

# Plot the ROC curves
plot_roc_curve(models, test_data)

# Step 6: Variable Importance Plot for Random Forest
rf_var_imp <- varImp(rf_model)
print(rf_var_imp)

# Plot the Variable Importance
plot(rf_var_imp, main = "Variable Importance - Random Forest")


# Generate validation data
set.seed(123)
normal <- generate_data(10, 100, "normal", "Normal")
lognormal <- generate_data(10, 100,  "LogNormal", "Non_Normal")
Uniform <- generate_data(10, 100, "Uniform", "Non_Normal")
Weibull <- generate_data(10, 100, "Weibull", "Non_Normal")

validation_data <- rbind(normal, lognormal)
validation_data <- validation_data[sample(nrow(validation_data)), ]
validation_data <- normalize_data(validation_data)

validation_data$Label <- factor(validation_data$Label, levels = levels(train_data$Label))

# Predict on validation data using trained models
log_pred_val <- predict(log_model, newdata = validation_data)
rf_pred_val <- predict(rf_model, newdata = validation_data)
ann_pred_val <- predict(ann_model, newdata = validation_data)

# Convert predictions to factors with the same levels as the true labels
log_pred_val <- factor(log_pred_val, levels = levels(validation_data$Label))
rf_pred_val <- factor(rf_pred_val, levels = levels(validation_data$Label))
ann_pred_val <- factor(ann_pred_val, levels = levels(validation_data$Label))

# Confusion Matrix for validation data
log_conf_matrix_val <- confusionMatrix(log_pred_val, validation_data$Label)
rf_conf_matrix_val <- confusionMatrix(rf_pred_val, validation_data$Label)
ann_conf_matrix_val <- confusionMatrix(ann_pred_val, validation_data$Label)

# Print the validation confusion matrices
print(log_conf_matrix_val)
print(rf_conf_matrix_val)
print(ann_conf_matrix_val)

