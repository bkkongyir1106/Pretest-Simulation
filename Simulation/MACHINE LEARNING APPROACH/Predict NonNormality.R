# Required Libraries
rm(list = ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest, infotheo, ineq, caret, pROC, ROCR, randomForest, evd, discretization, nnet, ggplot2)

# Step 1: Generate Random Samples
generate_nonnormal_samples <- function(n, dist = "exponential") {
  if (dist == "exponential") {
    samples <- rexp(n, rate = 1)
  } else if (dist == "lognormal") {
    samples <- rlnorm(n)
  } else if (dist == "chi_squared") {
    samples <- rchisq(n, df = 2)
  } else if (dist == "gamma") {
    samples <- rgamma(n, shape = 2, rate = 1)
  } else {
    stop("Unknown distribution")
  }
  return(samples)
}

calculate_zero_crossing_rate <- function(samples) {
  zero_crossings <- sum(diff(samples > 0))
  return(zero_crossings / length(samples))
}

calculate_gini_coefficient <- function(samples) {
  gini_value <- ineq(samples, type = "Gini")
  return(gini_value)
}

calculate_features <- function(samples) {
  skewness <- e1071::skewness(samples)
  kurtosis <- e1071::kurtosis(samples)
  jb_stat <- tseries::jarque.bera.test(samples)$statistic
  ad_stat <- nortest::ad.test(samples)$statistic
  zero_crossing_rate <- calculate_zero_crossing_rate(samples)
  gini_coefficient <- calculate_gini_coefficient(samples)
  
  features <- data.frame(
    Skewness = skewness,
    Kurtosis = kurtosis,
    JB_Statistic = as.numeric(jb_stat),
    AD_Statistic = as.numeric(ad_stat),
    Zero_Crossing_Rate = zero_crossing_rate,
    Gini_Coefficient = gini_coefficient
  )
  
  return(features)
}

# Step 2: Prepare Data for Machine Learning Models
generate_data <- function(n_samples, n_per_sample, label, dist = "normal") {
  data <- do.call(rbind, lapply(1:n_samples, function(x) {
    samples <- if (dist == "normal") {
      rnorm(n_per_sample)
    } else if (dist == "normal_25_12") {
      rnorm(n_per_sample, mean = 25, sd = sqrt(12))
    } else {
      generate_nonnormal_samples(n_per_sample, dist)
    }
    features <- calculate_features(samples)
    features$Label <- label
    return(features)
  }))
  return(data)
}

# Generate data
set.seed(123)
normal_data <- rbind(
  generate_data(1000, 8, "Normal", "normal"),
  generate_data(1000, 8, "Normal", "normal_25_12")
)

lognorma <- generate_data(500, 8, "Non_Normal", "lognormal")
chisq_data <- generate_data(500, 8, "Non_Normal", "chi_squared")
exp_data <- generate_data(500, 8, "Non_Normal", "exponential")
gamma_data <- generate_data(500, 8, "Non_Normal", "gamma")
non_normal_data <- rbind(lognorma, chisq_data, exp_data, gamma_data)

# Combine and shuffle data, excluding entropy and tail index columns
data <- rbind(normal_data, non_normal_data)
data <- data[sample(nrow(data)), ]

# Step 3: Split Data into Training and Test Sets
set.seed(123)
trainIndex <- createDataPartition(data$Label, p = .8, list = FALSE, times = 1)
train_data <- data[trainIndex, ]
test_data  <- data[-trainIndex, ]

# Ensure that Label is a factor in both train and test data
train_data$Label <- as.factor(train_data$Label)
test_data$Label <- as.factor(test_data$Label)

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
    pred_prob <- tryCatch({
      predict(model, newdata = test_data, type = "prob")[,2]
    }, error = function(e) {
      warning(paste("Probability prediction failed for model:", model_name))
      return(NULL)
    })
    
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


