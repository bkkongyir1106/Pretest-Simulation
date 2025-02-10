# Required Libraries
if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071,tseries, nortest,infotheo, ineq, caret, pROC,ROCR, randomForest,evd, discretization)

# Step 1: Generate Random Samples
generate_nonnormal_samples <- function(n, dist = "exponential") {
  if (dist == "exponential") {
    samples <- rexp(n, rate = 1)
  } else if (dist == "lognormal") {
    samples <- rlnorm(n)
  } else if (dist == "chi_squared") {
    samples <- rchisq(n, df = 2)
  } else {
    stop("Unknown distribution")
  }
  return(samples)
}

calculate_entropy <- function(samples) {
  discretized <- discretize(samples, disc = "equalfreq", nbins = 10)
  entropy_value <- entropy(discretized$X)  # Extract the relevant vector
  return(entropy_value)
}

calculate_zero_crossing_rate <- function(samples) {
  zero_crossings <- sum(diff(samples > 0))
  return(zero_crossings / length(samples))
}

calculate_gini_coefficient <- function(samples) {
  gini_value <- ineq(samples, type = "Gini")
  return(gini_value)
}

calculate_tail_index <- function(samples) {
  tail_index <- gpd.fit(samples)$par.ests[1]  # Shape parameter is the tail index
  return(tail_index)
}

calculate_features <- function(samples) {
  skewness <- e1071::skewness(samples)
  kurtosis <- e1071::kurtosis(samples)
  jb_stat <- tseries::jarque.bera.test(samples)$statistic
  ad_stat <- nortest::ad.test(samples)$statistic
  entropy_value <- calculate_entropy(samples)
  zero_crossing_rate <- calculate_zero_crossing_rate(samples)
  gini_coefficient <- calculate_gini_coefficient(samples)
  tail_index <- calculate_tail_index(samples)
  
  features <- data.frame(
    Skewness = skewness,
    Kurtosis = kurtosis,
    JB_Statistic = as.numeric(jb_stat),
    AD_Statistic = as.numeric(ad_stat),
    Entropy = entropy_value,
    Zero_Crossing_Rate = zero_crossing_rate,
    Gini_Coefficient = gini_coefficient,
    Tail_Index = tail_index
  )
  
  return(features)
}

# Step 2: Prepare Data for Machine Learning Models
generate_data <- function(n_samples, n_per_sample, label, dist = "normal") {
  data <- do.call(rbind, lapply(1:n_samples, function(x) {
    samples <- if (dist == "normal") {
      rnorm(n_per_sample)
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
normal_data <- generate_data(10000, 10, "Normal", "normal")
non_normal_data <- generate_data(1000, 10, "Non-Normal", "exponential")  # Modify distribution as needed

# Combine and shuffle data
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

# SVM Model
svm_model <- train(Label ~ ., data = train_data, method = "svmRadial")
svm_pred <- predict(svm_model, newdata = test_data)
svm_pred <- factor(svm_pred, levels = levels(test_data$Label))
svm_conf_matrix <- confusionMatrix(svm_pred, test_data$Label)
print(svm_conf_matrix)

# ANN Model
ann_model <- train(Label ~ ., data = train_data, method = "nnet", trace = FALSE)
ann_pred <- predict(ann_model, newdata = test_data)
ann_pred <- factor(ann_pred, levels = levels(test_data$Label))
ann_conf_matrix <- confusionMatrix(ann_pred, test_data$Label)
print(ann_conf_matrix)

# Step 5: Create Combined ROC Curve
plot_roc_curve <- function(models, test_data) {
  roc_curves <- list()
  
  for (model_name in names(models)) {
    model <- models[[model_name]]
    pred <- predict(model, newdata = test_data, type = "prob")[,2]
    roc_curve <- roc(test_data$Label, pred, levels = rev(levels(test_data$Label)))
    roc_curves[[model_name]] <- roc_curve
  }
  
  # Plot ROC curves
  ggroc(roc_curves) +
    ggtitle("ROC Curves for Different Models") +
    theme_minimal()
}

# Store models in a list
models <- list(
  "Logistic Regression" = log_model,
  "Random Forest" = rf_model,
  "SVM" = svm_model,
  "ANN" = ann_model
)

# Plot the ROC curves
plot_roc_curve(models, test_data)