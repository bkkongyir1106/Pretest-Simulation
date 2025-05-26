# --------------------------------------------------------------------------------
# Project: Machine Learning Approach for Simulation Data (With File Browser Validation)
# --------------------------------------------------------------------------------

# ---- Setup Environment ----
rm(list = ls())
if (!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest, gbm, lawstat, infotheo, ineq,
               caret, pROC, ROCR, randomForest, evd, discretization, nnet, 
               ggplot2)
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/MACHINE LEARNING APPROACH")
source("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/ML Models/fun.R")

# ---- Feature Engineering Functions ----
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
  peak <- max(samples)
  trough <- min(samples)
  return(ifelse(trough != 0, peak/abs(trough), NA))
}

calculate_range <- function(samples) {
  return(max(samples) - min(samples))
}

calculate_cv <- function(samples) {
  m <- mean(samples)
  return(ifelse(m != 0, sd(samples)/m, NA))
}

calculate_energy <- function(samples) {
  return(sum(samples^2))
}

calculate_features <- function(samples) {
  features <- data.frame(
    Skewness = e1071::skewness(samples, na.rm = TRUE),
    Kurtosis = e1071::kurtosis(samples, na.rm = TRUE),
    JB_Statistic = as.numeric(tseries::jarque.bera.test(samples)$statistic),
    AD_Statistic = as.numeric(nortest::ad.test(samples)$statistic),
    Zero_Crossing_Rate = calculate_zero_crossing_rate(samples),
    Gini = calculate_gini_coefficient(samples),
    Outliers = calculate_outliers(samples),
    Shapiro_Wilk = calculate_shapiro_wilk(samples),
    Range = calculate_range(samples),
    CV = calculate_cv(samples),
    Energy = calculate_energy(samples)
  )
  return(features)
}
# ---- Data Processing Functions ----
compute_norm_params <- function(train_data) {
  numeric_cols <- sapply(train_data, is.numeric)
  lapply(train_data[, numeric_cols, drop = FALSE], function(x) {
    list(min = min(x, na.rm = TRUE), 
         max = max(x, na.rm = TRUE))
  })
}

apply_normalization <- function(data, params) {
  normalized <- data
  for(col in names(params)) {
    if(col %in% colnames(normalized)) {
      normalized[[col]] <- (data[[col]] - params[[col]]$min) / 
        (params[[col]]$max - params[[col]]$min)
    }
  }
  normalized
}

# ---- File Browser Integration ----
browse_file <- function() {
  if (requireNamespace("tcltk", quietly = TRUE)) {
    file_path <- tcltk::tk_choose.files(
      caption = "Select Validation Dataset",
      filters = matrix(c("CSV files", ".csv", "All files", "*"), 2, 2, byrow = TRUE)
    )
  } else {
    file_path <- file.choose()
  }
  
  if (length(file_path) == 0 || file_path == "") {
    message("File selection cancelled by user")
    return(NULL)
  }
  return(file_path)
}

load_external_data <- function() {
  file_path <- browse_file()
  if (is.null(file_path)) return(NULL)
  
  raw_data <- tryCatch(
    read.csv(file_path, header = FALSE, stringsAsFactors = FALSE),
    error = function(e) {
      message("Error reading file: ", e$message)
      return(NULL)
    }
  )
  
  # Extract labels (last column) and samples (all other columns)
  labels <- raw_data[, ncol(raw_data)]
  samples_data <- raw_data[, -ncol(raw_data)]
  
  # Validate structure
  if (ncol(samples_data) < 1) {
    stop("CSV must contain at least 1 observation column + 1 label column")
  }
  
  if (!all(sapply(samples_data, is.numeric))) {
    stop("All observation columns must be numeric")
  }
  
  labels <- trimws(labels)
  if (!all(labels %in% c("Normal", "Non_Normal"))) {
    stop("Labels must be 'Normal' or 'Non_Normal'")
  }
  
  # Process each row as a sample
  features_list <- lapply(1:nrow(raw_data), function(i) {
    samples <- as.numeric(samples_data[i, ])
    if (length(samples) < 8) {
      stop("Sample ", i, " has fewer than 8 observations")
    }
    features <- calculate_features(samples)
    features$Label <- labels[i]
    features
  })
  
  do.call(rbind, features_list)
}

# ---- Core Analysis Workflow ----
set.seed(12345)

# Generate training data (10 observations per sample)
generate_data <- function(samples_size, N, dist = "normal", label) {
  do.call(rbind, lapply(1:N, function(x) {
    samples <- generate_samples(samples_size, dist)
    features <- calculate_features(samples)
    features$Label <- label
    features
  }))  # Added missing closing parentheses
}

normal_data <- generate_data(10, 1000, "normal", "Normal")
non_normal_data <- rbind(
  generate_data(10, 200, "LogNormal", "Non_Normal"),
  generate_data(10, 200, "Chi-Square", "Non_Normal"),
  generate_data(10, 200, "Exponential", "Non_Normal"),
  generate_data(10, 200, "Weibull", "Non_Normal"),
  generate_data(10, 200, "Laplace", "Non_Normal"),
  generate_data(10, 200, "Gamma", "Non_Normal")  # Fixed missing parentheses
)

combined_data <- rbind(normal_data, non_normal_data)
combined_data$Label <- factor(combined_data$Label, levels = c("Normal", "Non_Normal"))

# Split and normalize data
trainIndex <- createDataPartition(combined_data$Label, p = 0.7, list = FALSE)
train_data <- combined_data[trainIndex, ]
test_data <- combined_data[-trainIndex, ]

train_features <- train_data[, sapply(train_data, is.numeric)]
norm_params <- compute_norm_params(train_features)

Normalized_train_data <- cbind(
  apply_normalization(train_features, norm_params),
  Label = factor(train_data$Label)
)

Normalized_test_data <- cbind(
  apply_normalization(test_data[, sapply(test_data, is.numeric)], norm_params),
  Label = factor(test_data$Label)
)

# ---- Model Training ----
train_control <- trainControl(method = "cv", number = 5, classProbs = TRUE)

log_model <- train(Label ~ ., data = Normalized_train_data, 
                   method = "glm", family = "binomial",
                   trControl = train_control)

rf_model <- train(Label ~ ., data = Normalized_train_data,
                  method = "rf", trControl = train_control)

ann_model <- train(Label ~ ., data = Normalized_train_data,
                   method = "nnet", trace = FALSE,
                   trControl = train_control)

gbm_model <- train(Label ~ ., data = Normalized_train_data,
                   method = "gbm", verbose = FALSE,
                   trControl = train_control)

svm_model <- train(Label ~ ., data = Normalized_train_data,
                   method = "svmRadial",
                   trControl = train_control)

knn_model <- train(Label ~ ., data = Normalized_train_data,
                   method = "knn", trControl = train_control)

# ---- Validation Functions ----
validate_external_dataset <- function(model) {
  cat("\n=== Starting External Validation ===\n")
  cat("Please select validation dataset via file browser...\n")
  raw_features <- load_external_data()
  if (is.null(raw_features)) return(invisible(NULL))
  
  tryCatch({
    normalized_features <- cbind(
      apply_normalization(
        raw_features[, sapply(raw_features, is.numeric)], 
        norm_params
      ),
      Label = factor(raw_features$Label, levels = c("Normal", "Non_Normal"))
    )
    
    predictions <- predict(model, newdata = normalized_features)
    results <- confusionMatrix(predictions, normalized_features$Label)
    
    cat("\nValidation Results:\n")
    print(results$table)
    cat(sprintf("\nAccuracy: %.3f (95%% CI: %.3f-%.3f)\n",
                results$overall["Accuracy"],
                results$overall["AccuracyLower"],
                results$overall["AccuracyUpper"]))
    return(invisible(results))
    
  }, error = function(e) {
    cat("\nValidation failed:", e$message, "\n")
    return(invisible(NULL))
  })
}

# ---- Evaluation and Visualization ----
models <- list(
  "Logistic Regression" = log_model,
  "Random Forest" = rf_model,
  "ANN" = ann_model,
  "GBM" = gbm_model,
  "SVM" = svm_model,
  "KNN" = knn_model
)

# Test set performance
cat("\n=== Test Set Performance ===\n")
lapply(names(models), function(name) {
  cat("\n", name, ":\n", sep = "")
  pred <- predict(models[[name]], Normalized_test_data)
  print(confusionMatrix(pred, Normalized_test_data$Label))
})

# ROC Analysis
plot_roc_curves <- function(models, test_data) {
  roc_data <- lapply(names(models), function(name) {
    pred_prob <- predict(models[[name]], newdata = test_data, type = "prob")[,2]
    roc_obj <- roc(test_data$Label ~ pred_prob)
    data.frame(
      Model = name,
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities
    )
  })
  
  ggplot(do.call(rbind, roc_data), aes(FPR, TPR, color = Model)) +
    geom_line(size = 1) + 
    geom_abline(slope = 1, linetype = "dashed", alpha = 0.5) +
    labs(title = "ROC Curves", x = "False Positive Rate", y = "True Positive Rate") +
    theme_minimal() +
    theme(legend.position = "bottom")
}

plot_roc_curves(models, Normalized_test_data)

# Variable Importance
plot(varImp(rf_model), main = "Random Forest Feature Importance")

# ---- Execution ----
validate_external_dataset(rf_model)  # Triggers file browser

