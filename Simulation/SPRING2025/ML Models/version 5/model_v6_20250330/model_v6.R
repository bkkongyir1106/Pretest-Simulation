rm(list = ls())
# ---------------------------
# Set Directories & Load Functions
# ---------------------------
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/ML Models/version 5")
source("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/ML Models/fun.R")

if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest, gbm, lawstat, infotheo, ineq, caret, pROC, ROCR, randomForest, evd, discretization, nnet, ggplot2)

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

calculate_entropy <- function(samples) {
  return(infotheo::entropy(discretize(samples, nbins = 10)))
}

calculate_features <- function(samples) {
  skewness <- e1071::skewness(samples)
  kurtosis <- e1071::kurtosis(samples)
  mean_val <-  mean(samples)
  median_val  <- median(samples)
  variance <- var(samples)
  calculate_iqr  <- IQR(samples)
  calculate_mad   <- mad(samples)
  calculate_rms   <- sqrt(mean(samples^2))
  range <- max(samples) - min(samples)
  cv <- sd(samples) / mean(samples)
  energy <- sum(samples^2)
  
  # Statistical tests for normality
  shapiro_wilk  <- as.numeric(shapiro.test(samples)$statistic)
  shapiro_francia <- as.numeric(nortest::sf.test(samples)$statistic)
  lilliefors_stat <- as.numeric(nortest::lillie.test(samples)$statistic)
  cramer_von_mises <-  as.numeric(nortest::cvm.test(samples)$statistic)
  jb_stat <- as.numeric(tseries::jarque.bera.test(samples)$statistic)
  ad_stat <- as.numeric(nortest::ad.test(samples)$statistic)

  zero_crossing_rate <- calculate_zero_crossing_rate(samples)
  gini_coefficient <- calculate_gini_coefficient(samples)
  outliers <- calculate_outliers(samples)
  peak_to_trough <- max(samples) / abs(min(samples))
  entropy_rate <- calculate_entropy(samples)
  
  features <- data.frame(
    #JB_Statistics incorporate both skewness and kurtosis, so we can exclude them
    Skewness = skewness,
    Kurtosis = kurtosis,
    mean_val = mean_val,
    median_val = median_val,
    variance = variance,
    calculate_iqr = calculate_iqr,
    calculate_mad = calculate_mad,
    calculate_rms = calculate_rms,
    range = range,
    Energy = energy,
    sd = sd(samples),
    Coefficient_of_Variation = cv,
    Shapiro_Wilk = shapiro_wilk,
    shapiro_francia = shapiro_francia,
    Liliefors = lilliefors_stat,
    Cramer_Von_Mises = cramer_von_mises,
    JB_Statistic = jb_stat,
    AD_Statistic = ad_stat,
    Zero_Crossing_Rate = zero_crossing_rate,
    Outliers = outliers,
    peak_to_trough = peak_to_trough,
    gini_coefficient = gini_coefficient,
    entropy_rate = entropy_rate
  
  )
  return(features)
}

# ------------------------------------------
# Standardization & Normalization Function
# ------------------------------------------
preprocess_data <- function(train_data) {
  numeric_train <- train_data[, sapply(train_data, is.numeric)]
  
  # Step 1: Standardization (Z-score scaling: mean = 0, std = 1)
  preProcStandard <- preProcess(numeric_train, method = c("center", "scale"))
  train_std <- predict(preProcStandard, numeric_train)
  
  # Step 2: Normalization (Rescale to [0,1])
  preProcNorm <- preProcess(train_std, method = "range")
  train_norm <- predict(preProcNorm, train_std)
  # add labels back
  train_norm$Label <- train_data$Label
  
  return(list(train = train_norm, preProcStandard = preProcStandard, preProcNorm = preProcNorm))
}

# ---------------------------
# Data Generation Function 
# ---------------------------
generate_data <- function(sample_size, num.sim, dist = "normal", label) {
  data <- do.call(rbind, lapply(1:num.sim, function(x) {
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
sample_size <- 10 
num.sim <- 100            

# Normal data
normal_data1 <- generate_data(sample_size, 2*num.sim, "normal", "Normal")
normal_data2 <- generate_data(sample_size, 2*num.sim, "normal_5", "Normal")
normal_data3 <- generate_data(sample_size, 2*num.sim, "normal_15", "Normal")
normal_data4 <- generate_data(sample_size, 2*num.sim, "normal_25", "Normal")
normal_data5 <- generate_data(sample_size, 2*num.sim, "normal_50", "Normal")
normal_data6 <- generate_data(sample_size, 2*num.sim, "normal_100", "Normal")
# non-normal data 
lognormal <- generate_data(sample_size, num.sim, "LogNormal", "Non_Normal")
chisq_data   <- generate_data(sample_size, num.sim, "Chi_Square", "Non_Normal")
exp_data     <- generate_data(sample_size, num.sim, "Exponential", "Non_Normal")
Weibull      <- generate_data(sample_size, num.sim, "Weibull", "Non_Normal")
Pareto      <- generate_data(sample_size, num.sim, "Pareto", "Non_Normal")
Laplace      <- generate_data(sample_size, num.sim, "Laplace", "Non_Normal")
Gamma        <- generate_data(sample_size, num.sim, "Gamma", "Non_Normal")
Uniform      <- generate_data(sample_size, num.sim, "Uniform", "Non_Normal")
t       <- generate_data(sample_size, num.sim, "t", "Non_Normal")
t_15       <- generate_data(sample_size, num.sim, "t_5", "Non_Normal")
t_25       <- generate_data(sample_size, num.sim, "t_15", "Non_Normal")
beta       <- generate_data(sample_size, num.sim, "beta", "Non_Normal")

non_normal_data <- rbind(lognormal, chisq_data,  exp_data, Weibull, Pareto, Laplace, Gamma, Uniform, t, t_15)
normal_data <- rbind(normal_data1, normal_data3, normal_data4, normal_data5, normal_data6)

data_all <- rbind(normal_data, non_normal_data)
data_all$Label <- as.factor(data_all$Label)

# ---------------------------
# Standardize & Normalize Data
# ---------------------------
norm_result <- preprocess_data(data_all)
train_norm <- norm_result$train

# define appropriate class reference
train_norm$Label <- relevel(train_norm$Label, ref = "Non_Normal")

### Train the models
# ---------------------------------------------------
# Define Common Training Control for Cross-Validation
set.seed(12345)
# ---------------------------------------------------
ctrl <- trainControl(method = "cv", number = 10, 
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE, 
                     savePredictions = "final")


# Logistic Regression
log_model <- train(Label ~ .,
                   data = train_norm,
                   method = "glm", 
                   family = "binomial", 
                   trControl = ctrl, 
                   metric = "ROC")

# Random Forest
rf_model <- train(Label ~ ., 
                  data = train_norm, 
                  method = "rf", 
                  trControl = ctrl, 
                  metric = "ROC")

# Artificial Neural Network
ann_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "nnet", 
                   trControl = ctrl, 
                   metric = "ROC", 
                   trace = FALSE)

# Gradient Boosting
gbm_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "gbm", 
                   trControl = ctrl, 
                   metric = "ROC", 
                   verbose = FALSE)

# Support Vector Machines
svm_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "svmRadial", 
                   trControl = ctrl, 
                   metric = "ROC")

# K-Nearest Neighbors
knn_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "knn", 
                   trControl = ctrl, 
                   metric = "ROC")

# Collect all models into a list for later use
models_list <- list(
  "Logistic Regression" = log_model,
  "Random Forest" = rf_model,
  "ANN" = ann_model,
  "GBM" = gbm_model,
  "SVM" = svm_model,
  "KNN" = knn_model
)

### Predictions on unseen data row - by - row
simulate_predictions <- function(distributions ,n_iter = 1000, n = sample_size,
                                 trained_models, preProcStandard, preProcNorm) {
  
  # Initialize storage for all model predictions
  all_predictions <- lapply(names(trained_models), function(x) {
    data.frame(True_Class = character(),
               Predicted_Class = character(),
               stringsAsFactors = FALSE)
  })
  names(all_predictions) <- names(trained_models)
  # get true class
  get_true_class <- function(dist) {
    if (startsWith(dist, "normal")) return("Normal") else return("Non_Normal")
  }
  # predict for each sample data individually
  for (dist in distributions) {
    for (i in 1:n_iter) {
      # generate data and identify the class
      samples <- generate_samples(n, dist)
      true_class <- get_true_class(dist)
      
      # Calculate and normalize features
      user_features <- calculate_features(samples)
      user_features_std <- predict(preProcStandard, user_features)
      user_features_norm <- predict(preProcNorm, user_features_std)
      # make predictions on all models
      for (model_name in names(trained_models)) {
        # Predicted class
        pred_class <- as.character(predict(trained_models[[model_name]], newdata = user_features_norm))
        
        # Predicted probabilities
        pred_probs <- predict(trained_models[[model_name]], newdata = user_features_norm, type = "prob")
        
        # Probability of class "Non_Normal"
        prob_non_normal <- pred_probs[, "Non_Normal"]
        
        # Store predictions along with probability
        all_predictions[[model_name]] <- rbind(
          all_predictions[[model_name]],
          data.frame(
            True_Class = true_class,
            Predicted_Class = pred_class,
            Prob_Non_Normal = prob_non_normal,
            stringsAsFactors = FALSE
          )
        )
      }
      
    }
  }
  
  # Compute confusion matrices and metrics for each model
  evaluation_results <- list()
  for (model_name in names(all_predictions)) {
    
    pred_df <- all_predictions[[model_name]]
    
    # Convert to factors with consistent levels
    pred_df$True_Class <- factor(pred_df$True_Class, levels = c("Non_Normal", "Normal"))
    pred_df$Predicted_Class <- factor(pred_df$Predicted_Class, levels = c("Non_Normal", "Normal"))
    
    cm <- confusionMatrix(pred_df$Predicted_Class, pred_df$True_Class, positive = "Non_Normal")
    
    # store metrics and the actual + predicted values
    evaluation_results[[model_name]] <- list(
      ConfusionMatrix = cm,
      Accuracy = cm$overall['Accuracy'],
      Sensitivity = cm$byClass['Sensitivity'],
      Specificity = cm$byClass['Specificity'],
      Precision = cm$byClass['Precision'],
      F1 = cm$byClass['F1'],
      Predictions = pred_df  
    )
  }
  
  return(evaluation_results)
}
# --------------------------------------------
#             Example predictions
# --------------------------------------------
set.seed(12345)

distribution_set <- c("normal", "normal_15", "t_10", "beta")
eval_results <- simulate_predictions(
  distributions = distribution_set,
  n_iter = 1000,
  n = sample_size,
  trained_models = models_list,
  preProcStandard = norm_result$preProcStandard,
  preProcNorm = norm_result$preProcNorm
)

# View results
for (model in names(eval_results)) {
  cat("\n==== Model:", model, "====\n")
  print(eval_results[[model]]$ConfusionMatrix)
  cat("Accuracy:", round(eval_results[[model]]$Accuracy, 4), "\n")
  cat("Sensitivity:", round(eval_results[[model]]$Sensitivity, 4), "\n")
  cat("Specificity:", round(eval_results[[model]]$Specificity, 4), "\n")
  cat("Precision:", round(eval_results[[model]]$Precision, 4), "\n")
  cat("F1 Score:", round(eval_results[[model]]$F1, 4), "\n")
}


# Export predictions for each model to CSV
export_predictions <- function(eval_results, output_dir = "predictions_output") {
  if (!dir.exists(output_dir)) dir.create(output_dir)
  
  for (model_name in names(eval_results)) {
    pred_df <- eval_results[[model_name]]$Predictions
    write.csv(pred_df,
              file = file.path(output_dir, paste0(model_name, "_predictions.csv")),
              row.names = FALSE)
  }
}

# export results
export_predictions(eval_results)

# Extract variable importance
rf_var_imp <- varImp(models_list$`Random Forest`)

# Convert to data frame for ggplot
rf_var_imp_df <- data.frame(Variable = rownames(rf_var_imp$importance), 
                            Importance = rf_var_imp$importance$Overall)

# Plot as a proper bar chart
ggplot(rf_var_imp_df, aes(x = reorder(Variable, Importance), y = Importance)) +
  geom_bar(stat = "identity", fill = "steelblue") +
  coord_flip() +  # Flip for better readability
  labs(title = "Variable Importance - Random Forest",
       x = "Variables",
       y = "Importance") +
  theme_minimal()

### Plot ROC Curves
plot_combined_roc_ggplot <- function(eval_results) {
  roc_data <- data.frame()

  for (model_name in names(eval_results)) {
    pred_df <- eval_results[[model_name]]$Predictions

    # Convert actual class to binary: 1 = Non_Normal, 0 = Normal
    actual_bin <- ifelse(pred_df$True_Class == "Non_Normal", 1, 0)

    # Predicted probability for class "Non_Normal"
    probs <- pred_df$Prob_Non_Normal

    # Create ROC object with explicit positive class
    roc_obj <- roc(actual_bin, probs, levels = c(0,1), direction = "<", quiet = TRUE)

    # Extract ROC curve points
    roc_points <- data.frame(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Model = model_name,
      AUC = round(auc(roc_obj), 3)
    )

    roc_data <- rbind(roc_data, roc_points)
  }

  # Plot all ROC curves together
  ggplot(roc_data, aes(x = FPR, y = TPR, color = Model)) +
    geom_line(size = 1.2) +
    geom_abline(linetype = "dashed", color = "gray") +
    labs(
      title = "Combined ROC Curves (Positive Class: Non_Normal)",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.title = element_blank()) +
    scale_color_brewer(palette = "Set1") +
    facet_wrap(~ Model + AUC, ncol = 2, scales = "free_y")
}

plot_combined_roc_ggplot(eval_results)


### Summary Table of All Metrics
get_metrics_summary <- function(eval_results) {
  summary_df <- data.frame(
    Model = character(),
    Accuracy = numeric(),
    Sensitivity = numeric(),
    Specificity = numeric(),
    Precision = numeric(),
    F1 = numeric(),
    stringsAsFactors = FALSE
  )
  
  for (model_name in names(eval_results)) {
    metrics <- eval_results[[model_name]]
    summary_df <- rbind(summary_df, data.frame(
      Model = model_name,
      Accuracy = round(metrics$Accuracy, 4),
      Sensitivity = round(metrics$Sensitivity, 4),
      Specificity = round(metrics$Specificity, 4),
      Precision = round(metrics$Precision, 4),
      F1 = round(metrics$F1, 4)
    ))
  }
  return(summary_df)
}

# usage
get_metrics_summary(eval_results)

### 
plot_combined_roc_ggplot <- function(eval_results) {
 
  roc_data <- data.frame()
  
  for (model_name in names(eval_results)) {
    pred_df <- eval_results[[model_name]]
    
    # Check and skip if Prob_Non_Normal is missing
    if (!"Prob_Non_Normal" %in% colnames(pred_df$Predictions)) next
    
    actual_bin <- ifelse(pred_df$Predictions$True_Class == "Non_Normal", 1, 0)
    probs <- pred_df$Predictions$Prob_Non_Normal
    
    roc_obj <- roc(actual_bin, probs, levels = c(0, 1), direction = "<", quiet = TRUE)
    auc_value <- round(auc(roc_obj), 3)
    
    roc_points <- data.frame(
      FPR = 1 - roc_obj$specificities,
      TPR = roc_obj$sensitivities,
      Model = paste0(model_name, " (AUC=", auc_value, ")")
    )
    
    roc_data <- rbind(roc_data, roc_points)
  }
  
  ggplot(roc_data, aes(x = FPR, y = TPR, color = Model, linetype = Model)) +
    geom_line(size = 1.2) +
    geom_abline(linetype = "dashed", color = "gray") +
    labs(
      title = "Combined ROC Curves (Positive Class: Non_Normal)",
      x = "False Positive Rate (1 - Specificity)",
      y = "True Positive Rate (Sensitivity)"
    ) +
    theme_minimal(base_size = 14) +
    theme(legend.title = element_blank())
}
plot_combined_roc_ggplot(eval_results)

