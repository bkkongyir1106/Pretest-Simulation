rm(list = ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest, gbm, lawstat, infotheo, ineq, caret, pROC, ROCR, randomForest, evd, discretization, nnet, ggplot2)

# ---------------------------
# Set Directories & Load Functions
# ---------------------------
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/ML Models/version 5")
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
    Energy = energy,
    mean_val = mean(samples),
    sd = sd(samples)
  )
  return(features)
}

# ---------------------------
# Standardization & Normalization Function
# ---------------------------
preprocess_data <- function(train_data) {
  numeric_train <- train_data[, sapply(train_data, is.numeric)]
  
  # Step 1: Standardization (Z-score scaling: mean = 0, std = 1)
  preProcStandard <- preProcess(numeric_train, method = c("center", "scale"))
  train_std <- predict(preProcStandard, numeric_train)
  
  # Step 2: Normalization (Rescale to [0,1])
  preProcNorm <- preProcess(train_std, method = "range")
  train_norm <- predict(preProcNorm, train_std)
  
  train_norm$Label <- train_data$Label
  
  return(list(train = train_norm, preProcStandard = preProcStandard, preProcNorm = preProcNorm))
}

# ---------------------------
# Data Generation Function (balanced sample size across groups)
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
sample_size <- 8  
N <- 100            

normal_data1 <- generate_data(sample_size, 4*N, "normal", "Normal")
normal_data2 <- generate_data(sample_size, 3*N, "normal_5", "Normal")
normal_data3 <- generate_data(sample_size, 3*N, "normal_100", "Normal")
# Generate non-normal data from several distributions
lognormal <- generate_data(sample_size, N, "LogNormal", "Non_Normal")
chisq_data   <- generate_data(sample_size, N, "Chi-Square", "Non_Normal")
exp_data     <- generate_data(sample_size, N, "Exponential", "Non_Normal")
Weibull      <- generate_data(sample_size, N, "Weibull", "Non_Normal")
Pareto      <- generate_data(sample_size, N, "Pareto", "Non_Normal")
Laplace      <- generate_data(sample_size, N, "Laplace", "Non_Normal")
Gamma        <- generate_data(sample_size, N, "Gamma", "Non_Normal")
Uniform      <- generate_data(sample_size, N, "Uniform", "Non_Normal")
t       <- generate_data(sample_size, N, "t", "Non_Normal")
t_15       <- generate_data(sample_size, N, "t_15", "Non_Normal")
t_25       <- generate_data(sample_size, N, "t_25", "Non_Normal")
beta       <- generate_data(sample_size, N, "beta", "Non_Normal")

non_normal_data <- rbind(lognormal, chisq_data,  exp_data, Weibull, Pareto, Laplace, Gamma, Uniform, t, t_25)
normal_data <- rbind(normal_data1, normal_data2, normal_data3)

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
# ---------------------------
# Define Common Training Control for Cross-Validation
set.seed(12345)
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


# Extract and display the tuning parameters for each model
list(
  Logistic_Regression = log_model$bestTune,
  Random_Forest = rf_model$bestTune,
  Neural_Network = ann_model$bestTune,
  Gradient_Boosting = gbm_model$bestTune,
  Support_Vector_Machine = svm_model$bestTune,
  K_Nearest_Neighbors = knn_model$bestTune
)

# Extract accuracy for each fold
list(
  Logistic_Regression = log_model$results$Accuracy,
  Random_Forest = rf_model$results$Accuracy,
  Neural_Network = ann_model$results$Accuracy,
  Gradient_Boosting = gbm_model$results$Accuracy,
  Support_Vector_Machine = svm_model$results$Accuracy,
  K_Nearest_Neighbors = knn_model$results$Accuracy
)



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
simulate_predictions <- function(distributions ,n_iter = 1000, n = 8,
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

set.seed(12345)

distribution_set <- c("normal", "normal_15", "beta", "Chi-Square")
eval_results <- simulate_predictions(
  distributions = distribution_set,
  n_iter = 1000,
  n = 8,
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
