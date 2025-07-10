source("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/Summer 2025/ML Approach/gen_data_fun.R")

if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest, gbm, lawstat, infotheo, ineq, caret, pROC, ROCR, randomForest,
               evd, discretization, nnet, ggplot2, mlbench, infotheo, dplyr, fractaldim)

# ---------------------------
# Define Optimized Feature Extraction Functions
# ---------------------------

# Core Statistical Features
calculate_core_stats <- function(samples) {
  list(
    Mean = mean(samples),
    Variance = var(samples),
    Skewness = e1071::skewness(samples),
    Kurtosis = e1071::kurtosis(samples),
    IQR = IQR(samples),
    MAD = mad(samples)
  )
}

# Normality Tests (keeping all classical tests)
calculate_normality_tests <- function(samples) {
  list(
    Jarque_Bera = as.numeric(tseries::jarque.bera.test(samples)$statistic),
    Anderson_Darling = as.numeric(nortest::ad.test(samples)$statistic),
    Shapiro_Wilk = shapiro.test(samples)$statistic,
    Shapiro_Francia = nortest::sf.test(samples)$statistic,
    Lilliefors = nortest::lillie.test(samples)$statistic,
    Cramer_Von_Misse = nortest::cvm.test(samples)$statistic
  )
}

# Information-Theoretic Features
calculate_information_features <- function(samples) {
  list(
    Entropy = infotheo::entropy(discretize(samples, nbins = 10)),
    Spectral_Entropy = {
      spec <- stats::spec.pgram(samples, plot = FALSE)
      p <- spec$spec / sum(spec$spec)
      p <- pmax(p, 1e-12)  # Avoid log(0)
      -sum(p * log(p))
    }
  )
}

# Fractal and Complexity Features
calculate_complexity_features <- function(samples) {
  list(
    Fractal_Dimension = as.numeric(fractaldim::fd.estimate(samples, methods = "madogram")$fd),
    Hjorth_Mobility = {
      activity <- var(samples)
      sqrt(var(diff(samples)) / activity)
    },
    Hjorth_Complexity = {
      activity <- var(samples)
      mobility <- sqrt(var(diff(samples)) / activity)
      (sqrt(var(diff(diff(samples))) / var(diff(samples)))) / mobility
    }
  )
}

# Tail Behavior Features
calculate_tail_features <- function(samples) {
  list(
    Gini_Coefficient = ineq(abs(samples - min(samples)), type = "Gini"),
    Energy = sum(samples^2)
  )
}

# ---------------------------
# Optimized Feature Extraction
# ---------------------------
calculate_features <- function(samples) {
  # Calculate all feature groups
  core_stats <- calculate_core_stats(samples)
  normality_tests <- calculate_normality_tests(samples)
  information_features <- calculate_information_features(samples)
  complexity_features <- calculate_complexity_features(samples)
  tail_features <- calculate_tail_features(samples)
  
  # Combine into single dataframe
  features <- data.frame(
    # Core statistics
    Mean = core_stats$Mean,
    Variance = core_stats$Variance,
    Skewness = core_stats$Skewness,
    Kurtosis = core_stats$Kurtosis,
    IQR = core_stats$IQR,
    MAD = core_stats$MAD,
    
    # Normality tests
    Jarque_Bera = normality_tests$Jarque_Bera,
    Anderson_Darling = normality_tests$Anderson_Darling,
    Shapiro_Wilk = normality_tests$Shapiro_Wilk,
    Shapiro_Francia = normality_tests$Shapiro_Francia,
    Lilliefors = normality_tests$Lilliefors,
    Cramer_Von_Misse = normality_tests$Cramer_Von_Misse,
    
    # Information and complexity
    #Entropy = information_features$Entropy,
    #Spectral_Entropy = information_features$Spectral_Entropy,
    #Fractal_Dimension = complexity_features$Fractal_Dimension,
    #Hjorth_Mobility = complexity_features$Hjorth_Mobility,
    #Hjorth_Complexity = complexity_features$Hjorth_Complexity,
    
    # Tail behavior
    Gini_Coefficient = tail_features$Gini_Coefficient,
    Energy = tail_features$Energy
  )
  
  return(features)
}


# ---------------------------------
# preprocessing function
# ---------------------------------
preprocess_data <- function(train_data) {
  numeric_train <- train_data[, sapply(train_data, is.numeric)]
  
  # Step 1: Standardization (Z-score scaling: mean = 0, std = 1)
  preProcStandard <- preProcess(numeric_train, method = c("medianImpute","center", "scale"))
  train_std <- predict(preProcStandard, numeric_train)
  
  # Step 2: Normalization (Rescale to [0,1])
  preProcNorm <- preProcess(train_std, method = "range")
  train_norm <- predict(preProcNorm, train_std)
  
  train_norm$Label <- train_data$Label
  
  return(list(train = train_norm, preProcStandard = preProcStandard, preProcNorm = preProcNorm))
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

set.seed(12345)
sample_size <- 30 
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
normal_data <- rbind(normal_data1, normal_data2, normal_data3, normal_data5, normal_data6)

data_all <- rbind(normal_data, non_normal_data)
data_all$Label <- as.factor(data_all$Label)

# -----------------------------
# Standardize & Normalize Data
# -----------------------------
norm_result <- preprocess_data(data_all)
train_norm <- norm_result$train

# define appropriate class reference
train_norm$Label <- relevel(train_norm$Label, ref = "Non_Normal")

# ---------------------------------------------------
# Define Common Training Control for Cross-Validation
set.seed(12345)
# ---------------------------------------------------
ctrl <- trainControl(method = "cv", number = 10, 
                     summaryFunction = defaultSummary,
                     classProbs = TRUE, 
                     search = "grid",
                     savePredictions = "final")

# Random Forest
rf_model <- train(Label ~ ., 
                  data = train_norm, 
                  method = "rf", 
                  trControl = ctrl, 
                  metric = "Accuracy")

# Artificial Neural Network
ann_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "nnet", 
                   trControl = ctrl, 
                   metric = "Accuracy", 
                   trace = FALSE)

# Gradient Boosting
gbm_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "gbm", 
                   trControl = ctrl, 
                   metric = "Accuracy", 
                   verbose = FALSE)

# Support Vector Machines
svm_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "svmRadial", 
                   trControl = ctrl, 
                   metric = "Accuracy")

# Collect all models into a list for later use
models_list <- list(
  "RF" = rf_model,
  "ANN" = ann_model,
  "GBM" = gbm_model,
  "SVM" = svm_model
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
  
  # assign appropriate names to predictions
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
distribution_set <- c("normal_5", "normal_25" ,"Gumbel", "triangular")

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



### Variable Importance Chart
# Extract variable importance
rf_var_imp <- varImp(models_list$`RF`, importance = TRUE)

# Convert to data frame
rf_var_imp_df <- data.frame(Variable = rownames(rf_var_imp$importance), 
                            Importance = rf_var_imp$importance$Overall)

# Sort by importance
rf_var_imp_df <- rf_var_imp_df[order(rf_var_imp_df$Importance), ]

# Save current par settings to restore later
old_par <- par(no.readonly = TRUE)
# Adjust margins
par(mar = c(5, 10, 4, 2))  

# Create the barplot
barplot(rf_var_imp_df$Importance,
        names.arg = rf_var_imp_df$Variable,
        horiz = TRUE,
        col = "steelblue",
        las = 1,  # Horizontal axis labels
        main = "Variable Importance - Random Forest",
        xlab = "Importance")

# Restore original par settings
par(old_par)


plot_combined_roc_base <- function(eval_results) {
  # Prepare plotting symbols and colors
  model_names <- names(eval_results)
  n_models <- length(model_names)
  colors <- rainbow(n_models)
  
  # Plot frame
  plot(0, 0, type = "n",
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate (1 - Specificity)",
       ylab = "True Positive Rate (Sensitivity)",
       main = "Combined ROC Curves (Positive Class: Non_Normal)")
  
  abline(a = 0, b = 1, col = "gray", lty = 2)  # Diagonal reference line
  
  legend_labels <- c()
  legend_colors <- c()
  
  for (i in seq_along(model_names)) {
    model_name <- model_names[i]
    pred_df <- eval_results[[model_name]]
    
    # Skip if prediction probabilities are missing
    if (!"Prob_Non_Normal" %in% colnames(pred_df$Predictions)) next
    
    actual_bin <- ifelse(pred_df$Predictions$True_Class == "Non_Normal", 1, 0)
    probs <- pred_df$Predictions$Prob_Non_Normal
    
    roc_obj <- pROC::roc(actual_bin, probs, levels = c(0, 1), direction = "<", quiet = TRUE)
    auc_val <- round(pROC::auc(roc_obj), 3)
    
    lines(1 - roc_obj$specificities, roc_obj$sensitivities,
          col = colors[i],
          lwd = 2)  # line thickness
    
    legend_labels <- c(legend_labels, paste0(model_name, " (AUC=", auc_val, ")"))
    legend_colors <- c(legend_colors, colors[i])
  }
  
  legend("bottomright",
         legend = legend_labels,
         col = legend_colors,
         lty = 1, lwd = 2,
         bty = "o",  # 'o' adds a box
         cex = 0.85,
         title = "Models")
  
  grid(col = "lightgray", lty = "dotted")
  
}

plot_combined_roc_base(eval_results)

