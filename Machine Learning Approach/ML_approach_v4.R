setwd("~/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach")
source("~/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach/gen_data_fun.R")

if(!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest, gbm, lawstat, infotheo, ineq, caret, pROC, ROCR, randomForest,
               evd, discretization, nnet, ggplot2, mlbench, infotheo, dplyr)

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

calculate_peak_to_trough <- function(samples) max(samples) / abs(min(samples))

calculate_box_test  <- function(samples) {
  as.numeric(Box.test(samples, lag = 1, type = "Ljung-Box")$statistic)
}


calculate_spectral_entropy <- function(samples) {
  spec <- stats::spec.pgram(samples, plot = FALSE)
  p <- spec$spec / sum(spec$spec)
  -sum(p * log(p))
}

calculate_spectral_centroid <- function(samples) {
  spec <- stats::spec.pgram(samples, plot = FALSE)
  p <- spec$spec / sum(spec$spec)
  freq <- spec$freq
  sum(freq * p) / sum(p)
}

# Fractal & complexity measures
calculate_hjorth <- function(samples) {
  activity   <- var(samples)
  mobility   <- sqrt(var(diff(samples)) / activity)
  complexity <- (sqrt(var(diff(diff(samples))) / var(diff(samples)))) / mobility
  c(Activity = activity, Mobility = mobility, Complexity = complexity)
}


calculate_fractal_dimension <- function(samples) {
  res <- fractaldim::fd.estimate(samples, methods = "madogram")
  as.numeric(res$fd)
}

calculate_energy <- function(samples) {
  return(sum(samples^2))
}

# Dispersion & variation
calculate_cv             <- function(samples) sd(samples) / mean(samples)
calculate_range          <- function(samples) max(samples) - min(samples)
calculate_entropy        <- function(samples) {
  infotheo::entropy(discretize(samples, nbins = 10))
}

# ---------------------------
# extract all features
# ---------------------------
calculate_features <- function(samples) {
# Central & spread
mean_val        <- mean(samples)
median_val      <- median(samples)
var_val         <- var(samples)
iqr_val         <- IQR(samples)
mad_val         <- mad(samples)
range_val       <- max(samples) - min(samples)
cv_val          <- calculate_cv(samples)
rms_val         <- sqrt(mean(samples^2))

# Normality & shape
skewness        <- e1071::skewness(samples)
kurtosis        <- e1071::kurtosis(samples)
jb_stat         <- as.numeric(tseries::jarque.bera.test(samples)$statistic)
ad_stat         <- as.numeric(nortest::ad.test(samples)$statistic)
sw_stat         <- shapiro.test(samples)$statistic
sf_stat         <- nortest::sf.test(samples)$statistic
lf_stat         <- nortest::lillie.test(samples)$statistic
cvm_stat        <- nortest::cvm.test(samples)$statistic

# Time & distributional
zcr             <- calculate_zero_crossing_rate(samples)
gini            <- calculate_gini_coefficient(samples)
outliers        <- calculate_outliers(samples)

# Others
entropy_val     <- calculate_entropy(samples)
pt_ratio        <- calculate_peak_to_trough(samples)
box_val         <- calculate_box_test(samples)
spec_entropy    <- calculate_spectral_entropy(samples)
spec_centroid   <- calculate_spectral_centroid(samples)
fd_val          <- calculate_fractal_dimension(samples)
hjorth_vals     <- calculate_hjorth(samples)
energy          <- calculate_energy(samples)

# create a dataframe
features <- data.frame(
  Skewness                = skewness,
  Kurtosis                = kurtosis,
  Jarque_Bera             = jb_stat,
  Anderson_Darling        = ad_stat,
  Shapiro_Wilk            = sw_stat,
  Shapiro_Francia         = sf_stat,
  Lilliefors              = lf_stat,
  Cramer_Von_Misse        = cvm_stat,
  Zero_Cross_Rate         = zcr,
  Gini_Coefficient        = gini,
  #Outliers                = outliers,
  Mean                     = mean_val,
  Median                   = median_val,
  Variance                 = var_val,
  IQR                      = iqr_val,
  MAD                      = mad_val,
  Range                    = range_val,
  CV                       = cv_val,
  Root_Mean_Square         = rms_val,
  energy                   = energy,
  Peak_to_Trough           = pt_ratio,
  Enropy                   = entropy_val,
  Spectral_Entropy         = spec_entropy,
  Spectral_Centroid        = spec_centroid,
  Box_Ljung_Stat           = box_val,
  Fractal_Dimension        = fd_val,
  Hjorth_Activity          = hjorth_vals["Activity"],
  Hjorth_Mobility          = hjorth_vals["Mobility"],
  Hjorth_Complexity        = hjorth_vals["Complexity"]
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
sample_size <- 8 
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

# ---------------------------
# Standardize & Normalize Data
# ---------------------------
norm_result <- preprocess_data(data_all)
train_norm <- norm_result$train

# define appropriate class reference
train_norm$Label <- relevel(train_norm$Label, ref = "Non_Normal")


# ---------------------------
# Define Common Training Control for Cross-Validation
set.seed(12345)
# ---------------------------
ctrl <- trainControl(method = "cv", number = 10, 
                     summaryFunction = defaultSummary,
                     classProbs = TRUE, 
                     search = "grid",
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
                   metric = "Accuracy")

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

# K-Nearest Neighbors
knn_model <- train(Label ~ .,
                   data = train_norm,
                   method = "knn",
                   trControl = ctrl,
                   metric = "Accuracy")


# Collect all models into a list for later use
models_list <- list(
  "LR" = log_model,
  "RF" = rf_model,
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
distribution_set <- c("normal_5", "normal_100" ,"Gumbel", "beta")

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
varImplot <- function(rf_model, top_n = NULL){
  # Extract variable importance from RF model
  rf_var_imp <- varImp(rf_model)
  
  # Convert to data frame
  rf_var_imp_df <- data.frame(Variable = rownames(rf_var_imp$importance), 
                              Importance = rf_var_imp$importance$Overall)
  
  # Sort by importance
  rf_var_imp_df <- rf_var_imp_df[order(rf_var_imp_df$Importance), ]
  
  # set up margins & restore on exit
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  # Adjust margins
  par(mar = c(5, 10, 4, 2))  
  
  # Create the bar plot
  barplot(rf_var_imp_df$Importance,
          names.arg = rf_var_imp_df$Variable,
          horiz = TRUE,
          col = "steelblue",
          las = 1,  # Horizontal axis labels
          main = "Variable Importance - Random Forest",
          xlab = "Importance")
}


plot_combined_roc_base <- function(eval_results) {
  # Prepare plotting symbols and colors
  model_names <- names(eval_results)
  n_models <- length(model_names)
  colors <- rainbow(n_models)
  
  # Plot frame
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate (1 - Specificity)",
       ylab = "True Positive Rate (Sensitivity)",
       main = "Combined ROC Curves (Positive Class: Non_Normal)")
  # Diagonal reference line
  abline(a = 0, b = 1, col = "gray", lty = 2) 
  
  # add legend
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
    
    lines(1 - roc_obj$specificities, roc_obj$sensitivities,col = colors[i], lwd = 2)  
    
    legend_labels <- c(legend_labels, paste0(model_name, " (AUC=", auc_val, ")"))
    legend_colors <- c(legend_colors, colors[i])
  }
  
  legend("bottomright",
         legend = legend_labels,
         col = legend_colors,
         lty = 1, lwd = 2,
         bty = "o",  
         cex = 0.85,
         title = "Models")
  
  grid(col = "lightgray", lty = "dotted")
  
}

# generate plots
varImplot(models_list$RF)

plot_combined_roc_base(eval_results)

# save results
save(models_list, norm_result, eval_results, file = "trained_models.RData")

# --------------------------------
# classify a single data set
# --------------------------------
classify_sample <- function(sample_data, trained_models, preProcStandard, preProcNorm) {
  # 1. Calculate features from sample data
  sample_features <- calculate_features(sample_data)
  
  # 2. Preprocess using saved scalers
  sample_std <- predict(preProcStandard, sample_features)
  sample_norm <- predict(preProcNorm, sample_std)
  
  # 3. Make predictions with all models
  results <- lapply(names(trained_models), function(model_name) {
    model <- trained_models[[model_name]]
    pred_class <- as.character(predict(model, newdata = sample_norm))
    pred_prob <- predict(model, newdata = sample_norm, type = "prob")$Non_Normal
    
    data.frame(
      Model = model_name,
      Prediction = pred_class,
      Non_Normal_Probability = round(pred_prob, 4),
      stringsAsFactors = FALSE
    )
  })
  
  # 4. Combine and return results
  do.call(rbind, results)
}

# ========================
# EXAMPLE USAGE
# ========================

# Generate test data (replace with user's actual data)
test_sample <- generate_samples(n = sample_size, dist = "LogNormal")

# Classify the sample
classification_results <- classify_sample(
  sample_data = test_sample,
  trained_models = models_list,
  preProcStandard = norm_result$preProcStandard,
  preProcNorm = norm_result$preProcNorm
)

# Print formatted results
print(classification_results)
