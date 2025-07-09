if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  e1071, tseries, nortest, gbm, lawstat, infotheo, ineq,
  caret, pROC, ROCR, randomForest, evd, discretization,
  nnet, ggplot2, pracma, fractaldim
)

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
  zero_crossings / (length(samples) - 1)
}

calculate_gini_coefficient <- function(samples) {
  samples_abs <- abs(samples - min(samples))
  ineq(samples_abs, type = "Gini")
}

calculate_outliers <- function(samples) {
  qnt <- quantile(samples, probs = c(0.25, 0.75))
  H <- 1.5 * IQR(samples)
  sum(samples < (qnt[1] - H) | samples > (qnt[2] + H))
}

# Statistical tests for normality
calculate_shapiro_wilk    <- function(samples) as.numeric(shapiro.test(samples)$statistic)
calculate_shapiro_francia <- function(samples) as.numeric(nortest::sf.test(samples)$statistic)
calculate_lilliefors      <- function(samples) as.numeric(nortest::lillie.test(samples)$statistic)
calculate_cramer_von_mises<- function(samples) as.numeric(nortest::cvm.test(samples)$statistic)

# Central tendency, spread & shape
calculate_mean           <- function(samples) mean(samples)
calculate_median         <- function(samples) median(samples)
calculate_variance       <- function(samples) var(samples)
calculate_iqr            <- function(samples) IQR(samples)
calculate_mad            <- function(samples) mad(samples)
calculate_rms            <- function(samples) sqrt(mean(samples^2))

# Dispersion & variation
calculate_cv             <- function(samples) sd(samples) / mean(samples)
calculate_range          <- function(samples) max(samples) - min(samples)
calculate_entropy        <- function(samples) {
  infotheo::entropy(discretize(samples, nbins = 10))
}

# Time-domain signal characteristics
calculate_peak_to_trough <- function(samples) max(samples) / abs(min(samples))
calculate_sample_size    <- function(samples) length(samples)

# Autocorrelation & whiteness
calculate_acf1           <- function(samples) {
  stats::acf(samples, plot=FALSE, lag.max=1)$acf[2]
}
calculate_box_test  <- function(samples) {
  as.numeric(Box.test(samples, lag = 1, type = "Ljung-Box")$statistic)
}

# Frequency-domain features
calculate_sample_entropy <- function(samples) {
  return(infotheo::entropy(discretize(samples, nbins = 10)))
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

# ---------------------------
# Master function: extract all features
# ---------------------------
calculate_features <- function(samples) {
  # Basic moments & tests
  skewness        <- e1071::skewness(samples)
  kurtosis        <- e1071::kurtosis(samples)
  jb_stat         <- as.numeric(tseries::jarque.bera.test(samples)$statistic)
  ad_stat         <- as.numeric(nortest::ad.test(samples)$statistic)
  # Normality & shape
  sw_stat         <- calculate_shapiro_wilk(samples)
  sf_stat         <- calculate_shapiro_francia(samples)
  lf_stat         <- calculate_lilliefors(samples)
  cvm_stat        <- calculate_cramer_von_mises(samples)
  # Time & distributional
  zcr             <- calculate_zero_crossing_rate(samples)
  gini            <- calculate_gini_coefficient(samples)
  outliers        <- calculate_outliers(samples)
  # Central & spread
  mean_val        <- calculate_mean(samples)
  median_val      <- calculate_median(samples)
  var_val         <- calculate_variance(samples)
  iqr_val         <- calculate_iqr(samples)
  mad_val         <- calculate_mad(samples)
  range_val       <- calculate_range(samples)
  cv_val          <- calculate_cv(samples)
  rms_val         <- calculate_rms(samples)
  # Entropy & size
  entropy_val     <- calculate_entropy(samples)
  samp_size       <- calculate_sample_size(samples)
  samp_entropy    <- calculate_sample_entropy(samples)
  # Time-domain signal features
  pt_ratio        <- calculate_peak_to_trough(samples)
  # Frequency-domain
  spec_entropy    <- calculate_spectral_entropy(samples)
  spec_centroid   <- calculate_spectral_centroid(samples)
  # Autocorrelation & whiteness
  acf1_val        <- calculate_acf1(samples)
  box_val         <- calculate_box_test(samples)
  # Fractal & complexity
  fd_val          <- calculate_fractal_dimension(samples)
  hjorth_vals     <- calculate_hjorth(samples)
  energy <- calculate_energy(samples)
  
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
    Outliers                = outliers,
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
    Spectral_Entropy         = spec_entropy,
    Spectral_Centroid        = spec_centroid,
    ACF1                     = acf1_val,
    Box_Ljung_Stat           = box_val,
    Fractal_Dimension        = fd_val,
    Hjorth_Activity          = hjorth_vals["Activity"],
    Hjorth_Mobility          = hjorth_vals["Mobility"],
    Hjorth_Complexity        = hjorth_vals["Complexity"]
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
  # add the labels back
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

# ---------------------------
# Generate the Data
# ---------------------------
set.seed(12345)
sample_size <- 10  
N <- 100            

normal_data1 <- generate_data(sample_size, 4*N, "normal", "Normal")
normal_data2 <- generate_data(sample_size, 3*N, "normal_50", "Normal")
normal_data3 <- generate_data(sample_size, 3*N, "normal_5", "Normal")
# non-normal data 
lognormal <- generate_data(sample_size, N, "LogNormal", "Non_Normal")
chisq_data   <- generate_data(sample_size, N, "Chi_Square", "Non_Normal")
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

# Define Common Training Control for Cross-Validation
set.seed(12345)
# ---------------------------
ctrl <- trainControl(method = "cv", number = 10,
                     search = "grid",
                     summaryFunction = twoClassSummary,
                     classProbs = TRUE,
                     savePredictions = "final")

# ctrl <- trainControl(
#   method = "cv",
#   number = 10,
#   summaryFunction = defaultSummary,  # <-- this returns Accuracy & Kappa
#   classProbs      = FALSE            # accuracy doesn’t need probabilities
# )


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


# # Extract and display the tuning parameters for each model
# list(
#   Logistic_Regression = log_model$bestTune,
#   Random_Forest = rf_model$bestTune,
#   Neural_Network = ann_model$bestTune,
#   Gradient_Boosting = gbm_model$bestTune,
#   Support_Vector_Machine = svm_model$bestTune,
#   K_Nearest_Neighbors = knn_model$bestTune
# )


# Collect all models into a list for later use
models_list <- list(
  "Logistic Regression" = log_model,
  "Random Forest" = rf_model,
  "ANN" = ann_model,
  "GBM" = gbm_model,
  "SVM" = svm_model,
  "KNN" = knn_model
)

save.image(paste0("machine_learning_models",".RData"))

# Harder Testing

### Predictions on unseen data row - by - row

simulate_predictions <- function(distributions ,n_iter = 1000, n = 10,
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

# example
distribution_set <- c("normal_15", "normal_100", "beta", "extremeskew")

eval_results <- simulate_predictions(
  distributions = distribution_set,
  n_iter = 1000,
  n = 10,
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

# -------------- Plots ---------
# 1. Extract and prepare data
rf_var_imp <- varImp(models_list$`Random Forest`)
rf_var_imp_df <- data.frame(
  Variable = rownames(rf_var_imp$importance),
  Importance = rf_var_imp$importance$Overall
)

# 2. Order by importance (ascending for proper bar order)
rf_var_imp_df <- rf_var_imp_df[order(rf_var_imp_df$Importance), ]

# 3. Calculate dynamic margin (prevents cutoff)
max_name_length <- max(nchar(as.character(rf_var_imp_df$Variable)))
left_margin <- min(20, max(8, max_name_length * 0.5))  # Auto-adjusting margin
par(mar = c(5, left_margin, 4, 4) + 0.1)  # bottom, left, top, right

# 4. Create horizontal barplot
barplot(
  height = rf_var_imp_df$Importance,
  horiz = TRUE,
  names.arg = rf_var_imp_df$Variable,
  las = 1,  # Horizontal labels
  cex.names = 0.85, 
  col = "steelblue",
  border = NA,
  space = 0.2,
  main = "Variable Importance - Random Forest",
  xlab = "Importance",
  xlim = c(0, max(rf_var_imp_df$Importance) * 1.1),
  panel.first = grid(ny = NA, col = "gray90", lty = "dotted")  # Vertical grid
)

# 5. Add importance values on right side
text(
  x = rf_var_imp_df$Importance * 1.05,  # Position right of bars
  y = seq_along(rf_var_imp_df$Importance) * 1.2 - 0.5,
  labels = round(rf_var_imp_df$Importance, 1),
  cex = 0.8,
  xpd = TRUE  # Allow writing in margin
)

### 2. Combined ROC Curves (base R)

plot_combined_roc_base <- function(eval_results) {
  cols <- rainbow(length(eval_results))
  i <- 1
  aucs <- character()
  
  for (model_name in names(eval_results)) {
    pred_df <- eval_results[[model_name]]
    if (!"Prob_Non_Normal" %in% colnames(pred_df$Predictions)) next
    
    # prepare ROC
    actual_bin <- ifelse(pred_df$Predictions$True_Class == "Non_Normal", 1, 0)
    probs      <- pred_df$Predictions$Prob_Non_Normal
    roc_obj    <- pROC::roc(actual_bin, probs, levels = c(0,1), direction = "<", quiet = TRUE)
    auc_val    <- round(pROC::auc(roc_obj), 3)
    aucs[i]    <- paste0(model_name, " (AUC=", auc_val, ")")
    
    # base R plot/lines
    if (i == 1) {
      plot(
        roc_obj$specificities,
        roc_obj$sensitivities,
        type = "l",
        col  = cols[i],
        lwd  = 2,
        xlim = c(1,0),           # flip x-axis to go from 0→1
        ylim = c(0,1),
        xlab = "False Positive Rate (1 - Specificity)",
        ylab = "True Positive Rate (Sensitivity)",
        main = "Combined ROC Curves\n(Positive Class: Non_Normal)"
      )
      # instead of abline(0, 1)
      abline(a = 1, b = -1, lty = "dashed", col = "gray")
      
    } else {
      lines(
        roc_obj$specificities,
        roc_obj$sensitivities,
        col = cols[i],
        lwd = 2
      )
    }
    i <- i + 1
  }
  
  legend(
    "bottomright",
    legend = aucs,
    col    = cols[1:(i-1)],
    lwd    = 2,
    bty    = "n"
  )
}

# Call it:
plot_combined_roc_base(eval_results)

