setwd("~/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach_to_NT/version 1")
# -----------------------------------
# 1.        INSTALL PACKAGES          #
# -----------------------------------
if (!require("pacman")) install.packages("pacman")
pacman::p_load(e1071, tseries, nortest, lawstat, moments, ineq, caret,
               randomForest, gbm, pROC, ROCR, nnet, ggplot2, mlbench, dplyr, kernlab,fractaldim
)

# -----------------------------------
# 2.  DATA GENERATION FUNCTION      #
# -----------------------------------
gen_data <- function(n, dist, par = NULL) {
  # Input validation
  if (!is.numeric(n) || n <= 0) stop("n must be a positive integer")
  dist <- tolower(dist)
  
  if (dist == "normal") {
    if (is.null(par)) par <- c(0, 1)
    samples <- rnorm(n, mean = par[1], sd = par[2])
    
  } else if (dist == "chi_square") {
    if (is.null(par)) par <- 3
    samples <- rchisq(n, df = par)
    
  }else if (dist == "gamma") {
    if (is.null(par)) par <- c(3, 0.1)
    samples <- rgamma(n, shape = par[1], rate = par[2])
    
  } else if (dist == "exponential") {
    if (is.null(par)) par <- 1
    samples <- rexp(n, rate = par)
    
  } else if (dist == "f") {
    if (is.null(par)) par <- c(6, 15)
    samples <- rf(n, df1 = par[1], df2 = par[2])
    
  } else if (dist == "t") {
    if (is.null(par)) par <- 3
    samples <- rt(n, df = par) 
    
  }else if (dist == "uniform") {
    if (is.null(par)) par <- c(0, 1)
    samples <- runif(n, min = par[1], max = par[2])
    
  } else if (dist == "laplace") {
    if (is.null(par)) par <- c(2, 7)  
    samples <- LaplacesDemon::rlaplace(n, location = par[1], scale = par[2])
    
  } else if (dist == "cauchy") {
    if (is.null(par)) par <- c(0, 1)
    samples <- rcauchy(n, location = par[1], scale = par[2])
    
  } else if (dist == "gumbel") {
    if (is.null(par)) par <- c(0, 1)
    samples <- evd::rgumbel(n, loc = par[1], scale = par[2])
    
  } else if (dist == "weibull") {
    if (is.null(par)) par <- c(1, 2)
    samples <- rweibull(n, shape = par[1], scale = par[2]) 
    
  } else if (dist == "lognormal") {
    if (is.null(par)) par <- c(0, 1)
    samples <- rlnorm(n, meanlog = par[1], sdlog = par[2])
    
  } else if (dist == "contaminated") {
    if (is.null(par)) par <- c(0.75, 0, 1, 5)
    br <- rbinom(n, size = 1, prob = par[1])
    sd_br <- ifelse(br == 1, par[3], par[4])
    samples <- rnorm(n, mean = par[2], sd = sd_br)
    
  } else if (dist == "pareto") {
    if (is.null(par)) par <- c(1, 3)
    samples <- VGAM::rpareto(n, scale = par[1], shape = par[2])
    
  }else if (dist == "beta") {
    if (is.null(par)) par <- c(2, 5)
    if (any(par <= 0)) stop("Beta parameters must be > 0")
    samples <- rbeta(n, shape1 = par[1], shape2 = par[2])
    
  }else if (dist == "logistic") {
    if (is.null(par)) par <- c(0, 1)
    if (par[2] <= 0) stop("Logistic scale must be > 0")
    samples <- rlogis(n, location = par[1], scale = par[2])
    
  } else {
    stop("Unsupported distribution: ", dist)
  }
  
  return(samples)
}

# ========================================================================
# -----------------------------------
# 3. FEATURE EXTRACTION FUNCTIONS
# -----------------------------------
# Zero Crossing Rate
calculate_zero_crossing_rate <- function(samples) {
  signs <- samples > 0
  zero_crossings <- sum(abs(diff(signs)))
  return(zero_crossings / (length(samples) - 1))
}

# Outlier Detection
calculate_outliers <- function(samples) {
  qnt <- quantile(samples, probs = c(0.25, 0.75))
  H <- 1.5 * IQR(samples)
  return(sum(samples < (qnt[1] - H) | samples > (qnt[2] + H)))
}

# Peak to Trough Ratio
calculate_peak_to_trough <- function(samples) {
  max_val <- max(samples)
  min_val <- min(samples)
  if(min_val == 0) return(NA)
  return(max_val / abs(min_val))
}

# Spectral Entropy
calculate_spectral_entropy <- function(samples) {
  spec <- spec.pgram(samples, plot = FALSE, na.action = na.omit)
  p <- spec$spec / sum(spec$spec)
  p <- p[p > 0]
  return(-sum(p * log(p)))
}

# Lin & Mudholkar Zp statistic
calculate_zp_statistic <- function(samples) {
  x <- sort(na.omit(samples))
  n <- length(x)
  n1 <- n - 1
  return(atanh(cor(x, ((sum(x^2) - x^2)/n1 - ((sum(x) - x)/n1)^2)^(1/3))))
}

# Vasicek K_{m,n} statistic
calculate_vasicek_kmn <- function(samples) {
  m = floor(sqrt(sum(!is.na(samples))))
  x <- sort(na.omit(samples))
  n <- length(x) 
  s <- sqrt(mean((x - mean(x))^2))
  return((n/(2*m*s)) * exp(mean(log(x[pmin(n, 1:n + m)] - x[pmax(1L, 1:n - m)]))))
}

# Safe Feature Calculation Wrapper
safe_calculate <- function(expr, default = NA) {
  tryCatch(expr, error = function(e) default)
}

# ----------------------------
#  calculate features
# ----------------------------
calculate_features <- function(samples) {
  # Remove NA values 
  samples <- na.omit(samples)
  
  # Basic statistics
  mean_val <- mean(samples)
  median_val <- median(samples)
  var_val <- var(samples)
  iqr_val <- IQR(samples)
  mad_val <- mad(samples)
  range_val <- max(samples) - min(samples)
  cv_val <- if(mean_val != 0) sd(samples)/mean_val else NA
  rms_val <- sqrt(mean(samples^2))
  
  # Shape statistics
  skewness <- safe_calculate(as.numeric(e1071::skewness(samples)))
  kurtosis <- safe_calculate(as.numeric(e1071::kurtosis(samples)))
  
  # Normality tests
  jb_stat <- safe_calculate(as.numeric(tseries::jarque.bera.test(samples)$statistic))
  ad_stat <- safe_calculate(as.numeric(nortest::ad.test(samples)$statistic))
  sw_stat <- safe_calculate(as.numeric(shapiro.test(samples)$statistic))
  sf_stat <- safe_calculate(as.numeric(nortest::sf.test(samples)$statistic))
  lf_stat <- safe_calculate(as.numeric(nortest::lillie.test(samples)$statistic))
  
  # Additional features
  zcr <- safe_calculate(as.numeric(calculate_zero_crossing_rate(samples)))
  gini <- safe_calculate(as.numeric(ineq::ineq(abs(samples - min(samples)), type = "Gini")))
  pt_ratio <- safe_calculate(as.numeric(calculate_peak_to_trough(samples)))
  spec_entropy <- safe_calculate(as.numeric(calculate_spectral_entropy(samples)))
  energy <- sum(samples^2)
  
  # Zp and Vasicek’s statistic
  zp_stat <- safe_calculate(as.numeric(calculate_zp_statistic(samples)))
  vasicek_stat <- safe_calculate(as.numeric(calculate_vasicek_kmn(samples)))
  
  # Statistical features
  Tail_Weight_Ratio <- safe_calculate(as.numeric(quantile(samples, 0.95) / quantile(samples, 0.05)))
  Moment_Ratio <- safe_calculate(as.numeric(moments::moment(samples, 4) / (moments::moment(samples, 2)^2)))
  Q_Q_Correlation = cor(quantile(samples, ppoints(100)), qnorm(ppoints(100)))
  
  # Create feature data frame
  features <- data.frame(
    Median = median_val, 
    Variance = var_val,
    IQR = iqr_val, 
    MAD = mad_val, 
    Range = range_val, 
    CV = cv_val,
    Root_Mean_Square = rms_val, 
    Skewness = skewness, 
    Kurtosis = kurtosis,
    Jarque_Bera = jb_stat, 
    Anderson_Darling = ad_stat,
    Shapiro_Wilk = sw_stat, 
    Shapiro_Francia = sf_stat,
    Lilliefors = lf_stat, 
    Zero_Cross_Rate = zcr, 
    Gini_Coefficient = gini,
    Peak_to_Trough = pt_ratio,
    Spectral_Entropy = spec_entropy,
    Q_Q_Correlation  = Q_Q_Correlation,
    Energy = energy, 
    Tail_Weight_Ratio = Tail_Weight_Ratio,
    Moment_Ratio = Moment_Ratio,
    Zp = zp_stat,
    Vasicek_Kmn = vasicek_stat,
    stringsAsFactors = FALSE
  )
  
  return(features)
}


# ------------------------------------------
# Standardization & Normalization Function
# ------------------------------------------
preprocess_data <- function(train_data) {
  # choose only numeric columns(exclude label column)
  numeric_train <- train_data[, sapply(train_data, is.numeric)]
  
  # Step 1: Standardization (Z-score scaling: mean = 0, std = 1)
  preProcStandard <- preProcess(numeric_train, method = c("center", "scale"))
  train_std <- predict(preProcStandard, numeric_train)
  
  # Step 2: Normalization (Rescale to [0,1])
  preProcNorm <- preProcess(train_std, method = "range")
  train_norm <- predict(preProcNorm, train_std)
  
  # add label comuln back
  train_norm$Label <- train_data$Label
  
  return(list(train = train_norm, 
              preProcStandard = preProcStandard, 
              preProcNorm = preProcNorm))
}

# --------------------------------------------
# Function to Generation and Extract Features 
# --------------------------------------------
generate_data <- function(sample_size, N, dist = "normal", par = NULL, label) {
  # row-by-row generation
  data <- do.call(rbind, lapply(1:N, function(x) {
    samples <- gen_data(sample_size, dist, par = par) 
    features <- calculate_features(samples)
    features$Label <- label
    return(features)
  }))
  return(data)
}


# generate the data and extract the features
set.seed(12345)
sample_size <- 8 
num.sim <- 100            

# Normal data
normal_data1 <- generate_data(sample_size, 2*num.sim, "normal", par = c(0, 1), "Normal")
normal_data2 <- generate_data(sample_size, 2*num.sim, "normal", par = c(5, 2) , "Normal")
normal_data3 <- generate_data(sample_size, 2*num.sim, "normal" , par = c(15, 8) ,"Normal")
normal_data4 <- generate_data(sample_size, 2*num.sim, "normal", par = c(30, 5),  "Normal")
normal_data5 <- generate_data(sample_size, 2*num.sim, "normal", par = c(100, 25), "Normal")

# non-normal data 
lognormal <- generate_data(sample_size, num.sim, "lognormal", par = c(0, 1), "Non_Normal")
chisq_data   <- generate_data(sample_size, num.sim, "chi_square", par = 7,"Non_Normal")
exp_data     <- generate_data(sample_size, num.sim, "exponential", par = 3, "Non_Normal")
Weibull      <- generate_data(sample_size, num.sim, "Weibull", par = c(2, 1), "Non_Normal")
Pareto      <- generate_data(sample_size, num.sim, "Pareto", par = c(1, 2), "Non_Normal")
Gamma        <- generate_data(sample_size, num.sim, "Gamma", par = c(2, 1), "Non_Normal")
Uniform      <- generate_data(sample_size, num.sim, "Uniform", par = c(0, 1), "Non_Normal")
t       <- generate_data(sample_size, num.sim, "t", par = 3, "Non_Normal")
f       <- generate_data(sample_size, num.sim, "f", par = c(6, 15), "Non_Normal")
Cauchy       <- generate_data(sample_size, num.sim, "Cauchy", par = c(0, 1), "Non_Normal")

non_normal_data <- rbind(lognormal, chisq_data, exp_data, Weibull, Pareto, Gamma, Uniform, t, f, Cauchy)
normal_data <- rbind(normal_data1, normal_data2, normal_data3, normal_data4, normal_data5)

train_data <- rbind(normal_data, non_normal_data)
train_data$Label <- as.factor(train_data$Label)

# ---------------------------
# Standardize & Normalize Data
# ---------------------------
norm_result <- preprocess_data(train_data)
train_norm <- norm_result$train

# define appropriate class reference
train_norm$Label <- relevel(train_norm$Label, ref = "Non_Normal")


# ------------------------------------
# 4. MODEL TRAINING FUNCTIONS
# ------------------------------------
# Training control with cross-validation
ctrl <- caret::trainControl(
  method = "cv", 
  number = 10,
  summaryFunction = twoClassSummary,
  classProbs = TRUE, 
  search = "grid",
  savePredictions = "final"
)

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
                  tuneGrid = expand.grid(mtry = c(2, 5, 10, 15)),
                  ntree = 500, 
                  importance = TRUE,
                  metric = "ROC")

# Artificial Neural Network
ann_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "nnet",
                   trControl = ctrl, 
                   tuneGrid = expand.grid(
                     size = c(5, 10, 15),
                     decay = c(0.001, 0.01, 0.1)),
                   MaxNWts = 2000, 
                   maxit = 200,
                   metric = "ROC", 
                   trace = FALSE)

# Gradient Boosting
gbm_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "gbm", 
                   trControl = ctrl, 
                   metric = "ROC", 
                   tuneGrid = expand.grid(
                     n.trees = c(100, 300, 500),
                     interaction.depth = c(1, 3, 5),
                     shrinkage = c(0.01, 0.1),
                     n.minobsinnode = c(5, 10)
                   ),
                   verbose = FALSE)

# Support Vector Machines
svm_model <- train(Label ~ ., 
                   data = train_norm, 
                   method = "svmRadial", 
                   trControl = ctrl, 
                   tuneGrid = expand.grid(
                     C = c(0.1, 1, 10, 100),
                     sigma = c(0.01, 0.1, 0.5, 1)
                   ),
                   metric = "ROC")

# K-Nearest Neighbors
knn_model <- train(Label ~ .,
                   data = train_norm,
                   method = "knn", 
                   tuneGrid = expand.grid(k = c(3, 5, 7, 9, 11)),
                   trControl = ctrl,
                   metric = "ROC")

# Train the XGBoost model
capture.output({
  xgb_model <- train(
    Label ~ ., 
    data = train_norm, 
    method = "xgbTree", 
    trControl = ctrl, 
    metric = "ROC", 
    tuneGrid = expand.grid(
      nrounds = c(100, 300, 500),
      max_depth = c(3, 6, 9),
      eta = c(0.01, 0.1),
      gamma = 0,
      colsample_bytree = 1,
      min_child_weight = c(1, 5),
      subsample = 1
    ),
    verbose = FALSE
  )
})


# Collect all models into a list for later use
models_list <- list(
  "LR" = log_model,
  "RF" = rf_model,
  "ANN" = ann_model,
  "GBM" = gbm_model,
  "SVM" = svm_model,
  "KNN" = knn_model,
  "XGB" = xgb_model
) 

# -------------------------------------------------------
#              5. Variable Importance Chart
# -------------------------------------------------------
varImplot <- function(rf_model, top_n = NULL){
  
  rf_var_imp <- varImp(rf_model)
  imp <- rf_var_imp$importance
  
  # Pick an importance vector
  if ("Overall" %in% colnames(imp)) {
    imp_vec <- imp$Overall
  } else if (ncol(imp) == 1) {
    imp_vec <- imp[, 1]
  } else {
    # two-class case: average the class-specific importances
    imp_vec <- rowMeans(imp)
  }
  
  # create a vip dataframe 
  rf_var_imp_df <- data.frame(
    Variable = rownames(imp),
    Importance = as.numeric(imp_vec),
    stringsAsFactors = FALSE
  )
  
  # Order features and Optionally keep only top_n
  rf_var_imp_df <- rf_var_imp_df[order(rf_var_imp_df$Importance, decreasing = TRUE), ]
  if (!is.null(top_n)) rf_var_imp_df <- head(rf_var_imp_df, top_n)
  
  # plot(reverse order so biggest on top in horizontal barplot)
  rf_var_imp_df <- rf_var_imp_df[order(rf_var_imp_df$Importance), ]
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  par(mar = c(5, 10, 4, 2))
  
  # create the bar plot
  barplot(rf_var_imp_df$Importance,
          names.arg = rf_var_imp_df$Variable,
          horiz = TRUE,
          col = "steelblue",
          las = 1,
          main = "Variable Importance - Random Forest",
          xlab = "Importance")
}

# generate plots
pdf("results/vip.pdf", width=8, height=8)
varImplot(models_list$RF)
dev.off()

# ---------------------------------------------------------- 
#    6. Predictions on unseen data row - by - row
# ----------------------------------------------------------
simulate_predictions <- function(distributions ,
                                 n_iter = 1000, 
                                 n = sample_size,
                                 trained_models, 
                                 preProcStandard, 
                                 preProcNorm) {
  
  # Initialize storage for all model predictions
  all_predictions <- lapply(names(trained_models), function(x) {
    data.frame(True_Class = character(),
               Predicted_Class = character(),
               Prob_Non_Normal = numeric(),
               stringsAsFactors = FALSE)})
  
  # assign appropriate names to predictions
  names(all_predictions) <- names(trained_models)
  # get true class
  get_true_class <- function(dist_name) {
    if (startsWith(dist_name, "normal")) "Normal" else "Non_Normal"
  }
  
  # predict for each sample data individually
  for (dist_obj in distributions) {
    
    dist_name <- dist_obj$dist
    dist_par  <- dist_obj$par
    
    for (i in 1:n_iter) {
      
      # generate data and identify the class
      samples <- gen_data(n, dist_name, par = dist_par)
      true_class <- get_true_class(dist_name)
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
    
    # create the confusion matrix
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

# ----------------------------------------------
#          Make the prediction
# ----------------------------------------------

# create test set distributions
set.seed(12345)
distribution_set <- list(
  list(dist = "normal", par = c(0, 5)),          
  list(dist = "normal", par = c(50, 15)),        
  list(dist = "normal", par = c(25, 10)), 
  list(dist = "normal", par = c(50, 15)),          
  list(dist = "normal", par = c(50, 15)),        
  list(dist = "normal", par = c(25, 10)), 
  list(dist = "gumbel", par = c(0, 1)),         
  list(dist = "beta", par = c(2, 5)),
  list(dist = "laplace", par = c(0, 4)),
  list(dist = "cauchy", par = c(0, 1)),
  list(dist = "uniform", par = c(0, 1)),
  list(dist = "contaminated", par = c(0.75, 0, 1, 5)),
  list(dist = "chi_square", par = 3)
)


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

# ----------------------------------------------------
#        7.  compare models through ROC curves
# ----------------------------------------------------
plot_combined_roc_base <- function(eval_results) {
  # Prepare plotting symbols and colors
  model_names <- names(eval_results)
  n_models <- length(model_names)
  colors <- c(
    "black","red","blue","forestgreen","orange",
    "purple","brown","magenta","pink","darkcyan"
  )
  
  # Plot frame
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate (1 - Specificity)",
       ylab = "True Positive Rate (Sensitivity)",
       main = "Model AUROC Comparison (Positive Class: Non_Normal)", cex.main = 0.75)
  
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
pdf("results/roc_curve.pdf", width=5, height=5)
plot_combined_roc_base(eval_results)
dev.off()

# -----------------------------------------------
#    8. classify a single data set
# -----------------------------------------------
classify_sample <- function(sample_data, trained_models, preProcStandard, preProcNorm,
                            indiv_threshold = 0.50, final_threshold = 0.60) {
  
  sample_features <- calculate_features(sample_data)
  sample_std  <- predict(preProcStandard, sample_features)
  sample_norm <- predict(preProcNorm, sample_std)
  
  results <- lapply(names(trained_models), function(model_name) {
    model <- trained_models[[model_name]]
    
    p_non_normal <- predict(model, newdata = sample_norm, type = "prob")[, "Non_Normal"]
    pred_class <- ifelse(p_non_normal >= indiv_threshold, "Non_Normal", "Normal")
    
    data.frame(
      Model = model_name,
      Non_Normal_Probability = round(p_non_normal, 4),
      Normal_Probability = round(1 - p_non_normal, 4),
      Prediction = pred_class,
      stringsAsFactors = FALSE
    )
  })
  
  df <- do.call(rbind, results)
  # ---- Majority rule (>= 75% Non_Normal) ----
  prop_non_normal <- mean(df$Prediction == "Non_Normal")
  final_prediction <- ifelse(prop_non_normal >= final_threshold, "Non_Normal", "Normal")
  
  list(
    per_model = df,
    final_prediction = final_prediction,
    prop_non_normal = prop_non_normal
  )
}


# save results
save.image(file = "results/trained_ml_model_n10.RData")


# ------------------------------------------------------
#         9. EXAMPLE USAGE
# ------------------------------------------------------
# Generate test data (replace with user's actual data)
test_sample <- gen_data(n = sample_size, dist = "t", par = 3)

classification_results <- classify_sample(
  sample_data = test_sample,
  trained_models = models_list,
  preProcStandard = norm_result$preProcStandard,
  preProcNorm = norm_result$preProcNorm,
  indiv_threshold = 0.50,   # per-model cutoff
  final_threshold = 0.60    # final vote cutoff
)

# Print formatted results
print(classification_results$per_model)
cat("\nProp Non_Normal votes:", round(classification_results$prop_non_normal, 3),
    "\nFinal prediction:", classification_results$final_prediction, "\n")


# Demo to chapter 1
N <- 1000
alpha <- 0.05
effect_size <- 0.75
threshold_vec <- c(0.0, 0.2, 0.3, 0.5, 0.75, 0.8, 1)
type1_error <- power <- numeric(length(threshold_vec))

# fixed vote proportion from your classification step
prop_non_normal <- classification_results$prop_non_normal

for (j in seq_along(threshold_vec)) {
  decision_threshold <- threshold_vec[j]
  
  # choose test based on final decision threshold 
  final_prediction <- ifelse(prop_non_normal >= decision_threshold, "Non_Normal", "Normal")
  
  pval_error <- pval_power <- numeric(N)
  
  for (i in 1:N) {
    x <- gen_data(n = 10, dist = "normal", par = c(0,1))
    
    if (final_prediction == "Normal") {
      pval_error[i] <- t.test(x)$p.value
      pval_power[i] <- t.test(x + effect_size)$p.value
    } else {
      pval_error[i] <- wilcox.test(x )$p.value
      pval_power[i] <- wilcox.test(x + effect_size)$p.value
    }
  }
  
  type1_error[j] <- mean(pval_error < alpha)
  power[j] <- mean(pval_power < alpha)
}

type1_error
power


