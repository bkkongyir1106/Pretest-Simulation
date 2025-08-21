# Load helper source files
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/user_framework/ShinyApp")

# load user defined functions

FilePath <- "~/Desktop/OSU/Research/Pretest-Simulation/user_framework/ShinyApp/sample_test_functions/anova_functions/"
source(paste0(FilePath, "gen_data.R"))
source(paste0(FilePath, "get_parameters.R"))
source(paste0(FilePath, "fn_get_norm_obj.R"))
source(paste0(FilePath, "fn_for_ds_test_1.R"))
source(paste0(FilePath, "fn_for_ds_test_2.R"))

# -------------------------------------------------------------------
#-------------- Define all other functions here ---------------------
# -------------------------------------------------------------------
#----------   function for computing auc ----------------------------
compute_area <- function(sample_sizes, y) {
  if (is.null(y) || length(y) != length(sample_sizes)) return(NA_real_)
  sum(diff(sample_sizes) * (head(y, -1) + tail(y, -1)) / 2) /
    (max(sample_sizes) - min(sample_sizes))
}

# --------------------------------------------------------
# ---------- General Normality check function ------------
# --------------------------------------------------------
normality_test <- function(data, test = "SW", alpha = 0.05) {
  pvals <- NULL
  all_normality_satisfied <- NULL
  
  # Case 1: Single numeric vector
  if (is.numeric(data) && is.atomic(data) && is.null(dim(data))) {
    pval <- generate_tests(data, test = test)$p.value
    pvals <- pval
    all_normality_satisfied <- pval > alpha
  }
  
  # Case 2: List of numeric vectors
  else if (is.list(data) && !is.data.frame(data)) {
    pvals <- sapply(data, function(sample) {
      generate_tests(sample, test = test)$p.value
    })
    names(pvals) <- names(data) %||% paste0("Sample", seq_along(pvals))
    all_normality_satisfied <- all(pvals > alpha)
  }
  
  # Case 3: Wide-format data frame
  else if (is.data.frame(data) && all(sapply(data, is.numeric))) {
    pvals <- sapply(as.list(data), function(sample) generate_tests(sample, test = test)$p.value)
    names(pvals) <- names(data)
    all_normality_satisfied <- all(pvals > alpha)
  }
  
  # Case 4: Long-format with group labels
  else if ((is.data.frame(data) || is.matrix(data)) && ncol(data) >= 2) {
    # We assume first column is group, second is value
    grouped_samples <- split(data[[2]], data[[1]])
    pvals <- sapply(grouped_samples, function(sample) {
      generate_tests(sample, test = test)$p.value
    })
    all_normality_satisfied <- all(pvals > alpha)
  }
  else {
    stop("Unsupported input type: must be numeric vector, list of vectors, or grouped data.")
  }
  
  return(list(p_values = pvals, normality_satisfied = all_normality_satisfied))
}
# ------------------------------------------------------------------------------ 
# ------------------ Downstream Test function -------------------------------
# ------------------------------------------------------------------------------
ds_test_function <- function(
    gen_data = gen_data,
    get_parameters = get_parameters,
    fn_to_get_norm_obj = fn_to_get_norm_obj,
    fn_for_norm_test = normality_test,
    fn_for_ds_test_1 = fn_for_ds_test_1,
    fn_for_ds_test_2 = fn_for_ds_test_2,
    paras            = NULL,
    alpha = 0.05,
    norm_test_method = "SW",
    ...
) {
  # generate dataset
  data <- if (!is.null(paras)) do.call(gen_data, paras) else gen_data()
  
  # get normality test object
  normality_test_object     <- fn_to_get_norm_obj(data)
  normality_test_pval_ds_test  <- fn_for_norm_test(data = normality_test_object, test = norm_test_method, alpha = alpha)
  
  # choose test
  if (isTRUE(normality_test_pval_ds_test$normality_satisfied)) {
    ds_norm_test_method <- fn_for_ds_test_1(data)
  } else {
    ds_norm_test_method <- fn_for_ds_test_2(data, ...)
  }
  
  return(ds_norm_test_method$p.value)
}

# ------------------------------------------------------------------------------
perform_ds_func <- function(
    sample_sizes       = c(10, 20, 30, 40, 50),
    Nsim               = 1e3,
    alpha              = 0.05,
    n_boot             = NULL,
    gen_data           = gen_data,
    get_parameters     = function(n) list(n = n),
    fn_to_get_norm_obj = fn_to_get_norm_obj,
    fn_for_norm_test   = fn_for_norm_test,
    fn_for_ds_test_1   = fn_for_ds_test_1,
    fn_for_ds_test_2   = fn_for_ds_test_2,
    norm_test_method   = "SW",
    ds_test_methods    = c("parametric", "nonparametric", "adaptive"),
    ...
) {
  ds_test_methods <- match.arg(ds_test_methods, several.ok = TRUE)
  results <- list()
  
  for (method in ds_test_methods) {
    cat("Running method:", method, "\n")
    ds_test_results <- numeric(length(sample_sizes))
    names(ds_test_results) <- paste0("n=", sample_sizes)
    
    for (i in seq_along(sample_sizes)) {
      n <- sample_sizes[i]
      rejections <- 0
      
      for (sim in seq_len(Nsim)) {
        paras <- get_parameters(n)
        data <- if (!is.null(paras)) do.call(gen_data, paras) else gen_data()
        
        p_value <- tryCatch({
          if (method == "adaptive") {
            ds_test_function(
              gen_data           = gen_data,
              get_parameters     = get_parameters,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_norm_test   = fn_for_norm_test,
              fn_for_ds_test_1   = fn_for_ds_test_1,
              fn_for_ds_test_2   = fn_for_ds_test_2,
              paras              = paras,
              alpha              = alpha,
              norm_test_method   = norm_test_method,
              ...
            )
          } else if (method == "parametric") {
            fn_for_ds_test_1(data)$p.value
          } else {
            fn_for_ds_test_2(data, ...)$p.value
          }
        }, error = function(e) NA)
        
        if (!is.na(p_value) && p_value < alpha) {
          rejections <- rejections + 1
        }
      }
      
      ds_test_results[i] <- rejections / Nsim
      cat("Completed n =", n, "| Method =", method, "| Power =", ds_test_results[i], "\n")
    }
    
    results[[method]] <- ds_test_results
  }
  
  return(results)
}

# ------------------------------------------------------------------------------
# ----------------- Main simulation: power/Type I error analysis ---------------
sizes <- c(10, 20, 30, 40, 50)
Nsim <- 1e3
alpha <- 0.05
norm_test <- "SW"

# Give the functions names
gen_data        <- gen_data
get_parameters  <- get_parameters
fn_get_norm_obj <- fn_get_norm_obj
fn_ds_test_1    <- fn_for_ds_test_1
fn_ds_test_2    <- fn_for_ds_test_2

ds_test_methods <- c("parametric", "nonparametric", "adaptive")

ds_test_results <- perform_ds_func(
  sample_sizes       = sizes,
  Nsim               = Nsim,
  alpha              = alpha,
  gen_data           = gen_data,
  get_parameters     = get_parameters,
  fn_to_get_norm_obj = fn_get_norm_obj,
  fn_for_norm_test   = normality_test,
  fn_for_ds_test_1   = fn_ds_test_1,
  fn_for_ds_test_2   = fn_ds_test_2,
  norm_test_method   = norm_test,
  ds_test_methods    = ds_test_methods
)

# Access individual test results
ds_test1_results  <- ds_test_results$parametric
ds_test2_results <- ds_test_results$nonparametric
ds_adaptive_results  <- ds_test_results$adaptive


auc_df <- data.frame(
  Method = c("Test 1", "Test 2", "Test 3(Adaptive)"),
  AUC = c(
    compute_area(sizes, ds_test1_results),
    compute_area(sizes, ds_test2_results),
    compute_area(sizes, ds_adaptive_results)
  )
)

results <- list(
  sizes = sizes,
  ds_test1_results = ds_test1_results,
  ds_test2_results = ds_test2_results,
  ds_adaptive_results = ds_adaptive_results,
  auc = auc_df
)


# ------------------------------------------------------------------------------
# ----------------- Plot Power/Type I error -----------------------
# Plot all three curves
plot(sizes, ds_test1_results,    
     type = "b", 
     pch = 19,
     ylim = c(0,1), 
     xlab = "Sample Size", 
     ylab = "Power", 
     col = "red",
     main = "ANOVA Test Power Comparison")
lines(sizes, ds_test2_results, type = "b", pch = 17, col = "blue")
lines(sizes, ds_adaptive_results,   type = "b", pch = 15, col = "green")

legend("bottomright",
       legend = c("Parametric (Test 1)", "Nonparametric (Test 2)", "Adaptive (Test 3)"),
       col = c("red", "blue", "green"),
       pch = c(19, 17, 15),
       lty = 1,
       title = "Method")

# =============================================================================
# -----------------------------------------------------------------------------
# Function to compute FPR and TPR for Normality Test Methods
# -----------------------------------------------------------------------------
fn_for_roc_curve_for_norm_test <- function(alpha_pretest, tests, Nsim = 1e3) {
  FPR <- TPR <- matrix(0, nrow = length(tests), ncol = length(alpha_pretest))
  rownames(FPR) <- rownames(TPR) <- tests
  colnames(FPR) <- colnames(TPR) <- paste0("alpha_", alpha_pretest)
  
  for (i in seq_along(tests)) {
    test_name <- tests[i]
    
    for (j in seq_along(alpha_pretest)) {
      alpha <- alpha_pretest[j]
      reject_H0 <- reject_H1 <- numeric(Nsim)
      
      for (k in 1:Nsim) {
        # Under H0 
        paras_H0 <- get_parameters(n, dist = "Normal") 
        normal_data <- do.call(gen_data, paras_H0)
        # Under H1 
        paras_H1 <- get_parameters(n, dist = "Uniform") 
        non_normal_data <- do.call(gen_data, paras_H1)
        # Get normality test object
        normal_data <- fn_get_norm_obj(normal_data)
        non_normal_data <- fn_get_norm_obj(non_normal_data)
        # perform normality test for H0 and H1
        pvals_H0 <- normality_test(normal_data, test = test_name, alpha = alpha)$p_values
        reject_H0[k] <- any(pvals_H0 < alpha, na.rm = TRUE)
        
        pvals_H1 <- normality_test(non_normal_data, test = test_name, alpha = alpha)$p_values
        reject_H1[k] <- any(pvals_H1 < alpha, na.rm = TRUE)
        
      }
      
      # calculate FPR & TPR
      FPR[i, j] <- mean(reject_H0, na.rm = TRUE)
      TPR[i, j] <- mean(reject_H1, na.rm = TRUE)
    }
  }
  
  return(list(FPR = FPR, TPR = TPR, alpha = alpha_pretest))
}

# --------------- Run the simulation -------------------
# 
roc_pval_ds_test <- fn_for_roc_curve_for_norm_test(
  alpha_pretest = seq(from = 0.01, to = 1, by = 0.05),
  tests = c("SW", "SF", "JB"),
  Nsim = 1e3
)

# ------------------------------------------------------------------------------
#              Function to plot ROC curves from FPR and TPR matrices
# ------------------------------------------------------------------------------
roc_curve_for_norm_test <- function(FPR, TPR, tests_to_plot = rownames(FPR), alpha = NULL,
                                    title = "ROC Curves for Different Normality Tests") {
  # define line colors and plot characters
  colors <- 1:length(tests_to_plot)
  plot_chars <- 1:length(tests_to_plot)
  
  # create an empty plot
  plot(0, 0, type = "n", 
       xlim = c(0, 1), 
       ylim = c(0, 1),
       xlab = "False Positive Rate (FPR)", 
       ylab = "True Positive Rate (TPR)",
       main = title)
  # add lines for each test
  for (i in seq_along(tests_to_plot)) {
    test <- tests_to_plot[i]
    lines(FPR[test, ], TPR[test, ], col = colors[i], lwd = 2)
    points(FPR[test, ], TPR[test, ], col = colors[i], pch = plot_chars[i], cex = 0.5)
  }
  
  # add reference line
  abline(0, 1, lty = 2, col = "gray")
  legend("bottomright", 
         legend = tests_to_plot, 
         col = colors, 
         pch = plot_chars,
         lwd = 2, 
         title = "Normality Tests")
}

# Plot ROC using selected tests
selected_tests <- c("SW", "SF", "JB")
roc_curve_for_norm_test(FPR = roc_pval_ds_test$FPR,
         TPR = roc_pval_ds_test$TPR,
         tests_to_plot = selected_tests,
         alpha = roc_pval_ds_test$alpha)


# =============================================================================
# -----------------------------------------------------------------------------
# Function to compute p-values for normality test & downstream tests
# -----------------------------------------------------------------------------
generate_pval<- function(N, test = "SW", ...) {
  # Initialize storage 
  pval_t.test_H0 <- pval_wilcox.test_H0 <- numeric(N)
  pval_t.test_H1 <- pval_wilcox.test_H1 <- numeric(N)
  norm_pvals_H0 <- norm_pvals_H1 <- vector("list", N)  
  
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  
  for (i in 1:N) {
    # Under null hypothesis 
    paras_H0 <- get_parameters(n) 
    data_H0 <- do.call(gen_data, paras_H0)
    
    # Under alternative hypothesis 
    paras_H1 <- get_parameters(n, means = c(0.0, 0.2, 0.5, 0.0, 0.8))  
    data_H1 <- do.call(gen_data, paras_H1)
    
    # Get normality test objects
    normality_test_object_H0 <- fn_get_norm_obj(data_H0)
    normality_test_object_H1 <- fn_get_norm_obj(data_H1)
    
    # perform normality test
    normality_test_H0 <- normality_test(normality_test_object_H0, test = test, alpha = 0.05)
    normality_test_H1 <- normality_test(normality_test_object_H1, test = test, alpha = 0.05)
    
    # Store normality p-values
    norm_pvals_H0[[i]] <- normality_test_H0$p_values
    norm_pvals_H1[[i]] <- normality_test_H1$p_values
    
    # Get test p-values under null and alternative
    pval_t.test_H0[i] <- fn_ds_test_1(data_H0)$p.value
    pval_wilcox.test_H0[i] <- fn_ds_test_2(data_H0)$p.value
    
    pval_t.test_H1[i] <- fn_ds_test_1(data_H1)$p.value
    pval_wilcox.test_H1[i] <- fn_ds_test_2(data_H1)$p.value
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  
  return(list(
    pval_t.test_H0 = pval_t.test_H0,
    pval_wilcox.test_H0 = pval_wilcox.test_H0,
    pval_t.test_H1 = pval_t.test_H1,
    pval_wilcox.test_H1 = pval_wilcox.test_H1,
    norm_pvals_H0 = norm_pvals_H0,
    norm_pvals_H1 = norm_pvals_H1
  ))
}

# --------------- Calculate power/Type I errors each test method --------------
perform_analysis <- function(N, distributions, alpha_pretest, test_alpha) {
  ds_test_results <- list()
  error_ds_test <- list()
  power_ds_test <- list()
  
  for(dist in distributions) {
    cat("Processing distribution:", dist, "\n")
    ds_test_results[[dist]] <- generate_pval(N, test = "SW")
    
    # Calculate Type I error rates (under H0)
    error_ds_test[[dist]] <- list(
      error_t.test = mean(ds_test_results[[dist]]$pval_t.test_H0 < test_alpha),
      error_wilcox.test = mean(ds_test_results[[dist]]$pval_wilcox.test_H0 < test_alpha),
      adaptive_wilcox_error = numeric(length(alpha_pretest)))
    
    # Calculate Power (under H1)
    power_ds_test[[dist]] <- list(
      power_t.test = mean(ds_test_results[[dist]]$pval_t.test_H1 < test_alpha),
      power_wilcox.test = mean(ds_test_results[[dist]]$pval_wilcox.test_H1 < test_alpha),
      adaptive_wilcox_power = numeric(length(alpha_pretest)))
    
    # Pre-compute decisions for each alpha level
    for(j in seq_along(alpha_pretest)) {
      alpha <- alpha_pretest[j]
      
      # For Type I error (H0)
      use_t_test_H0 <- sapply(ds_test_results[[dist]]$norm_pvals_H0, function(x) all(x > alpha))
      adaptive_pvals_H0 <- ifelse(use_t_test_H0,
                                  ds_test_results[[dist]]$pval_t.test_H0,
                                  ds_test_results[[dist]]$pval_wilcox.test_H0)
      error_ds_test[[dist]]$adaptive_wilcox_error[j] <- mean(adaptive_pvals_H0 < test_alpha)
      
      # For Power (H1)
      use_t_test_H1 <- sapply(ds_test_results[[dist]]$norm_pvals_H1, function(x) all(x > alpha))
      adaptive_pvals_H1 <- ifelse(use_t_test_H1,
                                  ds_test_results[[dist]]$pval_t.test_H1,
                                  ds_test_results[[dist]]$pval_wilcox.test_H1)
      power_ds_test[[dist]]$adaptive_wilcox_power[j] <- mean(adaptive_pvals_H1 < test_alpha)
    }
  }
  
  return(list(
    error_ds_test = error_ds_test,
    power_ds_test = power_ds_test
  ))
}

# run analysis to get power and error 
analysis_ds_tests <- perform_analysis(
  N = 1e4,
  distributions = c("Normal", "LogNormal"),
  alpha_pretest = seq(from = 0.001, to = 0.1, by = 0.01),
  test_alpha = 0.05
)

# ------------------------------------------------------------------------------
#                         Function to compute ROC-like metrics
# ------------------------------------------------------------------------------
compute_roc_metrics <- function(error_ds_test, power_ds_test, test_alpha) {
 
  # Non-normal case
  power_t_test_nonnormal <- power_ds_test[[ 2 ]]$power_t.test
  adaptive_power_nonnormal <- power_ds_test[[ 2 ]]$adaptive_wilcox_power
  adaptive_error_nonnormal <- error_ds_test[[ 2 ]]$adaptive_wilcox_error
  EPG <- adaptive_power_nonnormal - power_t_test_nonnormal
  EDE <- adaptive_error_nonnormal - test_alpha
  
  # normal case
  power_t_test_normal <- power_ds_test[[ 1 ]]$power_t.test
  adaptive_power_normal <- power_ds_test[[ 1 ]]$adaptive_wilcox_power
  adaptive_error_normal <- error_ds_test[[ 1 ]]$adaptive_wilcox_error
  EPL <- power_t_test_normal - adaptive_power_normal
  EIE <- adaptive_error_normal - test_alpha
  
  # Point estimates for benchmark comparison
  power_gain <- power_ds_test[[ 2 ]]$power_wilcox.test - power_ds_test[[ 2 ]]$power_t.test
  power_loss <- power_ds_test[[ 1 ]]$power_t.test - power_ds_test[[ 1 ]]$power_wilcox.test
  
  return(list(EPL = EPL,
              EPG = EPG,
              EIE = EIE,
              EDE = EDE,
              power_gain = power_gain,
              power_loss = power_loss)
         )
}

# Compute metrics after running analysis
metrics <- compute_roc_metrics(
  error_ds_test = analysis_ds_tests$error_ds_test,
  power_ds_test = analysis_ds_tests$power_ds_test,
  test_alpha    = 0.05
)

# --------------------------------------------------------------------
# Power and Type I error trade-off plots
plot_power_error_tradeoff <- function(alpha_pretest, metrics, file_name = NULL) {
 # if (!is.null(file_name)) pdf(file_name, width = 12, height = 10)
  par(mfrow = c(2, 2))
  plot(alpha_pretest, 
       metrics$EPL, 
       type = "l", 
       col = "blue", 
       lwd = 2,
       ylab = "Expected Power Loss (EPL)", 
       xlab = expression(alpha),
       main = "Power Loss (Normal)")
  
  plot(alpha_pretest, 
       metrics$EPG, 
       type = "l", 
       col = "red", 
       lwd = 2,
       ylab = "Expected Power Gain (EPG)", 
       xlab = expression(alpha),
       main = "Power Gain (LogNormal)")
  
  plot(alpha_pretest, 
       metrics$EIE, 
       type = "l", 
       col = "orange", 
       lwd = 2,
       ylab = "Type I Error Inflation", 
       xlab = expression(alpha),
       main = "Inflation (Normal)")
  
  plot(alpha_pretest, 
       metrics$EDE, 
       type = "l", 
       col = "green", 
       lwd = 2,
       ylab = "Type I Error Inflation", 
       xlab = expression(alpha),
       main = "Inflation (LogNormal)")
  if (!is.null(file_name)) dev.off()
}

# Generate all plots
plot_power_error_tradeoff(
  alpha_pretest = seq(from = 0.001, to = 0.1, by = 0.01),
  metrics = metrics
  #file_name = "comparison_power_error.pdf"
)

# --------------------------------------------------------------------
# ----------- Function to ROC curve data for downstream test ---------
generate_roc_tables <- function(N, distributions, sig_levels) {
  # data frames for Type I error and Power
  roc_data <- data.frame()
  
  for (dist in distributions) {
    cat("Generating ROC data for:", dist, "\n")
    res <- generate_pval(N, test = "SW", dist)
    
    for (alpha in sig_levels) {
      # Compute Type I error (under H0)
      error_t <- mean(res$pval_t.test_H0 < alpha)
      error_w <- mean(res$pval_wilcox.test_H0 < alpha)
      
      # Compute Power (under H1)
      power_t <- mean(res$pval_t.test_H1 < alpha)
      power_w <- mean(res$pval_wilcox.test_H1 < alpha)
      
      # Append results
      roc_data <- rbind(roc_data, 
                        data.frame(
                          Distribution = dist,
                          Method = "test 1",
                          Alpha = alpha,
                          Power = power_t,
                          TypeIError = error_t
                        ),
                        data.frame(
                          Distribution = dist,
                          Method = "test 2",
                          Alpha = alpha,
                          Power = power_w,
                          TypeIError = error_w
                        ))
    }
  }
  
  return(roc_data)
}

# Example usage:
roc_data <- generate_roc_tables(
  N = 1e3, 
  distributions = c("Normal", "LogNormal"), 
  sig_levels = seq(from = 0.01, to = 1, by = 0.05)
)

# --------------------- create the plots -------------------------------------
plot_power_vs_error_ROC_curve <- function(roc_data, file_name) {
  pdf(file_name, width = 10, height = 6)
  par(mfrow = c(1, length(unique(roc_data$Distribution))))
  
  methods <- unique(roc_data$Method)
  colors  <- c("blue", "red")
  pch_vals <- c(19, 17)
  
  for (dist in unique(roc_data$Distribution)) {
    plot(NA, xlim = c(0, 1), 
         ylim = c(0, 1),
         xlab = "Type I Error", 
         ylab = "Power",
         main = paste("ROC-like Curve (", dist, ")", sep = ""))
    
    for (m in seq_along(methods)) {
      method <- methods[m]
      data_subset <- subset(roc_data, Distribution == dist & Method == method)
      if (nrow(data_subset) > 0) {
        lines(data_subset$TypeIError, data_subset$Power, type = "l", col = colors[m], lwd = 2, pch = pch_vals[m])
      }
    }
    
    legend("bottomright", legend = methods, col = colors,
           lwd = 2, pch = pch_vals, title = "Method")
  }
  dev.off()
}


plot_power_vs_error_ROC_curve(
  roc_data = roc_data,
 file_name = "Power_vs_TypeIError_ROC_By_Distribution.pdf"
)

