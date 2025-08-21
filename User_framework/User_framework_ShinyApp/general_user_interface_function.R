# Load helper source files
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
setwd("~/Desktop/OSU/Research/Pretest-Simulation/User_interface/User_framework_ShinyApp")
# load user defined functions

FilePath <- "~/Desktop/OSU/Research/Research/user_framework/ShinyApp/sample_test_functions/t_test_functions/"
source(paste0(FilePath, "gen_data.R"))
source(paste0(FilePath, "get_parameters.R"))
source(paste0(FilePath, "fn_get_norm_obj.R"))
source(paste0(FilePath, "fn_for_ds_test_1.R"))
source(paste0(FilePath, "fn_for_ds_test_2.R"))

#---------- function for computing auc ------------
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
  # generate data set. Accept parameters if provided
  data <- if (!is.null(paras)) do.call(gen_data, paras) else gen_data()
  
  # get normality test object
  normality_test_object     <- fn_to_get_norm_obj(data)
  normality_test_pval  <- fn_for_norm_test(data = normality_test_object, test = norm_test_method, alpha = alpha)
  
  # choose ds test method based normality test results
  if (isTRUE(normality_test_pval$normality_satisfied)) {
    ds_test <- fn_for_ds_test_1(data)
  } else {
    ds_test<- fn_for_ds_test_2(data, ...)
  }
  
  return(ds_test$p.value)
}

# ----------------- Power/Type I error Analysis Test Function ---------------
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
    effect_size        = NULL,
    ...
) {
  ds_test_methods <- match.arg(ds_test_methods, several.ok = TRUE)
  results <- list()
  timing <- list()  # To store time information
  pval_storage <- list()   # initialize p-value storage
  
  for (method in ds_test_methods) {
    cat("Running method:", method, "\n")
    # store ds test results for each sample size
    ds_test_results <- numeric(length(sample_sizes))
    # store p_vlaues for each method for each sample size
    pval_storage[[method]] <- vector("list", length(sample_sizes))
    names(ds_test_results) <- paste0("n=", sample_sizes)
    # track run progress with progress bar
    pb <- txtProgressBar(min = 0, max = length(sample_sizes), style = 3)
    
    start_time <- Sys.time()  # Start timer for this method
    
    for (i in seq_along(sample_sizes)) {
      n <- sample_sizes[i]
      rejections <- 0
      pvals <- numeric(Nsim) # store p-values for this n
      for (sim in seq_len(Nsim)) {
        # get user-supplied parameters
        paras <- get_parameters(n, par = 7, effect_size = effect_size)
        data <- if (!is.null(paras)) do.call(gen_data, paras) else gen_data()
        # calculate ds p-values
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
        
        pvals[sim] <- p_value   # store each simulated p-value
        
        if (!is.na(p_value) && p_value < alpha) {
          rejections <- rejections + 1
        }
      }
      
      ds_test_results[i] <- rejections / Nsim
      pval_storage[[method]][[i]] <- pvals   # store all p-values
     
      setTxtProgressBar(pb, i)
  
    }
    close(pb)
    end_time <- Sys.time()  # End timer for this method
    timing[[method]] <- as.numeric(difftime(end_time, start_time, units = "secs"))
    
    results[[method]] <- ds_test_results
    cat("Method", method, "completed in", round(timing[[method]], 2), "seconds\n")
  }
  
  return(list(results = results, pvalues = pval_storage, timing = timing))
}

# -------------------------------------------------------
# ------- Main simulation: power/Type I error analysis --
Nsim <- 1e4
sizes <- c(10, 20, 30, 40, 50)
Nsim <- Nsim
alpha <- 0.05
norm_test <- "SW"
effect_size =  0.75

# test methods
ds_test_methods <- c("parametric", "nonparametric", "adaptive")

# run test
sim_output <- perform_ds_func(  
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
  ds_test_methods    = ds_test_methods,
  effect_size        = effect_size 
)

# Extract results and timing
ds_test_results <- sim_output$results
timing_results <- sim_output$timing  
pvals_all       <- sim_output$pvalues   

# Access individual test results
ds_test1_results  <- ds_test_results$parametric
ds_test2_results <- ds_test_results$nonparametric
ds_adaptive_results  <- ds_test_results$adaptive


auc_df <- data.frame(
  Method = c("Test 1(t_test)", "Test 2(Wilcoxon)", "Test 3(Adaptive)"),
  AUC = c(
    compute_area(sizes, ds_test1_results),
    compute_area(sizes, ds_test2_results),
    compute_area(sizes, ds_adaptive_results)
  ),
  Runtime = c(
    timing_results$parametric,
    timing_results$nonparametric,
    timing_results$adaptive
  )
)

results <- list(
  sizes = sizes,
  ds_test1_results = ds_test1_results,
  ds_test2_results = ds_test2_results,
  ds_adaptive_results = ds_adaptive_results,
  auc = auc_df,
  timing = timing_results,
  pvalues   = pvals_all    
)


# ------------------------------------------------------------------------------
# ----------------- Plot Power/Type I error -----------------------
# Plot all three curves
pdf("t_test_power_chisq.pdf", width=12, height=8)
plot(sizes, ds_test1_results,    
     type = "b", 
     pch = 19,
     ylim = c(0,1), 
     xlab = "Sample Size", 
     ylab = "Power", 
     col = "red",
     main = "Power Comparison for two sample location tests(chisq)")
lines(sizes, ds_test2_results, type = "b", pch = 17, col = "blue")
lines(sizes, ds_adaptive_results,   type = "b", pch = 15, col = "green")

legend("bottomright",
       legend = paste0(c("Test1(Parametric)","Test2(Nonparametric)","Test3(Adaptive)"),
                       " (AUC=", formatC(auc_df$AUC, format = "f", digits = 4),")"),
       col = c("red", "blue", "green"),
       pch = c(19, 17, 15),
       lty = 1,
       title = "Method")

dev.off()
# =============================================================================
# -----------------------------------------------------------------------------
# Function to compute FPR and TPR for Normality Test Methods
# -----------------------------------------------------------------------------
fn_for_roc_curve_for_norm_test <- function(n, alpha_pretest, H1_dist, tests, Nsim = 1e3) {
  FPR <- TPR <- matrix(0, nrow = length(tests), ncol = length(alpha_pretest))
  rownames(FPR) <- rownames(TPR) <- tests
  colnames(FPR) <- colnames(TPR) <- paste0("alpha_", alpha_pretest)
  
  pb <- txtProgressBar(min = 0, max = Nsim * length(tests) * length(alpha_pretest), style = 3)
  counter <- 0
  
  for (i in seq_along(tests)) {
    test_name <- tests[i]
    
    for (j in seq_along(alpha_pretest)) {
      alpha <- alpha_pretest[j]
      reject_H0 <- reject_H1 <- numeric(Nsim)
      
      for (k in 1:Nsim) {
        # Under H0 
        paras_H0 <- get_parameters(n, dist = "Normal", par = NULL) 
        normal_data <- do.call(gen_data, paras_H0)
        # Under H1 
        paras_H1 <- get_parameters(n, dist = H1_dist, par = 3) 
        non_normal_data <- do.call(gen_data, paras_H1)
        # Get normality test object
        normal_data <- fn_get_norm_obj(normal_data)
        non_normal_data <- fn_get_norm_obj(non_normal_data)
        # perform normality test for H0 and H1
        pvals_H0 <- normality_test(normal_data, test = test_name, alpha = alpha)$p_values
        reject_H0[k] <- any(pvals_H0 < alpha, na.rm = TRUE)
        
        pvals_H1 <- normality_test(non_normal_data, test = test_name, alpha = alpha)$p_values
        reject_H1[k] <- any(pvals_H1 < alpha, na.rm = TRUE)
        
        counter <- counter + 1
        setTxtProgressBar(pb, counter)
      }
      
      # calculate FPR & TPR
      FPR[i, j] <- mean(reject_H0, na.rm = TRUE)
      TPR[i, j] <- mean(reject_H1, na.rm = TRUE)
    }
  }
  close(pb)
  
  return(list(
            pvals_H0 = pvals_H0,
            pvals_H1 = pvals_H1,
            FPR = FPR, 
            TPR = TPR, 
            alpha = alpha_pretest))
}

# --------------- Run the simulation -------------------
Nsim = 1e4
alpha_pretest = seq(from = 0.0, to = 1, by = 0.05)
roc_pval_ds_test <- fn_for_roc_curve_for_norm_test(
  n = 10,
  alpha_pretest = alpha_pretest,
  H1_dist = "chi_square",
  tests = c("SW", "SF", "JB","SKEW") #, "DAP", "AD", "CVM"),
  ,
  Nsim = Nsim
)

# ------------------------------------------------------------------------------
#              Function to plot ROC curves from FPR and TPR matrices
# ------------------------------------------------------------------------------
roc_curve_for_norm_test <- function(FPR, TPR, tests_to_plot = rownames(FPR), alpha = NULL,
                                    title = "ROC Curves for Different Normality Tests(Chisq_3)") {
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
pdf("t_test_norm.pdf", width=12, height=8)
selected_tests <- c("SW", "SF", "JB", "SKEW")
roc_curve_for_norm_test(FPR = roc_pval_ds_test$FPR,
                        TPR = roc_pval_ds_test$TPR,
                        tests_to_plot = selected_tests,
                        alpha = roc_pval_ds_test$alpha)
dev.off()

# =============================================================================
# -----------------------------------------------------------------------------
# Function to compute p-values for normality test & downstream tests
# -----------------------------------------------------------------------------
generate_pval<- function(N, n, effect_size, test = "SW", dist = "Normal", ...) {
  # Initialize storage 
  pval_t.test_H0 <- pval_wilcox.test_H0 <- numeric(N)
  pval_t.test_H1 <- pval_wilcox.test_H1 <- numeric(N)
  norm_pvals_H0 <- norm_pvals_H1 <- vector("list", N)  
  
  pb <- txtProgressBar(min = 0, max = N, style = 3)
  
  for (i in 1:N) {
    # Under null hypothesis 
    paras_H0 <- get_parameters(n, dist = dist) 
    data_H0 <- do.call(gen_data, paras_H0)
    
    # Under alternative hypothesis 
    paras_H1 <- get_parameters(n, dist = dist, effect_size = effect_size)  
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

# --------------- Calculate power/Type I errors for each test method --------------
perform_analysis <- function(N, n, distributions, effect_size, test, alpha_pretest, test_alpha) {
  ds_test_results <- list()
  error_ds_test <- list()
  power_ds_test <- list()
  
  pb_dist <- txtProgressBar(min = 0, max = length(distributions), style = 3)
  
  for(dist_idx in seq_along(distributions)) {
    dist <- distributions[dist_idx]
    cat("Processing distribution:", dist, "\n")
    
    # Store all p-values from generate_pval
    ds_test_results[[dist]] <- generate_pval(N, n, effect_size = effect_size, test = "SW" , dist = dist)
    
    # Calculate Type I error rates (under H0)
    error_ds_test[[dist]] <- list(
      error_t.test = mean(ds_test_results[[dist]]$pval_t.test_H0 < test_alpha),
      error_wilcox.test = mean(ds_test_results[[dist]]$pval_wilcox.test_H0 < test_alpha),
      # storage for adaptive test for error
      adaptive_wilcox_error = numeric(length(alpha_pretest)))
    
    # Calculate Power (under H1)
    power_ds_test[[dist]] <- list(
      power_t.test = mean(ds_test_results[[dist]]$pval_t.test_H1 < test_alpha),
      power_wilcox.test = mean(ds_test_results[[dist]]$pval_wilcox.test_H1 < test_alpha),
      # storage for adaptive test for power
      adaptive_wilcox_power = numeric(length(alpha_pretest)))
    
    pb_alpha <- txtProgressBar(min = 0, max = length(alpha_pretest), style = 3)
    
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
      
      setTxtProgressBar(pb_alpha, j)
    }
    close(pb_alpha)
    setTxtProgressBar(pb_dist, dist_idx)
  }
  close(pb_dist)
  
  return(list(
    error_ds_test = error_ds_test,
    power_ds_test = power_ds_test,
    all_pvalues     = ds_test_results   # raw p-values saved here
  ))
}

# run analysis to get power and error
Nsim <- 1e6
alpha_pretest = seq(from = 0.005, to = 1, by = 0.005)
analysis_ds_tests <- perform_analysis(
  N = Nsim,
  n = 10,
  distributions = c("Normal", "lognormal"),
  effect_size = effect_size,
  test = "SW",
  alpha_pretest = alpha_pretest,
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
       main = "EPL for samples from (Normal)")
  
  plot(alpha_pretest, 
       metrics$EPG, 
       type = "l", 
       col = "red", 
       lwd = 2,
       ylab = "Expected Power Gain (EPG)", 
       xlab = expression(alpha),
       main = "EPG for samples from Chi-square")
  
  plot(alpha_pretest, 
       metrics$EIE, 
       type = "l", 
       col = "orange", 
       lwd = 2,
       ylab = "Inflation of Type I Error", 
       xlab = expression(alpha),
       main = "EI/D of Type I error (Normal)")
  
  plot(alpha_pretest, 
       metrics$EDE, 
       type = "l", 
       col = "green", 
       lwd = 2,
       ylab = "Type I Error Inflation", 
       xlab = expression(alpha),
       main = "EI/D of Type I error (Chi-square)")
  if (!is.null(file_name)) dev.off()
}

# Generate all plots
pdf("t_test_power_tradeoff.pdf", width=12, height=8)
plot_power_error_tradeoff(
  alpha_pretest = alpha_pretest,
  metrics = metrics
  #file_name = "comparison_power_error.pdf"
)
dev.off()

EPL_vs_EPG_ROC_like_curve <- function(metrics){
  par(mfrow = c(1, 2))
  plot(metrics$EPL,
       metrics$EPG, 
       type = "l", 
       col = "blue", 
       lwd = 2,
       ylab = "EPG", 
       xlab = "EPL",
       main = "ROC like curve: EPG vs EPL")
  
  plot(metrics$EDE,
       metrics$EIE, 
       type = "l", 
       col = "blue", 
       lwd = 2,
       ylab = "EIE", 
       xlab = "EDE",
       main = "ROC like curve: EIE vs EDE")
}

pdf("t_test_EPL_vs_EPG.pdf", width=12, height=8)
EPL_vs_EPG_ROC_like_curve(metrics)

dev.off()

# --------------------------------------------------------------------
# ----------- Function to ROC curve data for downstream test ---------
generate_roc_tables <- function(N, n, distributions, effect_size, sig_levels) {
  # data frames for Type I error and Power
  roc_data <- data.frame()
  
  pb <- txtProgressBar(min = 0, max = length(distributions) * length(sig_levels), style = 3)
  counter <- 0
  
  for (dist in distributions) {
    cat("Generating ROC data for:", dist, "\n")
    res <- generate_pval(N, n,  effect_size, test = "SW", dist = dist)
    
    # Convert norm_pvals_H0 from list to numeric vector
    norm_pvals_H0_vec <- sapply(res$norm_pvals_H0, function(x) if(is.null(x)) NA else x)
    
    for (alpha in sig_levels) {
      # Compute Type I error (under H0)
      error_t <- mean(res$pval_t.test_H0 < alpha, na.rm = TRUE)
      error_w <- mean(res$pval_wilcox.test_H0 < alpha, na.rm = TRUE)
      
      # For adaptive test 
      use_t_test <- !is.na(norm_pvals_H0_vec) & norm_pvals_H0_vec > 0.05
      adaptive_pvals_H0 <- ifelse(use_t_test, 
                                  res$pval_t.test_H0, 
                                  res$pval_wilcox.test_H0)
      error_adaptive <- mean(adaptive_pvals_H0 < alpha, na.rm = TRUE)
      
      # Compute Power (under H1)
      power_t <- mean(res$pval_t.test_H1 < alpha, na.rm = TRUE)
      power_w <- mean(res$pval_wilcox.test_H1 < alpha, na.rm = TRUE)
      
      adaptive_pvals_H1 <- ifelse(use_t_test,
                                  res$pval_t.test_H1,
                                  res$pval_wilcox.test_H1)
      power_adaptive <- mean(adaptive_pvals_H1 < alpha, na.rm = TRUE)
      
      # Append results
      roc_data <- rbind(roc_data, 
                        data.frame(
                          Distribution = dist,
                          Method = "test 1(t-test)",
                          Alpha = alpha,
                          Power = power_t,
                          TypeIError = error_t
                        ),
                        data.frame(
                          Distribution = dist,
                          Method = "test 2(Wilcoxon test)",
                          Alpha = alpha,
                          Power = power_w,
                          TypeIError = error_w
                        )
                        ,
                        data.frame(
                          Distribution = dist,
                          Method = "adaptive test",
                          Alpha = alpha,
                          Power = power_adaptive,
                          TypeIError = error_adaptive
                        ))
      
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
    }
  }
  close(pb)
  
  return(roc_data)
}

# Example usage:
Nsim <- 1e5
roc_data <- generate_roc_tables(
  N = Nsim, 
  n = 10,
  distributions = c("Normal", "chi_square"), 
  effect_size = effect_size,
  sig_levels = seq(from = 0.001, to = 1, by = 0.005)
)

# --------------------- create the plots -------------------------------------
plot_power_vs_error_ROC_curve <- function(roc_data, file_name = NULL) {
  #if (!is.null(file_name)) pdf(file_name, width = 12, height = 10)
  par(mfrow = c(1, length(unique(roc_data$Distribution))))
  
  methods <- unique(roc_data$Method)
  colors  <- c("blue", "red", "green")
  pch_vals <- c(19, 17, 18)
  
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
  #   dev.off()
}

pdf("t_test_power_error.pdf", width=12, height=8)
plot_power_vs_error_ROC_curve(
  roc_data = roc_data
  #file_name = "Power_vs_TypeIError_ROC_By_Distribution.pdf"
)

dev.off()

# Save workspace with timestamp
save.image(paste0("t_test_results_", format(Sys.time(), "%Y%m%d_%H%M%S"), ".RData"))
