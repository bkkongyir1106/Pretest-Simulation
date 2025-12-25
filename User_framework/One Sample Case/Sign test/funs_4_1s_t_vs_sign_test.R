# ----------------------------------------------------
#   Generate standardized data from various distributions
#                 Center around median 
# ----------------------------------------------------
generate_data <- function(n, dist, par = NULL) {
  # Input validation
  if (!is.numeric(n) || n <= 0) stop("n must be a positive integer")
  dist <- tolower(dist)
  
  if (dist == "normal") {
    if (is.null(par)) par <- c(0, 1)
    x <- rnorm(n, mean = par[1], sd = par[2])
    # Center around median (same as mean for normal)
    x <- x - median(x)
    
  } else if (dist == "chi_square") {
    if (is.null(par)) par <- 3
    x <- rchisq(n, df = par)
    median_chi <- par *(1-2/(9*par))^3
    sd_chi <- sqrt(2*par)
    x <- (x - median_chi)/sd_chi 
    
  } else if (dist == "gamma") {
    if (is.null(par)) par <- c(3, 0.1)
    x <- rgamma(n, shape = par[1], rate = par[2])
    median_g <- qgamma(0.5, shape = par[1], rate = par[2])
    sd_g <- sqrt(par[1] / par[2]^2)
    x <- (x - median_g)/sd_g 
    
  } else if (dist == "exponential") {
    if (is.null(par)) par <- 1
    x <- rexp(n, rate = par)
    median_e <- log(2)/par
    sd_e <- 1 / par
    x <- (x - median_e)/sd_e
    
  } else if (dist == "t") {
    if (is.null(par)) par <- 3
    if (par <= 2) stop("Degrees of freedom must be > 1")
    x <- rt(n, df = par)
    x <- x / sqrt(par / (par - 2))
    
  } else if (dist == "uniform") {
    if (is.null(par)) par <- c(0, 1)
    x <- runif(n, min = par[1], max = par[2])
    median_u <- (par[1] + par[2]) / 2
    sd_u <- sqrt((par[2] - par[1])^2 / 12)
    x <- (x - median_u)/sd_u
    
  } else if (dist == "laplace") {
    if (is.null(par)) par <- c(3, 1)
    x <- extraDistr::rlaplace(n, mu = par[1], sigma = par[2])
    x <- (x - par[1])/sqrt(2 * par[2]^2)
  
  } else if (dist == "weibull") {
    if (is.null(par)) par <- c(1, 2)
    x <- rweibull(n, shape = par[1], scale = par[2])
    median_w <- par[2] * (log(2))^(1/par[1])
    var_w <- par[2]^2 * (gamma(1 + 2 / par[1]) - gamma(1 + 1 / par[1])^2)
    x <- (x - median_w) / sqrt(var_w)
    
  } else if (dist == "lognormal") {
    if (is.null(par)) par <- c(0, 1)
    x <- rlnorm(n, meanlog = par[1], sdlog = par[2])
    median_ln <- exp(par[1])
    var_ln <- (exp(par[2]^2) - 1) * exp(2 * par[1] + par[2]^2)
    x <- (x - median_ln)/sqrt(var_ln) 
    
   }  else {
    stop("Unsupported distribution: ", dist)
  }
  
  return(x)
}

# ---------------------------------------------------------------------------
# Define function to perform different Normality tests
# ---------------------------------------------------------------------------
generate_tests <- function(x, test){
  if(test == "KS"){  # Kolmogorov-Smirnov test (Lilliefors)
    output <- nortest::lillie.test(x)
    
  } else if(test == "SW"){  # Shapiro-Wilk
    output <- shapiro.test(x)
    
  } else if(test == "JB"){  # Jarque-Bera
    output <- tseries::jarque.bera.test(x)
    
  } else if(test == "DAP"){  # D'Agostino test (skewness)
    output <- moments::agostino.test(x)
    
  } else if(test == "ANS"){  # Anscombe-Glynn test (kurtosis)
    output <- moments::anscombe.test(x)
    
  } else if(test == "AD"){  # Anderson-Darling
    output <- nortest::ad.test(x)
    
  } else if(test == "AD2"){  # Anderson-Darling from DescTools (more options)
    output <- DescTools::AndersonDarlingTest(x)
    
  } else if(test == "SF"){  # Shapiro-Francia
    output <- nortest::sf.test(x)
    
  } else if(test == "CVM"){  # Cramer-von Mises
    output <- nortest::cvm.test(x)
    
  } else if(test == "SKEW"){  # Skewness test (z-test)
    s <- moments::skewness(x)
    se <- sqrt(6/length(x))
    z <- s / se
    p <- 2 * (1 - pnorm(abs(z)))
    output <- list(statistic = z, p.value = p, method = "Skewness z-test")
    
  } else if(test == "KURT"){  # Kurtosis test (z-test)
    k <- moments::kurtosis(x)
    se <- sqrt(24/length(x))
    z <- (k - 3) / se
    p <- 2 * (1 - pnorm(abs(z)))
    output <- list(statistic = z, p.value = p, method = "Kurtosis z-test")
    
  } else {
    stop("Unknown test name. Choose from: KS, SW, JB, DAP, ANS, AD, AD2, SF, CVM, CHISQ, SKEW, KURT.")
  }
  
  return(output)
}

# ------------------------------------
# Define simple calculation functions 
# -------------------------------------

# calculate one-sample  test statistic 
OneSample_test_statistic <- function(x) {
  return((mean(x)* sqrt(length(x))) / sd(x))
}

# compute the area under the curve 
compute_area <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / (max(x) - min(x))
}

# -----------------------------------------------------------------------------#
# Define functions for performing various downstream tests 
# -----------------------------------------------------------------------------#

# -----------One Sample Permutation test ------
one_sample_permutation_test <- function(x, B){
  observe_stat <- OneSample_test_statistic(x)
  permuted_stat <- numeric(B)
  for (j in 1:B) {
    index <- sample(c(-1, 1), length(x), replace = TRUE)
    sample_data <- index * abs(x)
    permuted_stat[j] <- OneSample_test_statistic(sample_data)
  }
  return(mean(abs(permuted_stat) >= abs(observe_stat)))
}

## ---------- Sign test for the median ---------
sign_test <- function(data, mu0 = 0) {
  x0 <- data - mu0
  n_valid <- sum(x0 != 0)
  if (n_valid == 0) return(list(p.value = 1))
  signs <- sum(x0 > 0)
  return(list(p.value = binom.test(signs, n_valid, p = 0.5)$p.value))
}

# -----------One Sample Bootstrap test ------
one_sample_bootstrap_test <- function(n, dist, test_stat, n_bootstrap, effect_size = 0){
  x <- generate_data(n, dist) + effect_size
  # observe test stats
  observe_test_stat <- test_stat(x)
  # bootstrap test stats
  x <- x - mean(x)
  bootstrap_test_stat <- numeric(n_bootstrap)
  for(b in 1 : n_bootstrap){
    sample1 <- sample(x, n, replace = TRUE)
    bootstrap_test_stat[b] <- test_stat(sample1)
    bootstrap_p_val <- mean(abs(bootstrap_test_stat) >= abs(observe_test_stat))
  }
  return(bootstrap_p_val)
}

# ------------ framework functions -----------------------------------------
# generate data 
gen_data <- function(n = 20, effect_size = 0.0, sd = 1, dist = "Exponential", par = NULL) {
  dist <- tolower(dist)
  x <- effect_size + sd * generate_data(n, dist, par = par)
  return(x)
}

# get parameters
get_parameters <- function(n, effect_size = 0.5, sd = 1, dist = "exponential", par = NULL, ...) {
  list(
    n = n,
    effect_size = effect_size,
    sd = sd,
    dist = dist,
    par = par
  )
}

# get normality object
fn_to_get_norm_obj <- function(data) {
  return(data)
}

# ds test function 1: one-sample t-test
fn_for_ds_test_1 <- function(data) {
  test_result <- t.test(data, mu = 0)
  return(list(p.value = test_result$p.value))
}


# # ds test function 2: Wilcoxon Sign Rank test 
fn_for_ds_test_2 <- function(data, mu0 = 0) {
  x0 <- data - mu0
  n_valid <- sum(x0 != 0)
  if (n_valid == 0) return(list(p.value = 1))
  signs <- sum(x0 > 0)
  return(list(p.value = binom.test(signs, n_valid, p = 0.5)$p.value))
}

# -----------------------------------------------------------------------------#
# ---------- General Normality check function ---------------------------------#
# -----------------------------------------------------------------------------#
normality_test <- function(data, test = "SW", alpha = 0.05) {
  pvals <- NULL
  all_normality_satisfied <- NULL
  
  # Case 1: Single numeric vector
  if (is.numeric(data) && is.atomic(data) && is.null(dim(data))) {
    pvals <- generate_tests(data, test = test)$p.value
    all_normality_satisfied <- pvals > alpha
  }
  
  # Case 2: List of numeric vectors
  else if (is.list(data) && !is.data.frame(data)) {
    pvals <- sapply(data, function(sample) {
      generate_tests(sample, test = test)$p.value
    })
    # name results
    names(pvals) <- if (is.null(names(data))) paste0("Sample", seq_along(pvals)) else names(data)
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

# =============================================================================
# -----------------------------------------------------------------------------
# Function to compute FPR and TPR for Normality Test Methods
# -----------------------------------------------------------------------------
fn_for_roc_curve_for_norm_test <- function(n, alpha_pretest, H1_dist, tests, Nsim = 1e3) {
  FPR <- TPR <- matrix(0, nrow = length(tests), ncol = length(alpha_pretest))
  rownames(FPR) <- rownames(TPR) <- tests
  colnames(FPR) <- colnames(TPR) <- paste0("alpha_", alpha_pretest)
  
  cat("Generating normality ROC data","\n")
  pb <- txtProgressBar(min = 0, max = Nsim * length(tests) * length(alpha_pretest), style = 3)
  counter <- 0
  
  for (i in seq_along(tests)) {
    test_name <- tests[i]
    
    for (j in seq_along(alpha_pretest)) {
      alpha <- alpha_pretest[j]
      reject_H0 <- reject_H1 <- numeric(Nsim)
      
      for (k in 1:Nsim) {
        # Under H0 
        paras_H0 <- get_parameters(n, dist = "Normal", par = NULL) # pick parameters
        normal_data <- do.call(gen_data, paras_H0)  # generate normal data
        # Under H1 
        paras_H1 <- get_parameters(n, dist = H1_dist)  # get parameters
        non_normal_data <- do.call(gen_data, paras_H1) # generate data 
        # Get normality test objects
        normal_data <- fn_to_get_norm_obj(normal_data)
        non_normal_data <- fn_to_get_norm_obj(non_normal_data)
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

# -------------------------------------------------------------------------
# ----------------
# specialized trapezoidal auc function for Normality tests
# ---------------
compute_auc <- function(fpr, tpr, ensure_endpoints = TRUE, normalize = FALSE) {
  # sanity checks
  ok <- is.finite(fpr) & is.finite(tpr)
  x <- fpr[ok]
  y <- tpr[ok]
  
  # sort by x 
  o <- order(x, y)
  x <- x[o]
  y <- y[o]
  
  if (length(x) < 2) return(NA_real_)
  
  # include (0,0) and (1,1)
  if (ensure_endpoints) {
    if (x[1] > 0 || y[1] > 0) { x <- c(0, x); y <- c(0, y) }
    if (x[length(x)] < 1 || y[length(y)] < 1) { x <- c(x, 1); y <- c(y, 1) }
  }
  
  # trapezoidal rule 
  area <- sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)
  
  # normalization to [0,1] width
  if (normalize) {
    rng <- max(x) - min(x)
    if (rng > 0) area <- area / rng else area <- NA_real_
  }
  
  return(area)
}

# ------------------------------------------------------------------
# Plot Normality ROC Curves with AUCs in legend 
# ------------------------------------------------------------------
plot_norm_roc_curve <- function(FPR, TPR, tests_to_plot = rownames(FPR), 
                                alpha = NULL, title = NULL, dist_name = NULL, 
                                title_cex = 1, legend_digits = 3) {
  # add alternative dist to title
  if (is.null(title)) {
    suffix <- if (!is.null(dist_name)) paste(" -", dist_name) else ""
    title <- paste("ROC Curves for Different Normality Tests", suffix)
  }
  # define plot characters
  colors <- seq_along(tests_to_plot)
  plot_chars <- seq_along(tests_to_plot) + 14
  
  # create empty canvas
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate (FPR)", 
       ylab = "True Positive Rate (TPR)",
       main = title, cex.main = title_cex)
  # container for auc values
  aucs <- numeric(length(tests_to_plot))
  
  # add curves + AUC values
  for (i in seq_along(tests_to_plot)) {
    test <- tests_to_plot[i]
    x <- as.numeric(FPR[test, ])
    y <- as.numeric(TPR[test, ])
    lines(x, y, col = colors[i], lwd = 2)
    points(x, y, col = colors[i], pch = plot_chars[i], cex = 0.5)
    # compute aucs
    aucs[i] <- compute_auc(x, y, ensure_endpoints = TRUE, normalize = FALSE)
  }
  # add the 45 degrees line
  abline(0, 1, lty = 2, col = "gray")
  # create legend
  legend_labels <- sprintf("%s (AUC = %.*f)", tests_to_plot, legend_digits, aucs)
  legend("bottomright",legend = legend_labels,col = colors, 
         pch = plot_chars, lwd = 2,title = "Normality Tests", cex = 0.75)
  
  invisible(setNames(aucs, tests_to_plot))
}

# -----------------------------------------------------------------------------
# Function to compute p-values for normality test & downstream tests
# -----------------------------------------------------------------------------
generate_pval<- function(Nsim, n, effect_size,  alpha_norm = 0.05,  test = "SW", dist = "Normal", ...) {
  # Initialize p-values storage 
  pval_ds_test1_H0 <- pval_ds_test2_H0 <- numeric(Nsim)
  pval_ds_test1_H1 <- pval_ds_test2_H1 <- numeric(Nsim)
  norm_pvals_H0 <- norm_pvals_H1 <- vector("list", Nsim)  
  
  # set progress bar
  pb <- txtProgressBar(min = 0, max = Nsim, style = 3)
  
  for (i in 1:Nsim) {
    # Under H0: get parameters & generate data
    paras_H0 <- get_parameters(n, dist = dist, effect_size = 0) 
    data_H0 <- do.call(gen_data, paras_H0)
    
    # Under H1: get parameters & generate data 
    paras_H1 <- get_parameters(n, dist = dist, effect_size = effect_size)  
    data_H1 <- do.call(gen_data, paras_H1)
    
    # Get normality test objects
    normality_test_object_H0 <- fn_to_get_norm_obj(data_H0)
    normality_test_object_H1 <- fn_to_get_norm_obj(data_H1)
    
    # perform normality test
    normality_test_H0 <- normality_test(normality_test_object_H0, test = test, alpha = alpha_norm)
    normality_test_H1 <- normality_test(normality_test_object_H1, test = test, alpha = alpha_norm)
    
    # Store normality p-values
    norm_pvals_H0[[i]] <- normality_test_H0$p_values
    norm_pvals_H1[[i]] <- normality_test_H1$p_values
    
    # Get ds test p-values under H0
    pval_ds_test1_H0[i] <- fn_for_ds_test_1(data_H0)$p.value
    pval_ds_test2_H0[i] <- fn_for_ds_test_2(data_H0)$p.value
    
    # Get ds test p-values under H1
    pval_ds_test1_H1[i] <- fn_for_ds_test_1(data_H1)$p.value
    pval_ds_test2_H1[i] <- fn_for_ds_test_2(data_H1)$p.value
    
    setTxtProgressBar(pb, i)
  }
  close(pb)
  # return results
  return(list(
    pval_ds_test1_H0 = pval_ds_test1_H0,
    pval_ds_test2_H0 = pval_ds_test2_H0,
    pval_ds_test1_H1 = pval_ds_test1_H1,
    pval_ds_test2_H1 = pval_ds_test2_H1,
    norm_pvals_H0 = norm_pvals_H0,
    norm_pvals_H1 = norm_pvals_H1
  ))
}

# --------------------------------------------------------------------
#  Function to calculate power/Type I errors for each test method
# --------------------------------------------------------------------
perform_ds_analysis <- function(Nsim, n, distributions = c("exponential", "normal"), 
                             effect_size, test, alpha_pretest, test_alpha) {
  # define storage list
  ds_test_results <- list()
  error_ds_test <- list()
  power_ds_test <- list()
  
  # set progress bar
  pb_dist <- txtProgressBar(min = 0, max = length(distributions), style = 3)
  
  for(dist_idx in seq_along(distributions)) {
    dist <- distributions[dist_idx]
    
    cat("Processing distribution:", dist, "\n")
    
    # Store all p-values from generate_pval for each distribution
    ds_test_results[[dist]] <- generate_pval(Nsim, n, effect_size = effect_size, test = test, dist = dist)
    
    # Calculate Type I error rates (under H0) for test 1 & test 2
    error_ds_test[[dist]] <- list(
      error_ds_test1 = mean(ds_test_results[[dist]]$pval_ds_test1_H0 < test_alpha),
      error_ds_test2 = mean(ds_test_results[[dist]]$pval_ds_test2_H0 < test_alpha),
      
      # storage for adaptive test for error
      error_adaptive_test = numeric(length(alpha_pretest)))
    
    # Calculate Power (under H1) for for test 1 & test 2
    power_ds_test[[dist]] <- list(
      power_ds_test1 = mean(ds_test_results[[dist]]$pval_ds_test1_H1 < test_alpha),
      power_ds_test2 = mean(ds_test_results[[dist]]$pval_ds_test2_H1 < test_alpha),
      
      # storage for adaptive test for power
      power_adaptive_test = numeric(length(alpha_pretest)))
    
    # progress bar by alpha
    pb_alpha <- txtProgressBar(min = 0, max = length(alpha_pretest), style = 3)
    
    # perform adaptive test: compute decision rule for each alpha pretest level
    for(j in seq_along(alpha_pretest)) {
      alpha <- alpha_pretest[j]
      
      # For Type I error (H0): use test 1 if all norm_pvals_H0 > alpha, else test 2
      use_ds_test1_H0 <- sapply(ds_test_results[[dist]]$norm_pvals_H0, function(x) all(x > alpha))
      # adaptive test p-value
      adaptive_pvals_H0 <- ifelse(use_ds_test1_H0, ds_test_results[[dist]]$pval_ds_test1_H0, ds_test_results[[dist]]$pval_ds_test2_H0)
      # calculate adaptive type I error
      error_ds_test[[dist]]$error_adaptive_test[j] <- mean(adaptive_pvals_H0 < test_alpha)
      
      # For Power (H1): use test 1 if all norm_pvals_H1 > alpha, else test 2
      use_ds_test1_H1 <- sapply(ds_test_results[[dist]]$norm_pvals_H1, function(x) all(x > alpha))
      adaptive_pvals_H1 <- ifelse(use_ds_test1_H1, ds_test_results[[dist]]$pval_ds_test1_H1, ds_test_results[[dist]]$pval_ds_test2_H1)
      power_ds_test[[dist]]$power_adaptive_test[j] <- mean(adaptive_pvals_H1 < test_alpha)
      
      setTxtProgressBar(pb_alpha, j)
    }
    close(pb_alpha)
    setTxtProgressBar(pb_dist, dist_idx)
  }
  close(pb_dist)
  
  return(list(
    error_ds_test = error_ds_test,
    power_ds_test = power_ds_test,
    all_pvalues     = ds_test_results   
  ))
}

# -----------------------------------------------------------------------------#
# Function to compute Normality test trade-off metrics: EPL EPG, & ETIE        #
# -----------------------------------------------------------------------------#
compute_roc_metrics <- function(error_ds_test, power_ds_test, test_alpha) {
  # Non-normal case (index 1)
  power_ds_test_non_normal   <- power_ds_test[[1]]$power_ds_test1
  power_adaptive_ds_test_non_normal <- power_ds_test[[1]]$power_adaptive_test
  error_adaptive_ds_test_non_normal <- error_ds_test[[1]]$error_adaptive_test
  Expected_power_gain <- power_adaptive_ds_test_non_normal - power_ds_test_non_normal
  Expected_type1_error_inflation_non_normal <- error_adaptive_ds_test_non_normal - test_alpha
  
  # Normal case (index 2)
  power_ds_test_normal   <- power_ds_test[[2]]$power_ds_test1
  power_adaptive_ds_test_normal <- power_ds_test[[2]]$power_adaptive_test
  adaptive_error_normal <- error_ds_test[[2]]$error_adaptive_test
  Expected_power_loss      <- power_ds_test_normal - power_adaptive_ds_test_normal
  Expected_type1_error_inflation_normal <- adaptive_error_normal - test_alpha
  
  # Benchmarks (test_1 vs test_2)
  power_gain <- power_ds_test[[1]]$power_ds_test2 - power_ds_test[[1]]$power_ds_test1
  power_loss <- power_ds_test[[2]]$power_ds_test1 - power_ds_test[[2]]$power_ds_test2
  
  # store results in a list
  list(
    Expected_power_loss      = as.numeric(Expected_power_loss),
    Expected_power_gain      = as.numeric(Expected_power_gain),
    Expected_type1_error_inflation_normal = as.numeric(Expected_type1_error_inflation_normal),
    Expected_type1_error_inflation_non_normal = as.numeric(Expected_type1_error_inflation_non_normal),
    power_gain               = as.numeric(power_gain),
    power_loss               = as.numeric(power_loss)
  )
}

# =============================================================================
# Power and Type I Error Trade-off Analysis with Positive Tolerance Constraint
# 
# Features:
# - Tolerance applies ONLY to positive Type I error inflation
# - Negative inflation (conservative tests) are allowed
# - Excludes points where both distributions show negative inflation
# - Marks selected optimal point with clamped dot
# =============================================================================

plot_power_error_tradeoff <- function(
    distributions = c("exponential", "normal"),
    alpha_pretest,
    metrics,
    outer_title = "Power and Type I error trade-off for pre-testing for normality",
    tol_pos = 0,    # Tolerance for positive Type I error inflation only
    text_size_main = 0.75, 
    text_size_labels = 0.75, 
    point_size = 0.8, 
    line_padding = 0.001
) {
  # Save original graphical parameters and restore on exit
  original_par <- par(no.readonly = TRUE)
  on.exit(par(original_par))
  
  # Helper function to plot points clamped to plot boundaries
  plot_clamped_point <- function(x, y, color = "black", size = 1) {
    # Get current plot boundaries
    # par("usr") returns: c(x_min, x_max, y_min, y_max)
    plot_region <- par("usr")
    
    # Calculate tiny padding to stay inside borders
    # This prevents points from touching the very edge
    x_padding <- 1e-9 * (plot_region[2] - plot_region[1])  # Very small % of x-range
    y_padding <- 1e-9 * (plot_region[4] - plot_region[3])  # Very small % of y-range
    
    # Clamp x-coordinate: force it between left and right edges
    # Take x, but don't let it go below left edge or above right edge
    clamped_x <- min(max(x, plot_region[1] + x_padding), plot_region[2] - x_padding)
    
    
    # Clamp y-coordinate: force it between bottom and top edges  
    # Take y, but don't let it go below bottom edge or above top edge
    clamped_y <- min(max(y, plot_region[3] + y_padding), plot_region[4] - y_padding)
    
    # Plot the clamped point
    points(clamped_x, clamped_y, pch = 19, cex = size, col = color)
  }
  
  # Extract metrics from input data
  power_loss_normal <- as.numeric(metrics$Expected_power_loss)        # Normal: test1 - adaptive
  power_gain_non_normal <- as.numeric(metrics$Expected_power_gain)    # Non-normal: adaptive - test1
  type1_inflation_normal <- as.numeric(metrics$Expected_type1_error_inflation_normal)    # Normal: adaptive - alpha
  type1_inflation_non_normal <- as.numeric(metrics$Expected_type1_error_inflation_non_normal)  # Non-normal: adaptive - alpha
  
  # Filter valid data points (finite values)
  # Removes any missing, infinite, or NaN values that could break calculations
  valid_indices <- is.finite(alpha_pretest) & is.finite(power_loss_normal) & is.finite(power_gain_non_normal) & is.finite(type1_inflation_normal) & is.finite(type1_inflation_non_normal)
  
  alpha_values <- alpha_pretest[valid_indices]
  power_loss <- power_loss_normal[valid_indices]
  power_gain <- power_gain_non_normal[valid_indices]
  type1_infl_normal <- type1_inflation_normal[valid_indices]
  type1_infl_non_normal <- type1_inflation_non_normal[valid_indices]
  
  # Calculate Type I error inflation metrics
  numerical_precision <- .Machine$double.eps^0.5
  
  # Positive inflation (bad - we want to control this)
  positive_infl_normal <- pmax(type1_infl_normal, 0)           # Set negative values to 0
  positive_infl_non_normal <- pmax(type1_infl_non_normal, 0)   # Keep only positive inflation
  worst_case_positive_inflation <- pmax(positive_infl_normal, positive_infl_non_normal)
  # Identifies the maximum positive Type I error inflation across both distributions.
  # Negative inflation (conservative - acceptable but not ideal)
  both_negative <- (type1_infl_normal < 0 & type1_infl_non_normal < 0)
  not_both_negative <- !both_negative
  
  negative_infl_normal <- pmax(-type1_infl_normal, 0)
  negative_infl_non_normal <- pmax(-type1_infl_non_normal, 0)
  worst_case_negative_inflation <- pmax(negative_infl_normal, negative_infl_non_normal)
  
  # Define feasible set: within positive tolerance & not both negative
  feasible_indices <- which(
    worst_case_positive_inflation <= tol_pos + numerical_precision & not_both_negative
  )
  
  # Selection algorithm: Type I error control first, then power considerations
  select_optimal_alpha <- function(candidate_indices) {
    # Step 1: Minimize worst-case positive inflation (primary constraint)
    # Among all candidate alpha values, find the smallest possible positive Type I error inflation
    # Keep only those alpha values that achieve this minimum (within numerical tolerance)
    positive_inflation_values <- worst_case_positive_inflation[candidate_indices]
    min_positive_inflation <- min(positive_inflation_values, na.rm = TRUE)
    # Find alpha values with minimum worst-case positive inflation
    step1_indices <- candidate_indices[
      abs(positive_inflation_values - min_positive_inflation) <= numerical_precision * (1 + min_positive_inflation)
    ]
    
    # Step 2: Minimize worst-case negative inflation (avoid over-conservatism)
    # From the alpha values that minimize positive inflation, now minimize how conservative they are
    negative_inflation_values <- worst_case_negative_inflation[step1_indices]
    min_negative_inflation <- min(negative_inflation_values, na.rm = TRUE)
    # Among those, find alpha values with minimum negative inflation
    step2_indices <- step1_indices[
      abs(negative_inflation_values - min_negative_inflation) <= numerical_precision * (1 + min_negative_inflation)
    ]
    
    # Step 3-4: Optimize power trade-offs
    has_negative_power_loss <- any(power_loss[step2_indices] < -numerical_precision)
    
    if (has_negative_power_loss) {
      # When adaptive test is better for normal distribution
      # Maximize power gain for non-normal distributions
      gain_values <- power_gain[step2_indices]
      max_gain <- max(gain_values, na.rm = TRUE)
      step3_indices <- step2_indices[
        abs(gain_values - max_gain) <= numerical_precision * (1 + abs(max_gain))
      ]
      
      # Minimize power loss (more negative is better)
      loss_values <- power_loss[step3_indices]
      min_loss <- min(loss_values, na.rm = TRUE)
      step4_indices <- step3_indices[
        abs(loss_values - min_loss) <= numerical_precision * (1 + abs(min_loss))
      ]
    } else {
      # Standard case: minimize power loss, then maximize power gain
      loss_values <- power_loss[step2_indices]
      min_loss <- min(loss_values, na.rm = TRUE)
      step3_indices <- step2_indices[
        abs(loss_values - min_loss) <= numerical_precision * (1 + abs(min_loss))
      ]
      
      gain_values <- power_gain[step3_indices]
      max_gain <- max(gain_values, na.rm = TRUE)
      step4_indices <- step3_indices[
        abs(gain_values - max_gain) <= numerical_precision * (1 + abs(max_gain))
      ]
    }
    
    # Final tie-breaker: smallest alpha value
    if (length(step4_indices) > 1) {
      # If multiple alpha values still tied, pick the smallest alpha
      step4_indices[which.min(alpha_values[step4_indices])]
    } else {
      step4_indices
    }
  }
  # ------------------------------------------------------------------------
  # Apply selection algorithm with fallback strategies
  # If we have alpha values that satisfy both constraints
  # Use the full selection algorithm on these ideal candidates
  # when no alpha values meet the ideal constraints, use fallback strategies. 
  # -------------------------------------------------------------------------
  if (length(feasible_indices) > 0) {
    optimal_index <- select_optimal_alpha(feasible_indices)
    selection_note <- sprintf(" (within tolerance & not-both-negative: max +inflation ≤ %.5g)", tol_pos)
  } else {
    # Fallback 1: Allow not-both-negative, minimize positive inflation globally
    not_both_neg_indices <- which(not_both_negative)
    if (length(not_both_neg_indices) > 0) {
      min_positive_infl <- min(worst_case_positive_inflation[not_both_neg_indices], na.rm = TRUE)
      near_optimal_indices <- not_both_neg_indices[
        abs(worst_case_positive_inflation[not_both_neg_indices] - min_positive_infl) <= numerical_precision * (1 + min_positive_infl)
      ]
      optimal_index <- select_optimal_alpha(near_optimal_indices)
      selection_note <- sprintf(" (fallback not-both-negative: min max +inflation = %.5g)", min_positive_infl)
    } else {
      # Fallback 2: All points have both negative inflation - pick least negative
      min_negative_infl <- min(worst_case_negative_inflation, na.rm = TRUE)
      near_optimal_indices <- which(
        abs(worst_case_negative_inflation - min_negative_infl) <= numerical_precision * (1 + min_negative_infl)
      )
      optimal_index <- select_optimal_alpha(near_optimal_indices)
      selection_note <- sprintf(" (fallback both-negative: least negative = %.5g)", min_negative_infl)
    }
  }
  
  # Extract optimal values
  opt_alpha <- alpha_values[optimal_index]
  optimal_power_loss <- power_loss[optimal_index]
  optimal_power_gain <- power_gain[optimal_index]
  optimal_infl_normal <- type1_infl_normal[optimal_index]
  optimal_infl_non_normal <- type1_infl_non_normal[optimal_index]
  
  ## layout for plots 
  par(mfrow = c(2, 2), oma  = c(0.8, 0.2, 2.5, 0.2), mar  = c(4, 3.5, 1.5, 0.8), mgp = c(2.0, 0.6, 0)) 
  
  # plot function 
  create_plot <- function(x, y, color, ylab, main, optimal_index) {
    optimal_y <- y[optimal_index]  
    y_range <- range(y, na.rm = TRUE) + c(-1, 1) * line_padding
    
    plot(x, y, type = "l", col = color, lwd = 2, 
         ylim = y_range, ylab = ylab, xlab = expression(alpha[pre]), 
         cex.lab = text_size_labels, main = main, cex.main = text_size_main)
    abline(v = x[optimal_index], lty = 3, col = "grey40")
    plot_clamped_point(x[optimal_index], optimal_y, color, point_size)
  }
  
  # Expected Power Loss
  create_plot(alpha_values, power_loss, "blue", "Expected Power Loss", sprintf("Power Loss: %s", distributions[2]), optimal_index)
  # Expected Power Gain
  create_plot(alpha_values, power_gain, "red", "Expected Power Gain", sprintf("Power Gain: %s", distributions[1]), optimal_index)
  # Type I Inflation, normal
  create_plot(alpha_values, type1_infl_normal, "orange", "Type I Inflation", sprintf("Type I Inflation: %s", distributions[2]), optimal_index)
  abline(h = 0, lty = 2, col = "grey60")
  # Type I Inflation, non-normal
  create_plot(alpha_values, type1_infl_non_normal, "green4", "Type I Inflation", sprintf("Type I Inflation: %s", distributions[1]), optimal_index)
  
  abline(h = 0, lty = 2, col = "grey60")
  # Add titles 
  mtext(outer_title, side = 3, outer = TRUE, line = 1.2, cex = 1, font = 2)
  mtext(sprintf("Selected alpha = %.5f | Gain = %.5f | Loss = %.5f | Inflation (Normal, Non-normal) = (%.5f, %.5f) | tol = %.5g", opt_alpha, optimal_power_gain, optimal_power_loss, optimal_infl_normal, optimal_infl_non_normal, tol_pos),side = 3, outer = TRUE, line = 0.2, cex = 0.6)
  
  # Return results invisibly so they can be captured if needed
  # This allows users to access the selection details programmatically
  invisible(list(
    alpha_star = opt_alpha,
    power_loss = optimal_power_loss,
    power_gain = optimal_power_gain,
    inflation_normal = optimal_infl_normal,
    inflation_non_normal = optimal_infl_non_normal,
    selection_note = selection_note,
    feasible = length(feasible_indices) > 0,
    worst_case_positive_inflation = max(positive_infl_normal[optimal_index], positive_infl_non_normal[optimal_index]),
    worst_case_negative_inflation = max(negative_infl_normal[optimal_index], negative_infl_non_normal[optimal_index])
  ))
}

# -----------------------------------------------------------------------
#             Function for Power vs Type I error ROC curve 
# -----------------------------------------------------------------------
power_vs_error_roc_data <- function(N, n, distributions = c("exponential", "normal"), test = "SW", effect_size, alpha_pretest = 0.05, sig_levels) {
  
  # empty data frame for results
  roc_data <- data.frame()
  
  pb <- txtProgressBar(min = 0, max = length(distributions) * length(sig_levels), style = 3)
  counter <- 0
  
  for (dist in distributions) {
    cat("Generating ROC data for:", dist, "\n")
    # generate p-values
    res <- generate_pval(N, n, effect_size, test = test, dist = dist)
    
    # calculate power and type I error for each significance level
    for (alpha in sig_levels) {
      # H0 adaptive decision
      use_ds_test1_H0 <- sapply(res$norm_pvals_H0, function(v) all(as.numeric(v) > alpha_pretest))
      # adaptive test p-value
      adaptive_pvals_H0 <- ifelse(use_ds_test1_H0, res$pval_ds_test1_H0, res$pval_ds_test2_H0)
      # calculate Type I error for test 1 test 2 and adaptive test
      type1_error_test1 <- mean(res$pval_ds_test1_H0 <= alpha)
      type1_error_test2 <- mean(res$pval_ds_test2_H0 <= alpha)
      type1_error_adaptive <- mean(adaptive_pvals_H0 <= alpha)
      
      # H1 adaptive decision
      use_ds_test1_H1 <- sapply(res$norm_pvals_H1, function(v) all(as.numeric(v) > alpha_pretest))
      # adaptive test p-value
      adaptive_pvals_H1 <- ifelse(use_ds_test1_H1, res$pval_ds_test1_H1, res$pval_ds_test2_H1)
      # calculate power for test 1 test 2 and adaptive test
      power_test1 <- mean(res$pval_ds_test1_H1 <= alpha)
      power_test2 <- mean(res$pval_ds_test2_H1 <= alpha)
      power_adaptive <- mean(adaptive_pvals_H1 <= alpha)
      
      # SIMPLIFIED: Create consistent column names
      roc_data <- rbind(roc_data,
        data.frame(Distribution = dist, Method = "ds_test 1", Alpha = alpha, Power = power_test1, TypeI_error = type1_error_test1),
        data.frame(Distribution = dist, Method = "ds_test 2", Alpha = alpha, Power = power_test2, TypeI_error = type1_error_test2),
        data.frame(Distribution = dist, Method = "adaptive test", Alpha = alpha, Power = power_adaptive,TypeI_error = type1_error_adaptive)
      )
      counter <- counter + 1 
      setTxtProgressBar(pb, counter)
    }
  }
  close(pb)
  return(roc_data)
}
# ------------------------------------------------------------------------------
#             power vs type I error rate ROC like curve function
# ------------------------------------------------------------------------------
power_vs_error_ROC_curve <- function(roc_data, alpha_star = 0.05, zoom_xlim = 0.10) {
  # Set up layout
  n_dist <- length(unique(roc_data$Distribution))
  par(mfrow = c(2, n_dist), mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
  
  # Colors for methods
  colors <- c("ds_test 1" = "red", "ds_test 2" = "blue", "adaptive test" = "green")
  
  for (dist in unique(roc_data$Distribution)) {
    data_sub <- subset(roc_data, Distribution == dist)
    
    # Full view
    plot(NA, xlim = c(0, 1), ylim = c(0, 1), xlab = "Type I Error", ylab = "Power",main = paste(dist, "- Full View"))
    
    # Add reference lines
    abline(0, 1, col = "gray80")
    abline(v = alpha_star, col = "red", lty = 2)
    
    # Plot each method
    for (method in unique(data_sub$Method)) {
      method_data <- subset(data_sub, Method == method)
      lines(method_data$TypeI_error, method_data$Power, col = colors[method], lwd = 2)
      
      # Mark alpha_star point
      idx <- which.min(abs(method_data$Alpha - alpha_star))
      points(method_data$TypeI_error[idx], method_data$Power[idx], col = colors[method], pch = 16)
    }
    
    # Zoom view
    plot(NA, xlim = c(0, zoom_xlim), ylim = c(0, 1), xlab = "Type I Error", ylab = "Power",main = paste(dist, "- Zoom View"))
    
    abline(0, 1, col = "gray80")
    abline(v = alpha_star, col = "red", lty = 2)
    
    for (method in unique(data_sub$Method)) {
      method_data <- subset(data_sub, Method == method)
      lines(method_data$TypeI_error, method_data$Power, col = colors[method], lwd = 2)
      
      # Mark alpha_star point
      idx <- which.min(abs(method_data$Alpha - alpha_star))
      points(method_data$TypeI_error[idx], method_data$Power[idx], col = colors[method], pch = 16)
    }
  }
  
  # Add overall title
  mtext("Power vs Type I Error ROC Curves", outer = TRUE, cex = 1.2)
  
  # Add legend in a separate plot
  par(fig = c(0, 1, 0, 0.3), new = TRUE, mar = c(0, 0, 0, 0))
  plot(0, 0, type = "n", axes = FALSE, xlab = "", ylab = "")
  legend("center", legend = names(colors), col = colors, lwd = 2, pch = 16, bty = "n", horiz = TRUE, cex = 0.9)
}


# ============================================================================
# ================ Step 4: Main Downstream Test Function =====================
# --------- Adaptive test function: Decide test type(test 1 or test 2) -------
ds_test_function <- function(gen_data = gen_data, get_parameters = get_parameters, 
                             fn_to_get_norm_obj = fn_to_get_norm_obj, 
                             fn_for_norm_test = normality_test, 
                             fn_for_ds_test_1 = fn_for_ds_test_1, 
                             fn_for_ds_test_2 = fn_for_ds_test_2, 
                             paras = NULL, alpha = 0.05, norm_test_method = "SW", 
                             alpha_pre = 0.05, data = NULL, ...) {
  # generate data
  data <- if (!is.null(data)) data else if (!is.null(paras)) do.call(gen_data, paras) else gen_data()
  
  # get normality test object
  normality_test_object <- fn_to_get_norm_obj(data)
  
  # perform normality test
  normality_test_pval   <- fn_for_norm_test(data = normality_test_object, test = norm_test_method, alpha = alpha_pre)
  
  # decide downstream test based on normality test results
  if (isTRUE(normality_test_pval$normality_satisfied)) {
    ds_test <- fn_for_ds_test_1(data)
  } else {
    ds_test <- fn_for_ds_test_2(data, ...)
  }
  return(ds_test$p.value)
}

# --------- perform test 1, test 2 & adaptive test -------------
# --------- perform test 1, test 2 & adaptive test -------------
perform_ds_func <- function(sample_sizes = sample_sizes, 
                            distributions = c("exponential", "normal"), 
                            N = 1e3, alpha = 0.05, gen_data = gen_data, 
                            get_parameters = get_parameters, 
                            fn_to_get_norm_obj = fn_to_get_norm_obj, 
                            fn_for_norm_test = normality_test, 
                            fn_for_ds_test_1 = fn_for_ds_test_1, 
                            fn_for_ds_test_2 = fn_for_ds_test_2, 
                            norm_test_method = "SW", 
                            ds_test_methods = ds_test_methods, 
                            effect_size = effect_size, 
                            alpha_pre = alpha_pre, ...) {
  
  # pick correct ds test method
  ds_test_methods <- match.arg(ds_test_methods, several.ok = TRUE)
  # storage for all ds tests results
  results_by_dist <- list()
  
  # perform ds tests for each distribution
  for (dist in distributions) {
    cat("\n=== Running simulation for distribution:", dist, "===\n")
    
    # Initialize storage for results
    results_power <- list()
    results_type1 <- list()
    pval_storage <- list(H0 = list(), H1 = list())
    timing <- list()
    
    # Initialize storage for each method
    for (method in ds_test_methods) {
      results_type1[[method]] <- numeric(length(sample_sizes))
      results_power[[method]] <- numeric(length(sample_sizes))
      names(results_type1[[method]]) <- names(results_power[[method]]) <- paste0("n=", sample_sizes)
      pval_storage$H0[[method]] <- vector("list", length(sample_sizes))
      pval_storage$H1[[method]] <- vector("list", length(sample_sizes))
    }
    
    # set progress bar for sample sizes
    pb <- txtProgressBar(min = 0, max = length(sample_sizes), style = 3)
    
    for (i in seq_along(sample_sizes)) {
      n <- sample_sizes[i]
      
      # Track rejections for each method
      rej_H0_vec <- rej_H1_vec <- numeric(length(ds_test_methods))
      names(rej_H0_vec) <- names(rej_H1_vec) <- ds_test_methods
      
      # Initialize matrices to store p-values for all simulations and all methods
      p_H0_mat <- matrix(NA, nrow = N, ncol = length(ds_test_methods))
      p_H1_mat <- matrix(NA, nrow = N, ncol = length(ds_test_methods))
      colnames(p_H0_mat) <- colnames(p_H1_mat) <- ds_test_methods
      
      # Start timing for this sample size
      t0 <- Sys.time()
      
      for (sim in seq_len(N)) {
        # Generate data ONCE per simulation (common for all tests)
        paras_H0 <- get_parameters(n, effect_size = 0.0, dist = dist)
        paras_H1 <- get_parameters(n, effect_size = effect_size, dist = dist)
        dat_H0 <- do.call(gen_data, paras_H0)
        dat_H1 <- do.call(gen_data, paras_H1)
        
        # Run ALL tests on the SAME data
        for (method_idx in seq_along(ds_test_methods)) {
          method <- ds_test_methods[method_idx]
          
          if (method == "test_1") {
            p0 <- fn_for_ds_test_1(dat_H0)$p.value
            p1 <- fn_for_ds_test_1(dat_H1)$p.value
          } else if (method == "test_2") {
            p0 <- fn_for_ds_test_2(dat_H0)$p.value
            p1 <- fn_for_ds_test_2(dat_H1)$p.value
          } else {  # adaptive test
            p0 <- ds_test_function(
              gen_data = gen_data, 
              get_parameters = get_parameters,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_norm_test = fn_for_norm_test,
              fn_for_ds_test_1 = fn_for_ds_test_1,
              fn_for_ds_test_2 = fn_for_ds_test_2,
              paras = paras_H0, 
              alpha = alpha,
              norm_test_method = norm_test_method, 
              alpha_pre = alpha_pre,
              data = dat_H0,  # Pass the pre-generated data
              ...
            )
            p1 <- ds_test_function(
              gen_data = gen_data, 
              get_parameters = get_parameters,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_norm_test = fn_for_norm_test,
              fn_for_ds_test_1 = fn_for_ds_test_1,
              fn_for_ds_test_2 = fn_for_ds_test_2,
              paras = paras_H1, 
              alpha = alpha,
              norm_test_method = norm_test_method, 
              alpha_pre = alpha_pre,
              data = dat_H1,  # Pass the pre-generated data
              ...
            )
          }
          
          # Store p-values
          p_H0_mat[sim, method_idx] <- p0
          p_H1_mat[sim, method_idx] <- p1
          
          # Count rejections
          if (p0 < alpha) rej_H0_vec[method_idx] <- rej_H0_vec[method_idx] + 1
          if (p1 < alpha) rej_H1_vec[method_idx] <- rej_H1_vec[method_idx] + 1
        }
      }
      
      # Store results for each method
      for (method_idx in seq_along(ds_test_methods)) {
        method <- ds_test_methods[method_idx]
        results_type1[[method]][i] <- rej_H0_vec[method_idx] / N
        results_power[[method]][i] <- rej_H1_vec[method_idx] / N
        pval_storage$H0[[method]][[i]] <- p_H0_mat[, method_idx]
        pval_storage$H1[[method]][[i]] <- p_H1_mat[, method_idx]
      }
      
      # Record timing
      timing[[paste0("n=", n)]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      
      setTxtProgressBar(pb, i)
    }
    close(pb)
    
    # Collect results by distribution
    results_by_dist[[dist]] <- list(
      power = results_power, 
      type1 = results_type1,
      pvalues = pval_storage,
      timing = timing
    )
    
    cat("Distribution", dist, "completed\n")
  }
  
  return(results_by_dist)
}


# plot power versus type i error
plot_power_type1_results <- function(
    combined_power, combined_type1,
    ds_test_methods = c("test_1", "test_2", "adaptive"),
    distributions   = c("Non-normal", "Normal"),
    sizes, 
    alpha, 
    alpha_pre,
    filename = NULL, width = 7, height = 6
){
  plot_colors   <- c("red", "blue", "green")
  method_shapes <- c(19, 17, 15)
  method_names  <- ds_test_methods
  
  opened <- FALSE
  if (!is.null(filename)) { pdf(filename, width = width, height = height); opened <- TRUE }
  
  layout(matrix(c(1,2, 3,4, 5,5), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.25))
  
  op <- par(mar = c(2.6, 3.2, 1.2, 0.8), oma = c(0, 0, 1.2, 0),
            mgp = c(1.6, 0.45, 0), tcl = -0.2, xaxs = "r", yaxs = "r",
            cex.axis = 0.75, cex.lab = 0.75)
  on.exit({ par(op); if (opened) dev.off() }, add = TRUE)
  
  cex.main <- 0.7
  xr <- range(sizes, finite = TRUE); xpad <- 0.02 * diff(xr)
  xlim_pad <- c(xr[1] - xpad, xr[2] + xpad)
  
  for (dist in distributions) {
    dist_data <- subset(combined_power, Distribution == dist)
    plot(NA, xlim = xlim_pad, ylim = c(0, 1.01),
         xlab = "Sample Size", ylab = "Power",
         main = paste("Power -", dist), cex.main = cex.main, font.main = 1)
    for (j in seq_along(method_names)) {
      method <- method_names[j]
      lines(dist_data$n, dist_data[[method]], type = "b",
            col = plot_colors[j], pch = method_shapes[j], lwd = 2)
    }
    abline(h = 0.8, col = "gray50", lty = 2, lwd = 1.2)
  }
  
  for (dist in distributions) {
    dist_data <- subset(combined_type1, Distribution == dist)
    ymax <- max(dist_data[, method_names], na.rm = TRUE)
    ypad <- 0.02 * ymax
    plot(NA, xlim = xlim_pad, ylim = c(0, ymax + max(ypad, 1e-6)),
         xlab = "Sample Size", ylab = "P(Type I Error)",
         main = paste("P(Type I Error) -", dist), cex.main = cex.main, font.main = 1)
    for (j in seq_along(method_names)) {
      method <- method_names[j]
      lines(dist_data$n, dist_data[[method]], type = "b",
            col = plot_colors[j], pch = method_shapes[j], lwd = 2)
    }
    abline(h = alpha, col = "gray50", lty = 3, lwd = 1.2)
  }
  
  par(mar = c(0.2, 0.2, 0.2, 0.2)); plot.new()
  legend("center", legend = method_names, col = plot_colors, title = "Methods",
         pch = method_shapes, lwd = 2, horiz = TRUE, bty = "n", cex = 0.75)
  
  mtext(substitute(paste("Comparison of Test Methods Across Sample Sizes (pre-test ", alpha, " = ", v, ")"),
                   list(v = sprintf("%.3f", alpha_pre))),
        side = 3, outer = TRUE, line = 0.35, cex = 0.75)
}


perform_ds_power_by_effect <- function(
    fixed_n            = 10,
    effect_sizes       = c(0.1, 0.2, 0.3, 0.5, 1.0),
    distributions      = c("Non-normal", "Normal"),
    Nsim               = Nsim,
    alpha              = alpha,
    gen_data           = gen_data,
    get_parameters     = get_parameters,
    fn_to_get_norm_obj = fn_to_get_norm_obj,
    fn_for_norm_test   = normality_test,
    fn_for_ds_test_1   = fn_for_ds_test_1,
    fn_for_ds_test_2   = fn_for_ds_test_2,
    norm_test_method   = norm_test_method,
    ds_test_methods    = ds_test_methods,
    alpha_pre          = alpha_pre,
    ...
) {
  ds_test_methods <- match.arg(ds_test_methods, several.ok = TRUE)
  results_by_dist <- list()
  
  for (dist in distributions) {
    cat("\n=== Running POWER simulation for distribution:", dist,
        "with fixed n =", fixed_n, "===\n")
    
    results_power <- list(); timing <- list()
    pval_storage  <- list(H1 = list())
    
    for (method in ds_test_methods) {
      cat("Running method:", method, "\n")
      pow_vec <- numeric(length(effect_sizes))
      names(pow_vec) <- paste0("d=", effect_sizes)
      pval_storage$H1[[method]] <- vector("list", length(effect_sizes))
      pb <- txtProgressBar(min = 0, max = length(effect_sizes), style = 3)
      t0 <- Sys.time()
      
      for (i in seq_along(effect_sizes)) {
        d <- effect_sizes[i]; rej_H1 <- 0L; p_H1 <- numeric(Nsim)
        for (sim in seq_len(Nsim)) {
          paras_H1 <- get_parameters(fixed_n, effect_size = d, dist = dist)
          if (method == "test_1") {
            p1 <- fn_for_ds_test_1(do.call(gen_data, paras_H1))$p.value
          } else if (method == "test_2") {
            p1 <- fn_for_ds_test_2(do.call(gen_data, paras_H1))$p.value
          } else {
            p1 <- ds_test_function(
              gen_data = gen_data, get_parameters = get_parameters,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_norm_test   = fn_for_norm_test,
              fn_for_ds_test_1   = fn_for_ds_test_1,
              fn_for_ds_test_2   = fn_for_ds_test_2,
              paras = paras_H1, alpha = alpha,
              norm_test_method = norm_test_method, alpha_pre = alpha_pre, ...
            )
          }
          p_H1[sim] <- p1; if (p1 < alpha) rej_H1 <- rej_H1 + 1L
        }
        pow_vec[i] <- rej_H1 / Nsim
        pval_storage$H1[[method]][[i]] <- p_H1
        setTxtProgressBar(pb, i)
      }
      close(pb)
      
      timing[[method]]        <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      results_power[[method]] <- pow_vec
      cat("Method", method, "completed in", round(timing[[method]], 2), "seconds\n")
    }
    
    results_by_dist[[dist]] <- list(power = results_power,
                                    pvalues = pval_storage,
                                    timing  = timing)
  }
  results_by_dist
}

plot_power_by_effect_size <- function(
    combined_power, distributions = c("Non-normal", "Normal"),
    ds_test_methods = c("test_1", "test_2", "adaptive"),
    sample_sizes, effect_sizes, alpha_pre = 0.05,
    filename = NULL, width = 7, height = 4
) {
  plot_colors   <- c("red", "blue", "green")
  method_shapes <- c(19, 17, 15)
  method_names  <- ds_test_methods
  
  # pdf(filename, width = width, height = height)
  # op <- par(no.readonly = TRUE); on.exit({ par(op); dev.off() }, add = TRUE)
  
  # Reduce the bottom outer margin significantly
  par(mfrow = c(1, 2), mar = c(3.5, 3.5, 2, 1), oma = c(1, 0, 1, 0))
  panel_idx <- 0
  
  for (dist in distributions) for (n_now in sample_sizes) {
    panel_idx <- panel_idx + 1
    dist_data <- subset(combined_power, Distribution == dist & n == n_now)
    if (!nrow(dist_data)) {
      plot(NA, xlim = range(effect_sizes), ylim = c(0, 1),
           xlab = "", ylab = "",
           main = sprintf("%s (n = %d)", dist, n_now),
           cex.main = 0.7, cex.axis = 0.6)
      title(xlab = "Effect Size", ylab = "Power", line = 2, cex.lab = 0.65)
      next
    }
    plot(NA, xlim = range(effect_sizes, na.rm = TRUE), ylim = c(0, 1),
         xlab = "", ylab = "",
         main = sprintf("%s (n = %d)", dist, n_now),
         cex.main = 0.7, cex.axis = 0.6)
    title(xlab = "Effect Size", ylab = "Power", line = 2, cex.lab = 0.65)
    
    for (j in seq_along(method_names)) {
      method <- method_names[j]
      if (!method %in% names(dist_data)) next
      ok <- is.finite(dist_data$d) & is.finite(dist_data[[method]])
      if (!any(ok)) next
      x <- dist_data$d[ok]; y <- dist_data[[method]][ok]
      ord <- order(x); x <- x[ord]; y <- y[ord]
      lines(x, y, type = "l", col = plot_colors[j], lwd = 2)
      y_at_05 <- approx(x = x, y = y, xout = 0.5, rule = 2, ties = mean)$y
      points(0.5, y_at_05, pch = method_shapes[j], col = plot_colors[j], cex = 1.0)
    }
    abline(v = 0.5, col = "gray50", lty = 2, lwd = 1.2)
  }
  
  # Add main title closer to plots
  mtext(substitute(paste("Power vs Effect Size by Distribution | pre-test ", alpha, " = ", v), list(v = sprintf("%.3f", alpha_pre))),
        side = 3, outer = TRUE, line = 0, cex = 0.75, font = 2)
  
  # Add single legend at bottom center 
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", 
         legend = method_names, 
         title = "Method",
         col = plot_colors, 
         lwd = 2, 
         pch = method_shapes, 
         bty = "n",
         cex = 0.75,
         xpd = TRUE,
         horiz = TRUE,
         inset = c(0, -0.01))
}
# ===================== Core DS functions & plotting =====================

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
    alpha_pre = 0.05,
    ...
) {
  data <- if (!is.null(paras)) do.call(gen_data, paras) else gen_data()
  normality_test_object <- fn_to_get_norm_obj(data)
  normality_test_pval   <- fn_for_norm_test(data = normality_test_object,
                                            test = norm_test_method, alpha = alpha_pre)
  if (isTRUE(normality_test_pval$normality_satisfied)) {
    ds_test <- fn_for_ds_test_1(data)
  } else {
    ds_test <- fn_for_ds_test_2(data, ...)
  }
  ds_test$p.value
}

perform_ds_func <- function(
    sample_sizes       = sample_sizes,
    distributions      = c("Non-normal", "Normal"),
    Nsim               = 1e3,
    alpha              = 0.05,
    gen_data           = gen_data,
    get_parameters     = get_parameters,
    fn_to_get_norm_obj = fn_to_get_norm_obj,
    fn_for_norm_test   = normality_test,
    fn_for_ds_test_1   = fn_for_ds_test_1,
    fn_for_ds_test_2   = fn_for_ds_test_2,
    norm_test_method   = "SW",
    ds_test_methods    = ds_test_methods,
    effect_size        = effect_size,
    alpha_pre          = alpha_pre,
    ...
) {
  ds_test_methods <- match.arg(ds_test_methods, several.ok = TRUE)
  results_by_dist <- list()
  
  for (dist in distributions) {
    cat("\n=== Running simulation for distribution:", dist, "===\n")
    results_power <- list(); results_type1 <- list()
    timing <- list(); pval_storage <- list(H0 = list(), H1 = list())
    
    for (method in ds_test_methods) {
      cat("Running method:", method, "\n")
      pow_vec <- typ1_vec <- numeric(length(sample_sizes))
      names(pow_vec) <- names(typ1_vec) <- paste0("n=", sample_sizes)
      
      pval_storage$H0[[method]] <- vector("list", length(sample_sizes))
      pval_storage$H1[[method]] <- vector("list", length(sample_sizes))
      pb <- txtProgressBar(min = 0, max = length(sample_sizes), style = 3)
      t0 <- Sys.time()
      
      for (i in seq_along(sample_sizes)) {
        n <- sample_sizes[i]; rej_H0 <- rej_H1 <- 0L
        p_H0 <- p_H1 <- numeric(Nsim)
        
        for (sim in seq_len(Nsim)) {
          paras_H0 <- get_parameters(n, effect_size = 0.0, dist = dist)
          paras_H1 <- get_parameters(n, effect_size = effect_size, dist = dist)
          
          dat_H0 <- do.call(gen_data, paras_H0)
          dat_H1 <- do.call(gen_data, paras_H1)
          
          if (method == "test_1") {
            p0 <- fn_for_ds_test_1(dat_H0)$p.value
            p1 <- fn_for_ds_test_1(dat_H1)$p.value
          } else if (method == "test_2") {
            p0 <- fn_for_ds_test_2(dat_H0)$p.value
            p1 <- fn_for_ds_test_2(dat_H1)$p.value
          } else {
            # Get the alpha_pre for current sample size
            current_alpha_pre <- if(length(alpha_pre) == 1) alpha_pre else alpha_pre[i]
            
            p0 <- ds_test_function(
              gen_data = gen_data, get_parameters = get_parameters,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_norm_test   = fn_for_norm_test,
              fn_for_ds_test_1   = fn_for_ds_test_1,
              fn_for_ds_test_2   = fn_for_ds_test_2,
              paras = paras_H0, alpha = alpha,
              norm_test_method = norm_test_method, alpha_pre = current_alpha_pre, ...
            )
            p1 <- ds_test_function(
              gen_data = gen_data, get_parameters = get_parameters,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_norm_test   = fn_for_norm_test,
              fn_for_ds_test_1   = fn_for_ds_test_1,
              fn_for_ds_test_2   = fn_for_ds_test_2,
              paras = paras_H1, alpha = alpha,
              norm_test_method = norm_test_method, alpha_pre = current_alpha_pre, ...
            )
          }
          p_H0[sim] <- p0; p_H1[sim] <- p1
          if (p0 < alpha) rej_H0 <- rej_H0 + 1
          if (p1 < alpha) rej_H1 <- rej_H1 + 1
        }
        
        typ1_vec[i] <- rej_H0 / Nsim
        pow_vec[i]  <- rej_H1 / Nsim
        pval_storage$H0[[method]][[i]] <- p_H0
        pval_storage$H1[[method]][[i]] <- p_H1
        setTxtProgressBar(pb, i)
      }
      close(pb)
      
      timing[[method]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      results_type1[[method]] <- typ1_vec
      results_power[[method]] <- pow_vec
      cat("Method", method, "completed in", round(timing[[method]], 2), "seconds\n")
    }
    
    results_by_dist[[dist]] <- list(
      power = results_power, type1 = results_type1,
      pvalues = pval_storage, timing = timing
    )
  }
  results_by_dist
}

plot_power_type1_results <- function(
    combined_power, combined_type1,
    ds_test_methods = c("test_1", "test_2", "adaptive"),
    distributions   = c("Non-normal", "Normal"),
    sizes, alpha, alpha_pre,
    filename = NULL, width = 7, height = 6
){
  plot_colors   <- c("red", "blue", "green")
  method_shapes <- c(19, 17, 15)
  method_names  <- ds_test_methods
  
  opened <- FALSE
  if (!is.null(filename)) { pdf(filename, width = width, height = height); opened <- TRUE }
  
  layout(matrix(c(1,2, 3,4, 5,5), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.25))
  
  op <- par(mar = c(2.6, 3.2, 1.2, 0.8), oma = c(0, 0, 1.2, 0),
            mgp = c(1.6, 0.45, 0), tcl = -0.2, xaxs = "r", yaxs = "r",
            cex.axis = 0.75, cex.lab = 0.75)
  on.exit({ par(op); if (opened) dev.off() }, add = TRUE)
  
  cex.main <- 0.7
  xr <- range(sizes, finite = TRUE); xpad <- 0.02 * diff(xr)
  xlim_pad <- c(xr[1] - xpad, xr[2] + xpad)
  
  for (dist in distributions) {
    dist_data <- subset(combined_power, Distribution == dist)
    plot(NA, xlim = xlim_pad, ylim = c(0, 1.01),
         xlab = "Sample Size", ylab = "Power",
         main = paste("Power -", dist), cex.main = cex.main, font.main = 1)
    for (j in seq_along(method_names)) {
      method <- method_names[j]
      lines(dist_data$n, dist_data[[method]], type = "b",
            col = plot_colors[j], pch = method_shapes[j], lwd = 2)
    }
    abline(h = 0.8, col = "gray50", lty = 2, lwd = 1.2)
  }
  
  for (dist in distributions) {
    dist_data <- subset(combined_type1, Distribution == dist)
    ymax <- max(dist_data[, method_names], na.rm = TRUE)
    ypad <- 0.02 * ymax
    plot(NA, xlim = xlim_pad, ylim = c(0, ymax + max(ypad, 1e-6)),
         xlab = "Sample Size", ylab = "P(Type I Error)",
         main = paste("P(Type I Error) -", dist), cex.main = cex.main, font.main = 1)
    for (j in seq_along(method_names)) {
      method <- method_names[j]
      lines(dist_data$n, dist_data[[method]], type = "b",
            col = plot_colors[j], pch = method_shapes[j], lwd = 2)
    }
    abline(h = alpha, col = "gray50", lty = 3, lwd = 1.2)
  }
  
  par(mar = c(0.2, 0.2, 0.2, 0.2)); plot.new()
  legend("center", legend = method_names, col = plot_colors, title = "Methods",
         pch = method_shapes, lwd = 2, horiz = TRUE, bty = "n", cex = 0.75)
  
  # Create appropriate title based on alpha_pre type
  if(length(alpha_pre) == 1) {
    main_title <- sprintf("Comparison of Test Methods Across Sample Sizes (pre-test alpha = %.3f)", alpha_pre)
  } else {
    alpha_pre_str <- paste(sprintf("%.3f", alpha_pre), collapse = ", ")
    main_title <- sprintf("Comparison of Test Methods | alpha_pre by n: [%s]", alpha_pre_str)
  }
  mtext(main_title, side = 3, outer = TRUE, line = 0.35, cex = 0.75)
}

perform_ds_power_by_effect <- function(
    fixed_n            = 10,
    effect_sizes       = c(0.1, 0.2, 0.3, 0.5, 1.0),
    distributions      = c("Non-normal", "Normal"),
    Nsim               = Nsim,
    alpha              = alpha,
    gen_data           = gen_data,
    get_parameters     = get_parameters,
    fn_to_get_norm_obj = fn_to_get_norm_obj,
    fn_for_norm_test   = normality_test,
    fn_for_ds_test_1   = fn_for_ds_test_1,
    fn_for_ds_test_2   = fn_for_ds_test_2,
    norm_test_method   = norm_test_method,
    ds_test_methods    = ds_test_methods,
    alpha_pre          = alpha_pre,
    ...
) {
  ds_test_methods <- match.arg(ds_test_methods, several.ok = TRUE)
  results_by_dist <- list()
  
  for (dist in distributions) {
    cat("\n=== Running POWER simulation for distribution:", dist,
        "with fixed n =", fixed_n, "===\n")
    
    results_power <- list(); timing <- list()
    pval_storage  <- list(H1 = list())
    
    for (method in ds_test_methods) {
      cat("Running method:", method, "\n")
      pow_vec <- numeric(length(effect_sizes))
      names(pow_vec) <- paste0("d=", effect_sizes)
      pval_storage$H1[[method]] <- vector("list", length(effect_sizes))
      pb <- txtProgressBar(min = 0, max = length(effect_sizes), style = 3)
      t0 <- Sys.time()
      
      for (i in seq_along(effect_sizes)) {
        d <- effect_sizes[i]; rej_H1 <- 0L; p_H1 <- numeric(Nsim)
        for (sim in seq_len(Nsim)) {
          paras_H1 <- get_parameters(fixed_n, effect_size = d, dist = dist)
          if (method == "test_1") {
            p1 <- fn_for_ds_test_1(do.call(gen_data, paras_H1))$p.value
          } else if (method == "test_2") {
            p1 <- fn_for_ds_test_2(do.call(gen_data, paras_H1))$p.value
          } else {
            p1 <- ds_test_function(
              gen_data = gen_data, get_parameters = get_parameters,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_norm_test   = fn_for_norm_test,
              fn_for_ds_test_1   = fn_for_ds_test_1,
              fn_for_ds_test_2   = fn_for_ds_test_2,
              paras = paras_H1, alpha = alpha,
              norm_test_method = norm_test_method, alpha_pre = alpha_pre, ...
            )
          }
          p_H1[sim] <- p1; if (p1 < alpha) rej_H1 <- rej_H1 + 1L
        }
        pow_vec[i] <- rej_H1 / Nsim
        pval_storage$H1[[method]][[i]] <- p_H1
        setTxtProgressBar(pb, i)
      }
      close(pb)
      
      timing[[method]]        <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      results_power[[method]] <- pow_vec
      cat("Method", method, "completed in", round(timing[[method]], 2), "seconds\n")
    }
    
    results_by_dist[[dist]] <- list(power = results_power,
                                    pvalues = pval_storage,
                                    timing  = timing)
  }
  results_by_dist
}

plot_power_by_effect_size <- function(
    combined_power, distributions = c("Non-normal", "Normal"),
    ds_test_methods = c("test_1", "test_2", "adaptive"),
    sample_sizes, effect_sizes, alpha_pre = 0.05,
    filename = NULL, width = 7, height = 4
) {
  plot_colors   <- c("red", "blue", "green")
  method_shapes <- c(19, 17, 15)
  method_names  <- ds_test_methods
  
  # pdf(filename, width = width, height = height)
  # op <- par(no.readonly = TRUE); on.exit({ par(op); dev.off() }, add = TRUE)
  
  # Reduce the bottom outer margin significantly
  par(mfrow = c(1, 2), mar = c(3.5, 3.5, 2, 1), oma = c(1, 0, 1, 0))
  panel_idx <- 0
  
  for (dist in distributions) for (n_now in sample_sizes) {
    panel_idx <- panel_idx + 1
    dist_data <- subset(combined_power, Distribution == dist & n == n_now)
    if (!nrow(dist_data)) {
      plot(NA, xlim = range(effect_sizes), ylim = c(0, 1),
           xlab = "", ylab = "",
           main = sprintf("%s (n = %d)", dist, n_now),
           cex.main = 0.7, cex.axis = 0.6)
      title(xlab = "Effect Size", ylab = "Power", line = 2, cex.lab = 0.65)
      next
    }
    plot(NA, xlim = range(effect_sizes, na.rm = TRUE), ylim = c(0, 1),
         xlab = "", ylab = "",
         main = sprintf("%s (n = %d)", dist, n_now),
         cex.main = 0.7, cex.axis = 0.6)
    title(xlab = "Effect Size", ylab = "Power", line = 2, cex.lab = 0.65)
    
    for (j in seq_along(method_names)) {
      method <- method_names[j]
      if (!method %in% names(dist_data)) next
      ok <- is.finite(dist_data$d) & is.finite(dist_data[[method]])
      if (!any(ok)) next
      x <- dist_data$d[ok]; y <- dist_data[[method]][ok]
      ord <- order(x); x <- x[ord]; y <- y[ord]
      lines(x, y, type = "l", col = plot_colors[j], lwd = 2)
      y_at_05 <- approx(x = x, y = y, xout = 0.5, rule = 2, ties = mean)$y
      points(0.5, y_at_05, pch = method_shapes[j], col = plot_colors[j], cex = 1.0)
    }
    abline(v = 0.5, col = "gray50", lty = 2, lwd = 1.2)
  }
  
  # Add main title closer to plots
  mtext(sprintf("Power vs Effect Size by Distribution | pre-test = %.3f", alpha_pre),
        outer = TRUE, side = 3, line = 0, cex = 0.75, font = 2)
  
  # Add single legend at bottom center 
  par(fig = c(0, 1, 0, 1), oma = c(0, 0, 0, 0), mar = c(0, 0, 0, 0), new = TRUE)
  plot(0, 0, type = "n", bty = "n", xaxt = "n", yaxt = "n")
  legend("bottom", 
         legend = method_names, 
         title = "Method",
         col = plot_colors, 
         lwd = 2, 
         pch = method_shapes, 
         bty = "n",
         cex = 0.75,
         xpd = TRUE,
         horiz = TRUE,
         inset = c(0, -0.01))
}