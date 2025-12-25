setwd("~/Desktop/OSU/Research/Pretest-Simulation/User_framework/two_sample/wilcoxon/Wilcoxon_v2")
# Add Dependencies Check
required_packages <- c("LaplacesDemon", "VGAM")
check_packages <- function() {
  missing <- setdiff(required_packages, rownames(installed.packages()))
  if (length(missing) > 0) {
    stop("Please install required packages: ", paste(missing, collapse = ", "))
  }
}

# ----------------------------------------------------
#   Generate standardized data from various distributions
# ----------------------------------------------------
generate_data <- function(n, dist, par = NULL, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  # Input validation
  if (!is.numeric(n) || n <= 0) stop("n must be a positive integer")
  dist <- tolower(dist)
  
  if (dist == "normal") {
    if (is.null(par)) par <- c(0, 1)
    x <- rnorm(n, mean = par[1], sd = par[2])
    
  } else if (dist == "chi_square") {
    if (is.null(par)) par <- 3
    x <- (rchisq(n, df = par) - par) / sqrt(2 * par)
    
  } else if (dist == "gamma") {
    if (is.null(par)) par <- c(3, 0.1)
    mean_g <- par[1] / par[2]
    sd_g <- sqrt(par[1] / par[2]^2)
    x <- (rgamma(n, shape = par[1], rate = par[2]) - mean_g) / sd_g
    
  } else if (dist == "exponential") {
    if (is.null(par)) par <- 1
    mean_e <- 1 / par
    sd_e <- 1 / par
    x <- (rexp(n, rate = par) - mean_e) / sd_e
    
  } else if (dist == "t") {
    if (is.null(par)) par <- 3
    if (par <= 1) stop("Degrees of freedom must be > 1")
    x <- rt(n, df = par) / sqrt(par / (par - 2))
    
  } else if (dist == "uniform") {
    if (is.null(par)) par <- c(0, 1)
    mean_u <- (par[1] + par[2]) / 2
    sd_u <- sqrt((par[2] - par[1])^2 / 12)
    x <- (runif(n, min = par[1], max = par[2]) - mean_u) / sd_u
    
  } else if (dist == "laplace") {
    if (is.null(par)) par <- c(2, 7)  
    x <- LaplacesDemon::rlaplace(n, location = par[1], scale = par[2])
    x <- (x - par[1]) / sqrt(2 * par[2]^2)
    
  } else if (dist == "weibull") {
    if (is.null(par)) par <- c(1, 2)
    mean_w <- par[2] * gamma(1 + 1 / par[1])
    var_w <- par[2]^2 * (gamma(1 + 2 / par[1]) - gamma(1 + 1 / par[1])^2)
    x <- (rweibull(n, shape = par[1], scale = par[2]) - mean_w) / sqrt(var_w)
    
  } else if (dist == "lognormal") {
    if (is.null(par)) par <- c(0, 1)
    mean_ln <- exp(par[1] + par[2]^2 / 2)
    var_ln <- (exp(par[2]^2) - 1) * exp(2 * par[1] + par[2]^2)
    x <- (rlnorm(n, meanlog = par[1], sdlog = par[2]) - mean_ln) / sqrt(var_ln)
    
  }  else if (dist == "cauchy") {
    # Note: Cauchy has no finite variance!
    if (is.null(par)) par <- c(0, 1)
    x <- rcauchy(n, location = par[1], scale = par[2])
    # Cannot standardize - maybe median center and IQR scale?
    x <- (x - median(x)) / IQR(x)
    
  } else if (dist == "poisson") {
    if (is.null(par)) par <- 5
    x <- (rpois(n, lambda = par) - par) / sqrt(par)  # Poisson standardization
    
  }else if (dist == "contaminated") {
    # par = c(prob_good, mean, sd_good, sd_bad)
    if (is.null(par)) par <- c(0.75, 0, 1, 5)
    br <- rbinom(n, size = 1, prob = par[1])
    sd_br <- ifelse(br == 1, par[3], par[4])
    x_raw <- rnorm(n, mean = par[2], sd = sd_br)
    # Variance: weighted avg of variances (mean is constant)
    var_c <- par[1] * par[3]^2 + (1 - par[1]) * par[4]^2
    x <- (x_raw - par[2]) / sqrt(var_c)
    
  } else if (dist == "pareto") {
    if (is.null(par)) par <- c(1, 3)
    if (par[2] <= 2) stop("Shape parameter for Pareto must be > 2 for variance to exist")
    mean_p <- par[1] * par[2] / (par[2] - 1)
    var_p <- (par[1]^2 * par[2]) / ((par[2] - 1)^2 * (par[2] - 2))
    x <- (VGAM::rpareto(n, scale = par[1], shape = par[2]) - mean_p) / sqrt(var_p)
    
  }else if (dist == "beta") {
    # par = c(alpha, beta), both > 0
    if (is.null(par)) par <- c(2, 5)
    if (any(par <= 0)) stop("Beta parameters must be > 0")
    mean_b <- par[1] / (par[1] + par[2])
    var_b  <- (par[1] * par[2]) / ((par[1] + par[2])^2 * (par[1] + par[2] + 1))
    x <- (rbeta(n, shape1 = par[1], shape2 = par[2]) - mean_b) / sqrt(var_b)
    
  }else if (dist == "logistic") {
    # par = c(location, scale>0)
    if (is.null(par)) par <- c(0, 1)
    if (par[2] <= 0) stop("Logistic scale must be > 0")
    mean_l <- par[1]
    var_l  <- (pi^2 / 3) * par[2]^2
    x <- (rlogis(n, location = par[1], scale = par[2]) - mean_l) / sqrt(var_l)
    
  } else {
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
    se <- sqrt(6/length(x)) # for large n
    z <- s / se 
    p <- 2 * (1 - pnorm(abs(z)))
    output <- list(statistic = z, p.value = p, method = "Skewness z-test")
    
  } else if(test == "KURT"){  # Kurtosis test (z-test)
    k <- moments::kurtosis(x)
    se <- sqrt(24/length(x)) # for large n
    z <- (k - 3) / se
    p <- 2 * (1 - pnorm(abs(z)))
    output <- list(statistic = z, p.value = p, method = "Kurtosis z-test")
    
  } else {
    stop("Unknown test name. Choose from: KS, SW, JB, DAP, ANS, AD, AD2, SF, CVM, SKEW, KURT.")
  }
  
  return(output)
}

# ------------------------------------
# Define simple calculation functions 
# -------------------------------------

# Generic test statistic function
TwoSample_test_statistic <- function(x, y, type = "welch") {
  if(type == "welch") {
    # Welch t-statistic
    return((mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))
  } else if(type == "pooled") {
    # Pooled t-statistic
    n1 <- length(x); n2 <- length(y)
    s_pooled <- sqrt(((n1-1)*var(x) + (n2-1)*var(y))/(n1+n2-2))
    return((mean(x) - mean(y))/(s_pooled*sqrt(1/n1 + 1/n2)))
  } else if(type == "median") {
    # Difference of medians
    return(median(x) - median(y))
  } else if(type == "mood") {
    # Mood's test statistic (median test)
    combined_median <- median(c(x, y))
    m1 <- sum(x > combined_median)
    m2 <- sum(y > combined_median)
    n1 <- length(x); n2 <- length(y)
    return((m1/n1 - m2/n2)/sqrt((m1+m2)*(n1+n2-m1-m2)/(n1*n2*(n1+n2))))
  }
}


# compute the area under the curve 
compute_area <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / (max(x) - min(x))
}


# Two Sample Permutation 
two_sample_permutation_test <- function(x, y, B, test_stat = TwoSample_test_statistic) {
  observed <- test_stat(x, y)
  combined <- c(x, y)
  n_x <- length(x)
  n_total <- length(combined)
  
  perm_stats <- replicate(B, {
    perm_indices <- sample(n_total)
    x_star <- combined[perm_indices[1:n_x]]
    y_star <- combined[perm_indices[(n_x+1):n_total]]
    test_stat(x_star, y_star)
  })
  
  # Two-sided p-value
  p_value <- mean(abs(perm_stats) >= abs(observed))
  return(list(p.value = p_value, statistic = observed))
}

# ------------ framework functions -----------------------------------------
# Data generation function
gen_data <- function(
    n1 = n,         
    n2 = n,         
    mean1 = 0.0,       
    effect_size = 0.0,     
    sd1 = 1,         
    sd2 = 1,         
    dist = "chi_square",
    par = NULL,
    ...
) {
  # generate both samples
  group1 <- mean1 + sd1 * generate_data(n1, dist, par = par)
  group2 <- effect_size + sd2 * generate_data(n2, dist, par = par)
  # return dataframe
  return(data.frame(
    group = factor(rep(c("x", "y"), c(n1, n2))),
    value = c(group1, group2)
  ))
}

# parameters function
get_parameters <- function(n, ...) {
  defaults <- list(
    n1 = n,
    n2 = n,
    mean1 = 0.0,
    effect_size = 0.0,
    sd1 = 1,
    sd2 = 1,
    dist = "chi_square",
    par = NULL
  )
  
  modifyList(defaults, list(...))
}

# get normality test object function
fn_to_get_norm_obj <- function(data) {
  return(data)
}

# downstream test 1(t-test) function
fn_for_ds_test_1 <- function(data) {
  x_data <- data$value[data$group == "x"]
  y_data <- data$value[data$group == "y"]
  test_result <- t.test(x = x_data, y = y_data, var.equal = FALSE)
  return(list(p.value = test_result$p.value))
}

# downstream test 2(Mann-Whitney U Test) function
fn_for_ds_test_2 <- function(data) {
  x_data <- data$value[data$group == "x"]
  y_data <- data$value[data$group == "y"]
  test_result <- wilcox.test(x_data, y_data)
  return(list(p.value = test_result$p.value))
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

# -----------------------------------------------------------------------------#
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
plot_norm_roc_curve <- function(FPR, TPR, 
                                tests_to_plot = rownames(FPR), 
                                alpha = NULL, 
                                title = NULL, 
                                dist_name = NULL, 
                                title_cex = 1, 
                                legend_digits = 3) {
  # add alternative dist to title
  if (is.null(title)) {
    suffix <- if (!is.null(dist_name)) paste(" -", dist_name) else ""
    title <- paste("ROC Curves for Different Normality Tests", suffix)
  }
  # define plot characters
  colors <- seq_along(tests_to_plot)
  plot_chars <- seq_along(tests_to_plot) + 14
  
  # create empty canvas
  plot(0, 0, type = "n", xlim = c(0, 1), 
       ylim = c(0, 1),
       xlab = "False Positive Rate (FPR)", 
       ylab = "True Positive Rate (TPR)",
       main = title, cex.main = title_cex)
  
  aucs <- numeric(length(tests_to_plot))
  
  # add curves + AUC values
  for (i in seq_along(tests_to_plot)) {
    test <- tests_to_plot[i]
    x <- as.numeric(FPR[test, ])
    y <- as.numeric(TPR[test, ])
    lines(x, y, col = colors[i], lwd = 2)
    points(x, y, col = colors[i], 
           pch = plot_chars[i], cex = 0.5)
    aucs[i] <- compute_auc(x, y, ensure_endpoints = TRUE, normalize = FALSE)
  }
  
  abline(0, 1, lty = 2, col = "gray")
  # create legend
  legend_labels <- sprintf("%s (AUC = %.*f)", tests_to_plot, legend_digits, aucs)
  legend("bottomright",
         legend = legend_labels,
         col = colors, 
         pch = plot_chars, 
         lwd = 2,
         title = "Normality Tests", 
         cex = 0.75)
  
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
perform_analysis <- function(Nsim, n, distributions = c("Non-normal", "Normal"), 
                             effect_size, test, alpha_pretest, test_alpha) {
  # define storage list
  ds_test_results <- list()
  error_ds_test <- list()
  power_ds_test <- list()
  
  pb_dist <- txtProgressBar(min = 0, max = length(distributions), style = 3)
  
  for(dist_idx in seq_along(distributions)) {
    dist <- distributions[dist_idx]
    cat("Processing distribution:", dist, "\n")
    
    # Store all p-values from generate_pval for each distribution
    ds_test_results[[dist]] <- generate_pval(Nsim, n, effect_size = effect_size, test = "SW", dist = dist)
    
    # Calculate Type I error rates (under H0)
    error_ds_test[[dist]] <- list(
      error_ds_test1 = mean(ds_test_results[[dist]]$pval_ds_test1_H0 < test_alpha),
      error_ds_test2 = mean(ds_test_results[[dist]]$pval_ds_test2_H0 < test_alpha),
      
      # storage for adaptive test for error
      error_adaptive_test = numeric(length(alpha_pretest)))
    
    # Calculate Power (under H1)
    power_ds_test[[dist]] <- list(
      power_ds_test1 = mean(ds_test_results[[dist]]$pval_ds_test1_H1 < test_alpha),
      power_ds_test2 = mean(ds_test_results[[dist]]$pval_ds_test2_H1 < test_alpha),
      
      # storage for adaptive test for power
      power_adaptive_test = numeric(length(alpha_pretest)))
    
    pb_alpha <- txtProgressBar(min = 0, max = length(alpha_pretest), style = 3)
    
    # perform adaptive test: compute decisions for each alpha pretest level
    for(j in seq_along(alpha_pretest)) {
      alpha <- alpha_pretest[j]
      
      # For Type I error (H0): use test 1 if all norm_pvals_H0 > alpha, else test 2
      use_ds_test1_H0 <- sapply(ds_test_results[[dist]]$norm_pvals_H0, function(x) all(x > alpha))
      # adaptive test p-value
      adaptive_pvals_H0 <- ifelse(use_ds_test1_H0, ds_test_results[[dist]]$pval_ds_test1_H0, ds_test_results[[dist]]$pval_ds_test2_H0)
      # calculate adaptive type I error
      error_ds_test[[dist]]$error_adaptive_test[j] <- mean(adaptive_pvals_H0 < test_alpha)
      
      # For Power (H1): use test 1 if all norm_pvals_H0 > alpha, else test 2
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
    distributions = c("Non-normal", "Normal"),
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
    worst_case_positive_inflation <= tol_pos + numerical_precision & 
      not_both_negative
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
  
  ## ----- 2×2 layout for plots -----
  par(mfrow = c(2, 2), oma  = c(0.8, 0.2, 2.5, 0.2), mar  = c(4, 3.5, 1.5, 0.8), mgp = c(2.0, 0.6, 0)) 
  
  # helper plot function 
  create_plot <- function(x, y, color, ylab, main, optimal_index) {
    optimal_y <- y[optimal_index]  
    y_range <- range(y, na.rm = TRUE) + c(-1, 1) * line_padding
    
    plot(x, y, type = "l", col = color, lwd = 2, ylim = y_range, ylab = ylab, xlab = expression(alpha[pre]), cex.lab = text_size_labels, main = main, cex.main = text_size_main)
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
power_vs_error_roc_data <- function(N, n, distributions = c("Non-normal", "Normal"), effect_size, alpha_pretest = 0.05, sig_levels) {
  
  # empty data frame for results
  roc_data <- data.frame()
  
  pb <- txtProgressBar(min = 0, max = length(distributions) * length(sig_levels), style = 3)
  counter <- 0
  
  for (dist in distributions) {
    cat("Generating ROC data for:", dist, "\n")
    # generate p-values
    res <- generate_pval(N, n, effect_size, test = "SW", dist = dist)
    
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
      
      # organize results in data frames
      roc_data <- rbind(
        roc_data,
        data.frame(Distribution = dist, Method = "ds_test 1", Alpha = alpha, Power = power_test1, type1_error = type1_error_test1),
        data.frame(Distribution = dist, Method = "ds_test 2", Alpha = alpha, Power = power_test2, type1_error = type1_error_test2),
        data.frame(Distribution = dist, Method = "adaptive test", Alpha = alpha, Power = power_adaptive,type1_error = type1_error_adaptive)
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
power_vs_error_ROC_curve <- function(roc_data, file = NULL, width = 8, height = 6.8, alpha_star = 0.05, zoom_xlim = 0.10, subtitle = "Power vs Type I error: full view (first row) and zoom view (second row)"){
  
  # Define plot parameters
  dists <- unique(roc_data$Distribution)
  methods <- unique(roc_data$Method)
  colors <- c("red", "blue", "green")
  point_char <- 16  # Solid circle
  
  # Layout
  layout(matrix(1:6, nrow = 3, byrow = TRUE), heights = c(1, 1, 0.30))
  
  # Set graphical parameters
  op <- par(oma = c(0.5, 0.7, 1.5, 0.6), mgp = c(1.6, 0.50, 0), tcl = -0.25, xaxs = "r", yaxs = "r", cex.axis = 0.9, cex.lab = 0.9, cex.main = 0.75)
  on.exit(par(op))
  
  # Plot panel function
  plot_panel <- function(dist, xlim, bottom_mar = 2.0) {
    par(mar = c(bottom_mar, 3.2, 1.0, 0.8))
    
    # Create empty plot
    plot(NA, xlim = xlim, ylim = c(0, 1), xlab = "P(Type I Error)", ylab = "Power", main = paste0("ROC-like Curve (", dist, ")"), cex.lab = 0.65)
    
    # Add lines for each method
    for (i in seq_along(methods)) {
      data_sub <- subset(roc_data, Distribution == dist & Method == methods[i])
      if (nrow(data_sub) == 0) next
      
      lines(data_sub$type1_error, data_sub$Power, col = colors[i], lwd = 2)
      
      # Mark alpha_star point
      idx <- which.min(abs(data_sub$Alpha - alpha_star))
      if (length(idx) == 1 && is.finite(idx)) {
        points(data_sub$type1_error[idx], data_sub$Power[idx], pch = point_char, col = colors[i])
      }
    }
    
    # Reference lines
    abline(a = 0, b = 1, col = "grey80")
    abline(v = alpha_star, col = "red", lty = 2)  
  }
  
  # Create plots
  plot_panel(dists[1], c(0, 1), 2.0)           # Full view - dist 1
  plot_panel(dists[2], c(0, 1), 2.0)           # Full view - dist 2
  plot_panel(dists[1], c(0, zoom_xlim), 2.5)   # Zoom view - dist 1
  plot_panel(dists[2], c(0, zoom_xlim), 2.5)   # Zoom view - dist 2
  
  # Legend function
  add_legend <- function(dist) {
    par(mar = c(0, 0.3, 0.6, 0.3))
    plot.new()
    
    legend_labels <- sapply(seq_along(methods), function(i) {
      data_sub <- subset(roc_data, Distribution == dist & Method == methods[i])
      if (nrow(data_sub) == 0) return(methods[i])
      # add alpha_star point
      idx <- which.min(abs(data_sub$Alpha - alpha_star))
      x <- data_sub$type1_error[idx]
      y <- data_sub$Power[idx]
      sprintf("%s (%.3f, %.3f)", methods[i], x, y)
    })
    
    legend("top", legend = legend_labels, title = "Methods", col = colors, lwd = 2, pch = point_char, cex = 0.75, bty = "n")
  }
  
  # Add legends
  add_legend(dists[1])
  add_legend(dists[2])
  
  # Main title
  mtext(subtitle, side = 3, outer = TRUE, line = 0.2, cex = 0.65)
}
# --------------------------------------------------------------------
# ================ Main Downstream Test Function =====================
# --------- Decide test type(test 1 or test 2) -----------------
ds_test_function <- function(gen_data = gen_data, get_parameters = get_parameters, fn_to_get_norm_obj = fn_to_get_norm_obj, fn_for_norm_test = normality_test, fn_for_ds_test_1 = fn_for_ds_test_1, fn_for_ds_test_2 = fn_for_ds_test_2, paras = NULL, alpha = 0.05, norm_test_method = "SW", alpha_pre = 0.05, ...) {
  # generate data
  data <- if (!is.null(paras)) do.call(gen_data, paras) else gen_data()
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

# perform test 1, test 2 & adaptive test
perform_ds_func <- function(sample_sizes = sample_sizes, distributions = c("Non-normal", "Normal"), N = 1e3, alpha = 0.05, gen_data = gen_data, get_parameters = get_parameters, fn_to_get_norm_obj = fn_to_get_norm_obj, fn_for_norm_test   = normality_test, fn_for_ds_test_1   = fn_for_ds_test_1, fn_for_ds_test_2 = fn_for_ds_test_2, norm_test_method = "SW", ds_test_methods = ds_test_methods, effect_size = effect_size, alpha_pre = alpha_pre,...) {
  # pick correct ds test method
  ds_test_methods <- match.arg(ds_test_methods, several.ok = TRUE)
  # storage for all ds tests results
  results_by_dist <- list()
  
  # perform ds tests for each distribution
  for (dist in distributions) {
    cat("\n=== Running simulation for distribution:", dist, "===\n")
    # storage for power and type I error
    results_power <- list()
    results_type1 <- list()
    # track runtime
    timing <- list() 
    # store all p-values
    pval_storage <- list(H0 = list(), H1 = list())
    
    # for each ds test method
    for (method in ds_test_methods) {
      cat("Running method:", method, "\n")
      # storage for each sample size
      pow_vec <- typ1_vec <- numeric(length(sample_sizes))
      names(pow_vec) <- names(typ1_vec) <- paste0("n=", sample_sizes)
      # store all p-values
      pval_storage$H0[[method]] <- vector("list", length(sample_sizes))
      pval_storage$H1[[method]] <- vector("list", length(sample_sizes))
      
      # set progress bar
      pb <- txtProgressBar(min = 0, max = length(sample_sizes), style = 3)
      # start time
      t0 <- Sys.time()
      
      for (i in seq_along(sample_sizes)) {
        n <- sample_sizes[i]
        #alpha_pre <- alpha_optimum
        
        rej_H0 <- rej_H1 <- 0L
        p_H0 <- p_H1 <- numeric(N)
        
        for (sim in seq_len(N)) {
          # get parameters
          paras_H0 <- get_parameters(n, effect_size = 0.0, dist = dist)
          paras_H1 <- get_parameters(n, effect_size = effect_size, dist = dist)
          # generate data
          dat_H0 <- do.call(gen_data, paras_H0)
          dat_H1 <- do.call(gen_data, paras_H1)
          
          if (method == "test_1") { # power & type 1 error p-values for test 1
            p0 <- fn_for_ds_test_1(dat_H0)$p.value
            p1 <- fn_for_ds_test_1(dat_H1)$p.value
            
          } else if (method == "test_2") {# power & type 1 error p-values for test 2
            p0 <- fn_for_ds_test_2(dat_H0)$p.value
            p1 <- fn_for_ds_test_2(dat_H1)$p.value
            
          } else {# type I error for adaptive test
            p0 <- ds_test_function(
              gen_data = gen_data, get_parameters = get_parameters,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_norm_test   = fn_for_norm_test,
              fn_for_ds_test_1   = fn_for_ds_test_1,
              fn_for_ds_test_2   = fn_for_ds_test_2,
              paras = paras_H0, alpha = alpha,
              norm_test_method = norm_test_method, alpha_pre = alpha_pre, ...
            )# power for adaptive test
            p1 <- ds_test_function(
              gen_data = gen_data, get_parameters = get_parameters,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_norm_test   = fn_for_norm_test,
              fn_for_ds_test_1   = fn_for_ds_test_1,
              fn_for_ds_test_2   = fn_for_ds_test_2,
              paras = paras_H1, alpha = alpha,
              norm_test_method = norm_test_method, alpha_pre = alpha_pre, ...
            )
          }# store p-values 
          p_H0[sim] <- p0
          p_H1[sim] <- p1
          if (p0 < alpha) rej_H0 <- rej_H0 + 1
          if (p1 < alpha) rej_H1 <- rej_H1 + 1
        }
        # calculate power & type I error rates
        typ1_vec[i] <- rej_H0 / N
        pow_vec[i]  <- rej_H1 / N
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
    # collect all results by distribution
    results_by_dist[[dist]] <- list(
      power = results_power, 
      type1 = results_type1,
      pvalues = pval_storage, 
      timing = timing
    )
  }
  return(results_by_dist)
}

# -----------------------------------------------------------------------------
#                     power & type I error plot function 
# -----------------------------------------------------------------------------
plot_power_type1 <- function(combined_power, combined_type1, optimal_alphas, methods = c("test_1", "test_2", "adaptive"), distributions = c("Non-normal", "Normal"), sizes, alpha, filename = NULL, width = 7, height = 6) {
  
  # Setup
  colors <- c("red", "blue", "green")
  shapes <- c(19, 17, 15)
  
  # Open PDF if file name provided
  if (!is.null(filename)) {
    pdf(filename, width = width, height = height)
    on.exit(dev.off())
  }
  
  # Plot layout: 2 plots, 2 plots, 1 legend
  layout(matrix(c(1,2,3,4,5,5), nrow = 3, byrow = TRUE), heights = c(1,1,0.25))
  op <- par(mar = c(2.6, 3.2, 1.2, 0.8), oma = c(0, 0,1.6,0), mgp = c(1.6, 0.45, 0), tcl = -0.2, cex.axis = 0.75, cex.lab = 0.75)
  on.exit(par(op), add = TRUE)
  
  # Plot function
  create_plot <- function(data, ylab, main, hline = NULL) {
    xlim <- range(sizes) + c(-1, 1) * 0.02 * diff(range(sizes))
    
    plot(NA, xlim = xlim, ylim = c(0, 1.01), xlab = "Sample Size", ylab = ylab, main = main)
    
    for (i in seq_along(methods)) {
      lines(data$n, data[[methods[i]]], type = "b", col = colors[i], pch = shapes[i], lwd = 2)
    }
    
    if (!is.null(hline)) abline(h = hline, col = "gray50", lty = 2)
  }
  # Create plots
  for (dist in distributions) {
    data_power <- subset(combined_power, Distribution == dist)
    create_plot(data_power, "Power", paste("Power -", dist), 0.75)
  }
  
  for (dist in distributions) {
    data_type1 <- subset(combined_type1, Distribution == dist)
    ymax <- max(data_type1[methods], na.rm = TRUE)
    plot(NA, xlim = range(sizes), ylim = c(0, ymax * 1.02), xlab = "Sample Size", ylab = "P(Type I Error)", main = paste("P(Type I Error) -", dist))
    
    for (i in seq_along(methods)) {
      lines(data_type1$n, data_type1[[methods[i]]], type = "b", col = colors[i], pch = shapes[i], lwd = 2)
    }
    abline(h = alpha, col = "gray50", lty = 3)
  }
  
  # Legend
  par(mar = c(0.2, 0.2, 0.2, 0.2))
  plot.new()
  legend("center", legend = methods, col = colors, pch = shapes, lwd = 2, horiz = TRUE, bty = "n", cex = 0.75, title = "Methods")
  
  # Title with optimal alphas
  mtext(sprintf("Test Methods Comparison (Optimal alpha: %s)", paste(round(optimal_alphas, 3), collapse = ", ")), side = 3, outer = TRUE, line = 0.35, cex = 0.75)
}
# ----------------------------------------------------------------------------
# calculate power for each effect size
# ----------------------------------------------------------------------------
perform_ds_power_by_effect <- function(fixed_n = 10, effect_sizes = NULL, distributions = c("Non-normal", "Normal"), N = N, alpha = alpha, gen_data = gen_data,
                                       get_parameters = get_parameters, fn_to_get_norm_obj = fn_to_get_norm_obj, fn_for_norm_test = normality_test,
                                       fn_for_ds_test_1 = fn_for_ds_test_1, fn_for_ds_test_2 = fn_for_ds_test_2, norm_test_method = norm_test_method,
                                       ds_test_methods = ds_test_methods, alpha_pre = alpha_pre, ...) {
  
  results_by_dist <- list()
  
  for (dist in distributions) {
    cat("Power simulation:", dist, "n =", fixed_n, "\n")
    # store results
    power_results <- list()
    timing <- list()
    pvalues <- list(H1 = list())
    
    for (method in ds_test_methods) {
      cat("Method:", method, "\n")
      
      power_vec <- numeric(length(effect_sizes))
      names(power_vec) <- paste0("d=", effect_sizes)
      pvalues$H1[[method]] <- vector("list", length(effect_sizes))
      
      pb <- txtProgressBar(0, length(effect_sizes), style = 3)
      start_time <- Sys.time()
      
      for (i in seq_along(effect_sizes)) {
        effect <- effect_sizes[i]
        rejections <- 0
        p_vals <- numeric(N)
        
        for (sim in 1:N) {
          paras <- get_parameters(fixed_n, effect_size = effect, dist = dist)
          data <- do.call(gen_data, paras)
          
          if (method == "test_1") { # test 1
            p_val <- fn_for_ds_test_1(data)$p.value
          } else if (method == "test_2") { # test 2
            p_val <- fn_for_ds_test_2(data)$p.value
          } else {# adaptive test
            p_val <- ds_test_function(
              gen_data = gen_data, get_parameters = get_parameters,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_norm_test = fn_for_norm_test,
              fn_for_ds_test_1 = fn_for_ds_test_1,
              fn_for_ds_test_2 = fn_for_ds_test_2,
              paras = paras, alpha = alpha,
              norm_test_method = norm_test_method, alpha_pre = alpha_pre, ...
            )
          }
          
          p_vals[sim] <- p_val
          if (p_val < alpha) rejections <- rejections + 1
        }
        # calculate power and store all p-values
        power_vec[i] <- rejections / N
        pvalues$H1[[method]][[i]] <- p_vals
        setTxtProgressBar(pb, i)
      }
      
      close(pb)
      timing[[method]] <- as.numeric(Sys.time() - start_time)
      power_results[[method]] <- power_vec
      cat("Completed in", round(timing[[method]], 2), "s\n")
    }
    
    results_by_dist[[dist]] <- list(
      power = power_results,
      pvalues = pvalues,
      timing = timing
    )
  }
  
  results_by_dist
}

# -----------------------------------------------------------------------
# plot power against effect sizes
# ----------------------------------------------------------------------
plot_power_by_effect_size <- function(power_results, 
                                      distributions = c("Non-normal", "Normal"),
                                      ds_test_methods = c("test_1", "test_2", "adaptive"),
                                      effect_sizes,
                                      alpha_pre = 0.05) {
  
  colors <- c("red", "blue", "green")
  shapes <- c(19, 17, 15)
  
  # Save and restore graphics settings
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  # Layout: 2 plots on top, legend row at bottom
  layout(
    matrix(c(1, 2,
             3, 3), nrow = 2, byrow = TRUE),
    heights = c(0.8, 0.2)
  )
  
  # Margins for the two main plots
  par(mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
  
  for (dist in distributions) {
    dist_data <- subset(power_results, Distribution == dist)
    
    plot(NA,
         xlim = range(effect_sizes),
         ylim = c(0, 1),
         xlab = "Effect Size",
         ylab = "Power",
         main = dist,
         cex.main = 0.9)
    
    for (i in seq_along(ds_test_methods)) {
      method <- ds_test_methods[i]
      if (method %in% names(dist_data)) {
        y <- dist_data[[method]]
        lines(effect_sizes, y, col = colors[i], lwd = 2)
        points(effect_sizes, y, pch = shapes[i], col = colors[i])
      }
    }
    
    abline(v = 0.5, col = "gray50", lty = 2)
  }
  
  # Overall title
  mtext(sprintf("Power vs Effect Size | pre-test = %.3f", alpha_pre),
        outer = TRUE, side = 3, line = 0, cex = 0.9)
  
  # Legend panel (bottom row)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
         legend = ds_test_methods,
         title = "Method",
         col = colors,
         lwd = 2,
         pch = shapes,
         bty = "n",
         horiz = TRUE,
         cex = 0.9)
}


# =============================================================================
#                             RUN SIMULATIONS                                 #
# =============================================================================
run_simulation <- function(n, N, Nsim, test_type, distributions, tol_pos , alpha_pre) {
  
  # create results folder
  if (!dir.exists("results")) dir.create("results")
  
  ## Keep a small results 
  results <- list(
    params  = list(N = N, Nsim = Nsim, test_type = test_type, distributions = distributions , alpha_pre = alpha_pre),
    objects = list(),
    files   = list()
  )
  # define run control parameters
  sample_size       <- c(10, 20, 30, 40, 50)
  test_alpha        <- 0.05
  alpha             <- 0.05
  select_norm_test  <- "SW"
  effect_size       <- 0.5
  norm_alpha_pre    <- seq(from = 0.0, to = 1, by = 0.05)
  alpha_pretest     <- seq(from = 0.009, to = 1, by = 0.0025)
  effect_sizes      <- c(0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0)
  sig_levels        <- seq(from = 0.005, to = 1, by = 0.005)
  ds_test_methods   <- c("test_1", "test_2", "adaptive")
  norm_test         <- c("SW", "SF", "KS", "JB","SKEW", "DAP", "AD", "CVM")
  
  ## ROC data for normality tests 
  roc_pval_ds_test <- fn_for_roc_curve_for_norm_test(
    n = 10,
    alpha_pretest = norm_alpha_pre,
    H1_dist = distributions[1],
    tests = norm_test,
    Nsim = Nsim
  )
  # store results
  results$objects$roc_pval_ds_test <- roc_pval_ds_test
  
  ## Plot ROC using selected tests 
  pdf(paste0("results/", test_type, "_test_norm_roc_curve.pdf"), width = 6, height = 5)
  selected_tests <- c("SW", "SF", "KS", "JB","SKEW", "DAP", "AD", "CVM")
  plot_norm_roc_curve(
    FPR = roc_pval_ds_test$FPR,
    TPR = roc_pval_ds_test$TPR,
    tests_to_plot = selected_tests,
    alpha = roc_pval_ds_test$alpha,
    dist_name = distributions[1]
  )
  dev.off()
  
  # Trade-off analysis (power & Type I error) 
  analysis_ds_tests <- perform_analysis(
    Nsim = Nsim,
    n = 10,
    distributions = distributions,
    effect_size = effect_size,
    test = select_norm_test,
    alpha_pretest = alpha_pretest,
    test_alpha = test_alpha
  )
  
  # store results
  results$objects$analysis_ds_tests <- analysis_ds_tests
  
  # Calculate metrics
  metrics <- compute_roc_metrics(
    error_ds_test = analysis_ds_tests$error_ds_test,
    power_ds_test = analysis_ds_tests$power_ds_test,
    test_alpha    = test_alpha
  )
  
  # store results
  results$objects$metrics <- metrics
  
  # Capture the trade off result to get alpha_star 
  pdf(paste0("results/", test_type, "_test_tradeoff.pdf"), width = 6, height = 5)
  tradeoff_result <- plot_power_error_tradeoff(
    distributions = distributions,
    alpha_pretest = alpha_pretest,
    metrics = metrics,
    tol_pos = tol_pos
  )
  dev.off()
  
  # Save the trade off result and extract alpha_star
  results$objects$tradeoff_result <- tradeoff_result
  alpha_star <- tradeoff_result$alpha_star
  results$objects$alpha_star <- alpha_star
  
  cat("=== Selected alpha_star from tradeoff analysis:", round(alpha_star, 4), "===\n")
  
  ## Generate data for Power vs Type I error like plot 
  roc_data <- power_vs_error_roc_data(
    N = N, 
    n = 10,
    distributions = distributions, 
    effect_size = effect_size,
    alpha_pretest = alpha_star,  
    sig_levels = sig_levels
  )
  
  # store results
  results$objects$roc_data <- roc_data
  
  #  Power vs Type I error ROC-like plots 
  pdf(paste0("results/", test_type, "_test_power_error_roc.pdf"), width = 6, height = 5)
  power_vs_error_ROC_curve(roc_data = roc_data)
  dev.off()
  
  # -------------------------------------------------------
  # ------- Main simulation: power/Type I error analysis --
  # -------------------------------------------------------
  sim_output <- perform_ds_func(
    sample_sizes       = sample_size,
    distributions      = distributions,
    N                  = N,
    alpha              = alpha,
    gen_data           = gen_data,
    get_parameters     = get_parameters,
    fn_to_get_norm_obj = fn_to_get_norm_obj,
    fn_for_norm_test   = normality_test,
    fn_for_ds_test_1   = fn_for_ds_test_1,
    fn_for_ds_test_2   = fn_for_ds_test_2,
    norm_test_method   = select_norm_test,
    ds_test_methods    = ds_test_methods,
    effect_size        = effect_size,
    alpha_pre          = alpha_star  
  )
  
  # store results
  results$objects$sim_output <- sim_output
  
  ## Per-distribution AUCs (printed and saved) 
  all_auc_tables <- list()
  for (dist in names(sim_output)) {
    cat("\n\n=== Results for", dist, "distribution ===\n")
    dist_results  <- sim_output[[dist]]
    power_results <- dist_results$power
    type1_results <- dist_results$type1
    
    # power
    power_df <- list(
      test_1    = power_results$test_1,
      test_2 = power_results$test_2,
      adaptive      = power_results$adaptive
    )
    
    # type I error
    TypeI_error_df <- list(
      test_1    = type1_results$test_1,
      test_2 = type1_results$test_2,
      adaptive      = type1_results$adaptive
    )
    
    # calculate aucs
    auc_power <- sapply(power_df, function(y) compute_area(sample_size, y))
    auc_type1 <- sapply(TypeI_error_df, function(y) compute_area(sample_size, y))
    
    # store in a dataframe
    auc_table <- data.frame(
      Method    = names(auc_power),
      AUC_Power = unname(auc_power),
      AUC_TypeI = unname(auc_type1),
      row.names = NULL
    )
    all_auc_tables[[dist]] <- auc_table
    print(auc_table, digits = 4)
  }
  
  # store results
  results$objects$auc_tables <- all_auc_tables
  
  ## Combine across distributions and plot 
  all_power_results <- list()
  all_type1_results <- list()
  
  for (dist in names(sim_output)) {
    dist_results <- sim_output[[dist]]
    
    power_df <- data.frame(
      n = sample_size,
      test_1    = dist_results$power$test_1,
      test_2 = dist_results$power$test_2,
      adaptive      = dist_results$power$adaptive,
      Distribution  = dist
    )
    type1_df <- data.frame(
      n = sample_size,
      test_1    = dist_results$type1$test_1,
      test_2 = dist_results$type1$test_2,
      adaptive      = dist_results$type1$adaptive,
      Distribution  = dist
    )
    
    all_power_results[[dist]] <- power_df
    all_type1_results[[dist]] <- type1_df
  }
  
  combined_power <- do.call(rbind, all_power_results)
  combined_type1 <- do.call(rbind, all_type1_results)
  results$objects$combined_power_n <- combined_power
  results$objects$combined_type1_n <- combined_type1
  
  # 
  pdf(paste0("results/", test_type, "_power_and_error_test.pdf"), width = 6, height = 5)
  plot_power_type1(
    combined_power  = combined_power,
    combined_type1  = combined_type1,
    optimal_alphas  = alpha_star,        
    methods         = ds_test_methods,   
    distributions   = distributions,
    sizes           = sample_size,
    alpha           = alpha
  )
  dev.off()
# -----------------------------------------------------------
  ## Effect-size power curves 
# ----------------------------------------------------------
  sample_size_es <- 10
  sim_power_list <- list()
  
  for (n in sample_size_es) {
    sim_power_list[[as.character(n)]] <- perform_ds_power_by_effect(
      fixed_n            = n,
      effect_sizes       = effect_sizes,
      distributions      = distributions,
      N                  = N,
      alpha              = alpha,
      gen_data           = gen_data,
      get_parameters     = get_parameters,
      fn_to_get_norm_obj = fn_to_get_norm_obj,
      fn_for_norm_test   = normality_test,
      fn_for_ds_test_1   = fn_for_ds_test_1,
      fn_for_ds_test_2   = fn_for_ds_test_2,
      norm_test_method   = select_norm_test,
      ds_test_methods    = ds_test_methods,
      alpha_pre          = alpha_star  
    )
  }
  
  # store results
  results$objects$sim_power_list <- sim_power_list
  
  all_power_df <- list()
  for (n in names(sim_power_list)) {
    for (dist in names(sim_power_list[[n]])) {
      pr <- sim_power_list[[n]][[dist]]$power
      all_power_df[[paste(dist, n, sep = "_")]] <- data.frame(
        d    = effect_sizes,
        test_1    = pr$test_1,
        test_2 = pr$test_2,
        adaptive      = pr$adaptive,
        Distribution  = dist,
        n             = as.integer(n)
      )
    }
  }
  
  # store results
  combined_power_effect <- do.call(rbind, all_power_df)
  results$objects$combined_power_effect <- combined_power_effect
  
  # plot power vs effect sizes
  pdf(paste0("results/", test_type, "_power_by_effect.pdf"), 6, 4)
  plot_power_by_effect_size(
    power_results   = combined_power_effect,
    distributions   = distributions,
    ds_test_methods = ds_test_methods,
    effect_sizes    = effect_sizes,
    alpha_pre       = alpha_star
  )
  dev.off()
  
  ## ---------- Save compact results + full workspace ----------
  save(results, file = paste0("results/", test_type, "_results.RData"))
  save.image(paste0("results/", test_type, "_workspace.RData"))
  
  # Final summary
  cat("\n=== SIMULATION COMPLETE ===\n")
  cat("Selected alpha_star used throughout:", round(alpha_star, 4), "\n")
  cat("Tradeoff feasibility:", if(tradeoff_result$feasible) "FEASIBLE" else "FALLBACK", "\n")
  cat("Worst-case positive inflation:", round(tradeoff_result$worst_case_positive_inflation, 6), "\n")
  cat("Inflation values (Normal, Non-normal):",
      round(tradeoff_result$inflation_normal, 6), ",",
      round(tradeoff_result$inflation_non_normal, 6), "\n")
  
  invisible(results)
}
## Example call
run_simulation(
  N = 1e2,
  Nsim = 1e2,
  test_type = "two_sample_t_vs_wilcoxon",
  distributions = c("exponential", "normal"),
  tol_pos = 1e-2,
  alpha_pre = 0.05
)

