source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
# ------------------------------------------------
#----------- Two-Sample Independent t-test -------
# ------------------------------------------------
generate_ind_ttest_data <- function(
    n1 = 20,         
    n2 = 20,         
    mean1 = 0,       
    mean2 = 0.5,     
    sd1 = 1,         
    sd2 = 1,         
    dist = "Exponential"  
) {
  
  group1 <- mean1 + sd1 * generate_data(n1, dist)
  group2 <- mean2 + sd2 * generate_data(n2, dist)
  
  return(data.frame(
    group = factor(rep(c("x", "y"), c(n1, n2))),
    value = c(group1, group2)
  ))
}


# ------ Define parameter generator for t-test ----
get_ttest_params <- function(n) {
  list(
    n1 = n,
    n2 = n,
    mean1 = 0,
    mean2 = 0.5,
    sd1 = 1,
    sd2 = 1,
    dist = "Exponential"
  )
}

# --------extract groups for normality test -------
extract_ttest_groups <- function(data) {
  split(data$value, data$group)
}


# --------Parametric: Two-sample t-test --------
two_ind_t_test <- function(data) {
  x_data <- data$value[data$group == "x"]
  y_data <- data$value[data$group == "y"]
  test_result <- t.test(x = x_data, y = y_data)
  return(list(p.value = test_result$p.value))
}


# -------Nonparametric: Mann-Whitney U Test-----
Mann_Whitney_U_test <- function(data) {
  x_data <- data$value[data$group == "x"]
  y_data <- data$value[data$group == "y"]
  test_result <- wilcox.test(x_data, y_data)
  return(list(p.value = test_result$p.value))
}

# Paired t-test
generate_paired_data <- function(
    n = 20,            
    mean_diff = 0.3,   
    sd_diff = 1,       
    dist = "Normal"    
) {
  # Generate standardized differences and scale
  diff <- mean_diff + sd_diff * generate_data(n, dist)
  
  return(diff)
}
# ------------------------------------------------
#            Simple Linear Regression
# ------------------------------------------------
# Simple Linear Regression
generate_regression_data <- function(
    n = 30,           
    beta0 = 0,        
    beta1 = 0.5, 
    x_dist = "Exponential",
    error_sd = 1,     
    error_dist = "Normal"  
) {
  # predictor
  x <- generate_data(n, dist = x_dist)
  # error
  error <- error_sd * generate_data(n, error_dist)
  # independent variable
  y <- beta0 + beta1 * x + error
  return(data.frame(x = x, y = y))
}

#-----get residuals from regression model-------
regression_residuals <- function(data){
  model <- lm(y ~ x , data = data)
  return(residuals(model))
}

# 1. Define parameter generator function
get_regression_params <- function(n) {
  list(
    n = n,
    beta0 = 0,
    beta1 = 0.5,
    x_dist = "Exponential",
    error_sd = 1,
    error_dist = "Normal"
  )
}

#------ perform Simple linear regression Analysis ----
simple_linear_regression <- function(data) {
  model <- lm(y ~ x, data = data)
  tidy_model <- broom::tidy(model)
  p_value <- tidy_model$p.value[tidy_model$term == "x"]
  return(list(p.value = p_value))
}

#---------- Bootstrap regression(nonparametric) -----
bootstrap_regression <- function(data, n_boot = 1000) {
  if (!require(boot)) install.packages("boot")
  library(boot)
  
  # 1. Fit H0: y = μ + error  → get residuals from intercept‐only model
  mu_hat   <- mean(data$y)
  resid0   <- data$y - mu_hat
  df_boot  <- data.frame(x = data$x, resid0 = resid0)
  
  # 2. Observed slope
  observed_slope <- coef(lm(y ~ x, data = data))[2]
  
  # 3. Statistic: sample residuals with replacement, add to μ, refit slope
  boot_stat <- function(df, idx) {
    y_star <- mu_hat + df$resid0[idx]    # bootstrap y under H0
    model  <- lm(y_star ~ df$x)         # refit slope
    coef(model)[2]
  }
  
  # 4. Run bootstrap
  boot_results <- boot::boot(df_boot, statistic = boot_stat, R = n_boot)
  
  # 5. Two‐sided p‐value
  p_value <- mean(abs(boot_results$t) >= abs(observed_slope))
  
  return(list(p.value = p_value))
}

# ----------------------------------------------------
#                   One-Way ANOVA
# ----------------------------------------------------

generate_anova_data <- function(
    n_per_group = 20,  
    means = c(0, 0.2, 0.4),  
    sd = 1,            
    dist = "Normal"     
) {
  k <- length(means)
  group_labels <- LETTERS[1:k]
  
  # Generate data for each group
  values <- unlist(lapply(means, function(m) {
    m + sd * generate_data(n_per_group, dist)
  }))
  
  return(data.frame(
    group = factor(rep(group_labels, each = n_per_group)),
    value = values
  ))
}

#--------------- ANOVA residuals ----------------
anova_residuals <- function(data){
  return(residuals(aov(formula = value ~ group, data = data)))
}

# -------- Define ANOVA parameter generator -----
get_anova_params <- function(n) {
  list(
    n_per_group = n,
    means = c(0.0, 0.5, 2.5),
    sd = 1,
    dist = "Exponential"
  )
}

#------------- ANOVA Test ---------------------
anova_main_test <- function(data) {
  aov_model <- aov(value ~ group, data = data)
  p_value <- summary(aov_model)[[1]]$"Pr(>F)"[1]
  return(list(p.value = p_value))
}

#------- Kruskal-Wallis Test(Nonparameteric)----
kruskal_wallis_test <- function(data) {
  test_result <- kruskal.test(value ~ group, data = data)
  return(list(p.value = test_result$p.value))
}

# --------------------------------------------------
#                 Logistic Regression
# --------------------------------------------------
 
generate_logistic_data <- function(
    n = 100,               
    beta0 = -1,            
    beta1 = 0.8,           
    x_sd = 1               
) {
  x <- rnorm(n, sd = x_sd)
  log_odds <- beta0 + beta1 * x
  prob <- plogis(log_odds)  
  y <- rbinom(n, size = 1, prob = prob)
  
  return(data.frame(x = x, y = factor(y)))
}
# --------------------------------------------------------
# ---------- Normality check function --------------------
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


# compute the area under the curve 
compute_area <- function(sample_sizes, y) {
  sum(diff(sample_sizes) * (head(y, -1) + tail(y, -1)) / 2) /
    (max(sample_sizes) - min(sample_sizes))
}

# ========================================================
user_interface <- function(
    gen_data = generate_regression_data,
    fn_to_get_norm_obj = regression_residuals,
    fn_for_norm_test = normality_test,
    fn_for_ds_test_1 = simple_linear_regression,
    fn_for_ds_test_2 = bootstrap_regression,
    paras = NULL,
    alpha = 0.05,
    test_method = "SW",
    ...
) {
  # generate dataset
  data <- if (!is.null(paras)) do.call(gen_data, paras) else gen_data()
  
  # residuals & normality check
  resid     <- fn_to_get_norm_obj(data)
  norm_res  <- fn_for_norm_test(data = resid, test = test_method, alpha = alpha)
  
  # choose test
  if (isTRUE(norm_res$normality_satisfied)) {
    result <- fn_for_ds_test_1(data)        
  } else {
    result <- fn_for_ds_test_2(data, ...)   
  }
  
  return(result$p.value)
}

calculate_power <- function(
    sample_sizes    = c(10, 20, 30, 40, 50),
    n_sim           = 1e3,
    alpha           = 0.05,
    n_boot          = NULL,
    gen_data        = generate_regression_data,
    get_parameters  = function(n) list(n = n),
    fn_to_get_norm_obj   = regression_residuals,
    fn_for_norm_test     = normality_test,
    fn_for_ds_test_1     = simple_linear_regression,
    fn_for_ds_test_2     = bootstrap_regression,
    test_method    = "SW",
    mode           = c("adaptive", "parametric", "nonparametric"),
    ...
) {
  mode <- match.arg(mode)
  power_results <- numeric(length(sample_sizes))
  names(power_results) <- paste0("n=", sample_sizes)
  
  for (i in seq_along(sample_sizes)) {
    n          <- sample_sizes[i]
    rejections <- 0
    
    for (sim in seq_len(n_sim)) {
      paras <- get_parameters(n)
      # generate data
      data <- if (!is.null(paras)) do.call(gen_data, paras) else gen_data()
      
      p_value <- tryCatch({
        if (mode == "adaptive") {
          user_interface(
            gen_data           = gen_data,
            fn_to_get_norm_obj = fn_to_get_norm_obj,
            fn_for_norm_test   = fn_for_norm_test,
            fn_for_ds_test_1   = fn_for_ds_test_1,
            fn_for_ds_test_2   = fn_for_ds_test_2,
            paras              = paras,
            alpha              = alpha,
            test_method        = test_method,
            ...
          )
        } else if (mode == "parametric") {
          fn_for_ds_test_1(data)$p.value
        } else {
          fn_for_ds_test_2(data, ...)$p.value
        }
      }, error = function(e) NA)
      
      if (!is.na(p_value) && p_value < alpha) {
        rejections <- rejections + 1
      }
    }
    
    power_results[i] <- rejections / n_sim
    cat("Completed n =", n, "| mode =", mode, "| Power =", power_results[i], "\n")
  }
  
  return(power_results)
}

# ------------ Power Curve Plotting -----------
plot_power_curve <- function(power_results, title = "Power Analysis") {
  sample_sizes <- as.numeric(sub("n=", "", names(power_results)))
  plot(sample_sizes, power_results, type = "b",
       xlab = "Sample Size", ylab = "Power",
       main = title, ylim = c(0, 1), pch = 19)
  abline(h = 0.8, lty = 2)
  grid()
}

# ===================================================