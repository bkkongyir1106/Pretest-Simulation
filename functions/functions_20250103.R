# Load necessary libraries
library("nortest")
library("dgof")
library("dplyr")
library("moments")
library("tseries")
library("LaplacesDemon")
library("VGAM")

# Generate data from different distributions
generate_data <- function(n, dist) {
  if (!is.numeric(n) || n <= 0) stop("Sample size 'n' must be a positive integer.")
  
  x <- switch(
    dist,
    "Standard Normal" = rnorm(n, mean = 0, sd = 1),
    "Chi-Square" = (rchisq(n, df = 3) - 3) / sqrt(6),
    "Gamma" = (rgamma(n, shape = 3, rate = 0.1) - 30) / sqrt(300),
    "Exponential" = rexp(n, rate = 1) - 1,
    "t" = rt(n, df = 7) / sqrt(7 / 5),
    "Uniform" = (runif(n, min = 0, max = 1) - 0.5) * sqrt(12),
    "Laplace" = rlaplace(n, location = 0, scale = 4) / sqrt(8),
    "Weibull" = {
      shape <- 1
      scale <- 2
      variance <- scale^2 * (gamma(1 + 2 / shape) - (gamma(1 + 1 / shape))^2)
      (rweibull(n, shape = shape, scale = scale) - scale * gamma(1 + 1 / shape)) / sqrt(variance)
    },
    "LogNormal" = {
      meanlog <- 0
      sdlog <- 1
      (rlnorm(n, meanlog = meanlog, sdlog = sdlog) - exp(meanlog + sdlog^2 / 2)) /
        sqrt((exp(sdlog^2) - 1) * exp(2 * meanlog + sdlog^2))
    },
    "Contaminated" = {
      br <- rbinom(n, size = 1, prob = 0.75)
      sd_br <- sqrt(1 + br * 24)
      rnorm(n, mean = 0, sd = sd_br) / sqrt(0.25 + 0.75 * 24)
    },
    "Pareto" = {
      shape <- 3
      if (shape <= 2) stop("Shape parameter for Pareto must be greater than 2.")
      mean_pareto <- shape / (shape - 1)
      sd_pareto <- sqrt(shape / ((shape - 1)^2 * (shape - 2)))
      (rpareto(n, shape = shape) - mean_pareto) / sd_pareto
    },
    stop("Invalid distribution type.")
  )
  
  return(x)
}

# Apply different normality tests
generate_tests <- function(x, test) {
  valid_tests <- c("KS", "SW", "JB", "DAP", "AD", "SF", "CVM", "CHISQ")
  if (!test %in% valid_tests) stop("Invalid test type.")
  
  output <- switch(
    test,
    "KS" = lillie.test(x),
    "SW" = shapiro.test(x),
    "JB" = jarque.test(x),
    "DAP" = agostino.test(x),
    "AD" = ad.test(x),
    "SF" = sf.test(x),
    "CVM" = cvm.test(x),
    "CHISQ" = chisq.test(x),
    stop("Invalid test type.")
  )
  
  return(output)
}

# One-sample test statistic
OneSample_test_statistic <- function(x) {
  return((mean(x) * sqrt(length(x))) / sd(x))
}

# Two-sample test statistic
TwoSample_test_statistic <- function(x, y) {
  pooled_var <- ((length(x) - 1) * var(x) + (length(y) - 1) * var(y)) / (length(x) + length(y) - 2)
  return((mean(x) - mean(y)) / sqrt(pooled_var / length(x) + pooled_var / length(y)))
}

# Compute area under the curve
compute_area <- function(x, y) {
  if (length(x) != length(y)) stop("Vectors x and y must have the same length.")
  return(sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2))
}

# One-sample test
OneSample_test <- function(x, test, alpha = 0.05, B = 1000) {
  valid_tests <- c("t", "Wilcox", "t_Wilcox", "perm")
  if (!test %in% valid_tests) stop("Invalid test type.")
  
  pval <- switch(
    test,
    "t" = t.test(x)$p.value,
    "Wilcox" = wilcox.test(x)$p.value,
    "t_Wilcox" = {
      if (shapiro.test(x)$p.value > alpha) t.test(x)$p.value else wilcox.test(x)$p.value
    },
    "perm" = {
      observe_stat <- OneSample_test_statistic(x)
      permuted_stat <- replicate(B, {
        index <- sample(c(-1, 1), length(x), replace = TRUE)
        sample_data <- index * abs(x)
        OneSample_test_statistic(sample_data)
      })
      mean(abs(permuted_stat) >= abs(observe_stat))
    },
    stop("Invalid test type.")
  )
  
  return(pval)
}

# Two-sample test
TwoSample_test <- function(x, y, test, alpha = 0.05, B = 1000) {
  valid_tests <- c("t", "Wilcox", "t_Wilcox", "perm")
  if (!test %in% valid_tests) stop("Invalid test type.")
  
  pval <- switch(
    test,
    "t" = t.test(x, y)$p.value,
    "Wilcox" = wilcox.test(x, y)$p.value,
    "t_Wilcox" = {
      if (shapiro.test(x)$p.value > alpha && shapiro.test(y)$p.value > alpha) {
        t.test(x, y)$p.value
      } else {
        wilcox.test(x, y)$p.value
      }
    },
    "perm" = {
      observe_stat <- TwoSample_test_statistic(x, y)
      data <- c(x, y)
      permuted_stat <- replicate(B, {
        sample_data <- sample(data)
        sample_x <- sample_data[1:length(x)]
        sample_y <- sample_data[(length(x) + 1):(length(x) + length(y))]
        TwoSample_test_statistic(sample_x, sample_y)
      })
      mean(abs(permuted_stat) >= abs(observe_stat))
    },
    stop("Invalid test type.")
  )
  
  return(pval)
}

# Bootstrap Power Calculation Function
bootstrap_power <- function(data, null_value, num_iterations = 1000, alpha = 0.05) {
  sample_size <- length(data)
  bootstrap_means <- numeric(num_iterations)
 
  for (i in 1:num_iterations) {
    resampled_data <- sample(data, size = sample_size, replace = TRUE)
    bootstrap_means[i] <- OneSample_test_statistic(resampled_data)
  }

  bootstrap_ci <- quantile(bootstrap_means, c(alpha / 2, 1 - alpha / 2))
  ci_excludes_null <- bootstrap_ci[1] > null_value || bootstrap_ci[2] < null_value
 
  p_value <- mean(bootstrap_means < null_value | bootstrap_means > null_value)
  
  list(
    bootstrap_confidence_interval = bootstrap_ci,
    p_value = p_value,
  )
}

set.seed(123)  
data <- data$x
null_value <- 0.5 

bootstrap_power(data, null_value = 0.5, num_iterations = 1000, alpha = 0.05)
results <- bootstrap_power(data, null_value)
print(results)
