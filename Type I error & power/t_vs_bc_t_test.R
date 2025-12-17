library(MASS)
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/user_defined_functions_center_median.R")

# Box-Cox transformation 
boxcox_simple <- function(data) {
  # Shift data to be positive
  shift <- 0
  if (min(data) <= 0) {
    shift <- abs(min(data)) + 0.001
    data <- data + shift
  }
  
  # Find optimal lambda 
  bc <- suppressMessages(boxcox(data ~ 1, lambda = seq(-2, 2, 0.1), plotit = FALSE))
  lambda <- bc$x[which.max(bc$y)]
  
  # Transform data
  if (lambda == 0) {
    transformed_data <- log(data)
  } else {
    transformed_data <- (data^lambda - 1) / lambda
  }
  
  return(list(
    transformed_data = transformed_data,
    lambda = lambda,
    shift = shift
  ))
}

# T-test on median using Box-Cox transformed data
test_t_on_median <- function(data) {
  # Since data is already centered around median, true_median = 0
  true_median <- 0
  
  # Transform data
  bc <- boxcox_simple(data)
  
  # Transform the true median 
  true_median_shifted <- true_median + bc$shift
  if (bc$lambda == 0) {
    transformed_median <- log(true_median_shifted)
  } else {
    transformed_median <- (true_median_shifted^bc$lambda - 1) / bc$lambda
  }
  
  # T-test on transformed data against transformed median
  t.test(bc$transformed_data, mu = transformed_median)$p.value
}

# Calculate Type I error rate for different distributions
type1_error_all_distributions <- function(n_sim = 100) {
  sizes <- c(10, 20, 30, 50, 100)
  distributions <- c("normal", "lognormal", "exponential", "chi_square")
  
  # Initialize results
  results <- data.frame(
    distribution = rep(distributions, each = length(sizes)),
    n = rep(sizes, times = length(distributions)),
    error_rate = NA,
    avg_lambda = NA
  )
  
  for (dist in distributions) {
    for (i in 1:length(sizes)) {
      n <- sizes[i]
      errors <- 0
      lambdas <- numeric(n_sim)
      
      for (j in 1:n_sim) {
        # Generate data centered around median using your function
        data <- generate_data(n, dist)
        
        # T-test on median using Box-Cox
        p_value <- test_t_on_median(data)
        if (p_value < 0.05) errors <- errors + 1
        
        lambdas[j] <- boxcox_simple(data)$lambda
      }
      
      # Store results
      idx <- which(results$distribution == dist & results$n == n)
      results$error_rate[idx] <- errors / n_sim
      results$avg_lambda[idx] <- mean(lambdas)
    }
  }
  
  return(results)
}

# Run it for all distributions
set.seed(123)
all_results <- type1_error_all_distributions(n_sim = 1000)
print(all_results)

# Combined plot for all distributions
plot_combined <- function(results) {
  distributions <- unique(results$distribution)
  colors <- c("red", "blue", "green", "purple")
  # Get the sample sizes from the results
  sizes <- unique(results$n)
  # Set up the plot
  plot(0, 0, type = "n", 
       xlim = range(sizes), ylim = c(0, 1),
       xlab = "Sample Size", ylab = "Type I Error Rate",
       main = "Type I Error Rates - Box-Cox Transformed T-test")
  
  # Add line for each distribution
  for (i in 1:length(distributions)) {
    dist <- distributions[i]
    dist_data <- results[results$distribution == dist, ]
    lines(dist_data$n, dist_data$error_rate, 
          type = "b", col = colors[i], lwd = 2, pch = 16)
  }
  
  # Add horizontal line at nominal level
  abline(h = 0.05, lty = 2, col = "gray", lwd = 2)
  
  # Add legend
  legend("topleft", legend = distributions, 
         col = colors, lwd = 2, pch = 16, bty = "b")
}

# Create combined plot
pdf(file = "one_sample_t_test_vs_bc_t_test.pdf",width = 6, height = 5 )
plot_combined(all_results)
dev.off()
# =============================================================================

## ==================== Two sample Case ====================

# Two-sample t-test on medians using Box-Cox transformed data
test_twosample_t_on_median <- function(data1, data2) {
  # Transform both groups
  bc1 <- boxcox_simple(data1)
  bc2 <- boxcox_simple(data2)
  
  # transform median 0 for both
  median_shifted1 <- 0 + bc1$shift
  median_shifted2 <- 0 + bc2$shift
  
  if (bc1$lambda == 0) {
    transformed_median1 <- log(median_shifted1)
  } else {
    transformed_median1 <- (median_shifted1^bc1$lambda - 1) / bc1$lambda
  }
  
  if (bc2$lambda == 0) {
    transformed_median2 <- log(median_shifted2)
  } else {
    transformed_median2 <- (median_shifted2^bc2$lambda - 1) / bc2$lambda
  }
  
  # Test if the transformed medians are equal
  t.test(bc1$transformed_data, bc2$transformed_data)$p.value
}

# Calculate Type I error rate 
type1_error_twosample <- function(n_sim = 1000) {
  sizes <- c(10, 20, 30, 50)
  distributions <- c("normal", "lognormal", "exponential", "chi_square")
  
  # Initialize results
  results <- data.frame(
    distribution = rep(distributions, each = length(sizes)),
    n = rep(sizes, times = length(distributions)),
    error_rate = NA,
    avg_lambda1 = NA,
    avg_lambda2 = NA
  )
  
  for (dist in distributions) {
    for (i in 1:length(sizes)) {
      n <- sizes[i]
      errors <- 0
      lambdas1 <- numeric(n_sim)
      lambdas2 <- numeric(n_sim)
      
      for (j in 1:n_sim) {
        # Generate two samples from the same distribution
        data1 <- generate_data(n, dist)
        data2 <- generate_data(n, dist)
        
        # Two-sample t-test on medians using Box-Cox
        p_value <- test_twosample_t_on_median(data1, data2)
        if (p_value < 0.05) errors <- errors + 1
        
        lambdas1[j] <- boxcox_simple(data1)$lambda
        lambdas2[j] <- boxcox_simple(data2)$lambda
      }
      
      # Store results
      idx <- which(results$distribution == dist & results$n == n)
      results$error_rate[idx] <- errors / n_sim
      results$avg_lambda1[idx] <- mean(lambdas1)
      results$avg_lambda2[idx] <- mean(lambdas2)
    }
  }
  
  return(results)
}

# Run it for all distributions
set.seed(123)
twosample_results <- type1_error_twosample(n_sim = 1000)
print(twosample_results)



# Combined plot for all distributions
plot_combined <- function(results) {
  distributions <- unique(results$distribution)
  colors <- c("red", "blue", "green", "purple")
  # Get the sample sizes from the results
  sizes <- unique(results$n)
  # Set up the plot
  plot(0, 0, type = "n", 
       xlim = range(sizes), ylim = c(0, 1),
       xlab = "Sample Size", ylab = "Type I Error Rate",
       main = "Two-Samples:Type I Error Rates - Box-Cox Transformed T-test")
  
  # Add line for each distribution
  for (i in 1:length(distributions)) {
    dist <- distributions[i]
    dist_data <- results[results$distribution == dist, ]
    lines(dist_data$n, dist_data$error_rate, 
          type = "b", col = colors[i], lwd = 2, pch = 16)
  }
  
  # Add horizontal line at nominal level
  abline(h = 0.05, lty = 2, col = "gray", lwd = 2)
  
  # Add legend
  legend("topleft", legend = distributions, 
         col = colors, lwd = 2, pch = 16, bty = "b")
}

# Create combined plot
pdf(file = "two_sample_t_test_vs_bc_t_test.pdf",width = 6, height = 5 )
plot_combined(twosample_results)
dev.off()
