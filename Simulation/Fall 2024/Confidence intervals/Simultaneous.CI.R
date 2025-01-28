# Load required libraries
library(ggplot2)

# Function to generate data, create QQ plot with confidence intervals, and count points outside bounds
generate_qq_plot_with_ci <- function(sample_size, dist_name, dist_params, confidence_level = 0.95) {
  
  # Generate data based on specified distribution
  data <- switch(
    dist_name,
    "exponential" = rexp(sample_size, rate = dist_params$rate),
    "normal" = rnorm(sample_size, mean = dist_params$mean, sd = dist_params$sd),
    "uniform" = runif(sample_size, min = dist_params$min, max = dist_params$max),
    "gamma" = rgamma(sample_size, shape = dist_params$shape, rate = dist_params$rate),
    stop("Unsupported distribution. Add your distribution to the 'switch' statement.")
  )
  
  # Calculate theoretical quantiles
  theoretical_quantiles <- switch(
    dist_name,
    "exponential" = qexp(ppoints(sample_size), rate = dist_params$rate),
    "normal" = qnorm(ppoints(sample_size), mean = dist_params$mean, sd = dist_params$sd),
    "uniform" = qunif(ppoints(sample_size), min = dist_params$min, max = dist_params$max),
    "gamma" = qgamma(ppoints(sample_size), shape = dist_params$shape, rate = dist_params$rate),
    stop("Unsupported distribution. Add your distribution to the 'switch' statement.")
  )
  
  # Sort the data
  sorted_data <- sort(data)
  
  # Calculate simultaneous confidence intervals using the Kolmogorov-Smirnov method
  alpha <- 1 - confidence_level
  ks_constant <- sqrt(-0.5 * log(alpha / 2))
  ci_upper <- theoretical_quantiles + ks_constant / sqrt(sample_size)
  ci_lower <- theoretical_quantiles - ks_constant / sqrt(sample_size)
  
  # Count points outside the confidence intervals
  outside_bounds <- sum(sorted_data < ci_lower | sorted_data > ci_upper)
  
  # Create a data frame for plotting
  qq_data <- data.frame(
    Theoretical = theoretical_quantiles,
    Sample = sorted_data,
    CI_Lower = ci_lower,
    CI_Upper = ci_upper
  )
  
  # Create QQ plot with confidence intervals
  qq_plot <- ggplot(qq_data, aes(x = Theoretical, y = Sample)) +
    geom_point(color = "blue", size = 2) +
    geom_line(aes(x = Theoretical, y = Theoretical), color = "red", linetype = "dashed") +  # Line y=x
    geom_ribbon(aes(ymin = CI_Lower, ymax = CI_Upper), alpha = 0.2, fill = "grey") +  # Confidence intervals
    labs(
      title = paste("QQ Plot with Simultaneous Confidence Intervals\n", 
                    "Distribution:", dist_name),
      x = "Theoretical Quantiles",
      y = "Sample Quantiles"
    ) +
    theme_minimal()
  
  # Print the plot
  print(qq_plot)
  
  # Return the number of points outside the confidence intervals
  return(outside_bounds)
}

# Example Usage
set.seed(123)  # For reproducibility

# Exponential distribution example
n_outside <- generate_qq_plot_with_ci(
  sample_size = 30, 
  dist_name = "exponential", 
  dist_params = list(rate = 1), 
  confidence_level = 0.95
)
cat("Number of points outside the confidence bounds (Exponential):", n_outside, "\n")

# Normal distribution example
n_outside <- generate_qq_plot_with_ci(
  sample_size = 10, 
  dist_name = "normal", 
  dist_params = list(mean = 0, sd = 1), 
  confidence_level = 0.95
)
cat("Number of points outside the confidence bounds (Normal):", n_outside, "\n")

# Uniform distribution example
n_outside <- generate_qq_plot_with_ci(
  sample_size = 100, 
  dist_name = "uniform", 
  dist_params = list(min = 0, max = 1), 
  confidence_level = 0.95
)
cat("Number of points outside the confidence bounds (Normal):", n_outside, "\n")

