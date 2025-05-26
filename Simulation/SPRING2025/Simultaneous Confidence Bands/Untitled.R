# Load necessary libraries
library(ggplot2)
library(qqplotr)
#install.packages("qqplotr")  # Install the package if not already installed

# Function to compute simultaneous confidence bands for Q-Q plot
qq_plot_with_simultaneous_ci <- function(data, title = "Q-Q Plot with Simultaneous Confidence Bands", B = 1000, conf_level = 0.95) {
  n <- length(data)
  sorted_data <- sort(data)  # Sort the sample data
  theoretical_quantiles <- qnorm(ppoints(n))  # Theoretical quantiles under normality
  
  # Bootstrap to construct confidence bands
  boot_samples <- replicate(B, sort(rnorm(n)))
  lower_band <- apply(boot_samples, 1, function(x) quantile(x, (1 - conf_level) / 2))
  upper_band <- apply(boot_samples, 1, function(x) quantile(x, 1 - (1 - conf_level) / 2))
  
  # Create a data frame
  df <- data.frame(Theoretical = theoretical_quantiles, Sample = sorted_data,
                   Lower = lower_band, Upper = upper_band)
  
  # Create the Q-Q plot with confidence bands
  p <- ggplot(df, aes(x = Theoretical, y = Sample)) +
    geom_point(color = "black") +  # Q-Q points
    geom_abline(slope = 1, intercept = 0, color = "red", linetype = "dashed") +  # Reference line
    geom_ribbon(aes(ymin = Lower, ymax = Upper), fill = "blue", alpha = 0.2) +  # Simultaneous CI bands
    theme_minimal() +
    labs(title = title, x = "Theoretical Quantiles (Normal)", y = "Sample Quantiles")
  
  print(p)
}

# Example usage with a normally distributed sample
set.seed(123)
sample_data <- rnorm(100)  # Normal data
qq_plot_with_simultaneous_ci(sample_data)

# Example usage with a non-normal sample
sample_data_non_normal <- rt(100, df = 3)  # Heavy-tailed t-distribution
qq_plot_with_simultaneous_ci(sample_data_non_normal)

# Example usage with a non-normal sample
sample_data_non_normal <- rexp(100, rate = 1)  # Heavy-tailed t-distribution
qq_plot_with_simultaneous_ci(sample_data_non_normal)
