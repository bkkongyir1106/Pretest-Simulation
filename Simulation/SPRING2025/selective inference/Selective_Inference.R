# Load required libraries
library(rstanarm)
library(loo)
library(bayesplot)
library(ggplot2)
library(tidyr)
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
# Simulate example data
set.seed(123)
n <- 10
data <- data.frame(
  x = rnorm(n),      # Normal group
  y = rt(n, df = 3)  # Non-normal group (t-distribution)
)

# 1. Fit Bayesian Models with Correct Specifications
# -------------------------------------------------
# Gaussian model
normal_model <- stan_glm(
  x ~ y,
  data = data,
  family = gaussian(),  # Normal distribution
  prior = normal(0, 2.5),
  prior_intercept = normal(0, 2.5),
  seed = 123
)

# Robust Student-t regression (using Gaussian family with Student-t priors)
robust_model <- stan_glm(
  x ~ y,
  data = data,
  family = gaussian(),  # Keep Gaussian, but use robust priors
  prior = student_t(4, 0, 2.5),
  prior_intercept = student_t(4, 0, 2.5),
  seed = 123
)

# 2. Calculate Model Weights Using LOO-CV
# -------------------------------------------------
# Compute LOO results
loo_normal <- loo(normal_model)
loo_robust <- loo(robust_model)

# Compute model weights
model_weights <- loo_model_weights(list(loo_normal, loo_robust))

cat("Model weights:\n",
    "- Normal model:", model_weights[1], "\n",
    "- Robust model:", model_weights[2], "\n")

# 3. Model-Averaged Posterior Distribution
# -------------------------------------------------
# Extract posterior samples for y coefficient
posterior_normal <- as.matrix(normal_model)[, "y"]
posterior_robust <- as.matrix(robust_model)[, "y"]

# Calculate weighted average
posterior_avg <- model_weights[1] * posterior_normal + 
  model_weights[2] * posterior_robust

# 4. Visualization
# -------------------------------------------------
# Reshape the data for ggplot
plot_data <- data.frame(
  Value = c(posterior_normal, posterior_robust, posterior_avg),
  Model = rep(c("Normal", "Robust", "Average"), each = length(posterior_normal))
)

# Plot the densities
ggplot(plot_data, aes(x = Value, color = Model)) +
  geom_density(size = 1) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  scale_color_manual(
    values = c("Normal" = "#1f77b4",
               "Robust" = "#ff7f0e",
               "Average" = "#2ca02c")
  ) +
  labs(title = "Model-Averaged Posterior Distribution of Effect Size",
       x = "Effect Size (y coefficient)",
       y = "Density") +
  theme_minimal() +
  theme(legend.position = "bottom")

# -----------------------------------------------------------------------
# -----------------------------Bootstrap Procedure ----------------------
# -----------------------------------------------------------------------
# Required Libraries
library(boot)
library(ggplot2)

# Bootstrap-Adjusted Selective Inference
adjusted_t_test <- function(data, B = as.integer(1e4), alpha_pre = 0.05) {
  
  original_test <- function(d) {
    # Pre-test for normality
    if (shapiro_p <- shapiro.test(d)$p.value > alpha_pre) {
      test_res <- t.test(d)
      return(test_res$statistic)
    } else {
      test_res <- wilcox.test(d)
      return(test_res$statistic)  
    }
  }
  
  # Observed test statistic
  T_obs <- original_test(data)
  
  # Bootstrap under H0 (centered data)
  centered_data <- data - mean(data)
  
  boot_dist <- replicate(B, {
    # Resample under H0
    bs_sample <- sample(centered_data, size = length(centered_data), replace = TRUE)
    
    if (shapiro_p_bs <- shapiro.test(bs_sample)$p.value > alpha_pre) {
      t.test(bs_sample)$statistic
    } else {
      wilcox.test(bs_sample)$statistic
    }
  })
  
  # Calculate adjusted p-value
  p_adj <- mean(abs(boot_dist) >= abs(T_obs))
  
  return(list(
    original_stat = T_obs,
    adjusted_p = p_adj,
    bootstrap_dist = boot_dist
  ))
}

# Simulation Example
set.seed(12345)
sample_size <- 30
dist = "Normal"
data <- generate_data(n, dist)

# Run analysis
results <- adjusted_t_test(data, B = 1e3)

# View results
cat("Original test statistic:", results$original_stat, "\n")
cat("Adjusted p-value:", results$adjusted_p, "\n")

# Plot bootstrap distribution
ggplot(data.frame(stat = results$bootstrap_dist), aes(stat)) +
  geom_histogram(aes(y = ..density..), bins = 30, fill = "steelblue", alpha = 0.7) +
  geom_vline(xintercept = results$original_stat, color = "red", linetype = "dashed") +
  labs(title = "Bootstrap Distribution of Test Statistics",
       x = "Test Statistic", y = "Density") +
  theme_minimal()

# -------------------------- Type I Error Simulation ------------------------

sim_typeI <- function(n_sims = as.integer(1000), n = sample_size) {
  false_positives <- 0
  for(i in 1:n_sims) {
    data <- generate_data(n, dist)
    res <- adjusted_t_test(data, B = as.integer(1000))
    if(res$adjusted_p < 0.05) false_positives <- false_positives + 1
  }
  cat("Empirical Type I Error Rate:", false_positives / n_sims, "\n")
}
sim_typeI()

# ---------------------------- Power Analysis -------------------------------

sim_power <- function(n_sims = 1e3, n = sample_size , effect_size = 0.5) {
  true_positives <- 0
  for(i in 1:n_sims) {
    data <- generate_data(n, dist) + effect_size
    res <- adjusted_t_test(data, B = 1e3)
    if(res$adjusted_p < 0.05) true_positives <- true_positives + 1
  }
  cat("Empirical Power:", true_positives / n_sims, "\n")
}
sim_power()

