# Two-Sample Independent t-test
generate_ind_ttest_data <- function(
    n1 = 20,         
    n2 = 20,         
    mean1 = 0,       
    mean2 = 0.5,     
    sd1 = 1,         
    sd2 = 1,         
    dist = "Normal"  
) {
  # Generate standardized data and scale to desired parameters
  group1 <- mean1 + sd1 * generate_data(n1, dist)
  group2 <- mean2 + sd2 * generate_data(n2, dist)
  
  return(data.frame(
    group = factor(rep(c("x", "y"), c(n1, n2))),
    value = c(group1, group2)
  ))
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

# Simple Linear Regression
generate_regression_data <- function(
    n = 30,           
    beta0 = 0,        
    beta1 = 0.5, 
    x_dist = "Exponential",
    error_sd = 1,     
    error_dist = "Normal"  
) {
  # Generate predictor
  x <- generate_data(n, dist = x_dist)
  
  # Generate standardized errors and scale
  error <- error_sd * generate_data(n, error_dist)
  
  y <- beta0 + beta1 * x + error
  return(data.frame(x = x, y = y))
}


# One-Way ANOVA
generate_anova_data <- function(
    n_per_group = 20,  # Sample size per group
    means = c(0, 0.2, 0.4),  # Group means
    sd = 1,            # Common SD
    dist = "Normal"     # Error distribution
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


# Logistic Regression (Unchanged - binary response)
generate_logistic_data <- function(
    n = 100,               # Sample size
    beta0 = -1,            # Intercept
    beta1 = 0.8,           # Slope (effect size)
    x_sd = 1               # Predictor SD
) {
  x <- rnorm(n, sd = x_sd)
  log_odds <- beta0 + beta1 * x
  prob <- plogis(log_odds)  # More efficient than 1/(1+exp(-x))
  y <- rbinom(n, size = 1, prob = prob)
  
  return(data.frame(x = x, y = factor(y)))
}

#--------functions to generate data for normality test ---------
regression_residuals <- function(data){
  model <- lm(y ~ x , data = data)
  return(residuals(model))
}

# ANOVA
anova_residuals <- function(data){
  return(residuals(aov(formula = value ~ group, data = data)))
}


# -----------------Downstream test ---------------
# Two-sample t-test (fixed)
two_ind_t_test <- function(data) {
  x_data <- data$value[data$group == "x"]
  y_data <- data$value[data$group == "y"]
  return(t.test(x = x_data, y = y_data))
}

# One-sample t-test (fixed)
one_sample_t_test <- function(data) {
  return(t.test(data))
}

# Simple linear regression
simple_linear_regression <- function(data) {
  model <- lm(y ~ x, data = data)
  return(broom::tidy(model))
}

# Mann-Whitney U Test (fixed)
Mann_Whitney_U_test <- function(data){
  x_data <- data$value[data$group == "x"]
  y_data <- data$value[data$group == "y"]
  return(wilcox.test(x_data, y_data))
}

# ANOVA Test (fixed to return actual test results)
anova_main_test <- function(data) {
  aov_model <- aov(value ~ group, data = data)
  return(summary(aov_model))
}



