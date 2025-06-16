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
  
  # Create correlated pre/post measurements
  pre <- rnorm(n)  # Baseline from standard normal
  post <- pre + diff
  
  return(data.frame(
    subject = factor(rep(1:n, 2)),
    time = rep(c("pre", "post"), each = n),
    value = c(pre, post)
  ))
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

# Chi-Square Test of Independence 
generate_chisq_data <- function(
    n = 100,               # Total sample size
    p_row = c(0.3, 0.7),   # Row marginal probabilities
    effect_size = 0.3      # Difference in conditional probabilities
) {
  if (length(p_row) != 2) stop("Only 2x2 tables supported")
  
  # Calculate cell probabilities based on effect size
  p_col1 <- c(p_row[1] * (0.5 + effect_size/2), 
              p_row[2] * (0.5 - effect_size/2))
  p_col2 <- c(p_row[1] * (0.5 - effect_size/2), 
              p_row[2] * (0.5 + effect_size/2))
  
  # Generate counts
  counts <- rmultinom(1, n, c(p_col1, p_col2))[ ,1]
  
  return(matrix(counts, nrow = 2, byrow = TRUE,
                dimnames = list(Row = c("A", "B"), Col = c("X", "Y"))))
}

# Correlation Test
generate_correlation_data <- function(
    n = 30,               # Sample size
    rho = 0.5,            # True correlation (effect size)
    dist = "Normal"        # Bivariate distribution
) {
  # Generate two independent standardized variables
  z1 <- generate_data(n, dist)
  z2 <- generate_data(n, dist)
  
  # Induce correlation
  x <- z1
  y <- rho * z1 + sqrt(1 - rho^2) * z2
  
  return(data.frame(x = x, y = y))
}

# Mann-Whitney U Test
generate_mannwhitney_data <- function(
    n1 = 20,              # Sample size group 1
    n2 = 20,              # Sample size group 2
    shift = 0.5,          # Location shift (effect size)
    dist = "Normal"       # Base distribution
) {
  # Generate standardized data and apply shift
  group1 <- generate_data(n1, dist)
  group2 <- shift + generate_data(n2, dist)
  
  return(data.frame(
    group = factor(rep(c("A", "B"), c(n1, n2))),
    value = c(group1, group2)
  ))
}

# Kruskal-Wallis Test
generate_kruskal_data <- function(
    n_per_group = 20,     # Sample size per group
    shifts = c(0, 0.3, 0.6), # Location shifts for groups
    dist = "Normal"       # Base distribution
) {
  k <- length(shifts)
  values <- vector("list", k)
  
  for (i in 1:k) {
    # Generate standardized data and apply shift
    values[[i]] <- shifts[i] + generate_data(n_per_group, dist)
  }
  
  return(data.frame(
    group = factor(rep(LETTERS[1:k], each = n_per_group)),
    value = unlist(values)
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

# -------- example ---------
# Generate data for independent t-test with Gamma distribution
ttest_data <- generate_ind_ttest_data(
  n1 = 30, n2 = 30, mean1 = 0, mean2 = 0.8, sd1 = 1, sd2 = 1, dist = "Gamma"
)

# Generate data for correlation test with Laplace distribution
corr_data <- generate_correlation_data(n = 50, rho = 0.7, dist = "Laplace")

# Generate data for Kruskal-Wallis with Weibull distribution
kruskal_data <- generate_kruskal_data(
  n_per_group = 25, shifts = c(0, 0.5, 1.0), dist = "Weibull"
)
