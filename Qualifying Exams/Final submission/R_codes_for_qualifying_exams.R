# --------------------------------------------------------
#   data generation function
# --------------------------------------------------------
generate_data <- function(n, dist) {
  dist <- tolower(dist)
  if (dist == "normal") {
    x <- rnorm(n, mean = 1, sd = 1)
  } else if (dist == "exponential") {
    x <- rexp(n, rate = 1) 
  }
  return(x)
}
# --------------------------------------------------------
# a) R function for adaptive one-sample location test   #
# ---------------------------------------------------------
adaptive_test <- function(x, alpha = 0.05, train_frac = 0.10){
  # randomly choose training set
  train_index <- sample(1:length(x), floor(train_frac * length(x)), replace = FALSE)
  # Training data for normality check
  x_train <- x[train_index]
  # Test data for location test
  x_test <- x[-train_index]
  # SW-test on training data
  if(shapiro.test(x_train)$p.value < alpha){ 
    # Non-normal: Use Wilcoxon sign rank test
    pval <- wilcox.test(x_test, mu = 1)$p.value 
  }else{# Normal: use t-test
    pval <- t.test(x_test, mu = 1)$p.value 
  }
  return(pval)
}

# ----------------------------------------------
# b)    10-fold cross-validated adaptive test  #
# -----------------------------------------------
adaptive_test_cv <- function(x, alpha = 0.05, K = 10){
  # Create K folds
  folds <- sample(rep(1:K, length.out = length(x)))
  # Store p-values from each fold
  pvals_fold <- numeric(K)
  for(k in 1:K){
    # Current test fold
    test_index <- which(folds ==k)
    # Remaining as training
    train_index <- setdiff(1:length(x), test_index)
    # Define Training and test sets
    x_train <- x[train_index]
    x_test <- x[test_index]
    if(shapiro.test(x_train)$p.value < alpha){ # Non-normal
      pvals_fold[k] <- wilcox.test(x_test, mu = 1)$p.value
    }else{# Normal
      pvals_fold[k] <- t.test(x_test, mu = 1)$p.value
    }
  }
  # ## Combine the K p-values with Minimum P-value method
  p_min <- min(pvals_fold)
  p_combined <- 1 - (1 - p_min)^K
  return(p_combined)
}

# -----------------------------------------------------------
# d) Implementation:Simulation study for Exp(1) and N(1,1)
# -----------------------------------------------------------
# Define run controls
set.seed(123)
Nsim         <- 1e4 # Number of simulations
alpha        <- 0.05 # Significance level
sample_size  <- c(30, 50, 100) # Sample sizes to test
distributions <- c("normal", "exponential") # Distributions to test
# Initialize result matrices
type1_adaptive    <- matrix(NA,
                            nrow = length(sample_size),
                            ncol = length(distributions),
                            dimnames = list(sample_size, distributions))
type1_adaptive_cv <- type1_adaptive  

# loop through each sample size and each distribution
for (i in seq_along(sample_size)) {
  n <- sample_size[i]  
  for (j in seq_along(distributions)) {
    dist <- distributions[j]  
    # store all p-vales
    pvals1 <- numeric(Nsim)  
    pvals2 <- numeric(Nsim) 
    # calculate p-values for each procedure
    for (s in 1:Nsim) {
      x <- generate_data(n, dist = dist)  
      pvals1[s] <- adaptive_test(x, alpha = alpha, train_frac = 0.10)
      pvals2[s] <- adaptive_test_cv(x, alpha = alpha, K = 10)
    }
    # calculate type I error rates
    type1_adaptive[i, j]  <- mean(pvals1 < alpha)  
    type1_adaptive_cv[i, j] <- mean(pvals2 < alpha)
  }
}

# print results
type1_adaptive        
type1_adaptive_cv     

# -------------------------------------------------
# 3 compare wilcoxon test to adaptive procedures 
# ------------------------------------------------
set.seed(123)
type1_wilcox <- matrix(NA, 
                       nrow = length(sample_size), 
                       ncol = length(distributions),
                       dimnames = list(sample_size, distributions))

for (i in seq_along(sample_size)) {
  n <- sample_size[i]
  for (j in seq_along(distributions)) {
    dist <- distributions[j]
    pvals3 <- numeric(Nsim)
    for (s in 1:Nsim) {
      x <- generate_data(n, dist)
      pvals3[s] <- wilcox.test(x, mu = 1)$p.value
    }
    type1_wilcox[i, j] <- mean(pvals3 < alpha)
  }
}

# organize results in tables
compare_methods <- list(
  normal = cbind(
    Adaptive    = type1_adaptive[ , "normal"],
    Adaptive_cv = type1_adaptive_cv[ , "normal"],
    Wilcoxon    = type1_wilcox[ , "normal"]
  ),
  exponential = cbind(
    Adaptive    = type1_adaptive[ , "exponential"],
    Adaptive_cv = type1_adaptive_cv[ , "exponential"],
    Wilcoxon    = type1_wilcox[ , "exponential"]
  )
)
# print results
compare_methods$normal
compare_methods$exponential










