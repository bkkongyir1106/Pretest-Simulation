## Data Splitting
# # --------------------Addressing Selection Effects ------------------
# # -------------------- Data Splitting -------------------------------
# # -------------------------------------------------------------------
# set.seed(123)
# dist_sum <- c("Normal", "Exponential", "LogNormal")
# sample_size <- c(8, 10, 15, 20, 25, 30)
# N <- 1e3
# alpha <- 0.05
# 
# type_I_error <- matrix(nrow = length(sample_size), ncol = length(dist_sum))
# 
# rownames(type_I_error) <- sample_size
# colnames(type_I_error) <- dist_sum
# 
# for (j in seq_along(dist_sum)) {
#   for (i in seq_along(sample_size)) {
#     dist <- dist_sum[j]
#     n <- sample_size[i]
# 
#     iter <- 0
#     pvals <- numeric(N)
#     count <- 0
# 
#     while (count < N) {
#       iter <- iter + 1
#       x <- generate_data(n, dist)
# 
#       if (shapiro.test(x[1:floor(n/2)])$p.value > alpha) {
#         count <- count + 1
#         pvals[count] <- t.test(x[(floor(n/2) + 1):n])$p.value
#       }
# 
#     }
# 
#     type_I_error[i, j] <- mean(pvals < alpha)
#   }
# }
# 
# # Print results
# print("Type I Error Rates:")
# print(type_I_error)


### Data splitting with two-stage in One Sample Tests
# -------------------- SETUP -------------------------------
set.seed(123)
distributions <- c("Normal", "Exponential", "Chi-Square", "LogNormal")
sample_size <- c(8, 10, 15, 20, 25, 30)
N_sim <- 1e3
alpha <- 0.05

# Initialize results matrices
initialize_matrix <- function(sample_size, distributions) {
  matrix(nrow = length(sample_size), ncol = length(distributions),
         dimnames = list(sample_size, distributions))
}

type_I_error_data_splitting <- initialize_matrix(sample_size, distributions)
type_I_error_two_stage <- initialize_matrix(sample_size, distributions)
bootstrap_type_I_error <- initialize_matrix(sample_size, distributions)

# -------------------- bootstrap function -------------------------------
bootstrap_one_sample <- function(x, effect_size, alpha, n_bootstrap) {
  t_observed <- (mean(x) - effect_size)/(sd(x)/sqrt(length(x)))
  x_centered <- x - mean(x) + effect_size
  
  boot_stats <- replicate(n_bootstrap, {
    boot_sample <- sample(x_centered, size = length(x), replace = TRUE)
    (mean(boot_sample) - effect_size)/(sd(boot_sample)/sqrt(length(x)))
  })
  
  mean(abs(boot_stats) >= abs(t_observed))
}

simulate_method <- function(n, dist, N_sim, alpha, method) {
  foreach(k = 1:N_sim, .combine = rbind, .packages = "stats") %dopar% {
    x <- generate_data(n, dist)
    pval_data_split <- NA
    pval_two_stage <- NA
    pval_boot <- NA
    
    # Data Splitting Procedure
    if(method == "data_split" || method == "all") {
      split_point <- floor(n/2)
      part1 <- x[1:split_point]
      part2 <- x[(split_point+1):n]
      
      if(shapiro.test(part1)$p.value > alpha) {
        pval_data_split <- t.test(part2, mu = 0)$p.value
      } else {
        pval_data_split <- bootstrap_one_sample(part2, 0, alpha, 1e3)
      }
    }
    
    # Two-Stage Procedure
    if(method == "two_stage" || method == "all") {
      if(shapiro.test(x)$p.value > alpha) {
        pval_two_stage <- t.test(x, mu = 0)$p.value
      } else {
        pval_two_stage <- bootstrap_one_sample(x, 0, alpha, 1e3)
      }
    }
    
    # Pure Bootstrap
    if(method == "bootstrap" || method == "all") {
      pval_boot <- bootstrap_one_sample(x, 0, alpha, 1e3)
    }
    
    c(data_split = pval_data_split, 
      two_stage = pval_two_stage, 
      bootstrap = pval_boot)
  }
}

# -------------------- PARALLEL EXECUTION -------------------------------
run_analysis <- function() {
  # Set up parallel backend
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  clusterExport(cl, c("generate_data", "bootstrap_one_sample"))
  
  tryCatch({
    for (j in seq_along(distributions)) {
      for (i in seq_along(sample_size)) {
        n <- sample_size[i]
        dist <- distributions[j]
        
        # Run all methods simultaneously
        results <- simulate_method(n, dist, N_sim, alpha, "all")
        
        # Update all matrices
        if("data_split" %in% colnames(results)) {
          type_I_error_data_splitting[i,j] <- mean(results[,"data_split"] < alpha, na.rm = TRUE)
        }
        if("two_stage" %in% colnames(results)) {
          type_I_error_two_stage[i,j] <- mean(results[,"two_stage"] < alpha, na.rm = TRUE)
        }
        if("bootstrap" %in% colnames(results)) {
          bootstrap_type_I_error[i,j] <- mean(results[,"bootstrap"] < alpha, na.rm = TRUE)
        }
        
        cat(sprintf("Completed: %s, n=%d\n", dist, n))
      }
    }
  }, finally = {
    stopCluster(cl)
  })
  
  list(
    data_splitting = type_I_error_data_splitting,
    two_stage = type_I_error_two_stage,
    bootstrap = bootstrap_type_I_error
  )
}

# -------------------- EXECUTE & SAVE -------------------------------
results <- run_analysis()

# Save results
save(results, file = "onesample_adaptive_datasplit_type_I_errors.RData")

# Print results
print("Data Splitting Results:")
print(results$data_splitting)
print("\nTwo-Stage Results:")
print(results$two_stage)
print("\nBootstrap Results:")
print(results$bootstrap)


# --------------------Addressing Selection Effects ------------------
# ------------------- Conditional Type I error ----------------------
# -------------------------------------------------------------------
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
set.seed(123)
distributions <- c("Normal", "Exponential", "Chi-Square", "LogNormal")
sample_size <- c(8, 10, 15, 20, 25, 30)
N_sim <- 1e3
alpha <- 0.05
R <- 1e3

# Initialize results matrices
initialize_matrix <- function(sample_size, distributions) {
  matrix(nrow = length(sample_size), ncol = length(distributions),
         dimnames = list(sample_size, distributions))
}

one_sample.data_split_boot_error <- initialize_matrix(sample_size, distributions)
one_sample.data_split_perm_error <- initialize_matrix(sample_size, distributions)
one_sample.adaptive_boot_error <- initialize_matrix(sample_size, distributions)
one_sample.adaptive_perm_error <- initialize_matrix(sample_size, distributions)
one_sample_boot_error <- initialize_matrix(sample_size, distributions)
one_sample_perm_error <- initialize_matrix(sample_size, distributions)

# ---------------------------- bootstrap function -----------------------------------
bootstrap_one_sample <- function(x, effect_size = 0, alpha = 0.05, n_bootstrap = R) {
  t_obs <- (mean(x) - effect_size) / (sd(x)/sqrt(length(x)))
  x_centered <- x - mean(x) + effect_size
  
  null_dist <- replicate(n_bootstrap, {
    bs_sample <- sample(x_centered, replace = TRUE)
    (mean(bs_sample) - effect_size) / (sd(bs_sample)/sqrt(length(x)))
  })
  
  mean(abs(null_dist) >= abs(t_obs))
}

# Parallel simulation function
run_simulations <- function(method) {
  cl <- makeCluster(detectCores() - 1)
  registerDoParallel(cl)
  clusterExport(cl, c("generate_data", "bootstrap_one_sample", "OneSample.test", "OneSample_test_statistic", "alpha", "R"))
  
  result_matrix <- initialize_matrix(sample_size, distributions)
  
  for (j in seq_along(distributions)) {
    for (i in seq_along(sample_size)) {
      n <- sample_size[i]
      dist <- distributions[j]
      
      pvals <- foreach(k = 1:N_sim, .combine = c, .packages = "stats") %dopar% {
        x <- generate_data(n, dist)
        
        switch(method,
               "data_split_boot" = {
                 split_point <- floor(n/2)
                 part1 <- x[1:split_point]
                 part2 <- x[(split_point+1):n]
                 
                 if (shapiro.test(part1)$p.value > alpha) {
                   t.test(part2, mu = 0)$p.value
                 } else {
                   bootstrap_one_sample(part2)
                 }
               },
               "data_split_perm" = {
                 split_point <- floor(n/2)
                 part1 <- x[1:split_point]
                 part2 <- x[(split_point+1):n]
                 
                 if (shapiro.test(part1)$p.value > alpha) {
                   t.test(part2, mu = 0)$p.value
                 } else {
                   OneSample.test(part2, test = "perm", alpha = 0.05, effect_size = 0.0, B = R)
                 }
               },
               "two_stage_boot" = {
                 if (shapiro.test(x)$p.value > alpha) {
                   t.test(x, mu = 0)$p.value
                 } else {
                   bootstrap_one_sample(x)
                 }
               },
               "two_stage_perm" = {
                 if (shapiro.test(x)$p.value > alpha) {
                   t.test(x, mu = 0)$p.value
                 } else {
                   OneSample.test(x, test = "perm", alpha = 0.05, effect_size = 0.0, B = R)
                 }
               },
               "bootstrap" = {
                 bootstrap_one_sample(x)
               },
      
               "perm" = {
                 OneSample.test(x, "perm", alpha = 0.05, effect_size = 0.0, B = R)
               }
        )
      }
      
      result_matrix[i, j] <- mean(pvals < alpha)
      cat(sprintf("%s: %s, n=%d completed\n", method, dist, n))
    }
  }
  
  stopCluster(cl)
  return(result_matrix)
}

# ---------------------------- EXECUTION ----------------------------------------
# Run all methods in parallel
one_sample.data_split_boot_error <- run_simulations("data_split_boot")
one_sample.data_split_perm_error <- run_simulations("data_split_perm")
one_sample.adaptive_boot_error <- run_simulations("two_stage_boot")
one_sample.adaptive_perm_error <- run_simulations("two_stage_perm")
one_sample_boot_error <- run_simulations("bootstrap")
one_sample_perm_error <- run_simulations("perm")

# ---------------------------- SAVE RESULTS -------------------------------------
save(one_sample.data_split_boot_error,  one_sample.data_split_perm_error, one_sample.adaptive_boot_error,
one_sample.adaptive_perm_error, one_sample_boot_error, one_sample_perm_error, file = "one_sample.adaptive_type_I_error_results.RData")

# ======================================= Two Sample Case ======================
### Addressing Selection Effects
# # --------------------Addressing Selection Effects ------------------
# # -------------------- Conditional Type I error -------------------------------
# # -------------------------------------------------------------------
# # Set simulation parameters
# set.seed(12345)
# distributions <- c("Normal","Exponential","Chi-Square", "LogNormal")
# sample_size <- c(8, 10, 15, 20, 25, 30)
# N_sim <- 1e4
# alpha <- 0.05
#
# # Initialize matrices to store results
# initialize_matrix <- function(sample_size, distributions) {
#   matrix(nrow = length(sample_size), ncol = length(distributions),
#          dimnames = list(sample_size, distributions))
# }
# twosample.type_I_error_data_splitting <- initialize_matrix(sample_size, distributions)
# twosample.type_I_error_two_stage <- initialize_matrix(sample_size, distributions)
#
# # -------------------- Data Splitting Procedure -------------------------
# for (j in seq_along(distributions)) {
#   for (i in seq_along(sample_size)) {
#     dist <- distributions[j]
#     n <- sample_size[i]
#     pvals <- numeric(N_sim)
#     count <- 0
#     iter <- 0
#
#     while (count < N_sim) {
#       iter <- iter + 1
#       x <- generate_data(n, dist)
#       y <- generate_data(n, dist)
#       split_point <- floor(n/2)
#
#       # First stage: Shapiro test on first half
#       if (shapiro.test(x[1:split_point])$p.value > alpha & shapiro.test(y[1:split_point])$p.value > alpha) {
#         count <- count + 1
#         # Second stage: t-test on second half
#         pvals[count] <- t.test(x[(split_point+1):n],y[(split_point+1):n] )$p.value
#       }
#     }
#     twosample.type_I_error_data_splitting[i, j] <- mean(pvals < alpha)
#   }
# }
#
# # -------------------- Two-Stage Procedure ------------------------------
# for (j in seq_along(distributions)) {
#   for (i in seq_along(sample_size)) {
#     dist <- distributions[j]
#     n <- sample_size[i]
#     pvals <- numeric(N_sim)
#     count <- 0
#     iter <- 0
#
#     while (count < N_sim) {
#       iter <- iter + 1
#       x <- generate_data(n, dist)
#       y <- generate_data(n, dist)
#       # First stage: Shapiro test on full data
#       if (shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha) {
#         count <- count + 1
#         # Second stage: t-test on full data
#         pvals[count] <- t.test(x, y)$p.value
#       }
#     }
#     twosample.type_I_error_two_stage[i, j] <- mean(pvals < alpha)
#   }
# }
#
# save(twosample.type_I_error_data_splitting,twosample.type_I_error_two_stage, file = "twosample.data_splitting.RData")

# ==============================================================================
effect = 0.5
N = 1000
R = 1000
x = rnorm(100)
pval = c()
for(b in 1 : N){
  pval[b] = bootstrap_one_sample(x , effect_size = effect, B = R)
}
power = mean(pval < 0.05)
power

# ==========================
# --------------------------
# Bootstrap Test Function
# --------------------------
bootstrap_test <- function(data, mu0, n_boot = 1000, alpha = 0.05, two_tailed = TRUE) {
  # Perform a bootstrap one-sample location test
  #
  # Args:
  #   data: Numeric vector of sample data
  #   mu0: Null hypothesis mean
  #   n_boot: Number of bootstrap replicates
  #   alpha: Significance level
  #   two_tailed: Whether to perform a two-tailed test
  #
  # Returns:
  #   List containing p-value and rejection decision
  
  n <- length(data)
  observed_stat <- mean(data) - mu0
  
  # Shift data to enforce H0: mu = mu0
  shifted_data <- data - mean(data) + mu0
  
  # Generate bootstrap distribution under H0
  boot_means <- replicate(n_boot, {
    mean(sample(shifted_data, size = length(data), replace = TRUE))
  })
  
  # Calculate p-value
  if (two_tailed) {
    p_value <- mean(abs(boot_means - mu0) >= abs(observed_stat))
  } else {
    # One-tailed (right-tail) test
    p_value <- mean((boot_means - mu0) >= observed_stat)
  }
  
  list(
    reject = (p_value < alpha),
    p_value = p_value
  )
}

# --------------------------
# Simulation Function
# --------------------------
simulate_power <- function(n_sim, n_sample, mu0, mu_true, alpha = 0.05, ...) {
  # Simulate Type I error or power for the bootstrap test
  #
  # Args:
  #   n_sim: Number of simulation iterations
  #   n_sample: Sample size
  #   mu0: Null hypothesis mean
  #   mu_true: True population mean (for data generation)
  #   alpha: Significance level
  #   ...: Additional arguments passed to bootstrap_test (e.g., n_boot)
  #
  # Returns:
  #   Rejection rate (Type I error or power)
  
  rejections <- c()
  
  for (i in 1:n_sim) {
    # Generate data under the specified true mean
    data <- rnorm(n_sample, mean = mu_true)
    data <- rnorm(n_sample, mean = mu_true)
    # Run bootstrap test
    test_result <- bootstrap_test(data, mu0 = mu0, alpha = alpha, ...)
    
    rejections[i] <- test_result$reject
  }
  
  mean(rejections)  # Return rejection rate
}

# --------------------------
# Example Usage
# --------------------------
set.seed(123)

# Type I Error (mu_true = mu0 = 0)
type_I_error <- simulate_power(
  n_sim = 1000,
  n_sample = 30,
  mu0 = 0,
  mu_true = 0  # Simulate under H0
)

# Power (mu_true = 0.5)
power <- simulate_power(
  n_sim = 1000,
  n_sample = 30,
  mu0 = 0,
  mu_true = 0.5  # Simulate under H1
)

cat("Type I error rate:", type_I_error, "\n")
cat("Power:", power)


# ============== using lapply functiomn =============================
effect_sizes <- seq(0, 1, by = 0.1)
power_results <- sapply(effect_sizes, function(mu1) {
  simulate_power(n_sim = 1000, n_sample = 30, mu0 = 0, mu_true = mu1)
})

sample_sizes <- c(20, 30, 50, 100)
power_by_n <- sapply(sample_sizes, function(n) {
  simulate_power(n_sim = 1000, n_sample = n, mu0 = 0, mu_true = 0.5)
})
 
# ====================== different sample sizes ====================
# --------------------------
# Bootstrap Test Function (unchanged)
# --------------------------
bootstrap_test <- function(data, mu0, n_boot = 10000, alpha = 0.05, two_tailed = TRUE) {
  n <- length(data)
  observed_stat <- mean(data) - mu0
  
  # Shift data to enforce H0: mu = mu0
  shifted_data <- data - mean(data) + mu0
  
  # Generate bootstrap distribution under H0
  boot_means <- replicate(n_boot, {
    mean(sample(shifted_data, size = n, replace = TRUE))
  })
  
  # Calculate p-value
  if (two_tailed) {
    p_value <- mean(abs(boot_means - mu0) >= abs(observed_stat))
  } else {
    p_value <- mean((boot_means - mu0) >= observed_stat)
  }
  
  list(reject = (p_value < alpha), p_value = p_value)
}

# --------------------------
# Simulation Function (modified for exponential data)
# --------------------------
simulate_power <- function(n_sim, n_sample, mu0, mu_true, alpha = 0.05, ...) {
  rejections <- logical(n_sim)
  
  for (i in 1:n_sim) {
    # Generate exponential data with mean = mu_true
    data <- rexp(n_sample, rate = 1 / mu_true)
    
    # Run bootstrap test
    test_result <- bootstrap_test(data, mu0 = mu0, alpha = alpha, ...)
    rejections[i] <- test_result$reject
  }
  
  mean(rejections)  # Return rejection rate (power)
}

# --------------------------
# Example: Power vs. Sample Size for Effect Size = 0.5
# --------------------------
set.seed(123)
mu0 <- 1            # Null hypothesis mean (H0: mu = 1)
effect_size <- 0.5  # True mean is mu0 + 0.5 = 1.5 (H1)
sample_sizes <- c(20, 30, 50, 100, 200)  # Sample sizes to test

# Simulate power for each sample size
power_results <- sapply(sample_sizes, function(n) {
  simulate_power(
    n_sim = 1000,
    n_sample = n,
    mu0 = mu0,
    mu_true = mu0 + effect_size  # True mean under H1
  )
})

# Create results table
results <- data.frame(
  Sample_Size = sample_sizes,
  Power = power_results
)

print(results)

# ========================
# --------------------------
# Bootstrap Test Function
# --------------------------
bootstrap_test <- function(data, mu0, n_boot = 10000, alpha = 0.05, two_tailed = TRUE) {
  n <- length(data)
  observed_stat <- mean(data) - mu0
  
  # Shift data to enforce H0: mu = mu0
  shifted_data <- data - mean(data) + mu0
  
  # Generate bootstrap distribution under H0
  boot_means <- replicate(n_boot, {
    mean(sample(shifted_data, size = n, replace = TRUE))
  })
  
  # Calculate p-value
  if (two_tailed) {
    p_value <- mean(abs(boot_means - mu0) >= abs(observed_stat))
  } else {
    p_value <- mean((boot_means - mu0) >= observed_stat)
  }
  
  list(reject = (p_value < alpha), p_value = p_value)
}

# --------------------------
# Unified Simulation Function
# --------------------------
simulate_power <- function(
    n_sim, 
    n_sample, 
    mu0, 
    effect_size, 
    alpha = 0.05, 
    distribution = "normal",
    ...
) {
  # Simulate power for a given distribution and sample size
  #
  # Args:
  #   n_sim: Number of simulation iterations
  #   n_sample: Sample size
  #   mu0: Null hypothesis mean
  #   effect_size: True difference from mu0 (H1: mu = mu0 + effect_size)
  #   distribution: "normal", "exponential", "chi-squared", "lognormal"
  #   ...: Additional arguments passed to bootstrap_test
  #
  # Returns:
  #   Power estimate and metadata
  
  rejections <- logical(n_sim)
  mu_true <- mu0 + effect_size  # True population mean
  
  for (i in 1:n_sim) {
    # Generate data under H1 for the specified distribution
    if (distribution == "normal") {
      data <- rnorm(n_sample, mean = mu_true)
    } else if (distribution == "exponential") {
      data <- rexp(n_sample, rate = 1 / mu_true)
    } else if (distribution == "chi-squared") {
      data <- rchisq(n_sample, df = mu_true)
    } else if (distribution == "lognormal") {
      # Adjust lognormal parameters to achieve mean = mu_true
      sigma <- 1  # Fixed shape parameter (sdlog)
      mu <- log(mu_true) - sigma^2 / 2
      data <- rlnorm(n_sample, meanlog = mu, sdlog = sigma)
    } else {
      stop("Invalid distribution specified")
    }
    
    # Run bootstrap test
    test_result <- bootstrap_test(data, mu0 = mu0, alpha = alpha, ...)
    rejections[i] <- test_result$reject
  }
  
  list(
    distribution = distribution,
    sample_size = n_sample,
    power = mean(rejections),
    mu0 = mu0,
    effect_size = effect_size
  )
}

# --------------------------
# Example: Power Across Distributions and Sample Sizes
# --------------------------
set.seed(123)
mu0 <- 1          # Null hypothesis mean
effect_size <- 0.5  # True effect size (H1: mu = 1.5)
sample_sizes <- c(20, 30, 50, 100)
distributions <- c("normal", "exponential", "chi-squared", "lognormal")

# Simulate power for all combinations
results <- list()
for (dist in distributions) {
  for (n in sample_sizes) {
    result <- simulate_power(
      n_sim = 1000,
      n_sample = n,
      mu0 = mu0,
      effect_size = effect_size,
      distribution = dist
    )
    results <- append(results, list(result))
  }
}

# Convert results to a data frame
results_df <- do.call(rbind, lapply(results, function(x) {
  data.frame(
    Distribution = x$distribution,
    Sample_Size = x$sample_size,
    Power = round(x$power, 3)
  )
}))

print(results_df)


# ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
# --------------------------
# One-Sample Bootstrap Test
# --------------------------
bootstrap_one_sample <- function(data, effect_size = 0, B = 1000, two_tailed = TRUE) {
  # Perform a bootstrap test for H0: mean = (observed mean - effect_size)
  #
  # Args:
  #   data: Numeric vector of sample data
  #   effect_size: Difference to test under H0 (default = 0)
  #   B: Number of bootstrap replicates
  #   two_tailed: Whether to perform a two-tailed test
  #
  # Returns:
  #   P-value for the test
  
  n <- length(data)
  observed_mean <- mean(data)
  mu0 <- observed_mean - effect_size  # Null hypothesis mean
  
  # Shift data to enforce H0: mu = mu0
  shifted_data <- data - observed_mean + mu0
  
  # Generate bootstrap distribution under H0
  boot_means <- replicate(B, {
    boot_sample <- sample(shifted_data, size = n, replace = TRUE)
    mean(boot_sample)
  })
  
  # Calculate p-value
  observed_stat <- observed_mean - mu0  # = effect_size
  if (two_tailed) {
    p_value <- mean(abs(boot_means - mu0) >= abs(observed_stat))
  } else {
    p_value <- mean((boot_means - mu0) >= observed_stat)
  }
  
  p_value
}


# --------------------------
# Generic Simulation Wrapper
# --------------------------
simulate_test <- function(
    n_sim = 1000, 
    n_sample = 30, 
    mu_true = 0, 
    effect_size = 0, 
    distribution = "normal",
    B = 1000
) {
  # Simulate Type I error (effect_size = 0) or power (effect_size ≠ 0)
  #
  # Args:
  #   n_sim: Number of simulation iterations
  #   n_sample: Sample size
  #   mu_true: True population mean
  #   effect_size: Null hypothesis effect size (H0: mean = mu_true - effect_size)
  #   distribution: "normal", "exponential", "chi-squared", "lognormal"
  #   B: Bootstrap replicates
  #
  # Returns:
  #   Rejection rate (Type I error or power)
  
  rejections <- logical(n_sim)
  
  for (i in 1:n_sim) {
    # Generate data from the specified distribution
    if (distribution == "normal") {
      data <- rnorm(n_sample, mean = mu_true)
    } else if (distribution == "exponential") {
      data <- rexp(n_sample, rate = 1 / mu_true)
    } else if (distribution == "chi-squared") {
      data <- rchisq(n_sample, df = mu_true)
    } else if (distribution == "lognormal") {
      sigma <- 1  # Fixed shape parameter
      mu <- log(mu_true) - sigma^2 / 2
      data <- rlnorm(n_sample, meanlog = mu, sdlog = sigma)
    } else {
      stop("Invalid distribution")
    }
    
    # Run bootstrap test
    p_value <- bootstrap_one_sample(data, effect_size = effect_size, B = B)
    rejections[i] <- (p_value < 0.05)  # Using alpha = 0.05
  }
  
  mean(rejections)
}


# --------------------------
# Example 1: Type I Error (Normal Distribution)
# --------------------------
set.seed(123)
type_I_error <- simulate_test(
  n_sim = 1000,
  mu_true = 0,      # True population mean
  effect_size = 0,  # Test H0: mean = 0
  distribution = "exponential"
)
cat("Type I Error (exponential):", type_I_error, "\n")

# --------------------------
# Example 2: Power (Exponential Distribution)
# --------------------------
power <- simulate_test(
  n_sim = 1000,
  mu_true = 1.5,    # True population mean
  effect_size = 0.5, # Test H0: mean = 1.5 - 0.5 = 1.0
  distribution = "exponential"
)
cat("Power (Exponential):", power, "\n")



# ======================= Bootstrap ===========================================
# ----------------------One sample bootstrap --------------------------------
bootstrap_one_sample <- function(data, effect_size, B ) {
  # Handle missing values
  data <- na.omit(data)
  n <- length(data)
  
  # Observed test statistic
  observed_stat <- mean(data) - effect_size
  
  # Shift data to enforce H0: mu = mu0
  shifted_data <- data - mean(data) + effect_size
  
  # Generate bootstrap distribution under H0
  boot_means <- replicate(B, {
    mean(sample(shifted_data, size = length(data), replace = TRUE))
  })
  
  # Calculate p-value
  pval <- mean(abs(boot_means - effect_size) >= abs(observed_stat))
  
  return(pval)
}

# ---------------------two sample bootstrap ---------------------------
bootstrap_two_sample <- function(data1, data2, effect_size , B) {
  
  # Observed test statistic (difference in means)
  observed_stat <- mean(data1) - mean(data2) - effect_size
  
  # Shift data to enforce H0: mu1 - mu2 = effect_size
  shifted_data1 <- data1 - mean(data1) + (mean(data2) + effect_size)
  shifted_data2 <- data2 - mean(data2) + mean(data2)  # No shift for group 2
  
  # Generate bootstrap distribution under H0
  boot_diffs <- replicate(B, {
    boot_sample1 <- sample(shifted_data1, size = length(data1), replace = TRUE)
    boot_sample2 <- sample(shifted_data2, size = length(data2), replace = TRUE)
    mean(boot_sample1) - mean(boot_sample2)
  })
  
  # Calculate p-value
  pval <- mean(abs(boot_diffs - effect_size) >= abs(observed_stat))
  
  return(pval)
}


# ----------------------One sample bootstrap --------------------------------

bootstrap_one_sample_a <- function(x, effect_size, B) {
  # Calculate observed test statistic
  t_observed <- OneSample_test_statistic(x)
  
  # Generate bootstrap statistics under H0
  boot_stats <- replicate(B, {
    boot_sample <- sample(x + effect_size, size = length(x), replace = TRUE) 
    OneSample_test_statistic(boot_sample)
  })
  
  # Compute two-tailed p-value
  p_value <- mean(abs(boot_stats) >= abs(t_observed))
  
  return(p_value)
}

# ------------------------------Two sample bootstrap -----------------------
bootstrap_two_sample_a <- function(x, y, effect_size, B) {
  # Calculate observed test statistic
  t_observed <- TwoSample_test_statistic(x, y + effect_size)
  
  # Generate bootstrap statistics under H0
  boot_stats <- replicate(B, {
    x_sample <- sample(x, size = length(x), replace = TRUE)
    y_sample <- sample(y, size = length(y), replace = TRUE)
    TwoSample_test_statistic(x_sample, y_sample + effect_size)
  })
  
  # Compute two-tailed p-value
  p_value <- mean(abs(boot_stats) >= abs(t_observed))
  
  return(p_value)
}
# ===========================================================================
#                            bootstrap p-value function
# ==========================================================================
bootstrap_one_sample_pvalue_CI <- function(x, null_effect, n_bootstrap) {
  x_centered <- x - mean(x) + null_effect
  T_obs <- OneSample_test_statistic(x, effect_size = null_effect)
  null_stats <- replicate(n_bootstrap, {
    resampled_x <- sample(x_centered, replace = TRUE)
    OneSample_test_statistic(resampled_x, effect_size = null_effect)
  })
  p_value <- mean(abs(null_stats) >= abs(T_obs))
  return(p_value)
}

bootstrap_one_sample_pvalue_CI(rnorm(10), null_effect = 0, n_bootstrap = 1000)
# ----------------------------------------------------
#      Simulate Type I Error (H₀ is true)
# ---------------------------------------------------
simulate_typeI_error <- function(null_effect, alpha, n, n_sim, n_bootstrap, dist) {
  rejections <- replicate(n_sim, {
    # Generate data under H₀ (already centered at 0)
    data <- generate_data(n, dist)
    # Compute p-value
    p_value <- bootstrap_one_sample_pvalue_CI(data, null_effect, n_bootstrap)
    # Check rejection
    p_value < alpha
  })
  mean(rejections)  # Type I Error Rate
}

# ----------------------------------------------------
#      Simulate Power (H₁ is true)
# ---------------------------------------------------
simulate_power <- function(null_effect, true_effect, alpha, n, n_sim, n_bootstrap, dist) {
  rejections <- replicate(n_sim, {
    # Generate normalized data under H₀, then shift to H₁
    data <- generate_data(n, dist) + true_effect
    # Compute p-value
    p_value <- bootstrap_one_sample_pvalue_CI(data, null_effect, n_bootstrap)
    # Check rejection
    p_value < alpha
  })
  mean(rejections)  # Power
}

# Example 1: Type I Error for Normal Data (Target α = 0.05)
set.seed(123)
typeI_error_normal <- simulate_typeI_error(
  null_effect = 0,
  alpha = 0.05,
  n = 10,
  n_sim = 1000,
  n_bootstrap = 1000,
  dist = "Exponential"
)
print(paste("Type I Error (Normal):", round(typeI_error_normal, 3)))

# Example 2: Power for Normal Data (True Effect = 0.5)
set.seed(123)
power_normal <- simulate_power(
  null_effect = 0,
  true_effect = 0.0,
  alpha = 0.05,
  n = 10,
  n_sim = 1000,
  n_bootstrap = 1000,
  dist = "LogNormal"
)
print(paste("Power (Normal):", round(power_normal, 3)))

# Example 3: Type I Error for Heavy-Tailed Data (e.g., Laplace)
typeI_error_laplace <- simulate_typeI_error(
  null_effect = 0.0,
  alpha = 0.05,
  n = 10,
  n_sim = 1000,
  n_bootstrap = 1000,
  dist = "LogNormal"
)
print(paste("Type I Error (Laplace):", round(typeI_error_laplace, 3)))


# ========================================================================
# ----------------------------- Bootstrapping -----------------------------------
# test = c("t", "Wilcoxon","t-Wilcoxn", "perm" ), B_perm= number of permutations 
# -------------------------------------------------------------------------------
{
  # bootstrap_one_sample_test <- function(x, test, effect_size, alpha, n_bootstrap, sample_size) {
  #   purrr::map_dbl(sample_size, ~{
  #     mean(replicate(n_bootstrap, {
  #       sample_data <- sample(x, .x, replace = TRUE) + effect_size
  #       OneSample.test(sample_data, test, alpha, B = 100) < alpha
  #     }))
  #   })
  # }
  # 
  # bootstrap_two_sample_test <- function(x, y, test, effect_size, alpha, n_bootstrap, sample_size) {
  #   purrr::map_dbl(sample_size, ~{
  #     mean(replicate(n_bootstrap, {
  #       sample_x <- sample(x, .x, replace = TRUE)
  #       sample_y <- sample(y, .x, replace = TRUE) + effect_size
  #       TwoSample.test(sample_x, sample_y, test, alpha, B = 100) < alpha
  #     }))
  #   })
  # }
}

