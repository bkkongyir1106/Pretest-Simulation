# Set working directory and source functions
setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/SPRING2025/DOWNSTREAM TESTS/compare test methods_20250310/twosamples/twosamples")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

# Set up simulation parameters
{
  Nsim <- 1e5
  sample_size <- c(10, 15, 20, 25)
  distribution <- "LogNormal"
  alpha <- 0.05
  B <- 1e3
}

# Initialize arrays for Type I error rates
power_sw <- Type_I_error_t.test <- Type_I_error_perm.test <- Type_I_error_adaptive.test <- 
  Type_I_error_conditional.test <- array(NA, dim = c(length(sample_size), length(distribution)), 
                                         dimnames = list(sample_size, distribution))

power_sw_split <- Type_I_error_t.test_split <- Type_I_error_perm.test_split <- Type_I_error_adaptive.test_split <- 
  Type_I_error_conditional.test_split <- array(NA, dim = c(length(sample_size), length(distribution)), 
                                               dimnames = list(sample_size, distribution))

# Initialize lists to store p-values
pvalues_list <- list()
pvalues_split_list <- list()

# Progress bar setup
pb <- txtProgressBar(min = 0, max = length(sample_size) * Nsim, style = 3)
progress_counter <- 0

# Run simulation for each sample size
for (j in seq_along(sample_size)) {
  n <- sample_size[j]
  
  # Initialize vectors to store p-values
  pval_sw_x <- numeric(Nsim)
  pval_sw_y <- numeric(Nsim)
  pval_t <- numeric(Nsim)
  pval_perm <- numeric(Nsim)
  
  pval_sw_split_x <- numeric(Nsim)
  pval_sw_split_y <- numeric(Nsim)
  pval_split_t <- numeric(Nsim)
  pval_split_perm <- numeric(Nsim)
  
  for (k in 1:Nsim) {
    x <- generate_data(n, distribution)
    y <- generate_data(n, distribution)
    
    # Perform Shapiro-Wilk test
    pval_sw_x[k] <- shapiro.test(x)$p.value
    pval_sw_y[k] <- shapiro.test(y)$p.value
    
    # Perform t-test
    pval_t[k] <- t.test(x, y, mu = 0)$p.value
    
    # Perform permutation test
    observe_stat <- (mean(x) - mean(y)) / sqrt(var(x)/length(x) + var(y)/length(y))
    data <- c(x, y)
    permuted_stat <- numeric(B)
    
    for (p in 1:B) {
      sample_data <- sample(data)
      sample_x <- sample_data[1:length(x)]
      sample_y <- sample_data[(length(x) + 1):(length(x) + length(y))]
      permuted_stat[p] <- (mean(sample_x) - mean(sample_y)) / sqrt(var(sample_x)/length(sample_x) + var(sample_y)/length(sample_y))
    }
    pval_perm[k] <- mean(abs(permuted_stat) >= abs(observe_stat))
    
    # Split method
    split_point <- floor(length(x)/2)
    first_half.x <- x[1:split_point]
    second_half.x <- x[(split_point + 1):length(x)]
    first_half.y <- y[1:split_point]
    second_half.y <- y[(split_point + 1):length(y)]
    
    pval_sw_split_x[k] <- shapiro.test(first_half.x)$p.value
    pval_sw_split_y[k] <- shapiro.test(first_half.y)$p.value
    
    # Perform t-test (split)
    pval_split_t[k] <- t.test(second_half.x, second_half.y, mu = 0)$p.value
    
    # Perform permutation test (split)
    observe_stat <- (mean(second_half.x) - mean(second_half.y)) / sqrt(var(second_half.x)/length(second_half.x) + var(first_half.y)/length(first_half.y))
    data <- c(second_half.x, second_half.y)
    permuted_stat <- numeric(B)
    
    for (p in 1:B) {
      sample_data <- sample(data)
      sample_x <- sample_data[1:length(second_half.x)]
      sample_y <- sample_data[(length(second_half.x) + 1):(length(second_half.x) + length(second_half.y))]
      permuted_stat[p] <- (mean(sample_x) - mean(sample_y)) / sqrt(var(sample_x)/length(sample_x) + var(sample_y)/length(sample_y))
    }
    pval_split_perm[k] <- mean(abs(permuted_stat) >= abs(observe_stat))
    
    progress_counter <- progress_counter + 1
    setTxtProgressBar(pb, progress_counter)
  }
  
  # Store p-values
  pvalues_list[[paste0("n_", n)]] <- list(
    pval_sw_x = pval_sw_x,
    pval_sw_y = pval_sw_y,
    pval_t = pval_t,
    pval_perm = pval_perm
  )
  
  pvalues_split_list[[paste0("n_", n)]] <- list(
    pval_sw_split_x = pval_sw_split_x,
    pval_sw_split_y = pval_sw_split_y,
    pval_split_t = pval_split_t,
    pval_split_perm = pval_split_perm
  )
  
  # Compute Type I errors (traditional method)
  power_sw[j, 1] <- mean(pval_sw_x < alpha & pval_sw_y < alpha)
  Type_I_error_t.test[j, 1] <- mean(pval_t < alpha)
  Type_I_error_perm.test[j, 1] <- mean(pval_perm < alpha)
  Type_I_error_adaptive.test[j, 1] <- mean((pval_sw_x > alpha & pval_sw_y > alpha) * (pval_t < alpha)) + 
    mean((pval_sw_x <= alpha & pval_sw_y <= alpha) * (pval_perm < alpha))
  Type_I_error_conditional.test[j, 1] <- round(mean(pval_t[pval_sw_x > alpha & pval_sw_y > alpha] < alpha), 3)
  
  # Compute Type I errors (split method)
  power_sw_split[j, 1] <- mean(pval_sw_split_x < alpha & pval_sw_split_y < alpha)
  Type_I_error_t.test_split[j, 1] <- mean(pval_split_t < alpha)
  Type_I_error_perm.test_split[j, 1] <- mean(pval_split_perm < alpha)
  Type_I_error_adaptive.test_split[j, 1] <- mean((pval_sw_split_x > alpha & pval_sw_split_y > alpha) * (pval_split_t < alpha)) + 
    mean((pval_sw_split_x <= alpha & pval_sw_split_y <= alpha) * (pval_split_perm < alpha))
  Type_I_error_conditional.test_split[j, 1] <- round(mean(pval_split_t[pval_sw_split_x > alpha & pval_sw_split_y > alpha] < alpha), 3)
}

close(pb)

# Save results
save(pvalues_list, pvalues_split_list, power_sw, Type_I_error_t.test, Type_I_error_perm.test, 
     Type_I_error_adaptive.test, Type_I_error_conditional.test, power_sw_split, Type_I_error_t.test_split, 
     Type_I_error_perm.test_split, Type_I_error_adaptive.test_split, Type_I_error_conditional.test_split, 
     file = "simulation_results.RData")

# Save p-values as CSV files
for (n in sample_size) {
  write.csv(pvalues_list[[paste0("n_", n)]], file = paste0("pvalues_n", n, ".csv"), row.names = FALSE)
  write.csv(pvalues_split_list[[paste0("n_", n)]], file = paste0("pvalues_split_n", n, ".csv"), row.names = FALSE)
}

print("Simulation completed. P-values and Type I errors have been stored and saved.")
