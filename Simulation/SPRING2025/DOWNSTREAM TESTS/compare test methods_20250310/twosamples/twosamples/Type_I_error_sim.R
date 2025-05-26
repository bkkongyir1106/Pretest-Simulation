setwd("~/Library/Mobile Documents/com~apple~CloudDocs/Research/Type I error")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

## Set up the simulation parameters
{
  Nsim <- 1e4
  sample_size <- c(8, 10, 15, 20, 25, 30)
  distribution <- c("LogNormal", "Exponential")
  alpha <- 0.05
  B <- 1e3
  effect_size <- 0.0
}

# Pre-allocate matrices for Type I error rates
power_sw <- Type_I_error_t.test  <- Type_I_error_perm.test <- Type_I_error_adaptive.test <- 
  Type_I_error_conditional.test <- array(data = NA, dim = c(length(sample_size), length(distribution)),
                                         dimnames = list(sample_size, distribution))
# split method
power_sw_split <- Type_I_error_t.test_split  <- Type_I_error_perm.test_split <- Type_I_error_adaptive.test_split <- 
  Type_I_error_conditional.test_split <- array(data = NA, dim = c(length(sample_size), length(distribution)),
                                         dimnames = list(sample_size, distribution))

# Create an empty list to store p-values 
pvalues_list <- list()

pb <- txtProgressBar(min = 0, max = length(distribution) * length(sample_size) * Nsim, style = 3)
progress_counter <- 0

for(i in seq_along(distribution)){
  dist <- distribution[i]
  for(j in seq_along(sample_size)){
    n <- sample_size[j]
    
    # create vectors to store p-values 
    pval_sw <- numeric(Nsim)
    pval_t <- numeric(Nsim)
    pval_perm <- numeric(Nsim)
    pval_adaptive <- numeric(Nsim)
    pval_conditional <- numeric(Nsim)
    
    # data splitting method
    pval_sw_split <- numeric(Nsim)
    pval_t_split <- numeric(Nsim)
    pval_perm_split <- numeric(Nsim)
    pval_adaptive_split <- numeric(Nsim)
    pval_conditional_split <- numeric(Nsim)
    
    for (k in 1:Nsim) {
      x <- generate_data(n, dist)
      y <- generate_data(n, dist)
      
      # Perform Shapiro-Wilk test on x
      pval_sw[k] <- shapiro.test(x)$p.value
      
      # Perform t-test
      pval_t[k] <- t.test(x, y, mu = 0)$p.value
      
      # Perform permutation test
      observe_stat <- median(x) - median(y)
      permuted_stat <- numeric(B)
      combined_data <- c(x, y)  
      
      for (p in 1:B) {
        permuted_data <- sample(combined_data)
        sample_x <- permuted_data[1:length(x)]
        sample_y <- permuted_data[(length(x) + 1):(length(x) + length(y))]
        permuted_stat[p] <- median(sample_x) - median(sample_y)
      }
      pval_perm[k] <- mean(abs(permuted_stat) >= abs(observe_stat))
      
      # Adaptive t-perm test
      if (shapiro.test(x)$p.value > alpha & shapiro.test(y)$p.value > alpha) {
        pval_adaptive[k] <- pval_t[k]
      } else {
        pval_adaptive[k] <- pval_perm[k]
      }
      
      # Conditional t-test
      if (pval_sw[k] > alpha & shapiro.test(y)$p.value > alpha) {
        pval_conditional[k] <- pval_t[k]
      } else {
        pval_conditional[k] <- NA
      }
    
      # Split method
      split_point <- floor(length(x)/2)
      first_half.x <- x[1:split_point]
      second_half.x <- x[(split_point + 1):length(x)]
      first_half.y <- y[1:split_point]
      second_half.y <- y[(split_point + 1):length(y)]
      
      pval_sw_split[k] <- shapiro.test(first_half.x)$p.value
      
      # Perform t-test (split)
      pval_t_split[k] <- t.test(second_half.x, second_half.y, mu = 0)$p.value
      
      # Perform permutation test
      observe_stat_split <- median(second_half.x) - median(second_half.y)
      permuted_stat_split <- numeric(B)
      combined_data_split <- c(second_half.x, second_half.y)  
      for (p in 1:B) {
        permuted_data_split <- sample(combined_data_split)
        sample_x_split <- permuted_data_split[1:length(second_half.x)]
        sample_y_split <- permuted_data_split[(length(second_half.x) + 1):(length(second_half.x) + length(second_half.y))]
        permuted_stat_split[p] <- median(sample_x_split) - median(sample_y_split)
      }
      pval_perm_split[k] <- mean(abs(permuted_stat_split) >= abs(observe_stat_split))
      
      # Adaptive t-perm test
      if (shapiro.test(first_half.x)$p.value > alpha & shapiro.test(first_half.y)$p.value > alpha) {
        pval_adaptive_split[k] <- pval_t_split[k]
      } else {
        pval_adaptive_split[k] <- pval_perm_split[k]
      }
      
      # Conditional t-test
      if (pval_sw_split[k] > alpha & shapiro.test(first_half.y)$p.value > alpha) {
        pval_conditional_split[k] <- pval_t_split[k]
      } else {
        pval_conditional_split[k] <- NA
      }
      
      progress_counter <- progress_counter + 1
      setTxtProgressBar(pb, progress_counter)
    }
    
    # Save the p-values 
    setting_name <- paste0("n", n, "_", dist)
    pvalues_list[[setting_name]] <- list(
      pval_sw = pval_sw,
      pval_t = pval_t,
      pval_perm = pval_perm,
      pval_adaptive = pval_adaptive,
      pval_conditional = pval_conditional,
      pval_sw_split = pval_sw_split,
      pval_t_split = pval_t_split,
      pval_perm_split = pval_perm_split,
      pval_adaptive_split = pval_adaptive_split,
      pval_conditional_split = pval_conditional_split
    )
    
    # Calculate power and Type I error rates
    power_sw[j, i] <- mean(pval_sw < alpha)
    Type_I_error_t.test[j, i] <- mean(pval_t < alpha)
    Type_I_error_perm.test[j, i] <- mean(pval_perm < alpha)
    Type_I_error_adaptive.test[j, i] <- mean(pval_adaptive < alpha)
    Type_I_error_conditional.test[j, i] <- round(mean(pval_conditional < alpha, na.rm = TRUE), 3)
    
    # split method
    power_sw_split[j, i] <- mean(pval_sw_split < alpha)
    Type_I_error_t.test_split[j, i] <- mean(pval_t_split < alpha)
    Type_I_error_perm.test_split[j, i] <- mean(pval_perm_split < alpha)
    Type_I_error_adaptive.test_split[j, i] <- mean(pval_adaptive_split < alpha)
    Type_I_error_conditional.test_split[j, i] <- round(mean(pval_conditional_split < alpha, na.rm = TRUE), 3)
  }
}

close(pb)

# Print results
print(power_sw)
print(Type_I_error_t.test)
print(Type_I_error_perm.test)
print(Type_I_error_adaptive.test)
print(Type_I_error_conditional.test)

#split method
print(power_sw_split)
print(Type_I_error_t.test_split)
print(Type_I_error_perm.test_split)
print(Type_I_error_adaptive.test_split)
print(Type_I_error_conditional.test_split)

# save RData
save.image(file = "TypeI_error_two_sample_test_RData")

save(sample_size, distribution,  Nsim, power_sw, Type_I_error_t.test, Type_I_error_adaptive.test, Type_I_error_conditional.test, Type_I_error_t.test_split,  Type_I_error_adaptive.test_split, file = "twosample_Type_I_error.RData")









