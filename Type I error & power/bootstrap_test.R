source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")

one_sample_bootstrap_test <- function(n, dist, test_stat, n_bootstrap, effect_size = 0){
  x <- generate_data(n, dist) + effect_size
  # observe test stats
  observe_test_stat <- test_stat(x)
  # bootstrap test stats
  x <- x - mean(x)
  bootstrap_test_stat <- numeric(n_bootstrap)
  for(b in 1 : n_bootstrap){
    sample1 <- sample(x, n, replace = TRUE)
    bootstrap_test_stat[b] <- test_stat(sample1)
    bootstrap_p_val <- mean(abs(bootstrap_test_stat) >= abs(observe_test_stat))
  }
  return(bootstrap_p_val)
}

set.seed(12345)
N <- 1000         
alpha <- 0.05
sample_sizes <- c(10, 20, 30, 40, 50)
distributions <- c("normal", "exponential")
n_bootstrap <- 1000
effect_size <- 0.5

# test_stat <- function(x){
#   sqrt(n) * (mean(x) - 0)/sd(x)
# }

test_stat <- function(x){
  mean(x)
}

run_simulation <- function(test_stat, effect_size = 0.0){
  results <- array(data = NA, dim = c(length(sample_sizes), length(distributions)), dimnames = list(sample_sizes = as.character(sample_sizes), distributions = distributions))
  for (i in seq_along(sample_sizes)) {
    n <- sample_sizes[i]
    for (j in seq_along(distributions)) {
      dist <- distributions[j]
      pval <- replicate(N, one_sample_bootstrap_test(n = n, dist = dist, test_stat = test_stat,  n_bootstrap = n_bootstrap, effect_size = effect_size))
      results[i, j] <- round(mean(pval < alpha, na.rm = TRUE), 3)
    }
  }
return(results)
}

# run simulation
run_simulation(test_stat = test_stat, effect_size = effect_size)


