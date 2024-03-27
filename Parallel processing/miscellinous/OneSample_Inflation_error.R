source("/home/kongyir/spring2024/User_defined_functions.R")
N <- 1e5; P <- 1e5; alpha <- 0.05
dist_sum <- c("Standard Normal",  "Uniform", "t", "Laplace", "Contaminated")
nvec <- c(8, 10, 15, 20, 30, 50)
sig.level <- c(0.01, 0.05, 0.07)

set.seed(33)

powervec <- numeric(length(nvec) * length(dist_sum) * length(sig.level))
error_t.test <- error_perm.test <- array(powervec, dim = c(6, 5, 3), dimnames = list(nvec, dist_sum, sig.level))

calculate_test_statistic <- function(x) {
  return((mean(x)*sqrt(length(x))) / sqrt(var(x)))
}

system.time({
  for (i in 1:length(nvec)) {
    n <- nvec[i]
    print(n)
    for (m in 1:length(sig.level)) {
      alpha <- sig.level[m]
      for (j in 1:length(dist_sum)) {
        dist <- dist_sum[j]
        pval <- pval_perm <- numeric(N)
        for (k in 1:N) {
          x <- generate_data(n, dist) 
          pval[k] <- t.test(x)$p.value
          observed_statistic <- calculate_test_statistic(x)
          permuted_statistics <- rep(0, P)
          for (l in 1:P) {
            myIndex <- sample(c(-1, 1), length(x), replace = TRUE)
            sample_data <- myIndex * abs(x)
            permuted_statistics[l] <- calculate_test_statistic(sample_data)
          }
          pval_perm[k] <- round(mean(abs(permuted_statistics) >= abs(observed_statistic)), 5)
        }
        error_t.test[i, j, m] <- mean(pval < alpha)
        error_perm.test[i, j, m] <- mean(pval_perm < alpha)
      }
    }
  }
  Inflation_error <- error_t.test - error_perm.test
})

error_t.test
error_perm.test
Inflation_error
save.image(paste0("OneSample_Inflation_error",".RData"))