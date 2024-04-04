
set.seed(33)
nsim <- 1e4; alpha <- 0.05; sam_size <- c(5, 10, 20) ; msam = c(10, 25, 50)

error_norm <- error_exp <- error_logn <- numeric(length(sam_size))

for (i in seq_along(sam_size)) {
  n <- sam_size[i]
  m = msam[i]
  pval_norm <- pval_exp <- pval_logn <- numeric(nsim)
  for (j in 1:nsim) {
    # standard normal
    n1 <- rnorm(n, mean = 0, sd = 1)
    n2 <- rnorm(m, mean = 0, sd = 1)
    pval_norm[j] <- t.test(n1, n2)$p.value
    test_stat = t.test(n1, n2)$statistic
    # exponential
    e1 = rexp(n, rate = 0.5) - 2
    e2 = rexp(m, rate = 0.5) - 2
    pval_exp[j] = t.test(e1, e2)$p.value
  }
  
  # Calculate errors and store them
  error_norm[i] <- mean(pval_norm < alpha)
  error_exp[i] <- mean(pval_exp < alpha)
  
}
#print errors
error_norm
error_exp
