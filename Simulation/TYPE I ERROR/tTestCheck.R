N = 1e3
samples = c(10, 20, 30)
alpha = 0.05
# t test
error = NULL
for(j in seq_along(samples)){
  n = samples[j]
  pval = numeric(N)
  for(i in 1 : N){
    x = rexp(n, rate = 3) - 1/3
    #x = rnorm(n, 0, 1)
    pval[i] = t.test(x)$p.value
  }
  error[j] = mean(pval < alpha)
}
error

# Two-stage procedure
error1 = NULL
for(j in seq_along(samples)){
  n = samples[j]
  pval1 = numeric(N)
  for(i in 1 : N){
    x = rexp(n, rate = 3) - 1/3
    #x = rnorm(n, 0, 1)
    if(shapiro.test(x)$p.value < alpha){
      pval1[i] = t.test(x)$p.value
    }else{
      pval1[i] = wilcox.test(x)$p.value
    }
  }
  error1[j] = mean(pval1 < alpha)
}
error1

# Conditional Type I error rate
error2 = NULL
for(j in seq_along(samples)){
  n = samples[j]
  nsim = 0
  samplepast = 0
  pval2 = NULL
  while(samplepast < N){
    x = rexp(n, rate = 3) - 1/3
    #x = rnorm(n, 0, 1)
    nsim = nsim + 1
    if(shapiro.test(x)$p.value > alpha){
      pval2[samplepast] = t.test(x)$p.value
      samplepast = samplepast + 1
    }
  }
  error2[j] = mean(pval2 < alpha)
}
error2

