N = 1000; alpha = 0.05
samp = c(10, 20, 50)
errort = errorw = NULL
for(j in 1 : 3){
  n = samp[j]
  pvalt = pvalw = numeric(N)
  for( i in 1: N){
    x = rexp(n, rate = 1) - 1
    #x = rnorm(10, mean = 0, sd = 1)
    pvalw[i] = wilcox.test(x)$p.value
    pvalt[i] = t.test(x)$p.value
  }
  errort[j] = mean(pvalt < alpha)
  errorw[j] = mean(pvalw < alpha)
}
errort
errorw