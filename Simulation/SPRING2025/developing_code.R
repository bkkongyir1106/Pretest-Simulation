
N = 100
alpha = 0.05
B = 100
dist = "Chi-Square"
effect_size = 0.0
sample_size <- c(10, 20, 30)
traditionl.Type1error = numeric(length(sample_size))
adaptive_t_perm.Type1error = numeric(length(sample_size))
adaptive_t_boot.Type1error = numeric(length(sample_size))
for (j in seq_along(sample_size)) {
  n <- sample_size[j]
  number_t_perm.test<- 1
  number.perm.test <- 1
  
  number_t_boot.test<- 1
  number.boot.test <- 1
  
  pval_t.test = c()
  pval.perm = c()
  pval_t_perm.test = c()
  pval_t_boot.test = c()
  pval.boot = c()
  for(i in 1: N){
    x = generate_data(n, dist)
    y = generate_data(n, dist)
    pval_t.test[i] <- TwoSample.test(x, y, test = "t", alpha = alpha, effect_size = effect_size, B = NULL)
    # adaptive permutation
    if(shapiro.test(x)$p.value > alpha){
      pval_t_perm.test[number_t_perm.test] <- TwoSample.test(x, y, test = "t", alpha = alpha, effect_size = effect_size, B = NULL)
      number_t_perm.test <- number_t_perm.test + 1
    }else{
      pval.perm[number.perm.test] <- TwoSample.test(x, y, test = "perm", alpha = alpha, effect_size = effect_size, B = B)
      number.perm.test = number.perm.test + 1
    }
    # adaptive bootstrap
    if(shapiro.test(x)$p.value > alpha){
      pval_t_boot.test[number_t_boot.test] <- TwoSample.test(x, y, test = "t", alpha = alpha, effect_size = effect_size, B = NULL)
      number_t_boot.test <- number_t_boot.test + 1
    }else{
      pval.boot[number.boot.test] <- bootstrap_two_sample_test(x, y, effect_size = effect_size, alpha = alpha, n_bootstrap = B, sample_size = n)
      number.boot.test = number.boot.test + 1
    }
  }
  traditionl.Type1error[j] <- mean(pval_t.test < alpha)
  adaptive_t_perm.Type1error[j] <- ((number_t_perm.test - 1)/N) * (if (length(pval_t_perm.test) == 0) 0 else mean(pval_t_perm.test < alpha)) +
    ((number.perm.test-1)/N) * (if (length(pval.perm) == 0) 0 else mean(pval.perm < alpha))
  adaptive_t_boot.Type1error[j] <- ((number_t_boot.test - 1)/N) * (if (length(pval_t_boot.test) == 0) 0 else mean(pval_t_boot.test < alpha)) +
    ((number.boot.test-1)/N) * (if (length(pval.boot) == 0) 0 else mean(pval.boot < alpha))
}
traditionl.Type1error
adaptive_t_perm.Type1error
adaptive_t_boot.Type1error














dist = "Normal"
alpha = 0.05
n = 10
N = 1000
pvalw = numeric(N)
for( i in 1 : N){
  y = generate_data(n, dist)
  pvalw[i] = t.test(y)$p.value 
}
error = mean(pvalw < alpha)
error
