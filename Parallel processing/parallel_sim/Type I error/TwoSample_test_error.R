N <- 1e4; alpha <- 0.05
sample_size <- c(10, 50) #c(10,20,30,40)
set.seed(33)

# t test method
TypeI_error_norm <- TypeI_error_exp <- TypeI_error_logn <- numeric(length(sample_size))
for( k in 1 : length(sample_size)){
  n <- sample_size[k]
  print(n)
  
  pval_norm <- pval_exp <- pval_logn<- numeric(N)
  
  for(j in 1 : N){
    
    # ===============================================================================================
    # e ~exp(theta),  E(X) = theta and VAR(X) = theta^2
    # ===============================================================================================
    n1 <- rnorm(n, mean = 0, sd = 1)
    n2 <- rnorm(n + 5, mean = 0, sd = 1)
    
    pval_norm[j] <- t.test(n1, n2)$p.value
    
    # ===============================================================================================
    # e ~exp(theta),  E(X) = theta and VAR(X) = theta^2
    # ===============================================================================================
    e1 <- rexp(n, rate = 1) - 1 
    e2 <- rexp(n + 5, rate = 1) - 1 
    
    pval_exp[j] <- t.test(e1, e2)$p.value
    
    # ===================================================================================================
    # l1 ~lognormal(mu, sigma^2)  E(X) = exp(mu + sigma^2/2) VAR(X) = [exp(sigma^2)-1]exp(2mu + sigma^2)
    # ===================================================================================================
    l1 <- (rlnorm(n, meanlog = 0, sdlog = 1) - exp(0 + 1/2))/sqrt((exp(1)-1)*exp(2*0 + 1))
    l2 <- (rlnorm(n + 5, meanlog = 0, sdlog = 1) - exp(0 + 1/2))/sqrt((exp(1)-1)*exp(2*0 + 1))
    
    pval_logn[j] <- t.test(l1, l2)$p.value
    
  }
  # Power of t test
  TypeI_error_norm[k] <- round(mean(pval_norm < alpha), 3)
  TypeI_error_exp[k] <- round(mean(pval_exp < alpha), 3)
  TypeI_error_logn[k] <- round(mean(pval_logn < alpha), 3)
}

TypeI_error_norm
TypeI_error_exp
TypeI_error_logn
