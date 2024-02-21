N <- 1e4; alpha <- 0.05; d <- 0.5; P= 1e4
sample_size <- c(10,20,30,40,50)
set.seed(33)

# t test method
power_t_test_normal <- power_t_test_exp <- power_t_test_t <- power_t_test_chisq  <- power_t_test_gamma <- 
  power_t_test_weibul <- power_t_test_logn<- numeric(length(sample_size))
for( k in 1 : length(sample_size)){
  n <- sample_size[k]
  print(n)
  
  pval_norm <- pval_exp <- pval_t <- pval_chisq <- pval_gamma <- pval_weibul <- pval_logn<- numeric(N)
  
  for(j in 1 : N){
    
    # ===============================================================================================
    # e ~exp(theta),  E(X) = theta and VAR(X) = theta^2
    # ===============================================================================================
    n1 <- rnorm(n, mean = 0, sd = 1)
    n2 <- rnorm(n, mean = 0, sd = 1) + d
    
    pval_norm[j] <- t.test(n1, n2, var.equal = F)$p.value
    
    # ===============================================================================================
    # e ~exp(theta),  E(X) = theta and VAR(X) = theta^2
    # ===============================================================================================
    e1 <- rexp(n, rate = 1) - 1 
    e2 <- rexp(n, rate = 1) -1 + d
    
    pval_exp[j] <- t.test(e1, e2, var.equal = F)$p.value
    
    # ===============================================================================================
    # t1 ~t_v,  E(X) = 0 for v >1 and VAR(X) = v/(v-2)
    # ===============================================================================================
    t1 <-rt(n, df = 3)/sqrt(3)
    t2 <-rt(n, df = 3)/sqrt(3) + d
    
    pval_t[j] <- t.test(t1, t2, var.equal = F)$p.value
    
    # ===============================================================================================
    # k1 ~X_k,  E(X) = k and VAR(X) = 2k
    # ===============================================================================================
    k1 <-(rchisq(n, df = 3) - 3)/sqrt(6)
    k2 <-(rchisq(n, df = 3) - 3)/sqrt(6) + d
    
    pval_chisq[j] <- t.test(k1, k2, var.equal = F)$p.value
    
    # ===============================================================================================
    # g1 ~gamma(alpha, beta), where beta = 1/theta so E(X) = alpha/theta and VAR(X) = alpha/(theta)^2
    # ===============================================================================================
    g1 <- (rgamma(n, shape = 3, rate = 0.1) -300)/sqrt(300)
    g2 <- (rgamma(n, shape = 3, rate = 0.1) -300)/sqrt(300) +d
    
    pval_gamma[j] <- t.test(g1, g2, var.equal = F)$p.value
    
    # ===============================================================================================
    # w1 ~Weibull(shape, scale)  E(X) = scale * Gamma(1 + 1/shape) and 
    # VAR(X) = scale^2[Gamma(1 + 2/shape) - (Gamma(1 + 1/shape))^2]
    # ===============================================================================================
    w1 <- (rweibull(n, shape = 1, scale = 2) -2*gamma(51/50))/sqrt(4*(gamma(3) - gamma(2)))
    w2 <- (rweibull(n, shape = 1, scale = 2) -2*gamma(51/50))/sqrt(4*(gamma(3) - gamma(2))) + d
    
    pval_weibul[j] <- t.test(w1, w2, var.equal = F)$p.value  
    
    # ===================================================================================================
    # l1 ~lognormal(mu, sigma^2)  E(X) = exp(mu + sigma^2/2) VAR(X) = [exp(sigma^2)-1]exp(2mu + sigma^2)
    # ===================================================================================================
    l1 <- (rlnorm(n, meanlog = 0, sdlog = 1) - exp(0 + 1/2))/sqrt((exp(1)-1)*exp(2*0 + 1))
    l2 <- (rlnorm(n, meanlog = 0, sdlog = 1) - exp(0 + 1/2))/sqrt((exp(1)-1)*exp(2*0 + 1)) + d
    
    pval_logn[j] <- t.test(l1, l2, var.equal = F)$p.value
    
  }
  # Power of t test
  power_t_test_normal[k] <- round(mean(pval_norm < alpha), 3)
  power_t_test_exp[k] <- round(mean(pval_exp < alpha), 3)
  power_t_test_t[k] <- round(mean(pval_t < alpha), 3)
  power_t_test_chisq[k] <- round(mean(pval_chisq < alpha), 3)
  power_t_test_gamma[k] <- round(mean(pval_gamma < alpha), 3)
  power_t_test_weibul[k] <- round(mean(pval_weibul < alpha), 3)
  power_t_test_logn[k] <- round(mean(pval_logn < alpha), 3)
}

# Function to calculate the test statistic (difference of means)
calculate_test_statistic <- function(x, y) {
  return((mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))
}

# Permutation test approach
power_perm_test_norm <- power_perm_test_exp <- power_perm_test_t <- power_perm_test_chisq <- 
power_perm_test_gamma <- power_perm_test_weibul <- power_perm_test_logn<- numeric(length(sample_size))

for(k in 1 : length(sample_size)){
  n <- sample_size[k]
  print(n)
  pval_perm_norm <- pval_perm_exp <- pval_perm_t <- pval_perm_logn <- pval_perm_chisq <- 
    pval_perm_gamma <- pval_perm_weibul<- numeric(N)
  for(j in 1 : N){
    # ===============================================================================================
    # e ~exp(theta),  E(X) = theta and VAR(X) = theta^2
    # ===============================================================================================
    n1 <- rnorm(n, mean = 0, sd = 1)
    n2 <- rnorm(n, mean = 0, sd = 1) + d
    permuted_data_norm <- c(n1, n2)
    # ===============================================================================================
    # e ~exp(theta),  E(X) = theta and VAR(X) = theta^2
    # ===============================================================================================
    e1 <- rexp(n, rate = 1) - 1 
    e2 <- rexp(n, rate = 1) - 1 + d
    permuted_data_exp <- c(e1, e2)
    
    # ===============================================================================================
    # t1 ~t_v,  E(X) = 0 for v >1 and VAR(X) = v/(v-2)
    # ===============================================================================================
    t1 <-rt(n, df = 3)/sqrt(3)
    t2 <-rt(n, df = 3)/sqrt(3) + d
    permuted_data_t <- c(t1, t2)
    
    # ===============================================================================================
    # k1 ~X_k,  E(X) = k and VAR(X) = 2k
    # ===============================================================================================
    k1 <-(rchisq(n, df = 3) - 3)/sqrt(6)
    k2 <-(rchisq(n, df = 3) - 3)/sqrt(6) + d
    permuted_data_chisq <- c(k1, k2)
    
    # ===============================================================================================
    # g1 ~gamma(alpha, beta), where beta = 1/theta so E(X) = alpha/theta and VAR(X) = alpha/(theta)^2
    # ===============================================================================================
    g1 <- (rgamma(n, shape = 3, rate = 0.1) -300)/sqrt(300)
    g2 <- (rgamma(n, shape = 3, rate = 0.1) -300)/sqrt(300) +d
    permuted_data_gamma <- c(g1, g2)
    # ===============================================================================================
    # w1 ~Weibull(shape, scale)  E(X) = scale * Gamma(1 + 1/shape) and 
    # VAR(X) = scale^2[Gamma(1 + 2/shape) - (Gamma(1 + 1/shape))^2]
    # ===============================================================================================
    w1 <- (rweibull(n, shape = 1, scale = 2) -2*gamma(51/50))/sqrt(4*(gamma(3) - gamma(2)))
    w2 <- (rweibull(n, shape = 1, scale = 2) -2*gamma(51/50))/sqrt(4*(gamma(3) - gamma(2))) + d
    permuted_data_weibull <- c(w1, w2)
    
    # ===================================================================================================
    # l1 ~lognormal(mu, sigma^2)  E(X) = exp(mu + sigma^2/2) VAR(X) = [exp(sigma^2)-1]exp(2mu + sigma^2)
    # ===================================================================================================
    l1 <- (rlnorm(n, meanlog = 0, sdlog = 1) - exp(0 + 1/2))/sqrt((exp(1)-1)*exp(2*0 + 1))
    l2 <- (rlnorm(n, meanlog = 0, sdlog = 1) - exp(0 + 1/2))/sqrt((exp(1)-1)*exp(2*0 + 1)) + d
    permuted_data_logn <- c(l1, l2)
    
    # Observed test statistic
    observed_statistic_norm <- calculate_test_statistic(n1, n2)
    observed_statistic_exp <- calculate_test_statistic(e1, e2)
    observed_statistic_t <- calculate_test_statistic(t1, t2)
    observed_statistic_chisq <- calculate_test_statistic(k1, k2)
    observed_statistic_gamma <- calculate_test_statistic(g1, g2)
    observed_statistic_weibul <- calculate_test_statistic(w1, w2)
    observed_statistic_logn <- calculate_test_statistic(l1, l2)
    
    permuted_statistics_norm <- permuted_statistics_exp <- permuted_statistics_t <-permuted_statistics_logn <- 
      permuted_statistics_chisq <- permuted_statistics_gamma <- permuted_statistics_weibul<- rep(0, P)
    
    for(i in 1 : P){
      permuted_data_norm <- sample(permuted_data_norm)
      permuted_data_exp <- sample(permuted_data_exp)
      permuted_data_t <- sample(permuted_data_t)
      permuted_data_chisq <- sample(permuted_data_chisq)
      permuted_data_gamma <- sample(permuted_data_gamma)
      permuted_data_weibull <- sample(permuted_data_weibull)
      permuted_data_logn <- sample(permuted_data_logn)
      
      # Calculate the test statistic for the permuted data
      # normal
      permuted_norm1 <- permuted_data_norm[1:length(n1)]
      permuted_norm2 <- permuted_data_norm[(length(n1) + 1):(length(n1) + length(n2))]
      permuted_statistics_norm[i] <- calculate_test_statistic(permuted_norm1, permuted_norm2)
      
      # exp
      permuted_exp1 <- permuted_data_exp[1:length(e1)]
      permuted_exp2 <- permuted_data_exp[(length(e1) + 1):(length(e1) + length(e2))]
      permuted_statistics_exp[i] <- calculate_test_statistic(permuted_exp1, permuted_exp2)
      # t distn
      permuted_t1 <- permuted_data_t[1:length(t1)]
      permuted_t2 <- permuted_data_t[(length(t1) + 1):(length(t1) + length(t2))]
      permuted_statistics_t[i] <- calculate_test_statistic(permuted_t1, permuted_t2)
      # chisq distn
      permuted_chisq1 <- permuted_data_chisq[1:length(k1)]
      permuted_chisq2 <- permuted_data_chisq[(length(k1) + 1):(length(k1) + length(k2))]
      permuted_statistics_chisq[i] <- calculate_test_statistic(permuted_chisq1, permuted_chisq2)
      
      # gamma
      permuted_gamma1 <- permuted_data_gamma[1:length(g1)]
      permuted_gamma2 <- permuted_data_gamma[(length(g1) + 1):(length(g1) + length(g2))]
      permuted_statistics_gamma[i] <- calculate_test_statistic(permuted_gamma1, permuted_gamma2)
      # weibul distn
      permuted_weibul1 <- permuted_data_weibull[1:length(w1)]
      permuted_weibul2 <- permuted_data_weibull[(length(w1) + 1):(length(w1) + length(w2))]
      permuted_statistics_weibul[i] <- calculate_test_statistic(permuted_weibul1, permuted_weibul2)
      # logn distn
      permuted_logn1 <- permuted_data_logn[1:length(l1)]
      permuted_logn2 <- permuted_data_chisq[(length(l1) + 1):(length(l1) + length(l2))]
      permuted_statistics_logn[i] <- calculate_test_statistic(permuted_logn1, permuted_logn2)
      
    }
    # Calculate the pvalue of test
    pval_perm_norm[j] <- round(mean(abs(permuted_statistics_norm) >= abs(observed_statistic_norm)), 5)
    pval_perm_exp[j] <- round(mean(abs(permuted_statistics_exp) >= abs(observed_statistic_exp)), 5)
    pval_perm_t[j] <- round(mean(abs(permuted_statistics_t) >= abs(observed_statistic_t)), 5) 
    pval_perm_chisq[j] <- round(mean(abs(permuted_statistics_chisq) >= abs(observed_statistic_chisq)), 5) 
    pval_perm_gamma[j] <- round(mean(abs(permuted_statistics_gamma) >= abs(observed_statistic_gamma)), 5)
    pval_perm_weibul[j] <- round(mean(abs(permuted_statistics_weibul) >= abs(observed_statistic_weibul)), 5) 
    pval_perm_logn[j] <- round(mean(abs(permuted_statistics_logn) >= abs(observed_statistic_logn)), 5) 
    
    
  }
  
  power_perm_test_norm[k] <- round(mean(pval_perm_norm < alpha), 5)
  power_perm_test_exp[k] <- round(mean(pval_perm_exp < alpha), 5)
  power_perm_test_t[k] <- round(mean(pval_perm_t < alpha), 5)
  power_perm_test_chisq[k] <- round(mean(pval_perm_chisq < alpha), 5)
  power_perm_test_gamma[k] <- round(mean(pval_perm_gamma < alpha), 5)
  power_perm_test_weibul[k] <- round(mean(pval_perm_weibul < alpha), 5)
  power_perm_test_logn[k] <- round(mean(pval_perm_logn < alpha), 5)
  
}
# power loss
power_loss_norm <- power_perm_test_norm - power_t_test_normal
power_loss_exp <- power_perm_test_exp - power_t_test_exp
power_loss_t <- power_perm_test_t - power_t_test_t
power_loss_chisq <- power_perm_test_chisq - power_t_test_chisq
power_loss_gamma <- power_perm_test_gamma - power_t_test_gamma
power_loss_weibul <- power_perm_test_weibul - power_t_test_weibul
power_loss_logn<- power_perm_test_logn - power_t_test_logn


# power under t test
power_t<-c(power_t_test_normal, power_t_test_exp, power_t_test_t, power_t_test_chisq, 
           power_t_test_gamma, power_t_test_weibul, power_t_test_logn)
power_table_t<-array(power_t, dim = c(5, 7), dimnames = list(sample_size,c("norm_t", 
                    "exp_t", "t distn_t","chisq", "Gamma", "Weibull", "Lognormal")))
print(power_table_t)

# power under permutation test
power_perm<-c(power_perm_test_norm, power_perm_test_exp, power_perm_test_t,  power_perm_test_chisq, 
              power_perm_test_gamma, power_perm_test_weibul, power_perm_test_logn)
power_table_perm<-array(power_perm, dim = c(5, 7), dimnames = list(sample_size, c("Normal", 
                          "exp_perm", "t_perm", "chisq", "Gamma", "Weibull", "Lognormal")))
print(power_table_perm)

# power loss
powerloss_data<-c(power_loss_norm, power_loss_exp,power_loss_t, power_loss_chisq, 
                  power_loss_gamma, power_loss_weibul, power_loss_logn)
powerloss_table<-array(powerloss_data, dim = c(5, 7), dimnames = list(sample_size, c("Normal", "exp_perm", 
                                                "t_perm", "chisq", "Gamma", "Weibull", "Lognormal")))
print(powerloss_table)


save.image(paste0("result2.0",".RData"))
