rm(list = ls())
setwd("~/Desktop/OSU/Research/Pretest-Simulation/Power and Type I error Rate/Rmarkdown/test statistics plots/")
set.seed(33); N <- 1e5; alpha <- 0.05; d <- 0.5; n <- 10; df = n - 1

norm_statistic <- log_statistic <-t_statistic <- u_statistic<-e_statistic <- k_statistic <-g_statistic <- w_statistic <- NULL
for(j in 1 : N){
  # ===============================================================================================
  # n1 ~N(0,1),  E(X) = 0  and VAR(X) = 1
  # ===============================================================================================
  n1 <- rnorm(n, mean = 0, sd = 1)
  
  norm_statistic[j] <- t.test(n1)$statistic
  
  # ===============================================================================================
  # e ~exp(theta),  E(X) = theta and VAR(X) = theta^2
  # ===============================================================================================
  
  e1 <- rexp(n, rate = 1)
  e_statistic[j] <- t.test(e1 - 1)$statistic
  
  # ===============================================================================================
  # u1 ~U(a, b),  E(X) = (a + b)/2 and VAR(X) = [(b-a)^2]/12
  # ===============================================================================================
  u1 <- runif(n, min = 0, max = 1)
  
  u_statistic[j] <- t.test((u1-0.5)/sqrt(1/12))$statistic
  
  # ===============================================================================================
  # t1 ~t_v,  E(X) = 0 for v >1 and VAR(X) = v/(v-2)
  # ===============================================================================================
  t1 <-rt(n, df = 3) 
  
  t_statistic[j] <- t.test((t1)/sqrt(3))$statistic
  
  # ===============================================================================================
  # k1 ~X_k,  E(X) = k and VAR(X) = 2k
  # ===============================================================================================
  k1 <-rchisq(n, df = 7) 
  
  k_statistic[j] <- t.test((k1 - 7)/sqrt(14))$statistic
  
  # ===============================================================================================
  # g1 ~gamma(alpha, beta), where beta = 1/theta so E(X) = shape/rate and VAR(X) = shape/(rate)^2
  # ===============================================================================================
  g1 <- rgamma(n, shape = 3, rate = 0.1) 
  
  g_statistic[j] <- t.test((g1-30)/sqrt(300))$statistic
  
  # ===============================================================================================
  # w1 ~Weibull(shape, scale)  E(X) = scale Gamma(1 + 1/shape) and 
  # VAR(X) = scale^2[Gamma(1 + 2/shape) - (Gamma(1 + 1/shape))^2]
  # ===============================================================================================
  w1 <- rweibull(n, shape = 1, scale = 2) 
  u <- (w1-2*gamma(2))/sqrt(4*(gamma(3) - gamma(2)))
  
  
  w_statistic[j] <- t.test(u)$statistic
  
  # ===================================================================================================
  # l1 ~lognormal(mu, sigma^2)  E(X) = exp(mu + sigma^2/2) VAR(X) = [exp(sigma^2)-1]exp(2mu + sigma^2)
  # ===================================================================================================
  l1 <- rlnorm(n, meanlog = 0, sdlog = 1)
  z1 <- (l1- exp(0+1/2))/sqrt(sqrt((exp(1)-1)*exp(2*0 + 1)))
  
  log_statistic[j] <- t.test(z1)$statistic
  
}

# Calculation of Type 1 error

#normal Distn
left_norm <- mean(norm_statistic < qt(0.025, df = df))
right_norm <- mean(norm_statistic > qt(0.975, df = df))
Type1_error_norm <- sum(left_norm + right_norm)

#exponential Distn
left_e <- mean(e_statistic < qt(0.025, df = df))
right_e <- mean(e_statistic > qt(0.975, df = df))
Type1_error_e <- sum(left_e + right_e)
#t Distn
left_t <- mean(t_statistic < qt(0.025, df = df))
right_t <- mean(t_statistic > qt(0.975, df = df))
Type1_error_t <- sum(left_t + right_t)
#uuniform Distn
left_u <- mean(u_statistic < qt(0.025, df = df))
right_u <- mean(u_statistic > qt(0.975, df = df))
Type1_error_u <- sum(left_u + right_u)
#Chi-Squared Distn
left_k <- mean(k_statistic < qt(0.025, df = df))
right_k <- mean(k_statistic > qt(0.975, df = df))
Type1_error_k <- sum(left_k + right_k)
#Weilbul Distn
left_w <- mean(w_statistic < qt(0.025, df = df))
right_w <- mean(w_statistic > qt(0.975, df = df))
Type1_error_w <- sum(left_w + right_w)
#Gamma Distn
left_g <- mean(g_statistic < qt(0.025, df = df))
right_g <- mean(g_statistic > qt(0.975, df = df))
Type1_error_g <- sum(left_g + right_g)
#Lognormal Distn
left_l <- mean(log_statistic < qt(0.025, df = df))
right_l <- mean(log_statistic > qt(0.975, df = df))
Type1_error_l <- sum(left_l + right_l)

save.image(paste0("one_samples_plots",".RData"))