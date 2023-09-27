rm(list = ls())
setwd("C:/Users/kongy/OneDrive/Desktop/Research/Pretest-Simulation/power loss")
set.seed(1)
N <- 1000000
alpha <- 0.05
df <- 0.5
sample_size <- c(10,20,30,40,50)
Type1_error_Norm <- Type1_error_exp <- Type1_error_unif <- Type1_error_beta <- Type1_error_gamma <- Type1_error_chisq <- Type1_error_weibull <- Type1_error_logn<- numeric(length(sample_size))
inflation_Type1_error_exp <- inflation_Type1_error_unif <- inflation_Type1_error_beta <- inflation_Type1_error_gamma <- inflation_Type1_error_chisq <- inflation_Type1_error_weibull <- inflation_Type1_error_logn <- numeric(length(sample_size))
Power_SW_test_exp <- Power_SW_test_unif <- Power_SW_test_beta <- Power_SW_test_gamma <- Power_SW_test_chisq <- Power_SW_test_weibul<- Power_SW_test_logn <- numeric(length(sample_size))
power_t_test_Norm <- power_t_test_exp <- power_t_test_unif <- power_t_test_gamma <- power_t_test_beta <- power_t_test_chisq <- power_t_test_weibull <- power_t_test_logn <- numeric(length(sample_size))
powrloss_exp <-powrloss_unif  <-powrloss_gamma <-powrloss_beta <-powrloss_chisq <- powrloss_weibull <- powrloss_logn <- numeric(length(sample_size))
for (i in 1 : length(sample_size)) {
  n <- sample_size[i]
  print(n)
  
  rejectH0_SW_test_exp <- rejectH0_SW_test_unif <- rejectH0_SW_test_beta <- rejectH0_SW_test_gamma <- rejectH0_SW_test_chisq <- rejectH0_SW_test_weibull<- rejectH0_SW_test_logn<- numeric(N)
  rejectH0_t_test_exp <- rejectH0_t_test_unif <- rejectH0_t_test_beta <- rejectH0_t_test_norm <- rejectH0_t_test_gamma <- rejectH0_t_test_chisq <- rejectH0_t_test_weibull <- rejectH0_t_test_logn<- numeric(N)
  powr_t_test_norm <- powr_t_test_exp <- powr_t_test_unif <- powr_t_test_beta <- powr_t_test_gamma <- powr_t_test_chisq <- powr_t_test_weibull<- powr_t_test_logn <-numeric(N)
  for(j in 1 : N){
    
    #========================================
    #       Normal distn                   #
    # =======================================
    x <- rnorm(n, mean =0, sd = 1)
    if(t.test(x)$p.value <= alpha){
      rejectH0_t_test_norm[j] <- 1
    }#Type I error
    if(t.test(x + df)$p.value <= alpha){
      powr_t_test_norm[j] <- 1
    }#power of t test
    #========================================
    #       exponential distn               #
    # =======================================
    e <- rexp(n, rate = 1) - 1
    if(t.test(e)$p.value <= alpha){
      rejectH0_t_test_exp[j] <- 1
    }#Type I error
    if(shapiro.test(e)$p.value <= alpha){
      rejectH0_SW_test_exp[j] <- 1
    }#power of SW test
    if(t.test(e + df)$p.value <= alpha){
      powr_t_test_exp[j] <- 1
    }#power of t test
    #========================================
    #       uniform distn                   #
    # =======================================
    u <- runif(n, min = 0, max = 1) - 0.5
    y <- u/sqrt(1/12)
    if(t.test(y)$p.value <= alpha){
      rejectH0_t_test_unif[j] <- 1
    }# Type I error 
    if(shapiro.test(y)$p.value <= alpha){
      rejectH0_SW_test_unif[j] <- 1
    }# power of SW test
    if(t.test(y + df)$p.value <= alpha){
      powr_t_test_unif[j] <- 1
    }#power of t test
    #========================================
    #       beta distn                      #
    # =======================================
    b0 <- rbeta(n, shape1 = 2, shape2 = 9) - 2/11
    b <- b0/(sqrt(3/242))
    if(t.test(b)$p.value <= alpha){
      rejectH0_t_test_beta[j] <- 1
    }# Type I error 
    if(shapiro.test(b)$p.value <= alpha){
      rejectH0_SW_test_beta[j] <- 1
    }#power of SW test
    if(t.test(b + df)$p.value <= alpha){
      powr_t_test_beta[j] <- 1
    }#power of t test
    
    #========================================
    #       gamma distn                     #
    # =======================================
    g <- rgamma(n, shape = 3, rate = 1) - 3
    z <- g/sqrt(3)
    if(t.test(z)$p.value <= alpha){
      rejectH0_t_test_gamma[j] <- 1
    }# Type I error 
    if(shapiro.test(z)$p.value <= alpha){
      rejectH0_SW_test_gamma[j] <- 1
    }#power of SW test
    if(t.test(z + df)$p.value <= alpha){
      powr_t_test_gamma[j] <- 1
    }#power of t test
    #========================================
    #       chisq distn                     #
    # =======================================
    c <-rchisq(n, df = 7) - 7
    k <- c/sqrt(14)
    if(t.test(k)$p.value <= alpha){
      rejectH0_t_test_chisq[j] <- 1
    }
    if(shapiro.test(k)$p.value <= alpha){
      rejectH0_SW_test_chisq[j] <- 1
    }
    if(t.test(k + df)$p.value <= alpha){
      powr_t_test_chisq[j] <- 1
    }#power of t test
    #========================================
    #       weibull distn                     #
    # =======================================
    w1 <- rweibull(n, shape = 7, scale = 1) 
    w <- (w1-gamma(8/7))/sqrt(gamma(9/7)-(gamma(8/7))^2)
    if(t.test(w)$p.value <= alpha){
      rejectH0_t_test_weibull[j] <- 1
    }
    if(shapiro.test(w)$p.value <= alpha){
      rejectH0_SW_test_weibull[j] <- 1
    }
    if(t.test(w + df)$p.value <= alpha){
      powr_t_test_weibull[j] <- 1
    }#power of t test
    
    #========================================
    #       Lognormal distn                     #
    # =======================================
    l <- rlnorm(n, meanlog = 0, sdlog = 0.25)
   sqrt((exp(0.25^2)-1)*exp(2*0 + 0.25^2))
    l1 <- (l- exp(0+0.25^2/2))/sqrt((exp(0.25^2)-1)*exp(2*0 + 0.25^2))
    if(t.test(l1)$p.value <= alpha){
      rejectH0_t_test_logn[j] <- 1
    }
    if(shapiro.test(l1)$p.value <= alpha){
      rejectH0_SW_test_logn[j] <- 1
    }
    if(t.test(l1 + df)$p.value <= alpha){
      powr_t_test_logn[j] <- 1
    }#power of t test
  }
  # Type I error
  Type1_error_Norm[i] <- round(mean(rejectH0_t_test_norm), 3)
  Type1_error_exp[i] <- round(mean(rejectH0_t_test_exp), 3)
  Type1_error_unif[i] <- round(mean(rejectH0_t_test_unif),3)
  Type1_error_beta[i] <- round(mean(rejectH0_t_test_beta), 3)
  Type1_error_gamma[i] <- round(mean(rejectH0_t_test_gamma), 3)
  Type1_error_chisq[i] <- round(mean(rejectH0_t_test_chisq), 3)
  Type1_error_weibull[i] <- round(mean(rejectH0_t_test_weibull), 3)
  Type1_error_logn[i] <- round(mean(rejectH0_t_test_logn), 3)
  #inflation of Type I error
  inflation_Type1_error_exp[i] <- Type1_error_exp[i] - Type1_error_Norm[i]
  inflation_Type1_error_unif[i] <- Type1_error_unif[i] - Type1_error_Norm[i]
  inflation_Type1_error_beta[i] <- Type1_error_beta[i] - Type1_error_Norm[i]
  inflation_Type1_error_gamma[i] <- Type1_error_gamma[i] - Type1_error_Norm[i]
  inflation_Type1_error_chisq[i] <- Type1_error_chisq[i] - Type1_error_Norm[i]
  inflation_Type1_error_weibull[i] <- Type1_error_weibull[i] - Type1_error_Norm[i]
  inflation_Type1_error_logn[i] <- Type1_error_logn[i] - Type1_error_Norm[i]
  #power of SW test
  Power_SW_test_exp[i] <- round(mean(rejectH0_SW_test_exp), 3)
  Power_SW_test_unif[i] <- round(mean(rejectH0_SW_test_unif), 3)
  Power_SW_test_beta[i] <- round(mean(rejectH0_SW_test_beta), 3)
  Power_SW_test_gamma[i] <- round(mean(rejectH0_SW_test_gamma), 3)
  Power_SW_test_chisq[i] <- round(mean(rejectH0_SW_test_chisq), 3)
  Power_SW_test_weibul[i] <- round(mean(rejectH0_SW_test_weibull), 3)
  Power_SW_test_logn[i] <- round(mean(rejectH0_SW_test_logn), 3)
  # Power of t test
  power_t_test_Norm[i] <- round(mean(powr_t_test_norm), 3)
  power_t_test_exp[i] <- round(mean(powr_t_test_exp), 3)
  power_t_test_unif[i] <- round(mean(powr_t_test_unif),3)
  power_t_test_beta[i] <- round(mean(powr_t_test_beta), 3)
  power_t_test_gamma[i] <- round(mean(powr_t_test_gamma), 3)
  power_t_test_chisq[i] <- round(mean(powr_t_test_chisq), 3)
  power_t_test_weibull[i] <- round(mean(powr_t_test_weibull), 3)
  power_t_test_logn[i] <- round(mean(powr_t_test_logn), 3)
  #power loss
  powrloss_exp[i] <- power_t_test_Norm[i] - power_t_test_exp[i]
  powrloss_unif[i] <- power_t_test_Norm[i] - power_t_test_unif[i]
  powrloss_gamma[i] <- power_t_test_Norm[i] - power_t_test_beta[i]
  powrloss_beta[i] <- power_t_test_Norm[i] - power_t_test_gamma[i]
  powrloss_chisq[i] <- power_t_test_Norm[i] - power_t_test_chisq[i]
  powrloss_weibull[i] <- power_t_test_Norm[i] - power_t_test_weibull[i]
  powrloss_logn[i] <- power_t_test_Norm[i] - power_t_test_logn[i]
}

par(mfrow=c(2,2))

plot(sample_size, Power_SW_test_exp, type="l", lwd= 2, col = "red", 
     ylim = c(0, 1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, Power_SW_test_unif, lwd= 2, col = "blue")
lines(sample_size, Power_SW_test_beta,lwd= 2, col = "purple")
lines(sample_size, Power_SW_test_gamma, lwd= 2,col = "green")
lines(sample_size, Power_SW_test_chisq, lwd= 2,col = "pink")
lines(sample_size, Power_SW_test_weibul, lwd= 2, col = "gray")
title(main = "Power SW test.")
legend("topleft", legend=c("exponential", "uniform", "beta", "gamma", "chisq","weibull"),
       col=c("red","blue", "purple", "green", "pink", "gray"), lty =1,cex = 0.5)


plot(sample_size, inflation_Type1_error_unif, type="l",lwd= 2, col = "red", 
     ylim = c(-0.05, 0.1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, inflation_Type1_error_beta, lwd= 2, col = "blue")
lines(sample_size, inflation_Type1_error_gamma,lwd= 2, col = "green")
lines(sample_size, inflation_Type1_error_chisq,lwd= 2, col = "pink")
lines(sample_size, inflation_Type1_error_weibull,lwd= 2, col = "gray")
abline(h=0, lwd=2)
title(main = "Inflation of Type 1 error rate")
legend("topleft", legend=c("expontial", "uniform", "beta", "chisq", "weibull"),
       col=c("red", "blue", "green", "pink", "gray"), lty =1,cex = 0.5)

#Power of t test
plot(sample_size, power_t_test_Norm, type="l", lwd= 2, col = "red", 
     ylim = c(0, 1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, power_t_test_exp, lwd= 2, col = "blue")
lines(sample_size, power_t_test_unif,lwd= 2, col = "gold")
lines(sample_size, power_t_test_beta,lwd= 2, col = "green")
lines(sample_size, power_t_test_gamma, col = "brown")
lines(sample_size, power_t_test_chisq, col = "pink")
lines(sample_size, power_t_test_weibull, col = "gray")
title(main = "Power of t test")
legend("topleft", legend=c("normal", "expontial", "uniform", "beta", "gamma", "weibull"),
       col=c( "red","blue", "gold", "green", "brown", "pink", "gray"), lty =1,cex = 0.5)

#power loss of t test
plot(sample_size, powrloss_exp, type="l", lwd= 2, col = "red", 
     ylim = c(-1, 1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, powrloss_unif,lwd= 2, col = "blue")
lines(sample_size, powrloss_gamma, lwd= 2, col = "green")
lines(sample_size, powrloss_beta, lwd= 2,col = "purple")
lines(sample_size, powrloss_chisq, lwd= 2,col = "pink")
lines(sample_size, powrloss_weibull, lwd= 2,col = "gray")
abline(h = 0)
title(main = "Power loss of t test")
legend("topleft", legend=c("expontial", "uniform", "gamma", "beta", "chisq", "weibull"),
       col=c("red", "blue", "green", "purple", "pink", "gray"), lty =1,cex = 0.5)

dev.off()

save.image(paste0("one_sample",".RData"))










