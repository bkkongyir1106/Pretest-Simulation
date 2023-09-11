rm(list = ls())
setwd("D:/OSU/Research_Fall2023/One sample t test")
set.seed(1)
N <- 10000
alpha <- 0.05
df <- 0.5
sample_size <- c(10,20,30,40,50)
Type1_error_Norm <- Type1_error_exp <- Type1_error_unif <- Type1_error_beta <- Type1_error_gamma <- Type1_error_chisq <-Type1_error_weibull <- numeric(length(sample_size))
inflation_Type1_error_exp <- inflation_Type1_error_unif <- inflation_Type1_error_beta <- inflation_Type1_error_gamma <- inflation_Type1_error_chisq <-inflation_Type1_error_weibull <- numeric(length(sample_size))
Power_SW_test_exp <- Power_SW_test_unif <- Power_SW_test_beta <- Power_SW_test_gamma <- Power_SW_test_chisq <-Power_SW_test_weibul <- numeric(length(sample_size))
power_t_test_Norm <- power_t_test_exp <- power_t_test_unif <- power_t_test_gamma <- power_t_test_beta <- power_t_test_chisq <-power_t_test_weibull <- numeric(length(sample_size))
powrloss_exp <-powrloss_unif  <-powrloss_gamma <-powrloss_beta <-powrloss_chisq <-powrloss_weibull <- numeric(length(sample_size))
for (i in 1 : length(sample_size)) {
  n <- sample_size[i]
  print(n)
  
  rejectH0_SW_test_exp <- rejectH0_SW_test_unif <- rejectH0_SW_test_beta <- rejectH0_SW_test_gamma <- rejectH0_SW_test_chisq <- rejectH0_SW_test_weibull<- numeric(N)
  rejectH0_t_test_exp <- rejectH0_t_test_unif <- rejectH0_t_test_beta <- rejectH0_t_test_norm <- rejectH0_t_test_gamma <- rejectH0_t_test_chisq <- rejectH0_t_test_weibull <- numeric(N)
  powr_t_test_norm <- powr_t_test_exp <- powr_t_test_unif <- powr_t_test_beta <- powr_t_test_gamma <- powr_t_test_chisq <- powr_t_test_weibull <-numeric(N)
  for(j in 1 : N){
    #========================================
    #       Normal distn                   #
    # =======================================
    x1 <- rnorm(n, mean =0, sd = 1)
    x2 <- rnorm(n, mean = 0, sd = 1)
    if(t.test(x1, x2)$p.value <= alpha){
      rejectH0_t_test_norm[j] <- 1
    }#Type I error
    if(t.test(x1, x2 + df)$p.value <= alpha){
      powr_t_test_norm[j] <- 1
    }#power of t test
    #========================================
    #       exponential distn               #
    # =======================================
    e1 <- rexp(n, rate = 1)
    e2 <- rexp(n, rate = 1)
    if(t.test(e1, e2)$p.value <= alpha){
      rejectH0_t_test_exp[j] <- 1
    }#Type I error
    if(shapiro.test(e1)$p.value <= alpha/2 | shapiro.test(e2)$p.value <= alpha/2){
      rejectH0_SW_test_exp[j] <- 1
    }#power of SW test
    if(t.test(e1 -1, e2 - 1 + df)$p.value <= alpha){
      powr_t_test_exp[j] <- 1
    }#power of t test
    #========================================
    #       uniform distn                   #
    # =======================================
    y1 <- runif(n, min = 0, max = 1)
    y2 <- runif(n, min = 0, max = 1)
    if(t.test(y1, y2)$p.value <= alpha){
      rejectH0_t_test_unif[j] <- 1
    }# Type I error 
    if(shapiro.test(y1)$p.value <= alpha/2 | shapiro.test(y1)$p.value <= alpha/2){
      rejectH0_SW_test_unif[j] <- 1
    }# power of SW test
    if(t.test((y1-0.5)/sqrt(1/12), (y2-0.5)/sqrt(1/12) + df)$p.value <= alpha){
      powr_t_test_unif[j] <- 1
    }#power of t test
    #========================================
    #       beta distn                      #
    # =======================================
    b1 <- rbeta(n, shape1 = 2, shape2 = 1) 
    b2 <- rbeta(n, shape1 = 2, shape2 = 1) 
    if(t.test(b1,b2)$p.value <= alpha){
      rejectH0_t_test_beta[j] <- 1
    }# Type I error 
    if(shapiro.test(b1)$p.value <= alpha/2 | shapiro.test(b2)$p.value <= alpha/2){
      rejectH0_SW_test_beta[j] <- 1
    }#power of SW test
    if(t.test((b1-2/3)/sqrt(1/18), (b2-2/3)/sqrt(1/18) + df)$p.value <= alpha){
      powr_t_test_beta[j] <- 1
    }#power of t test
    
    #========================================
    #       gamma distn                     #
    # =======================================
    g1 <- rgamma(n, shape = 3, rate = 1) 
    g2 <- rgamma(n, shape = 3, rate = 1) 
    if(t.test(g1, g2)$p.value <= alpha){
      rejectH0_t_test_gamma[j] <- 1
    }# Type I error 
    if(shapiro.test(g1)$p.value <= alpha/2 | shapiro.test(g2)$p.value <= alpha/2){
      rejectH0_SW_test_gamma[j] <- 1
    }#power of SW test
    if(t.test((g1-3)/sqrt(3), (g2-3)/sqrt(3) + df)$p.value <= alpha){
      powr_t_test_gamma[j] <- 1
    }#power of t test
    #========================================
    #       chisq distn                     #
    # =======================================
    k1 <-rchisq(n, df = 4) 
    k2 <-rchisq(n, df = 4) 
    if(t.test(k1, k2)$p.value <= alpha){
      rejectH0_t_test_chisq[j] <- 1
    }
    if(shapiro.test(k1)$p.value <= alpha/2 | shapiro.test(k2)$p.value <= alpha/2){
      rejectH0_SW_test_chisq[j] <- 1
    }
    if(t.test((k1 - 4)/sqrt(8), (k2 - 4)/sqrt(8) + df)$p.value <= alpha){
      powr_t_test_chisq[j] <- 1
    }#power of t test
    #========================================
    #       Weibull distn                     #
    # =======================================
    w1 <- rweibull(n, shape = 7, scale = 1) 
    w2 <- rweibull(n, shape = 7, scale = 1) 
    if(t.test(w1, w2)$p.value <= alpha){
      rejectH0_t_test_weibull[j] <- 1
    }# Type I error 
    if(shapiro.test(w1)$p.value <= alpha/2 | shapiro.test(w2)$p.value <= alpha/2){
      rejectH0_SW_test_weibull[j] <- 1
    }#power of SW test
    if(t.test((w1-gamma(8/7))/sqrt(gamma(9/7)-(gamma(8/7))^2), (w2-gamma(8/7))/sqrt(gamma(9/7)-(gamma(8/7))^2) + df)$p.value <= alpha){
      powr_t_test_weibull[j] <- 1
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
  #inflation of Type I error
  inflation_Type1_error_exp[i] <- Type1_error_exp[i] - Type1_error_Norm[i]
  inflation_Type1_error_unif[i] <- Type1_error_unif[i] - Type1_error_Norm[i]
  inflation_Type1_error_beta[i] <- Type1_error_beta[i] - Type1_error_Norm[i]
  inflation_Type1_error_gamma[i] <- Type1_error_gamma[i] - Type1_error_Norm[i]
  inflation_Type1_error_chisq[i] <- Type1_error_chisq[i] - Type1_error_Norm[i]
  inflation_Type1_error_weibull[i] <- Type1_error_weibull[i] - Type1_error_Norm[i]
  #power of SW test
  Power_SW_test_exp[i] <- round(mean(rejectH0_SW_test_exp), 3)
  Power_SW_test_unif[i] <- round(mean(rejectH0_SW_test_unif), 3)
  Power_SW_test_beta[i] <- round(mean(rejectH0_SW_test_beta), 3)
  Power_SW_test_gamma[i] <- round(mean(rejectH0_SW_test_gamma), 3)
  Power_SW_test_chisq[i] <- round(mean(rejectH0_SW_test_chisq), 3)
  Power_SW_test_weibul[i] <- round(mean(rejectH0_SW_test_weibull), 3)
  # Power of t test
  power_t_test_Norm[i] <- round(mean(powr_t_test_norm), 3)
  power_t_test_exp[i] <- round(mean(powr_t_test_exp), 3)
  power_t_test_unif[i] <- round(mean(powr_t_test_unif),3)
  power_t_test_beta[i] <- round(mean(powr_t_test_beta), 3)
  power_t_test_gamma[i] <- round(mean(powr_t_test_gamma), 3)
  power_t_test_chisq[i] <- round(mean(powr_t_test_chisq), 3)
  power_t_test_weibull[i] <- round(mean(powr_t_test_weibull), 3)
  #power loss
  powrloss_exp[i] <- power_t_test_Norm[i] - power_t_test_exp[i]
  powrloss_unif[i] <- power_t_test_Norm[i] - power_t_test_unif[i]
  powrloss_gamma[i] <- power_t_test_Norm[i] - power_t_test_beta[i]
  powrloss_beta[i] <- power_t_test_Norm[i] - power_t_test_gamma[i]
  powrloss_chisq[i] <- power_t_test_Norm[i] - power_t_test_chisq[i]
  powrloss_weibull[i] <- power_t_test_Norm[i] - power_t_test_weibull[i]
}

par(mfrow=c(2,2))
#power of SW test
plot(sample_size, Power_SW_test_exp, type="l", lwd= 2, col = "red", 
     ylim = c(0, 1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, Power_SW_test_unif, lwd= 2, col = "blue")
lines(sample_size, Power_SW_test_beta,lwd= 2, col = "purple")
lines(sample_size, Power_SW_test_gamma, lwd= 2,col = "green")
lines(sample_size, Power_SW_test_chisq, lwd= 2,col = "pink")
lines(sample_size, Power_SW_test_weibul, lwd= 2,col = "violet")
title(main = "Power SW test.")
legend("topleft", legend=c("exponential", "uniform", "beta", "gamma", "chisq", "Weibull"),
       col=c("red","blue", "purple", "green", "pink", "violet"), lty =1,cex = 0.5)
#Type I error rate
plot(sample_size, Type1_error_Norm, type="l", lwd= 2, col = "red", 
     ylim = c(0, 0.1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, Type1_error_exp, lwd= 2, col = "blue")
lines(sample_size, Type1_error_unif,lwd= 2, col = "purple")
lines(sample_size, Type1_error_beta, lwd= 2,col = "green")
lines(sample_size, Type1_error_gamma, lwd= 2,col = "pink")
lines(sample_size, Type1_error_chisq, lwd= 2,col = "cyan")
lines(sample_size, Type1_error_weibull, lwd= 2,col = "violet")
abline(h = 0.05)
title(main = "Type I error of t test.")
legend("topleft", legend=c("normal", "exponential", "uniform", "beta", "gamma", "chisq", "weibull"),
       col=c("red","blue", "purple", "green", "pink", "cyan", "violet"), lty =1,cex = 0.5)

# inflation of Type I error rate
plot(sample_size, inflation_Type1_error_unif, type="l",lwd= 2, col = "red", 
     ylim = c(-0.05, 0.1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, inflation_Type1_error_beta, lwd= 2, col = "blue")
lines(sample_size, inflation_Type1_error_gamma,lwd= 2, col = "green")
lines(sample_size, inflation_Type1_error_chisq,lwd= 2, col = "pink")
lines(sample_size, inflation_Type1_error_weibull,lwd= 2, col = "violet")
abline(h=0)
title(main = "Inflation of Type 1 error rate")
legend("topleft", legend=c("expontial", "uniform", "beta", "chisq", "weibull"),
       col=c("red", "blue", "green", "pink", "violet"), lty =1,cex = 0.5)

#Power of t test
plot(sample_size, power_t_test_Norm, type="l", lwd= 2, col = "red", 
     ylim = c(0, 1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, power_t_test_exp, lwd= 2, col = "blue")
lines(sample_size, power_t_test_unif,lwd= 2, col = "gold")
lines(sample_size, power_t_test_beta,lwd= 2, col = "green")
lines(sample_size, power_t_test_gamma, col = "brown")
lines(sample_size, power_t_test_chisq, col = "pink")
lines(sample_size, power_t_test_weibull, col = "violet")
title(main = "Power of t test")
legend("topleft", legend=c("normal", "expontial", "uniform", "beta", "gamma"),
       col=c( "red","blue", "gold", "green", "brown", "pink", "violet"), lty =1,cex = 0.5)

#power loss of t test
plot(sample_size, powrloss_exp, type="l", lwd= 2, col = "red", 
     ylim = c(-1, 1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, powrloss_unif,lwd= 2, col = "blue")
lines(sample_size, powrloss_gamma, lwd= 2, col = "green")
lines(sample_size, powrloss_beta, lwd= 2,col = "purple")
lines(sample_size, powrloss_chisq, lwd= 2,col = "pink")
lines(sample_size, powrloss_weibull, lwd= 2,col = "violet")
abline(h = 0)
title(main = "Power loss of t test")
legend("topleft", legend=c("expontial", "uniform", "gamma", "beta", "chisq", "weibull"),
       col=c("red", "blue", "green", "purple", "violet"), lty =1,cex = 0.5)

dev.off()

save.image(paste0("two_samples",".RData"))




