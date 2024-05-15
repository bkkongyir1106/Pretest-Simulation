rm(list = ls())
setwd("D:/OSU/Research_Fall2023/One sample t test")
set.seed(1)
N <- 10000
alpha <- 0.05
sample_size <- c(10,20,30,40,50)
Type1_error_Norm <- Type1_error_exp <- Type1_error_unif <- numeric(length(sample_size))
inflation_Type1_error_exp <- inflation_Type1_error_unif <- inflation_Type1_error_beta <- numeric(length(sample_size))
Power_SW_test_exp <- Power_SW_test_unif <- Power_SW_test_beta <- numeric(length(sample_size))

for (i in 1 : length(sample_size)) {
  n <- sample_size[i]
  print(n)
  
  rejectH0_SW_test_exp <- rejectH0_SW_test_unif <- rejectH0_SW_test_beta <- numeric(N)
  rejectH0_t_test_exp <- rejectH0_t_test_unif <- rejectH0_t_test_beta <- rejectH0_t_test_norm <- numeric(N)
  for(j in 1 : N){
    x <- rnorm(n, mean =0, sd = 1)
    if(t.test(x)$p.value <= alpha){
      rejectH0_t_test_norm[j] <- 1
    }
    e <- rexp(n, rate = 1) - 1
    if(t.test(e)$p.value <= alpha){
      rejectH0_t_test_exp[j] <- 1
    }
    if(shapiro.test(e)$p.value <= alpha){
      rejectH0_SW_test_exp[j] <- 1
    }
    u <- runif(n, min = 0, max = 1) - 0.5
    y <- u/sqrt(1/12)
    if(t.test(y)$p.value <= alpha){
      rejectH0_t_test_unif[j] <- 1
    }
    if(shapiro.test(y)$p.value <= alpha){
      rejectH0_SW_test_unif[j] <- 1
    }
  }
  # Type I error
  Type1_error_Norm[i] <- round(mean(rejectH0_t_test_norm), 3)
  Type1_error_exp[i] <- round(mean(rejectH0_t_test_exp), 3)
  Type1_error_unif[i] <- round(mean(rejectH0_t_test_unif),3)
  inflation_Type1_error_exp[i] <- Type1_error_exp[i] - Type1_error_Norm[i]
  inflation_Type1_error_unif[i] <- Type1_error_unif[i] - Type1_error_Norm[i]
  
  #power of SW test
  Power_SW_test_exp[i] <- round(mean(rejectH0_SW_test_exp), 3)
  Power_SW_test_unif[i] <- round(mean(rejectH0_SW_test_unif), 3)
}

par(mfrow=c(2,2))

plot(sample_size, Power_SW_test_exp, type="l", col = "red", 
     ylim = c(0, 1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, Power_SW_test_unif, col = "blue")
title(main = "Power SW test.")
legend("topleft", legend=c("power SW for exp", "power SW for unif"),
       col=c("red","blue"), lty =1:2,cex = 0.5)


plot(sample_size, inflation_Type1_error_exp, type="l", col = "red", 
     ylim = c(-0.02, 0.1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, inflation_Type1_error_unif, col = "blue")
title(main = "Inflation of Type 1 error rate")
legend("topleft", legend=c("expontial", "uniform"),
       col=c("blue", "red"), lty =1:2,cex = 0.5)

dev.off()

save.image(paste0("inflation_error",".RData"))










