setwd("D:/OSU/Research_Fall2023/power loss")
set.seed(1)
N <- 100000
alpha <- 0.05
sample_size <- c(10,20,30,40,50)
Type1_error_Norm <- Type1_error_unif <- Type1_error_exp <- Inflation_TypeI_error_unif <-Inflation_TypeI_error_exp <- numeric(length(sample_size))
power_SW_exp1 <- power_SW_unif1 <- power_SW_exp_1n2 <- power_SW_unif_1n2 <- numeric(length(sample_size))
for(i in 1: length(sample_size)){
  n <- sample_size[i]
  print(n)
  #simulate from normal distn
  error_Norm <- 0
  for( j in 1 : N){
    x1 <- rnorm(n, mean = 0, sd = 1)
    x2 <- rnorm(n, mean = 0, sd = 1)
    if(t.test(x1, x2)$p.value <= alpha){
      error_Norm <- error_Norm + 1
    }
  }
  #simulate from exp distn
  error_exp <- 0
  powr_SW_exp1 <- 0
  powr_SW_exp2 <- 0
  for(k in 1 : N) {
    y1 <- rexp(n, 1)
    y2 <- rexp(n, 1)
    if(t.test(y1,y2)$p.value <= alpha){
      error_exp <- error_exp + 1
    }
    #power of SW test for exp distn
    if(shapiro.test(y1)$p.value <= alpha){
      powr_SW_exp1 <- powr_SW_exp1 + 1
    }
    if(shapiro.test(y2)$p.value <= alpha){
      powr_SW_exp2 <- powr_SW_exp2 + 1
    }
  }
  #simulate from uniform distn
  
  error_unif <- 0
  powr_SW_unif1 <- 0
  powr_SW_unif2 <- 0
  for( j in 1 : N){
    u1 <- runif(n, min = 0, max = 1)
    u2 <- runif(n, min = 0, max = 1)
    if(t.test(u1, u2)$p.value <= alpha){
      error_unif <- error_unif + 1
    }
    #power of SW test for unif distn
    if(shapiro.test(u1)$p.value <= alpha){
      powr_SW_unif1 <- powr_SW_unif1 + 1
    }
    if(shapiro.test(u2)$p.value <= alpha){
      powr_SW_unif2 <- powr_SW_unif2 + 1
    }
  }
  #power of SW for sample 1
  power_SW_exp1[i] <- round(powr_SW_exp1/N, 3)
  power_SW_unif1[i] <- round(powr_SW_unif1/N, 3)
  
  #power of SW for sample 1 & 2
  power_SW_exp_1n2[i] <- round((powr_SW_exp1 + powr_SW_exp1)/(2*N), 3)
  power_SW_unif_1n2[i] <- round((powr_SW_unif1 + powr_SW_unif2)/(2*N), 3)
  
  Type1_error_Norm[i] <- round(error_Norm/N, 3)
  Type1_error_unif[i] <- round(error_exp/N, 3)
  Type1_error_exp[i] <- round(error_unif/N, 3)
  
  Inflation_TypeI_error_unif[i] <- Type1_error_unif[i] - Type1_error_Norm[i]
  Inflation_TypeI_error_exp[i] <- Type1_error_exp[i] - Type1_error_Norm[i]
}

par(mfrow=c(2,2))

plot(sample_size, power_SW_exp_1n2, type="l", col = "red", 
     ylim = c(0, 1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, power_SW_unif_1n2, col = "blue")
title(main = "Power SW test.")
legend("topleft", legend=c("power SW for exp", "power SW for unif"),
       col=c("red","blue"), lty =1:2,cex = 0.5)

plot(sample_size, Inflation_TypeI_error_unif, type="l", col = "red", 
     ylim = c(-0.05, 0.05),xlab = "Sample Size", ylab = "Inflation of Type I error")
lines(sample_size, Inflation_TypeI_error_exp, col = "blue")
title(main = "Inflation of Type I error")
legend("topleft", legend=c("inflation of Type I error--unif", "inflation of Type I error--exp"),
       col=c("red","blue"), lty =1:2,cex = 0.5)

dev.off()

save.image(paste0("Inflation_Type1_error",".RData"))
