setwd("D:/OSU/Research_Fall2023/power loss")
set.seed(1)
N <- 10000
alpha <- 0.05
sample_size <- c(10,20,30,40,50)
Type1_error_Norm <- Type1_error_unif <- Type1_error_exp <- numeric(length(sample_size))
power_SW_exp <- power_SW_unif <- Inflation_TypeI_error_unif <-Inflation_TypeI_error_exp <- numeric(length(sample_size))

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
  powr_SW_exp <- 0
  for(k in 1 : N) {
    y1 <- rexp(n, 1)
    y2 <- rexp(n, 1)
    if(t.test(y1,y2)$p.value <= alpha){
      error_exp <- error_exp + 1
    }
    #power of SW test for exp distn
    if(shapiro.test(y1)$p.value <= alpha){
      powr_SW_exp <- powr_SW_exp + 1
    }
  }
  #simulate from uniform distn
  error_unif <- 0
  powr_SW_unif <- 0
  for( j in 1 : N){
    u1 <- runif(n, min = 0, max = 1)
    u2 <- runif(n, min = 0, max = 1)
    if(t.test(u1, u2)$p.value <= alpha){
      error_unif <- error_unif + 1
    }
    #power of SW test for unif distn
    if(shapiro.test(u1)$p.value <= alpha){
      powr_SW_unif <- powr_SW_unif + 1
    }
  }
  
  power_SW_exp[i] <- round(powr_SW_exp/N, 3)
  power_SW_unif[i] <- round(powr_SW_unif/N, 3)
  
  Type1_error_Norm[i] <- round(error_Norm/N, 3)
  Type1_error_unif[i] <- round(error_exp/N, 3)
  Type1_error_exp[i] <- round(error_unif/N, 3)

  Inflation_TypeI_error_unif[i] <- Type1_error_unif[i] - Type1_error_Norm[i]
  Inflation_TypeI_error_exp[i] <- Type1_error_exp[i] - Type1_error_Norm[i]
}


plot(sample_size, power_SW_exp, type="l",lwd=2, col = "red", 
     ylim = c(0, 1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, power_SW_unif,lwd=2, col = "blue")
title(main = "Power SW test.")
legend("topleft", legend=c("power SW for exp", "power SW for unif"),
       col=c("red","blue"), lty =1:2,cex = 0.5)

plot(sample_size, Type1_error_Norm, type="l", lwd=2, col = "red", 
     ylim = c(0, 0.1),xlab = "Sample Size", ylab = "Inflation of Type I error")
lines(sample_size, Type1_error_unif,lwd=2, col = "blue")
lines(sample_size, Type1_error_exp,lwd=2, col = "brown")
abline(h=0.05, lwd=2,col="black")
title(main = "Type I error Rates.")
legend("topleft", legend=c("Type1_error_Norm", "Type1_error_unif", "Type1_error_exp"),
       col=c("red","blue", "brown"), lty =1,cex = 0.5)

plot(sample_size, Inflation_TypeI_error_unif, type="l",lwd=2, col = "red", 
     ylim = c(-0.05, 0.05),xlab = "Sample Size", ylab = "Inflation of Type I error")
lines(sample_size, Inflation_TypeI_error_exp,lwd=2, col = "blue")
abline(h=0.0, lwd=2,col="black")
title(main = "Inflation of Type I error.")
legend("topleft", legend=c("inflation of Type I error--unif", "inflation of Type I error--exp"),
       col=c("red","blue"), lty =1:2,cex = 0.5)

save.image(paste0("Type1_error",".RData"))
