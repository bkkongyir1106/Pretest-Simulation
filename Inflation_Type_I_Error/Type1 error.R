setwd("D:/OSU/Research_Fall2023/power loss")
set.seed(1)
N <- 100000
alpha <- 0.05
sample_size <- c(5,10,20,30,40,50,100)
Type1_error_Norm <- Type1_error_unif <- Type1_error_exp <- Inflation_TypeI_error_unif <-Inflation_TypeI_error_exp <- numeric(length(sample_size))

for(i in 1: length(sample_size)){
  n <- sample_size[i]
  print(n)
  #simulate from normal distn
  error_Norm <- 0
  for( j in 1 : N){
    X1 <- rnorm(n, mean = 0, sd = 1)
    x2 <- rnorm(n, mean = 0, sd = 1)
    if(t.test(x1, x2)$p.value <= alpha){
      error_Norm <- error_Norm + 1
    }
  }
  #simulate from exp distn
  error_exp <- 0
  for(k in 1 : N) {
    y1 <- rexp(n, 1)
    y2 <- rexp(n, 1)
    if(t.test(y1,y2)$p.value <= alpha){
      error_exp <- error_exp + 1
    }
  }
  #simulate from uniform distn
  
  error_unif <- 0
  for( j in 1 : N){
    u1 <- runif(n, min = 0, max = 1)
    u2 <- runif(n, min = 0, max = 1)
    if(t.test(u1, u2)$p.value <= alpha){
      error_unif <- error_unif + 1
    }
  }
  Type1_error_Norm[i] <- round(error_Norm/N, 3)
  Type1_error_unif[i] <- round(error_exp/N, 3)
  Type1_error_exp[i] <- round(error_unif/N, 3)

  Inflation_TypeI_error_unif[i] <- Type1_error_unif[i] - Type1_error_Norm[i]
  Inflation_TypeI_error_exp[i] <- Type1_error_exp[i] - Type1_error_Norm[i]
}
Inflation_TypeI_error_unif
Inflation_TypeI_error_exp
png(file="error_inflation.jpeg", width=600, height=600)
plot(sample_size, Type1_error_Norm, type="l", col = 1, 
     ylim = c(-0.03, 0.08),xlab = "Sample Size", ylab = "Type I error")
lines(sample_size, Type1_error_unif, col = 2)
lines(sample_size, Type1_error_exp, col = 3)
lines(sample_size, Inflation_TypeI_error_unif, col = 4)
lines(sample_size, Inflation_TypeI_error_unif, col = 5)
title(main = "Type I Error Rate.")
legend("topleft", legend=c("Type1_error_Norm", "Type1_error_unif", "Type1_error_exp", "Inflation_TypeI_error_unif", "Inflation_TypeI_error_unif"),col=1:5, 
       lty =1:5,cex = 0.5)
dev.off()

save.image(paste0("Type1_error",".RData"))
