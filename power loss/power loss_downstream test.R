rm(list = ls())
setwd("D:/OSU/Research_Fall2023/power loss")
set.seed(1)
N <- 100000
alpha <- 0.05
sample_size <- c(10,20,30,40,50)
df <-0.5
power_Norm <- power_unif <- power_exp <- powerloss_unif <-powerloss_exp <- numeric(length(sample_size))
powr_SW.test_unif <- powr_SW.test_exp<- numeric(length(sample_size))
for(i in 1: length(sample_size)){
  n <- sample_size[i]
  print(n)
  #simulate from normal distn
  powr_Norm <- 0
  for( j in 1 : N){
    x1 <- rnorm(n, mean = 0, sd = 1)
    x2 <- rnorm(n, mean = 0, sd = 1) +df
    if(t.test(x1, x2)$p.value <= alpha){
      powr_Norm <- powr_Norm + 1
    }
  }
  #simulate from exp distn
  powr_exp <- 0
  powr_SW_exp <- 0
  for(k in 1 : N) {
    y1 <- rexp(n, 1)
    y2 <- rexp(n, 1)
    if(shapiro.test(y1)$p.value <= alpha){
      powr_SW_exp <- powr_SW_exp + 1
    }#power of SW test
    e1 <- y1 - 1
    e2 <- y2 - 1 + df
    if(t.test(e1,e2)$p.value <= alpha){
      powr_exp <- powr_exp + 1
    }# power of t test
  }
  #simulate from uniform distn
  powr_unif <- 0
  powr_SW_unif <- 0
  for( j in 1 : N){
    w1 <- runif(n, min = 0, max = 1) 
    w2 <- rnorm(n, mean = 0, sd = 1) 
    if(shapiro.test(w1)$p.value <= alpha){
      powr_SW_unif <- powr_SW_unif + 1
    }#power of SW test
    
    u1 <- (w1 - 0.5)/sqrt(1/12) 
    u2 <- (w2 - 0.5)/sqrt(1/12) + df
    if(t.test(u1, u2)$p.value <= alpha){
      powr_unif <- powr_unif + 1
    }# power of t test
    
  }
  #power of SW test
  powr_SW.test_unif[i] <- round(powr_SW_unif/N,3)
  powr_SW.test_exp[i] <- round(powr_SW_exp/N,3)
  #power of t test
  power_Norm[i] <- round(powr_Norm/N, 3)
  power_unif[i] <- round(powr_exp/N, 3)
  power_exp[i] <- round(powr_unif/N, 3)
  #power loss
  powerloss_unif[i] <- power_Norm[i] - power_unif[i]
  powerloss_exp[i] <- power_Norm[i] - power_exp[i]
}

png(file="power_SW.jpeg", width=600, height=600)

# par(mfrow=c(2,2))

plot(sample_size, powr_SW.test_exp, type="l", lwd=2, col = "red", 
 ylim = c(0, 1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, powr_SW.test_unif, lwd=2,col = "blue")
title(main = "Power SW test.")
legend("topleft", legend=c("power SW for unif", "power SW for exp"),
       col=c("red","blue"), lty =1,cex = 0.5)

dev.off()
plot(sample_size, power_Norm, type="l", lwd=2,col = "red", 
     ylim = c(0, 1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, power_unif,lwd=2, col = "blue")
lines(sample_size, power_exp,lwd=2, col = "black")
title(main = "Power of two samples t test.")
legend("topleft", legend=c("power_Norm", "power_unif", "power_exp"),
       col=c("red","blue", "black"), lty =1,cex = 0.5)


plot(sample_size, powerloss_unif, type="l",lwd=2, col = "red", 
     ylim = c(-0.1, 0.1),xlab = "Sample Size", ylab = "Power loss")
lines(sample_size, powerloss_exp, lwd=2,col = "blue")
abline(h=0.0, col="black", lwd=2)
title(main = "Power loss of uniform and exponential distns")
legend("topleft", legend=c("powerloss_unif", "powerloss_exp"),
       col=c("red","blue"), lty =1,cex = 0.5)

dev.off()

save.image(paste0("power_loss",".RData"))

