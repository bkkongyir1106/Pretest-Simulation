setwd("D:/OSU/Research_Fall2023/power loss")
set.seed(1)
N <- 100000
alpha <- 0.05
sample_size <- c(5,10,20,30,40,50,100)
df <-0.5
power_Norm <- power_unif <- power_exp <- powerloss_unif <-powerloss_exp <- numeric(length(sample_size))

for(i in 1: length(sample_size)){
  n <- sample_size[i]
  print(n)
  #simulate from normal distn
  powr_Norm <- 0
  for( j in 1 : N){
    X1 <- rnorm(n, mean = 0, sd = 1)
    x2 <- rnorm(n, mean = 0, sd = 1) +df
    if(t.test(x1, x2)$p.value <= alpha){
      powr_Norm <- powr_Norm + 1
    }
  }
  #simulate from exp distn
  powr_exp <- 0
  for(k in 1 : N) {
    y1 <- rexp(n, 1) - 1
    y2 <- rexp(n, 1)-1 + df
    if(t.test(y1,y2)$p.value <= alpha){
      powr_exp <- powr_exp + 1
    }
  }
  #simulate from uniform distn
  
  powr_unif <- 0
  for( j in 1 : N){
    u1 <- runif(n, min = 0, max = 1) - 0.5
    u3 <- rnorm(n, mean = 0, sd = 1) - 0.5
    u2 <- u3/sqrt(1/12) + df
    if(t.test(u1, u2)$p.value <= alpha){
      powr_unif <- powr_unif + 1
    }
  }
  power_Norm[i] <- round(powr_Norm/N, 3)
  power_unif[i] <- round(powr_exp/N, 3)
  power_exp[i] <- round(powr_unif/N, 3)
  
  powerloss_unif[i] <- power_Norm[i] - power_unif[i]
  powerloss_exp[i] <- power_Norm[i] - power_exp[i]
}
powerloss_unif
powerloss_exp
png(file="power_loss.jpeg", width=600, height=600)
plot(sample_size, power_Norm, type="l", col = 1, 
     ylim = c(-0.5, 1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, power_unif, col = 2)
lines(sample_size, power_exp, col = 3)
lines(sample_size, powerloss_unif, col = 4)
lines(sample_size, powerloss_exp, col = 5)
title(main = "Power of downstream t test.")
legend("topleft", legend=c("power_Norm", "power_unif", "power_exp", "powerloss_unif", "powerloss_exp"),col=1:5, 
       lty =1:5,cex = 0.5)
dev.off()

save.image(paste0("power_loss",".RData"))

