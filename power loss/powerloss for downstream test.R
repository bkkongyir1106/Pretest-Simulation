setwd("D:/OSU/Research_Fall2023/power loss")
set.seed(24)
N <- 100000
alpha <- 0.05
sample_size <- c(10,20,30,40,50)
df <-0.5
Type1_error<- power_SW_unif <- power_SW_exp <-power_norm <- power_unif <- power_exp <- power_loss_unif <- power_loss_exp <-  numeric(length(sample_size))

# simulate from normal dist
for(i in 1 : length(sample_size)){
  n <- sample_size[i]
  print(n)
  error1 <- 0
  error2 <- 0
  Norm_rejecth0 <- 0
  for(k in 1 : N) {
    x1 <- rnorm(n, mean = 0, sd = 1)
    x2 <- rnorm(n, mean = 0, sd = 1) + df
    if(shapiro.test(x1)$p.value <= alpha){
      error1 <- error1 + 1
    }
    if(shapiro.test(x2)$p.value <= alpha){
      error2 <- error2 + 1
    }
    error <- error1 + error2
    
    if(t.test(x1,x2)$p.value <= alpha){
      Norm_rejecth0 <- Norm_rejecth0 + 1
    }
  }
  # simulate from exp dist
  power_SW_test_exp1 <- 0
  power_SW_test_exp2 <- 0
  exp_rejecth0 <- 0
  for(k in 1 : N) {
    y1 <- rexp(n, 1)
    y2 <- rexp(n, 1)-1 + df
    if(shapiro.test(y1)$p.value <= alpha){
      power_SW_test_exp1 <- power_SW_test_exp1 + 1
    }
    if(shapiro.test(y2)$p.value <= alpha){
      power_SW_test_exp2 <- power_SW_test_exp2 + 1
    }
    power_SW_test_exp <- power_SW_test_exp1 + power_SW_test_exp2
    
    if(t.test(y1,y2)$p.value <= alpha){
      exp_rejecth0 <- exp_rejecth0 + 1
    }
  }
  # simulate from unif dist
  power_S.test1 <- 0
  power_S.test2 <- 0
  unif_rejecth0 <-0
  for(k in 1 : N) {
    u01 <- runif(n, min = 0, max = 1)-0.5
    u02 <- runif(n, min = 0, max = 1) -0.5 + df
    u1 <-u01/sqrt(1/12)
    u2 <- u02/sqrt(1/12)
    if(shapiro.test(u1)$p.value <= alpha){
      power_S.test1 <- power_S.test1 + 1
    }
    if(shapiro.test(u2)$p.value <= alpha){
      power_S.test2 <- power_S.test2 + 1
    }
    power_SW_test <- power_S.test1 + power_S.test2
    
    if(t.test(u1,u2)$p.value <= alpha){
      unif_rejecth0 <- unif_rejecth0 + 1
    }
  }
  Type1_error[i] <- round(error/(2*N), 3) 
  power_SW_unif[i] <- round(power_SW_test/(2*N), 3) 
  power_SW_exp[i] <- round(power_SW_test_exp/(2*N), 3) 
  power_norm[i] <- round(Norm_rejecth0/N, 3) 
  power_unif[i] <- round(unif_rejecth0/N, 3)
  power_exp[i] <- round(exp_rejecth0/N, 3)
  power_loss_unif[i] <- abs(power_norm[i] - power_unif[i])
  power_loss_exp[i] <- abs(power_norm[i] - power_exp[i])
}

png(file="powerloss.jpeg", width=600, height=600)
plot(sample_size, power_SW, type="l", col = 1, 
     ylim = c(0, 1),xlab = "Sample Size", ylab = "Power")
lines(sample_size, power_loss_unif, col = 2)
lines(sample_size, power_loss_exp, col = 3)
title(main = "Power of pretest and downstream test.")
legend("topleft", legend=c("power of SW test","Power loss for unif", "power loss for exp"),col=1:3, 
       lty =c(1:3),cex = 0.5)
dev.off()

save.image(paste0("powerloss",".RData"))
