setwd("D:/OSU/Research_Fall2023/power loss")
set.seed(1)
N <- 100000
alpha <- 0.05
sample_size <- c(10,20,30,40,50)
df <-1
Type1_error1 <- Type1_error2 <- power_pretest1 <- power_pretest2 <- numeric(length(sample_size))
power_norm <- power_unif <- power_loss <- numeric(length(sample_size))
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
      
      if(t.test(x1,x2)$p.value <= alpha){
        Norm_rejecth0 <- Norm_rejecth0 + 1
      }
    }
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
      if(t.test(u1,u2)$p.value <= alpha){
        unif_rejecth0 <- unif_rejecth0 + 1
      }
    }
    Type1_error1[i] <- round(error1/N, 3) 
    Type1_error2[i] <- round(error2/N, 3)
    power_pretest1[i] <- round(power_S.test1/N, 3) 
    power_pretest2[i] <- round(power_S.test2/N, 3) 
    power_norm[i] <- round(Norm_rejecth0/N, 3) 
    power_unif[i] <- round(unif_rejecth0/N, 3)
    power_loss[i] <-power_norm[i] - power_unif[i]
  }

png(file="Type I error.jpeg", width=600, height=600)
plot(sample_size, Type1_error1, type="l", col = 1, 
     ylim = c(-0.01, 1),xlab = "Sample Size", ylab = "Type I error")
lines(sample_size, Type1_error2, col = 2)
title(main = "Type I error of Pretest.")
legend("bottomright", legend=c("Type1_error1","Type1_error2"),col=1:2, 
       lty =c(1:2),cex = 0.5)
dev.off()
png(file="Power of SW test.jpeg", width=600, height=600)
plot(sample_size, power_pretest1, type="l", col = 1, 
     ylim = c(-0.01, 1),xlab = "Sample Size", ylab = "Power of SW test")
lines(sample_size, power_pretest2, col = 2)
title(main = "Power of SW test.")
legend("bottomright", legend=c("power_pretest1", "power_pretest2"),col=1:2, 
       lty =c(1:2),cex = 0.5)
dev.off()
png(file="Power of downstream test.jpeg", width=600, height=600)
plot(sample_size, power_norm, type="l", col = 1, 
     ylim = c(-0.5, 1),xlab = "Sample Size", ylab = "power of t test")
lines(sample_size, power_unif, col = 2)
lines(sample_size, power_loss, col= 3)
title(main = "Power of downstream test.")
legend("bottomright", legend=c("power_norm", "power_unif", "power_loss"),col=1:3, 
       lty =c(1:3),cex = 0.5)

# plot(sample_size, Type1_error1, type="l", col = 1, 
#      ylim = c(0, 1),xlab = "effect Size", ylab = "Power loss")
# lines(sample_size, Type1_error2, col = 2)
# lines(sample_size, power_pretest1, col = 3)
# lines(sample_size, power_pretest2, col = 4)
# lines(sample_size, power_norm,  col = 5)
# lines(sample_size, power_unif,  col = 6)
# lines(sample_size, power_loss,  col = 7)
# title(main = "Power loss under t test for uniform dist.")
# legend("bottomright", legend=c("Type1_error1","Type1_error2","power_pretest1", "power_pretest2","power_norm", "power_unif", "power_loss" ),col=1:7, 
#        lty =c(1:6),cex = 0.5)
dev.off()
#install.packages("writexl")
library("writexl")

powerloss_unif <-as.data.frame.matrix(power_loss)

write_xlsx(powerloss_unif, 'D:/OSU/Research_Fall2023/power loss/powerloss_unif.xlsx')
save.image(paste0("New_powerloss_unif",".RData"))
