setwd("D:/OSU/Research_Fall2023/power loss")
set.seed(1)
N <- 10000
alpha <- 0.05
sample_size <- c(10,20,30,40,50)
#effect_size <-c(0.05,0.25, 0.5, 0.75, 1)
effect_size <-seq(from = -1, to = 1, by = 0.15)
pval <-numeric(length(sample_size)*length(effect_size))
power_loss <- power_norm <-power_unif <-array(pval, dim = c(length(sample_size), 
          length(effect_size)), dimnames = list(sample_size, effect_size))

for(i in 1 : length(sample_size)){
  n <- sample_size[i]
  print(n)
  for (j in 1 : length(effect_size)) {
    df <- effect_size[j]
    print(df)
    Norm_rejecth0_t <- 0
    Norm_rejecth0_w <- 0
    for(k in 1 : N) {
      n1 <- rnorm(n, mean = 0, sd = 1)
      n2 <- rnorm(n, mean = 0, sd = 1) + df
      if(shapiro.test(n1)$p.value > alpha & shapiro.test(n2)$p.value > alpha){
        if(t.test(n1,n2)$p.value <= alpha){
          Norm_rejecth0_t <- Norm_rejecth0_t + 1
        }
      }
      else if(shapiro.test(n1)$p.value <= alpha | shapiro.test(n2)$p.value <= alpha){
        if(wilcox.test(n1,n2)$p.value <= alpha){
          Norm_rejecth0_w <- Norm_rejecth0_w + 1
        }
      }
    }
    
    unif_rejecth0_t <- 0
    unif_rejecth0_w <- 0
    for(k in 1 : N) {
      u01 <- runif(n, min = 0, max = 1)-0.5
      u02 <- runif(n, min = 0, max = 1) -0.5 + df
      u1 <-u01/sqrt(1/12)
      u2 <- u02/sqrt(1/12)
      if(shapiro.test(u1)$p.value > alpha & shapiro.test(u2)$p.value > alpha){
        if(t.test(u1,u2)$p.value <= alpha){
          unif_rejecth0_t <- unif_rejecth0_t + 1
        }
      }
      else if(shapiro.test(u1)$p.value <= alpha | shapiro.test(u2)$p.value <= alpha){
        if(wilcox.test(u1,u2)$p.value <= alpha){
          unif_rejecth0_w <- unif_rejecth0_w + 1
        }
      }
    }
    
    power_norm[i,j] <- round((Norm_rejecth0_t + Norm_rejecth0_w)/N, 3) 
    power_unif[i,j] <- round((unif_rejecth0_t + unif_rejecth0_w)/N, 3)
    power_loss[i,j] <- power_norm[i,j]-power_unif[i,j]
  }
}
power_loss

png(file="unif_two_stage_powerloss.jpeg", width=1000, height=600)
plot(effect_size, power_loss[1,], type="l", col = "black", 
     ylim = c(-1, 0.5),xlab = "effect Size", ylab = "Power loss")
lines(effect_size, power_loss[2,], col = "blue")
lines(effect_size, power_loss[3,],  col = "red")
lines(effect_size, power_loss[4,],  col = "pink")
lines(effect_size, power_loss[5,],  col = "green")
title(main = "Power loss under t test for exponential dist.")
legend("topright", legend=sample_size,col=c("black", "blue","red", "pink", "green"), 
       lty =c(1:4),cex = 0.5)
dev.off()
#install.packages("writexl")
library("writexl")

powerloss_unif_two_stage <-as.data.frame.matrix(power_loss)

write_xlsx(powerloss_unif_two_stage, 'D:/OSU/Research_Fall2023/power loss/powerloss_unif_two_stage.xlsx')
save.image(paste0("powerloss_unif_two stage",".RData"))
