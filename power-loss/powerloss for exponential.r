setwd("D:/OSU/Research_Fall2023/power loss")
set.seed(1)
N <- 10000
alpha <- 0.05
sample_size <- c(10,20,30,40,50)
effect_size <-c(0.05,0.25, 0.5, 0.75, 1)

pval <-numeric(length(sample_size)*length(effect_size))
power_loss <- power_norm <-power_unif <-array(pval, dim = c(length(sample_size), 
                  length(effect_size)), dimnames = list(sample_size, effect_size))

for(i in 1 : length(sample_size)){
  n <- sample_size[i]
  print(n)
  for (j in 1 : length(effect_size)) {
    df <- effect_size[j]
    print(df)
    
    Norm_rejecth0 <- 0
    for(k in 1 : N) {
      n1 <- rnorm(n, mean = 0, sd = 1)
      n2 <- rnorm(n, mean = 0, sd = 1) + df
      if(t.test(n1,n2)$p.value <= alpha){
        Norm_rejecth0 <- Norm_rejecth0 + 1
      }
    }
    
    unif_rejectho <-0
    for(k in 1 : N) {
      e1 <- rexp(n, rate= 1)-1
      e2 <- rexp(n, rate= 1)-1 + df
      
      if(t.test(e1,e2)$p.value <= alpha){
        unif_rejectho <- unif_rejectho + 1
      }
    }
    
    power_norm[i,j] <- round(Norm_rejecth0/N, 3) 
    power_unif[i,j] <- round(unif_rejectho/N, 3)
    power_loss[i,j] <- power_norm[i,j]-power_unif[i,j]
  }
}

myplot <- plot(sample_size, power_loss[,1], type="l", col = "black", pch = 3, 
          ylim = c(-0.07, 0.07),xlab = "Sample Size", ylab = "Power loss")
lines(sample_size, power_loss[,2], type="b", col = "blue", pch = 4)
lines(sample_size, power_loss[,3], type="b", col = "red", pch = 5)
lines(sample_size, power_loss[,4], type="b", col = "pink", pch = 6)
lines(sample_size, power_loss[,5], type="b", col = "green", pch = 7)
title(main = "Power loss under t test for exponential dist.")
legend("topleft", legend=effect_size,col=c("black", "blue","red", "pink", "green"), 
       lty =c(1:4),cex = 0.5)

#install.packages("writexl")
library("writexl")

powerloss_exp <-as.data.frame.matrix(power_loss)

write_xlsx(powerloss_exp, 'D:/OSU/Research_Fall2023/power loss/powerloss_exp.xlsx')

