setwd("D:/OSU/Research_Fall2023/two stage strategy I")
set.seed(1)
sample_size <- c(10,20,30,40,50)
shapiro_alpha <- c(0.1,0.05,0.01, 0.005)
alpha <- 0.05
N <- 100000

null_vec <-numeric(length(sample_size)*length(shapiro_alpha))
error <- array(data = null_vec, dim = c(length(sample_size),
    length(shapiro_alpha)), dimnames = list(sample_size, shapiro_alpha))

for( i in 1 : length(sample_size)){
  n <-sample_size[i]
  print(n)
  
  for( j in 1 : length(shapiro_alpha)){
    
    pre_alpha <- shapiro_alpha[j]
    print(pre_alpha)
    
    t_reject_h0 <-0
    w_reject_h0 <-0
    sample_passed <- 0
    sample_fail <- 0
    for(k in 1 : N){
      
      x1 <- rexp(n, rate = 1)
      x2 <- rexp(n, rate = 1)
      
      if(shapiro.test(x1)$p.value > pre_alpha & shapiro.test(x2)$p.value > pre_alpha){
        sample_passed <- sample_passed +1
        if(t.test(x1,x2)$p.value < alpha){
          t_reject_h0 <- t_reject_h0 + 1 
          
        }
      }
      else if((shapiro.test(x1)$p.value <=  pre_alpha | shapiro.test(x2)$p.value <= pre_alpha)){
        sample_fail <- sample_fail + 1
        if(wilcox.test(x1,x2)$p.value <= alpha){
          w_reject_h0 <- w_reject_h0 + 1
        }
      }
    }
    
    error[i,j]<- round((t_reject_h0 + w_reject_h0)/N, 3)
  }
}
error
myplot<-plot(sample_size,error[,1], type = "l", pch=24,lwd=3,ylim = c(0.00,0.25), 
             ylab = "Type I Error", xlab = "Sample Size")
title(main = "Type I Error Rate--Two stage strategy I exponential")
points(sample_size, error[,2],type="b", pch=0,  lty=2, lwd=3,col="red")
points(sample_size, error[,3],  type="b", pch=23,  lty=3, lwd=3,col="blue")
points(sample_size, error[,4],  type="b", pch=21,  lty=4, lwd=3,col="green")
points(sample_size, error[,5],  type="b", pch=25,  lty=5, lwd=3,col="orange")
legend("topleft", legend=c("0.100", "0.050", '0.010',"0.005","w/o s-w test"),
       col=c("black","red","blue", "green", "orange"), lty =c(1:5),cex = 0.5)

Error_two_stage_exp_I<- as.data.frame(error)

#export data into excel
library("writexl")
write_xlsx(Error_two_stage_exp_I,"D:/OSU/Research_Fall2023/two stage strategy I/Error_two_stage_exp_I.xlsx")

save.image(paste0("Error_two_stage_exp_I",".RData"))

