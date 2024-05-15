setwd("D:/OSU/Research_Fall2023/type I error/t test strategy II")
sample_size <- c(10,20,30,40,50)
shapiro_alpha <- c(0.1,0.05,0.01, 0.005,0)
alpha <- 0.05
N <- 10000
set.seed(1)

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
    sample_passed <- 0
    while(sample_passed < N){
      
      x1 <- rnorm(n, mean = 0, sd = 1)
      x2 <- rnorm(n, mean = 0, sd = 1)
      
      if(shapiro.test(c(x1-mean(x1), x2-mean(x2)))$p.value > pre_alpha ){
        sample_passed <- sample_passed +1
        if(t.test(x1,x2)$p.value < alpha){
          t_reject_h0 <- t_reject_h0 + 1 
          
        }
        
      }
      
    }
    error[i,j]<- round(t_reject_h0/N, 3)
  }
}
error
png(file="uniform_plot.jpeg", width=600, height=600)
myplot<-plot(sample_size,error[,1], type = "l", pch=24,lwd=3,ylim = c(0.00,0.25), 
             ylab = "Type I Error", xlab = "Sample Size")
title(main = "Type I Error Rate--strategy II -Normal...t.test")
points(sample_size, error[,2],type="b", pch=0,  lty=2, lwd=3,col="red")
points(sample_size, error[,3],  type="b", pch=23,  lty=3, lwd=3,col="blue")
points(sample_size, error[,4],  type="b", pch=21,  lty=4, lwd=3,col="green")
points(sample_size, error[,5],  type="b", pch=25,  lty=5, lwd=3,col="orange")
legend("topleft", legend=c("0.100", "0.050", '0.010',"0.005","w/o s-w test"),
       col=c("black","red","blue", "green", "orange"), lty =c(1:5),cex = 0.5)
dev.off()
Error_norm_II<- as.data.frame(error)

#export data into excel
library("writexl")
write_xlsx(Error_norm_II,"D:/OSU/Research_Fall2023/type I error/t test strategy II/Error.Normal_II.xlsx")

save.image(paste0("Normal.II_test",".RData"))

