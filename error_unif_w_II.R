setwd("D:/OSU/Research/Pval/strategy II")
set.seed(1)
alpha<-0.05
shapiro_alpha <- c(0.1,0.05,0.01, 0.005,0)
sample_size<-c(10,20,30,40,50)
pval<- numeric(length(shapiro_alpha)*length(sample_size))
pvalue<-array(pval, dim = c(length(sample_size), length((shapiro_alpha))),
              dimnames = list(sample_size,c(0.100, 0.050, 0.010, 0.005, "w/o pretest")))

num_sim<-0
for(i in 1: length(sample_size)){
  n<-sample_size[i]
  print(n)
  for(j in 1: length(shapiro_alpha)){
    s.alpha <- shapiro_alpha[j]
    print(s.alpha)
    sample_passed<-0 
    reject_h0 <-0
    while (sample_passed<10000){
      x1 <- runif(n, 0, 1)
      x2 <- runif(n, 0, 1)
      if(shapiro.test(c(x1 - mean(x1),x2 - mean(x2)))$p.value <= alpha){
        sample_passed <- sample_passed + 1
        
        if(wilcox.test(x1,x2)$p.value < alpha){
          reject_h0 <- reject_h0 + 1
        }
        
        num_sim <- num_sim + 1
      }
      
    }
    pvalue[i,j]<-round(reject_h0/sample_passed, 3)
    
  }
  
}

pvalue
myplot<-plot(sample_size,pvalue[,1], type = "l", pch=24,lwd=3,ylim = c(0.00,0.25), 
             ylab = "Type I Error", xlab = "Sample Size")
title(main = "Type I Error Rate -Normal_w.test")
points(sample_size, pvalue[,2],type="b", pch=0,  lty=2, lwd=3,col="red")
points(sample_size, pvalue[,3],  type="b", pch=23,  lty=3, lwd=3,col="blue")
points(sample_size, pvalue[,4],  type="b", pch=21,  lty=4, lwd=3,col="green")
points(sample_size, pvalue[,5],  type="b", pch=25,  lty=5, lwd=3,col="orange")
legend("topleft", legend=c("0.100", "0.050", '0.010',"0.005","w/o s-w test"),
       col=c("black","red","blue", "green", "orange"), lty =c(1:5),cex = 0.5)

error.unif_w_II <- as.data.frame(pvalue)

#export data into excel
library("writexl")
write_xlsx(error.unif_w_II, "D:/OSU/Research/Pval/strategy II/Error.unif.w_II.test.xlsx")

save.image(paste0("error.unif.w_II","myplot",.RData"))