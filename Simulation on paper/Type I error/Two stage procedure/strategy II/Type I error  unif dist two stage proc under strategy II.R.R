setwd("D:/OSU/Research/Pval/strategy II")
set.seed(1)
alpha<-0.05
N<-100000
shapiro_alpha <- c(0.1,0.05,0.01, 0.005,0.000)
sample_size<-c(10,20,30,40,50)
pval<- numeric(length(shapiro_alpha)*length(sample_size))
pvalue<-array(pval, dim = c(length(sample_size), length((shapiro_alpha))),
              dimnames = list(sample_size,shapiro_alpha))

for(i in 1: length(sample_size)){
  n<-sample_size[i]
  print(n)
  for(j in 1: length(shapiro_alpha)){
    s.alpha<-shapiro_alpha[j]
    print(s.alpha)
    reject_h0.t <- 0
    reject_h0.w <- 0
    for(k in 1:N){
      x1<-runif(n, 0, 1)
      x2<-runif(n, 0, 1)
      if(shapiro.test(c(x1 - mean(x1),x2 - mean(x2)))$p.value > s.alpha){
        if(t.test(x1,x2)$p.value < alpha){
          reject_h0.t <- reject_h0.t + 1
        }
      }
      else if(shapiro.test(c(x1 - mean(x1),x2 - mean(x2)))$p.value <= s.alpha){
        if(wilcox.test(x1,x2)$p.value < alpha){
          reject_h0.w <- reject_h0.w + 1
        }
      }
      
    }
    pvalue[i,j]<-round((reject_h0.t+reject_h0.w)/N, 3)
    
  }
  
  
  
}

pvalue
myplot<-plot(sample_size,pvalue[,1], type = "l", pch=24,lwd=3,ylim = c(0.00,0.25), 
             ylab = "Type I Error", xlab = "Sample Size")
title(main = "Type I Error Rate--strategy II -unif...test")
points(sample_size, pvalue[,2],type="b", pch=0,  lty=2, lwd=3,col="red")
points(sample_size, pvalue[,3],  type="b", pch=23,  lty=3, lwd=3,col="blue")
points(sample_size, pvalue[,4],  type="b", pch=21,  lty=4, lwd=3,col="green")
points(sample_size, pvalue[,5],  type="b", pch=25,  lty=5, lwd=3,col="orange")
legend("topleft", legend=c("0.100", "0.050", '0.010',"0.005","w/o s-w test"),
       col=c("black","red","blue", "green", "orange"), lty =c(1:5),cex = 0.5)

Error_unif_II<- as.data.frame(pvalue)

#export data into excel
library("writexl")
write_xlsx(Error_unif_II,"D:/OSU/Research/Pval/strategy II/Error.unif_II.xlsx")

save.image(paste0("unif_II.test",".RData"))


