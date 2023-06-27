setwd("D:/OSU/Research/Pval/power")
set.seed(1)
alpha<-0.05
shapiro_alpha <- c(0.1,0.05,0.01, 0.005,0)
sample_size<-c(10,20,30,40,50)
pval<- numeric(length(shapiro_alpha)*length(sample_size))
power<-array(pval, dim = c(length(sample_size), length((shapiro_alpha))),
             dimnames = list(sample_size,c(0.100, 0.050, 0.010, 0.005, "w/o pretest")))

num_sim<-0
for(i in 1: length(sample_size)){
  n<-sample_size[i]
  print(n)
  for(j in 1: length(shapiro_alpha)){
    s.alpha<-shapiro_alpha[j]
    print(s.alpha)
    reject_h0.t <- 0
    reject_h0.w <- 0
    passed <- 0
    fail <- 0
    for (k in 1:100000){
      x1 <-rnorm(n, 0.0, 1)
      x2 <- rnorm(n, 0.6, 1)
      if(shapiro.test(x1)$p.value > s.alpha & shapiro.test(x2)$p.value > s.alpha){
        passed <- passed + 1
        if(t.test(x1,x2)$p.value < alpha){
          reject_h0.t <- reject_h0.t + 1
        }
      }
        else if(wilcox.test(x1,x2)$p.value<= alpha){
            reject_h0.w <- reject_h0.w+1
        }
      
      error <-reject_h0.t+reject_h0.w
    }
    power[i,j]<-error/100000
  }
  
  
}
power
Power_norm<- as.data.frame(power)

#export data into excel
library("writexl")
write_xlsx(Power_norm,"D:/OSU/Research/Pval/power/Power_norm.xlsx")

save.image(paste0("Power_norm",".RData"))

