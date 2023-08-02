setwd("D:/OSU/Research/Pval/power")
set.seed(1)
alpha<-0.05
sample_size<-c(10,20,30,40,50)
power_w.test<- numeric(length(sample_size))


for( i in 1:length(sample_size)){
  n <- sample_size[i]
  print(i)
  reject_h0<-0
  for (k in 1:100000){
    x1 <-rnorm(n, 0.0, 1.0)
    x2 <- rnorm(n, 0.6, 1.0)
    if(wilcox.test(x1,x2)$p.value<alpha){
      reject_h0 <- reject_h0+1
    }
    
  }
  power_w.test[i]<-reject_h0/100000
  
}
power_w.test

norm_power_w.test_only<- as.data.frame(power_w.test)

#export data into excel
library("writexl")
write_xlsx(norm_power_w.test_only,"D:/OSU/Research/Pval/power/norm_power_w.test_only.xlsx")

save.image(paste0("norm_power_w.test_only",".RData"))

