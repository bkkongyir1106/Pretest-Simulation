setwd("D:/OSU/Research_Fall2023")
sample_size <- c(10,20,30,40,50)
shapiro_alpha <- c(0.1,0.05,0.01, 0.005,0)
alpha <- 0.05
N <- 10000

null_vec <-numeric(length(sample_size)*length(shapiro_alpha))
error <- array(data = null_vec, dim = c(length(sample_size),
    length(shapiro_alpha)), dimnames = list(sample_size, shapiro_alpha))

set.seed(1)
for( i in 1 : length(sample_size)){
  n <-sample_size[i]
  print(n)
  
  for( j in 1 : length(shapiro_alpha)){
    
    pre_alpha <- shapiro_alpha[j]
    print(pre_alpha)
    
    t_reject_h0 <-0
    sample_passed <- 0
    while(sample_passed < N){
      
      x1 <- runif(n, min = 0, max = 1)
      x2 <- runif(n, min = 0, max = 1)
      
      if(shapiro.test(x1)$p.value > pre_alpha & shapiro.test(x2)$p.value > pre_alpha){
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
