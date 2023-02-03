set.seed(110)
num.sim <- 50
alpha <- 0.05
mydist <- c("Expotential","Uniform", "Standard Normal")
mysum <- c(10,20,30,40,50)

mypval<-numeric(length(mydist)*length(mysum))
mypvalues <- array(mypval, dim = c(5, 3), dimnames = list(mysum,mydist))

for(i in 1:length(mysum)){
  n<-mysum[i]
  pval <- NULL
  for (k in 1:length(mydist)){ 
    dist<-mydist[k]
    for (j in 1: num.sim) {
      
      if(k==1){x<-rexp(n,1)} else if (k==2){ x<-runif(n,0,1)} else {x<-rnorm(n,0,1)}
      st<-shapiro.test(x)
      p[j]<-st$p.value
    }
    pval[k]<-mean(p>alpha)
  }
  mypvalues[i,] <- t(as.matrix(pval))
}

mypvalues
