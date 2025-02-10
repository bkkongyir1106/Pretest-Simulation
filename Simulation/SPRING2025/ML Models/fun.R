# load necessary libraries
rm(list = ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(dgof, tseries, nortest,dplyr, moments, LaplacesDemon, VGAM)

#Generate data from different distribution but located similarly
generate_samples <- function(n, dist){
  if(dist == "normal"){ 
     x<- rnorm(n, mean = 100, sd = 25)
     #samples <- (x-mean(x))/sd(x)
     samples <- x
  }
  if(dist == "Chi-Square"){
    x <- rchisq(n, df = 3)
    #samples <- (x -mean(x))/sd(x)
    samples <- x
  }
  if(dist == "Gamma"){
    x <- rgamma(n, shape = 3, rate = 0.1)
    #samples <- (x -mean(x))/sd(x)
    samples <- x
  }
  if(dist == "Exponential"){
    x <- rexp(n, rate = 1) 
    #samples <- (x -mean(x))/sd(x)
    samples <- x
  }
  if(dist == "t"){
    x <- rt(n, df = 7)
    #samples <- (x -mean(x))/sd(x)
    samples <- x
  }
  if(dist == "Uniform"){
    x <- runif(n, min = 0, max = 1)
    #samples <- (x -mean(x))/sd(x)
    samples <- x
  }
  if(dist == "Laplace"){
    x <- rlaplace(n , location = 0, scale = 4)
    #samples <- (x -mean(x))/sd(x)
    samples <- x
  }
  if(dist == "Weibull"){
    x <- rweibull(n, shape = 1, scale = 2) 
    #samples <- (x -mean(x))/sd(x)
    samples <- x
  }
  if(dist == "LogNormal"){
    x <- rlnorm(n, meanlog = 0, sdlog = 1)
    #samples <- (x -mean(x))/sd(x)
    samples <- x
  }
  
  if(dist == "Contaminated"){
    br <- rbinom(n , size = 1 , prob = 0.75)
    sd_br <- sqrt(1 + br * 24)
    x <- rnorm(n, mean = 0, sd = sd_br)
    #samples <- (x -mean(x))/sd(x)
    samples <- x
  }
  if(dist == "Pareto"){
    shape = 3
    x <- rpareto(n, shape = shape)
    #samples <- (x -mean(x))/sd(x)
    samples <- x
  }
  return(samples)
}
