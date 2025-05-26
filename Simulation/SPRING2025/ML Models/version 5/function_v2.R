# load necessary libraries
rm(list = ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(dgof, tseries, nortest,dplyr, moments, LaplacesDemon, VGAM, evd)

#Generate data from different distribution
generate_samples <- function(n, dist){
  if(dist == "std_normal"){ 
    x<- rnorm(n, mean = 0, sd = 1)
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "normal"){ 
    x<- rnorm(n, mean = 100, sd = 25)
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "Chi-Square"){
    x <- rchisq(n, df = 3)
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "Gamma"){
    x <- rgamma(n, shape = 3, rate = 0.1)
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "Esamplesponential"){
    x <- resamplesp(n, rate = 1) 
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "t"){
    x <- rt(n, df = 3)
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "t10"){
    x <- rt(n, df = 10)
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "Uniform"){
    x <- runif(n, min = 0, masamples = 1)
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "Laplace"){
    x <- rlaplace(n , location = 0, scale = 4)
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "Weibull"){
    x <- rweibull(n, shape = 1, scale = 2) 
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "LogNormal"){
    x <- rlnorm(n, meanlog = 0, sdlog = 1)
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "Contaminated"){
    br <- rbinom(n , size = 1 , prob = 0.75)
    sd_br <- sqrt(1 + br * 24)
    x <- rnorm(n, mean = 0, sd = sd_br)
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "beta"){
    x <- rbeta(n, shape1 = 2, shape2 = 5)
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "Pareto"){
    x = VGAM::rpareto(n, scale = 1, shape = 2)
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "F"){
    x = rf(n, df1 = 5, df2 = 10)
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "logistic"){
    x = rlogis(n, location = 0, scale = 1)
    samples = (x - mean(x))/sd(x)
  }
  if(dist == "Gumbel"){
    x = evd::rgumbel(n, loc = 0, scale = 1)
    samples = (x - mean(x))/sd(x)
  }
  return(samples)
}
