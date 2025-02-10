# load necessary libraries
rm(list = ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(dgof, tseries, nortest,dplyr, moments, LaplacesDemon, VGAM, evd)

#Generate data from different distribution
generate_samples <- function(n, dist){
  if(dist == "normal"){ 
    samples<- rnorm(n, mean = 100, sd = 25)
  }
  if(dist == "Chi-Square"){
    samples <- rchisq(n, df = 3)
  }
  if(dist == "Gamma"){
    samples <- rgamma(n, shape = 3, rate = 0.1)
  }
  if(dist == "Exponential"){
    samples <- rexp(n, rate = 1) 
  }
  if(dist == "t"){
    samples <- rt(n, df = 7)
  }
  if(dist == "Uniform"){
    samples <- runif(n, min = 0, max = 1)
  }
  if(dist == "Laplace"){
    samples <- rlaplace(n , location = 0, scale = 4)
  }
  if(dist == "Weibull"){
    samples <- rweibull(n, shape = 1, scale = 2) 
  }
  if(dist == "LogNormal"){
    samples <- rlnorm(n, meanlog = 0, sdlog = 1)
  }
  if(dist == "Contaminated"){
    br <- rbinom(n , size = 1 , prob = 0.75)
    sd_br <- sqrt(1 + br * 24)
    samples <- rnorm(n, mean = 0, sd = sd_br)
  }
  if(dist == "beta"){
    samples <- rbeta(n, shape1 = 2, shape2 = 5)
  }
  if(dist == "Pareto"){
    samples = VGAM::rpareto(n, scale = 1, shape = 2)
  }
  if(dist == "F"){
    samples = rf(n, df1 = 5, df2 = 10)
  }
  if(dist == "logistic"){
    samples = rlogis(n, location = 0, scale = 1)
  }
  if(dist == "Gumbel"){
    samples = evd::rgumbel(n, loc = 0, scale = 1)
  }
  return(samples)
}

generate_samples(10, "Gumbel")
