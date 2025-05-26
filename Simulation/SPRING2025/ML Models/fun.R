# load necessary libraries
rm(list = ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(dgof, tseries, nortest,dplyr, moments, LaplacesDemon, VGAM)

#Generate data from different distribution but located similarly
generate_samples <- function(n, dist){
  if(dist == "normal"){ 
    samples<- rnorm(n, mean = 0, sd = 1)
  }
  if(dist == "normal_100"){ 
    samples <- rnorm(n, mean = 100, sd = 25)
  }
  if(dist == "normal_50"){ 
    samples <- rnorm(n, mean = 50, sd = 15)
  }
if(dist == "normal_15"){ 
  samples <- rnorm(n, mean = 15, sd = 8)
}
  if(dist == "normal_5"){ 
    samples <- rnorm(n, mean = 5, sd = 2)
  }
  if(dist == "Chi-Square"){
    samples <- rchisq(n, df = 3)
  }
  if(dist == "Gamma"){
    samples <- rgamma(n, shape = 1, rate = 1)
  }
  if(dist == "Exponential"){
    samples <- rexp(n, rate = 1) 
  }
  if(dist == "t"){
    samples <- rt(n, df = 3)
  }
  if(dist == "t_7"){
    samples <- rt(n, df = 7)
  }
  if(dist == "t_15"){
    samples <- rt(n, df = 15)
  }
  if(dist == "Uniform"){
    samples <- runif(n, min = 0, max = 1)
  }
  if(dist == "Uniform_12.5"){
    samples <- runif(n, min = 25, max = 50)
  }
  if(dist == "Laplace"){
    samples <- rlaplace(n , location = 0, scale = 4)
  }
  if(dist == "Weibull"){
    samples <- rweibull(n, shape = 1, scale = 1) 
  }
  if(dist == "LogNormal"){
    samples <- rlnorm(n, meanlog = 0, sdlog = 1)
  }
  if(dist == "Contaminated"){
    br <- rbinom(n , size = 1 , prob = 0.75)
    sd_br <- sqrt(1 + br * 24)
    samples <- rnorm(n, mean = 0, sd = sd_br)
  }
  if(dist == "Pareto"){
    shape = 3
    samples <- rpareto(n, shape = shape)
  }
  if(dist == "beta"){
    samples <- rbeta(n, shape1 = 2, shape2 = 5)
  }
  if(dist == "cauchy"){
    samples <- rcauchy(n, location = 0, scale = 1)
  }
  if(dist == "spike"){
    spike_loc <- 0
    spike_prop <- 0.1
    samples <- c(rep(spike_loc, n*spike_prop), rnorm(n*(1-spike_prop)))
  }
  if(dist == "extremeskew"){
    y <- rexp(n)
    samples <- y^3 * sample(c(-1, 1), n, replace = TRUE)
  }
  return(samples)
}

