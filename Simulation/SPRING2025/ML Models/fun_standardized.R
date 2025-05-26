# load necessary libraries
rm(list = ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(dgof, tseries, nortest,dplyr, moments, LaplacesDemon, VGAM)

#Generate data from different distribution but located similarly
generate_samples <- function(n, dist){
  if(dist == "normal"){ 
    x<- rnorm(n, mean = 0, sd = 1)
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "normal_25"){ 
    x <- rnorm(n, mean = 100, sd = 25)
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "Chi-Square"){
    x <- rchisq(n, df = 3)
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "Gamma"){
    x <- rgamma(n, shape = 1, rate = 1)
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "Exponential"){
    x <- rexp(n, rate = 1) 
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "t"){
    x <- rt(n, df = 3)
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "t_15"){
    x <- rt(n, df = 15)
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "t_25"){
    x <- rt(n, df = 25)
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "Uniform"){
    x <- runif(n, min = 0, max = 1)
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "Laplace"){
    x <- rlaplace(n , location = 0, scale = 4)
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "Weibull"){
    x <- rweibull(n, shape = 1, scale = 1) 
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "LogNormal"){
    x <- rlnorm(n, meanlog = 0, sdlog = 1)
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  
  if(dist == "Contaminated"){
    br <- rbinom(n , size = 1 , prob = 0.75)
    sd_br <- sqrt(1 + br * 24)
    x <- rnorm(n, mean = 0, sd = sd_br)
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "Pareto"){
    shape = 3
    x <- rpareto(n, shape = shape)
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "beta"){
    x <- rbeta(n, shape1 = 2, shape2 = 5)
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "cauchy"){
    x <- rcauchy(n, location = 0, scale = 1)
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "spike"){
    spike_loc <- 0
    spike_prop <- 0.1
    x <- c(rep(spike_loc, n*spike_prop), rnorm(n*(1-spike_prop)))
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  if(dist == "extremeskew"){
    y <- rexp(n)
    x <- y^3 * sample(c(-1, 1), n, replace = TRUE)
    samples <- scale(x, center = TRUE, scale = TRUE)
  }
  return(samples)
}

