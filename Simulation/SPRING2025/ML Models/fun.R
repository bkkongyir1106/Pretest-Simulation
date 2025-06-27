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
  if(dist == "normal_25"){
    samples <- rnorm(n, mean = 50, sd = 15)
  }
  if(dist == "normal_15"){
    samples <- rnorm(n, mean = 15, sd = 8)
  }
  if(dist == "normal_5"){
    samples <- rnorm(n, mean = 5, sd = 2)
  }
  if(dist == "Chi_Square"){
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
  if(dist == "t_5"){
    samples <- rt(n, df = 5)
  }
  if(dist == "t_10"){
    samples <- rt(n, df = 10)
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
    samples <- VGAM::rpareto(n = 1000, scale = 1, shape = 2)
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

# 
# generate_samples <- function(n, dist) {
#   samples <- switch(dist,
#                     normal        = rnorm(n, mean = 0,   sd = 1),
#                     normal_100    = rnorm(n, mean = 100, sd = 25),
#                     normal_50     = rnorm(n, mean = 50,  sd = 15),
#                     normal_15     = rnorm(n, mean = 15,  sd = 8),
#                     normal_5      = rnorm(n, mean = 5,   sd = 2),
#                     Chi_Square    = rchisq(n, df = 3),
#                     Gamma         = rgamma(n, shape = 1, rate = 1),
#                     Exponential   = rexp(n, rate = 1),
#                     t             = rt(n, df = 3),
#                     t_7           = rt(n, df = 7),
#                     t_15          = rt(n, df = 15),
#                     Uniform       = runif(n, min = 0,  max = 1),
#                     Uniform_12_5  = runif(n, min = 25, max = 50),
#                     Laplace       = rlaplace(n, location = 0, scale = 4),
#                     Weibull       = rweibull(n, shape = 1, scale = 1),
#                     LogNormal     = rlnorm(n, meanlog = 0, sdlog = 1),
#                     Contaminated  = { br <- rbinom(n,1,0.75); rnorm(n,0,sqrt(1+br*24)) },
#                     Pareto        = rpareto(n,  shape = 3),
#                     beta          = rbeta(n, shape1 = 2, shape2 = 5),
#                     cauchy        = rcauchy(n, location = 0, scale = 1),
#                     spike         = {
#                       spike_prop <- 0.1
#                       k <- round(n * spike_prop)
#                       c(rep(0, k), rnorm(n - k))
#                     },
#                     extremeskew   = { y <- rexp(n); y^3 * sample(c(-1,1), n, replace = TRUE) },
#                     stop("Unknown distribution: ", dist)
#   )
#   return(samples)
# }
# 

