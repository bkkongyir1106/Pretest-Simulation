# Load necessary libraries
#rm(list = ls())
if(!require("pacman")) install.packages("pacman")
pacman::p_load(dgof, tseries, nortest, dplyr, moments, LaplacesDemon, VGAM, extraDistr, statmod, evd)

# Generate data from different distributions but located similarly
generate_samples <- function(n, dist){
  if(dist == "normal"){
    samples <- rnorm(n, mean = 0, sd = 1)
  }
  if(dist == "normal_100"){
    samples <- rnorm(n, mean = 100, sd = 25)
  }
  if(dist == "normal_50"){
    samples <- rnorm(n, mean = 50, sd = 15)
  }
  if(dist == "normal_25"){
    samples <- rnorm(n, mean = 25, sd = 10)
  }
  if(dist == "normal_15"){
    samples <- rnorm(n, mean = 15, sd = 8)
  }
  if(dist == "normal_5"){
    samples <- rnorm(n, mean = 5, sd = 2)
  }
  if(dist == "normal_a"){
    samples <- rnorm(n, mean = 2, sd = 1)
  }
  if(dist == "Chi_Square"){
    samples <- rchisq(n, df = 3)
  }
  if(dist == "Gamma"){
    samples <- rgamma(n, shape = 2, rate = 1)
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
  if(dist == "Uniform_100"){
    samples <- runif(n, min = 20, max = 100)
  }
  if(dist == "Uniform_10"){
    samples <- runif(n, min = 10, max = 50)
  }
  if(dist == "Laplace"){
    samples <- LaplacesDemon::rlaplace(n, 0, 4)
  }
  if(dist == "Weibull"){
    samples <- rweibull(n, shape = 2, scale = 1)
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
    samples <- VGAM::rpareto(n = n, scale = 1, shape = 2)
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
  if(dist == "binomial"){
    samples <- rbinom(n, size = 10, prob = 0.5)
  }
  if(dist == "poisson"){
    samples <- rpois(n, lambda = 3)
  }
  if(dist == "geometric"){
    samples <- rgeom(n, prob = 0.3)
  }
  if(dist == "negative_binomial"){
    samples <- rnbinom(n, size = 5, prob = 0.5)
  }
  if(dist == "inverse_gamma"){
    samples <- LaplacesDemon::rinvgamma(n, shape = 3, scale = 2)
  }
  if(dist == "inverse_gaussian"){
    samples <- statmod::rinvgauss(n, mean = 1, shape = 1)
  }
  if(dist == "bimodal"){
    samples <- c(rnorm(n/2, mean = -2, sd = 1), rnorm(n/2, mean = 2, sd = 1))
  }
  if(dist == "triangular"){
    samples <- extraDistr::rtriang(n, a = 0, b = 1, c = 0.5)
  }
  if(dist == "logistic"){
    samples <- rlogis(n, location = 0, scale = 1)
  }
  if(dist == "Gumbel"){
    samples <- evd::rgumbel(n, loc = 0, scale = 1)
  }
  if(dist == "Rayleigh"){
    samples <- extraDistr::rrayleigh(n, sigma = 1)
  }
  return(samples)
}
