# load necessary libraries
library("nortest")
library("dgof")
library("dplyr")
library(moments)
library(tseries)
library(LaplacesDemon)
library(VGAM)
#Generate data from different distribution but located similarly
generate_data <- function(n, dist){
    if(dist == "Standard Normal"){ 
      x <- rnorm(n, mean = 0, sd = 1)
    }
    if(dist == "Chi-Square"){
      x <- (rchisq(n, df = 3) - 3)/sqrt(6) 
    }
    if(dist == "Gamma"){
      x <- (rgamma(n, shape = 3, rate = 0.1) - 30)/sqrt(300)
    }
    if(dist == "Exponential"){
      x <- rexp(n, rate = 1) - 1 
    }
    if(dist == "t"){
      x <- (rt(n, df = 7))/sqrt(7/5) 
    }
    if(dist == "Uniform"){
      x <- (runif(n, min = 0, max = 1) - 0.5)*sqrt(12) 
    }
    if(dist == "Laplace"){
      x <- rlaplace(n , location = 0, scale = 4)/sqrt(8)
    }
    if(dist == "Weibull"){
      x <- (rweibull(n, shape = 1, scale = 2) - 2*gamma(51/50))/sqrt(4*(gamma(3) - gamma(2))) 
    }
    if(dist == "LogNormal"){
      x <- (rlnorm(n, meanlog = 0, sdlog = 1) - exp(0 + 1/2))/sqrt((exp(1)-1)*exp(2*0 + 1)) 
    }
   
    if(dist == "Contaminated"){
      br <- rbinom(n , size = 1 , prob = 0.75)
      sd_br <- sqrt(1 + br * 24)
      x <- rnorm(n, mean = 0, sd = sd_br)/sqrt(0.25 + 0.75*24)
    }
  if(dist == "Pareto"){
    shape = 3
    mean_pareto <- shape / (shape - 1)
    sd_pareto <- shape/(((shape - 1)^2)*(shape - 2))
    x <- (rpareto(n, shape = shape) - mean_pareto)/sd_pareto
  }
  return(x)
}

#Apply different tests
generate_tests <- function(x, test){
  if(test == "KS"){
    output <- lillie.test(x)
  }
  if(test == "SW"){
    output <- shapiro.test(x)
  }
  if(test == "JB"){
    output <- jarque.test(x)
  }
  if(test == "DAP"){
    output <- agostino.test(x)
  }
  if(test == "AD"){
    output <- ad.test(x)
  }
  if(test == "SF"){
    output <- sf.test(x)
  }
  if(test == "CVM"){
    output <- cvm.test(x)
  }
  if(test == "CHISQ"){
    output <- chisq.test(x)
  }
  
return(output)
}
