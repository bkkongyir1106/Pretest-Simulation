# Load necessary libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  e1071, tseries, dgof, nortest, ggplot2, dplyr,
  tidyverse, LaplacesDemon, VGAM, patchwork, reshape2
)

# ----------------------------------------------------
#   Generate standardized data from various distributions
# ----------------------------------------------------
generate_data <- function(n, dist, par = NULL) {
  # Input validation
  if (!is.numeric(n) || n <= 0) stop("n must be a positive integer")
  dist <- tolower(dist)
  
  if (dist == "normal") {
    if (is.null(par)) par <- c(0, 1)
    x <- rnorm(n, mean = par[1], sd = par[2])
    
  } else if (dist == "chi_square") {
    if (is.null(par)) par <- 3
    x <- (rchisq(n, df = par) - par) / sqrt(2 * par)
    
  }else if (dist == "chi_sq7") {
    if (is.null(par)) par <- 7
    x <- (rchisq(n, df = par) - par) / sqrt(2 * par)
    
  }else if (dist == "gamma") {
    if (is.null(par)) par <- c(3, 0.1)
    mean_g <- par[1] / par[2]
    sd_g <- sqrt(par[1] / par[2]^2)
    x <- (rgamma(n, shape = par[1], rate = par[2]) - mean_g) / sd_g
    
  } else if (dist == "exponential") {
    if (is.null(par)) par <- 1
    mean_e <- 1 / par
    sd_e <- 1 / par
    x <- (rexp(n, rate = par) - mean_e) / sd_e
    
  } else if (dist == "t") {
    if (is.null(par)) par <- 3
    x <- rt(n, df = par) / sqrt(par / (par - 2))
    
  } else if (dist == "t_5") {
  if (is.null(par)) par <- 5
  x <- rt(n, df = par) / sqrt(par / (par - 2))
  
}else if (dist == "uniform") {
    if (is.null(par)) par <- c(0, 1)
    mean_u <- (par[1] + par[2]) / 2
    sd_u <- sqrt((par[2] - par[1])^2 / 12)
    x <- (runif(n, min = par[1], max = par[2]) - mean_u) / sd_u
    
  } else if (dist == "laplace") {
    if (is.null(par)) par <- c(2, 7)  
    x <- LaplacesDemon::rlaplace(n, location = par[1], scale = par[2])
    x <- (x - par[1]) / sqrt(2 * par[2]^2)
    
  } else if (dist == "weibull") {
    if (is.null(par)) par <- c(1, 2)
    mean_w <- par[2] * gamma(1 + 1 / par[1])
    var_w <- par[2]^2 * (gamma(1 + 2 / par[1]) - gamma(1 + 1 / par[1])^2)
    x <- (rweibull(n, shape = par[1], scale = par[2]) - mean_w) / sqrt(var_w)
    
  } else if (dist == "lognormal") {
    if (is.null(par)) par <- c(0, 1)
    mean_ln <- exp(par[1] + par[2]^2 / 2)
    var_ln <- (exp(par[2]^2) - 1) * exp(2 * par[1] + par[2]^2)
    x <- (rlnorm(n, meanlog = par[1], sdlog = par[2]) - mean_ln) / sqrt(var_ln)
    
  } else if (dist == "contaminated") {
    # par = c(prob_good, mean, sd_good, sd_bad)
    if (is.null(par)) par <- c(0.75, 0, 1, 5)
    br <- rbinom(n, size = 1, prob = par[1])
    sd_br <- ifelse(br == 1, par[3], par[4])
    x_raw <- rnorm(n, mean = par[2], sd = sd_br)
    # Variance: weighted avg of variances (mean is constant)
    var_c <- par[1] * par[3]^2 + (1 - par[1]) * par[4]^2
    x <- (x_raw - par[2]) / sqrt(var_c)
    
  } else if (dist == "pareto") {
    if (is.null(par)) par <- c(1, 3)
    if (par[2] <= 2) stop("Shape parameter for Pareto must be > 2 for variance to exist")
    mean_p <- par[1] * par[2] / (par[2] - 1)
    var_p <- (par[1]^2 * par[2]) / ((par[2] - 1)^2 * (par[2] - 2))
    x <- (VGAM::rpareto(n, scale = par[1], shape = par[2]) - mean_p) / sqrt(var_p)
    
  }else if (dist == "beta") {
  # par = c(alpha, beta), both > 0
  if (is.null(par)) par <- c(2, 5)
  if (any(par <= 0)) stop("Beta parameters must be > 0")
  mean_b <- par[1] / (par[1] + par[2])
  var_b  <- (par[1] * par[2]) / ((par[1] + par[2])^2 * (par[1] + par[2] + 1))
  x <- (rbeta(n, shape1 = par[1], shape2 = par[2]) - mean_b) / sqrt(var_b)
  
}else if (dist == "logistic") {
  # par = c(location, scale>0)
  if (is.null(par)) par <- c(0, 1)
  if (par[2] <= 0) stop("Logistic scale must be > 0")
  mean_l <- par[1]
  var_l  <- (pi^2 / 3) * par[2]^2
  x <- (rlogis(n, location = par[1], scale = par[2]) - mean_l) / sqrt(var_l)
  
} else {
    stop("Unsupported distribution: ", dist)
  }
  
  return(x)
}


# ---------------------------------------
# Define the functions to generate data 
# ---------------------------------------
Generate_data <- function(datagen.type, n = NULL, dist = NULL, two_samples = FALSE, 
                          priors = NULL, x_weights = NULL, y_weights = NULL, ...) {
  
  # generation via function
  if (datagen.type == 1) {
    if (!is.null(priors) && length(priors) == 2 && !is.null(x_weights) && dist == "mixture") {
      x_weights <- x_weights / sum(x_weights)
      # Generate data for x 
      x <- vector(length = n)
      for (i in 1:n) {
        component <- sample(1:2, size = 1, prob = x_weights)
        x[i] <- generate_data(1, priors[component])
      }
      if (two_samples) {
        y_weights <- y_weights / sum(y_weights)
        # Generate data for y
        y <- vector(length = n)
        for (i in 1:n) {
          component <- sample(1:2, size = 1, prob = y_weights)
          y[i] <- generate_data(1, priors[component])
        }
      } else {
        y <- NULL  
      }
    } else {
      x <- generate_data(n, dist)  
      if (two_samples) {
        y <- generate_data(n, dist)  
      } else {
        y <- NULL  
      }
    }
  }
  # loading from a CSV file
  if (datagen.type == 2) {
    # Prompt user to select a CSV file
    file_path <- file.choose()  
    data <- read.csv(file_path, header = TRUE)  
    if (two_samples) {
      if (ncol(data) < 2) {
        stop("The CSV file should contain at least two columns for two samples (x and y).")
      }
      x <- data[[1]]  
      y <- data[[2]]  
    } else {
      x <- data[[1]]  
      y <- NULL       
    }
  }
  
  return(list(x = x, y = y))
}

# ---------------------------------------------------------------------------
# Define function to perform different Normality tests
# ---------------------------------------------------------------------------
generate_tests <- function(x, test){
  if(test == "KS"){  # Kolmogorov-Smirnov test (Lilliefors)
    output <- nortest::lillie.test(x)
    
  } else if(test == "SW"){  # Shapiro-Wilk
    output <- shapiro.test(x)
    
  } else if(test == "JB"){  # Jarque-Bera
    output <- tseries::jarque.bera.test(x)
    
  } else if(test == "DAP"){  # D'Agostino test (skewness)
    output <- moments::agostino.test(x)
    
  } else if(test == "ANS"){  # Anscombe-Glynn test (kurtosis)
    output <- moments::anscombe.test(x)
    
  } else if(test == "AD"){  # Anderson-Darling
    output <- nortest::ad.test(x)
  
  } else if(test == "SF"){  # Shapiro-Francia
    output <- nortest::sf.test(x)
    
  } else if(test == "CVM"){  # Cramer-von Mises
    output <- nortest::cvm.test(x)
    
  } else if(test == "CHISQ"){  # Pearson Chi-Square
    # Need to bin data since chisq.test is for frequencies
    breaks <- pretty(x, n = sqrt(length(x)))  # Sturges rule
    observed <- table(cut(x, breaks = breaks))
    expected <- dnorm(breaks[-length(breaks)], mean(x), sd(x))
    expected <- expected / sum(expected) * sum(observed)
    output <- chisq.test(x = observed, p = expected, rescale.p = TRUE, simulate.p.value = TRUE)
    
  } else if(test == "SKEW"){  # Skewness test (z-test)
    s <- moments::skewness(x)
    se <- sqrt(6/length(x))
    z <- s / se
    p <- 2 * (1 - pnorm(abs(z)))
    output <- list(statistic = z, p.value = p, method = "Skewness z-test")
    
  } else if(test == "KURT"){  # Kurtosis test (z-test)
    k <- moments::kurtosis(x)
    se <- sqrt(24/length(x))
    z <- (k - 3) / se
    p <- 2 * (1 - pnorm(abs(z)))
    output <- list(statistic = z, p.value = p, method = "Kurtosis z-test")
    
  } else {
    stop("Unknown test name. Choose from: KS, SW, JB, DAP, ANS, AD, AD2, SF, CVM, CHISQ, SKEW, KURT.")
  }
  
  return(output)
}


# ------------------------------------
# Define simple calculation functions 
# -------------------------------------
# calculate one-sample  test statistic 
OneSample_test_statistic <- function(x) {
  return((mean(x)* sqrt(length(x))) / sd(x))
}
# calculate two-sample  test statistic 
TwoSample_test_statistic <- function(x, y) {
  return((mean(x) - mean(y))/sqrt(var(x)/length(x) + var(y)/length(y)))
}

# compute the area under the curve 
compute_area <- function(x, y) {
  (sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2)) / (max(sample_size) - min(sample_size))
}

# ---------------------------------------------------------
# Define functions for performing various downstream tests 
# ---------------------------------------------------------

# -----------One Sample Permutation test ------
one_sample_permutation_test <- function(x, B){
  observe_stat <- OneSample_test_statistic(x)
  permuted_stat <- numeric(B)
  for (j in 1:B) {
    index <- sample(c(-1, 1), length(x), replace = TRUE)
    sample_data <- index * abs(x)
    permuted_stat[j] <- OneSample_test_statistic(sample_data)
  }
  
  return(mean(abs(permuted_stat) >= abs(observe_stat)))
}

# --------- Two Sample Permutation -----------
two_sample_permutation_test <- function(x, y, B) {
  observed <- TwoSample_test_statistic(x, y)
  combined <- c(x, y)
  perm_stats <- replicate(B, {
    permuted <- sample(combined)
    x_star <- permuted[1:length(x)]
    y_star <- permuted[(length(x)+1):(2*length(x))]
    TwoSample_test_statistic(x_star, y_star)
  })
  
  return(mean(abs(perm_stats) >= abs(observed)))
}

