# Load necessary libraries
if (!require("pacman")) install.packages("pacman")
pacman::p_load(
  e1071, tseries, dgof, nortest, ggplot2, dplyr,
  tidyverse, LaplacesDemon, VGAM, patchwork, reshape2, extraDistr
)

# ----------------------------------------------------
#   Generate standardized data from various distributions
#   Centered around median for sign test purposes
# ----------------------------------------------------
generate_data <- function(n, dist, par = NULL) {
  # Input validation
  if (!is.numeric(n) || n <= 0) stop("n must be a positive integer")
  dist <- tolower(dist)
  
  if (dist == "normal") {
    if (is.null(par)) par <- c(0, 1)
    x <- rnorm(n, mean = par[1], sd = par[2])
    # Center around median (same as mean for normal)
    x <- x - median(x)
    
  } else if (dist == "chi_square") {
    if (is.null(par)) par <- 3
    x <- rchisq(n, df = par)
    # Center around median
    #median_chi <- qchisq(0.5, df = par)
    median_chi <- par *(1-2/(9*par))^3
    sd_chi <- sqrt(2*par)
    x <- (x - median_chi) 
    
  } else if (dist == "gamma") {
    if (is.null(par)) par <- c(3, 0.1)
    x <- rgamma(n, shape = par[1], rate = par[2])
    # Center around median
    median_g <- qgamma(0.5, shape = par[1], rate = par[2])
    sd_g <- sqrt(par[1] / par[2]^2)
    x <- (x - median_g)/sd_g 
    
  } else if (dist == "exponential") {
    if (is.null(par)) par <- 1
    x <- rexp(n, rate = par)
    # Center around median
    median_e <- qexp(0.5, rate = par)
    median_e <- log(2)/par
    sd_e <- 1 / par
    x <- (x - median_e)/sd_e
    
  } else if (dist == "t") {
    if (is.null(par)) par <- 3
    if (par <= 1) stop("Degrees of freedom must be > 1")
    x <- rt(n, df = par)
    # Center around median (0 for symmetric t-distribution)
    x <- x / sqrt(par / (par - 2))
    
  } else if (dist == "uniform") {
    if (is.null(par)) par <- c(0, 1)
    x <- runif(n, min = par[1], max = par[2])
    # Center around median
    median_u <- (par[1] + par[2]) / 2
    sd_u <- sqrt((par[2] - par[1])^2 / 12)
    x <- (x - median_u)/sd_u
    
  } else if (dist == "laplace") {
    if (is.null(par)) par <- c(3, 1)
    x <- extraDistr::rlaplace(n, mu = par[1], sigma = par[2])
    x <- x - par[1]
    sd_l <- sqrt(2 * par[2]^2)
    x <- x / sd_l
    
  } else if (dist == "weibull") {
    if (is.null(par)) par <- c(1, 2)
    x <- rweibull(n, shape = par[1], scale = par[2])
    # Center around median
    median_w <- par[2] * (log(2))^(1/par[1])
    var_w <- par[2]^2 * (gamma(1 + 2 / par[1]) - gamma(1 + 1 / par[1])^2)
    x <- (x - median_w) / sqrt(var_w)
    
  } else if (dist == "lognormal") {
    if (is.null(par)) par <- c(0, 1)
    x <- rlnorm(n, meanlog = par[1], sdlog = par[2])
    # Center around median
    median_ln <- exp(par[1])
    var_ln <- (exp(par[2]^2) - 1) * exp(2 * par[1] + par[2]^2)
    x <- (x - median_ln) 
    
  } else if (dist == "contaminated") {
    # par = c(prob_good, mean, sd_good, sd_bad)
    if (is.null(par)) par <- c(0.75, 0, 1, 5)
    br <- rbinom(n, size = 1, prob = par[1])
    sd_br <- ifelse(br == 1, par[3], par[4])
    x_raw <- rnorm(n, mean = par[2], sd = sd_br)
    # Center around median (approximate for contaminated normal)
    x <- x_raw - median(x_raw)
    var_c <- par[1] * par[3]^2 + (1 - par[1]) * par[4]^2
    x <- x / sqrt(var_c)
    
  } else if (dist == "pareto") {
    if (is.null(par)) par <- c(1, 3)
    if (par[2] <= 2) stop("Shape parameter for Pareto must be > 2 for variance to exist")
    x <- VGAM::rpareto(n, scale = par[1], shape = par[2])
    # Center around median
    median_p <- par[1] * 2^(1/par[2])
    var_p <- (par[1]^2 * par[2]) / ((par[2] - 1)^2 * (par[2] - 2))
    x <- (x - median_p) / sqrt(var_p)
    
  } else {
    stop("Unsupported distribution: ", dist)
  }
  
  return(x)
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
    
  } else if(test == "AD2"){  # Anderson-Darling from DescTools (more options)
    output <- DescTools::AndersonDarlingTest(x)
    
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

