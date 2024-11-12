setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/NORMALITY TEST METHODS")

# Install necessary packages if not already installed
if (!require("nortest")) install.packages("nortest", dependencies = TRUE)
if (!require("tseries")) install.packages("tseries", dependencies = TRUE)
if (!require("goftest")) install.packages("goftest", dependencies = TRUE)

# Load libraries
library(nortest)     # For Anderson-Darling, Cramer-Von Mises
library(tseries)     # For Jarque-Bera
library(goftest)     # For Kolmogorov-Smirnov

# Parameters
N = 1e4       # Number of repetitions for each alpha level to average TPR and FPR
n = 10        # Sample size for each test
alpha.level = seq(from = 0, to = 1, by = 0.01)  # Range of significance levels

# Define Laplace function (if not available)
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n) - 0.5
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

# Initialize vectors to store False Positive Rates (FPR) and True Positive Rates (TPR)
FPR = numeric(length(alpha.level))  # FPR for normal distribution
TPR_exp = TPR_chisq = TPR_lognorm = TPR_unif = TPR_laplace = numeric(length(alpha.level))

# Loop over each significance level
for(i in seq_along(alpha.level)) {
  alpha = alpha.level[i]

  # Vectors to store p-values for each repetition for each distribution
  pval_norm = pval_exp = pval_chisq = pval_lognorm = pval_unif = pval_laplace = numeric(N)

  # Repeat the test N times for each alpha level
  for(j in 1:N) {
    # Generate new samples for each test
    x_norm = rnorm(n = n, mean = 0, sd = 1)      # Normal sample
    x_exp = rexp(n = n, rate = 1)                # Exponential sample
    x_chisq = rchisq(n = n, df = 2)              # Chi-squared sample
    x_lognorm = rlnorm(n = n, meanlog = 0, sdlog = 1)  # Lognormal sample
    x_unif = runif(n = n, min = 0, max = 1)      # Uniform sample
    x_laplace = rlaplace(n = n, location = 0, scale = 1) # Laplace sample

    # Perform the Shapiro-Wilk test and store p-values
    pval_norm[j] = shapiro.test(x_norm)$p.value
    pval_exp[j] = shapiro.test(x_exp)$p.value
    pval_chisq[j] = shapiro.test(x_chisq)$p.value
    pval_lognorm[j] = shapiro.test(x_lognorm)$p.value
    pval_unif[j] = shapiro.test(x_unif)$p.value
    pval_laplace[j] = shapiro.test(x_laplace)$p.value
  }

  # Calculate False Positive Rate (FPR) for normal samples and TPRs for other distributions
  FPR[i] = mean(pval_norm < alpha)         # Proportion of normal samples rejected (false positives)
  TPR_exp[i] = mean(pval_exp < alpha)      # TPR for exponential samples
  TPR_chisq[i] = mean(pval_chisq < alpha)  # TPR for chi-squared samples
  TPR_lognorm[i] = mean(pval_lognorm < alpha)  # TPR for lognormal samples
  TPR_unif[i] = mean(pval_unif < alpha)    # TPR for uniform samples
  TPR_laplace[i] = mean(pval_laplace < alpha)  # TPR for Laplace samples
}
# save data
#save(FPR, TPR_exp, TPR_chisq,TPR_lognorm, TPR_unif, TPR_laplace, file = "ROC_SW.test.RData")

load("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/NORMALITY TEST METHODS/ROC_SW.test.RData")
# Plot the ROC curves
plot(FPR, TPR_exp, type = "l", col = "blue", lwd = 2, xlab = "False Positive Rate (FPR)", ylab = "True Positive Rate (TPR)", main = "ROC Curve for Shapiro-Wilk Test")
lines(FPR, TPR_chisq, col = "green", lwd = 2)
lines(FPR, TPR_lognorm, col = "purple", lwd = 2)
lines(FPR, TPR_unif, col = "orange", lwd = 2)
lines(FPR, TPR_laplace, col = "brown", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "red")  # Diagonal reference line

# Add a legend
legend("bottomright", legend = c("Exponential", "Chi-Squared", "Lognormal", "Uniform", "Laplace"), col = c("blue", "green", "purple", "orange", "brown"), lwd = 2)

# %%%%%%%%%%%%%% Compare Different Normality Test Methods %%%%%%%%%%%%%%%%%% #

# Define Laplace function (if not available)
rlaplace <- function(n, location = 0, scale = 1) {
  u <- runif(n) - 0.5
  location - scale * sign(u) * log(1 - 2 * abs(u))
}

# Parameters
N = 1e4       # Number of repetitions for each alpha level to average TPR and FPR
n = 10        # Sample size for each test
alpha.level = seq(from = 0, to = 1, by = 0.01)  # Range of significance levels

# Initialize lists to store False Positive Rates (FPR) and True Positive Rates (TPR) for each test
tests = c("SW", "KS", "AD", "DAP", "SF", "JB", "CVM")
FPR_list = TPR_list = list()
for (test in tests) {
  FPR_list[[test]] = numeric(length(alpha.level))
  TPR_list[[test]] = numeric(length(alpha.level))
}

# Function to run each normality test and return p-value
get_p_value <- function(x, test) {
  switch(test,
         "SW" = shapiro.test(x)$p.value,
         "KS" = ks.test(x, "pnorm", mean(x), sd(x))$p.value,
         "AD" = ad.test(x)$p.value,
         "DAP" = agostino.test(x)$p.value,  # D'Agostino test
         "SF" = nortest::sf.test(x)$p.value,
         "JB" = jarque.bera.test(x)$p.value,
         "CVM" = cvm.test(x)$p.value)
}

# Loop over each significance level
for(i in seq_along(alpha.level)) {
  alpha = alpha.level[i]

  # Repeat the test N times for each alpha level
  for (test in tests) {
    pval_norm = pval_exp = pval_chisq = pval_lognorm = pval_unif = pval_laplace = numeric(N)

    for (j in 1:N) {
      # Generate samples
      x_norm = rnorm(n = n, mean = 0, sd = 1)            # Normal sample
      #x_exp = rexp(n = n, rate = 1)                      # Exponential sample
      #x_chisq = rchisq(n = n, df = 2)                    # Chi-squared sample
      x_lognorm = rlnorm(n = n, meanlog = 0, sdlog = 1)  # Lognormal sample
      # x_unif = runif(n = n, min = 0, max = 1)            # Uniform sample
      # x_laplace = rlaplace(n = n, location = 0, scale = 1)  # Laplace sample

      # Get p-values for each sample for the current test
      pval_norm[j] = get_p_value(x_norm, test)
      #pval_exp[j] = get_p_value(x_exp, test)
      #pval_chisq[j] = get_p_value(x_chisq, test)
      pval_lognorm[j] = get_p_value(x_lognorm, test)
      # pval_unif[j] = get_p_value(x_unif, test)
      # pval_laplace[j] = get_p_value(x_laplace, test)
    }

    # Calculate FPR and TPR for the current test and alpha level
    FPR_list[[test]][i] = mean(pval_norm < alpha)                    # FPR (False Positive Rate)
    TPR_list[[test]][i] = mean(pval_lognorm < alpha)  # TPR (True Positive Rate)
  }
}
#save data
#save(FPR_list, TPR_list, file = "ROC_Normality.test.RData")

load("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/NORMALITY TEST METHODS/ROC_Normality.test.RData")
# Plot the ROC curves for each test
tests = c("SW", "KS", "AD", "DAP", "SF", "JB", "CVM")
plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1), xlab = "False Positive Rate (FPR)", ylab = "True Positive Rate (TPR)", main = "ROC Curves for Different Normality Tests")
colors = rainbow(length(tests))
for (i in seq_along(tests)) {
  lines(FPR_list[[tests[i]]], TPR_list[[tests[i]]], col = colors[i], lwd = 2, type = "l")
}
abline(a = 0, b = 1, lty = 2, col = "black")  # Diagonal reference line
legend("bottomright", legend = tests, col = colors, lwd = 2, title = "Normality Tests")


### ROC Curves for t-test vs two-stage method data generation function

# # Parameters
N <- 1e3             # Increased number of repetitions for smoother ROC curve
n <- 10              # Sample size for each test
d <- 0.5             # Effect size
sw.alpha <- 0.05     # Significance level for the Shapiro-Wilk test

# Range of significance levels for ROC curve
alpha.level <- seq(0.001, 0.8, 0.1)

# Initialize vectors to store False Positive Rates (FPR) and True Positive Rates (TPR)
FPR_t <- TPR_t <- FPR_w_t <- TPR_w_t <- numeric(length(alpha.level))

# Loop over each significance level
for (i in seq_along(alpha.level)) {
  alpha <- alpha.level[i]

  # Initialize vectors for storing p-values
  pval_error_t <- pval_error_w_t <- numeric(N)
  pval_powr_t <- pval_powr_w_t <- numeric(N)

  # Generate data for N samples for each alpha level
  for (j in 1:N) {
    x <- rexp(n, rate = 1)
    y <- rexp(n, rate = 1)

    # Calculate p-values under null hypothesis (no effect)
    pval_error_t[j] <- t.test(x - 1, y - 1)$p.value

    # Calculate p-values under alternative hypothesis (with effect)
    pval_powr_t[j] <- t.test(x - 1, y - 1 + d)$p.value

    # Test for normality and use appropriate test
    if (shapiro.test(x)$p.value > sw.alpha && shapiro.test(y)$p.value > sw.alpha) {
      # Use t-test if both samples pass Shapiro-Wilk normality test
      pval_error_w_t[j] <- t.test(x - 1, y - 1)$p.value
      pval_powr_w_t[j] <- t.test(x - 1, y - 1 + d)$p.value
    } else {
      # Use Wilcoxon test if normality assumption is violated
      pval_error_w_t[j] <- wilcox.test(x - 1, y - 1)$p.value
      pval_powr_w_t[j] <- wilcox.test(x - 1, y - 1 + d)$p.value
    }
  }

  # Calculate FPR and TPR for each test at the current alpha level
  FPR_t[i] <- mean(pval_error_t < alpha)
  TPR_t[i] <- mean(pval_powr_t < alpha)
  FPR_w_t[i] <- mean(pval_error_w_t < alpha)
  TPR_w_t[i] <- mean(pval_powr_w_t < alpha)
}

save(FPR_t, TPR_t, FPR_w_t, TPR_w_t, file = "t.test_vs_two_stage_ROC.RData")

# load data 
load("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Fall 2024/ROC Curves/t.test_vs_two_stage_ROC.RData")
# Plot the ROC curves with smoothing
plot(FPR_t, TPR_t, type = "l", col = "blue", lwd = 2, ylim = c(0, 1),
     xlab = "False Positive Rate (FPR)", ylab = "True Positive Rate (TPR)",
     main = "ROC Curve for t vs two-stage method for Exponential Data")
lines(FPR_w_t, TPR_w_t, type = "l", col = "red", lwd = 2)
abline(a = 0, b = 1, lty = 2, col = "grey")  # Diagonal reference line for random guessing

# Add a legend
legend("bottomright", title = "test type", legend = c("t-test", "t/Wilcoxon test"), col = c("blue", "red"), lwd = 2)