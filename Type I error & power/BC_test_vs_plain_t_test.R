# ==========================================================
# Compare Type I error: Corrected Box–Cox test vs Plain t-test
# Families: Gamma, Lognormal, Weibull, Normal
# Panel plot: 2 x 2 (one panel per family)
# ==========================================================
setwd("~/Desktop/OSU/Research/Pretest-Simulation/Type I error & power")
# --- deps ---
library(MASS)   # boxcox
library(stats)

set.seed(123)

# --- global controls ---
alpha  <- 0.05
reps   <- 1e5
ns     <- c(10, 20, 30, 40, 50)

# --- scientific null on ORIGINAL scale ---
m0 <- 1.0

# --- family parameters (fixed under H0) ---
# Gamma(shape=k_g, scale=theta_g) with mean m0 = k_g * theta_g
k_g       <- 1          # Exponential(m0)
# Lognormal: mean m0 => mu = log(m0) - sigma2/2
sigma2_ln <- 0.5
# Weibull(shape=k_w, scale=theta_w) with mean m0 = theta_w * Gamma(1 + 1/k_w)
k_w       <- 1.5
# Normal(mean=m0, sd=sd_norm)
sd_norm   <- 1.0

# --- numerical safety for Box–Cox ---
eps_clip <- 1e-4

# --- Box–Cox helper ---
g_boxcox <- function(x, lambda) {
  if (abs(lambda) < 1e-8) log(x) else (x^lambda - 1) / lambda
}

# ===============================
# Correct transformed null E[g_lambda(X)] under H0
# ===============================
# Gamma(k, theta=m0/k)
mu0_gamma <- function(lambda, m0, k) {
  if (abs(lambda) < 1e-8) {
    return(digamma(k) + log(m0 / k))
  } else {
    if (lambda <= -k) return(NA_real_)  # existence: lambda > -k
    theta <- m0 / k
    term <- (theta^lambda) * exp(lgamma(k + lambda) - lgamma(k))
    return((term - 1) / lambda)
  }
}

# Lognormal(mu, sigma2), mu = log(m0) - sigma2/2
mu0_lognormal <- function(lambda, m0, sigma2) {
  if (abs(lambda) < 1e-8) {
    return(log(m0) - sigma2 / 2)
  } else {
    mu <- log(m0) - sigma2 / 2
    term <- exp(lambda * mu + 0.5 * lambda^2 * sigma2)  # E[X^lambda]
    return((term - 1) / lambda)
  }
}

# Weibull(k, theta) with mean m0 = theta * Gamma(1 + 1/k)
mu0_weibull <- function(lambda, m0, k) {
  theta <- m0 / gamma(1 + 1 / k)
  if (abs(lambda) < 1e-8) {
    gamma_const <- -digamma(1)  # Euler's constant
    return(log(theta) - gamma_const / k)
  } else {
    # existence: 1 + lambda/k > 0  => lambda > -k
    if (lambda <= -k) return(NA_real_)
    term <- (theta^lambda) * gamma(1 + lambda / k)      # E[X^lambda]
    return((term - 1) / lambda)
  }
}

# ===============================
# Data generators under H0
# ===============================
rH0_gamma     <- function(n, m0, k) {
  theta <- m0 / k
  rgamma(n, shape = k, scale = theta)
}

rH0_lognormal <- function(n, m0, sigma2) {
  mu <- log(m0) - sigma2 / 2
  rlnorm(n, meanlog = mu, sdlog = sqrt(sigma2))
}

rH0_weibull   <- function(n, m0, k) {
  theta <- m0 / gamma(1 + 1 / k)
  rweibull(n, shape = k, scale = theta)
}

rH0_normal    <- function(n, m0, sd_norm) {
  rnorm(n, mean = m0, sd = sd_norm)
}

# ===============================
# Generic Box–Cox mean-null test (for Gamma/Lognormal/Weibull)
# ===============================
test_boxcox_mean_generic <- function(x, mu0_fun, mu0_args, lambda_boundary = -Inf,
                                     alpha = 0.05, eps_clip = 1e-4) {
  stopifnot(all(is.finite(x)), all(x > 0))  # Box–Cox needs positivity
  
  # Estimate lambda by profile likelihood
  bc <- boxcox(x ~ 1, plotit = FALSE)
  lambda_hat <- bc$x[which.max(bc$y)]
  
  # Enforce boundary if needed (Gamma/Weibull: lambda > -k)
  lambda_used <- max(lambda_hat, lambda_boundary + eps_clip)
  
  gx <- g_boxcox(x, lambda_used)
  
  # Correct transformed null
  mu0_t <- do.call(mu0_fun, c(list(lambda = lambda_used), mu0_args))
  
  # Fallback to log transform if mapping failed
  if (!is.finite(mu0_t)) {
    lambda_used <- 0
    gx <- log(x)
    mu0_t <- do.call(mu0_fun, c(list(lambda = 0), mu0_args))
  }
  
  tt <- t.test(gx, mu = mu0_t)
  list(lambda_hat = lambda_hat, lambda_used = lambda_used, p = tt$p.value)
}

# ===============================
# Comparison per family: corrected vs plain t-test
# ===============================
estimate_type1_comparison <- function(family = c("Gamma","Lognormal","Weibull","Normal"),
                                      ns, reps, m0, k_g, sigma2_ln, k_w, sd_norm,
                                      alpha = 0.05, eps_clip = 1e-4) {
  family <- match.arg(family)
  out <- vector("list", length(ns))
  
  for (i in seq_along(ns)) {
    n <- ns[i]
    rej_corr <- logical(reps)
    rej_plain <- logical(reps)
    
    for (b in seq_len(reps)) {
      if (family == "Gamma") {
        x <- rH0_gamma(n, m0, k_g)
        # corrected (Box–Cox mapped mean-null)
        res <- test_boxcox_mean_generic(
          x,
          mu0_fun   = mu0_gamma,
          mu0_args  = list(m0 = m0, k = k_g),
          lambda_boundary = -k_g,
          alpha = alpha, eps_clip = eps_clip
        )
        rej_corr[b] <- (res$p < alpha)
        # plain t-test on original scale
        rej_plain[b] <- (t.test(x, mu = m0)$p.value < alpha)
        
      } else if (family == "Lognormal") {
        x <- rH0_lognormal(n, m0, sigma2_ln)
        res <- test_boxcox_mean_generic(
          x,
          mu0_fun   = mu0_lognormal,
          mu0_args  = list(m0 = m0, sigma2 = sigma2_ln),
          lambda_boundary = -Inf,
          alpha = alpha, eps_clip = eps_clip
        )
        rej_corr[b] <- (res$p < alpha)
        rej_plain[b] <- (t.test(x, mu = m0)$p.value < alpha)
        
      } else if (family == "Weibull") {
        x <- rH0_weibull(n, m0, k_w)
        res <- test_boxcox_mean_generic(
          x,
          mu0_fun   = mu0_weibull,
          mu0_args  = list(m0 = m0, k = k_w),
          lambda_boundary = -k_w,
          alpha = alpha, eps_clip = eps_clip
        )
        rej_corr[b] <- (res$p < alpha)
        rej_plain[b] <- (t.test(x, mu = m0)$p.value < alpha)
        
      } else if (family == "Normal") {
        x <- rH0_normal(n, m0, sd_norm)
        # For Normal, the “above” test is the plain t-test; both coincide
        rej_corr[b]  <- (t.test(x, mu = m0)$p.value < alpha)
        rej_plain[b] <- (t.test(x, mu = m0)$p.value < alpha)
      }
    }
    
    out[[i]] <- data.frame(
      family = family,
      n = n,
      reps = reps,
      type1_corrected = mean(rej_corr),
      type1_plain     = mean(rej_plain)
    )
  }
  
  do.call(rbind, out)
}

# ===============================
# Run comparisons for all families
# ===============================
cmp_gamma   <- estimate_type1_comparison("Gamma",    ns, reps, m0, k_g, sigma2_ln, k_w, sd_norm, alpha, eps_clip)
cmp_lognorm <- estimate_type1_comparison("Lognormal",ns, reps, m0, k_g, sigma2_ln, k_w, sd_norm, alpha, eps_clip)
cmp_weibull <- estimate_type1_comparison("Weibull",  ns, reps, m0, k_g, sigma2_ln, k_w, sd_norm, alpha, eps_clip)
cmp_normal  <- estimate_type1_comparison("Normal",   ns, reps, m0, k_g, sigma2_ln, k_w, sd_norm, alpha, eps_clip)

# create the plots
pdf("bc_transform_vs_plain_t_test_typeI_err.pdf", width = 6, height = 6)
plot_family_panel_noleg <- function(df, family_name, alpha, ylim_pad = 0.02) {
  df <- df[order(df$n), ]
  y_all <- c(df$type1_corrected, df$type1_plain, alpha)
  y_min <- max(0, min(y_all) - ylim_pad)
  y_max <- min(1, max(y_all) + ylim_pad)
  
  plot(df$n, df$type1_corrected, type = "o", lwd = 2, pch = 16,
       col = "blue",
       xlab = "Sample size (n)", ylab = "Type I error rate",
       ylim = c(y_min, y_max), main = family_name)
  lines(df$n, df$type1_plain, type = "o", lwd = 2, pch = 17, col = "red", lty = 2)
  abline(h = alpha, lty = 3)
}

# ---- 2x2 panels + bottom legend + common title ----
op <- par(no.readonly = TRUE)

# Give outer margins so the common title has space (top = 3 lines here)
par(oma = c(0, 0, 3, 0))

# layout: 3 rows (top two rows are panels; bottom row is the legend spanning both columns)
layout(matrix(c(1,2,3,4,5,5), nrow = 3, byrow = TRUE), heights = c(1, 1, 0.28))

# Panel margins
par(mar = c(4,4,2,1))
plot_family_panel_noleg(cmp_gamma,   "Gamma",    alpha)
plot_family_panel_noleg(cmp_lognorm, "Lognormal",alpha)
plot_family_panel_noleg(cmp_weibull, "Weibull",  alpha)
plot_family_panel_noleg(cmp_normal,  "Normal",   alpha)

# Common title (suptitle)
mtext("Type I Error: Corrected Box–Cox vs Plain t-test", side = 3, outer = TRUE, line = 1, cex = 1.2, font = 2)

# Bottom legend panel
par(mar = c(0,0,0,0))
plot.new()
legend("center",
       legend = c("Corrected Box–Cox test", "Plain t-test", sprintf("alpha = %.2f", alpha)),
       col    = c("blue","red","black"),
       lty    = c(1,2,3),
       lwd    = c(2,2,2),
       pch    = c(16,17,NA),
       horiz  = TRUE,
       bty    = "n",
       xpd    = NA,
       cex    = 1)

par(op)  # restore
dev.off()

