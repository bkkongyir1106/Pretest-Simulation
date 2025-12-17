# ---------- General Box–Cox + t-test Type I error (Exp, Gamma, Chi-square) ----------

# Box-Cox transform
bc_transform <- function(x, lambda) {
  if (abs(lambda) < 1e-12) log(x) else (x^lambda - 1) / lambda
}

# Numerical null mean: E[ BC(X + delta; lambda) ] for chosen distribution
# dist_type: "exp", "gamma", "chisq"
# params: list of parameters:
#   - exp:   list(rate)
#   - gamma: list(shape, rate)
#   - chisq: list(df)
bc_null_mean_numeric <- function(dist_type, params, delta, lambda) {
  # pdf over [0, ∞)
  pdf_fun <- switch(
    dist_type,
    "exp"   = function(x) dexp(x, rate = params$rate),
    "gamma" = function(x) dgamma(x, shape = params$shape, rate = params$rate),
    "chisq" = function(x) dchisq(x, df = params$df),
    stop("Unsupported dist_type")
  )
  
  if (abs(lambda) < 1e-12) {
    # E[ log(X + delta) ]
    integrand <- function(x) log(x + delta) * pdf_fun(x)
  } else {
    # E[ ((X + delta)^lambda - 1) / lambda ]
    integrand <- function(x) (( (x + delta)^lambda - 1 ) / lambda) * pdf_fun(x)
  }
  
  val <- integrate(integrand, lower = 0, upper = Inf,
                   rel.tol = 1e-8, abs.tol = 0)$value
  return(val)
}

# One replicate: generate data, shift if min<0, ensure positivity for log case,
# apply Box-Cox, compute correct null mean (with same pre-transform offset),
# one-sample t-test vs that null, return reject indicator.
replicate_once <- function(n, dist_type, params, lambda, alpha,
                           tiny_eps = 1e-8) {
  # RNG for X
  x <- switch(
    dist_type,
    "exp"   = rexp(n, rate = params$rate),
    "gamma" = rgamma(n, shape = params$shape, rate = params$rate),
    "chisq" = rchisq(n, df = params$df),
    stop("Unsupported dist_type")
  )
  
  # Per instruction: only shift if min is negative (these are ≥0 in theory,
  # but we follow the rule literally).
  shift_c <- if (min(x) < 0) -min(x) else 0
  x_shifted <- x + shift_c
  
  # For lambda=0 (log), ensure strictly positive arguments; add tiny eps if needed.
  delta_eps <- if (abs(lambda) < 1e-12 && min(x_shifted) <= 0) tiny_eps else 0
  
  # Total pre-transform offset included in null mean
  delta_total <- shift_c + delta_eps
  
  # Apply Box–Cox
  y <- bc_transform(x_shifted + delta_eps, lambda)
  
  # Correct null mean under the DGP with the *same* offset
  mu0 <- bc_null_mean_numeric(dist_type, params, delta_total, lambda)
  
  # One-sample t-test (Welch is equivalent here)
  pval <- t.test(y, mu = mu0)$p.value
  as.numeric(pval < alpha)
}

# Master wrapper to estimate Type I error by Monte Carlo
type1_error_boxcox <- function(n = 30,
                               dist_type = c("exp","gamma","chisq"),
                               params = list(rate = 1),
                               lambda = 0.5,
                               alpha = 0.05,
                               nrep = 5000,
                               seed = 123) {
  dist_type <- match.arg(dist_type)
  set.seed(seed)
  rejs <- replicate(nrep, replicate_once(n, dist_type, params, lambda, alpha))
  mean(rejs)
}

# ---------- Examples ----------

# 1) Exponential(rate=1)
t1_exp_lam05 <- type1_error_boxcox(n = 10,
                                   dist_type = "exp",
                                   params = list(rate = 1),
                                   lambda = 0.5,
                                   alpha = 0.05,
                                   nrep = 3000)

# 2) Gamma(shape=2, rate=1)  (mean = shape/rate = 2)
t1_gamma_lam0 <- type1_error_boxcox(n = 10,
                                    dist_type = "gamma",
                                    params = list(shape = 2, rate = 1),
                                    lambda = 0.0,     # log transform
                                    alpha = 0.05,
                                    nrep = 3000)

# 3) Chi-square(df=5)  (Gamma(shape=df/2, rate=1/2))
t1_chisq_lam1 <- type1_error_boxcox(n = 10,
                                    dist_type = "chisq",
                                    params = list(df = 5),
                                    lambda = 1.0,     # BC(x;1)=x-1
                                    alpha = 0.05,
                                    nrep = 3000)

cat(sprintf("Type I error | Exp(rate=1), lambda=0.5:  %.3f\n", t1_exp_lam05))
cat(sprintf("Type I error | Gamma(2,1),  lambda=0.0:  %.3f\n", t1_gamma_lam0))
cat(sprintf("Type I error | ChiSq(df=5), lambda=1.0:  %.3f\n", t1_chisq_lam1))
