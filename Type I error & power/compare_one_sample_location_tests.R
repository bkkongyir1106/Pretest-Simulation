## =========================================================================
## Compare Type I Error & Power for Five One-Sample Tests (Base R)
## - Tests: t-test, Wilcoxon, Bootstrap (mean), Sign (median), Permutation
## - One effect_size controls scenario:
##     effect_size == 0  -> Type I error
##     effect_size != 0  -> Power
## - For each distribution, a panel overlays Type I error (dashed) & Power (solid)
## =========================================================================

## ---- Load dependencies and user functions --------------------------------
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/user_defined_functions_center_median.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

setwd("~/Desktop/OSU/Research/Pretest-Simulation/Type I error & power")

suppressPackageStartupMessages({
  require(foreach)
  require(snow)
  require(doSNOW)
  require(parallel)
})

## ---- Parallel setup helpers ----------------------------------------------
par_set <- function(cores_reserve = 2) {
  cores <- parallel::detectCores()
  cores_use <- max(1, cores - cores_reserve)
  if (Sys.info()[["sysname"]] == "Windows") {
    cl <- parallel::makeCluster(cores_use)
    doSNOW::registerDoSNOW(cl)  
  } else {
    cl <- snow::makeSOCKcluster(cores_use)
    doSNOW::registerDoSNOW(cl)
  }
  cl
}
close_cluster <- function(cl) if (!is.null(cl)) parallel::stopCluster(cl)

## ---- Test functions -------------------------------------------------------
# Bootstrap test for the mean (non-studentized)
one_sample_bootstrap_test <- function(x, mu0 = 0, B = 1000) {
  T_obs <- mean(x) - mu0
  x_centered <- x - mean(x) + mu0
  boot_stats <- replicate(B, mean(sample(x_centered, replace = TRUE)) - mu0)
  mean(abs(boot_stats) >= abs(T_obs))
}

# Sign test for the median
sign_test <- function(x, mu0 = 0) {
  x0 <- x - mu0
  n_valid <- sum(x0 != 0)
  if (n_valid == 0) return(1)  # all ties -> no evidence
  signs <- sum(x0 > 0)
  binom.test(signs, n_valid, p = 0.5)$p.value
}

# Permutation (sign-flip) test on mean
perm_test <- function(x, B = 1000) {
  n <- length(x)
  sx <- sd(x)
  if (!is.finite(sx) || sx == 0) return(1)
  t_obs <- mean(x) / (sx / sqrt(n))
  t_perm <- replicate(B, {
    s <- sample(c(-1, 1), n, replace = TRUE)
    xp <- x * s
    sp <- sd(xp)
    if (!is.finite(sp) || sp == 0) return(0)
    mean(xp) / (sp / sqrt(n))
  })
  mean(abs(t_perm) >= abs(t_obs))
}

## ---- Parameters -----------------------------------------------------------
Nsim         <- 1e3        # repetitions per (n, dist) under H0 and H1
alpha        <- 0.05
B            <- 1e3        # bootstrap / permutation resamples
effect_size  <- 0.50       # 0 => Type I error; nonzero => Power
distributions <- c("Normal", "t", "Chi_Square", "LogNormal")
sample_size   <- c(8, 10, 15, 20, 25, 30, 50)

## ---- Progress bar + cluster ----------------------------------------------
my_cl <- par_set(cores_reserve = 2)
on.exit(close_cluster(my_cl), add = TRUE)

ntasks <- length(sample_size) * length(distributions)
pb <- txtProgressBar(max = ntasks, style = 3)
on.exit(close(pb), add = TRUE)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)

## ---- Simulation under H0 and H1 ------------------------------------------
system.time({
  sim_out <- foreach(i = seq_along(sample_size),
                     .packages = c("stats"),
                     .export = c("generate_data",
                                 "one_sample_bootstrap_test",
                                 "sign_test",
                                 "perm_test",
                                 "alpha", "B", "Nsim", "effect_size",
                                 "distributions", "sample_size")) %:%
    
    foreach(j = seq_along(distributions),
            .options.snow = opts) %dopar% {
              
              n    <- sample_size[i]
              dist <- distributions[j]
              
              ## task-specific RNG seed for reproducibility across workers
              set.seed(12345 + i * 1000 + j)
              
              pval_t_H0 <- pval_wilcox_H0 <- pval_boot_H0 <- pval_sign_H0 <- pval_perm_H0 <- numeric(Nsim)
              pval_t_H1 <- pval_wilcox_H1 <- pval_boot_H1 <- pval_sign_H1 <- pval_perm_H1 <- numeric(Nsim)
              
              for (rep in seq_len(Nsim)) {
                x <- generate_data(n, dist)  # H0: centered at 0 by your generator
                
                # --- Type I error (H0) ---
                pval_t_H0[rep]      <- tryCatch(t.test(x, mu = 0)$p.value, error = function(e) NA_real_)
                pval_wilcox_H0[rep] <- tryCatch(wilcox.test(x, mu = 0)$p.value, error = function(e) NA_real_)
                pval_boot_H0[rep]   <- tryCatch(one_sample_bootstrap_test(x, mu0 = 0, B = B), error = function(e) NA_real_)
                pval_sign_H0[rep]   <- tryCatch(sign_test(x, mu0 = 0), error = function(e) NA_real_)
                pval_perm_H0[rep]   <- tryCatch(perm_test(x, B = B), error = function(e) NA_real_)
                
                # --- Power (H1) ---
                x1 <- x + effect_size
                pval_t_H1[rep]      <- tryCatch(t.test(x1, mu = 0)$p.value, error = function(e) NA_real_)
                pval_wilcox_H1[rep] <- tryCatch(wilcox.test(x1, mu = 0)$p.value, error = function(e) NA_real_)
                pval_boot_H1[rep]   <- tryCatch(one_sample_bootstrap_test(x1, mu0 = 0, B = B), error = function(e) NA_real_)
                pval_sign_H1[rep]   <- tryCatch(sign_test(x1, mu0 = 0), error = function(e) NA_real_)
                pval_perm_H1[rep]   <- tryCatch(perm_test(x1, B = B), error = function(e) NA_real_)
              }
              
              list(
                error_t_test         = mean(pval_t_H0      < alpha, na.rm = TRUE),
                error_wilcox_test    = mean(pval_wilcox_H0 < alpha, na.rm = TRUE),
                error_bootstrap_test = mean(pval_boot_H0   < alpha, na.rm = TRUE),
                error_sign_test      = mean(pval_sign_H0   < alpha, na.rm = TRUE),
                error_perm_test      = mean(pval_perm_H0   < alpha, na.rm = TRUE),
                
                power_t_test         = mean(pval_t_H1      < alpha, na.rm = TRUE),
                power_wilcox_test    = mean(pval_wilcox_H1 < alpha, na.rm = TRUE),
                power_bootstrap_test = mean(pval_boot_H1   < alpha, na.rm = TRUE),
                power_sign_test      = mean(pval_sign_H1   < alpha, na.rm = TRUE),
                power_perm_test      = mean(pval_perm_H1   < alpha, na.rm = TRUE)
              )
            }
  
  ## close cluster once computation finishes
  close_cluster(my_cl)
})

## ---- Collect results into arrays -----------------------------------------
make_arr <- function() array(NA_real_,
                             dim = c(length(sample_size), length(distributions)),
                             dimnames = list(sample_size, distributions))

TypeI_error_t.test         <- make_arr()
TypeI_error_wilcox.test    <- make_arr()
TypeI_error_bootstrap.test <- make_arr()
TypeI_error_sign.test      <- make_arr()
TypeI_error_perm.test      <- make_arr()

power_t.test               <- make_arr()
power_wilcox.test          <- make_arr()
power_bootstrap.test       <- make_arr()
power_sign.test            <- make_arr()
power_perm.test            <- make_arr()

for (i in seq_along(sample_size)) {
  for (j in seq_along(distributions)) {
    TypeI_error_t.test[i, j]         <- sim_out[[i]][[j]]$error_t_test
    TypeI_error_wilcox.test[i, j]    <- sim_out[[i]][[j]]$error_wilcox_test
    TypeI_error_bootstrap.test[i, j] <- sim_out[[i]][[j]]$error_bootstrap_test
    TypeI_error_sign.test[i, j]      <- sim_out[[i]][[j]]$error_sign_test
    TypeI_error_perm.test[i, j]      <- sim_out[[i]][[j]]$error_perm_test
    
    power_t.test[i, j]               <- sim_out[[i]][[j]]$power_t_test
    power_wilcox.test[i, j]          <- sim_out[[i]][[j]]$power_wilcox_test
    power_bootstrap.test[i, j]       <- sim_out[[i]][[j]]$power_bootstrap_test
    power_sign.test[i, j]            <- sim_out[[i]][[j]]$power_sign_test
    power_perm.test[i, j]            <- sim_out[[i]][[j]]$power_perm_test
  }
}

## ---- Optional: print tables ----------------------------------------------
cat("Type I error Rates for t-test\n");         print(TypeI_error_t.test)
cat("Type I error Rates for Wilcoxon-test\n");  print(TypeI_error_wilcox.test)
cat("Type I error Rates for Bootstrap test\n"); print(TypeI_error_bootstrap.test)
cat("Type I error Rates for Sign test\n");      print(TypeI_error_sign.test)
cat("Type I error Rates for Permutation test\n"); print(TypeI_error_perm.test)

cat("Power for t-test\n");         print(power_t.test)
cat("Power for Wilcoxon-test\n");  print(power_wilcox.test)
cat("Power for Bootstrap test\n"); print(power_bootstrap.test)
cat("Power for Sign test\n");      print(power_sign.test)
cat("Power for Permutation test\n"); print(power_perm.test)

## ---- Plot: Type I error (dashed) and Power (solid) in the same panel -----
plot_type1_power_samepanel <- function(sample_size, distributions, alpha,
                                       TypeI_error_t.test, TypeI_error_wilcox.test,
                                       TypeI_error_bootstrap.test, TypeI_error_sign.test, TypeI_error_perm.test,
                                       power_t.test, power_wilcox.test,
                                       power_bootstrap.test, power_sign.test, power_perm.test,
                                       main_title = NULL,
                                       legend_height = 0.30, spacer_height = 0.04) {
  
  methods <- c("t-test", "Wilcoxon", "Bootstrap", "Sign", "Permutation")
  cols <- c("steelblue4", "seagreen4", "darkorange3", "purple4", "firebrick3")
  pchs <- c(16, 17, 15, 18, 8)
  
  type1_list <- list(TypeI_error_t.test, TypeI_error_wilcox.test, TypeI_error_bootstrap.test,
                     TypeI_error_sign.test, TypeI_error_perm.test)
  power_list <- list(power_t.test, power_wilcox.test, power_bootstrap.test,
                     power_sign.test, power_perm.test)
  
  # layout: grid of distribution panels + legend row
  nD <- length(distributions)
  nrow <- ceiling(nD / 2); ncol <- min(2, nD)
  layout_mat <- rbind(
    matrix(seq_len(nrow * ncol), nrow = nrow, byrow = TRUE),
    rep(nrow * ncol + 1, ncol)  # legend row index
  )
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  layout(layout_mat, heights = c(rep(1, nrow), legend_height))
  par(mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  
  for (j in seq_along(distributions)) {
    plot(range(sample_size), c(0, 1), type = "n",
         xlab = "Sample size (n)", ylab = "Rate", main = distributions[j])
    abline(h = alpha, col = "gray50", lty = 3, lwd = 1.5)
    grid(col = "lightgray", lty = "dotted")
    
    for (m in seq_along(methods)) {
      ti <- type1_list[[m]][, j]
      pw <- power_list[[m]][, j]
      lines(sample_size, ti, type = "b", lty = 2, lwd = 2, col = cols[m], pch = pchs[m])
      lines(sample_size, pw, type = "b", lty = 1, lwd = 2, col = cols[m], pch = pchs[m])
    }
  }
  
  # Legend panel (with some spacer)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  # spacer area
  usr <- par("usr")
  rect(usr[1], usr[3], usr[2], usr[3] + (usr[4]-usr[3]) * spacer_height, border = NA, col = NA)
  legend("center",
         legend = c(paste0(methods, " (Type I)"),
                    paste0(methods, " (Power)")),
         col    = c(cols, cols),
         lty    = c(rep(2, length(methods)), rep(1, length(methods))),
         lwd    = 2,
         pch    = c(pchs, pchs),
         ncol   = 2,
         title  = "Five Methods",
         bty    = "n", cex = 0.95)
  
  if (is.null(main_title)) {
    main_title <- if (effect_size == 0)
      "Type I Error vs Sample Size (five methods)"
    else
      sprintf("Type I Error (dashed) & Power (solid), effect_size = %.2f", effect_size)
  }
  mtext(main_title, outer = TRUE, cex = 1.2, font = 2, line = 0.2)
}

## ---- Draw the figure ------------------------------------------------------
plot_type1_power_samepanel(
  sample_size, distributions, alpha,
  TypeI_error_t.test, TypeI_error_wilcox.test, TypeI_error_bootstrap.test, TypeI_error_sign.test, TypeI_error_perm.test,
  power_t.test, power_wilcox.test,  power_bootstrap.test, power_sign.test, power_perm.test
)

## ---- Optional: save to PDF -----------------------------------------------
pdf("TypeI_Power_FiveMethods.pdf", width = 10, height = 8)
plot_type1_power_samepanel(
  sample_size, distributions, alpha,
  TypeI_error_t.test, TypeI_error_wilcox.test, TypeI_error_bootstrap.test, TypeI_error_sign.test, TypeI_error_perm.test,
  power_t.test,       power_wilcox.test,       power_bootstrap.test,       power_sign.test,       power_perm.test
)
dev.off()

# ================ power and ERROR =====
## =========================
## Separate plots: POWER and TYPE I ERROR (base R)
## =========================

# Shared styling
.methods <- c("t-test", "Wilcoxon", "Bootstrap", "Sign", "Permutation")
.cols    <- c("steelblue4", "seagreen4", "green", "purple4", "firebrick3")
.pchs    <- c(16, 17, 15, 18, 8)

# Generic panel plotter (one metric at a time)
plot_metric_panels <- function(sample_size, distributions, matrices_list, methods = .methods,
                               cols = .cols, pchs = .pchs,
                               main_title = "Metric vs Sample Size",
                               ylab = "Rate", ylim = c(0, 1),
                               alpha_line = NULL, legend_title = "Five Methods",
                               legend_height = 0.28) {
  
  stopifnot(length(matrices_list) == length(methods),
            length(cols) == length(methods),
            length(pchs) == length(methods))
  
  # Layout: grid of distribution panels + dedicated legend row
  nD   <- length(distributions)
  nrow <- ceiling(nD / 2)
  ncol <- min(2, nD)
  
  layout_mat <- rbind(
    matrix(seq_len(nrow * ncol), nrow = nrow, byrow = TRUE),
    rep(nrow * ncol + 1, ncol) # legend row
  )
  
  old_par <- par(no.readonly = TRUE)
  on.exit(par(old_par), add = TRUE)
  
  layout(layout_mat, heights = c(rep(1, nrow), legend_height))
  par(mar = c(4, 4, 3, 1), oma = c(0, 0, 2, 0))
  
  # Panels
  for (j in seq_along(distributions)) {
    plot(range(sample_size), ylim, type = "n",
         xlab = "Sample size (n)", ylab = ylab, main = distributions[j])
    if (!is.null(alpha_line)) abline(h = alpha_line, col = "gray50", lty = 3, lwd = 1.5)
    grid(col = "lightgray", lty = "dotted")
    
    for (m in seq_along(methods)) {
      y <- matrices_list[[m]][, j]
      lines(sample_size, y, type = "b", lwd = 2, col = cols[m], pch = pchs[m])
    }
  }
  
  # Legend row
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
         legend = methods,
         col    = cols,
         pch    = pchs,
         lty    = 1,
         lwd    = 2,
         ncol   = 5,
         bty    = "n",
         title  = legend_title,
         cex    = 0.95)
  
  mtext(main_title, outer = TRUE, cex = 1.2, font = 2, line = 0.2)
}

# ---------- Power plot ----------
plot_power_panels <- function(sample_size, distributions,
                              power_t.test, power_wilcox.test, power_bootstrap.test,
                              power_sign.test, power_perm.test,
                              effect_size = NA_real_) {
  
  power_list <- list(power_t.test, power_wilcox.test, power_bootstrap.test,
                     power_sign.test, power_perm.test)
  
  ttl <- if (is.na(effect_size)) "Power vs Sample Size"
  else sprintf("Power vs Sample Size (effect_size = %.2f)", effect_size)
  
  plot_metric_panels(sample_size, distributions, power_list,
                     methods = .methods, cols = .cols, pchs = .pchs,
                     main_title = ttl, ylab = "Power", ylim = c(0, 1),
                     alpha_line = NULL, legend_title = "Five Methods")
}

# ---------- Type I error plot ----------
plot_type1_panels <- function(sample_size, distributions, alpha,
                              TypeI_error_t.test, TypeI_error_wilcox.test,
                              TypeI_error_bootstrap.test, TypeI_error_sign.test,
                              TypeI_error_perm.test) {
  
  type1_list <- list(TypeI_error_t.test, TypeI_error_wilcox.test, TypeI_error_bootstrap.test,
                     TypeI_error_sign.test, TypeI_error_perm.test)
  
  plot_metric_panels(sample_size, distributions, type1_list,
                     methods = .methods, cols = .cols, pchs = .pchs,
                     main_title = "Type I Error vs Sample Size",
                     ylab = "Type I error rate", ylim = c(0, 0.75),
                     alpha_line = alpha, legend_title = "Five Methods")
}

## =========================
## Draw the two figures
## =========================

# POWER
# dev.new(width = 10, height = 8)  # optional: open a second device
plot_power_panels(
  sample_size, distributions,
  power_t.test, power_wilcox.test, power_bootstrap.test, power_sign.test, power_perm.test,
  effect_size = effect_size
)

# TYPE I ERROR
# dev.new(width = 10, height = 8)  # optional: open another device
plot_type1_panels(
  sample_size, distributions, alpha,
  TypeI_error_t.test, TypeI_error_wilcox.test, TypeI_error_bootstrap.test, TypeI_error_sign.test, TypeI_error_perm.test
)

## =========================
## Optional: save to PDF
## =========================
pdf("Power_FiveMethods.pdf", width = 10, height = 8)
plot_power_panels(sample_size, distributions,
                  power_t.test, power_wilcox.test, power_bootstrap.test,
                  power_sign.test, power_perm.test,
                  effect_size = effect_size)
dev.off()

pdf("TypeI_FiveMethods.pdf", width = 10, height = 8)
plot_type1_panels(sample_size, distributions, alpha,
                  TypeI_error_t.test, TypeI_error_wilcox.test,
                  TypeI_error_bootstrap.test, TypeI_error_sign.test, TypeI_error_perm.test)
dev.off()

