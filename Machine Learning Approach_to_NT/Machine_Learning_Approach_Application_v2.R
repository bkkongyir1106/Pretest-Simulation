setwd("~/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach_to_NT")
# Add Dependencies Check
required_packages <- c("LaplacesDemon", "VGAM")
check_packages <- function() {
  missing <- setdiff(required_packages, rownames(installed.packages()))
  if (length(missing) > 0) {
    stop("Please install required packages: ", paste(missing, collapse = ", "))
  }
}

# loard ML models
load("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach_to_NT/results/trained_ml_model_n10.RData")

# ----------------------------------------------------
#   Generate standardized data from various distributions
#                 Center around median 
# ----------------------------------------------------
generate_data <- function(n, dist, par = NULL) {
  # Input validation
  if (!is.numeric(n) || n <= 0) stop("n must be a positive integer")
  dist <- tolower(dist)
  
  if (dist == "normal") {
    if (is.null(par)) par <- c(0, 1)
    x <- rnorm(n, mean = par[1], sd = par[2])
    x <- (x - par[1]) / par[2]
    
  } else if (dist == "chi_square") {
    if (is.null(par)) par <- 3
    x <- rchisq(n, df = par)
    median_chi <- par *(1-2/(9*par))^3
    sd_chi <- sqrt(2*par)
    x <- (x - median_chi)/sd_chi 
    
  } else if (dist == "gamma") {
    if (is.null(par)) par <- c(3, 0.1)
    x <- rgamma(n, shape = par[1], rate = par[2])
    median_g <- qgamma(0.5, shape = par[1], rate = par[2])
    sd_g <- sqrt(par[1] / par[2]^2)
    x <- (x - median_g)/sd_g 
    
  } else if (dist == "exponential") {
    if (is.null(par)) par <- 1
    x <- rexp(n, rate = par)
    median_e <- log(2)/par
    sd_e <- 1 / par
    x <- (x - median_e)/sd_e
    
  } else if (dist == "t") {
    if (is.null(par)) par <- 3
    if (par <= 2) stop("Degrees of freedom must be > 1")
    x <- rt(n, df = par)
    x <- x / sqrt(par / (par - 2))
    
  } else if (dist == "uniform") {
    if (is.null(par)) par <- c(0, 1)
    x <- runif(n, min = par[1], max = par[2])
    median_u <- (par[1] + par[2]) / 2
    sd_u <- sqrt((par[2] - par[1])^2 / 12)
    x <- (x - median_u)/sd_u
    
  } else if (dist == "laplace") {
    if (is.null(par)) par <- c(3, 1)
    x <- extraDistr::rlaplace(n, mu = par[1], sigma = par[2])
    x <- (x - par[1])/sqrt(2 * par[2]^2)
    
  } else if (dist == "weibull") {
    if (is.null(par)) par <- c(1, 2)
    x <- rweibull(n, shape = par[1], scale = par[2])
    median_w <- par[2] * (log(2))^(1/par[1])
    var_w <- par[2]^2 * (gamma(1 + 2 / par[1]) - gamma(1 + 1 / par[1])^2)
    x <- (x - median_w) / sqrt(var_w)
    
  } else if (dist == "lognormal") {
    if (is.null(par)) par <- c(0, 1)
    x <- rlnorm(n, meanlog = par[1], sdlog = par[2])
    median_ln <- exp(par[1])
    var_ln <- (exp(par[2]^2) - 1) * exp(2 * par[1] + par[2]^2)
    x <- (x - median_ln)/sqrt(var_ln) 
    
  }  else {
    stop("Unsupported distribution: ", dist)
  }
  
  return(x)
}


# ------------------------------------
# Define simple calculation functions 
# -------------------------------------

# compute the area under the curve 
compute_area <- function(x, y) {
  sum(diff(x) * (head(y, -1) + tail(y, -1)) / 2) / (max(x) - min(x))
}

# ------------ framework functions -----------------------------------------
# generate data 
gen_data <- function(n = 20, effect_size = 0.0, sd = 1, dist = "Exponential", par = NULL) {
  dist <- tolower(dist)
  x <- effect_size + sd * generate_data(n, dist, par = par)
  return(x)
}

# get parameters
get_parameters <- function(n, effect_size = 0.5, sd = 1, dist = "exponential", par = NULL, ...) {
  list(
    n = n,
    effect_size = effect_size,
    sd = sd,
    dist = dist,
    par = par
  )
}

# get normality object
fn_to_get_norm_obj <- function(data) {
  return(data)
}

# ds test function 1: one-sample t-test
fn_for_ds_test_1 <- function(data) {
  test_result <- t.test(data, mu = 0)
  return(list(p.value = test_result$p.value))
}


# # ds test function 2: Wilcoxon Sign Rank test 
fn_for_ds_test_2 <- function(data, mu0 = 0) {
  x0 <- data - mu0
  n_valid <- sum(x0 != 0)
  if (n_valid == 0) return(list(p.value = 1))
  signs <- sum(x0 > 0)
  return(list(p.value = binom.test(signs, n_valid, p = 0.5)$p.value))
}


# --------------------------------
# classify a single data set
# --------------------------------
classify_sample <- function(sample_data, trained_models, preProcStandard, preProcNorm,
                            indiv_threshold = 0.50, final_threshold = 0.60) {
  
  sample_features <- calculate_features(sample_data)
  sample_std  <- predict(preProcStandard, sample_features)
  sample_norm <- predict(preProcNorm, sample_std)
  
  results <- lapply(names(trained_models), function(model_name) {
    model <- trained_models[[model_name]]
    
    p_non_normal <- predict(model, newdata = sample_norm, type = "prob")[, "Non_Normal"]
    pred_class <- ifelse(p_non_normal >= indiv_threshold, "Non_Normal", "Normal")
    
    data.frame(
      Model = model_name,
      Non_Normal_Probability = round(p_non_normal, 4),
      Normal_Probability = round(1 - p_non_normal, 4),
      Prediction = pred_class,
      stringsAsFactors = FALSE
    )
  })
  
  df <- do.call(rbind, results)
  # ---- Majority rule (>= 75% Non_Normal) ----
  prop_non_normal <- mean(df$Prediction == "Non_Normal")
  final_prediction <- ifelse(prop_non_normal >= final_threshold, "Non_Normal", "Normal")
  
  list(
    per_model = df,
    final_prediction = final_prediction,
    prop_non_normal = prop_non_normal
  )
}


# -----------------------------------------------------------------------------
# Function to compute p-values for normality test & downstream tests
# -----------------------------------------------------------------------------
generate_pval <- function(Nsim, n, effect_size,
                          indiv_threshold = 0.50,
                          final_threshold = 0.60,
                          dist = "Normal", 
                          show_progress = TRUE,
                          progress_cb = NULL, ...) {
  
  pval_ds_test1_H0 <- pval_ds_test2_H0 <- numeric(Nsim)
  pval_ds_test1_H1 <- pval_ds_test2_H1 <- numeric(Nsim)
  
  class_H0 <- character(Nsim)
  class_H1 <- character(Nsim)
  
  if (show_progress && is.null(progress_cb)) {
    pb <- txtProgressBar(min = 0, max = Nsim, style = 3)
    on.exit(close(pb), add = TRUE)
  }
  
  for (i in 1:Nsim) {
    
    data_H0 <- do.call(gen_data, get_parameters(n, dist = dist, effect_size = 0))
    data_H1 <- do.call(gen_data, get_parameters(n, dist = dist, effect_size = effect_size))
    
    class_H0[i] <- classify_sample(
      sample_data = data_H0,
      trained_models = models_list,
      preProcStandard = norm_result$preProcStandard,
      preProcNorm = norm_result$preProcNorm,
      indiv_threshold = indiv_threshold,
      final_threshold = final_threshold
    )$final_prediction
    
    class_H1[i] <- classify_sample(
      sample_data = data_H1,
      trained_models = models_list,
      preProcStandard = norm_result$preProcStandard,
      preProcNorm = norm_result$preProcNorm,
      indiv_threshold = indiv_threshold,
      final_threshold = final_threshold
    )$final_prediction
    
    pval_ds_test1_H0[i] <- fn_for_ds_test_1(data_H0)$p.value
    pval_ds_test2_H0[i] <- fn_for_ds_test_2(data_H0)$p.value
    
    pval_ds_test1_H1[i] <- fn_for_ds_test_1(data_H1)$p.value
    pval_ds_test2_H1[i] <- fn_for_ds_test_2(data_H1)$p.value
    
    # update either external callback OR internal bar
    if (!is.null(progress_cb)) progress_cb()
    if (show_progress && is.null(progress_cb)) setTxtProgressBar(pb, i)
  }
  
  list(
    pval_ds_test1_H0 = pval_ds_test1_H0,
    pval_ds_test2_H0 = pval_ds_test2_H0,
    pval_ds_test1_H1 = pval_ds_test1_H1,
    pval_ds_test2_H1 = pval_ds_test2_H1,
    class_H0 = class_H0,
    class_H1 = class_H1
  )
}


# --------------------------------------------------------------------
#  Function to calculate power/Type I errors for each test method
# --------------------------------------------------------------------
perform_analysis <- function(Nsim, n, distributions,
                             effect_size,
                             indiv_threshold = 0.50,
                             final_threshold = 0.60,
                             test_alpha = 0.05,
                             progress_cb = NULL) {
  
  ds_test_results <- list()
  error_ds_test <- list()
  power_ds_test <- list()
  
  for (dist in distributions) {
    
    ds_test_results[[dist]] <- generate_pval(
      Nsim = Nsim, n = n, effect_size = effect_size,
      indiv_threshold = indiv_threshold,
      final_threshold = final_threshold,
      dist = dist,
      show_progress = FALSE,          # disable internal bars
      progress_cb = progress_cb 
    )
    
    # base tests
    error_ds_test[[dist]] <- list(
      error_ds_test1 = mean(ds_test_results[[dist]]$pval_ds_test1_H0 < test_alpha),
      error_ds_test2 = mean(ds_test_results[[dist]]$pval_ds_test2_H0 < test_alpha)
    )
    
    power_ds_test[[dist]] <- list(
      power_ds_test1 = mean(ds_test_results[[dist]]$pval_ds_test1_H1 < test_alpha),
      power_ds_test2 = mean(ds_test_results[[dist]]$pval_ds_test2_H1 < test_alpha)
    )
    
    # adaptive (use H0 classification for type I, H1 classification for power)
    adaptive_pvals_H0 <- ifelse(ds_test_results[[dist]]$class_H0 == "Normal",
                                ds_test_results[[dist]]$pval_ds_test1_H0,
                                ds_test_results[[dist]]$pval_ds_test2_H0)
    
    adaptive_pvals_H1 <- ifelse(ds_test_results[[dist]]$class_H1 == "Normal",
                                ds_test_results[[dist]]$pval_ds_test1_H1,
                                ds_test_results[[dist]]$pval_ds_test2_H1)
    
    error_ds_test[[dist]]$error_adaptive_test <- mean(adaptive_pvals_H0 < test_alpha)
    power_ds_test[[dist]]$power_adaptive_test <- mean(adaptive_pvals_H1 < test_alpha)
  }
  
  list(error_ds_test = error_ds_test,
       power_ds_test = power_ds_test,
       all_pvalues = ds_test_results)
}

# -----------------------------------------------------------------------------#
# Function to compute Normality test trade-off metrics: EPL EPG, & ETIE        #
# -----------------------------------------------------------------------------#
compute_roc_metrics <- function(error_ds_test, power_ds_test, test_alpha) {
  # Non-normal case (index 1)
  power_ds_test_non_normal   <- power_ds_test[[1]]$power_ds_test1
  power_adaptive_ds_test_non_normal <- power_ds_test[[1]]$power_adaptive_test
  error_adaptive_ds_test_non_normal <- error_ds_test[[1]]$error_adaptive_test
  Expected_power_gain <- power_adaptive_ds_test_non_normal - power_ds_test_non_normal
  Expected_type1_error_inflation_non_normal <- error_adaptive_ds_test_non_normal - test_alpha
  
  # Normal case (index 2)
  power_ds_test_normal   <- power_ds_test[[2]]$power_ds_test1
  power_adaptive_ds_test_normal <- power_ds_test[[2]]$power_adaptive_test
  adaptive_error_normal <- error_ds_test[[2]]$error_adaptive_test
  Expected_power_loss      <- power_ds_test_normal - power_adaptive_ds_test_normal
  Expected_type1_error_inflation_normal <- adaptive_error_normal - test_alpha
  
  # Benchmarks (test_1 vs test_2)
  power_gain <- power_ds_test[[1]]$power_ds_test2 - power_ds_test[[1]]$power_ds_test1
  power_loss <- power_ds_test[[2]]$power_ds_test1 - power_ds_test[[2]]$power_ds_test2
  
  list(
    Expected_power_loss      = as.numeric(Expected_power_loss),
    Expected_power_gain      = as.numeric(Expected_power_gain),
    Expected_type1_error_inflation_normal = as.numeric(Expected_type1_error_inflation_normal),
    Expected_type1_error_inflation_non_normal = as.numeric(Expected_type1_error_inflation_non_normal),
    power_gain               = as.numeric(power_gain),
    power_loss               = as.numeric(power_loss)
  )
}

# =============================================================================
# Power and Type I Error Trade-off Analysis with Positive Tolerance Constraint
# 
# Features:
# - Tolerance applies ONLY to positive Type I error inflation
# - Negative inflation (conservative tests) are allowed
# - Excludes points where both distributions show negative inflation
# - Marks selected optimal point with clamped dot
# =============================================================================
plot_power_error_tradeoff <- function(
    distributions = c("Non-normal", "Normal"),
    threshold_grid,   
    metrics,
    outer_title = "Power and Type I error trade-off for pre-testing for normality",
    tol_pos = 0,   # Tolerance for positive Type I error inflation only
    text_size_main = 0.75,
    text_size_labels = 0.75,
    point_size = 0.8,
    line_padding = 0.001
) {
  # Save original graphical parameters and restore on exit
  original_par <- par(no.readonly = TRUE)
  on.exit(par(original_par))
  
  plot_clamped_point <- function(x, y, color = "black", size = 1) {
    plot_region <- par("usr")
    x_padding <- 1e-9 * (plot_region[2] - plot_region[1])
    y_padding <- 1e-9 * (plot_region[4] - plot_region[3])
    clamped_x <- min(max(x, plot_region[1] + x_padding), plot_region[2] - x_padding)
    clamped_y <- min(max(y, plot_region[3] + y_padding), plot_region[4] - y_padding)
    points(clamped_x, clamped_y, pch = 19, cex = size, col = color)
  }
  
  power_loss_normal <- as.numeric(metrics$Expected_power_loss)
  power_gain_non_normal <- as.numeric(metrics$Expected_power_gain)
  type1_inflation_normal <- as.numeric(metrics$Expected_type1_error_inflation_normal)
  type1_inflation_non_normal <- as.numeric(metrics$Expected_type1_error_inflation_non_normal)
  
  valid_indices <- is.finite(threshold_grid) & is.finite(power_loss_normal) &
    is.finite(power_gain_non_normal) & is.finite(type1_inflation_normal) &
    is.finite(type1_inflation_non_normal)
  
  alpha_values <- threshold_grid[valid_indices]  # now: indiv_threshold values
  power_loss <- power_loss_normal[valid_indices]
  power_gain <- power_gain_non_normal[valid_indices]
  type1_infl_normal <- type1_inflation_normal[valid_indices]
  type1_infl_non_normal <- type1_inflation_non_normal[valid_indices]
  
  ord <- order(alpha_values)
  alpha_values <- alpha_values[ord]
  power_loss <- power_loss[ord]
  power_gain <- power_gain[ord]
  type1_infl_normal <- type1_infl_normal[ord]
  type1_infl_non_normal <- type1_infl_non_normal[ord]
  
  
  numerical_precision <- .Machine$double.eps^0.5
  
  positive_infl_normal <- pmax(type1_infl_normal, 0)
  positive_infl_non_normal <- pmax(type1_infl_non_normal, 0)
  worst_case_positive_inflation <- pmax(positive_infl_normal, positive_infl_non_normal)
  
  both_negative <- (type1_infl_normal < 0 & type1_infl_non_normal < 0)
  not_both_negative <- !both_negative
  
  negative_infl_normal <- pmax(-type1_infl_normal, 0)
  negative_infl_non_normal <- pmax(-type1_infl_non_normal, 0)
  worst_case_negative_inflation <- pmax(negative_infl_normal, negative_infl_non_normal)
  
  feasible_indices <- which(
    worst_case_positive_inflation <= tol_pos + numerical_precision &
      not_both_negative
  )
  
  select_optimal_alpha <- function(candidate_indices) {
    positive_inflation_values <- worst_case_positive_inflation[candidate_indices]
    min_positive_inflation <- min(positive_inflation_values, na.rm = TRUE)
    step1_indices <- candidate_indices[
      abs(positive_inflation_values - min_positive_inflation) <= numerical_precision * (1 + min_positive_inflation)
    ]
    
    negative_inflation_values <- worst_case_negative_inflation[step1_indices]
    min_negative_inflation <- min(negative_inflation_values, na.rm = TRUE)
    step2_indices <- step1_indices[
      abs(negative_inflation_values - min_negative_inflation) <= numerical_precision * (1 + min_negative_inflation)
    ]
    
    has_negative_power_loss <- any(power_loss[step2_indices] < -numerical_precision)
    
    if (has_negative_power_loss) {
      gain_values <- power_gain[step2_indices]
      max_gain <- max(gain_values, na.rm = TRUE)
      step3_indices <- step2_indices[
        abs(gain_values - max_gain) <= numerical_precision * (1 + abs(max_gain))
      ]
      
      loss_values <- power_loss[step3_indices]
      min_loss <- min(loss_values, na.rm = TRUE)
      step4_indices <- step3_indices[
        abs(loss_values - min_loss) <= numerical_precision * (1 + abs(min_loss))
      ]
    } else {
      loss_values <- power_loss[step2_indices]
      min_loss <- min(loss_values, na.rm = TRUE)
      step3_indices <- step2_indices[
        abs(loss_values - min_loss) <= numerical_precision * (1 + abs(min_loss))
      ]
      
      gain_values <- power_gain[step3_indices]
      max_gain <- max(gain_values, na.rm = TRUE)
      step4_indices <- step3_indices[
        abs(gain_values - max_gain) <= numerical_precision * (1 + abs(max_gain))
      ]
    }
    
    if (length(step4_indices) > 1) {
      step4_indices[which.min(alpha_values[step4_indices])]
    } else {
      step4_indices
    }
  }
  
  if (length(feasible_indices) > 0) {
    optimal_index <- select_optimal_alpha(feasible_indices)
    selection_note <- sprintf(" (within tolerance & not-both-negative: max +inflation ≤ %.5g)", tol_pos)
  } else {
    not_both_neg_indices <- which(not_both_negative)
    if (length(not_both_neg_indices) > 0) {
      min_positive_infl <- min(worst_case_positive_inflation[not_both_neg_indices], na.rm = TRUE)
      near_optimal_indices <- not_both_neg_indices[
        abs(worst_case_positive_inflation[not_both_neg_indices] - min_positive_infl) <= numerical_precision * (1 + min_positive_infl)
      ]
      optimal_index <- select_optimal_alpha(near_optimal_indices)
      selection_note <- sprintf(" (fallback not-both-negative: min max +inflation = %.5g)", min_positive_infl)
    } else {
      min_negative_infl <- min(worst_case_negative_inflation, na.rm = TRUE)
      near_optimal_indices <- which(
        abs(worst_case_negative_inflation - min_negative_infl) <= numerical_precision * (1 + min_negative_infl)
      )
      optimal_index <- select_optimal_alpha(near_optimal_indices)
      selection_note <- sprintf(" (fallback both-negative: least negative = %.5g)", min_negative_infl)
    }
  }
  
  opt_alpha <- alpha_values[optimal_index]
  optimal_power_loss <- power_loss[optimal_index]
  optimal_power_gain <- power_gain[optimal_index]
  optimal_infl_normal <- type1_infl_normal[optimal_index]
  optimal_infl_non_normal <- type1_infl_non_normal[optimal_index]
  
  par(mfrow = c(2, 2),
      oma  = c(0.8, 0.2, 2.5, 0.2),
      mar  = c(4, 3.5, 1.5, 0.8),
      mgp = c(2.0, 0.6, 0))
  
  create_plot <- function(x, y, color, ylab, main, optimal_index) {
    optimal_y <- y[optimal_index]
    y_range <- range(y, na.rm = TRUE) + c(-1, 1) * line_padding
    
    plot(x, y, type = "l", col = color, lwd = 2, ylim = y_range,
         ylab = ylab,
         xlab = "indiv_threshold",  # was expression(alpha[pre])
         cex.lab = text_size_labels,
         main = main, cex.main = text_size_main)
    
    abline(v = x[optimal_index], lty = 3, col = "grey40")
    plot_clamped_point(x[optimal_index], optimal_y, color, point_size)
  }
  
  create_plot(alpha_values, power_loss, "blue", "Expected Power Loss", sprintf("Power Loss: %s", distributions[2]), optimal_index)
  create_plot(alpha_values, power_gain, "red", "Expected Power Gain", sprintf("Power Gain: %s", distributions[1]), optimal_index)
  create_plot(alpha_values, type1_infl_normal, "orange", "Type I Inflation", sprintf("Type I Inflation: %s", distributions[2]), optimal_index)
  abline(h = 0, lty = 2, col = "grey60")
  create_plot(alpha_values, type1_infl_non_normal, "green4", "Type I Inflation", sprintf("Type I Inflation: %s", distributions[1]), optimal_index)
  abline(h = 0, lty = 2, col = "grey60")
  
  mtext(outer_title, side = 3, outer = TRUE, line = 1.2, cex = 1, font = 2)
  mtext(sprintf("Selected threshold = %.5f | Gain = %.5f | Loss = %.5f | Inflation (Normal, Non-normal) = (%.5f, %.5f) | tol = %.5g",
                opt_alpha, optimal_power_gain, optimal_power_loss, optimal_infl_normal, optimal_infl_non_normal, tol_pos),
        side = 3, outer = TRUE, line = 0.2, cex = 0.6)
  
  invisible(list(
    alpha_star = opt_alpha,  
    power_loss = optimal_power_loss,
    power_gain = optimal_power_gain,
    inflation_normal = optimal_infl_normal,
    inflation_non_normal = optimal_infl_non_normal,
    selection_note = selection_note,
    feasible = length(feasible_indices) > 0,
    worst_case_positive_inflation = max(positive_infl_normal[optimal_index], positive_infl_non_normal[optimal_index]),
    worst_case_negative_inflation = max(negative_infl_normal[optimal_index], negative_infl_non_normal[optimal_index])
  ))
}


# -----------------------------------------------------------------------
#             Function for Power vs Type I error ROC curve
# -----------------------------------------------------------------------
power_vs_error_roc_data <- function(
    N, n,
    distributions = c("Non-normal", "Normal"),
    effect_size,
    threshold_grid = 0.05,      # now acts as indiv_threshold
    sig_levels,
    final_threshold = 0.60     # majority rule threshold 
) {
  
  # empty data frame for results
  roc_data <- data.frame()
  
  pb <- txtProgressBar(min = 0, max = length(distributions) * length(sig_levels), style = 3)
  counter <- 0
  
  for (dist in distributions) {
    cat("Generating ROC data for:", dist, "\n")
    
    # generate p-values + ML classifications
    res <- generate_pval(
      Nsim = N,
      n = n,
      effect_size = effect_size,
      indiv_threshold = threshold_grid,     
      final_threshold = final_threshold,   
      dist = dist
    )
    
    # calculate power and type I error for each significance level
    for (alpha in sig_levels) {
      
      # H0 adaptive decision (classification-based)
      use_ds_test1_H0 <- (res$class_H0 == "Normal")
      adaptive_pvals_H0 <- ifelse(use_ds_test1_H0, res$pval_ds_test1_H0, res$pval_ds_test2_H0)
      
      # Type I error
      type1_error_test1    <- mean(res$pval_ds_test1_H0 <= alpha)
      type1_error_test2    <- mean(res$pval_ds_test2_H0 <= alpha)
      type1_error_adaptive <- mean(adaptive_pvals_H0 <= alpha)
      
      # H1 adaptive decision (classification-based)
      use_ds_test1_H1 <- (res$class_H1 == "Normal")
      adaptive_pvals_H1 <- ifelse(use_ds_test1_H1, res$pval_ds_test1_H1, res$pval_ds_test2_H1)
      
      # Power
      power_test1    <- mean(res$pval_ds_test1_H1 <= alpha)
      power_test2    <- mean(res$pval_ds_test2_H1 <= alpha)
      power_adaptive <- mean(adaptive_pvals_H1 <= alpha)
      
      # organize results
      roc_data <- rbind(
        roc_data,
        data.frame(Distribution = dist, Method = "ds_test 1", Alpha = alpha, Power = power_test1, type1_error = type1_error_test1),
        data.frame(Distribution = dist, Method = "ds_test 2", Alpha = alpha, Power = power_test2, type1_error = type1_error_test2),
        data.frame(Distribution = dist, Method = "adaptive test", Alpha = alpha, Power = power_adaptive, type1_error = type1_error_adaptive)
      )
      
      counter <- counter + 1
      setTxtProgressBar(pb, counter)
    }
  }
  
  close(pb)
  return(roc_data)
}


# ------------------------------------------------------------------------------
#             power vs type I error rate ROC like curve function
# ------------------------------------------------------------------------------
power_vs_error_ROC_curve <- function(roc_data, file = NULL, width = 8, height = 6.8, alpha_star = 0.05, zoom_xlim = 0.10, subtitle = "Power vs Type I error: full view (first row) and zoom view (second row)"){
  
  # Define plot parameters
  dists <- unique(roc_data$Distribution)
  methods <- unique(roc_data$Method)
  colors <- c("red", "blue", "green")
  point_char <- 16  # Solid circle
  
  # Layout
  layout(matrix(1:6, nrow = 3, byrow = TRUE), heights = c(1, 1, 0.30))
  
  # Set graphical parameters
  op <- par(oma = c(0.5, 0.7, 1.5, 0.6), mgp = c(1.6, 0.50, 0), tcl = -0.25, xaxs = "r", yaxs = "r", cex.axis = 0.9, cex.lab = 0.9, cex.main = 0.75)
  on.exit(par(op))
  
  # Plot panel function
  plot_panel <- function(dist, xlim, bottom_mar = 2.0) {
    par(mar = c(bottom_mar, 3.2, 1.0, 0.8))
    
    # Create empty plot
    plot(NA, xlim = xlim, ylim = c(0, 1), xlab = "P(Type I Error)", ylab = "Power", main = paste0("ROC-like Curve (", dist, ")"), cex.lab = 0.65)
    
    # Add lines for each method
    for (i in seq_along(methods)) {
      data_sub <- subset(roc_data, Distribution == dist & Method == methods[i])
      if (nrow(data_sub) == 0) next
      
      lines(data_sub$type1_error, data_sub$Power, col = colors[i], lwd = 2)
      
      # Mark alpha_star point
      idx <- which.min(abs(data_sub$Alpha - alpha_star))
      if (length(idx) == 1 && is.finite(idx)) {
        points(data_sub$type1_error[idx], data_sub$Power[idx], pch = point_char, col = colors[i])
      }
    }
    
    # Reference lines
    abline(a = 0, b = 1, col = "grey80")
    abline(v = alpha_star, col = "red", lty = 2)  
  }
  
  # Create plots
  plot_panel(dists[1], c(0, 1), 2.0)           # Full view - dist 1
  plot_panel(dists[2], c(0, 1), 2.0)           # Full view - dist 2
  plot_panel(dists[1], c(0, zoom_xlim), 2.5)   # Zoom view - dist 1
  plot_panel(dists[2], c(0, zoom_xlim), 2.5)   # Zoom view - dist 2
  
  # Legend function
  add_legend <- function(dist) {
    par(mar = c(0, 0.3, 0.6, 0.3))
    plot.new()
    
    legend_labels <- sapply(seq_along(methods), function(i) {
      data_sub <- subset(roc_data, Distribution == dist & Method == methods[i])
      if (nrow(data_sub) == 0) return(methods[i])
      # add alpha_star point
      idx <- which.min(abs(data_sub$Alpha - alpha_star))
      x <- data_sub$type1_error[idx]
      y <- data_sub$Power[idx]
      sprintf("%s (%.3f, %.3f)", methods[i], x, y)
    })
    
    legend("top", legend = legend_labels, title = "Methods", col = colors, lwd = 2, pch = point_char, cex = 0.75, bty = "n")
  }
  
  # Add legends
  add_legend(dists[1])
  add_legend(dists[2])
  
  # Main title
  mtext(subtitle, side = 3, outer = TRUE, line = 0.2, cex = 0.65)
}
# --------------------------------------------------------------------
# =============================================================================
# ================ Main Downstream Test Function (ML adaptive) ================
# --------- Decide test type(test 1 or test 2) using ML classifier ------------
ds_test_function <- function(
    gen_data = gen_data,
    get_parameters = get_parameters,
    fn_to_get_norm_obj = fn_to_get_norm_obj,
    fn_for_ds_test_1 = fn_for_ds_test_1,
    fn_for_ds_test_2 = fn_for_ds_test_2,
    paras = NULL,
    alpha = 0.05,
    
    # --- ML adaptive inputs ---
    trained_models,
    preProcStandard,
    preProcNorm,
    indiv_threshold = 0.50,
    final_threshold = 0.60,
    ...
) {
  # generate data
  data <- if (!is.null(paras)) do.call(gen_data, paras) else gen_data()
  
  # keep this (you said we still need fn_to_get_norm_obj)
  normality_obj <- fn_to_get_norm_obj(data)
  
  # ML classification replaces normality test
  class_label <- classify_sample(
    sample_data     = normality_obj,
    trained_models  = trained_models,
    preProcStandard = preProcStandard,
    preProcNorm     = preProcNorm,
    indiv_threshold = indiv_threshold,
    final_threshold = final_threshold
  )$final_prediction
  
  # choose downstream test
  if (identical(class_label, "Normal")) {
    fn_for_ds_test_1(data)$p.value
  } else {
    fn_for_ds_test_2(data, ...)$p.value
  }
}



# perform test 1, test 2 & adaptive test
perform_ds_func <- function(
    sample_sizes = sample_sizes,
    distributions = c("Non-normal", "Normal"),
    N = 1e3,
    alpha = 0.05,
    gen_data = gen_data,
    get_parameters = get_parameters,
    fn_to_get_norm_obj = fn_to_get_norm_obj,
    fn_for_ds_test_1 = fn_for_ds_test_1,
    fn_for_ds_test_2 = fn_for_ds_test_2,
    ds_test_methods = ds_test_methods,
    effect_size = effect_size,
    
    # --- ML adaptive args (NEW) ---
    trained_models,
    preProcStandard,
    preProcNorm,
    indiv_threshold = 0.50,
    final_threshold = 0.60,
    ...
) {
  # pick correct ds test method
  ds_test_methods <- match.arg(ds_test_methods, several.ok = TRUE)
  # storage for all ds tests results
  results_by_dist <- list()
  
  # perform ds tests for each distribution
  for (dist in distributions) {
    cat("\n=== Running simulation for distribution:", dist, "===\n")
    # storage for power and type I error
    results_power <- list()
    results_type1 <- list()
    # track runtime
    timing <- list()
    # store all p-values
    pval_storage <- list(H0 = list(), H1 = list())
    
    # for each ds test method
    for (method in ds_test_methods) {
      cat("Running method:", method, "\n")
      # storage for each sample size
      pow_vec <- typ1_vec <- numeric(length(sample_sizes))
      names(pow_vec) <- names(typ1_vec) <- paste0("n=", sample_sizes)
      # store all p-values
      pval_storage$H0[[method]] <- vector("list", length(sample_sizes))
      pval_storage$H1[[method]] <- vector("list", length(sample_sizes))
      
      # set progress bar
      pb <- txtProgressBar(min = 0, max = length(sample_sizes), style = 3)
      # start time
      t0 <- Sys.time()
      
      for (i in seq_along(sample_sizes)) {
        n <- sample_sizes[i]
        
        rej_H0 <- rej_H1 <- 0L
        p_H0 <- p_H1 <- numeric(N)
        
        for (sim in seq_len(N)) {
          # get parameters
          paras_H0 <- get_parameters(n, effect_size = 0.0, dist = dist)
          paras_H1 <- get_parameters(n, effect_size = effect_size, dist = dist)
          # generate data
          dat_H0 <- do.call(gen_data, paras_H0)
          dat_H1 <- do.call(gen_data, paras_H1)
          
          if (method == "test_1") {
            p0 <- fn_for_ds_test_1(dat_H0)$p.value
            p1 <- fn_for_ds_test_1(dat_H1)$p.value
            
          } else if (method == "test_2") {
            p0 <- fn_for_ds_test_2(dat_H0)$p.value
            p1 <- fn_for_ds_test_2(dat_H1)$p.value
            
          } else {
            # adaptive test (ML classifier decides)
            p0 <- ds_test_function(
              gen_data = gen_data,
              get_parameters = get_parameters,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_ds_test_1 = fn_for_ds_test_1,
              fn_for_ds_test_2 = fn_for_ds_test_2,
              paras = paras_H0,
              alpha = alpha,
              
              trained_models  = trained_models,
              preProcStandard = preProcStandard,
              preProcNorm     = preProcNorm,
              indiv_threshold = indiv_threshold,
              final_threshold = final_threshold
            )
            
            p1 <- ds_test_function(
              gen_data = gen_data,
              get_parameters = get_parameters,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_ds_test_1 = fn_for_ds_test_1,
              fn_for_ds_test_2 = fn_for_ds_test_2,
              paras = paras_H1,
              alpha = alpha,
              
              trained_models  = trained_models,
              preProcStandard = preProcStandard,
              preProcNorm     = preProcNorm,
              indiv_threshold = indiv_threshold,
              final_threshold = final_threshold
            )
          }
          
          # store p-values
          p_H0[sim] <- p0
          p_H1[sim] <- p1
          if (p0 < alpha) rej_H0 <- rej_H0 + 1
          if (p1 < alpha) rej_H1 <- rej_H1 + 1
        }
        
        # calculate power & type I error rates
        typ1_vec[i] <- rej_H0 / N
        pow_vec[i]  <- rej_H1 / N
        pval_storage$H0[[method]][[i]] <- p_H0
        pval_storage$H1[[method]][[i]] <- p_H1
        setTxtProgressBar(pb, i)
      }
      close(pb)
      
      timing[[method]] <- as.numeric(difftime(Sys.time(), t0, units = "secs"))
      results_type1[[method]] <- typ1_vec
      results_power[[method]] <- pow_vec
      cat("Method", method, "completed in", round(timing[[method]], 2), "seconds\n")
    }
    
    # collect all results by distribution
    results_by_dist[[dist]] <- list(
      power = results_power,
      type1 = results_type1,
      pvalues = pval_storage,
      timing = timing
    )
  }
  
  return(results_by_dist)
}

# -----------------------------------------------------------------------------
#                     power & type I error plot function 
# -----------------------------------------------------------------------------
plot_power_type1 <- function(combined_power, combined_type1, optimal_alphas, methods = c("test_1", "test_2", "adaptive"), distributions = c("Non-normal", "Normal"), sizes, alpha, filename = NULL, width = 7, height = 6) {
  
  # Setup
  colors <- c("red", "blue", "green")
  shapes <- c(19, 17, 15)
  
  # Open PDF if file name provided
  if (!is.null(filename)) {
    pdf(filename, width = width, height = height)
    on.exit(dev.off())
  }
  
  # Plot layout: 2 plots, 2 plots, 1 legend
  layout(matrix(c(1,2,3,4,5,5), nrow = 3, byrow = TRUE), heights = c(1,1,0.25))
  op <- par(mar = c(2.6, 3.2, 1.2, 0.8), oma = c(0, 0,1.6,0), mgp = c(1.6, 0.45, 0), tcl = -0.2, cex.axis = 0.75, cex.lab = 0.75)
  on.exit(par(op), add = TRUE)
  
  # Plot function
  create_plot <- function(data, ylab, main, hline = NULL) {
    xlim <- range(sizes) + c(-1, 1) * 0.02 * diff(range(sizes))
    
    plot(NA, xlim = xlim, ylim = c(0, 1.01), xlab = "Sample Size", ylab = ylab, main = main)
    
    for (i in seq_along(methods)) {
      lines(data$n, data[[methods[i]]], type = "b", col = colors[i], pch = shapes[i], lwd = 2)
    }
    
    if (!is.null(hline)) abline(h = hline, col = "gray50", lty = 2)
  }
  # Create plots
  for (dist in distributions) {
    data_power <- subset(combined_power, Distribution == dist)
    create_plot(data_power, "Power", paste("Power -", dist), 0.75)
  }
  
  for (dist in distributions) {
    data_type1 <- subset(combined_type1, Distribution == dist)
    ymax <- max(data_type1[methods], na.rm = TRUE)
    plot(NA, xlim = range(sizes), ylim = c(0, ymax * 1.02), xlab = "Sample Size", ylab = "P(Type I Error)", main = paste("P(Type I Error) -", dist))
    
    for (i in seq_along(methods)) {
      lines(data_type1$n, data_type1[[methods[i]]], type = "b", col = colors[i], pch = shapes[i], lwd = 2)
    }
    abline(h = alpha, col = "gray50", lty = 3)
  }
  
  # Legend
  par(mar = c(0.2, 0.2, 0.2, 0.2))
  plot.new()
  legend("center", legend = methods, col = colors, pch = shapes, lwd = 2, horiz = TRUE, bty = "n", cex = 0.75, title = "Methods")
  
  # Title with optimal alphas
  mtext(sprintf("Test Methods Comparison (optimal individual threshold: %s)", paste(round(optimal_alphas, 3), collapse = ", ")), side = 3, outer = TRUE, line = 0.35, cex = 0.75)
}

# ----------------------------------------------------------------------------
# calculate power for each effect size
# ----------------------------------------------------------------------------
perform_ds_power_by_effect <- function(
    fixed_n = 10, effect_sizes = NULL,
    distributions = c("Non-normal", "Normal"),
    N = N, alpha = alpha,
    gen_data = gen_data,
    get_parameters = get_parameters,
    fn_to_get_norm_obj = fn_to_get_norm_obj,
    fn_for_ds_test_1 = fn_for_ds_test_1,
    fn_for_ds_test_2 = fn_for_ds_test_2,
    ds_test_methods = ds_test_methods,
    
    # --- ML adaptive args ---
    trained_models,
    preProcStandard,
    preProcNorm,
    indiv_threshold = 0.50,
    final_threshold = 0.60,
    ...
) {
  
  results_by_dist <- list()
  
  for (dist in distributions) {
    cat("Power simulation:", dist, "n =", fixed_n, "\n")
    power_results <- list()
    timing <- list()
    pvalues <- list(H1 = list())
    
    for (method in ds_test_methods) {
      cat("Method:", method, "\n")
      
      power_vec <- numeric(length(effect_sizes))
      names(power_vec) <- paste0("d=", effect_sizes)
      pvalues$H1[[method]] <- vector("list", length(effect_sizes))
      
      pb <- txtProgressBar(0, length(effect_sizes), style = 3)
      start_time <- Sys.time()
      
      for (i in seq_along(effect_sizes)) {
        effect <- effect_sizes[i]
        rejections <- 0
        p_vals <- numeric(N)
        
        for (sim in 1:N) {
          paras <- get_parameters(fixed_n, effect_size = effect, dist = dist)
          data <- do.call(gen_data, paras)
          
          if (method == "test_1") {
            p_val <- fn_for_ds_test_1(data)$p.value
            
          } else if (method == "test_2") {
            p_val <- fn_for_ds_test_2(data)$p.value
            
          } else {
            # adaptive test: ML classification replaces normality p-values 
            p_val <- ds_test_function(
              gen_data = gen_data,
              get_parameters = get_parameters,
              fn_to_get_norm_obj = fn_to_get_norm_obj,
              fn_for_ds_test_1 = fn_for_ds_test_1,
              fn_for_ds_test_2 = fn_for_ds_test_2,
              paras = paras,
              alpha = alpha,
              
              trained_models  = trained_models,
              preProcStandard = preProcStandard,
              preProcNorm     = preProcNorm,
              indiv_threshold = indiv_threshold,
              final_threshold = final_threshold
            )
            
            
          }
          
          p_vals[sim] <- p_val
          if (p_val < alpha) rejections <- rejections + 1
        }
        
        power_vec[i] <- rejections / N
        pvalues$H1[[method]][[i]] <- p_vals
        setTxtProgressBar(pb, i)
      }
      
      close(pb)
      timing[[method]] <- as.numeric(Sys.time() - start_time)
      power_results[[method]] <- power_vec
      cat("Completed in", round(timing[[method]], 2), "s\n")
    }
    
    results_by_dist[[dist]] <- list(
      power = power_results,
      pvalues = pvalues,
      timing = timing
    )
  }
  
  results_by_dist
}

# -----------------------------------------------------------------------
# plot power against effect sizes
# ----------------------------------------------------------------------
plot_power_by_effect_size <- function(power_results, 
                                      distributions = c("Non-normal", "Normal"),
                                      ds_test_methods = c("test_1", "test_2", "adaptive"),
                                      effect_sizes,
                                      indiv_threshold = 0.05) {
  
  colors <- c("red", "blue", "green")
  shapes <- c(19, 17, 15)
  
  # Save and restore graphics settings
  op <- par(no.readonly = TRUE)
  on.exit(par(op))
  
  # Layout: 2 plots on top, legend row at bottom
  layout(
    matrix(c(1, 2,
             3, 3), nrow = 2, byrow = TRUE),
    heights = c(0.8, 0.2)
  )
  
  # Margins for the two main plots
  par(mar = c(4, 4, 2, 1), oma = c(0, 0, 2, 0))
  
  for (dist in distributions) {
    dist_data <- subset(power_results, Distribution == dist)
    
    plot(NA,
         xlim = range(effect_sizes),
         ylim = c(0, 1),
         xlab = "Effect Size",
         ylab = "Power",
         main = dist,
         cex.main = 0.9)
    
    for (i in seq_along(ds_test_methods)) {
      method <- ds_test_methods[i]
      if (method %in% names(dist_data)) {
        y <- dist_data[[method]]
        lines(effect_sizes, y, col = colors[i], lwd = 2)
        points(effect_sizes, y, pch = shapes[i], col = colors[i])
      }
    }
    
    abline(v = 0.5, col = "gray50", lty = 2)
  }
  
  # Overall title
  mtext(sprintf("Power vs Effect Size | optimal individual threshold = %.3f", indiv_threshold),
        outer = TRUE, side = 3, line = 0, cex = 0.9)
  
  # Legend panel (bottom row)
  par(mar = c(0, 0, 0, 0))
  plot.new()
  legend("center",
         legend = ds_test_methods,
         title = "Method",
         col = colors,
         lwd = 2,
         pch = shapes,
         bty = "n",
         horiz = TRUE,
         cex = 0.9)
}


# =============================================================================
#                             RUN SIMULATIONS                                 #
# =============================================================================
# ------------------------------------------------------------
# Wrapper: compute tradeoff metrics over a grid of indiv_threshold
# ------------------------------------------------------------
compute_tradeoff_metrics_over_grid <- function(
    threshold_grid,
    Nsim, n,
    distributions,
    effect_size,
    final_threshold = 0.60,
    test_alpha = 0.05
) {
  out <- list(
    Expected_power_loss = rep(NA_real_, length(threshold_grid)),
    Expected_power_gain = rep(NA_real_, length(threshold_grid)),
    Expected_type1_error_inflation_normal = rep(NA_real_, length(threshold_grid)),
    Expected_type1_error_inflation_non_normal = rep(NA_real_, length(threshold_grid)),
    power_gain = rep(NA_real_, length(threshold_grid)),
    power_loss = rep(NA_real_, length(threshold_grid))
  )
  
  # set progress bar
  total_steps <- length(threshold_grid) * length(distributions) * Nsim
  pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
  on.exit(close(pb), add = TRUE)
  
  step <- 0L
  progress_cb <- function() {
    step <<- step + 1L # super-assignment
    if (step %% 10L == 0L || step == total_steps) setTxtProgressBar(pb, step)
  }
  
  for (k in seq_along(threshold_grid)) {
    thr <- threshold_grid[k]
    
    analysis_k <- perform_analysis(
      Nsim = Nsim,
      n = n,
      distributions = distributions,
      effect_size = effect_size,
      indiv_threshold = thr,
      final_threshold = final_threshold,
      test_alpha = test_alpha,
      progress_cb = progress_cb
    )
    
    m_k <- compute_roc_metrics(
      error_ds_test = analysis_k$error_ds_test,
      power_ds_test = analysis_k$power_ds_test,
      test_alpha = test_alpha
    )
    
    out$Expected_power_loss[k] <- m_k$Expected_power_loss
    out$Expected_power_gain[k] <- m_k$Expected_power_gain
    out$Expected_type1_error_inflation_normal[k] <- m_k$Expected_type1_error_inflation_normal
    out$Expected_type1_error_inflation_non_normal[k] <- m_k$Expected_type1_error_inflation_non_normal
    out$power_gain[k] <- m_k$power_gain
    out$power_loss[k] <- m_k$power_loss
    
  }
  
  out
}

# run the entire simulation
run_simulation <- function(n, N, Nsim, test_type, distributions, tol_pos ) {
  
  # create results folder
  if (!dir.exists("results")) dir.create("results")
  
  ## Keep a small results 
  results <- list(
    params  = list(N = N, Nsim = Nsim, test_type = test_type, distributions = distributions ),
    objects = list(),
    files   = list()
  )
  # define run control parameters
  sample_size       <- c(10, 20, 30, 40, 50)
  test_alpha        <- 0.05
  alpha             <- 0.05
  select_norm_test  <- "SW"
  effect_size       <- 0.5
  norm_alpha_pre    <- seq(from = 0.0, to = 1, by = 0.05)
  threshold_grid     <- seq(from = 0.1, to = 1, by = 0.1)
  effect_sizes      <- c(0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0)
  sig_levels        <- seq(from = 0.005, to = 1, by = 0.005)
  ds_test_methods   <- c("test_1", "test_2", "adaptive")
  norm_test         <- c("SW", "SF", "KS", "JB","SKEW", "DAP", "AD", "CVM")
  
  
  # ----------------------------------------------------
  # Trade-off analysis (compute curves over threshold_grid)
  # ----------------------------------------------------
  metrics_grid <- compute_tradeoff_metrics_over_grid(
    threshold_grid = threshold_grid,
    Nsim = Nsim,
    n = 10,
    distributions = distributions,
    effect_size = effect_size,
    final_threshold = 0.60,
    test_alpha = test_alpha
  )
  
  results$objects$metrics_grid <- metrics_grid
  
  pdf(paste0("results/", test_type, "_test_tradeoff.pdf"), width = 6, height = 5)
  tradeoff_result <- plot_power_error_tradeoff(
    distributions = distributions,
    threshold_grid = threshold_grid,
    metrics = metrics_grid,
    tol_pos = tol_pos
  )
  dev.off()
  
  
  # Save the trade off result and extract alpha_star
  results$objects$tradeoff_result <- tradeoff_result
  alpha_star <- tradeoff_result$alpha_star
  results$objects$alpha_star <- alpha_star
  
  cat("=== Selected alpha_star from tradeoff analysis:", round(alpha_star, 4), "===\n")
  
  ## Generate data for Power vs Type I error like plot 
  roc_data <- power_vs_error_roc_data(
    N = N, 
    n = 10,
    distributions = distributions, 
    effect_size = effect_size,
    threshold_grid = alpha_star,  
    sig_levels = sig_levels
  )
  
  # store results
  results$objects$roc_data <- roc_data
  
  #  Power vs Type I error ROC-like plots 
  pdf(paste0("results/", test_type, "_test_power_error_roc.pdf"), width = 6, height = 5)
  power_vs_error_ROC_curve(roc_data = roc_data)
  dev.off()
  
  # -------------------------------------------------------
  # ------- Main simulation: power/Type I error analysis --
  # -------------------------------------------------------
  sim_output <- perform_ds_func(
    sample_sizes       = sample_size,
    distributions      = distributions,
    N                  = N,
    alpha              = alpha,
    gen_data           = gen_data,
    get_parameters     = get_parameters,
    fn_to_get_norm_obj = fn_to_get_norm_obj,
    fn_for_norm_test   = normality_test,   
    fn_for_ds_test_1   = fn_for_ds_test_1,
    fn_for_ds_test_2   = fn_for_ds_test_2,
    norm_test_method   = select_norm_test, 
    ds_test_methods    = ds_test_methods,
    effect_size        = effect_size,
    indiv_threshold          = alpha_star,       
    
    # --- REQUIRED for ML adaptive path ---
    trained_models     = models_list,
    preProcStandard    = norm_result$preProcStandard,
    preProcNorm        = norm_result$preProcNorm,
    final_threshold    = 0.60
  )
  
  
  # store results
  results$objects$sim_output <- sim_output
  
  ## Per-distribution AUCs (printed and saved) 
  all_auc_tables <- list()
  for (dist in names(sim_output)) {
    cat("\n\n=== Results for", dist, "distribution ===\n")
    dist_results  <- sim_output[[dist]]
    power_results <- dist_results$power
    type1_results <- dist_results$type1
    
    # power
    power_df <- list(
      test_1    = power_results$test_1,
      test_2 = power_results$test_2,
      adaptive      = power_results$adaptive
    )
    
    # type I error
    TypeI_error_df <- list(
      test_1    = type1_results$test_1,
      test_2 = type1_results$test_2,
      adaptive      = type1_results$adaptive
    )
    
    # calculate aucs
    auc_power <- sapply(power_df, function(y) compute_area(sample_size, y))
    auc_type1 <- sapply(TypeI_error_df, function(y) compute_area(sample_size, y))
    
    # store in a dataframe
    auc_table <- data.frame(
      Method    = names(auc_power),
      AUC_Power = unname(auc_power),
      AUC_TypeI = unname(auc_type1),
      row.names = NULL
    )
    all_auc_tables[[dist]] <- auc_table
    print(auc_table, digits = 4)
  }
  
  # store results
  results$objects$auc_tables <- all_auc_tables
  
  ## Combine across distributions and plot 
  all_power_results <- list()
  all_type1_results <- list()
  
  for (dist in names(sim_output)) {
    dist_results <- sim_output[[dist]]
    
    power_df <- data.frame(
      n = sample_size,
      test_1    = dist_results$power$test_1,
      test_2 = dist_results$power$test_2,
      adaptive      = dist_results$power$adaptive,
      Distribution  = dist
    )
    type1_df <- data.frame(
      n = sample_size,
      test_1    = dist_results$type1$test_1,
      test_2 = dist_results$type1$test_2,
      adaptive      = dist_results$type1$adaptive,
      Distribution  = dist
    )
    
    all_power_results[[dist]] <- power_df
    all_type1_results[[dist]] <- type1_df
  }
  
  combined_power <- do.call(rbind, all_power_results)
  combined_type1 <- do.call(rbind, all_type1_results)
  results$objects$combined_power_n <- combined_power
  results$objects$combined_type1_n <- combined_type1
  
  # 
  pdf(paste0("results/", test_type, "_power_and_error_test.pdf"), width = 6, height = 5)
  plot_power_type1(
    combined_power  = combined_power,
    combined_type1  = combined_type1,
    optimal_alphas  = alpha_star,        
    methods         = ds_test_methods,   
    distributions   = distributions,
    sizes           = sample_size,
    alpha           = alpha
  )
  dev.off()
  # -----------------------------------------------------------
  ## Effect-size power curves 
  # ----------------------------------------------------------
  sample_size_es <- 10
  sim_power_list <- list()
  
  for (n in sample_size_es) {
    sim_power_list[[as.character(n)]] <- perform_ds_power_by_effect(
      fixed_n            = n,
      effect_sizes       = effect_sizes,
      distributions      = distributions,
      N                  = N,
      alpha              = alpha,
      gen_data           = gen_data,
      get_parameters     = get_parameters,
      fn_to_get_norm_obj = fn_to_get_norm_obj,
      fn_for_norm_test   = normality_test,   
      fn_for_ds_test_1   = fn_for_ds_test_1,
      fn_for_ds_test_2   = fn_for_ds_test_2,
      norm_test_method   = select_norm_test, 
      ds_test_methods    = ds_test_methods,
      indiv_threshold          = alpha_star,       
      
      # --- REQUIRED for ML adaptive path ---
      trained_models     = models_list,
      preProcStandard    = norm_result$preProcStandard,
      preProcNorm        = norm_result$preProcNorm,
      final_threshold    = 0.60
    )
    
  }
  
  # store results
  results$objects$sim_power_list <- sim_power_list
  
  all_power_df <- list()
  for (n in names(sim_power_list)) {
    for (dist in names(sim_power_list[[n]])) {
      pr <- sim_power_list[[n]][[dist]]$power
      all_power_df[[paste(dist, n, sep = "_")]] <- data.frame(
        d    = effect_sizes,
        test_1    = pr$test_1,
        test_2 = pr$test_2,
        adaptive      = pr$adaptive,
        Distribution  = dist,
        n             = as.integer(n)
      )
    }
  }
  
  # store results
  combined_power_effect <- do.call(rbind, all_power_df)
  results$objects$combined_power_effect <- combined_power_effect
  
  # plot power vs effect sizes
  pdf(paste0("results/", test_type, "_power_by_effect.pdf"), 6, 4)
  plot_power_by_effect_size(
    power_results   = combined_power_effect,
    distributions   = distributions,
    ds_test_methods = ds_test_methods,
    effect_sizes    = effect_sizes,
    indiv_threshold       = alpha_star
  )
  dev.off()
  
  ## ---------- Save compact results + full workspace ----------
  save(results, file = paste0("results/", test_type, "_results.RData"))
  save.image(paste0("results/", test_type, "_workspace.RData"))
  
  # Final summary
  cat("\n=== SIMULATION COMPLETE ===\n")
  cat("=== Selected indiv_threshold_star from tradeoff analysis:", round(alpha_star, 4), "===\n")
  cat("Tradeoff feasibility:", if(tradeoff_result$feasible) "FEASIBLE" else "FALLBACK", "\n")
  cat("Worst-case positive inflation:", round(tradeoff_result$worst_case_positive_inflation, 6), "\n")
  cat("Inflation values (Normal, Non-normal):",
      round(tradeoff_result$inflation_normal, 6), ",",
      round(tradeoff_result$inflation_non_normal, 6), "\n")
  
  invisible(results)
}
## Example call
run_simulation(
  N = 1e4,
  Nsim = 1e5,
  test_type = "one_sample_t_vs_sign_ml",
  distributions = c("exponential", "normal"),
  tol_pos = 1e-2
)

