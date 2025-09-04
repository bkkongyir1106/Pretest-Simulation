gen_data <- function(
    n = n,           
    beta0 = 0,        
    beta1 = effect_size, 
    x_dist = "Exponential",
    error_sd = 1,     
    dist = "Normal"  
) {
  # predictor
  x <- generate_data(n, dist = x_dist)
  # error
  error <- error_sd * generate_data(n, dist)
  # independent variable
  y <- beta0 + beta1 * x + error
  return(data.frame(x = x, y = y))
}

# function for parameters
get_parameters <- function(n, effect_size = NULL, ...) {
  defaults <- list(
    n = n,
    beta0 = 2,
    beta1 = 0,  # Default no effect
    x_dist = "Exponential",
    error_sd = 1,
    dist = "lognormal"
  )
  
  args <- list(...)
  
  # If effect_size provided, use it for beta1
  if (!is.null(effect_size)) {
    args$beta1 <- effect_size
  }
  
  modifyList(defaults, args)
}

# function for norm object
fn_to_get_norm_obj <- function(data){
  model <- lm(y ~ x , data = data)
  return(residuals(model))
}

# ds function 1
fn_for_ds_test_1<- function(data) {
  model <- lm(y ~ x, data = data)
  tidy_model <- broom::tidy(model)
  p_value <- tidy_model$p.value[tidy_model$term == "x"]
  return(list(p.value = p_value))
}

fn_for_ds_test_2 <- function(data, n_boot = 1000) {
  # Data validation
  if (!all(c("x", "y") %in% names(data))) {
    stop("Data must contain 'x' and 'y' columns")
  }
  
  if (nrow(data) < 2) {
    stop("Insufficient data for bootstrap")
  }
  
  # Observed slope from original data
  observed_slope <- coef(lm(y ~ x, data = data))[2]
  
  # Paired bootstrap procedure (Case 1 from lecture)
  boot_stat <- function(data, idx) {
    # Resample pairs (y_i, x_i) with replacement
    boot_data <- data[idx, ]
    model_boot <- lm(y ~ x, data = boot_data)
    coef(model_boot)[2]  # Return bootstrap slope estimate
  }
  
  # Perform bootstrap
  boot_results <- boot::boot(data, statistic = boot_stat, R = n_boot)
  
  # Calculate p-value using the guideline: compare observed statistic 
  # with bootstrap distribution under H0 (using θ* - θ̂)
  # For H0: β₁ = 0, we compare |β₁_obs - 0| with |β₁* - β₁_obs|
  bootstrap_diffs <- abs(boot_results$t - observed_slope)
  observed_diff <- abs(observed_slope - 0)  # Since H0: β₁ = 0
  
  p_value <- mean(bootstrap_diffs >= observed_diff)
  
  return(list(
    p.value = p_value,
    observed.slope = observed_slope,
    boot.results = boot_results
  ))
}



# # ds function 2
# fn_for_ds_test_2 <- function(data, n_boot = 1000) {
#   # if (!require(boot)) install.packages("boot")
#   # library(boot)
#   mu_hat   <- mean(data$y)
#   resid0   <- data$y - mu_hat
#   df_boot  <- data.frame(x = data$x, resid0 = resid0)
#   
#   # Observed slope
#   observed_slope <- coef(lm(y ~ x, data = data))[2]
#   
#   boot_stat <- function(df, idx) {
#     y_star <- mu_hat + df$resid0[idx]    
#     model  <- lm(y_star ~ df$x)         
#     coef(model)[2]
#   }
#   
#   boot_results <- boot::boot(df_boot, statistic = boot_stat, R = n_boot)
#   
#   p_value <- mean(abs(boot_results$t) >= abs(observed_slope))
#   
#   return(list(p.value = p_value))
# }
# 
