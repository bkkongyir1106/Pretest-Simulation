setwd("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/Summer 2025/user_framework")
source("~/Desktop/OSU/Research/Pretest-Simulation/Simulation/Summer 2025/user_framework/my_functions_v2.R")

calculate_power <- function(
    sample_sizes    = c(10, 20, 30, 40, 50),
    n_sim           = 1e3,
    alpha           = 0.05,
    n_boot          = NULL,
    gen_data,
    get_parameters,
    fn_to_get_norm_obj,
    fn_for_norm_test,
    fn_for_ds_test_1,
    fn_for_ds_test_2,
    test_method    = "SW",
    mode           = c("adaptive", "parametric", "nonparametric"),
    update_progress = NULL,   # <-- add this
    ...
) {
  mode <- match.arg(mode)
  power_results <- numeric(length(sample_sizes))
  names(power_results) <- paste0("n=", sample_sizes)
  
  total_steps <- length(sample_sizes) * n_sim
  step_counter <- 0
  
  for (i in seq_along(sample_sizes)) {
    n <- sample_sizes[i]
    rejections <- 0
    
    for (sim in seq_len(n_sim)) {
      step_counter <- step_counter + 1
      if (!is.null(update_progress)) update_progress(step_counter / total_steps)
      
      paras <- get_parameters(n)
      data <- do.call(gen_data, paras)
      
      p_value <- tryCatch({
        if (mode == "adaptive") {
          user_interface(
            gen_data           = gen_data,
            fn_to_get_norm_obj = fn_to_get_norm_obj,
            fn_for_norm_test   = fn_for_norm_test,
            fn_for_ds_test_1   = fn_for_ds_test_1,
            fn_for_ds_test_2   = fn_for_ds_test_2,
            paras              = paras,
            alpha              = alpha,
            test_method        = test_method,
            ...
          )
        } else if (mode == "parametric") {
          fn_for_ds_test_1(data)$p.value
        } else {
          fn_for_ds_test_2(data, ...)$p.value
        }
      }, error = function(e) NA)
      
      if (!is.na(p_value) && p_value < alpha) {
        rejections <- rejections + 1
      }
    }
    
    power_results[i] <- rejections / n_sim
    cat("Completed n =", n, "| mode =", mode, "| Power =", power_results[i], "\n")
  }
  
  return(power_results)
}


# ------------ Power Curve Plotting -----------
plot_power_curve <- function(power_results, title = "Power Analysis") {
  sample_sizes <- as.numeric(sub("n=", "", names(power_results)))
  plot(sample_sizes, power_results, type = "b",
       xlab = "Sample Size", ylab = "Power",
       main = title, ylim = c(0, 1), pch = 19)
  abline(h = 0.8, lty = 2)
  grid()
}

# ----------------------------------------------------
#                     ---- Example ----
# ----------------------------------------------------
# sample_sizes <- c(10, 20, 30, 40, 50)
# n_sim <- 1e4
# 
# # Parametric only 
# power_anova_param    <- calculate_power(
#   sample_sizes   = sample_sizes,
#   n_sim          = n_sim,
#   get_parameters = get_anova_params,
#   gen_data       = generate_anova_data,
#   fn_to_get_norm_obj = anova_residuals,
#   fn_for_ds_test_1   = anova_main_test,
#   fn_for_ds_test_2   = kruskal_wallis_test,
#   mode           = "parametric"
# )
# 
# # Nonparametric only 
# power_anova_nonpar    <- calculate_power(
#   sample_sizes   = sample_sizes,
#   n_sim          = n_sim,
#   get_parameters = get_anova_params,
#   gen_data       = generate_anova_data,
#   fn_to_get_norm_obj = anova_residuals,
#   fn_for_ds_test_1   = anova_main_test,
#   fn_for_ds_test_2   = kruskal_wallis_test,
#   mode           = "nonparametric"
#   #n_boot = 1000
# )
# 
# # Adaptive (normality-based)
# power_anova_adapt    <- calculate_power(
#   sample_sizes   = sample_sizes,
#   n_sim          = n_sim,
#   get_parameters = get_anova_params,
#   gen_data       = generate_anova_data,
#   fn_to_get_norm_obj = anova_residuals,
#   fn_for_ds_test_1   = anova_main_test,
#   fn_for_ds_test_2   = kruskal_wallis_test,
#   mode           = "adaptive"
#   #n_boot = 1000
# )
# 
# # AUC
# # --- Compute AUC for each ---
# auc_param  <- compute_area(sample_sizes, power_anova_param)
# auc_nonpar <- compute_area(sample_sizes, power_anova_nonpar)
# auc_adapt  <- compute_area(sample_sizes, power_anova_adapt)
# # create a data frame of AUC results
# auc_table <- data.frame(
#   Test = c("Parametric", "Nonparametric", "Adaptive"),
#   AUC  = c(auc_param, auc_nonpar, auc_adapt)
# )
# 
# # display the table
# print(auc_table)
# 
# #-------------------------------------------------------
# # save RData
# # save(
# #   sample_sizes,
# #   power_anova_param,
# #   power_anova_nonpar,
# #   power_anova_adapt,
# #   auc_table,
# #   file = "power_compare_anova_test.RData"
# # )
# #------------------------------------------------------
# 
# 
# # Plot all three curves
# plot(sample_sizes, power_anova_param,    type = "b", pch = 19, 
#      ylim = c(0,1), xlab = "Sample Size", ylab = "Power", col = "red",
#      main = "ANOVA Test Power Comparison")
# lines(sample_sizes, power_anova_nonpar, type = "b", pch = 17, col = "blue")
# lines(sample_sizes, power_anova_adapt,   type = "b", pch = 15, col = "green")
# legend("bottomright",
#        legend_labels <- c(
#          paste0("Parametric (AUC=", round(auc_param, 4), ")"),
#          paste0("Nonparametric (AUC=", round(auc_nonpar, 4), ")"),
#          paste0("Adaptive (AUC=", round(auc_adapt, 4), ")")
#        ),
#        pch    = c(19,17,15),col = c("red", "blue", "green"), lty = 1)
