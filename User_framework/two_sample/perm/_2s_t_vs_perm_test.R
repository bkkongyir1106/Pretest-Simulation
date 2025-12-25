setwd("~/Desktop/OSU/Research/Pretest-Simulation/User_framework_v1/two_sample/perm")
## keep startup clean — no extra packages auto-attached
options(defaultPackages = c("datasets","utils","grDevices","graphics","stats","methods"))
Sys.setenv(R_DEFAULT_PACKAGES = "NULL")

## ====== Cluster-friendly args, working dir, and sourcing ======

args <- commandArgs(trailingOnly = TRUE)
get_arg <- function(name, default = NULL) {
  i <- grep(paste0("^--", name, "="), args)
  if (length(i)) sub(paste0("^--", name, "="), "", args[i[1]]) else default
}

# Working directory where the scripts live (defaults to current dir if not passed)
wd  <- get_arg("wd",  getwd())
# Output directory for all results (PDF/RData/etc.)
out <- get_arg("out", file.path(wd, "results"))

dir.create(out, showWarnings = FALSE, recursive = TRUE)
setwd(wd)
cat("Working directory:", wd,  "\n")
cat("Output directory: ", out, "\n")

# Source helper script sitting next to this file
helper <- file.path(wd, "funs_4_2s_t_vs_perm_test.R")
if (!file.exists(helper)) stop("Helper script not found: ", helper)
source(helper)

# ====================================================================

## ====== Minimal, safe PDF wrapper (closes even on error) ======
save_pdf <- function(file, width, height, plot_fun) {
  grDevices::pdf(file, width = width, height = height)
  on.exit(try(grDevices::dev.off(), silent = TRUE), add = TRUE)
  plot_fun()
  invisible(TRUE)
}
# ====================================================================
# --------------- Run the simulation -------------------
run_simulation <- function(n, N,  Nsim, test_type, distributions, tol_pos , alpha_pre) {
  
  ## Keep a small results 
  results <- list(
    params  = list(N = N, Nsim = Nsim, test_type = test_type, distributions = distributions , alpha_pre = alpha_pre),
    objects = list(),
    files   = list()
  )
  
  sample_size       <- c(10, 20, 30, 40, 50)
  test_alpha        <- 0.05
  alpha             <- 0.05
  select_norm_test  <- "SW"
  effect_size       <- 0.5
  norm_alpha_pre    <- seq(from = 0.0, to = 1, by = 0.05)
  alpha_pretest     <- seq(from = 0.009, to = 1, by = 0.0025)
  effect_sizes      <- c(0.0, 0.1, 0.2, 0.3, 0.5, 0.8, 1.0)
  sig_levels        <- seq(from = 0.005, to = 1, by = 0.005)
  ds_test_methods   <- c("test_1", "test_2", "adaptive")
  norm_test         <- c("SW", "SF", "KS", "JB","SKEW", "DAP", "AD", "CVM")
  
  ## ROC data for normality tests 
  roc_pval_ds_test <- fn_for_roc_curve_for_norm_test(
    n = 10,
    alpha_pretest = norm_alpha_pre,
    H1_dist = distributions[1],
    tests = norm_test,
    Nsim = Nsim
  )
  results$objects$roc_pval_ds_test <- roc_pval_ds_test
  
  ## Plot ROC using selected tests
  results$files$norm_roc <- file.path(out, paste(test_type, "_test_norm_roc_curve.pdf"))
  save_pdf(results$files$norm_roc, 7, 6, function() {
    selected_tests <- c("SW", "SF", "KS", "JB","SKEW", "DAP", "AD", "CVM")
    plot_norm_roc_curve(
      FPR = roc_pval_ds_test$FPR,
      TPR = roc_pval_ds_test$TPR,
      tests_to_plot = selected_tests,
      alpha = roc_pval_ds_test$alpha,
      dist_name = distributions[1]
    )
  })
  
  # Trade-off analysis (power & Type I error) 
  analysis_ds_tests <- perform_analysis(
    Nsim = Nsim,
    n = 10,
    distributions = distributions,
    effect_size = effect_size,
    test = select_norm_test,
    alpha_pretest = alpha_pretest,
    test_alpha = test_alpha
  )
  results$objects$analysis_ds_tests <- analysis_ds_tests
  
  metrics <- compute_roc_metrics(
    error_ds_test = analysis_ds_tests$error_ds_test,
    power_ds_test = analysis_ds_tests$power_ds_test,
    test_alpha    = test_alpha
  )
  results$objects$metrics <- metrics
  
  # Capture the trade off result to get alpha_star
  results$files$tradeoff <- file.path(out, paste(test_type, "_test_tradeoff.pdf"))
  tradeoff_result <- NULL
  save_pdf(results$files$tradeoff, 7, 6, function() {
    tradeoff_result <<- plot_power_error_tradeoff(
      distributions = distributions,
      alpha_pretest = alpha_pretest,
      metrics = metrics,
      tol_pos = tol_pos
    )
  })
  
  # Save the trade off result and extract alpha_star
  results$objects$tradeoff_result <- tradeoff_result
  alpha_star <- tradeoff_result$alpha_star
  results$objects$alpha_star <- alpha_star
  
  cat("=== Selected alpha_star from tradeoff analysis:", round(alpha_star, 4), "===\n")
  
  ## Power vs Type I error ROC-style plots 
  roc_data <- power_vs_error_roc_data(
    N = N, 
    n = 10,
    distributions = distributions, 
    effect_size = effect_size,
    alpha_pretest = alpha_star,  
    sig_levels = sig_levels
  )
  results$objects$roc_data <- roc_data
  
  results$files$power_error_roc <- file.path(out, paste(test_type, "_test_power_error_roc.pdf"))
  save_pdf(results$files$power_error_roc, 6, 5, function() {
    power_vs_error_ROC_curve(roc_data = roc_data)
  })
  
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
    alpha_pre          = alpha_star  # Use alpha_star here
  )
  results$objects$sim_output <- sim_output
  
  ## Per-distribution AUCs (printed and saved) 
  all_auc_tables <- list()
  for (dist in names(sim_output)) {
    cat("\n\n=== Results for", dist, "distribution ===\n")
    dist_results  <- sim_output[[dist]]
    power_results <- dist_results$power
    type1_results <- dist_results$type1
    
    power_df <- list(
      test_1    = power_results$test_1,
      test_2 = power_results$test_2,
      adaptive      = power_results$adaptive
    )
    TypeI_error_df <- list(
      test_1    = type1_results$test_1,
      test_2 = type1_results$test_2,
      adaptive      = type1_results$adaptive
    )
    
    auc_power <- sapply(power_df, function(y) compute_area(sample_size, y))
    auc_type1 <- sapply(TypeI_error_df, function(y) compute_area(sample_size, y))
    
    auc_table <- data.frame(
      Method    = names(auc_power),
      AUC_Power = unname(auc_power),
      AUC_TypeI = unname(auc_type1),
      row.names = NULL
    )
    all_auc_tables[[dist]] <- auc_table
    print(auc_table, digits = 4)
  }
  results$objects$auc_tables <- all_auc_tables
  
  ## Combine across distributions and plot -
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
  
  results$files$power_type1 <- file.path(out, paste(test_type, "_power_and_error_test.pdf"))
  
  save_pdf(results$files$power_type1, 6, 5, function() {
    plot_power_type1(
      combined_power  = combined_power,
      combined_type1  = combined_type1,
      optimal_alphas  = alpha_star,        
      methods         = ds_test_methods,   
      distributions   = distributions,
      sizes           = sample_size,
      alpha           = alpha
    )
  })
  
  
  
  ## Effect-size power curves 
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
      alpha_pre          = alpha_star  
    )
  }
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
  combined_power_effect <- do.call(rbind, all_power_df)
  results$objects$combined_power_effect <- combined_power_effect
  
  # Where to save
  results$files$power_by_effect <- file.path(out, paste0(test_type, "_power_by_effect.pdf"))
  # Open device with desired size, plot, then close
  pdf(results$files$power_by_effect, width = 6, height = 4)
  plot_power_by_effect_size(
    power_results   = combined_power_effect,
    distributions   = distributions,
    ds_test_methods = ds_test_methods,
    effect_sizes    = effect_sizes,
    alpha_pre       = alpha_star
  )
  
  dev.off()
  
  
  ## ---------- Save compact results + full workspace ----------
  rdata_results <- file.path(out, paste0(test_type, "_results.RData"))
  save(results, file = rdata_results)
  
  rdata_workspace <- file.path(out, paste0(test_type, "_workspace.RData"))
  save.image(rdata_workspace)
  
  # Final summary
  cat("\n=== SIMULATION COMPLETE ===\n")
  cat("Selected alpha_star used throughout:", round(alpha_star, 4), "\n")
  cat("Tradeoff feasibility:", if(tradeoff_result$feasible) "FEASIBLE" else "FALLBACK", "\n")
  cat("Worst-case positive inflation:", round(tradeoff_result$worst_case_positive_inflation, 6), "\n")
  cat("Inflation values (Normal, Non-normal):",
      round(tradeoff_result$inflation_normal, 6), ",",
      round(tradeoff_result$inflation_non_normal, 6), "\n")
  
  
  invisible(results)
}
## Example call
run_simulation(
  Nsim = 1e3, N = 1e4,
  test_type = "two_sample_t_vs_perm",
  distributions = c("exponential", "normal"),
  tol_pos = 1e-2,
  alpha_pre = 0.05
)

