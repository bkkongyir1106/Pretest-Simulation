# Set directories in local computer
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/utility.R")

setwd("/Users/benedictkongyir/Library/Mobile Documents/com~apple~CloudDocs/PhD Thesis/Expected Power loss & Expected Inflation of Error")

{
  Nsim <- 1e4
  B <- 1e3
  alpha <- 0.05
  distributions <- c("Normal", "Uniform","Exponential",  "LogNormal")
  sample_sizes <- c(8, 10, 15, 20, 25, 30, 50)
  effect_size <- 0.5
}

calculate_test_statistic <- function(x, y) {
  (mean(x) - mean(y)) / sqrt(var(x)/length(x) + var(y)/length(y))
}

permutation_test <- function(x, y, B = 1000) {
  observed <- calculate_test_statistic(x, y)
  combined <- c(x, y)
  perm_stats <- replicate(B, {
    permuted <- sample(combined)
    x_star <- permuted[1:length(x)]
    y_star <- permuted[(length(x)+1):(2*length(x))]
    calculate_test_statistic(x_star, y_star)
  })
  mean(abs(perm_stats) >= abs(observed))
}


run_simulation <- function(sample_sizes, distributions, Nsim, effect_size, alpha, B) {
  
  total_steps <- length(sample_sizes) * length(distributions) * Nsim
  current_step <- 0
  pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
  
  results <- vector("list", length(sample_sizes))
  names(results) <- as.character(sample_sizes)
  
  for (i in seq_along(sample_sizes)) {
    n <- sample_sizes[i]
    results[[i]] <- vector("list", length(distributions))
    names(results[[i]]) <- distributions
    
    for (j in seq_along(distributions)) {
      dist <- distributions[j]
      
      pval_t_H0 <- pval_perm_H0 <- pval_adaptive_H0 <- numeric(Nsim)
      pval_t_H1 <- pval_perm_H1 <- pval_adaptive_H1 <- numeric(Nsim)
      pval_sw_x <- pval_sw_y <- numeric(Nsim)
    
      for (s in 1:Nsim) {
        x <- generate_data(n, dist)
        y <- generate_data(n, dist)
        
        # Shapiro-Wilk normality tests
        pval_sw_x[s] <- shapiro.test(x)$p.value
        pval_sw_y[s] <- shapiro.test(y)$p.value
        
        # under H0
        pval_t_H0[s] <- t.test(x, y)$p.value
        pval_perm_H0[s] <- permutation_test(x, y, B)
        
        # under H1
        pval_t_H1[s] <- t.test(x, y + effect_size)$p.value
        pval_perm_H1[s] <- permutation_test(x, y + effect_size, B)
        
        # Adaptive procedure
        decision <- (pval_sw_x[s] <= alpha) | (pval_sw_y[s] <= alpha)
        pval_adaptive_H0[s] <- ifelse(decision, pval_t_H0[s], pval_perm_H0[s])
        pval_adaptive_H1[s] <- ifelse(decision, pval_t_H1[s], pval_perm_H1[s])
        
        # Update progress
        current_step <- current_step + 1
        setTxtProgressBar(pb, current_step)
      }
      
      results[[i]][[j]] <- list(
        error_t_test = mean(pval_t_H0 < alpha),
        error_perm_test = mean(pval_perm_H0 < alpha),
        error_adaptive_test = mean(pval_adaptive_H0 < alpha),
        power_t_test = mean(pval_t_H1 < alpha),
        power_perm_test = mean(pval_perm_H1 < alpha),
        power_adaptive_test = mean(pval_adaptive_H1 < alpha)
      )
    }
  }
  close(pb)
  return(results)
}

# run simulation
results <- run_simulation(
  sample_sizes = sample_sizes,
  distributions = distributions,
  Nsim = Nsim,
  effect_size = effect_size,
  alpha = alpha,
  B = B
)

# organize output
dims <- c(length(sample_sizes), length(distributions))
dimnames_list <- list(as.character(sample_sizes), distributions)

TypeIerror_t_test <- TypeIerror_perm_test <- TypeIerror_Adaptive_test <- array(NA, dim = dims, dimnames = dimnames_list
)

power_t_test <- power_perm_test <- power_Adaptive_test <- array(NA, dim = dims, dimnames = dimnames_list
)

for(i in seq_along(sample_sizes)){
  for(j in seq_along(distributions)){
    # Type I error
    TypeIerror_t_test[i, j] <- results[[i]][[j]]$error_t_test
    TypeIerror_perm_test[i, j] <- results[[i]][[j]]$error_perm_test
    TypeIerror_Adaptive_test[i, j] <- results[[i]][[j]]$error_adaptive_test
    # Power
    power_t_test[i, j] <- results[[i]][[j]]$power_t_test
    power_perm_test[i, j] <- results[[i]][[j]]$power_perm_test
    power_Adaptive_test[i, j] <- results[[i]][[j]]$power_adaptive_test
  }
}

# Power loss & Inflation of Type I error
powerloss <- power_perm_test - power_t_test
Inflation_TypeI_error <- TypeIerror_t_test - alpha

# Expected power loss & Expected Inflation of Type I error
Expected_power_loss <- power_Adaptive_test -  power_t_test
Expected_Inflation_TypeI_error <- TypeIerror_Adaptive_test - alpha

# save RData
save(
  Nsim,
  B,
  sample_sizes,
  distributions,
  effect_size,
  TypeIerror_t_test,
  TypeIerror_perm_test,
  TypeIerror_Adaptive_test,
  power_t_test,
  power_perm_test,
  power_Adaptive_test,
  powerloss,
  Inflation_TypeI_error,
  Expected_power_loss,
  Expected_Inflation_TypeI_error,
  file = "Powerloss_inflation_of_error.RData"
)

# print results
cat("Type I error rates for each test\n")
print(TypeIerror_t_test)
print(TypeIerror_perm_test)
print(TypeIerror_Adaptive_test)

cat("Power for each test\n")
print(power_t_test)
print(power_perm_test)
print(power_Adaptive_test)

cat("Power loss & Expected Power loss\n")
print(powerloss)
print(Expected_power_loss)

cat("Inflation of Type I error & Expected Inflation of Type I error\n")
print(Inflation_TypeI_error)
print(Expected_Inflation_TypeI_error)

# ---------- function to create plots -----------------
plot_test_comparison <- function(plot_data, 
                                 sample_sizes, 
                                 distributions, 
                                 y_label = "Power",
                                 main_title_prefix = "Power Comparison for",
                                 ylim = NULL,
                                 file_name = NULL) {
  
  test_colors <- c("blue", "red", "green")
  test_lty <- c(1, 2, 3)
  test_names <- c("t-test", "perm-test", "adaptive-test")
  
  # dynamically set y limits & abline
  if (is.null(ylim)) {
    if (tolower(y_label) == "power") {
      ylim <- c(0, 1)
    } else if (tolower(y_label) == "type i error") {
      ylim <- c(0, 0.1)
    } else {
      ylim <- c(0, 1)  # fallback
    }
  }
  
  # Save to PDF if file_name is provided
  if (!is.null(file_name)) {
    pdf(file_name, width = 10, height = 8)
    on.exit(dev.off())  # Ensure it closes on exit
  }
  par(mfrow = c(2, 2))
  
  for (j in seq_along(distributions)) {
    dist <- distributions[j]
    
    # Extract values for current distribution 
    values <- rbind(
      plot_data[[1]][, j],
      plot_data[[2]][, j],
      plot_data[[3]][, j]
    )
    
    # Plot the first method
    plot(
      sample_sizes, 
      values[1, ], 
      type = "o", 
      col = test_colors[1], 
      lty = test_lty[1],
      ylim = ylim, 
      xlab = "Sample Size (n)", 
      ylab = y_label,
      main = paste(main_title_prefix, dist), pch = 16
    )
    
    # Add reference line 
    if (tolower(y_label) == "power") {
      abline(h = 0.8, col = "grey", lty = 2)
    } else if (tolower(y_label) == "type i error") {
      abline(h = 0.05, col = "grey", lty = 2)
    }
    
    # Overlay remaining methods
    for (k in 2:3) {
      lines(sample_sizes, 
            values[k, ], 
            type = "o", 
            col = test_colors[k],
            lty = test_lty[k], 
            pch = 16)
    }
    # add legend
    legend("bottomright", 
           legend = test_names, 
           title = "Test Method",
           col = test_colors, 
           lty = test_lty, 
           pch = 16, 
           bty = "o")
  }
}

# power comparison
plot_test_comparison(
  plot_data = list(power_t_test, 
                   power_perm_test, 
                   power_Adaptive_test),
  sample_sizes = sample_sizes,
  distributions = distributions,
  y_label = "Power",
  file_name = "Power_Comparison.pdf"
)

# Type I error
plot_test_comparison(
  plot_data = list(TypeIerror_t_test, 
                   TypeIerror_perm_test, 
                   TypeIerror_Adaptive_test),
  sample_sizes = sample_sizes,
  distributions = distributions,
  y_label = "Type I Error",
  file_name = "TypeIError_Comparison.pdf"
)

# --------------------------------------------
plot_difference_metrics <- function(metric_matrix_list, 
                                    sample_sizes, 
                                    distributions, 
                                    y_label = "Power Loss",
                                    main_title_prefix = "Difference Plot for",
                                    colors = c("blue", "red"),
                                    lty = c(1, 2),
                                    legend_labels = c("Observed", "Expected"),
                                    ref_line = 0,
                                    file_name = NULL,
                                    ylim = NULL) {
  
  # Determine y-axis range if not provided
  if (is.null(ylim)) {
    # Automatically range from both matrices
    combined <- do.call(rbind, metric_matrix_list)
    ylim <- range(combined, na.rm = TRUE)
  }
  
  # Save to PDF if filename is given
  if (!is.null(file_name)) {
    pdf(file_name, width = 10, height = 8)
    on.exit(dev.off())
  }
  
  par(mfrow = c(2, 2))
  
  for (j in seq_along(distributions)) {
    dist <- distributions[j]
    
    # Extract the two lines for this distribution
    observed <- metric_matrix_list[[1]][, j]
    expected <- metric_matrix_list[[2]][, j]
    
    plot(sample_sizes, observed, type = "o", col = colors[1], lty = lty[1], pch = 16,
         ylim = ylim, xlab = "Sample Size (n)", ylab = y_label,
         main = paste(main_title_prefix, dist))
    
    lines(sample_sizes, expected, type = "o", col = colors[2], lty = lty[2], pch = 16)
    
    abline(h = ref_line, col = "grey", lty = 2)
    
    legend("topright", legend = legend_labels, col = colors, lty = lty, pch = 16, bty = "o")
  }
}

plot_difference_metrics(
  metric_matrix_list = list(powerloss, Expected_power_loss),
  sample_sizes = sample_sizes,
  distributions = distributions,
  y_label = "Power Loss",
  main_title_prefix = "Power Loss vs Expected Power Loss for",
  file_name = "PowerLoss_vs_Expected.pdf"
)


plot_difference_metrics(
  metric_matrix_list = list(Inflation_TypeI_error, Expected_Inflation_TypeI_error),
  sample_sizes = sample_sizes,
  distributions = distributions,
  y_label = "Inflation of Type I Error",
  main_title_prefix = "Inflation vs Expected Inflation for",
  file_name = "TypeIErrorInflation_vs_Expected.pdf"
)
