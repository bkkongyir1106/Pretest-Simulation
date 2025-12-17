# Load custom test and data functions
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")

setwd("~/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach/Machine Learning Pipeline")
# Parameters
Nsim <- 1e5
sample_size <- 10
alpha_pretest <- seq(from = 0.0, to = 1, by = 0.05)
tests <- c("SW", "KS", "AD", "DAP", "SF", "JB", "CVM", "SKEW") 

# -----------------------------------------------------------------------------
# Function to compute FPR and TPR for a list of tests over alpha thresholds
# -----------------------------------------------------------------------------
ROC_curve_function <- function(sample_size, alpha_pretest, tests, Nsim = 1e3) {
  FPR <- matrix(0, nrow = length(tests), ncol = length(alpha_pretest))
  TPR <- matrix(0, nrow = length(tests), ncol = length(alpha_pretest))
  rownames(FPR) <- rownames(TPR) <- tests
  colnames(FPR) <- colnames(TPR) <- paste0("alpha_", alpha_pretest)
  # set progress bar
  total_steps <- length(tests) * length(alpha_pretest) * Nsim
  current_step <- 0
  pb <- txtProgressBar(min = 0, max = total_steps, style = 3)
  
  for (i in seq_along(tests)) {
    test_name <- tests[i]
    
    for (j in seq_along(alpha_pretest)) {
      alpha <- alpha_pretest[j]
      pval_norm <- pval_non_normal <- numeric(Nsim)
      
      for (k in 1:Nsim) {
        # generate data
        normal_data <- generate_data(sample_size, dist = "Normal", par = NULL)
        non_normal_data <- generate_data(sample_size, dist = "chi_square", par = NULL)
        # calculate p-values
        pval_norm[k] <- generate_tests(normal_data, test_name)$p.value
        pval_non_normal[k] <- generate_tests(non_normal_data, test_name)$p.value
        
        # Update progress
        current_step <- current_step + 1
        setTxtProgressBar(pb, current_step)
      }
      # calculate FPR & TPR
      FPR[i, j] <- mean(pval_norm < alpha)
      TPR[i, j] <- mean(pval_non_normal < alpha)
      
    }
  }
  close(pb)
  return(list(FPR = FPR, 
              TPR = TPR, 
              alpha = alpha_pretest
              ))
}

# Run the simulation
roc_results <- ROC_curve_function(
  sample_size = sample_size,
  alpha_pretest = alpha_pretest,
  tests = tests,
  Nsim = Nsim
)


# -------------------------------------------------------
# Function to plot ROC curves from FPR and TPR matrices
# -------------------------------------------------------
plot_ROC <- function(FPR, 
                     TPR, 
                     tests_to_plot = rownames(FPR), 
                     alpha = NULL,
                     title = "ROC Curves for Different Normality Tests") {
  
  colors <- 1: length(tests_to_plot)
  plot_chars <- 1:length(tests_to_plot)
  
  #Set margins and text sizes
  par(mar = c(5, 5, 4, 2))  
  plot(0, 0, type = "n",
       xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate (FPR)", 
       ylab = "True Positive Rate (TPR)",
       main = title)
  
  for (i in seq_along(tests_to_plot)) {
    test <- tests_to_plot[i]
    lines(FPR[test, ], TPR[test, ], col = colors[i], lwd = 2)
    points(FPR[test, ], TPR[test, ], col = colors[i], pch = plot_chars[i], cex = 0.75)
  }
  
  abline(0, 1, lty = 2, col = "gray")
  legend("bottomright", 
         legend = tests_to_plot, 
         col = colors, 
         pch = plot_chars,
         lwd = 2, 
         title = "Normality Tests")
}

# Plot ROC using selected tests
pdf("ROC_curves_for_normality_tests_10_chi_sq.pdf", width = 8, height = 6)
selected_tests <- c("SW", "KS", "AD", "DAP", "SF", "JB", "CVM", "SKEW")
plot_ROC(FPR = roc_results$FPR,
         TPR = roc_results$TPR,
         tests_to_plot = selected_tests,
         alpha = roc_results$alpha)

dev.off()
# -------- save RData -----------
save(
  Nsim,
  sample_size,
  alpha_pretest,
  tests,
  roc_results,
  file = "ROC_curve_data_10_chi_sq.RData"
)

