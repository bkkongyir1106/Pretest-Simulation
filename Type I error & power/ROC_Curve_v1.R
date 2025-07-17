# Load your custom functions
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/User_defined_functions.R")

# Parameters
Nsim <- 1e5
sample_size <- 10
alpha.level <- seq(from = 0, to = 1, by = 0.01)
tests <- c("SW", "KS", "AD", "DAP", "SF", "JB", "CVM", "ANS", "AD2", "SKEW", "KURT")

# Initialize FPR and TPR matrices
FPR <- matrix(0, nrow = length(tests), ncol = length(alpha.level))
TPR <- matrix(0, nrow = length(tests), ncol = length(alpha.level))
rownames(FPR) <- rownames(TPR) <- tests
colnames(FPR) <- colnames(TPR) <- paste0("alpha_", alpha.level)

# Main simulation loop
for (i in seq_along(tests)) {
  test <- tests[i]
  
  for (j in seq_along(alpha.level)) {
    alpha <- alpha.level[j]
    
    pval_norm <- pval_non_normal <- numeric(Nsim)
    
    for (k in 1:Nsim) {
      normal_data <- generate_data(sample_size, dist = "Normal", par = NULL)
      non_normal_data <- generate_data(sample_size, dist = "LogNormal", par = NULL)  
      
      pval_norm[k] <- generate_tests(normal_data, test)$p.value
      pval_non_normal[k] <- generate_tests(non_normal_data, test)$p.value
    }
    
    FPR[i, j] <- mean(pval_norm < alpha)
    TPR[i, j] <- mean(pval_non_normal < alpha)
  }
}


# plot ROC curve
selected_tests <- c("SW", "KS", "AD", "DAP", "SF", "JB", "CVM",  "AD2", "SKEW") #, "ANS", "KURT")
plot_character <- 1:length(selected_tests)
colors <- rainbow(length(selected_tests))

plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
     xlab = "False Positive Rate (FPR)", ylab = "True Positive Rate (TPR)",
     main = "ROC Curves for Different Normality Tests")

for (i in seq_along(selected_tests)) {
  test <- selected_tests[i]
  lines(FPR[test, ], TPR[test, ], col = colors[i], lwd = 2)
  points(FPR[test, ], TPR[test, ], col = colors[i], pch = plot_character[i], cex = 0.5)
}

abline(0, 1, lty = 2, col = "gray")

legend("bottomright", legend = selected_tests,
       col = colors, pch = plot_character, lwd = 2, title = "Normality Tests")

