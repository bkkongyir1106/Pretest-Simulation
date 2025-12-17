library(MASS)
source("~/Desktop/OSU/Research/Pretest-Simulation/functions/user_defined_functions_center_median.R")

# Box-Cox transformation 
boxcox_simple <- function(data) {
  # Shift data to be positive
  shift <- 0
  if (min(data) <= 0) {
    shift <- abs(min(data)) + 0.001
    data <- data + shift
  }
  
  # Find optimal lambda 
  bc <- suppressMessages(boxcox(data ~ 1, lambda = seq(-2, 2, 0.1), plotit = FALSE))
  lambda <- bc$x[which.max(bc$y)]
  
  # Transform data
  if (lambda == 0) {
    transformed_data <- log(data)
  } else {
    transformed_data <- (data^lambda - 1) / lambda
  }
  
  return(list(
    transformed_data = transformed_data,
    lambda = lambda,
    shift = shift
  ))
}

# T-test on median using Box-Cox transformed data
test_t_on_median <- function(data) {
  # Since data is already centered around median, true_median = 0
  true_median <- 0
  
  # Transform data
  bc <- boxcox_simple(data)
  
  # Transform the true median 
  true_median_shifted <- true_median + bc$shift
  if (bc$lambda == 0) {
    transformed_median <- log(true_median_shifted)
  } else {
    transformed_median <- (true_median_shifted^bc$lambda - 1) / bc$lambda
  }
  
  # T-test on transformed data against transformed median
  t.test(bc$transformed_data, mu = transformed_median)$p.value
}

# Calculate Type I error rate for different distributions
type1_error_all_distributions <- function(n_sim = 1000) {
  sizes <- c(10, 20, 30, 50)
  distributions <- c("normal", "lognormal", "exponential")
  
  # Initialize results
  results <- data.frame(
    distribution = rep(distributions, each = length(sizes)),
    n = rep(sizes, times = length(distributions)),
    error_rate = NA,
    avg_lambda = NA
  )
  
  # Create a list to store all lambda values for histogram plotting
  lambda_storage <- list()
  
  for (dist in distributions) {
    lambda_storage[[dist]] <- list()
    
    for (i in 1:length(sizes)) {
      n <- sizes[i]
      errors <- 0
      lambdas <- numeric(n_sim)
      
      for (j in 1:n_sim) {
        # Generate data centered around median using your function
        data <- generate_data(n, dist)
        
        # T-test on median using Box-Cox
        p_value <- test_t_on_median(data)
        if (p_value < 0.05) errors <- errors + 1
        
        lambdas[j] <- boxcox_simple(data)$lambda
      }
      
      # Store results
      idx <- which(results$distribution == dist & results$n == n)
      results$error_rate[idx] <- errors / n_sim
      results$avg_lambda[idx] <- mean(lambdas)
      
      # Store lambdas for histogram
      lambda_storage[[dist]][[as.character(n)]] <- lambdas
    }
  }
  
  return(list(results = results, lambdas = lambda_storage))
}

# Combined plot for all distributions (MISSING FUNCTION ADDED)
plot_combined <- function(results) {
  distributions <- unique(results$distribution)
  colors <- c("red", "blue", "green", "purple")
  # Get the sample sizes from the results
  sizes <- unique(results$n)
  # Set up the plot
  plot(0, 0, type = "n", 
       xlim = range(sizes), ylim = c(0, 1),
       xlab = "Sample Size", ylab = "Type I Error Rate",
       main = "Type I Error Rates - Box-Cox Transformed T-test")
  
  # Add line for each distribution
  for (i in 1:length(distributions)) {
    dist <- distributions[i]
    dist_data <- results[results$distribution == dist, ]
    lines(dist_data$n, dist_data$error_rate, 
          type = "b", col = colors[i], lwd = 2, pch = 16)
  }
  
  # Add horizontal line at nominal level
  abline(h = 0.05, lty = 2, col = "gray", lwd = 2)
  
  # Add legend
  legend("topleft", legend = distributions, 
         col = colors, lwd = 2, pch = 16, bty = "n")
}

# Function to create lambda histograms
plot_lambda_histograms <- function(lambda_storage) {
  distributions <- names(lambda_storage)
  sizes <- c(10, 20, 30, 50)
  
  # Set up plotting parameters
  par(mfrow = c(length(distributions), length(sizes)), 
      mar = c(2, 2, 2, 1), oma = c(4, 4, 4, 1))
  
  for (dist in distributions) {
    for (n in sizes) {
      lambdas <- lambda_storage[[dist]][[as.character(n)]]
      
      # Create histogram
      hist(lambdas, breaks = 20, main = "", 
           xlab = "", ylab = "", col = "lightblue",
           xlim = c(-2, 2))
      
      # Add title for first row only
      if (dist == distributions[1]) {
        title(main = paste("n =", n), line = 0.5, cex.main = 0.8)
      }
      
      # Add y-axis label for first column only
      if (n == sizes[1]) {
        mtext(dist, side = 2, line = 2.5, cex = 0.7)
      }
      
      # Add mean lambda line
      abline(v = mean(lambdas), col = "red", lwd = 2, lty = 2)
      
      # Add lambda = 1 line (no transformation)
      abline(v = 1, col = "darkgreen", lwd = 1, lty = 3)
      
      # Add lambda = 0 line (log transformation)
      abline(v = 0, col = "blue", lwd = 1, lty = 3)
    }
  }
  
  # Add overall titles
  mtext("Lambda Distribution from Box-Cox Transformation", 
        side = 3, outer = TRUE, line = 1, cex = 1.2)
  mtext("Lambda Value", side = 1, outer = TRUE, line = 2, cex = 0.8)
  mtext("Frequency", side = 2, outer = TRUE, line = 2, cex = 0.8)
  
  # Reset plotting parameters
  par(mfrow = c(1, 1))
}

# Alternative: Individual histogram plots for each distribution
plot_lambda_histograms_by_distribution <- function(lambda_storage) {
  distributions <- names(lambda_storage)
  sizes <- c(10, 20, 30, 50)
  colors <- c("red", "blue", "green", "purple", "orange")
  
  for (dist in distributions) {
    # Create a PDF for each distribution
    pdf(file = paste0("lambda_histograms_", dist, ".pdf"), width = 10, height = 6)
    
    par(mfrow = c(2, 3), mar = c(4, 4, 3, 1))
    
    for (i in 1:length(sizes)) {
      n <- sizes[i]
      lambdas <- lambda_storage[[dist]][[as.character(n)]]
      
      hist(lambdas, breaks = 20, 
           main = paste(dist, "- n =", n),
           xlab = "Lambda", ylab = "Frequency",
           col = colors[i],
           xlim = c(-2, 2))
      
      abline(v = mean(lambdas), col = "red", lwd = 2, lty = 2)
      abline(v = 1, col = "darkgreen", lwd = 1, lty = 3)
      abline(v = 0, col = "blue", lwd = 1, lty = 3)
      
      legend("topright", 
             legend = c(paste("Mean =", round(mean(lambdas), 3)),
                        "Lambda = 1", "Lambda = 0"),
             col = c("red", "darkgreen", "blue"),
             lty = c(2, 3, 3), lwd = c(2, 1, 1),
             cex = 0.7, bty = "n")
    }
    
    # Add one more panel with summary statistics
    plot(1, type = "n", xaxt = "n", yaxt = "n", 
         xlab = "", ylab = "", bty = "n")
    
    summary_text <- c(paste("Distribution:", dist),
                      "Summary Statistics:",
                      paste("Sample sizes:", paste(sizes, collapse = ", ")),
                      "")
    
    # Add mean lambdas for each sample size
    for (n in sizes) {
      lambdas <- lambda_storage[[dist]][[as.character(n)]]
      summary_text <- c(summary_text, 
                        paste("n =", n, ": mean lambda =", round(mean(lambdas), 3)))
    }
    
    text(1, 1, paste(summary_text, collapse = "\n"), cex = 0.9)
    
    dev.off()
  }
}

# Run the simulation
set.seed(123)
simulation_results <- type1_error_all_distributions(n_sim = 1000)
all_results <- simulation_results$results
lambda_storage <- simulation_results$lambdas

print(all_results)

# Create the combined Type I error plot
pdf(file = "one_sample_t_test_vs_bc_t_test.pdf", width = 6, height = 5)
plot_combined(all_results)
dev.off()

# Create lambda histograms
pdf(file = "lambda_distributions_all.pdf", width = 12, height = 8)
plot_lambda_histograms(lambda_storage)
dev.off()

# Create individual histogram plots for each distribution
plot_lambda_histograms_by_distribution(lambda_storage)

# Also create a summary table of lambda statistics
lambda_summary <- data.frame()
for (dist in names(lambda_storage)) {
  for (n in c(10, 20, 30, 50)) {
    lambdas <- lambda_storage[[dist]][[as.character(n)]]
    lambda_summary <- rbind(lambda_summary, 
                            data.frame(
                              Distribution = dist,
                              SampleSize = n,
                              MeanLambda = mean(lambdas),
                              SdLambda = sd(lambdas),
                              MinLambda = min(lambdas),
                              MaxLambda = max(lambdas)
                            ))
  }
}

print("Lambda Summary Statistics:")
print(lambda_summary)

# Save lambda summary to CSV
write.csv(lambda_summary, "lambda_summary_statistics.csv", row.names = FALSE)