load("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach/ML vs SW test/ROC_SW_n10.RData")
load("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach/ML vs SW test/trained_models_n10.RData")
ml_n10 <- eval_results
sw_n10 <- roc_results

load("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach/ML vs SW test/trained_models_n30.RData")
load("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach/ML vs SW test/ROC_SW_n30.RData")
ml_n30 <- eval_results
sw_n30 <- roc_results

setwd("~/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach/ML vs SW test")

plot_combined_roc_base <- function(eval_results, roc_results, main_title = "ROC Curves (Positive Class: Non_Normal)") {
  # Prepare plotting symbols, colors, and line types
  model_names <- names(eval_results)
  colors <- c("cyan", "brown", "blueviolet", "red", "green", "blue", "orange")
  line_types <- c(1, 1, 1, 1, 1, 1, 2)  # Solid for ML models, dashed for normality tests
  point_chars <- c(16, 17, 18, 15, 1, 2, 3)  # Different point characters for each
  
  # Plot frame
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate (1 - Specificity)",
       ylab = "True Positive Rate (Sensitivity)",
       main = main_title)
  # Diagonal reference line
  abline(a = 0, b = 1, col = "gray", lty = 2) 
  
  # add legend
  legend_labels <- c()
  legend_colors <- c()
  legend_lty <- c()
  legend_pch <- c()
  
  # Plot machine learning models
  for (i in seq_along(model_names)) {
    model_name <- model_names[i]
    pred_df <- eval_results[[model_name]]
    
    # Skip if prediction probabilities are missing
    if (!"Prob_Non_Normal" %in% colnames(pred_df$Predictions)) next
    
    actual_bin <- ifelse(pred_df$Predictions$True_Class == "Non_Normal", 1, 0)
    probs <- pred_df$Predictions$Prob_Non_Normal
    
    roc_obj <- pROC::roc(actual_bin, probs, levels = c(0, 1), direction = "<", quiet = TRUE)
    auc_val <- round(pROC::auc(roc_obj), 3)
    
    # Plot line
    lines(1 - roc_obj$specificities, roc_obj$sensitivities, 
          col = colors[i], lwd = 2, lty = line_types[i])  
    
    # Add points at regular intervals
    n_points <- 10
    idx <- seq(1, length(roc_obj$sensitivities), length.out = n_points)
    points(1 - roc_obj$specificities[idx], roc_obj$sensitivities[idx], 
           col = colors[i], pch = point_chars[i], cex = 1.2)
    
    legend_labels <- c(legend_labels, paste0(model_name, " (AUC=", auc_val, ")"))
    legend_colors <- c(legend_colors, colors[i])
    legend_lty <- c(legend_lty, line_types[i])
    legend_pch <- c(legend_pch, point_chars[i])
  }
  
  # Function to calculate AUC using trapezoidal rule (base R)
  calculate_auc <- function(fpr, tpr) {
    # Order by FPR
    order_idx <- order(fpr)
    fpr_ordered <- fpr[order_idx]
    tpr_ordered <- tpr[order_idx]
    
    # Calculate AUC using trapezoidal rule
    auc <- 0
    for (i in 2:length(fpr_ordered)) {
      width <- fpr_ordered[i] - fpr_ordered[i-1]
      avg_height <- (tpr_ordered[i] + tpr_ordered[i-1]) / 2
      auc <- auc + width * avg_height
    }
    return(auc)
  }
  
  # Add SW test ROC curve
  sw_idx <- length(model_names) + 1
  sw_fpr <- roc_results$FPR["SW", ]
  sw_tpr <- roc_results$TPR["SW", ]
  sw_auc <- round(calculate_auc(sw_fpr, sw_tpr), 3)
  lines(sw_fpr, sw_tpr, col = colors[sw_idx], lwd = 2, lty = line_types[sw_idx])
  # Add points for SW test
  n_points <- 10
  idx <- seq(1, length(sw_fpr), length.out = n_points)
  points(sw_fpr[idx], sw_tpr[idx], col = colors[sw_idx], pch = point_chars[sw_idx], cex = 1.2)
  legend_labels <- c(legend_labels, paste0("SW Test (AUC=", sw_auc, ")"))
  legend_colors <- c(legend_colors, colors[sw_idx])
  legend_lty <- c(legend_lty, line_types[sw_idx])
  legend_pch <- c(legend_pch, point_chars[sw_idx])
  
  # add legend
  legend("bottomright",
         legend = legend_labels,
         col = legend_colors,
         lty = legend_lty,
         pch = legend_pch,
         lwd = 2,
         bty = "o",  
         cex = 0.85,
         title = "Models & Tests")
  
  grid(col = "lightgray", lty = "dotted")
}

# Now run the function
pdf("roc_curve_sw.pdf", width=10, height=6)
par(mfrow = c(1, 2))
plot_combined_roc_base(eval_results = ml_n10, roc_results= sw_n10, main_title = "Sample Size = 10" )
plot_combined_roc_base(eval_results = ml_n30, roc_results= sw_n30 , main_title = "Sample Size = 30")

dev.off()
