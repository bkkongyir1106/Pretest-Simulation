## Compare SW ROC vs. ML ROC curves (base R)
## - Expects: FPR_list, TPR_list in env (from your RData), and eval_results from your ML code
## - Only uses SW from the RData; ignores KS/AD/DAP/SF/JB/CVM
## - Positive class is "Non_Normal" (matches your pipeline)

setwd("~/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach")
load("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/NORMALITY TEST METHODS/ROC_Normality.test.RData")
load("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach/trained_models.RData")

plot_ml_vs_sw_roc <- function(eval_results,
                              FPR_list, TPR_list,
                              sw_key = "SW",
                              legend_loc = "bottomright",
                              file = "SW_test_ML_models.pdf") {
  
  # Optional PDF export
  if (!is.null(file)) {
    pdf(file, width = 8, height = 6)
    on.exit(dev.off(), add = TRUE)
  }
  
  # compute AUC from an ROC curve (FPR, TPR) via trapezoids ---
  roc_auc_from_curve <- function(fpr, tpr) {
    ok <- is.finite(fpr) & is.finite(tpr)
    fpr <- fpr[ok]
    tpr <- tpr[ok]
    ord <- order(fpr)
    f <- fpr[ord]
    t <- tpr[ord]
    # ensure endpoints are present
    if (length(f) == 0) return(NA_real_)
    if (f[1] > 0 || t[1] > 0) { f <- c(0, f); t <- c(0, t) }
    if (tail(f, 1) < 1 || tail(t, 1) < 1) { f <- c(f, 1); t <- c(t, 1) }
    sum(diff(f) * (head(t, -1) + tail(t, -1)) / 2)
  }
  
  # --- Get SW curve from lists ---
  if (!sw_key %in% names(FPR_list) || !sw_key %in% names(TPR_list)) {
    stop(sprintf("Couldn't find '%s' in FPR_list/TPR_list.", sw_key))
  }
  sw_fpr <- as.numeric(FPR_list[[sw_key]])
  sw_tpr <- as.numeric(TPR_list[[sw_key]])
  sw_auc <- roc_auc_from_curve(sw_fpr, sw_tpr)
  
  # --- Prepare ML model names (exclude any normality-test placeholders) ---
  drop_names <- c("Shapiro-Wilk", "SW", "KS", "AD", "DAP", "SF", "JB", "CVM")
  ml_names <- setdiff(names(eval_results), drop_names)
  if (length(ml_names) == 0) {
    warning("No ML models detected in eval_results after filtering; plotting SW only.")
  }
  
  # --- Base frame ---
  plot(0, 0, type = "n", xlim = c(0, 1), ylim = c(0, 1),
       xlab = "False Positive Rate (1 - Specificity)",
       ylab = "True Positive Rate (Sensitivity)",
       main = "ROC: Shapiro–Wilk vs. Machine Learning (Positive: Non_Normal)")
  abline(a = 0, b = 1, col = "gray", lty = 2)
  grid(col = "lightgray", lty = "dotted")
  
  # --- Plot ML curves ---
  if (length(ml_names) > 0) {
    cols <- setNames(rainbow(length(ml_names)), ml_names)
    legend_labels <- character(0)
    legend_cols   <- character(0)
    legend_lty    <- integer(0)
    legend_lwd    <- numeric(0)
    
    for (model_name in ml_names) {
      pred_df <- eval_results[[model_name]]$Predictions
      if (!("Prob_Non_Normal" %in% colnames(pred_df))) next
      y_true <- ifelse(pred_df$True_Class == "Non_Normal", 1, 0)
      scores <- pred_df$Prob_Non_Normal
      keep <- is.finite(scores) & !is.na(scores)
      if (!any(keep)) next
      
      roc_obj <- pROC::roc(y_true[keep], scores[keep], levels = c(0, 1),
                           direction = "<", quiet = TRUE)
      auc_val <- as.numeric(round(pROC::auc(roc_obj), 3))
      
      lines(1 - roc_obj$specificities, roc_obj$sensitivities,
            col = cols[model_name], lty = 1, lwd = 2)
      
      legend_labels <- c(legend_labels, sprintf("%s (AUC=%.3f)", model_name, auc_val))
      legend_cols   <- c(legend_cols, cols[model_name])
      legend_lty    <- c(legend_lty, 1)
      legend_lwd    <- c(legend_lwd, 2)
    }
  } else {
    legend_labels <- legend_cols <- character(0)
    legend_lty <- legend_lwd <- numeric(0)
  }
  
  # --- Plot SW curve (black dashed) ---
  lines(sw_fpr, sw_tpr, col = "black", lwd = 3, lty = 2)
  sw_label <- sprintf("Shapiro–Wilk (AUC=%.3f)", sw_auc)
  
  # --- Legend (SW + ML) ---
  legend(legend_loc,
         legend = c(sw_label, legend_labels),
         col    = c("black", legend_cols),
         lty    = c(2, legend_lty),
         lwd    = c(3, legend_lwd),
         bty    = "o", cex = 0.9,
         title  = "Models / Test")
}

## ---- Usage ---------------------------------------------------------------
## 1) Make sure you've run your ML pipeline and have `eval_results`
## 2) Load the ROC points for tests (already in your snippet):
##    load("/Users/benedictkongyir/Desktop/OSU/Research/Pretest-Simulation/Simulation/NORMALITY TEST METHODS/ROC_Normality.test.RData")
##    (provides FPR_list, TPR_list with an entry "SW")
## 3) Call:
plot_ml_vs_sw_roc(eval_results, FPR_list, TPR_list, sw_key = "SW")
# Optional to save as PDF:
plot_ml_vs_sw_roc(eval_results, FPR_list, TPR_list, sw_key = "SW",
                  file = "ROC_SW_vs_ML.pdf")

