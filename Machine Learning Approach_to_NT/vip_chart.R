rm(list = ls())
pacman::p_load(caret, pROC, ROCR, randomForest)

setwd("~/Desktop/OSU/Research/Pretest-Simulation/Machine Learning Approach_to_NT")
# Function to load and process a model
process_model <- function(filename, suffix) {
  load(filename)
  var_imp <- varImp(models_list$RF)
  imp_values <- var_imp$importance$Overall
  
  # Get the order for sorting
  ord <- order(imp_values, decreasing = FALSE)
  
  return(list(imp = imp_values, ord = ord, var_imp = var_imp))
}

# Process both models
model_10 <- process_model("trained_models_vip_8.RData", "_8")
model_50 <- process_model("trained_models_vip_50.RData", "_50")

# Create plots
pdf("vip.pdf", width = 10, height = 6)
par(mfrow = c(1, 2), mar = c(5, 8, 4, 2))

# Plot for n = 10
barplot(
  model_10$imp[model_10$ord],
  names.arg = rownames(model_10$var_imp$importance)[model_10$ord],
  horiz = TRUE,
  col = "steelblue",
  main = "Variable Importance\nRandom Forest (n = 8)",
  xlab = "Importance",
  las = 1,
  cex.names = 0.8
)

# Plot for n = 50
barplot(
  model_50$imp[model_50$ord],
  names.arg = rownames(model_50$var_imp$importance)[model_50$ord],
  horiz = TRUE,
  col = "steelblue",
  main = "Variable Importance\nRandom Forest (n = 50)",
  xlab = "Importance",
  las = 1,
  cex.names = 0.8
)

dev.off()
