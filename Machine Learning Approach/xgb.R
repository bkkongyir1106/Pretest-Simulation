# Install once if needed:
install.packages("xgboost")
install.packages("Matrix")      # for sparse model matrices
install.packages("caret")       # for train/test split & metrics

library(xgboost)
library(Matrix)
library(caret)

set.seed(42)

# ----- Data -----
data(iris)
iris$Species <- as.integer(iris$Species) - 1   # xgboost wants 0..K-1 labels

# Train/test split
idx <- createDataPartition(iris$Species, p = 0.8, list = FALSE)
train <- iris[idx, ]
test  <- iris[-idx, ]

# Model matrices (no intercept)
X_train <- model.matrix(~ . - 1 - Species, data = train)
X_test  <- model.matrix(~ . - 1 - Species, data = test)
y_train <- train$Species
y_test  <- test$Species

dtrain <- xgb.DMatrix(data = X_train, label = y_train)
dtest  <- xgb.DMatrix(data = X_test,  label = y_test)

# ----- Train with early stopping -----
params <- list(
  objective = "multi:softprob",  # probabilities per class
  num_class = length(unique(y_train)),
  eval_metric = "mlogloss",
  max_depth = 4,
  eta = 0.1,
  subsample = 0.8,
  colsample_bytree = 0.8,
  min_child_weight = 1
)

watchlist <- list(train = dtrain, eval = dtest)

xgb_model <- xgb.train(
  params = params,
  data = dtrain,
  nrounds = 500,
  watchlist = watchlist,
  early_stopping_rounds = 20,
  print_every_n = 20
)

# ----- Evaluate -----
# Predicted class = argmax over class probabilities
prob <- predict(xgb_model, dtest)                  # vector
prob <- matrix(prob, ncol = params$num_class, byrow = TRUE)
pred_class <- max.col(prob) - 1

confusionMatrix(
  factor(pred_class, levels = 0:2),
  factor(y_test,    levels = 0:2)
)

# ----- Feature importance -----
imp <- xgb.importance(model = xgb_model, feature_names = colnames(X_train))
print(imp)
xgb.plot.importance(imp)
