# ==============================================================================
# 0. Environment Setup and Data Loading
# ==============================================================================
library(data.table)
library(caret)
library(glmnet)
library(Boruta)
library(randomForest)
library(pROC)
library(ggplot2)
library(ggpubr)
library(ggprism)
library(reshape2)
library(dplyr)
# ---------------- Load dataset ----------------
data <- fread("./tree.txt", data.table = FALSE)
# First column: outcome variable
colnames(data)[1] <- "type"
data$type <- as.factor(data$type)
class_levels <- levels(data$type)
neg_class <- class_levels[1]
pos_class <- class_levels[2]

cat(sprintf(
  "Dataset loaded: %d samples, %d metabolite features\n",
  nrow(data),
  ncol(data)-1
))
# ==============================================================================
# 1. Stratified Training/Test Split
# ==============================================================================

set.seed(123)
# 70% training, 30% held-out test
train_idx <- createDataPartition(
  data$type,
  p = 0.7,
  list = FALSE
)

train_data <- data[train_idx, ]
test_data  <- data[-train_idx, ]
cat("Training samples:", nrow(train_data), "\n")
cat("Held-out test samples:", nrow(test_data), "\n")
# ==============================================================================
# 2. Consensus Feature Selection
#    LASSO + SVM-RFE + Boruta
#    Training Set Only
# ==============================================================================
features_x <- as.matrix(train_data[, -1])
labels_y <- train_data$type

# 10 resampled training subsets
set.seed(123)

cv_folds <- createFolds(
  labels_y,
  k = 10,
  returnTrain = TRUE
)

feature_counts <- setNames(
  rep(0, ncol(features_x)),
  colnames(features_x)
)

cat("Running consensus feature selection...\n")

for(i in seq_along(cv_folds)){

  cat("Iteration:", i, "/10\n")
  subset_x <- features_x[cv_folds[[i]], ]
  subset_y <- labels_y[cv_folds[[i]]]

  # -----------------------------
  # A. LASSO regression
  # -----------------------------
  cv_lasso <- cv.glmnet(
    subset_x,
    subset_y,
    family = "binomial",
    alpha = 1
  )
  lasso_coef <- coef(
    cv_lasso,
    s = "lambda.min"
  )
  lasso_selected <- rownames(lasso_coef)[
    which(lasso_coef != 0)
  ]
  # remove intercept
  lasso_selected <- setdiff(
    lasso_selected,
    "(Intercept)"
  )

  # -----------------------------
  # B. SVM-RFE
  # -----------------------------

  svm_control <- rfeControl(
    functions = caretFuncs,
    method = "cv",
    number = 5
  )
  svm_rfe <- rfe(
    subset_x,
    subset_y,
    sizes = c(5,10,15,20),
    rfeControl = svm_control,
    method = "svmLinear"
  )
  svm_selected <- predictors(svm_rfe)

  # -----------------------------
  # C. Boruta
  # -----------------------------
  boruta_result <- Boruta(
    x = subset_x,
    y = subset_y,
    maxRuns = 300,
    doTrace = 0
  )
  boruta_selected <-
    getSelectedAttributes(
      boruta_result,
      withTentative = FALSE
    )
  # -----------------------------
  # Consensus features
  # -----------------------------
  consensus_features <-
    Reduce(
      intersect,
      list(
        lasso_selected,
        svm_selected,
        boruta_selected
      )
    )
  if(length(consensus_features)>0){

    feature_counts[consensus_features] <-
      feature_counts[consensus_features]+1
  }
}
# Retain metabolites selected in >=80% iterations
final_features <- names(
  feature_counts[
    feature_counts >= 8
  ]
)
if(length(final_features)==0){
  stop(
    "No metabolites fulfilled the 80% consensus criterion."
  )
}
cat(
  "Selected features:",
  paste(final_features, collapse=", "),
  "\n"
)
# ==============================================================================
# 3. Random Forest Model Training
# ==============================================================================
train_subset <- train_data[
  ,
  c("type", final_features)
]
# mtry determined by square-root rule
mtry_value <-
  max(
    1,
    floor(
      sqrt(length(final_features))
    )
  )
rf_control <- trainControl(
  method = "cv",
  number = 5,
  classProbs = TRUE
)
set.seed(123)
rf_model <- train(
  type ~ .,
  data = train_subset,
  method = "rf",
  tuneGrid = data.frame(
    mtry = mtry_value
  ),
  ntree = 500,
  nodesize = 5,
  trControl = rf_control,
  importance = TRUE
)
final_model <- rf_model$finalModel
# ==============================================================================
# 4. Held-out Test Set Evaluation
# ==============================================================================
test_prob <-
  predict(
    rf_model,
    newdata = test_data,
    type="prob"
  )[,pos_class]
test_class <-
  predict(
    rf_model,
    newdata=test_data
  )
# Confusion matrix
cm <- confusionMatrix(
  test_class,
  test_data$type,
  positive = pos_class
)
# ROC
roc_obj <- roc(
  test_data$type,
  test_prob,
  levels=c(
    neg_class,
    pos_class
  )
)
# Bootstrap CI
auc_ci <-
  ci.auc(
    roc_obj,
    method="bootstrap",
    boot.n=1000
  )
cat("\n========== Held-out Test Performance ==========\n")
cat(
  sprintf(
    "AUC %.3f (95%% CI %.3f-%.3f)\n",
    auc_ci[2],
    auc_ci[1],
    auc_ci[3]
  )
)
cat(
  sprintf(
    "Accuracy %.3f\n",
    cm$overall["Accuracy"]
  )
)
cat(
  sprintf(
    "Sensitivity %.3f\n",
    cm$byClass["Sensitivity"]
  )
)

cat(
  sprintf(
    "Specificity %.3f\n",
    cm$byClass["Specificity"]
  )
)
# ==============================================================================
# 5. Feature Importance
# ==============================================================================

importance_df <-
  as.data.frame(
    importance(final_model)
  )

importance_df$Feature <-
  rownames(importance_df)

importance_df <-
  importance_df[
    order(
      importance_df$MeanDecreaseGini,
      decreasing=TRUE
    ),
  ]
ggplot(
  importance_df,
  aes(
    x=reorder(
      Feature,
      MeanDecreaseGini
    ),
    y=MeanDecreaseGini
  )
)+
geom_col()+
coord_flip()+
theme_prism()+
labs(
  x="Metabolites",
  y="Mean Decrease Gini"
)
# ==============================================================================
# Save model
# ==============================================================================

saveRDS(
  rf_model,
  "NMIBC_metabolomic_classification_RF_model.rds"
)
