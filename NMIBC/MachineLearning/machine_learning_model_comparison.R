
################ Install Dependencies ################
setwd("D:\\ZYD\\R\\meta\\RF")
install.packages("mlr3extralearners")
BiocManager::install("distr6")
install.packages("mlr3tuning")
install.packages("mlr3pipelines")
install.packages("mlr3verse")
install.packages("compareGroups")
install.packages("abess")
install.packages("paradox")
install.packages("future")
install.packages("mlr3extralearners", repos = "https://cran.r-project.org")
install.packages("RWeka")
remotes::install_github("mlr-org/mlr3extralearners@*release")
lrn("regr.gbm")
install.packages("gitcreds")
library(gitcreds)
gitcreds_set()

library(mlr3)              # Core ML framework
library(mlr3learners)      # Base learners
library(mlr3filters)       # Feature filters
library(mlr3fselect)       # Feature selection
library(mlr3extralearners) # Extra learners
library(mlr3tuning)        # Hyperparameter tuning
library(mlr3pipelines)     # Pipelines
library(mlr3verse)         # ML universe
library(compareGroups)     # Baseline tables
library(tidyverse)
library(data.table)
library(cli)
library(parallelly)
library(Rcpp)
library(rlang)
library(stringi)
library(RWeka)
library(usethis)
library(xgboost)

install.packages("xgboost")
install.packages("distr6")
install.packages("mlr3proba")
remotes::install_github("mlr-org/mlr3proba")
library(mlr3proba)

# Other utilities
library(abess)   # Adaptive best subset selection
library(paradox) # Parameter space definition
library(future)  # Parallel computing support

data <- fread("./cohort1.txt", data.table = FALSE)

# Convert target variable 'type' to factor
data$type <- as.factor(data$type)

# Create classification task
task_train <- TaskClassif$new(
  id = 'buty',          # Task ID
  backend = data,       # Dataset
  target = 'type',      # Target variable
  positive = "BCA"      # Positive class
)

################ Explore and Filter Learners ################
# View all available learners
learners_all = as.data.table(list_mlr3learners())
# Can view more columns: select = c("id", "mlr3_package", "required_packages")

# Filter learners based on criteria
# [Adjustable] Filter criteria
learners_to_try <- learners_all %>%
  filter(class == "classif") %>%                 # Classification
  filter(grepl("twoclass", properties)) %>%      # Two-class support
  filter(grepl("integer", feature_types))        # Integer feature support

# Keep only name and ID
learners_to_try <- learners_to_try %>%
  dplyr::select(c(name, id)) %>%
  as.data.frame()

# [Adjustable] Select desired learner IDs
selected_ids <- c(
  "classif.AdaBoostM1",   # AdaBoost
  "classif.gbm",          # GBM
  "classif.log_reg",      # Logistic Regression
  "classif.randomForest", # Random Forest
  "classif.rpart",        # Decision Tree
  "classif.svm",          # SVM
  "classif.xgboost"       # XGBoost
)

# Final list of filtered learners
filtered_learners <- learners_to_try %>%
  filter(id %in% selected_ids)

################ Build Learner Objects ################
# Create learner objects for selected models
myvar <- list()

for (i in filtered_learners$name) {
  # Build learner, set predict type to prob
  myvar[[i]] <- lrn(
    learners_to_try[learners_to_try$name == i,][[2]],
    id = i,
    predict_type = "prob"
  )
}

# Save learner objects (optional)
# [Adjustable] Save path
save(myvar, file = "./out/all.learner.Rdata")
# Load saved learner objects
# load(file = "./out/all.learner.Rdata")

################ Cross-Validation Setup ################
# [Adjustable] CV folds
# 3-fold CV
mycv = rsmp("cv", folds = 3)
# 10-fold CV
# mycv = rsmp("cv", folds = 10)

# [Adjustable] Max evaluations for tuning
# Small value for quick testing
myterm_evals_test = 20
# Large value for formal tuning
myterm_evals = 2000

################ AdaBoostM1 Configuration ################
# Parameter space for tuning
AdaBoostM1_pars <- paradox::ps(
  P = paradox::p_int(lower = 90, upper = 100), # Confidence percentage
  S = paradox::p_int(lower = 1, upper = 5),    # Random seed
  I = paradox::p_int(lower = 1, upper = 15)    # Iterations
)

# Base learner config
AdaBoostM1_learner <- lrn("classif.AdaBoostM1", predict_type = "prob")

# Create AutoTuner
lrn_AdaBoostM1 = mlr3tuning::AutoTuner$new(
  learner = AdaBoostM1_learner,
  resampling = mycv,
  measure = msr("classif.auc"),                                  # Measure: AUC
  terminator = mlr3tuning::trm("evals", n_evals = myterm_evals), # Max evals
  tuner = mlr3tuning::tnr("random_search"),                      # Random search
  search_space = AdaBoostM1_pars,                                # Search space
  store_models = TRUE                                            # Store models
)

################ GBM Configuration ################
# Parameter space
gbm_pars <- paradox::ps(
  n.trees = paradox::p_int(lower = 50, upper = 800),     # Number of trees
  interaction.depth = paradox::p_int(lower = 1, upper = 5) # Interaction depth
)

# Base learner config
gbm_learner <- lrn("classif.gbm", predict_type = "prob")

# Create AutoTuner
lrn_gbm = mlr3tuning::AutoTuner$new(
  learner = gbm_learner,
  resampling = mycv,
  measure = msr("classif.auc"),
  terminator = mlr3tuning::trm("evals", n_evals = myterm_evals),
  tuner = mlr3tuning::tnr("random_search"),
  search_space = gbm_pars,
  store_models = TRUE
)

################ Logistic Regression Configuration ################
# Base learner config
log_learner <- lrn("classif.log_reg", predict_type = "prob")

# Create AutoTuner (empty search space, default params)
lrn_log = mlr3tuning::AutoTuner$new(
  learner = log_learner,
  resampling = mycv,
  measure = msr("classif.auc"),
  terminator = mlr3tuning::trm("evals", n_evals = myterm_evals),
  tuner = mlr3tuning::tnr("random_search"),
  search_space = paradox::ps(), # Empty search space
  store_models = TRUE
)

################ SVM Configuration ################
# Parameter space
svm_pars <- paradox::ps(
  gamma = paradox::p_dbl(lower = 0.1, upper = 10),                         # Gamma
  kernel = paradox::p_fct(levels = c("polynomial", "radial", "sigmoid")), # Kernel
  type = paradox::p_fct(levels = c("C-classification"))                   # SVM type
)

# Base learner config
svm_learner <- lrn("classif.svm", predict_type = "prob")

# Create AutoTuner
lrn_svm = mlr3tuning::AutoTuner$new(
  learner = svm_learner,
  resampling = mycv,
  measure = msr("classif.auc"),
  terminator = mlr3tuning::trm("evals", n_evals = myterm_evals),
  tuner = mlr3tuning::tnr("random_search"),
  search_space = svm_pars,
  store_models = TRUE
)

################ Random Forest Configuration ################
# Parameter space
rf_pars <- paradox::ps(
  ntree = paradox::p_int(lower = 10, upper = 500),   # Number of trees
  mtry = paradox::p_int(lower = 6, upper = 12),      # mtry
  nodesize = paradox::p_int(lower = 1, upper = 20),  # nodesize
  maxnodes = paradox::p_int(lower = 20, upper = 80)  # maxnodes
)

# Base learner config
rf_learner <- lrn("classif.randomForest", predict_type = "prob")

# Create AutoTuner
lrn_rf = mlr3tuning::AutoTuner$new(
  learner = rf_learner,
  resampling = mycv,
  measure = msr("classif.auc"),
  terminator = mlr3tuning::trm("evals", n_evals = myterm_evals),
  tuner = mlr3tuning::tnr("random_search"),
  search_space = rf_pars,
  store_models = TRUE
)

################ Decision Tree Configuration ################
# Parameter space
dt_pars <- paradox::ps(
  minsplit = paradox::p_int(lower = 1, upper = 20),   # minsplit
  maxdepth = paradox::p_int(lower = 3, upper = 10),   # maxdepth
  cp = paradox::p_dbl(lower = 0.001, upper = 0.1)     # complexity parameter (cp)
)

# Base learner config
dt_learner <- lrn("classif.rpart", predict_type = "prob")

# Create AutoTuner
lrn_dt = mlr3tuning::AutoTuner$new(
  learner = dt_learner,
  resampling = mycv,
  measure = msr("classif.auc"),
  terminator = mlr3tuning::trm("evals", n_evals = myterm_evals),
  tuner = mlr3tuning::tnr("random_search"),
  search_space = dt_pars,
  store_models = TRUE
)

################ XGBoost Configuration ################
# Parameter space
XGBoost_pars <- paradox::ps(
  eta = paradox::p_dbl(lower = 0, upper = 1),                          # Learning rate
  max_depth = paradox::p_int(lower = 3, upper = 20, default = 3),      # max_depth
  min_child_weight = paradox::p_int(lower = 1, upper = 10),            # min_child_weight
  subsample = paradox::p_dbl(lower = 0.3, upper = 1, default = 0.6),   # subsample
  colsample_bytree = paradox::p_dbl(lower = 0.3, upper = 1),           # colsample_bytree
  nrounds = paradox::p_int(lower = 1, upper = 100)                     # nrounds
)

# Base learner config
XGBoost_learner <- lrn("classif.xgboost", predict_type = "prob")

# Create AutoTuner
lrn_XGBoost = mlr3tuning::AutoTuner$new(
  learner = XGBoost_learner,
  resampling = mycv,
  measure = msr("classif.auc"),
  terminator = mlr3tuning::trm("evals", n_evals = myterm_evals),
  tuner = mlr3tuning::tnr("random_search"),
  search_space = XGBoost_pars,
  store_models = TRUE
)

################ Aggregate Learners ################
# Combine all configured learners into a list
all_learners <- list(
  lrn_AdaBoostM1,
  lrn_gbm,
  lrn_log,
  lrn_svm,
  lrn_rf,
  lrn_dt,
  lrn_XGBoost
)

################ Benchmarking and Model Comparison ################
# Load required packages
library(mlr3verse)
library(ggplot2)

# Set seed for reproducibility
set.seed(123)

# Define measures
measures <- list(
  msr("classif.auc"),         # AUC
  msr("classif.acc"),         # Accuracy
  msr("classif.sensitivity"), # Sensitivity
  msr("classif.specificity")  # Specificity
)

# Create benchmark design
design <- benchmark_grid(
  tasks = task_train,
  learners = all_learners,
  resamplings = rsmp("cv", folds = 3) # 3-fold CV
)

# Run benchmark (takes time)
bmr <- benchmark(design, store_models = TRUE)
save(bmr, file = "bmr2.Rdata")

# Load bmr.RData
load("bmr.RData")
tryCatch({
  load("bmr.RData")
}, error = function(e) {
  print(e)
})

# View aggregated results
bmr_results <- bmr$aggregate(measures)
print(bmr_results[, .(learner_id, classif.auc, classif.acc)])

# Assuming bmr_results is a list of dataframes to be saved
bmr_results_df <- do.call(rbind, bmr_results)
write.csv(bmr_results_df, "bmr_results2.csv", row.names = FALSE)

# Simple method: use mlr3viz directly
autoplot(bmr, type = "roc") +
  ggtitle("Models ROC Curves Comparison") +
  theme_minimal() +
  theme(legend.position = "bottom") +
  scale_color_discrete(name = "Model") +
  labs(x = "False Positive Rate", y = "True Positive Rate")

# Visualize model comparison (AUC)
autoplot(bmr, measure = msr("classif.auc")) +
  ggtitle("Models AUC Performance Comparison") +
  theme_minimal()

################ Best Model Training and Validation ################
# Select model with highest AUC
best_learner_id <- bmr_results[which.max(classif.auc), learner_id]
best_learner <- all_learners[[which(sapply(all_learners, function(x) x$id == best_learner_id))]]

# Train final model on full dataset
best_learner$train(task_train)

# Predict (using training data for demonstration)
prediction <- best_learner$predict(task_train)

# View performance metrics
cat("Final model performance on training set:\n")
prediction$score(measures)

# Plot ROC curve
autoplot(prediction, type = "roc") +
  ggtitle("ROC Curve") +
  theme_minimal()

install.packages("precrec")
library(precrec)

################ Feature Importance Analysis ################
if ("importance" %in% best_learner$properties) {
  # Extract feature importance
  importance <- best_learner$importance()

  # Visualization
  importance_df <- data.frame(
    Feature = names(importance),
    Importance = importance
  )

  ggplot(importance_df, aes(x = reorder(Feature, Importance), y = Importance)) +
    geom_col(fill = "steelblue") +
    coord_flip() +
    labs(x = "Metabolite Features", y = "Importance Score", title = "Feature Importance Analysis") +
    theme_minimal()
} else {
  message("Current model does not support feature importance analysis")
}

################ Model Saving and Deployment ################
# Save best model
saveRDS(best_learner, file = "./out/best_model.rds")
