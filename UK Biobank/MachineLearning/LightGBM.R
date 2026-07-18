
cat("Starting pipeline...\n")
required_packages <- c("tidyverse", "lightgbm", "caret", "pROC", "janitor", "isotone")
new_packages <- required_packages[!(required_packages %in% installed.packages()[,"Package"])]
if(length(new_packages)) {
  cat("Installing missing packages: ", paste(new_packages, collapse = ", "), "\n")
  install.packages(new_packages, repos = "https://cloud.r-project.org")
}

suppressPackageStartupMessages({
  library(tidyverse)
  library(lightgbm)
  library(caret)     
  library(pROC)      
  library(janitor)   
  library(isotone)   
})

# ==============================================================================
# 1 Helper Functions
# ==============================================================================


select_params_combo <- function(params_dict, nb_params = 100) {
  all_combos <- expand.grid(params_dict, stringsAsFactors = FALSE)
  
  if (nrow(all_combos) < nb_params) {
    return(all_combos)
  }
  

  selected_combos <- all_combos %>% sample_n(nb_params)
  
  return(selected_combos)
}

normal_imp <- function(imp_df) {
  if (nrow(imp_df) == 0) return(imp_df)
  
  total_gain <- sum(imp_df$Gain, na.rm = TRUE)
  if (total_gain == 0) {
    imp_df$Gain_Norm <- 0
  } else {
    imp_df$Gain_Norm <- imp_df$Gain / total_gain
  }
  return(imp_df)
}

clean_feature_names <- function(names) {
  names <- gsub("\\|", "_", names)
  names <- gsub(":", "_", names)
  names <- gsub(",", "_", names)
  names <- gsub("\"", "", names)
  names <- gsub("'", "", names)
  names <- gsub(" ", "_", names)
  names <- gsub("_+", "_", names)
  return(names)
}


get_cov_f_lst <- function(sex_id, all_covs) {

  if (!is.na(sex_id) && (sex_id == 1 || sex_id == 2)) {
    return(setdiff(all_covs, "Sex"))
  } else {
    return(all_covs)
  }
}

read_target <- function(disease_name, dpath) {

  fpath <- file.path(dpath, paste0(disease_name, ".csv"))
  if (!file.exists(fpath)) {
    fpath <- file.path(dpath, paste0(disease_name, "2.csv"))
  }
  
  if (!file.exists(fpath)) {
    stop(paste("Target file not found for:", disease_name))
  }
  
  df <- read_csv(fpath, show_col_types = FALSE)
  

  req_cols <- c("eid", "target_y", "BL2Target_yrs")
  if (!all(req_cols %in% names(df))) {
    stop(paste("File", fpath, "missing required columns:", paste(setdiff(req_cols, names(df)), collapse=", ")))
  }

  df_filtered <- df %>%
    filter(BL2Target_yrs > 0) %>%
    select(all_of(req_cols))
  
  return(df_filtered)
}

get_top_nmr_lst <- function(data_df, target_col, feature_cols, nmr_features, n_repeats = 10) {
  
  cat(sprintf("    [Feature Selection] Running %d iterations to select stable features...\n", n_repeats))
  
  agg_imp <- tibble(Feature = nmr_features, TotalGain = 0, Frequency = 0)
  
  params <- list(
    objective = "binary",
    metric = "auc",
    verbosity = -1,
    seed = 2024
  )
  
  for (i in 1:n_repeats) {
    set.seed(2024 + i)
    idx <- createDataPartition(data_df[[target_col]], p = 0.8, list = FALSE)
    sub_train <- data_df[idx, ]
    
    dtrain <- lgb.Dataset(
      data = as.matrix(sub_train[, feature_cols]),
      label = sub_train[[target_col]]
    )
    
    model <- lgb.train(
      params = params,
      data = dtrain,
      nrounds = 100,
      verbose = -1
    )
    
    imp <- lgb.importance(model, percentage = FALSE)
    
    if (nrow(imp) > 0) {
      agg_imp <- agg_imp %>%
        left_join(imp %>% select(Feature, Gain), by = "Feature") %>%
        mutate(
          Gain = replace_na(Gain, 0),
          TotalGain = TotalGain + Gain,
          Frequency = Frequency + ifelse(Gain > 0, 1, 0)
        ) %>%
        select(-Gain)
    }
    
    rm(dtrain, model, sub_train, idx)
    gc()
  }
  
  top_nmr <- agg_imp %>%
    mutate(AvgGain = TotalGain / n_repeats) %>%
    arrange(desc(AvgGain)) %>%
    head(30) %>%
    pull(Feature)
  
  return(top_nmr)
}

get_best_params <- function(train_df, features, target_col, candidate_params_df, n_folds = 3) {
  
  best_auc <- -1
  best_params <- list()
  
  folds <- createFolds(train_df[[target_col]], k = n_folds, list = TRUE)
  
  for (i in 1:nrow(candidate_params_df)) {
    current_params <- as.list(candidate_params_df[i, ])
    cat(sprintf("    CV Round: %d/%d | Leaves: %d | LR: %.3f\n", 
                i, nrow(candidate_params_df), current_params$num_leaves, current_params$learning_rate))
    
    params <- c(current_params, list(
      objective = "binary",
      metric = "auc",
      verbosity = -1,
      seed = 2024,
      num_threads = 1
    ))
    
    fold_aucs <- numeric(n_folds)
    best_iters <- numeric(n_folds)
    
    for (k in 1:n_folds) {
      gc() 
      val_idx <- folds[[k]]
      train_idx <- setdiff(seq_len(nrow(train_df)), val_idx)
      
      X_tr <- as.matrix(train_df[train_idx, features])
      y_tr <- train_df[[target_col]][train_idx]
      X_val <- as.matrix(train_df[val_idx, features])
      y_val <- train_df[[target_col]][val_idx]
      
      d_tr <- lgb.Dataset(data = X_tr, label = y_tr)
      d_val <- lgb.Dataset(data = X_val, label = y_val)
      
      m <- tryCatch({
        lgb.train(
          params = params,
          data = d_tr,
          nrounds = 1000,
          valids = list(val = d_val),
          early_stopping_rounds = 50,
          verbose = -1
        )
      }, error = function(e) {
        return(NULL)
      })
      
      if (!is.null(m)) {
        fold_aucs[k] <- m$best_score
        best_iters[k] <- m$best_iter
      } else {
        fold_aucs[k] <- 0
      }
      
      rm(d_tr, d_val, m, X_tr, X_val, y_tr, y_val)
      gc()
    }
    
    current_auc <- mean(fold_aucs)
    
    if (current_auc > best_auc) {
      best_auc <- current_auc
      best_params <- params
      best_params$nrounds <- ceiling(mean(best_iters))
    }
  }
  
  return(best_params)
}

run_cv_evaluation <- function(train_df, features, target_col, params) {
  
  folds <- createFolds(train_df[[target_col]], k = 5, list = TRUE)
  cv_aucs <- numeric(5)
  
  cat("    [Outer CV] Running 5-Fold Stratified CV on Training Set...\n")
  
  for (k in 1:5) {
    gc()
    val_idx <- folds[[k]]
    tr_idx <- setdiff(seq_len(nrow(train_df)), val_idx)
    
    d_tr <- lgb.Dataset(
      data = as.matrix(train_df[tr_idx, features]),
      label = train_df[[target_col]][tr_idx]
    )
    
    X_val <- as.matrix(train_df[val_idx, features])
    y_val <- train_df[[target_col]][val_idx]

    m <- lgb.train(
      params = params,
      data = d_tr,
      nrounds = if(!is.null(params$nrounds)) params$nrounds else 100,
      verbose = -1
    )
    
    preds <- predict(m, X_val)
    
    cv_aucs[k] <- as.numeric(pROC::roc(y_val, preds, quiet = TRUE)$auc)
  }
  
  mean_auc <- mean(cv_aucs)
  sd_auc <- sd(cv_aucs)
  cat(sprintf("      CV AUC: %.4f (+/- %.4f)\n", mean_auc, sd_auc))
  
  return(mean_auc)
}

model_train_pred <- function(train_df, test_df, features, target_col, best_params, n_calib_folds = 5) {
  
  X_train <- as.matrix(train_df[, features])
  y_train <- train_df[[target_col]]
  X_test <- as.matrix(test_df[, features])
  
  dtrain <- lgb.Dataset(data = X_train, label = y_train)
  
  if (n_calib_folds < 2) {
      cat("      [Info] Skipping calibration (n_calib_folds < 2).\n")
      final_model <- lgb.train(
        params = best_params,
        data = dtrain,
        nrounds = if(!is.null(best_params$nrounds)) best_params$nrounds else 100,
        verbose = -1
      )
      
      raw_preds <- predict(final_model, X_test)
      calibrated_preds <- raw_preds

      imp_df <- lgb.importance(final_model, percentage = FALSE)
      
      return(list(
        y_pred = calibrated_preds,
        y_raw = raw_preds,
        importance = imp_df
      ))
  }

  folds <- createFolds(y_train, k = n_calib_folds, list = TRUE)
  oof_preds <- numeric(length(y_train))
  
  for (f in 1:n_calib_folds) {
    gc()
    idx_val <- folds[[f]]
    idx_tr <- setdiff(seq_along(y_train), idx_val)
    
    X_tr_fold <- X_train[idx_tr, , drop=FALSE]
    y_tr_fold <- y_train[idx_tr]
    
    d_tr_fold <- lgb.Dataset(data = X_tr_fold, label = y_tr_fold)
    
    m_fold <- tryCatch({
        lgb.train(
          params = best_params,
          data = d_tr_fold,
          nrounds = if(!is.null(best_params$nrounds)) best_params$nrounds else 100,
          verbose = -1
        )
    }, error = function(e) {
        cat(paste("      [Error] OOF Fold", f, "failed:", e$message, "\n"))
        return(NULL)
    })
    
    if (!is.null(m_fold)) {
        oof_preds[idx_val] <- predict(m_fold, X_train[idx_val, , drop=FALSE])
    }
    
    rm(d_tr_fold, m_fold, X_tr_fold, y_tr_fold)
    gc()
  }
  
  oof_auc <- tryCatch({
    as.numeric(pROC::roc(y_train, oof_preds, quiet = TRUE)$auc)
  }, error = function(e) 0.5)
  
  cat(sprintf("      [Calibration] OOF AUC: %.4f\n", oof_auc))
  
  iso_reg <- isoreg(x = oof_preds, y = y_train)
  
  calibrate_prob <- as.stepfun(iso_reg)
  
  final_model <- lgb.train(
    params = best_params,
    data = dtrain,
    nrounds = if(!is.null(best_params$nrounds)) best_params$nrounds else 100,
    verbose = -1
  )
  
  raw_preds <- predict(final_model, X_test)
  
  cat(sprintf("      [Debug] Raw Preds: Mean=%.4f, SD=%.4f, Min=%.4f, Max=%.4f\n", 
              mean(raw_preds), sd(raw_preds), min(raw_preds), max(raw_preds)))
  
  calibrated_preds <- calibrate_prob(raw_preds)
  
  calibrated_preds <- pmin(pmax(calibrated_preds, 0), 1)
  
  cat(sprintf("      [Debug] Calibrated Preds: Mean=%.4f, SD=%.4f\n", 
              mean(calibrated_preds), sd(calibrated_preds)))
  
  if (var(calibrated_preds) < 1e-6) {
    cat("      [Warning] Calibrated predictions have near-zero variance. Check model performance.\n")
  }
  
  imp_df <- lgb.importance(final_model, percentage = FALSE)
  
  cat("      [Explanation] Calculating SHAP values...\n")
  # Fix: predcontrib is deprecated, use type='contrib'
  shap_mat <- predict(final_model, X_test, type = "contrib")
  
  return(list(
    y_pred = calibrated_preds,
    y_raw = raw_preds,
    importance = imp_df,
    shap_values = shap_mat,
    test_eid = test_df$eid 
  ))
}

evaluate_performance_bootstrap <- function(y_true, y_pred, n_boot = 1000) {
  
  cat(sprintf("    [Eval] Running Bootstrap (n=%d) for CI...\n", n_boot))
  
  aucs <- numeric(n_boot)
  n <- length(y_true)
  
  set.seed(2024)
  for (i in 1:n_boot) {
    idx <- sample(n, n, replace = TRUE)
    
    aucs[i] <- tryCatch({
      as.numeric(pROC::roc(y_true[idx], y_pred[idx], quiet = TRUE)$auc)
    }, error = function(e) NA)
  }
  
  aucs <- aucs[!is.na(aucs)]
  
  if (length(aucs) == 0) return(list(median = NA, lower = NA, upper = NA))
  
  res <- list(
    median = median(aucs),
    lower = quantile(aucs, 0.025),
    upper = quantile(aucs, 0.975)
  )
  
  cat(sprintf("      AUC: %.3f (95%% CI: %.3f - %.3f)\n", res$median, res$lower, res$upper))
  return(res)
}

analyze_shap_values <- function(shap_mat, feature_df) {
  common_feats <- intersect(colnames(shap_mat), colnames(feature_df))
  if (length(common_feats) == 0) return(NULL)
  
  shap_mat <- shap_mat[, common_feats, drop=FALSE]
  feature_df <- feature_df[, common_feats, drop=FALSE]
  
  mean_abs_shap <- colMeans(abs(shap_mat))
  
  total_shap <- sum(mean_abs_shap)
  norm_shap <- if (total_shap > 0) mean_abs_shap / total_shap else 0

  cor_vals <- numeric(length(common_feats))
  for (i in seq_along(common_feats)) {
    feat <- common_feats[i]
    if (sd(feature_df[[feat]]) > 0 && sd(shap_mat[, feat]) > 0) {
      cor_vals[i] <- cor(feature_df[[feat]], shap_mat[, feat])
    } else {
      cor_vals[i] <- 0
    }
  }
  
  res_df <- data.frame(
    Feature = common_feats,
    MeanAbsSHAP = mean_abs_shap,
    NormalizedSHAP = norm_shap,
    DirectionCor = cor_vals
  ) %>%
    arrange(desc(MeanAbsSHAP))
  
  return(res_df)
}

visualize_results <- function(pred_df, imp_df, dpath) {
  
  plot_dir <- dpath
  # if (!dir.exists(plot_dir)) dir.create(plot_dir)
  
  diseases <- unique(pred_df$disease)
  
  model_colors <- c(
    "Covariates" = "#1f77b4",       # Blue
    "NMR" = "#2ca02c",              # Green
    "NMR+Covariates" = "#d62728"    # Red
  )
  
  for (dz in diseases) {
    cat(paste0("  Plotting for disease: ", dz, "\n"))
    
    # --- 1. ROC Curve with CI Ribbon ---
    dz_preds <- pred_df %>% filter(disease == dz)
    models <- unique(dz_preds$model)
    
    plot_data <- data.frame()
    label_map <- c()
    
    for (m in models) {
      m_data <- dz_preds %>% filter(model == m)
      
      if (nrow(m_data) < 10) next
      
      r <- pROC::roc(m_data$target_y, m_data$y_raw, quiet = TRUE)
    
      cat(paste0("    Calculating CI for ", m, "...\n"))
      ci_res <- pROC::ci.se(r, specificities = seq(0, 1, 0.01), conf.level = 0.95, progress = "none")
      
      specs <- as.numeric(colnames(ci_res))
      lowers <- ci_res[1, ]
      medians <- ci_res[2, ]
      uppers <- ci_res[3, ]
      
      tmp_df <- data.frame(
        specificity = specs,
        sensitivity = medians,
        ymin = lowers,
        ymax = uppers,
        model = m
      )
      plot_data <- bind_rows(plot_data, tmp_df)
      
      auc_med <- m_data$auc_median[1]
      auc_low <- m_data$auc_lower[1]
      auc_up  <- m_data$auc_upper[1]
      
      label_str <- sprintf("%s   %.3f (%.3f-%.3f)", m, auc_med, auc_low, auc_up)
      label_map[m] <- label_str
    }
    
    if (nrow(plot_data) > 0) {
      p_roc <- ggplot(plot_data, aes(x = 1 - specificity, y = sensitivity, group = model)) +
        geom_ribbon(aes(ymin = ymin, ymax = ymax, fill = model), alpha = 0.2) +
        geom_line(aes(color = model), linewidth = 1) +
        geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#d62728") +
        scale_x_continuous(name = "False positive rate", limits = c(0, 1), expand = c(0, 0)) +
        scale_y_continuous(name = "True positive rate", limits = c(0, 1), expand = c(0, 0)) +
        scale_color_manual(values = model_colors, labels = label_map) +
        scale_fill_manual(values = model_colors, labels = label_map) +
        theme_bw() +
        theme(
          panel.grid.major = element_line(linetype = "dashed", color = "gray80"),
          panel.grid.minor = element_blank(),
          legend.position = c(0.95, 0.05), 
          legend.justification = c(1, 0),
          legend.background = element_rect(fill = "gray95", color = NA),
          legend.title = element_blank(),
          legend.text = element_text(size = 9),
          legend.key = element_blank(),
          plot.title = element_text(hjust = 0.5, face = "bold")
        ) +
        labs(title = dz)
      
      ggsave(file.path(plot_dir, paste0(dz, "_ROC_Ribbon.png")), p_roc, width = 6, height = 6)
    }
    
    # --- 2. Feature Importance (Top 20) ---
    dz_imp <- imp_df %>% 
      filter(disease == dz, model == "NMR+Covariates")
    
    if (nrow(dz_imp) > 0) {
      top_feats <- dz_imp %>%
        group_by(Feature) %>%
        summarise(Gain = mean(Gain_Norm), .groups = "drop") %>%
        arrange(desc(Gain)) %>%
        head(20)
      
      p_imp <- ggplot(top_feats, aes(x = reorder(Feature, Gain), y = Gain)) +
        geom_bar(stat = "identity", fill = "steelblue") +
        coord_flip() +
        labs(title = paste("Top 20 Feature Importance -", dz),
             subtitle = "Model: NMR+Covariates",
             x = "Feature",
             y = "Normalized Gain") +
        theme_minimal()
      
      ggsave(file.path(plot_dir, paste0(dz, "_Importance.png")), p_imp, width = 8, height = 10)
    }
  }
}

# ==============================================================================

main <- function() {
  nb_cpus <- 8
  my_seed <- 2024
  
  N_FEAT_SELECT_ITERS <- 2 
  INNER_CV_FOLDS <- 2     
  N_PARAM_COMBOS <- 1      
  CALIB_FOLDS <- 3       
  N_BOOTSTRAP <- 200      
  
  set.seed(my_seed)
  dpath <- "d:/ZYD/R/meta/biobank/machine/prevelant" 
  
  params_dict <- list(
    learning_rate = c(0.01, 0.05, 0.1),
    num_leaves = c(31, 63, 127),
    max_depth = c(-1, 7, 10),
    feature_fraction = c(0.8, 0.9, 1.0),
    bagging_fraction = c(0.8, 0.9, 1.0),
    bagging_freq = c(1, 5),
    lambda_l1 = c(0, 0.1, 1, 10),
    lambda_l2 = c(0, 0.1, 1, 10)
  )
  
  candidate_params_lst <- select_params_combo(params_dict, nb_params = N_PARAM_COMBOS)
  
  cat("Loading data...\n")

  disease_list_df <- read_csv(file.path(dpath, "IncidentDiseaseTable.csv"), show_col_types = FALSE)
  
  cov_df <- read_csv(file.path(dpath, "Covariates.csv"), show_col_types = FALSE)

  if ("target_y" %in% names(cov_df)) {
    cov_df <- cov_df %>% select(-target_y)
  }
  
  names(cov_df) <- clean_feature_names(names(cov_df))
  
  nmr_df <- read_csv(file.path(dpath, "NMR_Preprocessed.csv"), show_col_types = FALSE)
  
  if ("Participant.ID" %in% names(nmr_df)) {
    nmr_df <- nmr_df %>% rename(eid = Participant.ID)
  }
  
  names(nmr_df) <- clean_feature_names(names(nmr_df))
  
  if (any(duplicated(names(cov_df)))) {
    cat("Warning: Duplicated column names in Covariates after cleaning. Fixing...\n")
    names(cov_df) <- make.unique(names(cov_df))
  }
  if (any(duplicated(names(nmr_df)))) {
    cat("Warning: Duplicated column names in NMR after cleaning. Fixing...\n")
    names(nmr_df) <- make.unique(names(nmr_df))
  }
  
  nmr_cov_df <- inner_join(nmr_df, cov_df, by = "eid")
  
  all_cols <- names(nmr_cov_df)
  exclude_cols <- c("eid") # 排除 ID 列
  nmr_features <- setdiff(names(nmr_df), "eid")
  cov_features <- setdiff(names(cov_df), "eid")
  
  final_pred_results <- list()
  final_imp_results <- list()
  
  diseases <- disease_list_df$NAME
  
  for (tgt_disease in diseases) {
    cat(paste0("\nProcessing disease: ", tgt_disease, "\n"))
    
    disease_info <- disease_list_df %>% filter(NAME == tgt_disease)
    sex_id <- disease_info$SEX[1]
    
    target_df <- tryCatch({
      read_target(tgt_disease, dpath)
    }, error = function(e) {
      cat(paste("  Error reading target:", e$message, "\n"))
      return(NULL)
    })
    
    if (is.null(target_df)) next
    
    if ("target_y" %in% names(nmr_cov_df)) {
      nmr_cov_df <- nmr_cov_df %>% select(-target_y)
    }
    
    full_data <- inner_join(nmr_cov_df, target_df, by = "eid")
    
    if (!is.na(sex_id)) {
      
      if (sex_id == 1) { # Male only
        cat("  Filtering for Males (Sex=1)...\n")
        full_data <- full_data %>% filter(Sex == 1)
      } else if (sex_id == 2) { # Female only
        cat("  Filtering for Females (Sex=0)...\n") 
        full_data <- full_data %>% filter(Sex == 0)
      } else {
        cat("  No sex filtering (Analyzing both sexes)...\n")
      }
    }
    
    if (nrow(full_data) < 50) {
      cat("  Not enough samples after filtering. Skipping.\n")
      next
    }
    
    all_feature_cols <- c(nmr_features, cov_features)
    all_feature_cols <- intersect(all_feature_cols, names(full_data))
    
    non_num_cols <- all_feature_cols[!sapply(full_data[, all_feature_cols], is.numeric)]
    if (length(non_num_cols) > 0) {
      cat(paste("  Converting non-numeric features to numeric:", paste(non_num_cols, collapse=", "), "\n"))
      for (col in non_num_cols) {
        full_data[[col]] <- as.numeric(as.factor(full_data[[col]]))
      }
    }
    
    cat(paste("  Rows in full_data:", nrow(full_data), "\n"))
    
    curr_cov_lst <- get_cov_f_lst(sex_id, cov_features)
    all_features <- colnames(full_data)[(grep("^eid$|^target_y$|^Sex$", colnames(full_data), invert = TRUE))]
    all_features <- setdiff(all_features, curr_cov_lst)
    
    nmr_features <- all_features[grep("_0$", all_features)]
    
    cat(paste("  Total NMR features (Instance 0 only):", length(nmr_features), "\n"))

    if (!"target_y" %in% names(full_data)) {
      cat("  [DEBUG] target_y is MISSING from full_data!\n")
      cat("  [DEBUG] Columns in nmr_cov_df:", paste(head(names(nmr_cov_df), 10), collapse=", "), "...\n")
      cat("  [DEBUG] Columns in target_df:", paste(names(target_df), collapse=", "), "\n")
      stop("Critical error: target_y missing")
    } else {
      cat("  [DEBUG] target_y is present in full_data.\n")
    }
    
    if(length(nmr_features) == 0) {
      stop("No Instance 0 features found! Please check feature names.")
    }
    
    # 5. 10-Fold External Cross-Validation
    outer_folds <- createFolds(full_data$target_y, k = 10, list = TRUE)
    
    fold_pred_results <- list()
    fold_imp_results <- list()
    fold_shap_results <- list()
    
    (Feature -> Count)
    feature_stability_counter <- list()
    
    pb <- txtProgressBar(min = 0, max = 10, style = 3)
    
    for (f_idx in 1:10) {
      setTxtProgressBar(pb, f_idx)
      cat(sprintf("\n  [Outer Fold %d/10] Processing...\n", f_idx))
      
      test_idx <- outer_folds[[f_idx]]
      train_df <- full_data[-test_idx, ]
      test_df <- full_data[test_idx, ]
    
      model_configs <- list(
        list(name = "NMR", feats = nmr_features),
        list(name = "Covariates", feats = curr_cov_lst),
        list(name = "NMR+Covariates", feats = c(nmr_features, curr_cov_lst))
      )
      
      top_nmr_30 <- tryCatch({
        get_top_nmr_lst(train_df, "target_y", c(nmr_features, curr_cov_lst), nmr_features, n_repeats = N_FEAT_SELECT_ITERS)
      }, error = function(e) {
        cat(paste("    [Error] Feature selection failed:", e$message, "\n"))
        return(NULL)
      })
      
      if (is.null(top_nmr_30)) {
          cat("    Skipping fold due to feature selection failure.\n")
          next
      }
      
      for (ft in top_nmr_30) {
        if (is.null(feature_stability_counter[[ft]])) {
          feature_stability_counter[[ft]] <- 1
        } else {
          feature_stability_counter[[ft]] <- feature_stability_counter[[ft]] + 1
        }
      }
      
      model_configs[[1]]$feats <- top_nmr_30
      model_configs[[3]]$feats <- c(top_nmr_30, curr_cov_lst)
      
      for (config in model_configs) {
        model_name <- config$name
        model_feats <- config$feats
        
        # cat(paste("    Training model:", model_name, "\n"))
        
        if (length(model_feats) == 0) next
        
        best_params <- tryCatch({
            get_best_params(train_df, model_feats, "target_y", candidate_params_lst, n_folds = INNER_CV_FOLDS)
        }, error = function(e) {
            cat(paste("    [Error] Tuning failed for", model_name, ":", e$message, "\n"))
            return(NULL)
        })
        
        if (is.null(best_params)) next
        
        res <- tryCatch({
          model_train_pred(train_df, test_df, model_feats, "target_y", best_params, n_calib_folds = CALIB_FOLDS)
        }, error = function(e) {
          cat(paste("    [Error] Fold", f_idx, "Model", model_name, "failed:", e$message, "\n"))
          return(NULL)
        })
        
        if (is.null(res)) next
      
        pred_res <- test_df %>%
          select(eid, target_y) %>%
          mutate(
            y_pred = res$y_pred,
            y_raw = res$y_raw,
            model = model_name,
            disease = tgt_disease,
            fold = f_idx
          )
        fold_pred_results[[length(fold_pred_results) + 1]] <- pred_res

        imp_res <- res$importance %>%
          mutate(
            model = model_name,
            disease = tgt_disease,
            fold = f_idx
          ) %>%
          normal_imp()
        
        fold_imp_results[[length(fold_imp_results) + 1]] <- imp_res
        

        if (model_name == "NMR+Covariates") {
             fold_shap_results[[length(fold_shap_results) + 1]] <- list(
               shap = res$shap_values,
               data = test_df[, model_feats],
               eid = res$test_eid,
               fold = f_idx
             )
        }
      } # End Model Loop
      
      rm(train_df, test_df)
      gc()
    } # End Fold Loop
    
    close(pb)
    
    final_pred_results <- bind_rows(fold_pred_results)
    
    models_list <- unique(final_pred_results$model)
    boot_results_list <- list()
    
    for (m in models_list) {
      cat(paste0("    [Eval] Calculating Bootstrap AUC for model: ", m, "\n"))
      m_data <- final_pred_results %>% filter(model == m)
      
      boot_res <- evaluate_performance_bootstrap(m_data$target_y, m_data$y_raw, n_boot = 1000)
      
      boot_results_list[[m]] <- data.frame(
        model = m,
        auc_median = boot_res$median,
        auc_lower = boot_res$lower,
        auc_upper = boot_res$upper
      )
    }
    
    boot_results_df <- bind_rows(boot_results_list)
    
    final_pred_results <- final_pred_results %>%
      left_join(boot_results_df, by = "model")
  
    final_imp_results <- bind_rows(fold_imp_results)
    stability_df <- tibble(
      Feature = names(feature_stability_counter),
      SelectionCount = unlist(feature_stability_counter),
      Frequency = SelectionCount / 10
    ) %>%
      arrange(desc(Frequency))
    
    write_csv(stability_df, file.path(dpath, paste0("feature_stability_", tgt_disease, ".csv")))
    cat(paste("  Feature stability saved to feature_stability_", tgt_disease, ".csv\n", sep=""))


    if (length(fold_shap_results) > 0) {
      all_shap_mat <- do.call(rbind, lapply(fold_shap_results, function(x) x$shap))
      all_shap_data <- do.call(bind_rows, lapply(fold_shap_results, function(x) as.data.frame(x$data)))
      
      shap_file <- file.path(dpath, paste0("shap_values_", tgt_disease, "_10fold.rds"))
      saveRDS(list(shap = all_shap_mat, data = all_shap_data), shap_file)
      
      shap_summary <- analyze_shap_values(all_shap_mat, all_shap_data)
      if (!is.null(shap_summary)) {
          shap_sum_file <- file.path(dpath, paste0("shap_summary_", tgt_disease, ".csv"))
          write_csv(shap_summary, shap_sum_file)
      }
    }

    write_csv(final_pred_results, file.path(dpath, paste0("final_predictions_10fold_", tgt_disease, ".csv")))
    write_csv(final_imp_results, file.path(dpath, paste0("final_feature_importance_10fold_", tgt_disease, ".csv")))
    
    cat("  Generating plots for this disease...\n")
    visualize_results(final_pred_results, final_imp_results, dpath)
    
  } # End Disease Loop
  
  cat("\nAll diseases processed.\n")
}

main()
