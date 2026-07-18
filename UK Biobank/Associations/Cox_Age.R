# Load necessary libraries
options(warn = -1)
suppressPackageStartupMessages({
  library(survival)
  library(parallel)
})

# Function to fit model for one metabolite
fit_cox <- function(metabolite, data, covs) {
  # Extract x_nmr
  x_nmr <- data[[metabolite]]
  
  # Remove NAs in x_nmr
  idx <- !is.na(x_nmr)
  sub_data <- data[idx, ]
  sub_x <- x_nmr[idx]
  
  # Check events >= 50
  if (sum(sub_data$target_y) < 50) {
    return(NULL)
  }
  
  # Construct formula
  # We use data frame subsetting or direct variable passing to avoid environment issues in parallel
  # But coxph needs a data frame or environment.
  # Let's create a temporary dataframe for the model
  model_df <- sub_data[, c("BL2Target_yrs", "target_y", covs)]
  model_df$x_nmr <- sub_x
  
  # Formula string
  f_str <- paste("Surv(BL2Target_yrs, target_y) ~", paste(c(covs, "x_nmr"), collapse = " + "))
  
  # Fit model with retry logic for stability
  # Try standard fit with Efron ties (default but explicit)
  fit <- tryCatch({
    survival::coxph(as.formula(f_str), data = model_df, ties = "efron")
  }, error = function(e) return(NULL), warning = function(w) return(NULL))
  
  # Retry with more iterations if failed
  if (is.null(fit)) {
    fit <- tryCatch({
      ctl <- survival::coxph.control(iter.max = 100, outer.max = 20)
      survival::coxph(as.formula(f_str), data = model_df, ties = "efron", control = ctl)
    }, error = function(e) return(NULL))
  }
  
  if (is.null(fit)) return(NULL)
  
  # Extract results for x_nmr
  # summary(fit)$coefficients gives matrix with cols: coef, exp(coef), se(coef), z, Pr(>|z|)
  summ <- summary(fit)
  coefs <- summ$coefficients
  cis <- summ$conf.int
  
  # Find x_nmr row
  row_idx <- grep("x_nmr", rownames(coefs))
  if (length(row_idx) == 0) return(NULL)
  
  hr <- coefs[row_idx, "exp(coef)"]
  pval <- coefs[row_idx, "Pr(>|z|)"]
  # conf.int columns: exp(coef), exp(-coef), lower .95, upper .95
  # indices might vary by version, usually 3 and 4 are CIs
  lower <- cis[row_idx, 3]
  upper <- cis[row_idx, 4]
  
  return(c(
    metabolite = metabolite,
    hr = hr,
    lower = lower,
    upper = upper,
    pval = pval,
    n = fit$n,
    nevent = fit$nevent
  ))
}

run_analysis_pipeline <- function(df, label, out_file, nmr_cols, covariates) {
  cat(sprintf("\n--- Running Analysis for %s ---\n", label))
  
  # Filter: BL2Target_yrs > 0
  n_before <- nrow(df)
  df <- df[df$BL2Target_yrs > 0, ]
  cat(sprintf("Filtered BL2Target_yrs > 0: %d -> %d\n", n_before, nrow(df)))
  
  if (nrow(df) == 0) {
    cat("No samples left after filtering. Skipping.\n")
    return()
  }
  
  # Subset DF for Parallel Efficiency
  all_cols_needed <- unique(c("eid", "BL2Target_yrs", "target_y", covariates, nmr_cols))
  df <- df[, intersect(all_cols_needed, names(df))]
  cat(sprintf("Subsetted dataframe to %d columns for efficiency.\n", ncol(df)))
  
  # Run in parallel
  n_cores <- parallel::detectCores(logical = FALSE)
  n_cores <- max(1, n_cores - 1)
  # n_cores <- 1 # Force sequential for debugging
  
  cat(sprintf("Running with %d cores...\n", n_cores))
  
  if (n_cores > 1) {
    cl <- parallel::makeCluster(n_cores)
    parallel::clusterEvalQ(cl, library(survival))
    results_list <- parallel::parLapply(cl, nmr_cols, fit_cox, data = df, covs = covariates)
    parallel::stopCluster(cl)
  } else {
    results_list <- lapply(nmr_cols, function(m) fit_cox(m, df, covariates))
  }
  
  # Combine results
  results_list <- results_list[!sapply(results_list, is.null)]
  
  if (length(results_list) == 0) {
    cat("No successful models found.\n")
    return()
  }
  
  # Convert to data frame
  res_df <- do.call(rbind, results_list)
  res_df <- as.data.frame(res_df, stringsAsFactors = FALSE)
  
  # Convert numeric columns
  num_cols <- c("hr", "lower", "upper", "pval", "n", "nevent")
  res_df[num_cols] <- lapply(res_df[num_cols], as.numeric)
  
  # Multiple Testing Correction
  res_df$pval_bfi <- p.adjust(res_df$pval, method = "bonferroni")
  
  # Formatting
  format_res <- function(r) {
    stars <- ""
    pval <- as.numeric(r['pval_bfi'])
    if (pval < 0.001) stars <- "***"
    else if (pval < 0.01) stars <- "**"
    else if (pval < 0.05) stars <- "*"
    
    sprintf("%.2f [%.2f-%.2f]%s", 
            as.numeric(r['hr']), 
            as.numeric(r['lower']), 
            as.numeric(r['upper']), 
            stars)
  }
  
  res_df$Formatted_Result <- apply(res_df, 1, format_res)
  
  # Clean metabolite names
  res_df$Metabolite_Clean <- gsub(" \\| Instance 0", "", res_df$metabolite)
  
  # Select columns
  final_df <- res_df[, c("metabolite", "Metabolite_Clean", "n", "nevent", "hr", "lower", "upper", "pval", "pval_bfi", "Formatted_Result")]
  colnames(final_df) <- c("Metabolite_Raw", "Metabolite", "N", "Events", "HR", "CI_Lower", "CI_Upper", "P_Value", "P_Bonferroni", "Output")
  
  # Sort by P-value
  final_df <- final_df[order(final_df$P_Value), ]
  
  # Save
  write.csv(final_df, out_file, row.names = FALSE)
  cat(sprintf("Saved results to %s\n", out_file))
  
  # Print top 5
  print(head(final_df[, c("Metabolite", "Output", "P_Bonferroni")], 5))
}

main <- function() {
  cat("Loading data...\n")
  
  # 1. Load Data
  if (!file.exists("NMR_Preprocessed.csv")) stop("NMR_Preprocessed.csv not found")
  if (!file.exists("Covariates.csv")) stop("Covariates.csv not found")
  if (!file.exists("bladdercancer2.csv")) stop("bladdercancer2.csv not found")
  
  nmr_df <- read.csv("NMR_Preprocessed.csv", check.names = FALSE)
  cov_df <- read.csv("Covariates.csv")
  tgt_df <- read.csv("bladdercancer2.csv")
  
  # 2. Preprocessing
  cat("Preprocessing...\n")
  
  # Identify NMR columns (containing "| Instance 0")
  nmr_cols <- grep("\\| Instance 0", names(nmr_df), value = TRUE)
  if (length(nmr_cols) == 0) {
    stop("No NMR columns found with '| Instance 0' suffix.")
  }
  
  # Subset NMR data
  nmr_df <- nmr_df[, c("eid", nmr_cols)]
  
  # Merge data
  df <- merge(nmr_df, cov_df, by = "eid")
  df <- merge(df, tgt_df, by = "eid")
  
  # --- Fix: Ensure Categorical Variables are Factors ---
  # Assuming 'Smoke' and 'Statin' are categorical covariates
  if ("Smoke" %in% names(df)) df$Smoke <- as.factor(df$Smoke)
  if ("Statin" %in% names(df)) df$Statin <- as.factor(df$Statin)
  if ("Sex" %in% names(df)) df$Sex <- as.factor(df$Sex) # Crucial: Sex is now a covariate
  # --------------------------------------------------
  
  # Define base covariates
  base_covariates <- c("Age", "TDI", "BMI", "Smoke", "Statin", "FastingTime")
  
  # Define full covariates (with Sex)
  full_covariates <- c("Age", "Sex", "TDI", "BMI", "Smoke", "Statin", "FastingTime")
  
  # Check if covariates exist
  missing_covs <- setdiff(full_covariates, names(df))
  if (length(missing_covs) > 0) {
    stop(paste("Missing covariates:", paste(missing_covs, collapse = ", ")))
  }
  
  # 3. Analysis Pipeline Wrapper
  run_subset <- function(sub_df, label, filename) {
      if (nrow(sub_df) == 0) {
          cat(sprintf("No samples for %s. Skipping.\n", label))
          return()
      }
      
      n_sexes <- length(unique(sub_df$Sex))
      if (n_sexes > 1) {
          cat(sprintf("[%s] Both sexes present. Using full model.\n", label))
          covs <- full_covariates
      } else {
          cat(sprintf("[%s] Only one sex present. Using base model.\n", label))
          covs <- base_covariates
      }
      
      run_analysis_pipeline(sub_df, label, filename, nmr_cols, covs)
  }

  # --- Run Analysis for Age >= 65 ---
  if ("Age" %in% names(df)) {
      df_old <- df[df$Age >= 65, ]
      cat(sprintf("Subset Age >= 65: %d samples\n", nrow(df_old)))
      run_subset(df_old, "Age >= 65", "Cox_Analysis_Results_Age65Plus_R.csv")
      
      # --- Run Analysis for Age < 65 ---
      df_young <- df[df$Age < 65, ]
      cat(sprintf("Subset Age < 65: %d samples\n", nrow(df_young)))
      run_subset(df_young, "Age < 65", "Cox_Analysis_Results_AgeUnder65_R.csv")
  } else {
      stop("Age column not found")
  }
}

main()
