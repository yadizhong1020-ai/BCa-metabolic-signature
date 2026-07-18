
rm(list=ls()) 
library(data.table) 
library(dplyr) 

disease_list <- c("BladderCancer") 
bins <- c(-15, -10, -5, -3, 0) 
bins_level <- levels(cut(c(min(bins), max(bins)), breaks = bins, include.lowest = TRUE)) 

COVAR_FILE <- "Covariates2.csv" 
MATCH_FILE <- "2matchresult.csv" 
METABOLITE_FILE <- "Finalmeta.csv" 

covariates_df <- as.data.frame(fread(COVAR_FILE, sep = ",", header = T)) 
metabolite_data_raw <- as.data.frame(fread(METABOLITE_FILE, sep = ",", header = T)) 

# --- CRITICAL FIX: Ensure correct data types ---
# Convert Smoke to factor
covariates_df$Smoke <- as.factor(covariates_df$Smoke)
# Convert numeric covariates explicitly to numeric to avoid them being treated as factors
covariates_df$FastingTime <- as.numeric(covariates_df$FastingTime)
covariates_df$Statin <- as.numeric(covariates_df$Statin)

metabolite_list <- names(metabolite_data_raw)[-1] 

metabolite_data <- merge(covariates_df, metabolite_data_raw, by = "eid") 
 
base_covariates_list <- c("Smoke", "FastingTime", "Statin") 

cat("Starting residual calculation for", length(metabolite_list), "metabolites...\n")

for (metabolite in metabolite_list) { 
  temp_model_vars <- c(metabolite, base_covariates_list)
  rows_to_keep <- complete.cases(metabolite_data[, temp_model_vars])
  
  if (sum(rows_to_keep) < 2) {
      cat("Skipping", metabolite, "- not enough data\n")
      next
  }
  subset_data <- metabolite_data[rows_to_keep, ]
  current_covariates <- c()
  
  for (covar in base_covariates_list) {
      vals <- subset_data[[covar]]
      
      if (is.factor(vals) || is.character(vals)) {
          if (length(unique(vals)) >= 2) {
              current_covariates <- c(current_covariates, covar)
          }
      } else {
          if (var(vals, na.rm=TRUE) > 0) {
              current_covariates <- c(current_covariates, covar)
          }
      }
  }
  if (length(current_covariates) == 0) {
      formula <- as.formula(paste0(metabolite, " ~ 1"))
  } else {
      formula <- as.formula(paste0(metabolite, " ~ ", paste(current_covariates, collapse = " + ")))
  }
  model_resid <- resid(lm(formula, data = metabolite_data, subset = rows_to_keep)) 
  metabolite_data[rows_to_keep, metabolite] <- model_resid 
} 

cat("Residual calculation completed.\n")

cat("Preparing data for Z-score calculation...\n")

metabolite_data <- metabolite_data[, c("eid", metabolite_list)]

result <- data.frame(Metabolite = metabolite_list) 
for (bin in bins_level) { 
  result[[bin]] <- NA 
} 

for (i in 1:length(disease_list)) { 
  disease <- disease_list[i] 
  cat("Processing disease:", disease, "\n")
  
  disease_result <- result 
  
  match_file_path <- paste0("matched_result/", disease, "_matched_result.txt")
  if (!file.exists(match_file_path)) {
      cat("  Error: Matched file not found:", match_file_path, "\n")
      # Try searching in current directory as fallback
      fallback_path <- paste0(disease, "_matched_result.txt")
      if (file.exists(fallback_path)) {
          match_file_path <- fallback_path
          cat("  Found in current directory:", match_file_path, "\n")
      } else {
          next
      }
  }
  
  matched_data <- fread(match_file_path, sep = "\t", header = T) 
  
  if (!"target_y" %in% names(matched_data)) {
      cat("  Error: 'target_y' column missing in matched data.\n")
      next
  }
  
  matched_data <- matched_data %>% 
    group_by(subclass) %>% 
    mutate(BL2Target_yrs_for_1 = ifelse(target_y == 1, BL2Target_yrs, 0)) %>% 
    mutate(BL2Target_yrs = sum(BL2Target_yrs_for_1)) %>% 
    ungroup() %>% 
    select(-BL2Target_yrs_for_1) 
    # bind_rows() removed as it's not needed after ungroup() on a dataframe
  
  matched_data <- as.data.frame(matched_data) 
  
  matched_data <- matched_data[matched_data$BL2Target_yrs <= 15 & matched_data$BL2Target_yrs >= 0, ] 
  
  matched_data$BL2Target_yrs <- -matched_data$BL2Target_yrs 
  matched_data$bin <- cut(matched_data$BL2Target_yrs, bins) 
  
  for (bin in bins_level) { 
    data <- matched_data[matched_data$bin == bin, ] 
    # remove NA bins
    data <- data[!is.na(data$bin), ]
    
    if (nrow(data) == 0) {
        next
    }
    
    case_ids <- data[data$target_y == 1,]$eid 
    control_ids <- data[data$target_y == 0,]$eid 
    
    if (length(case_ids) == 0 || length(control_ids) == 0) {
        next
    }
    
    case_metabolite_data <- metabolite_data[metabolite_data$eid %in% case_ids, -1, drop=FALSE] 
    control_metabolite_data <- metabolite_data[metabolite_data$eid %in% control_ids, -1, drop=FALSE] 
    
    case_metabolite_means <- colMeans(case_metabolite_data, na.rm = TRUE) 
    control_metabolite_means <- colMeans(control_metabolite_data, na.rm = TRUE) 
    control_metabolite_sd <- apply(control_metabolite_data, 2, sd, na.rm = TRUE) 
    
    case_control_meta_zscore <- (case_metabolite_means - control_metabolite_means) / control_metabolite_sd 
    
    disease_result[[bin]] <- case_control_meta_zscore 
  } 
  
  writepath <- paste0("matrix_result/", disease, "_metabolite_matrix_result.txt") 
  fwrite(disease_result, writepath, sep = " ", row.names = F, col.names = T, quote = F, na = "NA") 
  cat("  Result written to:", writepath, "\n")
}
cat("Analysis completed.\n")
