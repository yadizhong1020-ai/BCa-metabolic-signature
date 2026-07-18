
rm(list=ls())
library(data.table)
library(dplyr)
library(pheatmap)
library(grid)

METABOLITE_FILE <- "Finalmeta.csv" 
COVAR_FILE <- "Covariates2.csv"
MATCH_FILE <- "matched_result/BladderCancer_matched_result.txt"
if (!file.exists(MATCH_FILE)) MATCH_FILE <- "BladderCancer_matched_result.txt"
if (!file.exists(MATCH_FILE)) MATCH_FILE <- "matched_result.csv" 

OUTPUT_PLOT <- "Metabolite_Cluster_Heatmap_Fixed.pdf"
time_points <- seq(0, 15, by = 0.5) 

cat("Loading data...\n")
meta_df <- as.data.frame(fread(METABOLITE_FILE))
cov_df <- as.data.frame(fread(COVAR_FILE))

if (file.exists(MATCH_FILE)) {
    match_df <- fread(MATCH_FILE)
} else {
    match_files <- list.files(pattern = "matched_result.*\\.txt|matched_result.*\\.csv")
    if(length(match_files) > 0) {
        match_df <- fread(match_files[1])
        cat("Using match file:", match_files[1], "\n")
    } else {
        stop("No matched result file found.")
    }
}

match_df <- match_df %>%
  group_by(subclass) %>%
  mutate(Time_to_Diagnosis = sum(ifelse(target_y == 1, BL2Target_yrs, 0))) %>%
  ungroup() %>%
  filter(Time_to_Diagnosis >= 0 & Time_to_Diagnosis <= 15) %>%
  mutate(Years_Before_Diagnosis = -Time_to_Diagnosis)

full_data <- merge(cov_df, meta_df, by = "eid")

base_covariates <- c("Smoke", "FastingTime", "Statin") 
full_data$Smoke <- as.factor(full_data$Smoke)
full_data$FastingTime <- as.numeric(full_data$FastingTime)
full_data$Statin <- as.numeric(full_data$Statin)

metabolite_list <- names(meta_df)[-1]
# Debug: metabolite_list <- metabolite_list[1:20]


trajectory_matrix <- matrix(NA, nrow = length(metabolite_list), ncol = length(time_points))
rownames(trajectory_matrix) <- metabolite_list
colnames(trajectory_matrix) <- paste0("Y_", time_points)

cat("Calculating trajectories for", length(metabolite_list), "metabolites...\n")

for (i in seq_along(metabolite_list)) {
    metabolite <- metabolite_list[i]
    
    model_vars <- c(metabolite, base_covariates)
    rows_to_keep <- complete.cases(full_data[, model_vars])
    if (sum(rows_to_keep) < 50) next
    
    subset_data <- full_data[rows_to_keep, ]
    
    current_covs <- c()
    for (cv in base_covariates) {
        if (length(unique(subset_data[[cv]])) >= 2) current_covs <- c(current_covs, cv)
    }
    
    if (length(current_covs) == 0) form <- as.formula(paste(metabolite, "~ 1"))
    else form <- as.formula(paste(metabolite, "~", paste(current_covs, collapse="+")))
    
    fit <- lm(form, data = subset_data)
    subset_data$Resid <- resid(fit)
    
    clean_meta <- subset_data[, c("eid", "Resid")]
    plot_data <- merge(match_df, clean_meta, by = "eid")
    
    if (nrow(plot_data) < 20) next
    
    diff_data <- plot_data %>%
      group_by(subclass, Years_Before_Diagnosis) %>%
      summarise(
        Diff = mean(Resid[target_y == 1], na.rm=T) - mean(Resid[target_y == 0], na.rm=T),
        .groups = "drop"
      ) %>%
      filter(!is.na(Diff))
    
    if (nrow(diff_data) < 10) next

    loess_fit <- try(loess(Diff ~ Years_Before_Diagnosis, data = diff_data, span = 0.75), silent=TRUE)
    
    if (!inherits(loess_fit, "try-error")) {
        pred_x <- -time_points 
        preds <- predict(loess_fit, newdata = data.frame(Years_Before_Diagnosis = pred_x))
        

        if (any(is.na(preds))) {
            valid_idx <- which(!is.na(preds))
            if (length(valid_idx) > 0) {
                if (valid_idx[1] > 1) {
                    preds[1:(valid_idx[1]-1)] <- preds[valid_idx[1]]
                }
                if (valid_idx[length(valid_idx)] < length(preds)) {
                    preds[(valid_idx[length(valid_idx)]+1):length(preds)] <- preds[valid_idx[length(valid_idx)]]
                }
            }
        }
        # -------------------------
        
        trajectory_matrix[i, ] <- preds
    }
    
    if (i %% 20 == 0) cat(".")
}
cat("\n")

valid_rows <- apply(trajectory_matrix, 1, function(x) !all(is.na(x)))
trajectory_matrix <- trajectory_matrix[valid_rows, ]

z_matrix <- t(scale(t(trajectory_matrix)))
z_matrix[is.na(z_matrix)] <- 0

# dist_mat <- dist(z_matrix, method = "euclidean")
# hclust_res <- hclust(dist_mat, method = "ward.D2")

cat("Generating heatmap...\n")

col_anno <- data.frame(Years_Before = time_points)
rownames(col_anno) <- colnames(z_matrix)

pdf(OUTPUT_PLOT, width = 10, height = 15)
pheatmap(z_matrix, 
         cluster_rows = TRUE, 
         cluster_cols = FALSE, 
         show_rownames = TRUE, 
         show_colnames = FALSE,
         clustering_method = "ward.D2",
         clustering_distance_rows = "euclidean",
         fontsize_row = 5,
         main = "Metabolite Trajectory Clustering (Bladder Cancer)",
         annotation_col = col_anno,
         color = colorRampPalette(c("navy", "white", "firebrick3"))(100))
dev.off()

cat("Heatmap saved to:", OUTPUT_PLOT, "\n")
