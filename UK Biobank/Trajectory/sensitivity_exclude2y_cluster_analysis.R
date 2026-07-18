rm(list = ls())

library(data.table)
library(dplyr)
library(ggplot2)

METABOLITE_FILE <- "Finalmeta.csv"
COVAR_FILE <- "Covariates2.csv"
MATCH_FILE <- "matched_result/BladderCancer_matched_result.txt"
if (!file.exists(MATCH_FILE)) MATCH_FILE <- "BladderCancer_matched_result.txt"

OUTPUT_PLOT <- "Cluster_Trajectory_Hierarchical_exclude2y.pdf"
OUTPUT_CLUSTER <- "Metabolite_Cluster_Assignments_exclude2y.csv"
OUTPUT_TRAJ <- "Metabolite_Cluster_Trajectory_Values_exclude2y.csv"
OUTPUT_REPORT <- "cluster_trends_report_exclude2y.txt"
OUTPUT_CASE_SUMMARY <- "Sensitivity_Exclude2Y_Sample_Summary.csv"

pred_times <- seq(-15, 0, by = 1)
base_covariates <- c("Smoke", "FastingTime", "Statin")

cat("Loading data...\n")
meta_df <- as.data.frame(fread(METABOLITE_FILE))
cov_df <- as.data.frame(fread(COVAR_FILE))
match_df <- fread(MATCH_FILE)

original_case_n <- match_df %>%
  filter(target_y == 1, BL2Target_yrs >= 0, BL2Target_yrs <= 15) %>%
  nrow()


match_df <- match_df %>%
  group_by(subclass) %>%
  mutate(Time_to_Diagnosis = sum(ifelse(target_y == 1, BL2Target_yrs, 0), na.rm = TRUE)) %>%
  ungroup() %>%
  filter(Time_to_Diagnosis > 2, Time_to_Diagnosis <= 15) %>%
  mutate(Years_Before_Diagnosis = -Time_to_Diagnosis)

filtered_case_n <- match_df %>%
  filter(target_y == 1) %>%
  nrow()

sample_summary <- data.frame(
  Metric = c("Original_incident_cases_0_to_15y", "After_excluding_cases_within_2y"),
  Count = c(original_case_n, filtered_case_n)
)
write.csv(sample_summary, OUTPUT_CASE_SUMMARY, row.names = FALSE)
cat("Sample summary saved to:", OUTPUT_CASE_SUMMARY, "\n")

full_data <- merge(cov_df, meta_df, by = "eid")
full_data$Smoke <- as.factor(full_data$Smoke)
full_data$FastingTime <- as.numeric(full_data$FastingTime)
full_data$Statin <- as.numeric(full_data$Statin)

metabolite_list <- names(meta_df)[-1]

plot_data_list <- list()

cat("Calculating trajectories after excluding first 2 years...\n")

for (metabolite in metabolite_list) {
  model_vars <- c(metabolite, base_covariates)
  subset_data <- full_data[complete.cases(full_data[, model_vars]), ]
  if (nrow(subset_data) < 50) next

  current_covs <- c()
  for (cv in base_covariates) {
    if (length(unique(subset_data[[cv]])) >= 2) current_covs <- c(current_covs, cv)
  }

  if (length(current_covs) == 0) {
    fit <- lm(as.formula(paste(metabolite, "~ 1")), data = subset_data)
  } else {
    fit <- lm(as.formula(paste(metabolite, "~", paste(current_covs, collapse = "+"))), data = subset_data)
  }

  subset_data$Resid <- resid(fit)
  clean_meta <- subset_data[, c("eid", "Resid")]
  merged <- merge(match_df, clean_meta, by = "eid")
  if (nrow(merged) < 20) next

  diff_df <- merged %>%
    group_by(subclass, Years_Before_Diagnosis) %>%
    summarise(
      Diff = mean(Resid[target_y == 1], na.rm = TRUE) -
        mean(Resid[target_y == 0], na.rm = TRUE),
      .groups = "drop"
    ) %>%
    filter(!is.na(Diff))

  if (nrow(diff_df) < 10) next

  loess_fit <- try(loess(Diff ~ Years_Before_Diagnosis, data = diff_df, span = 0.8), silent = TRUE)
  if (inherits(loess_fit, "try-error")) next

  preds <- predict(loess_fit, newdata = data.frame(Years_Before_Diagnosis = pred_times))

  if (any(is.na(preds))) {
    valid <- which(!is.na(preds))
    if (length(valid) == 0) next
    if (valid[1] > 1) preds[1:(valid[1] - 1)] <- preds[valid[1]]
    if (valid[length(valid)] < length(preds)) {
      preds[(valid[length(valid)] + 1):length(preds)] <- preds[valid[length(valid)]]
    }
  }

  if (any(is.na(preds))) next

  plot_data_list[[length(plot_data_list) + 1]] <- data.frame(
    Metabolite = metabolite,
    Time = pred_times,
    Value = preds
  )
}

all_traj_df <- do.call(rbind, plot_data_list)

if (is.null(all_traj_df) || nrow(all_traj_df) == 0) {
  stop("No trajectories were generated after excluding the first 2 years.")
}

cat("Performing hierarchical clustering...\n")

setDT(all_traj_df)
wide_mat <- dcast(all_traj_df, Metabolite ~ Time, value.var = "Value")
mat_vals <- as.matrix(wide_mat[, -1])
rownames(mat_vals) <- wide_mat$Metabolite

z_mat <- t(scale(t(mat_vals)))
z_mat[is.na(z_mat)] <- 0

dist_mat <- dist(z_mat, method = "euclidean")
hclust_res <- hclust(dist_mat, method = "ward.D2")
cluster_ids <- cutree(hclust_res, k = 4)

clusters <- data.frame(
  Metabolite = names(cluster_ids),
  Cluster = cluster_ids
)
write.csv(clusters, OUTPUT_CLUSTER, row.names = FALSE)
cat("Cluster assignments saved to:", OUTPUT_CLUSTER, "\n")

time_cols <- paste0("Year_", abs(pred_times), "_Before")
traj_list <- list()

for (metabolite in rownames(z_mat)) {
  vals <- as.numeric(z_mat[metabolite, ])
  row_df <- data.frame(t(vals))
  names(row_df) <- time_cols
  row_df$Metabolite <- metabolite
  traj_list[[metabolite]] <- row_df
}

traj_df <- do.call(rbind, traj_list)
traj_df <- traj_df[, c("Metabolite", time_cols)]
traj_df <- merge(clusters, traj_df, by = "Metabolite")
traj_df <- traj_df %>% arrange(Cluster, Metabolite)
write.csv(traj_df, OUTPUT_TRAJ, row.names = FALSE)
cat("Trajectory values saved to:", OUTPUT_TRAJ, "\n")

cluster_means <- traj_df %>%
  group_by(Cluster) %>%
  summarise(across(all_of(time_cols), mean, na.rm = TRUE), .groups = "drop")

report_lines <- c(
  "=== Sensitivity Analysis: Excluding Cases Diagnosed Within 2 Years ===",
  paste0("Original incident cases (0-15y): ", original_case_n),
  paste0("After exclusion (>2-15y): ", filtered_case_n),
  " "
)

for (i in seq_len(nrow(cluster_means))) {
  cid <- cluster_means$Cluster[i]
  vals <- as.numeric(cluster_means[i, time_cols])
  start_val <- vals[1]
  mid_val <- vals[8]
  end_val <- vals[length(vals)]

  change_total <- end_val - start_val
  change_early <- mid_val - start_val
  change_late <- end_val - mid_val

  trend <- if (change_total > 0.5) {
    "Strongly Rising"
  } else if (change_total < -0.5) {
    "Strongly Falling"
  } else if (abs(change_total) < 0.2 && sd(vals) < 0.2) {
    "Stable"
  } else if (change_early > 0.2 && change_late < -0.2) {
    "Rise then Fall (Inverted U)"
  } else if (change_early < -0.2 && change_late > 0.2) {
    "Fall then Rise (U-shape)"
  } else if (change_total > 0) {
    "Slowly Rising"
  } else {
    "Slowly Falling"
  }

  report_lines <- c(
    report_lines,
    paste0("Cluster ", cid, ":"),
    sprintf("  Start (Year -15): %.3f", start_val),
    sprintf("  Mid   (Year -8):  %.3f", mid_val),
    sprintf("  End   (Year 0):   %.3f", end_val),
    paste0("  Overall Trend: ", trend),
    "------------------------------------------------"
  )
}

writeLines(report_lines, OUTPUT_REPORT)
cat("Trend report saved to:", OUTPUT_REPORT, "\n")

final_plot_df <- merge(all_traj_df, clusters, by = "Metabolite")
final_plot_df <- final_plot_df %>%
  group_by(Metabolite) %>%
  mutate(Z_Value = (Value - mean(Value, na.rm = TRUE)) / sd(Value, na.rm = TRUE)) %>%
  ungroup()

final_plot_df$Z_Value[is.na(final_plot_df$Z_Value)] <- 0
final_plot_df$Cluster_Label <- paste("Cluster", final_plot_df$Cluster)

mean_lines <- final_plot_df %>%
  group_by(Cluster_Label, Time) %>%
  summarise(Z_Value = mean(Z_Value, na.rm = TRUE), .groups = "drop")

p <- ggplot() +
  geom_point(
    data = final_plot_df,
    aes(x = Time, y = Z_Value, color = Cluster_Label),
    size = 1.2,
    alpha = 0.4
  ) +
  geom_line(
    data = final_plot_df,
    aes(x = Time, y = Z_Value, group = Metabolite, color = Cluster_Label),
    alpha = 0.3,
    linewidth = 0.4
  ) +
  geom_line(
    data = mean_lines,
    aes(x = Time, y = Z_Value, color = Cluster_Label),
    linewidth = 1.5
  ) +
  facet_wrap(~Cluster_Label, scales = "free") +
  scale_color_manual(values = c("#D53E4F", "#3288BD", "#66C2A5", "#5E4FA2")) +
  theme_classic() +
  labs(
    x = "Years before diagnosis",
    y = "Relative level (Z-score)",
    title = "Metabolite Trajectory Clusters After Excluding First 2 Years",
    subtitle = "Sensitivity analysis for reverse causation"
  ) +
  theme(
    legend.position = "none",
    strip.background = element_blank(),
    strip.text = element_text(size = 14, face = "bold"),
    axis.title = element_text(size = 12),
    axis.text = element_text(size = 10)
  )

ggsave(OUTPUT_PLOT, p, width = 10, height = 8)
cat("Plot saved to:", OUTPUT_PLOT, "\n")
cat("Sensitivity analysis completed successfully.\n")
