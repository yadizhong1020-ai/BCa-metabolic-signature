
rm(list=ls()) 
library(data.table) 
library(dplyr)
library(car)

# Source local DEswan function
source("DEswan_source.R")

# ============================================================================== 
# 1. Load Data
# ============================================================================== 

metabolite_data_raw <- as.data.frame(fread("Finalmeta.csv", sep = ",", header = T)) 
covariates_df <- as.data.frame(fread("Covariates2.csv", sep = ",", header = T)) 

metabolite_data <- merge(covariates_df, metabolite_data_raw, by = "eid") 

# ============================================================================== 
# 2. Setup Parameters
# ============================================================================== 

buckets_size <- 5  
qt_variable <- as.numeric(metabolite_data$Age)
min_age <- floor(min(qt_variable, na.rm=T)) + 5
max_age <- ceiling(max(qt_variable, na.rm=T)) - 5
window_center <- seq(min_age, max_age, by = 1) 

covariates_input <- metabolite_data[, c("Sex", "Smoke", "BMI", "target_y")]

metabolite_list <- names(metabolite_data_raw)[-1]

# ============================================================================== 
# 3. Run DEswan in Chunks
# ============================================================================== 

chunk_size <- 50
n_vars <- length(metabolite_list)
n_chunks <- ceiling(n_vars / chunk_size)

all_p <- list()
all_coeff <- list()

print(paste0("Processing buckets size: ", buckets_size)) 
cat("Total variables:", n_vars, "in", n_chunks, "chunks.\n")

for (i in 1:n_chunks) {
    start_idx <- (i - 1) * chunk_size + 1
    end_idx <- min(i * chunk_size, n_vars)
    
    current_vars <- metabolite_list[start_idx:end_idx]
    cat("\n--- Processing Chunk", i, "/", n_chunks, "(Vars", start_idx, "-", end_idx, ") ---\n")
    
    data_input_chunk <- metabolite_data[, current_vars, drop=FALSE]
    
    # Run DEswan for this chunk safely
    res_chunk <- tryCatch({
        DEswan(
          data.df = data_input_chunk, 
          qt = qt_variable, 
          window.center = window_center, 
          buckets.size = buckets_size, 
          covariates = covariates_input
        )
    }, error = function(e) {
        cat("Error in chunk", i, ":", e$message, "\n")
        return(NULL)
    })
    
    if (!is.null(res_chunk)) {
        # Save temp files
        writepath_p_chunk <- paste0("temp_p_chunk_", i, ".txt")
        writepath_coeff_chunk <- paste0("temp_coeff_chunk_", i, ".txt")
        
        fwrite(res_chunk$p, writepath_p_chunk, sep="\t", row.names=F, col.names=T, quote=F)
        fwrite(res_chunk$coeff, writepath_coeff_chunk, sep="\t", row.names=F, col.names=T, quote=F)
        
        all_p[[i]] <- res_chunk$p
        all_coeff[[i]] <- res_chunk$coeff
        cat("Chunk", i, "completed and saved.\n")
    } else {
        cat("Chunk", i, "failed. Skipping.\n")
    }
}

# ============================================================================== 
# 4. Merge and Save Final Results
# ============================================================================== 

cat("\nMerging results...\n")
final_p <- do.call(rbind, all_p)
final_coeff <- do.call(rbind, all_coeff)

if (!is.null(final_p)) {
    writepath_p <- paste0("deswan_buckets_size_", buckets_size, "_p_df.txt") 
    fwrite(final_p, writepath_p, sep="\t", row.names=F, col.names=T, quote=F) 
    cat("- Saved P-values:", writepath_p, "\n")
}

if (!is.null(final_coeff)) {
    writepath_coeff <- paste0("deswan_buckets_size_", buckets_size, "_coeff_df.txt") 
    fwrite(final_coeff, writepath_coeff, sep="\t", row.names=F, col.names=T, quote=F)
    cat("- Saved Coefficients:", writepath_coeff, "\n")
}

# Clean up temp files
temp_files <- list.files(pattern = "temp_.*_chunk_.*\\.txt")
if (length(temp_files) > 0) {
    file.remove(temp_files)
    cat("Cleaned up temp files.\n")
}

cat("DEswan analysis complete.\n")
