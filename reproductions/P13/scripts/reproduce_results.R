# ~~~ reproductions/P13/scripts/reproduce_results.R ~~~
# - reproduce Miller et al. 2022 figure 1c enrichment
# * make sure your CWD is set to .../null-model-sensitivity/reproductions/P13
library(RLSeq)

data_dir <- "data"
results_dir <- "results"

# data downloading
sample_ids <- c("SRX1025894", "SRX1025896", "SRX2683605", "SRX2675009")

#base_url <- "https://rlbase-data.s3.amazonaws.com"
#peak_urls <- paste0(base_url, "/peaks/", sample_ids, "_hg38.broadPeak")
#for (url in peak_urls) {
#  dest <- paste0(data_dir, "/", strsplit(url, split = "/")[[1]][5])
#  download.file(url, dest)
#}

# define results dataframe
results <- data.frame(
  sample_id = character(),
  p_value = numeric(),
  global_z_score = numeric(),
  shift_bp = numeric(),
  local_z_score = numeric()
)

# reproduce local z-score results of all 4 samples
cat("Running Figure 1c permutation tests...\n")
for (i in seq_along(sample_ids)) {
  sample_id <- sample_ids[i]
  cat(sprintf("[%d/%d] %s ... ", i, length(sample_ids), sample_id))
  
  # construct the RLRanges object direct from the broadPeak file
  peak_file <- file.path(data_dir, paste0(sample_id, "_hg38.broadPeak"))
  rlr <- RLRanges(peaks = peak_file, genome = "hg38", sampleName = sample_id)
  
  # run the RLFS enrichment
  rlr <- suppressWarnings(analyzeRLFS(rlr, quiet = TRUE))
  rlfs_result <- rlresult(rlr, "rlfsRes")
  
  # extract results
  perm_result <- rlfs_result$perTestResults[[1]]
  local_result <- rlfs_result$`Z-scores`[[1]]
  
  # store global and local results together in one data frame
  results <- rbind(
    results,
    data.frame(
      sample_id = sample_id,
      p_value = perm_result$pval,
      global_z_score = perm_result$zscore,
      shift_bp = local_result$shifts,
      local_z_score = local_result$shifted.z.scores
    )
  )
  
  cat(
    sprintf(
      "p = %.6f, z = %.4f\n",
      perm_result$pval, perm_result$zscore
    )
  )
}

cat("\nPermutation testing complete!\n")

# write the results to csv
write.csv(
  results, file.path(results_dir, "fig1c_results.csv"), row.names = FALSE
)