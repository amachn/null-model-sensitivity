# ~~~ reproductions/P03/scripts/reproduce_results.R ~~~
# - reproduce Chan et al. 2022 figure 5-4c genomic annotation enrichment
# * make sure your CWD is set to .../null-model-sensitivity/reproductions/P03
library(rtracklayer)
library(regioneR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

data_dir <- "data"
results_dir <- "results"

# define file paths
chromatin_file <- file.path(data_dir, "hg38_im90_chromatin_states.bed")
nad_file <- file.path(data_dir, "hg38_im90_nads.bed")
lad_file <- file.path(data_dir, "hg38_im90_lads.bed")

# load genomic intervals
chromatin_states <- import.bed(chromatin_file)
nads <- import.bed(nad_file)
lads <- import.bed(lad_file)

# construct the canonical hg38 genome (chr1-22, chrX, chrY)
human_canon <- getGenomeAndMask("hg38", mask = NA)$genome |>
  filterChromosomes(chr.type = "canonical", organism = "hg")

# create hg38 gene annotations and restrict to canonical genes
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
all_genes <- suppressMessages(genes(txdb))
canonical_genes <- all_genes[seqnames(all_genes) %in% seqnames(human_canon)]

# define genomic annotation targets
# - genes are taken from the hg38 UCSC knownGene TxDb and restricted
# - NAD/LAD intervals are from the IMR-90 files found in the original repo
chromatin_state_names <- unique(chromatin_states$name)

targets <- list(
  Genes = canonical_genes,
  NADs  = nads,
  LADs  = lads
)

results <- data.frame(
  chromatin_state = character(),
  target = character(),
  p_value = numeric(),
  z_score = numeric()
)

# permutations
# - for each of the 25 chromatin states, test overlap w/ genes, NADs, and LADs
# - regioneR::permTest parameters match the author's original 
# - no seed is set because the original analysis doesn't use one, 
#   and thus the results can vary quantitatively
total_tests <- length(chromatin_state_names) * length(targets)
test_num <- 0

cat("Running Figure 5-4c permutation tests...\n")
cat("Chromatin states:", length(chromatin_state_names), "\n")
cat("Targets:", paste(names(targets), collapse = ", "), "\n")
cat("Total tests:", total_tests, "\n\n")

for (state in chromatin_state_names) {
  # get all intervals of the current chromatin state
  state_regions <- subset(chromatin_states, name == state)
  
  for (target in names(targets)) {
    test_num <- test_num + 1
    
    cat(
      sprintf(
        "[%d/%d] %s vs %s | state intervals: %d ... ", 
        test_num, total_tests, state, target, length(state_regions)
      )
    )
    
    test <- permTest(
      A = state_regions, B = targets[[target]], ntimes = 100,
      genome = human_canon, universe = chromatin_states,
      randomize.function = resampleRegions, evaluate.function = numOverlaps,
      per.chromosome = FALSE, verbose = FALSE
    )
    
    # hide summary output
    capture.output(
      test_summary <- summary(test)
    )
    
    cat(
      sprintf(
        "p = %.6f, z = %.4f\n", 
        test_summary$pvalue, test_summary$zscore
      )
    )
    
    # append results to the dataframe
    results <- rbind(
      results,
      data.frame(
        chromatin_state = state,
        target = target,
        p_value = test_summary$pvalue,
        z_score = test_summary$zscore
      )
    )
  }
}

cat("\nPermutation testing complete!\n")

# write the raw results to csv
write.csv(
  results, file.path(results_dir, "figure5-4c_results.csv"), row.names = FALSE
)

# convert the results into a z-score matrix usable for the figure
zscore_mat <- matrix(
  results$z_score,
  nrow = length(chromatin_state_names),
  ncol = length(targets),
  byrow = TRUE,
  dimnames = list(
    chromatin_state_names,
    names(targets)
  )
)
write.csv(zscore_mat, file.path(results_dir, "figure5-4c_zscores.csv"))

cat("\nZ-score matrix:\n")
print(zscore_mat)