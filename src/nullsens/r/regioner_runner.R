# ~~~ src/nullsens/r/regioner_runner.R ~~~
# - R backend entry point that runs the full process of:
#   - input -> data loading -> randomization -> permutation -> output
suppressPackageStartupMessages({
  library(jsonlite)
  library(regioneR)
  library(rtracklayer)
  library(Biostrings)
  library(GenomeInfoDb)
  library(GenomicRanges)
  library(IRanges)
  library(S4Vectors)
})

source("io.R")

source("nulls/chromatin.R")
source("nulls/gc.R")
source("nulls/mappability.R")
source("nulls/build_randomizer.R")

# read the JSON request from stdin
input <- paste(readLines(file("stdin"), warn = FALSE), collapse = "\n")
config <- jsonlite::fromJSON(input)

# load region sets
genome_gr <- as_genome_gr(config$genome)
regions_a <- conform_regions(import_regions(config$regions_a), genome_gr)
regions_b <- conform_regions(import_regions(config$regions_b), genome_gr)

# configure reproducibility
if (!is.null(config$seed)) {
  set.seed(config$seed)
}

# timing - start
analysis_start <- proc.time()[["elapsed"]]

# build requested null model
randomize_function <- build_randomizer(
  config       = config$null,
  refs         = config$references,
  regions      = regions_a,
  genome       = config$genome,
  max_attempts = config$max_attempts
)

# run one enrichment test
backend_output <- capture.output({
  test <- regioneR::permTest(
    A                  = regions_a,
    B                  = regions_b,
    ntimes             = config$permutations,
    genome             = genome_gr,
    randomize.function = randomize_function,
    evaluate.function  = regioneR::numOverlaps,
    force.parallel     = config$cores > 1L,
    mc.cores           = config$cores,
    mc.set.seed        = FALSE,
    verbose            = FALSE
  )
})

# timing - end
runtime_seconds <- proc.time()[["elapsed"]] - analysis_start

# write errors to stderr
if (length(backend_output) > 0) { writeLines(backend_output, con = stderr()) }

# write standardized output 
result <- test[[1]]

null_mean <- mean(result$permuted, na.rm = TRUE)

# - this format allows the Python backend to unpack the JSON result with **
output <- list(
  n_peaks            = length(regions_a),
  observed           = as.numeric(result$observed),
  permutation_values = as.numeric(result$permuted),
  null_mean          = null_mean,
  null_sd            = sd(result$permuted, na.rm = TRUE),
  log2_fold_change   = log2(as.numeric(result$observed) / null_mean),
  z_score            = as.numeric(result$zscore),
  p_value            = as.numeric(result$pval),
  runtime_seconds    = as.numeric(runtime_seconds)
)

cat(jsonlite::toJSON(output, auto_unbox = TRUE, digits = NA))