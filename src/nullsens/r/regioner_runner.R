# ~~~ src/nullsens/r/regioner_runner.R ~~~
# - R backend entry point that runs the full process of:
#   - input -> data loading -> randomization -> permutation -> output
suppressPackageStartupMessages({
  library(jsonlite)
  library(regioneR)
  library(rtracklayer)
})

source("nulls/build_randomizer.R")
source("nulls/length.R")

# read the JSON request from stdin
input <- paste(readLines(file("stdin"), warn = FALSE), collapse = "\n")
config <- jsonlite::fromJSON(input)

# load region sets
import_regions <- function(path) {
  tryCatch(
    rtracklayer::import(path),
    error = function(err) {
      stop("Unable to import region file: ", path, "\n", conditionMessage(err))
    }
  )
}

regions_a <- import_regions(config$regions_a)
regions_b <- import_regions(config$regions_b)

# configure reproducibility
if (!is.null(config$seed)) {
  set.seed(config$seed)
}

# build requested null model
randomize_function <- build_randomizer(
  config = config$null,
  refs = config$references,
  regions = regions_a,
  genome = config$genome
)

# run one enrichment test
backend_output <- capture.output({
  test <- regioneR::permTest(
    A = regions_a,
    B = regions_b,
    ntimes = config$permutations,
    genome = config$genome,
    randomize.function = randomize_function,
    evaluate.function = regioneR::numOverlaps,
    verbose = FALSE
  )
})

if (length(backend_output) > 0) { writeLines(backend_output, con = stderr()) }

result <- test[[1]]

# write standardized output 
# - this format allows the Python backend to unpack the JSON result with **
output <- list(
  observed = as.numeric(result$observed),
  permutation_values = as.numeric(result$permuted),
  null_mean = mean(result$permuted, na.rm = TRUE),
  null_sd = sd(result$permuted, na.rm = TRUE),
  z_score = as.numeric(result$zscore),
  p_value = as.numeric(result$pval)
)

cat(jsonlite::toJSON(output, auto_unbox = TRUE, digits = NA))