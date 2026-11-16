# run_nulls.R
source("pre-processing.R")     # as_genome_gr, conform, read_any
source("mappability.R")        # prepare_mappability, validate_mappability
source("gc.R")                 # prepare_gc, validate_gc
source("chromatin.R")          # prepare_chromatin, state_of, propose_chromatin, validate_chromatin
source("Anna_and_sheridan.R")  # validate_candidates, build_randomizer

library(regioneR)

# ---- settings (defined once, used everywhere)
gc_tol     <- 0.10
gc_maxn    <- 0.1
min_accept <- 0.02    # drop peaks accepted in < 2% of single-pass proposals
n_draws    <- 100     # draws for the acceptance estimate
ntimes     <- 100      # permutations per row (GC rows are slow: ~45 s each)

# ---- inputs
genome    <- "mm10"                # name: needed by prepare_gc
genome_gr <- as_genome_gr(genome)  # GRanges: for permTest and conforming

A <- reduce(conform(read_any("Mouse-Genome-Data/H3K4me3_Striata_intervals.merged.any_3_samples.final.annotated_with_description.txt"), genome_gr))
B <- reduce(conform(read_any("Mouse-Genome-Data/Htt_ChIPSeq.downsample_30000000.filtered_p1e3.merged_peaks.any_2_samples.final.conditionally_reproducible.txt", header = FALSE), genome_gr))

refs <- list(
  mappability = "Mouse-Genome-Data/mm10-blacklist.v2.bed",
  chromatin   = "C:/Users/sheri/Desktop/Capstone/all_combined/Mouse-Genome-Data/ENCFF252YVM.bed"
)
stopifnot(file.exists(unlist(refs)))

# ---- configs
off <- list(
  mappability_masked = FALSE, mappability = list(threshold = NULL, invert = FALSE),
  gc_matched = FALSE,         gc = list(tolerance = gc_tol, flank = 0, max_n = gc_maxn),
  chromatin_matched = FALSE,  chromatin = list(pad = 0, label_col = NULL,
                                               priority = NULL, strict = FALSE)
)

chrom_cfg <- list(pad = 0, label_col = "class",
                  priority = c("PLS", "pELS", "dELS", "CTCF-only",
                               "DNase-H3K4me3", "DNase-only", "Low-DNase"),
                  strict = FALSE)

configs <- list(
  baseline      = off,
  mappability   = modifyList(off, list(mappability_masked = TRUE)),
  gc            = modifyList(off, list(gc_matched = TRUE)),
  chromatin     = modifyList(off, list(chromatin_matched = TRUE, chromatin = chrom_cfg)),
  map_gc        = modifyList(off, list(mappability_masked = TRUE, gc_matched = TRUE)),
  map_chromatin = modifyList(off, list(mappability_masked = TRUE,
                                       chromatin_matched = TRUE, chromatin = chrom_cfg)),
  gc_chromatin  = modifyList(off, list(gc_matched = TRUE,
                                       chromatin_matched = TRUE, chromatin = chrom_cfg)),
  all_three     = modifyList(off, list(mappability_masked = TRUE, gc_matched = TRUE,
                                       chromatin_matched = TRUE, chromatin = chrom_cfg))
)

# ---- step 1: drop peaks that fail their own GC check (e.g. too many Ns)
gc_ctx <- prepare_gc(genome, tolerance = gc_tol, flank = 0, max_n = gc_maxn)
A_gc   <- A[validate_gc(gc_ctx, A, A)]
message("Dropped ", length(A) - length(A_gc), " peaks failing their own GC check")

# ---- step 2: per-region acceptance under each GC-containing config
# max_attempts = 1 -> one propose+validate pass; diag$attempts > 1 marks failures
region_acceptance <- function(config, A_in, draws = n_draws) {
  rnd <- build_randomizer(config, refs, regions = A_in, genome = genome,
                          max_attempts = 1L)
  fails <- integer(length(A_in))
  for (i in seq_len(draws)) {
    tryCatch(rnd(A_in, genome = genome_gr), error = function(e) NULL)
    fails <- fails + (environment(rnd)$diag$attempts > 1L)
  }
  1 - fails / draws
}

gc_cfgs  <- configs[sapply(configs, function(c) isTRUE(c$gc_matched))]
set.seed(1)
acc_list <- lapply(gc_cfgs, region_acceptance, A_in = A_gc)

print(sapply(acc_list, function(a) sum(a < min_accept)))   # peaks each config would drop
keep_idx <- Reduce(`&`, lapply(acc_list, function(a) a >= min_accept))
A_common <- A_gc[keep_idx]
print(c(full = length(A), after_gc_check = length(A_gc), common = length(A_common)))

# ---- step 3: run every row on the same peak set
run_null <- function(label, config, A_use = A_common, ntimes = 100,
                     max_attempts = 1000L, stat = numOverlaps) {
  tryCatch({
    rnd <- build_randomizer(config, refs, regions = A_use, genome = genome,
                            max_attempts = max_attempts)
    set.seed(1)
    secs <- system.time(
      pt <- permTest(A = A_use, B = B, ntimes = ntimes, genome = genome_gr,
                     randomize.function = rnd, evaluate.function = stat,
                     verbose = FALSE))[["elapsed"]]
    r <- pt[[1]]; e <- mean(r$permuted)
    data.frame(null = label, n_peaks = length(A_use), observed = r$observed,
               expected = e, log2fold = log2(r$observed / e),
               z = r$zscore, p = r$pval, secs = secs, note = "")
  }, error = function(err)
    data.frame(null = label, n_peaks = length(A_use), observed = NA, expected = NA,
               log2fold = NA, z = NA, p = NA, secs = NA, note = conditionMessage(err)))
}

res <- do.call(rbind, Map(run_null, names(configs), configs, ntimes = ntimes))
res
