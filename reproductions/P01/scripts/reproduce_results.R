# ~~~ reproductions/P01/scripts/reproduce_results.R ~~~
# - reproduce Vaz et al. 2025 figure 3e genomic annotation enrichment
# * make sure your CWD is set to .../null-model-sensitivity/reproductions/P01
library(regioneR)

data_dir <- "data"
results_dir <- "results"

# baseline sncRNA genomic coords
sncRNAs <- c(
  miRNA = "BL-mirs-all-matches-hg38-coord-6c-2608.bed",
  piRNA = "Bl-pirs-all-matches-hg18-6c-lifted-orig-hg38-2408-rem-random.bed",
  tRF   = "Bl-tRFs-all-matches-hg19-6c-lifted-orig-hg38.bed"
)

# genomic annotation coords
annotations <- c(
  Gene       = "chr1-22.gene.bed",
  Intron     = "chr1-22.intron.bed",
  Exon       = "chr1-22.exonplus.bed", # extra annotations, not used in fig3e
  CodingExon = "chr1-22.codingexon.bed",
  CpG        = "chr1-22.CpG.bed",
  UTR5       = "chr1-22.5utr.bed",
  UTR3       = "chr1-22.3utr.bed",
  Downstream = "chr1-22.Gene-2000Down.bed",
  Upstream   = "chr1-22.Gene-2000UP.bed"
)

# reproduction results dataframe
results <- data.frame(
  sncRNA = character(),
  annotation = character(),
  p_value = numeric(),
  z_score = numeric()
)

cat("Starting reproduction...\n")
for (sncRNA_name in names(sncRNAs)) {
  # create a container for the .bed file (chromosome, start, end, name, strand)
  A <- toGRanges(file.path(data_dir, sncRNAs[[sncRNA_name]]))
  
  cat("\nLoaded", sncRNA_name, "container.\n")
  for (annotation_name in names(annotations)) {
    B <- toGRanges(file.path(data_dir, annotations[[annotation_name]]))
    
    cat("\nRunning:", sncRNA_name, "vs", annotation_name, "\n")
    
    # done to mirror the author's method of scripting, wherein this is run
    # each time a comparison was done
    set.seed(12345)
    
    # NOTE: mask has been swapped out to produce two different sets of results
    # - after further analysis, unmasked appears to match the original figure
    start <- Sys.time()
    pt <- overlapPermTest(
      A, B, ntimes = 100, genome = "hg38", mask = NA, # NA or getMask("hg38")
      count.once = TRUE, alternative = "auto"
    )
    end <- Sys.time()
    
    cat(sncRNA_name, "vs", annotation_name, "completed in", end - start, "s\n")
    
    # append results to the dataframe
    results <- rbind(
      results,
      data.frame(
        sncRNA = sncRNA_name,
        annotation = annotation_name,
        p_value = pt$numOverlaps$pval,
        z_score = pt$numOverlaps$zscore
      )
    )
  }
  cat("\nCompleted", sncRNA_name, "container.\n")
}

# write the results to our csv
write.csv( # change file name based on masked/unmasked
  results, file.path(results_dir, "figure3e_unmasked.csv"), row.names = FALSE
)

cat("\nCompleted reproduction!\nResults:\n")
print(results)