# ~~~ reproductions/P03/scripts/reproduce_results.R ~~~
# - reproduce Chan et al. 2022 figure 5-4c genomic annotation enrichment
# * make sure your CWD is set to .../null-model-sensitivity/reproductions/P03
library(rtracklayer)
library(regioneR)
library(TxDb.Hsapiens.UCSC.hg38.knownGene)

data_dir <- "data"

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
all_genes <- genes(txdb)
canonical_genes <- all_genes[seqnames(all_genes) %in% seqnames(human_canon)]

# define targets
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
for (state_name in chromatin_state_names) {
  state_regions <- subset(chromatin_states, name == state_name)
  
  for (target_name in names(targets)) {
    test <- permTest(
      A = state_regions, B = targets[[target_name]], ntimes = 100,
      genome = human_canon, universe = chromatin_states,
      randomize.function = resampleRegions, evaluate.function = numOverlaps,
      per.chromosome = FALSE, verbose = TRUE
    )
    
    test_summary <- summary(test)
    
    results <- rbind(
      results,
      data.frame(
        chromatin_state = state_name,
        target = target_name,
        p_value = test_summary$pvalue,
        z_score = test_summary$zscore
      )
    )
  }
}
