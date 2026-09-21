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

# hg38 gene annotations
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
genes_hg38 <- genes(txdb)
