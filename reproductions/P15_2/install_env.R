if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager", repos = "https://cloud.r-project.org")
stopifnot(startsWith(as.character(getRversion()), "4.3"))
BiocManager::install(version = "3.18", ask = FALSE, update = FALSE)
BiocManager::install(c("annotatr", "regioneR", "GenomicRanges", "rtracklayer", "AnnotationHub",
                       "TxDb.Hsapiens.UCSC.hg38.knownGene", "org.Hs.eg.db"), ask = FALSE, update = FALSE)
install.packages(c("data.table", "R.utils", "ggplot2"), repos = "https://cloud.r-project.org")
for (p in c("annotatr", "regioneR", "TxDb.Hsapiens.UCSC.hg38.knownGene"))
  cat(sprintf("%-36s %s\n", p, as.character(packageVersion(p))))
