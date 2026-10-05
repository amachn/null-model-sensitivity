# Recreate Figure S7 of Wenz et al. 2026 (AJHG 113:260) directly from the authors'
# published enrichment table (Data S1 = mmc2.xlsx, Table S3).
# Run from P03/:  Rscript code/figS7_from_table.R
suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
  library(tidyr)
})

xlsx <- "data/mmc2.xlsx"
dir.create("results", showWarnings = FALSE)

# colours used in the published figure (RColorBrewer Dark2 #4 and #5)
set_cols <- c("Random Region Median" = "#E7298A", "ATAC Peak Regions" = "#66A61E")

bar_plot <- function(df, xlab, title = NULL) {
  long <- pivot_longer(df, c("Random Region Median", "ATAC Peak Regions"),
                       names_to = "Set", values_to = "n")
  long$Set <- factor(long$Set, levels = names(set_cols))
  ggplot(long, aes(x = label, y = n, fill = Set)) +
    geom_col(position = "dodge") +
    scale_fill_manual(values = set_cols) +
    labs(x = xlab, y = "Number of Peaks", title = title) +
    theme_bw() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold", size = 11),
          axis.text.y = element_text(face = "bold", size = 11),
          axis.title.y = element_text(face = "bold", size = 14),
          axis.title.x = element_text(size = 14),
          plot.title = element_text(face = "bold.italic", size = 16),
          plot.margin = margin(5.5, 5.5, 5.5, 40),
          legend.position = "bottom")
}

# ---- Figure S7: Table S3 (ENCODE/Roadmap liver histone marks, 1000 shuffles) ----
s3 <- read_excel(xlsx, sheet = "S3 - AllLiverPeaks NIH Roadmap", skip = 2)
s7 <- data.frame(label = s3$EpigenomicMark,
                 `Random Region Median` = s3$Shuffled_Median,
                 `ATAC Peak Regions` = s3$caQTL_Count, check.names = FALSE)
p7 <- bar_plot(s7, "Epigenomic Mark")
ggsave("results/FigS7_from_TableS3.png", p7, width = 7, height = 6.5, dpi = 300)
ggsave("results/FigS7_from_TableS3.pdf", p7, width = 7, height = 6.5)

# consistency check: recompute the median and enrichment from the stored shuffled counts
r3 <- as.matrix(s3[, grep("^Shuffled_.*_Count$", names(s3))])
chk <- data.frame(feature = s3$EpigenomicMark, n_shuffles = ncol(r3),
                  median = apply(r3, 1, median),
                  enrichment = s3$caQTL_Count / apply(r3, 1, median),
                  paper_median = s3$Shuffled_Median, paper_enrichment = s3$Enrichment)
print(chk, digits = 4)
