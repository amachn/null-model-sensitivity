# Plot Figure S7 from the re-run results (figS7_roadmap_rerun.sh) in the style of the
# published figure, and compare with Table S3.
# Run from P03/:  Rscript code/figS7_plot_rerun.R
suppressPackageStartupMessages({
  library(ggplot2)
  library(readxl)
  library(tidyr)
  library(data.table)
})

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

obs <- fread("results/S7_rerun/observed.csv")
shf <- fread("results/S7_rerun/shuffled_counts.csv")
marks <- setdiff(names(obs), "mark")
o <- unlist(obs[1, ..marks]); r <- as.matrix(shf[, ..marks])
s7 <- data.frame(EpigenomicMark = marks, ATAC_Count = o, Shuffled_Median = apply(r, 2, median),
                 pval_enriched = pmax(colSums(sweep(r, 2, o, ">=")), 1) / nrow(r),
                 pval_depleted = pmax(colSums(sweep(r, 2, o, "<=")), 1) / nrow(r),
                 n_shuffles = nrow(r))
s7$Enrichment <- s7$ATAC_Count / s7$Shuffled_Median
fwrite(s7, "results/TableS3_rerun.csv")
p7 <- bar_plot(data.frame(label = s7$EpigenomicMark, `Random Region Median` = s7$Shuffled_Median,
                          `ATAC Peak Regions` = s7$ATAC_Count, check.names = FALSE),
               "Epigenomic Mark")
ggsave("results/FigS7_rerun.png", p7, width = 7, height = 6.5, dpi = 300)
ggsave("results/FigS7_rerun.pdf", p7, width = 7, height = 6.5)

# ---- comparison with the paper's Table S3 ----
t3 <- read_excel("data/mmc2.xlsx", sheet = "S3 - AllLiverPeaks NIH Roadmap", skip = 2)
cmp <- merge(data.frame(feature = t3$EpigenomicMark, paper_obs = t3$caQTL_Count,
                        paper_random_median = t3$Shuffled_Median, paper_enrichment = t3$Enrichment),
             data.frame(feature = s7$EpigenomicMark, rerun_obs = s7$ATAC_Count,
                        rerun_random_median = s7$Shuffled_Median, rerun_enrichment = s7$Enrichment))
cmp$same_direction <- sign(cmp$paper_enrichment - 1) == sign(cmp$rerun_enrichment - 1)
fwrite(cmp, "results/S7_rerun_vs_paper.csv")
print(cmp, digits = 4)
