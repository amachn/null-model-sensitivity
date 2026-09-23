# ~~~ reproductions/P03/scripts/build_figure.R ~~~
# - build a figure similar to Chan et al. figure 5-4c showing reproduced results
# * make sure your CWD is set to .../null-model-sensitivity/reproductions/P03
library(pheatmap)

# order chromatin states to match the original figure
figure_order <- c(
  # misc.
  "21_Het", "24_ReprPC", "20_ZNF/Rpts", "19_DNase", "25_Quies",
  # promoter
  "4_PromD2", "23_PromBiv", "22_PromP", "1_TssA", "2_PromU", "3_PromD1",
  # transcription
  "5_Tx5'", "7_Tx3'", "6_Tx", "8_TxWk",
  # enhancer
  "11_TxEnh3'", "18_EnhAc", "12_TxEnhW", "13_EnhA1", "10_TxEnh5'",
  "9_TxReg", "15_EnhAF", "14_EnhA2", "17_EnhW2", "16_EnhW1"
)

# assign readable labels
figure_labels <- c(
  # misc.
  "heterochromatin", "repressed polycomb", "ZNF repeats", 
  "DNase hypersensitive", "undefined",
  # promoter
  "prom. downstream 2", "prom. bivalent", "prom. poised", "active TSS",
  "prom. upstream TSS", "prom. downstream 1",
  # transcription
  "transcribed - 5'", "transcribed - 3'", 
  "strong transcription", "weak transcription",
  # enhancer
  "3' enhancer", "H3K27ac enhancer", "trans/weak enhancer", "active enhancer 1",
  "5' enhancer", "txn regulartion", "active enhancer flank", "active enhancer 2", 
  "weak enhancer 2", "weak enhancer 1"
)

# read in zscore results and organize by original figure order
results_dir <- "results"
zscore_file <- file.path(results_dir, "figure5-4c_zscores.csv")

zscore_mat <- as.matrix(read.csv(zscore_file, row.names = 1))
plot_mat <- zscore_mat[figure_order, ]

# define heatmap colors and color breaks
heatmap_colors <- colorRampPalette(c("darkgreen", "white", "tomato"))(100)
heatmap_breaks <- seq(-60, 60, length.out = 101)

# plot heatmap
pheatmap(
  plot_mat,
  
  # coloring
  color = heatmap_colors,
  border_color = "grey80",
  breaks = heatmap_breaks,
  
  # labels
  labels_row = figure_labels,
  fontsize_row = 10,
  angle_col = "45",
  
  # disable clustering
  cluster_rows = FALSE,
  cluster_cols = FALSE,
  
  # resize cells and add gaps
  cellwidth = 18, 
  cellheight = 18,
  gaps_row = c(5, 11, 15),
  
  # output file
  filename = file.path(results_dir, "figure5-4c_reproduction.png"),
  width = 5,
  height = 9
)
