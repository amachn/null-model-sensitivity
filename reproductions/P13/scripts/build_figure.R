# ~~~ reproductions/P13/scripts/build_figure.R ~~~
# - build a figure similar to Miller et al. figure 1c showing reproduced results
# * make sure your CWD is set to .../null-model-sensitivity/reproductions/P13
library(ggplot2)
library(ggprism)

# load data
results_dir <- "results"
results_file <- file.path(results_dir, "fig1c_results.csv")
results <- read.csv(results_file)

# define samples for plotting
sample_colors <- c(
  "SRX1025894" = "#43abe1",
  "SRX1025896" = "#79919b",
  "SRX2683605" = "#e662a4",
  "SRX2675009" = "#a08a96"
)

sample_labels <- c(
  "SRX1025894" = "NT2 [DRIP]\nSRX1025894",
  "SRX1025896" = "NT2 (RNH1) [DRIP]\nSRX1025896",
  "SRX2683605" = "HEK293 [R-ChIP]\nSRX2683605",
  "SRX2675009" = "HEK293 (WKKD) [R-ChIP]\nSRX2675009"
)

# factor sample id and scale z-scores to 0-1
results$sample_id <- factor(results$sample_id, levels = names(sample_colors))
results$scaled_z_score <- ave(
  results$local_z_score,
  results$sample_id,
  FUN = function(zscore) {
    (zscore - min(zscore)) / (max(zscore) - min(zscore))
  }
)

# plot results to match Figure 1c
ggplot(results, aes(x = shift_bp, y = scaled_z_score, color = sample_id)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey", alpha = 0.5) +
  geom_line(linewidth = 0.6) +
  facet_wrap(~ sample_id, nrow = 1, labeller = as_labeller(sample_labels)) +
  scale_color_manual(values = sample_colors) +
  scale_y_continuous(limits = c(0, 1.1), expand = c(0, 0)) +
  labs(
    title = "R-loop peak enrichment around RLFS",
    x = "Distance from RLFS (bp)",
    y = "Z-score (scaled)"
  ) +
  theme_prism(border = TRUE, base_size = 10) +
  theme(
    legend.position = "none", 
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  file.path(results_dir, "fig1c_reproduction.png"),
  width = 13, height = 4, dpi = 300
)
