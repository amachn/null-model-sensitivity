# ~~~ reproductions/P01/scripts/build_figure.R ~~~
# - build a figure similar to Vaz et al. figure 3e showing reproduced results
# * make sure your CWD is set to .../null-model-sensitivity/reproductions/P01
library(dplyr)
library(ggplot2)

results_dir <- "results"
results <- read.csv(file.path(results_dir, "figure3e_unmasked.csv"))

plot_df <- results |>
  filter(annotation != "Exon") |> # exclude ExonPlus
  mutate(
    annotation = recode( # rename to follow the figure
      annotation,
      CodingExon = "Exon",
      UTR5 = "5'UTR",
      UTR3 = "3'UTR",
      Downstream = "Down",
      Upstream = "Up"
    ),
    annotation = factor( # keep order consistent with the figure
      annotation,
      levels = c(
        "Gene", "Intron", "Exon", "CpG",
        "5'UTR", "3'UTR", "Down", "Up"
      )
    ),
    sig = p_value < 0.01, # determine if the result is significant
    star_y = ifelse(z_score >= 0, z_score + 1, z_score - 3),
    star_label = ifelse(sig, "*", "") # do we add a "*" to this results bar?
  )
  
ggplot(plot_df, aes(x = annotation, y = z_score, fill = sncRNA)) +
  geom_col(position = position_dodge(width = 0.9), width = 0.8) +
  geom_text(
    aes(label = star_label, y = star_y),
    position = position_dodge(width = 0.9),
    size = 6
  ) +
  geom_hline(yintercept = 0, linewidth = 0.4) +
  geom_hline(yintercept = c(-2, 2), linetype = "dashed", linewidth = 0.4) +
  labs(
    title = "Reproduced Figure 3e (Vaz et al. 2025)",
    x = "Genomic Regions",
    y = "Z-score"
  ) +
  theme_classic(base_size = 14) +
  theme(
    legend.position = "inside",
    legend.position.inside = c(0.1, 0.9),
    plot.title = element_text(hjust = 0.5)
  )

ggsave(
  file.path(results_dir, "figure3e_reproduction.png"), 
  width = 10, height = 6, dpi = 300
)
