suppressPackageStartupMessages({ library(data.table); library(ggplot2) })

summarise_null <- function(obs, null) {
  null <- null[annotation %in% names(obs)]
  s <- null[, .(null_median = median(n), null_mean = mean(n), null_sd = sd(n),
                null_lo = quantile(n, 0.025), null_hi = quantile(n, 0.975), n_random = .N), by = annotation]
  s[, observed := as.numeric(obs[annotation])]
  p <- null[, .(p_enrich = (1 + sum(n >= obs[annotation[1]])) / (.N + 1),
                p_deplete = (1 + sum(n <= obs[annotation[1]])) / (.N + 1)), by = annotation]
  s <- merge(s, p, by = "annotation")
  s[, `:=`(fold = observed / null_median, z = (observed - null_mean) / null_sd,
           p_two_sided = pmin(1, 2 * pmin(p_enrich, p_deplete)))]
  s[, direction := fifelse(observed > null_median, "enriched", fifelse(observed < null_median, "depleted", "none"))]
  s[, significant := p_two_sided < 0.05]
  s[]
}

compare_to_published <- function(res, fig, figure_id) {
  f <- fig[figure == figure_id]
  m <- merge(res, f[, .(annotation, obs_approx, null_median_approx, direction_in_figure, caption_claim)],
             by = "annotation", all = TRUE)
  m[, `:=`(obs_pct_diff = round(100 * (observed - obs_approx) / obs_approx, 1),
           null_pct_diff = round(100 * (null_median - null_median_approx) / null_median_approx, 1),
           direction_matches_figure = direction == direction_in_figure)]
  m[, call := fifelse(!significant, "not_significant", direction)]
  m[, call_matches_caption := fifelse(is.na(caption_claim), NA, call == caption_claim)]
  m[]
}

plot_paper_style <- function(res, title, file_stem, query_label = "caQTL peaks") {
  d <- rbind(res[, .(annotation, set = "Random Region Median", n = null_median, lo = null_lo, hi = null_hi)],
             res[, .(annotation, set = query_label, n = observed, lo = NA_real_, hi = NA_real_)])
  d[, set := factor(set, levels = c("Random Region Median", query_label))]
  star <- res[significant == TRUE, .(annotation, y = pmax(observed, null_hi) * 1.04)]
  g <- ggplot(d, aes(annotation, n, fill = set)) +
    geom_col(position = position_dodge(0.9), width = 0.9) +
    geom_errorbar(aes(ymin = lo, ymax = hi), position = position_dodge(0.9), width = 0.3, na.rm = TRUE) +
    geom_text(data = star, aes(annotation, y, label = "*"), inherit.aes = FALSE, size = 6) +
    scale_fill_manual(values = c("#E7298A", "#66A61E"), name = "Set") +
    labs(x = "Annotation", y = "Number of Peaks", title = title,
         subtitle = sprintf("error bars: 95%% range of %d random sets; * empirical two-sided p < 0.05", res$n_random[1])) +
    theme_bw(base_size = 11) +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "bold"), legend.position = "bottom",
          plot.title = element_text(face = "bold.italic"), plot.margin = margin(5, 5, 5, 40))
  ggsave(paste0(file_stem, ".pdf"), g, width = 7, height = 6)
  ggsave(paste0(file_stem, ".png"), g, width = 7, height = 6, dpi = 150)
  invisible(g)
}

plot_vs_published <- function(cmp, title, file_stem) {
  d <- rbind(cmp[, .(annotation, what = "observed", reproduced = observed, published = obs_approx)],
             cmp[, .(annotation, what = "random median", reproduced = null_median, published = null_median_approx)])
  g <- ggplot(d, aes(published, reproduced, colour = what)) +
    geom_abline(linetype = 2, colour = "grey50") + geom_point(size = 2.5) +
    geom_text(aes(label = sub("^hg38_(genes_|enhancers_)?", "", annotation)), size = 2.6, vjust = -0.8, show.legend = FALSE) +
    scale_x_log10() + scale_y_log10() +
    labs(x = "Published (read from figure, approx.)", y = "Reproduced", colour = NULL, title = title) +
    theme_bw() + theme(legend.position = "bottom")
  ggsave(paste0(file_stem, ".pdf"), g, width = 6.5, height = 6)
  ggsave(paste0(file_stem, ".png"), g, width = 6.5, height = 6, dpi = 150)
  invisible(g)
}
