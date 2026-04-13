#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

option_list <- list(
  make_option("--input-root", type = "character", dest = "input_root",
              help = "Root directory containing phenotype-specific LOO outputs"),
  make_option("--outdir", type = "character", dest = "outdir",
              help = "Output directory"),
  make_option("--pattern", type = "character", dest = "pattern",
              default = "^loo_summary_.*\\.tsv$",
              help = "Regex pattern for LOO summary files [default: %default]"),
  make_option("--width", type = "double", dest = "width",
              default = 14,
              help = "Figure width in inches [default: %default]"),
  make_option("--height-per-facet", type = "double", dest = "height_per_facet",
              default = 2.8,
              help = "Figure height per facet row in inches [default: %default]"),
  make_option("--top-n-influential", type = "integer", dest = "top_n_influential",
              default = 5L,
              help = "Top N most influential studies per phenotype/contrast [default: %default]"),
  make_option("--drop-skipped", action = "store_true", dest = "drop_skipped",
              default = FALSE,
              help = "Drop rows with status != ok before plotting")
)

opt <- parse_args(OptionParser(option_list = option_list))

req <- c("input_root", "outdir")
miss <- req[vapply(req, function(x) is.null(opt[[x]]) || is.na(opt[[x]]) || !nzchar(opt[[x]]), logical(1))]
if (length(miss) > 0) {
  stop("Missing required arguments: ", paste(miss, collapse = ", "))
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

summary_files <- list.files(
  path = opt$input_root,
  pattern = opt$pattern,
  recursive = TRUE,
  full.names = TRUE
)

if (length(summary_files) == 0) {
  stop("No LOO summary files found under: ", opt$input_root)
}

read_one_summary <- function(f) {
  x <- fread(f, data.table = FALSE)

  rel_dir <- dirname(sub(paste0("^", normalizePath(opt$input_root, winslash = "/"), "/?"), "", normalizePath(f, winslash = "/")))
  rel_parts <- strsplit(rel_dir, "/", fixed = TRUE)[[1]]

  phenotype <- if (length(rel_parts) >= 1) rel_parts[1] else NA_character_
  contrast_dir <- if (length(rel_parts) >= 2) rel_parts[2] else NA_character_

  x$phenotype <- phenotype

  if (!"contrast" %in% colnames(x) || all(is.na(x$contrast) | x$contrast == "")) {
    x$contrast <- contrast_dir
  }

  x$source_file <- basename(f)
  x$source_dir <- dirname(f)
  x
}

df <- do.call(rbind, lapply(summary_files, read_one_summary))

if (opt$drop_skipped && "status" %in% colnames(df)) {
  df <- df[df$status == "ok", , drop = FALSE]
}

if (nrow(df) == 0) {
  stop("No rows remaining after filtering")
}

num_cols <- intersect(
  c(
    "pearson_log2fc", "spearman_stat", "jaccard_sig",
    "prop_direction_concordant", "median_abs_dlog2fc",
    "mean_abs_dlog2fc", "p90_abs_dlog2fc", "p95_abs_dlog2fc",
    "prop_full_sig_retained", "top_n_abs_stat_overlap",
    "betaConv_rate_loo", "median_dispersion",
    "median_dispersion_shrinkage", "median_dispersion_fit_error",
    "dispersion_slope", "prop_dispersion_outliers"
  ),
  colnames(df)
)

for (nm in num_cols) {
  df[[nm]] <- as.numeric(df[[nm]])
}

if (!all(c("phenotype", "contrast", "left_out_study") %in% colnames(df))) {
  stop("Combined table is missing one or more required columns: phenotype, contrast, left_out_study")
}

df$facet_label <- paste(df$phenotype, df$contrast, sep = " | ")

write.table(
  df,
  file = file.path(opt$outdir, "loo_summary_all.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

agg_mean <- aggregate(
  df[, num_cols, drop = FALSE],
  by = df[, c("phenotype", "contrast"), drop = FALSE],
  FUN = function(x) mean(x, na.rm = TRUE)
)

write.table(
  agg_mean,
  file = file.path(opt$outdir, "loo_summary_by_contrast_mean.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

agg_median <- aggregate(
  df[, num_cols, drop = FALSE],
  by = df[, c("phenotype", "contrast"), drop = FALSE],
  FUN = function(x) median(x, na.rm = TRUE)
)

write.table(
  agg_median,
  file = file.path(opt$outdir, "loo_summary_by_contrast_median.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

df_influence <- df

scale01 <- function(x) {
  rng <- range(x, na.rm = TRUE)
  if (!all(is.finite(rng)) || diff(rng) == 0) {
    return(rep(0, length(x)))
  }
  (x - rng[1]) / diff(rng)
}

if (all(c("pearson_log2fc", "jaccard_sig", "median_abs_dlog2fc") %in% colnames(df_influence))) {
  df_influence$influence_score <- 
    scale01(1 - df_influence$pearson_log2fc) +
    scale01(1 - df_influence$jaccard_sig) +
    scale01(df_influence$median_abs_dlog2fc)
} else {
  df_influence$influence_score <- NA_real_
}

top_influential <- do.call(
  rbind,
  lapply(split(df_influence, interaction(df_influence$phenotype, df_influence$contrast, drop = TRUE)), function(z) {
    z <- z[order(-z$influence_score, z$pearson_log2fc, z$jaccard_sig), , drop = FALSE]
    head(z, opt$top_n_influential)
  })
)

write.table(
  top_influential,
  file = file.path(opt$outdir, "top_influential_studies.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)

facet_levels <- unique(df$facet_label)
facet_levels <- facet_levels[order(facet_levels)]
df$facet_label <- factor(df$facet_label, levels = facet_levels)

order_within_facet <- do.call(
  rbind,
  lapply(split(df, df$facet_label), function(z) {
    z <- z[order(z$pearson_log2fc, z$median_abs_dlog2fc, decreasing = c(FALSE, TRUE)), , drop = FALSE]
    z$left_out_study_ordered <- factor(z$left_out_study, levels = z$left_out_study)
    z
  })
)

plot_df <- order_within_facet
n_facets <- length(unique(plot_df$facet_label))
fig_height <- max(6, n_facets * opt$height_per_facet)

theme_main <- theme_bw(base_size = 11) +
  theme(
    strip.background = element_rect(fill = "grey95", colour = "grey80"),
    panel.grid.major.y = element_blank(),
    panel.grid.minor = element_blank(),
    legend.position = "none"
  )

p1 <- ggplot(plot_df, aes(x = pearson_log2fc, y = left_out_study_ordered)) +
  geom_point(size = 2.5) +
  geom_segment(aes(x = 0, xend = pearson_log2fc, y = left_out_study_ordered, yend = left_out_study_ordered),
               linewidth = 0.4) +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 1) +
  labs(
    title = "Leave-one-out sensitivity of pooled RNA-seq differential-expression models",
    subtitle = "Panel A shows the primary robustness metric: Pearson correlation of gene-wise log2 fold changes between full and leave-one-out models",
    x = "Pearson correlation of log2FC (LOO vs full)",
    y = "Left-out study"
  ) +
  theme_main

p2 <- ggplot(plot_df, aes(x = median_abs_dlog2fc, y = left_out_study_ordered)) +
  geom_point(size = 2.5) +
  geom_segment(aes(x = 0, xend = median_abs_dlog2fc, y = left_out_study_ordered, yend = left_out_study_ordered),
               linewidth = 0.4) +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 1) +
  labs(
    x = "Median |Δlog2FC|",
    y = "Left-out study"
  ) +
  theme_main

p3 <- ggplot(plot_df, aes(x = jaccard_sig, y = left_out_study_ordered)) +
  geom_point(size = 2.5) +
  geom_segment(aes(x = 0, xend = jaccard_sig, y = left_out_study_ordered, yend = left_out_study_ordered),
               linewidth = 0.4) +
  facet_wrap(~ facet_label, scales = "free_y", ncol = 1) +
  labs(
    x = "Jaccard similarity of significant genes",
    y = "Left-out study",
    caption = "Significant genes defined using FDR < 0.05. Lower Pearson correlation and higher median |Δlog2FC| indicate stronger study influence."
  ) +
  theme_main

p_main <- p1 / p2 / p3 + plot_layout(heights = c(1.15, 1, 1))

ggsave(
  filename = file.path(opt$outdir, "main_loo_comparison_plot.pdf"),
  plot = p_main,
  width = opt$width,
  height = fig_height * 2.2
)

ggsave(
  filename = file.path(opt$outdir, "main_loo_comparison_plot.png"),
  plot = p_main,
  width = opt$width,
  height = fig_height * 2.2,
  dpi = 300
)

if ("prop_direction_concordant" %in% colnames(plot_df)) {
  p4 <- ggplot(plot_df, aes(x = prop_direction_concordant, y = left_out_study_ordered)) +
    geom_point(size = 2.5) +
    geom_segment(aes(x = 0, xend = prop_direction_concordant, y = left_out_study_ordered, yend = left_out_study_ordered),
                 linewidth = 0.4) +
    facet_wrap(~ facet_label, scales = "free_y", ncol = 1) +
    labs(
      x = "Proportion of concordant effect direction",
      y = "Left-out study"
    ) +
    theme_main

  ggsave(
    filename = file.path(opt$outdir, "supplement_direction_concordance.pdf"),
    plot = p4,
    width = opt$width,
    height = fig_height
  )
}

if ("median_dispersion_fit_error" %in% colnames(plot_df)) {
  p5 <- ggplot(plot_df, aes(x = median_dispersion_fit_error, y = left_out_study_ordered)) +
    geom_point(size = 2.5) +
    geom_segment(aes(x = 0, xend = median_dispersion_fit_error, y = left_out_study_ordered, yend = left_out_study_ordered),
                 linewidth = 0.4) +
    facet_wrap(~ facet_label, scales = "free_y", ncol = 1) +
    labs(
      x = "Median dispersion fit error",
      y = "Left-out study"
    ) +
    theme_main

  ggsave(
    filename = file.path(opt$outdir, "supplement_dispersion_fit_error.pdf"),
    plot = p5,
    width = opt$width,
    height = fig_height
  )
}