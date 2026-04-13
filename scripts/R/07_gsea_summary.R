#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(readr)
  library(dplyr)
  library(stringr)
  library(ggplot2)
  library(forcats)
})

option_list <- list(
  make_option("--go", type = "character", help = "Path to GOBP readable GSEA TSV"),
  make_option("--kegg", type = "character", help = "Path to KEGG readable GSEA TSV"),
  make_option("--outdir", type = "character", default = "gsea_figures", help = "Output directory [default: %default]"),
  make_option("--fdr", type = "double", default = 0.05, help = "FDR threshold [default: %default]"),
  make_option("--top-kegg", dest = "top_kegg", type = "integer", default = 10, help = "Number of top KEGG terms per direction [default: %default]"),
  make_option("--top-go", dest = "top_go", type = "integer", default = 12, help = "Number of top GO terms per direction [default: %default]"),
  make_option("--title", type = "character", default = "acute exercise (Post vs Pre)",
            help = "Main title suffix [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

opt$top_kegg <- as.integer(if (is.null(opt$top_kegg) || is.na(opt$top_kegg)) 10 else opt$top_kegg)
opt$top_go   <- as.integer(if (is.null(opt$top_go) || is.na(opt$top_go)) 12 else opt$top_go)

if (is.null(opt$go) || is.null(opt$kegg)) {
  stop("Both --go and --kegg must be supplied.")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

theme_set(theme_bw(base_size = 11))

plot_theme <- theme(
  panel.grid.minor = element_blank(),
  panel.grid.major.y = element_line(linewidth = 0.2),
  plot.title = element_text(face = "bold"),
  plot.subtitle = element_text(size = 10),
  axis.title = element_text(face = "bold"),
  strip.text = element_text(face = "bold"),
  legend.title = element_text(face = "bold"),
  legend.position = "right"
)

read_gsea <- function(path, db_name, fdr_cutoff = 0.05) {
  read_tsv(path, show_col_types = FALSE) %>%
    mutate(
      database = db_name,
      regulation = if_else(NES >= 0, "Activated (up-regulated)", "Suppressed (down-regulated)"),
      regulation = factor(
        regulation,
        levels = c("Suppressed (down-regulated)", "Activated (up-regulated)")
      ),
      log10_fdr = -log10(p.adjust),
      core_n = if_else(is.na(core_enrichment) | core_enrichment == "", NA_integer_, str_count(core_enrichment, "/") + 1L),
      significant = p.adjust < fdr_cutoff,
      Description = str_squish(Description)
    )
}

top_terms <- function(df, n, regulation_keep) {
  n <- as.integer(if (is.null(n) || length(n) == 0 || is.na(n)) 10 else n)

  df %>%
    filter(significant, regulation == regulation_keep) %>%
    arrange(
      if (regulation_keep == "Activated (up-regulated)") desc(NES) else NES,
      p.adjust
    ) %>%
    slice_head(n = n) %>%
    mutate(Description = str_wrap(Description, width = 44))
}

make_dotplot <- function(df, top_n, title_txt, subtitle_txt) {
  plot_df <- bind_rows(
    top_terms(df, top_n, "Activated (up-regulated)"),
    top_terms(df, top_n, "Suppressed (down-regulated)")
  ) %>%
    mutate(
      Description = fct_reorder(
        Description,
        if_else(regulation == "Activated (up-regulated)", NES, -abs(NES))
      )
    )

  ggplot(plot_df, aes(x = NES, y = Description)) +
    geom_vline(xintercept = 0, linetype = 2, linewidth = 0.4) +
    geom_point(aes(size = log10_fdr, colour = regulation), alpha = 0.9) +
    scale_colour_manual(
      values = c(
        "Suppressed (down-regulated)" = "blue",
        "Activated (up-regulated)" = "red"
      )
    ) +
    labs(
      title = title_txt,
      subtitle = subtitle_txt,
      x = "GSEA NES",
      y = NULL,
      size = expression(-log[10]("FDR")),
      colour = NULL
    ) +
    plot_theme
}

go <- read_gsea(opt$go, "GOBP", opt$fdr)
kegg <- read_gsea(opt$kegg, "KEGG", opt$fdr)

# Export summary tables
write_tsv(
  go %>% filter(significant) %>% arrange(desc(NES)),
  file.path(opt$outdir, "gobp_significant_terms.tsv")
)

write_tsv(
  kegg %>% filter(significant) %>% arrange(desc(NES)),
  file.path(opt$outdir, "kegg_significant_terms.tsv")
)

write_tsv(
  bind_rows(
    top_terms(kegg, opt$top_kegg, "Activated (up-regulated)"),
    top_terms(kegg, opt$top_kegg, "Suppressed (down-regulated)")
  ) %>% select(database, Description, NES, p.adjust, regulation, core_enrichment),
  file.path(opt$outdir, "top_kegg_dotplot_terms.tsv")
)

write_tsv(
  bind_rows(
    top_terms(go, opt$top_go, "Activated (up-regulated)"),
    top_terms(go, opt$top_go, "Suppressed (down-regulated)")
  ) %>% select(database, Description, NES, p.adjust, regulation, core_enrichment),
  file.path(opt$outdir, "top_gobp_dotplot_terms.tsv")
)

# Plot 1: KEGG dotplot
p_kegg <- make_dotplot(
  kegg,
  top_n = opt$top_kegg,
  title_txt = paste0("KEGG GSEA:", opt$title),
  subtitle_txt = "Activated pathways in red and suppressed pathways in blue"
)

ggsave(
  filename = file.path(opt$outdir, "figure1_kegg_dotplot.png"),
  plot = p_kegg,
  width = 9.5,
  height = 6.5,
  dpi = 300
)

# Plot 2: GOBP dotplot
p_go <- make_dotplot(
  go,
  top_n = opt$top_go,
  title_txt = paste0("GOBP GSEA: ", opt$title),
  subtitle_txt = "Activated terms in red and suppressed terms in blue"
)

ggsave(
  filename = file.path(opt$outdir, "figure2_gobp_dotplot.png"),
  plot = p_go,
  width = 10,
  height = 8,
  dpi = 300
)