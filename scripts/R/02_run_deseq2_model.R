#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(ggplot2)
  library(pheatmap)
  library(RColorBrewer)
})

option_list <- list(
  make_option("--dds-rds", type = "character", dest = "dds_rds", help = "Input DESeqDataSet RDS"),
  make_option("--outdir", type = "character", dest = "outdir", help = "Output directory"),
  make_option("--design", type = "character", dest = "design", default = "~ Study + Condition", help = "Design formula [default: %default]"),
  make_option("--contrast", , dest = "contrast", action = "append", type = "character", default = NULL,
              help = "Contrast in format factor,levelA,levelB. Repeat flag for multiple contrasts"),
  make_option("--fit-type", type = "character", dest = "fit_type", default = "parametric", help = "DESeq2 fitType [default: %default]"),
  make_option("--sf-type", type = "character", dest = "sf_type", default = "ratio", help = "DESeq2 sfType [default: %default]"),
  make_option("--test", type = "character", dest = "test", default = "Wald", help = "DESeq2 test [default: %default]"),
  make_option("--min-count-pca", type = "integer", dest = "min_count_pca", default = 10L, help = "Minimum count used for rlog PCA subset [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

req <- c("dds_rds", "outdir")
miss <- req[vapply(req, function(x) is.null(opt[[x]]) || is.na(opt[[x]]), logical(1))]
if (length(miss) > 0) stop("Missing required arguments: ", paste(miss, collapse = ", "))
if (is.null(opt$contrast) || length(opt$contrast) == 0) stop("Provide at least one --contrast")

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
dds <- readRDS(opt$dds_rds)
if (!inherits(dds, "DESeqDataSet")) stop("Input RDS is not a DESeqDataSet")

design(dds) <- as.formula(opt$design)
dds <- DESeq(dds, test = opt$test, fitType = opt$fit_type, sfType = opt$sf_type)
saveRDS(dds, file.path(opt$outdir, paste0(gsub("dds_|.rds", "", basename(opt$dds_rds)), "_dds_fitted.rds")))

capture.output({
  cat("design:\t", paste(deparse(design(dds)), collapse = " "), "\n", sep = "")
  cat("test:\t", opt$test, "\n", sep = "")
  cat("fitType:\t", opt$fit_type, "\n", sep = "")
  cat("sfType:\t", opt$sf_type, "\n", sep = "")
  print(resultsNames(dds))
  print(summary(sizeFactors(dds)))
  print(mcols(dds)$dispGeneEst[1:5])
}, file = file.path(opt$outdir, "model_fit_summary.txt"))

write.table(
  data.frame(sample = colnames(dds), sizeFactor = sizeFactors(dds)),
  file.path(opt$outdir, "size_factors.tsv"), sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  as.data.frame(mcols(dds)),
  file.path(opt$outdir, "row_metrics.tsv"), sep = "\t", quote = FALSE, row.names = TRUE
)

pdf(file.path(opt$outdir, "dispersion_plot.pdf"), width = 7, height = 6)
plotDispEsts(dds)
dev.off()

vsd <- vst(dds, blind = FALSE)
sd_mat <- dist(t(assay(vsd)))
sd_mat <- as.matrix(sd_mat)
ann <- as.data.frame(colData(dds))

pdf(file.path(opt$outdir, "sample_distance_heatmap.pdf"), width = 9, height = 8)
pheatmap(sd_mat,
         clustering_distance_rows = "correlation",
         clustering_distance_cols = "correlation")
dev.off()

pca_groups <- intersect(c("PHENOTYPE", "Study", "Condition"), colnames(colData(dds)))
if (length(pca_groups) == 0) pca_groups <- colnames(colData(dds))[1]
pca_df <- plotPCA(vsd, intgroup = pca_groups, returnData = TRUE)
percentVar <- round(100 * attr(pca_df, "percentVar"))
color_var <- if ("Condition" %in% colnames(pca_df)) "Condition" else pca_groups[1]
p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[color_var]])) +
  geom_point(size = 3) +
  xlab(sprintf("PC1: %s%% variance", percentVar[1])) +
  ylab(sprintf("PC2: %s%% variance", percentVar[2])) +
  theme_bw()
ggsave(file.path(opt$outdir, "pca_plot.pdf"), plot = p, width = 7, height = 6)
write.table(pca_df, file.path(opt$outdir, "pca_coordinates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

norm_counts <- counts(dds, normalized = TRUE)
con <- gzfile(file.path(opt$outdir, "normalized_counts.tsv.gz"), open = "wt")
write.table(as.data.frame(norm_counts), con, sep = "\t", quote = FALSE, row.names = TRUE)
close(con)

for (ct in opt$contrast) {
  parts <- strsplit(ct, ",", fixed = TRUE)[[1]]
  if (length(parts) != 3) stop("Invalid contrast: ", ct, " ; expected factor,levelA,levelB")
  factor_name <- parts[1]
  level_a <- parts[2]
  level_b <- parts[3]
  res <- results(dds, contrast = c(factor_name, level_a, level_b))
  res_df <- as.data.frame(res)
  res_df$feature_id <- rownames(res_df)
  res_df <- res_df[, c("feature_id", setdiff(colnames(res_df), "feature_id"))]

  label <- paste(factor_name, level_a, "vs", level_b, sep = "_")
  write.table(res_df, file.path(opt$outdir, sprintf("results_%s.tsv", label)), sep = "\t", quote = FALSE, row.names = FALSE)

  pdf(file.path(opt$outdir, sprintf("MA_%s.pdf", label)), width = 7, height = 6)
  plotMA(res, ylim = c(-5, 5))
  dev.off()

  pdf(file.path(opt$outdir, sprintf("pvalue_hist_%s.pdf", label)), width = 7, height = 6)
  hist(res$pvalue, breaks = 50, main = label, xlab = "p-value")
  dev.off()
}