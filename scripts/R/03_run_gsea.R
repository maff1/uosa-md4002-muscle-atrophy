#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(data.table)
  library(clusterProfiler)
  library(org.Hs.eg.db)
  library(enrichplot)
  library(ggplot2)
})

option_list <- list(
  make_option("--results-tsv", type = "character", dest = "results_tsv", help = "Input DESeq2 results TSV"),
  make_option("--outdir", type = "character", dest = "outdir", help = "Output directory"),
  make_option("--id-col", type = "character", dest = "id_col", default = "feature_id", help = "Identifier column [default: %default]"),
  make_option("--rank-col", type = "character", dest = "rank_col", default = "stat", help = "Ranking column [default: %default]"),
  make_option("--id-type", type = "character", dest = "id_type", default = "SYMBOL", help = "Input identifier type for bitr [default: %default]"),
  make_option("--min-gs-size", type = "integer", dest = "min_gs_size", default = 10L, help = "minGSSize [default: %default]"),
  make_option("--max-gs-size", type = "integer", dest = "max_gs_size", default = 500L, help = "maxGSSize [default: %default]"),
  make_option("--pvalue-cutoff", type = "double", dest = "pvalue_cutoff", default = 0.05, help = "pvalueCutoff [default: %default]"),
  make_option("--padj-cutoff", type = "double", dest = "padj_cutoff", default = 0.05, help = "Adjusted p-value cutoff used for plotting [default: %default]")
)
opt <- parse_args(OptionParser(option_list = option_list))

req <- c("results_tsv", "outdir")
miss <- req[vapply(req, function(x) is.null(opt[[x]]) || is.na(opt[[x]]), logical(1))]
if (length(miss) > 0) stop("Missing required arguments: ", paste(miss, collapse = ", "))

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)
res <- fread(opt$results_tsv, data.table = FALSE)
req_cols <- c(opt$id_col, opt$rank_col)
if (!all(req_cols %in% colnames(res))) stop("Missing required columns in results table: ", paste(setdiff(req_cols, colnames(res)), collapse = ", "))

res <- res[!is.na(res[[opt$rank_col]]) & !is.na(res[[opt$id_col]]), , drop = FALSE]
res[[opt$id_col]] <- as.character(res[[opt$id_col]])
res <- res[res[[opt$id_col]] != "", , drop = FALSE]

map <- bitr(unique(res[[opt$id_col]]), fromType = opt$id_type, toType = c("ENTREZID", "SYMBOL"), OrgDb = org.Hs.eg.db)
map <- unique(map)
res_map <- merge(res, map, by.x = opt$id_col, by.y = opt$id_type, all.x = FALSE, all.y = FALSE)
res_map <- res_map[!duplicated(res_map$ENTREZID), , drop = FALSE]

rank_df <- res_map[, c("ENTREZID", opt$rank_col)]
rank_df <- rank_df[order(rank_df[[opt$rank_col]], decreasing = TRUE), , drop = FALSE]
geneList <- rank_df[[opt$rank_col]]
names(geneList) <- rank_df$ENTREZID
geneList <- sort(geneList, decreasing = TRUE)

write.table(res_map, file.path(opt$outdir, "mapped_input_results.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

run_and_save <- function(obj, prefix) {
  if (is.null(obj) || nrow(as.data.frame(obj)) == 0) return(invisible(NULL))
  df <- as.data.frame(obj)
  write.table(df, file.path(opt$outdir, sprintf("%s.tsv", prefix)), sep = "\t", quote = FALSE, row.names = FALSE)

  df_sig <- df[df$p.adjust <= opt$padj_cutoff, , drop = FALSE]
  if (nrow(df_sig) == 0) return(invisible(NULL))

  pdf(file.path(opt$outdir, sprintf("dotplot_%s.pdf", prefix)), width = 9, height = 7)
  print(dotplot(obj, showCategory = min(20, nrow(df_sig)), split = NULL))
  dev.off()

  pdf(file.path(opt$outdir, sprintf("emapplot_%s.pdf", prefix)), width = 10, height = 8)
  print(emapplot(pairwise_termsim(obj), showCategory = min(30, nrow(df_sig))))
  dev.off()

  ledge <- setReadable(obj, OrgDb = org.Hs.eg.db, keyType = "ENTREZID")
  ledge_df <- as.data.frame(ledge)
  write.table(ledge_df, file.path(opt$outdir, sprintf("%s_readable.tsv", prefix)), sep = "\t", quote = FALSE, row.names = FALSE)
}

gsea_bp <- gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "BP", minGSSize = opt$min_gs_size, maxGSSize = opt$max_gs_size, pvalueCutoff = opt$pvalue_cutoff, verbose = FALSE)
gsea_mf <- gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "MF", minGSSize = opt$min_gs_size, maxGSSize = opt$max_gs_size, pvalueCutoff = opt$pvalue_cutoff, verbose = FALSE)
gsea_cc <- gseGO(geneList = geneList, OrgDb = org.Hs.eg.db, keyType = "ENTREZID", ont = "CC", minGSSize = opt$min_gs_size, maxGSSize = opt$max_gs_size, pvalueCutoff = opt$pvalue_cutoff, verbose = FALSE)
gsea_kegg <- gseKEGG(geneList = geneList, organism = "hsa", keyType = "ncbi-geneid", minGSSize = opt$min_gs_size, maxGSSize = opt$max_gs_size, pvalueCutoff = opt$pvalue_cutoff, verbose = FALSE)

run_and_save(gsea_bp, "gsea_go_bp")
run_and_save(gsea_mf, "gsea_go_mf")
run_and_save(gsea_cc, "gsea_go_cc")
run_and_save(gsea_kegg, "gsea_kegg")