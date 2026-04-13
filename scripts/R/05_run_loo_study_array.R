#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(DESeq2)
  library(ggplot2)
  library(S4Vectors)
  library(BiocParallel)
})

option_list <- list(
  make_option("--dds-rds", type = "character", dest = "dds_rds",
              help = "Input DESeqDataSet RDS"),
  make_option("--outdir", type = "character", dest = "outdir",
              help = "Output directory"),
  make_option("--design", type = "character", dest = "design",
              default = "~ Study + Condition",
              help = "Design formula [default: %default]"),
  make_option("--contrast", action = "append", type = "character", dest = "contrast",
              default = NULL,
              help = "Contrast in format factor,levelA,levelB. Repeat flag for multiple contrasts"),
  make_option("--study-col", type = "character", dest = "study_col",
              default = "Study",
              help = "Study column in colData used for leave-one-out [default: %default]"),
  make_option("--leave-out-study", type = "character", dest = "leave_out_study",
              help = "Study level to leave out for this run"),
  make_option("--fit-type", type = "character", dest = "fit_type",
              default = "parametric",
              help = "DESeq2 fitType [default: %default]"),
  make_option("--sf-type", type = "character", dest = "sf_type",
              default = "ratio",
              help = "DESeq2 sfType [default: %default]"),
  make_option("--test", type = "character", dest = "test",
              default = "Wald",
              help = "DESeq2 test [default: %default]"),
  make_option("--padj-cutoff", type = "double", dest = "padj_cutoff",
              default = 0.05,
              help = "Adjusted p-value cutoff [default: %default]"),
  make_option("--lfc-threshold", type = "double", dest = "lfc_threshold",
              default = 0,
              help = "Optional absolute LFC threshold for direction concordance filtering [default: %default]"),
  make_option("--top-n", type = "integer", dest = "top_n",
              default = 100L,
              help = "Top N genes by absolute stat for overlap metric [default: %default]"),
  make_option("--min-replicates-for-replace", type = "double", dest = "min_replicates_for_replace",
              default = Inf,
              help = "DESeq2 minReplicatesForReplace [default: %default]"),
  make_option("--n-cores", type = "integer", dest = "n_cores",
              default = 1L,
              help = "Number of CPU cores for DESeq2 parallel fitting [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

req <- c("dds_rds", "outdir", "leave_out_study")
miss <- req[vapply(req, function(x) is.null(opt[[x]]) || is.na(opt[[x]]) || !nzchar(opt[[x]]), logical(1))]
if (length(miss) > 0) {
  stop("Missing required arguments: ", paste(miss, collapse = ", "))
}
if (is.null(opt$contrast) || length(opt$contrast) == 0) {
  stop("Provide at least one --contrast")
}
if (is.na(opt$n_cores) || opt$n_cores < 1) {
  stop("--n-cores must be >= 1")
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

dds <- readRDS(opt$dds_rds)
if (!inherits(dds, "DESeqDataSet")) stop("Input RDS is not a DESeqDataSet")
if (!opt$study_col %in% colnames(colData(dds))) {
  stop("Study column not found in colData: ", opt$study_col)
}

bp_param <- if (opt$n_cores > 1) {
  BiocParallel::MulticoreParam(workers = opt$n_cores, progressbar = FALSE)
} else {
  BiocParallel::SerialParam(progressbar = FALSE)
}

drop_unused_levels <- function(dds_obj) {
  cd <- as.data.frame(colData(dds_obj))
  cd[] <- lapply(cd, function(x) {
    if (is.factor(x)) droplevels(x) else x
  })
  colData(dds_obj) <- S4Vectors::DataFrame(cd)
  dds_obj
}

design_estimable <- function(dds_obj, design_formula) {
  vars <- all.vars(as.formula(design_formula))
  for (v in vars) {
    if (!v %in% colnames(colData(dds_obj))) {
      return(list(ok = FALSE, reason = paste("design variable missing:", v)))
    }
    x <- colData(dds_obj)[[v]]
    if (is.factor(x) && nlevels(droplevels(x)) < 2) {
      return(list(ok = FALSE, reason = paste("factor has <2 levels after subsetting:", v)))
    }
  }
  list(ok = TRUE, reason = NA_character_)
}

parse_contrast <- function(ct) {
  parts <- trimws(strsplit(ct, ",", fixed = TRUE)[[1]])
  if (length(parts) != 3) {
    stop("Invalid contrast: ", ct, " ; expected factor,levelA,levelB")
  }
  list(
    factor_name = parts[1],
    level_a = parts[2],
    level_b = parts[3],
    label = paste(parts[1], parts[2], "vs", parts[3], sep = "_")
  )
}

fit_dds <- function(dds_obj, design_formula, test, fit_type, sf_type,
                    min_replicates_for_replace, bp_param, n_cores) {
  design(dds_obj) <- as.formula(design_formula)
  DESeq(
    dds_obj,
    test = test,
    fitType = fit_type,
    sfType = sf_type,
    minReplicatesForReplace = min_replicates_for_replace,
    parallel = (n_cores > 1),
    BPPARAM = bp_param
  )
}

extract_res_df <- function(dds_fit, factor_name, level_a, level_b) {
  res <- results(dds_fit, contrast = c(factor_name, level_a, level_b))
  res_df <- as.data.frame(res)
  res_df$feature_id <- rownames(res_df)
  res_df <- res_df[, c("feature_id", setdiff(colnames(res_df), "feature_id"))]
  res_df$contrast_factor <- factor_name
  res_df$contrast_numerator <- level_a
  res_df$contrast_denominator <- level_b
  res_df$contrast_description <- paste0(level_a, " vs ", level_b)
  res_df$log2fc_definition <- paste0(level_a, " - ", level_b)
  res_df
}

jaccard_index <- function(x, y) {
  ux <- unique(x)
  uy <- unique(y)
  u <- union(ux, uy)
  if (length(u) == 0) return(NA_real_)
  length(intersect(ux, uy)) / length(u)
}

beta_conv_rate <- function(dds_fit) {
  mm <- as.data.frame(mcols(dds_fit))
  if ("betaConv" %in% colnames(mm)) {
    return(mean(mm$betaConv %in% TRUE, na.rm = TRUE))
  }
  NA_real_
}

extract_dispersion_metrics <- function(dds_fit) {
  mm <- as.data.frame(mcols(dds_fit))

  if (!all(c("dispGeneEst", "dispFit", "baseMean") %in% colnames(mm))) {
    return(list(
      median_dispersion = NA_real_,
      median_dispersion_shrinkage = NA_real_,
      median_dispersion_fit_error = NA_real_,
      dispersion_slope = NA_real_,
      prop_dispersion_outliers = NA_real_
    ))
  }

  dg <- mm$dispGeneEst
  dfit <- mm$dispFit
  bm <- mm$baseMean
  dmap <- dispersions(dds_fit)

  valid_shrink <- !is.na(dg) & !is.na(dmap) & dg > 0 & dmap > 0
  valid_fit <- !is.na(dg) & !is.na(dfit) & dg > 0 & dfit > 0
  valid_slope <- !is.na(bm) & !is.na(dmap) & bm > 0 & dmap > 0

  median_dispersion <- if (any(!is.na(dmap) & dmap > 0)) median(dmap[dmap > 0], na.rm = TRUE) else NA_real_
  median_dispersion_shrinkage <- if (sum(valid_shrink) > 0) median(abs(log10(dg[valid_shrink]) - log10(dmap[valid_shrink])), na.rm = TRUE) else NA_real_
  median_dispersion_fit_error <- if (sum(valid_fit) > 0) median(abs(log10(dg[valid_fit]) - log10(dfit[valid_fit])), na.rm = TRUE) else NA_real_
  dispersion_slope <- if (sum(valid_slope) > 10) unname(coef(lm(log10(dmap[valid_slope]) ~ log10(bm[valid_slope])))[2]) else NA_real_
  prop_dispersion_outliers <- if (sum(valid_fit) > 0) mean(abs(log10(dg[valid_fit]) - log10(dfit[valid_fit])) > 1, na.rm = TRUE) else NA_real_

  list(
    median_dispersion = median_dispersion,
    median_dispersion_shrinkage = median_dispersion_shrinkage,
    median_dispersion_fit_error = median_dispersion_fit_error,
    dispersion_slope = dispersion_slope,
    prop_dispersion_outliers = prop_dispersion_outliers
  )
}

top_n_overlap_abs_stat <- function(res_a, res_b, n = 100L) {
  a <- res_a[!is.na(res_a$stat), c("feature_id", "stat"), drop = FALSE]
  b <- res_b[!is.na(res_b$stat), c("feature_id", "stat"), drop = FALSE]
  if (nrow(a) == 0 || nrow(b) == 0) return(NA_real_)
  a <- a[order(abs(a$stat), decreasing = TRUE), , drop = FALSE]
  b <- b[order(abs(b$stat), decreasing = TRUE), , drop = FALSE]
  top_a <- head(a$feature_id, n)
  top_b <- head(b$feature_id, n)
  if (length(top_a) == 0 || length(top_b) == 0) return(NA_real_)
  length(intersect(top_a, top_b)) / min(length(top_a), length(top_b))
}

safe_quantile <- function(x, prob) {
  x <- x[is.finite(x)]
  if (length(x) == 0) return(NA_real_)
  unname(quantile(x, probs = prob, na.rm = TRUE))
}

studies <- unique(as.character(colData(dds)[[opt$study_col]]))
studies <- studies[!is.na(studies) & studies != ""]
if (!opt$leave_out_study %in% studies) {
  stop("--leave-out-study not found in ", opt$study_col, ": ", opt$leave_out_study)
}

full_fit <- fit_dds(
  dds_obj = dds,
  design_formula = opt$design,
  test = opt$test,
  fit_type = opt$fit_type,
  sf_type = opt$sf_type,
  min_replicates_for_replace = opt$min_replicates_for_replace,
  bp_param = bp_param,
  n_cores = opt$n_cores
)

full_disp_metrics <- extract_dispersion_metrics(full_fit)
saveRDS(full_fit, file.path(opt$outdir, "dds_full_fitted.rds"))

capture.output({
  cat("leave_out_study:\t", opt$leave_out_study, "\n", sep = "")
  cat("study_col:\t", opt$study_col, "\n", sep = "")
  cat("design:\t", paste(deparse(design(full_fit)), collapse = " "), "\n", sep = "")
  cat("n_cores:\t", opt$n_cores, "\n", sep = "")
  cat("minReplicatesForReplace:\t", opt$min_replicates_for_replace, "\n", sep = "")
  cat("n_samples_full:\t", ncol(full_fit), "\n", sep = "")
  cat("n_genes_full:\t", nrow(full_fit), "\n", sep = "")
  cat("betaConv_rate_full:\t", beta_conv_rate(full_fit), "\n", sep = "")
  cat("median_dispersion_full:\t", full_disp_metrics$median_dispersion, "\n", sep = "")
  cat("median_dispersion_fit_error_full:\t", full_disp_metrics$median_dispersion_fit_error, "\n", sep = "")
}, file = file.path(opt$outdir, paste0("run_summary_", opt$leave_out_study, ".txt")))

keep <- as.character(colData(dds)[[opt$study_col]]) != opt$leave_out_study
dds_loo <- dds[, keep]
dds_loo <- drop_unused_levels(dds_loo)

chk <- design_estimable(dds_loo, opt$design)
if (!chk$ok) {
  stop("LOO design not estimable after leaving out ", opt$leave_out_study, ": ", chk$reason)
}

loo_fit <- fit_dds(
  dds_obj = dds_loo,
  design_formula = opt$design,
  test = opt$test,
  fit_type = opt$fit_type,
  sf_type = opt$sf_type,
  min_replicates_for_replace = opt$min_replicates_for_replace,
  bp_param = bp_param,
  n_cores = opt$n_cores
)

loo_disp_metrics <- extract_dispersion_metrics(loo_fit)
saveRDS(loo_fit, file.path(opt$outdir, paste0("dds_loo_", opt$leave_out_study, "_fitted.rds")))

summary_list <- list()

for (ct in opt$contrast) {
  cc <- parse_contrast(ct)
  factor_name <- cc$factor_name
  level_a <- cc$level_a
  level_b <- cc$level_b
  label <- cc$label

  message("Running contrast: ", label, " ; leaving out study: ", opt$leave_out_study)

  full_levels <- unique(as.character(colData(full_fit)[[factor_name]]))
  if (!all(c(level_a, level_b) %in% full_levels)) {
    stop("Contrast levels not present in full model for ", label)
  }

  loo_levels <- unique(as.character(colData(loo_fit)[[factor_name]]))
  if (!all(c(level_a, level_b) %in% loo_levels)) {
    warning("Contrast levels absent after leaving out study for ", label, ". Skipping.")
    next
  }

  contrast_dir <- file.path(opt$outdir, label)
  dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)

  full_res_df <- extract_res_df(full_fit, factor_name, level_a, level_b)
  loo_res_df  <- extract_res_df(loo_fit,  factor_name, level_a, level_b)

  write.table(full_res_df,
              file.path(contrast_dir, sprintf("results_full_%s.tsv", label)),
              sep = "\t", quote = FALSE, row.names = FALSE)
  write.table(loo_res_df,
              file.path(contrast_dir, sprintf("results_loo_%s_%s.tsv", opt$leave_out_study, label)),
              sep = "\t", quote = FALSE, row.names = FALSE)

  cmp <- merge(
    full_res_df[, c("feature_id", "log2FoldChange", "stat", "padj")],
    loo_res_df[,  c("feature_id", "log2FoldChange", "stat", "padj")],
    by = "feature_id",
    suffixes = c("_full", "_loo"),
    all = FALSE
  )

  cmp <- cmp[!is.na(cmp$log2FoldChange_full) & !is.na(cmp$log2FoldChange_loo), , drop = FALSE]
  cmp$abs_dlog2fc <- abs(cmp$log2FoldChange_loo - cmp$log2FoldChange_full)
  cmp$direction_full <- sign(cmp$log2FoldChange_full)
  cmp$direction_loo <- sign(cmp$log2FoldChange_loo)

  cmp_dir <- if (opt$lfc_threshold > 0) {
    cmp[abs(cmp$log2FoldChange_full) >= opt$lfc_threshold |
          abs(cmp$log2FoldChange_loo) >= opt$lfc_threshold, , drop = FALSE]
  } else {
    cmp
  }

  prop_concordant <- if (nrow(cmp_dir) > 0) mean(cmp_dir$direction_full == cmp_dir$direction_loo, na.rm = TRUE) else NA_real_
  pearson_lfc <- if (nrow(cmp) > 1) suppressWarnings(cor(cmp$log2FoldChange_full, cmp$log2FoldChange_loo, method = "pearson", use = "complete.obs")) else NA_real_

  stat_ok <- !is.na(cmp$stat_full) & !is.na(cmp$stat_loo)
  spearman_stat <- if (sum(stat_ok) > 1) suppressWarnings(cor(cmp$stat_full[stat_ok], cmp$stat_loo[stat_ok], method = "spearman", use = "complete.obs")) else NA_real_

  full_sig <- full_res_df$feature_id[!is.na(full_res_df$padj) & full_res_df$padj < opt$padj_cutoff]
  loo_sig  <- loo_res_df$feature_id[!is.na(loo_res_df$padj) & loo_res_df$padj < opt$padj_cutoff]

  jaccard_sig <- jaccard_index(full_sig, loo_sig)
  prop_full_sig_retained <- if (length(full_sig) > 0) mean(full_sig %in% loo_sig) else NA_real_
  top_overlap <- top_n_overlap_abs_stat(full_res_df, loo_res_df, n = opt$top_n)

  cmp$contrast <- label
  cmp$left_out_study <- opt$leave_out_study

  write.table(
    cmp,
    file.path(contrast_dir, sprintf("comparison_full_vs_loo_%s_%s.tsv", opt$leave_out_study, label)),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  pdf(file.path(contrast_dir, sprintf("MA_full_%s.pdf", label)), width = 7, height = 6)
  plotMA(results(full_fit, contrast = c(factor_name, level_a, level_b)),
         main = paste0("Full: ", level_a, " vs ", level_b),
         ylim = c(-5, 5))
  dev.off()

  pdf(file.path(contrast_dir, sprintf("MA_loo_%s_%s.pdf", opt$leave_out_study, label)), width = 7, height = 6)
  plotMA(results(loo_fit, contrast = c(factor_name, level_a, level_b)),
         main = paste0("LOO ", opt$leave_out_study, ": ", level_a, " vs ", level_b),
         ylim = c(-5, 5))
  dev.off()

  summary_list[[label]] <- data.frame(
    contrast = label,
    left_out_study = opt$leave_out_study,
    n_samples_loo = ncol(loo_fit),
    n_genes_loo = nrow(loo_fit),
    pearson_log2fc = pearson_lfc,
    spearman_stat = spearman_stat,
    jaccard_sig = jaccard_sig,
    prop_direction_concordant = prop_concordant,
    median_abs_dlog2fc = median(cmp$abs_dlog2fc, na.rm = TRUE),
    mean_abs_dlog2fc = mean(cmp$abs_dlog2fc, na.rm = TRUE),
    p90_abs_dlog2fc = safe_quantile(cmp$abs_dlog2fc, 0.90),
    p95_abs_dlog2fc = safe_quantile(cmp$abs_dlog2fc, 0.95),
    n_sig_full = length(full_sig),
    n_sig_loo = length(loo_sig),
    prop_full_sig_retained = prop_full_sig_retained,
    top_n_abs_stat_overlap = top_overlap,
    n_genes_compared = nrow(cmp),
    n_genes_tested_loo = sum(!is.na(loo_res_df$pvalue)),
    n_genes_padj_non_na_loo = sum(!is.na(loo_res_df$padj)),
    betaConv_rate_loo = beta_conv_rate(loo_fit),
    median_dispersion_loo = loo_disp_metrics$median_dispersion,
    median_dispersion_shrinkage_loo = loo_disp_metrics$median_dispersion_shrinkage,
    median_dispersion_fit_error_loo = loo_disp_metrics$median_dispersion_fit_error,
    dispersion_slope_loo = loo_disp_metrics$dispersion_slope,
    prop_dispersion_outliers_loo = loo_disp_metrics$prop_dispersion_outliers,
    stringsAsFactors = FALSE
  )
}

if (length(summary_list) > 0) {
  summary_df <- do.call(rbind, summary_list)
  write.table(
    summary_df,
    file.path(opt$outdir, paste0("loo_summary_", opt$leave_out_study, ".tsv")),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}