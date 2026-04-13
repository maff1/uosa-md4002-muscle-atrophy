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
              help = "DESeq2 minReplicatesForReplace. Use Inf to disable outlier replacement/refit [default: %default]"),
  make_option("--n-cores", type = "integer", dest = "n_cores",
              default = 1L,
              help = "Number of CPU cores for DESeq2 parallel fitting [default: %default]")
)

opt <- parse_args(OptionParser(option_list = option_list))

req <- c("dds_rds", "outdir")
miss <- req[vapply(req, function(x) is.null(opt[[x]]) || is.na(opt[[x]]), logical(1))]
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

  median_dispersion <- if (any(!is.na(dmap) & dmap > 0)) {
    median(dmap[dmap > 0], na.rm = TRUE)
  } else {
    NA_real_
  }

  median_dispersion_shrinkage <- if (sum(valid_shrink) > 0) {
    median(abs(log10(dg[valid_shrink]) - log10(dmap[valid_shrink])), na.rm = TRUE)
  } else {
    NA_real_
  }

  median_dispersion_fit_error <- if (sum(valid_fit) > 0) {
    median(abs(log10(dg[valid_fit]) - log10(dfit[valid_fit])), na.rm = TRUE)
  } else {
    NA_real_
  }

  dispersion_slope <- if (sum(valid_slope) > 10) {
    unname(coef(lm(log10(dmap[valid_slope]) ~ log10(bm[valid_slope])))[2])
  } else {
    NA_real_
  }

  prop_dispersion_outliers <- if (sum(valid_fit) > 0) {
    mean(abs(log10(dg[valid_fit]) - log10(dfit[valid_fit])) > 1, na.rm = TRUE)
  } else {
    NA_real_
  }

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

saveRDS(
  full_fit,
  file.path(opt$outdir, paste0(gsub("dds_|\\.rds$", "", basename(opt$dds_rds)), "_dds_full_fitted.rds"))
)

capture.output({
  cat("design:\t", paste(deparse(design(full_fit)), collapse = " "), "\n", sep = "")
  cat("test:\t", opt$test, "\n", sep = "")
  cat("fitType:\t", opt$fit_type, "\n", sep = "")
  cat("sfType:\t", opt$sf_type, "\n", sep = "")
  cat("study_col:\t", opt$study_col, "\n", sep = "")
  cat("minReplicatesForReplace:\t", opt$min_replicates_for_replace, "\n", sep = "")
  cat("n_cores:\t", opt$n_cores, "\n", sep = "")
  cat("n_samples:\t", ncol(full_fit), "\n", sep = "")
  cat("n_genes:\t", nrow(full_fit), "\n", sep = "")
  print(resultsNames(full_fit))
  print(summary(sizeFactors(full_fit)))
  cat("betaConv_rate:\t", beta_conv_rate(full_fit), "\n", sep = "")
  cat("median_dispersion:\t", full_disp_metrics$median_dispersion, "\n", sep = "")
  cat("median_dispersion_shrinkage:\t", full_disp_metrics$median_dispersion_shrinkage, "\n", sep = "")
  cat("median_dispersion_fit_error:\t", full_disp_metrics$median_dispersion_fit_error, "\n", sep = "")
  cat("dispersion_slope:\t", full_disp_metrics$dispersion_slope, "\n", sep = "")
  cat("prop_dispersion_outliers:\t", full_disp_metrics$prop_dispersion_outliers, "\n", sep = "")
}, file = file.path(opt$outdir, "full_model_summary.txt"))

pdf(file.path(opt$outdir, "dispersion_plot_full_model.pdf"), width = 7, height = 6)
plotDispEsts(full_fit)
dev.off()

studies <- unique(as.character(colData(dds)[[opt$study_col]]))
studies <- studies[!is.na(studies) & studies != ""]

for (ct in opt$contrast) {
  cc <- parse_contrast(ct)
  factor_name <- cc$factor_name
  level_a <- cc$level_a
  level_b <- cc$level_b
  label <- cc$label

  message("Running full model contrast: ", level_a, " vs ", level_b,
          " (log2FC = ", level_a, " - ", level_b, ")")

  if (!factor_name %in% colnames(colData(full_fit))) {
    stop("Contrast factor not found in colData: ", factor_name)
  }

  avail_levels <- unique(as.character(colData(full_fit)[[factor_name]]))
  if (!all(c(level_a, level_b) %in% avail_levels)) {
    stop("Contrast levels not present in full model for ", label)
  }

  contrast_dir <- file.path(opt$outdir, label)
  dir.create(contrast_dir, recursive = TRUE, showWarnings = FALSE)

  full_res_df <- extract_res_df(full_fit, factor_name, level_a, level_b)
  write.table(
    full_res_df,
    file.path(contrast_dir, sprintf("results_full_%s.tsv", label)),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  pdf(file.path(contrast_dir, sprintf("MA_full_%s.pdf", label)), width = 7, height = 6)
  plotMA(
    results(full_fit, contrast = c(factor_name, level_a, level_b)),
    main = paste0(level_a, " vs ", level_b, " (", level_a, " - ", level_b, ")"),
    ylim = c(-5, 5)
  )
  dev.off()

  pdf(file.path(contrast_dir, sprintf("pvalue_hist_full_%s.pdf", label)), width = 7, height = 6)
  hist(full_res_df$pvalue, breaks = 50,
       main = paste0("Full model: ", level_a, " vs ", level_b),
       xlab = "p-value")
  dev.off()

  full_sig <- full_res_df$feature_id[!is.na(full_res_df$padj) & full_res_df$padj < opt$padj_cutoff]

  loo_summary_list <- list()
  delta_long_list <- list()

  for (st in studies) {
    message("  Leave out study: ", st)

    keep <- as.character(colData(dds)[[opt$study_col]]) != st
    dds_loo <- dds[, keep]
    dds_loo <- drop_unused_levels(dds_loo)

    chk <- design_estimable(dds_loo, opt$design)
    if (!chk$ok) {
      loo_summary_list[[st]] <- data.frame(
        contrast = label,
        left_out_study = st,
        status = "skipped",
        reason = chk$reason,
        n_samples = ncol(dds_loo),
        n_genes = nrow(dds_loo),
        pearson_log2fc = NA_real_,
        spearman_stat = NA_real_,
        jaccard_sig = NA_real_,
        prop_direction_concordant = NA_real_,
        median_abs_dlog2fc = NA_real_,
        mean_abs_dlog2fc = NA_real_,
        p90_abs_dlog2fc = NA_real_,
        p95_abs_dlog2fc = NA_real_,
        n_sig_full = length(full_sig),
        n_sig_loo = NA_integer_,
        prop_full_sig_retained = NA_real_,
        top_n_abs_stat_overlap = NA_real_,
        n_genes_compared = NA_integer_,
        n_genes_tested_loo = NA_integer_,
        n_genes_padj_non_na_loo = NA_integer_,
        betaConv_rate_loo = NA_real_,
        median_dispersion = NA_real_,
        median_dispersion_shrinkage = NA_real_,
        median_dispersion_fit_error = NA_real_,
        dispersion_slope = NA_real_,
        prop_dispersion_outliers = NA_real_,
        stringsAsFactors = FALSE
      )
      next
    }

    loo_levels <- unique(as.character(colData(dds_loo)[[factor_name]]))
    if (!all(c(level_a, level_b) %in% loo_levels)) {
      loo_summary_list[[st]] <- data.frame(
        contrast = label,
        left_out_study = st,
        status = "skipped",
        reason = "contrast levels absent after leaving out study",
        n_samples = ncol(dds_loo),
        n_genes = nrow(dds_loo),
        pearson_log2fc = NA_real_,
        spearman_stat = NA_real_,
        jaccard_sig = NA_real_,
        prop_direction_concordant = NA_real_,
        median_abs_dlog2fc = NA_real_,
        mean_abs_dlog2fc = NA_real_,
        p90_abs_dlog2fc = NA_real_,
        p95_abs_dlog2fc = NA_real_,
        n_sig_full = length(full_sig),
        n_sig_loo = NA_integer_,
        prop_full_sig_retained = NA_real_,
        top_n_abs_stat_overlap = NA_real_,
        n_genes_compared = NA_integer_,
        n_genes_tested_loo = NA_integer_,
        n_genes_padj_non_na_loo = NA_integer_,
        betaConv_rate_loo = NA_real_,
        median_dispersion = NA_real_,
        median_dispersion_shrinkage = NA_real_,
        median_dispersion_fit_error = NA_real_,
        dispersion_slope = NA_real_,
        prop_dispersion_outliers = NA_real_,
        stringsAsFactors = FALSE
      )
      next
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

    saveRDS(
      loo_fit,
      file.path(contrast_dir, sprintf("dds_fitted_loo_%s_%s.rds", st, label))
    )

    capture.output({
      cat("left_out_study:\t", st, "\n", sep = "")
      cat("contrast:\t", label, "\n", sep = "")
      cat("minReplicatesForReplace:\t", opt$min_replicates_for_replace, "\n", sep = "")
      cat("n_cores:\t", opt$n_cores, "\n", sep = "")
      cat("n_samples:\t", ncol(loo_fit), "\n", sep = "")
      cat("n_genes:\t", nrow(loo_fit), "\n", sep = "")
      print(resultsNames(loo_fit))
      print(summary(sizeFactors(loo_fit)))
      cat("betaConv_rate:\t", beta_conv_rate(loo_fit), "\n", sep = "")
      cat("median_dispersion:\t", loo_disp_metrics$median_dispersion, "\n", sep = "")
      cat("median_dispersion_shrinkage:\t", loo_disp_metrics$median_dispersion_shrinkage, "\n", sep = "")
      cat("median_dispersion_fit_error:\t", loo_disp_metrics$median_dispersion_fit_error, "\n", sep = "")
      cat("dispersion_slope:\t", loo_disp_metrics$dispersion_slope, "\n", sep = "")
      cat("prop_dispersion_outliers:\t", loo_disp_metrics$prop_dispersion_outliers, "\n", sep = "")
    }, file = file.path(contrast_dir, sprintf("model_summary_loo_%s_%s.txt", st, label)))

    loo_res_df <- extract_res_df(loo_fit, factor_name, level_a, level_b)
    write.table(
      loo_res_df,
      file.path(contrast_dir, sprintf("results_loo_%s_%s.tsv", st, label)),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    pdf(file.path(contrast_dir, sprintf("MA_loo_%s_%s.pdf", st, label)), width = 7, height = 6)
    plotMA(
      results(loo_fit, contrast = c(factor_name, level_a, level_b)),
      main = paste0("LOO ", st, ": ", level_a, " vs ", level_b),
      ylim = c(-5, 5)
    )
    dev.off()

    cmp <- merge(
      full_res_df[, c("feature_id", "log2FoldChange", "stat", "padj")],
      loo_res_df[, c("feature_id", "log2FoldChange", "stat", "padj")],
      by = "feature_id",
      suffixes = c("_full", "_loo"),
      all = FALSE
    )

    cmp <- cmp[!is.na(cmp$log2FoldChange_full) & !is.na(cmp$log2FoldChange_loo), , drop = FALSE]
    cmp$abs_dlog2fc <- abs(cmp$log2FoldChange_loo - cmp$log2FoldChange_full)
    cmp$direction_full <- sign(cmp$log2FoldChange_full)
    cmp$direction_loo <- sign(cmp$log2FoldChange_loo)

    if (opt$lfc_threshold > 0) {
      keep_dir <- abs(cmp$log2FoldChange_full) >= opt$lfc_threshold |
        abs(cmp$log2FoldChange_loo) >= opt$lfc_threshold
      cmp_dir <- cmp[keep_dir, , drop = FALSE]
    } else {
      cmp_dir <- cmp
    }

    prop_concordant <- if (nrow(cmp_dir) > 0) {
      mean(cmp_dir$direction_full == cmp_dir$direction_loo, na.rm = TRUE)
    } else {
      NA_real_
    }

    pearson_lfc <- if (nrow(cmp) > 1) {
      suppressWarnings(cor(
        cmp$log2FoldChange_full,
        cmp$log2FoldChange_loo,
        method = "pearson",
        use = "complete.obs"
      ))
    } else {
      NA_real_
    }

    stat_ok <- !is.na(cmp$stat_full) & !is.na(cmp$stat_loo)
    spearman_stat <- if (sum(stat_ok) > 1) {
      suppressWarnings(cor(
        cmp$stat_full[stat_ok],
        cmp$stat_loo[stat_ok],
        method = "spearman",
        use = "complete.obs"
      ))
    } else {
      NA_real_
    }

    loo_sig <- loo_res_df$feature_id[!is.na(loo_res_df$padj) & loo_res_df$padj < opt$padj_cutoff]
    jaccard_sig <- jaccard_index(full_sig, loo_sig)

    prop_full_sig_retained <- if (length(full_sig) > 0) {
      mean(full_sig %in% loo_sig)
    } else {
      NA_real_
    }

    top_overlap <- top_n_overlap_abs_stat(full_res_df, loo_res_df, n = opt$top_n)

    cmp$contrast <- label
    cmp$left_out_study <- st
    write.table(
      cmp,
      file.path(contrast_dir, sprintf("comparison_full_vs_loo_%s_%s.tsv", st, label)),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    delta_long_list[[st]] <- cmp[, c(
      "feature_id", "contrast", "left_out_study",
      "log2FoldChange_full", "log2FoldChange_loo",
      "abs_dlog2fc", "direction_full", "direction_loo"
    )]

    loo_summary_list[[st]] <- data.frame(
      contrast = label,
      left_out_study = st,
      status = "ok",
      reason = NA_character_,
      n_samples = ncol(dds_loo),
      n_genes = nrow(dds_loo),
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
      median_dispersion = loo_disp_metrics$median_dispersion,
      median_dispersion_shrinkage = loo_disp_metrics$median_dispersion_shrinkage,
      median_dispersion_fit_error = loo_disp_metrics$median_dispersion_fit_error,
      dispersion_slope = loo_disp_metrics$dispersion_slope,
      prop_dispersion_outliers = loo_disp_metrics$prop_dispersion_outliers,
      stringsAsFactors = FALSE
    )
  }

  loo_summary_df <- do.call(rbind, loo_summary_list)
  write.table(
    loo_summary_df,
    file.path(contrast_dir, sprintf("loo_summary_%s.tsv", label)),
    sep = "\t", quote = FALSE, row.names = FALSE
  )

  if (length(delta_long_list) > 0) {
    delta_long_df <- do.call(rbind, delta_long_list)
    write.table(
      delta_long_df,
      file.path(contrast_dir, sprintf("abs_dlog2fc_long_%s.tsv", label)),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    ok_df <- subset(loo_summary_df, status == "ok")

    p1 <- ggplot(ok_df, aes(x = left_out_study, y = pearson_log2fc)) +
      geom_col() +
      coord_flip() +
      theme_bw() +
      ylab("Pearson correlation of log2FC") +
      xlab("Left-out study")
    ggsave(file.path(contrast_dir, sprintf("pearson_log2fc_%s.pdf", label)),
           plot = p1, width = 7, height = 5)

    p2 <- ggplot(ok_df, aes(x = left_out_study, y = jaccard_sig)) +
      geom_col() +
      coord_flip() +
      theme_bw() +
      ylab("Jaccard index of significant genes") +
      xlab("Left-out study")
    ggsave(file.path(contrast_dir, sprintf("jaccard_sig_%s.pdf", label)),
           plot = p2, width = 7, height = 5)

    p3 <- ggplot(subset(delta_long_df, !is.na(abs_dlog2fc)),
                 aes(x = left_out_study, y = abs_dlog2fc)) +
      geom_boxplot(outlier.size = 0.2) +
      coord_flip() +
      theme_bw() +
      ylab("|Δlog2FC|") +
      xlab("Left-out study")
    ggsave(file.path(contrast_dir, sprintf("abs_dlog2fc_boxplot_%s.pdf", label)),
           plot = p3, width = 8, height = 6)

    p4 <- ggplot(ok_df, aes(x = left_out_study, y = prop_direction_concordant)) +
      geom_col() +
      coord_flip() +
      theme_bw() +
      ylab("Proportion of concordant effect direction") +
      xlab("Left-out study")
    ggsave(file.path(contrast_dir, sprintf("direction_concordance_%s.pdf", label)),
           plot = p4, width = 7, height = 5)

    p5 <- ggplot(ok_df, aes(x = left_out_study, y = spearman_stat)) +
      geom_col() +
      coord_flip() +
      theme_bw() +
      ylab("Spearman correlation of stat") +
      xlab("Left-out study")
    ggsave(file.path(contrast_dir, sprintf("spearman_stat_%s.pdf", label)),
           plot = p5, width = 7, height = 5)

    p6 <- ggplot(ok_df, aes(x = left_out_study, y = top_n_abs_stat_overlap)) +
      geom_col() +
      coord_flip() +
      theme_bw() +
      ylab(sprintf("Top %d overlap by |stat|", opt$top_n)) +
      xlab("Left-out study")
    ggsave(file.path(contrast_dir, sprintf("topn_abs_stat_overlap_%s.pdf", label)),
           plot = p6, width = 7, height = 5)

    p7 <- ggplot(ok_df, aes(x = left_out_study, y = median_dispersion_fit_error)) +
      geom_col() +
      coord_flip() +
      theme_bw() +
      ylab("Median dispersion fit error") +
      xlab("Left-out study")
    ggsave(file.path(contrast_dir, sprintf("dispersion_fit_error_%s.pdf", label)),
           plot = p7, width = 7, height = 5)

    p8 <- ggplot(ok_df, aes(x = left_out_study, y = prop_dispersion_outliers)) +
      geom_col() +
      coord_flip() +
      theme_bw() +
      ylab("Proportion of dispersion outliers") +
      xlab("Left-out study")
    ggsave(file.path(contrast_dir, sprintf("dispersion_outliers_%s.pdf", label)),
           plot = p8, width = 7, height = 5)

    p9 <- ggplot(ok_df, aes(x = left_out_study, y = betaConv_rate_loo)) +
      geom_col() +
      coord_flip() +
      theme_bw() +
      ylab("betaConv rate") +
      xlab("Left-out study")
    ggsave(file.path(contrast_dir, sprintf("betaConv_rate_%s.pdf", label)),
           plot = p9, width = 7, height = 5)
  }
}