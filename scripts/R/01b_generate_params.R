#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(data.table)
  library(DESeq2)
  library(SummarizedExperiment)
})

x <- fread("contrasts.tsv", sep = "\t", header = TRUE, data.table = FALSE)

out <- list()

for (i in seq_len(nrow(x))) {
  dds <- readRDS(x$dds_rds[i])
  studies <- unique(as.character(colData(dds)[["Study"]]))
  studies <- studies[!is.na(studies) & studies != ""]

  out[[i]] <- data.frame(
    phenotype = x$phenotype[i],
    dds_rds = x$dds_rds[i],
    contrast = x$contrast[i],
    leave_out_study = studies,
    stringsAsFactors = FALSE
  )
}

y <- do.call(rbind, out)
write.table(y, file = "loo_params.tsv", sep = "\t", quote = FALSE, row.names = FALSE)