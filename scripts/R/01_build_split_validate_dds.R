#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(optparse)
  library(SummarizedExperiment)
  library(DESeq2)
  library(openxlsx2)
})

option_list <- list(
  make_option("--se-rds", dest = "se_rds", type = "character",
              help = "Input SummarizedExperiment RDS"),
  make_option("--metadata-xlsx", dest = "metadata_xlsx", type = "character",
              help = "Metadata Excel file"),
  make_option("--metadata-clean-csv", dest = "metadata_clean_csv", type = "character",
              help = "Clean phenotype CSV with Run and Condition"),
  make_option("--assay-name", dest = "assay_name", type = "character",
              default = "salmon.merged.gene_counts",
              help = "Assay name in SummarizedExperiment [default %default]"),
  make_option("--design", dest = "design", type = "character",
              default = "~ Study + Condition",
              help = "Design formula [default %default]"),
  make_option("--split-col", dest = "split_col", type = "character",
              default = "PHENOTYPE",
              help = "Column used to split metadata [default %default]"),
  make_option("--outdir", dest = "outdir", type = "character",
              help = "Output directory")
)

opt <- parse_args(OptionParser(option_list = option_list))

req <- c("se_rds", "metadata_xlsx", "metadata_clean_csv", "outdir")
miss <- req[!nzchar(vapply(req, function(x) ifelse(is.null(opt[[x]]), "", opt[[x]]), character(1)))]
if (length(miss) > 0) {
  stop("Missing required arguments: ", paste(paste0("--", gsub("_", "-", miss)), collapse = ", "))
}

dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

f_read_metadata <- function(path_metadata) {
  wb_1 <- read_xlsx(file = path_metadata, sheet = 1)
  wb_4 <- read_xlsx(file = path_metadata, sheet = 4)

  wb_1$ID_EXTERNAL <- trimws(wb_1$ID_EXTERNAL)
  wb_4$ID_EXTERNAL <- trimws(wb_4$ID_EXTERNAL)

  wb_1 <- wb_1[!duplicated(wb_1), , drop = FALSE]
  wb_4 <- wb_4[!duplicated(wb_4), , drop = FALSE]

  df <- merge(wb_1, wb_4, by = "ID_EXTERNAL")
  df <- df[, !colnames(df) %in% c("N_TOTAL", "NOTES"), drop = FALSE]

  if ("Run" %in% colnames(df)) {
    df$Run <- trimws(df$Run)
  }

  df
}

f_validate_metadata <- function(metadata, split_col) {
  needed <- c("Run", "Condition", "Study", split_col)
  missing_cols <- setdiff(needed, colnames(metadata))
  if (length(missing_cols) > 0) {
    stop("Metadata is missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  if (anyNA(metadata$Run) || any(metadata$Run == "")) {
    stop("Metadata contains missing/empty Run values")
  }

  if (anyDuplicated(metadata$Run) > 0) {
    dup <- unique(metadata$Run[duplicated(metadata$Run)])
    stop("Metadata contains duplicated Run values, e.g.: ", paste(head(dup, 10), collapse = ", "))
  }

  if (anyNA(metadata[[split_col]]) || any(metadata[[split_col]] == "")) {
    stop("Metadata contains missing/empty values in split column: ", split_col)
  }

  invisible(TRUE)
}

f_make_dds <- function(se_file, metadata, assay_name, design_formula) {
  se_obj <- readRDS(se_file)

  if (!inherits(se_obj, "SummarizedExperiment")) {
    stop("Input RDS is not a SummarizedExperiment")
  }

  if (!assay_name %in% assayNames(se_obj)) {
    stop("Assay not found: ", assay_name)
  }

  cts <- assay(se_obj, assay_name)
  all_samples <- colnames(cts)

  metadata$Run <- trimws(metadata$Run)
  metadata <- metadata[!is.na(metadata$Run) & metadata$Run != "", , drop = FALSE]

  keep_samples <- intersect(metadata$Run, all_samples)
  if (length(keep_samples) == 0) {
    stop("No metadata rows matched count matrix samples")
  }

  md <- metadata[match(keep_samples, metadata$Run), , drop = FALSE]
  cts_sub <- cts[, keep_samples, drop = FALSE]

  stopifnot(identical(colnames(cts_sub), md$Run))

  mat <- as.matrix(cts_sub)
  mode(mat) <- "numeric"
  mat <- round(mat)
  storage.mode(mat) <- "integer"

  rownames(md) <- md$Run
  md$Study <- factor(md$Study)
  md$Condition <- factor(md$Condition)

  if (nlevels(md$Condition) < 2) {
    stop("Less than two Condition levels present after subsetting")
  }

  dds <- DESeqDataSetFromMatrix(
    countData = mat,
    colData = md,
    design = as.formula(design_formula),
    ignoreRank = TRUE
  )

  annotations <- rowData(se_obj)
  ann_df <- as.data.frame(annotations)

  if ("gene_id" %in% colnames(ann_df)) {
    ann_idx <- match(rownames(dds), ann_df$gene_id)
    rowData(dds) <- cbind(rowData(dds), annotations[ann_idx, , drop = FALSE])
  }

  dds <- estimateSizeFactors(dds)
  dds
}

f_write_validation <- function(dds, phenotype_name, outdir) {
  pheno_dir <- file.path(outdir, phenotype_name)
  dir.create(pheno_dir, recursive = TRUE, showWarnings = FALSE)

  saveRDS(dds, file = file.path(pheno_dir, paste0("dds_", phenotype_name, ".rds")))

  coldata_df <- as.data.frame(colData(dds))
  write.csv(
    coldata_df,
    file = file.path(pheno_dir, paste0("coldata_", phenotype_name, ".csv")),
    row.names = TRUE,
    quote = TRUE
  )

  counts_df <- data.frame(
    metric = c(
      "n_genes",
      "n_samples",
      "size_factors_estimated",
      "all_size_factors_positive"
    ),
    value = c(
      nrow(dds),
      ncol(dds),
      !all(is.na(sizeFactors(dds))),
      all(sizeFactors(dds) > 0, na.rm = TRUE)
    ),
    stringsAsFactors = FALSE
  )
  write.csv(
    counts_df,
    file = file.path(pheno_dir, paste0("validation_summary_", phenotype_name, ".csv")),
    row.names = FALSE,
    quote = TRUE
  )

  tab_condition <- as.data.frame(table(colData(dds)$Condition), stringsAsFactors = FALSE)
  colnames(tab_condition) <- c("Condition", "N")
  write.csv(
    tab_condition,
    file = file.path(pheno_dir, paste0("condition_counts_", phenotype_name, ".csv")),
    row.names = FALSE,
    quote = TRUE
  )

  tab_study <- as.data.frame(table(colData(dds)$Study), stringsAsFactors = FALSE)
  colnames(tab_study) <- c("Study", "N")
  write.csv(
    tab_study,
    file = file.path(pheno_dir, paste0("study_counts_", phenotype_name, ".csv")),
    row.names = FALSE,
    quote = TRUE
  )

  tab_joint <- as.data.frame(table(colData(dds)$Study, colData(dds)$Condition), stringsAsFactors = FALSE)
  colnames(tab_joint) <- c("Study", "Condition", "N")
  write.csv(
    tab_joint,
    file = file.path(pheno_dir, paste0("study_condition_counts_", phenotype_name, ".csv")),
    row.names = FALSE,
    quote = TRUE
  )
}

df_metadata <- f_read_metadata(opt$metadata_xlsx)

pheno_file <- read.csv(opt$metadata_clean_csv, stringsAsFactors = FALSE, check.names = FALSE)
if (ncol(pheno_file) < 2) {
  stop("Clean phenotype CSV must contain at least two columns: Run and Condition")
}
colnames(pheno_file)[1:2] <- c("Run", "Condition")

df_metadata$Run <- trimws(df_metadata$Run)
pheno_file$Run <- trimws(pheno_file$Run)
pheno_file$Condition <- trimws(pheno_file$Condition)

pheno_file <- pheno_file[!is.na(pheno_file$Run) & pheno_file$Run != "", , drop = FALSE]
pheno_file <- pheno_file[!duplicated(pheno_file$Run), , drop = FALSE]

df_metadata <- merge(df_metadata, pheno_file[, c("Run", "Condition")], by = "Run")
df_metadata$Study <- df_metadata$ID_EXTERNAL

dup_runs <- unique(df_metadata$Run[duplicated(df_metadata$Run)])

if (length(dup_runs) > 0) {
  message("Keeping first occurrence for duplicated Run entries: ", length(dup_runs))
  message("Examples: ", paste(head(dup_runs, 10), collapse = ", "))
}

df_metadata <- df_metadata[!duplicated(df_metadata$Run), , drop = FALSE]

f_validate_metadata(df_metadata, opt$split_col)

phenotypes <- unique(df_metadata[[opt$split_col]])
phenotypes <- phenotypes[!is.na(phenotypes) & phenotypes != ""]

write.csv(
  data.frame(PHENOTYPE = phenotypes, stringsAsFactors = FALSE),
  file = file.path(opt$outdir, "phenotypes_detected.csv"),
  row.names = FALSE,
  quote = TRUE
)

for (ph in phenotypes) {
  message("Processing phenotype: ", ph)

  md_sub <- df_metadata[df_metadata[[opt$split_col]] == ph, , drop = FALSE]

  dds <- f_make_dds(
    se_file = opt$se_rds,
    metadata = md_sub,
    assay_name = opt$assay_name,
    design_formula = opt$design
  )

  f_write_validation(
    dds = dds,
    phenotype_name = ph,
    outdir = opt$outdir
  )
}