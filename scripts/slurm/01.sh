#!/bin/bash
#SBATCH --job-name=01_build
#SBATCH --mail-user=ew255@st-andrews.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --export=ALL
#SBATCH --cpus-per-task=4
#SBATCH --mem=32G
#SBATCH --partition=large-long
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --chdir=/sharedscratch/ew255/md4002-erik-muscle-atrophy

source /software/conda/ew255/conda/etc/profile.d/conda.sh
conda activate /software/conda/ew255/conda/envs/md4002_erik
source .env

INPDIR="${DATA_DIR}/datasets/hsa-rnaseq-muscle"
OUTDIR="/sharedscratch/ew255/md4002-erik-muscle-atrophy/results"
SCRIPTDIR="/sharedscratch/ew255/md4002-erik-muscle-atrophy/scripts/R"

Rscript "${SCRIPTDIR}/01_build_split_validate_dds.R" \
  --se-rds "${INPDIR}/null.merged.gene.SummarizedExperiment.rds" \
  --metadata-xlsx "${INPDIR}/datasets-muscle.xlsx" \
  --metadata-clean-csv "${INPDIR}/pheno-clean.csv" \
  --assay-name "salmon.merged.gene_counts" \
  --design "~ Study + Condition" \
  --split-col "PHENOTYPE" \
  --outdir "${OUTDIR}/dds_by_phenotype"
