#!/bin/bash
#SBATCH --job-name=04_gsea
#SBATCH --mail-user=ew255@st-andrews.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --export=ALL
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --partition=large-long,large-short
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --chdir=/sharedscratch/ew255/md4002-erik-muscle-atrophy

source /software/conda/ew255/conda/etc/profile.d/conda.sh
conda activate /software/conda/ew255/conda/envs/r-gsea
source .env

INPDIR="${DATA_DIR}/datasets/hsa-rnaseq-muscle"
OUTDIR="/sharedscratch/ew255/md4002-erik-muscle-atrophy/results"
SCRIPTDIR="/sharedscratch/ew255/md4002-erik-muscle-atrophy/scripts/R"

OUTDIR="./results/gsea"

Rscript "${SCRIPTDIR}/03_run_gsea.R" \
  --results-tsv "./results/deseq_results/Ageing/results_Condition_ageing_old_vs_ageing_young.tsv" \
  --id-type "ENSEMBL" \
  --id-col "feature_id" \
  --rank-col "stat" \
  --outdir "$OUTDIR"
