#!/bin/bash
#SBATCH --job-name=02_model
#SBATCH --mail-user=ew255@st-andrews.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --export=ALL
#SBATCH --cpus-per-task=12
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

Rscript "${SCRIPTDIR}/02_run_deseq2_model.R" \
  --dds-rds "./results/dds_by_phenotype/Atrophy/dds_Atrophy.rds" \
  --design "~ Study + Condition" \
  --contrast "Condition,atrophy_post,atrophy_pre" \
  --outdir "${OUTDIR}/models"
