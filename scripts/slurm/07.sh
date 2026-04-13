#!/bin/bash
#SBATCH --job-name=07_agg_loo
#SBATCH --mail-user=ew255@st-andrews.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --export=ALL
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --partition=large-long,large-short
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --chdir=/sharedscratch/ew255/md4002-erik-muscle-atrophy

source /software/conda/ew255/conda/etc/profile.d/conda.sh
conda activate /software/conda/ew255/conda/envs/md4002_erik

INPUT_ROOT="results/deseq_loo_array"
OUTDIR="results/loo_aggregated"
SCRIPTDIR="/sharedscratch/ew255/md4002-erik-muscle-atrophy/scripts/R"

Rscript "${SCRIPTDIR}/06_aggregate_loo.R" \
  --input-root "${INPUT_ROOT}" \
  --outdir "${OUTDIR}" \
  --width 14 \
  --height-per-facet 2.8 \
  --top-n-influential 5