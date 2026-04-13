#!/bin/bash
#SBATCH --job-name=06_model_loo
#SBATCH --mail-user=ew255@st-andrews.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --export=ALL
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --array=1-22
#SBATCH --partition=large-long,large-short
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --chdir=/sharedscratch/ew255/md4002-erik-muscle-atrophy

source /software/conda/ew255/conda/etc/profile.d/conda.sh
conda activate /software/conda/ew255/conda/envs/md4002_erik

SCRIPTDIR="/sharedscratch/ew255/md4002-erik-muscle-atrophy/scripts/R"
PARAMS="loo_params.tsv"

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" "$PARAMS")

PHENO=$(echo "$LINE" | cut -f1)
DDS=$(echo "$LINE" | cut -f2)
CONTRAST=$(echo "$LINE" | cut -f3)
LEAVE_OUT_STUDY=$(echo "$LINE" | cut -f4)

echo "Processing phenotype: $PHENO"
echo "DDS: $DDS"
echo "Contrast: $CONTRAST"
echo "Leave-out study: $LEAVE_OUT_STUDY"

OUTDIR="results/deseq_loo_array/${PHENO}"

Rscript "${SCRIPTDIR}/05_run_loo_study_array.R" \
  --dds-rds "$DDS" \
  --design "~ Study + Condition" \
  --study-col "Study" \
  --leave-out-study "$LEAVE_OUT_STUDY" \
  --contrast "$CONTRAST" \
  --outdir "$OUTDIR" \
  --padj-cutoff 0.05 \
  --top-n 100 \
  --min-replicates-for-replace Inf \
  --n-cores 1
