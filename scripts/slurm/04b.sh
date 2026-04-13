#!/bin/bash
#SBATCH --job-name=04b_gsea_array
#SBATCH --mail-user=ew255@st-andrews.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --export=ALL
#SBATCH --cpus-per-task=12
#SBATCH --mem=32G
#SBATCH --array=1-3
#SBATCH --partition=large-long
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --chdir=/sharedscratch/ew255/md4002-erik-muscle-atrophy

source /software/conda/ew255/conda/etc/profile.d/conda.sh
conda activate /software/conda/ew255/conda/envs/r-gsea
source .env

INPDIR="${DATA_DIR}/datasets/hsa-rnaseq-muscle"
SCRIPTDIR="/sharedscratch/ew255/md4002-erik-muscle-atrophy/scripts/R"

cat > gsea_input.tsv << 'EOF'
phenotype	file
Atrophy	./results/deseq_results/Atrophy/results_Condition_atrophy_post_vs_atrophy_pre.tsv
Exercise	./results/deseq_results/Exercise/results_Condition_exercise_post_vs_exercise_pre.tsv
Ageing	./results/deseq_results/Ageing/results_Condition_ageing_old_vs_ageing_young.tsv
EOF

PARAMS=gsea_input.tsv

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $PARAMS)

PHENO=$(echo "$LINE" | cut -f1)
FILE=$(echo "$LINE" | cut -f2)

echo "Processing phenotype: $PHENO"
echo "File: $FILE"

OUTDIR=results/gsea_array/${PHENO}

Rscript "${SCRIPTDIR}/03_run_gsea.R" \
  --results-tsv "$FILE" \
  --id-type "ENSEMBL" \
  --id-col "feature_id" \
  --rank-col "stat" \
  --outdir "$OUTDIR"
