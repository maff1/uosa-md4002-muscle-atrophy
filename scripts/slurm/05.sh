#!/bin/bash
#SBATCH --job-name=05_model_loo
#SBATCH --mail-user=ew255@st-andrews.ac.uk
#SBATCH --mail-type=END,FAIL
#SBATCH --export=ALL
#SBATCH --cpus-per-task=12
#SBATCH --mem=48G
#SBATCH --array=1-3
#SBATCH --partition=large-long,bigmem,vbigmem
#SBATCH --output=logs/%x_%j.out
#SBATCH --error=logs/%x_%j.err
#SBATCH --chdir=/sharedscratch/ew255/md4002-erik-muscle-atrophy

source /software/conda/ew255/conda/etc/profile.d/conda.sh
conda activate /software/conda/ew255/conda/envs/md4002_erik
source .env

INPDIR="${DATA_DIR}/datasets/hsa-rnaseq-muscle"
SCRIPTDIR="/sharedscratch/ew255/md4002-erik-muscle-atrophy/scripts/R"

cat > contrasts.tsv << 'EOF'
phenotype	dds_rds	contrast
Atrophy	./results/dds_by_phenotype/Atrophy/dds_Atrophy.rds	Condition,atrophy_post,atrophy_pre
Exercise	./results/dds_by_phenotype/Exercise/dds_Exercise.rds	Condition,exercise_post,exercise_pre
Ageing	./results/dds_by_phenotype/Ageing/dds_Ageing.rds	Condition,ageing_old,ageing_young
EOF

PARAMS=contrasts.tsv

LINE=$(sed -n "$((SLURM_ARRAY_TASK_ID+1))p" $PARAMS)

PHENO=$(echo "$LINE" | cut -f1)
DDS=$(echo "$LINE" | cut -f2)
CONTRAST=$(echo "$LINE" | cut -f3)

echo "Processing phenotype: $PHENO"
echo "DDS: $DDS"
echo "Contrast: $CONTRAST"

OUTDIR=results/deseq_loo/${PHENO}

Rscript "${SCRIPTDIR}/04_run_loo_study.R" \
  --dds-rds "$DDS" \
  --design "~ Study + Condition" \
  --study-col "Study" \
  --contrast "$CONTRAST" \
  --outdir "$OUTDIR" \
  --padj-cutoff 0.05 \
  --top-n 100 \
  --min-replicates-for-replace Inf \
  --n-cores "${SLURM_CPUS_PER_TASK}"
