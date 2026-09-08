#!/bin/bash
#SBATCH --job-name=16S_pipeline
#SBATCH --mail-user=your_email@example.com
#SBATCH --mail-type=end,fail
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=10
#SBATCH --mem-per-cpu=10G
#SBATCH --time=168:00:00
#SBATCH --qos=standard
#SBATCH --output=logs/pipeline_%j.out
#SBATCH --error=logs/pipeline_%j.err

# Create logs directory
mkdir -p logs

# Run the pipeline
./qiime2_picrust2_pipeline.sh \
    --work-dir /project/microbiome \
    --seqs-dir /project/data/reads \
    --forward-pattern _R1 \
    --reverse-pattern _R2 \
    --n-char 5 \
    --o-manifest results/manifest \
    --o-qiime results/qiime2 \
    --forward-primer CCTACGGGNGGCWGCAG \
    --reverse-primer GGACTACNVGGGTWTCTAAT \
    --d-denoise dada2 \
    --p-picrust true \
    --threads ${SLURM_CPUS_PER_TASK}
