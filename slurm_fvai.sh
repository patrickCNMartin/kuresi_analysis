#!/bin/bash
#SBATCH --job-name=fvai
#SBATCH --output=logs/%x-%j.out
#SBATCH --error=logs/%x-%j.err
#SBATCH --time=24:00:00
#SBATCH --cpus-per-task=1
#SBATCH --mem=1G
#SBATCH --partition=defq
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --mail-type=FAIL,BEGIN,END
#SBATCH --mail-user=<patrick.martin@cshs.org>

module load snakemake
module load apptainer

# making directories
echo "Creating output directories"
mkdir -p report logs

# Run pipeline
echo "Starting pipeline at $(date +%Y%m%d_%H%M%S)"

snakemake \
    --profile slurm \
    --use-singularity \
    --singularity-args "--bind ./data:/data" \
    --cores 10 \
    --max-jobs-per-second 1\
    --max-status-checks-per-second 0.1 \
    2>&1 | tee logs/snakemake_$(date +%Y%m%d_%H%M%S).log

echo "Pipeline completed at $(date +%Y%m%d_%H%M%S)"