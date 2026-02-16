#!/bin/bash
set -e  # Exit on any error

# Setup
echo "=========================================="
echo "Pipeline Submission Started"
echo "Start time: $(date +%Y%m%d_%H%M%S)"
echo "=========================================="
# Load Snakemake
module load snakemake
module load slurm
snakemake --version

# Create output directories
echo "Creating output directories..."
mkdir -p logs/slurm
mkdir -p results/{qc,vesalius,kuresi,plots}

snakemake \
    --profile config/snakemake_profile \
    --latency-wait 600 \
    --nolock \
    --verbose \
    --printshellcmds \
    2>&1 | tee logs/snakemake_$(date +%Y%m%d_%H%M%S).log