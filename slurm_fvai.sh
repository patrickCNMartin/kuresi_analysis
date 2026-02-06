#!/bin/bash
#SBATCH --job-name=snakemake_spatial
#SBATCH --partition=defq
#SBATCH --time=72:00:00
#SBATCH --mem=4G
#SBATCH --cpus-per-task=1
#SBATCH --output=logs/snakemake_submission_%j.out
#SBATCH --error=logs/snakemake_submission_%j.err
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=<patrick.martin@cshs.org>

# ============================================================================
# Snakemake Pipeline Submission Script
# Submit with: sbatch submit_pipeline.sh
# Check status: squeue -u $USER
# ============================================================================

set -e  # Exit on any error

# Setup
echo "=========================================="
echo "Pipeline Submission Started"
echo "Start time: $(date +%Y%m%d_%H%M%S)"
echo "=========================================="

# Load Snakemake
echo "Loading Snakemake module..."
module load snakemake
module load slurm
snakemake --version

# Create output directories
echo "Creating output directories..."
mkdir -p logs/slurm
mkdir -p results/{qc,vesalius,kuresi,plots}

# ============================================================================
# Run Snakemake Pipeline
# ============================================================================

echo ""
echo "Starting Snakemake pipeline at $(date +%Y%m%d_%H%M%S)"
echo "=========================================="

snakemake \
    --profile config/snakemake_profile \
    --verbose \
    2>&1 | tee logs/snakemake_$(date +%Y%m%d_%H%M%S).log

# Capture exit code
SNAKEMAKE_EXIT=$?

# ============================================================================
# Post-pipeline reporting
# ============================================================================

echo ""
echo "=========================================="
echo "Pipeline Completed"
echo "End time: $(date +%Y%m%d_%H%M%S)"
echo "Exit code: $SNAKEMAKE_EXIT"
echo "=========================================="

# Check if all outputs were created
echo ""
echo "Checking output files..."
echo "QC steps completed: $(find results/qc -name '*.rds' | wc -l) / 6"
echo "Vesalius steps completed: $(find results/vesalius -name '*.rds' | wc -l) / 6"
echo "Kuresi steps completed: $(find results/kuresi -name '*.rds' | wc -l) / 6"
echo "Final plots completed: $(find results/plots -name 'final_plot.pdf' | wc -l) / 6"

# Exit with Snakemake's exit code
exit $SNAKEMAKE_EXIT