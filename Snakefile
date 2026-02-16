# Snakemake workflow for multi-sample spatial analysis - HPC/SLURM Edition
# Pipeline: QC -> Vesalius -> Kuresi -> Plotting

import os
from pathlib import Path

# Load configuration
configfile: "config/config.yaml"

# Define sample names
SAMPLES = config["samples"]

# Define output directories
QC_DIR = "qc"
VESALIUS_DIR = "vesalius"
KURESI_DIR = "kuresi"
PLOTS_DIR = "plots"
RESULTS_BASE = "results"

# Create a rule to define all final outputs
rule all:
    input:
        expand(f"{RESULTS_BASE}/{{sample}}/{QC_DIR}/qc_data.rds", sample=SAMPLES),
        expand(f"{RESULTS_BASE}/{{sample}}/{VESALIUS_DIR}/vesalius_data.rds", sample=SAMPLES),
        expand(f"{RESULTS_BASE}/{{sample}}/{KURESI_DIR}/kuresi_competition_scores.rds", sample=SAMPLES),
        expand(f"{RESULTS_BASE}/{{sample}}/report.txt", sample=SAMPLES),



# ============================================================================
# RULE 1: QC - Quality Control
# ============================================================================
rule qc:
    input:
        coordinates = lambda wildcards: config["input_files"][wildcards.sample]["coordinates"],
        counts = lambda wildcards: config["input_files"][wildcards.sample]["counts"],
        scale = lambda wildcards: config["input_files"][wildcards.sample]["scale"],
        image = lambda wildcards: config["input_files"][wildcards.sample]["image"],
    output:
        qc_rds = f"{RESULTS_BASE}/{{sample}}/{QC_DIR}/qc_data.rds",
        qc_1 = f"{RESULTS_BASE}/{{sample}}/{QC_DIR}/feature_map.pdf",
        qc_2 = f"{RESULTS_BASE}/{{sample}}/{QC_DIR}/count_map.pdf",
        qc_3 = f"{RESULTS_BASE}/{{sample}}/{QC_DIR}/feature_map_win_loose.pdf",
        qc_4 = f"{RESULTS_BASE}/{{sample}}/{QC_DIR}/count_map_win_loose.pdf",
        qc_report = f"{RESULTS_BASE}/{{sample}}/{QC_DIR}/qc_report.txt",
    params:
        sample_name = "{sample}",
        output_dir = f"{RESULTS_BASE}/{{sample}}/{QC_DIR}",
        bin_size = config["qc"]["bin_size"],
        min_cells = config["qc"]["min_cells"],
        min_features = config["qc"]["min_features"]
    threads: config["resources"]["qc"]["cpus"]
    resources:
        mem_mb = config["resources"]["qc"]["mem_mb"],
        time_min = config["resources"]["qc"]["time_min"],
        slurm_partition = config["slurm"]["partition"],
        slurm_extra = f"--mail-type={config['slurm']['mail_type']} --mail-user={config['slurm']['mail_user']}",
    log:
        "logs/qc/{sample}.log"
    shell:
        """
        module load R/4.4.0 rlibs/4.4.0 hdf5
        mkdir -p {params.output_dir}
        
        Rscript scripts/qc.r \
            --sample_name {params.sample_name} \
            --coordinates {input.coordinates} \
            --counts {input.counts} \
            --scale {input.scale} \
            --min_cells {params.min_cells} \
            --min_features {params.min_features} \
            --bin_size {params.bin_size} \
            --output_dir {params.output_dir} \
            --report_file {output.qc_report} \
            2>&1 | tee {log}
        """


# ============================================================================
# RULE 2: Vesalius - Spatial analysis and dimensionality reduction
# ============================================================================
rule vesalius:
    input:
        qc_data = rules.qc.output.qc_rds,
    output:
        vesalius_rds = f"{RESULTS_BASE}/{{sample}}/{VESALIUS_DIR}/vesalius_data.rds",
        vesalius_territory = f"{RESULTS_BASE}/{{sample}}/{VESALIUS_DIR}/vesalius_territory_plot.pdf",
        vesalius_image = f"{RESULTS_BASE}/{{sample}}/{VESALIUS_DIR}/vesalius_image_plot.pdf",
        vesalius_report = f"{RESULTS_BASE}/{{sample}}/{VESALIUS_DIR}/vesalius_report.txt",
    params:
        sample_name = "{sample}",
        output_dir = f"{RESULTS_BASE}/{{sample}}/{VESALIUS_DIR}",
        feature_set = config["vesalius"]["feature_set"],
        resolution = config["vesalius"]["resolution"],
        dim_reduc = config["vesalius"]["dim_reduc"],
        dims = config["vesalius"]["dims"],
        equalize = config["vesalius"]["equalize"],
        sigma = config["vesalius"]["sigma"],
        iter = config["vesalius"]["iter"],
        col_resolution = config["vesalius"]["col_resolution"],
        distance = config["vesalius"]["distance"],
    threads: config["resources"]["vesalius"]["cpus"]
    resources:
        mem_mb = config["resources"]["vesalius"]["mem_mb"],
        time_min = config["resources"]["vesalius"]["time_min"],
        slurm_partition = config["slurm"]["partition"],
        slurm_extra = f"--mail-type={config['slurm']['mail_type']} --mail-user={config['slurm']['mail_user']}",
    log:
        "logs/vesalius/{sample}.log"
    shell:
        """
        module load R/4.4.0 rlibs/4.4.0 hdf5
        mkdir -p {params.output_dir}
        
        Rscript scripts/vesalius.r \
            --sample_name {params.sample_name} \
            --input_rds {input.qc_data} \
            --feature_set {params.feature_set} \
            --resolution {params.resolution} \
            --dim_reduc {params.dim_reduc} \
            --dims {params.dims} \
            --equalize {params.equalize} \
            --sigma {params.sigma} \
            --iter {params.iter} \
            --col_resolution {params.col_resolution} \
            --distance {params.distance} \
            --output_dir {params.output_dir} \
            --report_file {output.vesalius_report} \
            --cores {threads} \
            2>&1 | tee {log}
        """


# ============================================================================
# RULE 3: Kuresi - Statistical analysis
# ============================================================================
rule kuresi:
    input:
        vesalius_data = rules.vesalius.output.vesalius_rds,
    output:
        kuresi_rds = f"{RESULTS_BASE}/{{sample}}/{KURESI_DIR}/kuresi_competition_scores.rds",
        kuresi_plot = f"{RESULTS_BASE}/{{sample}}/{KURESI_DIR}/kuresi_score_plot.pdf",
        kuresi_report = f"{RESULTS_BASE}/{{sample}}/{KURESI_DIR}/kuresi_report.txt",
    params:
        sample_name = "{sample}",
        output_dir = f"{RESULTS_BASE}/{{sample}}/{KURESI_DIR}",
        method = config["kuresi"]["method"],
        scale = config["kuresi"]["scale"],
        center = config["kuresi"]["center"],
    threads: config["resources"]["kuresi"]["cpus"]
    resources:
        mem_mb = config["resources"]["kuresi"]["mem_mb"],
        time_min = config["resources"]["kuresi"]["time_min"],
        slurm_partition = config["slurm"]["partition"],
        slurm_extra = f"--mail-type={config['slurm']['mail_type']} --mail-user={config['slurm']['mail_user']}",
    log:
        "logs/kuresi/{sample}.log"
    shell:
        """
        module load R/4.4.0 rlibs/4.4.0 hdf5
        mkdir -p {params.output_dir}
        
        Rscript scripts/kuresi.r \
            --sample_name {params.sample_name} \
            --input_rds {input.vesalius_data} \
            --method {params.method} \
            --scale {params.scale} \
            --center {params.center} \
            --output_dir {params.output_dir} \
            --report_file {output.kuresi_report} \
            2>&1 | tee {log}
        """


# ============================================================================
# RULE 4: Concatenate all reports into unified report
# ============================================================================
rule concatenate_reports:
    input:
        qc_report = rules.qc.output.qc_report,
        vesalius_report = rules.vesalius.output.vesalius_report,
        kuresi_report = rules.kuresi.output.kuresi_report,
    output:
        unified_report = f"{RESULTS_BASE}/{{sample}}/report.txt",
    shell:
        """
        cat {input.qc_report} {input.vesalius_report} {input.kuresi_report} > {output.unified_report}
        """
