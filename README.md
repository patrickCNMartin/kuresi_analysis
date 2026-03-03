# kuresi_analysis

Cancer clone competition analysis pipeline using [Kuresi](https://github.com/patrickCNMartin/Kuresi).

## Overview

This repository contains a Snakemake pipeline for running end-to-end competitive fitness analysis on cancer single cell and spatial transcriptomics data. 

## Repository Structure

```
kuresi_analysis/
├── config/          # Pipeline configuration files
├── data/            # Placeholder for processed data (see below)
├── envs/            # Required R packages
└── Snakefile        # Pipeline entry point
```

## Data

Processed data is not distributed with this repository. Placeholders are provided in `data/` alongside md5 checksums for verification. 


## Usage

### Prerequisites

- [Snakemake](https://snakemake.readthedocs.io/) (v7+)
- R packages shown in `envs/install_packages.r`

Note: The packages will be installed in a `lib_cache` directory create in the root directory. 


### Running the pipeline

Edit `config/config.yaml` to point to your data and set any relevant parameters, then:

```bash
snakemake \
    --profile config/snakemake_profile \
    --latency-wait 600 \
    --nolock \
    --verbose \
    --printshellcmds
```

Note: The analysis was run on a HPC with SLURM. The `config/snakemake_profile` handles HPC specific configuration. 