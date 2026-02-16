#!/usr/bin/env Rscript

# This will manually install all the packages in the lib_cache directory
# Conda won't work - renv won't work - I am tired of fighting with the HPC
# Ideally this would have been a container 
# Will do if I have the time to make one - otherwise...
install.packages("devtools",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("Seurat",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("ggplot2",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("dplyr",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("tidyr",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("ggpubr",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("RColorBrewer",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("remotes",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("imager",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("imagerExtra",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("arrow",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("sf",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("argparser",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("future",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("future.apply",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("data.table",repo = "https://cloud.r-project.org", lib = "lib_cache/")
install.packages("Matrix",repo = "https://cloud.r-project.org", lib = "lib_cache/")
BiocManager::install("rhdf5", lib = "lib_cache/")

remotes::install_github("WonLab-CS/Vesalius@4c0b03c", lib = "lib_cache/") 
remotes::install_github("patrickCNMartin/Kuresi@fabfcbc", lib = "lib_cache/")