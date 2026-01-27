#!/usr/bin/env Rscript
#-----------------------------------------------------------------------------#
# LIBRARIES
#-----------------------------------------------------------------------------#
library(Kuresi)
library(argparser)

#-----------------------------------------------------------------------------#
# ARGS & OPTIONS
#-----------------------------------------------------------------------------#
p <- arg_parser("Kuresi Ratio Score")

p <- add_arguments(p, "--counts", short = "-c", help = "Path to RNA count Matrix", type = "character")
p <- add_arguments(p, "--coordinates", short = "-co", help = "Path to Sptial Coordinates", type = "character")
p <- add_arguments(p, "--method", short = "-m", help = "Ratio Method to use", type = "character")
p <- add_arguments(p, "--resolution", short = "-r", help = "Image resolution factor", type = "numeric")
p <- add_arguments(p, "--dim_reduc", short = "-d", help = "Dimensionality Reduction Method", type = "character")
p <- add_arguments(p, "--dims", short = "-di", help = "Number of Latent space dimensions", type = "numeric")
p <- add_arguments(p, "--sigma", short = "-s", help = "Sigma Factor for smoothing", type = "numeric")
p <- add_arguments(p, "--iter", short = "-i", help = "Number of smoothing iteration", type = "numeric")
p <- add_arguments(p, "--col_resolution", short = "-cr", help = "Number of color segments", type = "numeric")
p <- add_arguments(p, "--distance", short = "-df", help = "Distance factor for territory pooling", type = "numeric")

argv <- parse_args(p)
counts <- argv$input
coord <- argv$coordinates


max_size <- 100000 * 1024^2
options(future.globals.maxSize = max_size)