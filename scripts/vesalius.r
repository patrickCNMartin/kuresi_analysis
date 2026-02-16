#!/usr/bin/env Rscript
# cheat libs because the hpc configs sucks and I don't have time for anything else
library(ggplot2)
#-----------------------------------------------------------------------------#
# LIBRARIES
#-----------------------------------------------------------------------------#
library(Kuresi, lib.loc = "lib_cache/")
library(future, lib.loc = "lib_cache/")
library(future.apply, lib.loc = "lib_cache/")
library(Seurat, lib.loc = "lib_cache/")
library(vesalius, lib.loc = "lib_cache/")
library(ggpubr, lib.loc = "lib_cache/")
library(argparser, lib.loc = "lib_cache/")
set.seed(1547)
#-----------------------------------------------------------------------------#
# ARGS & OPTIONS
#-----------------------------------------------------------------------------#
p <- arg_parser("Vesalius Territory Detection")

p <- add_argument(p, "--sample_name", short = "-sm", help = "Sample Name", type = "character")

p <- add_argument(p, "--input_rds", short = "-i", help = "Path to rds", type = "character")

p <- add_argument(p, "--feature_set", short = "-fs", help = "Use full feature set", type = "logical")

p <- add_argument(p, "--resolution", short = "-r", help = "Image Resolution", type = "numeric")

p <- add_argument(p, "--dim_reduc", short = "-d", help = "Dim reduction method", type = "character")

p <- add_argument(p, "--dims", short = "-di", help = "Number of latent space dimensions", type = "numeric")

p <- add_argument(p, "--equalize", short = "-eq", help = "Piecewise eq percentage", type = "numeric")

p <- add_argument(p, "--sigma", short = "-s", help = "Sigma for smoothing", type = "numeric")

p <- add_argument(p, "--iter", short = "-it", help = "Smoothing iterations", type = "numeric")

p <- add_argument(p, "--col_resolution", short = "-cr", help = "Number of color segments", type = "numeric")

p <- add_argument(p, "--distance", short = "-ds", help = "Territory pooling distance", type = "numeric")

p <- add_argument(p, "--output_dir", short = "-od", help = "Output Directory for results and meta data", type = "character")

p <- add_argument(p, "--report_file", short = "-rf", help = "Report File to generate", type = "character")

p <- add_argument(p, "--cores", short = "-c", help = "Number of cores", type = "numeric")

argv <- parse_args(p)

sample_name <- argv$sample_name
rds <- argv$input_rds
feature_set <- argv$feature_set
resolution <- argv$resolution
dim_reduc <- argv$dim_reduc
dims <- argv$dims
eq <- argv$equalize
sigma <- argv$sigma
iter <- argv$iter
col_resolution <- argv$col_resolution
distance <- argv$distance
output_dir <- argv$output_dir
report_file <- argv$report_file
cores <- argv$cores

max_size <- 10000 * 1024^2
options(future.globals.maxSize = max_size)
plan(multicore, workers = cores)
#-----------------------------------------------------------------------------#
# INPUT
#-----------------------------------------------------------------------------#
cat("Vesalius Start \n")
if (feature_set) {
    gene_sets <- win_lose_genes()
    win_genes <- gene_sets$win
    lose_genes <- gene_sets$lose    
    features <- c(win_genes, lose_genes)
} else {
    features <- NULL  
}

data <- readRDS(rds)

coordinates <- as.data.frame(data$coord)
cat("Data Loading - Completed\n")
#-----------------------------------------------------------------------------#
# Vesalius Territories
#-----------------------------------------------------------------------------#
vesalius <- build_vesalius_assay(coordinates,
                                 data$counts,
                                 data$image,
                                 assay = sample_name,
                                 scale = data$scale,
                                 verbose = FALSE)
cat("Object Build - Completed\n")

vesalius <- generate_embeddings(vesalius,
                                dim_reduction = dim_reduc,
                                dimensions = dims,
                                tensor_resolution = resolution,
                                features = features,
                                verbose = FALSE)
cat("Generate Embeddings - Completed\n")

vesalius <- equalize_image(vesalius,
                           sleft = eq,
                           sright = eq,
                           verbose = FALSE)
cat("Image Equalization - Completed\n")

vesalius <- smooth_image(vesalius,
                        dimensions = seq(1,dims),
                        sigma = sigma,
                        iter = iter,
                        verbose = FALSE)
cat("Image Smoothing - Completed\n")

vesalius <- segment_image(vesalius,
                          dimensions = seq(1, dims),
                          col_resolution = col_resolution,
                          verbose = FALSE)
cat("Image Segmentation - Completed\n")

vesalius <- isolate_territories(vesalius,
                                capture_radius = distance,
                                verbose = FALSE)
cat("Territory Pooling - Completed\n")

#-----------------------------------------------------------------------------#
# Plotting Vesalius
#-----------------------------------------------------------------------------#
img <- image_plot(vesalius)
ter <- territory_plot(vesalius, cex_pt = 2)
ggsave(file.path(output_dir, "vesalius_image_plot.pdf"), plot = img, width = 12, height = 12, units = "in")
ggsave(file.path(output_dir, "vesalius_territory_plot.pdf"), plot = ter, width = 12, height = 12, units = "in")
#-----------------------------------------------------------------------------#
# Export Metrics
#-----------------------------------------------------------------------------#
n_cells <- nrow(vesalius@territories)
n_territories <- length(unique(vesalius@territories$Territory))
isolated <- "isolated" %in% vesalius@territories$Territory
# Write QC report
report_text <- sprintf(
    "Vesalius Report for %s\n%s\nMetrics:\n  Cells: %d\n  Territories: %d\n  Isolated Territories: %s\n",
    sample_name,
    strrep("=", 50),
    n_cells,
    n_territories,
    isolated
)

writeLines(report_text, con = file.path(report_file))
#-----------------------------------------------------------------------------#
# Export and Save
#-----------------------------------------------------------------------------#
saveRDS(vesalius, file = file.path(output_dir,"vesalius_data.rds"))
#-----------------------------------------------------------------------------#
# DONE
#-----------------------------------------------------------------------------#
cat("Vesalius Territory Detection - Completed\n")