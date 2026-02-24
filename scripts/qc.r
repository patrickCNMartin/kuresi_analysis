#!/usr/bin/env Rscript
#-----------------------------------------------------------------------------#
# LIBRARIES
#-----------------------------------------------------------------------------#
library(future, lib.loc = "lib_cache/")
library(future.apply, lib.loc = "lib_cache/")
library(ggplot2)
library(magick)
library(Kuresi, lib.loc = "lib_cache/")
library(ggpubr, lib.loc = "lib_cache/")
library(vesalius, lib.loc = "lib_cache/")
library(argparser, lib.loc = "lib_cache/")
set.seed(1547)

#-----------------------------------------------------------------------------#
# ARGS & OPTIONS
#-----------------------------------------------------------------------------#
p <- arg_parser("Data QC")

p <- add_argument(p, "--sample_name", short = "-sm", help = "Sample Name", type = "character")

p <- add_argument(p, "--counts", short = "-cs", help = "Path to counts", type = "character")

p <- add_argument(p, "--coordinates", short = "-co", help = "Path to coordinates", type = "character")

p <- add_argument(p, "--image", short = "-im", help = "Path to H&E image", type = "character")

p <- add_argument(p, "--scale", short = "-sc", help = "Path to scale factors", type = "character")

p <- add_argument(p, "--cancer_type", short = "-ct", help = "Cancer type", type = "character")

p <- add_argument(p, "--count_threshold", short = "-cot", help = "Gene Count threshold", type = "numeric")

p <- add_argument(p, "--bin_size", short = "-r", help = "Number of bins (assuming 2um base data)", type = "numeric")

p <- add_argument(p, "--min_features", short = "-mf", help = "Min Number of Features per cells", type = "numeric")

p <- add_argument(p, "--min_cells", short = "-mc", help = "Min number of cells per feature", type = "numeric")

p <- add_argument(p, "--output_dir", short = "-od", help = "Output Directory for results and meta data", type = "character")

p <- add_argument(p, "--report_file", short = "-rf", help = "Report File to generate", type = "character")

argv <- parse_args(p)

sample_name <- argv$sample_name
counts <- argv$counts
coordinates <- argv$coordinates
image <- argv$image
scale <- argv$scale
cancer_type <- argv$cancer_type
count_threshold <- argv$count_threshold
bin_size <- argv$bin_size
min_features <- argv$min_features
min_cells <- argv$min_cells
output_dir <- argv$output_dir
report_file <- argv$report_file

source("scripts/utils/io.r")
source("scripts/utils/utils.r")
source("scripts/utils/viz.r")
#-----------------------------------------------------------------------------#
# Load data
#-----------------------------------------------------------------------------#
cat("QC Start \n")
data <- load_visiumhd(counts, coordinates, image , scale)
# Win lose sets
gene_sets <- win_lose_genes()
win_genes <- gene_sets$win
lose_genes <- gene_sets$lose
features <- c(win_genes, lose_genes)

# cancer marker set
cancer <- cancer_maker_list(cancer_type)


cat("Data Loading - Completed\n")
#-----------------------------------------------------------------------------#
# Binning data
#-----------------------------------------------------------------------------#
aggregated_data <- bin_spatial_data(data$counts, data$coord, bin_size = bin_size)

#-----------------------------------------------------------------------------#
# Count Filtering
#-----------------------------------------------------------------------------#
cancer_counts <- filter_counts(aggregated_data$counts,
                        min_features,
                        min_cells,
                        cancer,
                        count_threshold,
                        return_cancer = TRUE)
cancer_coord <- aggregated_data$coords[aggregated_data$coords$barcodes %in% colnames(cancer_counts), ]
counts <- filter_counts(aggregated_data$counts,
                        min_features,
                        min_cells,
                        cancer,
                        count_threshold)
coord <- aggregated_data$coords[aggregated_data$coords$barcodes %in% colnames(counts), ]
#-----------------------------------------------------------------------------#
# Get Cancer cells
#-----------------------------------------------------------------------------#
qc_1 <- view_feature_map(cancer_counts,
                         cancer_coord,
                         img = data$image,
                         bin_size = bin_size)
qc_2 <- view_count_map(counts,
                       coord,
                       img = data$image,
                       bin_size = bin_size)
ggsave(file.path(output_dir, "cancer_feature_map.tiff"), plot = qc_1, width = 8, height = 8, units = "in")
ggsave(file.path(output_dir, "cancer_count_map.tiff"), plot = qc_2, width = 8, height = 8, units = "in")
#-----------------------------------------------------------------------------#
# Get QC plots
#-----------------------------------------------------------------------------#
qc_3 <- view_feature_map(counts,
                         coord,
                         img = data$image,
                         bin_size = bin_size)
qc_4 <- view_count_map(counts,
                       coord,
                       img = data$image,
                       bin_size = bin_size)
ggsave(file.path(output_dir, "feature_map.tiff"), plot = qc_3, width = 12, height = 12, units = "in")
ggsave(file.path(output_dir, "count_map.tiff"), plot = qc_4, width = 12, height = 12, units = "in")

#-----------------------------------------------------------------------------#
# Get QC plots of only features of interest
#-----------------------------------------------------------------------------#
qc_5 <- view_feature_map(counts,
                         coord,
                         features,
                         img = data$image,
                         bin_size = bin_size)
qc_6 <- view_count_map(counts,
                       coord,
                       features,
                       img = data$image,
                       bin_size = bin_size)
ggsave(file.path(output_dir, "feature_map_win_loose.tiff"), plot = qc_5, width = 8, height = 8, units = "in")
ggsave(file.path(output_dir, "count_map_win_loose.tiff"), plot = qc_6, width = 8, height = 8, units = "in")
#-----------------------------------------------------------------------------#
# Export metrics
#-----------------------------------------------------------------------------#
mean_counts <- mean(colSums(counts))
mean_genes  <- mean(colSums(counts > 0))
# Write QC report
report_text <- sprintf(
    "QC Report for %s\n%s\nMetrics:\n  Cells: %d\n  Genes: %d\n  Mean counts/cell: %.1f\n  Mean genes/cell: %.0f\n ",
    sample_name,
    strrep("=", 50),
    ncol(counts),
    nrow(counts),
    mean_counts,
    mean_genes
)

writeLines(report_text, con = file.path(report_file))
#-----------------------------------------------------------------------------#
# Export data
# Not exporting image since magick pointers don't persist across sessions
# Note: not scale here either since the scale factors here are pointless
#-----------------------------------------------------------------------------#
out <- list("counts" = counts , "coord" = coord)
saveRDS(out, file = file.path(output_dir,"qc_data.rds"))