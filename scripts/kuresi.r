#!/usr/bin/env Rscript
#-----------------------------------------------------------------------------#
# LIBRARIES
#-----------------------------------------------------------------------------#
library(ggplot2)
library(magick)
library(Kuresi, lib.loc = "lib_cache/")
library(future, lib.loc = "lib_cache/")
library(future.apply, lib.loc = "lib_cache/")
library(Seurat, lib.loc = "lib_cache/")
library(vesalius, lib.loc = "lib_cache/")
library(argparser, lib.loc = "lib_cache/")
set.seed(1547)
#-----------------------------------------------------------------------------#
# ARGS & OPTIONS
#-----------------------------------------------------------------------------#
p <- arg_parser("Kuresi Ratio Score\n")

p <- add_argument(p, "--sample_name", short = "-sm", help = "Sample Name", type = "character")
p <- add_argument(p, "--input_rds", short = "-i", help = "Vesalius Object location", type = "character")
p <- add_argument(p, "--image", short = "-im", help = "Path to image", type = "character")
p <- add_argument(p, "--method", short = "-m", help = "Kuresi Method", type = "character")
p <- add_argument(p, "--scale", short = "-s", help = "Scale Ratio Score", type = "logical")
p <- add_argument(p, "--center", short = "-c", help = "Center Ratio Score", type = "logical")
p <- add_argument(p, "--rank", short = "-r", help = "Rank Ratio Score", type = "logical")
p <- add_argument(p, "--bin_size", short = "-b", help = "Bin Size", type = "numeric")
p <- add_argument(p, "--output_dir", short = "-o", help = "Output Kuresi Scores", type = "character")
p <- add_argument(p, "--report_file", short = "-rf", help = "Report File to generate", type = "character")

argv <- parse_args(p)
sample_name <- argv$sample_name
input_rds <- argv$input_rds
image <- argv$image
method <- argv$method
scale <- argv$scale
center <- argv$center
rank <- argv$rank
bin_size <- argv$bin_size
output_dir <- argv$output_dir
report_file <- argv$report_file
source("scripts/utils/viz.r")
source("scripts/utils/utils.r")
#-----------------------------------------------------------------------------#
# INPUT
#-----------------------------------------------------------------------------#
vesalius <- readRDS(input_rds)
counts <- vesalius@counts$log_norm
territories <- vesalius@territories
image <- magick::image_read(image)
gene_sets <- win_lose_genes()
win_genes <- gene_sets$win
lose_genes <- gene_sets$lose
#-----------------------------------------------------------------------------#
# Compute Scores
#-----------------------------------------------------------------------------#
if (method == "ELO") {
    kuresi_score <- compute_competition_outcomes(counts,
                                                 territories,
                                                 win_genes,
                                                 group_name = "Territory",
                                                 gene_set2 = lose_genes,
                                                 log_fc = 0.2)
} else {
    kuresi_score <- compute_ratio_bygroup(counts,
                                          territories,
                                          group_name = "Territory",
                                          genes_1 = win_genes,
                                          genes_2 = lose_genes,
                                          method = method,
                                          scale = scale,
                                          center = center,
                                          rank = rank,
                                          add_name = method,
                                          verbose = FALSE)
}
#-----------------------------------------------------------------------------#
# Kuresi report
#-----------------------------------------------------------------------------#
max_score <- max(kuresi_score[,method])
min_score  <- min(kuresi_score[,method])
mean_score <- mean(kuresi_score[,method])
med_score <- median(kuresi_score[,method])
n_score <- length(unique(kuresi_score[,method]))
# Write QC report
report_text <- sprintf(
    "Kuresi Report for %s\n%s\nMetrics:\n  Max Score: %.0f\n  Min Score: %.0f\n  Mean Score: %.1f\n  Median Score: %.1f\n  Score Count: %.0f\n",
    sample_name,
    strrep("=", 50),
    max_score,
    min_score,
    mean_score,
    med_score,
    n_score
)

writeLines(report_text, con = file.path(report_file))
#-----------------------------------------------------------------------------#
# Plotting Kuresi
#-----------------------------------------------------------------------------#
if (method == "ELO") {
    # PLACE HOLDER
} else {
    kuresi_score <- kuresi_score[,c("x","y", method)]
    colnames(kuresi_score) <- gsub(method,"score",colnames(kuresi_score))
}
kuresi_score_orig <- set_original_coordinates(vesalius,kuresi_score)
kur <- view_scores(kuresi_score_orig,img = image, bin_size, palette = "kuresi")
ggsave(file.path(output_dir, "kuresi_score_plot.tiff"), plot = kur, width = 8, height = 8, units = "in")
#-----------------------------------------------------------------------------#
# Export and Save
#-----------------------------------------------------------------------------#
saveRDS(kuresi_score, file = file.path(output_dir, "kuresi_competition_scores.rds"))
#-----------------------------------------------------------------------------#
# DONE
#-----------------------------------------------------------------------------#
cat("Kuresi Competition Score - Completed")