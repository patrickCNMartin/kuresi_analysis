#!/usr/bin/env Rscript
#-----------------------------------------------------------------------------#
# LIBRARIES
#-----------------------------------------------------------------------------#
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
p <- add_argument(p, "--method", short = "-m", help = "Kuresi Method", type = "character")
p <- add_argument(p, "--scale", short = "-s", help = "Scale Ratio Score", type = "logical")
p <- add_argument(p, "--center", short = "-c", help = "Center Ratio Score", type = "logical")
p <- add_argument(p, "--output_dir", short = "-o", help = "Output Kuresi Scores", type = "character")

argv <- parse_args(p)
sample_name <- argv$sample_name
input_rds <- argv$input_rds
territories <- argv$territories
method <- argv$method
scale <- argv$scale
center <- argv$center
output_dir <- argv$output_dir

#-----------------------------------------------------------------------------#
# INPUT
#-----------------------------------------------------------------------------#
vesalius <- readRDS(input_rds)
counts <- vesalius@counts$log_norm
territories <- vesalius@territories
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
    "Kuresi Report for %s\n%s\nMetrics:\n  Max Score: %d\n  Min Score: %d\n  Mean Score: %.1f\n  Median Score: %.1f\n  Score Count: %d\n",
    sample_name,
    strrep("=", 50),
    max_score,
    min_score,
    mean_score,
    med_score,
    n_score
)

writeLines(report_text, con = file.path(output_dir, "kuresi_report.txt"))
#-----------------------------------------------------------------------------#
# Plotting Kuresi
#-----------------------------------------------------------------------------#
if (method == "ELO") {
    # PLACE HOLDER
} else {
    kuresi_score <- kuresi_score[,c("x","y", method)]
}

kur <- score_plot(kuresi_score)
ggsave(file.path(output_dir, "kuresi_score_plot.pdf"), plot = kur, width = 12, height = 12, units = "in")
#-----------------------------------------------------------------------------#
# Export and Save
#-----------------------------------------------------------------------------#
saveRDS(kuresi_score, file = file.path(output_dir, "kuresi_competition_scores.rds"))
#-----------------------------------------------------------------------------#
# DONE
#-----------------------------------------------------------------------------#
cat("Kuresi Competition Score - Completed")