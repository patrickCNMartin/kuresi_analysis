library(arrow, lib.loc = "lib_cache/")
library(Seurat, lib.loc = "lib_cache/")
library(imager, lib.loc = "lib_cache/")




load_visiumhd <- function(
    counts,
    coordinates,
    image,
    scale_factor,
    bin_size = 16,
    downsample = NULL)
{
    # Coordinate loading
    if (grepl("parquet", coordinates)) {
        coordinates <- arrow::read_parquet(coordinates)
        coordinates <- coordinates[coordinates$in_tissue == 1,]
        coordinates <- coordinates[order(coordinates$barcode, decreasing = FALSE),]
    } else {
        stop("Please Use Binned outputs")
    }
    # Counts
    if (grepl("h5", counts)) {
        counts <- Seurat::Read10X_h5(counts)
        counts <- counts[,order(colnames(counts), decreasing = FALSE)]
    } else {
        stop("Please use h5 files")
    }
    # Image - for now assume Null
    image <- NULL
    scale <- "auto"
    # Downsample for testing purposes
    if (!is.null(downsample) && is(downsample, "numeric")) {
        barcodes <- sort(sample(coordinates$barcodes, downsample, replace = FALSE))
        coordinates <- coordinates[coordinates$barcodes %in% barcodes, ]
        counts <- counts[, colnames(counts) %in% barcodes]
    }
    return(list("counts" = counts,
                "coordinates" = coordinates,
                "image" = image,
                "scale" = scale))
}