#library(ggplot2)
library(imager)
library(RColorBrewer)



view_feature_map <- function(counts, coord, features = NULL) {
    if (!is.null(features)) {
        counts <- counts[rownames(counts) %in% features, ]
    }
    feature_count <- colSums(counts > 0)
    locs <- match(colnames(counts), coord$barcodes)
    coord$feature_counts <- feature_count[locs]
    str(coord)
    g <- ggplot(coord, aes(x = x , y = y, col = feature_counts)) +
         geom_point(size = 0.2, stroke = 0) +
         scale_color_distiller(palette = "Spectral") +
         theme_bw() +
         labs(title = "Features per Cell")
    return(g)
}

view_count_map <- function(counts, coord, features = NULL) {
    if (!is.null(features)) {
        counts <- counts[rownames(counts) %in% features, ]
    }
    cell_count <- colSums(counts)
    locs <- match(colnames(counts), coord$barcodes)
    coord$cell_count <- cell_count[locs]
    str(coord)
    g <- ggplot(coord, aes(x = x , y = y, col = cell_count)) +
         geom_point(size = 0.2, stroke = 0) +
         scale_color_distiller(palette = "Spectral") +
         theme_bw() +
         labs(title = "Total Counts per Cell")
    return(g)
}