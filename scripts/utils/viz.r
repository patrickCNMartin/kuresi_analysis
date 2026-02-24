#library(ggplot2)
library(imager, lib.loc = "lib_cache/")
library(RColorBrewer, lib.loc = "lib_cache/")
library(ggnewscale, lib.loc = "lib_cache/")
library(magick)
library(grid)

view_feature_map <- function(counts,
                             coord,
                             features = NULL,
                             img = NULL,
                             bin_size = 16,
                             max_display_px = 4000) {
  
  # 1. Process Feature Counts
  if (!is.null(features)) {
    counts <- counts[rownames(counts) %in% features, , drop = FALSE]
  }
  
  feature_count <- colSums(counts > 0)
  locs <- match(colnames(counts), coord$barcodes)
  coord$feature_counts <- feature_count[locs]
  
  g <- ggplot()
  
  # 2. Handle Image
  if (!is.null(img) && is(img, "magick-image")) {
    
    info <- image_info(img)
    img_w <- as.numeric(info$width)
    img_h <- as.numeric(info$height)
    
    # Downsample for memory-safe rendering only
    display_scale <- min(1, max_display_px / max(img_w, img_h))
    if (display_scale < 1) {
      img <- image_scale(img, paste0(round(img_w * display_scale), "x"))
    }
    
    img_grob <- grid::rasterGrob(as.raster(img),
                                 interpolate = TRUE,
                                 width  = grid::unit(1, "npc"),
                                 height = grid::unit(1, "npc"))
    
    # Image and coordinates are in the same pixel space — no scale needed
    # Only flip y since ggplot y goes bottom-up, pixels go top-down
    g <- g + annotation_custom(img_grob,
                               xmin = 0, xmax = img_w,
                               ymin = 0, ymax = img_h)
    
    coord$x <- as.numeric(coord$pxl_col_in_fullres)
    coord$y <- img_h - as.numeric(coord$pxl_row_in_fullres)
    
    g <- g +
      scale_x_continuous(limits = c(0, img_w), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, img_h), expand = c(0, 0))
    
  } else {
    coord$x <- coord$array_col
    coord$y <- coord$array_row
  }
  
  # 3. Add Point Layer
  g <- g +
    geom_point(
      data = coord,
      aes(x = x, y = y, color = feature_counts),
      size = 0.05 * bin_size,
      alpha = 0.8
    ) +
    scale_color_distiller(palette = "Spectral") +
    theme_bw() +
    coord_fixed() +
    labs(title = "Features per Cell", color = "Counts", x = "x", y = "y")
  
  return(g)
}

view_count_map <- function(counts,
                           coord,
                           features = NULL, 
                           img = NULL,
                           bin_size = 16,
                           max_display_px = 4000) {
    
    if (!is.null(features)) {
        counts <- counts[rownames(counts) %in% features, , drop = FALSE]
    }
    cell_count <- colSums(counts)
    locs <- match(colnames(counts), coord$barcodes)
    coord$cell_count <- cell_count[locs]
    
    g <- ggplot()
    
    if (!is.null(img) && is(img, "magick-image")) {
        
        info <- image_info(img)
        img_w <- as.numeric(info$width)
        img_h <- as.numeric(info$height)
        
        display_scale <- min(1, max_display_px / max(img_w, img_h))
        if (display_scale < 1) {
            img <- image_scale(img, paste0(round(img_w * display_scale), "x"))
        }
        
        img_grob <- grid::rasterGrob(as.raster(img),
                                     interpolate = TRUE,
                                     width  = grid::unit(1, "npc"),
                                     height = grid::unit(1, "npc"))
        
        g <- g + annotation_custom(img_grob,
                                   xmin = 0, xmax = img_w,
                                   ymin = 0, ymax = img_h)
        
        coord$x <- as.numeric(coord$pxl_col_in_fullres)
        coord$y <- img_h - as.numeric(coord$pxl_row_in_fullres)
        
        g <- g +
            scale_x_continuous(limits = c(0, img_w), expand = c(0, 0)) +
            scale_y_continuous(limits = c(0, img_h), expand = c(0, 0))
        
    } else {
        coord$x <- coord$array_col
        coord$y <- coord$array_row
    }
    
    g <- g +
        geom_point(data = coord,
                   aes(x = x, y = y, col = cell_count),
                   size = 0.05 * bin_size,
                    alpha = 0.8) +
        scale_color_distiller(palette = "Spectral") +
        theme_bw() +
        coord_fixed() +
        labs(title = "Total Counts per Cell", color = "Counts", x = "x", y = "y")
         
    return(g)
}

view_scores <- function(coord,
                        img = NULL,
                        bin_size = 16,
                        max_display_px = 4000,
                        palette = "vesalius") {
  coord$score <- as.factor(coord$score)
  g <- ggplot()
  
  # 2. Handle Image
  if (!is.null(img) && is(img, "magick-image")) {
    
    info <- image_info(img)
    img_w <- as.numeric(info$width)
    img_h <- as.numeric(info$height)
    
    # Downsample for memory-safe rendering only
    display_scale <- min(1, max_display_px / max(img_w, img_h))
    if (display_scale < 1) {
      img <- image_scale(img, paste0(round(img_w * display_scale), "x"))
    }
    
    img_grob <- grid::rasterGrob(as.raster(img),
                                 interpolate = TRUE,
                                 width  = grid::unit(1, "npc"),
                                 height = grid::unit(1, "npc"))
    
    # Image and coordinates are in the same pixel space — no scale needed
    # Only flip y since ggplot y goes bottom-up, pixels go top-down
    g <- g + annotation_custom(img_grob,
                               xmin = 0, xmax = img_w,
                               ymin = 0, ymax = img_h)
    
    coord$x <- as.numeric(coord$x)
    coord$y <- img_h - as.numeric(coord$y)
    
    g <- g +
      scale_x_continuous(limits = c(0, img_w), expand = c(0, 0)) +
      scale_y_continuous(limits = c(0, img_h), expand = c(0, 0))
    
  } 
  
  g <- g +
    geom_point(
      data = coord,
      aes(x = x, y = y, color = score),
      size = 0.05 * bin_size,
      alpha = 0.8
    ) +
    theme_bw() +
    coord_fixed() +
    guides(colour = guide_legend(
        override.aes = list(size = bin_size * 0.8)))
  if (palette == "vesalius") {
    g <- g + scale_color_manual(values = get_ves_pal(coord))+
        labs(title = "Competition Territories", color = "Territory", x = "x", y = "y")
  } else if (palette == "kuresi") {
    g <- g + scale_color_manual(values = get_kur_pal(coord))+
        labs(title = "Competition Territories", color = "Rank", x = "x", y = "y")
  } else {
    stop("Unknown palette")
  }
  return(g)

}

get_ves_pal <- function(score) {
  ter_col <- length(levels(score$score))
  base_colours <- c(
      "#E69F00",
      "#56B4E9",
      "#009E73",
      "#F0E442",
      "#0072B2",
      "#D55E00",
      "#CC79A7",
      "#999999")
  if (ter_col < length(base_colours)) {
      ter_pal <- colorRampPalette(base_colours[seq(1, ter_col)])
  } else {
      ter_pal <- colorRampPalette(base_colours)
  }
  ter_col <- sample(ter_pal(ter_col), ter_col)
  return(ter_col)
}

get_kur_pal <- function(score) {
  ter_col <- length(levels(score$score))
  base_colours <- c(
    "#850101",
    "#cd8878",
    "#f1f1b1",
    "#9CAAC4",
    "#1F3B70"
  )
  if (ter_col < length(base_colours)) {
    ter_pal <- colorRampPalette(base_colours[seq(1, ter_col)])
  } else {
    ter_pal <- colorRampPalette(base_colours)
  }
  ter_col <- ter_pal(ter_col)
  return(ter_col)
}
