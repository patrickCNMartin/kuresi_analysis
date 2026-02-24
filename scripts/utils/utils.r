library(data.table, lib.loc = "lib_cache/")
library(Matrix, lib.loc = "lib_cache/")


filter_counts <- function(counts,
    min_features,
    min_cells,
    cancer = NULL,
    count_threshold = NULL,
    return_cancer = FALSE) {
  
  if(!is.null(cancer) && !is.null(count_threshold)) {
    cancer_cells <- counts[cancer, ]
    cancer_cells <- which(Matrix::colSums(cancer_cells) > count_threshold)
    counts <- counts[,cancer_cells]
    if (return_cancer)  {
      return(counts)
    }
  }
  genes_per_cell <- Matrix::colSums(counts > 0)
  keep_cells <- genes_per_cell >= min_features
  counts <- counts[, keep_cells]
  
  cells_per_gene <- Matrix::rowSums(counts > 0)
  keep_genes <- cells_per_gene >= min_cells
  counts <- counts[keep_genes, ]
  
  return(counts)
}


bin_spatial_data <- function(counts, coordinates, bin_size = 16) {
  # 1. Align data: Ensure coordinates match the matrix columns
  # Using 'barcode' (singular) as per your tibble structure
  dt <- as.data.table(coordinates)
  dt <- dt[barcode %in% colnames(counts)]
  
  # Ensure the count matrix columns are in the exact same order as our coordinate rows
  # This is critical for the matrix multiplication step
  counts <- counts[, dt$barcode]
  
  # 2. Assign Bin IDs based on array coordinates
  dt[, bin_x := (array_row %/% bin_size) * bin_size]
  dt[, bin_y := (array_col %/% bin_size) * bin_size]
  dt[, bin_id := .GRP, by = .(bin_x, bin_y)]
  
  # 3. Create the Mapping Matrix (Cells to Bins)
  # Rows = Cells, Cols = Bins. A '1' indicates a cell belongs to a bin.
  mapping <- sparseMatrix(
    i = seq_len(nrow(dt)),
    j = dt$bin_id,
    x = 1,
    dims = c(nrow(dt), max(dt$bin_id))
  )
  
  # 4. Matrix Multiplication for Speed
  # [Genes x Cells] %*% [Cells x Bins] = [Genes x Bins]
  # This sums the counts for all cells within each bin
  binned_counts <- counts %*% mapping
  
  # 5. Compute Mean Locations for all coordinate types
  binned_coords <- dt[, .(
    array_row = mean(array_row),
    array_col = mean(array_col),
    pxl_row_in_fullres = mean(pxl_row_in_fullres),
    pxl_col_in_fullres = mean(pxl_col_in_fullres),
    cell_count = .N
  ), by = .(bin_id)][order(bin_id)]
  
  # Create a 'barcodes' column using the bin_id to identify each new bin
  binned_coords[, barcodes := as.character(bin_id)]
  
  # Set matrix colnames to match the new binned coordinate identifiers
  colnames(binned_counts) <- binned_coords$barcodes
  
  # Final formatting of the coordinate table
  binned_coords <- binned_coords[, .(
    barcodes, 
    array_row, 
    array_col, 
    pxl_row_in_fullres, 
    pxl_col_in_fullres,
    cell_count
  )]
  
  return(list(
    counts = binned_counts,
    coords = binned_coords
  ))
}



get_orig_coordinates <- function(vesalius_assay) {
    trial <- vesalius_assay@territories[,c("barcodes","x","y","Territory")]
    orig_coord <- vesalius_assay@meta$orig_coord
    
    in_trial <- match(orig_coord$barcodes, trial$barcodes) 
    
    trial$x <- orig_coord$x_orig[in_trial]
    trial$y <- orig_coord$y_orig[in_trial]
    colnames(trial) <- gsub("Territory","score",colnames(trial))
    trial <- trial[,c("x","y","score")]
    return(trial)
}

set_original_coordinates <- function(vesalius_assay, score) {
    orig_coord <- vesalius_assay@meta$orig_coord
    in_trial <- match(orig_coord$barcodes, score$barcodes) 
    score$x <- orig_coord$x_orig[in_trial]
    score$y <- orig_coord$y_orig[in_trial]
    score <- score[,c("x","y","score")]
    return(score)

}