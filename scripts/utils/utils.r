library(data.table, lib.loc = "lib_cache/")
library(Matrix, lib.loc = "lib_cache/")


filter_counts <- function(counts, min_features, min_cells) {
  
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
  dt <- as.data.table(coordinates)
  dt <- dt[barcodes %in% colnames(counts)]
  
  # Ensure the count matrix columns are in the exact same order as our coordinate rows
  counts <- counts[, dt$barcodes]
  
  # 2. Assign Bin IDs based on spatial proximity
  dt[, bin_x := (x %/% bin_size) * bin_size]
  dt[, bin_y := (y %/% bin_size) * bin_size]
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
  binned_counts <- counts %*% mapping
  
  # Convert sums to averages by dividing by the number of cells per bin
  cells_per_bin <- colSums(mapping)
  binned_counts <- t(t(binned_counts) / cells_per_bin)
  
  # 5. Aggregate Coordinates (Average of original x/y)
  binned_coords <- dt[, .(
    avg_x = mean(x), 
    avg_y = mean(y), 
    cell_count = .N
  ), by = bin_id][order(bin_id)]
  
  # Final Correspondence: Set column names of matrix to the bin_id
  colnames(binned_counts) <- binned_coords$bin_id
  
  return(list(
    counts = binned_counts,
    coords = binned_coords
  ))
}