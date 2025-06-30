#' Add cell type hierarchy to BayesCNet object
#'
#' @param object A BayesCNet object
#' @param edges Character vector of edges in format c("parent1", "child1", "parent2", "child2", ...)
#' @param lambda Decay parameter for tree-based distance (default: 1.0)
#'
#' @return Updated BayesCNet object with hierarchy
#' @export
#'
#' @examples
#' \dontrun{
#' edges <- c("HSPC", "T Cells", "HSPC", "B Cells", "T Cells", "CD4 T", "T Cells", "CD8 T")
#' bcnet <- AddCellTypeHierarchy(bcnet, edges, lambda = 1.0)
#' }
AddCellTypeHierarchy <- function(object, edges, lambda = 1.0) {

  # Validate inputs
  if (!inherits(object, "BayesCNet")) {
    stop("object must be a BayesCNet object")
  }

  if (length(edges) %% 2 != 0) {
    stop("edges must have even number of elements (pairs of parent-child)")
  }

  if (lambda <= 0) {
    stop("lambda must be positive")
  }

  # Get cell types from data
  data_cell_types <- unique(object@cell_metadata$cell_type)

  # Create graph from edges
  message("Building cell type hierarchy...")
  g <- igraph::graph(edges, directed = TRUE)

  # Get all cell types in hierarchy
  hierarchy_cell_types <- unique(c(edges))

  # Check if all data cell types are in hierarchy
  missing_types <- setdiff(data_cell_types, hierarchy_cell_types)
  if (length(missing_types) > 0) {
    stop("Cell types in data not found in hierarchy: ",
         paste(missing_types, collapse = ", "))
  }

  # Compute distance matrix using shortest paths
  dist_matrix <- igraph::shortest.paths(g, mode = "all")

  # Apply exponential decay to create weighted adjacency
  adjacency_matrix <- exp(-lambda * dist_matrix)
  diag(adjacency_matrix) <- 1

  # Subset to only cell types present in data
  adjacency_matrix <- adjacency_matrix[data_cell_types, data_cell_types]

  # Compute degree matrix
  degree_matrix <- diag(rowSums(adjacency_matrix))

  # Compute graph Laplacian
  laplacian <- degree_matrix - adjacency_matrix

  # Compute covariance matrix (used as prior in inference)
  n_types <- length(data_cell_types)
  Sigma2 <- solve(laplacian + (1/n_types) * diag(n_types))

  # Standardize the covariance matrix
  diag_elements <- diag(Sigma2)
  Sigma2 <- Sigma2 / sqrt(diag_elements %*% t(diag_elements))

  # Store in object
  object@cell_type_hierarchy <- list(
    graph = g,
    adjacency = adjacency_matrix,
    laplacian = laplacian,
    lambda = lambda,
    sigma2 = Sigma2,
    cell_types = data_cell_types
  )

  object@parameters$hierarchy_lambda <- lambda

  message("Cell type hierarchy added successfully")
  message("  ", length(data_cell_types), " cell types in hierarchy")

  return(object)
}

#' Add gene annotation to BayesCNet object
#'
#' @param object A BayesCNet object
#' @param annotation Data frame with columns: gene, chr, tss
#' @param species Species for automatic annotation (not implemented yet)
#'
#' @return Updated BayesCNet object with gene annotation
#' @export
#'
#' @examples
#' \dontrun{
#' # annotation should have columns: gene, chr, tss
#' bcnet <- AddGeneAnnotation(bcnet, gene_annotation_df)
#' }
AddGeneAnnotation <- function(object, annotation, species = NULL) {

  # Validate object
  if (!inherits(object, "BayesCNet")) {
    stop("object must be a BayesCNet object")
  }

  # Validate annotation
  if (!is.data.frame(annotation)) {
    stop("annotation must be a data frame")
  }

  required_cols <- c("gene", "chr", "tss")
  missing_cols <- setdiff(required_cols, colnames(annotation))
  if (length(missing_cols) > 0) {
    stop("annotation missing required columns: ",
         paste(missing_cols, collapse = ", "))
  }

  # Check data types
  if (!is.character(annotation$gene)) {
    stop("gene column must be character")
  }
  if (!is.character(annotation$chr)) {
    stop("chr column must be character")
  }
  if (!is.numeric(annotation$tss)) {
    stop("tss column must be numeric")
  }

  # Get genes in RNA data
  rna_genes <- rownames(object@RNA)

  # Check overlap
  common_genes <- intersect(rna_genes, annotation$gene)
  if (length(common_genes) == 0) {
    stop("No genes found in both RNA data and annotation")
  }

  # Filter annotation to genes in data
  annotation <- annotation[annotation$gene %in% common_genes, ]

  # Remove duplicates (keep first occurrence)
  if (any(duplicated(annotation$gene))) {
    warning("Duplicate genes found in annotation. Keeping first occurrence.")
    annotation <- annotation[!duplicated(annotation$gene), ]
  }

  # Ensure chr format is consistent
  if (!all(grepl("^chr", annotation$chr))) {
    message("Adding 'chr' prefix to chromosome names")
    annotation$chr <- paste0("chr", annotation$chr)
  }

  # Store annotation
  object@gene_annotation <- annotation

  message("Gene annotation added successfully")
  message("  ", nrow(annotation), " genes annotated")
  message("  ", length(setdiff(rna_genes, annotation$gene)),
          " genes in RNA data without annotation")

  return(object)
}

#' Aggregate cells by KNN
#'
#' @param object A BayesCNet object
#' @param k Number of nearest neighbors (default: 50)
#' @param max_overlap Maximum overlap between metacells (default: 0.8)
#' @param seed Random seed (default: 123)
#' @param normalize Whether to normalize aggregated data (default: TRUE)
#'
#' @return Updated BayesCNet object with aggregated data
#' @export
AggregateByKNN <- function(object, k = 50, max_overlap = 0.8, seed = 123,
                           normalize = TRUE) {

  # Validate inputs
  if (!inherits(object, "BayesCNet")) {
    stop("object must be a BayesCNet object")
  }

  if (k < 2) {
    stop("k must be at least 2")
  }

  if (max_overlap < 0 || max_overlap > 1) {
    stop("max_overlap must be between 0 and 1")
  }

  # Check if we have enough cells per cell type
  cell_counts <- table(object@cell_metadata$cell_type)
  too_small <- names(cell_counts)[cell_counts < k]
  if (length(too_small) > 0) {
    stop("Cell types with fewer than k cells: ",
         paste(too_small, collapse = ", "))
  }

  message("Aggregating cells by KNN (k = ", k, ")...")

  # Check if embeddings exist
  if (nrow(object@cell_embeddings) == 0) {
    stop("No cell embeddings found. This should not happen - please report this bug.")
  }

  # Use the stored embeddings
  cell_coord <- object@cell_embeddings

  # Get cell type groupings from metadata
  cell_types <- object@cell_metadata$cell_type
  unique_cell_types <- unique(cell_types)

  # Initialize aggregated data matrices
  rna_new <- matrix(0, nrow = nrow(object@RNA), ncol = 1)
  atac_new <- matrix(0, nrow = nrow(object@ATAC), ncol = 1)
  cell_sample <- matrix(0, nrow = 1, ncol = k + 1)

  # Process each cell type separately
  for (i in 1:length(unique_cell_types)) {
    current_type <- unique_cell_types[i]
    message("  Aggregating cell type: ", current_type)

    # Get indices for current cell type
    type_indices <- which(cell_types == current_type)

    # Get coordinates for this cell type
    type_coords <- cell_coord[type_indices, , drop = FALSE]

    # Create subset data for this cell type
    rna_subset <- object@RNA[, type_indices, drop = FALSE]
    atac_subset <- object@ATAC[, type_indices, drop = FALSE]

    # Aggregate this cell type
    agg_result <- aggregate_single_celltype(
      rna_data = rna_subset,
      atac_data = atac_subset,
      cell_coord = type_coords,
      k_neigh = k,
      max_overlap = max_overlap,
      seed = seed + i,  # Different seed for each cell type
      verbose = FALSE
    )

    # Append results
    rna_new <- cbind(rna_new, agg_result$rna)
    atac_new <- cbind(atac_new, agg_result$atac)

    # Adjust cell indices to global indices
    cell_indices <- agg_result$cell_sample
    cell_indices[cell_indices > 0] <- type_indices[cell_indices[cell_indices > 0]]

    # Add cell type label
    cell_type_col <- matrix(current_type,
                            nrow = nrow(cell_indices),
                            ncol = 1)
    cell_indices <- cbind(cell_indices, cell_type_col)

    cell_sample <- rbind(cell_sample, cell_indices)
  }

  # Remove initialization columns/rows
  rna_new <- rna_new[, -1, drop = FALSE]
  atac_new <- atac_new[, -1, drop = FALSE]
  cell_sample <- cell_sample[-1, , drop = FALSE]

  # Normalize if requested
  if (normalize) {
    message("Normalizing aggregated data...")
    # CPM normalization for RNA
    rna_new <- sweep(rna_new, 2, colSums(rna_new), FUN = "/") * 1e6
    # ATAC is already binarized, no need for additional normalization
  }

  # Store aggregated data
  object@aggregated_data <- list(
    rna = rna_new,
    atac = atac_new,
    cell_assignments = cell_sample,
    k_neighbors = k,
    normalized = normalize
  )

  object@parameters$aggregation_k <- k
  object@parameters$aggregation_max_overlap <- max_overlap
  object@parameters$aggregation_normalized <- normalize

  message("Aggregation complete")
  message("  ", ncol(rna_new), " metacells created")

  return(object)
}

#' Select variable genes
#'
#' @param object A BayesCNet object
#' @param method Method for selection: "markers" or "custom"
#' @param genes Custom gene list (if method = "custom")
#' @param min.pct Minimum percentage of cells expressing gene (default: 0.1)
#' @param logfc.threshold Log fold change threshold (default: 0.1)
#'
#' @return Updated BayesCNet object with selected genes
#' @export
AddVariableGenes <- function(object, method = "markers", genes = NULL,
                             min.pct = 0.1, logfc.threshold = 0.1) {

  # Validate inputs
  if (!inherits(object, "BayesCNet")) {
    stop("object must be a BayesCNet object")
  }

  if (!method %in% c("markers", "custom")) {
    stop("method must be 'markers' or 'custom'")
  }

  # Check if gene annotation is added
  if (nrow(object@gene_annotation) == 0) {
    stop("Gene annotation must be added before selecting variable genes")
  }

  if (method == "custom") {
    if (is.null(genes)) {
      stop("genes must be provided when method = 'custom'")
    }

    # Validate custom genes
    annotated_genes <- object@gene_annotation$gene
    valid_genes <- intersect(genes, annotated_genes)

    if (length(valid_genes) == 0) {
      stop("No provided genes found in gene annotation")
    }

    if (length(valid_genes) < length(genes)) {
      warning(length(genes) - length(valid_genes),
              " genes not found in annotation")
    }

    object@variable_genes <- valid_genes

  } else if (method == "markers") {
    message("Finding cell type marker genes...")

    # Create temporary Seurat object
    temp_seurat <- CreateSeuratObject(counts = object@RNA)
    temp_seurat@meta.data$cell_type <- object@cell_metadata$cell_type
    Idents(temp_seurat) <- object@cell_metadata$cell_type

    # Normalize
    temp_seurat <- NormalizeData(temp_seurat, verbose = FALSE)

    # Find markers
    markers <- FindAllMarkers(temp_seurat,
                              only.pos = TRUE,
                              min.pct = min.pct,
                              logfc.threshold = logfc.threshold,
                              verbose = FALSE)

    # Get unique marker genes
    marker_genes <- unique(markers$gene)

    # Filter to annotated genes
    annotated_genes <- object@gene_annotation$gene
    valid_genes <- intersect(marker_genes, annotated_genes)

    if (length(valid_genes) == 0) {
      stop("No marker genes found in gene annotation")
    }

    object@variable_genes <- valid_genes
  }

  object@parameters$variable_gene_method <- method
  object@parameters$variable_gene_min_pct <- min.pct
  object@parameters$variable_gene_logfc <- logfc.threshold

  message("Variable genes selected: ", length(object@variable_genes))

  return(object)
}
