#' Run BayesCNet inference
#'
#' @param object A prepared BayesCNet object
#' @param window Regulatory window size around TSS (default: 250000)
#' @param regularization Regularization parameter (default: 1e-6)
#' @param cores Number of cores for parallel processing (default: 1)
#' @param verbose Print progress messages (default: TRUE)
#'
#' @return BayesCNet object with inference results
#' @export
#'
#' @examples
#' \dontrun{
#' bcnet <- RunBayesCNet(bcnet, window = 250000, cores = 4)
#' }
RunBayesCNet <- function(object, window = 250000, regularization = 1e-6,
                         cores = 1, verbose = TRUE) {

  # Validate object is ready
  if (!inherits(object, "BayesCNet")) {
    stop("object must be a BayesCNet object")
  }

  # Check all required components
  if (length(object@cell_type_hierarchy) == 0) {
    stop("Cell type hierarchy not defined. Run AddCellTypeHierarchy() first.")
  }

  if (nrow(object@gene_annotation) == 0) {
    stop("Gene annotation not added. Run AddGeneAnnotation() first.")
  }

  if (length(object@aggregated_data) == 0) {
    stop("Data not aggregated. Run AggregateByKNN() first.")
  }

  if (length(object@variable_genes) == 0) {
    stop("Variable genes not selected. Run AddVariableGenes() first.")
  }

  # Validate parameters
  if (window <= 0) {
    stop("window must be positive")
  }

  if (regularization <= 0) {
    stop("regularization must be positive")
  }

  if (cores < 1) {
    stop("cores must be at least 1")
  }

  message("Starting BayesCNet inference...")
  message("  Genes to analyze: ", length(object@variable_genes))
  message("  Regulatory window: ", format(window, big.mark = ","), " bp")
  message("  Using ", cores, " core(s)")

  # Extract required data
  data_rna <- object@aggregated_data$rna
  data_atac <- object@aggregated_data$atac
  cell_assignments <- object@aggregated_data$cell_assignments

  # Get cell type vector from assignments
  cell_types_vector <- cell_assignments[, ncol(cell_assignments)]
  celltypes_order <- unique(cell_types_vector)

  # Get Sigma2 (prior covariance)
  Sigma2_est <- object@cell_type_hierarchy$sigma2

  # Prepare gene information
  gene_annotation <- object@gene_annotation
  variable_genes <- object@variable_genes
  gene_info <- gene_annotation[gene_annotation$gene %in% variable_genes, ]

  # Store parameters
  object@parameters$inference_window <- window
  object@parameters$inference_regularization <- regularization
  object@parameters$inference_cores <- cores

  # Run inference
  if (cores > 1) {
    # Parallel processing
    if (verbose) message("Running parallel inference...")

    # Set up parallel backend
    cl <- parallel::makeCluster(cores)
    on.exit(parallel::stopCluster(cl))

    # Export required objects and functions
    parallel::clusterExport(cl, c("data_rna", "data_atac", "cell_types_vector",
                                  "Sigma2_est", "gene_info", "window",
                                  "regularization", "process_gene_inference",
                                  "find_overlapping_coordinates"),
                            envir = environment())

    # Load required packages on workers
    parallel::clusterEvalQ(cl, {
      library(Matrix)
      library(expm)
    })

    # Run parallel inference
    results_list <- parallel::parLapply(cl, 1:nrow(gene_info), function(i) {
      process_gene_inference(
        gene_idx = i,
        gene_info = gene_info,
        data_rna = data_rna,
        data_atac = data_atac,
        cell_types_vector = cell_types_vector,
        Sigma2_est = Sigma2_est,
        window = window,
        regularization = regularization
      )
    })

  } else {
    # Sequential processing with progress bar
    if (verbose) {
      pb <- txtProgressBar(min = 0, max = nrow(gene_info), style = 3)
    }

    results_list <- list()
    for (i in 1:nrow(gene_info)) {
      if (verbose) setTxtProgressBar(pb, i)

      results_list[[i]] <- process_gene_inference(
        gene_idx = i,
        gene_info = gene_info,
        data_rna = data_rna,
        data_atac = data_atac,
        cell_types_vector = cell_types_vector,
        Sigma2_est = Sigma2_est,
        window = window,
        regularization = regularization
      )
    }

    if (verbose) close(pb)
  }

  # Combine results
  results_list <- results_list[!sapply(results_list, is.null)]

  if (length(results_list) == 0) {
    stop("No valid results obtained from inference")
  }

  inference_results <- do.call(rbind, results_list)
  rownames(inference_results) <- NULL

  # Store results
  object@inference_results <- inference_results

  # Summary statistics
  n_connections <- nrow(inference_results)
  n_genes <- length(unique(inference_results$Gene))
  n_celltypes <- length(unique(inference_results$CellType))

  message("\nInference complete!")
  message("  Total connections: ", format(n_connections, big.mark = ","))
  message("  Genes analyzed: ", n_genes)
  message("  Cell types: ", n_celltypes)

  return(object)
}

#' Get cell-type-specific network
#'
#' @param object A BayesCNet object with completed inference
#' @param celltype Cell type name
#' @param min_importance Minimum importance score threshold (default: NULL)
#' @param top_n Return top N connections (default: NULL, return all)
#'
#' @return Data frame with network connections for specified cell type
#' @export
GetCellTypeSpecificNetwork <- function(object, celltype,
                                       min_importance = NULL,
                                       top_n = NULL) {

  # Validate inputs
  if (!inherits(object, "BayesCNet")) {
    stop("object must be a BayesCNet object")
  }

  if (nrow(object@inference_results) == 0) {
    stop("No inference results found. Run RunBayesCNet() first.")
  }

  # Check if celltype exists
  available_celltypes <- unique(object@inference_results$CellType)
  if (!celltype %in% available_celltypes) {
    stop("Cell type '", celltype, "' not found. Available cell types: ",
         paste(available_celltypes, collapse = ", "))
  }

  # Filter results
  results <- object@inference_results[object@inference_results$CellType == celltype, ]

  # Apply importance threshold if specified
  if (!is.null(min_importance)) {
    results <- results[results$Importance >= min_importance, ]
  }

  # Sort by importance
  results <- results[order(results$Importance, decreasing = TRUE), ]

  # Return top N if specified
  if (!is.null(top_n) && top_n < nrow(results)) {
    results <- results[1:top_n, ]
  }

  # Remove cell type column since it's redundant
  results$CellType <- NULL

  return(results)
}

#' Get all networks from BayesCNet inference
#'
#' @param object A BayesCNet object with completed inference
#' @param min_importance Minimum importance score threshold (default: NULL)
#'
#' @return Data frame with all network connections
#' @export
GetAllNetworks <- function(object, min_importance = NULL) {

  # Validate inputs
  if (!inherits(object, "BayesCNet")) {
    stop("object must be a BayesCNet object")
  }

  if (nrow(object@inference_results) == 0) {
    stop("No inference results found. Run RunBayesCNet() first.")
  }

  # Get results
  results <- object@inference_results

  # Apply importance threshold if specified
  if (!is.null(min_importance)) {
    results <- results[results$Importance >= min_importance, ]
  }

  # Sort by cell type and importance
  results <- results[order(results$CellType, -results$Importance), ]

  return(results)
}
