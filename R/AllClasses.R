#' The BayesCNet class
#'
#' @slot RNA RNA count matrix (genes x cells)
#' @slot ATAC ATAC count matrix (peaks x cells)
#' @slot cell_metadata DataFrame with cell metadata
#' @slot gene_annotation DataFrame with gene TSS information
#' @slot cell_type_hierarchy List containing hierarchy information
#' @slot aggregated_data List containing aggregated cell data
#' @slot variable_genes Character vector of selected genes
#' @slot cell_embeddings Matrix of cell embeddings (cells x 2)
#' @slot embedding_type Character string indicating embedding type
#' @slot inference_results DataFrame with inference results
#' @slot parameters List of parameters used
#' @slot project_name Character string for project name
#'
#' @exportClass BayesCNet
setClass(
  "BayesCNet",
  slots = c(
    RNA = "matrix",
    ATAC = "matrix",
    cell_metadata = "data.frame",
    gene_annotation = "data.frame",
    cell_type_hierarchy = "list",
    aggregated_data = "list",
    variable_genes = "character",
    cell_embeddings = "matrix",
    embedding_type = "character",
    inference_results = "data.frame",
    parameters = "list",
    project_name = "character"
  ),
  prototype = list(
    RNA = matrix(),
    ATAC = matrix(),
    cell_metadata = data.frame(),
    gene_annotation = data.frame(),
    cell_type_hierarchy = list(),
    aggregated_data = list(),
    variable_genes = character(),
    cell_embeddings = matrix(),
    embedding_type = character(),
    inference_results = data.frame(),
    parameters = list(),
    project_name = "BayesCNet_project"
  ),
  validity = function(object) {
    errors <- character()

    # Check if RNA and ATAC have same number of cells
    if (ncol(object@RNA) > 0 && ncol(object@ATAC) > 0) {
      if (ncol(object@RNA) != ncol(object@ATAC)) {
        errors <- c(errors, "RNA and ATAC must have same number of cells")
      }
    }

    # Check cell metadata matches
    if (nrow(object@cell_metadata) > 0 && ncol(object@RNA) > 0) {
      if (nrow(object@cell_metadata) != ncol(object@RNA)) {
        errors <- c(errors, "Cell metadata must have one row per cell")
      }
    }

    # Check required columns in cell metadata
    if (nrow(object@cell_metadata) > 0) {
      if (!"cell_type" %in% colnames(object@cell_metadata)) {
        errors <- c(errors, "Cell metadata must contain 'cell_type' column")
      }
    }

    # Check gene annotation format if provided
    if (nrow(object@gene_annotation) > 0) {
      required_cols <- c("gene", "chr", "tss")
      missing_cols <- setdiff(required_cols, colnames(object@gene_annotation))
      if (length(missing_cols) > 0) {
        errors <- c(errors, paste("Gene annotation missing columns:",
                                  paste(missing_cols, collapse = ", ")))
      }
    }

    if (length(errors) == 0) TRUE else errors
  }
)

#' Show method for BayesCNet
#' @param object A BayesCNet object
setMethod("show", "BayesCNet", function(object) {
  cat("An object of class BayesCNet\n")
  cat("Project:", object@project_name, "\n\n")

  # Data dimensions
  cat("Data:\n")
  if (ncol(object@RNA) > 0) {
    cat("  RNA:", nrow(object@RNA), "genes x", ncol(object@RNA), "cells\n")
  }
  if (ncol(object@ATAC) > 0) {
    cat("  ATAC:", nrow(object@ATAC), "peaks x", ncol(object@ATAC), "cells\n")
  }

  # Cell types
  if (nrow(object@cell_metadata) > 0 && "cell_type" %in% colnames(object@cell_metadata)) {
    cell_types <- unique(object@cell_metadata$cell_type)
    cat("  Cell types:", length(cell_types), "\n")
  }

  # Cell embeddings
  if (nrow(object@cell_embeddings) > 0) {
    cat("  Cell embeddings:", object@embedding_type,
        "(", nrow(object@cell_embeddings), "cells x",
        ncol(object@cell_embeddings), "dims )\n")
  }

  # Status of preparation steps
  cat("\nPreparation status:\n")

  # Hierarchy
  if (length(object@cell_type_hierarchy) > 0) {
    cat("  ✓ Cell type hierarchy defined (lambda =",
        object@cell_type_hierarchy$lambda, ")\n")
  } else {
    cat("  ✗ Cell type hierarchy not defined\n")
  }

  # Gene annotation
  if (nrow(object@gene_annotation) > 0) {
    cat("  ✓ Gene annotation added (", nrow(object@gene_annotation), " genes)\n")
  } else {
    cat("  ✗ Gene annotation not added\n")
  }

  # Aggregation
  if (length(object@aggregated_data) > 0) {
    cat("  ✓ Data aggregated (k =", object@aggregated_data$k_neighbors, ")\n")
  } else {
    cat("  ✗ Data not aggregated\n")
  }

  # Variable genes
  if (length(object@variable_genes) > 0) {
    cat("  ✓ Variable genes selected (", length(object@variable_genes), " genes)\n")
  } else {
    cat("  ✗ Variable genes not selected\n")
  }

  # Inference
  if (nrow(object@inference_results) > 0) {
    cat("\nInference:\n")
    cat("  ✓ Complete -", nrow(object@inference_results), "connections found\n")
    n_celltypes <- length(unique(object@inference_results$CellType))
    n_genes <- length(unique(object@inference_results$Gene))
    cat("  ", n_celltypes, "cell types,", n_genes, "genes analyzed\n")
  } else {
    cat("\nInference: Not run\n")
  }

  # Check if ready for inference
  if (length(object@cell_type_hierarchy) == 0 ||
      nrow(object@gene_annotation) == 0 ||
      length(object@aggregated_data) == 0 ||
      length(object@variable_genes) == 0) {
    cat("\nNote: Complete data preparation before running inference\n")
  }
})
