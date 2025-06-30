#' Accessor methods for BayesCNet objects
#'
#' Methods to access data from BayesCNet objects
#'
#' @param x A BayesCNet object
#' @return The requested data

#' @rdname accessors
#' @export
setMethod("GetRNA", "BayesCNet", function(x) {
  x@RNA
})

#' @rdname accessors
#' @export
setMethod("GetATAC", "BayesCNet", function(x) {
  x@ATAC
})

#' @rdname accessors
#' @export
setMethod("GetCellMetadata", "BayesCNet", function(x) {
  x@cell_metadata
})

#' @rdname accessors
#' @export
setMethod("GetGeneAnnotation", "BayesCNet", function(x) {
  x@gene_annotation
})

#' @rdname accessors
#' @export
setMethod("GetVariableGenes", "BayesCNet", function(x) {
  x@variable_genes
})

#' @rdname accessors
#' @export
setMethod("GetParameters", "BayesCNet", function(x) {
  x@parameters
})

#' @rdname accessors
#' @export
setMethod("GetCellEmbeddings", "BayesCNet", function(x) {
  x@cell_embeddings
})

#' Summary method for BayesCNet
#' @param object A BayesCNet object
#' @export
setMethod("summary", "BayesCNet", function(object) {
  cat("BayesCNet Summary\n")
  cat("================\n\n")

  # Basic info
  cat("Project:", object@project_name, "\n")
  cat("RNA data:", nrow(object@RNA), "genes x", ncol(object@RNA), "cells\n")
  cat("ATAC data:", nrow(object@ATAC), "peaks x", ncol(object@ATAC), "cells\n")

  # Cell types
  if ("cell_type" %in% colnames(object@cell_metadata)) {
    ct_table <- table(object@cell_metadata$cell_type)
    cat("\nCell type distribution:\n")
    print(ct_table)
  }

  # Analysis status
  cat("\nAnalysis status:\n")
  cat("- Cell type hierarchy:",
      ifelse(length(object@cell_type_hierarchy) > 0, "Yes", "No"), "\n")
  cat("- Gene annotation:",
      ifelse(nrow(object@gene_annotation) > 0,
             paste0("Yes (", nrow(object@gene_annotation), " genes)"), "No"), "\n")
  cat("- Data aggregated:",
      ifelse(length(object@aggregated_data) > 0,
             paste0("Yes (k=", object@aggregated_data$k_neighbors, ")"), "No"), "\n")
  cat("- Variable genes:",
      ifelse(length(object@variable_genes) > 0,
             paste0("Yes (", length(object@variable_genes), " genes)"), "No"), "\n")
  cat("- Inference complete:",
      ifelse(nrow(object@inference_results) > 0,
             paste0("Yes (", nrow(object@inference_results), " connections)"), "No"), "\n")

  # Parameters if inference is complete
  if (nrow(object@inference_results) > 0) {
    cat("\nInference parameters:\n")
    cat("- Window size:", object@parameters$inference_window, "bp\n")
    cat("- Regularization:", object@parameters$inference_regularization, "\n")
    cat("- Cores used:", object@parameters$inference_cores, "\n")
  }

  invisible(object)
})

#' Check if BayesCNet object is ready for inference
#' @param object A BayesCNet object
#' @return Logical indicating if object is ready
#' @export
IsReadyForInference <- function(object) {
  if (!inherits(object, "BayesCNet")) {
    stop("object must be a BayesCNet object")
  }

  ready <- TRUE
  messages <- character()

  if (length(object@cell_type_hierarchy) == 0) {
    ready <- FALSE
    messages <- c(messages, "Cell type hierarchy not defined")
  }

  if (nrow(object@gene_annotation) == 0) {
    ready <- FALSE
    messages <- c(messages, "Gene annotation not added")
  }

  if (length(object@aggregated_data) == 0) {
    ready <- FALSE
    messages <- c(messages, "Data not aggregated")
  }

  if (length(object@variable_genes) == 0) {
    ready <- FALSE
    messages <- c(messages, "Variable genes not selected")
  }

  if (!ready) {
    message("Object not ready for inference. Missing:")
    for (msg in messages) {
      message("  - ", msg)
    }
  } else {
    message("Object is ready for inference!")
  }

  return(ready)
}
