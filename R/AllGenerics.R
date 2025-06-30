# This file contains all generic function definitions for BayesCNet

#' Basic accessors for BayesCNet objects
#' @param x BayesCNet object
#' @rdname accessors
#' @export
setGeneric("GetRNA", function(x) standardGeneric("GetRNA"))

#' @rdname accessors
#' @export
setGeneric("GetATAC", function(x) standardGeneric("GetATAC"))

#' @rdname accessors
#' @export
setGeneric("GetCellMetadata", function(x) standardGeneric("GetCellMetadata"))

#' @rdname accessors
#' @export
setGeneric("GetGeneAnnotation", function(x) standardGeneric("GetGeneAnnotation"))

#' @rdname accessors
#' @export
setGeneric("GetVariableGenes", function(x) standardGeneric("GetVariableGenes"))

#' @rdname accessors
#' @export
setGeneric("GetParameters", function(x) standardGeneric("GetParameters"))

#' @rdname accessors
#' @export
setGeneric("GetCellEmbeddings", function(x) standardGeneric("GetCellEmbeddings"))
