#' Create a BayesCNet object from Seurat
#'
#' @param seurat_object A Seurat object with RNA and ATAC assays
#' @param rna_assay Name of RNA assay in Seurat object
#' @param atac_assay Name of ATAC assay in Seurat object
#' @param cell_type_col Column name in metadata for cell types
#' @param project_name Character string for project name
#'
#' @return A BayesCNet object
#' @export
#'
#' @examples
#' \dontrun{
#' bcnet <- CreateBayesCNet(pbmc_seurat)
#' }
CreateBayesCNet <- function(
    seurat_object,
    rna_assay = "RNA",
    atac_assay = "ATAC",
    cell_type_col = "celltype",
    project_name = NULL
) {

  # Validate Seurat object
  if (!inherits(seurat_object, "Seurat")) {
    stop("Input must be a Seurat object")
  }

  # Check required assays exist
  if (!rna_assay %in% names(seurat_object@assays)) {
    stop("RNA assay '", rna_assay, "' not found in Seurat object")
  }
  if (!atac_assay %in% names(seurat_object@assays)) {
    stop("ATAC assay '", atac_assay, "' not found in Seurat object")
  }

  # Extract RNA data
  message("Extracting RNA data...")
  rna_counts <- GetAssayData(seurat_object, assay = rna_assay, slot = "counts")
  rna_counts <- as.matrix(rna_counts)

  # Extract ATAC data
  message("Extracting ATAC data...")
  atac_counts <- GetAssayData(seurat_object, assay = atac_assay, slot = "counts")
  atac_counts <- as.matrix(atac_counts)

  # Validate dimensions
  if (ncol(rna_counts) != ncol(atac_counts)) {
    stop("RNA and ATAC assays must have the same number of cells")
  }

  # Extract cell metadata
  cell_metadata <- seurat_object@meta.data

  # Validate cell type column
  if (!cell_type_col %in% colnames(cell_metadata)) {
    # Try to use Idents if cell type column not found
    warning("Cell type column '", cell_type_col, "' not found. Using Idents()")
    cell_metadata$cell_type <- as.character(Idents(seurat_object))
  } else {
    cell_metadata$cell_type <- as.character(cell_metadata[[cell_type_col]])
  }

  # Check for NA cell types
  if (any(is.na(cell_metadata$cell_type))) {
    stop("Cell type information contains NA values")
  }

  # Set project name
  if (is.null(project_name)) {
    project_name <- seurat_object@project.name
    if (is.null(project_name) || project_name == "") {
      project_name <- "BayesCNet_project"
    }
  }

  # Extract cell embeddings
  message("Extracting cell embeddings...")

  # Check for available reductions
  available_reductions <- names(seurat_object@reductions)

  # Determine which reduction to use
  embedding_type <- NULL
  cell_embeddings <- matrix()

  # Priority: wnn.umap > umap
  if ("wnn.umap" %in% available_reductions) {
    cell_embeddings <- seurat_object@reductions$wnn.umap@cell.embeddings
    embedding_type <- "wnn.umap"
    message("  Using wnn.umap embeddings")
  } else if ("umap" %in% available_reductions) {
    cell_embeddings <- seurat_object@reductions$umap@cell.embeddings
    embedding_type <- "umap"
    message("  Using umap embeddings")
  } else {
    stop("No UMAP embeddings found in Seurat object. ",
         "Please run RunUMAP() or FindMultiModalNeighbors() + IntegrateEmbeddings() first.")
  }

  # Ensure embeddings are 2D
  if (ncol(cell_embeddings) != 2) {
    stop("UMAP embeddings must be 2-dimensional. Found ", ncol(cell_embeddings), " dimensions.")
  }

  # Create BayesCNet object
  object <- new("BayesCNet",
                RNA = rna_counts,
                ATAC = atac_counts,
                cell_metadata = cell_metadata,
                cell_embeddings = cell_embeddings,
                embedding_type = embedding_type,
                project_name = project_name)

  # Store initial parameters
  object@parameters$rna_assay <- rna_assay
  object@parameters$atac_assay <- atac_assay
  object@parameters$cell_type_col <- cell_type_col

  message("Created BayesCNet object")
  message("  RNA: ", nrow(rna_counts), " genes x ", ncol(rna_counts), " cells")
  message("  ATAC: ", nrow(atac_counts), " peaks x ", ncol(atac_counts), " cells")
  message("  Cell types: ", paste(unique(cell_metadata$cell_type), collapse = ", "))

  return(object)
}

#' Create BayesCNet from matrices
#'
#' @param rna_matrix RNA count matrix (genes x cells)
#' @param atac_matrix ATAC count matrix (peaks x cells)
#' @param cell_metadata Data frame with cell metadata (must include cell_type column)
#' @param project_name Character string for project name
#'
#' @return A BayesCNet object
#' @export
CreateBayesCNetFromMatrices <- function(
    rna_matrix,
    atac_matrix,
    cell_metadata,
    project_name = "BayesCNet_project"
) {

  # Validate inputs
  if (!is.matrix(rna_matrix) && !inherits(rna_matrix, "dgCMatrix")) {
    stop("RNA must be a matrix or sparse matrix")
  }
  if (!is.matrix(atac_matrix) && !inherits(atac_matrix, "dgCMatrix")) {
    stop("ATAC must be a matrix or sparse matrix")
  }

  # Convert to regular matrices
  rna_matrix <- as.matrix(rna_matrix)
  atac_matrix <- as.matrix(atac_matrix)

  # Validate dimensions
  if (ncol(rna_matrix) != ncol(atac_matrix)) {
    stop("RNA and ATAC must have the same number of cells")
  }

  if (nrow(cell_metadata) != ncol(rna_matrix)) {
    stop("Cell metadata must have one row per cell")
  }

  # Check for cell_type column
  if (!"cell_type" %in% colnames(cell_metadata)) {
    stop("Cell metadata must contain 'cell_type' column")
  }

  # Create object
  object <- new("BayesCNet",
                RNA = rna_matrix,
                ATAC = atac_matrix,
                cell_metadata = cell_metadata,
                project_name = project_name)

  # Compute cell embeddings
  message("Computing cell embeddings...")

  # Create temporary Seurat object
  temp_seurat <- CreateSeuratObject(counts = rna_matrix, project = project_name)
  temp_seurat[["ATAC"]] <- CreateAssayObject(counts = atac_matrix)
  temp_seurat@meta.data$cell_type <- cell_metadata$cell_type
  Idents(temp_seurat) <- cell_metadata$cell_type

  # Check modalities
  has_rna <- nrow(rna_matrix) > 0 && ncol(rna_matrix) > 0
  has_atac <- nrow(atac_matrix) > 0 && ncol(atac_matrix) > 0

  if (has_rna && has_atac) {
    message("  Running WNN workflow for RNA+ATAC...")

    # Process RNA
    DefaultAssay(temp_seurat) <- "RNA"
    temp_seurat <- NormalizeData(temp_seurat, verbose = FALSE)
    temp_seurat <- FindVariableFeatures(temp_seurat, verbose = FALSE)
    temp_seurat <- ScaleData(temp_seurat, verbose = FALSE)
    temp_seurat <- RunPCA(temp_seurat, verbose = FALSE)

    # Process ATAC
    DefaultAssay(temp_seurat) <- "ATAC"
    temp_seurat <- RunTFIDF(temp_seurat)
    temp_seurat <- FindTopFeatures(temp_seurat, min.cutoff = 'q0')
    temp_seurat <- RunSVD(temp_seurat)

    # Run WNN
    temp_seurat <- FindMultiModalNeighbors(
      temp_seurat,
      reduction.list = list("pca", "lsi"),
      dims.list = list(1:30, 2:30),
      verbose = FALSE
    )

    temp_seurat <- RunUMAP(
      temp_seurat,
      nn.name = "weighted.nn",
      reduction.name = "wnn.umap",
      reduction.key = "wnnUMAP_",
      verbose = FALSE
    )

    object@cell_embeddings <- temp_seurat@reductions$wnn.umap@cell.embeddings
    object@embedding_type <- "wnn.umap"

  } else {
    stop("Cannot compute embeddings without RNA and ATAC data")
  }

  message("  Cell embeddings computed successfully")

  return(object)
}

