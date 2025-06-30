# Unit tests for BayesCNet package

test_that("BayesCNet object creation works", {
  # Create small test data
  set.seed(123)
  n_genes <- 100
  n_peaks <- 200
  n_cells <- 50

  rna_mat <- matrix(rpois(n_genes * n_cells, 5),
                    nrow = n_genes, ncol = n_cells)
  rownames(rna_mat) <- paste0("Gene", 1:n_genes)
  colnames(rna_mat) <- paste0("Cell", 1:n_cells)

  atac_mat <- matrix(rbinom(n_peaks * n_cells, 1, 0.3),
                     nrow = n_peaks, ncol = n_cells)
  rownames(atac_mat) <- paste0("chr1_", 1:n_peaks, "_", 1:n_peaks + 1000)
  colnames(atac_mat) <- paste0("Cell", 1:n_cells)

  cell_meta <- data.frame(
    cell_type = rep(c("TypeA", "TypeB"), each = 25),
    row.names = paste0("Cell", 1:n_cells)
  )

  # Test object creation
  bcnet <- CreateBayesCNetFromMatrices(
    rna_matrix = rna_mat,
    atac_matrix = atac_mat,
    cell_metadata = cell_meta
  )

  expect_s4_class(bcnet, "BayesCNet")
  expect_equal(nrow(GetRNA(bcnet)), n_genes)
  expect_equal(nrow(GetATAC(bcnet)), n_peaks)
  expect_equal(nrow(GetCellMetadata(bcnet)), n_cells)
})

test_that("Cell type hierarchy works", {
  # Create test object
  set.seed(123)
  rna_mat <- matrix(rpois(100 * 60, 5), nrow = 100, ncol = 60)
  atac_mat <- matrix(rbinom(200 * 60, 1, 0.3), nrow = 200, ncol = 60)

  cell_meta <- data.frame(
    cell_type = rep(c("Root", "TypeA", "TypeB"), each = 20),
    row.names = paste0("Cell", 1:60)
  )

  bcnet <- CreateBayesCNetFromMatrices(rna_mat, atac_mat, cell_meta)

  # Add hierarchy
  edges <- c("Root", "TypeA", "Root", "TypeB")
  bcnet <- AddCellTypeHierarchy(bcnet, edges, lambda = 1.0)

  expect_true(length(bcnet@cell_type_hierarchy) > 0)
  expect_equal(bcnet@cell_type_hierarchy$lambda, 1.0)
  expect_s4_class(bcnet@cell_type_hierarchy$graph, "igraph")
})

test_that("Gene annotation validation works", {
  # Create test object
  set.seed(123)
  rna_mat <- matrix(rpois(100 * 50, 5), nrow = 100, ncol = 50)
  rownames(rna_mat) <- paste0("Gene", 1:100)
  atac_mat <- matrix(rbinom(200 * 50, 1, 0.3), nrow = 200, ncol = 50)

  cell_meta <- data.frame(
    cell_type = rep(c("TypeA", "TypeB"), each = 25)
  )

  bcnet <- CreateBayesCNetFromMatrices(rna_mat, atac_mat, cell_meta)

  # Test valid annotation
  gene_annot <- data.frame(
    gene = paste0("Gene", 1:50),
    chr = rep("chr1", 50),
    tss = seq(1000, 50000, 1000)
  )

  bcnet <- AddGeneAnnotation(bcnet, gene_annot)
  expect_equal(nrow(bcnet@gene_annotation), 50)

  # Test invalid annotation (missing columns)
  bad_annot <- data.frame(
    gene = paste0("Gene", 1:50),
    chromosome = rep("chr1", 50)  # Wrong column name
  )

  expect_error(AddGeneAnnotation(bcnet, bad_annot),
               "missing required columns")
})

test_that("Aggregation parameters are validated", {
  # Create test object
  set.seed(123)
  n_cells <- 30  # Small number of cells

  rna_mat <- matrix(rpois(100 * n_cells, 5), nrow = 100, ncol = n_cells)
  atac_mat <- matrix(rbinom(200 * n_cells, 1, 0.3), nrow = 200, ncol = n_cells)

  cell_meta <- data.frame(
    cell_type = rep("TypeA", n_cells)
  )

  bcnet <- CreateBayesCNetFromMatrices(rna_mat, atac_mat, cell_meta)

  # Test k too large
  expect_error(AggregateByKNN(bcnet, k = 50),
               "Cell types with fewer than k cells")

  # Test invalid max_overlap
  expect_error(AggregateByKNN(bcnet, k = 5, max_overlap = 1.5),
               "max_overlap must be between 0 and 1")
})

test_that("IsReadyForInference works correctly", {
  # Create minimal object
  set.seed(123)
  rna_mat <- matrix(rpois(100 * 50, 5), nrow = 100, ncol = 50)
  atac_mat <- matrix(rbinom(200 * 50, 1, 0.3), nrow = 200, ncol = 50)

  cell_meta <- data.frame(
    cell_type = rep(c("TypeA", "TypeB"), each = 25)
  )

  bcnet <- CreateBayesCNetFromMatrices(rna_mat, atac_mat, cell_meta)

  # Should not be ready initially
  expect_false(IsReadyForInference(bcnet))

  # Add required components one by one
  # Note: In practice, you would add real data here
  # This is just to test the validation logic
})
