# BayesCNet
BayesCNet: Hierarchical Bayesian inference of cell type-specific gene regulatory networks from single-cell multi-omics data. Leverages cell type relationships to improve network accuracy, especially for rare cell populations.
# BayesCNet

**B**ayesian **C**ell-type-specific **Net**work Inference from Single-cell Multi-omics Data

## Overview

BayesCNet is an R package for inferring cell-type-specific gene regulatory networks from paired single-cell RNA-seq and ATAC-seq data. It uses a hierarchical Bayesian framework that leverages cell type relationships to improve network inference accuracy.

## Key Features

- **Cell-type-specific inference**: Infers distinct regulatory networks for each cell type
- **Hierarchical modeling**: Uses cell type relationships to share information across related cell types
- **Multi-omics integration**: Combines scRNA-seq and scATAC-seq data
- **Scalable**: Supports parallel processing for large datasets
- **User-friendly**: Two-stage workflow with clear data preparation and inference steps

## Installation

```r
# Install from GitHub (once available)
# devtools::install_github("yourusername/BayesCNet")

# For development
devtools::load_all("path/to/BayesCNet")
```

## Quick Start

```r
library(BayesCNet)

# Stage 1: Data Preparation
# Create BayesCNet object from Seurat
bcnet <- CreateBayesCNet(seurat_object)

# Add cell type hierarchy
edges <- c("HSPC", "T Cells", "T Cells", "CD4 T", "T Cells", "CD8 T")
bcnet <- AddCellTypeHierarchy(bcnet, edges, lambda = 1.0)

# Add gene annotation
bcnet <- AddGeneAnnotation(bcnet, gene_annotation_df)

# Aggregate cells
bcnet <- AggregateByKNN(bcnet, k = 50)

# Select variable genes
bcnet <- AddVariableGenes(bcnet, method = "markers")

# Stage 2: Run Inference
bcnet <- RunBayesCNet(bcnet, window = 250000, cores = 4)

# Access Results
all_networks <- GetAllNetworks(bcnet)
cd4_network <- GetCellTypeSpecificNetwork(bcnet, "CD4 T Cells")
```

## Workflow

### Stage 1: Data Preparation

1. **Create BayesCNet object**: Import data from Seurat or raw matrices
2. **Define cell type hierarchy**: Specify relationships between cell types
3. **Add gene annotation**: Provide TSS information for regulatory window definition
4. **Aggregate cells**: Create metacells using KNN approach
5. **Select genes**: Choose cell-type marker genes or custom gene set

### Stage 2: Inference

1. **Run BayesCNet**: Execute the Bayesian inference model
2. **Access results**: Extract cell-type-specific networks

## Input Requirements

- **RNA data**: Gene expression counts (genes × cells)
- **ATAC data**: Peak accessibility counts (peaks × cells)
- **Cell metadata**: Must include cell type annotations
- **Gene annotation**: TSS coordinates for genes of interest

## Output

The main output is a data frame containing:
- `CellType`: Cell type name
- `Gene`: Target gene
- `Peak1`: Promoter peak
- `Peak2`: Enhancer peak
- `Pmean`: Posterior mean of regulatory coefficient
- `Pvar`: Posterior variance
- `Importance`: Importance score (|Pmean|/sqrt(Pvar))

## Method Details

BayesCNet implements a hierarchical Bayesian model that:
1. Models gene expression as a function of chromatin accessibility
2. Shares information across cell types based on their hierarchical relationships
3. Infers cell-type-specific regulatory coefficients
4. Provides uncertainty estimates for each connection

## Citation

If you use BayesCNet in your research, please cite:
[Citation information to be added]

## License

GPL-3

## Contact

[Your contact information]
