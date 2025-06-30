# BayesCNet Example Usage
# This script demonstrates the two-stage workflow

library(BayesCNet)
library(Seurat)

# Assuming you have a Seurat object with RNA and ATAC data
# load("pbmc_seurat.rda")

#############################################
# Stage 1: Data Preparation
#############################################

# 1. Create BayesCNet object from Seurat
bcnet <- CreateBayesCNet(
  seurat_object = pbmc,
  rna_assay = "RNA",
  atac_assay = "ATAC",
  cell_type_col = "celltype"
)

# Check the object
show(bcnet)

# 2. Add cell type hierarchy
# Define edges for PBMC cell type tree
edges <- c(
  "HSPC", "Monocytes",
  "HSPC", "T Cells",
  "HSPC", "B Cells",
  "HSPC", "NK",
  "HSPC", "Dendritic Cells",
  "Monocytes", "CD14 Mono",
  "Monocytes", "CD16 Mono",
  "T Cells", "CD4 T Cells",
  "T Cells", "CD8 T Cells",
  "CD4 T Cells", "CD4 Naive",
  "CD4 T Cells", "CD4 TCM",
  "CD8 T Cells", "CD8 Naive",
  "CD8 T Cells", "CD8 TEM_1",
  "B Cells", "Naive B",
  "B Cells", "Memory B"
)

bcnet <- AddCellTypeHierarchy(
  object = bcnet,
  edges = edges,
  lambda = 1.0  # Controls the strength of hierarchy influence
)

# 3. Add gene annotation
# You need a data frame with columns: gene, chr, tss
# For example:
gene_annotation <- data.frame(
  gene = c("IL2", "CD4", "CD8A", "MS4A1", "GNLY"),
  chr = c("chr4", "chr12", "chr2", "chr11", "chr2"),
  tss = c(123456789, 6934951, 87011728, 60223282, 85932145)
)
# In practice, load comprehensive annotation from biomaRt or similar

bcnet <- AddGeneAnnotation(
  object = bcnet,
  annotation = gene_annotation
)

# 4. Aggregate cells into metacells
bcnet <- AggregateByKNN(
  object = bcnet,
  k = 50,           # Number of neighbors
  max_overlap = 0.8, # Maximum overlap between metacells
  seed = 123,
  normalize = TRUE   # CPM normalization
)

# 5. Select variable genes
# Method 1: Find cell type markers automatically
bcnet <- AddVariableGenes(
  object = bcnet,
  method = "markers",
  min.pct = 0.1,
  logfc.threshold = 0.1
)

# Method 2: Use custom gene list
# bcnet <- AddVariableGenes(
#   object = bcnet,
#   method = "custom",
#   genes = c("IL2", "CD4", "CD8A", "MS4A1")
# )

# Check if ready for inference
IsReadyForInference(bcnet)

# View summary
summary(bcnet)

#############################################
# Stage 2: Run Inference
#############################################

# Run BayesCNet inference
bcnet <- RunBayesCNet(
  object = bcnet,
  window = 250000,      # 250kb regulatory window
  regularization = 1e-6,
  cores = 4,            # Use 4 cores for parallel processing
  verbose = TRUE
)

#############################################
# Access Results
#############################################

# Get all networks
all_networks <- GetAllNetworks(bcnet)
head(all_networks)

# Get cell-type-specific network
cd4_network <- GetCellTypeSpecificNetwork(
  object = bcnet,
  celltype = "CD4 T Cells",
  min_importance = 2.0,  # Filter by importance score
  top_n = 100           # Get top 100 connections
)

# View top connections
head(cd4_network[order(cd4_network$Importance, decreasing = TRUE), ])

# Save results
write.csv(all_networks, "bayescnet_all_networks.csv", row.names = FALSE)
write.csv(cd4_network, "bayescnet_cd4_network.csv", row.names = FALSE)

# Access raw data if needed
rna_data <- GetRNA(bcnet)
atac_data <- GetATAC(bcnet)
cell_meta <- GetCellMetadata(bcnet)
params <- GetParameters(bcnet)

#############################################
# Filtering and Analysis
#############################################

# Get high-confidence connections (high importance, low variance)
high_conf <- all_networks[
  all_networks$Importance > 3 &
    all_networks$Pvar < 0.1,
]

# Cell type comparison
library(dplyr)
network_stats <- all_networks %>%
  group_by(CellType) %>%
  summarise(
    n_connections = n(),
    mean_importance = mean(Importance),
    n_genes = n_distinct(Gene)
  )

print(network_stats)

# Gene-centric view
il2_connections <- all_networks[all_networks$Gene == "IL2", ]
table(il2_connections$CellType)
