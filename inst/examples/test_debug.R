# test_debug.R - Manual debugging script
library(devtools)
load_all(".")  # Load the package

# Enable verbose output
options(BayesCNet.verbose = TRUE)


###### 1.5. Test Object Creation from seurat obj ######
cat("\n=== Testing Object Creation ===\n")

load("/Users/fengdi/Library/CloudStorage/OneDrive-UniversityofFlorida/GraphModel/Rcode/DIRECT_NET/pbmc_seurat.rda")
pbmc
head(pbmc@meta.data)
seurat_object = pbmc

selected_cells = c("CD16 Mono", "CD4 TEM", "CD8 TEM_1", "Naive B")
pbmc_subset <- subset(pbmc, subset = celltype %in% selected_cells)
pbmc_subset$celltype <- droplevels(pbmc_subset$celltype)
pbmc_subset@assays$SCT = NULL

save(pbmc_subset, file = "pbmc_demo.r")


bcnet_seurat <- CreateBayesCNet(pbmc_subset, cell_type_col = "celltype")

####### create bayescnet obj pass from seurat object.!!!!!!

###### 1. Test Object Creation from matrices ######
cat("\n=== Testing Object Creation ===\n")

pbmc_subset

# Extract RNA data
rna_counts <- GetAssayData(pbmc_subset, assay = rna_assay, slot = "counts")
rna_mat <- as.matrix(rna_counts)

# Extract ATAC data
atac_counts <- GetAssayData(pbmc_subset, assay = atac_assay, slot = "counts")
atac_mat <- as.matrix(atac_counts)

cell_meta <- pbmc_subset@meta.data
cell_meta$cell_type = cell_meta$celltype

# Test creation
bcnet <- CreateBayesCNetFromMatrices(rna_mat, atac_mat, cell_meta)
print(bcnet)

# Check slots
cat("\nRNA dimensions:", dim(bcnet@RNA), "\n")
cat("ATAC dimensions:", dim(bcnet@ATAC), "\n")
cat("Cell types:", unique(bcnet@cell_metadata$cell_type), "\n")

## pass create bayescnet class from raw matrices; paired RNA and ATAC !!!!!!

###### 2. Test Cell Type Hierarchy ######
cat("\n=== Testing Cell Type Hierarchy ===\n")

# Simple hierarchy
edges <- c("HSPC", "Monocytes",
           "HSPC", "T Cells",
           "HSPC", "B Cells",
           "HSPC", "NK",
           "HSPC", "Dendritic Cells",
           "Monocytes", "CD14 Mono",
           "Monocytes", "CD16 Mono",
           "T Cells", "CD4 T Cells",
           "T Cells", "CD8 T Cells",
           "T Cells", "gdT",
           "T Cells", "MAIT",
           "CD4 T Cells", "CD4 Naive",
           "CD4 T Cells", "CD4 TCM",
           "CD4 T Cells", "CD4 TEM",
           "CD4 T Cells", "Treg",
           "CD8 T Cells", "CD8 Naive",
           "CD8 T Cells", "CD8 TEM_1",
           "CD8 T Cells", "CD8 TEM_2",
           "B Cells", "Naive B",
           "B Cells", "Memory B",
           "B Cells", "Intermediate B",
           "B Cells", "Plasma",
           "Dendritic Cells", "cDC",
           "Dendritic Cells", "pDC")
bcnet <- AddCellTypeHierarchy(bcnet, edges, lambda = 1.0)
bcnet_seurat <- AddCellTypeHierarchy(bcnet_seurat, edges, lambda = 1.0)

# Check hierarchy
cat("Hierarchy added successfully\n")
cat("Lambda:", bcnet@cell_type_hierarchy$lambda, "\n")
cat("Sigma2 dimensions:", dim(bcnet@cell_type_hierarchy$sigma2), "\n")

### both obj created from seurat and matrice pass AddCellTypeHierarchy

###### 3. Test Gene Annotation ######
cat("\n=== Testing Gene Annotation ===\n")

library(EnsDb.Hsapiens.v86)
## create genome info (using only protein-coding genes)
# hg38 gene annotation- get TSS of all transcripts of all protein coding genes
edb <- EnsDb.Hsapiens.v86
listTables(edb)
listColumns(edb)
# Retrieve all transcripts for protein-coding genes
transcripts_info <- transcripts(edb,filter = TxBiotypeFilter("protein_coding"),columns=c("gene_name","tx_seq_start","tx_cds_seq_start"))
genome(transcripts_info) <- "hg38"
# TSS is the start position for '+' strand and end position for '-' strand
# Calculate TSS based on strand
TSS_starts <- ifelse(strand(transcripts_info) == "+", start(transcripts_info), end(transcripts_info))
TSS_starts <- IRanges(start = TSS_starts, width = 2)
# TSS dataframe
tss_df <- data.frame(
  chr = paste0("chr", seqnames(transcripts_info)),
  start = start(TSS_starts),
  end = end(TSS_starts),
  genes = mcols(transcripts_info)$gene_name
)
unik <- !duplicated(tss_df$genes) # filter out different transcript
tss_df = tss_df[unik,]
genome.info = tss_df
colnames(genome.info) <- c("Chrom","Starts","Ends","genes")
# keep only regular chromosome
# Define regular chromosomes
regular_chromosomes <- paste0("chr", c(1:22, "X", "Y"))
# Filter genome.info to keep only rows with regular chromosomes
genome.info <- genome.info[genome.info$Chrom %in% regular_chromosomes, ]

tss_annotation <- data.frame(
  chr = genome.info$Chrom,
  tss = genome.info$Starts,
  gene = genome.info$genes,
  stringsAsFactors = FALSE
)

gene_annot = tss_annotation
####

bcnet <- AddGeneAnnotation(bcnet, gene_annot)

bcnet_seurat <- AddGeneAnnotation(bcnet_seurat, gene_annot)

cat("Genes annotated:", nrow(bcnet@gene_annotation), "\n")

# both pass adding gene annotation function

###### 4. Test Aggregation ######
cat("\n=== Testing Cell Aggregation ===\n")



# Now test aggregation
bcnet <- AggregateByKNN(bcnet, k = 50)

bcnet_seurat <- AggregateByKNN(bcnet_seurat, k = 50)

cat("Aggregation complete\n")
cat("Aggregated RNA dimensions:", dim(bcnet_seurat@aggregated_data$rna), "\n")
cat("Aggregated ATAC dimensions:", dim(bcnet_seurat@aggregated_data$atac), "\n")

###### 5. Test Variable Gene Selection ######
cat("\n=== Testing Variable Gene Selection ===\n")

# Method 1: Custom genes
custom_genes <- paste0("Gene", 1:20)
bcnet_test <- AddVariableGenes(bcnet_seurat, method = "custom", genes = custom_genes)
cat("Custom genes selected:", length(bcnet_test@variable_genes), "\n")

# Method 2: Find markers (requires normalization)
bcnet_seurat <- AddVariableGenes(bcnet_seurat, method = "markers",
                                 min.pct = 0.01, logfc.threshold = 0.01)

cat("Marker genes found:", length(bcnet_seurat@variable_genes), "\n")

###### 6. Check if Ready for Inference ######
cat("\n=== Checking Readiness ===\n")
IsReadyForInference(bcnet_seurat)

###### 7. Test Small Inference ######
cat("\n=== Testing Inference (Small) ===\n")

# Run on very few genes to test
if (length(bcnet_seurat@variable_genes) > 5) {
  bcnet_seurat@variable_genes <- bcnet_seurat@variable_genes[1:20]
}

# Run inference
bcnet_result <- RunBayesCNet(bcnet_seurat,
                             window = 100000,  # Smaller window for testing
                             cores = 1,       # Single core for debugging
                             verbose = TRUE)

# Check results

CD4Naive_net <- GetCellTypeSpecificNetwork(bcnet_result, "CD4 Naive")
all_nets <- GetAllNetworks(bcnet_result)






