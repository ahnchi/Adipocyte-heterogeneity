library(Seurat)
library(DESeq2)
library(msigdbr)
remotes::install_github("montilab/hypeR")
library(hypeR)
library(gtools)
library(limma)

# Pseudobulk per sample and cell type
pseudobulk_counts <- AggregateExpression(
  seurat,
  assays = "Pcoding",
  group.by = c("orig.ident", "CellType"),
  return.seurat = FALSE,
  slot = "counts"
)$Pcoding 
counts_mat <- as.matrix(pseudobulk_counts)
counts_mat <- round(counts_mat)

# Function to remove low expressed genes
keep_genes <- sapply(unique(celltypes), function(ct) {
  group_cols <- which(celltypes == ct)
  group_mat <- counts_mat[, group_cols, drop = FALSE]
  rowSums(group_mat >= 10) >= 3  # Keep genes with ≥10 counts in ≥3 samples for this cell type
})
keep_any <- rowSums(keep_genes) > 0
counts_mat_filtered <- counts_mat[keep_any, ] # 17399 protein coding genes

# Ensure all values are non-negative integers
storage.mode(counts_mat_filtered) <- "integer"

# Create a colData data.frame from sample names
sample_info <- data.frame(
  Sample = colnames(counts_mat_filtered),
  Subject = sub("^(SOMMA\\d+)_.*", "\\1", colnames(counts_mat_filtered)),
  CellType = sub("^SOMMA\\d+_", "", colnames(counts_mat_filtered))
)
rownames(sample_info) <- sample_info$Sample

# Create DESeq2 dataset
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat_filtered,
  colData = sample_info,
  design = ~ CellType
)

# Run DESeq2
dds <- DESeq(dds)

# Perform contrast (e.g., Adip2 vs Adip1)
res <- results(dds, contrast = c("CellType", "Adip2", "Adip1"))
res_ordered <- res[order(res$padj), ]
volcano_df <- as.data.frame(res_ordered) %>%
  mutate(
    gene = rownames(res_ordered),
    neglog10p = -log10(pvalue),
    direction = case_when(
      padj < 0.05 & log2FoldChange > 0 ~ "Adip2",
      padj < 0.05 & log2FoldChange < 0 ~ "Adip1",
      TRUE ~ "NS"
    )
  )

# CAMERA-PR: Adip2 vs. Adip1
stats_vec <- volcano_df$stat
names(stats_vec) <- rownames(volcano_df)  # Gene names

# Remove NAs
stats_vec <- stats_vec[!is.na(stats_vec)]

msigdb_info()
genesets <- msigdb_gsets("Homo sapiens", "C2", "CP:REACTOME", clean=TRUE)
gs_list <- as.list(genesets$list())

idx <- ids2indices(gs_list, names(stats_vec), remove.empty = TRUE)
camera_res <- cameraPR(stat = stats_vec, index = idx, inter.gene.cor = 0.01)
camera_res <- camera_res %>%
  tibble::rownames_to_column("Pathway")



### Pseudobulk in Xenium
setwd("/Users/ahn/Library/CloudStorage/OneDrive-AdventHealth/Desktop/2 - PROJECTS/1 - PEPPER/Spatial_Xenium/data")
# Load integrated seurat object (Xenium)
merge.obj <- readRDS("./pepper_xenium_seurat_qc.RDS") 

pseudobulk_counts <- AggregateExpression(
  merge.obj, 
  assays = "RNA",
  group.by = c("samples", "CellType"),
  return.seurat = FALSE,
  slot = "counts"
)$RNA 
counts_mat <- as.matrix(pseudobulk_counts)
counts_mat <- round(counts_mat)

celltypes <- sub(".*_", "", colnames(counts_mat))
keep_genes <- sapply(unique(celltypes), function(ct) {
  group_cols <- which(celltypes == ct)
  group_mat <- counts_mat[, group_cols, drop = FALSE]
  rowSums(group_mat >= 10) >= 3  # Keep genes with ≥10 counts in ≥3 samples for this cell type
})
keep_any <- rowSums(keep_genes) > 0
counts_mat_filtered <- counts_mat[keep_any, ] 
# Ensure all values are non-negative integers
storage.mode(counts_mat_filtered) <- "integer"
# Fix colData properly
sample_info <- data.frame(
  Sample = colnames(counts_mat_filtered),
  Subject = sub("_.*", "", colnames(counts_mat_filtered)),   # e.g., g10638
  CellType = sub(".*_", "", colnames(counts_mat_filtered))   # e.g., Adip1
)
rownames(sample_info) <- sample_info$Sample

# Ensure factors
sample_info$Subject <- factor(sample_info$Subject)
sample_info$CellType <- factor(sample_info$CellType)

# Build DESeq2 object
dds <- DESeqDataSetFromMatrix(
  countData = counts_mat_filtered,
  colData = sample_info,
  design = ~ Subject + CellType
)

# Run DESeq2
dds <- DESeq(dds)
