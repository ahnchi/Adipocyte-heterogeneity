library(granulator)
library(dplyr)
library(Seurat)

# Acquire signature matrix
# Import snRNAseq object
seurat <- readRDS("/Users/ahn/Library/CloudStorage/OneDrive-AdventHealth/Desktop/2 - PROJECTS/1 - PEPPER/MAC/somma at/seurat_Merged_SOMMA_Harmony.RDS")
Idents(seurat) <- seurat$CellType
DefaultAssay(object = seurat) <- 'RNA'

avg_exp <- AverageExpression(seurat, return.seurat = FALSE)$RNA
sig <- as.matrix(avg_exp) # sig is the full signature matrix

## Obtain HVG
hvg_counts <- seq(1000, 10000, by = 1000)

# Store HVG lists
hvg_list <- lapply(hvg_counts, function(n) {
  obj <- FindVariableFeatures(seurat, selection.method = "vst", nfeatures = n)
  VariableFeatures(obj)
})

# Name the list
names(hvg_list) <- paste0("HVG_", hvg_counts) # Sets of HVG stored in hvg_list

## signature matrix
sig1 <- sig[rownames(sig) %in% hvg_list$HVG_9000, ] # Using top 9000 HVG as signature matrix
sig1<-  as.matrix(sig1)

# Import bulk RNAseq data that will be deconvoluted. 
MS <- read.csv("/Users/ahn/Library/CloudStorage/OneDrive-AdventHealth/Desktop/2 - PROJECTS/1 - PEPPER/MAC/somma at/somma.at.bulk.matrix.csv", row.names = 1)
MS <- as.matrix(MS)


# Run deconvoltion
decon <- deconvolute(m = MS, sigMatrix = sig1, methods = "dtangle")

# Extract cell proportion data
dtangle_props = decon$proportions$dtangle_sig1 


## Deconvolution in MHO study
MS2 <- read.csv("/Users/ahn/Library/CloudStorage/OneDrive-AdventHealth/Desktop/2 - PROJECTS/1 - PEPPER/pepper_decon/Adipose_bulk_MHOstudy.csv", row.names = 1)
MS2 <- as.matrix(MS2)
decon_MS2 <- deconvolute(m = MS2, sigMatrix = sig1, methods = "dtangle")
# Extract cell proportion data
dtangle_props_MS2 = decon_MS2$proportions$dtangle_sig1 
