# Load packages
library(Seurat)
library(ggplot2)
library(RColorBrewer)
library(harmony)

folder_path <- "/Users/ahn/Library/CloudStorage/OneDrive-AdventHealth/Desktop/2 - PROJECTS/1 - PEPPER"
setwd(folder_path)

##load seurat object, set default assay and identification.
seurat <- readRDS('pepper_seurat_110225.RDS')
DefaultAssay(object = seurat) <- "RNA" 
Idents(seurat) <- seurat$CellType

# Identify marker genes per cell type
markers <- FindAllMarkers(
  seurat,
  assay = "filtered",
  only.pos = TRUE,       # Only keep upregulated genes
  min.pct = 0.25,
  logfc.threshold = 0.25
)
# Or compare two cell types directly. E.g., Adipocyte-2 genes over Adipocyte-1
markers_adip2_vs_adip1 <- FindMarkers(
  seurat,
  ident.1 = "Adip2",
  ident.2 = "Adip1",
  group.by = "CellType", 
  assay = "RNA",
  only.pos = TRUE,       
  min.pct = 0.5,
  logfc.threshold = 0.5)

# Integrated UMAP
DimPlot(seurat, reduction = "umap", group.by = "CellType") +
  theme_classic()

# Dot plot 
genes <- c('IL7R',
        "KLRD1",
        "PDGFRA",
        'NRXN3',
        'STEAP4',
        'CPA3',
        'MAFB',
        'PECAM1',
        "ADIPOQ",
        "LEP"
)
DotPlot(seurat, features = genes, assay = "SCT") + coord_flip() +
  scale_colour_gradientn(colours = rev(brewer.pal(n =11, name = 'PiYG'))) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),  # Remove x-axis title
        axis.title.y = element_blank())


## Macrophage subclustering
Macro <- subset(seurat, idents = 'Macrophage')
DefaultAssay(object = Macro) <- 'RNA'

##Recluster with harmony integration
merge.list <- SplitObject(Macro, split.by="orig.ident")
merge.list <- lapply(X = merge.list, 
                     FUN = SCTransform, 
                     method = "glmGamPoi", 
                     return.only.var.genes = FALSE)
var.features <- SelectIntegrationFeatures(object.list = merge.list, nfeatures = 2000)

merge.sct <- merge(x = merge.list[[1]], y = merge.list[2:length(merge.list)], merge.data=TRUE)
VariableFeatures(merge.sct) <- var.features
merge.sct <- RunPCA(merge.sct, verbose = FALSE)
merge.sct <- RunHarmony(merge.sct, assay.use="SCT", group.by.vars = "orig.ident")
merge.sct <- RunUMAP(merge.sct, reduction = "harmony", dims = 1:15)
merge.sct <- FindNeighbors(merge.sct, reduction = "harmony", dims = 1:15) %>% FindClusters(resolution = c(0.1, 0.2)) 

# Plot UMAP. 
DimPlot(merge.sct, pt.size = 1, label = FALSE, label.size = 7, group.by = "SCT_snn_res.0.1", reduction = "umap") +
  theme(text = element_text(size=20)) 

Macro_genes <- c("MRC1",
                  "PDGFC",
                  "F13A1",
                  "NAV2",
                  'CD163L1',
                  'SLC9A9',
                  "ARHGAP24",
                  "EDA",
                  "ABCA6",
                  'PID1',
                  'CD163',
                  'ITGAX', #CD11c
                  "CD86",
                  "TREM2",
                  "MMP9", 
                  "CHIT1",
                  "LAPTM5",
                  "LIPA",
                  "CTSD",
                  'GPNMB'
)
# Dot plot for macrophage genes
DotPlot(merge.sct, features = Macro_genes, group.by = "SCT_snn_res.0.1", scale.min = 0) + coord_flip() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  scale_colour_gradientn(colours = rev(brewer.pal(n = 11, name = "PiYG")))

# Based on marker gene expression, identity 0 is Resident M. Identity 1 is LAM. 
merge.sct@meta.data <- merge.sct@meta.data %>% mutate(CellType_macro = case_when(
  SCT_snn_res.0.1 == "0"  ~ "Resident", 
  SCT_snn_res.0.1 == "1"  ~ "LAM"
)) 
