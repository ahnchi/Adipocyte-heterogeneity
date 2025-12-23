library(CellChat)
library(patchwork)
options(stringsAsFactors = FALSE)
library(reticulate)
library(uwot)
# Load snRNAseq seurat object
folder_path <- "/Users/ahn/Library/CloudStorage/OneDrive-AdventHealth/Desktop/2 - PROJECTS/1 - PEPPER"
setwd(folder_path)

##load seurat object, set default assay and identification.
seurat <- readRDS('pepper_seurat_110225.RDS')

cellchat <- createCellChat(object = seurat, group.by = "CellType", assay = "RNA")
seurat@meta.data$CellType <- as.factor(seurat@meta.data$CellType)

cellchat <- setIdent(cellchat, ident.use = "CellType") # set "labels" as default cell identity

groupSize <- as.numeric(table(cellchat@meta$CellType)) # number of cells in each cell group

CellChatDB <- CellChatDB.human 

dplyr::glimpse(CellChatDB$interaction)
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") # use Secreted Signaling
CellChatDB.use <- subsetDB(CellChatDB, search = "ECM-Receptor") # use Secreted Signaling
CellChatDB.use <- subsetDB(CellChatDB, search = "Cell-Cell Contact") # use Secreted Signaling
CellChatDB.use <- subsetDB(CellChatDB, search = "Non-protein Signaling") # use Secreted Signaling

cellchat@DB <- CellChatDB.use

cellchat <- subsetData(cellchat) 
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

cellchat <- computeCommunProb(cellchat, nboot = 100)
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)

cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

groupSize <- as.numeric(table(cellchat@meta$CellType))

# Get cell type levels
cell_types <- levels(cellchat@idents)


# Use in netVisual_circle
netVisual_circle(
  cellchat@net$count,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Number of interactions"
)

netVisual_circle(
  cellchat@net$weight,
  vertex.weight = groupSize,
  weight.scale = TRUE,
  label.edge = FALSE,
  title.name = "Interaction weights/strength"
)
