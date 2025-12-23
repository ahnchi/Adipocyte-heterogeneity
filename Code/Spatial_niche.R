library(FNN)
library(dplyr)
meta <- merge.obj_squid@meta.data

# Keep only cells with valid coords & celltype
keep <- is.finite(meta$x) & is.finite(meta$y) & !is.na(meta$CellType)
cells_keep <- rownames(meta)[keep]

coords <- as.matrix(meta[cells_keep, c("x","y")])
celltype <- droplevels(factor(meta$CellType[keep]))
subject  <- meta$samples[keep]

# Use all observed celltype levels across kept cells
all_types <- levels(celltype)
if (is.null(all_types)) all_types <- sort(unique(as.character(celltype)))

neighbors.k <- 40  # you can tune

# Compute per-subject kNN neighborhood composition
res_list <- lapply(split(cells_keep, subject), function(cells_sub) {
  idx <- match(cells_sub, cells_keep)
  X   <- coords[idx, , drop = FALSE]
  Tsub <- droplevels(celltype[idx])
  
  n <- nrow(X)
  if (n < 2) return(NULL)
  
  k <- min(neighbors.k, n - 1)  # avoid asking for more neighbors than exist
  
  knn <- FNN::get.knn(X, k = k)
  nn_index <- knn$nn.index  # matrix n x k (indices relative to 'X')
  
  # Count neighbor types for each focal cell
  mat <- matrix(0, nrow = n, ncol = length(all_types),
                dimnames = list(cells_sub, all_types))
  for (i in seq_len(n)) {
    nb_types <- Tsub[nn_index[i, ]]
    tab <- table(nb_types)
    mat[i, names(tab)] <- as.numeric(tab)
  }
  
  # Convert to proportions; guard against 0 division
  prop <- sweep(mat, 1, rowSums(mat), `/`)
  prop[!is.finite(prop)] <- 0
  prop
})

# Bind all subjects
niche_df <- do.call(rbind, res_list)

# Align to the Seurat object’s cell order; fill missing with 0
miss_cells <- setdiff(colnames(merge.obj_squid), rownames(niche_df))
if (length(miss_cells) > 0) {
  add0 <- matrix(0, nrow = length(miss_cells), ncol = ncol(niche_df),
                 dimnames = list(miss_cells, colnames(niche_df)))
  niche_df <- rbind(niche_df, add0)
}
niche_df <- niche_df[colnames(merge.obj_squid), , drop = FALSE]

# Add as an assay
merge.obj_squid[["Niche"]] <- CreateAssayObject(counts = t(niche_df))
DefaultAssay(merge.obj_squid) <- "Niche"

# Quick sanity check
GetAssayData(merge.obj_squid, assay = "Niche")[1:5, 1:5]
meta <- merge.obj_squid@meta.data

cairo_ps("xenium_niche_cell.eps", width = 4, height = 12, pointsize = 12, fallback_resolution = 300)
ggplot(merge.obj_squid@meta.data, aes(x = x, y = y, color = CellType)) +
  geom_point(size = 0.4, alpha = 0.8) +
  scale_color_manual(values = celltype_colors) +
  facet_wrap(~ samples, scales = "free", ncol = 2) +
  theme_void() 
dev.off()

meta <- merge.obj_squid@meta.data
niche <- as.data.frame(t(GetAssayData(merge.obj_squid, assay = "Niche")))

# Choose one cell type (e.g., Macrophage1)
meta$Macrophage2_prop <- niche$Macrophage2

ggplot(meta, aes(x = x, y = y, color = Macrophage2_prop)) +
  geom_point(size = 0.5) +
  scale_color_viridis_c(option = "plasma") +
  facet_wrap(~ samples, scales = "free", ncol = 2) +
  theme_void() 


niche <- as.data.frame(t(GetAssayData(merge.obj_squid, assay = "Niche")))
meta <- merge.obj_squid@meta.data

for (ct in colnames(niche)) {
  meta[[ct]] <- niche[[ct]]
  
  p <- ggplot(meta, aes(x = x, y = y, color = .data[[ct]])) +
    geom_point(size = 0.4) +
    scale_color_viridis_c(option = "plasma") +
    facet_wrap(~ samples, scales = "free", ncol = 2) +
    theme_void() +
    ggtitle(paste("Local", ct, "proportion (niche map)"))
  
  ggsave(paste0("nichemap_", ct, ".png"), p, width = 8, height = 6, dpi = 300)
}


niche_mat <- t(GetAssayData(merge.obj_squid, assay = "Niche"))

# Remove NaNs and normalize
niche_mat[is.na(niche_mat)] <- 0

# Add to Seurat for easy clustering
merge.obj_squid[["Niche"]] <- CreateAssayObject(counts = t(niche_mat))
DefaultAssay(merge.obj_squid) <- "Niche"
VariableFeatures(merge.obj_squid) <- rownames(merge.obj_squid[["Niche"]])

merge.obj_squid <- ScaleData(merge.obj_squid, features = VariableFeatures(merge.obj_squid))
merge.obj_squid <- RunPCA(merge.obj_squid, features = VariableFeatures(merge.obj_squid))
ElbowPlot(merge.obj_squid)

merge.obj_squid <- RunUMAP(merge.obj_squid, dims = 1:10)
merge.obj_squid <- FindNeighbors(merge.obj_squid, dims = 1:10)
merge.obj_squid <- FindClusters(merge.obj_squid, resolution = c(0.05, 0.1, 0.2, 0.5, 1))

png("niche_umap_res0.1.png", width = 1600, height = 1400, res = 300) # or export in pdf and rasterize the dot part in illustrator
DimPlot(merge.obj_squid, group.by = "Niche_snn_res.0.1", reduction = "umap", cols = brewer.pal(10, "Set3")) 
dev.off()

table(merge.obj_squid$Niche_snn_res.0.05)
table(merge.obj_squid$Niche_snn_res.0.2)
table(merge.obj_squid$Niche_snn_res.0.1)
table(merge.obj_squid$Niche_snn_res.0.15)

table(merge.obj_squid$Niche_snn_res.0.2) #9 niches

meta <- merge.obj_squid@meta.data
Idents(merge.obj_squid) <- "Niche_snn_res.0.1"
niche_mat <- as.data.frame(t(GetAssayData(merge.obj_squid, assay = "Niche")))
niche_mat$niche_cluster <- Idents(merge.obj_squid)

cluster_means <- niche_mat %>%
  group_by(niche_cluster) %>%
  summarize(across(everything(), mean, na.rm = TRUE)) %>%
  as.data.frame()