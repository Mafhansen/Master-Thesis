library(Seurat); library(future); library(patchwork)
plan("multisession", workers = 10)
library(ggplot2); library(viridis); library(simspec); library(dplyr)
library(STvEA); library(SeuratDisk); library(Matrix); library(RColorBrewer)
library(gplots); library(plotly); library(hrbrthemes); library(leiden)

## Define paths to measurement .csv files.
path_5pcw <- "C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220420_CODEX_fetalheart_5pcw/CODEX_data_matrix.csv"
path_7pcw <- "C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220322_CODEX_fetalheart_7pcw/CODEX_data_matrix.csv"

## Load measurement .csv files and prepare for Seurat object construction.
data_5pcw <- read.csv(path_5pcw)
data_7pcw <- read.csv(path_7pcw)

for (i in 1:length(data_5pcw[1,])) {
  temp <- tail(unlist(strsplit(names(data_5pcw)[i],'\\.')),n=1)
  names(data_5pcw)[i] = temp
}
drops <- c("Empty","Blank1","Blank2")
data_5pcw <- data_5pcw[ , !(names(data_5pcw) %in% drops)]
rownames(data_5pcw) <- data_5pcw['cell_id'][,1]

for (i in 1:length(data_7pcw[1,])) {
  temp <- tail(unlist(strsplit(names(data_7pcw)[i],'\\.')),n=1)
  names(data_7pcw)[i] = temp
}
drops <- c("Empty","Blank1","Blank2")
data_7pcw <- data_7pcw[ , !(names(data_7pcw) %in% drops)]
rownames(data_7pcw) <- data_7pcw['cell_id'][,1]

## Create Seurat objects
codex.5pcw <- CreateSeuratObject(t(data_5pcw[,7:32]), assay='Akoya')
codex.7pcw <- CreateSeuratObject(t(data_7pcw[,7:32]), assay='Akoya')

codex.5pcw <- AddMetaData(codex.5pcw, metadata=data_5pcw['Area'],col.name='Area')
codex.5pcw <- AddMetaData(codex.5pcw, metadata=data_5pcw['x'],col.name='x')
codex.5pcw <- AddMetaData(codex.5pcw, metadata=data_5pcw['y'],col.name='y')

codex.7pcw <- AddMetaData(codex.7pcw, metadata=data_7pcw['Area'],col.name='Area')
codex.7pcw <- AddMetaData(codex.7pcw, metadata=data_7pcw['x'],col.name='x')
codex.7pcw <- AddMetaData(codex.7pcw, metadata=data_7pcw['y'],col.name='y')

sample_id_5pcw <- vector(mode='character',length=length(codex.5pcw@meta.data[["cell_id"]]))
sample_id_5pcw[1:length(codex.5pcw@meta.data[["cell_id"]])] = '5pcw'
sample_id_7pcw <- vector(mode='character',length=length(codex.7pcw@meta.data[["cell_id"]]))
sample_id_7pcw[1:length(codex.7pcw@meta.data[["cell_id"]])] = '7pcw'

codex.5pcw <- AddMetaData(codex.5pcw, metadata=sample_id_5pcw, col.name='sample')
codex.7pcw <- AddMetaData(codex.7pcw, metadata=sample_id_7pcw, col.name='sample')

## Filter.
codex.5pcw <- subset(codex.5pcw, subset = Area >= 94.67 & Area <= 3786.98)
codex.7pcw <- subset(codex.7pcw, subset = Area >= 94.67 & Area <= 3786.98)

## Merge objects.
codex.obj <- merge(codex.5pcw, codex.7pcw, add.cell.ids=c('5pcw','7pcw'))
rm(codex.5pcw); rm(codex.7pcw)

## Normalization and dimension reduction.
codex.obj <- NormalizeData(object = codex.obj, normalization.method = "CLR", margin = 2)
mean_obj <- apply(codex.obj@assays[["Akoya"]]@data, 1, mean)
sd_obj <- apply(codex.obj@assays[["Akoya"]]@data, 1, sd)
codex.obj@assays[["Akoya"]]@data <- (codex.obj@assays[["Akoya"]]@data - mean_obj) / sd_obj
codex.obj <- ScaleData(codex.obj)
VariableFeatures(codex.obj) <- rownames(codex.obj)  # since the panel is small, treat all features as variable.
codex.obj <- RunPCA(object = codex.obj, npcs = 26, verbose = FALSE, approx=FALSE)

ElbowPlot(codex.obj, ndims = ncol(Embeddings(codex.obj, "pca")))

## Clustering.
codex.obj <- cluster_sim_spectrum(codex.obj, use_scale=T, cluster_resolution=0.4, dims_use=1:20, label_tag='sample', 
                                  corr_method="spearman", spectrum_type="corr_ztransform", redo_pca=T, num_pcs_compute=20, num_pcs_use=20)
codex.obj <- RunUMAP(object = codex.obj, reduction="css", dims = 1:ncol(Embeddings(codex.obj, "css")), verbose = FALSE)
codex.obj <- FindNeighbors(object = codex.obj, reduction="css", dims = 1:ncol(Embeddings(codex.obj, "css")), verbose = FALSE)

membership <- leiden(codex.obj@graphs$Akoya_snn, resolution_parameter=0.6, seed=1)

levels(membership) <- c(1:length(unique(membership)))
codex.obj@meta.data[["seurat_clusters"]] <- membership
levels(codex.obj@meta.data[["seurat_clusters"]]) <- c(1:length(unique(membership)))
Idents(object = codex.obj) <- codex.obj@meta.data[["seurat_clusters"]]
Idents(codex.obj) <- factor(x = Idents(codex.obj), levels = paste(c(1:length(unique(membership)))))

## Plot UMAPs
UMAPPlot(codex.obj, label=TRUE, size=0.05)+
  ggplot2::ggtitle("CSS data integration")
UMAPPlot(codex.obj,label=TRUE, size=0.05, group.by="sample")+
  ggplot2::ggtitle("CSS data integration")

DimPlot(codex.obj, label = TRUE, label.box = FALSE, size=0.1, group.by="sample") + 
  ggplot2::theme(panel.background = element_rect(colour = "black", size=0.1))

## Plot Heatmaps
AvgExpObj_mat <- matrix(nrow=length(VariableFeatures(codex.obj)),ncol=length(unique(codex.obj@meta.data[["seurat_clusters"]])))
colnames(AvgExpObj_mat) <- sort(unique(codex.obj@meta.data[["seurat_clusters"]]))
rownames(AvgExpObj_mat) <- VariableFeatures(codex.obj)
for (i in 1:length(unique(codex.obj@meta.data[["seurat_clusters"]]))) {
  if (i <= 25) {
    AvgExpObj_mat[,i] <- unlist(apply(subset(x=codex.obj,idents=c(as.character(i)))@assays[["Akoya"]]@data,1,mean))
  } else {
    AvgExpObj_mat[,i] <- apply(as.data.frame(t(codex.obj@assays[['Akoya']]@data))[which(Idents(codex.obj)==26),],2,mean)
  }
}
par(oma=c(0,0,0,0),cex.main=1)
heatmap.2(AvgExpObj_mat, scale = "none", col = rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(n = 100)), 
          trace = "none", symbreaks = T, symkey=T, symm=F, key.xlab="", key.title="Expression Level", breaks=seq(-3, 3, length.out=101),
          lhei=c(1,3),lwid=c(1,3))

DoHeatmap(codex.obj) + NoLegend()

## Find differentially expressed markers.
cl_markers <- FindAllMarkers(codex.obj, only.pos = FALSE, min.pct = 0.25, logfc.threshold = 0.50, test.use = "bimod")
cl_markers %>% group_by(cluster) #%>% top_n(n = 10, wt = avg_log2FC)
cl_markers$p_val_adj = p.adjust(cl_markers$p_val, method='fdr')


#df_coord <- data.frame(x=codex.obj@meta.data[["x"]],y=codex.obj@meta.data[["y"]])
#
#test <- CreateFOV(df_coord,type='centroids',radius=80, assay='Akoya',key='fov_', name='centroids_')
#
#codex.obj <- DefaultFOV(codex.obj, assay='Akoya',value='test')
#
#codex.obj@images[["image"]] <- test
#
#ImageDimPlot(codex.obj, size=0.5, fov='test', axes=TRUE)+ 
#  NoGrid() + 
#  ggplot2::coord_flip()+
#  ggplot2::scale_x_reverse()+
#  ggplot2::theme(axis.title.x=element_blank(),
#                 axis.title.y=element_blank(),
#                 panel.background = element_rect(colour = "black", size=0.1),
#                 plot.title = element_text(hjust = 0.5))
