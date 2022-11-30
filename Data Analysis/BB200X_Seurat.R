library(Seurat); library(future); library(patchwork)
plan("multisession", workers = 10)
library(ggplot2); library(viridis); library(simspec); library(dplyr)
library(STvEA); library(SeuratDisk); library(Matrix); library(RColorBrewer)
library(gplots); library(plotly); library(hrbrthemes); library(leiden)
library(monocle3)

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
codex.5pcw$cell_id <- data_5pcw['cell_id']

codex.7pcw <- AddMetaData(codex.7pcw, metadata=data_7pcw['Area'],col.name='Area')
codex.7pcw <- AddMetaData(codex.7pcw, metadata=data_7pcw['x'],col.name='x')
codex.7pcw <- AddMetaData(codex.7pcw, metadata=data_7pcw['y'],col.name='y')
codex.7pcw$cell_id <- data_7pcw['cell_id']

sample_id_5pcw <- vector(mode='character',length=length(codex.5pcw@meta.data[["cell_id"]]))
sample_id_5pcw[1:length(codex.5pcw@meta.data[["cell_id"]])] = '5pcw'
sample_id_7pcw <- vector(mode='character',length=length(codex.7pcw@meta.data[["cell_id"]]))
sample_id_7pcw[1:length(codex.7pcw@meta.data[["cell_id"]])] = '7pcw'

codex.5pcw <- AddMetaData(codex.5pcw, metadata=sample_id_5pcw, col.name='sample')
codex.7pcw <- AddMetaData(codex.7pcw, metadata=sample_id_7pcw, col.name='sample')

## Filter
codex.5pcw <- subset(codex.5pcw, subset = Area >= 94.67 & Area <= 3786.98)
codex.7pcw <- subset(codex.7pcw, subset = Area >= 94.67 & Area <= 3786.98)

## Normalization and dimension reduction.
codex.5pcw <- NormalizeData(object = codex.5pcw, normalization.method = "CLR", margin = 2)
mean_5pcw <- apply(codex.5pcw@assays[["Akoya"]]@data, 1, mean)
sd_5pcw <- apply(codex.5pcw@assays[["Akoya"]]@data, 1, sd)
codex.5pcw@assays[["Akoya"]]@data <- (codex.5pcw@assays[["Akoya"]]@data - mean_5pcw) / sd_5pcw
codex.5pcw <- ScaleData(codex.5pcw)
VariableFeatures(codex.5pcw) <- rownames(codex.5pcw)  # since the panel is small, treat all features as variable.
codex.5pcw <- RunPCA(object = codex.5pcw, npcs = length(VariableFeatures(codex.5pcw)), verbose = FALSE, approx=FALSE)

codex.7pcw <- NormalizeData(object = codex.7pcw, normalization.method = "CLR", margin = 2)
mean_7pcw <- apply(codex.7pcw@assays[["Akoya"]]@data, 1, mean)
sd_7pcw <- apply(codex.7pcw@assays[["Akoya"]]@data, 1, sd)
codex.7pcw@assays[["Akoya"]]@data <- (codex.7pcw@assays[["Akoya"]]@data - mean_7pcw) / sd_7pcw
codex.7pcw <- ScaleData(codex.7pcw)
VariableFeatures(codex.7pcw) <- rownames(codex.7pcw)  # since the panel is small, treat all features as variable.
codex.7pcw <- RunPCA(object = codex.7pcw, npcs = length(VariableFeatures(codex.5pcw)), verbose = FALSE, approx=FALSE)


###   PCA plot (start) ----------------------------------------------------------------------------------
pca_5pcw <- codex.5pcw[["pca"]]
# Get the total variance:
eigValues_5pcw = (pca_5pcw@stdev)^2  # EigenValues
total_variance_5pcw <- sum(eigValues_5pcw)
varExplained_5pcw = cumsum(100*eigValues_5pcw / total_variance_5pcw)

pca_7pcw <- codex.7pcw[["pca"]]
# Get the total variance:
eigValues_7pcw = (pca_7pcw@stdev)^2  ## EigenValues
total_variance_7pcw <- sum(eigValues_7pcw)
varExplained_7pcw = cumsum(100*eigValues_7pcw / total_variance_7pcw)

PC_vec <- c(1:26)

par( mar= c(5,5,1,5) )
plot(PC_vec, eigValues_5pcw[1:26], pch=16, axes=FALSE, ylim=c(0,8), xlab="", ylab="", 
     type="b",col="#000000")
axis(1, xlim=c(1,26),col="black",las=1)  ## las=1 makes horizontal labels
axis(2, ylim=c(0,9),col="black",las=1)  ## las=1 makes horizontal labels
mtext("Eigenvalue",side=2,line=2.5)
mtext("Principal Component",side=1,line=2.5)
box()
par(new=TRUE)
plot(PC_vec, eigValues_7pcw[1:26], pch=16, axes=FALSE, ylim=c(0,9), xlab="", ylab="", 
     type="b",col="#CC0000")
par(new=TRUE)
## Plot the second plot and put axis scale on right
plot(PC_vec, varExplained_5pcw[1:26], pch=15,  xlab="", ylab="", ylim=c(0,100), 
     axes=FALSE, type="b", col="#999999")
## a little farther out (line=4) to make room for labels
mtext("% of Variance Explained",side=4,col="#000000",line=2.5) 
axis(4, ylim=c(0,100), col="#000000",col.axis="#000000",las=1)
par(new=TRUE)
plot(PC_vec, varExplained_7pcw[1:26], pch=15,  xlab="", ylab="", ylim=c(0,100), 
     axes=FALSE, type="b", col="#FF9933")
par(new=FALSE)
legend(15.5, 85, legend=c("Eigenval 5pcw", "Eigenval 7pcw", "VarExplained 5pcw", "VarExplained 7pcw"),
       col=c("#000000", "#CC0000", "#999999", "#FF9933"), lty=1, cex=0.8, box.lty=0)
title(main="Principal Component Variance")

###   PCA plot (end) ----------------------------------------------------------------------------------

## Clustering
codex.5pcw <- RunUMAP(object = codex.5pcw, reduction="pca", dims = 1:20, verbose = F)
codex.5pcw <- FindNeighbors(object = codex.5pcw, reduction="pca", dims = 1:20, verbose = F)

membership <- leiden(codex.5pcw@graphs$Akoya_snn, resolution_parameter=0.6, seed=1, n_iterations=-1)

levels(membership) <- c(1:length(unique(membership)))
codex.5pcw@meta.data[["seurat_clusters"]] <- membership
levels(codex.5pcw@meta.data[["seurat_clusters"]]) <- c(1:length(unique(membership)))
Idents(object = codex.5pcw) <- codex.5pcw@meta.data[["seurat_clusters"]]
Idents(codex.5pcw) <- factor(x = Idents(codex.5pcw), levels = paste(c(1:length(unique(membership)))))

codex.7pcw <- RunUMAP(object = codex.7pcw, reduction="pca", dims = 1:20, verbose = F)
codex.7pcw <- FindNeighbors(object = codex.7pcw, reduction="pca", dims = 1:20, verbose = F)

membership <- leiden(codex.7pcw@graphs$Akoya_snn, resolution_parameter=0.6, seed=1, n_iterations=-1)

levels(membership) <- c(1:length(unique(membership)))
codex.7pcw@meta.data[["seurat_clusters"]] <- membership
levels(codex.7pcw@meta.data[["seurat_clusters"]]) <- c(1:length(unique(membership)))
Idents(object = codex.7pcw) <- codex.7pcw@meta.data[["seurat_clusters"]]
Idents(codex.7pcw) <- factor(x = Idents(codex.7pcw), levels = paste(c(1:length(unique(membership)))))

## Plot UMAPs
DimPlot(codex.5pcw, label = TRUE, label.box = FALSE, size=0.1) + 
  ggplot2::ggtitle("UMAP Embedding 5pcw")
DimPlot(codex.7pcw, label = TRUE, label.box = FALSE, size=0.1) + 
  ggplot2::ggtitle("UMAP Embedding 7pcw")

###   Manual annotation (start) ----------------------------------------------------------------------------------
new.cluster.ids <- c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23)

new.cluster.ids <- c('Ukn', 'PDPN EpCs', 'EpC derived', 'V-CMs', 'FBs', 'p. Ukn', 'CD44 EpCs', 'EPCAM ECs', 
                     'ECs', 'A-CMs', 'Prolif. EpCs', 'PDPN EpCs', 'FBs', 'MPs', 'p. FBs', 'V-CMs', 'SMCs', 'EpCs',
                     'V-CMs', 'A-CMs/ECs', 'V-CMs', 'PDPN EpCs', 'GRHL2+ Ukn')
new.cluster.ids <- c('V-CMs', 'Ukn', 'SMCs', 'Immune Cells', 'Skl.-MCs', 'SMCs', 'EndoMT', 'PDPN EpCs', 
                     'ECs', 'A-CMs', 'FBs', 'ECs/FBs', 'p. Immune Cells', 'V-CMs', 'DCN+ MPs', 'V-CMs')

codex.temp <- codex.7pcw

names(new.cluster.ids) <- levels(codex.temp)
codex.temp <- RenameIdents(codex.temp, new.cluster.ids)

DimPlot(codex.temp, label = TRUE, label.box = FALSE, size=0.1) + 
  ggplot2::ggtitle("UMAP Embedding 5pcw")

DF <- data.frame(x = codex.temp@meta.data$cell_id, y = codex.temp@active.ident)

DF <- DF %>% group_by(y) %>% 
  summarise(Nb = n()) %>%
  mutate(C = sum(Nb)) %>%
  mutate(percent = Nb/C*100)

ggplot(DF, aes(x = y, y = percent, fill = y))+
  theme_classic()+
  geom_bar(stat = "identity")+
  geom_text(aes(label = paste(round(percent),"%")), position = position_stack(vjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1), axis.title.x=element_blank(),
        plot.margin = margin(10,10,30,10), legend.title=element_blank()) +
  ggtitle('Cell-type composition 5 pcw')+
  ylab('Percent % ')

DF %>%
  arrange(percent) %>%
  ggplot(DF, aes(x = y, y = percent, fill = y))+
  theme_classic()+
  geom_bar(stat = "identity")+
  geom_text(aes(label = paste(round(percent),"%")), position = position_stack(vjust = 0.5))+
  theme(axis.text.x = element_text(angle = 45, vjust = 0.5, hjust=1), axis.title.x=element_blank(),
        plot.margin = margin(10,10,30,10), legend.title=element_blank()) +
  ggtitle('Cell-type composition 5 pcw')+
  ylab('Percent % ')

###   Manual annotation (end) ----------------------------------------------------------------------------------


###   Heatmaps (start) ----------------------------------------------------------------------------------
AvgExp5pcw_mat <- matrix(nrow=length(VariableFeatures(codex.5pcw)),ncol=length(unique(codex.5pcw@meta.data[["seurat_clusters"]])))
colnames(AvgExp5pcw_mat) <- sort(unique(codex.5pcw@meta.data[["seurat_clusters"]]))
rownames(AvgExp5pcw_mat) <- VariableFeatures(codex.5pcw)
for (i in 1:length(unique(codex.5pcw@meta.data[["seurat_clusters"]]))) {
  AvgExp5pcw_mat[,i] <- unlist(apply(subset(x=codex.5pcw,idents=c(as.character(i)))@assays[["Akoya"]]@data,1,mean))
}
par(oma=c(0,0,0,0),cex.main=1)
heatmap.2(AvgExp5pcw_mat, scale = "none", col = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(n = 100), 
          trace = "none", symbreaks = T, symkey=T, symm=F, key.xlab="", key.title="Expression Level", breaks=seq(-3, 3, length.out=101),
          lhei=c(1,3),lwid=c(1,3),Rowv=F,density.info='none')

AvgExp7pcw_mat <- matrix(nrow=length(VariableFeatures(codex.7pcw)),ncol=length(unique(codex.7pcw@meta.data[["seurat_clusters"]])))
colnames(AvgExp7pcw_mat) <- sort(unique(codex.7pcw@meta.data[["seurat_clusters"]]))
rownames(AvgExp7pcw_mat) <- VariableFeatures(codex.7pcw)
for (i in 1:length(unique(codex.7pcw@meta.data[["seurat_clusters"]]))) {
  AvgExp7pcw_mat[,i] <- unlist(apply(subset(x=codex.7pcw,idents=c(as.character(i)))@assays[["Akoya"]]@data,1,mean))
}
par(oma=c(0,0,0,0),cex.main=1)
heatmap.2(AvgExp7pcw_mat, scale = "none", col = colorRampPalette(brewer.pal(n = 9, name = "YlOrRd"))(n = 100), 
          trace = "none", symbreaks = T, symkey=T, symm=F, key.xlab="", key.title="Expression Level", breaks=seq(-3, 3, length.out=101),
          lhei=c(1,3),lwid=c(1,3),Rowv=F,density.info='none')

###   Heatmaps (end) ----------------------------------------------------------------------------------


codex.test <- codex.7pcw
df_coord <- data.frame(x=codex.test@meta.data[["x"]],y=codex.test@meta.data[["y"]])
test <- CreateFOV(df_coord,type='centroids',radius=80, assay='Akoya',key='fov_', name='centroids_')
codex.test@images[["image"]] <- test
codex.test@meta.data[['seurat_clusters']] <- factor(x = Idents(codex.test), levels = paste(c(1:length(unique(codex.test@meta.data[["seurat_clusters"]])))))


path_temp <-"C:/Users/matti/Documents/KTH/BB200X Degree Project in Biotechnology, Second Cycle/Image Processing/20220322_CODEX_fetalheart_7pcw/CODEX_data_matrix_2.csv"
codex.temp <- LoadAkoya(filename = path_temp, type = "processor", fov = "fetal.heart.5pcw")
codex.temp <- NormalizeData(object = codex.temp, normalization.method = "CLR", margin = 2)
mean_temp <- apply(codex.temp@assays[["Akoya"]]@data, 1, mean)
sd_temp <- apply(codex.temp@assays[["Akoya"]]@data, 1, sd)
codex.temp@assays[["Akoya"]]@data <- (codex.temp@assays[["Akoya"]]@data - mean_temp) / sd_temp
codex.temp <- ScaleData(codex.temp)
VariableFeatures(codex.temp) <- rownames(codex.temp)  # since the panel is small, treat all features as variable.
codex.temp <- RunPCA(object = codex.temp, npcs = length(VariableFeatures(codex.5pcw)), verbose = FALSE, approx=FALSE)
codex.temp <- RunUMAP(object = codex.temp, reduction="pca", dims = 1:20, verbose = F)
codex.temp <- FindNeighbors(object = codex.temp, reduction="pca", dims = 1:20, verbose = F)
codex.temp <- FindClusters(object = codex.temp, verbose = F, resolution = 0.6, n.start = 1)


ggplotColours <- function(n = 6, h = c(0, 360) + 15){
  if ((diff(h) %% 360) < 1) h[2] <- h[2] - 360/n
  hcl(h = (seq(h[1], h[2], length = n)), c = 100, l = 65)
}

color_list <- ggplotColours(n=12)


ImageDimPlot(codex.test, size=0.5, fov='image', axes=TRUE)+ 
  theme_classic()+
  NoGrid() + 
  ggplot2::coord_flip()+
  ggplot2::scale_x_reverse()+
  ggplot2::theme(axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 panel.background = element_rect(colour = "black", size=0.1),
                 plot.title = element_text(hjust = 0.5))



ImageDimPlot(codex.test, size=0.5, fov="fetal.heart.5pcw", axes=TRUE)+ 
  NoGrid() + 
  ggplot2::theme_classic()+
  ggplot2::coord_flip()+
  ggplot2::ggtitle("Spatial 5 pcw")+
  ggplot2::scale_x_reverse()+
  ggplot2::theme(axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 plot.title = element_text(hjust = 0.5))

ImageDimPlot(codex.7pcw, size=0.5, fov="fetal.heart.7pcw", axes=TRUE)+ 
  NoGrid() + 
  ggplot2::theme_classic()+
  ggplot2::coord_flip()+
  ggplot2::ggtitle("Spatial 7 pcw")+
  ggplot2::scale_x_reverse()+
  ggplot2::theme(axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 plot.title = element_text(hjust = 0.5))


ImageFeaturePlot(codex.5pcw, fov = "fetal.heart.5pcw", features = c("CTNNB1"), min.cutoff = "q10", max.cutoff = "q90", size=2) + 
  NoGrid() + 
  ggplot2::theme_classic()+ ggplot2::scale_fill_gradientn(colors=rev(colorRampPalette(brewer.pal(n = 11, name = "RdBu"))(n = 100)))+
  ggplot2::coord_flip()+
  ggplot2::scale_x_reverse()
  ggplot2::theme(axis.line=element_blank(),axis.text.x=element_blank(),
                 axis.text.y=element_blank(),axis.ticks=element_blank(),
                 axis.title.x=element_blank(),
                 axis.title.y=element_blank(),
                 panel.background = element_rect(colour = "black", size=0.1),
                 plot.title = element_text(hjust = 0.5))



###   STVEA (Start) ------------------------------------------------------------------------------------
  # (Not used)
stvea.5pcw.protein <- as.data.frame(t(codex.5pcw@assays[["Akoya"]]@data))

stvea.5pcw.spatial <- data.frame('x' = c(1:length(codex.5pcw$cell_id))*0, 'y' = c(1:length(codex.5pcw$cell_id))*0)
stvea.5pcw.spatial$x <- codex.5pcw$x
stvea.5pcw.spatial$y <- codex.5pcw$y
stvea.5pcw.spatial$x <- stvea.5pcw.spatial$x - (min(stvea.5pcw.spatial$x) - 1)
stvea.5pcw.spatial$y <- stvea.5pcw.spatial$y - (min(stvea.5pcw.spatial$y) - 1)
stvea.5pcw.spatial_nm <- as.data.frame(cbind(x=stvea.5pcw.spatial$x*325, y=stvea.5pcw.spatial$y*325))

stvea.5pcw.area <- codex.5pcw$Area

stvea.5pcw.blanks <- data.frame('blank1_CH2' = c(1:length(codex.5pcw$cell_id))*0, 'blank1_CH3' = c(1:length(codex.5pcw$cell_id))*0+1, 
                               'blank1_CH4' = c(1:length(codex.5pcw$cell_id))*0, 'blank2_CH2' = c(1:length(codex.5pcw$cell_id))*0+1,
                               'blank2_CH3' = c(1:length(codex.5pcw$cell_id))*0, 'blank2_CH4' = c(1:length(codex.5pcw$cell_id))*0+1)

stvea.5pcw <- SetDataCODEX(codex_protein = stvea.5pcw.protein,
                          codex_size = stvea.5pcw.area,
                          codex_spatial = stvea.5pcw.spatial_nm,
                          codex_blanks = stvea.5pcw.blanks)

stvea.5pcw <- CleanCODEX(stvea.5pcw, model = "nb")

stvea.5pcw@codex_clusters <- codex.5pcw$seurat_clusters

protein_adj_5pcw <- AdjScoreProteins(stvea.5pcw, k=5, num_cores=1)
AdjScoreHeatmap(protein_adj_5pcw)

heatmap_matrix <- matrix(rep(0,length(unique(protein_adj_5pcw$f))*length(unique(protein_adj_5pcw$g))), ncol=length(unique(protein_adj_5pcw$g)))
row.names(heatmap_matrix) <- unique(protein_adj_5pcw$f)[order(unique(protein_adj_5pcw$f))]
colnames(heatmap_matrix) <- unique(protein_adj_5pcw$g)[order(unique(protein_adj_5pcw$g))]
for (i in 1:nrow(protein_adj_5pcw)) {
  heatmap_matrix[protein_adj_5pcw[i,"f"],protein_adj_5pcw[i,"g"]] <- log10(protein_adj_5pcw[i,"q"]+1e-15)
  heatmap_matrix[protein_adj_5pcw[i,"g"],protein_adj_5pcw[i,"f"]] <- log10(protein_adj_5pcw[i,"q"]+1e-15)
}

par(oma=c(0,0,0,0),cex.main=1)
heatmap.2(heatmap_matrix, scale = "none", col = colorRampPalette(c('#6666FF','#000066'))(n = 100),dendrogram = 'both', 
          trace = "none", symbreaks = F, symkey=F, symm=F, density.info="none", revC=T, key.xlab="",key.title="adj. p-val.",
          lhei=c(1,3),lwid=c(1,3),colsep=1:nrow(heatmap_matrix),rowsep=1:nrow(heatmap_matrix),sepcolor='#336633',
          sepwidth=c(0.001,0.001))



cluster_adj_5pcw <- AdjScoreClustersCODEX(stvea.5pcw, k=5)
AdjScoreHeatmap(cluster_adj_5pcw)

heatmap_matrix <- matrix(rep(0,length(unique(cluster_adj_5pcw$f))*length(unique(cluster_adj_5pcw$g))), ncol=length(unique(cluster_adj_5pcw$g)))
row.names(heatmap_matrix) <- unique(cluster_adj_5pcw$f)[order(unique(cluster_adj_5pcw$f))]
colnames(heatmap_matrix) <- unique(cluster_adj_5pcw$g)[order(unique(cluster_adj_5pcw$g))]
for (i in 1:nrow(cluster_adj_5pcw)) {
  heatmap_matrix[cluster_adj_5pcw[i,"f"],cluster_adj_5pcw[i,"g"]] <- log10(cluster_adj_5pcw[i,"q"]+1e-15)
  heatmap_matrix[cluster_adj_5pcw[i,"g"],cluster_adj_5pcw[i,"f"]] <- log10(cluster_adj_5pcw[i,"q"]+1e-15)
}
par(oma=c(0,0,0,0),cex.main=1)
heatmap.2(heatmap_matrix, scale = "none", col = colorRampPalette(c('#990033','#FFCC00'))(n = 100),dendrogram = 'both', 
          trace = "none", symbreaks = F, symkey=F, symm=F, density.info="none", revC=T, key.xlab="",key.title="adj. p-val.",
          lhei=c(1,3),lwid=c(1,3),colsep=1:nrow(heatmap_matrix),rowsep=1:nrow(heatmap_matrix),sepcolor='#996600',
          sepwidth=c(0.005,0.005))



stvea.7pcw.protein <- as.data.frame(t(codex.7pcw@assays[["Akoya"]]@data))

stvea.7pcw.spatial <- data.frame('x' = c(1:length(codex.7pcw$cell_id))*0, 'y' = c(1:length(codex.7pcw$cell_id))*0)
stvea.7pcw.spatial$x <- codex.7pcw$x
stvea.7pcw.spatial$y <- codex.7pcw$y
stvea.7pcw.spatial$x <- stvea.7pcw.spatial$x - (min(stvea.7pcw.spatial$x) - 1)
stvea.7pcw.spatial$y <- stvea.7pcw.spatial$y - (min(stvea.7pcw.spatial$y) - 1)
stvea.7pcw.spatial_nm <- as.data.frame(cbind(x=stvea.7pcw.spatial$x*325, y=stvea.7pcw.spatial$y*325))

stvea.7pcw.area <- codex.7pcw$Area

stvea.7pcw.blanks <- data.frame('blank1_CH2' = c(1:length(codex.7pcw$cell_id))*0, 'blank1_CH3' = c(1:length(codex.7pcw$cell_id))*0+1, 
                                'blank1_CH4' = c(1:length(codex.7pcw$cell_id))*0, 'blank2_CH2' = c(1:length(codex.7pcw$cell_id))*0+1,
                                'blank2_CH3' = c(1:length(codex.7pcw$cell_id))*0, 'blank2_CH4' = c(1:length(codex.7pcw$cell_id))*0+1)

stvea.7pcw <- SetDataCODEX(codex_protein = stvea.7pcw.protein,
                           codex_size = stvea.7pcw.area,
                           codex_spatial = stvea.7pcw.spatial_nm,
                           codex_blanks = stvea.7pcw.blanks)

stvea.7pcw <- CleanCODEX(stvea.7pcw, model = "nb")

stvea.7pcw@codex_clusters <- codex.7pcw$seurat_clusters

protein_adj_7pcw <- AdjScoreProteins(stvea.7pcw, k=5, num_cores=1)
AdjScoreHeatmap(protein_adj_7pcw)

heatmap_matrix <- matrix(rep(0,length(unique(protein_adj_7pcw$f))*length(unique(protein_adj_7pcw$g))), ncol=length(unique(protein_adj_7pcw$g)))
row.names(heatmap_matrix) <- unique(protein_adj_7pcw$f)[order(unique(protein_adj_7pcw$f))]
colnames(heatmap_matrix) <- unique(protein_adj_7pcw$g)[order(unique(protein_adj_7pcw$g))]
for (i in 1:nrow(protein_adj_7pcw)) {
  heatmap_matrix[protein_adj_7pcw[i,"f"],protein_adj_7pcw[i,"g"]] <- log10(protein_adj_7pcw[i,"q"]+1e-15)
  heatmap_matrix[protein_adj_7pcw[i,"g"],protein_adj_7pcw[i,"f"]] <- log10(protein_adj_7pcw[i,"q"]+1e-15)
}

par(oma=c(0,0,0,0),cex.main=1)
heatmap.2(heatmap_matrix, scale = "none", col = colorRampPalette(c('#6666FF','#000066'))(n = 100),dendrogram = 'both', 
          trace = "none", symbreaks = F, symkey=F, symm=F, density.info="none", revC=T, key.xlab="",key.title="adj. p-val.",
          lhei=c(1,3),lwid=c(1,3),colsep=1:nrow(heatmap_matrix),rowsep=1:nrow(heatmap_matrix),sepcolor='#336633',
          sepwidth=c(0.001,0.001))

cluster_adj_7pcw <- AdjScoreClustersCODEX(stvea.7pcw, k=5)
AdjScoreHeatmap(cluster_adj_7pcw)

heatmap_matrix <- matrix(rep(0,length(unique(cluster_adj_7pcw$f))*length(unique(cluster_adj_7pcw$g))), ncol=length(unique(cluster_adj_7pcw$g)))
row.names(heatmap_matrix) <- unique(cluster_adj_7pcw$f)[order(unique(cluster_adj_7pcw$f))]
colnames(heatmap_matrix) <- unique(cluster_adj_7pcw$g)[order(unique(cluster_adj_7pcw$g))]
for (i in 1:nrow(cluster_adj_7pcw)) {
  heatmap_matrix[cluster_adj_7pcw[i,"f"],cluster_adj_7pcw[i,"g"]] <- log10(cluster_adj_7pcw[i,"q"]+1e-15)
  heatmap_matrix[cluster_adj_7pcw[i,"g"],cluster_adj_7pcw[i,"f"]] <- log10(cluster_adj_7pcw[i,"q"]+1e-15)
}

par(oma=c(0,0,0,0),cex.main=1)
heatmap.2(heatmap_matrix, scale = "none", col = colorRampPalette(c('#990033','#FFCC00'))(n = 100),dendrogram = 'both', 
          trace = "none", symbreaks = F, symkey=F, symm=F, density.info="none", revC=T, key.xlab="",key.title="adj. p-val.",
          lhei=c(1,3),lwid=c(1,3),colsep=1:nrow(heatmap_matrix),rowsep=1:nrow(heatmap_matrix),sepcolor='#996600',
          sepwidth=c(0.005,0.005))
###   STVEA (End) ------------------------------------------------------------------------------------

###   Export data (start) -------------------------------------------------------------------------------------
ObjNo <- codex.5pcw@meta.data[["cell_id"]]
x <- codex.5pcw@meta.data[["x"]]
y <- codex.5pcw@meta.data[["y"]]
Area <- codex.5pcw@meta.data[["Area"]]
data <- as.data.frame(t(codex.5pcw@assays[["Akoya"]]@data))
cluster <- codex.5pcw@meta.data[["seurat_clusters"]]
umap_emb <- codex.5pcw@reductions[["umap"]]@cell.embeddings

csv_export <- data.frame(ObjNo,x,y,Area,data,cluster,umap_emb)
write.csv(csv_export, file = "C:/Users/matti/Desktop/utag_5pcw_test.csv")


ObjNo <- codex.7pcw@meta.data[["cell_id"]]
x <- codex.7pcw@meta.data[["x"]]
y <- codex.7pcw@meta.data[["y"]]
Area <- codex.7pcw@meta.data[["Area"]]
data <- as.data.frame(t(codex.7pcw@assays[["Akoya"]]@data))
cluster <- codex.7pcw@meta.data[["seurat_clusters"]]
umap_emb <- codex.7pcw@reductions[["umap"]]@cell.embeddings

csv_export <- data.frame(ObjNo,x,y,Area,data,cluster,umap_emb)
write.csv(csv_export, file = "C:/Users/matti/Desktop/utag_7pcw_test.csv")

###   Export data (end) -------------------------------------------------------------------------------------

## Compare expression levels between 5 and 7 pcw.
codex.obj <- merge(codex.5pcw, codex.7pcw, add.cell.ids=c('5pcw','7pcw'))
codex.obj <- ScaleData(codex.obj)

Idents(object = codex.obj) <- codex.obj@meta.data$'sample'
codex.obj@assays[["Akoya"]]@var.features <- codex.7pcw@assays[["Akoya"]]@var.features

markers_seurat1_vs_seurat2 <- FindMarkers(codex.obj, ident.1 = "7pcw", ident.2 = "5pcw", only.pos = FALSE, min.pct = 0, logfc.threshold = 0, method = "bimod")
markers_seurat1_vs_seurat2$p_val_adj = p.adjust(markers_seurat1_vs_seurat2$p_val, method='fdr')

df_pcw_expression <- data.frame(x=rownames(markers_seurat1_vs_seurat2), y=markers_seurat1_vs_seurat2$avg_log2FC)
df_pcw_expression <- df_pcw_expression[order(df_pcw_expression$y, decreasing = TRUE), ]
df_pcw_expression$x <- factor(df_pcw_expression$x, levels = df_pcw_expression$x)

ggplot(df_pcw_expression, aes(x, y)) +
  geom_bar(stat = "identity", fill=colorRampPalette(brewer.pal(n = 9, name = "Spectral"))(n = 20))+
  theme_light()+
  theme(legend.position = "None",axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1))+
  geom_text(aes(label = paste(round(y,1)),vjust = ifelse(y >= 0, 0, 1)))+
  geom_hline(yintercept=0)+geom_vline(xintercept=0)+
  ggtitle('Global DE between 7 pcw & 5 pcw')+
  theme(axis.title.x=element_blank())+
  ylab('Avg. Log2 FC')+ylim(c(-10,10))


df_pcw_expression <- data.frame(x=rownames(markers_seurat1_vs_seurat2), y=markers_seurat1_vs_seurat2$avg_log2FC, q=markers_seurat1_vs_seurat2$p_val_adj)
df_pcw_expression <- df_pcw_expression[order(df_pcw_expression$y, decreasing = TRUE), ]
df_pcw_expression <- df_pcw_expression[df_pcw_expression$q < 0.01, ]
df_pcw_expression$x <- factor(df_pcw_expression$x, levels = df_pcw_expression$x)



