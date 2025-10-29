setwd("/mnt/8w/data7/yiyuan/Bone/test1/NEW_MPs")
library(Seurat)
library(Matrix)
library(harmony)
library(ggplot2)
library(tidyverse)
sce = readRDS("/mnt/8w/data7/yiyuan/Bone/test1/MP/sce_second.1.rds")
head(sce@meta.data)

png(file = "second_cluster.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce, 
  group.by = "cluster", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
)
dev.off()

png(file = "ElbowPlot.png", width = 3000, height = 2500, res = 300)
ElbowPlot(sce, ndims = 40, reduction = "pca")
dev.off()

#看看前20个主成分
pdf(file="PCA.pdf",width=10,height=8)
DimHeatmap(sce, dims = 1:20, cells = 100, balanced = TRUE,nfeatures = 10)
dev.off()
#看看前20个主成分
pdf(file="PCA.pdf",width=10,height=8)
DimHeatmap(sce, dims = 1:20, cells = 100, balanced = TRUE,nfeatures = 10)
dev.off()

pct <- sce[["pca"]]@stdev / sum(sce[["pca"]]@stdev) * 100 ; 
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)

pc.use
#13
pcs=pc.use


set.seed(123)
sce <- RunUMAP(sce,reduction = "harmony", dims = 1:pcs)
sce <- RunTSNE(sce, dims.use = 1:pcs, reduction = "harmony",do.fast = TRUE)
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:pcs)
sce <- FindClusters(sce, algorithm = 1 , resolution = c(0.5),cluster.name="MP_0.5")
sce <- FindClusters(sce, algorithm = 1 , resolution = c(1),cluster.name="MP_1")



png(file = "second_MP1.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce,
  reduction = "umap", 
  group.by = "MP_1", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
)
dev.off()



sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25,
                              thresh.use = 0.25)
write.csv(sce.markers,file=paste0('MP_all.markers.csv'))

top5 <- sce.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
write.csv(top5,file=paste0('MP_top5.markers.csv'))
pdf(file="all.markers_top5_heatmap.pdf",width=30,height=10)
DoHeatmap(sce,top5$gene,size=3)
dev.off()
sce@meta.data %>% head()


Osteoclasts = c("MMP9","CTSK","ACP5","CKB")
MastCells = c("TPSB2","TPSAB1","CPA3","HPGD","CSF1")
BCells = c('IGKC','CD69','IGLC2')
Neutrophils = c('S100A8','S100A9','IFTM2')
Macrophages = c('RNASE1','C1QB','SELENOP')
Monocytes = c('S100A9','FCN1','EREG')
DCs = c("NAPSB","CPVL","HLA-DPA1","HLA-DPB1","HLA-DRB5","CD1C","LAMP3","CCL22","CCL17","PLD4","GZMB")
cDCs = c('FCER1A','CD1C','AREG')
pDCs = c('PTGDS','GZMB','JCHAIN')



pdf(file="Dotplot_seurat_clusters.pdf",width=30,height=8)
DotPlot(sce,group.by = 'seurat_clusters',cols = c('blue','red'),features=unique(cellMarker),cluster.idents = F)
dev.off()

cellMarker <- c(Osteoclasts,MastCells,BCells,Neutrophils,Macrophages,Monocytes,DCs,cDCs,pDCs)

pdf(file="Dotplot_seurat_clusters2.pdf",width=50,height=8)
DotPlot_fixed(sce,group.by = 'seurat_clusters',cols = c('blue','red'),features=unique(cellMarker),cluster.idents = F)+
 scale_fill_gsea()+
 theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 15,
                                   angle = 90,
                                   hjust = 1,
                                   vjust = .5),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
dev.off()


celltype=data.frame(ClusterID=0:16,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(7),2]='Osteoclasts'
celltype[celltype$ClusterID %in% c(15),2]='MastCells'
celltype[celltype$ClusterID %in% c(),2]='BCells'
celltype[celltype$ClusterID %in% c(),2]='Neutrophils'
celltype[celltype$ClusterID %in% c(0,1,2,3,4,5),2]='Macrophages'
celltype[celltype$ClusterID %in% c(),2]='Monocytes'
celltype[celltype$ClusterID %in% c(8,9,10,11,13,14,6,12,16),2]='DCs'
celltype[celltype$ClusterID %in% c(),2]='unknow'
head(celltype)


sce@meta.data$celltype_new = "unknow"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_new'] <- celltype$celltype[i]}
table(sce@meta.data$celltype_new)

Osteoclasts = c("MMP9","CTSK","ACP5","CKB")
MastCells = c("TPSB2","TPSAB1","CPA3","HPGD","CSF1")
Macrophages = c('RNASE1','C1QB','SELENOP')
DCs = c("CPVL","HLA-DPA1","HLA-DPB1","HLA-DRB5")
cellMarker1 <- c(Osteoclasts,MastCells,Macrophages,DCs)
sce$celltype_new <- factor(sce$celltype_new,levels=c('DCs','Macrophages','MastCells','Osteoclasts'))

pdf(file="Dotplot_type_new.pdf",width=20,height=8)
DotPlot_fixed(sce,group.by = 'celltype_new',cols = c('blue','red'),features=unique(cellMarker),cluster.idents = F)+
 scale_fill_gsea()+
 theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 15,
                                   angle = 90,
                                   hjust = 1,
                                   vjust = .5),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
dev.off()



png(file = "second_celltype_umap.png", width = 3300, height = 2500, res = 300)
DimPlot(
  sce,
  reduction = "umap", 
  group.by = "celltype_new", 
  label = TRUE,
  label.size = 5,  # 设置标签的字体大小
   cols = c(
   'DCs' = '#ef070f',  
   'MastCells' = '#aba3e7',  
   'Macrophages' = '#D4E157',  
   'Osteoclasts' = '#A8DADC'

))
dev.off()

cell.Patient <- aggregate(sce$sample, list(sce$sample, sce$gname, sce$celltype_new), length)
colnames(cell.Patient) <- c("sample","gname","cellType","number")
cell.Patient$Patient.type <- paste(cell.Patient$sample, cell.Patient$gname, sep = "-")
sce@meta.data$Patient.type <- paste(sce@meta.data$sample, sce@meta.data$gname, sep = "-")
patient <- unique(sce$Patient.type)
for(i in patient){
    if(i == patient[1]){
        cell.pro.Patient <- subset(cell.Patient, Patient.type == i)
        cell.pro.Patient$proportion <- cell.pro.Patient$number/sum(cell.pro.Patient$number)
    } else {
        single.patient <- subset(cell.Patient, Patient.type == i)
        single.patient$proportion <- single.patient$number/sum(single.patient$number)
        cell.pro.Patient <- rbind(cell.pro.Patient, single.patient)
    }    
}
cell.pro.Patient.OA <- cell.pro.Patient[cell.pro.Patient$gname == "OA",]
cell.pro.Patient.HA <- cell.pro.Patient[cell.pro.Patient$gname == "HA",]
cell.pro.Patient.OA$x_start <- as.numeric(factor(cell.pro.Patient.OA$Patient.type)) - 0.5
cell.pro.Patient.OA$x_end <- as.numeric(factor(cell.pro.Patient.OA$Patient.type)) + 0.5
cell.pro.Patient.HA$x_start <- as.numeric(factor(cell.pro.Patient.HA$Patient.type)) + length(unique(cell.pro.Patient.OA$Patient.type)) - 0.5
cell.pro.Patient.HA$x_end <- as.numeric(factor(cell.pro.Patient.HA$Patient.type)) + length(unique(cell.pro.Patient.OA$Patient.type)) + 0.5
cell.pro.Patient <- rbind(cell.pro.Patient.OA, cell.pro.Patient.HA)
cell.pro.Patient$Patient.type <- factor(cell.pro.Patient$Patient.type, level = c(cell.pro.Patient$Patient.type[grep("OA",cell.pro.Patient$Patient.type)]%>%unique(),
                                               cell.pro.Patient$Patient.type[grep("HA",cell.pro.Patient$Patient.type)]%>%unique()))
library("ggalluvial")
library("rlang")
library(ggplot2)



patient_col = c(
   'DCs' = '#ef070f',  
   'MastCells' = '#aba3e7',  
   'Macrophages' = '#D4E157', 
   'Osteoclasts' = '#A8DADC'
)

# 设置图形的高度和宽度
pdf(file = "proportion_second.pdf", width = 10, height = 12)

# 绘制堆积柱状图
ggplot(data = cell.pro.Patient, mapping = aes(x = Patient.type, fill = cellType, y = proportion)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = patient_col) +
  labs(x = "Samples", y = "Ratio of Cell Types (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_segment(mapping = aes(x = x_start, y = 1, xend = x_end, yend = 1, color = gname), size = 4)
dev.off()

library(tidyr)
library(reshape2)
sce$sample = factor(sce$sample,levels = c('OA1','OA2','OA3','HA1','HA2','HA3'))
tb=table(sce$sample, sce$celltype_new)
head(tb)
library (gplots)

png(file = "balloonplot_1_celltype.png", width = 4000, height = 2500, res = 300)
balloonplot(tb)
dev.off()

saveRDS(sce,file="/mnt/8w/data7/yiyuan/Bone/test1/NEW_MPs/sce_second.final.rds")

sce_subset1 <- sce[, sce$celltype_new %in% c("Macrophages")]
sce1=sce_subset1
saveRDS(sce1,file="/mnt/8w/data7/yiyuan/Bone/test1/NEW_MPs/sce_second.Macrophages.rds")

png(file = "second_celltype_umap1.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce1,
  reduction = "umap", 
  group.by = "celltype", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
   cols = c(
   'SC_M1' = '#ef070f',  
   'SC_M2' = '#FAD02E',  
   'SC_M3' = '#A8DADC',  
   'SC_M4' = '#aba3e7',  
   'SC_M5' = '#D4E157',  
   'SC_M6' = '#F2826C'
))
dev.off()

png(file = "VlnPlot_MPs_celltype.png", width = 2500, height = 5000, res = 300)
VlnPlot(sce, features = cellMarker,
        stack=T,pt.size=0,  group.by = "celltype_new", 
        flip = T,
        add.noise = T)+#横纵轴不标记任何东西
  theme(axis.text.y = element_blank(), #不显示坐标刻度
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'none')
dev.off()



####提取巨噬细胞二级分群
setwd("/mnt/8w/data7/yiyuan/Bone/test1/NEW_MPs/new")

sce1
head(sce1@meta.data)
table(sce1@meta.data$celltype_new,sce1@meta.data$sample)
sce1 <- NormalizeData(sce1, normalization.method =  "LogNormalize",
                     scale.factor = 10000)
#GetAssay(sce,assay = "RNA")
sce1@meta.data %>% head()
sce1 <- FindVariableFeatures(sce1,
                            selection.method = "vst", nfeatures = 2000)
# 中心化，为下一步PCA做准备
sce1 <- ScaleData(sce1, vars.to.regress = c("nCount_RNA", "nFeature_RNA"))
sce1 <- RunPCA(object = sce1, pc.genes = VariableFeatures(sce1))
library(harmony)
sce1 <- sce1%>% harmony::RunHarmony("rawname")

# 设置颜色：指定HA和OA组的颜色
group_colors <- c("HA" = "#1f77b4", "OA" = "#ff7f0e")  # 你可以替换为任何你想要的颜色代码
png(file = "sample_harmony.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce1, 
  group.by = "gname", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
  cols = group_colors  # 设置颜色
)
dev.off()

saveRDS(sce1,file="sce_Macrophages.1.rds")


png(file = "Macrophages_cluster.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce1, 
  group.by = "cluster", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
)
dev.off()

png(file = "ElbowPlot.png", width = 3000, height = 2500, res = 300)
ElbowPlot(sce1, ndims = 40, reduction = "pca")
dev.off()

#看看前20个主成分
pdf(file="PCA.pdf",width=10,height=8)
DimHeatmap(sce1, dims = 1:20, cells = 100, balanced = TRUE,nfeatures = 10)
dev.off()

pct <- sce1[["pca"]]@stdev / sum(sce1[["pca"]]@stdev) * 100 ; 
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)

pc.use
#13
pcs=pc.use

set.seed(123)
sce1 <- RunUMAP(sce1,reduction = "harmony", dims = 1:pcs)
sce1 <- RunTSNE(sce1, dims.use = 1:pcs, reduction = "harmony",do.fast = TRUE)
sce1 <- FindNeighbors(sce1, reduction = "harmony", dims = 1:pcs)
sce1 <- FindClusters(sce1, algorithm = 1 , resolution = c(0.5),cluster.name="Macrophages_0.5")

head(sce1@meta.data)
write.csv(sce1@meta.data, file = "sce_meta.data_Macrophages_0.5.csv", row.names = TRUE)

png(file = "second_Macrophages0.5.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce1,
  reduction = "umap", 
  group.by = "Macrophages_0.5", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
)
dev.off()
############################################################
sce.markers <- FindAllMarkers(object = sce1, only.pos = TRUE, min.pct = 0.25,
                              thresh.use = 0.25)
write.csv(sce.markers,file=paste0('Macrophages_all.markers.csv'))
table(sce1@meta.data$seurat_clusters,sce1@meta.data$sample)

top5 <- sce.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
write.csv(top5,file=paste0('Macrophages_top5.markers.csv'))
pdf(file="all.markers_top5_heatmap.pdf",width=30,height=10)
DoHeatmap(sce1,top5$gene,size=3)
dev.off()
sce1@meta.data %>% head()
##############################################################

SC_M1 = c('S100A8','S100A9','S100A10')
SC_M2 = c('LYVE1','MRC1','FOLR2')
SC_M3 = c('IL1B','CCL3','CCL4','CXCL2','CXCL8')
SC_M4 = c('APOC1','APOE','ACP5')
SC_M5 = c('FCN1','VCAN')
SC_M6 = c('ISG15','IFIT1','IFIT2','IFIT3','IFI6')

cellMarker <- c(SC_M1,SC_M2,SC_M3,SC_M4,SC_M5,SC_M6)
source("/mnt/3w/data2/yiyuan/script/2024.11.work/DotPlot_fixed.R")
library(ggplot2)
library(Seurat)
library(ggsci)
library(cowplot)
######使用修改后的dotplot_fixed画气泡图
pdf(file="Dotplot_Macrophages_0.5.pdf",width=20,height=8)
DotPlot_fixed(sce1,group.by = 'Macrophages_0.5',cols = c('blue','red'),features=unique(cellMarker),cluster.idents = F)+
 scale_fill_gsea()+
 theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 15,
                                   angle = 90,
                                   hjust = 1,
                                   vjust = .5),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
dev.off()

png(file = "VlnPlot_Macrophages_celltype.png", width = 2500, height = 5000, res = 300)
VlnPlot(sce1, features = cellMarker,
        stack=T,pt.size=0,  group.by = "Macrophages_0.5", 
        flip = T,
        add.noise = T)+#横纵轴不标记任何东西
  theme(axis.text.y = element_blank(), #不显示坐标刻度
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'none')
dev.off()


celltype=data.frame(ClusterID=0:7,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(1,3),2]='SC_M1'
celltype[celltype$ClusterID %in% c(2,5),2]='SC_M2'
celltype[celltype$ClusterID %in% c(0),2]='SC_M3'
celltype[celltype$ClusterID %in% c(4),2]='SC_M4'
celltype[celltype$ClusterID %in% c(6),2]='SC_M5'
celltype[celltype$ClusterID %in% c(7),2]='SC_M6'
head(celltype)

sce1@meta.data$celltype_Macrophages = "unknow"
for(i in 1:nrow(celltype)){
  sce1@meta.data[which(sce1@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype_Macrophages'] <- celltype$celltype[i]}
table(sce1@meta.data$celltype_Macrophages)

sce1$celltype_Macrophages <- factor(sce1$celltype_Macrophages,levels=c('SC_M1','SC_M2','SC_M3','SC_M4','SC_M5','SC_M6'))
pdf(file="Dotplot_type_Macrophages.pdf",width=30,height=8)
DotPlot_fixed(sce1,group.by = 'celltype_Macrophages',cols = c('blue','red'),features=unique(cellMarker),cluster.idents = F)+
 scale_fill_gsea()+
 theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 10,
                                   angle = 90,
                                   hjust = 1,
                                   vjust = .5),
        axis.text.y = element_text(size = 10),
        legend.text = element_text(size = 10),
        legend.title = element_text(size = 10))
dev.off()

png(file = "Macrophages_celltype_umap.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce1,
  reduction = "umap", 
  group.by = "celltype_Macrophages", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
   cols = c(
   'SC_M1' = '#ef070f',  
   'SC_M2' = '#FAD02E',  
   'SC_M3' = '#A8DADC',  
   'SC_M4' = '#aba3e7',  
   'SC_M5' = '#D4E157',  
   'SC_M6' = '#6c85f2'

))
dev.off()

cell.Patient <- aggregate(sce1$sample, list(sce1$sample, sce1$gname, sce1$celltype_Macrophages), length)
colnames(cell.Patient) <- c("sample","gname","cellType","number")
cell.Patient$Patient.type <- paste(cell.Patient$sample, cell.Patient$gname, sep = "-")
sce1@meta.data$Patient.type <- paste(sce1@meta.data$sample, sce1@meta.data$gname, sep = "-")
patient <- unique(sce1$Patient.type)
for(i in patient){
    if(i == patient[1]){
        cell.pro.Patient <- subset(cell.Patient, Patient.type == i)
        cell.pro.Patient$proportion <- cell.pro.Patient$number/sum(cell.pro.Patient$number)
    } else {
        single.patient <- subset(cell.Patient, Patient.type == i)
        single.patient$proportion <- single.patient$number/sum(single.patient$number)
        cell.pro.Patient <- rbind(cell.pro.Patient, single.patient)
    }    
}
cell.pro.Patient.OA <- cell.pro.Patient[cell.pro.Patient$gname == "OA",]
cell.pro.Patient.HA <- cell.pro.Patient[cell.pro.Patient$gname == "HA",]
cell.pro.Patient.OA$x_start <- as.numeric(factor(cell.pro.Patient.OA$Patient.type)) - 0.5
cell.pro.Patient.OA$x_end <- as.numeric(factor(cell.pro.Patient.OA$Patient.type)) + 0.5
cell.pro.Patient.HA$x_start <- as.numeric(factor(cell.pro.Patient.HA$Patient.type)) + length(unique(cell.pro.Patient.OA$Patient.type)) - 0.5
cell.pro.Patient.HA$x_end <- as.numeric(factor(cell.pro.Patient.HA$Patient.type)) + length(unique(cell.pro.Patient.OA$Patient.type)) + 0.5
cell.pro.Patient <- rbind(cell.pro.Patient.OA, cell.pro.Patient.HA)
cell.pro.Patient$Patient.type <- factor(cell.pro.Patient$Patient.type, level = c(cell.pro.Patient$Patient.type[grep("OA",cell.pro.Patient$Patient.type)]%>%unique(),
                                               cell.pro.Patient$Patient.type[grep("HA",cell.pro.Patient$Patient.type)]%>%unique()))
library("ggalluvial")
library("rlang")
library(ggplot2)
patient_col = c(
   'SC_M1' = '#ef070f',  
   'SC_M2' = '#FAD02E',  
   'SC_M3' = '#A8DADC',  
   'SC_M4' = '#aba3e7',  
   'SC_M5' = '#D4E157',  
   'SC_M6' = '#6c85f2'
)

# 设置图形的高度和宽度
pdf(file = "proportion_second.pdf", width = 10, height = 12)

# 绘制堆积柱状图
ggplot(data = cell.pro.Patient, mapping = aes(x = Patient.type, fill = cellType, y = proportion)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = patient_col) +
  labs(x = "Samples", y = "Ratio of Cell Types (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_segment(mapping = aes(x = x_start, y = 1, xend = x_end, yend = 1, color = gname), size = 4)
dev.off()


library(tidyr)
library(reshape2)
sce1$sample = factor(sce1$sample,levels = c('OA1','OA2','OA3','HA1','HA2','HA3'))
tb=table(sce1$sample, sce1$celltype_Macrophages)
head(tb)
library (gplots)

png(file = "balloonplot_2_celltype.png", width = 4000, height = 2500, res = 300)
balloonplot(tb)
dev.off()


saveRDS(sce1,file="sce_Macrophages.final.rds")

png(file = "VlnPlot_Macrophages_celltype2.png", width = 2500, height = 5000, res = 300)
VlnPlot(sce1, features = cellMarker,
        stack=T,pt.size=0,  group.by = "celltype_Macrophages", 
        flip = T,
        add.noise = T)+#横纵轴不标记任何东西
  theme(axis.text.y = element_blank(), #不显示坐标刻度
        axis.ticks.y = element_blank(),
        axis.title = element_blank(),
        axis.text.x = element_text(colour = 'black',size = 10,angle = 90),
        legend.position = 'none')
dev.off()
