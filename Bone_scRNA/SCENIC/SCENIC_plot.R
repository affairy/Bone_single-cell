###micromamba activate SCENIC
# ==============================================================================
# plot in R
# ==============================================================================
library(Seurat)
library(SCopeLoomR)
library(AUCell)
library(SCENIC)
library(dplyr)
library(KernSmooth)
library(RColorBrewer)
library(plotly)
library(BiocParallel)
library(grid)
library(ComplexHeatmap)
library(data.table)
library(scRNAseq)
library(patchwork)
library(ggplot2)
library(stringr)
library(circlize)

setwd("/mnt/data4/yiyuan/Bone/scRNA/scenic/Ch/output/plots")

# load data
loom <- open_loom('/mnt/data4/yiyuan/Bone/scRNA/scenic/Ch/output/aucell.loom')

regulons_incidMat <- SCopeLoomR::get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4]
##regulonsToGeneLists 函数返回每个转录因子调控的靶基因列表：
regulons <- SCENIC::regulonsToGeneLists(regulons_incidMat)
head(regulons,3)
#get_regulons_AUC 返回计算好的转录因子的活性分数矩阵，行为转录因子，列为细胞
regulonAUC <- SCopeLoomR::get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAUC
#get_regulon_thresholds 返回每个转录因子的 AUC 阈值：
regulonAucThresholds <- SCopeLoomR::get_regulon_thresholds(loom)
head(regulonAucThresholds,3)
embeddings <- SCopeLoomR::get_embeddings(loom)
embeddings
# list()
close_loom(loom)

###处理原始scRNA-Seq对象
seurat.data=readRDS("/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/sce_Ch_resubtype.final.rds")
#把每个转录因子的活性合并到 metadata 里面：
# precoess
sub_regulonAUC <- regulonAUC[,match(colnames(seurat.data),colnames(regulonAUC))]

# check
identical(colnames(sub_regulonAUC), colnames(seurat.data))

# combine
seurat.data@meta.data = cbind(seurat.data@meta.data ,t(assay(sub_regulonAUC)))
head(seurat.data@meta.data)
seurat.data$Cgroup <- paste0(seurat.data$gname, "_", seurat.data$celltype)


###小提琴图
#VlnPlot(seurat.data, features = c("TCF4(+)","GATA3(+)"),pt.size = 0)
#FeaturePlot(object = seurat.data,features = c("TCF4(+)","GATA3(+)"))
#RidgePlot(seurat.data, features = c("TCF4(+)","GATA3(+)") , ncol = 2)

####热图可视化
# ==============================================================================
# AVG expression
# ==============================================================================
cellClusters <- data.frame(row.names = colnames(seurat.data),
                           seurat_clusters = as.character(seurat.data$gname)) |>
  dplyr::mutate(seurat_clusters = ifelse(is.na(seurat_clusters),"unkown",seurat_clusters))

cellsPerGroup <- split(rownames(cellClusters),cellClusters$seurat_clusters)

# remove extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),]

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells)
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

# scale
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),center = T, scale=T))

# heatmap
png("Heatmap20.png", width = 15, height = 8, units = "in", res = 300)
Heatmap(t(regulonActivity_byGroup_Scaled[1:20, ]))
dev.off()


png("Heatmap_gname.png", width = 30, height = 8, units = "in", res = 300)
Heatmap(t(regulonActivity_byGroup_Scaled))
dev.off()


######Fb
setwd("/mnt/data4/yiyuan/Bone/scRNA/scenic/Fb/output/plots")

# load data
loom <- open_loom('/mnt/data4/yiyuan/Bone/scRNA/scenic/Fb/output/aucell.loom')

regulons_incidMat <- SCopeLoomR::get_regulons(loom, column.attr.name="Regulons")
regulons_incidMat[1:4,1:4]
##regulonsToGeneLists 函数返回每个转录因子调控的靶基因列表：
regulons <- SCENIC::regulonsToGeneLists(regulons_incidMat)
head(regulons,3)
## 转换为矩阵形式，按 regulon 分列
#regulons_matrix <- sapply(regulons, function(x) paste(x, collapse = ","))
#regulons_matrix <- as.data.frame(t(regulons_matrix))
#
## 保存为 CSV 文件
#write.csv(regulons_matrix, file = "regulons_matrix.csv", row.names = FALSE)

#get_regulons_AUC 返回计算好的转录因子的活性分数矩阵，行为转录因子，列为细胞
regulonAUC <- SCopeLoomR::get_regulons_AUC(loom,column.attr.name='RegulonsAUC')
regulonAUC
#get_regulon_thresholds 返回每个转录因子的 AUC 阈值：
regulonAucThresholds <- SCopeLoomR::get_regulon_thresholds(loom)
head(regulonAucThresholds,3)
embeddings <- SCopeLoomR::get_embeddings(loom)
embeddings
# list()
close_loom(loom)

###处理原始scRNA-Seq对象
seurat.data=readRDS("/mnt/data7/yiyuan/Bone/test1/Fibroblasts/new/sce_Fb_resubtype.final.rds")
head(seurat.data@meta.data)
colnames(seurat.data)=seurat.data$cellid
#把每个转录因子的活性合并到 metadata 里面：
# precoess
sub_regulonAUC <- regulonAUC[,match(colnames(seurat.data),colnames(regulonAUC))]

# check
identical(colnames(sub_regulonAUC), colnames(seurat.data))

# combine
seurat.data@meta.data = cbind(seurat.data@meta.data ,t(assay(sub_regulonAUC)))

head(seurat.data@meta.data)
seurat.data$Cgroup <- paste0(seurat.data$gname, "_", seurat.data$celltype)
###小提琴图
#VlnPlot(seurat.data, features = c("TCF4(+)","GATA3(+)"),pt.size = 0)
#FeaturePlot(object = seurat.data,features = c("TCF4(+)","GATA3(+)"))
#RidgePlot(seurat.data, features = c("TCF4(+)","GATA3(+)") , ncol = 2)

####热图可视化
# ==============================================================================
# AVG expression
# ==============================================================================
cellClusters <- data.frame(row.names = colnames(seurat.data),
                           seurat_clusters = as.character(seurat.data$gname)) |>
  dplyr::mutate(seurat_clusters = ifelse(is.na(seurat_clusters),"unkown",seurat_clusters))

cellsPerGroup <- split(rownames(cellClusters),cellClusters$seurat_clusters)

# remove extened regulons
sub_regulonAUC <- sub_regulonAUC[onlyNonDuplicatedExtended(rownames(sub_regulonAUC)),]

# Calculate average expression:
regulonActivity_byGroup <- sapply(cellsPerGroup,
                                  function(cells)
                                    rowMeans(getAUC(sub_regulonAUC)[,cells]))

# scale
regulonActivity_byGroup_Scaled <- t(scale(t(regulonActivity_byGroup),center = T, scale=T))

# heatmap
png("Heatmap20.png", width = 15, height = 8, units = "in", res = 300)
Heatmap(t(regulonActivity_byGroup_Scaled[1:20, ]))
dev.off()


png("Heatmap_gname.png", width = 30, height = 8, units = "in", res = 300)
Heatmap(t(regulonActivity_byGroup_Scaled))
dev.off()