#软骨导出矩阵
library(Seurat)
setwd("/mnt/data4/yiyuan/Bone/scRNA/scenic/Ch")
#remotes::install_github("satijalab/seurat-data")
sce = readRDS("/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/sce_Ch_resubtype.final.rds")
head(t(as.matrix(GetAssayData(sce, assay = "origin.counts", layer = "counts")))[1:3, 1:3])

head(t(as.matrix(sce@assays$origin.counts@counts))[1:3,1:3])
colnames(t(as.matrix(GetAssayData(sce, assay = "origin.counts", layer = "counts"))))

# output matrix
write.csv(t(as.matrix(GetAssayData(sce, assay = "origin.counts", layer = "counts"))),file="scenic.data.csv")


####成纤维导出矩阵
setwd("/mnt/data4/yiyuan/Bone/scRNA/scenic/Fb")
library(Seurat)
sce = readRDS("/mnt/data7/yiyuan/Bone/test1/Fibroblasts/new/sce_Fb_resubtype.final.rds")
head(t(as.matrix(GetAssayData(sce, assay = "origin.counts", layer = "counts")))[1:3, 1:3])

#head(t(as.matrix(sce@assays$origin.counts@counts))[1:3,1:3])
colnames(t(as.matrix(GetAssayData(sce, assay = "origin.counts", layer = "counts"))))

# output matrix
write.csv(t(as.matrix(GetAssayData(sce, assay = "origin.counts", layer = "counts"))),file="scenic.data.csv")

