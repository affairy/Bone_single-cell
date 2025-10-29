###将软骨细胞与滑膜成纤维和巨噬细胞进行合并通讯分析
##将seuratv5转seuratv4格式

library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
options(stringsAsFactors = FALSE)
##先分别读取软骨，滑膜成纤维，滑膜巨噬细胞区分HA和OA
setwd("/mnt/8w/data7/yiyuan/Bone/test1/Fibroblasts/new/monocle")
##软骨已经区分存放于/mnt/8w/data7/yiyuan/Bone/test/Ch/Re_subtype/monocle
###滑膜成纤维
sce <- readRDS("/mnt/8w/data7/yiyuan/Bone/test1/Fibroblasts/new/sce_Fb_resubtype.final.rds")

###分别提取HA和OA组
HA <- sce[, sce$gname %in% c("HA")]
HA
OA <- sce[, sce$gname %in% c("OA")]
OA
saveRDS(HA,file="sce_Sy_Fb_resubtype.final_HA.rds")
saveRDS(OA,file="sce_Sy_Fb_resubtype.final_OA.rds")

#####滑膜巨噬细胞
setwd("/mnt/8w/data7/yiyuan/Bone/test1/NEW_MPs/new/monocle")
sce1 <- readRDS("/mnt/8w/data7/yiyuan/Bone/test1/NEW_MPs/new/sce_Macrophages.final.rds")
###分别提取HA和OA组
HA <- sce1[, sce1$gname %in% c("HA")]
HA
OA <- sce1[, sce1$gname %in% c("OA")]
OA
saveRDS(HA,file="sce_Sy_Mo_resubtype.final_HA.rds")
saveRDS(OA,file="sce_Sy_Mo_resubtype.final_OA.rds")

#######根据分别提取的结果将三种细胞类型的细胞进行合并，区分HA组和OA组
##分别读取HA的三种细胞类型
setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat")
HA_Ch <- readRDS("/mnt/8w/data7/yiyuan/Bone/test/Ch/Re_subtype/monocle/sce_Ch_resubtype.final_HA.rds")
HA_Sy_Fb <- readRDS("/mnt/8w/data7/yiyuan/Bone/test1/Fibroblasts/new/monocle/sce_Sy_Fb_resubtype.final_HA.rds")
HA_Sy_Mo <- readRDS("/mnt/8w/data7/yiyuan/Bone/test1/NEW_MPs/new/monocle/sce_Sy_Mo_resubtype.final_HA.rds")

HA_Ch <- RenameCells(HA_Ch, new.names = HA_Ch@meta.data$cellid)
HA_Sy_Fb <- RenameCells(HA_Sy_Fb, new.names = HA_Sy_Fb@meta.data$cellid)
HA_Sy_Mo <- RenameCells(HA_Sy_Mo, new.names = HA_Sy_Mo@meta.data$cellid)

HA_intermediate <- merge(HA_Sy_Mo, HA_Sy_Fb, project = "HA")
HA_final <- merge(HA_intermediate, HA_Sy_Mo, project = "HA")
HA_final <- merge(HA_Ch, y = c(HA_Sy_Fb, HA_Sy_Mo), project = "HA")
HA_final

# merge all datasets, adding a cell ID to make sure cell names are unique
combined <- merge(
  x = HA_Ch,
  y = list(pbmc1k, pbmc5k, pbmc10k),
  add.cell.ids = c("500", "1k", "5k", "10k")
)

mat_mo <- GetAssayData(HA_Sy_Mo, assay = "origin.counts", slot = "counts")
mat_fb <- GetAssayData(HA_Sy_Fb, assay = "origin.counts", slot = "counts")

combined_mat <- cbind(mat_mo, mat_fb)

meta_mo <- HA_Sy_Mo@meta.data
meta_fb <- HA_Sy_Fb@meta.data

combined_meta <- rbind(meta_mo, meta_fb)

HA_intermediate <- CreateSeuratObject(counts = combined_mat, meta.data = combined_meta, project = "HA")

# 提取 HA_Ch 的表达矩阵和元数据
mat_ch <- GetAssayData(HA_Ch, assay = "origin.counts", slot = "counts")
meta_ch <- HA_Ch@meta.data

# 合并表达矩阵
combined_mat <- cbind(mat_mo, mat_fb, mat_ch)

# 合并元数据
combined_meta <- rbind(meta_mo, meta_fb, meta_ch)

# 创建新的 Seurat 对象
HA_intermediate <- CreateSeuratObject(counts = combined_mat, meta.data = combined_meta, project = "HA")

HA_Ch <- RenameCells(HA_Ch, new.names = paste0("Ch_", HA_Ch@meta.data$cellid))
HA_Sy_Fb <- RenameCells(HA_Sy_Fb, new.names = paste0("Fb_", HA_Sy_Fb@meta.data$cellid))
HA_Sy_Mo <- RenameCells(HA_Sy_Mo, new.names = paste0("Mo_", HA_Sy_Mo@meta.data$cellid))

mat_mo <- GetAssayData(HA_Sy_Mo, assay = "origin.counts", slot = "counts")
meta_mo <- HA_Sy_Mo@meta.data

mat_fb <- GetAssayData(HA_Sy_Fb, assay = "origin.counts", slot = "counts")
meta_fb <- HA_Sy_Fb@meta.data

mat_ch <- GetAssayData(HA_Ch, assay = "origin.counts", slot = "counts")
meta_ch <- HA_Ch@meta.data

common_genes <- intersect(rownames(mat_mo), intersect(rownames(mat_fb), rownames(mat_ch)))
mat_mo <- mat_mo[common_genes, , drop = FALSE]
mat_fb <- mat_fb[common_genes, , drop = FALSE]
mat_ch <- mat_ch[common_genes, , drop = FALSE]

common_meta <- intersect(intersect(colnames(meta_mo), colnames(meta_fb)), colnames(meta_ch))

meta_mo <- meta_mo[, common_meta, drop = FALSE]
meta_fb <- meta_fb[, common_meta, drop = FALSE]
meta_ch <- meta_ch[, common_meta, drop = FALSE]

# 合并表达矩阵
combined_mat <- cbind(mat_mo, mat_fb, mat_ch)

# 合并元数据
combined_meta <- rbind(meta_mo, meta_fb, meta_ch)
HA_intermediate <- CreateSeuratObject(counts = combined_mat, meta.data = combined_meta, project = "HA")
# 检查对象的基本信息
HA_intermediate

saveRDS(HA_intermediate,file="sce_merged.final.rds")
HA_intermediate[["RNA4"]] <- as(object = HA_intermediate[["RNA"]], Class = "Assay")

DefaultAssay(HA_intermediate) <- "RNA4"

saveRDS(HA_intermediate,file="sce_merged_final_seuratV4.rds")

###3w cellchat 环境
library(Seurat)
library(tidyverse)
library(CellChat)
library(NMF)
library(ggalluvial)
library(patchwork)
library(ggplot2)
library(svglite)
options(stringsAsFactors = FALSE)
setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat")
####转换成seuratV4格式进行运行 因为cellchat安装的版本为4.0的seurat
sce <- readRDS("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/sce_merged_final_seuratV4.rds")
head(sce@meta.data)

cellchat <- createCellChat(object=sce,meta = sce@meta.data,group.by = "celltype")

cellchat
summary(cellchat)

groupSize <- as.numeric(table(cellchat@idents))  

groupSize

CellChatDB <- CellChatDB.human
colnames(CellChatDB$interaction)
# [1] "interaction_name"   "pathway_name"       "ligand"            
# [4] "receptor"           "agonist"            "antagonist"        
# [7] "co_A_receptor"      "co_I_receptor"      "evidence"          
#[10] "annotation"         "interaction_name_2"
CellChatDB$interaction[1:4,1:4]

showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)
#[1] "Secreted Signaling" "ECM-Receptor"       "Cell-Cell Contact" 
#选择"Secreted Signaling"进行后续细胞互作分析
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
# use all CellChatDB for cell-cell communication analysis
# CellChatDB.use <- CellChatDB # simply use the default CellChatDB
# set the used database in the object
cellchat@DB <- CellChatDB.use # set the used database in the object
cellchat <- subsetData(cellchat)
#future::plan("multiprocess", workers = 4)
#相当于Seurat的FindMarkers，找每个细胞群中高表达的配体受体
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat) #Identify over-expressed ligand-receptor interactions (pairs) within the used CellChatDB
#上一步运行的结果储存在cellchat@LR$LRsig
cellchat <- projectData(cellchat, PPI.human) 
#找到配体受体关系后，projectData将配体受体对的表达值投射到PPI上，来对@data.signaling中的表达值进行校正。结果保存在@data.project
cellchat@data.project[1:4,1:4]

levels(cellchat@idents)
# [1] "SC_M1" "SC_M2" "SC_M3" "SC_M4" "SC_M5" "SC_M6" "Fb1"   "Fb2"   "Fb3"  
#[10] "Fb4"   "Fb5"   "Fb6"   "Fb7"   "Ch1"   "Ch3"   "Ch2"   "Ch6"   "Ch7"  
#[19] "Ch5"   "Ch4"   "Ch8"  
# 从因子变量中移除未使用的水平
cellchat@idents <- droplevels(cellchat@idents)

#根据表达值推测细胞互作的概率（cellphonedb是用平均表达值代表互作强度）。
cellchat <- computeCommunProb(cellchat, raw.use = FALSE, population.size = TRUE) #如果不想用上一步PPI矫正的结果，raw.use = TRUE即可。
# Filter out the cell-cell communication if there are only few number of cells in certain cell groups
cellchat <- filterCommunication(cellchat, min.cells = 10)
df.net <- subsetCommunication(cellchat)
write.csv(df.net, "net_lr.csv")

cellchat <- computeCommunProbPathway(cellchat)
df.netp <- subsetCommunication(cellchat, slot.name = "netP")
write.csv(df.netp, "net_pathway.csv")

cellchat <- aggregateNet(cellchat)
# 设置绘图为 1 行 2 列的布局
par(mfrow = c(1, 2), xpd = TRUE)
groupSize=21
# 打开 PDF 图形设备，指定保存的文件名和路径
pdf("TIL_net_number_strength.pdf", width = 10, height = 6)

# 绘制第一个图：Number of interactions
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge = FALSE, title.name = "Number of interactions")

# 绘制第二个图：Interaction weights/strength
netVisual_circle(cellchat@net$weight, vertex.weight = groupSize, weight.scale = TRUE, 
                 label.edge = FALSE, title.name = "Interaction weights/strength")

# 关闭 PDF 图形设备，保存文件
dev.off()

mat <- cellchat@net$count

# 设置图形布局为 3x3
par(mfrow = c(3, 3), xpd = TRUE)

# 打开 PDF 图形设备，指定保存路径和文件名
pdf("TIL_net_number_individual.pdf", width = 8, height = 8)

# 循环绘制每一行的网络图
for (i in 1:nrow(mat)) {
  # 创建一个新的矩阵，只有第 i 行的值
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  
  # 绘制单个网络图
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = TRUE, 
                   arrow.width = 0.2, arrow.size = 0.1, edge.weight.max = max(mat), 
                   title.name = rownames(mat)[i])
}

# 关闭图形设备并保存到 PDF 文件
dev.off()

mat <- cellchat@net$weight
par(mfrow = c(3,3), xpd=T)
pdf("TIL_net_strength_individual.pdf", width = 8, height = 8)
for (i in 1:nrow(mat)) {
  mat2 <- matrix(0, nrow = nrow(mat), ncol = ncol(mat), dimnames = dimnames(mat))
  mat2[i, ] <- mat[i, ]
  netVisual_circle(mat2, vertex.weight = groupSize, weight.scale = T, arrow.width = 0.2,
                   arrow.size = 0.1, edge.weight.max = max(mat), title.name = rownames(mat)[i])
}
dev.off()

save(cellchat, file = "cellchat_1.RData")

cellchat@netP$pathways
# [1] "TGFb"       "ANGPTL"     "MK"         "MIF"        "CXCL"      
# [6] "FGF"        "CCL"        "VISFATIN"   "ANNEXIN"    "PTN"       
#[11] "BMP10"      "RANKL"      "ACTIVIN"    "GDF"        "BMP"       
#[16] "GRN"        "IGF"        "GAS"        "SEMA3"      "PROS"      
#[21] "GALECTIN"   "LIFR"       "EGF"        "COMPLEMENT" "PDGF"      
#[26] "OSM"        "IL6"        "IL1"        "IL4"        "VEGF"      
#[31] "PERIOSTIN"  "CSF"        "CHEMERIN"   "IFN-II"     "TWEAK"     
#[36] "LT"         "PTH"        "PSAP"       "IL16"       "CD40"      
#[41] "TRAIL"      "TNF"        "IL2"        "LIGHT"      "BAFF"      
#[46] "CD30"       "BTLA"       "CD137"      "OPIOID"     "KIT"       
#[51] "FASLG"      "NPR2"       "NT"         "GH"         "EPO"       
#[56] "PRL"        "NGF"        "PARs"       "CSF3"       "FLT3"      
#[61] "FSH"  
pathways.show <- c("TGFb")
levels(cellchat@idents)
vertex.receiver = c(14,15,16,17,18,19,20,21,9,8,7)
setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/TGFb")
pdf("TIL_TGFb_hierarchy.pdf", width = 20, height = 10)
netVisual_aggregate(cellchat, signaling = pathways.show,  vertex.receiver = vertex.receiver, layout = "hierarchy")
dev.off()
par(mfrow=c(1,1))
pdf("TIL_TGFb_circle.pdf", width = 20, height = 10)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "circle")
dev.off()
par(mfrow=c(1,1))
pdf("TIL_TGFb_chord.pdf", width = 10, height = 10)
netVisual_aggregate(cellchat, signaling = pathways.show, layout = "chord")
dev.off()
par(mfrow=c(1,1))
pdf("TIL_TGFb_heatmap.pdf", width = 10, height = 10)
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()
pdf("TIL_TGFb_LR_contribution.pdf", width = 10, height = 10)
netAnalysis_contribution(cellchat, signaling = pathways.show)
dev.off()

# 获取所有信号通路
pathways.show.all <- cellchat@netP$pathways

# 查看细胞类型（以便选择目标细胞类型）
levels(cellchat@idents)

# 设置目标细胞类型的索引（选择适合的细胞类型）


# 创建文件夹保存图形
if (!dir.exists("all_pathways_com_circle")) {
  dir.create("all_pathways_com_circle")
}
setwd("all_pathways_com_circle")

# 循环遍历所有信号通路进行绘制
for (i in 1:length(pathways.show.all)) {

  # 绘制并保存每个信号通路的不同图
  # Circle layout
  pdf(paste0(pathways.show.all[i], "_circle.pdf"), width = 20, height = 10)
  netVisual_aggregate(cellchat, signaling = pathways.show.all[i], layout = "circle")
  dev.off()

  # Chord layout
  pdf(paste0(pathways.show.all[i], "_chord.pdf"), width = 10, height = 10)
  netVisual_aggregate(cellchat, signaling = pathways.show.all[i], layout = "chord")
  dev.off()

}

# 获取所有信号通路
pathways.show.all <- cellchat@netP$pathways

# 循环遍历所有信号通路进行绘制
for (i in 1:length(pathways.show.all)) {
  
  # 动态获取当前信号通路的名称
  pathway <- pathways.show.all[i]
  
  # 设置目标细胞类型的索引（根据实际需要选择细胞类型）
  vertex.receiver <- c(14, 15, 16, 17, 18, 19, 20, 21, 9, 8, 7)
  
  # 设置工作目录
  setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/TGFb/all_pathways_com_circle/all_pathways_com_circle/test")
  
  # 绘制并保存Hierarchy布局
  pdf(paste0("TIL_", pathway, "_hierarchy.pdf"), width = 20, height = 10)
  netVisual_aggregate(cellchat, signaling = pathway, vertex.receiver = vertex.receiver, layout = "hierarchy")
  dev.off()
  
  # 如果需要，您可以添加更多的图形绘制，例：
  # 绘制Circle布局
  # pdf(paste0("TIL_", pathway, "_circle.pdf"), width = 20, height = 10)
  # netVisual_aggregate(cellchat, signaling = pathway, layout = "circle")
  # dev.off()
}

# [1] "TGFb"       "ANGPTL"     "MK"         "MIF"        "CXCL"      
# [6] "FGF"        "CCL"        "VISFATIN"   "ANNEXIN"    "PTN"       
#[11] "BMP10"      "RANKL"      "ACTIVIN"    "GDF"        "BMP"       
#[16] "GRN"        "IGF"        "GAS"        "SEMA3"      "PROS"      
#[21] "GALECTIN"   "LIFR"       "EGF"        "COMPLEMENT" "PDGF"      
#[26] "OSM"        "IL6"        "IL1"        "IL4"        "VEGF"      
#[31] "PERIOSTIN"  "CSF"        "CHEMERIN"   "IFN-II"     "TWEAK"     
#[36] "LT"         "PTH"        "PSAP"       "IL16"       "CD40"      
#[41] "TRAIL"      "TNF"        "IL2"        "LIGHT"      "BAFF"      
#[46] "CD30"       "BTLA"       "CD137"      "OPIOID"     "KIT"       
#[51] "FASLG"      "NPR2"       "NT"         "GH"         "EPO"       
#[56] "PRL"        "NGF"        "PARs"       "CSF3"       "FLT3"      
#[61] "FSH"  
pathways.show <- c("CCL")
par(mfrow=c(1,1))
pdf("TIL_CCL_heatmap.pdf", width = 10, height = 10)
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()


levels(cellchat@idents)
#
# [1] "SC_M1" "SC_M2" "SC_M3" "SC_M4" "SC_M5" "SC_M6" "Fb1"   "Fb2"   "Fb3"  
#[10] "Fb4"   "Fb5"   "Fb6"   "Fb7"   "Ch1"   "Ch3"   "Ch2"   "Ch6"   "Ch7"  
#[19] "Ch5"   "Ch4"   "Ch8"  
# show all the significant interactions (L-R pairs)
#需要指定受体细胞和配体细胞
setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/dot")
p = netVisual_bubble(cellchat, sources.use = c(14), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch1.pdf", p, width = 10, height = 30) 
# save as TIL/Mye_Lymph_bubble.pdf
p = netVisual_bubble(cellchat, sources.use = c(15), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch3.pdf", p, width = 10, height = 30) 
# save as TIL/Mye_Lymph_bubble.pdf
p = netVisual_bubble(cellchat, sources.use = c(16), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch2.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(17), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch6.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(18), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch7.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(19), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch5.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(20), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch4.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(21), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch8.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(1), 
                     targets.use = c(21,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_SC_M1.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(2), 
                     targets.use = c(21,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_SC_M2.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(3), 
                     targets.use = c(21,1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_SC_M3.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(4), 
                     targets.use = c(21,1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_SC_M4.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(5), 
                     targets.use = c(21,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_SC_M5.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(6), 
                     targets.use = c(21,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_SC_M6.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(7), 
                     targets.use = c(21,1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb1.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(8), 
                     targets.use = c(21,1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb2.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(9), 
                     targets.use = c(21,1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb3.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(10), 
                     targets.use = c(21,1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb4.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(11), 
                     targets.use = c(21,1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb5.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(12), 
                     targets.use = c(21,1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb6.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(13), 
                     targets.use = c(21,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb7.pdf", p, width = 10, height = 30) 

setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/singlerole")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 8,width = 15, height = 25)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 8,width = 15, height = 25)
pdf("TIL_SNA_SignalingPattern.pdf", width = 18, height = 26)
ht1 + ht2
dev.off()  

save(cellchat, file = "cellchat_HA_final.RData")

setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/dot")
load("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/singlerole/cellchat_HA_final.RData")

p = netVisual_bubble(cellchat, targets.use = c(14), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch1.1.pdf", p, width = 10, height = 30) 



p = netVisual_bubble(cellchat, targets.use = c(15), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch3.1.pdf", p, width = 10, height = 30) 


# save as TIL/Mye_Lymph_bubble.pdf
p = netVisual_bubble(cellchat, targets.use = c(16), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch2.1.pdf", p, width = 10, height = 30) 



p = netVisual_bubble(cellchat, targets.use = c(17), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch6.1.pdf", p, width = 10, height = 30) 


p = netVisual_bubble(cellchat, targets.use = c(18), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch7.1.pdf", p, width = 10, height = 30) 


p = netVisual_bubble(cellchat, targets.use = c(19), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch5.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(20), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch4.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(21), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch8.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(1), 
                     sources.use = c(21,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_SC_M1.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(2), 
                     sources.use = c(21,1,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_SC_M2.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(3), 
                     sources.use = c(21,1,2,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_SC_M3.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(4), 
                     sources.use = c(21,1,2,3,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_SC_M4.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(5), 
                     sources.use = c(21,1,2,3,4,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_SC_M5.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(6), 
                     sources.use = c(21,1,2,3,4,5,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_SC_M6.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(7), 
                     sources.use = c(21,1,2,3,4,5,6,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb1.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(8), 
                     sources.use = c(21,1,2,3,4,5,6,7,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb2.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(9), 
                     sources.use = c(21,1,2,3,4,5,6,7,8,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb3.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(10), 
                     sources.use = c(21,1,2,3,4,5,6,7,8,9,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb4.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(11), 
                     sources.use = c(21,1,2,3,4,5,6,7,8,9,10,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb5.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(12), 
                     sources.use = c(21,1,2,3,4,5,6,7,8,9,10,11,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb6.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(13), 
                     sources.use = c(21,1,2,3,4,5,6,7,8,9,10,11,12,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb7.1.pdf", p, width = 10, height = 30) 



####细胞通讯图修改
setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/HA_NEW")
load("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/singlerole/cellchat_HA_final.RData")
levels(cellchat@idents)
p = netVisual_bubble(cellchat, sources.use = c(7,10), targets.use = c(14,16,15,20,19,17,18,21), signaling = c("MK"), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb1_4.pdf", p, width = 10, height = 10)

p1 = netVisual_bubble(cellchat, sources.use = c(14,16,19), targets.use = c(7,8,9,10,11,12,13), signaling = c("VISFATIN"), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch1_2_5.pdf", p1, width = 10, height = 5)

pathways.show <- c("TGFb","ANGPTL","MK","MIF","CXCL","FGF","CCL","VISFATIN","ANNEXIN","PTN","BMP10","RANKL","ACTIVIN","GDF","BMP")
setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/HA_NEW")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", signaling = pathways.show, font.size = 8,width = 15, height = 10)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", signaling = pathways.show, font.size = 8,width = 15, height = 10)
pdf("TIL_SNA_SignalingPattern_HA.pdf", width = 18, height = 10)
ht1 + ht2
dev.off()  
pdf("TIL_SNA_SignalingPattern_HA_outing.pdf", width = 18, height = 10)
ht1
dev.off()  





