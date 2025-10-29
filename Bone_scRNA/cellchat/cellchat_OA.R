##分别读取OA的三种细胞类型
library(Seurat)
library(tidyverse)
setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/OA")
OA_Ch <- readRDS("/mnt/8w/data7/yiyuan/Bone/test/Ch/Re_subtype/monocle/sce_Ch_resubtype.final_OA.rds")
OA_Sy_Fb <- readRDS("/mnt/8w/data7/yiyuan/Bone/test1/Fibroblasts/new/monocle/sce_Sy_Fb_resubtype.final_OA.rds")
OA_Sy_Mo <- readRDS("/mnt/8w/data7/yiyuan/Bone/test1/NEW_MPs/new/monocle/sce_Sy_Mo_resubtype.final_OA.rds")

# 为每个数据集的细胞添加标识
OA_Ch <- RenameCells(OA_Ch, new.names = paste0("Ch_", OA_Ch@meta.data$cellid))
OA_Sy_Fb <- RenameCells(OA_Sy_Fb, new.names = paste0("Fb_", OA_Sy_Fb@meta.data$cellid))
OA_Sy_Mo <- RenameCells(OA_Sy_Mo, new.names = paste0("Mo_", OA_Sy_Mo@meta.data$cellid))

OA_Sy_Mo$celltype=OA_Sy_Mo$celltype_Macrophages
# 提取 OA_Sy_Mo 的表达矩阵和元数据
mat_mo <- GetAssayData(OA_Sy_Mo, assay = "origin.counts", slot = "counts")
meta_mo <- OA_Sy_Mo@meta.data

# 提取 OA_Sy_Fb 的表达矩阵和元数据
mat_fb <- GetAssayData(OA_Sy_Fb, assay = "origin.counts", slot = "counts")
meta_fb <- OA_Sy_Fb@meta.data

# 提取 OA_Ch 的表达矩阵和元数据
mat_ch <- GetAssayData(OA_Ch, assay = "origin.counts", slot = "counts")
meta_ch <- OA_Ch@meta.data


common_genes <- intersect(rownames(mat_mo), intersect(rownames(mat_fb), rownames(mat_ch)))
mat_mo <- mat_mo[common_genes, , drop = FALSE]
mat_fb <- mat_fb[common_genes, , drop = FALSE]
mat_ch <- mat_ch[common_genes, , drop = FALSE]

common_meta <- intersect(intersect(colnames(meta_mo), colnames(meta_fb)), colnames(meta_ch))


meta_mo <- meta_mo[, common_meta, drop = FALSE]
meta_fb <- meta_fb[, common_meta, drop = FALSE]
meta_ch <- meta_ch[, common_meta, drop = FALSE]

combined_mat <- cbind(mat_mo, mat_fb, mat_ch)

combined_meta <- rbind(meta_mo, meta_fb, meta_ch)
OA_intermediate <- CreateSeuratObject(counts = combined_mat, meta.data = combined_meta, project = "OA")
# 检查对象的基本信息
OA_intermediate

saveRDS(OA_intermediate,file="sce_merged.final.rds")
OA_intermediate[["RNA4"]] <- as(object = OA_intermediate[["RNA"]], Class = "Assay")

DefaultAssay(OA_intermediate) <- "RNA4"

saveRDS(OA_intermediate,file="sce_merged_final_seuratV4.rds")


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
setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/OA")
####转换成seuratV4格式进行运行 因为cellchat安装的版本为4.0的seurat
sce <- readRDS("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/OA/sce_merged_final_seuratV4.rds")
head(sce@meta.data)

cellchat <- createCellChat(object=sce,meta = sce@meta.data,group.by = "celltype")
#cellchat <- createCellChat(ec@assays$RNA@data, meta = ec@meta.data, group.by = "cell_type")
cellchat
#An object of class CellChat created from a single dataset 
# 32959 genes.
# 57769 cells. 
#CellChat analysis of single cell RNA-seq data!
summary(cellchat)
#  Length    Class     Mode 
#       1 CellChat       S4 
#str(cellchat)
#levels(cellchat@idents)
#cellchat <- setIdent(cellchat, ident.use = "cell_type")
groupSize <- as.numeric(table(cellchat@idents))  
#查看每个cluster有多少个细胞，后面画图的时候需要用到这个值
groupSize
# [1] 2969 1548 1660  666  475  236 5499 6302 6048 9297 2785  985  350 1146 3475
#[16] 4988 3890  619 2394 1636  801
CellChatDB <- CellChatDB.human
colnames(CellChatDB$interaction)
# [1] "interaction_name"   "pathway_name"       "ligand"            
# [4] "receptor"           "agonist"            "antagonist"        
# [7] "co_A_receptor"      "co_I_receptor"      "evidence"          
#[10] "annotation"         "interaction_name_2"
CellChatDB$interaction[1:4,1:4]
#                       interaction_name pathway_name ligand      receptor
#TGFB1_TGFBR1_TGFBR2 TGFB1_TGFBR1_TGFBR2         TGFb  TGFB1     TGFbR1_R2
#TGFB2_TGFBR1_TGFBR2 TGFB2_TGFBR1_TGFBR2         TGFb  TGFB2     TGFbR1_R2
#TGFB3_TGFBR1_TGFBR2 TGFB3_TGFBR1_TGFBR2         TGFb  TGFB3     TGFbR1_R2
#TGFB1_ACVR1B_TGFBR2 TGFB1_ACVR1B_TGFBR2         TGFb  TGFB1 ACVR1B_TGFbR2

showDatabaseCategory(CellChatDB)

unique(CellChatDB$interaction$annotation)#查看可以选择的侧面，也就是上图左中的三种
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
#         Mo_OA1_ACTGCCTAGACACCAACGAACAAGTGG Mo_OA1_TGCCTGATCAGGCTGTTGAACAAGTGG
#TNFRSF18                          0.0000000                        0.000000000
#TNFRSF4                           0.0000000                        0.000000000
#TNFRSF14                          0.5744968                        0.005147793
#TNFRSF25                          0.0000000                        0.000000000
#         Mo_OA1_CGAGATAGTCCAGAAGATAACAAGTGG Mo_OA1_ACGAAGCTCCTGTGGTATAACAAGTGG
#TNFRSF18                        0.000000000                        0.000000000
#TNFRSF4                         0.000000000                        0.000000000
#TNFRSF14                        0.005738629                        0.003056466
#TNFRSF25                        0.000000000                        0.000000000
levels(cellchat@idents)
# [1] "SC_M1" "SC_M2" "SC_M3" "SC_M4" "SC_M5" "SC_M6" "Fb1"   "Fb2"   "Fb3"  
#[10] "Fb4"   "Fb5"   "Fb6"   "Fb7"   "Ch5"   "Ch2"   "Ch3"   "Ch1"   "Ch8"  
#[19] "Ch4"   "Ch6"   "Ch7" 
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
  setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/OA/test")
  
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

setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/OA/heatmap")
# [1] "ANGPTL"     "TGFb"       "MIF"        "CXCL"       "ACTIVIN"   
# [6] "FGF"        "CCL"        "MK"         "BMP"        "MSTN"      
#[11] "GDF"        "BMP10"      "VISFATIN"   "PDGF"       "CHEMERIN"  
#[16] "EGF"        "IGF"        "ANNEXIN"    "PTN"        "RANKL"     
#[21] "SEMA3"      "GRN"        "LIFR"       "PROS"       "IL6"       
#[26] "TWEAK"      "GALECTIN"   "CSF"        "GAS"        "TNF"       
#[31] "OSM"        "VEGF"       "COMPLEMENT" "IL4"        "LT"        
#[36] "IFN-II"     "IL1"        "CD40"       "PSAP"       "ANGPT"     
#[41] "IL2"        "BAFF"       "IL16"       "PERIOSTIN"  "KIT"       
#[46] "CD137"      "NT"         "NRG"        "TRAIL"      "LIGHT"     
#[51] "PARs"       "NGF"        "HGF"        "NPR2"       "AGT"       
#[56] "FASLG"      "CX3C"       "CD30"       "BTLA"       "EPO"       
#[61] "GDNF"       "GH"         "PRL"        "CSF3"       "FSH"  
pathways.show <- c("PRL")
par(mfrow=c(1,1))
pdf("TIL_PRL_heatmap.pdf", width = 10, height = 10)
netVisual_heatmap(cellchat, signaling = pathways.show, color.heatmap = "Reds")
dev.off()

levels(cellchat@idents)
# [1] "SC_M1" "SC_M2" "SC_M3" "SC_M4" "SC_M5" "SC_M6" "Fb1"   "Fb2"   "Fb3"  
#[10] "Fb4"   "Fb5"   "Fb6"   "Fb7"   "Ch5"   "Ch2"   "Ch3"   "Ch1"   "Ch8"  
#[19] "Ch4"   "Ch6"   "Ch7"

setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/OA/dot")
p = netVisual_bubble(cellchat, sources.use = c(14), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch5.pdf", p, width = 10, height = 30) 
# save as TIL/Mye_Lymph_bubble.pdf
p = netVisual_bubble(cellchat, sources.use = c(15), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch2.pdf", p, width = 10, height = 30) 
# save as TIL/Mye_Lymph_bubble.pdf
p = netVisual_bubble(cellchat, sources.use = c(16), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch3.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(17), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch1.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(18), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch8.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(19), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch4.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(20), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch6.pdf", p, width = 10, height = 30) 
p = netVisual_bubble(cellchat, sources.use = c(21), 
                     targets.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch7.pdf", p, width = 10, height = 30) 
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

setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/OA/singlerole")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", font.size = 8,width = 15, height = 25)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", font.size = 8,width = 15, height = 25)
pdf("TIL_SNA_SignalingPattern.pdf", width = 18, height = 26)
ht1 + ht2
dev.off()  

save(cellchat, file = "cellchat_OA_final.RData")
load("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/OA/singlerole/cellchat_OA_final.RData")
setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/OA/contribution")

pathways.show.all <- cellchat@netP$pathways
# check the order of cell identity to set suitable vertex.receiver
levels(cellchat@idents)

for (i in 1:length(pathways.show.all)) {

    # 计算并保存每个L-R对信号通路的贡献
    gg <- netAnalysis_contribution(cellchat, signaling = pathways.show.all[i])
    
    # 保存贡献图
    ggsave(filename = paste0(pathways.show.all[i], "_L-R_contribution.pdf"), 
           plot = gg, width = 12, height = 12, units = 'in', dpi = 300)
}


setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/OA/dot")


p = netVisual_bubble(cellchat, targets.use = c(14), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,15,16,17,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch5.1.pdf", p, width = 10, height = 30) 



p = netVisual_bubble(cellchat, targets.use = c(15), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,16,17,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch2.1.pdf", p, width = 10, height = 30) 


# save as TIL/Mye_Lymph_bubble.pdf
p = netVisual_bubble(cellchat, targets.use = c(16), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,17,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch3.1.pdf", p, width = 10, height = 30) 



p = netVisual_bubble(cellchat, targets.use = c(17), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,18,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch1.1.pdf", p, width = 10, height = 30) 


p = netVisual_bubble(cellchat, targets.use = c(18), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,19,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch8.1.pdf", p, width = 10, height = 30) 


p = netVisual_bubble(cellchat, targets.use = c(19), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,20,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch4.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(20), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,21), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch6.1.pdf", p, width = 10, height = 30) 

p = netVisual_bubble(cellchat, targets.use = c(21), 
                     sources.use = c(1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch7.1.pdf", p, width = 10, height = 30) 

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

#######细胞通讯
setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/OA_NEW")
load("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/OA/singlerole/cellchat_OA_final.RData")

levels(cellchat@idents)

p = netVisual_bubble(cellchat, sources.use = c(7,10), targets.use = c(17,15,16,19,14,20,21,18), signaling = c("MK"), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Fb1_4_OA.pdf", p, width = 10, height = 10)

p1 = netVisual_bubble(cellchat, sources.use = c(17,15,14), targets.use = c(7,8,9,10,11,12,13), signaling = c("VISFATIN"), remove.isolate = FALSE)
ggsave("Mye_Lymph_bubble_Ch1_2_5_OA.pdf", p1, width = 10, height = 5)


pathways.show <- c("ANGPTL","TGFb","MIF","CXCL","ACTIVIN","FGF","CCL","MK","BMP","MSTN","GDF","BMP10","VISFATIN","PDGF","CHEMERIN")
setwd("/mnt/8w/data7/yiyuan/Bone/merged_cellchat/OA_NEW")
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
ht1 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "outgoing", signaling = pathways.show, font.size = 8,width = 15, height = 10)
ht2 <- netAnalysis_signalingRole_heatmap(cellchat, pattern = "incoming", signaling = pathways.show, font.size = 8,width = 15, height = 10)
pdf("TIL_SNA_SignalingPattern_OA.pdf", width = 18, height = 10)
ht1 + ht2
dev.off()  
pdf("TIL_SNA_SignalingPattern_OA_outing.pdf", width = 18, height = 10)
ht1
dev.off()  
