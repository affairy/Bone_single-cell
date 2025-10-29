
##RAW读取
library(Seurat)
library(Matrix)
library(harmony)
library(ggplot2)
library(tidyverse)

setwd('/mnt/8w/data7/yiyuan/Bone/test') #单细胞矩阵存放路径
dir = c('/mnt/8w/data2/jiangr/Bone/OA1_1125/OA1_1125/outs/raw',
        '/mnt/8w/data2/jiangr/Bone/OA2_1127/OA2_1127/outs/raw',
        '/mnt/8w/data2/jiangr/Bone/OA3_1130/OA3_1130/outs/raw',
        '/mnt/8w/data2/jiangr/Bone/HAC1_1130/HAC1_1130/outs/raw', 
        '/mnt/8w/data2/jiangr/Bone/HAC2_1226/HAC2_1226/outs/raw',
        '/mnt/8w/data2/jiangr/Bone/HAC3_0228/HAC3_0228/outs/raw',
        '/mnt/8w/data7/yiyuan/Bone/HAC_0922/raw',
        '/mnt/8w/data7/yiyuan/Bone/HAC_0928/raw')
names(dir) = c('OA1_1125','OA2_1127','OA3_1130','HAC1_1130','HAC2_1226','HAC3_0228','HAC_0922','HAC_0928')
scRNAlist <- list()
#把每个样本数据创建一个seurat对象，并存放到列表scRNAlist里
for(i in 1:length(dir)){ # nolint
  counts <- Read10X(data.dir = dir[i],gene.column = 2)
  scRNAlist[[i]] <- CreateSeuratObject(counts,  min.cells=1, min.features = 150)
}
#使用merge将两个单细胞合并成一个seurat对象
scRNA <- merge(scRNAlist[[1]], y = scRNAlist[2:8])
exp=LayerData(scRNA, assay="RNA", layer='counts')
table(scRNA@meta.data$orig.ident)
#  HAC  HAC1  HAC2  HAC3   OA1   OA2   OA3 
#25947 10723  9107 36941 12390 11152 15886 
sce=scRNA

######添加raw_name
sample_id <- str_split(colnames(sce), "_", simplify = TRUE)[, 1:2]
sample_id <- apply(sample_id, 1, paste, collapse = "_")
sce@meta.data$raw_name <- sample_id

sce@meta.data %>% head()
table(sce@meta.data$raw_name)
#HAC_0922  HAC_0928 HAC1_1130 HAC2_1226 HAC3_0228  OA1_1125  OA2_1127  OA3_1130 
#    9476     16471     10723      9107     36941     12390     11152     15886 
##添加sample
library(dplyr)
sce$sample<- recode(sce$raw_name,
                        "OA1_1125" = "OA1",
                        "OA2_1127" = "OA2",
                        "OA3_1130" = "OA3",
                        "HAC1_1130" = "HA1",
                        "HAC2_1226" = "HA2",
                        "HAC3_0228" = "HA3",
                        "HAC_0922" = "HA4",
                        "HAC_0928" = "HA5",
                        .default = NA_character_)
sce@meta.data %>% head()

table(sce@meta.data$sample)
#  HA1   HA2   HA3   HA4   HA5   OA1   OA2   OA3 
#10723  9107 36941  9476 16471 12390 11152 15886

###创建cellid
library(stringr)
# 分割 rownames 并提取第3到第5列
barcodes_parts <- str_split(rownames(sce@meta.data), pattern = "_", simplify = TRUE)[, 3:5]
# 去除下划线并合并
sce@meta.data$barcodes <- apply(barcodes_parts, 1, function(x) paste(gsub("_", "", x), collapse = ""))
# 输出结果检查
head(sce@meta.data$barcodes)
sce@meta.data$cellid <- paste0(sce@meta.data$sample, "_", sce@meta.data$barcodes)
sce@meta.data %>% head()
#                                       orig.ident nCount_RNA nFeature_RNA
#OA1_1125_GATCCATGC_AACGCTAGT_AACAAGTGG        OA1       7080         2025
#OA1_1125_TACCGTCTG_AACGTCCAA_AACAAGTGG        OA1       5076         1440
#OA1_1125_TCACTGGAA_AACGTCCAA_AACAAGTGG        OA1       6050         1838
#OA1_1125_GAGCAGCTT_AAGGTGGTA_AACAAGTGG        OA1        355          215
#OA1_1125_CGATCGGTA_AATCCGGTG_AACAAGTGG        OA1       8491         2009
#OA1_1125_AGACGTTCA_ACCAGGTCA_AACAAGTGG        OA1       6404         2016
#                                       raw_name sample
#OA1_1125_GATCCATGC_AACGCTAGT_AACAAGTGG OA1_1125    OA1
#OA1_1125_TACCGTCTG_AACGTCCAA_AACAAGTGG OA1_1125    OA1
#OA1_1125_TCACTGGAA_AACGTCCAA_AACAAGTGG OA1_1125    OA1
#OA1_1125_GAGCAGCTT_AAGGTGGTA_AACAAGTGG OA1_1125    OA1
#OA1_1125_CGATCGGTA_AATCCGGTG_AACAAGTGG OA1_1125    OA1
#OA1_1125_AGACGTTCA_ACCAGGTCA_AACAAGTGG OA1_1125    OA1
#                                                          barcodes
#OA1_1125_GATCCATGC_AACGCTAGT_AACAAGTGG GATCCATGCAACGCTAGTAACAAGTGG
#OA1_1125_TACCGTCTG_AACGTCCAA_AACAAGTGG TACCGTCTGAACGTCCAAAACAAGTGG
#OA1_1125_TCACTGGAA_AACGTCCAA_AACAAGTGG TCACTGGAAAACGTCCAAAACAAGTGG
#OA1_1125_GAGCAGCTT_AAGGTGGTA_AACAAGTGG GAGCAGCTTAAGGTGGTAAACAAGTGG
#OA1_1125_CGATCGGTA_AATCCGGTG_AACAAGTGG CGATCGGTAAATCCGGTGAACAAGTGG
#OA1_1125_AGACGTTCA_ACCAGGTCA_AACAAGTGG AGACGTTCAACCAGGTCAAACAAGTGG
#                                                                cellid
#OA1_1125_GATCCATGC_AACGCTAGT_AACAAGTGG OA1_GATCCATGCAACGCTAGTAACAAGTGG
#OA1_1125_TACCGTCTG_AACGTCCAA_AACAAGTGG OA1_TACCGTCTGAACGTCCAAAACAAGTGG
#OA1_1125_TCACTGGAA_AACGTCCAA_AACAAGTGG OA1_TCACTGGAAAACGTCCAAAACAAGTGG
#OA1_1125_GAGCAGCTT_AAGGTGGTA_AACAAGTGG OA1_GAGCAGCTTAAGGTGGTAAACAAGTGG
#OA1_1125_CGATCGGTA_AATCCGGTG_AACAAGTGG OA1_CGATCGGTAAATCCGGTGAACAAGTGG
#OA1_1125_AGACGTTCA_ACCAGGTCA_AACAAGTGG OA1_AGACGTTCAACCAGGTCAAACAAGTGG

main.anno <- data.table::fread("/mnt/8w/data7/yiyuan/Bone/anno/25_P23092210_B1_血友病软骨_cellinfo/Main/P23092210_B1/P23092210_B1_cellinfo.xls", data.table=F)
dim(main.anno)
#[1] 53287    27

colnames(main.anno)[1] <- "cellid"
# 检查修改后的列名
head(main.anno)
subset_cols <- c("cellid", "rawname", "gname","sample_colors","cluster", "cluster_standard","cluster_colors")
# 提取指定列
main.anno_subset <- main.anno %>% select(all_of(subset_cols))
dim(main.anno_subset)
#[1] 53287     7
head(main.anno_subset)
#                           cellid  rawname gname sample_colors      cluster
#1 OA1_AACACACAGAACGCTAGTACGACGTAT OA1_1125    OA       #0067AA Chondrocytes
#2 OA1_AACACACAGAACGCTAGTACGCCACTT OA1_1125    OA       #0067AA Chondrocytes
#3 OA1_AACACACAGAAGGTGGTAAGAATGACC OA1_1125    OA       #0067AA Chondrocytes
#4 OA1_AACACACAGAATCCGGTGGAACACAGA OA1_1125    OA       #0067AA Chondrocytes
#5 OA1_AACACACAGACAACAGGTTGAAGCTGA OA1_1125    OA       #0067AA Chondrocytes
#6 OA1_AACACACAGACACCAACGCACACCTCA OA1_1125    OA       #0067AA  Fibroblasts
#  cluster_standard cluster_colors
#1     Chondrocytes        #0067AA
#2     Chondrocytes        #0067AA
#3     Chondrocytes        #0067AA
#4     Chondrocytes        #0067AA
#5     Chondrocytes        #0067AA
#6      Fibroblasts        #00A23F

sce$celltype_cluster <- main.anno_subset[match(sce@meta.data$cellid, main.anno_subset$cellid), "cluster"]
sce$celltype_cluster <- ifelse(is.na(sce$celltype_cluster), "Unknow", sce$celltype_cluster)
table(sce@meta.data$celltype_cluster)
#Chondrocytes          ECs  Fibroblasts  ImmuneCells       Unknow 
#       37791          351        13097         2048        68859 

sce$celltype_cluster_standard <- main.anno_subset[match(sce@meta.data$cellid, main.anno_subset$cellid), "cluster_standard"]
sce$celltype_cluster_standard <- ifelse(is.na(sce$celltype_cluster_standard), "Unknow", sce$celltype_cluster_standard)
table(sce@meta.data$celltype_cluster_standard)
#     Chondrocytes Endothelial cells       Fibroblasts      Immune cells      Unknow 
#            37791               351             13097              2048       68859 

sce$cluster_colors <- main.anno_subset[match(sce@meta.data$cellid, main.anno_subset$cellid), "cluster_colors"]
sce$cluster_colors <- ifelse(is.na(sce$cluster_colors), "Unknow", sce$cluster_colors)
table(sce@meta.data$cluster_colors)
# #0067AA #00A23F #FF1F1D #FF7F00  Unknow 
#  37791   13097    2048     351   68859 
sce$gname <- main.anno_subset[match(sce@meta.data$cellid, main.anno_subset$cellid), "gname"]
sce$gname <- ifelse(is.na(sce$cluster_colors), "Unknow", sce$gname)
table(sce@meta.data$gname)
#   HA    OA 
#29461 23826 
sce$rawname <- main.anno_subset[match(sce@meta.data$cellid, main.anno_subset$cellid), "rawname"]
sce$rawname <- ifelse(is.na(sce$cluster_colors), "Unknow", sce$rawname)
table(sce@meta.data$rawname)

# HAC_0922  HAC_0928 HAC1_1130 HAC2_1226 HAC3_0228  OA1_1125  OA2_1127  OA3_1130 
#     3074      3733      6488      6057     10109      8797      6407      8622 
# 将 sce@meta.data 导出为 CSV 文件
write.csv(sce@meta.data, file = "sce_meta.data.csv", row.names = FALSE)
table(sce@meta.data$celltype_cluster_standard,sce@meta.data$sample)
                   
#                      HA1   HA2   HA3   HA4   HA5   OA1   OA2   OA3
#  Chondrocytes       4050  5153  4203  2629  2807  8275  5244  5430
#  Endothelial cells   111     4   161     8    26     0     0    41
#  Fibroblasts        2020   845  4144   405   891   522  1158  3112
#  Immune cells        307    55  1601    32     9     0     5    39
#  Unknow             4235  3050 26832  6402 12738  3593  4745  7264

sce <- sce[, sce$celltype_cluster_standard %in% c("Chondrocytes","Endothelial cells","Fibroblasts","Immune cells")]
table(sce@meta.data$celltype_cluster_standard,sce@meta.data$sample)
#                     HA1  HA2  HA3  HA4  HA5  OA1  OA2  OA3
#  Chondrocytes      4050 5153 4203 2629 2807 8275 5244 5430
#  Endothelial cells  111    4  161    8   26    0    0   41
#  Fibroblasts       2020  845 4144  405  891  522 1158 3112
#  Immune cells       307   55 1601   32    9    0    5   39
saveRDS(sce,file="sce_first.rds")

sce <- PercentageFeatureSet(object = sce, pattern = "^MT-", col.name = "percent.mt")

write.csv(sce@meta.data, file = "sce_meta.data_4celltype.csv", row.names = FALSE)

sce <- subset(sce, subset = nFeature_RNA > 100 & nCount_RNA > 350 & percent.mt < 60)
sce
table(sce@meta.data$celltype_cluster_standard,sce@meta.data$sample)
#                     HA1  HA2  HA3  HA4  HA5  OA1  OA2  OA3
#  Chondrocytes      4050 5153 4203 2629 2807 8275 5244 5430
#  Endothelial cells  111    4  161    8   26    0    0   41
#  Fibroblasts       2020  845 4144  405  891  522 1158 3112
#  Immune cells       307   55 1601   32    9    0    5   39


sce <- NormalizeData(sce, normalization.method =  "LogNormalize",
                     scale.factor = 10000)
GetAssay(sce,assay = "RNA")
sce@meta.data %>% head()
#前2000个高变feature RNA
sce <- FindVariableFeatures(sce,
                            selection.method = "vst", nfeatures = 2000)
# 中心化，为下一步PCA做准备
sce <- ScaleData(sce, vars.to.regress = c("nCount_RNA", "nFeature_RNA"))
sce <- RunPCA(object = sce, pc.genes = VariableFeatures(sce))

sce1=sce
#sce1 = readRDS('/mnt/8w/data7/yiyuan/Bone/test/sce_first.1.rds')
library(harmony)
sce1<- sce1%>% harmony::RunHarmony("rawname")

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

png(file = "rawname_harmony.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce1, 
  group.by = "rawname", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
)
dev.off()

saveRDS(sce1,file="sce_first.1.rds")

####筛选

png(file = "celltype_cluster_standard.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce1, 
  group.by = "celltype_cluster_standard", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
)
dev.off()

#### PCA拐点定量识别
pct <- sce1[["pca"]]@stdev / sum(sce1[["pca"]]@stdev) * 100 ; 
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)

pc.use
# 13 

png(file = "ElbowPlot.png", width = 3000, height = 2500, res = 300)
ElbowPlot(sce, ndims = 40, reduction = "pca")
dev.off()

#看看前20个主成分
pdf(file="PCA.pdf",width=10,height=8)
DimHeatmap(sce, dims = 1:20, cells = 100, balanced = TRUE,nfeatures = 10)
dev.off()

pcs=pc.use

set.seed(123)
sce1 <- RunUMAP(sce1,reduction = "harmony", dims = 1:pcs)
sce1 <- RunTSNE(sce1, dims.use = 1:pcs, reduction = "harmony",do.fast = TRUE)
sce1 <- FindNeighbors(sce1, reduction = "harmony", dims = 1:pcs)
sce1 <- FindClusters(sce1, algorithm = 1 , resolution = c(0.4))

write.csv(sce1@meta.data, file = "sce_meta.data_4celltype_0.4.csv", row.names = FALSE)

png(file = "umap_0.4_seurat_clusters.png", width = 3000, height = 2500, res = 300)
DimPlot(sce1,reduction = "umap",  group.by = "seurat_clusters", label = TRUE)
dev.off()

png(file = "umap_0.4_celltype_cluster_standard.png", width = 3000, height = 2500, res = 300)
DimPlot(sce1,reduction = "umap",  group.by = "celltype_cluster_standard", label = TRUE)
dev.off()

png(file = "umap_0.4_FeaturePlot.png", width = 3500, height = 1500, res = 300)
FeaturePlot(sce1,features = c("PCNA","MKI67"),reduction = "umap")
dev.off()

saveRDS(sce1,file="sce_first.2.rds")
#sce2 = readRDS('/mnt/8w/data7/yiyuan/Bone/test/sce_first.2.rds')

sce1[["origin.counts"]] <- JoinLayers(sce1[["RNA"]])
DefaultAssay(sce1) <- "origin.counts"

s.genes=Seurat::cc.genes.updated.2019$s.genes
g2m.genes=Seurat::cc.genes.updated.2019$g2m.genes
sce1 <- CellCycleScoring(sce1, s.features = s.genes, g2m.features = g2m.genes, set.ident = TRUE)

write.csv(sce1@meta.data, file = "sce_meta.data_4celltype_cellcycle.csv", row.names = FALSE)

head(sce1@meta.data,2)

png(file = "cellcyclescore.png", width = 2000, height = 1500, res = 300)
sce1@meta.data  %>% ggplot(aes(S.Score,G2M.Score))+geom_point(aes(color=Phase))+
    theme_minimal()
dev.off()

png(file = "umap_cellcyclescore.png", width = 2000, height = 1500, res = 300)
DimPlot(sce1,reduction = "umap")
dev.off()

#png(file = "RidgePlot.png", width = 2000, height = 1500, res = 300)
#RidgePlot(sce1, 
#          features = c("PCNA", "TOP2A", "MCM6", "MKI67"), 
#          ncol = 2)
#dev.off()
sce1$CC.Difference <- sce1$S.Score - sce1$G2M.Score
sce1 <- ScaleData(sce1, 
                    vars.to.regress = "CC.Difference", 
                    features = rownames(sce1))
sce1 <- RunPCA(object = sce1, pc.genes = VariableFeatures(sce))

pct <- sce1[["pca"]]@stdev / sum(sce1[["pca"]]@stdev) * 100 ; 
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)

pc.use

set.seed(123)
sce1 <- RunUMAP(sce1,reduction = "harmony", dims = 1:pc.use)
sce1 <- RunTSNE(sce1, dims.use = 1:pc.use, reduction = "harmony",do.fast = TRUE)
sce1 <- FindNeighbors(sce1, reduction = "harmony", dims = 1:pc.use)
sce1 <- FindClusters(sce1, algorithm = 1, resolution = 0.4, graph.name = "RNA_snn")


png(file = "umap_Phase.png", width = 2000, height = 1500, res = 300)
DimPlot(sce1,reduction = "umap",group.by = "Phase")
dev.off()

write.csv(sce1@meta.data, file = "sce_meta.data_4celltype_new0.4.csv", row.names = FALSE)

png(file = "umap_0.4_seurat_clusters_new.png", width = 3000, height = 2500, res = 300)
DimPlot(sce1,reduction = "umap",  group.by = "seurat_clusters", label = TRUE)
dev.off()

png(file = "umap_0.4_celltype_cluster_standard_new.png", width = 3000, height = 2500, res = 300)
DimPlot(sce1,reduction = "umap",  group.by = "celltype_cluster_standard", label = TRUE)
dev.off()

saveRDS(sce1,file="sce_first.3.rds")

png(file = "umap_0.4_celltype_cluster_standard_new.png", width = 3000, height = 2500, res = 300)
DimPlot(sce1, 
        reduction = "umap", 
        group.by = "celltype_cluster_standard", 
        label = TRUE, 
        cols = sce1@meta.data$cluster_colors)   # 使用 cluster_colors 作为颜色
dev.off()

###先挑选软骨
sce_subset <- sce1[, sce1$celltype_cluster_standard %in% c("Chondrocytes")]
sce_subset
head(sce_subset)
table(sce_subset@meta.data$celltype_cluster_standard,sce_subset@mreta.data$sample)
#######从一级中提取二级

sce.all=readRDS("/mnt/8w/data7/yiyuan/Bone/test/sce_first.3.rds")

library(tidyr)
library(reshape2)
sce.all$sample = factor(sce.all$sample,levels = c('OA1','OA2','OA3','HA1','HA2','HA3','HA4','HA5'))
tb=table(sce.all$sample, sce.all$celltype_cluster_standard)
head(tb)
library (gplots)
png(file = "balloonplot_1_celltype_cluster_standard.png", width = 3000, height = 2500, res = 300)
balloonplot(tb)
dev.off()


cell.Patient <- aggregate(sce.all$sample, list(sce.all$sample, sce.all$gname, sce.all$celltype_cluster_standard), length)
colnames(cell.Patient) <- c("sample","gname","cellType","number")
cell.Patient$Patient.type <- paste(cell.Patient$sample, cell.Patient$gname, sep = "-")
sce.all@meta.data$Patient.type <- paste(sce.all@meta.data$sample, sce.all@meta.data$gname, sep = "-")
patient <- unique(sce.all$Patient.type)
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


#library(RColorBrewer)
#brewer.pal.info
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
##处理后有73种差异还比较明显的颜色，基本够用
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#set.seed(seed = 1234)
#patient_col <- sample(col_vector, 50)
patient_col = c(
   'Chondrocytes' = '#ef070f',  # 浅粉色
   'Endothelial cells' = '#A8DADC',  # 浅蓝绿色
   'Fibroblasts' = '#D4E157',  # 浅橄榄绿色
   'Immune cells' = '#F2826C' 
)

# 设置图形的高度和宽度
pdf(file = "proportion_first.pdf", width = 10, height = 12)

# 绘制堆积柱状图
ggplot(data = cell.pro.Patient, mapping = aes(x = Patient.type, fill = cellType, y = proportion)) +
  geom_col(position = "stack") +
  scale_fill_manual(values = patient_col) +
  labs(x = "Samples", y = "Ratio of Cell Types (%)") +
  theme_classic() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
  geom_segment(mapping = aes(x = x_start, y = 1, xend = x_end, yend = 1, color = gname), size = 4)
dev.off()