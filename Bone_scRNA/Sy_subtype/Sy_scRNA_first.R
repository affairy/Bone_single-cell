##RAW读取
library(Seurat)
library(Matrix)
library(harmony)
library(ggplot2)
library(tidyverse)

setwd('/mnt/8w/data7/yiyuan/Bone/test1') #单细胞矩阵存放路径
dir = c('/mnt/8w/data2/jiangr/Bone/OAS1_1209/OAS1_1209/outs/raw',
        '/mnt/8w/data2/jiangr/Bone/OAS2_1209/OAS2_1209/outs/raw',
        '/mnt/8w/data2/jiangr/Bone/OAS3_0108/OAS3_0108/outs/raw',
        '/mnt/8w/data2/jiangr/Bone/HAS1_1130/HAS1_1130/outs/raw', 
        '/mnt/8w/data2/jiangr/Bone/HAS2_1226/HAS2_1226/outs/raw',
        '/mnt/8w/data2/jiangr/Bone/HAS3_0228/HAS3_0228/outs/raw')
names(dir) = c('OAS1_1209','OAS2_1209','OAS3_0108','HAS1_1130','HAS2_1226','HAS3_0228')
scRNAlist <- list()
#把每个样本数据创建一个seurat对象，并存放到列表scRNAlist里
for(i in 1:length(dir)){ # nolint
  counts <- Read10X(data.dir = dir[i],gene.column = 2)
  scRNAlist[[i]] <- CreateSeuratObject(counts,  min.cells=1, min.features = 150)
}
#使用merge将两个单细胞合并成一个seurat对象
scRNA <- merge(scRNAlist[[1]], y = scRNAlist[2:6])

exp=LayerData(scRNA, assay="RNA", layer='counts')

table(scRNA@meta.data$orig.ident)
# HAS1  HAS2  HAS3  OAS1  OAS2  OAS3 
#34146 32196 51630 92126 64953 31481 
sce=scRNA

######添加rawname
sample_id <- str_split(colnames(sce), "_", simplify = TRUE)[, 1:2]
sample_id <- apply(sample_id, 1, paste, collapse = "_")
sce@meta.data$rawname <- sample_id

sce@meta.data %>% head()
table(sce@meta.data$rawname)
#HAS1_1130 HAS2_1226 HAS3_0228 OAS1_1209 OAS2_1209 OAS3_0108 
#    34146     32196     51630     92126     64953     31481 
##添加sample
library(dplyr)
sce$sample<- recode(sce$rawname,
                        "OAS1_1209" = "OA1",
                        "OAS2_1209" = "OA2",
                        "OAS3_0108" = "OA3",
                        "HAS1_1130" = "HA1",
                        "HAS2_1226" = "HA2",
                        "HAS3_0228" = "HA3",
                        .default = NA_character_)
sce@meta.data %>% head()

table(sce@meta.data$sample)
#  HA1   HA2   HA3   OA1   OA2   OA3 
#34146 32196 51630 92126 64953 31481 

###创建cellid
library(stringr)
# 分割 rownames 并提取第3到第5列
barcodes_parts <- str_split(rownames(sce@meta.data), pattern = "_", simplify = TRUE)[, 3:5]
# 去除下划线并合并
sce@meta.data$barcodes <- apply(barcodes_parts, 1, function(x) paste(gsub("_", "", x), collapse = ""))
# 输出结果检查
head(sce@meta.data$barcodes)
#[1] "AGACGTTCAAACGCTAGTAACAAGTGG" "AGCTCCTTGAACGCTAGTAACAAGTGG"
#[3] "GCGAGTAACAACGCTAGTAACAAGTGG" "GCCTACGATAACGTCCAAAACAAGTGG"
#[5] "TACCTCTCCAACGTCCAAAACAAGTGG" "TGTGGACACAACGTCCAAAACAAGTGG"

sce@meta.data$cellid <- paste0(sce@meta.data$sample, "_", sce@meta.data$barcodes)
sce@meta.data %>% head()
#                                        orig.ident nCount_RNA nFeature_RNA
#OAS1_1209_AGACGTTCA_AACGCTAGT_AACAAGTGG       OAS1       2189          957
#OAS1_1209_AGCTCCTTG_AACGCTAGT_AACAAGTGG       OAS1        312          239
#OAS1_1209_GCGAGTAAC_AACGCTAGT_AACAAGTGG       OAS1        216          177
#OAS1_1209_GCCTACGAT_AACGTCCAA_AACAAGTGG       OAS1        820          481
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG       OAS1       7711         2569
#OAS1_1209_TGTGGACAC_AACGTCCAA_AACAAGTGG       OAS1        424          325
#                                                           barcodes
#OAS1_1209_AGACGTTCA_AACGCTAGT_AACAAGTGG AGACGTTCAAACGCTAGTAACAAGTGG
#OAS1_1209_AGCTCCTTG_AACGCTAGT_AACAAGTGG AGCTCCTTGAACGCTAGTAACAAGTGG
#OAS1_1209_GCGAGTAAC_AACGCTAGT_AACAAGTGG GCGAGTAACAACGCTAGTAACAAGTGG
#OAS1_1209_GCCTACGAT_AACGTCCAA_AACAAGTGG GCCTACGATAACGTCCAAAACAAGTGG
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG TACCTCTCCAACGTCCAAAACAAGTGG
#OAS1_1209_TGTGGACAC_AACGTCCAA_AACAAGTGG TGTGGACACAACGTCCAAAACAAGTGG
#                                                                 cellid
#OAS1_1209_AGACGTTCA_AACGCTAGT_AACAAGTGG OA1_AGACGTTCAAACGCTAGTAACAAGTGG
#OAS1_1209_AGCTCCTTG_AACGCTAGT_AACAAGTGG OA1_AGCTCCTTGAACGCTAGTAACAAGTGG
#OAS1_1209_GCGAGTAAC_AACGCTAGT_AACAAGTGG OA1_GCGAGTAACAACGCTAGTAACAAGTGG
#OAS1_1209_GCCTACGAT_AACGTCCAA_AACAAGTGG OA1_GCCTACGATAACGTCCAAAACAAGTGG
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG OA1_TACCTCTCCAACGTCCAAAACAAGTGG
#OAS1_1209_TGTGGACAC_AACGTCCAA_AACAAGTGG OA1_TGTGGACACAACGTCCAAAACAAGTGG
#                                          rawname sample
#OAS1_1209_AGACGTTCA_AACGCTAGT_AACAAGTGG OAS1_1209    OA1
#OAS1_1209_AGCTCCTTG_AACGCTAGT_AACAAGTGG OAS1_1209    OA1
#OAS1_1209_GCGAGTAAC_AACGCTAGT_AACAAGTGG OAS1_1209    OA1
#OAS1_1209_GCCTACGAT_AACGTCCAA_AACAAGTGG OAS1_1209    OA1
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG OAS1_1209    OA1
#OAS1_1209_TGTGGACAC_AACGTCCAA_AACAAGTGG OAS1_1209    OA1
#

main.anno <- data.table::fread("/mnt/8w/data7/yiyuan/Bone/anno/31_P22112601_B1_滑膜_cellinfo/Main/P22112601_B1/P22112601_B1_cellinfo.xls", data.table=F)
dim(main.anno)
#[1] 90608    20
colnames(main.anno)[1] <- "cellid"
# 检查修改后的列名
head(main.anno)
#                           cellid   rawname sample gname platform nFeature_RNA
#1 OA1_AACACACAGAACGCTAGTATGAGTCCT OAS1_1209    OA1    OA      sgr         2452
#2 OA1_AACACACAGAACGCTAGTGAGTTCGTC OAS1_1209    OA1    OA      sgr         2262
#3 OA1_AACACACAGAACGTCCAACTCTGCTTA OAS1_1209    OA1    OA      sgr         1067
#4 OA1_AACACACAGAAGGTGGTACTCTGCTTA OAS1_1209    OA1    OA      sgr         1314
#5 OA1_AACACACAGAAGGTGGTATCCAGACGA OAS1_1209    OA1    OA      sgr         1757
#6 OA1_AACACACAGAATCCGGTGATTAGCCTG OAS1_1209    OA1    OA      sgr         1182
#  nCount_RNA total_counts_mt percent_mt n_genes_by_counts total_counts
#1       9747             923   9.469581              2452         9747
#2       6700             443   6.611940              2262         6700
#3       2587             412  15.925783              1067         2587
#4       2855             190   6.654991              1314         2855
#5       4806             295   6.138160              1757         4806
#6       2952             184   6.233062              1182         2952
#  pct_counts_mt louvain raw_cluster sample_colors gname_colors
#1      9.469581       3           4       #0067AA      #0067AA
#2      6.611940       0           1       #0067AA      #0067AA
#3     15.925783       7           8       #0067AA      #0067AA
#4      6.654991       0           1       #0067AA      #0067AA
#5      6.138160       0           1       #0067AA      #0067AA
#6      6.233062       3           4       #0067AA      #0067AA
#  raw_cluster_colors     cluster cluster_standard cluster_colors
#1            #FF1F1D Fibroblasts      Fibroblasts        #FF7F00
#2            #0067AA Fibroblasts      Fibroblasts        #FF7F00
#3            #B6B800 Fibroblasts      Fibroblasts        #FF7F00
#4            #0067AA Fibroblasts      Fibroblasts        #FF7F00
#5            #0067AA Fibroblasts      Fibroblasts        #FF7F00
#6            #FF1F1D Fibroblasts      Fibroblasts        #FF7F00
subset_cols <- c("cellid", "rawname", "gname","sample_colors","cluster", "cluster_standard","cluster_colors")
# 提取指定列
main.anno_subset <- main.anno %>% select(all_of(subset_cols))
dim(main.anno_subset)
#[1] 90608     7
head(main.anno_subset)
#                           cellid   rawname gname sample_colors     cluster
#1 OA1_AACACACAGAACGCTAGTATGAGTCCT OAS1_1209    OA       #0067AA Fibroblasts
#2 OA1_AACACACAGAACGCTAGTGAGTTCGTC OAS1_1209    OA       #0067AA Fibroblasts
#3 OA1_AACACACAGAACGTCCAACTCTGCTTA OAS1_1209    OA       #0067AA Fibroblasts
#4 OA1_AACACACAGAAGGTGGTACTCTGCTTA OAS1_1209    OA       #0067AA Fibroblasts
#5 OA1_AACACACAGAAGGTGGTATCCAGACGA OAS1_1209    OA       #0067AA Fibroblasts
#6 OA1_AACACACAGAATCCGGTGATTAGCCTG OAS1_1209    OA       #0067AA Fibroblasts
#  cluster_standard cluster_colors
#1      Fibroblasts        #FF7F00
#2      Fibroblasts        #FF7F00
#3      Fibroblasts        #FF7F00
#4      Fibroblasts        #FF7F00
#5      Fibroblasts        #FF7F00
#6      Fibroblasts        #FF7F00

sce$celltype_cluster <- main.anno_subset[match(sce@meta.data$cellid, main.anno_subset$cellid), "cluster"]
sce$celltype_cluster <- ifelse(is.na(sce$celltype_cluster), "Unknow", sce$celltype_cluster)
table(sce@meta.data$celltype_cluster)
#            BCells                ECs        Fibroblasts          MastCells 
#              1098               6568              50346               3010 
#               MPs         MuralCells               pDCs ProliferatingCells 
#             20936               4839                 83                242 
#            TandNK             Unknow 
#              3486             215924 
#
sce$celltype_cluster_standard <- main.anno_subset[match(sce@meta.data$cellid, main.anno_subset$cellid), "cluster_standard"]
sce$celltype_cluster_standard <- ifelse(is.na(sce$celltype_cluster_standard), "Unknow", sce$celltype_cluster_standard)
table(sce@meta.data$celltype_cluster_standard)
#                     B cells            Endothelial cells 
#                        1098                         6568 
#                 Fibroblasts                   Mast cells 
#                       50346                         3010 
#      Mononuclear phagocytes                  Mural cells 
#                       20936                         4839 
#Plasmacytoid dendritic cells          Proliferating cells 
#                          83                          242 
#              T and NK cells                       Unknow 
#                        3486                       215924 

sce$cluster_colors <- main.anno_subset[match(sce@meta.data$cellid, main.anno_subset$cellid), "cluster_colors"]
sce$cluster_colors <- ifelse(is.na(sce$cluster_colors), "Unknow", sce$cluster_colors)
table(sce@meta.data$cluster_colors)
# #0067AA #00A23F #01C1CC #A763AC #B45B5D #B6B800 #FF1F1D #FF7F00 #FF8AB6  Unknow 
#   6568    4839      83    1098    3486   20936     242   50346    3010  215924 
sce$gname <- main.anno_subset[match(sce@meta.data$cellid, main.anno_subset$cellid), "gname"]
sce$gname <- ifelse(is.na(sce$cluster_colors), "Unknow", sce$gname)
table(sce@meta.data$gname)
#   HA    OA 
#43027 47581 
sce$rawname <- main.anno_subset[match(sce@meta.data$cellid, main.anno_subset$cellid), "rawname"]
sce$rawname <- ifelse(is.na(sce$cluster_colors), "Unknow", sce$rawname)
table(sce@meta.data$rawname)
#HAS1_1130 HAS2_1226 HAS3_0228 OAS1_1209 OAS2_1209 OAS3_0108 
#    15726     10809     16492     17326     12723     17532 

# 将 sce@meta.data 导出为 CSV 文件
write.csv(sce@meta.data, file = "sce_meta.data.csv", row.names = FALSE)
table(sce@meta.data$celltype_cluster_standard,sce@meta.data$sample)
#                              
#                                 HA1   HA2   HA3   OA1   OA2   OA3
#  B cells                        462    19   552    36    14    15
#  Endothelial cells              660  1894  2392   666   892    64
#  Fibroblasts                   7895  6258  4927 11214  8851 11201
#  Mast cells                    1471    50   627    81   187   594
#  Mononuclear phagocytes        3662  1023  5157  4173  1562  5359
#  Mural cells                    164  1467  1431   937   820    20
#  Plasmacytoid dendritic cells    14     0    64     1     0     4
#  Proliferating cells             25    18    32    66    92     9
#  T and NK cells                1373    80  1310   152   305   266
#  Unknow                       18420 21387 35138 74800 52230 13949

sce <- sce[, sce$celltype_cluster_standard %in% c("B cells","Endothelial cells","Fibroblasts",
                                                  "Mast cells","Mononuclear phagocytes","Mural cells",
                                                  "Plasmacytoid dendritic cells","Proliferating cells",
                                                  "T and NK cells" )]
table(sce@meta.data$celltype_cluster_standard,sce@meta.data$sample)
#                              
#                                 HA1   HA2   HA3   OA1   OA2   OA3
#  B cells                        462    19   552    36    14    15
#  Endothelial cells              660  1894  2392   666   892    64
#  Fibroblasts                   7895  6258  4927 11214  8851 11201
#  Mast cells                    1471    50   627    81   187   594
#  Mononuclear phagocytes        3662  1023  5157  4173  1562  5359
#  Mural cells                    164  1467  1431   937   820    20
#  Plasmacytoid dendritic cells    14     0    64     1     0     4
#  Proliferating cells             25    18    32    66    92     9
#  T and NK cells                1373    80  1310   152   305   266

saveRDS(sce,file="sce_first.rds")

sce <- PercentageFeatureSet(object = sce, pattern = "^MT-", col.name = "percent.mt")

write.csv(sce@meta.data, file = "sce_meta.data_9celltype.csv", row.names = FALSE)

sce <- subset(sce, subset = nFeature_RNA > 200 & nCount_RNA > 450 & percent.mt < 30)
sce
table(sce@meta.data$celltype_cluster_standard,sce@meta.data$sample)
#                              
#                                 HA1   HA2   HA3   OA1   OA2   OA3
#  B cells                        462    19   552    36    14    15
#  Endothelial cells              660  1894  2392   666   892    64
#  Fibroblasts                   7895  6258  4927 11214  8851 11201
#  Mast cells                    1471    50   627    81   187   594
#  Mononuclear phagocytes        3662  1023  5157  4173  1562  5359
#  Mural cells                    164  1467  1431   937   820    20
#  Plasmacytoid dendritic cells    14     0    64     1     0     4
#  Proliferating cells             25    18    32    66    92     9
#  T and NK cells                1373    80  1310   152   305   266
sce[["origin.counts"]] <- JoinLayers(sce[["RNA"]])
DefaultAssay(sce) <- "origin.counts"

sce <- NormalizeData(sce, normalization.method =  "LogNormalize",
                     scale.factor = 10000)
#GetAssay(sce,assay = "RNA")
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
  label.size = 3,  # 设置标签的字体大小
)
dev.off()

#### PCA拐点定量识别
pct <- sce1[["pca"]]@stdev / sum(sce1[["pca"]]@stdev) * 100 ; 
cumu <- cumsum(pct)
pc.use <- min(which(cumu > 90 & pct < 5)[1],sort(which((pct[1:length(pct) - 1] - pct[2:length(pct)]) > 0.1),decreasing = T)[1] + 1)

pc.use
# 16

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
sce1 <- FindClusters(sce1, algorithm = 1 , resolution = c(1))

sce1@meta.data %>% head()


write.csv(sce1@meta.data, file = "sce_meta.data_4celltype_1.csv", row.names = FALSE)

png(file = "umap_1_seurat_clusters.png", width = 3000, height = 2500, res = 300)
DimPlot(sce1,reduction = "umap",  group.by = "seurat_clusters", label = TRUE)
dev.off()

png(file = "umap_1_celltype_cluster_standard.png", width = 3000, height = 2500, res = 300)
DimPlot(sce1,reduction = "umap",  group.by = "celltype_cluster_standard", label = TRUE)
dev.off()

saveRDS(sce1,file="sce_first.2.rds")

###先挑成纤维
sce_subset <- sce1[, sce1$celltype_cluster_standard %in% c("Fibroblasts")]
sce_subset
head(sce_subset)
table(sce_subset@meta.data$celltype_cluster_standard,sce_subset@meta.data$sample)

####挑选巨噬细胞
sce_subset1 <- sce1[, sce1$celltype_cluster_standard %in% c("Mononuclear phagocytes")]
sce_subset1
head(sce_subset1)
table(sce_subset1@meta.data$celltype_cluster_standard,sce_subset1@meta.data$sample)

#######从一级中提取二级

#rm(list=ls())
#options(stringsAsFactors = F)

library(tidyr)
library(reshape2)
sce$sample = factor(sce$sample,levels = c('OA1','OA2','OA3','HA1','HA2','HA3'))
tb=table(sce$sample, sce$celltype_cluster_standard)
head(tb)
library (gplots)
setwd("/mnt/8w/data7/yiyuan/Bone/test1")
png(file = "balloonplot_1_celltype_cluster_standard.png", width = 4000, height = 2500, res = 300)
balloonplot(tb)
dev.off()


cell.Patient <- aggregate(sce$sample, list(sce$sample, sce$gname, sce$celltype_cluster_standard), length)
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


#library(RColorBrewer)
#brewer.pal.info
#qual_col_pals = brewer.pal.info[brewer.pal.info$category == 'qual',]
##处理后有73种差异还比较明显的颜色，基本够用
#col_vector = unlist(mapply(brewer.pal, qual_col_pals$maxcolors, rownames(qual_col_pals)))
#set.seed(seed = 1234)
#patient_col <- sample(col_vector, 50)
patient_col = c(
   'B cells' = '#ef070f',  # 浅粉色
   'Endothelial cells' = '#FAD02E',  # 浅黄色
   'Fibroblasts' = '#A8DADC',  # 浅蓝绿色
   'Mast cells' = '#aba3e7',  
   'Mononuclear phagocytes' = '#712820',  # 浅蓝色
   'Mural cells' = '#B7D7A8',  # 浅绿色
   'Plasmacytoid dendritic cells' = '#D4E157',  # 浅橄榄绿色
   'Proliferating cells' = '#F2826C',
   'T and NK cells' = '#B45B5D'
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