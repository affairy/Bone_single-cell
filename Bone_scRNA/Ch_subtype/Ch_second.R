#################
setwd("/mnt/8w/data7/yiyuan/Bone/test/Ch")
library(Seurat)
library(Matrix)
library(harmony)
library(ggplot2)
library(tidyverse)
###先挑选软骨
sce_subset <- sce1[, sce1$celltype_cluster_standard %in% c("Chondrocytes")]
sce_subset
head(sce_subset@meta.data)
# 删除不需要的列
sce_subset@meta.data <- sce_subset@meta.data[, !(colnames(sce_subset@meta.data) %in% c("RNA_snn_res.0.5", "RNA_snn_res.0.1", "RNA_snn_res.0.2", "RNA_snn_res.0.3","raw_name"))]
# 查看删除后的元数据
head(sce_subset@meta.data)
#                                     orig.ident nCount_RNA nFeature_RNA sample
#HAC_0928_AACACACAGAATCCGGTGCCTCTTAAG        HAC       3338          907    HA5
#HAC_0928_AACACACAGACAACAGGTCTAGGTTGC        HAC       1215          526    HA5
#HAC_0928_AACACACAGACACCAACGAGCGTGTAT        HAC       6470         1677    HA5
#HAC_0928_AACACACAGACCTGCTACACCTGAAGG        HAC       7963         1782    HA5
#HAC_0928_AACACACAGACTCATCGTTGCCTTGGA        HAC       1925          617    HA5
#HAC_0928_AACACACAGAGACGAAGTGTCTCGAAG        HAC       3263         1033    HA5
#                                                        barcodes
#HAC_0928_AACACACAGAATCCGGTGCCTCTTAAG AACACACAGAATCCGGTGCCTCTTAAG
#HAC_0928_AACACACAGACAACAGGTCTAGGTTGC AACACACAGACAACAGGTCTAGGTTGC
#HAC_0928_AACACACAGACACCAACGAGCGTGTAT AACACACAGACACCAACGAGCGTGTAT
#HAC_0928_AACACACAGACCTGCTACACCTGAAGG AACACACAGACCTGCTACACCTGAAGG
#HAC_0928_AACACACAGACTCATCGTTGCCTTGGA AACACACAGACTCATCGTTGCCTTGGA
#HAC_0928_AACACACAGAGACGAAGTGTCTCGAAG AACACACAGAGACGAAGTGTCTCGAAG
#                                                              cellid
#HAC_0928_AACACACAGAATCCGGTGCCTCTTAAG HA5_AACACACAGAATCCGGTGCCTCTTAAG
#HAC_0928_AACACACAGACAACAGGTCTAGGTTGC HA5_AACACACAGACAACAGGTCTAGGTTGC
#HAC_0928_AACACACAGACACCAACGAGCGTGTAT HA5_AACACACAGACACCAACGAGCGTGTAT
#HAC_0928_AACACACAGACCTGCTACACCTGAAGG HA5_AACACACAGACCTGCTACACCTGAAGG
#HAC_0928_AACACACAGACTCATCGTTGCCTTGGA HA5_AACACACAGACTCATCGTTGCCTTGGA
#HAC_0928_AACACACAGAGACGAAGTGTCTCGAAG HA5_AACACACAGAGACGAAGTGTCTCGAAG
#                                     celltype_cluster celltype_cluster_standard
#HAC_0928_AACACACAGAATCCGGTGCCTCTTAAG     Chondrocytes              Chondrocytes
#HAC_0928_AACACACAGACAACAGGTCTAGGTTGC     Chondrocytes              Chondrocytes
#HAC_0928_AACACACAGACACCAACGAGCGTGTAT     Chondrocytes              Chondrocytes
#HAC_0928_AACACACAGACCTGCTACACCTGAAGG     Chondrocytes              Chondrocytes
#HAC_0928_AACACACAGACTCATCGTTGCCTTGGA     Chondrocytes              Chondrocytes
#HAC_0928_AACACACAGAGACGAAGTGTCTCGAAG     Chondrocytes              Chondrocytes
#                                     cluster_colors gname  rawname percent.mt
#HAC_0928_AACACACAGAATCCGGTGCCTCTTAAG        #0067AA    HA HAC_0928   2.007190
#HAC_0928_AACACACAGACAACAGGTCTAGGTTGC        #0067AA    HA HAC_0928   4.362140
#HAC_0928_AACACACAGACACCAACGAGCGTGTAT        #0067AA    HA HAC_0928   3.972179
#HAC_0928_AACACACAGACCTGCTACACCTGAAGG        #0067AA    HA HAC_0928   2.725104
#HAC_0928_AACACACAGACTCATCGTTGCCTTGGA        #0067AA    HA HAC_0928   5.454545
#HAC_0928_AACACACAGAGACGAAGTGTCTCGAAG        #0067AA    HA HAC_0928   5.516396
#                                     seurat_clusters RNA_snn_res.0.4
#HAC_0928_AACACACAGAATCCGGTGCCTCTTAAG               0               0
#HAC_0928_AACACACAGACAACAGGTCTAGGTTGC               2               2
#HAC_0928_AACACACAGACACCAACGAGCGTGTAT               1               1
#HAC_0928_AACACACAGACCTGCTACACCTGAAGG               2               2
#HAC_0928_AACACACAGACTCATCGTTGCCTTGGA               4               4
#HAC_0928_AACACACAGAGACGAAGTGTCTCGAAG               4               4
#                                         S.Score    G2M.Score Phase old.ident
#HAC_0928_AACACACAGAATCCGGTGCCTCTTAAG -0.01798674 -0.014401187    G1         0
#HAC_0928_AACACACAGACAACAGGTCTAGGTTGC  0.01891331 -0.009213498     S         1
#HAC_0928_AACACACAGACACCAACGAGCGTGTAT -0.03327628  0.042704271   G2M         2
#HAC_0928_AACACACAGACCTGCTACACCTGAAGG -0.01424363  0.046587328   G2M         1
#HAC_0928_AACACACAGACTCATCGTTGCCTTGGA -0.03524419 -0.024230199    G1         5
#HAC_0928_AACACACAGAGACGAAGTGTCTCGAAG -0.02374740  0.026601844   G2M         0
#                                     CC.Difference
#HAC_0928_AACACACAGAATCCGGTGCCTCTTAAG  -0.003585548
#HAC_0928_AACACACAGACAACAGGTCTAGGTTGC   0.028126810
#HAC_0928_AACACACAGACACCAACGAGCGTGTAT  -0.075980553
#HAC_0928_AACACACAGACCTGCTACACCTGAAGG  -0.060830958
#HAC_0928_AACACACAGACTCATCGTTGCCTTGGA  -0.011013994
#HAC_0928_AACACACAGAGACGAAGTGTCTCGAAG  -0.050349244

table(sce_subset@meta.data$celltype_cluster_standard,sce_subset@meta.data$sample)
#              
#                HA1  HA2  HA3  HA4  HA5  OA1  OA2  OA3
#  Chondrocytes 4050 5153 4203 2629 2807 8275 5244 5430   37791
####和公司注释多了些细胞 OA1,OA3,HA1,HA2,HA3
#	             OA1	 OA2	 OA3	 HA1	 HA2	 HA3	 HA4	 HA5
#Chondrocytes	8272	5244	5423	4046	5150	4016	2629	2807   37587
write.csv(sce_subset@meta.data, file = "sce_subset_meta.data_Chondrocytes.csv", row.names = FALSE)

saveRDS(sce_subset,file="/mnt/8w/data7/yiyuan/Bone/test/Ch/sce_first.Ch.rds")

sce_subset = readRDS('/mnt/8w/data7/yiyuan/Bone/test/Ch/sce_first.Ch.rds')

sce=sce_subset

Ch.anno <- data.table::fread("/mnt/8w/data7/yiyuan/Bone/anno/25_P23092210_B1_血友病软骨_cellinfo/Chondrocytes_NewColor/P23092210_B1/P23092210_B1_cellinfo.xls", data.table=F)
dim(Ch.anno)
#[1] 37587    28
colnames(Ch.anno)[1] <- "cellid"
# 检查修改后的列名
head(Ch.anno)
#                           cellid  rawname sample gname platform nFeature_RNA
#1 OA1_AACACACAGAACGCTAGTACGACGTAT OA1_1125    OA1    OA      10X         2598
#2 OA1_AACACACAGAACGCTAGTACGCCACTT OA1_1125    OA1    OA      10X          665
#3 OA1_AACACACAGAAGGTGGTAAGAATGACC OA1_1125    OA1    OA      10X         3430
#4 OA1_AACACACAGAATCCGGTGGAACACAGA OA1_1125    OA1    OA      10X         3085
#5 OA1_AACACACAGACAACAGGTTGAAGCTGA OA1_1125    OA1    OA      10X         2583
#6 OA1_AACACACAGACACCAACGCTTAGCGAC OA1_1125    OA1    OA      10X         3076
#  nCount_RNA total_counts_mt percent_mt total_counts_plat pct_counts_plat
#1      13681             238   1.739639                 0               0
#2       1725              47   2.724638                 0               0
#3      17969             365   2.031276                 0               0
#4      17737             300   1.691380                 0               0
#5      11008             250   2.271075                 0               0
#6      19121             460   2.405732                 0               0
#  total_counts_ribo pct_counts_ribo n_genes_by_counts total_counts
#1               878        6.417660              2598        13681
#2               288       16.695652               665         1725
#3              1351        7.518504              3430        17969
#4              1111        6.263742              3085        17737
#5               720        6.540698              2583        11008
#6              1288        6.736049              3076        19121
#  pct_counts_mt     S_score   G2M_score phase louvain raw_cluster sample_colors
#1      1.739639 -0.10156769 -0.50612745    G1       4           5       #0067AA
#2      2.724638 -0.04404145  0.03553922   G2M       7           8       #0067AA
#3      2.031276 -0.00006640 -0.34068627    G1       4           5       #0067AA
#4      1.691380 -0.07672379 -0.48529412    G1       0           1       #0067AA
#5      2.271075 -0.08369868 -0.16911765    G1       1           2       #0067AA
#6      2.405732 -0.14042779 -0.68382353    G1       4           5       #0067AA
#  gname_colors raw_cluster_colors cluster_standard cluster cluster_colors
#1      #0067AA            #A763AC              Ch5     Ch5        #54181A
#2      #0067AA            #B6B800              Ch8     Ch8        #F2826C
#3      #0067AA            #A763AC              Ch5     Ch5        #54181A
#4      #0067AA            #0067AA              Ch1     Ch1        #A25B1D
#5      #0067AA            #FF7F00              Ch2     Ch2        #9D1F20
#6      #0067AA            #A763AC              Ch5     Ch5        #54181A
#  annot_modified
#1            Ch5
#2            Ch8
#3            Ch5
#4            Ch1
#5            Ch2
#6            Ch5
subset_cols <- c("cellid", "rawname", "gname","sample_colors","cluster","cluster_colors")
# 提取指定列
Ch.anno_subset <- Ch.anno %>% select(all_of(subset_cols))
dim(Ch.anno_subset)
#[1] 37587     6
head(Ch.anno_subset)
#                           cellid  rawname gname sample_colors cluster
#1 OA1_AACACACAGAACGCTAGTACGACGTAT OA1_1125    OA       #0067AA     Ch5
#2 OA1_AACACACAGAACGCTAGTACGCCACTT OA1_1125    OA       #0067AA     Ch8
#3 OA1_AACACACAGAAGGTGGTAAGAATGACC OA1_1125    OA       #0067AA     Ch5
#4 OA1_AACACACAGAATCCGGTGGAACACAGA OA1_1125    OA       #0067AA     Ch1
#5 OA1_AACACACAGACAACAGGTTGAAGCTGA OA1_1125    OA       #0067AA     Ch2
#6 OA1_AACACACAGACACCAACGCTTAGCGAC OA1_1125    OA       #0067AA     Ch5
#  cluster_colors
#1        #54181A
#2        #F2826C
#3        #54181A
#4        #A25B1D
#5        #9D1F20
#6        #54181A
sce$cluster <- Ch.anno_subset[match(sce@meta.data$cellid, Ch.anno_subset$cellid), "cluster"]
sce$cluster <- ifelse(is.na(sce$cluster), "Unknow", sce$cluster)
table(sce@meta.data$cluster)
#   Ch1    Ch2    Ch3    Ch4    Ch5    Ch6    Ch7    Ch8 Unknow 
#  7346   6959   6943   6421   3491   2652   2438   1337    204

sce$cluster_colors <- Ch.anno_subset[match(sce@meta.data$cellid, Ch.anno_subset$cellid), "cluster_colors"]
sce$cluster_colors <- ifelse(is.na(sce$cluster_colors), "Unknow", sce$cluster_colors)
table(sce@meta.data$cluster_colors)
# #54181A #866262 #9D1F20 #A25B1D #B87A6D #CCB6B8 #E6773F #F2826C  Unknow 
#   3491    2438    6959    7346    6943    2652    6421    1337     204 

write.csv(sce@meta.data, file = "sce_meta.data_ch2.csv", row.names = FALSE)

table(sce@meta.data$cluster,sce@meta.data$sample)
   
#          HA1  HA2  HA3  HA4  HA5  OA1  OA2  OA3
#  Ch1     558 1179  351  732  381 2301  782 1062
#  Ch2     937 1063  725  680  852 1352  536  814
#  Ch3     723 1180  524  348  333 1897  975  963
#  Ch4     706  473 1246  361  630  596 1374 1035
#  Ch5     210  429  152  168  174 1292  690  376
#  Ch6     302  401  179  178  174  319  346  753
#  Ch7     443  239  628  106  151  264  316  291
#  Ch8     167  186  211   56  112  251  225  129
#  Unknow    4    3  187    0    0    3    0    7
#####unknow的细胞 公司删除？

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


library(harmony)
sce<- sce%>% harmony::RunHarmony("rawname")

# 设置颜色：指定HA和OA组的颜色
group_colors <- c("HA" = "#1f77b4", "OA" = "#ff7f0e")  # 你可以替换为任何你想要的颜色代码
png(file = "sample_harmony.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce, 
  group.by = "gname", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
  cols = group_colors  # 设置颜色
)
dev.off()

png(file = "rawname_harmony.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce, 
  group.by = "rawname", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
)
dev.off()

saveRDS(sce,file="sce_second.1.rds")

####筛选

png(file = "second_cluster.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce, 
  group.by = "cluster", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
)
dev.off()

#### PCA拐点定量识别
pct <- sce[["pca"]]@stdev / sum(sce[["pca"]]@stdev) * 100 ; 
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
sce <- RunUMAP(sce,reduction = "harmony", dims = 1:pcs)
sce <- RunTSNE(sce, dims.use = 1:pcs, reduction = "harmony",do.fast = TRUE)
sce <- FindNeighbors(sce, reduction = "harmony", dims = 1:pcs)
sce <- FindClusters(sce, algorithm = 1 , resolution = c(0.4),cluster.name="ch_0.4")

head(sce@meta.data)
write.csv(sce@meta.data, file = "sce_meta.data_ch_0.4.csv", row.names = FALSE)

png(file = "second_ch1.5.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce,
  reduction = "umap", 
  group.by = "ch_1.5", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
)
dev.off()
png(file = "second_ch0.4.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce,
  reduction = "umap", 
  group.by = "ch_0.4", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
)
dev.off()

png(file = "second_cluster_umap.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce,
  reduction = "umap", 
  group.by = "cluster", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
)
dev.off()

sce@meta.data <- sce@meta.data[, !(colnames(sce@meta.data) %in% c("ch_1", "ch_0.5"))]
# 查看删除后的元数据
head(sce@meta.data)

sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25,
                              thresh.use = 0.25)
write.csv(sce.markers,file=paste0('Ch_all.markers.csv'))
table(sce@meta.data$seurat_clusters,sce@meta.data$orig.ident)

top5 <- sce.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
write.csv(top5,file=paste0('Ch_top5.markers.csv'))
pdf(file="all.markers_top5_heatmap.pdf",width=30,height=10)
DoHeatmap(sce,top5$gene,size=3)
dev.off()
sce@meta.data %>% head()


########根据新的marker重新进行聚类

Ch1 = c("TF","FRZB","CYTL1","CNTFR","FOXA3")
Ch2 = c("CHI3L1", "CHI3L2","HMOX1","IFITM2","THBD","VCAM1","ICAM1")
Ch3 = c("LOX","LOXL2","LOXL3","CILP2","BMP2","SRPX2","SMOC1")
Ch4 = c("PDE4D","MT-CO1","MT-CO2","MT-CO3","MT-ND1","GLIS3")
Ch5 = c("IBSP","CNMD","COL10A1","SPP1","WWP2")
Ch6 = c("ID3","FOS","JUN","RGS16","ID1","TXNIP","KLF2","IER2")
Ch7 = c("COL1A1","COL1A2","ACAN","MMP2","MATN4","PTN")
Ch8 = c("TPPP3","S100A4","TMSB4X","CLEC3B","TXN","PRG4","IGFBP5")

cellMarker <- c(Ch1,Ch2,Ch3,Ch4,Ch5,Ch6,Ch7,Ch8)
source("/mnt/3w/data2/yiyuan/script/2024.11.work/DotPlot_fixed.R")
library(ggplot2)
library(Seurat)
library(ggsci)
library(cowplot)
######使用修改后的dotplot_fixed画气泡图
pdf(file="Dotplot_seurat_clusters_Ch0.4.pdf",width=60,height=8)
DotPlot_fixed(sce,group.by = 'seurat_clusters',cols = c('blue','red'),features=unique(cellMarker),cluster.idents = F)+
 scale_fill_gsea()+
 theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 15,
                                   angle = 90,
                                   hjust = 1,
                                   vjust = .5),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))
dev.off()
celltype=data.frame(ClusterID=0:7,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0),2]='Ch1'
celltype[celltype$ClusterID %in% c(1),2]='Ch2'
celltype[celltype$ClusterID %in% c(2),2]='Ch3'
celltype[celltype$ClusterID %in% c(3),2]='Ch4'
celltype[celltype$ClusterID %in% c(4),2]='Ch5'
celltype[celltype$ClusterID %in% c(5),2]='Ch6'
celltype[celltype$ClusterID %in% c(6),2]='Ch7'
celltype[celltype$ClusterID %in% c(7),2]='Ch8'
head(celltype)

sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)
# Ch1  Ch2  Ch3  Ch4  Ch5  Ch6  Ch7  Ch8 
#8042 6743 6109 4887 4824 3507 2070 1609 
sce$celltype <- factor(sce$celltype,levels=c("Ch1","Ch2","Ch3","Ch4","Ch5","Ch6","Ch7","Ch8"))
pdf(file="Dotplot_type_Ch.pdf",width=30,height=8)
DotPlot_fixed(sce,group.by = 'celltype',cols = c('blue','red'),features=unique(cellMarker),cluster.idents = F)+
 scale_fill_gsea()+
 theme(axis.title = element_blank(),
        axis.text.x = element_text(size = 15,
                                   angle = 90,
                                   hjust = 1,
                                   vjust = .5),
        axis.text.y = element_text(size = 15),
        legend.text = element_text(size = 15),
        legend.title = element_text(size = 15))
dev.off()

png(file = "second_celltype_umap.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce,
  reduction = "umap", 
  group.by = "celltype", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
   cols = c(
   'Ch1' = '#ef070f',  # 浅粉色
   'Ch2' = '#FAD02E',  # 浅黄色
   'Ch3' = '#A8DADC',  # 浅蓝绿色
   'Ch4' = '#aba3e7',  
   'Ch5' = '#712820',  # 浅蓝色
   'Ch6' = '#B7D7A8',  # 浅绿色
   'Ch7' = '#D4E157',  # 浅橄榄绿色
   'Ch8' = '#F2826C' 
))
dev.off()

table(sce@meta.data$celltype,sce@meta.data$sample)

#       HA1  HA2  HA3  HA4  HA5  OA1  OA2  OA3
#  Ch1  590 1465  396 1050  651 2235  568 1087
#  Ch2  764 1325  398  384  397 1485  859 1131
#  Ch3  189  405  282  110  135 2885 1731  372
#  Ch4  678  334 1106  154  221  404 1019  971
#  Ch5  787  529  964  521  877  484  167  495
#  Ch6  478  670  171  287  265  375  322  939
#  Ch7  347  186  533   74  129  204  329  268
#  Ch8  217  239  353   49  132  203  249  167

cell.Patient <- aggregate(sce$sample, list(sce$sample, sce$gname, sce$celltype), length)
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
   'Ch1' = '#ef070f',  # 浅粉色
   'Ch2' = '#FAD02E',  # 浅黄色
   'Ch3' = '#A8DADC',  # 浅蓝绿色
   'Ch4' = '#aba3e7',  
   'Ch5' = '#712820',  # 浅蓝色
   'Ch6' = '#B7D7A8',  # 浅绿色
   'Ch7' = '#D4E157',  # 浅橄榄绿色
   'Ch8' = '#F2826C' 
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



#分组
group_table <- as.data.frame(table(sce@meta.data$gname,sce@meta.data$celltype))
names(group_table) <- c("group","celltype","CellNumber")

# 计算每个组内的总和
group_table <- group_table %>%
  group_by(group) %>%
  mutate(total_cells = sum(CellNumber)) %>%
  ungroup()

# 计算每个细胞类型的占比
group_table <- group_table %>%
  mutate(percentage = CellNumber / total_cells)

pdf(file = "proportion_merged_group_celltype.pdf", width = 10, height = 12)
# 绘制按组的占比图
ggplot(group_table, aes(x = group, y = percentage, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +  # 使用填充的占比绘制柱状图
  scale_fill_manual(values = patient_col) +         # 使用手动颜色映射
  theme_minimal() +                                 # 使用简单的主题
  labs(x = "Group", y = "Percentage", fill = "Cell Type") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16))
dev.off()

library(tidyr)
library(reshape2)
sce$sample = factor(sce$sample,levels = c('OA1','OA2','OA3','HA1','HA2','HA3','HA4','HA5'))
tb=table(sce$sample, sce$celltype)
head(tb)
library (gplots)
png(file = "balloonplot_1_celltype.png", width = 3000, height = 2500, res = 300)
balloonplot(tb)
dev.off()

saveRDS(sce,file="sce_Ch_resubtype.final.rds")


setwd("/mnt/8w/data7/yiyuan/Bone/test/Ch/Re_subtype")

library(scRNAtoolVis)
Idents(sce1) <- "celltype"
cellMarker=c("TF","FRZB","CYTL1","CNTFR","FOXA3",
"CHI3L1","CHI3L2","HMOX1","IFITM2","THBD","VCAM1","ICAM1",
"LOX","LOXL2","LOXL3","CILP2","BMP2","SRPX2","SMOC1",
"PDE4D","MT-CO1","MT-CO2","MT-CO3","MT-ND1","GLIS3",
"IBSP","CNMD","COL10A1","SPP1","WWP2",
"ID3","FOS","JUN","RGS16","ID1","TXNIP","KLF2","IER2",
"COL1A1","COL1A2","ACAN","MMP2","MATN4","PTN",
"TPPP3","S100A4","TMSB4X","CLEC3B","TXN","PRG4","IGFBP5")
sce1 = readRDS("/mnt/8w/data7/yiyuan/Bone/test/Ch/Re_subtype/sce_Ch_resubtype.final.rds")
cellMarker=c("TF","FRZB","CYTL1","CNTFR","FOXA3","CHI3L1","CHI3L2","HMOX1","IFITM2","THBD","VCAM1","ICAM1",
"LOX","LOXL2","LOXL3","CILP2","BMP2","SRPX2","SMOC1",
"PDE4D","MT-CO1","MT-CO2","MT-CO3","MT-ND1","GLIS3",
"IBSP","CNMD","COL10A1","SPP1","WWP2",
"ID3","FOS","JUN","RGS16","ID1","TXNIP","KLF2","IER2",
"COL1A1","COL1A2","ACAN","MMP2","MATN4","PTN","TPPP3","S100A4","TMSB4X","TXN","PRG4","IGFBP5")###CLEC3B基因报错
png(file = "heatmap_Ch_celltype_marker.png", width = 3700, height = 4200, res = 300)
averageHeatmap(object = sce,
               markerGene = cellMarker,
               group.by = 'celltype',
               gene.order = cellMarker,
               column_split = 1:8,#分割列
               border = T)
dev.off()
#[1] "Your cluster annotation color is:" "#B969E5FF"                        
#[3] "#BD94EDFF"                         "#E6C0FBFF"                        
#[5] "#E477F5FF"                         "#FFC8CEFF"                        
#[7] "#B8FD8AFF"                         "#7CF8E1FF"                        
#[9] "#74F6F1FF" 
png(file = "umap_Ch_celltype_TEST.png", width = 3700, height = 3500, res = 300)
# legend key size
# add cell type
clusterCornerAxes(object = sce,
                  reduction = 'umap',
                  clusterCol = "celltype",
                  noSplit = T,
                  cellLabel = T,
                  cellLabelSize = 10)
dev.off()


######绘制软骨HA和OA的umap
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
setwd("/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Ch/second")
sce1 = readRDS("/mnt/8w/data7/yiyuan/Bone/test/Ch/Re_subtype/sce_Ch_resubtype.final.rds")
head(sce1@meta.data)

HA <- sce1[, sce1$gname %in% c("HA")]
HA
OA <- sce1[, sce1$gname %in% c("OA")]
OA
umap = sce1@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = sce1@meta.data$celltype) # 注释后的label信息 ，改为cell_type
###"#F08080","#1E90FF","#7CFC00","#FFFF00", "#808000","#FF00FF","#FA8072","#7B68EE"
allcolour=c( "#a4cde1","#67a4cc","#277fb8","#549da3","#96cb8f","#8bc96d","#4dae47","#5c9e43") #注意：颜色数量一顶要大于细胞类型数量
colnames(umap)[1:2]=c("UMAP_1","UMAP_2")
umap_data=umap
p <- ggplot(umap_data,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  geom_point(size = 0.5 , alpha =1 )  +  scale_color_manual(values = allcolour)
p21= p  + theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))


p31 <- p21 +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) #设置legend中 点的大小 

p41 <- p31 + 
  geom_segment(aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,
                   xend = min(umap$UMAP_1) +3, yend = min(umap$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                   xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$UMAP_1) +1.5, y = min(umap$UMAP_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
png(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Ch/2_28/Ch_second_celltype_umap.png", width = 3000, height = 2500, res = 300)
p31
dev.off()
umap = HA@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = HA@meta.data$celltype) # 注释后的label信息 ，改为cell_type
###"#F08080","#1E90FF","#7CFC00","#FFFF00", "#808000","#FF00FF","#FA8072","#7B68EE"
allcolour=c( "#a4cde1","#67a4cc","#277fb8","#549da3","#96cb8f","#8bc96d","#4dae47","#5c9e43" ) #注意：颜色数量一顶要大于细胞类型数量
colnames(umap)[1:2]=c("UMAP_1","UMAP_2")
umap_data=umap
p <- ggplot(umap_data,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  geom_point(size = 0.5 , alpha =1 )  +  scale_color_manual(values = allcolour)
p21= p  + theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))


p31 <- p21 +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) #设置legend中 点的大小 

p41 <- p31 + 
  geom_segment(aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,
                   xend = min(umap$UMAP_1) +3, yend = min(umap$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                   xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$UMAP_1) +1.5, y = min(umap$UMAP_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
png(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Ch/2_28/Ch_second_celltype_umap_HA.png", width = 3000, height = 2500, res = 300)
p31
dev.off()

umap = OA@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = OA@meta.data$celltype) # 注释后的label信息 ，改为cell_type
###"#F08080","#1E90FF","#7CFC00","#FFFF00", "#808000","#FF00FF","#FA8072","#7B68EE"
allcolour=c(  "#a4cde1","#67a4cc","#277fb8","#549da3","#96cb8f","#8bc96d","#4dae47","#5c9e43"  ) #注意：颜色数量一顶要大于细胞类型数量
colnames(umap)[1:2]=c("UMAP_1","UMAP_2")
umap_data=umap
p <- ggplot(umap_data,aes(x= UMAP_1 , y = UMAP_2 ,color = cell_type)) +  geom_point(size = 0.5 , alpha =1 )  +  scale_color_manual(values = allcolour)
p21= p  + theme(panel.grid.major = element_blank(), #主网格线
        panel.grid.minor = element_blank(), #次网格线
        panel.border = element_blank(), #边框
        axis.title = element_blank(),  #轴标题
        axis.text = element_blank(), # 文本
        axis.ticks = element_blank(),
        panel.background = element_rect(fill = 'white'), #背景色
        plot.background=element_rect(fill="white"))


p31 <- p21 +         
  theme(
    legend.title = element_blank(), #去掉legend.title 
    legend.key=element_rect(fill='white'), #
    legend.text = element_text(size=20), #设置legend标签的大小
    legend.key.size=unit(1,'cm') ) +  # 设置legend标签之间的大小
  guides(color = guide_legend(override.aes = list(size=5))) #设置legend中 点的大小 

p41 <- p31 + 
  geom_segment(aes(x = min(umap$UMAP_1) , y = min(umap$UMAP_2) ,
                   xend = min(umap$UMAP_1) +3, yend = min(umap$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                   xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 3),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$UMAP_1) +1.5, y = min(umap$UMAP_2) -1, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap$UMAP_1) -1, y = min(umap$UMAP_2) + 1.5, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
png(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Ch/2_28/Ch_second_celltype_umap_OA.png", width = 3000, height = 2500, res = 300)
p31
dev.off()


patient_col = c(  "#a4cde1","#67a4cc","#277fb8","#549da3","#96cb8f","#8bc96d","#4dae47","#5c9e43" )
#分组
group_table <- as.data.frame(table(HA@meta.data$rawname,HA@meta.data$celltype))
names(group_table) <- c("group","celltype","CellNumber")
table(HA$rawname,HA$celltype)
#             Ch1  Ch2  Ch3  Ch4  Ch5  Ch6  Ch7  Ch8
#  HAC_0922  1050  384  110  154  521  287   74   49
#  HAC_0928   651  397  135  221  877  265  129  132
#  HAC1_1130  590  764  189  678  787  478  347  217
#  HAC2_1226 1465 1325  405  334  529  670  186  239
#  HAC3_0228  396  398  282 1106  964  171  533  353
# 计算每个组内的总和
group_table <- group_table %>%
  group_by(group) %>%
  mutate(total_cells = sum(CellNumber)) %>%
  ungroup()

# 计算每个细胞类型的占比
group_table <- group_table %>%
  mutate(percentage = CellNumber / total_cells)

pdf(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Ch/2_28/second_merged_Sample_celltype_HA.pdf", width = 10, height = 12)
# 绘制按组的占比图
ggplot(group_table, aes(x = group, y = percentage, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +  # 使用填充的占比绘制柱状图
  scale_fill_manual(values = patient_col) +         # 使用手动颜色映射
  theme_minimal() +                                 # 使用简单的主题
  labs(x = "Group", y = "Percentage", fill = "Cell Type") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16))
dev.off()
#分组
group_table <- as.data.frame(table(OA@meta.data$rawname,OA@meta.data$celltype))
names(group_table) <- c("group","celltype","CellNumber")
table(OA$rawname,OA$celltype)
#            Ch1  Ch2  Ch3  Ch4  Ch5  Ch6  Ch7  Ch8
#  OA1_1125 2235 1485 2885  404  484  375  204  203
#  OA2_1127  568  859 1731 1019  167  322  329  249
#  OA3_1130 1087 1131  372  971  495  939  268  167

# 计算每个组内的总和
group_table <- group_table %>%
  group_by(group) %>%
  mutate(total_cells = sum(CellNumber)) %>%
  ungroup()

# 计算每个细胞类型的占比
group_table <- group_table %>%
  mutate(percentage = CellNumber / total_cells)

pdf(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Ch/2_28/second_merged_Sample_celltype_OA.pdf", width = 10, height = 12)
# 绘制按组的占比图
ggplot(group_table, aes(x = group, y = percentage, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +  # 使用填充的占比绘制柱状图
  scale_fill_manual(values = patient_col) +         # 使用手动颜色映射
  theme_minimal() +                                 # 使用简单的主题
  labs(x = "Group", y = "Percentage", fill = "Cell Type") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16))
dev.off()

group_table <- as.data.frame(table(sce1@meta.data$rawname,sce1@meta.data$celltype))
names(group_table) <- c("group","celltype","CellNumber")
table(sce1$rawname,sce1$celltype)
#             Ch1  Ch2  Ch3  Ch4  Ch5  Ch6  Ch7  Ch8
#  HAC_0922  1050  384  110  154  521  287   74   49
#  HAC_0928   651  397  135  221  877  265  129  132
#  HAC1_1130  590  764  189  678  787  478  347  217
#  HAC2_1226 1465 1325  405  334  529  670  186  239
#  HAC3_0228  396  398  282 1106  964  171  533  353
#  OA1_1125  2235 1485 2885  404  484  375  204  203
#  OA2_1127   568  859 1731 1019  167  322  329  249
#  OA3_1130  1087 1131  372  971  495  939  268  167
group_table <- group_table %>%
  group_by(group) %>%
  mutate(total_cells = sum(CellNumber)) %>%
  ungroup()

# 计算每个细胞类型的占比
group_table <- group_table %>%
  mutate(percentage = CellNumber / total_cells)

pdf(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Ch/2_28/second_merged_Sample_celltype.pdf", width = 10, height = 12)
# 绘制按组的占比图
ggplot(group_table, aes(x = group, y = percentage, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +  # 使用填充的占比绘制柱状图
  scale_fill_manual(values = patient_col) +         # 使用手动颜色映射
  theme_minimal() +                                 # 使用简单的主题
  labs(x = "Group", y = "Percentage", fill = "Cell Type") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16))
dev.off()

group_table <- as.data.frame(table(sce1@meta.data$gname,sce1@meta.data$celltype))
names(group_table) <- c("group","celltype","CellNumber")
table(sce1$gname,sce1$celltype)
#      Ch1  Ch2  Ch3  Ch4  Ch5  Ch6  Ch7  Ch8
#  HA 4152 3268 1121 2493 3678 1871 1269  990
#  OA 3890 3475 4988 2394 1146 1636  801  619
group_table <- group_table %>%
  group_by(group) %>%
  mutate(total_cells = sum(CellNumber)) %>%
  ungroup()

# 计算每个细胞类型的占比
group_table <- group_table %>%
  mutate(percentage = CellNumber / total_cells)

pdf(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Ch/2_28/second_merged_group_celltype.pdf", width = 10, height = 12)
# 绘制按组的占比图
ggplot(group_table, aes(x = group, y = percentage, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +  # 使用填充的占比绘制柱状图
  scale_fill_manual(values = patient_col) +         # 使用手动颜色映射
  theme_minimal() +                                 # 使用简单的主题
  labs(x = "Group", y = "Percentage", fill = "Cell Type") +
  theme(panel.grid = element_blank(),
        panel.background = element_rect(fill = "transparent",colour = NA),
        axis.line.x = element_line(colour = "black") ,
        axis.line.y = element_line(colour = "black") ,
        plot.title = element_text(lineheight=.8, face="bold", hjust=0.5, size =16))
dev.off()


####
##########关键的基因绘制
setwd("/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Ch/new/marker_gene")
library(scRNAtoolVis)
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
sce = readRDS("/mnt/8w/data7/yiyuan/Bone/test/Ch/Re_subtype/sce_Ch_resubtype.final.rds")
###1，2 linging ####红色的
png(file = "Feature_Ch软骨.png", width = 3700, height = 3500, res = 300)
# no facet group
featureCornerAxes(object = sce,reduction = 'umap',
                  groupFacet = NULL,
                  relLength = 0.5,
                  relDist = 0.2,
                  features = c("JUN","IBSP","COL1A1","FOS","PRG4","CILP"),
                  aspect.ratio = 1,themebg = 'bwCorner',minExp = 0,maxExp = 6)
dev.off()