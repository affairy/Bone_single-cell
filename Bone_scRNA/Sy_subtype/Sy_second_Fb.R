setwd("/mnt/8w/data7/yiyuan/Bone/test1/Fibroblasts")
library(Seurat)
library(Matrix)
library(harmony)
library(ggplot2)
library(tidyverse)
###先挑选软骨
sce_subset <- sce1[, sce1$celltype_cluster_standard %in% c("Fibroblasts")]
sce_subset
#An object of class Seurat 
#69378 features across 50346 samples within 2 assays 
#Active assay: origin.counts (34689 features, 2000 variable features)
# 3 layers present: data, counts, scale.data
# 1 other assay present: RNA
# 4 dimensional reductions calculated: pca, harmony, umap, tsne
head(sce_subset@meta.data)
table(sce_subset@meta.data$celltype_cluster_standard,sce_subset@meta.data$sample)
#              
#                HA1   HA2   HA3   OA1   OA2   OA3
#  Fibroblasts  7895  6258  4927 11214  8851 11201   50346
####和公司二级注释多了些细胞 
#	OA1	OA2	OA3	HA1	HA2	HA3	all
#Fibroblasts	10441	8199	10680	7539	6007	4474	47340

write.csv(sce_subset@meta.data, file = "sce_subset_meta.data_Fibroblasts.csv", row.names = FALSE)

saveRDS(sce_subset,file="/mnt/8w/data7/yiyuan/Bone/test1/Fibroblasts/sce_first.Fb.rds")

sce_subset = readRDS('/mnt/8w/data7/yiyuan/Bone/test1/Fibroblasts/sce_first.Fb.rds')

sce=sce_subset


Fb.anno <- data.table::fread("/mnt/8w/data7/yiyuan/Bone/anno/31_P22112601_B1_滑膜_cellinfo/Fibroblasts_3_color/P22112601_B1/P22112601_B1_cellinfo.xls", data.table=F)
dim(Fb.anno)
#[1] 47340    24
colnames(Fb.anno)[1] <- "cellid"
# 检查修改后的列名
head(Fb.anno)
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
#2      6.611940       2           3       #0067AA      #0067AA
#3     15.925783       0           1       #0067AA      #0067AA
#4      6.654991       0           1       #0067AA      #0067AA
#5      6.138160       2           3       #0067AA      #0067AA
#6      6.233062       3           4       #0067AA      #0067AA
#  raw_cluster_colors cluster_standard     S_score   G2M_score phase     cluster
#1            #FF1F1D         Lining_2  0.14940149 -0.29044118     S    Lining_2
#2            #00A23F      Sublining_2 -0.06273063 -0.21813725    G1 Sublining_2
#3            #0067AA      Sublining_1  0.04857049 -0.07965686     S Sublining_1
#4            #0067AA      Sublining_1 -0.03936039 -0.08210784    G1 Sublining_1
#5            #00A23F      Sublining_2  0.03996040 -0.12377451     S Sublining_2
#6            #FF1F1D         Lining_2  0.04365044 -0.06740196     S    Lining_2
#  cluster_colors annot_modified
#1        #E8B975       Lining_2
#2        #889658    Sublining_2
#3        #EEDF7A    Sublining_1
#4        #EEDF7A    Sublining_1
#5        #889658    Sublining_2
#6        #E8B975       Lining_2
subset_cols <- c("cellid", "rawname", "gname","sample_colors","cluster","cluster_colors")
# 提取指定列
Fb.anno_subset <- Fb.anno %>% select(all_of(subset_cols))
dim(Fb.anno_subset)
#[1] 47340     6
head(Fb.anno_subset)
#                           cellid   rawname gname sample_colors     cluster
#1 OA1_AACACACAGAACGCTAGTATGAGTCCT OAS1_1209    OA       #0067AA    Lining_2
#2 OA1_AACACACAGAACGCTAGTGAGTTCGTC OAS1_1209    OA       #0067AA Sublining_2
#3 OA1_AACACACAGAACGTCCAACTCTGCTTA OAS1_1209    OA       #0067AA Sublining_1
#4 OA1_AACACACAGAAGGTGGTACTCTGCTTA OAS1_1209    OA       #0067AA Sublining_1
#5 OA1_AACACACAGAAGGTGGTATCCAGACGA OAS1_1209    OA       #0067AA Sublining_2
#6 OA1_AACACACAGAATCCGGTGATTAGCCTG OAS1_1209    OA       #0067AA    Lining_2
#  cluster_colors
#1        #E8B975
#2        #889658
#3        #EEDF7A
#4        #EEDF7A
#5        #889658
#6        #E8B975
sce$cluster <- Fb.anno_subset[match(sce@meta.data$cellid, Fb.anno_subset$cellid), "cluster"]
sce$cluster <- ifelse(is.na(sce$cluster), "Unknow", sce$cluster)
table(sce@meta.data$cluster)
#
#   Lining_1    Lining_2 Sublining_1 Sublining_2 Sublining_3 Sublining_4 
#      10103        7254       12492        7463        5520        4508 
#     Unknow 
#       3006

sce$cluster_colors <- Fb.anno_subset[match(sce@meta.data$cellid, Fb.anno_subset$cellid), "cluster_colors"]
sce$cluster_colors <- ifelse(is.na(sce$cluster_colors), "Unknow", sce$cluster_colors)
table(sce@meta.data$cluster_colors)
# #5A5E09 #889658 #A09623 #CAC6AB #E8B975 #EEDF7A  Unknow 
#   5520    7463   10103    4508    7254   12492    3006 

write.csv(sce@meta.data, file = "sce_meta.data_Fb2.csv", row.names = FALSE)

table(sce@meta.data$cluster,sce@meta.data$sample)
#             
#               HA1  HA2  HA3  OA1  OA2  OA3
#  Lining_1    1273  579  382 1738 1486 4645
#  Lining_2    1427  423  544 1210 1304 2346
#  Sublining_1 2704 1591 1755 2660 2482 1300
#  Sublining_2  712  984  747 1911 1313 1796
#  Sublining_3  847 1263  499 1674  795  442
#  Sublining_4  576 1167  547 1248  819  151
#  Unknow       356  251  453  773  652  521

sce <- NormalizeData(sce, normalization.method =  "LogNormalize",
                     scale.factor = 10000)
#GetAssay(sce,assay = "RNA")
sce@meta.data %>% head()
#                                        orig.ident nCount_RNA nFeature_RNA
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG       OAS1       7711         2569
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG       OAS1       6996         2175
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG       OAS1       9603         2334
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG       OAS1       7622         2164
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG       OAS1       5099         1640
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG       OAS1       6043         1968
#                                                           barcodes
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG TACCTCTCCAACGTCCAAAACAAGTGG
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG GCGAGTAACAAGGTGGTAAACAAGTGG
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG CTTACGCAGACAACAGGTAACAAGTGG
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG TGCCTGATCACCAGGTCAAACAAGTGG
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG ATGGTCTCAACCGTACTCAACAAGTGG
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG CAGTCTTCGACCGTACTCAACAAGTGG
#                                                                 cellid
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG OA1_TACCTCTCCAACGTCCAAAACAAGTGG
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG OA1_GCGAGTAACAAGGTGGTAAACAAGTGG
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG OA1_CTTACGCAGACAACAGGTAACAAGTGG
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG OA1_TGCCTGATCACCAGGTCAAACAAGTGG
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG OA1_ATGGTCTCAACCGTACTCAACAAGTGG
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG OA1_CAGTCTTCGACCGTACTCAACAAGTGG
#                                          rawname sample celltype_cluster
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG OAS1_1209    OA1      Fibroblasts
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG OAS1_1209    OA1      Fibroblasts
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG OAS1_1209    OA1      Fibroblasts
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG OAS1_1209    OA1      Fibroblasts
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG OAS1_1209    OA1      Fibroblasts
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG OAS1_1209    OA1      Fibroblasts
#                                        celltype_cluster_standard
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG               Fibroblasts
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG               Fibroblasts
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG               Fibroblasts
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG               Fibroblasts
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG               Fibroblasts
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG               Fibroblasts
#                                        cluster_colors gname percent.mt
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG        #5A5E09    OA   6.471275
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG        #889658    OA   4.016581
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG        #E8B975    OA  13.183380
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG        #E8B975    OA  11.886644
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG        #E8B975    OA   8.080016
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG        #889658    OA   5.096806
#                                        origin.counts_snn_res.1 seurat_clusters
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG                       7               7
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG                       2               2
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG                       3               3
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG                       3               3
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG                       3               3
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG                       2               2
#                                            cluster
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG Sublining_3
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG Sublining_2
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG    Lining_2
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG    Lining_2
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG    Lining_2
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG Sublining_2

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
sce <- FindClusters(sce, algorithm = 1 , resolution = c(0.5),cluster.name="Fb_0.5")

head(sce@meta.data)
write.csv(sce@meta.data, file = "sce_meta.data_Fb_0.4.csv", row.names = FALSE)
write.csv(sce@meta.data, file = "sce_meta.data_Fb_0.5.csv", row.names = FALSE)

png(file = "second_Fb0.4.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce,
  reduction = "umap", 
  group.by = "Fb_0.4", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
)
dev.off()

png(file = "second_Fb0.5.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce,
  reduction = "umap", 
  group.by = "Fb_0.5", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
)
dev.off()


png(file = "second_Fb0.4_cluster_umap.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce,
  reduction = "umap", 
  group.by = "cluster", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
)
dev.off()

sce.markers <- FindAllMarkers(object = sce, only.pos = TRUE, min.pct = 0.25,
                              thresh.use = 0.25)
write.csv(sce.markers,file=paste0('Fb_all.markers.csv'))
table(sce@meta.data$seurat_clusters,sce@meta.data$orig.ident)
#    HAS1 HAS2 HAS3 OAS1 OAS2 OAS3
#  0 2257 1353 1398 2565 2116  818
#  1 1772  498  526 1481 1401 3420
#  2  817  460  216 1491 1007 3550
#  3  565 1089  737 1859 1270 1880
#  4  618 1129  384 1473  794  518
#  5 1124  433  729  860  835  826
#  6  362  748  524  505  453   27
#  7  144  254  234  775  862  130
#  8  236  294  179  205  113   32
top5 <- sce.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
write.csv(top5,file=paste0('Fb_top5.markers.csv'))
pdf(file="all.markers_top5_heatmap.pdf",width=30,height=10)
DoHeatmap(sce,top5$gene,size=3)
dev.off()
sce@meta.data %>% head()
#                                        orig.ident nCount_RNA nFeature_RNA
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG       OAS1       7711         2569
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG       OAS1       6996         2175
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG       OAS1       9603         2334
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG       OAS1       7622         2164
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG       OAS1       5099         1640
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG       OAS1       6043         1968
#                                                           barcodes
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG TACCTCTCCAACGTCCAAAACAAGTGG
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG GCGAGTAACAAGGTGGTAAACAAGTGG
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG CTTACGCAGACAACAGGTAACAAGTGG
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG TGCCTGATCACCAGGTCAAACAAGTGG
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG ATGGTCTCAACCGTACTCAACAAGTGG
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG CAGTCTTCGACCGTACTCAACAAGTGG
#                                                                 cellid
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG OA1_TACCTCTCCAACGTCCAAAACAAGTGG
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG OA1_GCGAGTAACAAGGTGGTAAACAAGTGG
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG OA1_CTTACGCAGACAACAGGTAACAAGTGG
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG OA1_TGCCTGATCACCAGGTCAAACAAGTGG
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG OA1_ATGGTCTCAACCGTACTCAACAAGTGG
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG OA1_CAGTCTTCGACCGTACTCAACAAGTGG
#                                          rawname sample celltype_cluster
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG OAS1_1209    OA1      Fibroblasts
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG OAS1_1209    OA1      Fibroblasts
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG OAS1_1209    OA1      Fibroblasts
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG OAS1_1209    OA1      Fibroblasts
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG OAS1_1209    OA1      Fibroblasts
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG OAS1_1209    OA1      Fibroblasts
#                                        celltype_cluster_standard
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG               Fibroblasts
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG               Fibroblasts
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG               Fibroblasts
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG               Fibroblasts
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG               Fibroblasts
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG               Fibroblasts
#                                        cluster_colors gname percent.mt
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG        #5A5E09    OA   6.471275
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG        #889658    OA   4.016581
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG        #E8B975    OA  13.183380
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG        #E8B975    OA  11.886644
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG        #E8B975    OA   8.080016
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG        #889658    OA   5.096806
#                                        origin.counts_snn_res.1 seurat_clusters
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG                       7               4
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG                       2               3
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG                       3               1
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG                       3               1
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG                       3               1
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG                       2               3
#                                            cluster Fb_0.4 Fb_0.5
#OAS1_1209_TACCTCTCC_AACGTCCAA_AACAAGTGG Sublining_3      4      4
#OAS1_1209_GCGAGTAAC_AAGGTGGTA_AACAAGTGG Sublining_2      1      3
#OAS1_1209_CTTACGCAG_ACAACAGGT_AACAAGTGG    Lining_2      3      1
#OAS1_1209_TGCCTGATC_ACCAGGTCA_AACAAGTGG    Lining_2      3      1
#OAS1_1209_ATGGTCTCA_ACCGTACTC_AACAAGTGG    Lining_2      3      1
#OAS1_1209_CAGTCTTCG_ACCGTACTC_AACAAGTGG Sublining_2      1      3
table(sce@meta.data$seurat_clusters,sce@meta.data$sample)
#   
#     HA1  HA2  HA3  OA1  OA2  OA3
#  0 2257 1353 1398 2565 2116  818
#  1 1772  498  526 1481 1401 3420
#  2  817  460  216 1491 1007 3550
#  3  565 1089  737 1859 1270 1880
#  4  618 1129  384 1473  794  518
#  5 1124  433  729  860  835  826
#  6  362  748  524  505  453   27
#  7  144  254  234  775  862  130
#  8  236  294  179  205  113   32

pdf(file="FeaturePlot_PDPN.pdf",width = 10, height = 10)
FeaturePlot(sce,features = c("PDPN"),raster = FALSE,reduction = 'umap')
dev.off()

pdf(file="FeaturePlot_THY1.pdf",width = 10, height = 10)
FeaturePlot(sce,features = c("THY1"),raster = FALSE,reduction = 'umap')
dev.off()

####新marker的绘制
###3，5，7合并
Fb1 = c("C3","CD34","CDH11","THY1","OLFML3","SRPX","IGF1","FGF7","MARCKS","FBLN1","SFRP1")
Fb2 = c("CLIC5","DEFB1","CD55","HBEGF","GPR1","CD55","THBS4","LTBP4","PRG4")
Fb3 = c("RGS16","HAS1","BTG2","IRF1","CXCL1")
Fb4 = c("CCL2","CXCL12","FOSB","NR4A1","NR4A2","KLF2","KLF4","ICAM1","APOC1","APOE","LRP1B","KAZN","MAGI2","MALAT1","NEAT1")
Fb5 = c("DKK3","TREM1","CRIP1","CRIP2","TPPP3","S100A4")
Fb6 = c("SFRP2","THY1","CYGB","SFRP4","VCAN","NOTCH3")
Fb7 = c("COL1A1","CADM1","MDK","MMP2","COL1A2")

cellMarker <- c(Fb1,Fb2,Fb3,Fb4,Fb5,Fb6,Fb7)
source("/mnt/3w/data2/yiyuan/script/2024.11.work/DotPlot_fixed.R")
library(ggplot2)
library(Seurat)
library(ggsci)
library(cowplot)
######使用修改后的dotplot_fixed画气泡图
pdf(file="Dotplot_seurat_clusters_Fb0.5.pdf",width=60,height=8)
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

celltype=data.frame(ClusterID=0:8,
                    celltype='unkown')
celltype[celltype$ClusterID %in% c(0),2]='Fb1'
celltype[celltype$ClusterID %in% c(1),2]='Fb2'
celltype[celltype$ClusterID %in% c(2),2]='Fb3'
celltype[celltype$ClusterID %in% c(3,5,7),2]='Fb4'
celltype[celltype$ClusterID %in% c(4),2]='Fb5'
celltype[celltype$ClusterID %in% c(6),2]='Fb6'
celltype[celltype$ClusterID %in% c(8),2]='Fb7'
#celltype[celltype$ClusterID %in% c(7),2]='Unknow'
head(celltype)

sce@meta.data$celltype = "NA"
for(i in 1:nrow(celltype)){
  sce@meta.data[which(sce@meta.data$seurat_clusters == celltype$ClusterID[i]),'celltype'] <- celltype$celltype[i]}
table(sce@meta.data$celltype)
#   Fb1    Fb2    Fb3    Fb4    Fb5    Fb6    Fb7 
#10507  9098  7541 14606  4916  2619  1059 
sce$celltype <- factor(sce$celltype,levels=c("Fb1","Fb2","Fb3","Fb4","Fb5","Fb6","Fb7"))
pdf(file="Dotplot_type_Fb.pdf",width=30,height=8)
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


png(file = "second_celltype_umap_new.png", width = 3000, height = 2500, res = 300)
DimPlot(
  sce,
  reduction = "umap", 
  group.by = "celltype", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
   cols = c(
   'Fb1' = '#ef070f',  # 浅粉色
   'Fb2' = '#FAD02E',  # 浅黄色
   'Fb3' = '#A8DADC',  # 浅蓝绿色
   'Fb4' = '#aba3e7',  
   'Fb5' = '#D4E157',  # 浅橄榄绿色
   'Fb6' = '#F2826C',
   'Fb7' = '#33A02CFF'
))
dev.off()

table(sce@meta.data$celltype,sce@meta.data$sample)
#       OA1  OA2  OA3  HA1  HA2  HA3
#  Fb1 2565 2116  818 2257 1353 1398
#  Fb2 1481 1401 3420 1772  498  526
#  Fb3 1491 1007 3550  817  460  216
#  Fb4 3494 2967 2836 1833 1776 1700
#  Fb5 1473  794  518  618 1129  384
#  Fb6  505  453   27  362  748  524
#  Fb7  205  113   32  236  294  179

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
   'Fb1' = '#ef070f',  # 浅粉色
   'Fb2' = '#FAD02E',  # 浅黄色
   'Fb3' = '#A8DADC',  # 浅蓝绿色
   'Fb4' = '#aba3e7',  
   'Fb5' = '#D4E157',  # 浅橄榄绿色
   'Fb6' = '#F2826C',
   'Fb7' = '#33A02CFF'
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

patient_col = c(
   'Fb1' = '#ef070f',  # 浅粉色
   'Fb2' = '#FAD02E',  # 浅黄色
   'Fb3' = '#A8DADC',  # 浅蓝绿色
   'Fb4' = '#aba3e7',  
   'Fb5' = '#D4E157',  # 浅橄榄绿色
   'Fb6' = '#F2826C',
   'Fb7' = '#33A02CFF'
)

pdf(file = "proportion_gname.pdf", width = 10, height = 12)
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
sce$sample = factor(sce$sample,levels = c('OA1','OA2','OA3','HA1','HA2','HA3'))
tb=table(sce$sample, sce$celltype)
head(tb)
library (gplots)

png(file = "balloonplot_1_celltype.png", width = 4000, height = 2500, res = 300)
balloonplot(tb)
dev.off()

saveRDS(sce,file="sce_Fb_resubtype.final.rds")

#############
###8W运行
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
setwd("/mnt/data7/yiyuan/Bone/test1/Fibroblasts/new/DEGs")
sce1 = readRDS("/mnt/data7/yiyuan/Bone/test1/Fibroblasts/new/sce_Fb_resubtype.final.rds")
head(sce1@meta.data)
Idents(sce1)="celltype"

sce.markers <- FindAllMarkers(object = sce1, only.pos = TRUE, min.pct = 0.25,
                              thresh.use = 0.25)
write.csv(sce.markers,file=paste0('Fb滑膜_all.markers_celltype.csv'))
table(sce1@meta.data$seurat_clusters,sce1@meta.data$orig.ident)

top5 <- sce.markers %>% group_by(cluster) %>% top_n(5, avg_log2FC)
write.csv(top5,file=paste0('Fb滑膜_top5.markers_celltype.csv'))

pdf(file="all.markers_top5_heatmap_Fb滑膜_celltype.pdf",width=30,height=10)
DoHeatmap(sce1,top5$gene,size=6)
dev.off()
sce1@meta.data %>% head()


####绘制热图
setwd("/mnt/8w/data7/yiyuan/Bone/test1/Fibroblasts/new")
library(scRNAtoolVis)
sce = readRDS("/mnt/8w/data7/yiyuan/Bone/test1/Fibroblasts/new/sce_Fb_resubtype.final.rds")
Fb1 = c("C3","CD34","CDH11","THY1","OLFML3","SRPX","IGF1","FGF7","MARCKS","FBLN1","SFRP1")
Fb2 = c("CLIC5","DEFB1","CD55","HBEGF","GPR1","CD55","THBS4","LTBP4","PRG4")
Fb3 = c("RGS16","HAS1","BTG2","IRF1","CXCL1")
Fb4 = c("CCL2","CXCL12","FOSB","NR4A1","NR4A2","KLF2","KLF4","ICAM1","APOC1","APOE","LRP1B","KAZN","MAGI2","MALAT1","NEAT1")
Fb5 = c("DKK3","TREM1","CRIP1","CRIP2","TPPP3","S100A4")
Fb6 = c("SFRP2","THY1","CYGB","SFRP4","VCAN","NOTCH3")
Fb7 = c("COL1A1","CADM1","MDK","MMP2","COL1A2")

cellMarker <- c(Fb1,Fb2,Fb3,Fb4,Fb5,Fb6,Fb7)
cellMarker=c("C3","CD34","CDH11","THY1","OLFML3","SRPX","IGF1","FGF7","MARCKS","FBLN1","SFRP1",
"CLIC5","DEFB1","CD55","HBEGF","GPR1","CD55","THBS4","LTBP4","PRG4",
"RGS16","HAS1","BTG2","IRF1","CXCL1",
"FOSB","NR4A1","NR4A2","KLF2","KLF4","ICAM1","CCL2","CXCL12","APOC1","APOE","LRP1B","KAZN","MAGI2","MALAT1","NEAT1",
"DKK3","TREM1","CRIP1","CRIP2","TPPP3","S100A4",
"SFRP2","THY1","CYGB","SFRP4","VCAN","NOTCH3",
"COL1A1","CADM1","MDK","MMP2","COL1A2")

Idents(sce) <- "celltype"
png(file = "heatmap_Fb滑膜_celltype_marker1.png", width = 3700, height = 4200, res = 300)
averageHeatmap(object = sce,
               markerGene = cellMarker,
               gene.order = cellMarker,
               column_split = 1:7,#分割列
               border = T)
dev.off()

##########关键的基因绘制
###1，2 linging ####红色的
png(file = "Feature_Fb滑膜.png", width = 3700, height = 3500, res = 300)
# no facet group
featureCornerAxes(object = sce,reduction = 'umap',
                  groupFacet = NULL,
                  relLength = 0.5,
                  relDist = 0.2,
                  features = c("C3","CD34","CDH11","THY1","PRG4"),
                  aspect.ratio = 1,minExp = 0,maxExp = 6)
dev.off()
png(file = "Feature_Fb滑膜_sublining.png", width = 3700, height = 3500, res = 300)
# no facet group
featureCornerAxes(object = sce,reduction = 'umap',
                  groupFacet = NULL,
                  relLength = 0.5,
                  relDist = 0.2,
                  features = c("C3","CD34","CDH11","THY1"),
                  aspect.ratio = 1,minExp = 0,maxExp = 6)
dev.off()

png(file = "Feature_Fb滑膜_lining.png", width = 3700, height = 3500, res = 300)
# no facet group
featureCornerAxes(object = sce,reduction = 'umap',
                  groupFacet = NULL,
                  relLength = 0.5,
                  relDist = 0.2,
                  features = c("CLIC5","DEFB1","CD55","HBEGF","PRG4"),
                  aspect.ratio = 1,minExp = 0,maxExp = 6)
dev.off()


#####滑膜成纤维HA和OA
library(Seurat)
library(patchwork)
library(dplyr)
library(ggplot2)
setwd("/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Sy/second")
sce1 = readRDS("/mnt/8w/data7/yiyuan/Bone/test1/Fibroblasts/new/sce_Fb_resubtype.final.rds")
head(sce1@meta.data)

HA <- sce1[, sce1$gname %in% c("HA")]
HA
OA <- sce1[, sce1$gname %in% c("OA")]
OA
umap = sce1@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = sce1@meta.data$celltype) # 注释后的label信息 ，改为cell_type
###"#F08080","#1E90FF","#7CFC00","#FFFF00", "#808000","#FF00FF","#FA8072","#7B68EE"
allcolour=c('#DBC9B3','#EED0E0','#EBAEA9','#CB95BB','#EBCC96','#AED0DF','#CBE5DE')#注意：颜色数量一顶要大于细胞类型数量
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
                   xend = min(umap$UMAP_1) +2, yend = min(umap$UMAP_2) ),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm")))+ 
  geom_segment(aes(x = min(umap$UMAP_1)  , y = min(umap$UMAP_2)  ,
                   xend = min(umap$UMAP_1) , yend = min(umap$UMAP_2) + 2),
               colour = "black", size=1,arrow = arrow(length = unit(0.3,"cm"))) +
  annotate("text", x = min(umap$UMAP_1) +1, y = min(umap$UMAP_2) -0.5, label = "UMAP_1",
           color="black",size = 3, fontface="bold" ) + 
  annotate("text", x = min(umap$UMAP_1) -0.5, y = min(umap$UMAP_2) + 1, label = "UMAP_2",
           color="black",size = 3, fontface="bold" ,angle=90) 
png(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Sy/second/Sy_second_celltype_umap.png", width = 3000, height = 2500, res = 300)
p41
dev.off()
umap = HA@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = HA@meta.data$celltype) # 注释后的label信息 ，改为cell_type
###"#F08080","#1E90FF","#7CFC00","#FFFF00", "#808000","#FF00FF","#FA8072","#7B68EE"
allcolour=c('#DBC9B3','#EED0E0','#EBAEA9','#CB95BB','#EBCC96','#AED0DF','#CBE5DE') #注意：颜色数量一顶要大于细胞类型数量
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
png(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Sy/second/Sy_second_celltype_umap_HA.png", width = 3000, height = 2500, res = 300)
p41
dev.off()

umap = OA@reductions$umap@cell.embeddings %>%  #坐标信息
  as.data.frame() %>% 
  cbind(cell_type = OA@meta.data$celltype) # 注释后的label信息 ，改为cell_type
###"#F08080","#1E90FF","#7CFC00","#FFFF00", "#808000","#FF00FF","#FA8072","#7B68EE"
allcolour=c('#DBC9B3','#EED0E0','#EBAEA9','#CB95BB','#EBCC96','#AED0DF','#CBE5DE') #注意：颜色数量一顶要大于细胞类型数量
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
png(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Sy/second/Sy_second_celltype_umap_OA.png", width = 3000, height = 2500, res = 300)
p41
dev.off()

patient_col=c('#DBC9B3','#EED0E0','#EBAEA9','#CB95BB','#EBCC96','#AED0DF','#CBE5DE')
#分组
group_table <- as.data.frame(table(HA@meta.data$rawname,HA@meta.data$celltype))
names(group_table) <- c("group","celltype","CellNumber")
table(HA$rawname,HA$celltype)
#             Fb1  Fb2  Fb3  Fb4  Fb5  Fb6  Fb7
#  HAS1_1130 2257 1772  817 1833  618  362  236
#  HAS2_1226 1353  498  460 1776 1129  748  294
#  HAS3_0228 1398  526  216 1700  384  524  179
# 计算每个组内的总和
group_table <- group_table %>%
  group_by(group) %>%
  mutate(total_cells = sum(CellNumber)) %>%
  ungroup()

# 计算每个细胞类型的占比
group_table <- group_table %>%
  mutate(percentage = CellNumber / total_cells)

pdf(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Sy/second/proportion_merged_Sample_celltype_HA.pdf", width = 10, height = 12)
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
#             Fb1  Fb2  Fb3  Fb4  Fb5  Fb6  Fb7
#  OAS1_1209 2565 1481 1491 3494 1473  505  205
#  OAS2_1209 2116 1401 1007 2967  794  453  113
#  OAS3_0108  818 3420 3550 2836  518   27   32

# 计算每个组内的总和
group_table <- group_table %>%
  group_by(group) %>%
  mutate(total_cells = sum(CellNumber)) %>%
  ungroup()

# 计算每个细胞类型的占比
group_table <- group_table %>%
  mutate(percentage = CellNumber / total_cells)

pdf(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Sy/second/proportion_merged_Sample_celltype_OA.pdf", width = 10, height = 12)
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
#             Fb1  Fb2  Fb3  Fb4  Fb5  Fb6  Fb7
#  HAS1_1130 2257 1772  817 1833  618  362  236
#  HAS2_1226 1353  498  460 1776 1129  748  294
#  HAS3_0228 1398  526  216 1700  384  524  179
#  OAS1_1209 2565 1481 1491 3494 1473  505  205
#  OAS2_1209 2116 1401 1007 2967  794  453  113
#  OAS3_0108  818 3420 3550 2836  518   27   32
group_table <- group_table %>%
  group_by(group) %>%
  mutate(total_cells = sum(CellNumber)) %>%
  ungroup()

# 计算每个细胞类型的占比
group_table <- group_table %>%
  mutate(percentage = CellNumber / total_cells)

pdf(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Sy/second/proportion_merged_Sample_celltype.pdf", width = 10, height = 12)
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
#      Fb1  Fb2  Fb3  Fb4  Fb5  Fb6  Fb7
#  HA 5008 2796 1493 5309 2131 1634  709
#  OA 5499 6302 6048 9297 2785  985  350
group_table <- group_table %>%
  group_by(group) %>%
  mutate(total_cells = sum(CellNumber)) %>%
  ungroup()

# 计算每个细胞类型的占比
group_table <- group_table %>%
  mutate(percentage = CellNumber / total_cells)

pdf(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Sy/second/proportion_merged_group_celltype.pdf", width = 10, height = 12)
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


sce = readRDS("/mnt/8w/data7/yiyuan/Bone/test1/Fibroblasts/new/sce_Fb_resubtype.final.rds")
png(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Sy/new/marker_gene/Feature_Fb滑膜_sublining.png", width = 3700, height = 3500, res = 300)
# no facet group
featureCornerAxes(object = sce,reduction = 'umap',
                  groupFacet = NULL,
                  relLength = 0.5,
                  relDist = 0.2,
                  features = c("C3","CD34","CDH11","THY1"),
                  aspect.ratio = 1,themebg = 'bwCorner',minExp = 0,maxExp = 6)
dev.off()

png(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Sy/new/marker_gene/Feature_Fb滑膜_lining.png", width = 3700, height = 3500, res = 300)
# no facet group
featureCornerAxes(object = sce,reduction = 'umap',
                  groupFacet = NULL,
                  relLength = 0.5,
                  relDist = 0.2,
                  features = c("CLIC5","DEFB1","CD55","HBEGF","PRG4"),
                  aspect.ratio = 1,themebg = 'bwCorner',minExp = 0,maxExp = 6)
dev.off()