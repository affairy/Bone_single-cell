setwd("/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/monocle")
library(Seurat)
library(Matrix)
library(harmony)
library(ggplot2)
library(tidyverse)
sce = readRDS("/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/sce_Ch_resubtype.final.rds")
###查看基本软骨信息
head(sce@meta.data)
table(sce$rawname,sce$celltype)

###分别提取HA和OA组
HA <- sce[, sce$gname %in% c("HA")]
HA
OA <- sce[, sce$gname %in% c("OA")]
OA

saveRDS(HA,file="sce_Ch_resubtype.final_HA.rds")
saveRDS(OA,file="sce_Ch_resubtype.final_OA.rds")
library(Seurat)
library(monocle)
library(dplyr)

HA = readRDS("/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/monocle/sce_Ch_resubtype.final_HA.rds")
table(Idents(HA))

table(HA$sample)
head(HA@meta.data)


table(HA$celltype)
set.seed(1234)

# count矩阵，官方建议用count
ct <- HA@assays$origin.counts$counts
# 基因注释
gene_ann <- data.frame(
  gene_short_name = row.names(ct), 
  row.names = row.names(ct)
)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
# 临床信息
pd <- new("AnnotatedDataFrame",
          data=HA@meta.data)

#新建CellDataSet对象
sc_cds <- newCellDataSet(
  ct, 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds

# 用于估算每个细胞的大小因子,目的是归一化基因表达数据以纠正测序深度的差异
sc_cds <- estimateSizeFactors(sc_cds)
# 估算基因表达的离散度，帮助识别高变异的基因，用于差异基因分析和排序
sc_cds <- estimateDispersions(sc_cds)

diff.wilcox = FindAllMarkers(HA)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "top10_diff_genes_wilcox.csv", row.names = F)
mycds=sc_cds
save(mycds, file = "monocle_results.RData")

setwd("/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/monocle/HA")
load("/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/monocle/monocle_results.RData")
##使用clusters差异表达基因
diff.genes <- read.csv('/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/monocle/HA/diff_genes_wilcox.csv')
diff.genes <- subset(diff.genes,p_val_adj<0.01)$gene
mycds <- setOrderingFilter(mycds, diff.genes)
p1 <- plot_ordering_genes(mycds)
p1
pdf(file="Selectgene.pdf",width=10,height=8)
plot_ordering_genes(mycds)
dev.off()

##########################使用disp.genes开展后续分析
#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)
#State轨迹分布图

plot1 <- plot_cell_trajectory(mycds, color_by = "State")
ggsave("State.pdf", plot = plot1, width = 6, height = 5)
#ggsave("pseudotime/State.png", plot = plot1, width = 6, height = 5)
##Cluster轨迹分布图
plot2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
ggsave("Cluster.pdf", plot = plot2, width = 6, height = 5)
#ggsave("pseudotime/Cluster.png", plot = plot2, width = 6, height = 5)
##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("Pseudotime.pdf", plot = plot3, width = 6, height = 5)
#ggsave("pseudotime/Pseudotime.png", plot = plot3, width = 6, height = 5)
##合并作图
plot4 <- plot_cell_trajectory(mycds, color_by = "celltype")
ggsave("celltype.pdf", plot = plot4, width = 6, height = 5)
plotc <- plot1|plot2|plot3|plot4
ggsave("Combination.pdf", plot = plotc, width = 15, height = 5)
#ggsave("pseudotime/Combination.png", plot = plotc, width = 10, height = 3.5)
##保存结果
write.csv(pData(mycds), "pseudotime.csv")

p1 <- plot_cell_trajectory(mycds, color_by = "State") + facet_wrap(~State, nrow = 2)
p2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters") + facet_wrap(~seurat_clusters, nrow = 2)
plotc <- p1/p2
ggsave("trajectory_facet2.png", plot = plotc, width = 20, height = 17)

p3 <- plot_cell_trajectory(mycds, color_by = "celltype") +
  facet_wrap(~celltype, nrow = 2) +
  theme(
    axis.text = element_text(size = 14),           # 坐标轴刻度字体
    axis.title = element_text(size = 16),          # 坐标轴标题字体
    strip.text = element_text(size = 14),          # facet 标签字体
    plot.title = element_text(size = 18, face = "bold"),  # 图表标题字体
    legend.text = element_text(size = 14),         # 图例文本字体
    legend.title = element_text(size = 16)         # 图例标题字体
  )
# 保存图片
ggsave("trajectory_facet_celltype.png", plot = p3, width = 25, height = 17)


######拟时序差异基因热图
#cluster差异基因
diff.genes <- read.csv('/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/monocle/HA/diff_genes_wilcox.csv')
sig_diff.genes <- subset(diff.genes,p_val_adj<0.0001&abs(avg_log2FC)>0.75)$gene
sig_diff.genes <- unique(as.character(sig_diff.genes))
diff_test <- differentialGeneTest(mycds[sig_diff.genes,], cores = 1, 
                              fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(diff_test,"Time_diff_test.csv",row.names = F)
sig_gene_names <- row.names(subset(diff_test, qval < 0.01))
p1 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=8,
                             show_rownames=T, return_heatmap=T)
ggsave("pseudotime_heatmap1.png", plot = p1, width = 5, height = 8)
#p <- plot_pseudotime_heatmap(mouse_monocle[Time_genes,], 
#                             num_cluster = 4, 
#                             show_rownames = T, 
#                             return_heatmap = T,
#                             hmcols = colorRampPalette(c("navy","white","firebrick3")))

save(mycds, file = "monocle_results_HA_final.RData")

clusters <- cutree(p1$tree_row,k=8)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
write.csv(clustering,"Time_clustering_all.csv", row.names = F)

Time_genes <- diff_test %>%
  top_n(n = 100, wt = qval) %>%  # 按 qval 排序并选择前 100 个基因
  arrange(qval) %>%  # 按 qval 排序
  pull(gene_short_name) %>%  # 获取基因名称
  as.character()  # 转换为字符向量
write.csv(Time_genes,"Time_diff_test100.csv",row.names = F)
p <- plot_pseudotime_heatmap(mycds[Time_genes, ], num_clusters = 8, show_rownames = TRUE, return_heatmap = TRUE)

# 保存热图到 PDF 文件
ggsave("Time_heatmapTop100.png", p, width = 5, height = 10)

###### AMTN\CHAD\COL1A1\CRTAC1\ANGPTL4\CLEC3A\COL1A2\CYTL1

#s.genes <- c("ANXA1")
#p3 <- plot_genes_jitter(mycds[s.genes,], grouping = "State", color_by = "State")
#p4 <- plot_genes_violin(mycds[s.genes,], grouping = "State", color_by = "State")
#p5 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "State")
#plotc <- p3|p4|p5
#ggsave("genes_visual.png", plot = plotc, width = 15, height = 5)
#
#
#s.genes <- c("ANXA1")
#p3 <- plot_genes_jitter(mycds[s.genes,], grouping = "celltype", color_by = "celltype")
#p4 <- plot_genes_violin(mycds[s.genes,], grouping = "celltype", color_by = "celltype")
#p5 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "celltype")
#plotc <- p3|p4|p5
#ggsave("genes_visual_celltype_ANXA1.png", plot = plotc, width = 30, height = 8)
#
#s.genes <- c("FPR1")
#p3 <- plot_genes_jitter(mycds[s.genes,], grouping = "celltype", color_by = "celltype")
#p4 <- plot_genes_violin(mycds[s.genes,], grouping = "celltype", color_by = "celltype")
#p5 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "celltype")
#plotc <- p3|p4|p5
#ggsave("genes_visual_celltype_FPR1.png", plot = plotc, width = 30, height = 8)


keygenes <- head(diff.genes,5)
keygenes_vector <- keygenes$gene
cds_subset <- mycds[keygenes_vector,]

p1 <-plot_genes_in_pseudotime(cds_subset, color_by = "State")

p2 <-plot_genes_in_pseudotime(cds_subset, color_by = "celltype")

p3 <-plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")

plotc <-p1|p2|p3
ggsave("genes_visual_keygenes.png", plot = plotc, width = 15, height = 8)

#colnames(pData(mycds))
#
#pData(mycds)$ANXA1 = log2(exprs(mycds)['ANXA1',]+1)
#
#p1 = plot_cell_trajectory(mycds, color_by = "ANXA1") + 
#  scale_color_viridis_c()  # 使用 Viridis 配色方案
#ggsave("trajectory_facet_anax1.png", plot = p1, width = 15, height = 8)
#
#pData(mycds)$FPR1 = log2(exprs(mycds)['FPR1',]+1)
#
#p2 = plot_cell_trajectory(mycds, color_by = "FPR1") + 
#  scale_color_viridis_c()  # 使用 Viridis 配色方案
#ggsave("trajectory_facet_FPR1.png", plot = p2, width = 15, height = 8)
#
save(mycds, file = "monocle_results_HA_Ch.RData")

######OA组数据分析
library(Seurat)
library(Matrix)
library(harmony)
library(ggplot2)
library(tidyverse)
library(Seurat)
library(monocle)
library(dplyr)
setwd("/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/monocle/OA")
OA = readRDS("/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/monocle/sce_Ch_resubtype.final_OA.rds")
table(Idents(OA))

table(OA$sample)
head(OA@meta.data)


table(OA$celltype)
set.seed(1234)

# count矩阵，官方建议用count
ct <- OA@assays$origin.counts$counts
# 基因注释
gene_ann <- data.frame(
  gene_short_name = row.names(ct), 
  row.names = row.names(ct)
)
fd <- new("AnnotatedDataFrame",
          data=gene_ann)
# 临床信息
pd <- new("AnnotatedDataFrame",
          data=OA@meta.data)

#新建CellDataSet对象
sc_cds <- newCellDataSet(
  ct, 
  phenoData = pd,
  featureData =fd,
  expressionFamily = negbinomial.size(),
  lowerDetectionLimit=1)
sc_cds

# 用于估算每个细胞的大小因子,目的是归一化基因表达数据以纠正测序深度的差异
sc_cds <- estimateSizeFactors(sc_cds)
# 估算基因表达的离散度，帮助识别高变异的基因，用于差异基因分析和排序
sc_cds <- estimateDispersions(sc_cds)

diff.wilcox = FindAllMarkers(OA)
all.markers = diff.wilcox %>% select(gene, everything()) %>% subset(p_val<0.05)
top10 = all.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
write.csv(all.markers, "diff_genes_wilcox.csv", row.names = F)
write.csv(top10, "top10_diff_genes_wilcox.csv", row.names = F)
mycds=sc_cds
save(mycds, file = "monocle_results_OA.RData")

setwd("/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/monocle/OA")
load("/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/monocle/OA/monocle_results_OA.RData")
##使用clusters差异表达基因
diff.genes <- read.csv('/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/monocle/HA/diff_genes_wilcox.csv')
diff.genes <- subset(diff.genes,p_val_adj<0.01)$gene
mycds <- setOrderingFilter(mycds, diff.genes)
p1 <- plot_ordering_genes(mycds)
p1
pdf(file="Selectgene.pdf",width=10,height=8)
plot_ordering_genes(mycds)
dev.off()

##########################使用disp.genes开展后续分析
#降维
mycds <- reduceDimension(mycds, max_components = 2, method = 'DDRTree')
#排序
mycds <- orderCells(mycds)
#State轨迹分布图

plot1 <- plot_cell_trajectory(mycds, color_by = "State")
ggsave("State.pdf", plot = plot1, width = 6, height = 5)
#ggsave("pseudotime/State.png", plot = plot1, width = 6, height = 5)
##Cluster轨迹分布图
plot2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters")
ggsave("Cluster.pdf", plot = plot2, width = 6, height = 5)
#ggsave("pseudotime/Cluster.png", plot = plot2, width = 6, height = 5)
##Pseudotime轨迹图
plot3 <- plot_cell_trajectory(mycds, color_by = "Pseudotime")
ggsave("Pseudotime.pdf", plot = plot3, width = 6, height = 5)
#ggsave("pseudotime/Pseudotime.png", plot = plot3, width = 6, height = 5)
##合并作图
plot4 <- plot_cell_trajectory(mycds, color_by = "celltype")
ggsave("celltype.pdf", plot = plot4, width = 6, height = 5)
plotc <- plot1|plot2|plot3|plot4
ggsave("Combination.pdf", plot = plotc, width = 15, height = 5)
#ggsave("pseudotime/Combination.png", plot = plotc, width = 10, height = 3.5)
##保存结果
write.csv(pData(mycds), "pseudotime.csv")

p1 <- plot_cell_trajectory(mycds, color_by = "State") + facet_wrap(~State, nrow = 2)
p2 <- plot_cell_trajectory(mycds, color_by = "seurat_clusters") + facet_wrap(~seurat_clusters, nrow = 2)
plotc <- p1/p2
ggsave("trajectory_facet2.png", plot = plotc, width = 20, height = 17)

p3 <- plot_cell_trajectory(mycds, color_by = "celltype") +
  facet_wrap(~celltype, nrow = 2) +
  theme(
    axis.text = element_text(size = 14),           # 坐标轴刻度字体
    axis.title = element_text(size = 16),          # 坐标轴标题字体
    strip.text = element_text(size = 14),          # facet 标签字体
    plot.title = element_text(size = 18, face = "bold"),  # 图表标题字体
    legend.text = element_text(size = 14),         # 图例文本字体
    legend.title = element_text(size = 16)         # 图例标题字体
  )
# 保存图片
ggsave("trajectory_facet_celltype.png", plot = p3, width = 25, height = 17)


######拟时序差异基因热图
#cluster差异基因
diff.genes <- read.csv('/mnt/data7/yiyuan/Bone/test/Ch/Re_subtype/monocle/HA/diff_genes_wilcox.csv')
sig_diff.genes <- subset(diff.genes,p_val_adj<0.0001&abs(avg_log2FC)>0.75)$gene
sig_diff.genes <- unique(as.character(sig_diff.genes))
diff_test <- differentialGeneTest(mycds[sig_diff.genes,], cores = 1, 
                              fullModelFormulaStr = "~sm.ns(Pseudotime)")
write.csv(diff_test,"Time_diff_test.csv",row.names = F)
sig_gene_names <- row.names(subset(diff_test, qval < 0.01))
p1 = plot_pseudotime_heatmap(mycds[sig_gene_names,], num_clusters=8,
                             show_rownames=T, return_heatmap=T)
ggsave("pseudotime_heatmap1.png", plot = p1, width = 5, height = 8)
#p <- plot_pseudotime_heatmap(mouse_monocle[Time_genes,], 
#                             num_cluster = 4, 
#                             show_rownames = T, 
#                             return_heatmap = T,
#                             hmcols = colorRampPalette(c("navy","white","firebrick3")))

save(mycds, file = "monocle_results_HA_final.RData")

clusters <- cutree(p1$tree_row,k=8)
clustering <- data.frame(clusters)
clustering[,1] <- as.character(clustering[,1])
colnames(clustering) <- "Gene_Clusters"
table(clustering)
write.csv(clustering,"Time_clustering_all.csv", row.names = F)

Time_genes <- diff_test %>%
  top_n(n = 100, wt = qval) %>%  # 按 qval 排序并选择前 100 个基因
  arrange(qval) %>%  # 按 qval 排序
  pull(gene_short_name) %>%  # 获取基因名称
  as.character()  # 转换为字符向量
write.csv(Time_genes,"Time_diff_test100.csv",row.names = F)
p <- plot_pseudotime_heatmap(mycds[Time_genes, ], num_clusters = 8, show_rownames = TRUE, return_heatmap = TRUE)

# 保存热图到 PDF 文件
ggsave("Time_heatmapTop100.png", p, width = 5, height = 10)

###### AMTN\CHAD\COL1A1\CRTAC1\ANGPTL4\CLEC3A\COL1A2\CYTL1

#s.genes <- c("ANXA1")
#p3 <- plot_genes_jitter(mycds[s.genes,], grouping = "State", color_by = "State")
#p4 <- plot_genes_violin(mycds[s.genes,], grouping = "State", color_by = "State")
#p5 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "State")
#plotc <- p3|p4|p5
#ggsave("genes_visual.png", plot = plotc, width = 15, height = 5)
#
#
#s.genes <- c("ANXA1")
#p3 <- plot_genes_jitter(mycds[s.genes,], grouping = "celltype", color_by = "celltype")
#p4 <- plot_genes_violin(mycds[s.genes,], grouping = "celltype", color_by = "celltype")
#p5 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "celltype")
#plotc <- p3|p4|p5
#ggsave("genes_visual_celltype_ANXA1.png", plot = plotc, width = 30, height = 8)
#
#s.genes <- c("FPR1")
#p3 <- plot_genes_jitter(mycds[s.genes,], grouping = "celltype", color_by = "celltype")
#p4 <- plot_genes_violin(mycds[s.genes,], grouping = "celltype", color_by = "celltype")
#p5 <- plot_genes_in_pseudotime(mycds[s.genes,], color_by = "celltype")
#plotc <- p3|p4|p5
#ggsave("genes_visual_celltype_FPR1.png", plot = plotc, width = 30, height = 8)


keygenes <- head(diff.genes,5)
keygenes_vector <- keygenes$gene
cds_subset <- mycds[keygenes_vector,]

p1 <-plot_genes_in_pseudotime(cds_subset, color_by = "State")

p2 <-plot_genes_in_pseudotime(cds_subset, color_by = "celltype")

p3 <-plot_genes_in_pseudotime(cds_subset, color_by = "Pseudotime")

plotc <-p1|p2|p3
ggsave("genes_visual_keygenes.png", plot = plotc, width = 15, height = 8)

#colnames(pData(mycds))
#
#pData(mycds)$ANXA1 = log2(exprs(mycds)['ANXA1',]+1)
#
#p1 = plot_cell_trajectory(mycds, color_by = "ANXA1") + 
#  scale_color_viridis_c()  # 使用 Viridis 配色方案
#ggsave("trajectory_facet_anax1.png", plot = p1, width = 15, height = 8)
#
#pData(mycds)$FPR1 = log2(exprs(mycds)['FPR1',]+1)
#
#p2 = plot_cell_trajectory(mycds, color_by = "FPR1") + 
#  scale_color_viridis_c()  # 使用 Viridis 配色方案
#ggsave("trajectory_facet_FPR1.png", plot = p2, width = 15, height = 8)
#
save(mycds, file = "monocle_results_HA_Ch.RData")

#nohup Rscript /mnt/data1/yiyuan/script/2025.1/Bone_remake/软骨总_monocle.R > ch_mergerd.log 2>&1 &
#nohup Rscript /mnt/data1/yiyuan/script/2025.1/Bone_remake/软骨总monocle_new.R > 软骨总monocle_new.log 2>&1 &
#nohup Rscript /mnt/data1/yiyuan/script/2025.1/Bone_remake/ha_sy_fb.R > ha_sy_fb.log 2>&1 &
#nohup Rscript /mnt/data1/yiyuan/script/2025.1/Bone_remake/滑膜成纤维总_monocle.R > 滑膜成纤维总.log 2>&1 &
#nohup Rscript /mnt/data1/yiyuan/script/2025.1/Bone_remake/滑膜巨噬细胞总_monocle.R > 滑膜巨噬细胞总.log 2>&1 &
