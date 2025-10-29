library(dplyr)

# 读取文件
df1 <- read.csv("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme3_group/markerList_HA_anno_promoter.tsv", sep = "\t", stringsAsFactors = FALSE)
df2 <- read.csv("/mnt/8w/data7/yiyuan/Bone/test1/Fibroblasts/new/DEGs/Fb滑膜_all.markers_gname.csv", stringsAsFactors = FALSE)

# 去除前后空格并转换为大写
df1$symbol <- toupper(trimws(df1$SYMBOL))
df2$gene <- toupper(trimws(df2$gene))

# 筛选cluster值为HA的基因
df2_HA <- subset(df2, cluster == "HA")

# 找出共同的基因
common_genes <- intersect(df1$symbol, df2_HA$gene)
 summary(common_genes)
#   Length     Class      Mode 
#      276 character character 
# 检查ARF1是否在共同基因列表中
if ("ARF1" %in% common_genes) {
  cat("ARF1 在重合的基因列表中。\n")
} else {
  cat("ARF1 不在重合的基因列表中。\n")
}

# 假设 df1 中包含以下列：'chrom', 'start', 'end', 'symbol'
# 筛选出 symbol 在 common_genes 中的行
bed_data <- subset(df1, SYMBOL %in% common_genes)
# 保存为 BED 文件
write.table(bed_data, file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme3_group/common_genes.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 选择并重命名列为 BED 格式所需的列名
bed_data1 <- bed_data[, c('seqnames', 'start', 'end')]
colnames(bed_data1) <- c('seqnames', 'start', 'end')

# 保存为 BED 文件
write.table(bed_data1, file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme3_group/common_genes.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)



library(dplyr)

# 读取文件
df1 <- read.csv("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme5_group/markerList_HA_0.7_anno_promoter.tsv", sep = "\t", stringsAsFactors = FALSE)
df2 <- read.csv("/mnt/8w/data7/yiyuan/Bone/test/Ch/Re_subtype/DEGs/Ch_all.markers_group.csv", stringsAsFactors = FALSE)

# 去除前后空格并转换为大写
df1$symbol <- toupper(trimws(df1$SYMBOL))
df2$gene <- toupper(trimws(df2$gene))

# 筛选cluster值为HA的基因
df2_HA <- subset(df2, cluster == "HA")

# 找出共同的基因
common_genes <- intersect(df1$symbol, df2_HA$gene)
 summary(common_genes)
#   Length     Class      Mode 
#      276 character character 
# 检查ARF1是否在共同基因列表中
if ("ARF1" %in% common_genes) {
  cat("ARF1 在重合的基因列表中。\n")
} else {
  cat("ARF1 不在重合的基因列表中。\n")
}

# 假设 df1 中包含以下列：'chrom', 'start', 'end', 'symbol'
# 筛选出 symbol 在 common_genes 中的行
bed_data <- subset(df1, SYMBOL %in% common_genes)
# 保存为 BED 文件
write.table(bed_data, file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme5_group/common_genes_0.7.tsv", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)

# 选择并重命名列为 BED 格式所需的列名
bed_data1 <- bed_data[, c('seqnames', 'start', 'end')]
colnames(bed_data1) <- c('seqnames', 'start', 'end')

# 保存为 BED 文件
write.table(bed_data1, file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme5_group/common_genes_0.7.bed", sep = "\t", quote = FALSE, row.names = FALSE, col.names = FALSE)
