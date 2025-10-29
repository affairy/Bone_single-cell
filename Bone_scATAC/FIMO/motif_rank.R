# 加载必要的库
library(ggplot2)
library(ggrepel)

# 读取文件（去除多余空格）
file_path <- "/mnt/8w/data4/yiyuan/Bone/scATAC/result/fimo/up_fimo_Ch/motif_count.txt"
motif_data <- read.table(file_path, header = FALSE, sep = "", stringsAsFactors = FALSE, fill = TRUE)

# 删除空行或不完整行
motif_data <- motif_data[complete.cases(motif_data), ]

# 重命名列名
colnames(motif_data) <- c("Frequency", "Motif")

# 按照 Frequency 降序排序
motif_data <- motif_data[order(motif_data$Frequency, decreasing = TRUE), ]

# 添加 rank 列
motif_data$Rank <- seq_len(nrow(motif_data))

# 输出 PDF 文件
pdf("/mnt/8w/data4/yiyuan/Bone/scATAC/result/fimo/up_fimo_Ch/Motif_Rank_Plot.pdf", width = 8, height = 6)

# 绘制 rank plot
ggplot(motif_data, aes(x = Rank, y = Frequency)) +
  geom_point(size = 2, color = "blue") +
  geom_line(group = 1, color = "blue", linetype = "dashed") +
  # 使用 ggrepel 自动调整标签位置，防止重叠
  geom_text_repel(data = motif_data[1:20, ], 
                  aes(label = Motif), 
                  size = 3, 
                  color = "red",
                  box.padding = 0.5, # 增加标签周围的间距
                  max.overlaps = 20) + # 防止标签重叠
  theme_minimal() +
  labs(x = "Rank", y = "Frequency", title = "Motif Rank Plot") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# 关闭输出设备
dev.off()


###滑膜
# 加载必要的库
library(ggplot2)
library(ggrepel)

# 读取文件（去除多余空格）
file_path <- "/mnt/8w/data4/yiyuan/Bone/scATAC/result/fimo/up_fimo/motif_count.txt"
motif_data <- read.table(file_path, header = FALSE, sep = "", stringsAsFactors = FALSE, fill = TRUE)

# 删除空行或不完整行
motif_data <- motif_data[complete.cases(motif_data), ]

# 重命名列名
colnames(motif_data) <- c("Frequency", "Motif")

# 按照 Frequency 降序排序
motif_data <- motif_data[order(motif_data$Frequency, decreasing = TRUE), ]

# 添加 rank 列
motif_data$Rank <- seq_len(nrow(motif_data))

# 输出 PDF 文件
pdf("/mnt/8w/data4/yiyuan/Bone/scATAC/result/fimo/up_fimo/Motif_Rank_Plot.pdf", width = 8, height = 6)

# 绘制 rank plot
ggplot(motif_data, aes(x = Rank, y = Frequency)) +
  geom_point(size = 2, color = "blue") +
  geom_line(group = 1, color = "blue", linetype = "dashed") +
  # 使用 ggrepel 自动调整标签位置，防止重叠
  geom_text_repel(data = motif_data[1:20, ], 
                  aes(label = Motif), 
                  size = 3, 
                  color = "red",
                  box.padding = 0.5, # 增加标签周围的间距
                  max.overlaps = 20) + # 防止标签重叠
  theme_minimal() +
  labs(x = "Rank", y = "Frequency", title = "Motif Rank Plot") +
  theme(axis.text.x = element_text(angle = 90, hjust = 1))

# 关闭输出设备
dev.off()
