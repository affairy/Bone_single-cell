#加载包
library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(ggsci)
library(igraph)
library(ggraph)
setwd("/mnt/8w/data4/yiyuan/Bone/scATAC/result/fimo")
##读取数据
data <- read.csv('/mnt/8w/data4/yiyuan/Bone/scATAC/result/fimo/HA_up_motif_network_Ch_new.csv', header = T)
head(data)
#构建边和节点对象
highlight_edges <- data.frame(
  from = c("FOS", "FOS", "EGR1", "EGR1", "EGR1", "EGR1", "EGR1",
           "KLF10", "KLF10", "KLF10", "KLF10", "KLF10","JUNB"),
  to = c("FOXO3", "FTL", "FOXO3", "FTH1", "HMOX1", "NAMPT", "PRRX1",
         "FOXO3", "FTL", "HMOX1", "NAMPT", "PRRX1", "NAMPT")
)

# ==== Step 2: 构建节点对象 ====
edges <- data.frame(from = data$TFs, to = data$Genes)
vertices <- data.frame(name = unique(c(data$TFs, data$Genes)), type = 'Targets', show_name = NA, size = 1)
vertices$type[match(unique(data$TFs), vertices$name)] <- 'TFs'
vertices$show_name <- vertices$name
vertices$size[match(unique(data$TFs), vertices$name)] <- 2

# ==== Step 3: 高亮蓝色节点 ====
# 获取需要设为蓝色的节点集合
blue_nodes <- unique(c(highlight_edges$from, highlight_edges$to))
vertices$color <- ifelse(vertices$name %in% blue_nodes, "blue", "red")

# ==== Step 4: 绘图 ====
ggraph_data <- igraph::graph_from_data_frame(d = edges, vertices = vertices, directed = TRUE)

p <- ggraph(graph = ggraph_data, layout = 'kk', circular = FALSE) +
  geom_edge_link(color = 'grey', edge_alpha = 0.55,
                 arrow = arrow(type = "closed", angle = 15, ends = "last", length = unit(0.1, "inches"))) +
  geom_node_point(aes(shape = type, fill = color, size = size, color = color), alpha = 0.75) +
  geom_node_text(aes(label = show_name,
                     color = color),  # 特殊标红文本可根据需要调整
                 size = 4,
                 repel = TRUE, nudge_x = 0.05, nudge_y = 0.2) +
  theme_void() +
  scale_color_identity() +  # 使用实际颜色名
  scale_fill_identity() +
  scale_size_continuous(range = c(2, 4)) +
  scale_shape_manual(values = c(22, 21)) +
  coord_equal() +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = "white"))

# 保存
ggsave("graph_plot_Ch_highlight_blue.png", plot = p, width = 18, height = 12, dpi = 300, bg = "white", limitsize = FALSE)




library(ggplot2)
library(RColorBrewer)
library(tidyverse)
library(ggsci)
library(igraph)
library(ggraph)
setwd("/mnt/8w/data4/yiyuan/Bone/scATAC/result/fimo")
##读取数据
data <- read.csv('/mnt/8w/data4/yiyuan/Bone/scATAC/result/fimo/HA_up_motif_network_SY_new.csv', header = T)
head(data)
#构建边和节点对象
highlight_edges <- data.frame(
  from = c("STAT1"),
  to = c("PRRX1")
)

# ==== Step 2: 构建节点对象 ====
edges <- data.frame(from = data$TFs, to = data$Genes)
vertices <- data.frame(name = unique(c(data$TFs, data$Genes)), type = 'Targets', show_name = NA, size = 1)
vertices$type[match(unique(data$TFs), vertices$name)] <- 'TFs'
vertices$show_name <- vertices$name
vertices$size[match(unique(data$TFs), vertices$name)] <- 2

# ==== Step 3: 高亮蓝色节点 ====
# 获取需要设为蓝色的节点集合
blue_nodes <- unique(c(highlight_edges$from, highlight_edges$to))
vertices$color <- ifelse(vertices$name %in% blue_nodes, "blue", "red")

# ==== Step 4: 绘图 ====
ggraph_data <- igraph::graph_from_data_frame(d = edges, vertices = vertices, directed = TRUE)

p <- ggraph(graph = ggraph_data, layout = 'kk', circular = FALSE) +
  geom_edge_link(color = 'grey', edge_alpha = 0.55,
                 arrow = arrow(type = "closed", angle = 15, ends = "last", length = unit(0.1, "inches"))) +
  geom_node_point(aes(shape = type, fill = color, size = size, color = color), alpha = 0.75) +
  geom_node_text(aes(label = show_name,
                     color = color),  # 特殊标红文本可根据需要调整
                 size = 4,
                 repel = TRUE) +
  theme_void() +
  scale_color_identity() +  # 使用实际颜色名
  scale_fill_identity() +
  scale_size_continuous(range = c(2, 4)) +
  scale_shape_manual(values = c(22, 21)) +
  coord_equal() +
  theme(legend.position = 'none',
        panel.background = element_rect(fill = "white"))

# 保存
ggsave("graph_plot_SY_highlight_blue.png", plot = p, width = 18, height = 12, dpi = 300, bg = "white", limitsize = FALSE)