library(ArchR)
#ArchR::installExtraPackages()
library(pheatmap) 
library(Rsamtools)
library(scran) 
#library(scater) 
library(dplyr)
library(Seurat)
library(patchwork)
library(SingleCellExperiment) 
library(ComplexHeatmap)
library(ggplot2)
library(stringr)
#BiocManager::install('EnsDb.Hsapiens.v86')
library(EnsDb.Hsapiens.v86) 
library(viridis) 

#addArchRGenome("hg19")#设置基因组（人hg19），根据自己实际情况
addArchRGenome("hg38")#设置基因组（人hg38）
# addArchRGenome("mm10")#小鼠
addArchRThreads(threads = 10)#设置线程、根据自己电脑实际情况
###八例scATAC数据 
###软骨分析注释
setwd("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch")
#input文件路径,只需要样本的atac_fragments.tsv.gz文件
input.file.list = c('/mnt/8w/data4/yiyuan/Bone/scATAC/data/20250207/8例完整矩阵/HAC0416/outs/fragments_corrected_dedup_count.tsv.gz',
                    '/mnt/8w/data4/yiyuan/Bone/scATAC/data/20250207/8例完整矩阵/HAEC1/outs/fragments_corrected_dedup_count.tsv.gz',
                    '/mnt/8w/data4/yiyuan/Bone/scATAC/data/20250207/8例完整矩阵/OAC_ZXG/outs/fragments_corrected_dedup_count.tsv.gz',
                    '/mnt/8w/data4/yiyuan/Bone/scATAC/data/20250207/8例完整矩阵/OAC_0510/outs/fragments_corrected_dedup_count.tsv.gz')
#设置样本名，相当于seurat中orig.ident。
sampleNames = c("HAC0416","HAEC1","OAC_ZXG","OAC_0510")

#创建Arrow文件
ArrowFiles <- createArrowFiles(inputFiles = input.file.list,
                               sampleNames = sampleNames,
                               minTSS = 4, #默认值这里有过滤，自行设置，不宜太高
                               minFrags = 1000,
                               addTileMat = TRUE,
                               addGeneScoreMat = TRUE,
                               excludeChr = c("chrM", "chrY", "chrX"))
###去除双细胞，计算双细胞评分
doubScores <- addDoubletScores(input = ArrowFiles,
                               k = 20, #Refers to how many cells near a "pseudo-doublet" to count.
                               knnMethod = "UMAP",
                               useMatrix = "TileMatrix",
                               nTrials=5,
                               LSIMethod = 1,
                               scaleDims = F,
                               corCutOff = 0.75,
                               UMAPParams = list(n_neighbors =30, 
                                                 min_dist = 0.3, 
                                                 metric = "cosine", 
                                                 verbose =T))
#第四步：构建ArchRProject，这个就类似于seurat对象，是后续分析的基础。

#构建ArchRProject
projncov <- ArchRProject(ArrowFiles = ArrowFiles,
                         outputDirectory = "ATAC_out", #结果存储文件夹
                         copyArrows = TRUE) #建议使用T，保存arrows副本

table(projncov@cellColData$Sample)#查看每个样本细胞数
# HAC0416    HAEC1 OAC_0510  OAC_ZXG 
#     447     1780     3232     1093 
#####数据质控
#过滤双细胞------------------------------------------------------------
proj.filter <- filterDoublets(projncov)
#Filtering 146 cells from ArchRProject!
#        HAEC1 : 31 of 1780 (1.7%)
#        OAC_0510 : 104 of 3232 (3.2%)
#        OAC_ZXG : 11 of 1093 (1%)
#        HAC0416 : 0 of 447 (0%)
table(proj.filter@cellColData$Sample)#查看每个样本细胞数
# HAC0416    HAEC1 OAC_0510  OAC_ZXG 
#     447     1749     3128     1082 
#####质控结果绘制
p5 <- plotFragmentSizes(ArchRProj = proj.filter)+
  ggtitle("Fragment Size Histogram")
p6 <- plotTSSEnrichment(ArchRProj = proj.filter)+
  ggtitle("TSS Enrichment")
plotPDF(p5,p6, name = "QC-Sample-FragSizes-TSSProfile_filter.pdf", ArchRProj = proj.filter, addDOC = FALSE, width = 5, height = 5)
###密度图
filterSample_cellNum <- table(proj.filter$Sample)#过滤后样本数
sampleplot_list <- list()
for (i in 1:length(sampleNames)) {
  proj.i <- proj.filter[proj.filter$Sample == sampleNames[i]]#提取每个样本ArchRproject做密度图
  
  p <- ggPoint(
    x = log10(proj.i$nFrags),
    y = proj.i$TSSEnrichment,
    colorDensity = TRUE,
    continuousSet = "sambaNight",
    xlabel = "Log10(Unique Fragments)",
    ylabel = "TSS Enrichment")+ 
    geom_hline(yintercept = 4, lty = "dashed") + 
    geom_vline(xintercept = 3, lty = "dashed")+
    ggtitle(paste0(sampleNames[i],"\n","Cells Pass Filter =",filterSample_cellNum[[i]]))
  
  sampleplot_list[[i]] <- p
}
sampleplot_list
plotPDF(sampleplot_list[[1]], sampleplot_list[[2]], sampleplot_list[[3]], sampleplot_list[[4]],
        name = "QC-Sample-TSSenrich-Frag_filter.pdf", 
        ArchRProj = proj.filter, addDOC = FALSE, width = 5, height = 5)


######################
#降维---dimension reduction 
proj.filter <- addIterativeLSI(ArchRProj = proj.filter,
                               useMatrix = "TileMatrix",
                               name = "IterativeLSI",
                               iterations = 4,#要执行的LSI迭代次数
                               clusterParams = list(resolution = 4,
                                                    sampleCells =10000,
                                                    n.start = 10),
                               varFeatures = 50000,
                               dimsToUse = 1:20,
                               force = TRUE,
                               seed=10)
head(proj.filter@cellColData)
#使用harmony去除批次效应
 proj.filter <- addHarmony(ArchRProj = proj.filter,
                          reducedDims = "IterativeLSI",
                          name = "Harmony",#harmony去除批次
                          groupBy = "Sample",#因为是多样本，所以按照样本去除批次
                          force = T)
###聚类分析
proj.filter <- addClusters(input = proj.filter, 
                           reducedDims = "Harmony", 
                           method = "Seurat",
                           name = "Clusters", 
                           resolution = 1, ###4改为1
                           force=TRUE, 
                           seed = 11)
table(proj.filter@cellColData$Clusters)#查看每个样本细胞数

#  C1  C10  C11  C12   C2   C3   C4   C5   C6   C7   C8   C9 
#  32  879 1483  706  420  338  679  150  189  574  307  649 
###更好的看看cluster在样本中的分布。用热图的形式展示。
cM <- confusionMatrix(paste0(proj.filter$Clusters), paste0(proj.filter$Sample))
cM <- cM / Matrix::rowSums(cM)
#p <- pheatmap::pheatmap(mat = as.matrix(cM), 
#                        color = paletteContinuous("whiteBlue"), border_color = "black")
# 打开PNG设备并指定文件路径和图片大小
png(filename = "heatmap_cluster.png", width = 800, height = 600)
# 生成heatmap
pheatmap::pheatmap(
  mat = as.matrix(cM), 
  color = paletteContinuous("whiteBlue"), 
  border_color = "black"
)

# 关闭图形设备
dev.off()
table(proj.filter@cellColData$Clusters,proj.filter@cellColData$Sample)
#      HAC0416 HAEC1 OAC_0510 OAC_ZXG
#  C1        0     8       16       8
#  C10     111   155      453     160
#  C11      55   330      811     287
#  C12      66   298      256      86
#  C2       12   177      122     109
#  C3       43   129      136      30
#  C4        8    74      481     116
#  C5       15    87       35      13
#  C6        5   116       61       7
#  C7       53   221      229      71
#  C8       68    84      110      45
#  C9       11    70      418     150
proj.filter <- addUMAP(ArchRProj = proj.filter, 
                       reducedDims = "Harmony", 
                       name = "UMAP_harmony",
                       nNeighbors = 30,
                       minDist = 0.5, 
                       metric = "cosine",
                       force=TRUE,
                       seed=12)

proj.filter <- addUMAP(ArchRProj = proj.filter, 
                       reducedDims = "IterativeLSI", 
                       name = "UMAP_LSI",
                       nNeighbors = 30,
                       minDist = 0.5, 
                       metric = "cosine",
                       force=TRUE,
                       seed=12)
proj.filter <- addTSNE(
    ArchRProj = proj.filter, 
    reducedDims = "IterativeLSI", 
    name = "TSNE_LSI",
    perplexity = 30,force = TRUE
)

proj.filter <- addTSNE(
    ArchRProj = proj.filter, 
    reducedDims = "Harmony", 
    name = "TSNE_Harmony",
    perplexity = 30,force = TRUE
)
proj.filter@embeddings
#List of length 4
#names(4): UMAP_harmony UMAP_LSI TSNE_LSI TSNE_Harmony

#saveArchRProject(proj.filter,outputDirectory = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/ATAC_out",overwrite = TRUE)
#proj.filter <- loadArchRProject("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/ATAC_out")

#聚类后一些基本聚类图的可视化---UMAP图
p7 <- plotEmbedding(ArchRProj = proj.filter, 
                    colorBy = "cellColData", 
                    name = "Sample", 
                    embedding = "UMAP_harmony",
                    size = 0.5)

p8 <- plotEmbedding(ArchRProj = proj.filter, 
                    colorBy = "cellColData", 
                    name = "Clusters", 
                    embedding = "UMAP_harmony",
                    size = 0.5)

plotPDF(p7,p8, name = "Plot-UMAP-Sample-Clusters.pdf", 
        ArchRProj = proj.filter, 
        addDOC = FALSE, width = 8, height = 8)

##2、鉴定cluster marker基因
markersGS <- getMarkerFeatures(ArchRProj = proj.filter,
                               useMatrix = "GeneScoreMatrix",
                               groupBy = "Clusters",
                               bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")
for (i in names(markerList)) {
  write.table(markerList[[i]], file = paste0("markerList_", i, ".txt"), sep = "\t", row.names = TRUE, col.names = TRUE)
}
saveArchRProject(proj.filter,outputDirectory = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/ATAC_filtered",overwrite = TRUE)


proj.filter <- loadArchRProject("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/ATAC_filtered")
setwd("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch")
seRNA <- readRDS("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/scRNA_first.rds")
seRNA
colnames(colData(seRNA))
table(colData(seRNA)$celltype_cluster_standard)
proj.filter <- addGeneIntegrationMatrix(
    ArchRProj = proj.filter, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "celltype_cluster_standard",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)
proj.filter<- saveArchRProject(ArchRProj = proj.filter, outputDirectory = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/ATAC_filtered_First", load = TRUE)

head(proj.filter@cellColData)

p9 <- plotEmbedding(ArchRProj = proj.filter, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un", 
                    embedding = "UMAP_harmony",
                    size = 0.5)

plotPDF(p9, name = "Plot-UMAP-celltype.pdf", 
        ArchRProj = proj.filter, 
        addDOC = FALSE, width = 8, height = 8)

###提取软骨细胞
# 细胞类型提取
Chondrocytes<- BiocGenerics::which(proj.filter$predictedGroup_Un %in% "Chondrocytes")
cellsSample <- proj.filter$cellNames[Chondrocytes]
proj.filter_Ch=proj.filter[cellsSample, ]
saveArchRProject(proj.filter_Ch,outputDirectory = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/ATAC_filtered_Ch",overwrite = TRUE)

proj.filter_Ch


p10 <- plotEmbedding(ArchRProj = proj.filter, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un", 
                    embedding = "UMAP_LSI",
                    size = 0.5)

plotPDF(p10, name = "Plot-UMAP_LSI-celltype.pdf", 
        ArchRProj = proj.filter, 
        addDOC = FALSE, width = 8, height = 8)

p11 <- plotEmbedding(ArchRProj = proj.filter, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un", 
                    embedding = "TSNE_LSI",
                    size = 0.5)

plotPDF(p11, name = "Plot-TSNE_LSI-celltype.pdf", 
        ArchRProj = proj.filter, 
        addDOC = FALSE, width = 8, height = 8)
p12 <- plotEmbedding(ArchRProj = proj.filter, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un", 
                    embedding = "TSNE_Harmony",
                    size = 0.5)

plotPDF(p12, name = "Plot-TSNE_Harmony-celltype.pdf", 
        ArchRProj = proj.filter, 
        addDOC = FALSE, width = 8, height = 8)


proj.filter_Ch <- loadArchRProject("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/ATAC_filtered_Ch")

####新的降维聚类
######################
#降维---dimension reduction 
proj.filter_Ch <- addIterativeLSI(ArchRProj = proj.filter_Ch,
                               useMatrix = "TileMatrix",
                               name = "IterativeLSI",
                               iterations = 4,#要执行的LSI迭代次数
                               clusterParams = list(resolution = 4,
                                                    sampleCells =10000,
                                                    n.start = 10),
                               varFeatures = 50000,
                               dimsToUse = 1:20,
                               force = TRUE,
                               seed=10)
head(proj.filter_Ch@cellColData)
#使用harmony去除批次效应
 proj.filter_Ch <- addHarmony(ArchRProj = proj.filter_Ch,
                          reducedDims = "IterativeLSI",
                          name = "Harmony",#harmony去除批次
                          groupBy = "Sample",#因为是多样本，所以按照样本去除批次
                          force = T)
###聚类分析
proj.filter_Ch <- addClusters(input = proj.filter_Ch, 
                           reducedDims = "Harmony", 
                           method = "Seurat",
                           name = "Clusters", 
                           resolution = 2, 
                           force=TRUE, 
                           seed = 11)
table(proj.filter_Ch@cellColData$Clusters)#查看每个样本细胞数
# C1 C10 C11 C12  C2  C3  C4  C5  C6  C7  C8  C9 
#587 227 531 689 650 317 607 372 276 371 357 222
table(proj.filter_Ch@cellColData$Clusters_2)
# C1 C10 C11 C12 C13 C14  C2  C3  C4  C5  C6  C7  C8  C9 
#224 217 515 225 287 858 412 291 553 391 354 279 354 246 
proj.filter_Ch <- addUMAP(ArchRProj = proj.filter_Ch, 
                       reducedDims = "Harmony", 
                       name = "UMAP_harmony",
                       nNeighbors = 30,
                       minDist = 0.5, 
                       metric = "cosine",
                       force=TRUE,
                       seed=12)

proj.filter_Ch <- addUMAP(ArchRProj = proj.filter_Ch, 
                       reducedDims = "IterativeLSI", 
                       name = "UMAP_LSI",
                       nNeighbors = 30,
                       minDist = 0.5, 
                       metric = "cosine",
                       force=TRUE,
                       seed=12)
proj.filter_Ch <- addTSNE(
    ArchRProj = proj.filter_Ch, 
    reducedDims = "IterativeLSI", 
    name = "TSNE_LSI",
    perplexity = 30,force = TRUE
)

proj.filter_Ch <- addTSNE(
    ArchRProj = proj.filter_Ch, 
    reducedDims = "Harmony", 
    name = "TSNE_Harmony",
    perplexity = 30,force = TRUE
)
proj.filter_Ch@embeddings

#聚类后一些基本聚类图的可视化---UMAP图
p7 <- plotEmbedding(ArchRProj = proj.filter_Ch, 
                    colorBy = "cellColData", 
                    name = "Sample", 
                    embedding = "UMAP_LSI",
                    size = 0.5)

p8 <- plotEmbedding(ArchRProj = proj.filter_Ch, 
                    colorBy = "cellColData", 
                    name = "Clusters_2", 
                    embedding = "UMAP_LSI",
                    size = 0.5)

plotPDF(p7,p8, name = "Plot-UMAP_LSI-Sample-Clusters.new.pdf", 
        ArchRProj = proj.filter_Ch, 
        addDOC = FALSE, width = 8, height = 8)

p7 <- plotEmbedding(ArchRProj = proj.filter_Ch, 
                    colorBy = "cellColData", 
                    name = "Sample", 
                    embedding = "UMAP_harmony",
                    size = 0.5)

p8 <- plotEmbedding(ArchRProj = proj.filter_Ch, 
                    colorBy = "cellColData", 
                    name = "Clusters_2", 
                    embedding = "UMAP_harmony",
                    size = 0.5)

plotPDF(p7,p8, name = "Plot-UMAP-Sample-Clusters.new.pdf", 
        ArchRProj = proj.filter_Ch, 
        addDOC = FALSE, width = 8, height = 8)


setwd("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/ATAC_filtered_Ch")
seRNA <- readRDS("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/scRNA_Ch.rds")
seRNA
colnames(colData(seRNA))
table(colData(seRNA)$celltype)
library(devtools)
load_all("/home/yiyuan/software/ArchR-development")
proj.filter_Ch <- addGeneIntegrationMatrix(
    ArchRProj = proj.filter_Ch, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE,
    groupRNA = "celltype",
    nameCell = "predictedCell_Un_2",
    nameGroup = "predictedGroup_Un_2",
    nameScore = "predictedScore_Un_2"
)

head(proj.filter_Ch@cellColData)
proj.filter_Ch@embeddings
#List of length 4
#names(4): UMAP_harmony UMAP_LSI TSNE_LSI TSNE_Harmony
p9 <- plotEmbedding(ArchRProj = proj.filter_Ch, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un_2", 
                    embedding = "UMAP_harmony",
                    size = 0.5)

plotPDF(p9, name = "Plot-UMAP-celltype_2.pdf", 
        ArchRProj = proj.filter_Ch, 
        addDOC = FALSE, width = 8, height = 8)

p9 <- plotEmbedding(ArchRProj = proj.filter_Ch, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un_2", 
                    embedding = "UMAP_LSI",
                    size = 0.5)

plotPDF(p9, name = "Plot-UMAP_LSI-celltype_2.pdf", 
        ArchRProj = proj.filter_Ch, 
        addDOC = FALSE, width = 8, height = 8)

p9 <- plotEmbedding(ArchRProj = proj.filter_Ch, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un_2", 
                    embedding = "TSNE_LSI",
                    size = 0.5)

plotPDF(p9, name = "Plot-TSNE_LSI-celltype_2.pdf", 
        ArchRProj = proj.filter_Ch, 
        addDOC = FALSE, width = 8, height = 8)

p9 <- plotEmbedding(ArchRProj = proj.filter_Ch, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un_2", 
                    embedding = "TSNE_Harmony",
                    size = 0.5)

plotPDF(p9, name = "Plot-TSNE_Harmony-celltype_2.pdf", 
        ArchRProj = proj.filter_Ch, 
        addDOC = FALSE, width = 8, height = 8)

cM <- as.matrix(confusionMatrix(proj.filter_Ch$Clusters_2, proj.filter_Ch$predictedGroup_Un_2))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM))
unique(unique(proj.filter_Ch$predictedGroup_Un_2))

#From scRNA
cCh1 <- paste0(paste0(1), collapse="|")
cCh1
cCh2 <- paste0(paste0(2), collapse="|")
cCh2
cCh3 <- paste0(paste0(3), collapse="|")
cCh3
cCh4 <- paste0(paste0(4), collapse="|")
cCh4
cCh5 <- paste0(paste0(5), collapse="|")
cCh5
cCh6 <- paste0(paste0(6), collapse="|")
cCh6
cCh7 <- paste0(paste0(7), collapse="|")
cCh7
cCh8 <- paste0(paste0(8), collapse="|")
cCh8

clustCh1 <- rownames(cM)[grep(cCh1, preClust)]
clustCh1
clustCh2 <- rownames(cM)[grep(cCh2, preClust)]
clustCh2
clustCh3 <- rownames(cM)[grep(cCh3, preClust)]
clustCh3
clustCh4 <- rownames(cM)[grep(cCh4, preClust)]
clustCh4
clustCh5 <- rownames(cM)[grep(cCh5, preClust)]
clustCh5
clustCh6 <- rownames(cM)[grep(cCh6, preClust)]
clustCh6
clustCh7 <- rownames(cM)[grep(cCh7, preClust)]
clustCh7
clustCh8 <- rownames(cM)[grep(cCh8, preClust)]
clustCh8

rnaCh1<- colnames(seRNA)[grep(cCh1, colData(seRNA)$celltype)]
head(rnaCh1)
rnaCh2<- colnames(seRNA)[grep(cCh2, colData(seRNA)$celltype)]
head(rnaCh2)
rnaCh3<- colnames(seRNA)[grep(cCh3, colData(seRNA)$celltype)]
head(rnaCh3)
rnaCh4<- colnames(seRNA)[grep(cCh4, colData(seRNA)$celltype)]
head(rnaCh4)
rnaCh5<- colnames(seRNA)[grep(cCh5, colData(seRNA)$celltype)]
head(rnaCh5)
rnaCh6<- colnames(seRNA)[grep(cCh6, colData(seRNA)$celltype)]
head(rnaCh6)
rnaCh7<- colnames(seRNA)[grep(cCh7, colData(seRNA)$celltype)]
head(rnaCh7)
rnaCh8<- colnames(seRNA)[grep(cCh8, colData(seRNA)$celltype)]
head(rnaCh8)

groupList <- SimpleList(
    Ch1 = SimpleList(
        ATAC = proj.filter_Ch$cellNames[proj.filter_Ch$Clusters_2 %in% clustCh1],
        RNA = rnaCh1
    ),
    Ch2 = SimpleList(
        ATAC = proj.filter_Ch$cellNames[proj.filter_Ch$Clusters_2 %in% clustCh2],
        RNA = rnaCh2
    ),
    Ch3 = SimpleList(
        ATAC = proj.filter_Ch$cellNames[proj.filter_Ch$Clusters_2 %in% clustCh3],
        RNA = rnaCh3
    ),
    Ch4 = SimpleList(
        ATAC = proj.filter_Ch$cellNames[proj.filter_Ch$Clusters_2 %in% clustCh4],
        RNA = rnaCh4
    ),
    Ch5 = SimpleList(
        ATAC = proj.filter_Ch$cellNames[proj.filter_Ch$Clusters_2 %in% clustCh5],
        RNA = rnaCh5
    ),
    Ch7 = SimpleList(
        ATAC = proj.filter_Ch$cellNames[proj.filter_Ch$Clusters_2 %in% clustCh7],
        RNA = rnaCh7
    )
)
#~5 minutes
proj.filter_Ch <- addGeneIntegrationMatrix(
    ArchRProj = proj.filter_Ch, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = FALSE, 
    groupList = groupList,
    groupRNA = "celltype",
    nameCell = "predictedCell_Co",
    nameGroup = "predictedGroup_Co",
    nameScore = "predictedScore_Co"
)

pal <- paletteDiscrete(values = colData(seRNA)$celltype)
 pal
#      Ch1       Ch2       Ch3       Ch4       Ch5       Ch6       Ch7       Ch8 
#"#D51F26" "#272E6A" "#208A42" "#89288F" "#F47D2B" "#FEE500" "#8A9FD1" "#C06CAB" 

p1 <- plotEmbedding(
    proj.filter_Ch, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un_2", 
    embedding = "UMAP_harmony",
    size = 0.5,
    pal = pal
)

p2 <- plotEmbedding(
    proj.filter_Ch, 
    colorBy = "cellColData", 
    name = "predictedGroup_Co",
    embedding = "UMAP_harmony",
    size = 0.5,
    pal = pal
)

plotPDF(p1,p2, name = "Plot-UMAP_harmony-RNA-Integration.pdf", ArchRProj = proj.filter_Ch, addDOC = FALSE, width = 5, height = 5)
#UMAP_harmony,UMAP_LSI,TSNE_LSI,TSNE_Harmony
p1 <- plotEmbedding(
    proj.filter_Ch, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un_2", 
    embedding = "UMAP_LSI",
    size = 0.5,
    pal = pal
)

p2 <- plotEmbedding(
    proj.filter_Ch, 
    colorBy = "cellColData", 
    name = "predictedGroup_Co",
    embedding = "UMAP_LSI",
    size = 0.5,
    pal = pal
)

plotPDF(p1,p2, name = "Plot-UMAP_LSI-RNA-Integration.pdf", ArchRProj = proj.filter_Ch, addDOC = FALSE, width = 5, height = 5)

proj.filter_Ch <- saveArchRProject(ArchRProj = proj.filter_Ch, outputDirectory = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-proj.filter_Ch", load = TRUE)

projHeme3 <- addGeneIntegrationMatrix(
    ArchRProj = proj.filter_Ch, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    force= TRUE,
    groupList = groupList,
    groupRNA = "celltype",
    nameCell = "predictedCell",
    nameGroup = "predictedGroup",
    nameScore = "predictedScore"
)

getAvailableMatrices(projHeme3)
#[1] "GeneIntegrationMatrix" "GeneScoreMatrix"       "TileMatrix" 
projHeme3 <- addImputeWeights(projHeme3)
cM <- confusionMatrix(projHeme3$Clusters_2, projHeme3$predictedGroup)
labelOld <- rownames(cM)
labelOld
# [1] "C6"  "C13" "C4"  "C14" "C1"  "C12" "C8"  "C5"  "C2"  "C7"  "C11" "C9" 
#[13] "C10" "C3"
labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew
# [1] "Ch2" "Ch4" "Ch2" "Ch4" "Ch2" "Ch1" "Ch5" "Ch2" "Ch3" "Ch1" "Ch7" "Ch1"
#[13] "Ch5" "Ch1"
projHeme3$Clusters2 <- mapLabels(projHeme3$Clusters_2, newLabels = labelNew, oldLabels = labelOld)

p1 <- plotEmbedding(projHeme3, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP_LSI",size = 0.5)
plotPDF(p1, name = "Plot-UMAP_LSI-Remap-Clusters.pdf", ArchRProj = projHeme3, addDOC = FALSE, width = 5, height = 5)
projHeme3 <- saveArchRProject(ArchRProj = projHeme3, outputDirectory = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme3", load = TRUE)
head(projHeme3@cellColData)

setwd("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch")
####peaks 分组call peak
projHeme3 <- loadArchRProject("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme3")

head(projHeme3@cellColData)
getAvailableMatrices(projHeme3)
getOutputDirectory(ArchRProj = projHeme3)

library(dplyr)
projHeme3$group<- recode(projHeme3$Sample,
                        "HAC0416" = "HA",
                        "HAEC1" = "HA",
                        "OAC_0510" = "OA",
                        "OAC_ZXG" = "OA",
                        .default = NA_character_)
projHeme3@cellColData %>% head()

projHeme3$Cgroup <- paste0(projHeme3$group, "_", projHeme3$Clusters2)
table(projHeme3$Cgroup)
#HA_Ch1 HA_Ch2 HA_Ch3 HA_Ch4 HA_Ch5 HA_Ch7 OA_Ch1 OA_Ch2 OA_Ch3 OA_Ch4 OA_Ch5 
#   282    493     74    246    106     58    759   1029    338    899    465 
#OA_Ch7 
#   457 
projHeme3 <- saveArchRProject(ArchRProj = projHeme3, outputDirectory = "Save-ProjHeme3_new", load = TRUE)
library(BSgenome.Hsapiens.UCSC.hg38)
#Making Pseudo-bulk Replicates
projHeme4 <- addGroupCoverages(ArchRProj = projHeme3, groupBy = "Cgroup")
#Calling Peaks with ArchR
pathToMacs2 <- "/home/yiyuan/micromamba/envs/depipe/bin/macs2"

projHeme4 <- addReproduciblePeakSet(
    ArchRProj = projHeme4, 
    groupBy = "Cgroup", 
    pathToMacs2 = pathToMacs2
)
#        Group nCells nCellsUsed nReplicates nMin nMax maxPeaks
#HA_Ch1 HA_Ch1    282        282           2  100  182   141000
#HA_Ch2 HA_Ch2    493        493           2  124  369   150000
#HA_Ch3 HA_Ch3     74         63           2   40   40    31500
#HA_Ch4 HA_Ch4    246        246           2   44  202   123000
#HA_Ch5 HA_Ch5    106        106           2   40   66    53000
#HA_Ch7 HA_Ch7     58         56           2   40   40    28000
#OA_Ch1 OA_Ch1    759        679           2  179  500   150000
#OA_Ch2 OA_Ch2   1029        728           2  228  500   150000
#OA_Ch3 OA_Ch3    338        338           2  102  236   150000
#OA_Ch4 OA_Ch4    899        738           2  238  500   150000
#OA_Ch5 OA_Ch5    465        465           2  101  364   150000
#OA_Ch7 OA_Ch7    457        457           2  140  317   150000
getPeakSet(projHeme4)

projHeme4 <- saveArchRProject(ArchRProj = projHeme4, outputDirectory = "Save-ProjHeme4_new", load = TRUE)

projHeme5 <- addPeakMatrix(projHeme4)

getOutputDirectory(ArchRProj = projHeme5)

getAvailableMatrices(projHeme5)
#Identifying Marker Peaks
table(projHeme5$Clusters2)
table(projHeme5$Cgroup)

markerPeaks <- getMarkerFeatures(
    ArchRProj = projHeme5, 
    useMatrix = "PeakMatrix", 
    groupBy = "Cgroup",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerPeaks
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")
markerList

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)
#Identified 12 markers!
# [1] "chr1:222465392-222465892"  "chr10:122282062-122282562"
# [3] "chr19:18650501-18651001"   "chr19:41327928-41328428"  
# [5] "chr2:223613183-223613683"  "chr21:34795604-34796104"  
# [7] "chr21:46084970-46085470"   "chr22:30205633-30206133"  
# [9] "chr22:46754272-46754772"   "chr5:150298976-150299476" 
#[11] "chr8:1755545-1756045"      "chr8:8227791-8228291"    
plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)
## Plotting ComplexHeatmap!

#markerPeaks1 <- getMarkerFeatures(
#    ArchRProj = projHeme5, 
#    useMatrix = "PeakMatrix", 
#    groupBy = "group",
#  bias = c("TSSEnrichment", "log10(nFrags)"),
#  testMethod = "wilcoxon"
#)


#pSet_HA <- subset(pSet, group == "HA")
#pSet_OA <- subset(pSet, group == "OA")

###火山图
markerTest <- getMarkerFeatures(
  ArchRProj = projHeme5, 
  useMatrix = "PeakMatrix",
  groupBy = "group",
  testMethod = "wilcoxon",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  useGroups = "HA",
  bgdGroups = "OA"
)
pma <- plotMarkers(seMarker = markerTest, name = "HA", cutOff = "FDR <= 0.5 & abs(Log2FC) >= 1", plotAs = "MA")
pv <- plotMarkers(seMarker = markerTest, name = "HA", cutOff = "FDR <= 0.5 & abs(Log2FC) >= 1", plotAs = "Volcano")
plotPDF(pma, pv, name = "HA-vs-OA-Markers-MA-Volcano_25.5.29", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

#Motif and Feature Enrichment with ArchR
projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")

pSet <- getPeakSet(ArchRProj = projHeme5)
pSet
pSet$name <- paste(seqnames(pSet), start(pSet), end(pSet), sep = "_")
pSet$type=names(pSet)
pSet
library(GenomicRanges)
library(tidyr)

# 获取 metadata 并拆分 type 列
mcols(pSet) <- mcols(pSet) %>%
  as.data.frame() %>%
  separate(type, into = c("group", "Clusters2"), sep = "_", extra = "merge", fill = "right")

pSet$type=names(pSet)
pSet
write.table(pSet, 
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme4_new/pSet_peaks.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)
matches <- getMatches(ArchRProj = projHeme5, name = "Motif")
rownames(matches) <- paste(seqnames(matches), start(matches), end(matches), sep = "_")
matches <- matches[pSet$name]

gr <- GRanges(seqnames = c("chr11"), ranges = IRanges(start = c(46380184), end = c(46380684)))
queryHits <- queryHits(findOverlaps(query = pSet, subject = gr, type = "within"))
colnames(matches)[which(assay(matches[queryHits,]))]
# [1] "TFAP2C_3"     "AHRR_17"      "AHR_28"       "NEUROD1_63"   "ARNT2_76"    
# [6] "MGA_103"      "FOSL2_105"    "NFE2L2_115"   "NFE2_119"     "JUND_124"    
#[11] "FOS_137"      "JUNB_139"     "FOSL1_142"    "JUN_143"      "ZFX_158"     
#[16] "ZNF263_159"   "ZNF76_164"    "KLF6_165"     "PATZ1_172"    "KLF5_175"    
#[21] "MAZ_178"      "KLF1_179"     "KLF3_184"     "PLAGL1_190"   "EGR1_195"    
#[26] "EGR2_196"     "CTCFL_198"    "RREB1_200"    "KLF16_205"    "EGR4_207"    
#[31] "KLF4_208"     "ZNF740_209"   "ZIC5_210"     "ZIC1_213"     "ZBTB7B_216"  
#[36] "ZNF281_221"   "ZNF148_222"   "KLF15_223"    "ZNF219_229"   "SP2_232"     
#[41] "ZBTB49_235"   "SP3_247"      "INSM1_248"    "KLF14_251"    "ZNF354C_256" 
#[46] "ZBTB7A_258"   "ZBTB7C_265"   "WT1_266"      "SP1_267"      "SP6_275"     
#[51] "SP5_279"      "SPI1_322"     "FOXP3_348"    "FOXP1_353"    "FOXO3_354"   
#[56] "FOXP4_358"    "FOXO1_361"    "FOXK1_364"    "FOXG1_368"    "FOXE3_376"   
#[61] "FOXO6_379"    "KIAA0415_381" "FOXD1_382"    "MEIS3_411"    "NKX32_432"   
#[66] "MEIS1_486"    "IRF1_629"     "MBD2_644"     "MECP2_645"    "SMARCC1_651" 
#[71] "TP63_704"     "TP73_705"     "NFIC_740"     "TBX15_782"    "TBX19_787"   
#[76] "T_789"        "PRKRIR_799"   "NR0B1_811"    "PURA_813"     "ZNF350_817"  
#[81] "FOXD4_830"    "FOXD4L1_831"  "FOXD4L3_832"  "FOXD4L5_833"  "FOXD4L6_834" 
#[86] "KLF2_846"     "KLF11_847"    "FOXJ1_853"    "FOXD4L4_856"  "FOXD4L2_857" 
#[91] "TBX18_869"    "TBX22_870"

gr <- GRanges(seqnames = c("chr11"), ranges = IRanges(start = c(46378516), end = c(46379016)))
queryHits <- queryHits(findOverlaps(query = pSet, subject = gr, type = "within"))
colnames(matches)[which(assay(matches[queryHits,]))]
# [1] "TFAP2D_2"    "TFAP2C_3"    "SREBF1_22"   "TCF15_46"    "MYOD1_48"   
# [6] "FERD3L_58"   "EBF1_67"     "FIGLA_88"    "SCXB_93"     "SCXA_96"    
#[11] "SREBF2_99"   "ZEB1_157"    "ZFX_158"     "ZNF263_159"  "ZNF76_164"  
#[16] "ZFY_166"     "REST_168"    "PATZ1_172"   "MAZ_178"     "SP8_226"    
#[21] "SP2_232"     "SP7_241"     "KLF14_251"   "SP1_267"     "SP9_283"    
#[26] "ELF2_326"    "TFCP2_392"   "NKX32_432"   "TGIF2_441"   "MEIS2_471"  
#[31] "TGIF2LX_492" "PKNOX1_497"  "PKNOX2_513"  "TGIF2LY_538" "THRA_672"   
#[36] "NR1I3_680"   "THRB_684"    "TP63_704"    "TP73_705"    "RELA_722"   
#[41] "THAP1_798"   "LMO2_808"    "MYF5_842"    "KLF11_847"   "TFCP2L1_858"
gr <- GRanges(seqnames = c("chr11"), ranges = IRanges(start = c(46379143), end = c(46379643)))
queryHits <- queryHits(findOverlaps(query = pSet, subject = gr, type = "within"))
colnames(matches)[which(assay(matches[queryHits,]))]
# [1] "MYOD1_48"    "TAL1_62"     "NEUROD1_63"  "EBF1_67"     "ZFX_158"    
# [6] "ZNF263_159"  "ZIC2_162"    "KLF6_165"    "ZFY_166"     "MZF1_171"   
#[11] "PATZ1_172"   "KLF5_175"    "MAZ_178"     "GLI1_186"    "CTCFL_198"  
#[16] "KLF4_208"    "ZNF148_222"  "KLF15_223"   "SP2_232"     "ZNF354C_256"
#[21] "ZBTB7A_258"  "ZNF238_260"  "ZBTB42_261"  "WT1_266"     "SP1_267"    
#[26] "NFYA_288"    "DNMT1_301"   "E2F6_317"    "ELF4_323"    "EHF_333"    
#[31] "ELF5_334"    "ETS2_340"    "ELF3_342"    "HOXA13_418"  "TGIF2_441"  
#[36] "MEIS2_471"   "NKX21_476"   "HHEX_491"    "TGIF2LX_492" "HOXB13_494" 
#[41] "PKNOX1_497"  "PKNOX2_513"  "TGIF2LY_538" "ESR1_661"    "HNF4A_662"  
#[46] "VDR_663"     "NR3C1_666"   "NR5A2_667"   "ESR2_678"    "NR1I2_682"  
#[51] "NR2C2_693"   "NFATC4_715"  "RFX2_724"    "SMAD1_744"   "SMAD2_745"  
#[56] "PURA_813"    "TAL2_822"    "ATOH1_844"  
gr <- GRanges(seqnames = c("chr11"), ranges = IRanges(start = c(46381440), end = c(46381940)))
queryHits <- queryHits(findOverlaps(query = pSet, subject = gr, type = "within"))
colnames(matches)[which(assay(matches[queryHits,]))]
# [1] "TFAP2D_2"    "TFAP2C_3"    "USF2_26"     "USF1_60"     "EBF1_67"    
# [6] "ATF3_132"    "ZFX_158"     "ZNF263_159"  "KLF6_165"    "MZF1_171"   
#[11] "PATZ1_172"   "KLF5_175"    "CTCF_177"    "MAZ_178"     "KLF1_179"   
#[16] "SP4_180"     "KLF7_189"    "EGR1_195"    "EGR2_196"    "CTCFL_198"  
#[21] "RREB1_200"   "KLF16_205"   "EGR4_207"    "KLF4_208"    "ZNF740_209" 
#[26] "ZNF281_221"  "ZNF148_222"  "KLF15_223"   "SP8_226"     "ZNF219_229" 
#[31] "SP2_232"     "SP7_241"     "ZNF524_243"  "SP3_247"     "KLF14_251"  
#[36] "ZBTB7A_258"  "WT1_266"     "SP1_267"     "ZKSCAN4_273" "SP6_275"    
#[41] "ZKSCAN3_276" "SP5_279"     "SP9_283"     "CGBP_298"    "DNMT1_301"  
#[46] "E2F6_317"    "GRHL1_391"   "MBD2_644"    "TP63_704"    "TP73_705"   
#[51] "NR0B1_811"   "PURA_813"    "KLF2_846"    "KLF11_847"   "TFCP2L1_858"
#[56] "SMAD5_866

### Motif Enrichment in Differential Peaks
motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = projHeme5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC >= 0.5"
  )

df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)
write.table(df, 
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme4_new/motifsUp.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)
ggUp <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(P-adj) Motif Enrichment") + 
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
png(file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme4_new/Rank_Sorted_TFs_Enriched_ggUp.png", width = 3700, height = 4200, res = 300)
ggUp
dev.off()

motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = projHeme5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.1 & Log2FC <= -0.5"
  )

df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)
write.table(df, 
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme4_new/motifsDo.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)
ggDo <- ggplot(df, aes(rank, mlog10Padj, color = mlog10Padj)) + 
  geom_point(size = 1) +
  ggrepel::geom_label_repel(
        data = df[rev(seq_len(30)), ], aes(x = rank, y = mlog10Padj, label = TF), 
        size = 1.5,
        nudge_x = 2,
        color = "black"
  ) + theme_ArchR() + 
  ylab("-log10(FDR) Motif Enrichment") +
  xlab("Rank Sorted TFs Enriched") +
  scale_color_gradientn(colors = paletteContinuous(set = "comet"))
png(file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme4_new/Rank_Sorted_TFs_Enriched_ggDo.png", width = 4000, height = 4200, res = 300)
ggDo
dev.off()

plotPDF(ggUp, ggDo, name = "HA-vs-OA-Markers-Motifs-Enriched", width = 10, height = 10, ArchRProj = projHeme5, addDOC = FALSE)

###Motif Enrichment in Marker Peaks
enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerPeaks,
    ArchRProj = projHeme5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.4 & Log2FC >= 0.5"
  )
heatmapEM <- plotEnrichHeatmap(enrichMotifs,transpose = TRUE)
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-HeatmapFDR <= 0.4 & Log2FC >= 0.5", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)
#Motif enrichment in arbitrary regions   Plotting motif logos
pwm <- getPeakAnnotation(projHeme5, "Motif")$motifs[["KLF10_826"]]
pwm
PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}

ppm <- PWMatrixToProbMatrix(pwm)
ppm
colSums(ppm) %>% range
library(ggseqlogo)
png("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme5_new/Plots/logo_plot_KLF10_826.png", width = 1200, height = 800, res = 150)
ggseqlogo(ppm, method = "bits")
dev.off()
#> motifs <- c("JUNB","FOS","KLF10","EGR1")
#> markerMotifs <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
#markerMotifs
# [1] "z:KLF10_826"          "z:EGR1_195"           "z:FOSL1_142"         
# [4] "z:JUNB_139"           "z:FOS_137"            "z:FOSB_121"  
pwm <- getPeakAnnotation(projHeme5, "Motif")$motifs[["FOS_137"]]
pwm
PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}

ppm <- PWMatrixToProbMatrix(pwm)
ppm
colSums(ppm) %>% range
library(ggseqlogo)
png("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme5_new/Plots/logo_plot_FOS_137.png", width = 1200, height = 800, res = 150)
ggseqlogo(ppm, method = "bits")
dev.off()



# [1] "z:KLF10_826"          "z:EGR1_195"           "z:JUN_143"           
# [4] "z:FOSL1_142"          "z:JUNB_139"           "z:FOS_137"           
# [7] "z:JUND_124"           "z:FOSB_121"           "z:FOSL2_105" 
pwm <- getPeakAnnotation(projHeme5, "Motif")$motifs[["FOSB_121"]]
pwm
PWMatrixToProbMatrix <- function(x){
  if (class(x) != "PWMatrix") stop("x must be a TFBSTools::PWMatrix object")
  m <- (exp(as(x, "matrix"))) * TFBSTools::bg(x)/sum(TFBSTools::bg(x))
  m <- t(t(m)/colSums(m))
  m
}

ppm <- PWMatrixToProbMatrix(pwm)
ppm
colSums(ppm) %>% range
library(ggseqlogo)
png("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme5_new/Plots/logo_plot_FOSB_121.png", width = 1200, height = 800, res = 150)
ggseqlogo(ppm, method = "bits")
dev.off()
#####软骨细胞的HA特异性TF的偏差分析和足迹分析：JUNB\FOS\KLF10\EGR1\
####Motif Deviations

if("Motif" %ni% names(projHeme5@peakAnnotation)){
    projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")
}
projHeme5 <- addBgdPeaks(projHeme5)
## Identifying Background Peaks!
projHeme5 <- addDeviationsMatrix(
  ArchRProj = projHeme5, 
  peakAnnotation = "Motif",
  force = TRUE
)
plotVarDev <- getVarDeviations(projHeme5, name = "MotifMatrix", plot = TRUE)
plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 10, height = 8, ArchRProj = projHeme5, addDOC = FALSE)
## Plotting Ggplot!

motifs <- c("JUNB","FOS","KLF10","EGR1")
markerMotifs <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs

markerMotifs <- grep("z:", markerMotifs, value = TRUE)
#markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs

p <- plotGroups(ArchRProj = projHeme5, 
  groupBy = "group", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(projHeme5)
)
plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation_25.5.30", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)


####Footprint Analysis

library(BSgenome.Hsapiens.UCSC.hg38)
if(is.null(projHeme5@projectMetadata$GroupCoverages$Clusters2)){
  projHeme5 <- addGroupCoverages(ArchRProj = projHeme5, groupBy = "group")
}

motifPositions <- getPositions(projHeme5)
motifs <- c("JUNB","FOS","KLF10","EGR1","JUN","FOSL2","FOSL1","JUND","FOSB")
selectedMotifs <- grep(paste(motifs, collapse = "|"), names(motifPositions), value = TRUE)
seFoot <- getFootprints(
  ArchRProj = projHeme5,
  positions = motifPositions[selectedMotifs],
  groupBy = "group"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projHeme5, 
  normMethod = "None",
  plotName = "Footprints-No-Normalization",
  addDOC = FALSE,
  smoothWindow = 5
)

#####Co-accessibility with ArchR
projHeme5 <- addCoAccessibility(
    ArchRProj = projHeme5,
    reducedDims = "IterativeLSI"
)

cA <- getCoAccessibility(
    ArchRProj = projHeme5,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = FALSE
)

cA

write.table(cA, 
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme5_new/CoAccessibility_0.5.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

cB=metadata(cA)[[1]]

write.table(cB, 
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme5_new/CoAccessibility_0.5_peak.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

cA <- getCoAccessibility(
    ArchRProj = projHeme5,
    corCutOff = 0.5,
    resolution = 1,
    returnLoops = TRUE
)

write.table(cA, 
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme5_new/CoAccessibility_0.5_loops.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)


p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "Cgroup", 
    geneSymbol = c("MDK"), 
    upstream = 20000000,
    downstream = 50000,
    loops = getCoAccessibility(projHeme5, corCutOff = 0.2,resolution = 10000,
    returnLoops = TRUE)
)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility_mdk_FOSL1.pdf", 
    ArchRProj = projHeme5, 
    addDOC = FALSE, width = 50, height = 5)


projHeme5 <- loadArchRProject("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme5_new")
projHeme5 <- addPeak2GeneLinks(
    ArchRProj = projHeme5,
    reducedDims = "IterativeLSI"
)

p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "group", 
    geneSymbol = c("JUNB","FOS","KLF10","EGR1"), 
    upstream = 150000,
    downstream = 100000,
    loops = getPeak2GeneLinks(projHeme5, corCutOff = 0.4,resolution = 1000,
    returnLoops = TRUE)
)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility_2025.6.5.pdf", 
    ArchRProj = projHeme5, 
    addDOC = FALSE, width = 20, height = 8)


p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "group", 
    geneSymbol =  c("JUNB","FOS","KLF10","EGR1"), 
    upstream = 100000,
    downstream = 80000,
    loops = getPeak2GeneLinks(projHeme5, corCutOff = 0.4,resolution = 100,returnLoops = TRUE)
)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks_2025.6.7_klf10_fos.pdf", 
    ArchRProj = projHeme5, 
    addDOC = FALSE, width = 15, height = 8)


p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "group", 
    geneSymbol =  c("JUNB"), 
    upstream = 120000,
    downstream = 120000,
    loops = getPeak2GeneLinks(projHeme5, corCutOff = 0.4,resolution = 100,returnLoops = TRUE)
)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks_2025.6.7_JUNB.2.pdf", 
    ArchRProj = projHeme5, 
    addDOC = FALSE, width = 15, height = 8)



projHeme5 <- saveArchRProject(ArchRProj = projHeme5, outputDirectory = "Save-ProjHeme5_new", load = TRUE)

coverageFiles <- getGroupBW(ArchRProj = projHeme5, groupBy = "Cgroup")


pSet <- getPeakSet(ArchRProj = projHeme5)
pSet
pSet$name <- paste(seqnames(pSet), start(pSet), end(pSet), sep = "_")
pSet$type=names(pSet)
pSet
library(GenomicRanges)
library(tidyr)
# 获取 metadata 并拆分 type 列
mcols(pSet) <- mcols(pSet) %>%
  as.data.frame() %>%
  separate(type, into = c("group", "Clusters2"), sep = "_", extra = "merge", fill = "right")

pSet$type=names(pSet)
pSet
pSet_HA <- subset(pSet, group == "HA")
pSet_OA <- subset(pSet, group == "OA")


write.table(pSet_HA, 
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme5_new/pSet_HA_peaks.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

write.table(pSet_OA, 
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme5_new/pSet_OA_peaks.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

            

#####ArchR对象转换Seurat对象
GSM_se <- getMatrixFromProject(
    ArchRProj = projHeme5,
    useMatrix <- "GeneScoreMatrix"
)

GSM_se

library(SingleCellExperiment)
library(Seurat)
GSM_sce <- as(GSM_se, "SingleCellExperiment")
counts <- as.matrix(assay(GSM_sce, "GeneScoreMatrix"))
## Warning in asMethod(object): sparse->dense coercion: allocating vector of size
## 1.8 GiB
assays(GSM_sce)$counts <- counts
libSizes <- colSums(counts)
sizeFactors <- libSizes/mean(libSizes)
assays(GSM_sce)$logcounts <- log2(t(t(counts)/sizeFactors) + 1)
rownames(GSM_sce) <- rowData(GSM_sce)$name
seuratObj <- as.Seurat(GSM_sce, counts = "counts", data = "logcounts")
## Warning: Feature names cannot have underscores ('_'), replacing with dashes
## ('-')
## Warning: Feature names cannot have underscores ('_'), replacing with dashes
## ('-')
## Warning: Feature names cannot have underscores ('_'), replacing with dashes
## ('-')
seuratObj
seuratObj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(projHeme5@embeddings$UMAP_LSI$df), key = "UMAP_", assay = DefaultAssay(seuratObj))


png(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Ch/new/marker_gene/second_celltype_umap.png", width = 3000, height = 2500, res = 300)
DimPlot(
  seuratObj,
  reduction = "UMAP_LSI", 
  group.by = "predictedGroup", 
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
saveRDS(seuratObj,file="/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Ch/new/marker_gene/scATAC_Ch_resubtype.final.rds")


