library(ArchR)
library(pheatmap) 
library(Rsamtools)
library(scran) 
library(dplyr)
library(Seurat)
library(patchwork)
library(SingleCellExperiment) 
library(ComplexHeatmap)
library(ggplot2)
library(stringr)
library(EnsDb.Hsapiens.v86) 
library(viridis) 
addArchRGenome("hg38")#设置基因组（人hg38）
addArchRThreads(threads = 10)#设置线程
###八例scATAC数据 
###软骨分析注释
setwd("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch")
#input文件路径,只需要样本的atac_fragments.tsv.gz文件
input.file.list = c('/mnt/8w/data4/yiyuan/Bone/scATAC/data/20250207/8例完整矩阵/HAC0416/outs/fragments_corrected_dedup_count.tsv.gz',
                    '/mnt/8w/data4/yiyuan/Bone/scATAC/data/20250207/8例完整矩阵/HAEC1/outs/fragments_corrected_dedup_count.tsv.gz',
                    '/mnt/8w/data4/yiyuan/Bone/scATAC/data/20250207/8例完整矩阵/OAC_ZXG/outs/fragments_corrected_dedup_count.tsv.gz',
                    '/mnt/8w/data4/yiyuan/Bone/scATAC/data/20250207/8例完整矩阵/OAC_0510/outs/fragments_corrected_dedup_count.tsv.gz')
#设置样本名
sampleNames = c("HAC0416","HAEC1","OAC_ZXG","OAC_0510")
#创建Arrow文件
ArrowFiles <- createArrowFiles(inputFiles = input.file.list,
                               sampleNames = sampleNames,
                               minTSS = 4, #默认值
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
#构建ArchRProject
projncov <- ArchRProject(ArrowFiles = ArrowFiles,
                         outputDirectory = "ATAC_out", #结果存储
                         copyArrows = TRUE) #

table(projncov@cellColData$Sample)

#####数据质控
#过滤双细胞------------------------------------------------------------
proj.filter <- filterDoublets(projncov)

table(proj.filter@cellColData$Sample)

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
                           resolution = 1, 
                           force=TRUE, 
                           seed = 11)
table(proj.filter@cellColData$Clusters)#查看每个样本细胞数


table(proj.filter@cellColData$Clusters,proj.filter@cellColData$Sample)

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

####降维聚类

#降维---dimension reduction 
proj.filter_Ch <- addIterativeLSI(ArchRProj = proj.filter_Ch,
                               useMatrix = "TileMatrix",
                               name = "IterativeLSI",
                               iterations = 4,
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
                          groupBy = "Sample",
                          force = T)
###聚类分析
proj.filter_Ch <- addClusters(input = proj.filter_Ch, 
                           reducedDims = "Harmony", 
                           method = "Seurat",
                           name = "Clusters", 
                           resolution = 2, 
                           force=TRUE, 
                           seed = 11)
table(proj.filter_Ch@cellColData$Clusters)

table(proj.filter_Ch@cellColData$Clusters_2)

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

labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew

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

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)
## Plotting ComplexHeatmap!


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


gr <- GRanges(seqnames = c("chr11"), ranges = IRanges(start = c(46378516), end = c(46379016)))
queryHits <- queryHits(findOverlaps(query = pSet, subject = gr, type = "within"))
colnames(matches)[which(assay(matches[queryHits,]))]

gr <- GRanges(seqnames = c("chr11"), ranges = IRanges(start = c(46379143), end = c(46379643)))
queryHits <- queryHits(findOverlaps(query = pSet, subject = gr, type = "within"))
colnames(matches)[which(assay(matches[queryHits,]))]

gr <- GRanges(seqnames = c("chr11"), ranges = IRanges(start = c(46381440), end = c(46381940)))
queryHits <- queryHits(findOverlaps(query = pSet, subject = gr, type = "within"))
colnames(matches)[which(assay(matches[queryHits,]))]


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
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.pdf", 
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
    name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks_klf10_fos.pdf", 
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
    name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks_JUNB.2.pdf", 
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

assays(GSM_sce)$counts <- counts
libSizes <- colSums(counts)
sizeFactors <- libSizes/mean(libSizes)
assays(GSM_sce)$logcounts <- log2(t(t(counts)/sizeFactors) + 1)
rownames(GSM_sce) <- rowData(GSM_sce)$name
seuratObj <- as.Seurat(GSM_sce, counts = "counts", data = "logcounts")

seuratObj
seuratObj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(projHeme5@embeddings$UMAP_LSI$df), key = "UMAP_", assay = DefaultAssay(seuratObj))


png(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Ch/new/marker_gene/second_celltype_umap.png", width = 3000, height = 2500, res = 300)
DimPlot(
  seuratObj,
  reduction = "UMAP_LSI", 
  group.by = "predictedGroup", 
  label = TRUE,
  label.size = 6,  
   cols = c(
   'Ch1' = '#ef070f',  
   'Ch2' = '#FAD02E',  
   'Ch3' = '#A8DADC',  
   'Ch4' = '#aba3e7',  
   'Ch5' = '#712820',  
   'Ch6' = '#B7D7A8',  
   'Ch7' = '#D4E157',  
   'Ch8' = '#F2826C' 
))
dev.off()
saveRDS(seuratObj,file="/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Ch/new/marker_gene/scATAC_Ch_resubtype.final.rds")


