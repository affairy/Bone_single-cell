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
# addArchRGenome("mm10")#小鼠
addArchRThreads(threads = 10)#设置线程
###八例scATAC数据 
###滑膜分析注释
setwd("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy")
#input文件路径,只需要样本的atac_fragments.tsv.gz文件
input.file.list = c('/mnt/8w/data4/yiyuan/Bone/scATAC/data/20250207/8例完整矩阵/HAS0416/outs/fragments_corrected_dedup_count.tsv.gz',
                    '/mnt/8w/data4/yiyuan/Bone/scATAC/data/20250207/8例完整矩阵/HAES1/outs/fragments_corrected_dedup_count.tsv.gz',
                    '/mnt/8w/data4/yiyuan/Bone/scATAC/data/20250207/8例完整矩阵/OAS_ZXG/outs/fragments_corrected_dedup_count.tsv.gz',
                    '/mnt/8w/data4/yiyuan/Bone/scATAC/data/20250207/8例完整矩阵/OAS_0510/outs/fragments_corrected_dedup_count.tsv.gz')
#设置样本名，相当于seurat中orig.ident。
sampleNames = c("HAS0416","HAES1","OAS_ZXG","OAS_0510")
ArrowFiles <- createArrowFiles(inputFiles = input.file.list,
                               sampleNames = sampleNames,
                               minTSS = 4, #默认值
                               minFrags = 1000,
                               addTileMat = TRUE,
                               addGeneScoreMat = TRUE,
                               excludeChr = c("chrM", "chrY", "chrX"))
###去除双细胞，计算双细胞评分
doubScores <- addDoubletScores(input = ArrowFiles,
                               k = 20, 
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
                         outputDirectory = "ATAC_out", 
                         copyArrows = TRUE) 

table(projncov@cellColData$Sample)#查看每个样本细胞数

#####数据质控
#过滤双细胞------------------------------------------------------------
proj.filter <- filterDoublets(projncov)

table(proj.filter@cellColData$Sample)#查看每个样本细胞数

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
  proj.i <- proj.filter[proj.filter$Sample == sampleNames[i]]
  
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
                           resolution = 1.5, ###4改为1
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

p7 <- plotEmbedding(ArchRProj = proj.filter, 
                    colorBy = "cellColData", 
                    name = "Sample", 
                    embedding = "UMAP_LSI",
                    size = 0.5)

p8 <- plotEmbedding(ArchRProj = proj.filter, 
                    colorBy = "cellColData", 
                    name = "Clusters", 
                    embedding = "UMAP_LSI",
                    size = 0.5)

plotPDF(p7,p8, name = "Plot-UMAP_LSI-Sample-Clusters.pdf", 
        ArchRProj = proj.filter, 
        addDOC = FALSE, width = 8, height = 8)
markersGS <- getMarkerFeatures(ArchRProj = proj.filter,
                               useMatrix = "GeneScoreMatrix",
                               groupBy = "Clusters",
                               bias = c("TSSEnrichment", "log10(nFrags)"),
                               testMethod = "wilcoxon")
markerList <- getMarkers(markersGS, cutOff = "FDR <= 0.05 & Log2FC >= 1")

saveArchRProject(proj.filter,outputDirectory = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/ATAC_filtered",overwrite = TRUE)
proj.filter <- loadArchRProject("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/ATAC_filtered")
setwd("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy")
seRNA <- readRDS("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/scRNA_first.rds")
seRNA
colnames(colData(seRNA))
table(colData(seRNA)$celltype_cluster_standard)
library(devtools)
load_all("/home/yiyuan/software/ArchR-development")
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

head(proj.filter@cellColData)

p9 <- plotEmbedding(ArchRProj = proj.filter, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un", 
                    embedding = "UMAP_harmony",
                    size = 0.5)

plotPDF(p9, name = "Plot-UMAP-celltype.pdf", 
        ArchRProj = proj.filter, 
        addDOC = FALSE, width = 8, height = 8)

p9 <- plotEmbedding(ArchRProj = proj.filter, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un", 
                    embedding = "UMAP_LSI",
                    size = 0.5)

plotPDF(p9, name = "Plot-UMAP_LSI-celltype.pdf", 
        ArchRProj = proj.filter, 
        addDOC = FALSE, width = 8, height = 8)

p9 <- plotEmbedding(ArchRProj = proj.filter, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un", 
                    embedding = "TSNE_LSI",
                    size = 1)

plotPDF(p9, name = "Plot-TSNE_LSI-celltype.pdf", 
        ArchRProj = proj.filter, 
        addDOC = FALSE, width = 8, height = 8)

###提取成纤维细胞
# 细胞类型提取
Fibroblasts<- BiocGenerics::which(proj.filter$predictedGroup_Un %in% "Fibroblasts")
cellsSample <- proj.filter$cellNames[Fibroblasts]
proj.filter_Fb=proj.filter[cellsSample, ]

p10 <- plotEmbedding(ArchRProj = proj.filter_Fb, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un", 
                    embedding = "UMAP_LSI",
                    size = 0.5)

plotPDF(p10, name = "Plot-UMAP_LSI-celltype_Fb.pdf", 
        ArchRProj = proj.filter, 
        addDOC = FALSE, width = 8, height = 8)
saveArchRProject(proj.filter_Fb,outputDirectory = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/ATAC_filtered_Fb",overwrite = TRUE)

proj.filter <- loadArchRProject("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/ATAC_out")
table(proj.filter$predictedGroup_Un)

saveArchRProject(proj.filter,outputDirectory = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/ATAC_filtered_First",overwrite = TRUE)

proj.filter <- addGeneIntegrationMatrix(
    ArchRProj = proj.filter, 
    useMatrix = "GeneScoreMatrix",
    matrixName = "GeneIntegrationMatrix",
    reducedDims = "IterativeLSI",
    seRNA = seRNA,
    addToArrow = TRUE,
    groupRNA = "celltype_cluster_standard",
    nameCell = "predictedCell_Un",
    nameGroup = "predictedGroup_Un",
    nameScore = "predictedScore_Un"
)

proj.filter_Fb <- loadArchRProject("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/ATAC_filtered_Fb")
head(proj.filter_Fb@cellColData)

####新的降维聚类
######################
#降维---dimension reduction 
proj.filter_Fb <- addIterativeLSI(ArchRProj = proj.filter_Fb,
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
head(proj.filter_Fb@cellColData)
#使用harmony去除批次效应
 proj.filter_Fb <- addHarmony(ArchRProj = proj.filter_Fb,
                          reducedDims = "IterativeLSI",
                          name = "Harmony",#harmony去除批次
                          groupBy = "Sample",#因为是多样本，所以按照样本去除批次
                          force = T)
###聚类分析
proj.filter_Fb <- addClusters(input = proj.filter_Fb, 
                           reducedDims = "Harmony", 
                           method = "Seurat",
                           name = "Clusters_2", 
                           resolution = 2, 
                           force=TRUE, 
                           seed = 11)
table(proj.filter_Fb@cellColData$Clusters)#查看每个样本细胞数

table(proj.filter_Fb@cellColData$Clusters_2)

proj.filter_Fb <- addUMAP(ArchRProj = proj.filter_Fb, 
                       reducedDims = "Harmony", 
                       name = "UMAP_harmony",
                       nNeighbors = 30,
                       minDist = 0.5, 
                       metric = "cosine",
                       force=TRUE,
                       seed=12)

proj.filter_Fb <- addUMAP(ArchRProj = proj.filter_Fb, 
                       reducedDims = "IterativeLSI", 
                       name = "UMAP_LSI",
                       nNeighbors = 30,
                       minDist = 0.5, 
                       metric = "cosine",
                       force=TRUE,
                       seed=12)
proj.filter_Fb <- addTSNE(
    ArchRProj = proj.filter_Fb, 
    reducedDims = "IterativeLSI", 
    name = "TSNE_LSI",
    perplexity = 30,force = TRUE
)

proj.filter_Fb <- addTSNE(
    ArchRProj = proj.filter_Fb, 
    reducedDims = "Harmony", 
    name = "TSNE_Harmony",
    perplexity = 30,force = TRUE
)

proj.filter_Fb@embeddings

#聚类后一些基本聚类图的可视化---UMAP图
p7 <- plotEmbedding(ArchRProj = proj.filter_Fb, 
                    colorBy = "cellColData", 
                    name = "Sample", 
                    embedding = "UMAP_LSI",
                    size = 0.5)

p8 <- plotEmbedding(ArchRProj = proj.filter_Fb, 
                    colorBy = "cellColData", 
                    name = "Clusters_2", 
                    embedding = "UMAP_LSI",
                    size = 0.5)

plotPDF(p7,p8, name = "Plot-UMAP_LSI-Sample-Clusters.new.pdf", 
        ArchRProj = proj.filter_Fb, 
        addDOC = FALSE, width = 8, height = 8)

p7 <- plotEmbedding(ArchRProj = proj.filter_Fb, 
                    colorBy = "cellColData", 
                    name = "Sample", 
                    embedding = "UMAP_harmony",
                    size = 0.5)

p8 <- plotEmbedding(ArchRProj = proj.filter_Fb, 
                    colorBy = "cellColData", 
                    name = "Clusters_2", 
                    embedding = "UMAP_harmony",
                    size = 0.5)

plotPDF(p7,p8, name = "Plot-UMAP_harmony-Sample-Clusters.new.pdf", 
        ArchRProj = proj.filter_Fb, 
        addDOC = FALSE, width = 8, height = 8)


seRNA <- readRDS("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/scRNA_Sy_Fb.rds")
seRNA
colnames(colData(seRNA))
table(colData(seRNA)$celltype)
library(devtools)
load_all("/home/yiyuan/software/ArchR-development")
proj.filter_Fb <- addGeneIntegrationMatrix(
    ArchRProj = proj.filter_Fb, 
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

head(proj.filter_Fb@cellColData)
proj.filter_Fb@embeddings

p9 <- plotEmbedding(ArchRProj = proj.filter_Fb, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un_2", 
                    embedding = "UMAP_harmony",
                    size = 0.5)

plotPDF(p9, name = "Plot-UMAP_harmony-celltype_2.pdf", 
        ArchRProj = proj.filter_Fb, 
        addDOC = FALSE, width = 8, height = 8)

p9 <- plotEmbedding(ArchRProj = proj.filter_Fb, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un_2", 
                    embedding = "UMAP_LSI",
                    size = 0.5)

plotPDF(p9, name = "Plot-UMAP_LSI-celltype_2.pdf", 
        ArchRProj = proj.filter_Fb, 
        addDOC = FALSE, width = 8, height = 8)

p9 <- plotEmbedding(ArchRProj = proj.filter_Fb, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un_2", 
                    embedding = "TSNE_LSI",
                    size = 0.5)

plotPDF(p9, name = "Plot-TSNE_LSI-celltype_2.pdf", 
        ArchRProj = proj.filter_Fb, 
        addDOC = FALSE, width = 8, height = 8)

p9 <- plotEmbedding(ArchRProj = proj.filter_Fb, 
                    colorBy = "cellColData", 
                    name = "predictedGroup_Un_2", 
                    embedding = "TSNE_Harmony",
                    size = 0.5)

plotPDF(p9, name = "Plot-TSNE_Harmony-celltype_2.pdf", 
        ArchRProj = proj.filter_Fb, 
        addDOC = FALSE, width = 8, height = 8)

cM <- as.matrix(confusionMatrix(proj.filter_Fb$Clusters_2, proj.filter_Fb$predictedGroup_Un_2))
preClust <- colnames(cM)[apply(cM, 1 , which.max)]
cbind(preClust, rownames(cM))
unique(unique(proj.filter_Fb$predictedGroup_Un_2))
table(colData(seRNA)$celltype)

#From scRNA
cFb1 <- paste0(paste0(1), collapse="|")
cFb1
cFb2 <- paste0(paste0(2), collapse="|")
cFb2
cFb3 <- paste0(paste0(3), collapse="|")
cFb3
cFb4 <- paste0(paste0(4), collapse="|")
cFb4
cFb5 <- paste0(c(paste0(5),(7)), collapse="|")
cFb5
cFb6 <- paste0(paste0(6), collapse="|")
cFb6
cFb7 <- paste0(paste0(7), collapse="|")
cFb7
clustFb1 <- rownames(cM)[grep(cFb1, preClust)]
clustFb1
#[1] "C17" "C14" "C16" "C2" 
clustFb2 <- rownames(cM)[grep(cFb2, preClust)]
clustFb2
#[1] "C5"  "C6"  "C13" "C7" 
clustFb3 <- rownames(cM)[grep(cFb3, preClust)]
clustFb3
#[1] "C4" "C9"
clustFb4 <- rownames(cM)[grep(cFb4, preClust)]
clustFb4
#[1] "C8"  "C18" "C10" "C19" "C12" "C11" "C3" 
clustFb5 <- rownames(cM)[grep(cFb5, preClust)]
clustFb5
#[1] "C15"
clustFb6 <- rownames(cM)[grep(cFb6, preClust)]
clustFb6
#character(0)
clustFb7 <- rownames(cM)[grep(cFb7, preClust)]
clustFb7
#[1] "C1"

rnaFb1<- colnames(seRNA)[grep(cFb1, colData(seRNA)$celltype)]
head(rnaFb1)
rnaFb2<- colnames(seRNA)[grep(cFb2, colData(seRNA)$celltype)]
head(rnaFb2)
rnaFb3<- colnames(seRNA)[grep(cFb3, colData(seRNA)$celltype)]
head(rnaFb3)
rnaFb4<- colnames(seRNA)[grep(cFb4, colData(seRNA)$celltype)]
head(rnaFb4)
rnaFb5<- colnames(seRNA)[grep(cFb5, colData(seRNA)$celltype)]
head(rnaFb5)
rnaFb6<- colnames(seRNA)[grep(cFb6, colData(seRNA)$celltype)]
head(rnaFb6)
rnaFb7<- colnames(seRNA)[grep(cFb7, colData(seRNA)$celltype)]
head(rnaFb7)


groupList <- SimpleList(
    Fb1 = SimpleList(
        ATAC = proj.filter_Fb$cellNames[proj.filter_Fb$Clusters_2 %in% clustFb1],
        RNA = rnaFb1
    ),
    Fb2 = SimpleList(
        ATAC = proj.filter_Fb$cellNames[proj.filter_Fb$Clusters_2 %in% clustFb2],
        RNA = rnaFb2
    ),
    Fb3 = SimpleList(
        ATAC = proj.filter_Fb$cellNames[proj.filter_Fb$Clusters_2 %in% clustFb3],
        RNA = rnaFb3
    ),
    Fb4 = SimpleList(
        ATAC = proj.filter_Fb$cellNames[proj.filter_Fb$Clusters_2 %in% clustFb4],
        RNA = rnaFb4
    ),
    Fb5 = SimpleList(
        ATAC = proj.filter_Fb$cellNames[proj.filter_Fb$Clusters_2 %in% clustFb5],
        RNA = rnaFb5
    )

)

proj.filter_Fb <- addGeneIntegrationMatrix(
    ArchRProj = proj.filter_Fb, 
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
    proj.filter_Fb, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un_2", 
    embedding = "UMAP_harmony",
    size = 0.5,
    pal = pal
)

p2 <- plotEmbedding(
    proj.filter_Fb, 
    colorBy = "cellColData", 
    name = "predictedGroup_Co",
    embedding = "UMAP_harmony",
    size = 0.5,
    pal = pal
)

plotPDF(p1,p2, name = "Plot-UMAP_harmony-RNA-Integration.pdf", ArchRProj = proj.filter_Fb, addDOC = FALSE, width = 5, height = 5)

p1 <- plotEmbedding(
    proj.filter_Fb, 
    colorBy = "cellColData", 
    name = "predictedGroup_Un_2", 
    embedding = "UMAP_LSI",
    size = 0.5,
    pal = pal
)

p2 <- plotEmbedding(
    proj.filter_Fb, 
    colorBy = "cellColData", 
    name = "predictedGroup_Co",
    embedding = "UMAP_LSI",
    size = 0.5,
    pal = pal
)

plotPDF(p1,p2, name = "Plot-UMAP_LSI-RNA-Integration.pdf", ArchRProj = proj.filter_Fb, addDOC = FALSE, width = 5, height = 5)


projHeme3 <- addGeneIntegrationMatrix(
    ArchRProj = proj.filter_Fb, 
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

projHeme3 <- saveArchRProject(ArchRProj = projHeme3, outputDirectory = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-proj.filter_Fb", load = TRUE)

getAvailableMatrices(projHeme3)
#[1] "GeneIntegrationMatrix" "GeneScoreMatrix"       "TileMatrix" 
projHeme3 <- addImputeWeights(projHeme3)


p1 <- plotEmbedding(
    ArchRProj = projHeme3, 
    colorBy = "GeneIntegrationMatrix", 
    name = markerGenes, 
    continuousSet = "horizonExtra",
    embedding = "UMAP_harmony",
    imputeWeights = getImputeWeights(projHeme3)
)
p2 <- plotEmbedding(
    ArchRProj = projHeme3, 
    colorBy = "GeneScoreMatrix", 
    continuousSet = "horizonExtra",
    name = markerGenes, 
    embedding = "UMAP_harmony",
    imputeWeights = getImputeWeights(projHeme3)
)
plotPDF(plotList = p1, 
    name = "Plot-UMAP-Marker-Genes-RNA-W-Imputation.pdf", 
    ArchRProj = projHeme3, 
    addDOC = FALSE, width = 5, height = 5)

cM <- confusionMatrix(projHeme3$Clusters_2, projHeme3$predictedGroup)
labelOld <- rownames(cM)
labelOld

labelNew <- colnames(cM)[apply(cM, 1, which.max)]
labelNew

projHeme3$Clusters2 <- mapLabels(projHeme3$Clusters_2, newLabels = labelNew, oldLabels = labelOld)

p1 <- plotEmbedding(projHeme3, colorBy = "cellColData", name = "Clusters2", embedding = "UMAP_harmony",size = 0.5)
plotPDF(p1, name = "Plot-UMAP_harmony-Remap-Clusters.pdf", ArchRProj = projHeme3, addDOC = FALSE, width = 5, height = 5)
projHeme3 <- saveArchRProject(ArchRProj = projHeme3, outputDirectory = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme3", load = TRUE)
head(projHeme3@cellColData)


setwd("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy")

projHeme3 <- loadArchRProject("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme3")
head(projHeme3@cellColData)
table(projHeme3$Sample)
getAvailableMatrices(projHeme3)
getOutputDirectory(ArchRProj = projHeme3)

library(dplyr)
projHeme3$group<- recode(projHeme3$Sample,
                        "HAS0416" = "HA",
                        "HAES1" = "HA",
                        "OAS_0510" = "OA",
                        "OAS_ZXG" = "OA",
                        .default = NA_character_)
projHeme3@cellColData %>% head()

projHeme3$Cgroup <- paste0(projHeme3$group, "_", projHeme3$Clusters2)
table(projHeme3$Cgroup)

projHeme3 <- saveArchRProject(ArchRProj = projHeme3, outputDirectory = "Save-ProjHeme3_new", load = TRUE)

library(BSgenome.Hsapiens.UCSC.hg38)
projHeme4 <- addGroupCoverages(ArchRProj = projHeme3, groupBy = "Cgroup")

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
markerList <- getMarkers(markerPeaks, cutOff = "FDR <= 0.01 & Log2FC >= 1")
markerList

heatmapPeaks <- plotMarkerHeatmap(
  seMarker = markerPeaks, 
  cutOff = "FDR <= 0.1 & Log2FC >= 0.5",
  transpose = TRUE
)

plotPDF(heatmapPeaks, name = "Peak-Marker-Heatmap", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)

markerPeaks1 <- getMarkerFeatures(
    ArchRProj = projHeme5, 
    useMatrix = "PeakMatrix", 
    groupBy = "group",
  bias = c("TSSEnrichment", "log10(nFrags)"),
  testMethod = "wilcoxon"
)

markerPeaks1

markerList1 <- getMarkers(markerPeaks1, cutOff = "FDR <= 0.1 & Log2FC >= 0.5")

markerList1$HA

heatmapPeaks1 <- plotMarkerHeatmap(
  seMarker = markerPeaks1, 
  cutOff = "FDR <= 0.05 & Log2FC >= 1",
  transpose = TRUE
)
plotPDF(heatmapPeaks1, name = "Peak-Marker-Heatmap_25.3.4_group", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)




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
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme4_new/pSet_peaks.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

pSet_HA <- subset(pSet, group == "HA")
pSet_OA <- subset(pSet, group == "OA")

write.table(pSet_HA, 
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme5_new/pSet_HA_peaks.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

write.table(pSet_OA, 
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme5_new/pSet_OA_peaks.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

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
pma

pv <- plotMarkers(seMarker = markerTest, name = "HA", cutOff = "FDR <= 0.5 & abs(Log2FC) >= 1", plotAs = "Volcano")
pv
plotPDF(pma, pv, name = "HA-vs-OA-Markers-MA-Volcano", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)


projHeme5 <- addMotifAnnotations(ArchRProj = projHeme5, motifSet = "cisbp", name = "Motif")

matches <- getMatches(ArchRProj = projHeme5, name = "Motif")
rownames(matches) <- paste(seqnames(matches), start(matches), end(matches), sep = "_")
matches <- matches[pSet$name]

gr <- GRanges(seqnames = c("chr11"), ranges = IRanges(start = c(46381185), end = c(46381685)))
queryHits <- queryHits(findOverlaps(query = pSet, subject = gr, type = "within"))
colnames(matches)[which(assay(matches[queryHits,]))]


motifsUp <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = projHeme5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.8 & Log2FC >= 0.5"
  )
df <- data.frame(TF = rownames(motifsUp), mlog10Padj = assay(motifsUp)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)
write.table(df, 
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme4_new/motifsUp.tsv", 
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
png(file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme4_new/Rank_Sorted_TFs_Enriched_ggUp.png", width = 3700, height = 4200, res = 300)
ggUp
dev.off()

motifsDo <- peakAnnoEnrichment(
    seMarker = markerTest,
    ArchRProj = projHeme5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.8 & Log2FC <= -0.5"
  )

df <- data.frame(TF = rownames(motifsDo), mlog10Padj = assay(motifsDo)[,1])
df <- df[order(df$mlog10Padj, decreasing = TRUE),]
df$rank <- seq_len(nrow(df))
head(df)
write.table(df, 
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme4_new/motifsDo.tsv", 
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
png(file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme4_new/Rank_Sorted_TFs_Enriched_ggDo.png", width = 4000, height = 4200, res = 300)
ggDo
dev.off()

plotPDF(ggUp, ggDo, name = "HA-vs-OA-Markers-Motifs-Enriched", width = 10, height = 10, ArchRProj = projHeme5, addDOC = FALSE)

enrichMotifs <- peakAnnoEnrichment(
    seMarker = markerPeaks,
    ArchRProj = projHeme5,
    peakAnnotation = "Motif",
    cutOff = "FDR <= 0.8 & Log2FC >= 0.5"
  )
heatmapEM <- plotEnrichHeatmap(enrichMotifs,transpose = TRUE)
plotPDF(heatmapEM, name = "Motifs-Enriched-Marker-HeatmapFDR <= 0.8 & Log2FC >= 0.5", width = 8, height = 6, ArchRProj = projHeme5, addDOC = FALSE)


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


projHeme5 <- addBgdPeaks(projHeme5)

projHeme5 <- addDeviationsMatrix(
  ArchRProj = projHeme5, 
  peakAnnotation = "Motif",
  force = TRUE
)

plotVarDev <- getVarDeviations(projHeme5, name = "MotifMatrix", plot = TRUE)

plotPDF(plotVarDev, name = "Variable-Motif-Deviation-Scores", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)
## Plotting Ggplot!

motifs <- c("PRRX1", "CEBPD", "KLF6", "MITF", "SNAI2","STAT1")
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
plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)

saveArchRProject(projHeme5,outputDirectory = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme5_group",overwrite = TRUE)

p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "group", 
    geneSymbol = c("PRRX1", "CEBPD", "KLF6", "MITF", "SNAI2","STAT1"),
    upstream = 2000,
    downstream = 2000
)
plotPDF(p, name = "Plot-Tracks-With-up_motif", width = 10, height = 8, ArchRProj = projHeme5, addDOC = FALSE)

p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "group", 
    geneSymbol = c( "CEBPD","STAT1"),
    upstream = 200000,
    downstream = 200000
)
plotPDF(p, name = "Plot-Tracks-With-up_motif_5.30", width = 20, height = 8, ArchRProj = projHeme5, addDOC = FALSE)

p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "group", 
    geneSymbol = c( "CEBPD","STAT1"),
    upstream = 100000,
    downstream = 100000
)
plotPDF(p, name = "Plot-Tracks-With-up_motif_5.30_1", width = 20, height = 8, ArchRProj = projHeme5, addDOC = FALSE)
motifs <- c("CEBPD","STAT1")
markerMotifs <- getFeatures(projHeme5, select = paste(motifs, collapse="|"), useMatrix = "MotifMatrix")
markerMotifs
#[1] "z:STAT1_773"          "z:CEBPD_152"          "deviations:STAT1_773"
#[4] "deviations:CEBPD_152"
markerMotifs <- grep("z:", markerMotifs, value = TRUE)
#markerMotifs <- markerMotifs[markerMotifs %ni% "z:SREBF1_22"]
markerMotifs

p <- plotGroups(ArchRProj = projHeme5, 
  groupBy = "group", 
  colorBy = "MotifMatrix", 
  name = markerMotifs,
  imputeWeights = getImputeWeights(projHeme5)
)
plotPDF(p, name = "Plot-Groups-Deviations-w-Imputation_25.5.30_sy", width = 5, height = 5, ArchRProj = projHeme5, addDOC = FALSE)


projHeme5 <- loadArchRProject("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme5_new/")

motifPositions <- getPositions(projHeme5)

motifPositions

motifs <- c("PRRX1", "CEBPD", "KLF6", "MITF", "SNAI2","STAT1")
markerMotifs <- unlist(lapply(motifs, function(x) grep(x, names(motifPositions), value = TRUE)))
#markerMotifs <- markerMotifs[markerMotifs %ni% "SREBF1_22"]
markerMotifs
## [1] "PRRX1_440" "CEBPD_152" "KLF6_165"  "MITF_91"   "SNAI2_161" "STAT1_773"
library(BSgenome.Hsapiens.UCSC.hg38)
if(is.null(projHeme5@projectMetadata$GroupCoverages$Clusters2)){
  projHeme5 <- addGroupCoverages(ArchRProj = projHeme5, groupBy = "group")
}
seFoot <- getFootprints(
  ArchRProj = projHeme5, 
  positions = motifPositions[markerMotifs], 
  groupBy = "group"
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projHeme5, 
  normMethod = "Subtract",
  plotName = "Footprints-Subtract-Bias",
  addDOC = FALSE,
  smoothWindow = 15
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projHeme5, 
  normMethod = "Divide",
  plotName = "Footprints-Divide-Bias",
  addDOC = FALSE,
  smoothWindow = 10
)

plotFootprints(
  seFoot = seFoot,
  ArchRProj = projHeme5, 
  normMethod = "None",
  plotName = "Footprints-No-Normalization",
  addDOC = FALSE,
  smoothWindow = 15
)

p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "group", 
    geneSymbol = c( "KLF6"),
    upstream = 100000,
    downstream = 100000
)
plotPDF(p, name = "Plot-Tracks-With-up_motif_6.2_1", width = 20, height = 8, ArchRProj = projHeme5, addDOC = FALSE)

#markerMotifs
# [1] "FOSL2_105"   "FOSL2_105"   "FOSB_121"    "FOS_137"     "FOSL1_142"  
# [6] "FOSL1_142"   "SMARCC1_651" "JUND_124"    "BACH1_130"   "JUND_124"   
#[11] "JUNB_139"    "JUN_143"     "JUNB_139"    "RUNX1_733"  
pwm <- getPeakAnnotation(projHeme5, "Motif")$motifs[["JUNB_139"]]
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
png("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme5_new/Plots/logo_plot_JUNB_139.png", width = 1200, height = 800, res = 150)
ggseqlogo(ppm, method = "bits")
dev.off()


pwm <- getPeakAnnotation(projHeme5, "Motif")$motifs[["CEBPD_152"]]
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
png("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme5_group/Plots/logo_plot_CEBPD_152.png", width = 1200, height = 800, res = 150)
ggseqlogo(ppm, method = "bits")
dev.off()



projHeme5 <- saveArchRProject(ArchRProj = projHeme5, outputDirectory = "Save-ProjHeme5_new", load = TRUE)

coverageFiles <- getGroupBW(ArchRProj = projHeme5, groupBy = "Cgroup")


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
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme5_new/CoAccessibility_0.5.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)

cB=metadata(cA)[[1]]

write.table(cB, 
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme5_new/CoAccessibility_0.5_peak.tsv", 
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
            file = "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme5_new/CoAccessibility_0.5_loops.tsv", 
            sep = "\t", 
            row.names = FALSE, 
            col.names = TRUE, 
            quote = FALSE)
p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "group", 
    geneSymbol = c("FOSL2", "FOS", "FOSL1", "SMARCC1", "JUND","BACH1","JUN","JUNB","RUNX1"),
    upstream = 100000,
    downstream = 100000
)
plotPDF(p, name = "Plot-Tracks-With-up_motif_6.2", width = 20, height = 8, ArchRProj = projHeme5, addDOC = FALSE)


projHeme5 <- addPeak2GeneLinks(
    ArchRProj = projHeme5,
    reducedDims = "IterativeLSI"
)

p

p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "group", 
    geneSymbol =  c("JUN","FOS","JUNB","SMARCC1"), 
    upstream = 60000,
    downstream = 60000,
    loops = getPeak2GeneLinks(projHeme5, corCutOff = 0.3,resolution = 100,returnLoops = TRUE)
)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.pdf", 
    ArchRProj = projHeme5, 
    addDOC = FALSE, width = 15, height = 8)


###########
projHeme5 <- loadArchRProject("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme5_new")
projHeme5 <- addPeak2GeneLinks(
    ArchRProj = projHeme5,
    reducedDims = "IterativeLSI"
)
p <- plotPeak2GeneHeatmap(ArchRProj = projHeme5, groupBy = "Clusters2")
setwd("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme5_new")
png(file = "plotPeak2GeneHeatmap.png", width = 3000, height = 3000, res = 300)
p
dev.off()
####
projHeme5 <- loadArchRProject("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme5_new")

projHeme5 <- addCoAccessibility(
    ArchRProj = projHeme5,
    reducedDims = "IterativeLSI"
)
p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "group", 
    geneSymbol =  c("ANGPTL4"), 
    upstream = 100000,
    downstream = 100000,
    loops = getCoAccessibility(projHeme5)
)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-CoAccessibility.ANGPTL4.pdf", 
    ArchRProj = projHeme5, 
    addDOC = FALSE, width = 15, height = 8)


projHeme5 <- addPeak2GeneLinks(
    ArchRProj = projHeme5,
    reducedDims = "IterativeLSI"
)
p <- plotBrowserTrack(
    ArchRProj = projHeme5, 
    groupBy = "group", 
    geneSymbol =  c("ANGPTL4"), 
    upstream = 100000,
    downstream = 100000,
    loops = getPeak2GeneLinks(projHeme5)
)
plotPDF(plotList = p, 
    name = "Plot-Tracks-Marker-Genes-with-Peak2GeneLinks.ANGPTL4.pdf", 
    ArchRProj = projHeme5, 
    addDOC = FALSE, width = 15, height = 8)

coverageFiles <- getGroupBW(ArchRProj = projHeme5, groupBy = "group",normMethod="None")
coverageFiles1 <- getGroupBW(ArchRProj = projHeme5, groupBy = "Cgroup",normMethod="None")
coverageFiles2 <- getGroupBW(ArchRProj = projHeme5, groupBy = "Sample")

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
seuratObj[["umap"]] <- CreateDimReducObject(embeddings = as.matrix(projHeme5@embeddings$UMAP_harmony$df), key = "UMAP_", assay = DefaultAssay(seuratObj))


png(file = "/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Sy/new/marker_gene/second_celltype_umap_new.png", width = 3000, height = 2500, res = 300)
DimPlot(
  seuratObj,
  reduction = "umap", 
  group.by = "predictedGroup", 
  label = TRUE,
  label.size = 6,  # 设置标签的字体大小
   cols = c('#DBC9B3','#EED0E0','#EBAEA9','#CB95BB','#EBCC96','#AED0DF','#CBE5DE'
))
dev.off()

saveRDS(seuratObj,file="/mnt/8w/data4/yiyuan/Bone/scRNA/补充figure/Sy/new/marker_gene/scATAC_Sy_resubtype.final.rds")



