#Ch
library("ChIPseeker")
library("clusterProfiler")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peak <- readPeakFile("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme5_group/markerList_HA.bed")

peakAnno <- annotatePeak(
    peak,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db")

write.table(
    as.data.frame(peakAnno),
    "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Ch/Save-ProjHeme5_group/markerList_HA_anno.tsv",
    sep="\t",
    row.names = F,
    quote = F)

#Sy
library("ChIPseeker")
library("clusterProfiler")
library("org.Hs.eg.db")
library("TxDb.Hsapiens.UCSC.hg38.knownGene")
txdb <- TxDb.Hsapiens.UCSC.hg38.knownGene
peak <- readPeakFile("/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme5_group/markerList_HA.bed")

peakAnno <- annotatePeak(
    peak,
    tssRegion = c(-3000, 3000),
    TxDb = txdb,
    annoDb = "org.Hs.eg.db")

write.table(
    as.data.frame(peakAnno),
    "/mnt/8w/data4/yiyuan/Bone/scATAC/result/Sy/Save-ProjHeme5_group/markerList_HA_anno.tsv",
    sep="\t",
    row.names = F,
    quote = F)