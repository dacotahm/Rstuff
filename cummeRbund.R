library(cummeRbund)

cuffData <- readCufflinks(dir = getwd(), genome = 'osmiaFullStageDL.fasta', gtfFile = 'merged.gtf')

# Load annotation
#annot <- read.table("osmiaFullStageDL.gtf",sep="\t",header=T,na.string="-")
#addFeatures(cuffData,annot,level="genes")


cuffDiffIDs <- getSig(cuffData, level = 'genes', alpha = 0.05)
cuffDiffGenes<- getGenes(cuffData, cuffDiffIDs)
CuffFeatureNames <- featureNames(cuffDiffGenes)





# QC plots
Cuffdisp <- dispersionPlot(genes(cuffData))
jpeg(filename = 'DispersionPlot.jpg', width = 1500, height = 1000, units = 'px', pointsize = 14, bg = 'white')
Cuffdisp
dev.off()

CuffGenes <- fpkmSCVPlot(genes(cuffData))
jpeg(filename = 'fpkm_SCV.jpg', width = 1500, height = 1000, units = 'px', pointsize = 14, bg = 'white')
CuffGenes
dev.off()

CuffDensityPlot <- csDensity(genes(cuffData))
jpeg(filename = 'DensityPlot.jpg', width = 1500, height = 1000, units = 'px', pointsize = 14, bg = 'white')
CuffDensityPlot
dev.off()

CuffBoxPlot <- csBoxplot(genes(cuffData))
jpeg(filename = 'BoxPlot.jpg', width = 1500, height = 1000, units = 'px', pointsize = 14, bg = 'white')
CuffBoxPlot
dev.off()

CuffBoxPlotReps <- csBoxplot(genes(cuffData), replicates = TRUE)
jpeg(filename = 'BoxPlot_reps.jpg', width = 1500, height = 1000, units = 'px', pointsize = 14, bg = 'white')
CuffBoxPlotReps
dev.off()

#CuffDendroReps <- csDendro(genes(cuffData), replicates = TRUE)
jpeg(filename = 'Dendro_reps.jpg', width = 1500, height = 1000, units = 'px', pointsize = 14, bg = 'white')
csDendro(genes(cuffData), replicates = TRUE)
dev.off()

CuffDendro <- csDendro(genes(cuffData))
jpeg(filename = 'Dendro.jpg', width = 1500, height = 1000, units = 'px', pointsize = 14, bg = 'white')
csDendro(genes(cuffData))
dev.off()

library(dendextend)
library(ggdendro)
library(ggplot2)
plot(CuffDendro)
CuffDendro %>% set("branches_lwd", 2) %>% set("branches_k_color") %>% plot(main = "Clustering by gene expression")
CuffDendro %>% rect.dendrogram(k=3, border = 8, lty = 2, lwd = 2)


#CuffVolcano
CuffCano <- csVolcanoMatrix(genes(cuffData))
jpeg(filename = 'VolcanoMatrix.jpg', width = 1500, height = 1000, units = 'px', pointsize = 14, bg = 'white')
CuffCano
dev.off()

#K-means clustering
cuffDiffIDs <- getSig(cuffData, level = 'genes', alpha = 0.05)
sigGenes <- getGenes(cuffData, cuffDiffIDs)

ic9 <- csCluster(sigGenes, k=9)
icp9 <- csClusterPlot(ic9)
jpeg(filename = 'k-9Cluster.jpg', width = 1500, height = 1000, units = 'px', pointsize = 14, bg = 'white')
icp9
dev.off()

ic16 <- csCluster(sigGenes, k=16)
icp16 <- csClusterPlot(ic16)
jpeg(filename = 'k-16Cluster.jpg', width = 1500, height = 1000, units = 'px', pointsize = 14, bg = 'white')
icp16
dev.off()

ic25 <- csCluster(sigGenes, k=25)
icp25 <- csClusterPlot(ic25)
jpeg(filename = 'k-25Cluster.jpg', width = 1500, height = 1000, units = 'px', pointsize = 14, bg = 'white')
icp25
dev.off()

#PCA, MDS
genes.PCA <- PCAplot(genes(cuffData))
jpeg(filename = 'PCA.jpg', width = 1500, height = 1000, units = 'px', pointsize = 14, bg = 'white')
genes.PCA
dev.off()

genes.PCA.reps <- PCAplot(genes(cuffData), replicates=T)
jpeg(filename = 'PCA_reps.jpg', width = 1500, height = 1000, units = 'px', pointsize = 14, bg = 'white')
genes.PCA.reps
dev.off()

genes.MDS <- MDSplot(genes(cuffData))
jpeg(filename = 'MDS.jpg', width = 1500, height = 1000, units = 'px', pointsize = 14, bg = 'white')
genes.MDS
dev.off()

genes.MDS.reps <- MDSplot(genes(cuffData), replicates=T)

genes.MDS.reps + geom_vline(aes(xintercept = 0), size = .5, color = 'red') + 
  geom_hline(aes(yintercept = 0), size = .5, color = 'red') + 
  scale_x_continuous(limits = c(-0.2, 0.2)) +
  scale_y_continuous(limits = c(-.15, .25))
  
jpeg(filename = 'MDS_reps.jpg', width = 1500, height = 1000, units = 'px', pointsize = 14, bg = 'white')
genes.MDS.reps
dev.off()

