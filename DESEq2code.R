######################## Pesticide Queens - DESeq2 ###############################

#### Loading packages ####
library("DESeq2")
library("vsn")
library("ggplot2")
library("tibble")
library("plotly")
library("dplyr")
library("ggpubr")
library("RColorBrewer")
library("pheatmap")
theme_set(theme_pubr())


#### import the data ####

# read counts
readcounts <- read.csv("Data/readcount.csv", header = T, row.names = 1)
dim(readcounts)
head(readcounts)

# coldata
coldata <- read.csv("Data/coldata.csv", header = T, row.names = 1)
dim(coldata)
head(coldata)

# verify if the coldata's rownames is equal to readcounts' colnames
all(rownames(coldata)%in%colnames(readcounts))
all(rownames(coldata)==colnames(readcounts))

#### Construct a DESeqDataSet ####

ddsData <- DESeqDataSetFromMatrix(countData = readcounts,
                                  colData = coldata,
                                  design= ~ Pollen + ~ Pesticide + ~ Colony + ~ Pollen:Pesticide)

ddsData

keep <- rowSums(counts(ddsData)) >= 10
ddsData <- ddsData[keep,]
ddsData



#### PCA ####

# transform expression levels using the regularized log transformation
rld <- rlog(ddsData, blind = FALSE)

# transform expression levels using the variance stabilizing transformation (VST)
vsd <- vst(ddsData, blind = FALSE)

# do these meanSD plots on data before DESeq
notAllZero <- (rowSums(counts(ddsData))>0)
meanSdPlot(log2(counts(ddsData,normalized=FALSE)[notAllZero,]+1))
meanSdPlot(assay(rld[notAllZero,]))
meanSdPlot(assay(vsd[notAllZero,]))

## PCA - for rld-transformed data

# Simplest PCA, not easy to interpret
plotPCA(rld, intgroup=c("Pollen","Pesticide"))

# Gettin' fancier with the PCA
pca.dataRLD <- plotPCA(rld, intgroup=c("Pollen","Pesticide"), returnData=TRUE)
percentVar <- round(100 * attr(pca.dataRLD, "percentVar"))
ggplot(pca.dataRLD, aes(PC1, PC2, color=Pollen, shape=Pesticide)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed() +ggtitle("RLD-transformed data")

#format the graphic 
p <- ggplot(pca.dataRLD, aes(PC1, PC2, color=Pollen, shape=Pesticide)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()
p + theme_bw()
p + theme_classic()
p + theme_light() +theme(legend.text=element_text(size=12)) + scale_color_manual(values=c("orchid", "hotpink1", "deeppink3"))


## PCA - for vsd-transformed data

# Simplest PCA for VSD
plotPCA(vsd, intgroup=c("Pollen","Pesticide"))

# Gettin' fancier with the PCA
pca.dataVSD <- plotPCA(vsd, intgroup=c("Pollen","Pesticide"), returnData=TRUE)
percentVar <- round(100 * attr(pca.dataVSD, "percentVar"))
ggplot(pca.dataVSD, aes(PC1, PC2, color=Pollen, shape=Pesticide)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()+ggtitle("VSD-transformed data")

#format the graphic 
q<-ggplot(pca.dataVSD, aes(PC1, PC2, color=Pollen, shape=Pesticide)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()+ggtitle("VSD-transformed data")
q + theme_bw()
q + theme_light()
q + theme_light() +theme(legend.text=element_text(size=12)) + scale_color_manual(values=c("orchid", "hotpink1", "deeppink3"))


#### Pheatmap ####
top50 <- order(rowMeans(counts(ddsData,normalized=FALSE)),decreasing=TRUE)[1:50]
top100 <- order(rowMeans(counts(ddsData,normalized=FALSE)),decreasing=TRUE)[1:100]
top500 <- order(rowMeans(counts(ddsData,normalized=FALSE)),decreasing=TRUE)[1:500]

nt <- normTransform(ddsData)
log2.norm.counts <- assay(nt)[top50,]

df <- as.data.frame(colData(ddsData)[,c("Pollen","Pesticide")])
df$Pollen = factor(df$Pollen, levels = c("Cistus", "Erica", "Mixed"))
PollenCol <- c("thistle3", "hotpink1", "deeppink3")
names(PollenCol) <- levels(df$Pollen)
df$Pesticide = factor(df$Pesticide, levels = c("untreated", "IMD-A", "IMD-B"))
PesticideCol <- c("lavender", "lightblue1", "steelblue1")
names(PesticideCol) <- levels(df$Pesticide)
AnnColour <- list(
  Pollen = PollenCol,
  Pesticide = PesticideCol)

distmat <- dist(t(top100))
makeColorRampPalette <- function(colors, cutoff.fraction, num.colors.in.palette)
{
  stopifnot(length(colors) == 4)
  ramp1 <- colorRampPalette(colors[1:2])(num.colors.in.palette * cutoff.fraction)
  ramp2 <- colorRampPalette(colors[3:4])(num.colors.in.palette * (1 - cutoff.fraction))
  return(c(ramp1, ramp2))
}
cutoff.distance <- 3  
cols <- makeColorRampPalette(c("green", "purple",    # distances 0 to 3 colored from white to red
                               "green", "purple"), # distances 3 to max(distmat) colored from green to black
                             cutoff.distance / max(distmat),
                             30)

#changing show_rownames from FALSE to TRUE made gene names visible. More clutter, but easy to ID blue row.
pheatmap(log2.norm.counts, cluster_rows = FALSE, show_rownames = T, cluster_cols = T, annotation_col = df)
#the following heatmaps are less informative and interesting when using dds instead
#of rld data. Reconsider input, and if transformed first
pheatmap(assay(rld)[top50,], cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = T, annotation_col = df)
pheatmap(assay(rld)[top500,], cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = T, annotation_col = df)
pheatmap(assay(rld)[top100,], cluster_rows = FALSE, show_rownames = FALSE, cluster_cols = T, annotation_col = df)
#without indiv bee IDs listed
pheatmap(assay(rld)[top50,], color = cols, cluster_rows = FALSE, show_rownames = FALSE, show_colnames = FALSE, cluster_cols = F, annotation_col = df, annotation_colors = AnnColour)
pheatmap(assay(rld)[top500,], color = cols, cluster_rows = FALSE, show_rownames = FALSE, show_colnames = FALSE, cluster_cols = T, annotation_col = df, annotation_colors = AnnColour)



