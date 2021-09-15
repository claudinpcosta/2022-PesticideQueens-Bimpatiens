######################## Pesticide Queens - GOinputs ###############################

#Using table from Colgan et al. 2019 for GO analyses
#https://github.com/wurmlab/Bter_neonicotinoid_exposure_experiment/blob/master/07_gene_ontology_enrichment_analysis/input/dmel_vs_bter_biomart.input_for_converter.output.txt
BterGO <- read.table("Data/dmel_vs_bter_biomart.input_for_converter.output.txt")
colnames(BterGO) <- c("Bter", "GOterms") 
head(BterGO)
dim(BterGO)

#You can get gene orthologs list between Bombus impatiens and Bombus terrestris in the EnsemblMetazoa website by BioMart:
#https://metazoa.ensembl.org/biomart/martview/1ecde54e5164b2539a7e20381e990ada
#I downloaded the file and save as orthologs.csv in Data folder. 
orthologs <- read.csv("Data/BimpBterOrt.csv")
head(orthologs)
dim(orthologs)

#Merge lists above
BimpBterGO <- merge(BterGO,orthologs, by = "Bter", all = TRUE)
head(BimpBterGO)
dim(BimpBterGO)
BimpBterGO <- BimpBterGO[!is.na(BimpBterGO$GOterms), ]
head(BimpBterGO)
dim(BimpBterGO)
BimpBterGO <- BimpBterGO[!is.na(BimpBterGO$Bimp), ]
head(BimpBterGO)
dim(BimpBterGO)
BimpGO <- subset (BimpBterGO, select = -Bter)
head(BimpGO)
dim(BimpGO)
BimpGO <- BimpGO[c("Bimp", "GOterms")]
head(BimpGO)
dim(BimpGO)
BimpGO <- BimpGO[!duplicated(BimpGO$Bimp), ]
head(BimpGO)
dim(BimpGO)

#read the file with BIMP names and Locus (in this dataset we used Locus names)
BimpLocBimp <- read.csv("Data/BimpLocBimp.csv", row.names = 1)
head(BimpLocBimp)
dim(BimpLocBimp)

BimpLocusGO <- merge(BimpGO,BimpLocBimp, by = "Bimp", all = TRUE)
head(BimpLocusGO)
dim(BimpLocusGO)
BimpLocusGO <- subset (BimpLocusGO, select=c(5,2))
head(BimpLocusGO)
dim(BimpLocusGO)
BimpLocusGO <- na.omit(BimpLocusGO)
head(BimpLocusGO)
dim(BimpLocusGO)
BimpLocusGO <- BimpLocusGO[!duplicated(BimpLocusGO$Locus), ]
head(BimpLocusGO)
dim(BimpLocusGO)


write.table(as.data.frame(BimpLocusGO), col.names = F, row.names = F,sep="\t",
          file="~/Dropbox/2021-PesticideQueens-Bimpatiens/Data/BimpGO.txt")

#remove quotes in txt (open the file and find/replace function)









