######################## Pesticide Queens - DESeq2 ###############################

#### Loading basic packages ####
library("DESeq2")
library("DEGreport")
library("ggplot2")
library("ggrepel")
library("tibble")
library("plotly")
library("dplyr")
library("ggpubr")
library("RColorBrewer")
library("vsn")
library("zFPKM")
library("gplots")
library("pheatmap")
theme_set(theme_pubr())


#### import the data ####

#read counts
readcounts <- read.csv("RNAseq/Data/readcount.csv", header = T, row.names = 1)
head(readcounts)
dim(readcounts)

#coldata
coldata <- read.csv("RNAseq/Data/coldata.csv", header = T, row.names = 1)
head(coldata)
dim(coldata)

#verify if the coldata's rownames is equal to readcounts' colnames
all(rownames(coldata)%in%colnames(readcounts))
all(rownames(coldata)==colnames(readcounts))


#### Gene lists by LRT ####

##construct a DESeqDataSet 
ddsData <- DESeqDataSetFromMatrix(countData = readcounts,
                                  colData = coldata,
                                  design= ~ Pollen + ~ Pesticide + ~ Colony + ~ Pollen:Pesticide)
ddsData
ddsData$Pesticide <- relevel(ddsData$Pesticide, ref = "untreated")
ddsData$Pollen <- relevel(ddsData$Pollen, ref = "Mixed")
keep <- rowSums(counts(ddsData)) >= 10
ddsData <- ddsData[keep,]
ddsData


##Gene list for Colony 
ddsCol <- DESeq(ddsData, test = "LRT", reduced = ~ Colony)

resCol <- results(ddsCol, alpha=0.05); resCol
summary(resCol)
sum(resCol$padj < 0.05, na.rm=TRUE)
resCol <- resCol[order(resCol$pvalue),]; resCol

GLCol <- subset(resCol, padj < 0.05)
summary(GLCol)
dim(GLCol)


##Gene list for Pollen (diet)
ddsPollen <- DESeq(ddsData, test = "LRT", reduced = ~ Pollen)

resPollen <- results(ddsPollen, alpha=0.05); resPollen
summary(resPollen)
sum(resPollen$padj < 0.05, na.rm=TRUE)
resPollen <- resPollen[order(resPollen$pvalue),]; resPollen

GLPollen <- subset(resPollen, padj < 0.05)
summary(GLPollen)
dim(GLPollen)


##Gene list for Pesticide (IMD)
ddsPes <- DESeq(ddsData, test = "LRT", reduced = ~ Pesticide)

resPes <- results(ddsPes, alpha=0.05); resPes
summary(resPes)
sum(resPes$padj < 0.05, na.rm=TRUE)
resPes <- resPes[order(resPes$pvalue),]; resPes

GLPes <- subset(resPes, padj < 0.05)
summary(GLPes)
dim(GLPes)


##Gene list for interaction
ddsInteraction <- DESeq(ddsData, test = "LRT", reduced = ~ Pollen:Pesticide)

resInteraction <- results(ddsInteraction, alpha=0.05); resInteraction
summary(resInteraction)
sum(resInteraction$padj < 0.05, na.rm=TRUE)
resInteraction <- resInteraction[order(resInteraction$pvalue),]; resInteraction

GLInteraction <- subset(resInteraction, padj < 0.05)
summary(GLInteraction)
dim(GLInteraction)

#### Downloading the interested LRT Gene lists with information about genes for manuscript ####

#For this I need get the information from Biomart and NCBI, as accession number, gene name and others.
#download packages that I will use
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("Biostrings")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("biomaRt")
# 
# # Now install biomartr from CRAN
# install.packages("biomartr", dependencies = TRUE)

library(biomaRt) # Access gene databases through R
library(biomartr) # Extension of biomaRt that is a little more user friendly - I will not use this for the code below 
library(tidyverse)
library(dplyr)

listEnsemblGenomes() #what are biomarts available? 
ensembl_metazoa<-useEnsemblGenomes(biomart = "metazoa_mart")
searchDatasets(ensembl_metazoa, patter="Bombus")
#This will get the dataset I need from Bombus impatiens
ensembl_bimp<-useEnsemblGenomes(biomart = "metazoa_mart",
                                dataset="bimpatiens_eg_gene")
#Lists all of the attributes in Bombus dataset
listAttributes(ensembl_bimp) 
#Search if you would like specific term in your attributes
searchAttributes(mart = ensembl_bimp, pattern = "Drosophila melanogaster")
# Get information that I will use, ("ensembl_gene_id" = Bimp; "refseq_peptide" = accession number, the same information from NCBI wilh the name "Protein.product")
BimpsProteinIds<-getBM(attributes = c('ensembl_gene_id',
                                      'refseq_peptide'), 
                       mart = ensembl_bimp)
head(BimpsProteinIds)
dim(BimpsProteinIds)
#change the column names, it will be helpful when we merge with NCBI information.
names(BimpsProteinIds)[names(BimpsProteinIds) == "ensembl_gene_id"] <- "Bimp"
names(BimpsProteinIds)[names(BimpsProteinIds) == "refseq_peptide"] <- "Protein.product"
head(BimpsProteinIds)
dim(BimpsProteinIds)

#Once you have accession numbers (refseq_peptide from ensembl) can get NCBI information here:
#https://www.ncbi.nlm.nih.gov/genome/browse/#!/proteins/3415/468201%7CBombus%20impatiens/
#I downloaded the file and save as ProteinTableBimp.csv in Data folder. 

BimpLocs <- read.csv("RNAseq/Data/ProteinTableBimp.csv")
head(BimpLocs,30)
dim(BimpLocs)

#merge information from Biomart and NCBI
BimpLocBimp <- merge(BimpsProteinIds,BimpLocs, by = "Protein.product")
head(BimpLocBimp)
dim(BimpLocBimp)
#I am interested in only some information
BimpLocBimp <- subset(BimpLocBimp, select=c(1:2,8:9,12))
head(BimpLocBimp)
dim(BimpLocBimp)
#change the name for Accession.Number will be helpful for results
names(BimpLocBimp)[names(BimpLocBimp) == "Protein.product"] <- "Accession.Number"
head(BimpLocBimp)
dim(BimpLocBimp)

write.csv(as.data.frame(BimpLocBimp),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Data/BimpLocBimp.csv")

#Now create a folder to save results for this analysis
dir.create("RNAseq/Results") #create a folder all results
dir.create("RNAseq/Results/LRT") #sub folder to save this analysis 


#Now get the table for each condition:

#colony DEGs
ColDEG <- as.data.frame(GLCol)
head(ColDEG)
dim(ColDEG)
ColDEG <- subset(ColDEG, select=c(2, 6))
head(ColDEG)
dim(ColDEG)
ColDEG <- tibble::rownames_to_column(ColDEG, "Locus")
head(ColDEG)
dim(ColDEG)
ColDEG$direction <- ifelse(ColDEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(ColDEG)
dim(ColDEG)
ColDEG <- merge(ColDEG,BimpLocBimp, by = "Locus", all = TRUE)
head(ColDEG)
dim(ColDEG)
ColDEG <- ColDEG[!is.na(ColDEG$direction), ]
dim(ColDEG)
head(ColDEG, 30)
ColDEG <- ColDEG[!duplicated(ColDEG$Locus), ]
head(ColDEG, 30)
dim(ColDEG)

#final DEGs list for Colony!!!!!
write.csv(as.data.frame(ColDEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/LRT/DEGsColony.csv")


#Pollen DEGs filtered
PolDEG <- as.data.frame(GLPollen)
head(PolDEG)
dim(PolDEG)
PolDEG <- subset(PolDEG, select=c(2, 6))
head(PolDEG)
dim(PolDEG)
PolDEG <- tibble::rownames_to_column(PolDEG, "Locus")
head(PolDEG)
dim(PolDEG)
PolDEG$direction <- ifelse(PolDEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(PolDEG)
dim(PolDEG)
PolDEG <- merge(PolDEG,BimpLocBimp, by = "Locus", all = TRUE)
head(PolDEG)
dim(PolDEG)
PolDEG <- PolDEG[!is.na(PolDEG$direction), ]
dim(PolDEG)
head(PolDEG, 30)
PolDEG <- PolDEG[!duplicated(PolDEG$Locus), ]
head(PolDEG, 30)
dim(PolDEG)

#DEGs list for POllen no Filtered!!!!!
write.csv(as.data.frame(PolDEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/LRT/DEGsPollenNoFiltered.csv")


#Pesticide DEGs
PesDEG <- as.data.frame(GLPes)
head(PesDEG)
dim(PesDEG)
PesDEG <- subset(PesDEG, select=c(2, 6))
head(PesDEG)
dim(PesDEG)
PesDEG <- tibble::rownames_to_column(PesDEG, "Locus")
head(PesDEG)
dim(PesDEG)
PesDEG$direction <- ifelse(PesDEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(PesDEG)
dim(PesDEG)
PesDEG <- merge(PesDEG,BimpLocBimp, by = "Locus", all = TRUE)
head(PesDEG)
dim(PesDEG)
PesDEG <- PesDEG[!is.na(PesDEG$direction), ]
dim(PesDEG)
head(PesDEG, 30)
PesDEG <- PesDEG[!duplicated(PesDEG$Locus), ]
head(PesDEG, 30)
dim(PesDEG)

#DEGs list for POllen no Filtered!!!!!
write.csv(as.data.frame(PesDEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/LRT/DEGsPesticideNoFiltered.csv")


#Interaction DEGs
IntDEG <- as.data.frame(GLInteraction)
head(IntDEG)
dim(IntDEG)
IntDEG <- subset(IntDEG, select=c(2, 6))
head(IntDEG)
dim(IntDEG)
IntDEG <- tibble::rownames_to_column(IntDEG, "Locus")
head(IntDEG)
dim(IntDEG)
IntDEG$direction <- ifelse(IntDEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(IntDEG)
dim(IntDEG)
IntDEG <- merge(IntDEG,BimpLocBimp, by = "Locus", all = TRUE)
head(IntDEG)
dim(IntDEG)
IntDEG <- IntDEG[!is.na(IntDEG$direction), ]
dim(IntDEG)
head(IntDEG, 30)
IntDEG <- IntDEG[!duplicated(IntDEG$Locus), ]
head(IntDEG, 30)
dim(IntDEG)

#DEGs list for POllen no Filtered!!!!!
write.csv(as.data.frame(IntDEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/LRT/DEGsInteractionNoFiltered.csv")


##Filtering
###NOTES: DEGs related to natal colony and the interaction of diet and treatment with pesticide were subsequently excluded from the Pollen (Diet) and Pesticide DEG lists for enrichment analysis in order to identify genes whose expression is most impacted by these two factors, independently

#Filter by colony

Pollen_FilterByCol <- PolDEG[-which(PolDEG$Locus %in% ColDEG$Locus & PolDEG$direction %in% ColDEG$direction ),]
head(Pollen_FilterByCol)
dim(Pollen_FilterByCol)

Pesticide_FilterByCol <- PesDEG[-which(PesDEG$Locus %in% ColDEG$Locus & PesDEG$direction %in% ColDEG$direction ),]
head(Pesticide_FilterByCol)
dim(Pesticide_FilterByCol)

Interaction_FilterByCol <- IntDEG[-which(IntDEG$Locus %in% ColDEG$Locus & IntDEG$direction %in% ColDEG$direction ),]
head(Interaction_FilterByCol)
dim(Interaction_FilterByCol)

#Filter by Interaction

DEGsPollen <- Pollen_FilterByCol[-which(Pollen_FilterByCol$Locus %in% Interaction_FilterByCol$Locus & Pollen_FilterByCol$direction %in% Interaction_FilterByCol$direction ),]
summary(DEGsPollen)
head(DEGsPollen)
dim(DEGsPollen)

DEGsPesticide <- Pesticide_FilterByCol[-which(Pesticide_FilterByCol$Locus %in% Interaction_FilterByCol$Locus & Pesticide_FilterByCol$direction %in% Interaction_FilterByCol$direction ),]
summary(DEGsPesticide)
head(DEGsPesticide)
dim(DEGsPesticide)

#Filter POLLEN by Pesticide

DEGsPollen1 <- DEGsPollen[-which(DEGsPollen$Locus %in% DEGsPesticide$Locus & DEGsPollen$direction %in% DEGsPesticide$direction ),]
summary(DEGsPollen1)
head(DEGsPollen1)
dim(DEGsPollen1)

#final DEGs list for Pollen - Filtered!!!!!
write.csv(as.data.frame(DEGsPollen1),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/LRT/DEGsPollen.csv")

#Filter PESTICIDE by Pollen

DEGsPesticide1 <- DEGsPesticide[-which(DEGsPesticide$Locus %in% DEGsPollen$Locus & DEGsPesticide$direction %in% DEGsPollen$direction ),]
summary(DEGsPesticide1)
head(DEGsPesticide1)
dim(DEGsPesticide1)

#final DEGs list for Pesticide - Filtered!!!!!
write.csv(as.data.frame(DEGsPesticide1),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/LRT/DEGsPesticide.csv")

#DEGs excluisve for interaction

Interaction_FilterByPollen <- Interaction_FilterByCol[-which(Interaction_FilterByCol$Locus %in% Pollen_FilterByCol$Locus & Interaction_FilterByCol$direction %in% Pollen_FilterByCol$direction ),]
head(Interaction_FilterByPollen)
dim(Interaction_FilterByPollen)

DEGsInteraction <- Interaction_FilterByPollen[-which(Interaction_FilterByPollen$Locus %in% Pesticide_FilterByCol$Locus & Interaction_FilterByPollen$direction %in% Pesticide_FilterByCol$direction ),]
summary(DEGsInteraction)
head(DEGsInteraction)
dim(DEGsInteraction)

#final DEGs list for Interaction - Filtered!!!!!
write.csv(as.data.frame(DEGsInteraction),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/LRT/DEGsInteraction.csv")


#### GO enrichment analyses for LRT Gene lists ####

#This script for gene ontology (GO) enrichment analysis was based on Colgan et al 2019 (https://github.com/wurmlab/Bter_neonicotinoid_exposure_experiment/blob/master/07_gene_ontology_enrichment_analysis/go_enrichment_analysis.Rmd)
#This script will prepare data for GO analysis and create a 'TopGO object' from which enrichment tests can be performed to explore GO terms significantly enriched within the dataset. 
#This script outputs a results table of significantly enriched GO terms.
#It should be run once per list, and will produce 3 output files each time it is run

#install packages
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("topGO")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("lintr")
# 
# if (!requireNamespace("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# BiocManager::install("lattice")

#Two input files are required for the running of the test:  
# 1) A GO database file: Locus (gene ID or transcript) and GO terms (Table was generated using code in "GoInputCode.R")
# 2) A genelist file:  DEGs ranked by p-values. 

library("topGO")
library("lintr")

####DEGs Interaction List
####Step 1: Read in GO annotations:
#input (1) for GO database - use for all DEGs list the same GO database
gene_to_go_mapping <- readMappings("RNAseq/Data/BimpGO.txt", sep = "\t", IDsep = ",")

#input (2) DEGs list - Interaction list - one input different for each DEGs
IntGO <- as.data.frame(resInteraction)
head(IntGO)
dim(IntGO)
IntGO <- subset(IntGO, select=c(5)) #select raw p-value 
head(IntGO)
dim(IntGO)
IntGO <- tibble::rownames_to_column(IntGO, "Locus")
head(IntGO)
dim(IntGO)
IntGO <- IntGO[order(IntGO$pvalue), ]
IntGO <- subset(x = IntGO, subset = !is.na(pvalue))
#Convert into topgo's genelist format:
topgo_genelist <- IntGO$pvalue
names(topgo_genelist) <- IntGO$Locus
#Define a cut-off for running fisher's exact test:  
cutoff_for_top_fivepercent <- quantile(x = topgo_genelist, probs = 0.05)

#create output folder
#1st create a folder to save results for GO analysis
dir.create("RNAseq/Results.GO")
#now output directory for each DEGs list
output_directory <- "~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results.GO/output_LRT.Int"
if (file.exists(output_directory)) {
  stop("The output directory:", output_directory, ", already exists",
       "Let's avoid overwriting")
} else {
  dir.create(output_directory)
}

for (go_category in c("BP", "MF", "CC")) {
  # STEP 2
  ## Build the GOdata object in topGO
  my_go_data <- new("topGOdata",
                    description = paste("GOtest", go_category, sep = "_"),
                    ontology    = go_category,
                    geneSel     = function(x) {
                      # fails to run without thi
                      return(x <= cutoff_for_top_fivepercent)
                    },
                    allGenes    = topgo_genelist,
                    gene2GO     = gene_to_go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = 100) # Modify to reduce/increase stringency.
  # STEP THREE
  ## Calculate ks test using 'weight01' algorithm:
  result_weight_ks <- runTest(object    = my_go_data,
                              algorithm = "weight01",
                              statistic = "ks")
  ## Calculate fisher exact test using 'weight01' algorithm:
  result_weight_fisher <- runTest(object    = my_go_data,
                                  algorithm = "weight01",
                                  statistic = "fisher")
  ## Combine results from statistical tests:
  result_weight_output <- GenTable(object    = my_go_data,
                                   weight_ks = result_weight_ks,
                                   weight_fisher = result_weight_fisher,
                                   orderBy   = "weight_ks",
                                   topNodes  = length(score(result_weight_ks)))
  ## Correct ks test for multiple testing:
  result_weight_output$weight_ks <- as.numeric(result_weight_output$weight_ks)
  result_weight_output$weight_fisher <- as.numeric(result_weight_output$weight_fisher)
  result_weight_output$weight_ks_adjusted <- p.adjust(p = result_weight_output$weight_ks,
                                                      method = c("fdr"))
  result_weight_output$weight_fisher_adjusted <- p.adjust(p = result_weight_output$weight_fisher,
                                                          method = c("fdr"))
  ## Subset calls with significance higher than expected:
  result_weight_output_sig <- subset(x      = result_weight_output,
                                     subset = (Significant > Expected) &
                                       (weight_ks < 0.05))
  ## Update column names:
  colnames(result_weight_output_sig)[6] <- gsub(" ", "_",
                                                colnames(result_weight_output_sig)[6])
  ## For significant terms, add an additional field called 'category' which
  ## will be used for plotting of individual go category:
  result_weight_output_sig$category <- go_category
  ## Remove gaps between Terms, which can cause downstream problems:
  result_weight_output_sig$Term <- gsub(" ", "_", result_weight_output_sig$Term)
  ## Print to console one of the GO terms of interest to check the distribution of that GO term across ranked genes:
  print(showGroupDensity(object  = my_go_data,
                         whichGO = head(result_weight_output_sig$GO.ID,
                                        n = 1),
                         ranks   = TRUE,
                         rm.one  = FALSE))
  ## Write to output:
  write.table(x         = result_weight_output_sig,
              file      = file.path(output_directory,
                                    paste(go_category,
                                          "sig.tsv", sep = "_")),
              row.names = FALSE,
              quote = FALSE,
              sep       = "\t")
}

####DEGs Pesticide List
####Step 1: Read in GO annotations:
#input (1) for GO database - use for all DEGs list the same GO database, already done above.
#gene_to_go_mapping <- readMappings("Data/BimpGO.txt", sep = "\t", IDsep = ",")

#input (2) DEGs list - Pesticide list
PesGO <- as.data.frame(resPes)
head(PesGO)
dim(PesGO)
PesGO <- subset(PesGO, select=c(5)) #select raw p-value 
head(PesGO)
dim(PesGO)
PesGO <- tibble::rownames_to_column(PesGO, "Locus")
head(PesGO)
dim(PesGO)
PesGO <- PesGO[order(PesGO$pvalue), ]
PesGO <- subset(x = PesGO, subset = !is.na(pvalue))
#Convert into topgo's genelist format:
topgo_genelist <- PesGO$pvalue
names(topgo_genelist) <- PesGO$Locus
#Define a cut-off for running fisher's exact test:  
cutoff_for_top_fivepercent <- quantile(x = topgo_genelist, probs = 0.05)

#create output folder
#now output directory for each DEGs list
output_directory <- "~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results.GO/output_LRT.Pes"
if (file.exists(output_directory)) {
  stop("The output directory:", output_directory, ", already exists",
       "Let's avoid overwriting")
} else {
  dir.create(output_directory)
}

for (go_category in c("BP", "MF", "CC")) {
  # STEP 2
  ## Build the GOdata object in topGO
  my_go_data <- new("topGOdata",
                    description = paste("GOtest", go_category, sep = "_"),
                    ontology    = go_category,
                    geneSel     = function(x) {
                      # fails to run without thi
                      return(x <= cutoff_for_top_fivepercent)
                    },
                    allGenes    = topgo_genelist,
                    gene2GO     = gene_to_go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = 100) # Modify to reduce/increase stringency.
  # STEP THREE
  ## Calculate ks test using 'weight01' algorithm:
  result_weight_ks <- runTest(object    = my_go_data,
                              algorithm = "weight01",
                              statistic = "ks")
  ## Calculate fisher exact test using 'weight01' algorithm:
  result_weight_fisher <- runTest(object    = my_go_data,
                                  algorithm = "weight01",
                                  statistic = "fisher")
  ## Combine results from statistical tests:
  result_weight_output <- GenTable(object    = my_go_data,
                                   weight_ks = result_weight_ks,
                                   weight_fisher = result_weight_fisher,
                                   orderBy   = "weight_ks",
                                   topNodes  = length(score(result_weight_ks)))
  ## Correct ks test for multiple testing:
  result_weight_output$weight_ks <- as.numeric(result_weight_output$weight_ks)
  result_weight_output$weight_fisher <- as.numeric(result_weight_output$weight_fisher)
  result_weight_output$weight_ks_adjusted <- p.adjust(p = result_weight_output$weight_ks,
                                                      method = c("fdr"))
  result_weight_output$weight_fisher_adjusted <- p.adjust(p = result_weight_output$weight_fisher,
                                                          method = c("fdr"))
  ## Subset calls with significance higher than expected:
  result_weight_output_sig <- subset(x      = result_weight_output,
                                     subset = (Significant > Expected) &
                                       (weight_ks < 0.05))
  ## Update column names:
  colnames(result_weight_output_sig)[6] <- gsub(" ", "_",
                                                colnames(result_weight_output_sig)[6])
  ## For significant terms, add an additional field called 'category' which
  ## will be used for plotting of individual go category:
  result_weight_output_sig$category <- go_category
  ## Remove gaps between Terms, which can cause downstream problems:
  result_weight_output_sig$Term <- gsub(" ", "_", result_weight_output_sig$Term)
  ## Print to console one of the GO terms of interest to check the distribution of that GO term across ranked genes:
  print(showGroupDensity(object  = my_go_data,
                         whichGO = head(result_weight_output_sig$GO.ID,
                                        n = 1),
                         ranks   = TRUE,
                         rm.one  = FALSE))
  ## Write to output:
  write.table(x         = result_weight_output_sig,
              file      = file.path(output_directory,
                                    paste(go_category,
                                          "sig.tsv", sep = "_")),
              row.names = FALSE,
              quote = FALSE,
              sep       = "\t")
}


####DEGs Pesticide List
####Step 1: Read in GO annotations:
#input (1) for GO database - use for all DEGs list the same GO database, already done above.
#gene_to_go_mapping <- readMappings("Data/BimpGO.txt", sep = "\t", IDsep = ",")

#input (2) DEGs list - Pollen list
PolGO <- as.data.frame(resPollen)
head(PolGO)
dim(PolGO)
PolGO <- subset(PolGO, select=c(5)) #select raw p-value 
head(PolGO)
dim(PolGO)
PolGO <- tibble::rownames_to_column(PolGO, "Locus")
head(PolGO)
dim(PolGO)
PolGO <- PolGO[order(PolGO$pvalue), ]
PolGO <- subset(x = PolGO, subset = !is.na(pvalue))
#Convert into topgo's genelist format:
topgo_genelist <- PolGO$pvalue
names(topgo_genelist) <- PolGO$Locus
#Define a cut-off for running fisher's exact test:  
cutoff_for_top_fivepercent <- quantile(x = topgo_genelist, probs = 0.05)

#create output folder
#now output directory for each DEGs list
output_directory <- "~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results.GO/output_LRT.Pol"
if (file.exists(output_directory)) {
  stop("The output directory:", output_directory, ", already exists",
       "Let's avoid overwriting")
} else {
  dir.create(output_directory)
}

for (go_category in c("BP", "MF", "CC")) {
  # STEP 2
  ## Build the GOdata object in topGO
  my_go_data <- new("topGOdata",
                    description = paste("GOtest", go_category, sep = "_"),
                    ontology    = go_category,
                    geneSel     = function(x) {
                      # fails to run without thi
                      return(x <= cutoff_for_top_fivepercent)
                    },
                    allGenes    = topgo_genelist,
                    gene2GO     = gene_to_go_mapping,
                    annot       = annFUN.gene2GO,
                    nodeSize    = 100) # Modify to reduce/increase stringency.
  # STEP THREE
  ## Calculate ks test using 'weight01' algorithm:
  result_weight_ks <- runTest(object    = my_go_data,
                              algorithm = "weight01",
                              statistic = "ks")
  ## Calculate fisher exact test using 'weight01' algorithm:
  result_weight_fisher <- runTest(object    = my_go_data,
                                  algorithm = "weight01",
                                  statistic = "fisher")
  ## Combine results from statistical tests:
  result_weight_output <- GenTable(object    = my_go_data,
                                   weight_ks = result_weight_ks,
                                   weight_fisher = result_weight_fisher,
                                   orderBy   = "weight_ks",
                                   topNodes  = length(score(result_weight_ks)))
  ## Correct ks test for multiple testing:
  result_weight_output$weight_ks <- as.numeric(result_weight_output$weight_ks)
  result_weight_output$weight_fisher <- as.numeric(result_weight_output$weight_fisher)
  result_weight_output$weight_ks_adjusted <- p.adjust(p = result_weight_output$weight_ks,
                                                      method = c("fdr"))
  result_weight_output$weight_fisher_adjusted <- p.adjust(p = result_weight_output$weight_fisher,
                                                          method = c("fdr"))
  ## Subset calls with significance higher than expected:
  result_weight_output_sig <- subset(x      = result_weight_output,
                                     subset = (Significant > Expected) &
                                       (weight_ks < 0.05))
  ## Update column names:
  colnames(result_weight_output_sig)[6] <- gsub(" ", "_",
                                                colnames(result_weight_output_sig)[6])
  ## For significant terms, add an additional field called 'category' which
  ## will be used for plotting of individual go category:
  result_weight_output_sig$category <- go_category
  ## Remove gaps between Terms, which can cause downstream problems:
  result_weight_output_sig$Term <- gsub(" ", "_", result_weight_output_sig$Term)
  ## Print to console one of the GO terms of interest to check the distribution of that GO term across ranked genes:
  print(showGroupDensity(object  = my_go_data,
                         whichGO = head(result_weight_output_sig$GO.ID,
                                        n = 1),
                         ranks   = TRUE,
                         rm.one  = FALSE))
  ## Write to output:
  write.table(x         = result_weight_output_sig,
              file      = file.path(output_directory,
                                    paste(go_category,
                                          "sig.tsv", sep = "_")),
              row.names = FALSE,
              quote = FALSE,
              sep       = "\t")
}


#### Gene lists by Wald-Test - Pollen ####

##construct a DESeqDataSet 
ddsWaldP <- DESeqDataSetFromMatrix(countData = readcounts,
                                  colData = coldata,
                                  design= ~ Pollen)
ddsWaldP
ddsWaldP$Pollen <- relevel(ddsWaldP$Pollen, ref = "Mixed")
keep <- rowSums(counts(ddsWaldP)) >= 10
ddsWaldP <- ddsWaldP[keep,]
ddsWaldP

ddsPollenW <- DESeq(ddsWaldP)
resultsNames(ddsPollenW)

Erica_vs_Mixed <- results(ddsPollenW, name = "Pollen_Erica_vs_Mixed", alpha = 0.05)
summary(Erica_vs_Mixed)
Erica_vs_Mixed <- subset(Erica_vs_Mixed, padj < 0.05)
summary(Erica_vs_Mixed)
dim(Erica_vs_Mixed)

Cistus_vs_Mixed <- results(ddsPollenW, name = "Pollen_Cistus_vs_Mixed", alpha = 0.05)
summary(Cistus_vs_Mixed)

contrast_E.C <- c("Pollen", "Erica", "Cistus") #as we have by resultsNames only two comparison, we need to use this code to get other comparisons. Here, we use "Cistus" as baseline. 
Erica_vs_Cistus <- results(ddsPollenW, contrast=contrast_E.C, alpha = 0.05)
summary(Erica_vs_Cistus)
Erica_vs_Cistus <- subset(Erica_vs_Cistus, padj < 0.05)
summary(Erica_vs_Cistus)
dim(Erica_vs_Cistus)


#### Downloading the interested Pollen comparison - Gene lists with information about genes for manuscript ####

#We will use the object from the 1st part of LRT: BiomartR and NCBI, then after we start to work with Pollen DEGs lists:

#Now create a folder to save results for this analysis. We already created the Results folder
dir.create("RNAseq/Results/Pairwise-Pollen") #we create a sub folder to save this analysis

#Now get the table for each condition:

#Erica_vs_Cistus DEGs
ExCDEG <- as.data.frame(Erica_vs_Cistus)
head(ExCDEG)
dim(ExCDEG)
ExCDEG <- subset(ExCDEG, select=c(2, 6))
head(ExCDEG)
dim(ExCDEG)
ExCDEG <- tibble::rownames_to_column(ExCDEG, "Locus")
head(ExCDEG)
dim(ExCDEG)
ExCDEG$direction <- ifelse(ExCDEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(ExCDEG)
dim(ExCDEG)
ExCDEG <- merge(ExCDEG,BimpLocBimp, by = "Locus", all = TRUE)
head(ExCDEG)
dim(ExCDEG)
ExCDEG <- ExCDEG[!is.na(ExCDEG$direction), ]
dim(ExCDEG)
head(ExCDEG, 30)
ExCDEG <- ExCDEG[!duplicated(ExCDEG$Locus), ]
head(ExCDEG, 30)
dim(ExCDEG)

write.csv(as.data.frame(ExCDEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/Pairwise-Pollen/DEGsErica_vs_Cistus.csv")


#Erica_vs_Mixed DEGs
MxEDEG <- as.data.frame(Erica_vs_Mixed)
head(MxEDEG)
dim(MxEDEG)
MxEDEG <- subset(MxEDEG, select=c(2, 6))
head(MxEDEG)
dim(MxEDEG)
MxEDEG <- tibble::rownames_to_column(MxEDEG, "Locus")
head(MxEDEG)
dim(MxEDEG)
MxEDEG$direction <- ifelse(MxEDEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(MxEDEG)
dim(MxEDEG)
MxEDEG <- merge(MxEDEG,BimpLocBimp, by = "Locus", all = TRUE)
head(MxEDEG)
dim(MxEDEG)
MxEDEG <- MxEDEG[!is.na(MxEDEG$direction), ]
dim(MxEDEG)
head(MxEDEG, 30)
MxEDEG <- MxEDEG[!duplicated(MxEDEG$Locus), ]
head(MxEDEG, 30)
dim(MxEDEG)

write.csv(as.data.frame(MxEDEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/Pairwise-Pollen/DEGsMixed_vs_Erica.csv")


#### Gene lists by Wald-Test - Pesticide ####

ddsWaldPe <- DESeqDataSetFromMatrix(countData = readcounts,
                                   colData = coldata,
                                   design= ~ Pesticide)
ddsWaldPe
ddsWaldPe$Pesticide <- relevel(ddsWaldPe$Pesticide, ref = "untreated")
keep <- rowSums(counts(ddsWaldPe)) >= 10
ddsWaldPe <- ddsWaldPe[keep,]
ddsWaldPe

ddsPesticideW <- DESeq(ddsWaldPe)
resultsNames(ddsPesticideW)

IMD.A_vs_untreated <- results(ddsPesticideW, name = "Pesticide_IMD.A_vs_untreated", alpha = 0.05)
summary(IMD.A_vs_untreated)
IMD.A_vs_untreated <- subset(IMD.A_vs_untreated, padj < 0.05)
summary(IMD.A_vs_untreated)
dim(IMD.A_vs_untreated)

IMD.B_vs_untreated <- results(ddsPesticideW, name = "Pesticide_IMD.B_vs_untreated", alpha = 0.05)
summary(IMD.B_vs_untreated)
IMD.B_vs_untreated <- subset(IMD.B_vs_untreated, padj < 0.05)
summary(IMD.B_vs_untreated)
dim(IMD.B_vs_untreated)

contrast_B.A <- c("Pesticide", "IMD.B", "IMD.A") #as we have by resultsNames only two comparison, we need to use this code to get other comparisons. Here, we use "IMD.A" as baseline. 
IMD.B_vs_IMD.A <- results(ddsPesticideW, contrast=contrast_B.A, alpha = 0.05)
summary(IMD.B_vs_IMD.A)
IMD.B_vs_IMD.A <- subset(IMD.B_vs_IMD.A, padj < 0.05)
summary(IMD.B_vs_IMD.A)
dim(IMD.B_vs_IMD.A)


#### Downloading the interested Pesticide comparison - Gene lists with information about genes for manuscript ####

#We will use the object from the 1st part of LRT: BiomartR and NCBI, then after we start to work with Pollen DEGs lists:

#Now create a folder to save results for this analysis. We already created the Results folder
dir.create("RNAseq/Results/Pairwise-Pesticide") #we create a sub folder to save this analysis

#Now get the table for each condition:

#IMD.B_vs_IMD.A DEGs
BxADEG <- as.data.frame(IMD.B_vs_IMD.A)
head(BxADEG)
dim(BxADEG)
BxADEG <- subset(BxADEG, select=c(2, 6))
head(BxADEG)
dim(BxADEG)
BxADEG <- tibble::rownames_to_column(BxADEG, "Locus")
head(BxADEG)
dim(BxADEG)
BxADEG$direction <- ifelse(BxADEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(BxADEG)
dim(BxADEG)
BxADEG <- merge(BxADEG,BimpLocBimp, by = "Locus", all = TRUE)
head(BxADEG)
dim(BxADEG)
BxADEG <- BxADEG[!is.na(BxADEG$direction), ]
dim(BxADEG)
head(BxADEG, 30)
BxADEG <- BxADEG[!duplicated(BxADEG$Locus), ]
head(BxADEG, 30)
dim(BxADEG)

write.csv(as.data.frame(BxADEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/Pairwise-Pesticide/DEGsIMDB_vs_IMDA.csv")


#IMD.A_vs_untreated DEGs
CxADEG <- as.data.frame(IMD.A_vs_untreated)
head(CxADEG)
dim(CxADEG)
CxADEG <- subset(CxADEG, select=c(2, 6))
head(CxADEG)
dim(CxADEG)
CxADEG <- tibble::rownames_to_column(CxADEG, "Locus")
head(CxADEG)
dim(CxADEG)
CxADEG$direction <- ifelse(CxADEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(CxADEG)
dim(CxADEG)
CxADEG <- merge(CxADEG,BimpLocBimp, by = "Locus", all = TRUE)
head(CxADEG)
dim(CxADEG)
CxADEG <- CxADEG[!is.na(CxADEG$direction), ]
dim(CxADEG)
head(CxADEG, 30)
CxADEG <- CxADEG[!duplicated(CxADEG$Locus), ]
head(CxADEG, 30)
dim(CxADEG)

write.csv(as.data.frame(CxADEG),
     file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/Pairwise-Pesticide/DEGsControl_vs_IMDA.csv")


#IMD.B_vs_untreated DEGs
CxBDEG <- as.data.frame(IMD.B_vs_untreated)
head(CxBDEG)
dim(CxBDEG)
CxBDEG <- subset(CxBDEG, select=c(2, 6))
head(CxBDEG)
dim(CxBDEG)
CxBDEG <- tibble::rownames_to_column(CxBDEG, "Locus")
head(CxBDEG)
dim(CxBDEG)
CxBDEG$direction <- ifelse(CxBDEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(CxBDEG)
dim(CxBDEG)
CxBDEG <- merge(CxBDEG,BimpLocBimp, by = "Locus", all = TRUE)
head(CxBDEG)
dim(CxBDEG)
CxBDEG <- CxBDEG[!is.na(CxBDEG$direction), ]
dim(CxBDEG)
head(CxBDEG, 30)
CxBDEG <- CxBDEG[!duplicated(CxBDEG$Locus), ]
head(CxBDEG, 30)
dim(CxBDEG)

write.csv(as.data.frame(CxBDEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/Pairwise-Pesticide/DEGsControl_vs_IMDB.csv")


#### Gene lists by Wald-Test - Group ####

ddsWald <- DESeqDataSetFromMatrix(countData = readcounts,
                                  colData = coldata,
                                  design= ~ Group)

ddsWald

keep <- rowSums(counts(ddsWald)) >= 10
ddsWald <- ddsWald[keep,]
ddsWald

ddsGroup <- DESeq(ddsWald)
resultsNames(ddsGroup)

# 1. List Group by Diet-H ####
H_IMDA_vs_H_CTL <- results(ddsGroup, name = "Group_H_IMDA_vs_H_CTL", alpha = 0.05)
summary(H_IMDA_vs_H_CTL)

H_IMDB_vs_H_CTL <- results(ddsGroup, name = "Group_H_IMDB_vs_H_CTL", alpha = 0.05)
summary(H_IMDB_vs_H_CTL)
H_IMDB_vs_H_CTL <- subset(H_IMDB_vs_H_CTL, padj < 0.05)
summary(H_IMDB_vs_H_CTL)
dim(H_IMDB_vs_H_CTL)

contrast_H_IMD <- c("Group", "H_IMDB", "H_IMDA") #as we have by resultsNames only two comparison, we need to use this code to get other comparisons. Here, we use "H_IMDA" as baseline. 
H_IMDA_vs_H_IMDB <- results(ddsGroup, contrast=contrast_H_IMD, alpha = 0.05)
summary(H_IMDA_vs_H_IMDB)
H_IMDA_vs_H_IMDB <- subset(H_IMDA_vs_H_IMDB, padj < 0.05)
summary(H_IMDA_vs_H_IMDB)
dim(H_IMDA_vs_H_IMDB)

# 1. Downloading the interested Group by Diet-H - Gene lists with information about genes for manuscript ####

#We will use the object from the 1st part of LRT: BiomartR and NCBI, then after we start to work with Group by Diet-H DEGs lists:

#Now create a folder to save results for this analysis. We already created the Results folder
dir.create("RNAseq/Results/Group-Diet-H") #we create a sub folder to save this analysis

#Now get the table for each condition:

#H_IMDB_vs_H_CTL DEGs
HbxHcDEG <- as.data.frame(H_IMDB_vs_H_CTL)
head(HbxHcDEG)
dim(HbxHcDEG)
HbxHcDEG <- subset(HbxHcDEG, select=c(2, 6))
head(HbxHcDEG)
dim(HbxHcDEG)
HbxHcDEG <- tibble::rownames_to_column(HbxHcDEG, "Locus")
head(HbxHcDEG)
dim(HbxHcDEG)
HbxHcDEG$direction <- ifelse(HbxHcDEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(HbxHcDEG)
dim(HbxHcDEG)
HbxHcDEG <- merge(HbxHcDEG,BimpLocBimp, by = "Locus", all = TRUE)
head(HbxHcDEG)
dim(HbxHcDEG)
HbxHcDEG <- HbxHcDEG[!is.na(HbxHcDEG$direction), ]
dim(HbxHcDEG)
head(HbxHcDEG, 30)
HbxHcDEG <- HbxHcDEG[!duplicated(HbxHcDEG$Locus), ]
head(HbxHcDEG, 30)
dim(HbxHcDEG)

write.csv(as.data.frame(HbxHcDEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/Group-Diet-H/DEGsH_IMDB_vs_H_CTL.csv")


#H_IMDA_vs_H_IMDB DEGs
HaxHbDEG <- as.data.frame(H_IMDA_vs_H_IMDB)
head(HaxHbDEG)
dim(HaxHbDEG)
HaxHbDEG <- subset(HaxHbDEG, select=c(2, 6))
head(HaxHbDEG)
dim(HaxHbDEG)
HaxHbDEG <- tibble::rownames_to_column(HaxHbDEG, "Locus")
head(HaxHbDEG)
dim(HaxHbDEG)
HaxHbDEG$direction <- ifelse(HaxHbDEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(HaxHbDEG)
dim(HaxHbDEG)
HaxHbDEG <- merge(HaxHbDEG,BimpLocBimp, by = "Locus", all = TRUE)
head(HaxHbDEG)
dim(HaxHbDEG)
HaxHbDEG <- HaxHbDEG[!is.na(HaxHbDEG$direction), ]
dim(HaxHbDEG)
head(HaxHbDEG, 30)
HaxHbDEG <- HaxHbDEG[!duplicated(HaxHbDEG$Locus), ]
head(HaxHbDEG, 30)
dim(HaxHbDEG)

write.csv(as.data.frame(HaxHbDEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/Group-Diet-H/DEGsH_IMDA_vs_H_IMDB.csv")


# 2. List Group by Diet-M ####
contrast_M_IMDA <- c("Group", "M_IMDA", "M_CTL")
M_IMDA_vs_M_CTL <- results(ddsGroup, contrast=contrast_M_IMDA, alpha = 0.05)
summary(M_IMDA_vs_M_CTL)

contrast_M_IMDB <- c("Group", "M_IMDB", "M_CTL")
M_IMDB_vs_M_CTL <- results(ddsGroup, contrast=contrast_M_IMDB, alpha = 0.05)
summary(M_IMDB_vs_M_CTL)

contrast_M_IMD <- c("Group", "M_IMDB", "M_IMDA")
M_IMDA_vs_M_IMDB <- results(ddsGroup, contrast=contrast_M_IMD, alpha = 0.05)
summary(M_IMDA_vs_M_IMDB)


# 3. List Group by Diet-R ####
contrast_R_IMDA <- c("Group", "R_IMDA", "R_CTL")
R_IMDA_vs_R_CTL <- results(ddsGroup, contrast=contrast_R_IMDA, alpha = 0.05)
summary(R_IMDA_vs_R_CTL)
R_IMDA_vs_R_CTL <- subset(R_IMDA_vs_R_CTL, padj < 0.05)
summary(R_IMDA_vs_R_CTL)
dim(R_IMDA_vs_R_CTL)

contrast_R_IMDB <- c("Group", "R_IMDB", "R_CTL")
R_IMDB_vs_R_CTL <- results(ddsGroup, contrast=contrast_R_IMDB, alpha = 0.05)
summary(R_IMDB_vs_R_CTL)

contrast_R_IMD <- c("Group", "R_IMDB", "R_IMDA")
R_IMDA_vs_R_IMDB <- results(ddsGroup, contrast=contrast_R_IMD, alpha = 0.05)
summary(R_IMDA_vs_R_IMDB)


# 3. Downloading the interested Group by R - Gene lists with information about genes for manuscript ####

#We will use the object from the 1st part of LRT: BiomartR and NCBI, then after we start to work with Group by R DEGs lists:

#Now create a folder to save results for this analysis. We already created the Results folder
dir.create("RNAseq/Results/Group-Diet-R") #we create a sub folder to save this analysis

#Now get the table for each condition:

#R_IMDA_vs_R_CTL DEGs
RaxRcDEG <- as.data.frame(R_IMDA_vs_R_CTL)
head(RaxRcDEG)
dim(RaxRcDEG)
RaxRcDEG <- subset(RaxRcDEG, select=c(2, 6))
head(RaxRcDEG)
dim(RaxRcDEG)
RaxRcDEG <- tibble::rownames_to_column(RaxRcDEG, "Locus")
head(RaxRcDEG)
dim(RaxRcDEG)
RaxRcDEG$direction <- ifelse(RaxRcDEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(RaxRcDEG)
dim(RaxRcDEG)
RaxRcDEG <- merge(RaxRcDEG,BimpLocBimp, by = "Locus", all = TRUE)
head(RaxRcDEG)
dim(RaxRcDEG)
RaxRcDEG <- RaxRcDEG[!is.na(RaxRcDEG$direction), ]
dim(RaxRcDEG)
head(RaxRcDEG, 30)
RaxRcDEG <- RaxRcDEG[!duplicated(RaxRcDEG$Locus), ]
head(RaxRcDEG, 30)
dim(RaxRcDEG)

write.csv(as.data.frame(RaxRcDEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/Group-Diet-R/DEGsR_IMDA_vs_R_CTL.csv")


# 4. List Group by Pesticide-IMD-A ####
contrast_H.R_A <- c("Group", "H_IMDA", "R_IMDA")
H_IMDA_vs_R_IMDA <- results(ddsGroup, contrast=contrast_H.R_A, alpha = 0.05)
summary(H_IMDA_vs_R_IMDA)
H_IMDA_vs_R_IMDA <- subset(H_IMDA_vs_R_IMDA, padj < 0.05)
summary(H_IMDA_vs_R_IMDA)
dim(H_IMDA_vs_R_IMDA)

contrast_M.R_A <- c("Group", "R_IMDA", "M_IMDA")
M_IMDA_vs_R_IMDA <- results(ddsGroup, contrast=contrast_M.R_A, alpha = 0.05)
summary(M_IMDA_vs_R_IMDA)

contrast_M.H_A <- c("Group", "H_IMDA", "M_IMDA")
M_IMDA_vs_H_IMDA <- results(ddsGroup, contrast=contrast_M.H_A, alpha = 0.05)
summary(M_IMDA_vs_H_IMDA)
M_IMDA_vs_H_IMDA <- subset(M_IMDA_vs_H_IMDA, padj < 0.05)
summary(M_IMDA_vs_H_IMDA)
dim(M_IMDA_vs_H_IMDA)


# 4. Downloading the interested Group by Pesticide-IMD-A - Gene lists with information about genes for manuscript ####

#We will use the object from the 1st part of LRT: BiomartR and NCBI, then after we start to work with Group by Pesticide-IMD-A DEGs lists:

#Now create a folder to save results for this analysis. We already created the Results folder
dir.create("RNAseq/Results/Group-Pesticide-IMD-A") #we create a sub folder to save this analysis

#Now get the table for each condition:

#H_IMDA_vs_R_IMDA DEGs
aHxaRDEG <- as.data.frame(H_IMDA_vs_R_IMDA)
head(aHxaRDEG)
dim(aHxaRDEG)
aHxaRDEG <- subset(aHxaRDEG, select=c(2, 6))
head(aHxaRDEG)
dim(aHxaRDEG)
aHxaRDEG <- tibble::rownames_to_column(aHxaRDEG, "Locus")
head(aHxaRDEG)
dim(aHxaRDEG)
aHxaRDEG$direction <- ifelse(aHxaRDEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(aHxaRDEG)
dim(aHxaRDEG)
aHxaRDEG <- merge(aHxaRDEG,BimpLocBimp, by = "Locus", all = TRUE)
head(aHxaRDEG)
dim(aHxaRDEG)
aHxaRDEG <- aHxaRDEG[!is.na(aHxaRDEG$direction), ]
dim(aHxaRDEG)
head(aHxaRDEG, 30)
aHxaRDEG <- aHxaRDEG[!duplicated(aHxaRDEG$Locus), ]
head(aHxaRDEG, 30)
dim(aHxaRDEG)

write.csv(as.data.frame(aHxaRDEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/Group-Pesticide-IMD-A/DEGsH_IMDA_vs_R_IMDA.csv")


#M_IMDA_vs_H_IMDA DEGs
aMxaHDEG <- as.data.frame(M_IMDA_vs_H_IMDA)
head(aMxaHDEG)
dim(aMxaHDEG)
aMxaHDEG <- subset(aMxaHDEG, select=c(2, 6))
head(aMxaHDEG)
dim(aMxaHDEG)
aMxaHDEG <- tibble::rownames_to_column(aMxaHDEG, "Locus")
head(aMxaHDEG)
dim(aMxaHDEG)
aMxaHDEG$direction <- ifelse(aMxaHDEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(aMxaHDEG)
dim(aMxaHDEG)
aMxaHDEG <- merge(aMxaHDEG,BimpLocBimp, by = "Locus", all = TRUE)
head(aMxaHDEG)
dim(aMxaHDEG)
aMxaHDEG <- aMxaHDEG[!is.na(aMxaHDEG$direction), ]
dim(aMxaHDEG)
head(aMxaHDEG, 30)
aMxaHDEG <- aMxaHDEG[!duplicated(aMxaHDEG$Locus), ]
head(aMxaHDEG, 30)
dim(aMxaHDEG)

write.csv(as.data.frame(aMxaHDEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/Group-Pesticide-IMD-A/DEGsM_IMDA_vs_H_IMDA.csv")


# 5. List Group by Pesticide-IMD-B ####
contrast_H.R_B <- c("Group", "H_IMDB", "R_IMDB")
H_IMDB_vs_R_IMDB <- results(ddsGroup, contrast=contrast_H.R_B, alpha = 0.05)
summary(H_IMDB_vs_R_IMDB)
H_IMDB_vs_R_IMDB <- subset(H_IMDB_vs_R_IMDB, padj < 0.05)
summary(H_IMDB_vs_R_IMDB)
dim(H_IMDB_vs_R_IMDB)

contrast_M.R_B <- c("Group", "R_IMDB", "M_IMDB")
M_IMDB_vs_R_IMDB <- results(ddsGroup, contrast=contrast_M.R_B, alpha = 0.05)
summary(M_IMDB_vs_R_IMDB)

contrast_M.H_B <- c("Group", "H_IMDB", "M_IMDB")
M_IMDB_vs_H_IMDB <- results(ddsGroup, contrast=contrast_M.H_B, alpha = 0.05)
summary(M_IMDB_vs_H_IMDB)
M_IMDB_vs_H_IMDB <- subset(M_IMDB_vs_H_IMDB, padj < 0.05)
summary(M_IMDB_vs_H_IMDB)
dim(M_IMDB_vs_H_IMDB)


# 5. Downloading the interested Group by Pesticide-IMD-B - Gene lists with information about genes for manuscript ####

#We will use the object from the 1st part of LRT: BiomartR and NCBI, then after we start to work with Group by Pesticide-IMD-B DEGs lists:

#Now create a folder to save results for this analysis. We already created the Results folder
dir.create("RNAseq/Results/Group-Pesticide-IMD-B") #we create a sub folder to save this analysis

#Now get the table for each condition:

#H_IMDB_vs_R_IMDB DEGs
bHxbRDEG <- as.data.frame(H_IMDB_vs_R_IMDB)
head(bHxbRDEG)
dim(bHxbRDEG)
bHxbRDEG <- subset(bHxbRDEG, select=c(2, 6))
head(bHxbRDEG)
dim(bHxbRDEG)
bHxbRDEG <- tibble::rownames_to_column(bHxbRDEG, "Locus")
head(bHxbRDEG)
dim(bHxbRDEG)
bHxbRDEG$direction <- ifelse(bHxbRDEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(bHxbRDEG)
dim(bHxbRDEG)
bHxbRDEG <- merge(bHxbRDEG,BimpLocBimp, by = "Locus", all = TRUE)
head(bHxbRDEG)
dim(bHxbRDEG)
bHxbRDEG <- bHxbRDEG[!is.na(bHxbRDEG$direction), ]
dim(bHxbRDEG)
head(bHxbRDEG, 30)
bHxbRDEG <- bHxbRDEG[!duplicated(bHxbRDEG$Locus), ]
head(bHxbRDEG, 30)
dim(bHxbRDEG)

write.csv(as.data.frame(bHxbRDEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/Group-Pesticide-IMD-B/DEGsH_IMDB_vs_R_IMDB.csv")


#M_IMDB_vs_H_IMDB DEGs
bMxbHDEG <- as.data.frame(M_IMDB_vs_H_IMDB)
head(bMxbHDEG)
dim(bMxbHDEG)
bMxbHDEG <- subset(bMxbHDEG, select=c(2, 6))
head(bMxbHDEG)
dim(bMxbHDEG)
bMxbHDEG <- tibble::rownames_to_column(bMxbHDEG, "Locus")
head(bMxbHDEG)
dim(bMxbHDEG)
bMxbHDEG$direction <- ifelse(bMxbHDEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(bMxbHDEG)
dim(bMxbHDEG)
bMxbHDEG <- merge(bMxbHDEG,BimpLocBimp, by = "Locus", all = TRUE)
head(bMxbHDEG)
dim(bMxbHDEG)
bMxbHDEG <- bMxbHDEG[!is.na(bMxbHDEG$direction), ]
dim(bMxbHDEG)
head(bMxbHDEG, 30)
bMxbHDEG <- bMxbHDEG[!duplicated(bMxbHDEG$Locus), ]
head(bMxbHDEG, 30)
dim(bMxbHDEG)

write.csv(as.data.frame(bMxbHDEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/Group-Pesticide-IMD-B/DEGsM_IMDB_vs_H_IMDB.csv")


# 6. List Group by Pesticide-Control ####
contrast_H.R_C <- c("Group", "H_CTL", "R_CTL")
H_CTL_vs_R_CTL <- results(ddsGroup, contrast=contrast_H.R_C, alpha = 0.05)
summary(H_CTL_vs_R_CTL)

contrast_M.R_C <- c("Group", "R_CTL", "M_CTL")
M_CTL_vs_R_CTL <- results(ddsGroup, contrast=contrast_M.R_C, alpha = 0.05)
summary(M_CTL_vs_R_CTL)

contrast_M.H_C <- c("Group", "H_CTL", "M_CTL")
M_CTL_vs_H_CTL <- results(ddsGroup, contrast=contrast_M.H_C, alpha = 0.05)
summary(M_CTL_vs_H_CTL)
M_CTL_vs_H_CTL <- subset(M_CTL_vs_H_CTL, padj < 0.05)
summary(M_CTL_vs_H_CTL)
dim(M_CTL_vs_H_CTL)


# 6. Downloading the interested Group by Pesticide-Control - Gene lists with information about genes for manuscript ####

#We will use the object from the 1st part of LRT: BiomartR and NCBI, then after we start to work with Group by Pesticide-Control DEGs lists:

#Now create a folder to save results for this analysis. We already created the Results folder
dir.create("RNAseq/Results/Group-Pesticide-Control") #we create a sub folder to save this analysis

#Now get the table for each condition:

#M_CTL_vs_H_CTL DEGs
cMxcHDEG <- as.data.frame(M_CTL_vs_H_CTL)
head(cMxcHDEG)
dim(cMxcHDEG)
cMxcHDEG <- subset(cMxcHDEG, select=c(2, 6))
head(cMxcHDEG)
dim(cMxcHDEG)
cMxcHDEG <- tibble::rownames_to_column(cMxcHDEG, "Locus")
head(cMxcHDEG)
dim(cMxcHDEG)
cMxcHDEG$direction <- ifelse(cMxcHDEG$log2FoldChange > 0,"up-regulated", "down-regulated")
head(cMxcHDEG)
dim(cMxcHDEG)
cMxcHDEG <- merge(cMxcHDEG,BimpLocBimp, by = "Locus", all = TRUE)
head(cMxcHDEG)
dim(cMxcHDEG)
cMxcHDEG <- cMxcHDEG[!is.na(cMxcHDEG$direction), ]
dim(cMxcHDEG)
head(cMxcHDEG, 30)
cMxcHDEG <- cMxcHDEG[!duplicated(cMxcHDEG$Locus), ]
head(cMxcHDEG, 30)
dim(cMxcHDEG)

write.csv(as.data.frame(cMxcHDEG),
          file="~/Dropbox/2022-PesticideQueens-Bimpatiens/RNAseq/Results/Group-Pesticide-Control/DEGsM_CTL_vs_H_CTL.csv")


#### Figures ####

####Fig - Volcanos####

#volcano Pollen DEGs
PollenDEGs.dataframe <- as.data.frame(results(ddsPollen, alpha=0.05))
head(PollenDEGs.dataframe)
dim(PollenDEGs.dataframe)

PollenDEGs.dataframe$log10_p <- -log10(PollenDEGs.dataframe$padj)
head(PollenDEGs.dataframe)
dim(PollenDEGs.dataframe)

threshold_Pol <- PollenDEGs.dataframe$padj <= 0.05 & PollenDEGs.dataframe$log2FoldChange >0 | PollenDEGs.dataframe$padj <= 0.05 & PollenDEGs.dataframe$log2FoldChange < 0
length(which(threshold_Pol))

PollenDEGs.dataframe$threshold <- threshold_Pol 
head(PollenDEGs.dataframe)
dim(PollenDEGs.dataframe)

ggplot(PollenDEGs.dataframe) +
  geom_point(aes(x=log2FoldChange, y=log10_p, colour=threshold), size = 0.5) +
  scale_x_continuous(name="log2 fold change", limits=c(-10, 10)) + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values=c("black", "red")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1))) 


threshold_Pol1 <- row.names(ddsPollen) %in% DEGsPollen$Locus
length(which(threshold_Pol1))

PollenDEGs.dataframe$threshold1 <- threshold_Pol1 
head(PollenDEGs.dataframe)
dim(PollenDEGs.dataframe)

ggplot(PollenDEGs.dataframe) +
  geom_point(aes(x=log2FoldChange, y=log10_p, colour=threshold1), size = 0.5) +
  scale_x_continuous(name="log2 fold change", limits=c(-10, 10)) + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values=c("black", "red")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1))) 


#volcano Pesticide DEGs
PesticideDEGs.dataframe <- as.data.frame(results(ddsPes, alpha=0.05))
head(PesticideDEGs.dataframe)
dim(PesticideDEGs.dataframe)

PesticideDEGs.dataframe$log10_p <- -log10(PesticideDEGs.dataframe$padj)
head(PesticideDEGs.dataframe)
dim(PesticideDEGs.dataframe)

threshold_Pest <- PesticideDEGs.dataframe$padj <= 0.05 & PesticideDEGs.dataframe$log2FoldChange >0 | PesticideDEGs.dataframe$padj <= 0.05 & PesticideDEGs.dataframe$log2FoldChange < 0
length(which(threshold_Pest))

PesticideDEGs.dataframe$threshold <- threshold_Pest 
head(PesticideDEGs.dataframe)
dim(PesticideDEGs.dataframe)

ggplot(PesticideDEGs.dataframe) +
  geom_point(aes(x=log2FoldChange, y=log10_p, colour=threshold), size = 0.5) +
  scale_x_continuous(name="log2 fold change", limits=c(-10, 10)) + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values=c("black", "red")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1))) 

threshold_Pest1 <- row.names(ddsPes) %in% DEGsPesticide$Locus
length(which(threshold_Pest1))

PesticideDEGs.dataframe$threshold1 <- threshold_Pest1 
head(PesticideDEGs.dataframe)
dim(PesticideDEGs.dataframe)

ggplot(PesticideDEGs.dataframe) +
  geom_point(aes(x=log2FoldChange, y=log10_p, colour=threshold1), size = 0.5) +
  scale_x_continuous(name="log2 fold change", limits=c(-10, 10)) + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values=c("black", "red")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1))) 


#volcano Interaction DEGs
InteractionDEGs.dataframe <- as.data.frame(results(ddsInteraction, alpha=0.05))
head(InteractionDEGs.dataframe)
dim(InteractionDEGs.dataframe)

InteractionDEGs.dataframe$log10_p <- -log10(InteractionDEGs.dataframe$padj)
head(InteractionDEGs.dataframe)
dim(InteractionDEGs.dataframe)

threshold_Int <- InteractionDEGs.dataframe$padj <= 0.05 & InteractionDEGs.dataframe$log2FoldChange >0 | InteractionDEGs.dataframe$padj <= 0.05 & InteractionDEGs.dataframe$log2FoldChange < 0
length(which(threshold_Int))

InteractionDEGs.dataframe$threshold <- threshold_Int 
head(InteractionDEGs.dataframe)
dim(InteractionDEGs.dataframe)

ggplot(InteractionDEGs.dataframe) +
  geom_point(aes(x=log2FoldChange, y=log10_p, colour=threshold), size = 0.5) +
  scale_x_continuous(name="log2 fold change", limits=c(-10, 10)) + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values=c("black", "red")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1))) 

threshold_Int1 <- row.names(ddsInteraction) %in% DEGsInteraction$Locus
length(which(threshold_Int1))

InteractionDEGs.dataframe$threshold1 <- threshold_Int1 
head(InteractionDEGs.dataframe)
dim(InteractionDEGs.dataframe)

ggplot(InteractionDEGs.dataframe) +
  geom_point(aes(x=log2FoldChange, y=log10_p, colour=threshold1), size = 0.5) +
  scale_x_continuous(name="log2 fold change", limits=c(-10, 10)) + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values=c("black", "red")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1))) 


#volcano Erica_vs_Cistus
Erica_vs_Cistus.dataframe <- as.data.frame(results(ddsPollenW, contrast=contrast_E.C, alpha = 0.05))
head(Erica_vs_Cistus.dataframe)
dim(Erica_vs_Cistus.dataframe)

Erica_vs_Cistus.dataframe$log10_p <- -log10(Erica_vs_Cistus.dataframe$padj)
head(Erica_vs_Cistus.dataframe)
dim(Erica_vs_Cistus.dataframe)

threshold_EC <- Erica_vs_Cistus.dataframe$padj <= 0.05 & Erica_vs_Cistus.dataframe$log2FoldChange >0 | Erica_vs_Cistus.dataframe$padj <= 0.05 & Erica_vs_Cistus.dataframe$log2FoldChange < 0
length(which(threshold_EC))

Erica_vs_Cistus.dataframe$threshold <- threshold_EC 
head(Erica_vs_Cistus.dataframe)
dim(Erica_vs_Cistus.dataframe)

vP1 <- ggplot(Erica_vs_Cistus.dataframe) +
  geom_point(aes(x=log2FoldChange, y=log10_p, colour=threshold), size = 0.5) +
  scale_x_continuous(name="log2 fold change", limits=c(-10, 10)) + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values=c("black", "red")) +
  theme_bw()+ggtitle("Diet-1 versus Diet-2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1))) 
vP1


#volcano Erica_vs_Mixed
Mixed_vs_Erica.dataframe <- as.data.frame(results(ddsPollenW, name = "Pollen_Erica_vs_Mixed", alpha = 0.05))
head(Mixed_vs_Erica.dataframe)
dim(Mixed_vs_Erica.dataframe)

Mixed_vs_Erica.dataframe$log10_p <- -log10(Mixed_vs_Erica.dataframe$padj)
head(Mixed_vs_Erica.dataframe)
dim(Mixed_vs_Erica.dataframe)

threshold_ME <- Mixed_vs_Erica.dataframe$padj <= 0.05 & Mixed_vs_Erica.dataframe$log2FoldChange >0 | Mixed_vs_Erica.dataframe$padj <= 0.05 & Mixed_vs_Erica.dataframe$log2FoldChange < 0
length(which(threshold_ME))

Mixed_vs_Erica.dataframe$threshold <- threshold_ME 
head(Mixed_vs_Erica.dataframe)
dim(Mixed_vs_Erica.dataframe)

Mixed_vs_Erica.dataframe <- Mixed_vs_Erica.dataframe[order(Mixed_vs_Erica.dataframe$padj), ] 
head(Mixed_vs_Erica.dataframe)
dim(Mixed_vs_Erica.dataframe)

vP2 <- ggplot(Mixed_vs_Erica.dataframe) +
  geom_point(aes(x=log2FoldChange, y=log10_p, colour=threshold), size = 0.5) +
  scale_x_continuous(name="log2 fold change", limits=c(-10, 10)) + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values=c("black", "red")) +
  theme_bw()+ggtitle("Combined diet versus Diet-2")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1))) 
vP2

#volcano IMD.B_vs_IMD.A
IMD.B_vs_IMD.A.dataframe <- as.data.frame(results(ddsPesticideW, contrast=contrast_B.A, alpha = 0.05))
head(IMD.B_vs_IMD.A.dataframe)
dim(IMD.B_vs_IMD.A.dataframe)

IMD.B_vs_IMD.A.dataframe$log10_p <- -log10(IMD.B_vs_IMD.A.dataframe$padj)
head(IMD.B_vs_IMD.A.dataframe)
dim(IMD.B_vs_IMD.A.dataframe)

threshold_BA <- IMD.B_vs_IMD.A.dataframe$padj <= 0.05 & IMD.B_vs_IMD.A.dataframe$log2FoldChange >0 | IMD.B_vs_IMD.A.dataframe$padj <= 0.05 & IMD.B_vs_IMD.A.dataframe$log2FoldChange < 0
length(which(threshold_BA))

IMD.B_vs_IMD.A.dataframe$threshold <- threshold_BA 
head(IMD.B_vs_IMD.A.dataframe)
dim(IMD.B_vs_IMD.A.dataframe)

IMD.B_vs_IMD.A.dataframe <- IMD.B_vs_IMD.A.dataframe[order(IMD.B_vs_IMD.A.dataframe$padj), ] 
head(IMD.B_vs_IMD.A.dataframe)
dim(IMD.B_vs_IMD.A.dataframe)

vI1 <- ggplot(IMD.B_vs_IMD.A.dataframe) +
  geom_point(aes(x=log2FoldChange, y=log10_p, colour=threshold), size = 0.5) +
  scale_x_continuous(name="log2 fold change", limits=c(-10, 10)) + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values=c("black", "red")) +
  theme_bw()+ggtitle("IMD-A versus IMD-B")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1))) 
vI1

#volcano untreated_vs_IMD.A
untreated_vs_IMD.A.dataframe <- as.data.frame(results(ddsPesticideW, name = "Pesticide_IMD.A_vs_untreated", alpha = 0.05))
head(untreated_vs_IMD.A.dataframe)
dim(untreated_vs_IMD.A.dataframe)

untreated_vs_IMD.A.dataframe$log10_p <- -log10(untreated_vs_IMD.A.dataframe$padj)
head(untreated_vs_IMD.A.dataframe)
dim(untreated_vs_IMD.A.dataframe)

threshold_unA <- untreated_vs_IMD.A.dataframe$padj <= 0.05 & untreated_vs_IMD.A.dataframe$log2FoldChange >0 | untreated_vs_IMD.A.dataframe$padj <= 0.05 & untreated_vs_IMD.A.dataframe$log2FoldChange < 0
length(which(threshold_unA))

untreated_vs_IMD.A.dataframe$threshold <- threshold_unA 
head(untreated_vs_IMD.A.dataframe)
dim(untreated_vs_IMD.A.dataframe)

untreated_vs_IMD.A.dataframe <- untreated_vs_IMD.A.dataframe[order(untreated_vs_IMD.A.dataframe$padj), ] 
head(untreated_vs_IMD.A.dataframe)
dim(untreated_vs_IMD.A.dataframe)

vI2 <- ggplot(untreated_vs_IMD.A.dataframe) +
  geom_point(aes(x=log2FoldChange, y=log10_p, colour=threshold), size = 0.5) +
  scale_x_continuous(name="log2 fold change", limits=c(-10, 10)) + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values=c("black", "red")) +
  theme_bw()+ggtitle("not exposed versus IMD-A")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1))) 
vI2

#volcano untreated_vs_IMD.B
untreated_vs_IMD.B.dataframe <- as.data.frame(results(ddsPesticideW, name = "Pesticide_IMD.B_vs_untreated", alpha = 0.05))
head(untreated_vs_IMD.B.dataframe)
dim(untreated_vs_IMD.B.dataframe)

untreated_vs_IMD.B.dataframe$log10_p <- -log10(untreated_vs_IMD.B.dataframe$padj)
head(untreated_vs_IMD.B.dataframe)
dim(untreated_vs_IMD.B.dataframe)

threshold_unB <- untreated_vs_IMD.B.dataframe$padj <= 0.05 & untreated_vs_IMD.B.dataframe$log2FoldChange >0 | untreated_vs_IMD.B.dataframe$padj <= 0.05 & untreated_vs_IMD.B.dataframe$log2FoldChange < 0
length(which(threshold_unB))

untreated_vs_IMD.B.dataframe$threshold <- threshold_unB 
head(untreated_vs_IMD.B.dataframe)
dim(untreated_vs_IMD.B.dataframe)

untreated_vs_IMD.B.dataframe <- untreated_vs_IMD.B.dataframe[order(untreated_vs_IMD.B.dataframe$padj), ] 
head(untreated_vs_IMD.B.dataframe)
dim(untreated_vs_IMD.B.dataframe)

vI3 <- ggplot(untreated_vs_IMD.B.dataframe) +
  geom_point(aes(x=log2FoldChange, y=log10_p, colour=threshold), size = 0.5) +
  scale_x_continuous(name="log2 fold change", limits=c(-10, 10)) + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values=c("black", "red")) +
  theme_bw()+ggtitle("not exposed versus IMD-B")+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1))) 
vI3

    ##Fig to export 
    ggarrange(vI2, vI3, vI1, vP2, vP1, labels = c("A", "B", "C", "D", "E"), ncol = 3, nrow = 2)


#volcano H_IMDB_vs_H_CTL
H_IMDB_vs_H_CTL.dataframe <- as.data.frame(results(ddsGroup, name = "Group_H_IMDB_vs_H_CTL", alpha = 0.05))
head(H_IMDB_vs_H_CTL.dataframe)
dim(H_IMDB_vs_H_CTL.dataframe)

H_IMDB_vs_H_CTL.dataframe$log10_p <- -log10(H_IMDB_vs_H_CTL.dataframe$padj)
head(H_IMDB_vs_H_CTL.dataframe)
dim(H_IMDB_vs_H_CTL.dataframe)

threshold_H.Bun <- H_IMDB_vs_H_CTL.dataframe$padj <= 0.05 & H_IMDB_vs_H_CTL.dataframe$log2FoldChange >0 | H_IMDB_vs_H_CTL.dataframe$padj <= 0.05 & H_IMDB_vs_H_CTL.dataframe$log2FoldChange < 0
length(which(threshold_H.Bun))

H_IMDB_vs_H_CTL.dataframe$threshold <- threshold_H.Bun 
head(H_IMDB_vs_H_CTL.dataframe)
dim(H_IMDB_vs_H_CTL.dataframe)

H_IMDB_vs_H_CTL.dataframe <- H_IMDB_vs_H_CTL.dataframe[order(H_IMDB_vs_H_CTL.dataframe$padj), ] 
head(H_IMDB_vs_H_CTL.dataframe)
dim(H_IMDB_vs_H_CTL.dataframe)

ggplot(H_IMDB_vs_H_CTL.dataframe) +
  geom_point(aes(x=log2FoldChange, y=log10_p, colour=threshold), size = 0.5) +
  scale_x_continuous(name="log2 fold change", limits=c(-10, 10)) + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values=c("black", "red")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1))) 


#volcano H_IMDA_vs_H_IMDB
H_IMDA_vs_H_IMDB.dataframe <- as.data.frame(results(ddsGroup, contrast=contrast_H_IMD, alpha = 0.05))
head(H_IMDA_vs_H_IMDB.dataframe)
dim(H_IMDA_vs_H_IMDB.dataframe)

H_IMDA_vs_H_IMDB.dataframe$log10_p <- -log10(H_IMDA_vs_H_IMDB.dataframe$padj)
head(H_IMDA_vs_H_IMDB.dataframe)
dim(H_IMDA_vs_H_IMDB.dataframe)

threshold_H.AB <- H_IMDA_vs_H_IMDB.dataframe$padj <= 0.05 & H_IMDA_vs_H_IMDB.dataframe$log2FoldChange >0 | H_IMDA_vs_H_IMDB.dataframe$padj <= 0.05 & H_IMDA_vs_H_IMDB.dataframe$log2FoldChange < 0
length(which(threshold_H.AB))

H_IMDA_vs_H_IMDB.dataframe$threshold <- threshold_H.AB 
head(H_IMDA_vs_H_IMDB.dataframe)
dim(H_IMDA_vs_H_IMDB.dataframe)

H_IMDA_vs_H_IMDB.dataframe <- H_IMDA_vs_H_IMDB.dataframe[order(H_IMDA_vs_H_IMDB.dataframe$padj), ] 
head(H_IMDA_vs_H_IMDB.dataframe)
dim(H_IMDA_vs_H_IMDB.dataframe)

ggplot(H_IMDA_vs_H_IMDB.dataframe) +
  geom_point(aes(x=log2FoldChange, y=log10_p, colour=threshold), size = 0.5) +
  scale_x_continuous(name="log2 fold change", limits=c(-10, 10)) + 
  ylab("-log10 adjusted p-value") +
  scale_colour_manual(values=c("black", "red")) +
  theme_bw()+
  theme(panel.grid.major = element_blank(), panel.grid.minor = element_blank())+
  theme(legend.position = "none",
        axis.title = element_text(size = rel(1))) 

####Fig - PCA####

#transform expression levels using the regularized log transformation
rld <-rlog(ddsData, blind = T)
#why do we use blind = FALSE? If many of genes have large differences in counts due to the experimental design, it is important to set blind=FALSE for downstream analysis.
#PCA - for rld-transformed data
#Gettin' fancier with the PCA
pca.dataRLD <- plotPCA(rld, intgroup=c("Pollen", "Pesticide"), returnData=TRUE)
percentVar <- round(100 * attr(pca.dataRLD, "percentVar"))
ggplot(pca.dataRLD, aes(PC1, PC2, color=Pollen, shape=Pesticide)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed() +ggtitle("RLD-transformed data")
#format the graphic 
p <- ggplot(pca.dataRLD, aes(PC1, PC2, color=Pollen, shape=Pesticide)) + geom_point(size=3) + xlab(paste0("PC1: ",percentVar[1],"% variance")) + ylab(paste0("PC2: ",percentVar[2],"% variance")) + coord_fixed()
p <- p + theme_bw()
p <- p + theme_classic()
p <- p + theme_light() +theme(legend.text=element_text(size=10)) 
p

####Fig - Heatmap####

#I used a FKPM data per group to run the heatmap and I based in this tutorial: https://igordot.github.io/tutorials/heatmaps-2017-07.nb.html

#read counts - Per group
fpkmGroup <- read.csv("RNAseq/Data/fpkm_group.csv", header = T, row.names = 1)
head(fpkmGroup)
dim(fpkmGroup)
heat <- zFPKM(fpkmGroup)
genes.sig <- HbxHcDEG #using the subset for DEGs 
heat <- heat[which(rownames(heat) %in% genes.sig$Locus), ]
heat <- as.matrix(sapply(heat, as.numeric))  
heat <- na.omit(heat)
is.numeric(heat[1:9])
head(heat)
dim(heat)

df <- read.csv("RNAseq/Data/coldataGroup.csv", header = T, row.names = 1)
dim(df)
head(df)
df$Diet = factor(df$Diet, levels = c("Diet-1", "Diet-2", "Combined diet"))
dietCol <- c("darkgreen", "darkolivegreen1","chartreuse3")
names(dietCol) <- levels(df$Diet)
df$Pesticide = factor(df$Pesticide, levels = c("not exposed", "IMD-A", "IMD-B"))
pesticideCol <- c("rosybrown2", "darkorange4", "darkorange2")
names(pesticideCol) <- levels(df$Pesticide)
AnnColour <- list(
  Diet = dietCol,
  Pesticide = pesticideCol)

htmG <- pheatmap(heat, cluster_rows = TRUE, show_rownames = FALSE, show_colnames = FALSE, cluster_cols = T, scale = "row", treeheight_row = 0,  annotation_col = df, annotation_colors = AnnColour,
         clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete")
htmG 


#read counts - Per Sample
fpkmGroupS <- read.csv("RNAseq/Data/fpkm_sample.csv", header = T, row.names = 1)
head(fpkmGroupS)
dim(fpkmGroupS)
heatS <- zFPKM(fpkmGroupS)
genes.sig <- HbxHcDEG #using the subset for DEGs 
heatS <- heatS[which(rownames(heatS) %in% genes.sig$Locus), ]
heatS <- as.matrix(sapply(heatS, as.numeric))  
heatS <- na.omit(heatS)
is.numeric(heatS[1:54])
heatS <- heatS[!is.infinite(rowSums(heatS)),]
head(heatS)
dim(heatS)

dfS <- read.csv("RNAseq/Data/coldata.csv", header = T, row.names = 1)
dim(dfS)
head(dfS)
dfS = subset(dfS, select = c(Pollen.diet,Pesticide))
dim(dfS)
head(dfS)
dfS[dfS == "Combined.diet"] <- "Combined diet"
dfS[dfS == "untreated"] <- "not exposed"
names(dfS)[1] <- "Diet"
dim(dfS)
head(dfS)
dfS$Diet = factor(dfS$Diet, levels = c("Diet-1", "Diet-2", "Combined diet"))
dietColS <- c("darkgreen", "darkolivegreen1","chartreuse3")
names(dietColS) <- levels(dfS$Diet)
dfS$Pesticide = factor(dfS$Pesticide, levels = c("not exposed", "IMD-A", "IMD-B"))
pesticideColS <- c("rosybrown2", "darkorange4", "darkorange2")
names(pesticideColS) <- levels(dfS$Pesticide)
AnnColourS <- list(
  Diet = dietColS,
  Pesticide = pesticideColS)

htmS <- pheatmap(heatS, cluster_rows = TRUE, show_rownames = FALSE, show_colnames = FALSE, cluster_cols = T, scale = "row", treeheight_row = 0,  annotation_col = dfS, annotation_colors = AnnColourS,
                clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean", clustering_method = "complete")
htmS 


####Fig - Plots for specific genes by Wald-Test Group####

#Diet Rockrose - diet 1
LOC100744789 <- counts(ddsGroup['LOC100744789',], normalized = TRUE)
r1 <- list(counts = as.numeric(LOC100744789), group = as.factor(coldata$Group))
r1 <- as.tibble(r1)
r1$group <- factor(r1$group, levels=c("R_CTL", "R_IMDA", "R_IMDB"))
r1 <- na.omit(r1)
g1 <- ggplot(r1, aes(group,counts)) + geom_boxplot(aes(fill = group))
g1 <- g1 + labs(x = "Pesticide exposure", y = "Normalized Counts ", title = "glucose dehydrogenase")
g1 <- g1 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
g1 <- g1 + scale_fill_manual(values=c("rosybrown2", "darkorange4", "darkorange2"))
g1 <- g1 + scale_x_discrete(breaks=c("R_CTL", "R_IMDA", "R_IMDB"), labels=c("not exposed", "IMD-A", "IMD-B"))
g1


LOC100742002 <- counts(ddsGroup['LOC100742002',], normalized = TRUE)
r2 <- list(counts = as.numeric(LOC100742002), group = as.factor(coldata$Group))
r2 <- as.tibble(r2)
r2$group <- factor(r2$group, levels=c("R_CTL", "R_IMDA", "R_IMDB"))
r2 <- na.omit(r2)
g2 <- ggplot(r2, aes(group,counts)) + geom_boxplot(aes(fill = group))
g2 <- g2 + labs(x = "Pesticide exposure", y = "Normalized Counts ", title = "odorant binding protein")
g2 <- g2 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
g2 <- g2 + scale_fill_manual(values=c("rosybrown2", "darkorange4", "darkorange2"))
g2 <- g2 + scale_x_discrete(breaks=c("R_CTL", "R_IMDA", "R_IMDB"), labels=c("not exposed", "IMD-A", "IMD-B"))
g2 

LOC100747941 <- counts(ddsGroup['LOC100747941',], normalized = TRUE)
r3 <- list(counts = as.numeric(LOC100747941), group = as.factor(coldata$Group))
r3 <- as.tibble(r3)
r3$group <- factor(r3$group, levels=c("R_CTL", "R_IMDA", "R_IMDB"))
r3 <- na.omit(r3)
g3 <- ggplot(r3, aes(group,counts)) + geom_boxplot(aes(fill = group))
g3 <- g3 + labs(x = "Pesticide exposure", y = "Normalized Counts ", title = "hexamerin")
g3 <- g3 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
g3 <- g3 + scale_fill_manual(values=c("rosybrown2", "darkorange4", "darkorange2"))
g3 <- g3 + scale_x_discrete(breaks=c("R_CTL", "R_IMDA", "R_IMDB"), labels=c("not exposed", "IMD-A", "IMD-B"))
g3 

LOC100747001 <- counts(ddsGroup['LOC100747001',], normalized = TRUE)
r4 <- list(counts = as.numeric(LOC100747001), group = as.factor(coldata$Group))
r4 <- as.tibble(r4)
r4$group <- factor(r4$group, levels=c("R_CTL", "R_IMDA", "R_IMDB"))
r4 <- na.omit(r4)
g4 <- ggplot(r4, aes(group,counts)) + geom_boxplot(aes(fill = group))
g4 <- g4 + labs(x = "Pesticide exposure", y = "Normalized Counts ", title = "LOC100747001")
g4 <- g4 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
g4 <- g4 + scale_fill_manual(values=c("rosybrown2", "darkorange4", "darkorange2"))
g4 <- g4 + scale_x_discrete(breaks=c("R_CTL", "R_IMDA", "R_IMDB"), labels=c("not exposed", "IMD-A", "IMD-B"))
g4 

fig_diet1 <- ggarrange(g1, g2, g3, g4,
                         ncol = 2, nrow = 2)
fig_diet1


#Diet Heather - diet 2
LOC100744789 <- counts(ddsGroup['LOC100744789',], normalized = TRUE)
h1 <- list(counts = as.numeric(LOC100744789), group = as.factor(coldata$Group))
h1 <- as.tibble(h1)
h1$group <- factor(h1$group, levels=c("H_CTL", "H_IMDA", "H_IMDB"))
h1 <- na.omit(h1)
p1 <- ggplot(h1, aes(group,counts)) + geom_boxplot(aes(fill = group))
p1 <- p1 + labs(x = "Pesticide exposure", y = "Normalized Counts ", title = "glucose dehydrogenase")
p1 <- p1 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p1 <- p1 + scale_fill_manual(values=c("rosybrown2", "darkorange4", "darkorange2"))
p1 <- p1 + scale_x_discrete(breaks=c("H_CTL", "H_IMDA", "H_IMDB"), labels=c("not exposed", "IMD-A", "IMD-B"))
p1 <- p1 + geom_signif(comparisons = list(c("H_CTL", "H_IMDB")),
            map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.5), tip_length = 0, vjust = 0.8)
p1

LOC100742002 <- counts(ddsGroup['LOC100742002',], normalized = TRUE)
h2 <- list(counts = as.numeric(LOC100742002), group = as.factor(coldata$Group))
h2 <- as.tibble(h2)
h2$group <- factor(h2$group, levels=c("H_CTL", "H_IMDA", "H_IMDB"))
h2 <- na.omit(h2)
p2 <- ggplot(h2, aes(group,counts)) + geom_boxplot(aes(fill = group))
p2 <- p2 + labs(x = "Pesticide exposure", y = "Normalized Counts ", title = "odorant binding protein")
p2 <- p2 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p2 <- p2 + scale_fill_manual(values=c("rosybrown2", "darkorange4", "darkorange2"))
p2 <- p2 + scale_x_discrete(breaks=c("H_CTL", "H_IMDA", "H_IMDB"), labels=c("not exposed", "IMD-A", "IMD-B"))
p2 <- p2 + geom_signif(comparisons = list(c("H_CTL", "H_IMDB")),
                       map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.5), tip_length = 0, vjust = 0.8)
p2 

LOC100747941 <- counts(ddsGroup['LOC100747941',], normalized = TRUE)
h3 <- list(counts = as.numeric(LOC100747941), group = as.factor(coldata$Group))
h3 <- as.tibble(h3)
h3$group <- factor(h3$group, levels=c("H_CTL", "H_IMDA", "H_IMDB"))
h3 <- na.omit(h3)
p3 <- ggplot(h3, aes(group,counts)) + geom_boxplot(aes(fill = group))
p3 <- p3 + labs(x = "Pesticide exposure", y = "Normalized Counts ", title = "hexamerin")
p3 <- p3 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p3 <- p3 + scale_fill_manual(values=c("rosybrown2", "darkorange4", "darkorange2"))
p3 <- p3 + scale_x_discrete(breaks=c("H_CTL", "H_IMDA", "H_IMDB"), labels=c("not exposed", "IMD-A", "IMD-B"))
p3 <- p3 + geom_signif(comparisons = list(c("H_CTL", "H_IMDB")),
                       map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.5), y_position = 2200, tip_length = 0, vjust = 0.8)
p3 <- p3 + geom_signif(comparisons = list(c("H_IMDA", "H_IMDB")),
                       map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.5), tip_length = 0, vjust = 0.8)
p3 

LOC100747001 <- counts(ddsGroup['LOC100747001',], normalized = TRUE)
h4 <- list(counts = as.numeric(LOC100747001), group = as.factor(coldata$Group))
h4 <- as.tibble(h4)
h4$group <- factor(h4$group, levels=c("H_CTL", "H_IMDA", "H_IMDB"))
h4 <- na.omit(h4)
p4 <- ggplot(h4, aes(group,counts)) + geom_boxplot(aes(fill = group))
p4 <- p4 + labs(x = "Pesticide exposure", y = "Normalized Counts ", title = "LOC100747001")
p4 <- p4 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
p4 <- p4 + scale_fill_manual(values=c("rosybrown2", "darkorange4", "darkorange2"))
p4 <- p4 + scale_x_discrete(breaks=c("H_CTL", "H_IMDA", "H_IMDB"), labels=c("not exposed", "IMD-A", "IMD-B"))
p4 <- p4 + geom_signif(comparisons = list(c("H_CTL", "H_IMDB")),
                       map_signif_level=c("***"=0.001,"**"=0.01,"*"=0.5), tip_length = 0, vjust = 0.8)
p4 

fig_diet2 <- ggarrange(p1, p2, p3, p4,
                      ncol = 2, nrow = 2)
fig_diet2


#Diet combined
LOC100744789 <- counts(ddsGroup['LOC100744789',], normalized = TRUE)
c1 <- list(counts = as.numeric(LOC100744789), group = as.factor(coldata$Group))
c1 <- as.tibble(c1)
c1$group <- factor(c1$group, levels=c("M_CTL", "M_IMDA", "M_IMDB"))
c1 <- na.omit(c1)
q1 <- ggplot(c1, aes(group,counts)) + geom_boxplot(aes(fill = group))
q1 <- q1 + labs(x = "Pesticide exposure", y = "Normalized Counts ", title = "glucose dehydrogenase")
q1 <- q1 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
q1 <- q1 + scale_fill_manual(values=c("rosybrown2", "darkorange4", "darkorange2"))
q1 <- q1 + scale_x_discrete(breaks=c("M_CTL", "M_IMDA", "M_IMDB"), labels=c("not exposed", "IMD-A", "IMD-B"))
q1 

LOC100742002 <- counts(ddsGroup['LOC100742002',], normalized = TRUE)
c2 <- list(counts = as.numeric(LOC100742002), group = as.factor(coldata$Group))
c2 <- as.tibble(c2)
c2$group <- factor(c2$group, levels=c("M_CTL", "M_IMDA", "M_IMDB"))
c2 <- na.omit(c2)
q2 <- ggplot(c2, aes(group,counts)) + geom_boxplot(aes(fill = group))
q2 <- q2 + labs(x = "Pesticide exposure", y = "Normalized Counts ", title = "odorant binding protein")
q2 <- q2 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
q2 <- q2 + scale_fill_manual(values=c("rosybrown2", "darkorange4", "darkorange2"))
q2 <- q2 + scale_x_discrete(breaks=c("M_CTL", "M_IMDA", "M_IMDB"), labels=c("not exposed", "IMD-A", "IMD-B"))
q2 

LOC100747941 <- counts(ddsGroup['LOC100747941',], normalized = TRUE)
c3 <- list(counts = as.numeric(LOC100747941), group = as.factor(coldata$Group))
c3 <- as.tibble(c3)
c3$group <- factor(c3$group, levels=c("M_CTL", "M_IMDA", "M_IMDB"))
c3 <- na.omit(c3)
q3 <- ggplot(c3, aes(group,counts)) + geom_boxplot(aes(fill = group))
q3 <- q3 + labs(x = "Pesticide exposure", y = "Normalized Counts ", title = "hexamerin")
q3 <- q3 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
q3 <- q3 + scale_fill_manual(values=c("rosybrown2", "darkorange4", "darkorange2"))
q3 <- q3 + scale_x_discrete(breaks=c("M_CTL", "M_IMDA", "M_IMDB"), labels=c("not exposed", "IMD-A", "IMD-B"))
q3 

LOC100747001 <- counts(ddsGroup['LOC100747001',], normalized = TRUE)
c4 <- list(counts = as.numeric(LOC100747001), group = as.factor(coldata$Group))
c4 <- as.tibble(c4)
c4$group <- factor(c4$group, levels=c("M_CTL", "M_IMDA", "M_IMDB"))
c4 <- na.omit(c4)
q4 <- ggplot(c4, aes(group,counts)) + geom_boxplot(aes(fill = group))
q4 <- q4 + labs(x = "Pesticide exposure", y = "Normalized Counts ", title = "LOC100747001")
q4 <- q4 + theme_bw(base_size = 14) + theme(legend.title = element_blank()) + theme(legend.position = "none")
q4 <- q4 + scale_fill_manual(values=c("rosybrown2", "darkorange4", "darkorange2"))
q4 <- q4 + scale_x_discrete(breaks=c("M_CTL", "M_IMDA", "M_IMDB"), labels=c("not exposed", "IMD-A", "IMD-B"))
q4 

fig_combined <- ggarrange(q1, q2, q3, q4,
                         ncol = 2, nrow = 2)
fig_combined



