#plot box plot for each  tissue showing HML-2 expression
#options(warn=2)
#options(warn=0, error=NULL)
library(data.table)
library(Biobase)
library(edgeR)
library(org.Hs.eg.db)
library(DESeq2)
library(GenomicAlignments)
library(GenomicRanges)
library(Rsamtools)
library(refGenome)
library(plyr)
library(stringr)
library(ggplot2)
library(matrixStats)
library(dplyr)
library(pheatmap)
library(biomaRt)
library("affy")
library(affycoretools)
library(EnvStats)
library(tidyverse)
library(qusage)
library("Hmisc")
library(gplots)
library("ggpubr")
library(gridExtra)
library(grid)
library(lattice)

#variables
dir = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression"
TPM_counts = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression/counts_matrices/TPM"
figures = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression/graphs"
boxplot = paste(figures,"individual_tissue_boxplot",sep="/")
setwd(dir)
getwd()

HML2_filelist = list.files(path = TPM_counts, pattern = "*_TPM_HML2.csv") #save data as a list
HML2_filelist
HML2_datalist = lapply(paste(TPM_counts, HML2_filelist, sep ="/"), read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
HML2_datalist
HML2_filelist_edited = as.character(strsplit(HML2_filelist, "_TPM_HML2.csv"))
HML2_filelist_edited
names(HML2_datalist) = HML2_filelist_edited
head(HML2_datalist)
lapply(HML2_datalist,dim)

for (i in (1:length(HML2_datalist))) {
  #read in data for each tissue
  HML2_Tissue = stringr::str_replace((names(HML2_datalist)[[i]]), '_TPM_HML2.csv', '')
  HML2_Tissue
  HML2_df <- HML2_datalist[[i]]
  head(HML2_df)[1:5,1:5]
  HML2_DF = HML2_df[!grepl("GAPDH", HML2_df$X),]
  HML2_DF = HML2_DF[!grepl("ACTB", HML2_DF$X),]
  HML2_DF$X
  HML2_DF$difference = NULL
  colnames(HML2_DF)
  rownames(HML2_DF)=HML2_DF$X
  HML2_DF$X = NULL
  head(HML2_DF)[1:5,1:5]
  
  #remove rows with all 0s
  HML2_DF_edited  = HML2_DF[rowSums(HML2_DF > 0) != 0, ]
  
  #reshape data frame for ggplot
  HML2_melt = reshape2::melt(as.matrix(HML2_DF))
  colnames(HML2_melt) = c("provirus","SRR_ID","TPM")
  
  #plot
  p<-ggplot(HML2_melt, aes(x=provirus, y=TPM)) + geom_boxplot()+
    geom_dotplot(binaxis='y', stackdir='center', binwidth = 0.1, stackratio = 0.025) + 
    labs(title=paste("HML-2 expression within", HML2_Tissue,"from GTEx",sep=" "),x="Provirus", y = "TPM")
  p + theme(axis.text.x = element_text(size = 7, angle = 90, hjust = 1))
  
  ggsave(paste(boxplot,paste(HML2_Tissue,"provirus_expresssion.png",sep ="_"),sep ="/"),width=11,height=8,units="in")
}

