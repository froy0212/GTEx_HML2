#check for a correlation between hypoxia-related genes and HML-2 expression
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
library(purrr)

#variables
# dir = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression"
# TPM_counts = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression/counts_matrices/TPM"
# figures = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression/graphs"
# metadata = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression/metadata"
# figures = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression/graphs"
# Demographics = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression/Demographics"
# Hardy = paste(figures,"Demographics/Hardy",sep="/")
dir = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/"
raw_counts = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/V15 Telescope/"
TPM_counts = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/V15 Counts/"
figures = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/Figures"
metadata = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8"
tissuemeta = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/Full GTEx_HML2_Expression 3.nosync"
Demographics = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8"
Hardy = paste(figures,"Demographics/Hardy",sep="/")
setwd(dir)
getwd()


#downloaded GO terms from GSEA. Grabbed any list with the word "hypoxia" in it. Ended with the following:
#https://www.gsea-msigdb.org/gsea/msigdb/genesets.jsp?collection=C5 - website
#GO_HYPOXIA_INDUCIBLE_FACTOR_1ALPHA_SIGNALING_PATHWAY
#GO_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY_IN_RESPONSE_TO_HYPOXIA
#GO_NEGATIVE_REGULATION_OF_CELLULAR_RESPONSE_TO_HYPOXIA
#GO_NEGATIVE_REGULATION_OF_HYPOXIA_INDUCED_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY
#GO_NEGATIVE_REGULATION_OF_HYPOXIA_INDUCED_INTRINSIC_APOPTOTIC_SIGNALING_PATHWAY
#GO_REGULATION_OF_CELLULAR_RESPONSE_TO_HYPOXIA
#GO_REGULATION_OF_TRANSCRIPTION_FROM_RNA_POLYMERASE_II_PROMOTER_IN_RESPONSE_TO_HYPOXIA

#load in the gmt files and combine them into one file
GO_hypoxia_list = list.files(pattern = "geneset*")
GO_hypoxia_list_1 = lapply(GO_hypoxia_list,read.gmt)
GO_hypoxia = unlist(GO_hypoxia_list_1, use.names = FALSE)
GO_hypoxia_unique = unique(GO_hypoxia)
GO_hypoxia_unique

#create a vector of HML-2 provirus names
pro_list <- read.delim("../../HML2list.txt", sep ="\t", header = FALSE)
HML2_proviruses = pro_list$V1
HML2_names = paste("HML-2",HML2_proviruses,sep="_")
#merge to create one vector
gene_names = c(GO_hypoxia_unique,HML2_names)
gene_names

#filter expression data to get hypoxia-related genes and HML-2 proviruses only
counts_filelist = list.files(path=TPM_counts,pattern = "*TPM_total.csv") #save data as a list
counts_filelist
#HML2_filelist <- HML2_filelist[-23]
counts_datalist = lapply(paste(path=TPM_counts, counts_filelist, sep = "/"), read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
counts_datalist

#name each data frame within the list by the tissue specifier within the file name
counts_filelist_edited = as.character(strsplit(counts_filelist, "_TPM_total.csv"))
counts_filelist_edited
names(counts_datalist) = counts_filelist_edited
head(counts_datalist)

for (i in (1:length(counts_datalist))) {
  Tissue = stringr::str_replace((names(counts_datalist)[[i]]), '_TPM_total.csv', '')
  Tissue
  #read in counts data
  df <- counts_datalist[[i]]
  rownames(df) = df$X
  df$X = NULL
  df$difference = NULL
  head(df)
  hml2_hypoxia_df = df[row.names(df) %in% gene_names,]
  write.csv(hml2_hypoxia_df,file = paste(Tissue,"hypoxia_HML2_TPM.csv",sep="_"))
}        

#run pearson correlation on hml2_hypoxia data frame (all v all) for each tissue
counts_filelist = list.files(pattern = "*hypoxia_HML2_TPM.csv") #save data as a list
counts_filelist
counts_datalist = lapply(counts_filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
counts_datalist

#name each data frame within the list by the tissue specifier within the file name
counts_filelist_edited = as.character(strsplit(counts_filelist, "_hypoxia_HML2_TPM.csv"))
counts_filelist_edited
names(counts_datalist) = counts_filelist_edited
head(counts_datalist)

for (i in (1:length(counts_datalist))) {
  Tissue = stringr::str_replace((names(counts_datalist)[[i]]), '_V15_hypoxia_HML2_TPM.csv', '')
  Tissue
  #read in counts data
  df <- counts_datalist[[i]]
  rownames(df) = df$X
  df$X = NULL
  df$difference = NULL
  head(df)
  
  #remove genes where average TPM <1
  df$average = rowMeans(df)
  df_filtered = subset(df,df$average > 1)
  head(df_filtered)
  df_filtered$average = NULL
  head(t(df_filtered))
  
  #calculate correlation with p value
  res = rcorr(as.matrix(t(df_filtered)))
  write.csv(res$r,file = paste(Tissue,"pearson_HML2_TPM.csv",sep="_"))
  
  #to replace any potential NA with 0 for the heatmap
  res$r[is.na(res$r)] <- 0
  
  #generate heatmap of correlation plot
  my_palette <- colorRampPalette(c("red", "black", "green"))
  mypath <- file.path(paste(dir,paste("Pearson_correlation_HML2VHypoxia", Tissue, ".png", sep = ""), sep ="/"))
  png(file=mypath,width = 400, height = 400, units='mm', res = 300)
  heatmap.2(x = res$r, col=my_palette, labCol=as.expression(lapply(colnames(res$r), function(a) bquote(bold(.(a))))), labRow=as.expression(lapply(rownames(res$r), function(a) bquote(bold(.(a))))), dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',cexRow = 0.8, cexCol = 0.8,lhei = c(1,9), margins = c(7, 7),keysize=0.5,key.par = list(cex=0.5))
  dev.off()
}
