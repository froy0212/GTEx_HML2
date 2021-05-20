#table showing the number of times a provirus is expressed in each tissue
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
dir = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression"
TPM_counts = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression/counts_matrices/TPM"
figures = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression/graphs"
metadata = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression/metadata"
setwd(dir)
getwd()


HML2_filelist = list.files(path = TPM_counts, pattern = "*_TPM_HML2.csv") #save data as a list
HML2_filelist
HML2_datalist = lapply(paste(TPM_counts, HML2_filelist,sep="/"), read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
HML2_datalist
HML2_filelist_edited = as.character(strsplit(HML2_filelist, "_TPM_HML2.csv"))
HML2_filelist_edited
names(HML2_datalist) = HML2_filelist_edited
head(HML2_datalist)
lapply(HML2_datalist,dim)

HML2_datalist_T = list()

for (i in (1:length(HML2_datalist))) {
  #read in hml-2 counts matrices 
  HML2_Tissue = stringr::str_replace((names(HML2_datalist)[[i]]), '_TPM_HML2.csv', '')
  HML2_Tissue
  HML2_df <- HML2_datalist[[i]]
  head(HML2_df)[1:5,1:5]
  HML2_DF = HML2_df[!grepl("GAPDH", HML2_df$X),]
  HML2_DF = HML2_DF[!grepl("ACTB", HML2_DF$X),]
  HML2_DF$X
  HML2_DF$difference = NULL
  head(HML2_DF)
  colnames(HML2_DF)
  
  #count the number of samples in each data frame
  HML2_DF$sample_size = ncol(HML2_DF[, -1])
  HML2_DF  = dplyr::select(HML2_DF,"sample_size",everything())
  
  #count the number of times a cell for a given provirus is >0 and save as a new column
  HML2_DF$count <- rowSums(HML2_DF[, -c(1,2)]!=0 )
  names(HML2_DF)[names(HML2_DF) == 'count'] <- paste(HML2_Tissue,'count',sep="_")
  names(HML2_DF)[names(HML2_DF) == 'sample_size'] <- paste(HML2_Tissue,'sample_size',sep="_")
  HML2_df_occurence = HML2_DF[,c("X",paste(HML2_Tissue,'count',sep="_"),paste(HML2_Tissue,'sample_size',sep="_"))]
  
  #add in number of total samples per tissue as a new row
  HML2_df_occ_samp = rbind(HML2_df_occurence[1,3],HML2_df_occurence)
  HML2_df_occ_samp[3] = NULL
  HML2_df_occ_samp[1,1] = "sample_size"
  head(HML2_df_occ_samp)
  HML2_df_occ_samp2= HML2_df_occ_samp %>% 
    mutate_at(vars(-X), funs(HML2_df_occ_samp[,2]/HML2_df_occ_samp[1,2]))
  head(HML2_df_occ_samp2)
  HML2_df_occ_samp3 = HML2_df_occ_samp2 %>% 
    mutate_at(vars(-X), funs(HML2_df_occ_samp2[,2]*100))
  colnames(HML2_df_occ_samp3)=c("X",paste(HML2_Tissue,'percent',sep="_"))
  HML2_df_occ_samp2_merged = join(HML2_df_occ_samp,HML2_df_occ_samp3,by="X")
  
  #add to empty list
  HML2_datalist_T[[HML2_Tissue]] <- HML2_df_occ_samp2_merged # add it to your list
}

#save table of provirus names(rownames) and tissue (colnames) with occurrence counts in each col
HML2_datalist_T
big_data = purrr::reduce(HML2_datalist_T, full_join, by = "X")
head(big_data)
rownames(big_data) = big_data$X
big_data$X = NULL
write.csv(big_data,"Provirus_occurrences_within_tissues_07282020.csv")
big_data_percent = big_data[,colnames(big_data) %like% "percent"]
write.csv(big_data_percent,"Provirus_occurrences_percent_08142020.csv")
