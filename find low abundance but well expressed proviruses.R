#find low abundance but well expressed proviruses
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
heatmap = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression/heatmap"
setwd(dir)
getwd()



#read in filtered dataframe
for_VLP_analysis_final = read.csv(paste(heatmap, "average individual HML-2 expression in GTEx, average >= 1.csv",sep="/"))
rownames(for_VLP_analysis_final) = for_VLP_analysis_final$X
for_VLP_analysis_final$X = NULL

#assign list of known proviruses to variable for removal
toRemoveForLowAbund = rownames(for_VLP_analysis_final)

#read in expression data for each tissue and save as list
filelist = list.files(path = TPM_counts, pattern = "*_TPM_HML2.csv") #save data as a list
filelist
datalist = lapply(paste(TPM_counts, filelist, sep ="/"), read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
datalist

#name each data frame within the list by the tissue specifier within the file name
filelist_edited = as.character(strsplit(filelist, "_TPM_HML2.csv"))
filelist_edited
names(datalist) = filelist_edited

HML2_datalist_LowAbund = list()

#identify proviruses that are expressed in low abundance across donors for each tissue which are also well expressed in a few donors/tissue
for (i in (1:length(datalist))) {
  #read in df
  HML2_Tissue = stringr::str_replace((names(datalist)[[i]]), '_TPM_HML2.csv', '')
  HML2_Tissue
  HML2_df <- datalist[[i]]
  head(HML2_df)[1:5,1:5]
  HML2_df$X
  HML2_df$difference = NULL
  head(HML2_df)
  
  #remove rows with all 0s
  HML2_df_remove0rows = HML2_df[apply(HML2_df[,-1], 1, function(x) !all(x==0)),]
  
  #get average for all non-0 entries
  is.na(HML2_df_remove0rows) <- HML2_df_remove0rows==0
  HML2_df_remove0rows$Average = rowMeans(HML2_df_remove0rows[,c(-1)], na.rm=TRUE)
  HML2_df_remove0rows$StDev = rowSds(as.matrix(HML2_df_remove0rows[,c(-1)]), na.rm=TRUE)
  HML2_DF_AvgStDev <- HML2_df_remove0rows[, c("X", "Average", "StDev")]
  colnames(HML2_DF_AvgStDev) = c("provirus", "Average", "StDev")
  head(HML2_DF_AvgStDev)
  
  #assign provirus names as the row names and remove the old column containing just provirus names
  rownames(HML2_DF_AvgStDev) = HML2_DF_AvgStDev$provirus
  HML2_DF_AvgStDev$provirus = NULL
  head(HML2_DF_AvgStDev)
  
  HML2_datalist_LowAbund[[HML2_Tissue]] <- HML2_DF_AvgStDev # add it to your list
}
HML2_datalist_LowAbund
big_data = do.call(rbind, HML2_datalist_LowAbund)
head(big_data)
big_data$test = rownames(big_data)
big_data$provirus = gsub("^.*?\\.","", big_data$test)
big_data$tissue = gsub("[.][\\s\\S]*$", "", big_data$test, perl = T)
rownames(big_data)  = c()
big_data$test = NULL
head(big_data)
big_data_transf = big_data[,c(3,1,2,4)]
head(big_data_transf)
big_data_transf_mod = big_data_transf[big_data_transf$Average >= 1, ]
for_VLP_analysis_LowAbund = reshape2::dcast(data = big_data_transf_mod,formula = provirus~tissue,fun.aggregate = sum,value.var = "Average")

#subtract known proviruses from the low abundance df
LowAbund_df = for_VLP_analysis_LowAbund[!for_VLP_analysis_LowAbund$provirus %in% toRemoveForLowAbund,]
rownames(LowAbund_df) = LowAbund_df$provirus
remove = c("ACTB","GAPDH")
LowAbund_df_final=LowAbund_df[!rownames(LowAbund_df) %in% remove,]

#save for records
write.csv(LowAbund_df_final,file=paste(heatmap,"HML2_individual_expression_LowAbund_08212020.csv",sep="/"))

#plot heatmap
order = c("1p31.1a", "1p31.1b", "1p34.3", "1p36.21a", "1p36.21c", "1q21.3", "1q22", "1q23.3", "1q24.1", "1q32.2", "1q43", "2q21.1", "3p12.3", "3p25.3", "3q12.3", "3q13.2", "3q21.2", "3q24", "3q27.2", "4p16.1a", "4p16.1b", "4p16.3a", "4p16.3b", "4q13.2", "4q32.1", "4q32.3", "4q35.2", "5p12", "5p13.3", "5q33.2", "5q33.3", "6p11.2", "6p21.1", "6p22.1", "6q14.1", "6q25.1", "7p22.1a", "7p22.1b", "7q11.21", "7q22.2", "7q34", "8p22", "8p23.1a", "8p23.1b", "8p23.1c", "8p23.1d", "8q11.1", "8q24.3a", "8q24.3b", "8q24.3c", "9q34.11", "9q34.3", "10p12.1", "10p14", "10q24.2", "11p15.4", "11q12.1", "11q12.3", "11q22.1", "11q23.3", "12p11.1", "12q13.2", "12q14.1", "12q24.11", "12q24.33", "14q11.2", "14q32.33", "15q25.2", "16p11.2", "16p13.3", "17p13.1", "19p12a", "19p12b", "19p12c", "19p12d", "19p12e", "19p13.3", "19q11", "19q13.12b", "19q13.41", "19q13.42", "20q11.22", "21q21.1", "22q11.21", "22q11.23", "Xq12", "Xq21.33", "Xq28a", "Xq28b", "Yp11.2", "Yq11.23a", "Yq11.23b")
paletteLength <- 50
myColor <- colorRampPalette(c("white", "blue"))(paletteLength)

LowAbund_df_final_mod = LowAbund_df_final %>%
  dplyr::slice(match(order,provirus))
rownames(LowAbund_df_final_mod)=LowAbund_df_final_mod$provirus
LowAbund_df_final_mod$provirus = NULL
LowAbund_df_plot = LowAbund_df_final_mod[!row.names(LowAbund_df_final_mod)%in%remove,]

map = pheatmap(as.matrix(LowAbund_df_plot), color = myColor, cluster_rows = FALSE, cluster_cols = FALSE, 
               main = "Average individual HML-2 expression in GTEx,heatmap",breaks = seq(0,2,by=0.2),fontsize = 8, fontsize_row = 6, 
               fontsize_col = 6)

#since there is only one provirus for this heatmap, I'll add it to the above heatmap that we're using in the paper
final_heatmap_df = rbind(for_VLP_analysis_final,LowAbund_df_plot)
write.csv(final_heatmap_df,paste(heatmap,"HML2_final_df_heatmap_08212020.csv",sep="/"))

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.5, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

png(paste(figures,"Average individual HML-2 expression in GTEx,low abundance well expressed.png",sep="/"))
map = pheatmap(as.matrix(final_heatmap_df), color = myColor, cluster_rows = FALSE, cluster_cols = FALSE, cellheight=10,cellwidth=12,
               main = "Average individual HML-2 expression in GTEx",breaks = seq(0,8,by=0.2),fontsize = 8, fontsize_row = 6, 
               fontsize_col = 6)
dev.off()
