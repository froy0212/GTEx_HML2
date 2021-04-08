#make heatmap of individual HML-2 GTEx expression data
#load libraries
library(ggplot2)
library(readr)
library(stringr)
library(plyr)
library(reshape)
library(reshape2)
library(pheatmap)

setwd("/Users/far122/Desktop/GTEx_HML2_Expression")
getwd()

wd = "/Users/far122/Desktop/GTEx_HML2_Expression"

#load in individual expression files for each tissue and assign to a list
filelist = list.files(pattern = "*HML2_individual_expression.csv") #save data as a list
filelist
datalist = lapply(filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
datalist
class(datalist)

#name each data frame within the list by the tissue specifier within the file name
filelist_edited = as.character(strsplit(filelist, "_HML2_individual_expression.csv"))
filelist_edited
names(datalist) = filelist_edited

#rename the columns in each data frame to help with merging later
colnames = c("provirus", "Average", "StDev")
datalist = lapply(datalist, setNames, colnames)
head(datalist)

#load in all known provirus names to later assign as row names
gene_names = read.table("all_proviruses.txt", header=FALSE)

#remove the HML-2_ and _new designation from the naming sceme
gene_names[] <- lapply(gene_names, gsub, pattern='HML-2_', replacement='')
gene_names[] <- lapply(gene_names, gsub, pattern='_new', replacement='')
gene_names

#name column to help with merging later
colnames(gene_names) = c("provirus")
head(gene_names)
dim(gene_names)

#make empty data frame and append final counts data via matching to the first column
#this step also removes the HML-2 and _new designation within the for loop
HML2_individual_Counts = as.data.frame(gene_names)
head(HML2_individual_Counts)
for (i in (1:length(datalist))) {
  print(i)
  df <- datalist[[i]]
  df$provirus <- stringr::str_replace(df$provirus, 'HML-2_', '')
  df$provirus <- stringr::str_replace(df$provirus, '_new', '')
  HML2_individual_Counts <- join(HML2_individual_Counts, df[,c(1,2)], by = "provirus", type = "left")
} 
head(HML2_individual_Counts)

#assign provirus names as the row names and remove the old column containing just provirus names
rownames(HML2_individual_Counts) = HML2_individual_Counts$provirus
HML2_individual_Counts$provirus = NULL
head(HML2_individual_Counts)

#convert all NA designations to 0
HML2_individual_Counts[is.na(HML2_individual_Counts)] <- 0

#change the column names to match with the corresponding tissue file names
colnames(HML2_individual_Counts) <- filelist_edited
head(HML2_individual_Counts)

#remove any rows where the sum is = 0
greater_than_0 <- as.data.frame(data.matrix(HML2_individual_Counts[rowSums(HML2_individual_Counts[,])>0, ]))
head(greater_than_0)

#save data 
write.csv(greater_than_0, file = paste('HML2_individual_expression', '.csv', sep="_"))

# plot a basic heatmap
paletteLength <- 50
myColor <- colorRampPalette(c("white", "blue"))(paletteLength)
map = pheatmap(greater_than_0, color = myColor, cluster_rows = FALSE, cluster_cols = FALSE, 
              main = "Average individual HML-2 expression in GTEx consortium", fontsize = 14, fontsize_row = 8, 
              fontsize_col = 8, breaks = seq(0,8,by=0.2))