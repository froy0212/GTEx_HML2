#Covariates for each tissue and each provirus
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
figures = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression/graphs"
Demographics = "/Users/huginn/Desktop/HML-2/GTEx_HML2_Expression/Demographics"
Hardy = paste(figures,"Demographics/Hardy",sep="/")
age = paste(figures,"Demographics/Age",sep="/")
sex = paste(figures,"Demographics/Sex",sep="/")
ischemic = paste(figures, "Demographics/Ischemic",sep="/")
POGTSG=paste(figures,"Demographics/POG_TSG",sep="/")
setwd(dir)
getwd()


#load in expression data and save
counts_filelist = list.files(path = TPM_counts,pattern = "*TPM_total.csv") #save data as a list
counts_filelist
counts_datalist = lapply(paste(TPM_counts,counts_filelist,sep="/"), read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
counts_datalist

#name each data frame within the list by the tissue specifier within the file name
counts_filelist_edited = as.character(strsplit(counts_filelist, "_TPM_total.csv"))
counts_filelist_edited
names(counts_datalist) = counts_filelist_edited
head(counts_datalist)

#add outlier function
is_outlier <- function(x) {
  return(x < quantile(x, 0.25) - 1.5 * IQR(x) | x > quantile(x, 0.75) + 1.5 * IQR(x))
}

#meta data
meta = read.csv(paste(metadata,"SubjID_Pheno.csv",sep="/"),header = TRUE)
meta_mod = meta[,c("SUBJID","AGE","SEX","TRISCHD","TRCHSTIND","DTHHRDY")]
head(meta_mod)

total_df_list = list()

#cancer related genes
#downloaded list of tumor supppressors and proto-oncogenes from COSMIC (https://cancer.sanger.ac.uk/cosmic/census?tier=1)
#read  in list of cancer related genes
COSMIC = read.csv(paste(metadata,"Census_allFri Aug 28 23_09_14 2020.csv",sep="/"))
TSG_POG = COSMIC[,c("Gene.Symbol","Role.in.Cancer")]

for (i in (1:length(counts_datalist))) {
  #Tissue = stringr::str_replace((names(counts_datalist)[[i]]), '_TPM_total.csv', '')
  Tissue = "Artery_Coronary"
  i = "Artery_Coronary"
  
  #read in counts data
  df <- counts_datalist[[i]]
  rownames(df) = df$X
  df$X = NULL
  df$difference = NULL
  head(df)
  
  #extract meta data that match sample IDs
  SRA_run_ID = read.delim(paste(metadata,paste(Tissue,"SraRunTable.txt",sep="_"),sep="/"),header=TRUE,sep="\t")
  SRA_run_ID_to_merge = SRA_run_ID[,c("Run","Sample_Name")]
  string = SRA_run_ID_to_merge$Sample_Name
  subjid = sub("^(.*?-.*?)-.*", "\\1", string)
  head(subjid)
  length(subjid)
  SRA_run_ID_to_merge$SUBJID = subjid
  meta = read.csv(paste(metadata,"SubjID_Pheno.csv",sep="/"),header = TRUE)
  meta_mod = meta[,c("SUBJID","AGE","SEX","DTHHRDY","TRISCHD","TRCHSTIND")]
  colnames(meta_mod) = c("SUBJID","age","sex","DTHHRDY","TRISCHD","TRCHSTIND")
  head(meta_mod)
  meta_updated = join(SRA_run_ID_to_merge,meta_mod,by="SUBJID",type="inner")
  head(meta_updated)
  dim(meta_updated)
  
  #filter expression data to the list of genes and proviruses Aidan is interested in
  Aidan = c("BPTF","SMARC","MAZ","PRa","PRb","Usf1","SETD2","ASH1L","TP53","SMYD2",
            "SP90","CCDC101","C-MYC","E2F","KAT5","KAT2A","KAT2B","KAT6A","AT2B",
            "PGR", "HNF1A", "HNF1C", "SRY", "TCF4", "LEF1",
            "HML-2_3q12.3", "HML-2_1q21.3", "HML-2_1q22", "HML-2_12q24.33",
            "HML-2_19p13.12b", "HML-2_17p13.1")
  Aidan_df = df[rownames(df) %in% Aidan,,drop =FALSE]
  
  #filter expression data to average TPM >= 1
  Aidan_df_filtered = Aidan_df[rowMeans(Aidan_df) >= 1,,drop=FALSE]
  
  #add in COV to reduce the X axis length
  Aidan_df_filtered$mean = rowMeans(Aidan_df_filtered)
  Aidan_df_filtered$stdev = rowSds(as.matrix(Aidan_df_filtered))
  Aidan_df_filtered = transform(Aidan_df_filtered, COV = stdev/mean)
  
  #select all rows where COV < 1
  Aidan_df_filtered_COV = Aidan_df_filtered[Aidan_df_filtered$COV < 1, ]
  
  Aidan_df_filtered_COV$COV = NULL
  Aidan_df_filtered_COV$stdev = NULL
  Aidan_df_filtered_COV$mean = NULL
  
  #scatterplot
  Aidan_df_filtered_COV = as.data.frame(t(Aidan_df_filtered_COV))
  
  #plot based on r2 linear regression
  dir.create(paste(dir,"correlation_plots",Tissue,sep="/"))
  par(mar=c(5,5,5,5))
  R2_Aidan_df = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(R2_Aidan_df) = c("gene","provirus","r2")
  R2_Aidan_slope_df = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(R2_Aidan_slope_df) = c("gene","provirus","slope")
  
  
  HML2_list = names(Aidan_df_filtered_COV)[grep("HML-2",names(Aidan_df_filtered_COV))]
  for (i in HML2_list) {
    print(i)
    for(j in colnames(Aidan_df_filtered_COV)){
      print(j)
      
      #data for regression line
      mod = lm(Aidan_df_filtered_COV[,j]~Aidan_df_filtered_COV[,i])
      modsum = summary(mod)
      r2 = modsum$adj.r.squared
      my.p = modsum$coefficients[2,4]

      #make data frame or r2
      tempdf = data.frame("gene" = j, "provirus" = i, "r2" = r2, "Tissue" = Tissue)
      
      R2_Aidan_df = rbind(R2_Aidan_df,tempdf)
      
      if (r2 > 0.5 && r2 < 1){
        #plot and decorate plot
        mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
        name = glue::glue("Correlation between {i} and {j} expression.jpeg")
        jpeg(filename=paste(dir,"correlation_plots",Tissue,name,sep="/"))
        plot(Aidan_df_filtered_COV[,i],
             Aidan_df_filtered_COV[,j],
             xlab=i,
             ylab=j,
             main = glue::glue("Correlation between {i} and {j} expression"),
             sub = mylabel)
        vec = par("usr")
        maxx = nth(vec, 2, order_by = NULL, default = default_missing(x))
        maxy = nth(vec, 4, order_by = NULL, default = default_missing(x))
        # mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
        # text(x = (maxx - 1.5), y = (maxy - 1.5), labels = mylabel)
        abline(mod, col="red")
        dev.off()
        
        #make data frame of slope
        tempSlopeDF_Aidan = data.frame("gene" = j, "provirus" = i, "slope" = mod$coefficients[2])
        rownames(tempSlopeDF_Aidan) = j
        
        R2_Aidan_slope_df = rbind(R2_Aidan_slope_df,tempSlopeDF_Aidan)
      }
    }
  }
  write.csv(R2_Aidan_df, paste(dir,"correlation_plots",Tissue,"r2_Aidan.csv",sep="/"))
  write.csv(R2_Aidan_slope_df, paste(dir,"correlation_plots",Tissue,"slope_Aidan.csv",sep="/"))
  
  #filter expression data to get HML-2 expression
  HML2_df = df[rownames(df) %like% "HML-2",,drop =FALSE]
  
  #filter expression data to average TPM >= 1
  HML2_filtered = HML2_df[rowMeans(HML2_df) >= 1,,drop=FALSE]
  
  #reshape for ggplot
  HML2_filtered$provirus = rownames(HML2_filtered)
  HML2_filtered$provirus = gsub("HML-2_","",HML2_filtered$provirus)
  HML2_mod = reshape2::melt(HML2_filtered)
  colnames(HML2_mod)=c("provirus","Run","TPM")
  HML2_meta = join(HML2_mod,meta_updated,by="Run")
  HML2_meta$Sample_Name = NULL
  HML2_meta$SUBJID = NULL
  HML2_meta$SEX<-ifelse(HML2_meta$sex==1,'Male','Female')
  
  #convert age to range
  HML2_meta$range = ifelse(HML2_meta$age >= 20 & HML2_meta$age <= 35, "20-35", ifelse(HML2_meta$age >=36 & HML2_meta$age <= 51,"36-51",ifelse(HML2_meta$age >= 52 & HML2_meta$age <= 70, "52-70", "NA")))
  
  #convert hardy score to characters
  HML2_meta$HRDYRange = ifelse(HML2_meta$DTHHRDY == 0 , "Zero", ifelse(HML2_meta$DTHHRDY == 1 ,"One", ifelse(HML2_meta$DTHHRDY == 2, "Two", ifelse(HML2_meta$DTHHRDY == 3, "Three","Four"))))
  
  #plot sex
  p<-ggplot(HML2_meta, aes(x=provirus, y=TPM, fill = SEX)) + geom_boxplot(coef = 6, outlier.shape=NA)+
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), binwidth = 0.5, stackratio = 0.025) + 
    labs(title=paste("Differences between biological sex in regards to HML-2 expression in", Tissue, sep = " "),x="Provirus", y = "HML2 TPM")
  p + stat_n_text(angle = 90, size=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(Tissue,"provirus expresssion divided by biological sex 08302020.png",sep ="  "))
  
  #plot age
  p<-ggplot(HML2_meta, aes(x=provirus, y=TPM, fill = range)) + geom_boxplot(coef = 6, outlier.shape=NA)+
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), binwidth = 0.5, stackratio = 0.025) + 
    labs(title=paste("Differences between age range in regards to HML-2 expression in", Tissue, sep = " "),x="Provirus", y = "HML2 TPM")
  p + stat_n_text(angle = 90, size=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(Tissue,"provirus expresssion divided by age range 08302020.png",sep ="  "))
  
  #plot Hardy
  p<-ggplot(HML2_meta, aes(x=provirus, y=TPM, fill = HRDYRange)) + geom_boxplot(coef = 6, outlier.shape=NA)+
    geom_dotplot(binaxis='y', stackdir='center', position=position_dodge(1), binwidth = 0.5, stackratio = 0.025) + 
    labs(title=paste("Differences between Hardy score in regards to HML-2 expression in", Tissue, sep = " "),x="Provirus", y = "HML2 TPM")
  p + stat_n_text(angle = 90, size=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  ggsave(paste(Tissue,"provirus expresssion divided by Hardy score 08302020.png",sep ="  "))
  
  #plot Ischemic time
  #calculate correlation with p value
  mod_df = HML2_meta[,c("provirus","Run","TPM")]
  mod_dfT=reshape2::dcast(mod_df,formula = Run~provirus, fun.aggregate = sum,value.var = "TPM")
  merge = unique(HML2_meta[,c("Run","TRISCHD","TRCHSTIND")])
  merge_final = transform(merge, TRISCHD = as.numeric(TRISCHD))
  merge_final = transform(merge_final, TRCHSTIND = as.numeric(TRCHSTIND))
  mod_dfT_merged = join(mod_dfT,merge_final,by="Run")
  rownames(mod_dfT_merged) = mod_dfT_merged$Run
  mod_dfT_merged$Run = NULL
  
  #to remove any NAs since rcorr won't work with them
  mod_dfT_merged[is.na(mod_dfT_merged)] <- 0
  
  #pearson correlation
  res = rcorr(as.matrix(mod_dfT_merged))
  write.csv(res$r,file = paste(Tissue,"pearson_HML2_TPM_ischemia_08302020.csv",sep="_"))
  
  #to replace any potential NA with 0 for the heatmap
  res$r[is.na(res$r)] <- 0
  
  #generate heatmap of correlation plot
  my_palette <- colorRampPalette(c("red", "black", "green"))
  mypath <- file.path(paste(dir,paste("Pearson_correlation_HML2VIschemia", Tissue, "08302020.png", sep = ""), sep ="/"))
  png(file=mypath,width = 400, height = 400, units='mm', res = 300)
  par(mar=c(1,1,1,1))
  heatmap.2(x = res$r, main = paste("Ischemic time in",Tissue,sep=" "),col=my_palette, labCol=as.expression(lapply(colnames(res$r), function(a) bquote(bold(.(a))))), labRow=as.expression(lapply(rownames(res$r), function(a) bquote(bold(.(a))))), dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',cexRow = 0.7, cexCol = 0.7)
  dev.off()
  
  #plot tumor suppresors
  #filter TSG_POG  to get tumor suppresssor  genes only 
  TSG_list = TSG_POG[TSG_POG$Role.in.Cancer %in% "TSG",]
  
  #extract expression data for HML2s and TSG and combine
  TSG_df = df[rownames(df) %in% TSG_list$Gene.Symbol,,drop =FALSE]
  HML2_df = df[rownames(df) %like% "HML-2",,drop =FALSE]
  TSG_HML2_df = rbind(TSG_df,HML2_df)
  rownames(TSG_HML2_df) = gsub("_new","", rownames(TSG_HML2_df))
  
  #filter expression data to average TPM >= 1
  TSG_HML2_filtered = TSG_HML2_df[rowMeans(TSG_HML2_df) >= 1,,drop=FALSE]
  
  #correlation and plot
  res = rcorr(as.matrix(t(TSG_HML2_filtered)))
  write.csv(res$r,file = paste(Tissue,"pearson_HML2_TSG_TPM_08302020.csv",sep="_"))
  
  #to replace any potential NA with 0 for the heatmap
  res$r[is.na(res$r)] <- 0
  
  #generate heatmap of correlation plot
  my_palette <- colorRampPalette(c("red", "black", "green"))
  mypath <- file.path(paste(dir,paste("Pearson_correlation_HML2VTSG", Tissue, "08302020.png", sep = ""), sep ="/"))
  png(file=mypath,width = 400, height = 400, units='mm', res = 300)
  par(mar=c(1,1,1,1))
  heatmap_TSG_df = res$r[grep("HML-2",rownames(res$r)),]
  rownames(heatmap_TSG_df) = gsub("HML-2_","", rownames(heatmap_TSG_df))
  colnames(heatmap_TSG_df) = gsub("HML-2_","", colnames(heatmap_TSG_df))
  heatmap.2(x = heatmap_TSG_df, main = paste("TSG in",Tissue,sep=" "),col=my_palette, labCol=as.expression(lapply(colnames(heatmap_TSG_df), function(a) bquote(bold(.(a))))), labRow=as.expression(lapply(rownames(heatmap_TSG_df), function(a) bquote(bold(.(a))))), dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',cexRow = 0.8, cexCol = 0.6,lhei = c(1,5),margins = c(7, 7),keysize=0.5,key.par = list(cex=0.5))
  dev.off()
  
  #plot proto-oncogenes
  POG_list = TSG_POG[TSG_POG$Role.in.Cancer %in% "oncogene",]
  
  #extract expression data for HML2s and POGs and combine
  POG_df = df[rownames(df) %in% POG_list$Gene.Symbol,,drop =FALSE]
  HML2_df = df[rownames(df) %like% "HML-2",,drop =FALSE]
  POG_HML2_df = rbind(POG_df,HML2_df)
  rownames(POG_HML2_df) = gsub("_new","", rownames(POG_HML2_df))
  
  #filter expression data to average TPM >= 1
  POG_HML2_filtered = POG_HML2_df[rowMeans(POG_HML2_df) >= 1,,drop=FALSE]
  
  #correlation and plot
  res = rcorr(as.matrix(t(POG_HML2_filtered)))
  write.csv(res$r,file = paste(Tissue,"pearson_HML2_POG_TPM_08302020.csv",sep="_"))
  
  #to replace any potential NA with 0 for the heatmap
  res$r[is.na(res$r)] <- 0
  
  #generate heatmap of correlation plot
  my_palette <- colorRampPalette(c("red", "black", "green"))
  mypath <- file.path(paste(dir,paste("Pearson_correlation_HML2VSPOG", Tissue, "08302020.png", sep = ""), sep ="/"))
  png(file=mypath,width = 400, height = 400, units='mm', res = 300)
  par(mar=c(1,1,1,1))
  heatmap_POG_df = res$r[grep("HML-2",rownames(res$r)),]
  rownames(heatmap_POG_df) = gsub("HML-2_","", rownames(heatmap_POG_df))
  colnames(heatmap_POG_df) = gsub("HML-2_","", colnames(heatmap_POG_df))
  heatmap.2(x = heatmap_POG_df, main = paste("POG in",Tissue,sep=" "),col=my_palette, labCol=as.expression(lapply(colnames(heatmap_POG_df), function(a) bquote(bold(.(a))))), labRow=as.expression(lapply(rownames(heatmap_POG_df), function(a) bquote(bold(.(a))))), dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',cexRow = 0.8, cexCol = 0.6,lhei = c(1,5),margins = c(7, 7),keysize=0.5,key.par = list(cex=0.5))
  dev.off()
  
  #add in COV to reduce the X axis length
  TSG_HML2_filtered$mean = rowMeans(TSG_HML2_filtered)
  TSG_HML2_filtered$stdev = rowSds(as.matrix(TSG_HML2_filtered))
  TSG_HML2_filtered = transform(TSG_HML2_filtered, COV = stdev/mean)
  
  POG_HML2_filtered$mean = rowMeans(POG_HML2_filtered)
  POG_HML2_filtered$stdev = rowSds(as.matrix(POG_HML2_filtered))
  POG_HML2_filtered = transform(POG_HML2_filtered, COV = stdev/mean)
  
  #select all rows where COV < 1
  TSG_HML2_filtered_COV = TSG_HML2_filtered[TSG_HML2_filtered$COV < 1, ]
  
  POG_HML2_filtered_COV = POG_HML2_filtered[POG_HML2_filtered$COV < 1, ]
  
  TSG_HML2_filtered_COV$COV = NULL
  TSG_HML2_filtered_COV$stdev = NULL
  TSG_HML2_filtered_COV$mean = NULL
  
  POG_HML2_filtered_COV$COV = NULL
  POG_HML2_filtered_COV$stdev = NULL
  POG_HML2_filtered_COV$mean = NULL
  
  #repeat correlation plot and graph
  res_TSG = rcorr(as.matrix(t(TSG_HML2_filtered_COV)))
  write.csv(res_TSG$r,file = paste(Tissue,"pearson_HML2_TSG_TPM_COV_09102020.csv",sep="_"))
  
  res_POG = rcorr(as.matrix(t(POG_HML2_filtered_COV)))
  write.csv(res_POG$r,file = paste(Tissue,"pearson_HML2_POG_TPM_COV_09102020.csv",sep="_"))
  
  #to replace any potential NA with 0 for the heatmap
  res_TSG$r[is.na(res_TSG$r)] <- 0
  
  res_POG$r[is.na(res_POG$r)] <- 0
  
  #generate heatmap of correlation plot
  heatmap_TSG_df = res_TSG$r[grep("HML-2",rownames(res_TSG$r)),]
  if(dim(heatmap_TSG_df)[1] > 0) {
    my_palette <- colorRampPalette(c("red", "black", "green"))
    mypath <- file.path(paste(dir,paste("Pearson_correlation_HML2VTSG_COV", Tissue, "09102020.png", sep = ""), sep ="/"))
    png(file=mypath,width = 400, height = 400, units='mm', res = 300)
    par(mar=c(1,1,1,1))
    rownames(heatmap_TSG_df) = gsub("HML-2_","", rownames(heatmap_TSG_df))
    colnames(heatmap_TSG_df) = gsub("HML-2_","", colnames(heatmap_TSG_df))
    heatmap.2(x = heatmap_TSG_df, main = paste("TSG in",Tissue,"with COV filter",sep=" "),col=my_palette, labCol=as.expression(lapply(colnames(heatmap_TSG_df), function(a) bquote(bold(.(a))))), labRow=as.expression(lapply(rownames(heatmap_TSG_df), function(a) bquote(bold(.(a))))), dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',cexRow = 0.8, cexCol = 0.8,lhei = c(1,5),margins = c(7, 7),keysize=0.5,key.par = list(cex=0.5))
    dev.off()
  }
  
  heatmap_POG_df = res_POG$r[grep("HML-2",rownames(res_POG$r)),]
  if(dim(heatmap_POG_df)[1] > 0) {
    my_palette <- colorRampPalette(c("red", "black", "green"))
    mypath <- file.path(paste(dir,paste("Pearson_correlation_HML2VPOG_COV", Tissue, "09102020.png", sep = ""), sep ="/"))
    png(file=mypath,width = 400, height = 400, units='mm', res = 300)
    par(mar=c(1,1,1,1))
    rownames(heatmap_POG_df) = gsub("HML-2_","", rownames(heatmap_POG_df))
    colnames(heatmap_POG_df) = gsub("HML-2_","", colnames(heatmap_POG_df))
    heatmap.2(x = heatmap_POG_df, main = paste("POG in",Tissue,"with COV filter",sep=" "),col=my_palette, labCol=as.expression(lapply(colnames(heatmap_POG_df), function(a) bquote(bold(.(a))))), labRow=as.expression(lapply(rownames(heatmap_POG_df), function(a) bquote(bold(.(a))))), dendrogram='none', Rowv=FALSE, Colv=FALSE,trace='none',cexRow = 0.8, cexCol = 0.8,lhei = c(1,5),margins = c(7, 7),keysize=0.5,key.par = list(cex=0.5))
    dev.off()
  }
  
  #scatterplot
  TSG_HML2_filtered_COV = as.data.frame(t(TSG_HML2_filtered_COV))
  POG_HML2_filtered_COV = as.data.frame(t(POG_HML2_filtered_COV))
  
  #plot based on r2 linear regression
  dir.create(paste(dir,"correlation_plots",Tissue,sep="/"))
  par(mar=c(5,5,5,5))
  R2_POG_df = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(R2_POG_df) = c("gene","provirus","r2")
  R2_TSG_df = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(R2_TSG_df) = c("gene","provirus","r2")
  R2_POG_slope_df = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(R2_POG_slope_df) = c("gene","provirus","slope")
  R2_TSG_slope_df = data.frame(matrix(ncol = 3, nrow = 0))
  colnames(R2_TSG_slope_df) = c("gene","provirus","slope")
  
  HML2_list = names(POG_HML2_filtered_COV)[grep("HML-2",names(POG_HML2_filtered_COV))]
  for (i in HML2_list) {
    print(i)
    for(j in colnames(POG_HML2_filtered_COV)){
      print(j)
      
      #data for regression line
      mod = lm(POG_HML2_filtered_COV[,j]~POG_HML2_filtered_COV[,i])
      modsum = summary(mod)
      r2 = modsum$adj.r.squared
      my.p = modsum$coefficients[2,4]
      
      #make data frame or r2
      tempdf = data.frame("gene" = j, "provirus" = i, "r2" = r2, "Tissue" = Tissue)
      
      R2_POG_df = rbind(R2_POG_df,tempdf)
      
      if (r2 > 0.5 && r2 < 1){
        #plot and decorate plot
        mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
        name = glue::glue("POG,Correlation between {i} and {j} expression.jpeg")
        jpeg(filename=paste(dir,"correlation_plots",Tissue,name,sep="/"))
        plot(POG_HML2_filtered_COV[,i],
             POG_HML2_filtered_COV[,j],
             xlab=i,
             ylab=j,
             main = glue::glue("Correlation between {i} and {j} expression"),
             sub = mylabel)
        vec = par("usr")
        maxx = nth(vec, 2, order_by = NULL, default = default_missing(x))
        maxy = nth(vec, 4, order_by = NULL, default = default_missing(x))
        # mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
        # text(x = (maxx - 1.5), y = (maxy - 1.5), labels = mylabel)
        abline(mod, col="red")
        dev.off()
        
        #make data frame of slope
        tempSlopeDF_POG = data.frame("gene" = j, "provirus" = i, "slope" = mod$coefficients[2])
        rownames(tempSlopeDF_POG) = j
        
        R2_POG_slope_df = rbind(R2_POG_slope_df,tempSlopeDF_POG)
        
      }
    }
    
    
    for(j in colnames(TSG_HML2_filtered_COV)){
      print(j)
      
      #data for regression line
      mod = lm(TSG_HML2_filtered_COV[,j]~TSG_HML2_filtered_COV[,i])
      modsum = summary(mod)
      r2 = modsum$adj.r.squared
      my.p = modsum$coefficients[2,4]
      
      #make data frame or r2
      tempdf_TSG = data.frame("gene" = j, "provirus" = i, "r2" = r2, "Tissue" = Tissue)
      
      R2_TSG_df = rbind(R2_TSG_df,tempdf_TSG)
      
      
      if(r2 > 0.5 && r2 < 1){
        #plot and decorate plot
        mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
        name = glue::glue("TSG,Correlation between {i} and {j} expression.jpeg")
        jpeg(filename=paste(dir,"correlation_plots",Tissue,name,sep="/"))
        plot(TSG_HML2_filtered_COV[,i],
             TSG_HML2_filtered_COV[,j],
             xlab=i,
             ylab=j,
             main = glue::glue("Correlation between {i} and {j} expression"),
             sub = mylabel)
        vec = par("usr")
        maxx = nth(vec, 2, order_by = NULL, default = default_missing(x))
        maxy = nth(vec, 4, order_by = NULL, default = default_missing(x))
        # mylabel = bquote(italic(R)^2 == .(format(r2, digits = 3)))
        # text(x = (maxx - 1.5), y = (maxy - 1.5), labels = mylabel)
        abline(mod, col="red")
        dev.off()
        
        #make data frame of slope
        tempSlopeDF_TSG = data.frame("gene" = j, "provirus" = i, "slope" = mod$coefficients[2])
        rownames(tempSlopeDF_TSG) = j
        
        R2_TSG_slope_df = rbind(R2_TSG_slope_df,tempSlopeDF_TSG)
        
      }
    }
    
  }
  R2_POG_df_reshaped = reshape2::dcast(R2_POG_df, gene ~ provirus, value.var = "r2")
  R2_TSG_df_reshaped = reshape2::dcast(R2_TSG_df, gene ~ provirus, value.var = "r2")
  # R2_POG_slope_df_reshaped = reshape2::dcast(R2_POG_slope_df, gene ~ provirus, value.var = "slope")
  # R2_TSG_slope_df_reshaped = reshape2::dcast(R2_TSG_slope_df, gene ~ provirus, value.var = "slope")
  
  write.csv(R2_TSG_df_reshaped, paste(dir,"correlation_plots",Tissue,"r2_TSG.csv",sep="/"))
  write.csv(R2_POG_df_reshaped, paste(dir,"correlation_plots",Tissue,"r2_POG.csv",sep="/"))
  write.csv(R2_TSG_slope_df, paste(dir,"correlation_plots",Tissue,"slope_TSG.csv",sep="/"))
  write.csv(R2_POG_slope_df, paste(dir,"correlation_plots",Tissue,"slope_POG.csv",sep="/"))
}
