#normalize and graph HML-2 expression from GTEx
##WHAT IS GOING ON WITH THE PROSTATE SAMPLES
##WHY DOES LOG2 +1 REMOVE THIS DIFFERENCE
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
library(affy)
library(affycoretools)
library(EnvStats)
library(tidyverse)

#variables
dir = "/Users/aidanburn/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8"
setwd(dir)
getwd()


#TPM normalize expression data and extract HML-2, GAPDH, and ACTB counts
  #load in individual expression files for each tissue and assign to a list
  filelist = list.files(pattern = "*V15_Telescope_output.csv") #save data as a list
  filelist
  datalist = lapply(filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)

  #name each data frame within the list by the tissue specifier within the file name
  filelist_edited = as.character(strsplit(filelist, "_V15_Telescope_output.csv"))
  filelist_edited
  names(datalist) = filelist_edited
  head(datalist)
  lapply(datalist,dim)
  
  #confirm runs and tissues match
  #for (i in (1:length(datalist))) {
  #  Tissue = stringr::str_replace((names(datalist)[[i]]), '_Telescope_output.csv', '')
  #  Tissue
  #  df <- datalist[[i]]
  #  rownames(df)=df$X
  #  df$X=NULL
  #  head(df)
  #  #load in SRA table for each tissue
  #  SRARunTable = read.delim(paste(Tissue,"SraRunTable.txt",sep="_"),header=TRUE,sep="\t")
  #  head(SRARunTable)
    
  #  table(SRARunTable$body_site)
  #  SRA_cofirmation_test = SRARunTable[c("Run","body_site")]
  #  head(SRA_cofirmation_test)
  #  dim(SRA_cofirmation_test)
  #  test = as.character(SRA_cofirmation_test$Run)
  #  stopifnot(all(names(colnames(df) %in% names(test))))
  #  }
  #doesn't stop, so they must all match. will proceed.
  
  #convert raw counts to TPM for every data frame in the list
  for (i in (1:length(datalist))) {
    Donor = stringr::str_replace((names(datalist)[[i]]), 'V_15_Telescope_output.csv', '')
    Donor
    df <- datalist[[i]]
    head(df)
    
    #convert to TPM and filter poorly expressed genes (TPM < 1 in at least half of the samples)
    #read in hg38.gtf for transcript coordinates
    ens = ensemblGenome()
    read.gtf(ens, "../hg38.gtf")
    class(ens)
    Identifier = "exon"
    hg38_Annotations = extractFeature(ens, Identifier)
    hg38_Annotations
    hg38_Annotation_df = data.frame(start=getGtf(hg38_Annotations)$start,end=getGtf(hg38_Annotations)$end,gene_id=getGtf(hg38_Annotations)$gene_id, transcript_id=getGtf(hg38_Annotations)$transcript_id)
    hg38_Annotation_df[1:5,]
    
    #for each row, get the difference between start and end
    hg38_Annotation_df$difference = hg38_Annotation_df$end-hg38_Annotation_df$start
    hg38_Annotation_df[1:5,]
    
    #for each gene_id, sum the difference to get transcript length. This also converts length to kb.
    #hg38_Annotation_df_lengthSum = ddply(hg38_Annotation_df, .(transcript_id), summarise, difference=(sum(difference)/1000)) #FENRIR
    hg38_Annotation_df_lengthSum = ddply(hg38_Annotation_df, .(gene_id), summarise, difference=(sum(difference)/1000)) #with Telescope
    hg38_Annotation_df_lengthSum[1:5,]
    
    # Divide the read counts by the length of each gene (transcript?) in kilobases. This gives you reads per kilobase (RPK).
    #stuff = merge(counts, hg38_Annotation_df_lengthSum, by.x = "X", by.y = "transcript_id") #FENRIR
    stuff = merge(df, hg38_Annotation_df_lengthSum, by.x = "X", by.y = "gene_id") #with Telescoppe
    head(stuff)
    tail(stuff)
    stuff_mod = stuff[,-1]
    rownames(stuff_mod) = stuff[,1]
    head(stuff_mod)
    tail(stuff_mod)
    RPK = stuff_mod/(stuff_mod$difference)
    head(RPK)
    
    # Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
    Per_Million=(colSums(RPK))/1e6
    Per_Million
    
    # Divide the RPK values by the “per million” scaling factor. This gives you TPM.
    Counts_HML2_TPM <- t(t(RPK)/Per_Million)
    #head(Counts_HML2_TPM)[1:5,1:5]
    
    #filter counts file to extract HML-2 and some house keeping genes (GAPDH, ACTB)
    HML2_counts = Counts_HML2_TPM[grepl("^HML-2", rownames(Counts_HML2_TPM)),]
    head(HML2_counts)
    length(rownames(HML2_counts))
    rownames(HML2_counts) <- stringr::str_replace(rownames(HML2_counts), 'HML-2_', '')
    rownames(HML2_counts) <- stringr::str_replace(rownames(HML2_counts), '_new', '')
    head(HML2_counts)
    GAPDH = Counts_HML2_TPM[rownames(Counts_HML2_TPM)=="GAPDH",]
    head(GAPDH)
    ACTB = Counts_HML2_TPM[rownames(Counts_HML2_TPM)=="ACTB",]
    head(ACTB)
    HML2_HKG_counts = rbind(HML2_counts,GAPDH,ACTB)
    head(HML2_HKG_counts)
    rownames(HML2_HKG_counts)
    
    write.csv(HML2_HKG_counts, file=paste(Donor, "TPM_HML2.csv", sep="_"))
    
    #log2(counts+1) transform 
    Counts_HML2_TPM_Filtered_Log2Trans = log2(HML2_HKG_counts + 1)
    head(Counts_HML2_TPM_Filtered_Log2Trans)
    rownames(Counts_HML2_TPM_Filtered_Log2Trans)
    write.csv(Counts_HML2_TPM_Filtered_Log2Trans, file=paste(Donor, "TPM_Log2_HML2.csv", sep="_"))
  }

#graph HML-2 expression per tissue in a box plot
  HML2_filelist = list.files(pattern = "*_TPM_HML2.csv") #save data as a list
  HML2_filelist
  HML2_datalist = lapply(HML2_filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
  HML2_datalist
  HML2_filelist_edited = as.character(strsplit(HML2_filelist, "_TPM_HML2.csv"))
  HML2_filelist_edited
  names(HML2_datalist) = HML2_filelist_edited
  head(HML2_datalist)
  lapply(HML2_datalist,dim)
  
  HML2_datalist_T = list()
  
  for (i in (1:length(HML2_datalist))) {
    HML2_Tissue = stringr::str_replace((names(HML2_datalist)[[i]]), '_TPM_HML2.csv', '')
    HML2_Tissue
    HML2_df <- HML2_datalist[[i]]
    head(HML2_df)[1:5,1:5]
    HML2_DF = HML2_df[!grepl("GAPDH", HML2_df$X),]
    HML2_DF = HML2_DF[!grepl("ACTB", HML2_DF$X),]
    HML2_DF = HML2_DF[!grepl("17p13.1", HML2_DF$X),]
    HML2_DF = HML2_DF[!grepl("8q24.3c", HML2_DF$X),]
    HML2_DF = HML2_DF[!grepl("21q21.1", HML2_DF$X), ]
    HML2_DF$X
    HML2_DF$difference = NULL
    head(HML2_DF)
    colnames(HML2_DF)
    HML2_DF_SampleSum = rbind(HML2_DF, data.frame(X = "HML2_Sum", t(colSums(HML2_DF[, -1]))))
    HML2_DF_SampleSum$X
    HML2_DF_SampleSum[HML2_DF_SampleSum$X %in% "HML2_Sum", ]
    #grab sum row and colnames
    HML2_SampleSum = HML2_DF_SampleSum[grepl("HML2_Sum", HML2_DF_SampleSum$X),]
    rownames(HML2_SampleSum) = HML2_SampleSum$X
    HML2_SampleSum$X = NULL
    head(HML2_SampleSum)
    #transpose and add new column with tissue id
    ERV_Counts_modT = as.data.frame(t(HML2_SampleSum))
    ERV_Counts_modT$tissue = HML2_Tissue
    head(ERV_Counts_modT)
    HML2_datalist_T[[HML2_Tissue]] <- ERV_Counts_modT # add it to your list
  }
  HML2_datalist_T
  lapply(HML2_datalist_T,dim)
  
  big_data = do.call(rbind, HML2_datalist_T)
  head(big_data)
  big_data = cbind(big_data, read.table(text=row.names(big_data), sep="_", 
                       header=FALSE, col.names = paste0("col", 1:2), stringsAsFactors=FALSE))
  big_data$col1 = NULL
  rownames(big_data)  = c()
  colnames(big_data) = c("HML2_TPM", "tissue")
  head(big_data)
  test = big_data[c("HML2_TPM", "tissue")]
  head(test)
  table(test$tissue)
  write.csv(table(test$tissue), "Tissue_Distribution_V15_provirustrim.csv")
  test$tissue <- as.factor(test$tissue)
  head(test)
  
  # ggplot code
  p<-ggplot(test, aes(x=tissue, y=HML2_TPM)) + geom_boxplot()+
    geom_dotplot(binaxis='y', stackdir='center', binwidth = 1, stackratio = 0.025) + 
    labs(title="Total HML-2 expression in GTEx",x="Tissue", y = "HML2 TPM (sum/sample)")
  p + stat_n_text(angle = 90, size=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
  #filter TPM expression to >= 1 and graph
  TPM_filter = test[test$HML2_TPM >= 1, ]
  
  pdf("Total HML-2 expression in GTEx, TPM per sample per tissue >= 1.pdf", width = 10)
  p<-ggplot(TPM_filter, aes(x=tissue, y=HML2_TPM)) + geom_boxplot()+
    geom_dotplot(binaxis='y', stackdir='center', binwidth = 1, stackratio = 0.025) + 
    labs(title="Total HML-2 expression in GTEx, TPM/sample/tissue >= 1",x="Tissue", y = "HML2 TPM (sum/sample)")
  p + stat_n_text(angle = 90, size=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  dev.off()
  
  ggsave(filename = "Total HML-2 expression in GTEx, TPM per sample per tissue >= 1.png", width = 7.5, height = 7)
  
  #what about log2?
  HML2_filelist = list.files(pattern = "*_TPM_Log2_HML2.csv") #save data as a list
  HML2_filelist
  HML2_datalist = lapply(HML2_filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
  HML2_datalist
  HML2_filelist_edited = as.character(strsplit(HML2_filelist, "_TPM_Log2_HML2.csv"))
  HML2_filelist_edited
  names(HML2_datalist) = HML2_filelist_edited
  head(HML2_datalist)
  
  HML2_datalist_T = list()
  
  for (i in (1:length(HML2_datalist))) {
    HML2_Tissue = stringr::str_replace((names(HML2_datalist)[[i]]), '_TPM_Log2_HML2.csv', '')
    HML2_Tissue
    HML2_df <- HML2_datalist[[i]]
    head(HML2_df)[1:5,1:5]
    HML2_DF = HML2_df[!grepl("GAPDH", HML2_df$X),]
    HML2_DF = HML2_DF[!grepl("ACTB", HML2_DF$X),]
    HML2_DF = HML2_DF[!grepl("17p13.1", HML2_DF$X),]
    HML2_DF = HML2_DF[!grepl("8q24.3c", HML2_DF$X),]
    HML2_DF = HML2_DF[!grepl("21q21.1", HML2_DF$X), ]
    HML2_DF$X
    HML2_DF$difference = NULL
    head(HML2_DF)
    colnames(HML2_DF)
    HML2_DF_SampleSum = rbind(HML2_DF, data.frame(X = "HML2_Sum", t(colSums(HML2_DF[, -1]))))
    HML2_DF_SampleSum$X
    HML2_DF_SampleSum[HML2_DF_SampleSum$X %in% "HML2_Sum", ]
    #grab sum row and colnames
    HML2_SampleSum = HML2_DF_SampleSum[grepl("HML2_Sum", HML2_DF_SampleSum$X),]
    rownames(HML2_SampleSum) = HML2_SampleSum$X
    HML2_SampleSum$X = NULL
    head(HML2_SampleSum)
    #transpose and add new column with tissue id
    ERV_Counts_modT = as.data.frame(t(HML2_SampleSum))
    ERV_Counts_modT$tissue = HML2_Tissue
    head(ERV_Counts_modT)
    HML2_datalist_T[[HML2_Tissue]] <- ERV_Counts_modT # add it to your list
  }
  HML2_datalist_T
  big_data = do.call(rbind, HML2_datalist_T)
  head(big_data)
  big_data = cbind(big_data, read.table(text=row.names(big_data), sep=".", 
                                        header=FALSE, col.names = paste0("col", 1:2), stringsAsFactors=FALSE))
  big_data$col1 = NULL
  rownames(big_data)  = c()
  colnames(big_data) = c("HML2_TPM", "tissue")
  head(big_data)
  test = big_data[c("HML2_TPM", "tissue")]
  head(test)
  
  test$tissue <- as.factor(test$tissue)
  head(test)
  
  # ggplot code
  p<-ggplot(test, aes(x=tissue, y=HML2_TPM)) + geom_boxplot()+
    geom_dotplot(binaxis='y', stackdir='center', binwidth = 1, stackratio = 0.025) + 
    labs(title="Total HML-2 expression in GTEx",x="Tissue", y = "HML2 TPM Log2 (sum/sample)")
  p + stat_n_text(angle = 90, size=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1))
  
#graph heatmap of individual HML-2 proviruses
  #load in individual expression files for each tissue and assign to a list
  filelist = list.files(pattern = "*_TPM_HML2.csv") #save data as a list
  filelist
  datalist = lapply(filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
  datalist

  #name each data frame within the list by the tissue specifier within the file name
  filelist_edited = as.character(strsplit(filelist, "_TPM_HML2.csv"))
  filelist_edited
  names(datalist) = filelist_edited
  
  HML2_datalist_AvgStDev = list()
  
  #calculate average and stdev for each provirus
  for (i in (1:length(datalist))) {
    HML2_Tissue = stringr::str_replace((names(HML2_datalist)[[i]]), '_TPM_HML2.csv', '')
    HML2_Tissue
    HML2_df <- HML2_datalist[[i]]
    head(HML2_df)[1:5,1:5]
    HML2_df$X
    HML2_df$difference = NULL
    head(HML2_df)
    
    HML2_df$Average = rowMeans(HML2_df[,c(-1)])
    HML2_df$StDev = rowSds(as.matrix(HML2_df[,c(-1)]))
    HML2_DF_AvgStDev <- HML2_df[, c("X", "Average", "StDev")]
    colnames(HML2_DF_AvgStDev) = c("provirus", "Average", "StDev")
    head(HML2_DF_AvgStDev)

    #assign provirus names as the row names and remove the old column containing just provirus names
    rownames(HML2_DF_AvgStDev) = HML2_DF_AvgStDev$provirus
    HML2_DF_AvgStDev$provirus = NULL
    head(HML2_DF_AvgStDev)
    
    HML2_datalist_AvgStDev[[HML2_Tissue]] <- HML2_DF_AvgStDev # add it to your list
  }
  HML2_datalist_AvgStDev
  big_data = do.call(rbind, HML2_datalist_AvgStDev)
  head(big_data)
  big_data$test = rownames(big_data)
  big_data$provirus = gsub("^.*?\\.","", big_data$test)
  big_data$tissue = gsub("[.][\\s\\S]*$", "", big_data$test, perl = T)
  rownames(big_data)  = c()
  big_data$test = NULL
  head(big_data)
  big_data_transf = big_data[,c(3,1,2,4)]
  head(big_data_transf)
  
  #save data 
  write.csv(big_data_transf, file = 'HML2_individual_expression_06152020.csv')
  
  #get data for heatmap
  HML2_heatmap = big_data_transf[, c("provirus", "Average", "tissue")]
  head(HML2_heatmap)
  HML2_heatmap_1 = reshape2::dcast(HML2_heatmap, provirus ~ tissue, value.var = "Average", fun.aggregate = mean)
  rownames(HML2_heatmap_1) = HML2_heatmap_1$provirus
  HML2_heatmap_1$provirus = NULL
  rownames(HML2_heatmap_1)
  HML2_heatmap_2 = HML2_heatmap_1[c("ACTB","GAPDH"),]
  head(HML2_heatmap_2)
  rownames(HML2_heatmap_2)
  remove = c("ACTB","GAPDH")
  HML2_heatmap_3 = HML2_heatmap_1[!row.names(HML2_heatmap_1)%in%remove,]
  head(HML2_heatmap_3)
  rownames(HML2_heatmap_3)
  HML2_heatmap_final=rbind(HML2_heatmap_3,HML2_heatmap_2)
  head(HML2_heatmap_final)
  rownames(HML2_heatmap_final)
  write.csv(HML2_heatmap_final, file = 'HML2_individual_expression_heatmap.csv')
  
  
  # plot a basic heatmap
  paletteLength <- 50
  myColor <- colorRampPalette(c("white", "blue"))(paletteLength)
  map = pheatmap(as.matrix(HML2_heatmap_final), color = myColor, cluster_rows = FALSE, cluster_cols = FALSE, 
                 main = "Average individual HML-2 expression in GTEx, with housekeeping genes", fontsize = 8, fontsize_row = 6, 
                 fontsize_col = 6)
  
  map = pheatmap(as.matrix(HML2_heatmap_3), color = myColor, cluster_rows = FALSE, cluster_cols = FALSE, 
                 main = "Average individual HML-2 expression in GTEx", fontsize = 8, fontsize_row = 6, 
                 fontsize_col = 6)
  
#separate HML-2 expression based on biological sex
    #filter males and females based on y chrom
      #load in individual expression files for each tissue and assign to a list
      filelist = list.files(pattern = "*V15_Telescope_output.csv") #save data as a list
      filelist
      datalist = lapply(filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
      datalist

      #name each data frame within the list by the tissue specifier within the file name
      filelist_edited = as.character(strsplit(filelist, "_Telescope_output.csv"))
      filelist_edited
      names(datalist) = filelist_edited
      head(datalist)
      
      #convert raw counts to TPM for whole genome counts
      for (i in (1:length(datalist))) {
        Tissue = stringr::str_replace((names(datalist)[[i]]), '_Telescope_output.csv', '')
        Tissue 
        df <- datalist[[i]]
        head(df)
        
        #convert to TPM and filter poorly expressed genes (TPM < 1 in at least half of the samples)
        #read in hg38.gtf for transcript coordinates
        ens = ensemblGenome()
        read.gtf(ens, "../hg38.gtf")
        class(ens)
        Identifier = "exon"
        hg38_Annotations = extractFeature(ens, Identifier)
        hg38_Annotations
        hg38_Annotation_df = data.frame(start=getGtf(hg38_Annotations)$start,end=getGtf(hg38_Annotations)$end,gene_id=getGtf(hg38_Annotations)$gene_id, transcript_id=getGtf(hg38_Annotations)$transcript_id)
        hg38_Annotation_df[1:5,]
        
        #for each row, get the difference between start and end
        hg38_Annotation_df$difference = hg38_Annotation_df$end-hg38_Annotation_df$start
        hg38_Annotation_df[1:5,]
        
        #for each gene_id, sum the difference to get transcript length. This also converts length to kb.
        hg38_Annotation_df_lengthSum = ddply(hg38_Annotation_df, .(gene_id), summarise, difference=(sum(difference)/1000)) #with Telescope
        hg38_Annotation_df_lengthSum[1:5,]
        
        # Divide the read counts by the length of each gene (transcript?) in kilobases. This gives you reads per kilobase (RPK).
        stuff = merge(df, hg38_Annotation_df_lengthSum, by.x = "X", by.y = "gene_id") #with Telescoppe
        head(stuff)
        tail(stuff)
        stuff_mod = stuff[,-1]
        rownames(stuff_mod) = stuff[,1]
        head(stuff_mod)
        tail(stuff_mod)
        RPK = stuff_mod/(stuff_mod$difference)
        head(RPK)
        
        # Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
        Per_Million=(colSums(RPK))/1e6
        Per_Million
        
        # Divide the RPK values by the “per million” scaling factor. This gives you TPM.
        Counts_HML2_TPM <- t(t(RPK)/Per_Million)
        head(Counts_HML2_TPM)[1:5,1:5]
        write.csv(Counts_HML2_TPM, file=paste(Tissue, "TPM_total.csv", sep="_"))
        
        #log2(counts+1) transform 
        Counts_HML2_TPM_Filtered_Log2Trans = log2(Counts_HML2_TPM + 1)
        head(Counts_HML2_TPM_Filtered_Log2Trans)
        rownames(Counts_HML2_TPM_Filtered_Log2Trans)
        write.csv(Counts_HML2_TPM_Filtered_Log2Trans, file=paste(Tissue, "TPM_Log2_total.csv", sep="_"))
      }
    
    #plot PCA to make sure they separate properly
      counts_filelist = list.files(pattern = "*TPM_Log2_total.csv") #save data as a list
      counts_filelist
      counts_datalist = lapply(counts_filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
      counts_datalist

      #name each data frame within the list by the tissue specifier within the file name
      counts_filelist_edited = as.character(strsplit(counts_filelist, "_TPM_Log2_total.csv"))
      counts_filelist_edited
      names(counts_datalist) = counts_filelist_edited
      head(counts_datalist)
      
      #load in meta data from Aidan
      meta = read.csv("../SubjID_Pheno_V8.csv",header = TRUE)
      meta_mod = meta[,c("SUBJID","AGE","SEX")]
      colnames(meta_mod) = c("SUBJID","age","sex")
      head(meta_mod)
      
      for (i in (1:length(counts_datalist))) {
        Tissue = stringr::str_replace((names(counts_datalist)[[i]]), '_V15', '')
        Tissue
        #read in counts data
        df <- counts_datalist[[i]]
        rownames(df) = df$X
        df$X = NULL
        df$difference = NULL
        head(df)
        
        #read in run table to filter and merge with meta data
        #SRA_run_ID = read.delim(paste(Tissue,"SraRunTable.txt",sep="_"),header=TRUE,sep="\t")
        #SRA_run_ID_to_merge = SRA_run_ID[,c("Run","Sample_Name","sex")]
        GTEx_sample_ID = read.delim("../GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt", header = TRUE, sep = "\t")
        GTEx_sample_ID_merge = GTEx_sample_ID[,c("SAMPID", "SMTSD")]
        GTEx_sample_ID_merge$SMTSD <- gsub(" -","", GTEx_sample_ID_merge$SMTSD)
        GTEx_sample_ID_merge$SMTSD <- gsub(" ", "_", GTEx_sample_ID_merge$SMTSD)
        GTEx_sample_ID_merge$SMTSD <- sub("_\\(.*", "", GTEx_sample_ID_merge$SMTSD)
        GTEx_sample_ID_tissue <- GTEx_sample_ID_merge %>% filter(SMTSD == Tissue)
        string = GTEx_sample_ID_tissue$SAMPID
        #subjid = sub("^([^.]+.[^.]+).*", "\\1", string)
        #subjid = gsub("\\.", "-", string)
        #colnames(df) = subjid
        subjid = sub("^(.*?-.*?)-.*", "\\1", string)
        head(subjid)
        length(subjid)
        # 385
        GTEx_sample_ID_tissue$SUBJID = subjid
        #SRA_run_ID_to_merge$SUBJID = subjid
        meta_mod$SUBJID
        meta_updated = join(GTEx_sample_ID_tissue,meta_mod,by="SUBJID",type="inner")
        head(meta_updated)
        dim(meta_updated)

        #k <- which(colnames(df) %in% meta_updated$Run)
        k <- which(colnames(df) %in% meta_updated$SAMPID)
        df <- df[,k]
        meta_filtered <- meta_updated[match(colnames(df),meta_updated$SAMPID),]
        rownames(meta_filtered) <- colnames(df)
        dim(df)
        # 26085   385
        #dim(meta_updated)
        # 385   6
        dim(meta_filtered)
        # 385   6
        
        obj = ExpressionSet(as.matrix(df))
        phenoData(obj) = AnnotatedDataFrame(meta_filtered)
        obj
        
        #grab rows from counts matrices that correspond to the Y chrom
        ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",host = "useast.ensembl.org")
        t2g<-getBM(attributes=c('chromosome_name','start_position','end_position','hgnc_symbol'), mart = ensembl)
        table(t2g$chromosome_name)
        # y = 522
        my_ids <- data.frame(hgnc_symbol=rownames(assayData(obj)$exprs))
        head(my_ids)
        my_ids_final <- merge(my_ids, t2g, by= 'hgnc_symbol')
        head(my_ids_final)
        table(my_ids_final$chromosome_name)
        # y = 70
        
        #add featureData to eset object for Y-specific 
        storageMode(obj) <- "environment"
        featureData(obj) = new("AnnotatedDataFrame", data=my_ids_final)
        storageMode(obj) <- "lockedEnvironment"
        print(obj)
        
        #extract expression data from counts matrix that originate from the Y chromosome
        listOfIDs = fData(obj)[grep("Y", fData(obj)$chromosome_name),]
        head(listOfIDs)
        Male_exprs = subset(assayData(obj)$exprs, rownames(assayData(obj)$exprs) %in% listOfIDs$hgnc_symbol)
        head(Male_exprs)
        dim(listOfIDs)
        #70  4
        dim(Male_exprs)
        #70 385
        
        #plot PCA for y chromosome based gene expression
        mypath <- file.path(paste(dir,paste("PCA_1_2_", Tissue, ".png", sep = ""), sep ="/"))
        png(file=mypath)
        plotPCA(Male_exprs, groups = as.numeric(pData(obj)$sex), groupnames = levels(pData(obj)$sex), pcs = c(1, 2), main = paste(Tissue,"Male vs Female Y chrom only",sep = " "), )
        dev.off()
      }
      #the tissues that have outliers are: Adipose Subcu, Breast Mammary, Cells Fibroblasts, Colon Transverse,
      #  Esophagus Mucosa, Heart Atrial Appendage, Lung, Pancreas, Skin NSE, Skin SE, Whole Blood
      outlier_list = counts_datalist[c("Adipose_Subcutaneous","Breast_Mammary_Tissue","Cells_Transformed_Fibroblasts","Colon_Transverse","Esophagus_Mucosa",
                                     "Heart_Atrial_Appendage","Lung","Pancreas","Skin_NSE","Skin_SE","Whole_Blood")]
      outlier_list
      names(outlier_list)

      for (i in (1:length(outlier_list))) {
        Tissue = stringr::str_replace((names(outlier_list)[[i]]), '_TPM_Log2_total.csv', '')
        Tissue
        #read in counts data
        df <- outlier_list[[i]]
        rownames(df) = df$X
        df$X = NULL
        df$difference = NULL
        head(df)
        
        #read in run table to merge with meta data
        SRA_run_ID = read.delim(paste(Tissue,"SraRunTable.txt",sep="_"),header=TRUE,sep="\t")
        SRA_run_ID_to_merge = SRA_run_ID[,c("Run","Sample_Name","sex")]
        string = SRA_run_ID_to_merge$Sample_Name
        subjid = sub("^(.*?-.*?)-.*", "\\1", string)
        head(subjid)
        length(subjid)
        # 385
        SRA_run_ID_to_merge$SUBJID = subjid
        meta_updated = join(SRA_run_ID_to_merge,meta_mod,by="SUBJID",type="inner")
        head(meta_updated)
        dim(meta_updated)
        
        k <- which(colnames(df) %in% meta_updated$Run)
        df <- df[,k]
        meta_filtered <- meta_updated[match(colnames(df),meta_updated$Run),]
        rownames(meta_filtered) <- colnames(df)
        dim(df)
        # 26085   385
        dim(meta_updated)
        # 385   6
        dim(meta_filtered)
        # 385   6
        
        obj = ExpressionSet(as.matrix(df))
        phenoData(obj) = AnnotatedDataFrame(meta_filtered)
        obj
        
        #grab rows from counts matrices that correspond to the Y chrom
        ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",host = "useast.ensembl.org")
        t2g<-getBM(attributes=c('chromosome_name','start_position','end_position','hgnc_symbol'), mart = ensembl)
        table(t2g$chromosome_name)
        # y = 522
        my_ids <- data.frame(hgnc_symbol=rownames(assayData(obj)$exprs))
        head(my_ids)
        my_ids_final <- merge(my_ids, t2g, by= 'hgnc_symbol')
        head(my_ids_final)
        table(my_ids_final$chromosome_name)
        # y = 70
        
        #add featureData to eset object for Y-specific 
        storageMode(obj) <- "environment"
        featureData(obj) = new("AnnotatedDataFrame", data=my_ids_final)
        storageMode(obj) <- "lockedEnvironment"
        print(obj)
        
        #extract expression data from counts matrix that originate from the Y chromosome
        listOfIDs = fData(obj)[grep("Y", fData(obj)$chromosome_name),]
        head(listOfIDs)
        Male_exprs = subset(assayData(obj)$exprs, rownames(assayData(obj)$exprs) %in% listOfIDs$hgnc_symbol)
        head(Male_exprs)
        dim(listOfIDs)
        #70  4
        dim(Male_exprs)
        #70 385
        
        #get PC data
        stuff = prcomp(Male_exprs)
        test = stuff$rotation
        dim(test)
        #385   70
        head(test)
        PC = as.data.frame(test[, -c(3:70)])
        PC$Run = rownames(PC)
        PC_1 = PC[c("Run","PC1")]
        ID_Sex = pData(obj)[c("Run","sex")]
        PC_Sex = join(ID_Sex,PC_1,by = "Run")
        write.csv(PC_Sex, file = paste(Tissue,"PC1_matrix_for_outlier_removal.csv",sep="_"))
      }

    #remove outliers
      #remove outliers , then redo PCA plot
      #outlier run IDs: SRR1317387, SRR1311266, SRR1366790, SRR1436085, SRR1420413, SRR1470313, SRR1479263, SRR1361860, SRR1445064, SRR1436763, SRR1356577 
      drop = c("SRR1317387", "SRR1311266", "SRR1366790", "SRR1436085", "SRR1420413", "SRR1470313", "SRR1479263", "SRR1361860", "SRR1445064", "SRR1436763", "SRR1356577")
      # remove outliers and remake eset obj, then re-plot PCA
      counts_filelist = list.files(pattern = "*TPM_Log2_total.csv") #save data as a list
      counts_filelist
      counts_datalist = lapply(counts_filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
      counts_datalist
      
      #name each data frame within the list by the tissue specifier within the file name
      counts_filelist_edited = as.character(strsplit(counts_filelist, "_TPM_Log2_total.csv"))
      counts_filelist_edited
      names(counts_datalist) = counts_filelist_edited
      head(counts_datalist)
      
      for (i in (1:length(counts_datalist))) {
        Tissue = stringr::str_replace((names(counts_datalist)[[i]]), '_TPM_Log2_total.csv', '')
        Tissue
        #read in counts data
        df <- counts_datalist[[i]]
        rownames(df) = df$X
        df$X = NULL
        df$difference = NULL
        head(df)
        df_outliersRemoved = df[,!(names(df) %in% drop)]
        
        #read in run table to merge with meta data
        SRA_run_ID = read.delim(paste(Tissue,"SraRunTable.txt",sep="_"),header=TRUE,sep="\t")
        SRA_run_ID_to_merge = SRA_run_ID[,c("Run","Sample_Name","sex")]
        string = SRA_run_ID_to_merge$Sample_Name
        subjid = sub("^(.*?-.*?)-.*", "\\1", string)
        head(subjid)
        length(subjid)
        # 385
        SRA_run_ID_to_merge$SUBJID = subjid
        meta_updated = join(SRA_run_ID_to_merge,meta_mod,by="SUBJID",type="inner")
        head(meta_updated)
        dim(meta_updated)
        
        k <- which(colnames(df_outliersRemoved) %in% meta_updated$Run)
        df_outliersRemoved <- df_outliersRemoved[,k]
        meta_filtered <- meta_updated[match(colnames(df_outliersRemoved),meta_updated$Run),]
        rownames(meta_filtered) <- colnames(df_outliersRemoved)
        dim(df_outliersRemoved)
        # 26085   385
        dim(meta_updated)
        # 385   6
        dim(meta_filtered)
        # 385   6
        
        obj = ExpressionSet(as.matrix(df_outliersRemoved))
        phenoData(obj) = AnnotatedDataFrame(meta_filtered)
        obj
        
        #grab rows from counts matrices that correspond to the Y chrom
        ensembl = useMart(biomart = "ENSEMBL_MART_ENSEMBL", dataset = "hsapiens_gene_ensembl",host = "useast.ensembl.org")
        t2g<-getBM(attributes=c('chromosome_name','start_position','end_position','hgnc_symbol'), mart = ensembl)
        table(t2g$chromosome_name)
        # y = 522
        my_ids <- data.frame(hgnc_symbol=rownames(assayData(obj)$exprs))
        head(my_ids)
        my_ids_final <- merge(my_ids, t2g, by= 'hgnc_symbol')
        head(my_ids_final)
        table(my_ids_final$chromosome_name)
        # y = 70
        
        #add featureData to eset object for Y-specific 
        storageMode(obj) <- "environment"
        featureData(obj) = new("AnnotatedDataFrame", data=my_ids_final)
        storageMode(obj) <- "lockedEnvironment"
        print(obj)
        
        #extract expression data from counts matrix that originate from the Y chromosome
        listOfIDs = fData(obj)[grep("Y", fData(obj)$chromosome_name),]
        head(listOfIDs)
        Male_exprs = subset(assayData(obj)$exprs, rownames(assayData(obj)$exprs) %in% listOfIDs$hgnc_symbol)
        head(Male_exprs)
        dim(listOfIDs)
        #70  4
        dim(Male_exprs)
        #70 465
        
        #plot PCA for y chromosome based gene expression
        mypath <- file.path(paste(dir,paste("PCA_1_2_outliers_removed", Tissue, ".png", sep = ""), sep ="/"))
        png(file=mypath)
        plotPCA(Male_exprs, groups = as.numeric(pData(obj)$sex), groupnames = levels(pData(obj)$sex), pcs = c(1, 2), main = paste(Tissue,"Male vs Female Y chrom only",sep = " "))
        dev.off()
      }
     
    #plot expression
    #now that we know which runs to remove, I'll remove those runs and then plot the expression data
      HML2_filelist = list.files(pattern = "*_TPM_HML2.csv") #save data as a list
      HML2_filelist
      HML2_datalist = lapply(HML2_filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
      HML2_datalist
      HML2_filelist_edited = as.character(strsplit(HML2_filelist, "_TPM_HML2.csv"))
      HML2_filelist_edited
      names(HML2_datalist) = HML2_filelist_edited
      head(HML2_datalist)
      
      HML2_datalist_T = list()
      
      for (i in (1:length(HML2_datalist))) {
        HML2_Tissue = stringr::str_replace((names(HML2_datalist)[[i]]), '_TPM_HML2.csv', '')
        HML2_Tissue
        HML2_df <- HML2_datalist[[i]]
        head(HML2_df)[1:5,1:5]
        HML2_df_1 = HML2_df[,!(names(HML2_df) %in% drop)]
        HML2_DF = HML2_df_1[!grepl("GAPDH", HML2_df_1$X),]
        HML2_DF = HML2_DF[!grepl("ACTB", HML2_DF$X),]
        HML2_DF$X
        HML2_DF$difference = NULL
        head(HML2_DF)
        colnames(HML2_DF)
        HML2_DF_SampleSum = rbind(HML2_DF, data.frame(X = "HML2_Sum", t(colSums(HML2_DF[, -1]))))
        HML2_DF_SampleSum$X
        HML2_DF_SampleSum[HML2_DF_SampleSum$X %in% "HML2_Sum", ]
        
        #grab sum row and colnames
        HML2_DF_SampleSum = HML2_DF_SampleSum[grepl("HML2_Sum", HML2_DF_SampleSum$X),]
        rownames(HML2_DF_SampleSum) = HML2_DF_SampleSum$X
        HML2_DF_SampleSum$X = NULL
        head(HML2_DF_SampleSum)
        
        #transpose and add new column with tissue id
        ERV_Counts_modT = as.data.frame(t(HML2_DF_SampleSum))
        ERV_Counts_modT$tissue = HML2_Tissue
        ERV_Counts_modT$Run <- rownames(ERV_Counts_modT)
        head(ERV_Counts_modT)
        
        #read in run table to merge with meta data
        SRA_run_ID = read.delim(paste(HML2_Tissue,"SraRunTable.txt",sep="_"),header=TRUE,sep="\t")
        SRA_run_ID_to_merge = SRA_run_ID[,c("Run","Sample_Name","body_site","sex")]
        string = SRA_run_ID_to_merge$Sample_Name
        subjid = sub("^(.*?-.*?)-.*", "\\1", string)
        head(subjid)
        length(subjid)
        # 385
        SRA_run_ID_to_merge$SUBJID = subjid
        meta_updated = join(SRA_run_ID_to_merge,meta_mod,by="SUBJID",type="inner")
        head(meta_updated)
        dim(meta_updated)

        ERV_Counts_modT_sex = join(ERV_Counts_modT, meta_updated, by = "Run", type = "left")
        head(ERV_Counts_modT_sex)
        
        HML2_datalist_T[[HML2_Tissue]] <- ERV_Counts_modT_sex # add it to your list
      }
      
      HML2_datalist_T
      big_data = do.call(rbind, HML2_datalist_T)
      head(big_data)
      big_data = cbind(big_data, read.table(text=row.names(big_data), sep=".", 
                                            header=FALSE, col.names = paste0("col", 1:2), stringsAsFactors=FALSE))
      test = big_data[c("Run","HML2_Sum", "tissue", "sex")]
      head(test)
      HML2_outliersRemoved = test[!test$Run %in% drop, ]
      head(HML2_outliersRemoved)
      dim(HML2_outliersRemoved)

      write.csv(HML2_outliersRemoved,file="HML2_outliers_removed_for_sex_differences.csv")
      
      test=read.csv(file="HML2_outliers_removed_for_sex_differences.csv", header = TRUE,sep=",")
      head(test)
      test$X = NULL
      
      test$tissue <- as.factor(test$tissue)
      test$sex <- as.factor(test$sex)
      head(test)
      #      SRR_ID HML2_TPM               tissue    sex
      #1 SRR1069097 36.02131 Adipose_Subcutaneous female
      #2 SRR1070184 52.11286 Adipose_Subcutaneous female
      #3 SRR1070503 34.25826 Adipose_Subcutaneous   male
      #4 SRR1071453 67.78992 Adipose_Subcutaneous   male
      #5 SRR1072247 71.25833 Adipose_Subcutaneous female
      #6 SRR1073581 68.49832 Adipose_Subcutaneous   male
      
      write.csv(table(test$tissue), "Tissue_Distribution_outlierDropped.csv")
      #checked the number of samples per tissue to make sure the samples were dropped. 
      # looks like they were removed. continuing.
      table(test$sex)
      #female   male 
      #3677   6090   

      # ggplot code
      p<-ggplot(test, aes(x=tissue, y=HML2_Sum, fill = sex)) + geom_boxplot()+
        geom_dotplot(binaxis='y', stackdir='center', binwidth = 1, stackratio = 0.025) + 
        labs(title="Differences between biological sex in regards to total HML-2 expression in GTEx",x="Tissue", y = "HML2 TPM (sum/sample)")
      p + stat_n_text(angle = 90, size=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#graph heatmap of individual HML-2 proviruses in sex-specific data
      #load in big expression matrix generated above
      filelist = list.files(pattern = "*_TPM_HML2.csv") #save data as a list
      filelist
      datalist = lapply(filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
      datalist
      
      #name each data frame within the list by the tissue specifier within the file name
      filelist_edited = as.character(strsplit(filelist, "_TPM_HML2.csv"))
      filelist_edited
      names(datalist) = filelist_edited
      
      HML2_datalist_AvgStDev = list()
      
      #calculate average and stdev for each provirus
      for (i in (1:length(datalist))) {
        HML2_Tissue = stringr::str_replace((names(HML2_datalist)[[i]]), '_TPM_HML2.csv', '')
        HML2_Tissue
        HML2_df <- HML2_datalist[[i]]
        head(HML2_df)[1:5,1:5]
        HML2_df$X
        HML2_df$difference = NULL
        head(HML2_df)
        
        #drop runs that are outliers
        HML2_df_1 = HML2_df[,!(names(HML2_df) %in% drop)]
        rownames(HML2_df_1) = HML2_df_1$X
        HML2_df_1$X = NULL
        head(HML2_df_1)
        
        #separate data frame into males and females
          #read in run table to merge with meta data
          SRA_run_ID = read.delim(paste(HML2_Tissue,"SraRunTable.txt",sep="_"),header=TRUE,sep="\t")
          SRA_run_ID_to_merge = SRA_run_ID[,c("Run","Sample_Name","body_site","sex")]
          string = SRA_run_ID_to_merge$Sample_Name
          subjid = sub("^(.*?-.*?)-.*", "\\1", string)
          head(subjid)
          length(subjid)
          # 385
          SRA_run_ID_to_merge$SUBJID = subjid
          meta_updated = join(SRA_run_ID_to_merge,meta_mod,by="SUBJID",type="inner")
          head(meta_updated)
          dim(meta_updated)
          
          #collect run names for males and females and assign them to a variable
          male_IDs_df = meta_updated[meta_updated$sex %in% "male",]
          female_IDs_df = meta_updated[meta_updated$sex %in% "female",]
          
          male_IDs = male_IDs_df$Run
          female_IDs = female_IDs_df$Run
          
          #extract columns thata match the female/male variables and save them as their own dataframes
          male_exprs_df = HML2_df_1[,colnames(HML2_df_1) %in% male_IDs]
          female_exprs_df = HML2_df_1[,colnames(HML2_df_1) %in% female_IDs]
        
        #get average expression for both data sets  
        male_exprs_df$Average = rowMeans(male_exprs_df)
        male_exprs_df$StDev = rowSds(as.matrix(male_exprs_df))
        male_HML2_DF_AvgStDev <- male_exprs_df[, c("Average", "StDev")]
        colnames(male_HML2_DF_AvgStDev) = c("male_Average", "male_StDev")
        head(male_HML2_DF_AvgStDev)
        
        female_exprs_df$Average = rowMeans(female_exprs_df)
        female_exprs_df$StDev = rowSds(as.matrix(female_exprs_df))
        female_HML2_DF_AvgStDev <- female_exprs_df[, c("Average", "StDev")]
        colnames(female_HML2_DF_AvgStDev) = c("female_Average", "female_StDev")
        head(female_HML2_DF_AvgStDev)
        
        #merge the two data frames into one
        HML2_DF_AvgStDev = cbind(male_HML2_DF_AvgStDev,female_HML2_DF_AvgStDev)
        
        HML2_datalist_AvgStDev[[HML2_Tissue]] <- HML2_DF_AvgStDev # add it to your list
      }
      HML2_datalist_AvgStDev
      big_data = do.call(rbind, HML2_datalist_AvgStDev)
      head(big_data)
      big_data$test = rownames(big_data)
      big_data$provirus = gsub("^.*?\\.","", big_data$test)
      big_data$tissue = gsub("[.][\\s\\S]*$", "", big_data$test, perl = T)
      rownames(big_data)  = c()
      big_data$test = NULL
      head(big_data)
      #change NAs in females and males to 0
      big_data$female_Average[is.na(big_data$female_Average)] = 0
      big_data$male_Average[is.na(big_data$male_Average)] = 0
      #get difference between males and females (>0 = female, <0 = male)
      big_data$difference = big_data$female_Average - big_data$male_Average
      big_data_transf = big_data[,c(5,7,6)]
      head(big_data_transf)
      
      #save data 
      write.csv(big_data_transf, file = 'HML2_individual_expression_sex_07272020.csv')
      
      #get data for heatmap
      HML2_heatmap = big_data_transf
      head(HML2_heatmap)
      HML2_heatmap_1 = reshape2::dcast(HML2_heatmap, provirus ~ tissue, value.var = "difference", fun.aggregate = mean)
      rownames(HML2_heatmap_1) = HML2_heatmap_1$provirus
      HML2_heatmap_1$provirus = NULL
      rownames(HML2_heatmap_1)
      HML2_heatmap_2 = HML2_heatmap_1[c("ACTB","GAPDH"),]
      head(HML2_heatmap_2)
      rownames(HML2_heatmap_2)
      remove = c("ACTB","GAPDH")
      HML2_heatmap_3 = HML2_heatmap_1[!row.names(HML2_heatmap_1)%in%remove,]
      head(HML2_heatmap_3)
      rownames(HML2_heatmap_3)
      HML2_heatmap_final=rbind(HML2_heatmap_3,HML2_heatmap_2)
      head(HML2_heatmap_final)
      rownames(HML2_heatmap_final)
      write.csv(HML2_heatmap_3, file = 'HML2_individual_expression_sex_heatmap_07272020.csv')
      
      # plot a basic heatmap
      paletteLength <- 10
      myColor <- colorRampPalette(c("blue","white"))(paletteLength)
      map = pheatmap(as.matrix(HML2_heatmap_3), color = myColor, cluster_rows = FALSE, cluster_cols = FALSE, 
                     main = "Difference of average HML-2 expression in GTEx based on sex", fontsize = 8, fontsize_row = 5, 
                     fontsize_col = 6)
      

#plot HML-2 based on age
      #load in expression data and save as a list
      HML2_filelist = list.files(pattern = "*_TPM_HML2.csv")
      HML2_filelist
      HML2_datalist = lapply(HML2_filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
      HML2_datalist
      HML2_filelist_edited = as.character(strsplit(HML2_filelist, "_TPM_HML2.csv"))
      HML2_filelist_edited
      names(HML2_datalist) = HML2_filelist_edited
      head(HML2_datalist)
      
      HML2_datalist_T = list()
      
      for (i in (1:length(HML2_datalist))) {
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
        HML2_DF_SampleSum = rbind(HML2_DF, data.frame(X = "HML2_Sum", t(colSums(HML2_DF[, -1]))))
        HML2_DF_SampleSum$X
        HML2_DF_SampleSum[HML2_DF_SampleSum$X %in% "HML2_Sum", ]
        
        #grab sum row and colnames
        HML2_DF_SampleSum = HML2_DF_SampleSum[grepl("HML2_Sum", HML2_DF_SampleSum$X),]
        rownames(HML2_DF_SampleSum) = HML2_DF_SampleSum$X
        HML2_DF_SampleSum$X = NULL
        head(HML2_DF_SampleSum)
        
        #transpose and add new column with tissue id
        ERV_Counts_modT = as.data.frame(t(HML2_DF_SampleSum))
        ERV_Counts_modT$tissue = HML2_Tissue
        ERV_Counts_modT$Run <- rownames(ERV_Counts_modT)
        head(ERV_Counts_modT)
        
        #read in run table to merge with meta data
        SRA_run_ID = read.delim(paste(HML2_Tissue,"SraRunTable.txt",sep="_"),header=TRUE,sep="\t")
        SRA_run_ID_to_merge = SRA_run_ID[,c("Run","Sample_Name","body_site")]
        string = SRA_run_ID_to_merge$Sample_Name
        subjid = sub("^(.*?-.*?)-.*", "\\1", string)
        head(subjid)
        length(subjid)
        # 385
        SRA_run_ID_to_merge$SUBJID = subjid
        meta_updated = join(SRA_run_ID_to_merge,meta_mod,by="SUBJID",type="inner")
        head(meta_updated)
        dim(meta_updated)
        
        ERV_Counts_modT_age = join(ERV_Counts_modT, meta_updated, by = "Run", type = "left")
        head(ERV_Counts_modT_age)
        
        HML2_datalist_T[[HML2_Tissue]] <- ERV_Counts_modT_age # add it to your list
      }
      
      HML2_datalist_T
      big_data = do.call(rbind, HML2_datalist_T)
      head(big_data)
      big_data = cbind(big_data, read.table(text=row.names(big_data), sep=".", 
                                            header=FALSE, col.names = paste0("col", 1:2), stringsAsFactors=FALSE))
      test = big_data[c("Run","HML2_Sum", "tissue", "age")]
      head(test)
      range(test$age)
      #20 70
      test$range = ifelse(test$age >= 20 & test$age <= 35, "20-35", ifelse(test$age >=36 & test$age <= 51,"36-51",ifelse(test$age >= 52 & test$age <= 70, "52-70", "NA")))
      head(test)
      
      write.csv(test,file="HML2_expression_age_based.csv")
      
      test=read.csv(file="HML2_expression_age_based.csv", header = TRUE,sep=",")
      head(test)
      test$X = NULL
      
      test$tissue <- as.factor(test$tissue)
      test$range <- as.factor(test$range)
      head(test)
      #         Run HML2_Sum               tissue age range
      #1 SRR1069097 36.02131 Adipose_Subcutaneous  49 36-51
      #2 SRR1070184 52.11286 Adipose_Subcutaneous  59 52-70
      #3 SRR1070503 34.25826 Adipose_Subcutaneous  45 36-51
      #4 SRR1071453 67.78992 Adipose_Subcutaneous  51 36-51
      #5 SRR1072247 71.25833 Adipose_Subcutaneous  70 52-70
      #6 SRR1073581 68.49832 Adipose_Subcutaneous  54 52-70
      
      table(test$range)
      # 20-35  36-51   52-70 
      # 1222   2760    5796   
      
      # ggplot code
      p<-ggplot(test, aes(x=tissue, y=HML2_Sum, fill = range)) + geom_boxplot()+
        geom_dotplot(binaxis='y', stackdir='center', binwidth = 1, stackratio = 0.025) + 
        labs(title="Differences in HML-2 expression in regards to age in GTEx",x="Tissue", y = "HML2 TPM (sum/sample)")
      p + stat_n_text(angle = 90, size=3) + theme(axis.text.x = element_text(angle = 90, hjust = 1))

#table showing the number of times a provirus is expressed in each tissue
      HML2_filelist = list.files(pattern = "*_TPM_HML2.csv") #save data as a list
      HML2_filelist
      HML2_datalist = lapply(HML2_filelist, read.csv, header=TRUE, sep =",", stringsAsFactors=FALSE)
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
        #count the number of times a cell for a given provirus is >0 and save as a new column
        HML2_DF$count <- rowSums(HML2_DF[, -1]!=0 )
        names(HML2_DF)[names(HML2_DF) == 'count'] <- paste(HML2_Tissue,'count',sep="_")
        HML2_df_occurence = HML2_DF[,c("X",paste(HML2_Tissue,'count',sep="_"))]
        #add to empty list
        HML2_datalist_T[[HML2_Tissue]] <- HML2_df_occurence # add it to your list
      }
      #save table of provirus names(rownames) and tissue (colnames) with occurrence counts in each col
      HML2_datalist_T
      big_data = reduce(HML2_datalist_T, full_join, by = "X")
      head(big_data)
      rownames(big_data) = big_data$X
      write.csv(big_data,"Provirus_occurences_within_tissues_07282020.csv")
      