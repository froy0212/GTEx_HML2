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

setwd("/Users/far122/Desktop/GTEx_HML2_Expression")
Tissue = "Adipose_Subcutaneous"

#load in raw counts data and rename rows to reflect provirus name
counts=read.csv(paste(Tissue, "Telescope_output.csv", sep="_"), header=TRUE, sep=",")
head(counts[,1:5])
tail(counts[,1:5])
#counts_mod = counts[,-1] #FENRIR only
#rownames(counts_mod) = counts[,1] #FENRIR only
#head(counts_mod[,1:5]) #FENRIR only

#convert to TPM and filter low expressed genes (TPM < 1 in at least half of the samples)
#read in hg38.gtf for transcript coordinates
ens = ensemblGenome()
read.gtf(ens, "hg38.gtf")
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
stuff = merge(counts, hg38_Annotation_df_lengthSum, by.x = "X", by.y = "gene_id") #with Telescoppe
head(stuff)
tail(stuff)
stuff_mod = stuff[,-1]
rownames(stuff_mod) = stuff[,1]
head(stuff_mod)
tail(stuff_mod)
#featured_genes = subset(hg38_Annotation_df_lengthSum, gene_id %in% rownames(counts_mod))
#head(featured_genes)
RPK = stuff_mod/(stuff_mod$difference)
head(RPK)

# Count up all the RPK values in a sample and divide this number by 1,000,000. This is your “per million” scaling factor.
Per_Million=(colSums(RPK))/1e6
Per_Million

# Divide the RPK values by the “per million” scaling factor. This gives you TPM.
#Counts_HML2_TPM = RPK/Per_Million

Counts_HML2_TPM <- t(t(RPK)/Per_Million)
Counts_HML2_TPM

 #filter data by TPM < 1 for half the samples
Counts_HML2_TPM_Filtered = Counts_HML2_TPM[, -which(rowMeans(Counts_HML2_TPM < 1) > 0.5)]


keepgenes <- rowSums(Counts_HML2_TPM > 0.5) > ncol(Counts_HML2_TPM)/4
Counts_HML2_TPM_Filtered = Counts_HML2_TPM[keepgenes, ]
dim(Counts_HML2_TPM_Filtered)
Counts_HML2_TPM_Filtered
write.csv(Counts_HML2_TPM_Filtered, file=paste(Tissue, "TPM_HML2.csv", sep="_"))

#log2(counts+1) transform 
Counts_HML2_TPM_Filtered_Log2Trans = log2(Counts_HML2_TPM_Filtered + 1)
head(Counts_HML2_TPM_Filtered_Log2Trans)
write.csv(Counts_HML2_TPM_Filtered_Log2Trans, file=paste(Tissue, "TPM_Log2_HML2.csv", sep="_"))
