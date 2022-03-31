setwd("../V15 Counts")
library(edgeR)
counts <-read.csv("Breast_Mammary_Tissue_V15_TPM_Total.csv", sep = ",", header = TRUE, stringsAsFactors = FALSE)
head(counts)
#Clean up Row/Column names and remove unnecessary columns
rownames(counts) <- counts[,1]
counts <- counts[,-1]
counts <- counts[-c(26086:26087),]
colnames(counts) <- gsub("\\.", "-",colnames(counts))

#Remove Y chromosome outliers
colnum <- which(colnames(counts) == "GTEX-11ILO-2226-SM-5A5L1")
counts <- counts[,-colnum]

# #Generate Meta file
meta = read.csv(paste(metadata,"SubjID_Pheno_V8.csv",sep="/"),header = TRUE)
meta_mod = meta[,c("SUBJID","AGE", "SEX")]
colnames(meta_mod) = c("SUBJID","age", "sex")
head(meta_mod)


GTEx_sample_ID = read.delim(paste(dir,"GTEx_Analysis_2017-06-05_v8_Annotations_GTEx_Analysis_2017-06-05_v8_Annotations_SampleAttributesDS.txt", sep = "/"), header = TRUE, sep = "\t")
GTEx_sample_ID_merge = GTEx_sample_ID[,c("SAMPID", "SMTSD")]
GTEx_sample_ID_merge$SMTSD <- gsub(" -","", GTEx_sample_ID_merge$SMTSD)
GTEx_sample_ID_merge$SMTSD <- gsub(" ", "_", GTEx_sample_ID_merge$SMTSD)
GTEx_sample_ID_merge$SMTSD <- sub("_\\(.*", "", GTEx_sample_ID_merge$SMTSD)
GTEx_sample_ID_tissue <- GTEx_sample_ID_merge %>% filter(SMTSD == "Breast_Mammary_Tissue")
string = GTEx_sample_ID_tissue$SAMPID
subjid = sub("^([^.]+.[^.]+).*", "\\1", string)
subjid = gsub("\\.", "-", string)
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
meta_updated$sex[meta_updated$sex == 1] <- "Male"
meta_updated$sex[meta_updated$sex == 2] <- "Female"
meta_updated$squish <- paste(meta_updated$SAMPID, meta_updated$sex, sep = "_")

colnames(counts) <- meta_updated$squish[match(names(counts),meta_updated$SAMPID)]
write.csv(counts, "Breast_TPM_Sexlabel.csv")


xcell <- read.delim("xCell_Breast_TPM_Sexlabel_xCell_1217110921.txt", sep = "\t")

male_col <- xcell %>% select(ends_with("Male", ignore.case = FALSE))
female_col <- xcell %>% select(ends_with("Female"))

xcell_ave <- as.data.frame(rowMeans(male_col))
xcell_ave[,2] <- rowMeans(female_col)
