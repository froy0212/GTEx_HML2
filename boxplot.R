library(ggplot2)

ERV_Counts=read.csv(file="/Users/far122/Downloads/Farrah_Test_Table.csv", header=TRUE,sep=",")
head(ERV_Counts)
colnames(ERV_Counts)
class(ERV_Counts)
ERV_Counts = ERV_Counts[ , -which(names(ERV_Counts) %in% c("X","Run", "COHORT"))]
head(ERV_Counts)

#for Telescope analysis
#HML2only =ERV_Counts[grepl("^HML-2", ERV_Counts$X),, drop = FALSE]
#dim(HML2only)
#FENRIR: 50142 38
#FENRIR-Telescope: 92 39
#ERV_Counts_mod = HML2only[,-1]
#rownames(ERV_Counts_mod) = HML2only[,1]
#head(ERV_Counts_mod)

#for FENRIR only
ERV_Counts_mod = ERV_Counts[,-1]
rownames(ERV_Counts_mod) = ERV_Counts[,1]
head(ERV_Counts_mod)
dim(ERV_Counts_mod)

ERV_expression_subset = subset(ERV_Counts_mod, grepl("^HML-2_",row.names(ERV_Counts_mod)))
head(ERV_expression_subset)

greater_than_0 <- as.data.frame(data.matrix(ERV_expression_subset[rowSums(ERV_expression_subset[,-1])>10, ]))
head(greater_than_0)

ERV_Counts_modT = as.data.frame(t(greater_than_0))
head(ERV_Counts_modT)


# ggplot code
p = ggplot(stack(ERV_Counts_modT), aes(x = ind, y = values)) + geom_boxplot()
p + theme(axis.text.x = element_text(size=8, angle = 90)) + labs(title=paste("HML-2 expression in GTEx for", Tissue, analysis_type,  sep =" "), x ="provirus", y = "fragment counts (raw)")
dim(ERV_Counts_modT)
#FENRIR: 38 64
#FENRIR-Telescope: 38 69
write.csv(t(ERV_Counts_modT), file = paste(Tissue, "HML-2_expression_in_GTEx",analysis_type,".csv", sep="_" ))

# makaing a bar grpah
Benchmarking_Counts=read.csv(file="Nerve_Tibial_HML2.csv", sep=",")
head(Benchmarking_Counts)
Benchmarking_Counts_Mod = Benchmarking_Counts[1:4,]
head(Benchmarking_Counts_Mod)
Benchmarking_Counts_rownames=Benchmarking_Counts_Mod[,-1]
rownames(Benchmarking_Counts_rownames)=Benchmarking_Counts_Mod[,1]
head(Benchmarking_Counts_rownames)
colnames(Benchmarking_Counts_rownames) = c("SRR1106195","SRR1106195 two reads unpaired", "SRR1106195 two reads paired", "SRR1106195 10 reads 80% paired", "SRR1106195 10 reads paired", "SRR1106195 100 reads 80% paired", "SRR1106195 100 reads paired", "SRR1106195 1000 reads 80% paired", "SRR1106195 1000 reads paired")

class(Benchmarking_Counts_rownames)
condition=rep(colnames(Benchmarking_Counts_rownames))
plot = ggplot(Benchmarking_Counts_rownames, aes(x=rownames(Benchmarking_Counts_rownames), y= values))
plot
