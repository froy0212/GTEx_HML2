
#Filter and Normalize GTEx counts for WGCNA
library(edgeR)
library(WGCNA)

wd = "~/Documents/Coffin Lab/GTEx_HML2_Expression.nosync/GTEx_V7_V8/V15 Telescope/"
setwd(wd)
#Read in counts file
Tissue = "Breast_Mammary_Tissue"
Tissue

Raw_Counts<- read.csv(paste(Tissue, "V15_Telescope_output.csv", sep="_"), header = TRUE, sep = ",", stringsAsFactors = FALSE)
head(Raw_Counts)
dim(Raw_Counts)

#Filter out all rows with <10 in 90% of samples
Raw_Counts_filt <- Raw_Counts[rowSums(Raw_Counts > 10) >= (.9 * ncol(Raw_Counts)),]
length(grep("HML-2",Raw_Counts$X))
length(grep("HML-2",Raw_Counts_filt$X))

#Make a DGE object
Counts_obj <- DGEList(counts = Raw_Counts_filt[,2:356], genes = Raw_Counts_filt[,1])

#EdgeR normalize
#rownames(Raw_Counts_filt) <- Raw_Counts_filt$X
#Raw_Counts_filt$X = NULL
Counts_obj_norm <- calcNormFactors(Counts_obj, method = "TMM")
cpm_log_counts <- cpm(Counts_obj_norm, log = TRUE, normalized.lib.sizes= TRUE)
cpm_log_counts_frame <- as.data.frame(cpm_log_counts)
rownames(cpm_log_counts_frame) = Raw_Counts_filt$X
head(cpm_log_counts_frame)
Norm_counts <- t(cpm_log_counts_frame)
Norm_counts_frame <- as.data.frame(Norm_counts)

#WGCNA Step 1
#Check for missing values
names(Norm_counts_frame)
gsg = goodSamplesGenes(Norm_counts_frame, verbose = 3)
gsg$allOK

#Cluster Samples
sampleTree = hclust(dist(Norm_counts_frame), method = "average")
sizeGrWindow(12,9)
#pdf(file = "Plots/sampleClustering.pdf", width = 12, height = 9);
par(cex = 0.6);
par(mar = c(0,4,2,0))
plot(sampleTree, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

# Plot a line to show the cut
abline(h = 150, col = "red");
# Determine cluster under the line
clust = cutreeStatic(sampleTree, cutHeight = 150, minSize = 10)
table(clust)
# clust 1 contains the samples we want to keep.
keepSamples = (clust==1)
WCounts = Norm_counts_frame[keepSamples, ]
nGenes = ncol(WCounts)
nSamples = nrow(WCounts)
#Check cutoff
keepClust = hclust(dist(WCounts), method = "average")
plot(keepClust, main = "Sample clustering to detect outliers", sub="", xlab="", cex.lab = 1.5,
     cex.axis = 1.5, cex.main = 2)

#Create network
options(stringsAsFactors = FALSE)
#Find correct power
powers = c(c(1:10), seq(from = 12, to=20, by=2))
sft = pickSoftThreshold(WCounts, powerVector = powers, verbose = 5)
sizeGrWindow(9, 5)
par(mfrow = c(1,2));
cex1 = 0.9;
plot(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     xlab="Soft Threshold (power)",ylab="Scale Free Topology Model Fit,signed R^2",type="n",
     main = paste("Scale independence"));
text(sft$fitIndices[,1], -sign(sft$fitIndices[,3])*sft$fitIndices[,2],
     labels=powers,cex=cex1,col="red")
abline(h=0.90,col="red")
plot(sft$fitIndices[,1], sft$fitIndices[,5],
     xlab="Soft Threshold (power)",ylab="Mean Connectivity", type="n",
     main = paste("Mean connectivity"))
text(sft$fitIndices[,1], sft$fitIndices[,5], labels=powers, cex=cex1,col="red")

#Create network

net = blockwiseModules(WCounts, power = 9,
                       TOMType = "unsigned", minModuleSize = 30,
                       reassignThreshold = 0, mergeCutHeight = 0.25,
                       numericLabels = TRUE, pamRespectsDendro = FALSE,
                       saveTOMs = TRUE,
                       saveTOMFileBase = "GTExBreastTOM",
                       verbose = 3,
                       maxBlockSize = 20000 )

#Check num of modules
table(net$colors)

# open a graphics window
sizeGrWindow(12, 9)
# Convert labels to colors for plotting
mergedColors = labels2colors(net$colors)
# Plot the dendrogram and the module colors underneath
plotDendroAndColors(net$dendrograms[[1]], mergedColors[net$blockGenes[[1]]],
                    "Module colors",
                    dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05)

moduleLabels = net$colors
moduleColors = labels2colors(net$colors)
MEs = net$MEs;
geneTree = net$dendrograms[[1]];
save(MEs, moduleLabels, moduleColors, geneTree,
     file = "GTExBreast-02-networkConstruction-auto.RData")


# bwnet = blockwiseModules(WCounts, maxBlockSize = 160000,
#                          power = 12, TOMType = "unsigned", minModuleSize = 30,
#                          reassignThreshold = 0, mergeCutHeight = 0.25,
#                          numericLabels = TRUE,
#                          saveTOMs = TRUE,
#                          saveTOMFileBase = "GTExCerebellumTOM",
#                          verbose = 3)


length(grep("HML-2",names(WCounts)))
length(grep("HML-2",names(WCounts)[moduleColors=="grey"]))

probes=names(WCounts)
# Create the starting data frame
geneInfo0 = data.frame(substanceBXH = probes,
                       moduleColor = moduleColors)



#Export to Cytoscape
lnames = load(file = "GTExBreast-02-networkConstruction-auto.RData");
lnames
TOM = TOMsimilarityFromExpr(WCounts, power = 14)

modules = c("blue");
# Select module probes
probes = names(WCounts)
inModule = is.finite(match(moduleColors, modules));
modProbes = probes[inModule];
#modGenes = annot$gene_symbol[match(modProbes, annot$substanceBXH)];
# Select the corresponding Topological Overlap
modTOM = TOM[inModule, inModule];
dimnames(modTOM) = list(modProbes, modProbes)
# Export the network into edge and node list files Cytoscape can read
cyt = exportNetworkToCytoscape(modTOM,
                               edgeFile = paste("CytoscapeInput-edges-", paste(modules, collapse="-"), ".txt", sep=""),
                               nodeFile = paste("CytoscapeInput-nodes-", paste(modules, collapse="-"), ".txt", sep=""),
                               weighted = TRUE,
                               threshold = 0.00,
                               nodeNames = modProbes,
                               nodeAttr = moduleColors[inModule]);
HML2networkedge = dplyr::filter(cyt$edgeData,grepl("HML-2",fromNode) | grepl("HML-2",toNode))
write.csv(HML2networkedge, "BreastHML2networkedge.csv")
#Visualize network
# Calculate topological overlap anew: this could be done more efficiently by saving the TOM
# calculated during module detection, but let us do it again here.
dissTOM = 1-TOMsimilarityFromExpr(WCounts, power = 12);
# Transform dissTOM with a power to make moderately strong connections more visible in the heatmap
plotTOM = dissTOM^7;
# Set diagonal to NA for a nicer plot
diag(plotTOM) = NA;
# Call the plot function
sizeGrWindow(9,9)
TOMplot(plotTOM, geneTree, moduleColors, main = "Network heatmap plot, all genes")

