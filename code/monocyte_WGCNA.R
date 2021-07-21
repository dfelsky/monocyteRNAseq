# construct WGCNA network in steps
library(WGCNA)
library(rms)
library(limma)
library(flashClust)
library(cowplot)
library(gplots)
library(ggdendro)
ROSmaster <- readRDS("input/ROSmaster_TWAS_input.rds")

# set working directory
setwd("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/WGCNA/mono/run2")

# read in starting dataset
mono.resid <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/normalization/monocytes_techandagesex_residuals.rds")

# cluster subjects and plot dendrogram prior to removal
datExpr <- t(mono.resid)
sampleTree1 <- flashClust(dist(datExpr),method = "average")
dend1 <- ggdendrogram(sampleTree1)+
  labs(title=paste0("n=",nrow(datExpr)," subjects, prior to pruning"))

# remove outliers identified by hierarchical clustering
manual_remove_subjects <- c("00246264","35072859","91804757","21135680","69866926","50303145")
mono <- mono.resid[,which(colnames(mono.resid) %nin% manual_remove_subjects)]
datExpr <- t(mono)

# cluster subjects and plot dendrogram after removal
sampleTree2 <- flashClust(dist(datExpr),method = "average")
dend2 <- ggdendrogram(sampleTree2)+
  labs(title=paste0("n=",nrow(datExpr)," subjects, afterpruning"))

pdf(file="subjectDendrograms.pdf",width=12,height=6)
print(plot_grid(dend1,dend2,ncol=2,labels = c("A","B")))
dev.off()

# save input dataset for WGCNA
saveRDS(datExpr,file="input_datExpr.rds")

# determine optimal power visually - here set at 11
powers = c(1:30)
sft = pickSoftThreshold(datExpr,
                        powerVector=powers, 
                        verbose =5, 
                        networkType="signed",
                        blockSize = ncol(datExpr),
                        corFnc = "bicor",
                        corOptions = list(maxPOutliers=0.05)) 

p1 <- ggplot(data=sft$fitIndices,aes(y=SFT.R.sq,x=Power))+
  geom_hline(yintercept=c(0.85,0.9),lty=2,col="red")+
  geom_line(size=1.5,alpha = 0.2, stat = "smooth", method = "loess",col="blue")+
  geom_text(aes(label=Power),size=3)+
  geom_vline(xintercept=11,lty=1,col="darkgreen")+
  labs(y="Signed R2")+
  theme_minimal()
p2 <- ggplot(data=sft$fitIndices,aes(y=mean.k.,x=Power))+
  geom_line(size=1.5,alpha = 0.2, stat = "smooth", method = "loess",col="blue")+
  geom_text(aes(label=Power),size=3)+
  geom_vline(xintercept=11,lty=1,col="darkgreen")+
  labs(y="Mean connectivity")+
  theme_minimal()

pdf(file="softPower_study.pdf",width = 12,height=6)
print(plot_grid(p1,p2,ncol=2,labels = c("A","B")))
dev.off()

# calculate adjacency matrix
adjacency= adjacency(datExpr,
                     power = 11,
                     type = "signed",
                     distFnc = "bicor",
                     distOptions = list(maxPOutliers=0.05))

# calculate TOM
TOM = TOMsimilarity(adjacency,
                    TOMType="signed",
                    verbose = T)

rownames(TOM) <- colnames(datExpr)
colnames(TOM) <- colnames(datExpr)

saveRDS(TOM,file="TOM_signed_power11.rds")

# perform hierarchical clustering of TOM dissimilarity
dissTOM = 1 - TOM
geneTree = flashClust(as.dist(dissTOM), method = "average")

# define modules
dynamicMods = cutreeDynamic(dendro = geneTree,
                            distM = dissTOM,
                            deepSplit = 3, ### this is changed from 2-3 for run2
                            pamRespectsDendro = FALSE,
                            minClusterSize = 20)

dynamicColors = labels2colors(dynamicMods)

pdf(file="geneDendrogram_unmerged.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree,
                    dynamicColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.02,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

# calculate eigengenes and cluster
MEList = moduleEigengenes(datExpr, colors = dynamicColors)
MEs = MEList$eigengenes
MEDiss = 1-cor(MEs)
METree = flashClust(as.dist(MEDiss), method = "average")
merge = mergeCloseModules(datExpr,
                          dynamicColors,
                          cutHeight = 0.1,
                          verbose = 3)
mergedColors = merge$colors
mergedMEs = merge$newMEs

pdf(file="geneDendrogram_merged.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree,
                    cbind(dynamicColors, mergedColors),
                    c("Dynamic Tree Cut", "Merged dynamic"),
                    dendroLabels = FALSE,
                    hang = 0.02,
                    addGuide = TRUE,
                    guideHang = 0.05)
dev.off()

modlengths <- as.data.frame(do.call(rbind,as.list(table(mergedColors))))
modlengths$module <- rownames(modlengths)
names(modlengths)[1] <- "members"
modlengths$hex <- col2hex(modlengths$module)
modlengths$module <- factor(modlengths$module,levels=modlengths$module[order(modlengths$members)])

pdf(file="genesPerModule_barPlot.pdf",width = 12,height=6)
print(ggplot(data=modlengths,aes(y=members,x=module))+
        geom_bar(aes(fill=module),stat="identity",show.legend = F)+
        scale_fill_manual(values=modlengths$hex[order(modlengths$members)])+
        theme_minimal()+
        labs(y="# of genes",x="Module",title = "Number of genes per module")+
        theme(axis.text.x=element_text(angle = -45, hjust = 0)))
dev.off()

moduleColors = mergedColors
names(moduleCOlors) <- colnames(datExpr)
colorOrder = c("grey", standardColors(100))
moduleLabels = match(moduleColors, colorOrder)-1
MEs = mergedMEs
saveRDS(list(MEs=MEs,
             labels=moduleLabels,
             colors=moduleColors,
             geneTree=geneTree),
             file = "net.rds")
