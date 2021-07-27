# construct WGCNA network in steps
library(WGCNA)
library(rms)
library(limma)
library(flashClust)
library(cowplot)
library(gplots)
library(ggdendro)
ROSmaster <- readRDS("input/ROSmaster_TWAS_input.rds")

# read in starting dataset
mono.resid <- readRDS("input/monocytes_techandagesex_residuals.rds")

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

pdf(file="output/WGCNA/mono/subjectDendrograms.pdf",width=12,height=6)
print(plot_grid(dend1,dend2,ncol=2,labels = c("A","B")))
dev.off()

# save input dataset for WGCNA
saveRDS(datExpr,file="output/WGCNA/mono/input_datExpr.rds")

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

pdf(file="output/WGCNA/mono/softPower_study.pdf",width = 12,height=6)
print(plot_grid(p1,p2,ncol=2,labels = c("A","B")))
dev.off()

###############################################################################
##############################
##################### ONE STEP
##############################
datExpr <- readRDS("output/WGCNA/mono/input_datExpr.rds")

setwd("output/WGCNA/mono/")
net2 <- blockwiseModules(datExpr = datExpr,
                         power=11,
                         TOMType = "signed",
                         minModuleSize = 30,
                         minKMEtoStay = 0.3,
                         corType = "bicor",
                         maxBlockSize = 20000,
                         maxPOutliers =0.05,
                         deepSplit = 3,
                         detectCutHeight = 0.995,
                         saveTOMs = T,
                         saveTOMFileBase = "mono2_signed_signed_minKMEtoStay03",
                         pamStage = T,
                         pamRespectsDendro = T,
                         networkType = "signed",
                         verbose = 10)
setwd("/Users/dfelsky/Documents/monocyteRNAseq")

table(net1$colors)
plotDendroAndColors(net2$dendrograms[[1]],colors = net2$colors)

saveRDS(net2,"output/WGCNA/mono/net.rds")

#pdf(file="geneDendrogram_unmerged.pdf", width = 12, height = 9)
plotDendroAndColors(geneTree,
                    dynamicColors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.02,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
#dev.off()

#pdf(file="output/WGCNA/mono/genesPerModule_barPlot.pdf",width = 12,height=6)
print(ggplot(data=modlengths,aes(y=members,x=module))+
        geom_bar(aes(fill=module),stat="identity",show.legend = F)+
        scale_fill_manual(values=modlengths$hex[order(modlengths$members)])+
        theme_minimal()+
        labs(y="# of genes",x="Module",title = "Number of genes per module")+
        theme(axis.text.x=element_text(angle = -45, hjust = 0)))
#dev.off()



