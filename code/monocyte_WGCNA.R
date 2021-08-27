# construct WGCNA network in steps
library(WGCNA)
library(rms)
library(limma)
library(flashClust)
library(cowplot)
library(gplots)
library(ggdendro)
library(ggthemes)
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
net <- blockwiseModules(datExpr = datExpr,
                         power=8,
                         TOMType = "signed",
                         minModuleSize = 30,
                         minKMEtoStay = 0.3,
                         corType = "bicor",
                         maxBlockSize = 20000,
                         maxPOutliers =0.05,
                         deepSplit = 3,
                         detectCutHeight = 0.995,
                         saveTOMs = T,
                         saveTOMFileBase = "mono_signed",
                         pamStage = T,
                         pamRespectsDendro = F,
                         networkType = "signed",
                         verbose = 10)
setwd("/Users/dfelsky/Documents/monocyteRNAseq")

saveRDS(net,"output/WGCNA/mono/net_PAM.rds")
net <- readRDS("output/WGCNA/mono/net_PAM.rds")

pdf(file="output/WGCNA/mono/geneDendrogram_PAM.pdf", width = 12, height = 9)
plotDendroAndColors(net$dendrograms[[1]],
                    net$colors,
                    "Dynamic Tree Cut",
                    dendroLabels = FALSE,
                    hang = 0.02,
                    addGuide = TRUE,
                    guideHang = 0.05,
                    main = "Gene dendrogram and module colors")
dev.off()

modlengths <- data.frame(module=names(table(net$colors)),
                         n=as.numeric(table(net$colors)))
modlengths$hex <- col2hex(modlengths$module)
modlengths$module <- factor(modlengths$module,levels=modlengths$module[order(modlengths$n)])
modlengths <- subset(modlengths, module %nin% "grey")

pdf(file="output/WGCNA/mono/genesPerModule_barPlot_PAM.pdf",width = 12,height=6)
print(ggplot(data=modlengths,aes(y=n,x=module))+
        geom_bar(aes(fill=module),stat="identity",show.legend = F)+
        scale_fill_manual(values=modlengths$hex[order(modlengths$n)])+
        theme_minimal()+
        labs(y="# of genes",x="Module",title = "Number of genes per module")+
        theme(axis.text.x=element_text(angle = -45, hjust = 0)))
dev.off()

###### for figure 4
tiff("paper/figures/Figure4_dendrogram_mono_PAM.tif",w=4,h=3,units="in",res=600)
plotDendroAndColors(net$dendrograms[[1]],
                    net$colors,
                    dendroLabels = FALSE,
                    hang = 0.02,
                    addGuide = F,
                    guideHang = 0.05,
                    main = "Monocytes")
dev.off()


#### get hubgenes
datExpr <- readRDS("output/WGCNA/mono/input_datExpr.rds")
net <- readRDS("output/WGCNA/mono/net_PAM.rds")

hubgenesmono <- chooseTopHubInEachModule(datExpr = datExpr,
                         net$colors,
                         power = 8,
                         type = "signed",
                         options(corType = "bicor",maxPOutliers=0.05))


hubdf <- data.frame(module=names(hubgenesmono),hub=all_genes$external_gene_name[match(hubgenesmono,all_genes$ensembl_gene_id)])
hubdf$n <- as.numeric(table(net_mono$colors))
hubdf$hex <- col2hex(hubdf$module)
hubdf$value <- 1
hubdf <- subset(hubdf, module %nin% "grey")

hexindex <- hubdf$hex
names(hexindex) <- hubdf$module

hubdf$hub <- factor(hubdf$hub,levels=hubdf$hub[order(hubdf$n,decreasing = T)])

pdf("paper/figures/Fig4_monocyte_module_barchart_and_hubgenes.pdf",w=5.5,h=2.8)
ggplot(data=hubdf,aes(x=n,y=hub,fill=module))+
  geom_bar(stat="identity",show.legend = F)+
  geom_text(aes(label=module),nudge_x = 100)+
  scale_fill_manual(values=hexindex)+
  labs(y="Hub gene",x="Number of Genes")+
  theme_classic()
dev.off()