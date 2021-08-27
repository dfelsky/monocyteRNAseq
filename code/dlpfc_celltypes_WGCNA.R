# construct WGCNA network in steps
library(WGCNA)
library(rms)
library(limma)
library(flashClust)
library(cowplot)
library(gplots)
library(BRETIGEA)
library(ggdendro)

load("input/all_genes_ensembl.RData")

# save input dataset for WGCNA
preadjust <- readRDS("output/WGCNA/dlpfc/input_datExpr.rds")

preadjust <- t(preadjust)
genesymbols <- all_genes$external_gene_name[match(rownames(preadjust),all_genes$ensembl_gene_id)]
rownames(preadjust) <- genesymbols
preadjust <- preadjust[!is.na(genesymbols),]
postadjust <- adjustBrainCells(preadjust,species = "human")

datExpr <- t(postadjust$expression)

saveRDS(datExpr,file="output/WGCNA/dlpfc/input_adjustedforcells_datExpr.rds")
saveRDS(postadjust$SVPs,file="output/WGCNA/dlpfc/dlpfc_estimatedBrainCells.rds")

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
  geom_vline(xintercept=16,lty=1,col="darkgreen")+
  labs(y="Signed R2")+
  theme_minimal()
p2 <- ggplot(data=sft$fitIndices,aes(y=mean.k.,x=Power))+
  geom_line(size=1.5,alpha = 0.2, stat = "smooth", method = "loess",col="blue")+
  geom_text(aes(label=Power),size=3)+
  geom_vline(xintercept=16,lty=1,col="darkgreen")+
  labs(y="Mean connectivity")+
  theme_minimal()

pdf(file="output/WGCNA/dlpfc/softPower_study_celltypes.pdf",width = 12,height=6)
print(plot_grid(p1,p2,ncol=2,labels = c("A","B")))
dev.off()

###############################################################################
##############################
##################### ONE STEP 
##############################
datExpr <- readRDS("output/WGCNA/dlpfc/input_adjustedforcells_datExpr.rds")

setwd("output/WGCNA/dlpfc/")
net <- blockwiseModules(datExpr = datExpr,
                         power=12,
                         TOMType = "signed",
                         minModuleSize = 30,
                         minKMEtoStay = 0.3,
                         corType = "bicor",
                         maxBlockSize = 20000,
                         maxPOutliers =0.05,
                         deepSplit = 3,
                         detectCutHeight = 0.995,
                         saveTOMs = T,
                         saveTOMFileBase = "dlpfc_signed_signed_minKMEtoStay03_noPAM_celltypes",
                         pamStage = T,loadTOM = T,
                         pamRespectsDendro = F,
                         networkType = "signed",
                         verbose = 10)
setwd("/Users/dfelsky/Documents/monocyteRNAseq")

saveRDS(net,"output/WGCNA/dlpfc/net_PAM_celltypes.rds")
net <- readRDS("output/WGCNA/dlpfc/net_PAM_celltypes.rds")

pdf(file="output/WGCNA/dlpfc/geneDendrogram_PAM_celltypes.pdf", width = 12, height = 9)
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

pdf(file="output/WGCNA/dlpfc/genesPerModule_barPlot_PAM_celltypes.pdf",width = 12,height=6)
print(ggplot(data=modlengths,aes(y=n,x=module))+
        geom_bar(aes(fill=module),stat="identity",show.legend = F)+
        scale_fill_manual(values=modlengths$hex[order(modlengths$n)])+
        theme_minimal()+
        labs(y="# of genes",x="Module",title = "Number of genes per module")+
        theme(axis.text.x=element_text(angle = -45, hjust = 0)))
dev.off()


########## plot for figure 4
net <- readRDS("output/WGCNA/dlpfc/net_PAM_celltypes.rds")

tiff("paper/figures/Figure4_dendrogram_dlpfc_PAM.tif",w=4,h=3,units="in",res=600)
plotDendroAndColors(net$dendrograms[[1]],
                    net$colors,
                    dendroLabels = FALSE,
                    hang = 0.02,
                    addGuide = F,
                    guideHang = 0.05,
                    main = "DLPFC")
dev.off()
