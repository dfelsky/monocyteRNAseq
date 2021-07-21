# consensus WGCNA
library(WGCNA)
library(rms)
library(limma)
library(flashClust)
library(cowplot)
library(gplots)

# set working directory
load("input/all_genes_ensembl.RData")

dlpfc_datExpr <- readRDS("../dlpfc/input_adjustedforcells_datExpr.rds")
colnames(dlpfc_datExpr) <- all_genes$ensembl_gene_id[match(colnames(dlpfc_datExpr),all_genes$external_gene_name)]
mono_datExpr <- readRDS("../mono/input_datExpr.rds")

commongenes <- intersect(colnames(dlpfc_datExpr),colnames(mono_datExpr))

multiExpr <- multiData(dlpfc=dlpfc_datExpr[,commongenes],mono=mono_datExpr[,commongenes])
checkSets(multiExpr)

cnet <- blockwiseConsensusModules(multiExpr,
                          power = 12, 
                          minModuleSize = 30, 
                          deepSplit = 2,
                          networkType = "signed",
                          pamRespectsDendro = FALSE,
                          mergeCutHeight = 0.25, 
                          numericLabels = TRUE,
                          minKMEtoStay = 0.2,
                          corType = "bicor",
                          maxBlockSize = 20000,
                          maxPOutliers =0.05,
                          saveTOMs = TRUE, 
                          verbose = 5)

saveRDS(cnet, file="output/WGCNA/cnet.rds")



#
cnet$colors <- labels2colors(cnet$colors)
cnet <- readRDS("cnet.rds")

plotDendroAndColors(cnet$dendrograms[[1]],
                    colors = cnet$colors,
                    dendroLabels = F)



###################################################################
#### get trait associations for each module and each phenotype ####
###################################################################
library(broom)
library(ggrepel)
pheno <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/ROSmaster_TWAS_input.rds")

monoME <- as.data.frame(cnet$multiMEs$mono)
names(monoME) <- paste0(names(monoME),"_mono")
monoME$projid <- rownames(monoME)

dlpfcME <- as.data.frame(cnet$multiMEs$dlpfc)
names(dlpfcME) <- paste0(names(dlpfcME),"_dlpfc")
dlpfcME$projid <- rownames(dlpfcME)

allME <- merge(monoME,dlpfcME,by="projid",all=T)
amr <- merge(allME,pheno,by="projid",all=T)

indepvec.pathology <- c("plaq_n_sqrt","plaq_d_sqrt","amyloid_sqrt","tangles_sqrt","nft_sqrt","ci_num2_gct","ci_num2_mct","arteriol_scler","caa_4gp","cvda_4gp2","dlbdx","hspath_any","tdp_stage4","parkdx","pathoAD","vm3123","pput3123","it3123","mf3123") 
indepvec.cog <- c(grep("at_draw",names(amr),value=T),grep("at_lastvisit",names(amr),value=T))

alldepend <- c(indepvec.pathology,indepvec.cog)
allindep <- names(allME)[-c(1,19,37)]

index <- 1
pvalues <- NULL
bvalues <- NULL
nvalues <- NULL
phenovalues <- NULL
modvalues <- NULL
for (indep in allindep) {
  for (depend in alldepend) {
  
    #form <- formula(paste(indep,"~",depend,"+",paste0(covars,collapse = " + ")))
    form <- formula(paste(indep,"~",depend))
    mod <- lm(data=amr, form)
    pvalues[index] <- tidy(mod)$p.value[2]
    bvalues[index] <- coef(mod)[2]
    nvalues[index] <- dim(mod$model)[1]

    modvalues[index] <- depend
    phenovalues[index] <- indep
    index <- index + 1
  }
}

results <- data.frame(pheno=modvalues,
                      module=phenovalues,
                      b=bvalues,
                      p=pvalues,
                      n=nvalues)
results$tissue <- unlist(lapply(strsplit(results$module,"_"),function(x) { x[[2]] }))
results$modnum <- gsub("data.ME","",unlist(lapply(strsplit(results$module,"_"),function(x) { x[[1]] })))

bonfT <- 0.05/nrow(results)

results$module <- factor(results$module,levels=names(allME)[-1])

ggplot(data=subset(results,p<0.05), aes(x=pheno,y=-log10(p)*sign(b),group=modnum))+
        geom_bar(stat="identity",position="dodge")+
        geom_hline(yintercept = c(-log10(bonfT),log10(bonfT)),col="red",lty=2)+ 
        geom_hline(yintercept = c(-log10(0.05),log10(0.05)),col="blue",lty=3)+ 
        geom_text_repel(data=subset(results,p<0.05),aes(label=modnum))+
  facet_wrap(~tissue,scales="free_x")+
        theme_minimal()+
        theme(axis.text.x=element_text(angle = -45, hjust = 0))



#################################
### enrichment of modules #######
#################################
library(clusterProfiler)
library(DOSE)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(enrichplot)
library(gridExtra)
library(rms)
library(cowplot)

netcols <- cnet$colors
names(netcols) <- mapIds(org.Hs.eg.db,
                          keys=names(cnet$colors), 
                          column="ENTREZID",
                          keytype="ENSEMBL",
                          multiVals="first")

ncs1 <- netcols[which(!is.na(names(netcols)))]


ncslist1 <- lapply(unique(ncs1),function(module){
  names(ncs1)[which(ncs1==module)]
})
names(ncslist1) <- unique(ncs1)
ncslist1[["0"]] <- NULL

ck1 <- compareCluster(geneCluster = ncslist1,
                      fun = "enrichGO",
                      OrgDb='org.Hs.eg.db',
                      ont="BP")
dotplot(ck1)


############################################
### correlate DLPFC and mono modules #######
############################################
library(reshape2)
library(ggsci)
library(ggthemes)
library(ggpubr)

datmat <- as.matrix(allME[,-1])
cormat <- rcorr(datmat)

cormeltr <- melt(cormat$r)
cormeltp <- melt(cormat$P)
bothmelt <- merge(cormeltr,cormeltp,by=c("Var1","Var2"))
names(bothmelt) <- c("x","y","cor","p")

bms <- bothmelt[which(bothmelt$x!=bothmelt$y),]
bms <- bms[grep("mono",bms$x),]
bms <- bms[-grep("mono",bms$y),]

plotlabs <- ifelse(bms$p<0.01,
                   round(-log10(bms$p)),
                   NA)

ggplot(data=bms,aes(y=y,x=x,fill=cor))+
  geom_tile(color = "white")+
  scale_fill_gradient2_tableau()+
  geom_text(label=plotlabs)+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()

