# Figure 4
library(rms)
library(WGCNA)
library(gplots)
library(broom)
library(dplyr)
library(ggdendro)
library(reshape2)
library(ggalluvial)
load("input/all_genes_ensembl.RData")
ttall_unlabelled <- readRDS("output/all_TWAS_results.rds")

### Panel A. network plot of monocyte transcriptome
# source("code/network_plot_for_WGCNA_monocytes.R")

#### Panel A. Barchart
datExpr <- readRDS("output/WGCNA/mono/input_datExpr.rds")
net_mono <- readRDS("output/WGCNA/mono/net_PAM.rds")
net_dlpfc <- readRDS("output/WGCNA/dlpfc/net_PAM_celltypes.rds")

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


#### Panel B. Dendrograms
tiff("paper/figures/Figure4_dendrogram_mono_PAM.tif",w=4,h=2,units="in",res=600)
plotDendroAndColors(net_mono$dendrograms[[1]],
                    net_mono$colors,
                    dendroLabels = FALSE,
                    hang = 0.02,
                    addGuide = F,
                    guideHang = 0.05,
                    main = "Monocytes")
dev.off()

tiff("paper/figures/Figure4_dendrogram_dlpfc_PAM.tif",w=4,h=2,units="in",res=600)
plotDendroAndColors(net_dlpfc$dendrograms[[1]],
                    net_dlpfc$colors,
                    dendroLabels = FALSE,
                    hang = 0.02,
                    addGuide = F,
                    guideHang = 0.05,
                    main = "DLPFC")
dev.off()


#### Panel B. Module overlap & Alluvial Plot
# identify module overlaps
dat_mono <- readRDS("output/WGCNA/mono/input_datExpr.rds")
colnames(dat_mono) <- all_genes$external_gene_name[match(colnames(dat_mono),all_genes$ensembl_gene_id)]
dat_mono <- dat_mono[,which(is.na(colnames(dat_mono))==F)]
dat_dlpfc <- readRDS("output/WGCNA/dlpfc/input_adjustedforcells_datExpr.rds")

net_mono <- readRDS("output/WGCNA/mono/net_PAM.rds")
names(net_mono$colors) <- all_genes$external_gene_name[match(names(net_mono$colors),all_genes$ensembl_gene_id)]
net_mono$colors <- net_mono$colors[which(is.na(names(net_mono$colors))==F)]
net_dlpfc <- readRDS("output/WGCNA/dlpfc/net_PAM_celltypes.rds")

# get common genes:
commongenes <- intersect(colnames(dat_mono),colnames(dat_dlpfc))
dat1common <- dat_mono[,commongenes]
dat2common <- dat_dlpfc[,commongenes]

colorh1 <- net_mono$colors[commongenes]
colorh2 <- net_dlpfc$colors[commongenes]

ME1common <- net_mono$MEs[,which(gsub("ME","",colnames(net_mono$MEs)) %in% unique(colorh1))]
ME2common <- net_dlpfc$MEs[,which(gsub("ME","",colnames(net_dlpfc$MEs)) %in% unique(colorh2))]

# module overlaps
mod_overlap <- overlapTableUsingKME(dat1 = dat1common,
                                    dat2 = dat2common,
                                    MEs1 = ME1common,
                                    MEs2 = ME2common,
                                    name1 = "mono",
                                    name2 = "dlpfc",
                                    colorh1 = colorh1,
                                    colorh2 = colorh2,
                                    omitGrey = TRUE)


overlapYN <- apply(mod_overlap$PvaluesHypergeo,2,function(x) {
  y <- p.adjust(x)
  ifelse(y < 0.05,1,0)
})

mod_overlap2 <- overlapTable(colorh1,colorh2)

plotdat1 <- melt(overlapYN)
plotdat1$Var1 <- gsub(".*_","",plotdat1$Var1)
plotdat1$Var2 <- gsub(".*_","",plotdat1$Var2)
names(plotdat1) <- c("Monocyte","DLPFC","sig")
plotdat2 <- melt(mod_overlap$PvaluesHypergeo)
plotdat2$Var1 <- gsub(".*_","",plotdat2$Var1)
plotdat2$Var2 <- gsub(".*_","",plotdat2$Var2)
names(plotdat2) <- c("Monocyte","DLPFC","pvalue")
plotdat3 <- melt(mod_overlap2$countTable)
names(plotdat3) <- c("Monocyte","DLPFC","n")

plotdat <- merge(plotdat1,plotdat2)
plotdat <- merge(plotdat,plotdat3)

plotdat_sub <- subset(plotdat, Monocyte!="grey" & DLPFC!="grey")

phenodend <- as.dendrogram(hclust(dist(mod_overlap$PvaluesHypergeo)))
pheno.order <- order.dendrogram(phenodend)

egenedend <- as.dendrogram(hclust(dist(t(mod_overlap$PvaluesHypergeo))))
egene.order <- order.dendrogram(egenedend)

plotdat_sub$Monocyte <- factor(x = plotdat_sub$Monocyte,
                               levels = unique(colorh1)[pheno.order], 
                               ordered = TRUE)
plotdat_sub$DLPFC <- factor(x = plotdat_sub$DLPFC,
                            levels = unique(colorh2)[egene.order], 
                            ordered = TRUE)

# heatmap of overlap likely used for supplement
pdf("paper/figures/Fig4_overlap_heatmap.pdf",w=8,h=3)
ggplot(data=plotdat_sub,aes(y=Monocyte,x=DLPFC,fill=-log10(pvalue)))+
  geom_tile()+
  scale_fill_gradient(low = "white",high="red")+
  geom_text(data=subset(plotdat_sub,sig==1),aes(label=n),col="black",size=3)+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))
dev.off()

############ alluvial plot for module membership
modsum <- lapply(mod_overlap$OverlappingGenes,function(x) { length(x) })
monomods <- unlist(lapply(strsplit(names(modsum),split = "_"),function(x) { x[2]}))
dlpfcmods <- unlist(lapply(strsplit(names(modsum),split = "_"),function(x) { x[4]}))

alludat <- data.frame(mono=monomods,dlpfc=dlpfcmods,Freq=unlist(modsum))
alludat <- subset(alludat,mono %nin% "grey" & dlpfc %nin% "grey" & Freq>0)
rownames(alludat) <- NULL

allulode <- to_lodes_form(alludat,axes=1:2)

modlist <- unique(c(alludat$mono,alludat$dlpfc))
hexcols <- col2hex(modlist)
names(hexcols) <- modlist

pdf("paper/figures/Fig4_alluvial.pdf",w=5,h=5)
ggplot(allulode, aes(x=x,stratum=stratum,alluvium=alluvium,y=Freq)) +
  geom_alluvium(width=1/12) +
  geom_stratum(aes(fill=stratum),width = 1/12, color = "white") +
  scale_fill_manual(values=hexcols)+
  theme_void()+
  theme(legend.position = "none")
dev.off()

#### Panel C. Module enrichment
# source("code/enrichment_WGCNA.R")

#### Panel D. module - trait associations
ROSmaster <- readRDS("input/ROSmaster_TWAS_input.rds")
load("/Users/dfelsky/Documents/data/all_genes_ensembl.RData")

###############################
###############################
###############################
### trait-ME associations
MEm <- net_mono$MEs
names(MEm) <- paste0("mono_",names(MEm))
MEm$projid <- rownames(net_mono$MEs)
MEd <- net_dlpfc$MEs
names(MEd) <- paste0("dlpfc_",names(MEd))
MEd$projid <- rownames(net_dlpfc$MEs)

MEboth <- merge(MEm,MEd,by="projid",all=T)

## read in dge_filtered for age_draw
dge_filtered <- readRDS("input/monocytes_filtered_only.rds")
RM <- merge(ROSmaster,dge_filtered$samples[,c("projid","age_draw")],by="projid",all.x=T)

md <- merge(MEboth,RM,by="projid",all.x=T)

####################################
####################################
#### define outcomes and covariates
####################################
covars.mono.pathology <- c("msex","age_death","pmi","age_draw")
covars.mono.cognition <- c("msex","age_draw","educ")

#covars.mono.pathology.blood <- c("msex","age_death","pmi","age_draw")
#covars.mono.cognition.blood <- c("msex","age_draw","educ")
covars.mono.pathology.blood <- c("msex","age_death","pmi","age_draw","hemoglbn_at_draw","mchc_at_draw","mcv_at_draw","platelet_at_draw","wbc_at_draw","fasting.f","hemotologic_rx_at_draw")
covars.mono.cognition.blood <- c("msex","age_draw","educ","hemoglbn_at_draw","mchc_at_draw","mcv_at_draw","platelet_at_draw","wbc_at_draw","fasting.f","hemotologic_rx_at_draw")

covars.dlpfc.pathology <- c("msex","pmi","age_death")
covars.dlpfc.cognition <- c("educ","msex","age_death","age_at_visit_at_lastvisit")

indepvec.pathology <- c("plaq_n_sqrt","plaq_d_sqrt","amyloid_sqrt","tangles_sqrt","nft_sqrt","ci_num2_gct","ci_num2_mct","arteriol_scler","caa_4gp","cvda_4gp2","dlbdx","hspath_any","tdp_stage4","parkdx","pathoAD","vm3123","pput3123","it3123","mf3123") 

indepvec.mono.cognition <- grep("cogn|mmse30",grep("at_draw",names(ROSmaster),value=T),value=T)
indepvec.dlpfc.cognition <- grep("cogn|mmse30",grep("lastvisit",names(ROSmaster),value=T),value=T)

pathnames <- c("Neuritic plaques","Diffuse plaques","Total AB","PHF tau","NFT","Gross cerebral infarcts","Micro cerebral infarcts","Arteriolosclerosis","Cerebral AA","Cerebral atherosclerosis","Lewy body stage","Hippocampal sclerosis","TDP-43","PD Dx","Patho AD","PAM VM Caudate","PAM post. putamen","PAM IT","PAM MF")
cognames <- c("Episodic memory","Perceptual orientation","Perceptual speed","Semantic memory","Working memory","Global","MMSE")

varnameindex <- data.frame(var=c(indepvec.pathology,indepvec.mono.cognition,indepvec.dlpfc.cognition),
                           name=c(pathnames,cognames,cognames))

nocovars <- FALSE
#######################################

index <- 1
pvalues <- NULL
bvalues <- NULL
nvalues <- NULL
phenovalues <- NULL
egenevalues <- NULL
tissuevals <- NULL
for(vartype in c("path","cog")) {
  for (tissue in c("monocyte_blood","dlpfc")) {
    if (tissue=="monocyte") {
      if (vartype=="path") {
        indepvec <- indepvec.pathology
        covars <- covars.mono.pathology
      } else {
        indepvec <- indepvec.mono.cognition
        covars <- covars.mono.cognition
      }
      egenelist <- grep("MEgrey",grep("mono_",names(MEm),value=T),value=T,invert = T)
    } else if (tissue=="monocyte_blood") {
      if (vartype=="path") {
        indepvec <- indepvec.pathology
        covars <- covars.mono.pathology.blood
      } else {
        indepvec <- indepvec.mono.cognition
        covars <- covars.mono.cognition.blood
      }
      egenelist <- grep("MEgrey",grep("mono_",names(MEm),value=T),value=T,invert = T)
    } else {
      if (vartype=="path") {
        indepvec <- indepvec.pathology
        covars <- covars.dlpfc.pathology
      } else {
        indepvec <- indepvec.dlpfc.cognition
        covars <- covars.dlpfc.cognition
      }
      egenelist <- grep("MEgrey",grep("dlpfc_",names(MEd),value=T),value=T,invert = T)
    }
    
    tablist <- apply(md[,indepvec],2,table) %>%
      lapply(length) %>%
      as.numeric()
    cont.index <- tablist > 2 #TRUE if numeric indepvec
    phenindex <- 1
    
    for (pheno in indepvec) {
      if (cont.index[phenindex]==T) {
        for (egene in egenelist) {
          if (nocovars==TRUE) {
            form <- formula(paste(pheno,"~",egene))
          } else {
            form <- formula(paste(pheno,"~",egene,"+",paste0(covars,collapse = " + ")))
          }
          mod <- lm(data=md, form)
          pvalues[index] <- tidy(mod)$p.value[2]
          bvalues[index] <- coef(mod)[2]
          nvalues[index] <- dim(mod$model)[1]
          egenevalues[index] <- egene
          phenovalues[index] <- pheno
          tissuevals[index] <- tissue
          index <- index+1
        }
      } else {
        for (egene in egenelist) {
          if (nocovars==TRUE) {
            form <- formula(paste(pheno,"~",egene))
          } else {
            form <- formula(paste(pheno,"~",egene,"+",paste0(covars,collapse = " + ")))
          }
          mod <- lrm(data=md, form)
          pvalues[index] <- anova(mod)[1,3]
          bvalues[index] <- coef(mod)[2]
          nvalues[index] <- mod$stats[1]
          egenevalues[index] <- egene
          phenovalues[index] <- pheno
          tissuevals[index] <- tissue
          index <- index+1
        }
      }
      phenindex <- phenindex + 1
    }
  }
}

results <- data.frame(tissue=tissuevals,
                      pheno=phenovalues,
                      egene=egenevalues,
                      b=bvalues,
                      p=pvalues,
                      n=nvalues,
                      signedp=-log10(pvalues)*sign(bvalues))

######################################
###### heatplot
hplotlist <- list()
for (tissue in c("monocyte_blood","dlpfc")) {
  
  ressub <- results[which(results$tissue==tissue),]
  ressub$pheno <- varnameindex$name[match(ressub$pheno,varnameindex$var)]
  #ressub$pheno <- paste0(ressub$pheno," (",ressub$n,")")
  ressub$egene <- gsub(".*ME","",ressub$egene)
  
  if (tissue=="dlpfc") {
    moddefs <- data.frame(module=names(table(net_dlpfc$colors)),
                          ngenes=as.character(table(net_dlpfc$colors)))
  } else {
    moddefs <- data.frame(module=names(table(net_mono$colors)),
                          ngenes=as.character(table(net_mono$colors)))
  }
  
  ressub$modn <- moddefs$ngenes[match(ressub$egene,moddefs$module)]
  ressub$egene <- paste0(ressub$egene," (",ressub$modn,")")
  
  bonfT <- 0.05/nrow(ressub)
  rescast <- dcast(ressub, pheno ~ egene,value.var = "signedp")
  
  rownames(rescast) <- rescast$pheno
  rescast$pheno <- NULL
  rescast <- as.matrix(rescast)
  
  # phenodend <- as.dendrogram(hclust(dist(rescast)))
  # pheno.order <- order.dendrogram(phenodend)
  # 
  egenedend <- as.dendrogram(hclust(dist(t(rescast))))
  egene.order <- order.dendrogram(egenedend)
  
  # ressub$pheno.f <- factor(x = ressub$pheno,
  #                                levels = rownames(rescast)[pheno.order], 
  #                                ordered = TRUE)
  ressub$egene.f <- factor(x = ressub$egene,
                           levels = colnames(rescast)[egene.order], 
                           ordered = TRUE)
  
  ressub$pheno.f <- factor(x = ressub$pheno,
                           levels = c(pathnames,cognames), 
                           ordered = TRUE)
  
  
  hplotlist[[tissue]] <- ggplot(data=ressub,aes(y=pheno.f,x=egene.f,fill=-log10(p)*sign(b)))+
    geom_tile(show.legend = F)+
    scale_fill_gradient2(low="cornflowerblue",high="brown1",mid = "white")+
    geom_text(data=subset(ressub, p < 0.05 & p > bonfT),label="*",col="black",size=3)+
    geom_text(data=subset(ressub, p < bonfT),label="**",col="black",size=3)+
    labs(y="Phenotype",x="Module",title=tissue)+
    theme_minimal()+
    theme(axis.text.x=element_text(angle = -45, hjust = 0))
}

pdf("paper/figures/Figure4_module_effects_heatmap.pdf",h=5,w=14)
plot_grid(hplotlist$monocyte_blood,hplotlist$dlpfc,nrow=1,rel_widths = c(1,2))
dev.off()



#### plot for significant gene effects distribution in modules
moddist_siggenes <- lapply(ttall_unlabelled$monocyte_blood[c(2,8,9,16:19)],function(pheno) {
  sapply(c("up","down"), function(updown) {
    if (updown=="up") {
      sgenes <- pheno$hugo[which(pheno$adj.P.Val < 0.05 & pheno$t > 0 )]
      table(net_mono$colors[match(sgenes,names(net_mono$colors))])
    } else {
  sgenes <- pheno$hugo[which(pheno$adj.P.Val < 0.05 & pheno$t < 0)]
  table(net_mono$colors[match(sgenes,names(net_mono$colors))])
    }
  })
})






