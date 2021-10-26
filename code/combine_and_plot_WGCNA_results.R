### analyze WGCNA output
library(WGCNA)
library(rms)
library(broom)
library(dplyr)
library(ggdendro)
library(reshape2)
library(ggalluvial)
library(gplots)

load("/Users/dfelsky/Documents/data/all_genes_ensembl.RData")

### overlap
dat_mono <- readRDS("output/WGCNA/mono/input_datExpr.rds")
colnames(dat_mono) <- all_genes$external_gene_name[match(colnames(dat_mono),all_genes$ensembl_gene_id)]
dat_mono <- dat_mono[,which(is.na(colnames(dat_mono))==F)]
dat_dlpfc <- readRDS("output/WGCNA/dlpfc/input_adjustedforcells_datExpr.rds")

net_mono <- readRDS("output/WGCNA/mono/net_PAM.rds")
names(net_mono$colors) <- all_genes$external_gene_name[match(names(net_mono$colors),all_genes$ensembl_gene_id)]
net_mono$colors <- net_mono$colors[which(is.na(names(net_mono$colors))==F)]
net_dlpfc <- readRDS("output/WGCNA/dlpfc/net_PAM_celltypes.rds")

#### get common genes:
commongenes <- intersect(colnames(dat_mono),colnames(dat_dlpfc))
dat1common <- dat_mono[,commongenes]
dat2common <- dat_dlpfc[,commongenes]

colorh1 <- net_mono$colors[commongenes]
colorh2 <- net_dlpfc$colors[commongenes]

ME1common <- net_mono$MEs[,which(gsub("ME","",colnames(net_mono$MEs)) %in% unique(colorh1))]
ME2common <- net_dlpfc$MEs[,which(gsub("ME","",colnames(net_dlpfc$MEs)) %in% unique(colorh2))]

#### module overlaps
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

tiff("paper/supp_figures/SuppFig3_overlap_heatmap.tif",w=8,h=3,units="in",res=300)
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

###################################################
##### module - trait associations, monocytes & DLPFC
ROSmaster <- readRDS("input/ROSmaster_TWAS_input.rds")
load("/Users/dfelsky/Documents/data/all_genes_ensembl.RData")

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

#test
md$smoke <- md$smoking
md$smoke[which(md$smoke==2)] <- NA

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


########## WHERE ARE SIG GENES IN MODULES?
################ plot distribution of significant effect genes in gene modules (kME?)
alltt <- readRDS("output/all_TWAS_results_wideformat.rds")
net <- readRDS("output/WGCNA/mono/net_PAM.rds")
mono <- readRDS("output/WGCNA/mono/input_datExpr.rds")
modlist <- data.frame(module=net$colors,gene=names(net$colors))

allkme <- WGCNA::signedKME(mono,net$MEs,corFnc="bicor",corOptions="maxPOutliers=0.05")
allkme$gene <- rownames(allkme)

allinfo <- merge(alltt,allkme,by="gene",all.y=T)
allinfo <- merge(allinfo,modlist,by="gene")

subset(allinfo, module=="tan") %>%
ggplot(aes(y=kMEtan,x=-log10(MONO_BLOOD_P.Value_mf3123)*sign(MONO_BLOOD_t_mf3123)))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_minimal()

# try chi-squared test for yes/no significant for MF3123, stratified by module







mfgenes <- alltt$gene[which(alltt$MONO_BLOOD_adj.P.Val_mf3123<0.05)]

cormat <- rcorr(mono[,mfgenes])
rmat <- reshape2::melt(cormat$r)
names(rmat)[3] <- "r"
pmat <- reshape2::melt(cormat$P)
names(pmat)[3] <- "p"
allmat <- cbind(rmat,pmat)


#######
library(SuperExactTest)
library(reshape2)
library(ggpubr)
library(ggthemes)
library(ggsci)
library(cowplot)

siggenes <- lapply(ttall_unlabelled$monocyte_blood, function(patho) {
    patho$gene[which(patho$adj.P.Val < 0.05)]
})

moddefs <- lapply(unique(net$colors), function(mod) {
  names(net$colors[which(net$colors==mod)])
})
names(moddefs) <- unique(net$colors)


hyperp <- lapply(siggenes, function(patho){
  unlist(lapply(moddefs, function(module){
    testlist <- list(patho,module)
    supres <- supertest(testlist,n=length(net$colors))
    as.numeric(supres$P.value["11"])
  }))
})

hypero <- lapply(siggenes, function(patho){
  unlist(lapply(moddefs, function(module){
    testlist <- list(patho,module)
    supres <- supertest(testlist,n=length(net$colors))
    as.numeric(supres$overlap.sizes["11"])
  }))
})

hyperp2 <- do.call(rbind,hyperp)
hypero2 <- do.call(rbind,hypero)

hyperp3 <- reshape2::melt(hyperp2,value.name = "p")
hypero3 <- reshape2::melt(hypero2,value.name = "overlap")

hyperall <- merge(hyperp3,hypero3,by=c("Var1","Var2"))

ggplot(data=hyperall, aes(y=Var1,x=Var2,fill=-log10(p)))+
  geom_tile()+
  scale_fill_continuous_tableau()+
  geom_text(data=subset(hyperall, p<0.01),aes(label=overlap),size=2)+
  theme_clean()


############################ SUPPLEMENTARY FIGURE SHOWING PC1 LOADINGS OFR EACH MODULE VS. T-STATS FROM DE ANALYSIS
##### Plot of PC1 (eigengene) loadings per module against PAM effect Sizes
mono_expr_ensg <- readRDS("output/WGCNA/mono_blood/input_datExpr.rds")
mono_mf3123_genes <- ttm$mf3123$gene[which(ttm$mf3123$adj.P.Val<0.05)]
netm <- readRDS("output/WGCNA/mono/net_PAM.rds")
ttwide <- readRDS("output/all_TWAS_results_wideformat.rds")

# get PC1 loadings per gene
moddefsmono <- lapply(unique(netm$colors), function(mod) {
  modgenes <- names(netm$colors[which(netm$colors==mod)])
  mono_mfpca <- prcomp(scale(mono_expr_ensg[,modgenes]))
  return(mono_mfpca)
})
names(moddefsmono) <- unique(netm$colors)

alldatmono <- lapply(unique(netm$colors), function(module) {
  pc1 <- moddefsmono[[module]]$rotation[,"PC1"]
  tstats <- ttwide[match(names(pc1),ttwide$gene),c("gene",grep("MONO_BLOOD_t_.*3123",names(ttwide),value=T))]
  df1 <- data.frame(module=module,
                    pc1=pc1)
  df1 <- cbind(df1,tstats)
  return(df1)
})

alldatmono2 <- do.call(rbind,alldatmono)
alldatmono2$Significant <- as.factor(ifelse(alldatmono2$gene %in% ttwide$gene[which(ttwide$MONO_BLOOD_adj.P.Val_mf3123<0.05)],"yes","no"))

pmono <- ggplot(data=subset(alldatmono2,module!="grey"), aes(x=pc1,y=MONO_BLOOD_t_mf3123))+
  geom_point(aes(col=Significant))+
  geom_hline(yintercept = 0, col="red",lty=2)+
  geom_vline(xintercept=0,col="red",lty=2)+
  facet_wrap(~module)+
  scale_color_jco()+
  labs(x="PC1 Loadings",y="Differential Expression t-statistic",title="Monocyte module PC1 loadings vs. MF PAM effect sizes")+
  theme_minimal()


###### compare to DLPFC as example
##### Plot of PC1 (eigengene) loadings per module against PAM effect Sizes
ttwide <- readRDS("output/all_TWAS_results_wideformat.rds")
dlpfc_expr <- readRDS("output/WGCNA/dlpfc/input_adjustedforcells_datExpr.rds")
dlpfc_nft_genes <- ttwide$hugo[which(ttwide$DLPFC_cells_adj.P.Val_nft_sqrt<0.05)]
netd <- readRDS("output/WGCNA/dlpfc/net_PAM_celltypes.rds")

# get PC1 loadings per gene
moddefsdlpfc <- lapply(unique(netd$colors), function(mod) {
  modgenes <- names(netd$colors[which(netd$colors==mod)])
  dlpfc_nftpca <- prcomp(scale(dlpfc_expr[,modgenes]))
  return(dlpfc_nftpca)
})
names(moddefsdlpfc) <- unique(netd$colors)

alldatdlpfc <- lapply(unique(netd$colors), function(module) {
  pc1 <- moddefsdlpfc[[module]]$rotation[,"PC1"]
  tstats <- ttwide[match(names(pc1),ttwide$hugo),c("hugo",grep("DLPFC_cells_t_.*nft_sqrt",names(ttwide),value=T))]
  df1 <- data.frame(module=module,
                    pc1=pc1)
  df1 <- cbind(df1,tstats)
  return(df1)
})

alldatdlpfc2 <- do.call(rbind,alldatdlpfc)
alldatdlpfc2$Significant <- as.factor(ifelse(alldatdlpfc2$hugo %in% ttwide$hugo[which(ttwide$DLPFC_cells_adj.P.Val_nft_sqrt<0.05)],"yes","no"))

pdlpfc <- ggplot(data=subset(alldatdlpfc2,module!="grey"), aes(x=pc1,y=DLPFC_cells_t_nft_sqrt))+
  geom_point(aes(col=Significant))+
  geom_hline(yintercept = 0, col="red",lty=2)+
  geom_vline(xintercept=0,col="red",lty=2)+
  facet_wrap(~module)+
  scale_color_jco()+
  labs(x="PC1 Loadings",y="Differential Expression t-statistic",title="DLPFC module PC1 loadings vs. NFT effect sizes")+
  theme_minimal()


tiff(file="paper/supp_figures/SuppFig4_module-trait_lackofassociation_PC1_tstat_correlation_plots_MFPAM_NFT.tif",w=20,h=10,units="in",res=300)
print(plot_grid(pmono,pdlpfc))
dev.off()
