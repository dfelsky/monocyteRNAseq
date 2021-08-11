### analyze WGCNA output
library(WGCNA)

load("/Users/dfelsky/Documents/data/all_genes_ensembl.RData")

### overlap
dat_mono <- readRDS("output/WGCNA/mono/input_datExpr.rds")
colnames(dat_mono) <- all_genes$external_gene_name[match(colnames(dat_mono),all_genes$ensembl_gene_id)]
dat_mono <- dat_mono[,which(is.na(colnames(dat_mono))==F)]
dat_dlpfc <- readRDS("output/WGCNA/dlpfc/input_adjustedforcells_datExpr.rds")

net_mono <- readRDS("output/WGCNA/mono/net_noPAM.rds")
names(net_mono$colors) <- all_genes$external_gene_name[match(names(net_mono$colors),all_genes$ensembl_gene_id)]
net_mono$colors <- net_mono$colors[which(is.na(names(net_mono$colors))==F)]
net_dlpfc <- readRDS("output/WGCNA/dlpfc/net_noPAM_celltypes.rds")

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
  ifelse(x < 0.0005,1,0)
})

mod_overlap2 <- overlapTable(colorh1,colorh2)

plotdat1 <- melt(overlapYN)
names(plotdat1)[3] <- "sig"
plotdat2 <- melt(mod_overlap$PvaluesHypergeo)
names(plotdat2)[3] <- "pvalue"
plotdat3 <- melt(mod_overlap2$countTable)
plotdat3$Var1 <- paste0("mono_",plotdat3$Var1)
plotdat3$Var2 <- paste0("dlpfc_",plotdat3$Var2)
names(plotdat3)[3] <- "n"

plotdat <- merge(plotdat1,plotdat2)
plotdat <- merge(plotdat,plotdat3)

plotdat_sub <- subset(plotdat, Var1!="mono_grey" & Var2!="dlpfc_grey")

ggplot(data=plotdat_sub,aes(y=Var2,x=Var1,fill=-log10(pvalue)))+
  geom_tile()+
  scale_fill_gradient(low = "white",high="red")+
  geom_text(aes(label=n),col="black",size=3)+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))

############
###################################################
##### module - trait associations, monocytes
library(broom)

ROSmaster <- readRDS("input/ROSmaster_TWAS_input.rds")
net <- readRDS("output/WGCNA/mono/net_noPAM.rds")

### trait-ME associations
ME <- net$MEs
ME$projid <- rownames(net$MEs)
md <- merge(ME,ROSmaster,by="projid")

egenelist <- names(net$MEs)[which(names(net$MEs) %nin% "MEgrey")] 
cont.phenos <- c("plaq_n_sqrt","plaq_d_sqrt","amyloid_sqrt","tangles_sqrt","nft_sqrt","arteriol_scler","caa_4gp","cvda_4gp2","dlbdx","tdp_stage4","vm3123","pput3123","it3123","mf3123") 
covars <- c("age_death","pmi","msex")

index <- 1
pvalues <- NULL
bvalues <- NULL
nvalues <- NULL
phenovalues <- NULL
egenevalues <- NULL
for (pheno in cont.phenos) {
  for (egene in egenelist) {
    form <- formula(paste(pheno,"~",egene,"+",paste0(covars,collapse = " + ")))
    mod <- lm(data=md, form)
    pvalues[index] <- tidy(mod)$p.value[2]
    bvalues[index] <- coef(mod)[2]
    nvalues[index] <- dim(mod$model)[1]
    egenevalues[index] <- egene
    phenovalues[index] <- pheno
    index <- index + 1
  }
}

results <- data.frame(pheno=phenovalues,
                      egene=egenevalues,
                      b=bvalues,
                      p=pvalues,
                      n=nvalues,
                      signedp=-log10(pvalues)*sign(bvalues))

bonfT <- 0.05/nrow(results)

library(ggdendro)
library(reshape2)

rescast <- dcast(results,pheno ~ egene,value.var = "signedp")
rownames(rescast) <- rescast$pheno
rescast$pheno <- NULL
rescast <- as.matrix(rescast)

phenodend <- as.dendrogram(hclust(dist(rescast)))
pheno.order <- order.dendrogram(phenodend)

egenedend <- as.dendrogram(hclust(dist(t(rescast))))
egene.order <- order.dendrogram(egenedend)

results$pheno.f <- factor(x = results$pheno,
                               levels = rownames(rescast)[pheno.order], 
                               ordered = TRUE)
results$egene.f <- factor(x = results$egene,
                               levels = colnames(rescast)[egene.order], 
                               ordered = TRUE)

ggplot(data=results,aes(y=pheno.f,x=egene.f,fill=-log10(p)*sign(b)))+
  geom_tile()+
  scale_fill_gradient2()+
  #geom_text(aes(label=n),col="black",size=3)+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))
  


