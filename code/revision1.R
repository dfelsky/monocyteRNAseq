library(rms)
library(reshape2)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(ggsci)
library(ggpubr)
library(cowplot)

# reviewer comments
# Reviewer 2, comment 1. Microglialness of correlated genes between monocyte and brain
setwd("paper/FINAL/nature_communications/revision1/mglia_gene_sets/")
load("Mouse_Human_GenelistDatabaseAugust2021.RData")

hm1 <- subset(human.master3, groups %in% c("Microglia","Microglia Development"))
hm2 <- subset(hm1, Species=="human")

mgliasets <- lapply(unique(hm2$listname),function(x) {
  hm2$ensembl_gene_id[which(hm2$listname==x)]
})
names(mgliasets) <- unique(hm2$listname)

ct <- readRDS("../../../../../output/cross_tissue_correlation_results.rds")
ha <- read.csv("humi_aged.csv")

## custom enrichment
library(SuperExactTest)
m2g <- list(mono_posfdr=as.character(ct$gene[which(ct$cor_r>0 & ct$cor_fdr<0.05)]),
            mono_negfdr=as.character(ct$gene[which(ct$cor_r<0 & ct$cor_fdr<0.05)]),
            mono_posp=as.character(ct$gene[which(ct$cor_r>0 & ct$cor_p<0.05)]),
            mono_negp=as.character(ct$gene[which(ct$cor_r<0 & ct$cor_p<0.05)]),
            mono_allfdr=as.character(ct$gene[which(ct$cor_fdr<0.05)]),
            mono_allp=as.character(ct$gene[which(ct$cor_p<0.05)]),
            humi=ha$gene)
m2 <- c(m2g,mgliasets)

mset <- supertest(m2, n=nrow(ct),degree = 2)
sm <- summary(mset)
smt <- sm$Table[,-7]

#smt[grep("posfdr",smt$Intersections),]
ssmt <- smt[which(smt$P.value<0.05),]
ssmt[grep("mono_",ssmt$Intersections),]


###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
#### peripheral blood biomarker DGE analysis is in submit_bloodbiomarker_TWAS_forrevision1.R
# here are plots and proessing of the results

### overlapping genes for amyloid phenotypes
tb <- readRDS("output/TWAS_results_monocyte_cognition_blood_biomarkers.rds")
bio9 <- names(tb)[c(2,3,4,5,6,8,9,10,11)]

biosig <- list("CRP"=tb$log_hcrp$gene[which(tb$log_hcrp$adj.P.Val<0.05)],
                   "IL1B"=tb$log_hil1b$gene[which(tb$log_hil1b$adj.P.Val<0.05)],
                   "IL1RA"=tb$log_hil1ra$gene[which(tb$log_hil1ra$adj.P.Val<0.05)],
                   "IL10"=tb$log_hil10$gene[which(tb$log_hil10$adj.P.Val<0.05)],
                   "IL6"=tb$log_hil6$gene[which(tb$log_hil6$adj.P.Val<0.05)],
               "IL6R"=tb$log_hil6r$gene[which(tb$log_hil6r$adj.P.Val<0.05)],
               "MMP9"=tb$log_hmmp9$gene[which(tb$log_hmmp9$adj.P.Val<0.05)],
               "TNFA"=tb$log_htnfa$gene[which(tb$log_htnfa$adj.P.Val<0.05)],
               "VCAM"=tb$log_hvcam$gene[which(tb$log_hvcam$adj.P.Val<0.05)])

biostest <- supertest(biosig,n=nrow(tb$log_hcrp))
pdf(file="paper/figures/fig2_mset_biomarkers.pdf",w=5,h=4)
plot.msets(biostest,Layout = "landscape",
           keep.empty.intersections = F,
           degree=c(2,3,4),
           show.elements = F,
           margin = c(1,10,10,5),
           color.on = "black",
           sort.by = "size",
           overlap.size.cex = 0.5,
           legend.text.cex = 0.5,
           color.scale.cex = 0.5,
           cex=0.5,yfrac = 0.5)
dev.off()

library(cowplot)
p1 <- ggplot(tb$log_hcrp,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title = "CRP")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=tb$log_hcrp[1:5,],aes(label=hugo))+
  theme_minimal()
p2 <- ggplot(tb$log_hil1b,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title = "IL1B")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=tb$log_hil1b[1:5,],aes(label=hugo))+
  theme_minimal()
p3 <- ggplot(tb$log_hil1ra,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title = "IL1RA")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=tb$log_hil1ra[1:5,],aes(label=hugo))+
  theme_minimal()
p4 <- ggplot(tb$log_hil10,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title = "IL10")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=tb$log_hil10[1:5,],aes(label=hugo))+
  theme_minimal()
p5 <- ggplot(tb$log_hil6,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title = "IL6")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=tb$log_hil6[1:5,],aes(label=hugo))+
  theme_minimal()
p6 <- ggplot(tb$log_hil6r,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title = "IL6R")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=tb$log_hil6r[1:5,],aes(label=hugo))+
  theme_minimal()
p7 <- ggplot(tb$log_hmmp9,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title = "MMP9")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=tb$log_hmmp9[1:5,],aes(label=hugo))+
  theme_minimal()
p8 <- ggplot(tb$log_htnfa,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title = "TNFa")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=tb$log_htnfa[1:5,],aes(label=hugo))+
  theme_minimal()
p9 <- ggplot(tb$log_hvcam,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title = "VCAM")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=tb$log_hvcam[1:5,],aes(label=hugo))+
  theme_minimal()

pdf("paper/figures/biomarker_TWAS_volcanos_revision1.pdf",h=12,w=12)
print(plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=3))
dev.off()

png("paper/figures/biomarker_TWAS_volcanos_revision1.png",units = "in",h=12,w=12,res = 600)
print(plot_grid(p1,p2,p3,p4,p5,p6,p7,p8,p9,ncol=3))
dev.off()

# summarize results

genelists <- lapply(tb[bio9], function(x) {
  y <- x[1:20,c("hugo","t","adj.P.Val")]
  y$t <- signif(y$t,2)
  y$FDR <- signif(y$adj.P.Val,2)
  y$adj.P.Val <- NULL
  y$Symbol <- y$hugo
  y$hugo <- NULL
  y <- y[,c(3,1,2)]
  y
})

top20 <- do.call(cbind,genelists)
write.csv(top20,file="paper/FINAL/nature_communications/revision1/bloodbiomarkers/top20.csv",row.names = F)

###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
### Degnorm analysis

load("paper/FINAL/nature_communications/revision1/degnorm/res_DegNorm_monocyte_batch1.Rda")
load("paper/FINAL/nature_communications/revision1/degnorm/res_DegNorm_dlpfc_batch1.Rda")

deg.mono1 <- res_DegNorm_batch1$DI
deg.dlpfc1 <- res_DegNorm_batch3$DI

dge_mono <- readRDS("input/monocytes_filtered_only.rds")
mono.qc <- dge_mono$samples
mono.qc$projid <- as.numeric(mono.qc$projid)

degmonot <- as.data.frame(t(deg.mono1))
degmonot$projid <- as.numeric(gsub("_.*","",gsub(".bam","",rownames(degmonot))))

mono.m <- merge(mono.qc,degmonot,by="projid")

gene.names <- grep("ENSG",names(mono.m),value=T)

allcors <- apply(mono.m[,gene.names],2,function(gene){
  cor.test(gene,mono.m$MEDIAN_3PRIME_BIAS)
})


###########################################################################################
###########################################################################################
###########################################################################################
###########################################################################################
### snRNAseq CTP vs. PAM analyses
## first test PAM measures vs. each cell type and cell subtype
ctp <- read.table("paper/FINAL/nature_communications/revision1/cellproportions/subcelltype_proportion_n424.txt", header=T,check.names = F)
ctp$subcelltype2 <- NULL
ctp <- as.data.frame(t(ctp))
names(ctp) <- ctp[1,]
ctp <- as.data.frame(ctp[-1,])
celltypes <- names(ctp)
ctp$projid <- as.numeric(rownames(ctp))
ctp <- as.data.frame(apply(ctp,2,as.numeric))

# check celltype distributions
ctp_melt <- melt(ctp,id.vars = "projid")
ggplot(data=ctp_melt,aes(x=value))+
  geom_histogram()+
  facet_wrap(~variable,nrow=10,scales = "free")+
  theme_minimal()

mono <- readRDS("output/mono_dge_filtered_used_biomarkers.rds")
mp <- mono$samples

ROSmaster <- readRDS("input/ROSmaster_TWAS_input.rds")

mp_RM <- merge(ROSmaster,mp,by=c("projid","msex"))
mp_RM_ctp <- merge(mp_RM,ctp,by=c("projid"))
RM_ctp <- merge(ROSmaster,ctp,by="projid")

covs <- c("pmi","msex","age_death")
phenotypes <- c("mf3123","it3123","pput3123","vm3123")

reslist <- list()
i=1
for (phenotype in phenotypes) {
  for (celltype in celltypes) {
    form <- as.formula(paste0(phenotype,"~",celltype,"+",paste0(covs,collapse = "+")))
    form_noctp <- as.formula(paste0(phenotype,"~",paste0(covs,collapse = "+")))
    
    mod_full <- lm(data=RM_ctp, form)
    mod_full_noctp <- lm(data=RM_ctp, form_noctp)
    summod_full <- summary(mod_full)
    summod_full_noctp <- summary(mod_full_noctp)
    n_full <- nrow(mod_full$model)
    t_full <- summod_full$coefficients[celltype,"t value"]
    p_full <- summod_full$coefficients[celltype,"Pr(>|t|)"]
    r2_adj_full <- summod_full$adj.r.squared
    r2_adj_full_noctp <- summod_full_noctp$adj.r.squared
    
    mod_sub <- lm(data=mp_RM_ctp, form)
    mod_sub_noctp <- lm(data=mp_RM_ctp, form_noctp)
    summod_sub <- summary(mod_sub)
    summod_sub_noctp <- summary(mod_sub_noctp)
    n_sub <- nrow(mod_sub$model)
    t_sub <- tryCatch(summod_sub$coefficients[celltype,"t value"],error=function(e) { NA })
    p_sub <- tryCatch(summod_sub$coefficients[celltype,"Pr(>|t|)"],error=function(e) { NA })
    r2_adj_sub <- summod_sub$adj.r.squared
    r2_adj_sub_noctp <- summod_sub_noctp$adj.r.squared
    
    reslist[[i]] <- c(pheno=phenotype,
                      cell=celltype,
                      n_full=n_full,
                      t_full=t_full,
                      p_full=p_full,
                      r2_adj_full=r2_adj_full,
                      r2_adj_full_noctp=r2_adj_full_noctp,
                      n_sub=n_sub,
                      t_sub=t_sub,
                      p_sub=p_sub,
                      r2_adj_sub=r2_adj_sub,
                      r2_adj_sub_noctp=r2_adj_sub_noctp)
    
    i=i+1
  }
}

allres <- as.data.frame(do.call(rbind,reslist))
allres[,seq(3,ncol(allres))] <- apply(allres[,seq(3,ncol(allres))],2,as.numeric)
allres$cell <- as.factor(allres$cell)
allres$pheno <- as.factor(allres$pheno)

allres <- subset(allres, pheno=="mf3123")

allres$cellfull <- with(allres, reorder(cell, t_full, mean,order=T))

plotfull <- ggplot(data=allres, aes(x=-log10(p_full)*sign(t_full),y=cellfull,col=pheno))+
  geom_vline(xintercept = c(0,-log10(0.05),log10(0.05)),lty=c(1,2,2),col=c("blue","red","red"))+
  geom_text_repel(data=subset(allres,p_full<0.05),aes(label=cellfull),show.legend = F)+
  geom_point()+
  theme_minimal()+
  scale_color_jco()+
  scale_y_discrete(limits = levels(allres$cellfull))+
  labs(title = "Association of PAM Phenotypes with snRNAseq Cell Subtypes (n=74-80)",
       subtitle= "Maximal phenotypic overlap (covariates: PMI, age, sex)",
       y="Cell subtype proportion estimate (snRNAseq)",
       col="PAM Phenotype",
       x="Signed -log10(p-value)")
  
# test for sample that overlaps maximally with phenotype data 
allres$cellsub <- with(allres, reorder(cell, t_sub, mean,order=T))

plotsub <- ggplot(data=allres, aes(x=-log10(p_sub)*sign(t_sub),y=cellsub,col=pheno))+
  geom_vline(xintercept = c(0,-log10(0.05),log10(0.05)),lty=c(1,2,2),col=c("blue","red","red"))+
  geom_text_repel(data=subset(allres,p_sub<0.05),aes(label=cellsub),show.legend = F)+
  geom_point()+
  theme_minimal()+
  scale_color_jco()+
  scale_y_discrete(limits = levels(allres$cellsub))+
  labs(title = "Association of PAM Phenotypes with snRNAseq Cell Subtypes (n=28-31)",
       subtitle= "Overlapping sample with monocyte RNAseq (covariates: PMI, age, sex)",
       y="Cell subtype proportion estimate (snRNAseq)",
       col="PAM Phenotype",
       x="Signed -log10(p-value)")

pdf(file="paper/figures/revision_PAM_association_with_subtypes_snRNAseq_V1.pdf",h=12,w=12)
print(plot_grid(plotfull,plotsub,nrow=1))
dev.off()

## second, test top mono genes that are associated with PAM against cell types and subtypes

ttall_unlabelled <- readRDS("output/all_TWAS_results.rds")

### name files for plotting
pathnames <- c("Neuritic plaques","Diffuse plaques","Total AB","PHF tau","NFT","Gross cerebral infarcts","Micro cerebral infarcts","Arteriolosclerosis","Cerebral AA","Cerebral atherosclerosis","Lewy body stage","Hippocampal sclerosis","TDP-43","PD Dx","Patho AD","PAM VM Caudate","PAM post. putamen","PAM IT","PAM MF")
cognames <- c("Episodic memory","Perceptual orientation","Perceptual speed","Semantic memory","Working memory","Global","MMSE")

ttall_labelled <- lapply(ttall_unlabelled, function(x) {
  names(x) <- c(pathnames,cognames)
  x
})

pamsig <- list(MF=ttall_unlabelled$monocyte_blood$mf3123$gene[which(ttall_unlabelled$monocyte_blood$mf3123$adj.P.Val<0.05)],
               IT=ttall_unlabelled$monocyte_blood$it3123$gene[which(ttall_unlabelled$monocyte_blood$it3123$adj.P.Val<0.05)],
               "Post. putamen"=ttall_unlabelled$monocyte_blood$pput3123$gene[which(ttall_unlabelled$monocyte_blood$pput3123$adj.P.Val<0.05)],
               "VM caudate"=ttall_unlabelled$monocyte_blood$vm3123$gene[which(ttall_unlabelled$monocyte_blood$vm3123$adj.P.Val<0.05)])

mg <- as.data.frame(t(mono$counts))
mg$projid <- rownames(mg)
mg <- mg[,c("projid",pamsig$MF)]

mpg <- merge(mp,mg,by="projid")
mpg2 <- merge(ROSmaster,mpg,by=c("projid","msex"))
mpg_rm_ctp <- merge(mpg2,ctp,by="projid")

covs <- c("pmi","msex","age_death","batch","age_draw","PCT_RIBOSOMAL_BASES")
phenotypes <- grep("ENSG",names(mpg_rm_ctp),value=T)

reslist2 <- list()
i=1
for (phenotype in phenotypes) {
  for (celltype in celltypes) {
    form <- as.formula(paste0(phenotype,"~",celltype,"+",paste0(covs,collapse = "+")))
    form_noctp <- as.formula(paste0(phenotype,"~",paste0(covs,collapse = "+")))
    
    mod_full <- lm(data=mpg_rm_ctp, form)
    mod_full_noctp <- lm(data=mpg_rm_ctp, form_noctp)
    summod_full <- summary(mod_full)
    summod_full_noctp <- summary(mod_full_noctp)
    n_full <- nrow(mod_full$model)
    t_full <- summod_full$coefficients[celltype,"t value"]
    p_full <- summod_full$coefficients[celltype,"Pr(>|t|)"]
    r2_adj_full <- summod_full$adj.r.squared
    r2_adj_full_noctp <- summod_full_noctp$adj.r.squared
    
    reslist2[[i]] <- c(pheno=phenotype,
                      cell=celltype,
                      n_full=n_full,
                      t_full=t_full,
                      p_full=p_full,
                      r2_adj_full=r2_adj_full,
                      r2_adj_full_noctp=r2_adj_full_noctp)
    
    i=i+1
  }
}

load("input/all_genes_ensembl.RData")
allres2 <- as.data.frame(do.call(rbind,reslist2))
allres2[,seq(3,ncol(allres2))] <- apply(allres2[,seq(3,ncol(allres2))],2,as.numeric)
allres2$hugo <- all_genes$external_gene_name[match(allres2$pheno,all_genes$ensembl_gene_id)]

allres2 <- allres2[order(allres2$p_full,decreasing = F),]
allres2$fdr <- p.adjust(allres2$p_full,method="fdr")

allres2s <- subset(allres2,fdr < 0.05)

ggplot(data=allres2s, aes(y=hugo,x=cell,fill=t_full))+
  geom_tile()+
  geom_text(aes(label=signif(p_full,2)))+
  theme_minimal()

# merge effects on celltypes with effects on PAM
ar <- ttall_unlabelled$monocyte_blood$mf3123
ar2 <- allres2
names(ar) <- paste0(names(ar),"_mf3123")
names(ar2) <- paste0(names(ar2),"_monogene")
mres <- merge(ar,ar2,by.x="gene_mf3123",by.y="pheno_monogene")

plist <- NULL
rholist <- NULL
for (celltype in celltypes) {
  mressub <- subset(mres, cell_monogene==celltype)
  ct1 <- cor.test(mressub$t_full_monogene, mressub$t_mf3123, method="spearman")
  plist[celltype] <- ct1$p.value
  rholist[celltype] <- ct1$estimate
}
allcors <- data.frame(cell=celltypes,p=plist,r=rholist)
allcors$cell <- factor(allcors$cell,levels=allcors$cell[order(allcors$r)])

ggplot(data=allcors,aes(x=r,y=cell))+
  geom_bar(stat="identity")

ggplot(data=mres,aes(y=t_full_monogene,x=t_mf3123))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~cell_monogene,scales="free")+
  theme_minimal()

