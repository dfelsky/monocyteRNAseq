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

mono <- readRDS("output/mono_dge_filtered_used_biomarkers.rds")
mp <- mono$samples

ROSmaster <- readRDS("input/ROSmaster_TWAS_input.rds")

mp_RM <- merge(ROSmaster,mp,by="projid")
mp_RM_ctp <- merge(mp_RM,ctp,by="projid")
RM_ctp <- merge(ROSmaster,ctp,by="projid")

# test for sample that overlaps with obs from DE analysis
covs <- c()
phenotypes <- c()

for (phenotype in phenotypes) {
  for (celltype in celltypes) {
    form <- as.formula(paste0(phenotype,"~",celltype,"+",paste0(covs,collapse = "+")))
    lm(data=RM_ctp, mf_3123 ~ celltype)
  }
}

# test for sample that overlaps maximally with phenotype data 



## second, test top mono genes that are associated with PAM against cell types and subtypes






