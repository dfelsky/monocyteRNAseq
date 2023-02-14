# reviewer comments
# Reviewer 2, comment 1. Microglialness of correlated genes between monocyte and brain
setwd("paper/FINAL/nature_communications/revision1/mglia_gene_sets/")

hm1 <- subset(human.master3, groups %in% c("Microglia","Microglia Development"))
hm2 <- subset(hm1, Species=="human")

mgliasets <- lapply(unique(hm2$listname),function(x) {
  hm2$ensembl_gene_id[which(hm2$listname==x)]
})
names(mgliasets) <- unique(hm2$listname)

ct <- readRDS("output/cross_tissue_correlation_results.rds")
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

### degNorm
if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# The following initializes usage of Bioc devel
BiocManager::install(version='devel')

BiocManager::install("DegNorm")