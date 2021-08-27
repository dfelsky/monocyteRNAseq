#### module enrichment #####
library(clusterProfiler)
library(DOSE)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(enrichplot)
library(gridExtra)
library(rms)
library(data.table)
library(ggraph)

#### Monocytes
#### load module definitions
net <- readRDS("output/WGCNA/mono/net_PAM.rds")

####
modlist <- lapply(unique(net$colors), function(module) {
  ensg <- names(net$colors)[which(net$colors==module)]
  as.vector(na.omit(mapIds(org.Hs.eg.db,
         keys=ensg, 
         column="ENTREZID",
         keytype="ENSEMBL",
         multiVals="first")))
})
names(modlist) <- unique(net$colors)
modlist$grey <- NULL

ck1 <- compareCluster(geneCluster = modlist,
                      fun = "enrichGO",
                      OrgDb='org.Hs.eg.db',
                      ont="BP")
ckr <- setReadable(ck1, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
cks <- simplify(ckr)

pdf("paper/figures/Fig4_monocyte_module_enrichment_dotplot.pdf",w=8,h=6)
clusterProfiler::dotplot(cks,font.size=8)
dev.off()


