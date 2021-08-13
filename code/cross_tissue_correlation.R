library(limma)
library(edgeR)
library(BRETIGEA)
library(ggsci)
library(ggrepel)
library(ggthemes)
library(ggpubr)
library(gridExtra)
library(reshape2)
library(rms)
library(cowplot)

load("input/all_genes_ensembl.RData")

#### cross-tissue gene correlation
mono <- readRDS("input/monocytes_techandagesex_residuals.rds")
removesubs.mono <- c("00246264","35072859","91804757","21135680","69866926","50303145")
mono <- mono[,which(colnames(mono) %nin% removesubs.mono)]

mono.blood <- readRDS("input/monocytes_techandagesexblood_residuals.rds")
mono.blood <- mono.blood[,which(colnames(mono.blood) %nin% removesubs.mono)]

dlpfc <- readRDS("input/dlpfc_techandagesex_residuals.rds")
removesubs.dlpfc <- c("20634274","23690880","50103967","32697960","20886846","43596435","16322424","30544882","74522154","21246218","74284255","11390174","69924281","85578107","11475462","00402800","01797756","24141372")
dlpfc <- dlpfc[,which(colnames(dlpfc) %nin% removesubs.dlpfc)]

### correct for cell types in DLPFC
genesymbols <- all_genes$external_gene_name[match(rownames(dlpfc),all_genes$ensembl_gene_id)]
rownames(dlpfc) <- genesymbols
dlpfc <- dlpfc[!is.na(genesymbols),]
dlpfc.adj <- adjustBrainCells(dlpfc,species = "human")
dlpfc <- dlpfc.adj$expression
rownames(dlpfc) <- all_genes$ensembl_gene_id[match(rownames(dlpfc),all_genes$external_gene_name)]
dlpfc.celltypes <- dlpfc

# commonlists no blood
commonsubs <- intersect(colnames(mono),colnames(dlpfc.celltypes))
commongenes <- intersect(rownames(mono),rownames(dlpfc.celltypes))

# commonlists with blood
commonsubs.blood <- intersect(colnames(mono.blood),colnames(dlpfc.celltypes))
commongenes.blood <- intersect(rownames(mono.blood),rownames(dlpfc.celltypes))

# new matched dfs no blood
m2 <- t(mono[commongenes,commonsubs])
d2 <- t(dlpfc.celltypes[commongenes,commonsubs])
# new matched dfs with blood
m2.blood <- t(mono.blood[commongenes.blood,commonsubs.blood])
d2.blood <- t(dlpfc.celltypes[commongenes.blood,commonsubs.blood])

# correlations no blood
allcors <- lapply(commongenes, FUN=function(x) { cor.test(m2[,x],d2[,x],method = "spearman") })
allcor_r <- as.numeric(unlist(lapply(allcors,FUN=function(x) { x$estimate })))
allcor_p <- as.numeric(unlist(lapply(allcors,FUN=function(x) { x$p.value })))
# correlations with blood
allcors.blood <- lapply(commongenes.blood, FUN=function(x) { cor.test(m2.blood[,x],d2.blood[,x],method = "spearman") })
allcor_r.blood <- as.numeric(unlist(lapply(allcors.blood,FUN=function(x) { x$estimate })))
allcor_p.blood <- as.numeric(unlist(lapply(allcors.blood,FUN=function(x) { x$p.value })))

# allinfo no blood
allinfo <- data.frame(gene=commongenes,
                      cor_r=allcor_r,
                      cor_p=allcor_p)
allinfo$cor_fdr <- p.adjust(allinfo$cor_p,method="BH")
allinfo$hugo <- all_genes$external_gene_id[match(allinfo$gene,all_genes$ensembl_gene_id)]
allinfo$cor_sig <- ifelse(allinfo$cor_p<0.05,1,0)
allinfo$cor_sig <- as.factor(ifelse(allinfo$cor_fdr<0.05,2,allinfo$cor_sig))
# allinfo with blood
allinfo.blood <- data.frame(gene=commongenes.blood,
                            cor_r_blood=allcor_r.blood,
                            cor_p_blood=allcor_p.blood)
allinfo.blood$cor_fdr_blood <- p.adjust(allinfo.blood$cor_p_blood,method="BH")
allinfo.blood$hugo <- all_genes$external_gene_id[match(allinfo.blood$gene,all_genes$ensembl_gene_id)]
allinfo.blood$cor_sig_blood <- ifelse(allinfo.blood$cor_p_blood<0.05,1,0)
allinfo.blood$cor_sig_blood <- as.factor(ifelse(allinfo.blood$cor_fdr_blood<0.05,2,allinfo.blood$cor_sig_blood))

#merge
allinfo.both <- merge(allinfo,allinfo.blood,by=c("gene","hugo"))
#saveRDS(allinfo.both, file="output/cross_tissue_correlation_results.rds")

##### enrichment analysis ranked by cross-tissue correlation
library(org.Hs.eg.db)
library(gprofiler2)
library(rrvgo)

allinfo.both <- readRDS("output/cross_tissue_correlation_results.rds")

poscor_fdr <- allinfo.both[which(allinfo.both$cor_fdr < 0.05 & allinfo.both$cor_r > 0),c("gene","cor_r")]
poscor_raw <- allinfo.both[which(allinfo.both$cor_p < 0.05 & allinfo.both$cor_r > 0),c("gene","cor_r")]
negcor_fdr <- allinfo.both[which(allinfo.both$cor_fdr < 0.05 & allinfo.both$cor_r < 0),c("gene","cor_r")]
negcor_raw <- allinfo.both[which(allinfo.both$cor_p < 0.05 & allinfo.both$cor_r < 0),c("gene","cor_r")]

raw_r_enrich_pos_fdr <- gost(query = poscor_fdr$gene[order(poscor_fdr$cor_r,decreasing = T)],organism = "hsapiens",significant = T,ordered_query = T)
raw_r_enrich_pos_raw <- gost(query = poscor_raw$gene[order(poscor_raw$cor_r,decreasing = T)],organism = "hsapiens",significant = T,ordered_query = T)
raw_r_enrich_neg_fdr <- gost(query = negcor_fdr$gene[order(negcor_fdr$cor_r,decreasing = F)],organism = "hsapiens",significant = T,ordered_query = T)
raw_r_enrich_neg_raw <- gost(query = negcor_raw$gene[order(negcor_raw$cor_r,decreasing = F)],organism = "hsapiens",significant = T,ordered_query = T)

allenrich_results <- list(pos_fdr=raw_r_enrich_pos_fdr,
                          pos_p=raw_r_enrich_pos_raw,
                          neg_fdr=raw_r_enrich_neg_fdr,
                          neg_p=raw_r_enrich_neg_raw)

processed_allenrich_results <- lapply(allenrich_results[c(1,2,4)], function(res) {
go_r_enrich <- subset(res$result, res$result$source %in% c("GO:BP"))
simMatrix <- calculateSimMatrix(go_r_enrich$term_id,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")
scores <- setNames(-log10(go_r_enrich$p_value), go_r_enrich$term_id)

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=0.7,
                                orgdb="org.Hs.eg.db")
hplot <- heatmapPlot(simMatrix,
            reducedTerms,
            annotateParent=TRUE,
            annotationLabel="parentTerm",
            fontsize=6)
tplot <- treemapPlot(reducedTerms)
return(list(simMatrix=simMatrix,
            scores=scores,
            reducedTerms=reducedTerms,
            hplot=hplot,
            tplot=tplot))
})

### plots with highlighted enrichment
# positive correlation, raw p
prp <- gostplot(raw_r_enrich_pos_raw, capped=F, interactive=F) +
  labs(title="positively correlated (raw p < 0.05), n=587")
pospubplot <- publish_gostplot(prp,highlight_terms = c("GO:0042605","GO:0042611","GO:0002483","WP:WP4884"))

# negative correlation, raw p
nrp <- gostplot(raw_r_enrich_neg_raw, capped=F, interactive=F) +
  labs(title="negatively correlated (raw p < 0.05), n=219")
negpubplot <- publish_gostplot(nrp,highlight_terms = c("GO:0003735","GO:0005840","GO:0006614"))

pdf("paper/figures/enrichment_dotplot_for_crosstissue_correlation.pdf",w=8,h=10)
print(plot_grid(pospubplot,negpubplot,nrow=2))
dev.off()

treemapPlot(processed_allenrich_results$pos_p$reducedTerms)
treemapPlot(processed_allenrich_results$neg_p$reducedTerms)

##### perform permutation for cross-tissue correlation
library(doParallel)
detectCores()
registerDoParallel(3)

iterations <- 1000

# permutations without blood correction
stime <- proc.time()
set.seed(1234)
permut.cors <- foreach(icount(iterations), .combine=c) %dopar% {
  mnew <- m2[sample(seq(1,nrow(m2))),]
  dnew <- d2[sample(seq(1,nrow(d2))),]
  allcors <- sapply(seq.int(ncol(mnew)), FUN=function(x) { cor(mnew[,x],dnew[,x],method = "spearman") })
  mean(allcors)
}
proc.time()-stime

# permutations with blood correction
stime <- proc.time()
set.seed(1234)
permut.cors.blood <- foreach(icount(iterations), .combine=c) %dopar% {
  mnew <- m2.blood[sample(seq(1,nrow(m2.blood))),]
  dnew <- d2.blood[sample(seq(1,nrow(d2.blood))),]
  allcors <- sapply(seq.int(ncol(mnew)), FUN=function(x) { cor(mnew[,x],dnew[,x],method = "spearman") })
  mean(allcors)
}
proc.time()-stime

permut.results <- list(permut.cors.blood=permut.cors.blood,
                       permut.cors=permut.cors)
saveRDS(permut.results,"output/cross_tissue_correlation_permutation_results_n1000.rds")

#### plot permutation results
## no blood correction
permutplot1 <- data.frame(rho=permut.cors)
pdf(file="paper/figures/cross_tissue_correlation_noblood_permutations.pdf",w=2.5,h=2.5)
ggplot(data=permutplot1,aes(x=rho))+
  geom_histogram(bins=100,fill="black")+
  geom_vline(xintercept = mean(allinfo.both$cor_r),size=1,col="blue")+
  annotate(geom="text",x=-0.009,y=25,label="p=0",size=5)+
  theme_tufte(base_family = "")
dev.off()

# with blood correction
permutplot2 <- data.frame(rho=permut.cors.blood)
pdf(file="paper/figures/cross_tissue_correlation_withblood_permutations.pdf",w=2.5,h=2.5)
ggplot(data=permutplot2,aes(x=rho))+
  geom_histogram(bins=100,fill="black")+
  geom_vline(xintercept = mean(allinfo.both$cor_r_blood),size=1,col="blue")+
  annotate(geom="text",x=-0.009,y=25,label="p=0",size=5)+
  theme_tufte(base_family = "")
dev.off()

### plot cross-tissue correlation results
# no blood correction
pdf(file="paper/figures/cross_tissue_correlation_noblood.pdf",w=8,h=5)
ggplot(data=allinfo.both, aes(x=cor_r,fill=cor_sig))+
  geom_histogram(bins=500)+
  scale_fill_manual(values = colors()[c(12,30,54)])+
  geom_vline(xintercept = c(0,median(allinfo$cor_r)),col=c("black","yellow"),lty=c(2,4))+
  annotate(x=c(-0.25,-0.19,0.01,0.2,0.34),y=c(13,25,50,25,13),label=as.numeric(table(as.numeric(as.character(allinfo$cor_sig))*sign(allinfo$cor_r))),geom = "label",col=colors()[c(54,30,24,30,54)])+
  labs(y="Number of genes",x="Spearman rho",title = "Cross-tissue correlation without additional clinical corrections")+
  theme_tufte(base_family = "")
dev.off()

# with blood correction
pdf(file="paper/figures/cross_tissue_correlation_withblood.pdf",w=8,h=5)
ggplot(data=allinfo.both, aes(x=cor_r_blood,fill=cor_sig_blood))+
  geom_histogram(bins=500)+
  scale_fill_manual(values = colors()[c(12,30,54)])+
  geom_vline(xintercept = c(0,median(allinfo.both$cor_r_blood)),col=c("black","yellow"),lty=c(2,4))+
  annotate(x=c(-0.25,-0.19,0.01,0.2,0.34),y=c(13,25,50,25,13),label=as.numeric(table(as.numeric(as.character(allinfo.both$cor_sig_blood))*sign(allinfo.both$cor_r_blood))),geom = "label",col=colors()[c(54,30,24,30,54)])+
  theme_tufte(base_family = "")
dev.off()


#################################
#################################
##### GTEX and xQTLserve analysis
#################################
#################################
#### phrased: "to test if effects of monocyte genes on brain phenotypes are due to mutual regulation or separate processes"
# read in xQTL serve data from ROSMAP

teq <- list()
for (chr in seq(1,22)) {
  print(paste0("downloading chromosome ",chr))
  
  xfile <- paste0("http://mostafavilab.stat.ubc.ca/xqtl/xQTL_updated_data/eQTLs/eQTLchr",chr,"_1Mb.csv")
  xfileloc <- paste0("input/xQTLserve/inputfile_chr",chr,".csv")
  download.file(url = xfile,destfile = xfileloc)
  
  print(paste0("reading chromosome ",chr))
  xdat <- read.csv(xfileloc)
  
  print(paste0("finding eQTLs for chromosome ",chr))
  xdatsub <- subset(xdat, ENSG %in% allinfo.both$gene)
  topeqtls <- lapply(unique(xdatsub$ENSG), function(gene) {
    datsub <- xdatsub[which(xdatsub$ENSG==gene),]
    datsub[which(datsub$p==min(datsub$p)),][1,]
  })
  teq[[chr]] <- do.call(rbind,topeqtls)
  
  print(paste0("removing and deleting chromosome ",chr))
  rm(xdat)
  gc()
  unlink(xfileloc)
}


teqall <- do.call(rbind,teq)
#saveRDS(teqall, "input/xQTLserve/xQTL_all_commongene_egenes.rds")

teqall <- readRDS("input/xQTLserve/xQTL_all_commongene_egenes.rds")

xmerged <- merge(allinfo.both,teqall,by.x="gene",by.y="ENSG",all=T)

####### read in GTEX data for frontal cortex
gtex_allfiles <- list.files("input/GTEx",full.names = T)

tomatch <- c("Brain")
gtexfiles <- unique(grep(paste(tomatch,collapse="|"),gtex_allfiles, value=TRUE))

allegeneslist <- lapply(gtexfiles, function(x){
  y <- read.table(x,sep="\t",header=T)
  y$tissue <- gsub("\\..*","",basename(x))
  y
})

allegenes <- do.call(rbind,allegeneslist)
rownames(allegenes) <- NULL
allegenes$gene <- gsub("\\..*","",allegenes$gene_id)
allegenes2 <- allegenes[,c("gene","maf","slope","slope_se","pval_perm","pval_beta","qval","pval_nominal_threshold","log2_aFC","tissue")]


allegenescast <- reshape(allegenes2,
                         v.names = names(allegenes2)[2:9],
                         timevar = "tissue",
                         idvar = "gene",
                         direction="wide")

gtm <- merge(xmerged,allegenescast,by="gene",all.x=T)

#### Plot cross-tissue correlation with eQTL data

eqtl_plot_1 <- ggplot(data=gtm, aes(x=cor_r,y=abs(slope.Brain_Frontal_Cortex_BA9/slope_se.Brain_Frontal_Cortex_BA9),fill=cor_sig))+
  geom_point(aes(col=cor_sig))+
  scale_color_manual(values = colors()[c(12,30,54)])+
  geom_text_repel(data=subset(gtm,cor_fdr < 0.05),aes(label=hugo),size=2.5)+
  geom_vline(xintercept = c(0,median(gtm$cor_r)),col=c("black","orange"),lty=c(2,4))+
  labs(y="abs(t-value)",x="cross-tissue rho")+
  theme_minimal()

eqtl_plot_2 <- ggplot(data=gtm, aes(x=cor_r,y=abs(t),fill=cor_sig))+
  geom_point(aes(col=cor_sig))+
  scale_color_manual(values = colors()[c(12,30,54)])+
  geom_text_repel(data=subset(gtm,cor_fdr < 0.05),aes(label=hugo),size=2.5)+
  geom_vline(xintercept = c(0,median(gtm$cor_r)),col=c("black","orange"),lty=c(2,4))+
  labs(y="abs(t-value)",x="cross-tissue rho")+
  theme_minimal()

plot_grid(eqtl_plot_1,eqtl_plot_2,ncol=2,labels = c("GTEX","xQTLserve"))

### cross-tabs

qtlvars <- grep("slope",names(gtm),value=T)
qtlvars1 <- grep("slope_se",qtlvars,value=T)
qtlvars2 <- grep("slope_se",qtlvars,value=T,invert = T)

allres <- list()
for (tnum in seq(1,length(qtlvars1))) {
  tstat <- abs(gtm[,qtlvars2[tnum]]/gtm[,qtlvars1[tnum]])
  ctest <- cor.test(gtm$cor_r_celltype,tstat)
  allres[[tnum]] <- c(ctest$estimate,ctest$p.value)
}

pdat1 <- as.data.frame(do.call(rbind,allres))
names(pdat1)[2] <- "p"
pdat1$tissue <- gsub(".*\\.","",qtlvars1)

pdat1$tissue <- factor(pdat1$tissue,levels=pdat1$tissue[order(pdat1$p)])
ggplot(data=pdat1, aes(y=cor,x=tissue))+
  geom_bar(stat="identity")+
  coord_flip()+
  theme_minimal()


boxplot(data=ad, abs(slope.Brain_Frontal_Cortex_BA9/slope_se.Brain_Frontal_Cortex_BA9) ~ cor_sig_celltype)
summary(lm(data=ad, abs(slope.Brain_Frontal_Cortex_BA9/slope_se.Brain_Frontal_Cortex_BA9) ~ cor_sig_celltype))


ad$gtex_brainQTL <- ifelse(ad$pval_perm.Brain_Frontal_Cortex_BA9<0.05,1,0)

fisher.test(table(ad$gtex_brainQTL,ad$cor_sig_celltype))

ggplot(data=gtm,aes(y=abs(slope.Brain_Frontal_Cortex_BA9/slope_se.Brain_Frontal_Cortex_BA9),x=cor_sig))+
  geom_violin(draw_quantiles = T)+
  geom_jitter(width=0.1,alpha=0.2,size=0.5)+
  theme_minimal()


summary(lm(data=gtm, abs(slope.Brain_Frontal_Cortex_BA9/slope_se.Brain_Frontal_Cortex_BA9) ~ cor_sig))

