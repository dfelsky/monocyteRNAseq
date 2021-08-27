library(rms)
library(reshape2)
library(ggsci)
library(ggpubr)
library(ggthemes)
library(ggrepel)
library(SuperExactTest)

### read TWAS results
ttall_unlabelled <- readRDS("output/all_TWAS_results.rds")

### name files for plotting
pathnames <- c("Neuritic plaques","Diffuse plaques","Total AB","PHF tau","NFT","Gross cerebral infarcts","Micro cerebral infarcts","Arteriolosclerosis","Cerebral AA","Cerebral atherosclerosis","Lewy body stage","Hippocampal sclerosis","TDP-43","PD Dx","Patho AD","PAM VM Caudate","PAM post. putamen","PAM IT","PAM MF")
cognames <- c("Episodic memory","Perceptual orientation","Perceptual speed","Semantic memory","Working memory","Global","MMSE")

ttall_labelled <- lapply(ttall_unlabelled, function(x) {
  names(x) <- c(pathnames,cognames)
  x
})

### get summary numbers of significant DE genes
list_forsummary05 <- lapply(ttall_labelled, function(tissue) {
  unlist(lapply(tissue, function(pheno) {
    length(pheno$gene[which(pheno$adj.P.Val < 0.05)])
  }))
})
resmat05 <- as.matrix(do.call(rbind,list_forsummary05))
rm05 <- melt(resmat05)
names(rm05) <- c("tissue","pheno","FDR_05")

list_forsummary1 <- lapply(ttall_labelled, function(tissue) {
  unlist(lapply(tissue, function(pheno) {
    length(pheno$gene[which(pheno$adj.P.Val < 0.1)])
  }))
})
resmat1 <- as.matrix(do.call(rbind,list_forsummary1))
rm1 <- melt(resmat1)
names(rm1) <- c("tissue","pheno","FDR_10")

all_rm <- merge(rm05,rm1)
all_rm_m <- melt(all_rm,is.vars=c("tissue","pheno"))
names(all_rm_m)[3] <- "siglevel"

pdf("paper/figures/TWAS_results_summary.pdf",w=8,h=6)
subset(all_rm_m, tissue %in% c("dlpfc_cells","monocyte_blood") & siglevel=="FDR_05") %>%
ggplot(aes(y=value,x=pheno,fill=siglevel))+
  geom_bar(stat = "identity",position = "identity",show.legend = F)+
  facet_wrap(~tissue,nrow=1,scales="free")+
  geom_text(aes(label=value),hjust=-0.1,size=2.5)+
  coord_flip()+
  scale_fill_aaas()+
  labs(x="Phenotype",y="Number of genes (FDR<0.05)")+
  theme_tufte(base_family = "sans")
dev.off()


### overlapping genes for PAM phenotypes
pamsig <- list(MF=ttall_unlabelled$monocyte_blood$mf3123$hugo[which(ttall_unlabelled$monocyte_blood$mf3123$adj.P.Val<0.05)],
               IT=ttall_unlabelled$monocyte_blood$it3123$hugo[which(ttall_unlabelled$monocyte_blood$it3123$adj.P.Val<0.05)],
               "Post. putamen"=ttall_unlabelled$monocyte_blood$pput3123$hugo[which(ttall_unlabelled$monocyte_blood$pput3123$adj.P.Val<0.05)],
               "VM caudate"=ttall_unlabelled$monocyte_blood$vm3123$hugo[which(ttall_unlabelled$monocyte_blood$vm3123$adj.P.Val<0.05)])

pamstest <- supertest(pamsig,n=nrow(ttall_unlabelled$monocyte_blood$mf3123))
pdf(file="paper/figures/fig2_mset_pam.pdf",w=5,h=4)
plot.msets(pamstest,Layout = "landscape",
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




plot.msets(pamstest,Layout = "landscape",
           keep.empty.intersections = T,
           degree=4,
           show.elements = T)



### overlapping genes for amyloid phenotypes
amyloidsig <- list("Total AB"=ttall_unlabelled$monocyte_blood$amyloid_sqrt$gene[which(ttall_unlabelled$monocyte_blood$amyloid_sqrt$adj.P.Val<0.05)],
               "Neuritic plaques"=ttall_unlabelled$monocyte_blood$plaq_n_sqrt$gene[which(ttall_unlabelled$monocyte_blood$plaq_n_sqrt$adj.P.Val<0.05)],
               "Diffuse plaques"=ttall_unlabelled$monocyte_blood$plaq_d_sqrt$gene[which(ttall_unlabelled$monocyte_blood$plaq_d_sqrt$adj.P.Val<0.05)],
               "Cerebral AA"=ttall_unlabelled$monocyte_blood$caa_4gp$gene[which(ttall_unlabelled$monocyte_blood$caa_4gp$adj.P.Val<0.05)],
               "Patho AD"=ttall_unlabelled$monocyte_blood$pathoAD$gene[which(ttall_unlabelled$monocyte_blood$pathoAD$adj.P.Val<0.05)])

amyloidstest <- supertest(amyloidsig,n=nrow(ttall_unlabelled$monocyte_blood$mf3123))
pdf(file="paper/figures/fig2_mset_amyloid.pdf",w=5,h=4)
plot.msets(amyloidstest,Layout = "landscape",
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

### overlapping genes for cognitive phenotypes
cogsig <- list("MMSE"=ttall_unlabelled$monocyte_blood$cts_mmse30_at_draw$gene[which(ttall_unlabelled$monocyte_blood$cts_mmse30_at_draw$adj.P.Val<0.05)],
                   "Global"=ttall_unlabelled$monocyte_blood$cogn_global_at_draw$gene[which(ttall_unlabelled$monocyte_blood$cogn_global_at_draw$adj.P.Val<0.05)],
                   "Working memory"=ttall_unlabelled$monocyte_blood$cogn_wo_at_draw$gene[which(ttall_unlabelled$monocyte_blood$cogn_wo_at_draw$adj.P.Val<0.05)],
                   "Semantic memory"=ttall_unlabelled$monocyte_blood$cogn_se_at_draw$gene[which(ttall_unlabelled$monocyte_blood$cogn_se_at_draw$adj.P.Val<0.05)],
                   "Perceptual speed"=ttall_unlabelled$monocyte_blood$cogn_ps_at_draw$gene[which(ttall_unlabelled$monocyte_blood$cogn_ps_at_draw$adj.P.Val<0.05)],
                   "Perceptual orientation"=ttall_unlabelled$monocyte_blood$cogn_po_at_draw$gene[which(ttall_unlabelled$monocyte_blood$cogn_po_at_draw$adj.P.Val<0.05)],
                   "Episodic memory"=ttall_unlabelled$monocyte_blood$cogn_ep_at_draw$gene[which(ttall_unlabelled$monocyte_blood$cogn_ep_at_draw$adj.P.Val<0.05)])

cogstest <- supertest(cogsig,n=nrow(ttall_unlabelled$monocyte_blood$mf3123))
pdf(file="paper/figures/fig2_mset_cognition.pdf",w=6,h=4)
plot.msets(cogstest,Layout = "landscape",
           keep.empty.intersections = F,
           degree=c(2,3,4,5,6,7),
           show.elements = F,
           margin = c(1,10,10,5),
           color.on = "black",
           sort.by = "size",
           overlap.size.cex = 0.5,
           legend.text.cex = 0.5,
           color.scale.cex = 0.5,
           cex=0.5,yfrac = 0.5)
dev.off()



##############
#### plot volcano plots for PAM phenotypes and diffuse plaques and atherosclerosis

p1 <- ggplot(ttall_unlabelled$monocyte_blood$mf3123,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title = "A. PAM (MF)")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=ttall_unlabelled$monocyte_blood$mf3123[1:5,],aes(label=hugo))+
  theme_minimal()
p2 <- ggplot(ttm$it3123,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title = "B. PAM (IT)")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=ttm$it3123[1:5,],aes(label=hugo))+
  theme_minimal()
p3 <- ggplot(ttm$vm3123,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title="C. PAM (VM)")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=ttm$vm3123[1:5,],aes(label=hugo))+
  theme_minimal()
p4 <- ggplot(ttm$pput3123,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title="D. PAM (PPUT)")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=ttm$pput3123[1:5,],aes(label=hugo))+
  theme_minimal()
p5 <- ggplot(ttm$plaq_d_sqrt,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title="E. Diffuse plaques")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=ttm$plaq_d_sqrt[1:5,],aes(label=hugo))+
  theme_minimal()
p6 <- ggplot(ttm$cvda_4gp2,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title="F. Cereb. atheroscl.")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=ttm$cvda_4gp2[1:5,],aes(label=hugo))+
  theme_minimal()

plot_grid(p1,p2,p3,p4,p5,p6,ncol=3)