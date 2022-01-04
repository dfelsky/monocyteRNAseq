#### Microglial and diffuse plaque monocyte and DLPFC proxy
library(SuperExactTest)
library(factoextra)
library(reshape2)
library(ggthemes)
library(ggsci)
library(ggpubr)
library(cowplot)
library(rms)
library(ggrepel)
load("input/all_genes_ensembl.RData")
ROSmaster <- readRDS("input/ROSmaster_TWAS_input.rds")

mono_expr <- readRDS("output/WGCNA/mono_blood/input_datExpr.rds")
ttall_unlabelled <- readRDS("output/all_TWAS_results.rds")
ttm <- ttall_unlabelled$monocyte_blood
monov <- readRDS("input/monocytes_v_blood.rds")

####################################################
### overlap of monocyte genes for PAM phenotypes ###
####################################################
mono_mf3123_genes <- ttm$mf3123$gene[which(ttm$mf3123$adj.P.Val<0.05)]
mono_it3123_genes <- ttm$it3123$gene[which(ttm$it3123$adj.P.Val<0.05)]
mono_vm3123_genes <- ttm$vm3123$gene[which(ttm$vm3123$adj.P.Val<0.05)]
mono_pput3123_genes <- ttm$pput3123$gene[which(ttm$pput3123$adj.P.Val<0.05)]

mono_pam_overlap <- supertest(list(IT=mono_mf3123_genes,
               MF=mono_it3123_genes,
               VM=mono_vm3123_genes,
               PPUT=mono_pput3123_genes),n = 9129,degree=c(2,3,4))

plot.msets(mono_pam_overlap,
           Layout="landscape",
           keep.empty.intersections=F,
           sort.by="size",
           show.elements = F,
           margin = c(0.5,5,8,2))

mono_mfpca1 <- prcomp(scale(mono_expr[,mono_mf3123_genes]))
pambiplot1 <- fviz_pca_var(mono_mfpca1,label = "none",title = "MF PAM")
mono_mfpca2 <- prcomp(scale(mono_expr[,mono_it3123_genes]))
pambiplot2 <- fviz_pca_var(mono_mfpca2,label = "none",title = "IT PAM")
mono_mfpca3 <- prcomp(scale(mono_expr[,mono_vm3123_genes]))
pambiplot3 <- fviz_pca_var(mono_mfpca3,label = "none",title = "VM PAM")
mono_mfpca4 <- prcomp(scale(mono_expr[,mono_pput3123_genes]))
pambiplot4 <- fviz_pca_var(mono_mfpca4,label = "none",title = "PPUT PAM")

pdf("paper/figures/PCA_biplot_PAM_siggenes.pdf",w=8,h=8)
print(plot_grid(pambiplot1,pambiplot2,pambiplot3,pambiplot4))
dev.off()

fviz_screeplot(mono_mfpca)

groups <- as.factor(pheno_mono$pathoAD)
fviz_pca_var(mono_mfpca)


########## last section modelling
# run PCA for all phenotypes, note variance explained, and package into data frame
mono_exp2 <- as.data.frame(mono_expr)
mono_exp2$projid <- rownames(mono_expr)

var.explained <- list()
allpca <- list()
ngenes <- NULL
for (pheno in names(ttm)) {
  siggenes <- ttm[[pheno]]$gene[which(ttm[[pheno]]$adj.P.Val < 0.05)]
  if (length(siggenes)<2) { next }
  ngenes[[pheno]] <- length(siggenes)
  pca1 <- prcomp(scale(mono_expr[,siggenes]))
  var.explained[[pheno]] <- summary(pca1)$importance[,1]
  num.feats <- min(c(10,length(siggenes)))
  newpcs <- pca1$x[,1:num.feats]
  colnames(newpcs) <- paste0("PC",seq(1,num.feats),"_",pheno)
  mono_exp2 <- cbind(mono_exp2,newpcs)
  allpca[[pheno]] <- pca1
}

ve <- as.data.frame(do.call(rbind,var.explained))
ve$n <- as.numeric(ngenes)
nameref <- readRDS("output/variable_name_reference.rds")
ve$Phenotype <- factor(rownames(ve), levels=rownames(ve),labels=nameref$varnames[match(rownames(ve),nameref$mono.variable)])
ve$Phenotype <- paste0(ve$Phenotype," (",ve$n,")")
ve$Phenotype <- factor(ve$Phenotype,levels=ve$Phenotype[order(ve$`Proportion of Variance`)])

pdf("paper/figures/Fig5_var_explained_for_each_phenotype.pdf",w=6,h=5)
ggplot(data=ve, aes(y=`Proportion of Variance`,x=Phenotype))+
  geom_bar(stat="identity")+
  theme_minimal()+
  labs(y="% Variance explained by PC1")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1))
dev.off()
  
mm <- merge(mono_exp2,ROSmaster,by="projid",all=T)
mm2 <- merge(mm,monov$targets,by="projid",all=T)
mm2$APOE <- factor(mm2$apoe4d,levels=c(0,1),labels=c("non-e4","e4"))
mm2$Dx <- mm2$cogdx
mm2$Dx[which(mm2$Dx %in% c(3,5,6))] <- NA
mm2$Dx[which(mm2$Dx==4)] <- 3

pvals <- NULL
phens <- NULL
interfacts <- NULL
PCs <- NULL
index <- 1
for (phen in rownames(ve)) {
  for (interfact in c("APOE","Dx","msex.x")) {
    for (PC in seq(1,5)) {
      if (length(grep(paste0("PC",PC,"_",phen),names(mm2)))==0) { next }
      yphen <- phen
      if (yphen %in% c("ci_num2_mct","parkdx","pathoAD")) {
        form <- formula(paste0(yphen," ~ ",paste0("PC",PC,"_",phen),"*",interfact))
        mod <- glm(data=mm2,form,family=binomial())
        summod <- summary(mod)
        pvals[index] <- summod$coefficients[4,4]
      } else {
        form <- formula(paste0(yphen," ~ ",paste0("PC",PC,"_",phen),"*",interfact))
        mod <- lm(data=mm2,form)
        pvals[index] <- anova(mod)[3,5]
      }
      phens[index] <- phen
      PCs[index] <- PC
      interfacts[index] <- interfact
      index <- index + 1
    }
  }
}

resframe <- data.frame(y=phens,int=interfacts,PC=PCs,p=pvals)
resframe <- subset(resframe, PC==1)
resframe$fdr <- p.adjust(resframe$p)
resframe$Interaction <- factor(resframe$int,levels=c("APOE","Dx","msex.x"),labels=c("APOE genotype","AD Dx","Sex"))
nameref <- readRDS("output/variable_name_reference.rds")
resframe$Phenotype <- factor(resframe$y, levels=unique(resframe$y),labels=nameref$varnames[match(unique(resframe$y),nameref$mono.variable)])
resframe$Significance <- ifelse(resframe$fdr < 0.05, "FDR < 0.05","NS")

pdf("paper/figures/Fig5_all_interactions_dotplot.pdf",w=8,h=5)
ggplot(data=resframe, aes(y=-log10(p),x=Phenotype,col=Interaction))+
  geom_point(aes(shape=Significance),size=3)+
  scale_shape_manual(values=c(19,15))+
  scale_color_tableau()+
  geom_hline(yintercept = -log10(0.05),col="red",lty=2)+
  theme_minimal()+
  labs(y="Interaction term -log10(p-value)")+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, size = 10, hjust = 1))
dev.off()

mm3 <- subset(mm2, Dx %in% c(1,2,3))
mm3$Dx <- factor(mm3$Dx,levels=c(1,2,3),labels=c("CN","MCI","AD"))

# individual plots for sig effects
pdf("paper/figures/Fig5_mmse_cogdx.pdf",w=8,h=5)
ggplot(data=mm3, aes(y=cts_mmse30_at_draw ,x=PC1_cts_mmse30_at_draw ,col=Dx,group=Dx))+
  geom_point()+
  geom_smooth(method="lm")+
  #facet_wrap(~Dx,scales="fixed")+
  labs(y="MMSE at draw",x="MMSE-related monocyte gene expression (PC1)")+
  scale_color_tableau()+
  theme_minimal()
dev.off()

ols(data=mm2, cts_mmse30_at_draw ~ PC1_cts_mmse30_at_draw*Dx)

# individual plots for sig effects
pdf("paper/figures/Fig5_cogn_global_apoe.pdf",w=6,h=5)
subset(mm2, is.na(apoe4d)==F) %>%
ggplot(aes(y=cogn_global_at_draw ,x=PC1_cogn_global_at_draw ,col=APOE))+
  geom_point()+
  geom_smooth(method="lm")+
  labs(y="Global cognition at draw",x="Cognition-related monocyte gene expression (PC1)")+
  scale_color_tableau()+
  theme_minimal()
dev.off()

