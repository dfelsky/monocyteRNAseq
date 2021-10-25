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
#colnames(mono_expr) <- all_genes$external_gene_name[match(colnames(mono_expr),all_genes$ensembl_gene_id)]
#mono_expr <- mono_expr[,which(is.na(colnames(mono_expr))==F)]
dlpfc_expr <- readRDS("output/WGCNA/dlpfc/input_adjustedforcells_datExpr.rds")

ttall_unlabelled <- readRDS("output/all_TWAS_results.rds")

ttm <- ttall_unlabelled$monocyte_blood
ttd <- ttall_unlabelled$dlpfc_cells

pheno_mono <- ROSmaster[match(rownames(mono_expr),ROSmaster$projid),]
pheno_dlpfc <- ROSmaster[match(rownames(dlpfc_expr),ROSmaster$projid),]

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

mm3 <- subset(mm2, cogdx %in% c(1,2,3))
mm3$Dx <- factor(mm3$cogdx,levels=c(1,2,3),labels=c("CN","MCI","AD"))

# individual plots for sig effects
pdf("paper/figures/Fig5_mmse_cogdx.pdf",w=8,h=5)
ggplot(data=mm3, aes(y=cts_mmse30_at_draw ,x=PC1_cts_mmse30_at_draw ,col=Dx))+
  geom_point()+
  geom_smooth(method="lm")+
  facet_wrap(~Dx,scales="free_x")+
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


#################################################
### overlap of dlpfc genes for PAM phenotypes ###
#################################################
dlpfc_mf3123_genes <- ttd$mf3123$hugo[which(ttd$mf3123$adj.P.Val<0.05)]
dlpfc_it3123_genes <- ttd$it3123$hugo[which(ttd$it3123$adj.P.Val<0.05)]
dlpfc_vm3123_genes <- ttd$vm3123$hugo[which(ttd$vm3123$adj.P.Val<0.05)]
dlpfc_pput3123_genes <- ttd$pput3123$hugo[which(ttd$pput3123$adj.P.Val<0.05)]

dlpfc_pam_overlap <- supertest(list(IT=dlpfc_mf3123_genes,
                                   MF=dlpfc_it3123_genes,
                                   VM=dlpfc_vm3123_genes,
                                   PPUT=dlpfc_pput3123_genes),n = 17422,degree=c(2,3,4))

plot.msets(dlpfc_pam_overlap,
           Layout="landscape",
           keep.empty.intersections=F,
           sort.by="size",
           show.elements = F,
           margin = c(0.5,5,8,2))

dlpfc_mfpca <- prcomp(scale(dlpfc_expr[,dlpfc_mf3123_genes]))
fviz_pca_biplot(dlpfc_mfpca)
fviz_screeplot(dlpfc_mfpca)

groups <- as.factor(pheno_dlpfc$pathoAD)
fviz_pca_var(dlpfc_mfpca)

########################################
### correlate DLPFC and mono PCs #######
########################################
dlpfc_pcs <- as.data.frame(dlpfc_mfpca$x)
names(dlpfc_pcs) <- paste0("dlpfc_",names(dlpfc_pcs))
dlpfc_pcs$projid <- rownames(dlpfc_pcs)
mono_pcs <- as.data.frame(mono_mfpca$x)
names(mono_pcs) <- paste0("mono_",names(mono_pcs))
mono_pcs$projid <- rownames(mono_pcs)

both_pcs <- merge(mono_pcs,dlpfc_pcs,by="projid")
datmat <- as.matrix(both_pcs[,-1])
cormat <- rcorr(datmat)

cormeltr <- melt(cormat$r)
cormeltp <- melt(cormat$P)
bothmelt <- merge(cormeltr,cormeltp,by=c("Var1","Var2"))
names(bothmelt) <- c("x","y","cor","p")

bms <- bothmelt[which(bothmelt$x!=bothmelt$y),]
bms <- bms[which(grep("dlpfc",bms$x) %nin% grep("dlpfc",bms$y)),]

plotlabs <- ifelse(bms$p<0.01,
                   round(-log10(bms$p)),
                   NA)

ggplot(data=bms,aes(y=y,x=x,fill=cor))+
  geom_tile(color = "white")+
  scale_fill_gradient2_tableau()+
  geom_text(label=plotlabs)+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))+
  coord_fixed()


##### YIELDS PC30 (mono) --> PC2/PC3 (dlpfc), PC9 (mono) --> PC5 (dlpfc)
# block 1
fviz_contrib(mono_mfpca,"var",axes = 30)
fviz_contrib(dlpfc_mfpca,"var",axes = 2)
fviz_contrib(dlpfc_mfpca,"var",axes = 3)

# block 2
fviz_contrib(mono_mfpca,"var",axes = 9)
fviz_contrib(dlpfc_mfpca,"var",axes = 5)

############ correlate individual genes ###############
monogenes <- as.data.frame(mono_expr[,mono_mf3123_genes])
names(monogenes) <- paste0("mono_",names(monogenes))
monogenes$projid <- rownames(monogenes)
dlpfcgenes <- as.data.frame(dlpfc_expr[,dlpfc_mf3123_genes])
names(dlpfcgenes) <- paste0("dlpfc_",names(dlpfcgenes))
dlpfcgenes$projid <- rownames(dlpfcgenes)

bothgenes <- merge(monogenes,dlpfcgenes,by="projid")
datmat <- as.matrix(bothgenes[,-1])
cormat <- rcorr(datmat)

cormeltr <- melt(cormat$r)
cormeltp <- melt(cormat$P)
bothmelt <- merge(cormeltr,cormeltp,by=c("Var1","Var2"))
names(bothmelt) <- c("x","y","cor","p")

bms <- bothmelt[which(bothmelt$x!=bothmelt$y),]
bms <- bms[which(grep("dlpfc_",bms$x) %nin% grep("dlpfc_",bms$y)),]

plotlabs <- ifelse(bms$p<0.01,
                   round(-log10(bms$p)),
                   NA)

ggplot(data=bms,aes(y=y,x=x,fill=cor))+
  geom_tile(color = "white")+
  scale_fill_gradient2_tableau()+
  geom_text(label=plotlabs)+
  theme_minimal()+ 
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))


##### YIELDS FP236383.3-, FABP5-, CD82-, HMBS+ (mono) --> PELO (dlpfc), MYL5+ (mono) --> KIAA0408 (dlpfc)
### checked on genemania.org (coexpression across multiple human studies)

mania <- read.table("TWAS/Jan2021/genemania_FABP5_CD82_HMBS_PELO_MYL5_KIAA0408.txt",sep="\t",header=T)

mgenes <- unique(c(mania$Gene.1,mania$Gene.2))
mgenes <- mgenes[which(mgenes %nin% c("FABP5","CD82","HMBS","PELO","MYL5","KIAA0408"))]
match(mgenes,ttm$it3123$hugo)
match(mgenes,ttd$it3123$hugo)
#### NOTHING #####
### also checked https://www.inetbio.org/humannet/ (humanNet v2)
# some candidate genes but no significant enrichment for candidates in provided guides (provided all 52 mono genes)


######## overlap with ellis paper
ellis <- read.csv("patrick_m5_paper_supptable3.csv")
ellis$mf3123_t <- ttd$mf3123$t[match(ellis$Ensembl,ttd$mf3123$gene)]
ellis$mf3123_p <- ttd$mf3123$P.Value[match(ellis$Ensembl,ttd$mf3123$gene)]
ellis$it3123_p <- ttd$it3123$P.Value[match(ellis$Ensembl,ttd$it3123$gene)]
ellis$it3123_t <- ttd$it3123$t[match(ellis$Ensembl,ttd$it3123$gene)]

# all modules
ellissiggenes <- ellis$Gene[which(ellis$P.value.Activated.Microglia<0.05)]
newsiggenesmf <- ellis$Gene[which(ellis$mf3123_p<0.05)]
newsiggenesit <- ellis$Gene[which(ellis$it3123_p<0.05)]
ellistest <- supertest(list(ellis=ellissiggenes,
               newmf=newsiggenesmf,
               newit=newsiggenesit),n = nrow(ttd$mf3123))
plot.msets(ellistest,Layout = "landscape")

# individual modules
ellissiggenes5 <- ellis$Gene[which(ellis$P.value.Activated.Microglia<0.05 & ellis$Module=="5")]
ellissiggenes113 <- ellis$Gene[which(ellis$P.value.Activated.Microglia<0.05 & ellis$Module=="113")]
ellissiggenes114 <- ellis$Gene[which(ellis$P.value.Activated.Microglia<0.05 & ellis$Module=="114")]
ellissiggenes115 <- ellis$Gene[which(ellis$P.value.Activated.Microglia<0.05 & ellis$Module=="115")]
ellissiggenes116 <- ellis$Gene[which(ellis$P.value.Activated.Microglia<0.05 & ellis$Module=="116")]

ellistest5 <- supertest(list(ellism5=ellissiggenes5,
                            newmf=newsiggenesmf,
                            newit=newsiggenesit),n = nrow(ttd$mf3123))
ellistest113 <- supertest(list(ellism113=ellissiggenes113,
                             newmf=newsiggenesmf,
                             newit=newsiggenesit),n = nrow(ttd$mf3123))
ellistest114 <- supertest(list(ellism114=ellissiggenes114,
                             newmf=newsiggenesmf,
                             newit=newsiggenesit),n = nrow(ttd$mf3123))
ellistest115 <- supertest(list(ellism115=ellissiggenes115,
                             newmf=newsiggenesmf,
                             newit=newsiggenesit),n = nrow(ttd$mf3123))
ellistest116 <- supertest(list(ellism116=ellissiggenes116,
                             newmf=newsiggenesmf,
                             newit=newsiggenesit),n = nrow(ttd$mf3123))

plot.msets(ellistest5,Layout = "landscape")
plot.msets(ellistest113,Layout = "landscape")
plot.msets(ellistest114,Layout = "landscape")
plot.msets(ellistest115,Layout = "landscape")
plot.msets(ellistest116,Layout = "landscape")


################################################################
#### scatterplot for DLPFC and MONO PCs and diffuse plaques ####
################################################################
mono_dp_genes <- ttm$plaq_d_sqrt$hugo[which(ttm$plaq_d_sqrt$adj.P.Val<0.05)]
dlpfc_dp_genes <- ttd$plaq_d_sqrt$hugo[which(ttd$plaq_d_sqrt$adj.P.Val<0.05)]

DP_overlap <- supertest(list(MONO_dp=mono_dp_genes,
                             DLPFC_dp=dlpfc_dp_genes),n = 8762)

plot.msets(DP_overlap,
           Layout="circular",
           keep.empty.intersections=F,
           sort.by="size",
           show.elements = T,
           margin = c(0.5,5,8,2))

mono_dp_pca <- prcomp(scale(mono_expr[,mono_dp_genes]))
dlpfc_dp_pca <- prcomp(scale(dlpfc_expr[,dlpfc_dp_genes]))

fviz_pca_biplot(mono_dp_pca,label = c("quali","var"),axes = c(3,5))
fviz_screeplot(mono_dp_pca)

dlpfc_pcs <- as.data.frame(dlpfc_dp_pca$x)
names(dlpfc_pcs) <- paste0("dlpfc_",names(dlpfc_pcs))
dlpfc_pcs$projid <- rownames(dlpfc_pcs)
mono_pcs <- as.data.frame(mono_dp_pca$x)
names(mono_pcs) <- paste0("mono_",names(mono_pcs))
mono_pcs$projid <- rownames(mono_pcs)

mono_dp_merge <- merge(mono_pcs,dlpfc_pcs,by="projid")
mdpm <- merge(mono_dp_merge,ROSmaster,by="projid")
mdpm <- subset(mdpm, cogdx %in% c(1,2,4))

p1 <- ggplot(data=mdpm,aes(y=mono_PC1,x=plaq_d_sqrt))+
  geom_point()+
  geom_smooth(method="lm",se=F)+
  facet_wrap(~cogdx,scales="free",ncol=6)+
  labs(y="mono PC1 (DP genes)",x="Diffuse plaque count (sqrt)",title="PC1 stratified by clinical diagnosis")+
  theme_minimal()
p2 <- ggplot(data=mdpm,aes(y=mono_PC2,x=plaq_d_sqrt))+
  geom_point()+
  geom_smooth(method="lm",se=F)+
  facet_wrap(~cogdx,scales="free",ncol=6)+
  labs(y="mono PC2 (DP genes)",x="Diffuse plaque count (sqrt)",title="PC2 stratified by clinical diagnosis")+
  theme_minimal()
p3 <- ggplot(data=mdpm,aes(y=mono_PC3,x=plaq_d_sqrt))+
  geom_point()+
  geom_smooth(method="lm",se=F)+
  facet_wrap(~cogdx,scales="free",ncol=6)+
  labs(y="mono PC3 (DP genes)",x="Diffuse plaque count (sqrt)",title="PC3 stratified by clinical diagnosis")+
  theme_minimal()
p4 <- ggplot(data=mdpm,aes(y=mono_PC4,x=plaq_d_sqrt))+
  geom_point()+
  geom_smooth(method="lm",se=F)+
  facet_wrap(~cogdx,scales="free",ncol=6)+
  labs(y="mono PC4 (DP genes)",x="Diffuse plaque count (sqrt)",title="PC4 stratified by clinical diagnosis")+
  theme_minimal()
p5 <- ggplot(data=mdpm,aes(y=mono_PC5,x=plaq_d_sqrt))+
  geom_point()+
  geom_smooth(method="lm",se=F)+
  facet_wrap(~cogdx,scales="free",ncol=6)+
  labs(y="mono PC5 (DP genes)",x="Diffuse plaque count (sqrt)",title="PC5 stratified by clinical diagnosis")+
  theme_minimal()

plot_grid(p1,p2,p3,p4,p5,nrow=5)

summary(lm(data=mdpm, plaq_d_sqrt ~ mono_PC1*as.factor(cogdx)))

ggplot(data=mdpm,aes(y=plaq_d_sqrt,x=dlpfc_PC1))+
  geom_point()+
  geom_smooth(method="lm",se=F)+
  facet_wrap(~cogdx,scales="free",ncol=6)+
  theme_minimal()

fviz_pca_biplot(dlpfc_dp_pca)
fviz_screeplot(dlpfc_dp_pca)

##############################################################################################
### Build predictive models one gene and one PCA at a time to compare growth (adjusted R2) ###
##############################################################################################
####################### Monocytes #########################
genemods <- NULL
pcmods <- NULL
phenvec <- NULL
itervec <- NULL
validgenemods <- NULL
validpcmods <- NULL

phenolist <- c("mf3123","it3123","vm3123","pput3123","plaq_d_sqrt","arteriol_scler","caa_4gp")

tissue <- "mono"
feature.num.limit <- 20

full.index <- 1
for (yphen in phenolist) {
  if (tissue=="mono") {
    expdat <- mono_expr
    tt <- ttm
    phenodat <- pheno_mono 
    } else { 
      expdat <- dlpfc_expr
      tt <- ttd
      phenodat <- pheno_dlpfc
      }
  
  tt[[yphen]] <- tt[[yphen]][which(is.na(tt[[yphen]]$gene)==F),]
  yphen_sig_genes <- tt[[yphen]]$gene[which(tt[[yphen]]$adj.P.Val<0.05)]
  y <- phenodat[,yphen]
  x <- scale(expdat[,yphen_sig_genes])
  yphen_pca <- prcomp(x)
  maxnumfeat <- ifelse(length(yphen_sig_genes) <= feature.num.limit,
                       length(yphen_sig_genes),
                       feature.num.limit)
  
  iter.index <- 1
  for (featnum in seq(1,maxnumfeat)) {
    GM <- ols(y ~ x[,1:featnum],y=T,x=T)
    genemods[full.index] <- ifelse(summary.lm(GM)[["adj.r.squared"]] < 0,0,summary.lm(GM)[["adj.r.squared"]])
    
    set.seed(999)
    valGM <- validate(GM,B=100,method=".632")
    validgenemods[full.index] <- ifelse(valGM[1,5]<0,0,valGM[1,5])
    
    PM <- ols(y ~ yphen_pca$x[,1:featnum],x=T,y=T)
    pcmods[full.index] <- ifelse(summary.lm(PM)[["adj.r.squared"]] < 0,0,summary.lm(PM)[["adj.r.squared"]])
    
    set.seed(999)
    valPM <- validate(PM,B=100,method=".632")
    validpcmods[full.index] <- ifelse(valPM[1,5]<0,0,valPM[1,5])
    
    phenvec[full.index] <- yphen
    itervec[full.index] <- iter.index
    
    iter.index  <- iter.index + 1
    full.index <- full.index + 1
  }
  }

result <- data.frame(genes=genemods,
                     genes.boot=validgenemods,
                     PC=pcmods,
                     PC.boot=validpcmods,
                     iter=itervec,
                     phenotype=phenvec)


resmelt <- melt(result,id.vars = c("iter","phenotype"))
resmelt$phenotype <- factor(resmelt$phenotype,levels=phenolist)

pr2 <- ggplot(subset(resmelt,variable %in% c("genes","PC")),aes(y=value,x=iter,col=variable))+
  geom_line(size=1.5)+
  facet_wrap(~phenotype,nrow=1,scales = "free_x")+
  scale_color_aaas()+
  labs(title="Adjusted R2",y="adjusted R2",x="Number of genes/PCs in model")+
  theme_minimal()

pvalid <- ggplot(subset(resmelt,variable %nin% c("genes","PC")),aes(y=value,x=iter,col=variable))+
  geom_line(size=1.5)+
  facet_wrap(~phenotype,nrow=1,scales = "free_x")+
  scale_color_aaas()+
  labs(title="Bootstrapped R2 (b=100)",y="Boot R2",x="Number of genes/PCs in model")+
  theme_minimal()

pdf("paper/figures/Fig5_PC_and_gene_modelling_bootstrap_predictor.pdf",w=10,h=6)
print(plot_grid(pr2,pvalid,nrow=2))
dev.off()



######################
### run other modelling with covariates in model
ROSmaster <- readRDS("input/ROSmaster_TWAS_input.rds")
monov <- readRDS("input/monocytes_v_blood.rds")
e1 <- as.data.frame(t(monov$E))
monodf <- cbind(monov$targets,e1)
monodf$est.lib.size <- as.numeric(scale(monodf$ESTIMATED_LIBRARY_SIZE))

tt <- ttall_unlabelled$monocyte_blood

genemods <- NULL
pcmods <- NULL
phenvec <- NULL
itervec <- NULL
validgenemods <- NULL
validpcmods <- NULL

phenolist <- c("mf3123","it3123","vm3123","pput3123","plaq_d_sqrt","arteriol_scler","caa_4gp")
techvars <- c("batch","PCT_USABLE_BASES","PERCENT_DUPLICATION","MEDIAN_3PRIME_BIAS","study","PCT_PF_READS_ALIGNED","est.lib.size")
covs <- c(techvars,"msex","age_death","pmi","age_draw","hemoglbn_at_draw","mchc_at_draw","mcv_at_draw","platelet_at_draw","wbc_at_draw","fasting.f","hemotologic_rx_at_draw")

tissue <- "mono"
feature.num.limit <- 10

full.index <- 1
for (yphen in phenolist) {
  
  if(yphen %in% c("mf3123","it3123","vm3123","pput3123")) {
    covars <- covs[which(covs %nin% "study")]
  } else {
    covars <- covs
  }
  
  monodfiter <- monodf
  monodfiter[,yphen] <- ROSmaster[,yphen][match(monodfiter$projid,ROSmaster$projid)]
  monodfiter <- cbind(monodfiter,ROSmaster[,c("age_death","pmi")][match(monodfiter$projid,ROSmaster$projid),])

  yphen_sig_genes <- tt[[yphen]]$gene[which(tt[[yphen]]$adj.P.Val<0.05)]
  x <- scale(monodfiter[,yphen_sig_genes])
  yphen_pca <- prcomp(x)
  monodfiter <- cbind(monodfiter,yphen_pca$x)
  maxnumfeat <- ifelse(length(yphen_sig_genes) <= feature.num.limit,
                       length(yphen_sig_genes),
                       feature.num.limit)
  
  iter.index <- 1
  for (featnum in seq(1,maxnumfeat)) {
    
    form1 <- formula(paste(yphen," ~ ",paste(yphen_sig_genes[1:featnum],collapse = "+"),"+",paste(covars,collapse="+")))
    GM <- ols(data=monodfiter,form1,y=T,x=T)
    genemods[full.index] <- ifelse(summary.lm(GM)[["adj.r.squared"]] < 0,0,summary.lm(GM)[["adj.r.squared"]])
    
    set.seed(999)
    valGM <- validate(GM,B=100,method=".632")
    validgenemods[full.index] <- ifelse(valGM[1,5]<0,0,valGM[1,5])
    
    
    form2 <- formula(paste(yphen," ~ ",paste0("PC",seq(1:featnum),collapse = "+"),"+",paste(covars,collapse="+")))
    PM <- ols(data=monodfiter,form2,x=T,y=T)
    pcmods[full.index] <- ifelse(summary.lm(PM)[["adj.r.squared"]] < 0,0,summary.lm(PM)[["adj.r.squared"]])
    
    set.seed(999)
    valPM <- validate(PM,B=100,method=".632")
    validpcmods[full.index] <- ifelse(valPM[1,5]<0,0,valPM[1,5])
    
    phenvec[full.index] <- yphen
    itervec[full.index] <- iter.index
    
    iter.index  <- iter.index + 1
    full.index <- full.index + 1
  }
}

result <- data.frame(genes=genemods,
                     genes.boot=validgenemods,
                     PC=pcmods,
                     PC.boot=validpcmods,
                     iter=itervec,
                     phenotype=phenvec)

resmelt <- melt(result,id.vars = c("iter","phenotype"))
resmelt$phenotype <- factor(resmelt$phenotype,levels=phenolist)

pr2 <- ggplot(subset(resmelt,variable %in% c("genes","PC")),aes(y=value,x=iter,col=variable))+
  geom_line(size=1.5)+
  facet_wrap(~phenotype,nrow=1,scales = "free_x")+
  scale_color_aaas()+
  labs(title="Adjusted R2",y="adjusted R2",x="Number of genes/PCs in model")+
  theme_minimal()

pvalid <- ggplot(subset(resmelt,variable %nin% c("genes","PC")),aes(y=value,x=iter,col=variable))+
  geom_line(size=1.5)+
  facet_wrap(~phenotype,nrow=1,scales = "free_x")+
  scale_color_aaas()+
  labs(title="Bootstrapped R2 (b=100)",y="Boot R2",x="Number of genes/PCs in model")+
  theme_minimal()

#
