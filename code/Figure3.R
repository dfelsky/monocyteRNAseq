### figure 3
library(rms)
library(reshape2)
library(ggsci)
library(ggpubr)
library(ggthemes)
library(ggrepel)
library(naturalsort)

### read TWAS results
ttall_unlabelled <- readRDS("output/all_TWAS_results.rds")
ttwide <- readRDS("output/all_TWAS_results_wideformat_labelled.rds")

### name files for plotting
pathnames <- c("Neuritic plaques","Diffuse plaques","Total AB","PHF tau","NFT","Gross cerebral infarcts","Micro cerebral infarcts","Arteriolosclerosis","Cerebral AA","Cerebral atherosclerosis","Lewy body stage","Hippocampal sclerosis","TDP-43","PD Dx","Patho AD","PAM VM Caudate","PAM post. putamen","PAM IT","PAM MF")
cognames <- c("Episodic memory","Perceptual orientation","Perceptual speed","Semantic memory","Working memory","Global","MMSE")

ttall_labelled <- lapply(ttall_unlabelled, function(x) {
  names(x) <- c(pathnames,cognames)
  x
})


### make hybrid chicago plot
phenos_to_plot <- c("Global","Cerebral_AA","Diffuse_plaques","Arteriolosclerosis","PAM_MF","PAM_IT")

plotlist <- list()
for (pheno in phenos_to_plot) {
  
  vars_to_get <- grep("adj.P.Val|_t_", grep(paste0(pheno,collapse = "|"),grep("MONO_BLOOD",names(ttwide),value=T),value=T),value=T)
  widesub <- ttwide[,c("gene","hugo","chr","location",vars_to_get)]
  z <- widesub[complete.cases(widesub),]
  
  z$value <- -log10(z[,paste0("MONO_BLOOD_adj.P.Val_",pheno)])*sign(z[,paste0("MONO_BLOOD_t_",pheno)])

  z$chr <- ifelse(z$chr=="X",23,z$chr)
  z$chr <- ifelse(z$chr=="MT",24,z$chr)
  z$transpos <- z$location + as.numeric(z$chr)*1e12
  z$transpos <- rank(z$transpos)
  z$chr <- factor(z$chr, levels=naturalsort(unique(z$chr)))
  
  z$chrcol <- ifelse(as.numeric(z$chr)%%2==0,1,0)
  
  xmins <- NULL
  xmaxs <- NULL
  startloc <- 1
  for (chr in c(1:24)) {
    xmins[chr] <- startloc
    xmaxs[chr] <- startloc + table(z$chr)[chr]
    startloc <- xmaxs[chr]
  }
  
  xmids <- xmins + (xmaxs-xmins)/2
  
  genes_toplot_df <- z[which(z[,paste0("MONO_BLOOD_adj.P.Val_",pheno)]<0.05),]
  genes_toplot_df <- genes_toplot_df[order(abs(genes_toplot_df$value),decreasing = T),]
  genes_toplot <- genes_toplot_df$hugo[1:20]
  
  if (pheno==last(phenos_to_plot)) {
    extralayer <- annotate("text",x=xmids,y=-6,label=c(1:22,"X","MT"),size=3)
  } else {
    extralayer <- geom_blank()
  }
  
  plotlist[[pheno]] <- ggplot(data=z,aes(x=transpos,y=value))+
    geom_point(aes(col=chrcol),show.legend = F,size=0.4)+
    annotate("rect",
             xmin=xmins,
             xmax=xmaxs,
             ymin=-Inf, 
             ymax=Inf,
             fill=rep(c("grey","white"),12),
             alpha=0.2)+
    extralayer +
    theme_classic()+
    geom_hline(yintercept = -log(0.05,base=10),col="orange",lty=2,size=0.5)+
    geom_hline(yintercept = log(0.05,base=10),col="orange",lty=2,size=0.5)+
    geom_hline(yintercept =0,col="blue",size=1)+
    geom_text_repel(data=subset(z,hugo %in% genes_toplot),
                    aes(label=paste0("italic(",hugo,")")),
                        size=2.3,
                        parse=TRUE)+
    labs(title=pheno)+
    theme(axis.line=element_blank(),
          axis.text.x=element_blank(),
          axis.ticks.x = element_blank(),
          axis.title.x = element_blank(),
          axis.title.y = element_blank())
}


plot_grid(plotlist[[1]],
          plotlist[[2]],
          plotlist[[3]],
          plotlist[[4]],
          plotlist[[5]],
          plotlist[[6]],
          ncol=1)



##############
#### plot volcano plots for PAM phenotypes and diffuse plaques and atherosclerosis

p1 <- ggplot(ttall_unlabelled$monocyte_blood$mf3123,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title = "A. PAM (MF)")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=ttall_unlabelled$monocyte_blood$mf3123[1:10,],aes(label=hugo))+
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