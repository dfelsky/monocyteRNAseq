## overlap between results
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
cognames <- c("Global","Episodic memory","Perceptual orientation","Perceptual speed","Semantic memory","Working memory","MMSE")

ttall_labelled <- lapply(ttall_unlabelled, function(x) {
  names(x) <- c(pathnames,cognames)
  x
})

###################
###################
####### cross-sample overlap per phenotype
commongenes <- intersect(ttall_unlabelled$monocyte$plaq_n_sqrt$gene,
                         ttall_unlabelled$dlpfc_cells$plaq_n_sqrt$gene)

ttm <- ttall_labelled$monocyte_blood
ttdc <- ttall_labelled$dlpfc_cells

index <- 1
stestlist <- list()
for (sigthres in c(0.05,0.1,0.2)) {
  for (pheno in c(pathnames,cognames)) {
    ttmp <- ttm[[pheno]][commongenes,]
    ttdp <- ttdc[[pheno]][commongenes,]
    
    monosig <- ttmp$gene[which(ttmp$adj.P.Val < sigthres)]
    dlpfcsig <- ttdp$gene[which(ttdp$adj.P.Val < sigthres)]
    
    if(length(monosig)==0 & length(dlpfcsig)==0) { 
      stestlist[[index]] <- data.frame(pheno=pheno,
                                       FDR=as.character(sigthres),
                                       dlpfcsig=0,
                                       monosig=0,
                                       overlap=NA,
                                       expected.overlap=NA,
                                       p=NA)
    } else {
      
      sig_sets <- list(mono=monosig,
                       dlpfc=dlpfcsig)
      
      sumobj <- summary(supertest(sig_sets,n=length(commongenes)))
      
      stestlist[[index]] <- data.frame(pheno=pheno,
                                       FDR=as.character(sigthres),
                                       dlpfcsig=as.numeric(sumobj$set.sizes["dlpfc"]),
                                       monosig=as.numeric(sumobj$set.sizes["mono"]),
                                       overlap=as.numeric(sumobj$otab[3]),
                                       expected.overlap=as.numeric(sumobj$etab[3]),
                                       p=as.numeric(sumobj$P.value[3]))
    }
    index <- index + 1
  }
}

setdata <- as.data.frame(do.call(rbind,stestlist))
setdata

setdata$pheno <- factor(setdata$pheno,levels=c(pathnames,cognames))

pdf("paper/figures/TWAS_result_HGoverlap_plot_monocyte_blood.pdf",h=5,w=8)
ggplot(data=setdata,aes(y=-log10(p),x=pheno,col=FDR))+
  geom_point(size=5)+
  scale_color_aaas()+
  coord_flip()+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text(aes(label=overlap),col="white",size=2,fontface="bold")+
  labs(y="Hypergeometric overlap significance (-log10(p))",x="Phenotype")+
  theme_hc()
dev.off()

##### cross-tissue effect correlations 
### read in wide effects
allwide <- readRDS("output/all_TWAS_results_wideformat_labelled.rds")
pathnames2 <- c("Neuritic_plaques","Diffuse_plaques","Total_AB","PHF_tau","NFT","Gross_cerebral_infarcts","Micro_Cerebral_infarcts","Arteriolosclerosis","Cerebral_AA","Cerebral_atherosclerosis","Lewy_body_stage","Hippocampal_sclerosis","TDP_43","PD_Dx","Patho_AD","PAM_VM_Caudate","PAM_post_putamen","PAM_IT","PAM_MF")
cognames2 <- c("Global","Episodic_memory","Perceptual_orientation","Perceptual_speed","Semantic_memory","Working_memory","MMSE")

allcordats <- list()
index <- 1
for (datasub in c("full","onlysigmono10")) {
  for (pheno in c(cognames2,pathnames2)) {
    mvar <- paste0("MONO_t_",pheno)
    dvar <- paste0("DLPFC_cells_t_",pheno)
    
    if (datasub=="full") {
      cordat <- allwide
      } else {
        cordat <- allwide[which(abs(allwide[,mvar]) >= quantile(abs(allwide[,mvar]),na.rm=T,0.9)),]  # subset is 90th percentile of absolute t-values for monocyte association
        }
    
      if (nrow(cordat) < 2) {
        print(paste0("too few observations to calculate subset correlation for pheno: ",pheno))
        allcordats[[index]] <- c(NA,p=NA,pheno=pheno,abs=abso,subset=datasub)
        index <- index+1
        } else {
          
          for (abso in c(TRUE,FALSE)) {
            if (abso==TRUE) {
              mdcor <- cor.test(abs(cordat[,mvar]),abs(cordat[,dvar]), method="spearman")
              } else {
                mdcor <- cor.test(cordat[,mvar],cordat[,dvar],method="spearman")
                }
            
            allcordats[[index]] <- c(mdcor$estimate,
                          p=mdcor$p.value,
                          pheno=pheno,
                          abs=abso,
                          subset=datasub)
            index <- index + 1
          }
        }
  }
}



aplotdat <- as.data.frame(do.call(rbind,allcordats))
aplotdat <- within(aplotdat,{
  rho <- as.numeric(as.character(rho))
  p <- as.numeric(as.character(p))
  pheno <- as.factor(pheno)
  abs <- as.character(abs)
  subset <- as.character(subset)
})

aplotdat$pheno <- factor(aplotdat$pheno,levels=c(pathnames2,cognames2))

ggplot(data=aplotdat, aes(y=rho,x=pheno,col=abs))+
  geom_point(aes(size=-log10(p),shape=abs))+
  geom_text_repel(data=subset(aplotdat,p < 0.05/length(unique(aplotdat$pheno))), aes(label=signif(p,digits=3)),nudge_x = 0.2)+
  geom_hline(yintercept = 0)+
  facet_wrap(~subset,ncol=2,scales="free")+
  scale_color_aaas()+
  coord_flip()+
  theme_minimal()


#### inspect scatterplots of individual correlations
phenosofinterest <- as.character(unique(aplotdat$pheno[which(aplotdat$p< 0.05/nrow(aplotdat))]))

plotlist1 <- list()
for (pheno in phenosofinterest) {
  mvar <- paste0("MONO_t_",pheno)
  dvar <- paste0("DLPFC_cells_t_",pheno)
  
  aestext <- paste("aes(x=",mvar,",y=",dvar,")")
  aeseval <- eval(parse(text=aestext))
  
  plotlist1[[pheno]] <- ggplot(data=allwide, aeseval)+
    geom_point()+
    geom_smooth(method="lm")+
    theme_minimal()
  
}
