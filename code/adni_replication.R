### plot ADNI replication
library(tidyr)
library(ggplot2)
library(ggthemes)
library(ggpubr)
library(ggsci)
library(ggrepel)
library(rms)

### read in ADNI results from Peter
adni <- read.csv("input/adni/adni_phenotypes_vs_rosmap.csv")
names(adni) <- gsub(names(adni),pattern="^P_",replacement="P__")
names(adni) <- gsub(names(adni),pattern="^T_",replacement="T__")
adni <- as_tibble(adni)

### read in ROSMAP results summary
ros <- read.csv("input/adni/FDR01_allresults_monocytes.csv")

### melt ADNI dataframe
adnim <- adni %>% 
  pivot_longer(
    !Gene,
    names_to = c(".value","pheno"),
    names_sep="__"
  )
head(adnim)
adnim <- as.data.frame(adnim)

### subset rosmap genes to phenotypes
# AMYLOID
pthres <- 0.05
ross <- subset(ros, adj.P.Val < pthres)
ross_amyloid <- subset(ross, phenotype %in% c("amyloid_sqrt","plaq_d_sqrt","plaq_n_sqrt","pathoAD"))
ross_amyloid$sample <- "rosmap"
# get matching amyloid phenos from adni
adni_amyloid <- adni[,c("Gene",grep("SUMMARYSUVR_WHOLECEREBNORM",names(adni),value = T))]
# make long format
names(adni_amyloid) <- gsub(names(adni_amyloid),pattern="^P_",replacement="P__")
names(adni_amyloid) <- gsub(names(adni_amyloid),pattern="^T_",replacement="T__")
adni_amyloid <- as_tibble(adni_amyloid)
adni_amyloidm <- adni_amyloid %>% 
  pivot_longer(
    !Gene,
    names_to = c(".value","pheno"),
    names_sep="__"
  )
adnim <- as.data.frame(adni_amyloidm)
adnim$sample <- "adni"
adnims <- subset(adnim, Gene %in% ross_amyloid$hugo)

#rbind
r1 <- ross_amyloid[,c(2,1,3,4,5)]
names(r1) <- names(adnims)

amy_merged <- rbind(r1,adnims)

# sort genes
amy_merged$Gene <- factor(amy_merged$Gene,
                          levels=names(sort(sapply(unique(amy_merged$Gene),function(x) mean(amy_merged[which(amy_merged$Gene==x),"T"])))))

#individual pheno plot
ggplot(dat=amy_merged, aes(y=-log10(P)*sign(T),x=Gene,col=pheno))+
  geom_hline(yintercept=c(-log10(0.1),0,log10(0.1)),
             col=c("orange","blue","orange"),
             lty=c(2,1,2))+
  scale_color_nejm()+
  geom_point(aes(shape=sample),size=3)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1,
                                   size = 12, hjust = 1))

# barplot samplewise
ggplot(dat=amy_merged, aes(y=-log10(P)*sign(T),x=Gene,col=sample,fill=sample))+
  geom_hline(yintercept=c(-log10(pthres),0,log10(pthres)),
             col=c("orange","blue","orange"),
             lty=c(2,1,2))+
  scale_color_nejm()+
  scale_fill_nejm()+
  geom_boxplot(width=0.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 12, hjust = 1))



# INFLAMMATORY
pthres <- 0.05
ross <- subset(ros, adj.P.Val < pthres)
ross_inf <- subset(ross, phenotype %in% c("it3123","mf3123","pput3123"))
ross_inf$sample <- "rosmap"
# get matching amyloid phenos from adni
adni_inf <- adni[,c("Gene",grep("TOTAL_WMH|TNFALPHA",names(adni),value = T))]
# make long format
names(adni_inf) <- gsub(names(adni_inf),pattern="^P_",replacement="P__")
names(adni_inf) <- gsub(names(adni_inf),pattern="^T_",replacement="T__")
adni_inf <- as_tibble(adni_inf)
adni_infm <- adni_inf %>% 
  pivot_longer(
    !Gene,
    names_to = c(".value","pheno"),
    names_sep="__"
  )
adnim <- as.data.frame(adni_infm)
adnim$sample <- "adni"
adnims <- subset(adnim, Gene %in% ross_inf$hugo)

#rbind
r1 <- ross_inf[,c(2,1,3,4,5)]
names(r1) <- names(adnims)

inf_merged <- rbind(r1,adnims)

# sort genes
inf_merged$Gene <- factor(inf_merged$Gene,
                          levels=names(sort(sapply(unique(inf_merged$Gene),function(x) mean(inf_merged[which(inf_merged$Gene==x),"T"])))))

#individual pheno plot
# ggplot(dat=inf_merged, aes(y=-log10(P)*sign(T),x=Gene,col=pheno))+
#   geom_hline(yintercept=c(-log10(0.1),0,log10(0.1)),
#              col=c("orange","blue","orange"),
#              lty=c(2,1,2))+
#   scale_color_nejm()+
#   geom_point(aes(shape=sample),size=3)+
#   theme_minimal()+
#   theme(axis.text.x = element_text(angle = 45, vjust = 1, 
#                                    size = 7, hjust = 1))

# barplot samplewise
ggplot(dat=inf_merged, aes(y=-log10(P)*sign(T),x=Gene,col=sample,fill=sample))+
  geom_hline(yintercept=c(-log10(pthres),0,log10(pthres)),
             col=c("orange","blue","orange"),
             lty=c(2,1,2))+
  scale_color_nejm()+
  scale_fill_nejm()+
  geom_boxplot(width=0.2)+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                   size = 7, hjust = 1))

# tile plot
inf_merged$pheno <- factor(inf_merged$pheno,levels=unique(inf_merged$pheno))

inf_merged2 <- inf_merged[which(inf_merged$pheno %nin% gsub("^._","",grep("cov_mono",names(adni),value=T))),]

ggplot(data=inf_merged, aes(x=Gene,y=pheno))+
  geom_tile(data=subset(inf_merged, P < 0.05),aes(fill=-log10(P)*sign(T)))+
  #geom_text(data=subset(inf_merged, P < 0.05),aes(label=signif(P,2)))+
  scale_fill_gradient2_tableau()+
  theme_minimal()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, 
                                 size = 7, hjust = 1))


