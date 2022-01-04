### myeloid lanscape results, Supp Figure 1
library(rms)
library(ggsci)
library(ggthemes)
library(ggpubr)
library(dplyr)

fullres <- read.csv("input/myeloidlandscape2/merged_all_monocyte_genes_ensembl.csv")
nameref <- readRDS("output/variable_name_reference.rds")
ttall <- readRDS("output/all_TWAS_results.rds")
ttm <- ttall$monocyte_blood
ttml <- do.call(rbind,ttm)[,c("pheno","gene","adj.P.Val","t")]

################## too long to run?...
donttest <- c("niareagansc: High - Low |  {includeLowExpr}")

testlist <- unique(fullres$test)[which(unique(fullres$test) %nin% donttest)]

compres <- lapply(testlist,function(test) {
  print(paste0("test: ",test))
  allpheno <- lapply(unique(ttml$pheno), function(pheno) {
    print(paste0("phenotype: ",pheno))
    
    df.tt <- ttml[which(ttml$pheno==pheno),]
    df.res <- fullres[which(fullres$test==test & fullres$human_ensembl_id %in% df.tt$gene),]
    
    df <- merge(df.tt,df.res,by.x="gene",by.y="human_ensembl_id")
    
    ml_sig <- df$gene[which(df$pvalue<0.05)]
    mono_sig <- df$gene[which(df$adj.P.Val < 0.05)]
    allsets <- list(ml_sig=ml_sig,mono_sig=mono_sig)
    
    if(length(unlist(allsets))==0) {
      data.frame(Intersections=NA,
                 Degree=NA,
                 Observed.Overlap=NA,
                 Expected.Overlap=NA,
                 FE=NA,
                 P.value=NA,
                 Elements=NA,
                 pheno=pheno,
                 experiment=df$experimentID[1],
                 ngenes=nrow(df),
                 test=test,
                 corp=NA,
                 corr=NA,
                 mlup_monoup=NA,
                 mlup_monodown=NA,
                 mldown_monoup=NA,
                 mldown_monodown=NA,
                 fp=NA)
      
      } else {
        stest <- SuperExactTest::supertest(allsets,n=nrow(df))
        sumtest <- summary(stest)$Table
        sumtest$pheno <- pheno
        sumtest$experiment <- df$experimentID[1]
        sumtest$ngenes <- nrow(df)
        sumtest$test <- test
        
        ctest <- cor.test(df$coefficient,df$t,method = "pearson")
        sumtest$corp <- ctest$p.value
        sumtest$corr <- ctest$estimate
        
        if (sumtest[which(sumtest$Degree==2),"P.value"] < 0.05) {
          
          df$ml_sigdir <- df$coefficient
          df$ml_sigdir[which(df$pvalue >= 0.05)] <- NA
          df$ml_sigdir <- ifelse(df$ml_sigdir > 0,"up","down")
          df$mono_sigdir <- df$t
          df$mono_sigdir[which(df$adj.P.Val >= 0.05)] <- NA
          df$mono_sigdir <- ifelse(df$mono_sigdir > 0,"up","down")
          
          sumtest$mlup_monoup <- table(df$mono_sigdir,df$ml_sigdir)[2,2]
          sumtest$mlup_monodown <- table(df$mono_sigdir,df$ml_sigdir)[1,2]
          sumtest$mldown_monoup <- table(df$mono_sigdir,df$ml_sigdir)[2,1]
          sumtest$mldown_monodown <- table(df$mono_sigdir,df$ml_sigdir)[1,1]
          
          ftest <- fisher.test(table(df$mono_sigdir,df$ml_sigdir))
          sumtest$fp <- ftest$p.value
          
        } else {

          sumtest$mlup_monoup <- NA
          sumtest$mlup_monodown <- NA
          sumtest$mldown_monoup <- NA
          sumtest$mldown_monodown <- NA
          sumtest$fp <- NA
          
          }
        sumtest
        }
    })
  do.call(rbind,allpheno)
  })


saveRDS(compres, "output/myeloidlanscape_compresults_nodirection.rds")
compres <- readRDS("output/myeloidlanscape_compresults_nodirection.rds")
cres <- do.call(rbind,compres)
cres <- subset(cres,Degree == 2)
cres$Elements <- NULL

#cres <- subset(cres, Intersections %in% c("ml_sig_up & mono_sig_up","ml_sig_down & mono_sig_down"))

cres2 <- cres[complete.cases(cres),] # this subsets for p<0.05 automatically
cres2 <- cres
cres2$label <- paste0(cres2$experiment," ... ",cres2$test)
cres2 <- subset(cres2, P.value < 0.05)

library(scatterpie)
library(ggrepel)

# subset(cres2, P.value < 0.005) %>%
#   ggplot(aes(y=label,x=pheno,fill=FE))+
#   scale_fill_gradient_tableau()+
#   geom_tile()+
#   geom_text(aes(label=signif(P.value,3)))+
#   facet_wrap(~Intersections)+
#   theme_minimal()+
#   theme(axis.text.x=element_text(angle = -45, hjust = 0))
# 
# subset(cres2, corp < 0.0005) %>%
#   ggplot(aes(y=label,x=pheno,fill=FE))+
#   scale_fill_gradient_tableau()+
#   geom_tile()+
#   geom_text(aes(label=signif(P.value,3)))+
#   facet_wrap(~Intersections)+
#   theme_minimal()+
#   theme(axis.text.x=element_text(angle = -45, hjust = 0))

cres2$label_num <- as.numeric(as.factor(cres2$label))
cres2$pheno_num <- as.numeric(as.factor(cres2$pheno))

phenoindex <- cres2$pheno_num
names(phenoindex) <- nameref$varnames[match(cres2$pheno,nameref$mono.variable)]
phenoindex <- phenoindex[!duplicated(phenoindex)]

labelindex <- cres2$label_num
names(labelindex) <- cres2$label
labelindex <- labelindex[!duplicated(labelindex)]

cres2$pielab <- paste0(cres2$Observed.Overlap," \n (",signif(cres2$fp,2),")")

tiff(filename = "paper/supp_figures/SuppFig1_myeloidlandscape_scatterpie_fp.tif",res = 300,units = "in",w=15,h=15)
ggplot(data=cres2,aes(y=label_num,x=pheno_num))+
  geom_scatterpie(aes(y=label_num,x=pheno_num),data=cres2,cols=names(cres2)[13:16],pie_scale = 2,color=NA)+
  scale_fill_few()+
  coord_equal()+
  geom_text(aes(label=signif(P.value,3)),size=2,nudge_x = 0.5,nudge_y = -0.4)+
  scale_y_continuous(breaks=seq(1,length(unique(cres2$label_num))),labels=names(labelindex)[match(seq(1,length(unique(cres2$label_num))),labelindex)])+
  scale_x_continuous(breaks=seq(1,length(unique(cres2$pheno_num))),labels=names(phenoindex)[match(seq(1,length(unique(cres2$pheno_num))),phenoindex)])+
  theme_minimal()+
  labs(y="Experiment...contrast",x="Phenotype")+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))
dev.off()

pdf("paper/supp_figures/SuppFig1_myeloidlandscape_scatterpie_fp.pdf",w=15,h=15)
ggplot(data=cres2,aes(y=label_num,x=pheno_num))+
  geom_scatterpie(aes(y=label_num,x=pheno_num),data=cres2,cols=names(cres2)[13:16],pie_scale = 2,color=NA)+
  scale_fill_few()+
  coord_equal()+
  geom_text(aes(label=pielab),size=2,nudge_x = 0.5,nudge_y = -0.4)+
  scale_y_continuous(breaks=seq(1,length(unique(cres2$label_num))),labels=names(labelindex)[match(seq(1,length(unique(cres2$label_num))),labelindex)])+
  scale_x_continuous(breaks=seq(1,length(unique(cres2$pheno_num))),labels=names(phenoindex)[match(seq(1,length(unique(cres2$pheno_num))),phenoindex)])+
  theme_minimal()+
  labs(y="Experiment...contrast",x="Phenotype")+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))
dev.off()

##############################################################################
########################## TRY T-STAT CORRELATIONS  INSTEAD OF PVALUE OVERLAP
res <- read.csv("input/myeloidlandscape2/merged_df_mono_myeloid3.csv")
nameref <- readRDS("output/variable_name_reference.rds")

compres <- lapply(names(which(table(res$experimentID)>1)),function(experiment) {
  allpheno <- lapply(unique(res$phenotype), function(pheno) {
    df <- subset(res, phenotype==pheno & experimentID==experiment)
    if(nrow(df) < 10) {  
      data.frame(p=NA,
                 cor=NA,
                 pheno=pheno,
                 experiment=experiment,
                 ngenes=nrow(df))
    } else {
      
      ctest <- cor.test(df$coefficient,df$t_statistic_mono,method = "pearson")
      data.frame(p=ctest$p.value,
                 cor=ctest$estimate,
                 pheno=pheno,
                 experiment=experiment,
                 ngenes=nrow(df))
    }
  })
  do.call(rbind,allpheno)
})
cres <- do.call(rbind,compres)

cres2 <- cres[complete.cases(cres),]

cres2$phenof <- nameref$varnames[match(cres2$pheno,nameref$mono.variable)]
#cres2$phenof <- paste0(cres2$phenof," (",cres2$ngenes,")")

tiff(file="paper/supp_figures/SuppFig1_myeloidlandscape_correlation.tif",w=6,h=5,units="in",res=300)
subset(cres2, p < 0.01) %>%
  ggplot(aes(x=phenof,y=experiment,fill=cor))+
  scale_fill_gradient2_tableau()+
  geom_tile()+
  geom_text(aes(label=paste(signif(p,2),"(",ngenes,")")),size=2)+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))
dev.off()



#######
df <- res[which(res$experimentID=="GSE68376" & res$phenotype=="it3123"),]
plot(df$coefficient,df$t_statistic_mono)
cor.test(df$coefficient,df$t_statistic_mono)

df <- res[which(res$experimentID=="GSE73721" & res$phenotype=="arteriol_scler"),]
ggplot(data=df,aes(x=coefficient,y=t_statistic_mono))+
  geom_point()+
  geom_smooth(method="lm")+
  geom_hline(yintercept = 0)+
  geom_vline(xintercept = 0)+
  theme_minimal()
cor.test(df$coefficient,df$t_statistic_mono)

