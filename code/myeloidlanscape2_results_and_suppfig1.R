### myeloid lanscape results, Supp Figure 1
library(rms)
library(ggsci)
library(ggthemes)
library(ggpubr)
library(dplyr)

res <- read.csv("input/myeloidlandscape2/merged_df_mono_myeloid3.csv")
nameref <- readRDS("output/variable_name_reference.rds")


compres <- lapply(names(which(table(res$experimentID)>1)),function(experiment) {
   allpheno <- lapply(unique(res$phenotype), function(pheno) {
    df <- subset(res, phenotype==pheno & experimentID==experiment)
    if(nrow(df)==0) {  
      data.frame(Intersections=NA,
                 Degree=NA,
                 Observed.Overlap=NA,
                 Expected.Overlap=NA,
                 FE=experiment,
                 P.value=NA,
                 Elements=NA,
                 pheno=pheno,
                 experiment=experiment,
                 ngenes=nrow(df))
      } else {
    
    ml_sig_up <- df$geneID[which(df$pvalue<0.05 & df$coefficient > 0)]
    ml_sig_down <- df$geneID[which(df$pvalue<0.05 & df$coefficient < 0)]
    # ml_sig <- df$geneID[which(df$pvalue<0.05)]
    # mono_sig <- df$geneID[which(df$mono_blood_p_value < 0.05)]
    mono_sig_up <- df$geneID[which(df$mono_blood_p_value < 0.05 & df$t_statistic_mono > 0)]
    mono_sig_down <- df$geneID[which(df$mono_blood_p_value < 0.05 & df$t_statistic_mono < 0)]
    
    allsets <- list(ml_sig_up=ml_sig_up,
                    ml_sig_down=ml_sig_down,
                    mono_sig_up=mono_sig_up,
                    mono_sig_down=mono_sig_down)
    
    # allsets <- list(ml_sig=ml_sig,
    #                 mono_sig=mono_sig)
    
    stest <- SuperExactTest::supertest(allsets,n=nrow(df))
    sumtest <- summary(stest)$Table
    # data.frame(p=sumtest["11","P.value"],
    #            overlap=sumtest["11","Observed.Overlap"],
    #            FE=sumtest["11","FE"],
    #            pheno=pheno,
    #            experiment=experiment)
    sumtest$pheno <- pheno
    sumtest$experiment <- experiment
    sumtest$ngenes <- nrow(df)
    sumtest
    
      }
  })
   do.call(rbind,allpheno)
})
cres <- do.call(rbind,compres)
cres$Elements <- NULL
cres <- subset(cres,Degree > 1)
cres <- subset(cres, Intersections %in% c("ml_sig_up & mono_sig_up","ml_sig_down & mono_sig_down"))

cres2 <- cres[complete.cases(cres),]

subset(cres2, p < 0.001) %>%
ggplot(aes(y=pheno,x=experiment,fill=FE))+
  scale_fill_gradient_tableau()+
  geom_tile()+
  geom_text(aes(label=overlap))+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))


##############################################################################
########################## TRY T-STAT CORRELATIONS  INSTEAD OF PVALUE OVERLAP
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

subset(cres2, p < 0.05) %>%
  ggplot(aes(x=phenof,y=experiment,fill=cor))+
  scale_fill_gradient2_tableau()+
  geom_tile()+
  geom_text(aes(label=paste(signif(p,2),"(",ngenes,")")))+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))
