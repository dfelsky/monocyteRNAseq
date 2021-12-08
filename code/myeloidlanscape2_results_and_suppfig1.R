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
compres <- lapply(unique(fullres$test),function(test) {
  print(paste0("test: ",test))
  allpheno <- lapply(unique(ttml$pheno), function(pheno) {
    print(paste0("phenotype: ",pheno))
    
    df.tt <- ttml[which(ttml$pheno==pheno),]
    df.res <- fullres[which(fullres$test==test & fullres$human_ensembl_id %in% df.tt$gene),]
    
    df <- merge(df.tt,df.res,by.x="gene",by.y="human_ensembl_id")
    
      
      ml_sig_up <- df$gene[which(df$pvalue<0.05 & df$coefficient > 0)]
      ml_sig_down <- df$gene[which(df$pvalue<0.05 & df$coefficient < 0)]
      mono_sig_up <- df$gene[which(df$adj.P.Val < 0.05 & df$t > 0)]
      mono_sig_down <- df$gene[which(df$adj.P.Val < 0.05 & df$t < 0)]
      
      # allsets <- list(ml_sig_up=ml_sig_up,
      #                 ml_sig_down=ml_sig_down,
      #                 mono_sig_up=mono_sig_up,
      #                 mono_sig_down=mono_sig_down)
      
       allsets <- list(ml_sig=ml_sig,
                       mono_sig=mono_sig)
      
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
                   test=test)
      } else {

      
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
      sumtest$test <- test
      sumtest
      
      }
  })
  do.call(rbind,allpheno)
})

saveRDS(compres, "output/myeloidlanscape_compresults.rds")
compres <- readRDS("output/myeloidlanscape_compresults.rds")
cres <- do.call(rbind,compres)
cres$Elements <- NULL
cres <- subset(cres,Degree == 2)
#cres <- subset(cres, Intersections %in% c("ml_sig_up & mono_sig_up","ml_sig_down & mono_sig_down"))

cres2 <- cres[complete.cases(cres),]

subset(cres2, P.value < 0.05) %>%
  ggplot(aes(y=test,x=pheno,fill=FE))+
  scale_fill_gradient_tableau()+
  geom_tile()+
  geom_text(aes(label=Observed.Overlap))+
  facet_wrap(~Intersections)+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))

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

