library(corrplot)
library(rms)
library(BRETIGEA)

## pathology
ttd <- readRDS("output/TWAS_results_dlpfc_pathology.rds")
ttdc_raw <- readRDS("output/TWAS_results_dlpfc_pathology_celltypes.rds")
ttdc <- lapply(ttdc_raw,function(tt){
  tt <- tt[,c(2:16)]
  tt <- tt[!duplicated(tt[,"gene"]),]
  rownames(tt) <- tt[,"gene"]
  tt
})
ttm <- readRDS("output/TWAS_results_monocyte_pathology.rds")

ttdp <- ttd
ttdcp <- ttdc
ttmp <- ttm

nameset1 <- c("Neuritic plaques","Diffuse plaques","Total AB","PHF tau","NFT","Gross cerebral infarcts","Micro Cerebral infarcts","Arteriolosclerosis","Cerebral AA","Cerebral atherosclerosis","Lewy body stage","Hippocampal sclerosis","TDP-43","PD Dx","Patho AD","PAM VM Caudate","PAM post. putamen","PAM IT","PAM MF")


## cognition
ttd <- readRDS("output/TWAS_results_dlpfc_cognition.rds")
ttdc_raw <- readRDS("output/TWAS_results_dlpfc_cognition_celltypes.rds")
ttdc <- lapply(ttdc_raw,function(tt){
  tt <- tt[,c(2:16)]
  tt <- tt[!duplicated(tt[,"gene"]),]
  rownames(tt) <- tt[,"gene"]
  tt
})
ttm <- readRDS("output/TWAS_results_monocyte_cognition.rds")

##### MERGE pathology and cognition results

ttm <- c(ttmp,ttm)
ttd <- c(ttdp,ttd)
ttdc <- c(ttdcp,ttdc)

nameset2 <- c("Global","Episodic memory","Perceptual orientation","Perceptual speed","Semantic Memory","Working Memory","MMSE")

allnames <- c(nameset1,nameset2)

## only use for clarity of plotting
names(ttm) <- allnames
names(ttd) <- allnames
names(ttdc) <- allnames

# collect summary of results per tissue
monores <- lapply(ttm,function(x) { describe(x$adj.P.Val<0.05) })
monores.sug <- lapply(ttm,function(x) { describe(x$adj.P.Val<0.1) })
dlpfcres <- lapply(ttd,function(x) { describe(x$adj.P.Val<0.05) })
dlpfcres.sug <- lapply(ttd,function(x) { describe(x$adj.P.Val<0.1) })
dlpfcres.cells <- lapply(ttdc,function(x) { describe(x$adj.P.Val<0.05) })
dlpfcres.sug.cells <- lapply(ttdc,function(x) { describe(x$adj.P.Val<0.1) })

sigresdc <- lapply(dlpfcres.cells,function(x) { x$values$frequency[2]})
sugresdc <- lapply(dlpfcres.sug.cells,function(x) { x$values$frequency[2]})
sigresd <- lapply(dlpfcres,function(x) { x$values$frequency[2]})
sugresd <- lapply(dlpfcres.sug,function(x) { x$values$frequency[2]})
sigresm <- lapply(monores,function(x) { x$values$frequency[2]})
sugresm <- lapply(monores.sug,function(x) { x$values$frequency[2]})

sigresdc <- as.data.frame(do.call(rbind,sigresdc))
sigresdc$tissue <- "dlpfc.cells"
sigresdc$sig <- "FDR_05"
sigresdc$pheno <- rownames(sigresdc)
sigresd <- as.data.frame(do.call(rbind,sigresd))
sigresd$tissue <- "dlpfc"
sigresd$sig <- "FDR_05"
sigresd$pheno <- rownames(sigresd)
sigresm <- as.data.frame(do.call(rbind,sigresm))
sigresm$tissue <- "mono"
sigresm$sig <- "FDR_05"
sigresm$pheno <- rownames(sigresm)

sugresdc <- as.data.frame(do.call(rbind,sugresdc))
sugresdc$tissue <- "dlpfc.cells"
sugresdc$sig <- "FDR_10"
sugresdc$pheno <- rownames(sugresdc)
sugresd <- as.data.frame(do.call(rbind,sugresd))
sugresd$tissue <- "dlpfc"
sugresd$sig <- "FDR_10"
sugresd$pheno <- rownames(sugresd)
sugresm <- as.data.frame(do.call(rbind,sugresm))
sugresm$tissue <- "mono"
sugresm$sig <- "FDR_10"
sugresm$pheno <- rownames(sugresm)

sigres <- rbind(sigresd,sigresdc,sigresm,sugresd,sugresdc,sugresm)
sigres[is.na(sigres)] <- 0
names(sigres)[1] <- "nsig"

sigres$pheno <- factor(sigres$pheno,levels=names(ttm)[c(1:5,15,9,6:8,10,12:13,11,14,16:19,20:26)])

ggplot(data=subset(sigres,tissue!="dlpfc"),aes(y=nsig,x=pheno,fill=sig))+
  geom_bar(stat = "identity",position = "identity",col="black",alpha=0.6)+
  facet_wrap(~tissue,nrow=3,scales = "free")+
  scale_fill_aaas()+
  theme_classic()+
  coord_flip()+
  labs(y="Number of significant genes",x="Neuropathology")+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))