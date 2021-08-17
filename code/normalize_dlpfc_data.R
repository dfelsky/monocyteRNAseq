### generate normalized dlpfc gene expression data
library(limma)
library(edgeR)
library(SummarizedExperiment)
library(ggplot2)
library(cowplot)
library(ggstatsplot)
library(stringr)

### read in custom functions (Klein)
source("code/ggcorrplot.R")
source("code/identifyCovariates.R")

### ROSMAP phenotype data
ROSmaster <- readRDS("input/ROSmaster_TWAS_input.rds")

# load raw count and QC metric data from Hans and Annie's run
dlpfc_raw <- as.matrix(read.table("input/raw/dlpfc/geneCounts.txt",sep="",header=T,row.names = "gene_ID",check.names = F))
dlpfc_meta <- read.table("input/raw/dlpfc/qualityMetrics.csv",sep="",header=T,colClasses = c(Sample="character"))

# make DGElist object
dge <- DGEList(counts=dlpfc_raw, 
               samples=dlpfc_meta,
               group = dlpfc_meta$Batch,
               lib.size = colSums(dlpfc_raw))

# remove RIN from metadata so that complete.cases all equal TRUE
dge$samples$RinScore <- NULL

# rename Sample column as projid and, Batch as batch
names(dge$samples)[which(names(dge$samples)=="Sample")] <- "projid"
names(dge$samples)[which(names(dge$samples)=="Batch")] <- "batch"

# add msex, study, and pmi to dge object
dge$samples <- cbind(dge$samples,ROSmaster[match(dge$samples$projid,ROSmaster$projid),c("msex","pmi","study","age_death")])
dge$samples <- within(dge$samples,{
  msex <- as.factor(msex)
  study <- as.factor(study)
  batch <- as.factor(batch)
  LOG_ESTIMATED_LIBRARY_SIZE <- log(ESTIMATED_LIBRARY_SIZE)
  LOG_PF_READS_ALIGNED <- log(PF_READS_ALIGNED)
})


# take only complete.cases
dge <- dge[,complete.cases(dge$samples)]


### XIST sex check
xist <- "ENSG00000229807"
dge$samples$XIST <- cpm(dge$counts[xist, ], log=TRUE)
ggplot(dge$samples, aes(x=msex, y=XIST)) +
  geom_violin() + 
  geom_jitter(width=0.3, show.legend = FALSE)

outlier_xist <- dge$samples$projid[which((dge$samples$msex==0 & dge$samples$XIST < 5) | (dge$samples$msex==1 & dge$samples$XIST > 5))]

dge <- dge[,dge$samples$projid %nin% outlier_xist]

saveRDS(dge,file="input/dlpfc_raw.rds")

### QC and forward confounder selection
source("code/bulkRNA_QC.R")
dge <- readRDS("input/dlpfc_raw.rds")

qc <- bulk_qc(ydge=dge,
              min.count = 15,
              pc.bygroup = T,
              variance.prune = F, 
              high.prune = F,
              IQR_highprune_threshold = 4,
              IQR_obs_threshold = 8,
              IQR_sub_threshold = 4,
              coerce.obs = T)


fileConn<-file("input/dlpfc_QCresults.txt")
writeLines(
  c(paste(t(qc$returnparams)[,1],collapse=": "),
    paste(t(qc$returnparams)[,2],collapse=": "),
    paste(t(qc$returnparams)[,3],collapse=": "),
    paste(t(qc$returnparams)[,4],collapse=": "),
    paste(t(qc$returnparams)[,5],collapse=": "),
    paste(t(qc$returnparams)[,6],collapse=": "),
    paste(t(qc$returnparams)[,7],collapse=": "),
    paste(t(qc$returnparams)[,8],collapse=": "),
    paste(t(qc$returnparams)[,9],collapse=": "),
    paste(t(qc$returnparams)[,10],collapse=": "),
    paste(qc$returnprune),
    paste(qc$returnprunehigh),
    paste(qc$returnprunevar),
    paste(qc$returnsub),
    paste(qc$timestamp)),
  fileConn)
close(fileConn)

dge_filtered <- qc$YDGE

saveRDS(dge_filtered,file="input/dlpfc_filtered_only.rds")


### variable selection
dge_filtered <- readRDS("input/dlpfc_filtered_only.rds")
dge_filtered <- calcNormFactors(dge_filtered)

covars <- c("PF_READS","PF_READS_ALIGNED","PCT_PF_READS_ALIGNED","PCT_RIBOSOMAL_BASES","PCT_CODING_BASES","PCT_UTR_BASES","PCT_INTRONIC_BASES","PCT_USABLE_BASES","PCT_INTERGENIC_BASES","PCT_MRNA_BASES","MEDIAN_3PRIME_BIAS","MEDIAN_CV_COVERAGE","MEDIAN_5PRIME_BIAS","PERCENT_DUPLICATION","LOG_ESTIMATED_LIBRARY_SIZE","LOG_PF_READS_ALIGNED","ESTIMATED_LIBRARY_SIZE","batch","msex","study","pmi","age_death")
techvars <- c("PF_READS","PF_READS_ALIGNED","PCT_PF_READS_ALIGNED","PCT_RIBOSOMAL_BASES","PCT_CODING_BASES","PCT_UTR_BASES","PCT_INTRONIC_BASES","PCT_USABLE_BASES","PCT_INTERGENIC_BASES","PCT_MRNA_BASES","MEDIAN_3PRIME_BIAS","MEDIAN_CV_COVERAGE","MEDIAN_5PRIME_BIAS","PERCENT_DUPLICATION","LOG_PF_READS_ALIGNED","LOG_ESTIMATED_LIBRARY_SIZE","ESTIMATED_LIBRARY_SIZE","batch","study","pmi")

### plot correlations
corMatrix <- getCorrelations(dge_filtered$samples,
                             x=covars, 
                             y=covars)
pdf("paper/figures/DLPFC_covariate_correlations.pdf",w=8,h=8)
ggcorrplot(corr=corMatrix$c, p.mat=corMatrix$p, tl.cex=6,sig.level = 0.05,lab = T,lab_size = 2,title=paste0("DLPFC, n=",nrow(dge_filtered$samples)))
dev.off()

### Now PCA to check correlations of confounders with PCs
pcs <- getPCs.dge(dge_filtered, n = 20)
corMatrix <- getCorrelations(cbind(dge_filtered$samples, pcs),
                             x=covars, 
                             y=colnames(pcs))
pdf("paper/figures/DLPFC_covariate_PC_effects.pdf",w=8,h=8)
ggcorrplot(corr=corMatrix$c, p.mat=corMatrix$p, tl.cex=7,sig.level = 0.05,lab = T,lab_size = 2,title=paste0("DLPFC, n=",nrow(dge_filtered$samples)))
dev.off()

### perform variable selection on only techvars
inclVars <- c("batch")
swSelection <- detectSignCovariates(dge_filtered, varsFix=inclVars, varsOpt=techvars[techvars != "batch"], npc=20, maxIter=50, values="counts")
swSelection

### sanity check
covariates  <- as.character(swSelection$variable[swSelection$iteration == max(swSelection$iteration) & swSelection$selected])
model.agesex <- model.matrix(~ batch +
                               MEDIAN_CV_COVERAGE+
                               PCT_RIBOSOMAL_BASES+
                               PCT_CODING_BASES+
                               PCT_UTR_BASES+
                               LOG_ESTIMATED_LIBRARY_SIZE+
                               LOG_PF_READS_ALIGNED+
                               MEDIAN_5PRIME_TO_3PRIME_BIAS+
                               PCT_PF_READS_ALIGNED+
                               study+
                               PERCENT_DUPLICATION+
                               MEDIAN_3PRIME_BIAS+
                               PCT_INTERGENIC_BASES+
                               pmi+
                               age_death+
                               msex,
                             data=dge_filtered$samples)

v.agesex <- voom(dge_filtered, model.agesex, plot=TRUE)
fit.agesex <- lmFit(v.agesex, design=model.agesex, method="robust",maxit=10000)
resids.agesex <- residuals(fit.agesex, y=v.agesex)

pcs <- getPCs.matrix(resids.agesex, n=20)
corMatrix <- getCorrelations(cbind(dge_filtered$samples, pcs),
                             x=covars, 
                             y=colnames(pcs))
pdf("paper/figures/DLPFC_covariate_PC_effects_aftercorrection.pdf",w=8,h=8)
ggcorrplot(corr=corMatrix$c, p.mat=corMatrix$p, tl.cex=7,sig.level = 0.05,lab = T,lab_size = 2,title=paste0("DLPFC, n=",nrow(dge_filtered$samples)))
dev.off()


# covariates listed in covariates object from selection above

model <- model.matrix(~ batch +
                        MEDIAN_CV_COVERAGE+
                        PCT_RIBOSOMAL_BASES+
                        PCT_CODING_BASES+
                        PCT_UTR_BASES+
                        LOG_ESTIMATED_LIBRARY_SIZE+
                        LOG_PF_READS_ALIGNED+
                        MEDIAN_5PRIME_TO_3PRIME_BIAS+
                        PCT_PF_READS_ALIGNED+
                        study+
                        PERCENT_DUPLICATION+
                        MEDIAN_3PRIME_BIAS+
                        PCT_INTERGENIC_BASES+
                        pmi,
                      data=dge_filtered$samples)


# make voom object for modeling
v <- voom(dge_filtered, model, plot=TRUE)
fit <- lmFit(v, model, method="robust",maxit=10000)
resids <- residuals(fit, y=v)

### get residuals for tech and tech + age/sex
model.agesex <- model.matrix(~ batch +
                               MEDIAN_CV_COVERAGE+
                               PCT_RIBOSOMAL_BASES+
                               PCT_CODING_BASES+
                               PCT_UTR_BASES+
                               LOG_ESTIMATED_LIBRARY_SIZE+
                               LOG_PF_READS_ALIGNED+
                               MEDIAN_5PRIME_TO_3PRIME_BIAS+
                               PCT_PF_READS_ALIGNED+
                               study+
                               PERCENT_DUPLICATION+
                               MEDIAN_3PRIME_BIAS+
                               PCT_INTERGENIC_BASES+
                               pmi+
                               age_death+
                               msex,
                             data=dge_filtered$samples)

v.agesex <- voom(dge_filtered, model.agesex, plot=TRUE)
fit.agesex <- lmFit(v.agesex, design=model.agesex, method="robust",maxit=10000)
resids.agesex <- residuals(fit.agesex, y=v.agesex)

rbeffect <- removeBatchEffect(v$E,covariates = model[,-1])
rbeffect.agesex <- removeBatchEffect(v.agesex$E,covariates = model.agesex[,-1])

saveRDS(resids,file="input/dlpfc_techonlynorm_residuals.rds")
saveRDS(rbeffect,file="input/dlpfc_techonlynorm_removebatcheffectlogCPM.rds")
saveRDS(resids.agesex,file="input/dlpfc_techandagesex_residuals.rds")
saveRDS(rbeffect.agesex,file="input/dlpfc_techandagesex_removebatcheffectlogCPM.rds")
saveRDS(v,file="input/dlpfc_v.rds")