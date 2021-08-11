### generate normalized monocyte gene expression data
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

### read in data with blood draw dates to calculate age at draw
mp <- read.csv("input/b123_pheno.txt",sep="\t",header=T,colClasses = c(projid="character"))
mp$projid <- str_pad(mp$projid, width=8, side="left", pad="0")

### load raw count and QC metric data
b1g <- read.table("input/raw/mono_batch1/genes_expectedCounts.txt",
                  row.names = "gene_ID",
                  header=T, 
                  sep="",
                  check.names = F)

b1m <- read.table("input/raw/mono_batch1/qualityMetrics.csv", check.names = F, sep="", row.names = "Sample", header=T,colClasses=c("Sample"="character"))

b2g <- read.table("input/raw/mono_batch2/rsem_genes_raw/genes_expectedCounts.txt", 
                  row.names = "gene_ID", 
                  header=T, 
                  sep=" ", 
                  check.names = F)

b2m <- read.table("input/raw/mono_batch2/qualityMetrics.csv", check.names = F, sep="", row.names = "Sample", header=T,colClasses=c("Sample"="character"))


### Merge batches and create batch index
counts <- cbind(b1g,b2g)
meta <- rbind(b1m,b2m)
meta$batch <- c(rep("batch1",ncol(b1g)),rep("batch2",ncol(b2g)))
meta <- meta[match(colnames(counts),rownames(meta)),]

### read in IDkey and assign projids to meta data
IDkey <- read.csv("input/All List_20171106v4-corrected.csv")
meta$projid <- IDkey$projid..individual.ID.[match(rownames(meta),IDkey$Library.ID.External.ID)]
meta$sampleID <- rownames(meta)

### merge in biological and experimental covariates
# from phenotype dataset
covs <- ROSmaster[,c("projid","msex","study","age_bl")]

# age at blood draw from mp dataframe
MP <- subset(mp,monoRNA_batch<3) # necessary because mp includes some repeated measures
drawvisit <- MP[,c("projid","visit")]

meta1 <- merge(meta,covs,by="projid",all.x=T)
meta2 <- merge(meta1,drawvisit,by="projid",all.x=T)
meta2$age_draw <- meta2$age_bl + meta2$visit

### exclude poor quality subjects from multiqc html reports
exclude_subjects <- c("20504017","50108912","31_E1","34_G1","30_B1","40_E1","42_E1","42_D7","43_E1") 
meta2 <- meta2[complete.cases(meta2),]
meta2 <- meta2[-which(meta2$sampleID %in% exclude_subjects),]
counts <- counts[,meta2$sampleID]
colnames(counts) <- meta2$projid
rownames(meta2) <- meta2$projid

### make DGElist object and format covariates correctly
dge <- DGEList(counts=counts, 
               samples=meta2,
               group = meta2$batch,
               lib.size = colSums(counts))

dge$samples <- within(dge$samples,{
  msex <- as.factor(msex)
  study <- as.factor(study)
  batch <- as.factor(batch)
})

### XIST sex check
xist <- "ENSG00000229807"
dge$samples$XIST <- cpm(dge$counts[xist, ], log=TRUE)
ggplot(dge$samples, aes(x=msex, y=XIST)) +
  geom_violin() + 
  geom_jitter(width=0.3, show.legend = FALSE)

outlier_xist <- dge$samples$projid[which((dge$samples$msex==0 & dge$samples$XIST < 7.5) | (dge$samples$msex==1 & dge$samples$XIST > 7.5))]

dge <- dge[,dge$samples$projid %nin% outlier_xist]

saveRDS(dge,file="input/monocytes_raw.rds")

### QC and forward confounder selection
source("code/bulkRNA_QC.R")
dge <- readRDS("input/monocytes_raw.rds")

qc <- bulk_qc(ydge=dge,
             min.count = 15,
             pc.bygroup = T,
             variance.prune = F, 
             high.prune = F,
             IQR_highprune_threshold = 4,
             IQR_obs_threshold = 8,
             IQR_sub_threshold = 4,
             coerce.obs = T)


fileConn<-file("input/monocytes_QCresults.txt")
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

lcpm <- cpm(dge_filtered,log=T)
allmeds <- apply(lcpm,2,median)
bp1 <- boxplot(allmeds,range = 3)
distribution_outliers <- names(which(allmeds < max(bp1$out)))
dge_filtered <- dge_filtered[,which(dge_filtered$samples$projid %nin% distribution_outliers)]

#saveRDS(dge_filtered,file="input/monocytes_filtered_only.rds")

### variable selection
dge_filtered <- readRDS("input/monocytes_filtered_only.rds")
dge_filtered <- calcNormFactors(dge_filtered)

bloodvars <- c("hemacrit_at_draw","hemoglbn_at_draw","mch_at_draw","mchc_at_draw","mcv_at_draw","platelet_at_draw","rbc_at_draw", "rdw_at_draw","wbc_at_draw","fasting.f","hemotologic_rx_at_draw")

dge_filtered$samples <- cbind(dge_filtered$samples,ROSmaster[match(dge_filtered$samples$projid,ROSmaster$projid),bloodvars])

covars <- c("PF_READS","PF_READS_ALIGNED","PCT_PF_READS_ALIGNED","PCT_RIBOSOMAL_BASES","PCT_USABLE_BASES","MEDIAN_3PRIME_BIAS","PERCENT_DUPLICATION","ESTIMATED_LIBRARY_SIZE","batch","msex","study","age_draw",bloodvars)
techvars <- c("PF_READS","PF_READS_ALIGNED","PCT_PF_READS_ALIGNED","PCT_RIBOSOMAL_BASES","PCT_USABLE_BASES","MEDIAN_3PRIME_BIAS","PERCENT_DUPLICATION","ESTIMATED_LIBRARY_SIZE","batch","study")

### plot correlations
corMatrix <- getCorrelations(dge_filtered$samples,
                             x=covars, 
                             y=covars)
ggcorrplot(corr=corMatrix$c, p.mat=corMatrix$p, tl.cex=6)

### Now PCA to check correlations of confounders with PCs
pcs <- getPCs.dge(dge_filtered, n = 20)
corMatrix <- getCorrelations(cbind(dge_filtered$samples, pcs),
                             x=covars, 
                             y=colnames(pcs))
ggcorrplot(corr=corMatrix$c, p.mat=corMatrix$p, tl.cex=7,sig.level = 0.05)

### perform variable selection on only techvars
inclVars <- c("batch")
swSelection <- detectSignCovariates(dge_filtered, varsFix=inclVars, varsOpt=techvars[techvars != "batch"], npc=20, maxIter=15, values="counts")
swSelection

### sanity check
covariates  <- as.character(swSelection$variable[swSelection$iteration == max(swSelection$iteration) & swSelection$selected])
resids <- getResiduals(dge_filtered, vars=covariates, values="counts")

pcs <- getPCs.matrix(resids$resids, n=20)
corMatrix <- getCorrelations(cbind(dge_filtered$samples, pcs),
                                 x=covars, 
                                 y=colnames(pcs))

ggcorrplot(corr=corMatrix$c, p.mat=corMatrix$p, tl.cex=7)


# covariates listed in covariates object from selection above
model <- model.matrix(~ batch +
                        PCT_USABLE_BASES +
                        PCT_PF_READS_ALIGNED +
                        PERCENT_DUPLICATION +
                        MEDIAN_3PRIME_BIAS +
                        ESTIMATED_LIBRARY_SIZE +
                        study,	
                      data=dge_filtered$samples)


# make voom object for modeling
v <- voomWithQualityWeights(dge_filtered, model, plot=TRUE)
fit <- lmFit(v, model, method="robust",maxit=10000)
resids <- residuals(fit, y=v)

### get residuals for tech and tech + age/sex
model.agesex <- model.matrix(~ batch +
                        PCT_USABLE_BASES +
                        PCT_PF_READS_ALIGNED +
                        PERCENT_DUPLICATION +
                        MEDIAN_3PRIME_BIAS +
                        ESTIMATED_LIBRARY_SIZE +
                        study +
                        msex +
                        age_draw,	
                      data=dge_filtered$samples)

v.agesex <- voomWithQualityWeights(dge_filtered, model.agesex, plot=TRUE)
fit.agesex <- lmFit(v.agesex, design=model.agesex, method="robust",maxit=10000)
resids.agesex <- residuals(fit.agesex, y=v.agesex)

### get residuals for tech and tech + age/sex + blood counts and hematologic diagnosis
# need to treat this a little different because some of these vars have missing values in the DGE object

bloodvars_tocovary <- bloodvars[c(2,4:6,9:11)]
projidstokeep <- dge_filtered$samples$projid[complete.cases(dge_filtered$samples[,bloodvars_tocovary])]
dge_filtered_sub <- dge_filtered[,projidstokeep]

model.agesex.blood <- model.matrix(~ batch +
                               PCT_USABLE_BASES +
                               PCT_PF_READS_ALIGNED +
                               PERCENT_DUPLICATION +
                               MEDIAN_3PRIME_BIAS +
                               ESTIMATED_LIBRARY_SIZE +
                               study +
                               msex +
                               age_draw +
                               hemotologic_rx_at_draw +
                               fasting.f +
                               hemoglbn_at_draw +
                               mcv_at_draw +
                               wbc_at_draw,
                             data=dge_filtered_sub$samples)



v.agesex.blood <- voomWithQualityWeights(dge_filtered_sub, model.agesex.blood, plot=TRUE)
fit.agesex.blood <- lmFit(v.agesex.blood, design=model.agesex.blood, method="robust",maxit=10000)
resids.agesex.blood <- residuals(fit.agesex.blood, y=v.agesex.blood)

rbeffect <- removeBatchEffect(v$E,batch = v$targets$batch,covariates = model[,-c(1,2)])
rbeffect.agesex <- removeBatchEffect(v.agesex$E,batch = v.agesex$targets$batch,covariates = model.agesex[,-c(1,2)])
rbeffect.agesex.blood <- removeBatchEffect(v.agesex.blood$E,batch = v.agesex.blood$targets$batch,covariates = model.agesex.blood[,-c(1,2)])

saveRDS(resids,file="input/monocytes_techonlynorm_residuals.rds")
saveRDS(rbeffect,file="input/monocytes_techonlynorm_removebatcheffectlogCPM.rds")
saveRDS(resids.agesex,file="input/monocytes_techandagesex_residuals.rds")
saveRDS(rbeffect.agesex,file="input/monocytes_techandagesex_removebatcheffectlogCPM.rds")
saveRDS(resids.agesex.blood,file="input/monocytes_techandagesexblood_residuals.rds")
saveRDS(rbeffect.agesex.blood,file="input/monocytes_techandagesexblood_removebatcheffectlogCPM.rds")
saveRDS(v,file="input/monocytes_v.rds")
saveRDS(v.agesex.blood,file="input/monocytes_v_blood.rds")
