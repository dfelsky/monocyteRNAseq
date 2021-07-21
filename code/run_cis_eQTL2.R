# Prep for matrixeQTL (adapted from "http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html")
library(MatrixEQTL)
library(stringr)
library(rms)
library(data.table)
library(scales)
library(parallel)
library(BootstrapQTL)
#library(TreeQTL) ### problem

### set directory and load reference data
base.dir = "/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/eQTL"
load("/Users/dfelsky/Documents/data/all_genes_ensembl.RData")
source("/Users/dfelsky/Documents/scripts/make_ROSmaster_feb12019.R")

### static files - read in once and assess subject overlap in loop
# expression
norm_dlpfc <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/normalization/dlpfc_techonlynorm_residuals.rds")

norm_mono <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/normalization/monocytes_techonlynorm_residuals.rds")


sub_list_dlpfcexp <- colnames(norm_dlpfc)
gene_list_dlpfcexp <- rownames(norm_dlpfc)

sub_list_monoexp <- colnames(norm_mono)
gene_list_monoexp <- rownames(norm_mono)

# retrieve top 10 genomic PCs
ROSMAP_eigen <- read.table("/Users/dfelsky/Documents/Amin_PRS_study/Data_files/PCA_Geno/geno_qc.eigenvec.txt",header=F,sep="")



# covariates 
varsofinterest <- c("msex","age_death","pmi")
covfile <- ROSmaster[,c("projid",varsofinterest)]
covfile <- covfile[which(complete.cases(covfile)),]
names(covfile)[1] <- "id"
sub_list_covfile <- covfile$id

### gene location files per chromosome
for (chr in seq(1,22)) {
geneloc_dlpfc <- all_genes[which(all_genes$ensembl_gene_id %in% gene_list_dlpfcexp),]
geneloc_dlpfc <- geneloc_dlpfc[,c(1,5,3,4)]
names(geneloc_dlpfc) <- c("geneid","chr","left","right")
geneloc_dlpfc <- geneloc_dlpfc[which(geneloc_dlpfc$chr==chr),]
geneloc_dlpfc$chr <- paste("chr",geneloc_dlpfc$chr,sep="")
write.table(geneloc_dlpfc,file=paste(base.dir,"/matrix_eqtl_input_files/geneloc_chr",chr,"_dlpfc.txt",sep=""),quote=F,col.names = T,row.names = F,sep="\t")

geneloc_mono <- all_genes[which(all_genes$ensembl_gene_id %in% gene_list_monoexp),]
geneloc_mono <- geneloc_mono[,c(1,5,3,4)]
names(geneloc_mono) <- c("geneid","chr","left","right")
geneloc_mono <- geneloc_mono[which(geneloc_mono$chr==chr),]
geneloc_mono$chr <- paste("chr",geneloc_mono$chr,sep="")
write.table(geneloc_mono,file=paste(base.dir,"/matrix_eqtl_input_files/geneloc_chr",chr,"_mono.txt",sep=""),quote=F,col.names = T,row.names = F,sep="\t")
}

############################################################
############################################################
# LOOP (genetics, expression, and covariate cols must match)
############################################################

for (chr in seq(1,22)) {

rosmap_geno <- as.data.frame(fread(paste("/Users/dfelsky/Documents/monocyte_twas/eqtl/original_genetic_data/chr",chr,"_cleaned.raw",sep=""),header=T,check.names = F),check.names=F)
rosmap_geno$FID <- ROSmaster$projid[match(rosmap_geno$IID,ROSmaster$IID)]
names(rosmap_geno)[1] <- "id"
rosmap_geno <- rosmap_geno[complete.cases(rosmap_geno),]
sub_list_geno <- rosmap_geno$id

# generate snpinfo file
snpdata <- as.data.frame(fread(paste("/Users/dfelsky/Documents/monocyte_twas/eqtl/original_genetic_data/chr",chr,"_cleaned.bim",sep=""),header=T,check.names = F),check.names=F)

names(snpdata) <- c("chr","snp","cm","pos","a1","a2")
snpdata$snp <- paste(snpdata$snp,snpdata$a1,sep="_")
snpdata$chr <- paste("chr",snpdata$chr,sep="")
snpdata <- snpdata[,c("snp","chr","pos")]
write.table(snpdata,file=paste(base.dir,"/matrix_eqtl_input_files/chr",chr,"_SNPinfo.txt",sep=""),quote=F,col.names = T,row.names = F,sep="\t")

####### get overlapping subject list for dlpfc and mono datasets
overlapdlpfc <- intersect(intersect(sub_list_geno,sub_list_dlpfcexp),sub_list_covfile)
overlapmono <- intersect(intersect(sub_list_geno,sub_list_monoexp),sub_list_covfile)

##### WRITE genotype, expression, and covariate files
# genotype
rosmap_geno1 <- rosmap_geno[which(rosmap_geno$id %in% overlapdlpfc),]
rosmap_geno1 <- rosmap_geno1[match(rosmap_geno1$id,overlapdlpfc),]
rosmap_genot <- t(rosmap_geno1[,-c(2:6)])
colnames(rosmap_genot) <- rosmap_genot[1,]
rosmap_genot <- rosmap_genot[-1,]
write.table(data.frame("id"=rownames(rosmap_genot),rosmap_genot,check.names = F),file=paste(base.dir,"/matrix_eqtl_input_files/chr",chr,"_SNPdata_dlpfc.txt",sep=""),quote=F,col.names = T,row.names = F,sep="\t")

rosmap_geno1 <- rosmap_geno[which(rosmap_geno$id %in% overlapmono),]
rosmap_geno1 <- rosmap_geno1[match(rosmap_geno1$id,overlapmono),]
rosmap_genot <- t(rosmap_geno1[,-c(2:6)])
colnames(rosmap_genot) <- rosmap_genot[1,]
rosmap_genot <- rosmap_genot[-1,]
write.table(data.frame("id"=rownames(rosmap_genot),rosmap_genot,check.names = F),file=paste(base.dir,"/matrix_eqtl_input_files/chr",chr,"_SNPdata_mono.txt",sep=""),quote=F,col.names = T,row.names = F,sep="\t")

# expression
gld <- read.table(paste(base.dir,"/matrix_eqtl_input_files/geneloc_chr",chr,"_dlpfc.txt",sep=""),header=T)
chrgenesd <- gld$geneid[which(gld$chr==paste0("chr",chr))]

norm_dlpfc1 <- norm_dlpfc[,which(colnames(norm_dlpfc) %in% overlapdlpfc)]
norm_dlpfc1 <- norm_dlpfc1[,match(overlapdlpfc,colnames(norm_dlpfc1))]
norm_dlpfc1 <- norm_dlpfc1[which(rownames(norm_dlpfc1) %in% chrgenesd),]
write.table(data.frame("id"=rownames(norm_dlpfc1),norm_dlpfc1,check.names = F),file=paste(base.dir,"/matrix_eqtl_input_files/GE_chr",chr,"_dlpfc.txt",sep=""),quote=F,col.names = T,row.names = F,sep="\t")

gldmo <- read.table(paste(base.dir,"/matrix_eqtl_input_files/geneloc_chr",chr,"_mono.txt",sep=""),header=T)
chrgenesm <- gldmo$geneid[which(gldmo$chr==paste0("chr",chr))]

norm_mono1 <- norm_mono[,which(colnames(norm_mono) %in% overlapmono)]
norm_mono1 <- norm_mono1[,match(overlapmono,colnames(norm_mono1))]
norm_mono1 <- norm_mono1[which(rownames(norm_mono1) %in% chrgenesm),]
write.table(data.frame("id"=rownames(norm_mono1),norm_mono1,check.names = F),file=paste(base.dir,"/matrix_eqtl_input_files/GE_chr",chr,"_mono.txt",sep=""),quote=F,col.names = T,row.names = F,sep="\t")

# covariate
covdlpfc1 <- covfile[which(covfile$id %in% overlapdlpfc),]
covdlpfc1 <- covdlpfc1[match(overlapdlpfc,covdlpfc1$id),]
covfilet <- t(covdlpfc1)
write.table(covfilet,file=paste(base.dir,"/matrix_eqtl_input_files/cov_dlpfc.txt",sep=""),quote=F,col.names = F,row.names = T,sep="\t")

covmono1 <- covfile[which(covfile$id %in% overlapmono),]
covmono1 <- covmono1[match(overlapmono,covmono1$id),]
covfilet <- t(covmono1)
write.table(covfilet,file=paste(base.dir,"/matrix_eqtl_input_files/cov_mono.txt",sep=""),quote=F,col.names = F,row.names = T,sep="\t")

}

###########################
###########################
###########################
###########################
base.dir = "/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/eQTL"
for (tissue in c("dlpfc","mono")) {
for (chr in seq(1,22)) {

# run analyses
#outputfile = paste(base.dir,"/output/chr",chr,"_eqtl_",tissue,"_trans_boot.txt",sep="")
outputfile.cis = paste(base.dir,"/output/chr",chr,"_eqtl_",tissue,"_cis.txt",sep="")
outputfile.cis.boot = paste(base.dir,"/output/chr",chr,"_eqtl_",tissue,"_cis_boot.txt",sep="")
useModel = modelLINEAR
cisDist=1e6

# Genotype file name
SNP_file_name = paste0(base.dir, "/matrix_eqtl_input_files/chr",chr,"_SNPdata_",tissue,".txt")
snps_location_file_name = paste0(base.dir, "/matrix_eqtl_input_files/chr",chr,"_SNPinfo.txt")

# Gene expression file name
expression_file_name = paste0(base.dir, "/matrix_eqtl_input_files/GE_chr",chr,"_",tissue,".txt")
gene_location_file_name = paste0(base.dir, "/matrix_eqtl_input_files/geneloc_chr",chr,"_",tissue,".txt")

# Covariates file name
covariates_file_name = paste0(base.dir, "/matrix_eqtl_input_files/cov_",tissue,".txt")

## Load genotype data
snps = SlicedData$new()
snps$fileDelimiter = "\t"  # TAB delimiter   
snps$fileOmitCharacters = "NA" # denote missing values;
snps$fileSkipRows = 1          # one row of column labels
snps$fileSkipColumns = 1       # one column of row labels
snps$fileSliceSize = 2000      # read file in slices of 2,000 rows
snps$LoadFile(SNP_file_name)
## Load gene expression data
gene = SlicedData$new()
gene$fileDelimiter = "\t"      
gene$fileOmitCharacters = "NA" # denote missing values;
gene$fileSkipRows = 1          # one row of column labels
gene$fileSkipColumns = 1       # one column of row labels
gene$fileSliceSize = 2000      # read file in slices of 2,000 rows
gene$LoadFile(expression_file_name)
## Load covariates
cvrt = SlicedData$new()
cvrt$fileDelimiter = "\t"      
cvrt$fileOmitCharacters = "NA" # denote missing values;
cvrt$fileSkipRows = 1          # one row of column labels
cvrt$fileSkipColumns = 1       # one column of row labels
if(length(covariates_file_name)>0) {
  cvrt$LoadFile(covariates_file_name)
}
snpspos = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
genepos = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)

# BootstrapQTL modification. https://www.ncbi.nlm.nih.gov/pmc/articles/PMC6294523/
# throws error for chr18 because no sig eQTLs (same result with TreeQTL), so skips chr18 in monocytes
if (chr==18 & tissue=="mono") { next }
eQTLs <- BootstrapQTL(snps=snps, 
                      gene=gene, 
                      snpspos=snpspos, 
                      genepos=genepos, 
                      cvrt=cvrt,
                      useModel = useModel,
                      cisDist = cisDist,
                      eGene_detection_file_name=outputfile.cis.boot,
                      n_bootstraps=1000, 
                      n_cores=3)
saveRDS(eQTLs,file=paste0(base.dir,"/output/bootstrapQTL_chr",chr,"_",tissue,".Rdata"))

# vanilla matrix eQTL; run TreeQTL after for proper single and multi-tissue results
#me = Matrix_eQTL_main(
#  snps = snps,
#  gene = gene,
#  cvrt = cvrt,
#  output_file_name = outputfile,
#  output_file_name.cis = outputfile.cis,
#  pvOutputThreshold  = 0,
#  useModel = useModel,
#  errorCovariance = numeric(),
#  verbose = TRUE,
#  pvOutputThreshold.cis = 0.1,
#  snpspos = snpspos,
#  genepos = genepos,
#  cisDist = cisDist,
#  pvalue.hist = TRUE,
#  min.pv.by.genesnp = TRUE,
#  noFDRsaveMemory = FALSE)
#save(me,file=paste0(base.dir,"/output/me_chr",chr,"_",tissue,".Rdata"))
}
}


#### TreeQTL correction (single tissue)
for (tissue in c("mono","dlpfc")) {
  for (chr in seq(1,22)) {
    snps_location_file_name = paste0(base.dir, "/matrix_eqtl_input_files/chr",chr,"_SNPinfo.txt")
    gene_location_file_name = paste0(base.dir, "/matrix_eqtl_input_files/geneloc_chr",chr,"_",tissue,".txt")
    snp_map = read.table(snps_location_file_name, header = TRUE, stringsAsFactors = FALSE)
    gene_map = read.table(gene_location_file_name, header = TRUE, stringsAsFactors = FALSE)
    
    # Get number of genes nearby to each SNP
    n_tests_per_SNP <- get_n_tests_per_SNP(snp_map,gene_map, nearby = TRUE,dist = cisDist)
    eSNPs <- get_eSNPs(n_tests_per_SNP, paste0(base.dir,"/output/chr",chr,"_eqtl_",tissue,"_cis.txt"),level1=0.1)
    saveRDS(eSNPs,file=paste0(base.dir,"/output/eSNPs_chr",chr,"_",tissue,"_cis.rds"))
    # Generate txt output file with full set of eAssociations
    eAssoc_snps <- get_eAssociations(eSNPs, 
                                     n_tests_per_SNP, 
                                     paste0(base.dir,"/output/chr",chr,"_eqtl_",tissue,"_cis.txt"),
                                     paste0(base.dir,"/output/eAssoc_bysnp_chr",chr,"_eqtl_",tissue,"_cis.txt"),
                                     by_snp = TRUE)
    # Identify eGenes regulated by nearby SNPs
    # Get number of genes nearby to each SNP
    n_tests_per_gene <- get_n_tests_per_gene(snp_map,gene_map,nearby = TRUE,dist = cisDist)
    eGenes <- get_eGenes(n_tests_per_gene,paste0(base.dir,"/output/chr",chr,"_eqtl_",tissue,"_cis.txt"))
    saveRDS(eGenes,file=paste0(base.dir,"/output/eGenes_chr",chr,"_",tissue,"_cis.rds"))
    # Generate txt output file with full set of eAssociations
    eAssoc_genes <- get_eAssociations(eGenes, 
                                      n_tests_per_gene, 
                                      paste0(base.dir,"/output/chr",chr,"_eqtl_",tissue,"_cis.txt"),
                                      paste0(base.dir,"/output/eAssoc_bygene_chr",chr,"_eqtl_",tissue,"_cis.txt"),
                                      by_snp = FALSE)
  }
}

########### create SQLite databaset of all cis-eQTL results
library(RSQLite)
library(data.table)
library(naturalsort)

db <- dbConnect(SQLite(), dbname="/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/eQTL/database/dlpfc_mono_eqtls_boot.sqlite")

fdlpfc = naturalsort(grep("eAssoc_bygene.*dlpfc",list.files("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/eQTL/output/"),value=T))
fmono = naturalsort(grep("eAssoc_bygene.*mono",list.files("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/eQTL/output/"),value=T))

for (tissue in c("dlpfc","mono")) {
for (chr in seq(1,22)) {
  if (chr==1 & tissue=="dlpfc") {
    alleqtl <- as.data.frame(fread(paste0("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/eQTL/output/chr",chr,"_eqtl_",tissue,"_cis.txt"),header=T,check.names = F),check.names=F)
    alleqtl$chr <- chr
    alleqtl$tissue <- tissue
    names(alleqtl) <- c("snp","gene","beta","t","p","bbfdr","chr","tissue")
    alleqtl <- alleqtl[,c(8,7,1,2,3,4,5,6)]
  } else {
  toadd <- as.data.frame(fread(paste0("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/eQTL/output/chr",chr,"_eqtl_",tissue,"_cis.txt"),header=T,check.names = F),check.names=F)
  toadd$chr <- chr
  toadd$tissue <- tissue
  names(toadd) <- c("snp","gene","beta","t","p","bbfdr","chr","tissue")
  toadd <- toadd[,c(8,7,1,2,3,4,5,6)]
  alleqtl <- rbind(alleqtl,toadd)
  }
}
}

dbWriteTable(db,"eqtl",alleqtl)
saveRDS(alleqtl,file="/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/eQTL/database/alleqtls_2.rds")





###################################################
#### make all eGenes data frame
base.dir = "/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/eQTL"
for (tissue in c("mono","dlpfc")){
  for (chr in seq(1,22)) {
    if (chr==1) {
      alleGenes<-  readRDS(paste0(base.dir,"/output/eGenes_chr",chr,"_",tissue,"_cis.rds"))
      alleGenes$chr <- chr
      alleGenes$tissue <- tissue
      } else {
        toadd <- readRDS(paste0(base.dir,"/output/eGenes_chr",chr,"_",tissue,"_cis.rds"))
        toadd$chr <- chr
        toadd$tissue <- tissue
        alleGenes <- rbind(alleGenes,toadd)
      }
    }
}


