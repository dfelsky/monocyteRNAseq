# Prep for matrixeQTL (adapted from "http://www.bios.unc.edu/research/genomic_software/Matrix_eQTL/runit.html")
library(MatrixEQTL)
library(stringr)
library(rms)
library(data.table)
library(scales)
library(parallel)
library(BootstrapQTL)

### set directory and load reference data
load("input/all_genes_ensembl.RData")
ROSmaster <- readRDS("input/ROSmaster_TWAS_input.rds")

### read in dge object which has age at draw and also gives us the final set of subjects to use, and merge into ROSmaster
dge_filtered <- readRDS("output/mono_dge_filtered_used.rds")
agedraw <- dge_filtered$samples[,c("projid","age_draw")]
ROSmaster <- merge(ROSmaster,agedraw,by="projid")

### static files - read in once and assess subject overlap in loop
# expression
norm_mono <- readRDS("input/monocytes_techonlynorm_residuals.rds")

sub_list_monoexp <- colnames(norm_mono)
gene_list_monoexp <- rownames(norm_mono)

# retrieve top 10 genomic PCs
ROSMAP_eigen <- read.table("input/eqtl/PCA.eigenvec",header=T,sep="")
ROSmaster <- merge(ROSMAP_eigen,ROSmaster,by="IID")

# covariates 
varsofinterest <- c("msex","age_draw",paste0("PC",seq(1,10)))
covfile <- ROSmaster[,c("projid",varsofinterest)]
covfile <- covfile[which(complete.cases(covfile)),]
names(covfile)[1] <- "id"

## get master list of subjects with geno, exp, and covariate data
master_subs <- intersect(sub_list_monoexp,ROSmaster$projid)

### gene location files per chromosome
for (chr in seq(1,22)) {
geneloc_mono <- all_genes[which(all_genes$ensembl_gene_id %in% gene_list_monoexp),]
geneloc_mono <- geneloc_mono[,c(1,5,3,4)]
names(geneloc_mono) <- c("geneid","chr","left","right")
geneloc_mono <- geneloc_mono[which(geneloc_mono$chr==chr),]
geneloc_mono$chr <- paste("chr",geneloc_mono$chr,sep="")
write.table(geneloc_mono,file=paste("input/eqtl/matrix_eqtl_input_files/geneloc_chr",chr,"_mono.txt",sep=""),quote=F,col.names = T,row.names = F,sep="\t")
}

############################################################
############################################################
# LOOP (genetics, expression, and covariate cols must match)
############################################################
for (chr in seq(1,22)) {
  print(paste0("##### processing chromosome ",chr))
  
  print("## reading genotype data")
  rosmap_geno <- as.data.frame(fread(paste("input/eqtl/original_genetic_data/chr",chr,"_cleaned.raw",sep=""),header=T,check.names = F),check.names=F)
  rosmap_geno$FID <- ROSmaster$projid[match(rosmap_geno$IID,ROSmaster$IID)]
  names(rosmap_geno)[1] <- "id"
  rosmap_geno <- rosmap_geno[complete.cases(rosmap_geno),]
  
  # generate snpinfo file
  snpdata <- as.data.frame(fread(paste("input/eqtl/original_genetic_data/chr",chr,"_cleaned.bim",sep=""),header=F,check.names = F),check.names=F)
  
  names(snpdata) <- c("chr","snp","cm","pos","a1","a2")
  snpdata$snp <- paste(snpdata$snp,snpdata$a2,sep="_")
  snpdata$chr <- paste("chr",snpdata$chr,sep="")
  snpdata <- snpdata[,c("snp","chr","pos")]
  write.table(snpdata,file=paste("input/eqtl/matrix_eqtl_input_files/chr",chr,"_SNPinfo.txt",sep=""),quote=F,col.names = T,row.names = F,sep="\t")
  
  print("## writing genotype data")
  ##### WRITE genotype, expression, and covariate files
  # genotype
  rosmap_geno1 <- rosmap_geno[which(rosmap_geno$id %in% master_subs),]
  rownames(rosmap_geno1) <- rosmap_geno1$id
  rosmap_geno1 <- rosmap_geno1[master_subs,]
  rosmap_geno1 <- rosmap_geno1[,-c(1:6)]
  rosmap_genot <- t(rosmap_geno1)
  write.table(data.frame("id"=rownames(rosmap_genot),rosmap_genot,check.names = F),file=paste("input/eqtl/matrix_eqtl_input_files/chr",chr,"_SNPdata_mono.txt",sep=""),quote=F,col.names = T,row.names = F,sep="\t")
  
  print("## reading expression data")
  # expression
  gldmo <- read.table(paste("input/eqtl/matrix_eqtl_input_files/geneloc_chr",chr,"_mono.txt",sep=""),header=T)
  chrgenesm <- gldmo$geneid[which(gldmo$chr==paste0("chr",chr))]
  norm_mono1 <- norm_mono[which(rownames(norm_mono) %in% chrgenesm),master_subs]
  
  print("## writing expression data")
  write.table(data.frame("id"=rownames(norm_mono1),norm_mono1,check.names = F),file=paste("input/eqtl/matrix_eqtl_input_files/GE_chr",chr,"_mono.txt",sep=""),quote=F,col.names = T,row.names = F,sep="\t")
  
  print("## writing covariate data")
  # covariate
  covmono1 <- covfile
  rownames(covmono1) <- covmono1$id
  covmono1$id <- NULL
  covfilet <- t(covmono1)
  covfilet <- covfilet[,master_subs]
  write.table(data.frame("id"=rownames(covfilet),covfilet,check.names = F),file=paste("input/eqtl/matrix_eqtl_input_files/cov_mono.txt",sep=""),quote=F,col.names = T,row.names = F,sep="\t")
}

###########################
###########################
###########################
###########################
for (tissue in c("mono")) {
  for (chr in seq(1,22)) {
  
  # run analyses
  #outputfile = paste(base.dir,"/output/chr",chr,"_eqtl_",tissue,"_trans_boot.txt",sep="")
  outputfile.cis = paste("output/eqtl/chr",chr,"_eqtl_",tissue,"_cis.txt",sep="")
  outputfile.cis.boot = paste("output/eqtl/chr",chr,"_eqtl_",tissue,"_cis_boot.txt",sep="")
  useModel = modelLINEAR
  cisDist=1e6
  
  # Genotype file name
  SNP_file_name = paste0("input/eqtl/matrix_eqtl_input_files/chr",chr,"_SNPdata_",tissue,".txt")
  snps_location_file_name = paste0("input/eqtl/matrix_eqtl_input_files/chr",chr,"_SNPinfo.txt")
  
  # Gene expression file name
  expression_file_name = paste0("input/eqtl/matrix_eqtl_input_files/GE_chr",chr,"_",tissue,".txt")
  gene_location_file_name = paste0("input/eqtl/matrix_eqtl_input_files/geneloc_chr",chr,"_",tissue,".txt")
  
  # Covariates file name
  covariates_file_name = paste0("input/eqtl/matrix_eqtl_input_files/cov_",tissue,".txt")
  
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
  #if (chr==18 & tissue=="mono") { next }
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
  saveRDS(eQTLs,file=paste0("output/eqtl/bootstrapQTL_chr",chr,"_",tissue,".rds"))
  }
}

###################################################
#### make all eGenes data frame
for (chr in seq(1,22)) {
  if (chr==1) {
    alleGenes <-  readRDS(paste0("output/eqtl/BootstrapQTL_chr",chr,"_mono.rds"))
    alleGenes$chr <- chr
    alleGenes$tissue <- tissue
    } else {
      toadd <- readRDS(paste0("output/eqtl/BootstrapQTL_chr",chr,"_mono.rds"))
      toadd$chr <- chr
      toadd$tissue <- tissue
      alleGenes <- rbind(alleGenes,toadd)
    }
}

saveRDS(alleGenes, "output/eqtl/allegenes.rds")



