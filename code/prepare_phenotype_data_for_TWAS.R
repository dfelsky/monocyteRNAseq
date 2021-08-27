# prepare ROSmaster phenotype data for TWAS
### load libraries
library(rms)
library(limma)
library(edgeR)
library(reshape2)
library(naturalsort)
library(ggrepel)
library(RColorBrewer)
library(MASS)
library(dplyr)
library(stringr)

# load ROSMAP phenotype data and gene symbols reference, and consolidate over old and new updates
source("/Users/dfelsky/Documents/scripts/make_ROSmaster_feb12019.R")
load("input/all_genes_ensembl.RData")
all_genes$mid <- all_genes$start_position + abs((all_genes$end_position - all_genes$start_position)/2)

# format and transform key outcome variables
ROSmaster <- within(ROSmaster,{
  plaq_n_sqrt <- sqrt(plaq_n)
  plaq_d_sqrt <- sqrt(plaq_d)
  amyloid_sqrt <- sqrt(amyloid)
  tangles_sqrt <- sqrt(tangles)
  nft_sqrt <- sqrt(nft)
  
  tdp_stage4 <- as.numeric(tdp_stage4)
  dlbdx <- as.numeric(dlbdx)
  cvda_4gp2 <- as.numeric(cvda_4gp2)
  caa_4gp <- as.numeric(caa_4gp)
  age_death <- as.numeric(age_death)
  pmi <- as.numeric(pmi)
  
  ci_num2_gct <- as.factor(ci_num2_gct)
  ci_num2_mct <- as.factor(ci_num2_mct)
  hspath_any <- as.factor(hspath_any)
  parkdx <- as.factor(parkdx)
  pathoAD <- as.factor(pathoAD)
})


#### get cross-sectional cognition at time of monocyte draw and death separately
mp <- read.csv("input/b123_pheno.txt",sep="\t",header=T,colClasses = c(projid="character"))
mp$projid <- str_pad(mp$projid, width=8, side="left", pad="0")
MP <- subset(mp,monoRNA_batch<3) # necessary because mp includes some repeated measures
drawvisit <- MP[,c("projid","visit")]

roslong <- read.csv("input/dataset_978_long_10-13-2020.csv",colClasses = c(projid="character"))

# add cognitive variables at time of draw
longlist <- lapply(unique(drawvisit$projid), function(x){
  monovis <- drawvisit$visit[which(drawvisit$projid==x)]
  longsub <- subset(roslong, projid==x)
  minindex <- which(longsub$fu_year==monovis)
  longsub[minindex,]
})

newdat <- do.call(rbind,longlist)

newdatsub <- newdat #newdat[,c(1,42:48)] # keep all variables for "_at_draw"
names(newdatsub)[-1] <- paste0(names(newdatsub)[-1],"_at_draw")

ROSmaster <- merge(ROSmaster,newdatsub,by="projid",all.x=T)

# add cognitive variables at last visit
longlist2 <- lapply(unique(roslong$projid), function(x){
  pl <- roslong[which(roslong$projid==x & is.na(roslong$cogn_global)==F),]
  maxindex <- which.max(pl$fu_year)
  pl[maxindex,]
})

newdat2 <- do.call(rbind,longlist2)
newdatsub2 <- newdat2[,c(1,42:48,63)]
names(newdatsub2)[-1] <- paste0(names(newdatsub2)[-1],"_at_lastvisit")

ROSmaster <- merge(ROSmaster,newdatsub2,by="projid",all.x=T)

ROSmaster <- within(ROSmaster,{
fasting.f <- as.factor(ifelse(fasting_at_draw==9,3,fasting_at_draw))
})

saveRDS(ROSmaster,file="input/ROSmaster_TWAS_input.rds")