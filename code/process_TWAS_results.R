### monoctye base analysis
ttmpath <- readRDS("output/TWAS_results_monocyte_pathology.rds")
ttmcog <- readRDS("output/TWAS_results_monocyte_cognition.rds")
ttm <- c(ttmpath,ttmcog)

### monocyte with blood covariates analysis
ttmbpath <- readRDS("output/TWAS_results_monocyte_pathology_blood.rds")
ttmbcog <- readRDS("output/TWAS_results_monocyte_cognition_blood.rds")
ttmb <- c(ttmbpath,ttmbcog)

### monocyte without blood covariates but same reduced sample size (sensitivity)
ttmbspath <- readRDS("output/TWAS_results_monocyte_pathology_noblood_sens.rds")
ttmbscog <- readRDS("output/TWAS_results_monocyte_cognition_noblood_sens.rds")
ttmbs <- c(ttmbspath,ttmbscog)

### DLPFC base analysis
ttdpath <- readRDS("output/TWAS_results_dlpfc_pathology.rds")
ttdcog <- readRDS("output/TWAS_results_dlpfc_cognition.rds")
ttd <- c(ttdpath,ttdcog)

### DLPFC with cell type covariates analysis
ttdcpath_raw <- readRDS("output/TWAS_results_dlpfc_pathology_celltypes.rds")
ttdcpath <- lapply(ttdcpath_raw,function(tt){
  tt <- tt[,c(2:16)]
  tt <- tt[!duplicated(tt[,"gene"]),]
  rownames(tt) <- tt[,"gene"]
  tt
})
ttdccog_raw <- readRDS("output/TWAS_results_dlpfc_cognition_celltypes.rds")
ttdccog <- lapply(ttdccog_raw,function(tt){
  tt <- tt[,c(2:16)]
  tt <- tt[!duplicated(tt[,"gene"]),]
  rownames(tt) <- tt[,"gene"]
  tt
})
ttdc <- c(ttdcpath,ttdccog)

####### all TWAS results
ttall_unlabelled <- list(monocyte=ttm,
              monocyte_blood=ttmb,
              monocyte_sensitivity=ttmbs,
              dlpfc=ttd,
              dlpfc_cells=ttdc)

saveRDS(ttall_unlabelled,"output/all_TWAS_results.rds")

####################################
### generate wide-format version ###
####################################
####################################
index <- 1
for (pheno in names(ttm)) {
  if (index==1) {
    y <- ttm[[pheno]]
    names(y)[c(1:9,14:15)] <- paste0("MONO_",names(y)[c(1:9,14:15)],"_",pheno)
  } else {
    z <- ttm[[pheno]]
    names(z)[c(1:9,14:15)] <- paste0("MONO_",names(z)[c(1:9,14:15)],"_",pheno)
    y <- merge(y,z,by=c("gene","hugo","location","chr"))
  }
  index <- index+1
}
ttmwide <- y

index <- 1
for (pheno in names(ttmb)) {
  if (index==1) {
    y <- ttmb[[pheno]]
    names(y)[c(1:9,14:15)] <- paste0("MONO_BLOOD_",names(y)[c(1:9,14:15)],"_",pheno)
  } else {
    z <- ttmb[[pheno]]
    names(z)[c(1:9,14:15)] <- paste0("MONO_BLOOD_",names(z)[c(1:9,14:15)],"_",pheno)
    y <- merge(y,z,by=c("gene","hugo","location","chr"))
  }
  index <- index+1
}
ttmbwide <- y

index <- 1
for (pheno in names(ttdc)) {
  if (index==1) {
    y <- ttdc[[pheno]]
    names(y)[c(1:9,14:15)] <- paste0("DLPFC_cells_",names(y)[c(1:9,14:15)],"_",pheno)
  } else {
    z <- ttdc[[pheno]]
    names(z)[c(1:9,14:15)] <- paste0("DLPFC_cells_",names(z)[c(1:9,14:15)],"_",pheno)
    y <- merge(y,z,by=c("gene","hugo","location","chr"))
  }
  index <- index+1
}
ttdcwide <- y

ttall_wide <- merge(ttmwide,ttmbwide,by=c("gene","hugo","location","chr"))
ttall_wide <- merge(ttall_wide,ttdcwide,by=c("gene","hugo","location","chr"),all=T)

saveRDS(ttall_wide,"output/all_TWAS_results_wideformat.rds")


#### version of ttall_wide with common variable names between dlpfc and monocyte cog vars
pathnames2 <- c("Neuritic_plaques","Diffuse_plaques","Total_AB","PHF_tau","NFT","Gross_cerebral_infarcts","Micro_Cerebral_infarcts","Arteriolosclerosis","Cerebral_AA","Cerebral_atherosclerosis","Lewy_body_stage","Hippocampal_sclerosis","TDP_43","PD_Dx","Patho_AD","PAM_VM_Caudate","PAM_post_putamen","PAM_IT","PAM_MF")
cognames2 <- c("Global","Episodic_memory","Perceptual_orientation","Perceptual_speed","Semantic_memory","Working_memory","MMSE")

ttm_labelled <- ttm
names(ttm_labelled) <- c(pathnames2,cognames2)
index <- 1
for (pheno in names(ttm_labelled)) {
  if (index==1) {
    y <- ttm_labelled[[pheno]]
    names(y)[c(1:9,14:15)] <- paste0("MONO_",names(y)[c(1:9,14:15)],"_",pheno)
  } else {
    z <- ttm_labelled[[pheno]]
    names(z)[c(1:9,14:15)] <- paste0("MONO_",names(z)[c(1:9,14:15)],"_",pheno)
    y <- merge(y,z,by=c("gene","hugo","location","chr"))
  }
  index <- index+1
}
ttmwidel <- y

ttmb_labelled <- ttmb
names(ttmb_labelled) <- c(pathnames2,cognames2)
index <- 1
for (pheno in names(ttmb_labelled)) {
  if (index==1) {
    y <- ttmb_labelled[[pheno]]
    names(y)[c(1:9,14:15)] <- paste0("MONO_BLOOD_",names(y)[c(1:9,14:15)],"_",pheno)
  } else {
    z <- ttmb_labelled[[pheno]]
    names(z)[c(1:9,14:15)] <- paste0("MONO_BLOOD_",names(z)[c(1:9,14:15)],"_",pheno)
    y <- merge(y,z,by=c("gene","hugo","location","chr"))
  }
  index <- index+1
}
ttmbwidel <- y

ttdc_labelled <- ttdc
names(ttdc_labelled) <- c(pathnames2,cognames2)
index <- 1
for (pheno in names(ttdc_labelled)) {
  if (index==1) {
    y <- ttdc_labelled[[pheno]]
    names(y)[c(1:9,14:15)] <- paste0("DLPFC_cells_",names(y)[c(1:9,14:15)],"_",pheno)
  } else {
    z <- ttdc_labelled[[pheno]]
    names(z)[c(1:9,14:15)] <- paste0("DLPFC_cells_",names(z)[c(1:9,14:15)],"_",pheno)
    y <- merge(y,z,by=c("gene","hugo","location","chr"))
  }
  index <- index+1
}
ttdcwidel <- y

ttall_widel <- merge(ttmwidel,ttmbwidel,by=c("gene","hugo","location","chr"))
ttall_widel <- merge(ttall_widel,ttdcwidel,by=c("gene","hugo","location","chr"),all=T)

saveRDS(ttall_widel,"output/all_TWAS_results_wideformat_labelled.rds")

