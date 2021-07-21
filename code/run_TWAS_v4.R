### Monocyte TWAS for paper
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
library(BRETIGEA)

# load ROSMAP phenotype data and gene symbols reference, and consolidate over old and new updates
source("/Users/dfelsky/Documents/scripts/make_ROSmaster_feb12019.R")
load("/Users/dfelsky/Documents/data/all_genes_ensembl.RData")
all_genes$mid <- all_genes$start_position + abs((all_genes$end_position - all_genes$start_position)/2)

# read in blood cell count phenotype file
mp <- read.csv("/Users/dfelsky/Documents/data/ROSMAP_RNAseq/monocyte/bacth1-3/b123_pheno.txt",sep="\t",header=T,colClasses = c(projid="character"))
mp$projid <- str_pad(mp$projid, width=8, side="left", pad="0")
MP <- subset(mp,monoRNA_batch<3) # necessary because mp includes some repeated measures
drawvisit <- MP[,c("projid","visit")]

roslong <- read.csv("/Users/dfelsky/Documents/data/ROSMAP_phenotype_data/update_10162020/dataset_978_long_10-13-2020.csv",colClasses = c(projid="character"))

longlist <- lapply(unique(drawvisit$projid), function(x){
  monovis <- drawvisit$visit[which(drawvisit$projid==x)]
  longsub <- subset(roslong, projid==x)
  minindex <- which.min(abs(longsub$fu_year - monovis))
  longsub[minindex,]
})

newdat <- do.call(rbind,longlist)


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


##################################
#### set parameters for TWAS  ####
##################################
outdir <- "/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/"
setwd(outdir)

outcome.category.list <- c("pathology","cognition")
tissue.list <- c("monocyte","dlpfc")
celltype.adj <- c("TRUE","FALSE")

techvars.mono <- c("batch","PCT_USABLE_BASES","PERCENT_DUPLICATION","MEDIAN_3PRIME_BIAS","study","PCT_PF_READS_ALIGNED","ESTIMATED_LIBRARY_SIZE")
techvars.dlpfc <- c("batch", "MEDIAN_CV_COVERAGE", "PCT_RIBOSOMAL_BASES", "PCT_CODING_BASES", "PCT_UTR_BASES", "LOG_ESTIMATED_LIBRARY_SIZE", "LOG_PF_READS_ALIGNED", "MEDIAN_5PRIME_TO_3PRIME_BIAS", "PCT_PF_READS_ALIGNED", "study", "PERCENT_DUPLICATION", "MEDIAN_3PRIME_BIAS", "PCT_INTERGENIC_BASES")
  
covars.mono.pathology <- c("msex","age_death","pmi","age_draw")
covars.mono.cognition <- c("msex","age_draw","educ")

covars.dlpfc.pathology <- c("msex","pmi","age_death")
covars.dlpfc.cognition <- c("educ","msex","age_death")

indepvec.pathology <- c("plaq_n_sqrt","plaq_d_sqrt","amyloid_sqrt","tangles_sqrt","nft_sqrt","ci_num2_gct","ci_num2_mct","arteriol_scler","caa_4gp","cvda_4gp2","dlbdx","hspath_any","tdp_stage4","parkdx","pathoAD","vm3123","pput3123","it3123","mf3123") 

indepvec.cognition <- c("cogn_global_random_slope","cogn_ep_random_slope","cogn_po_random_slope","cogn_ps_random_slope","cogn_se_random_slope","cogn_wo_random_slope")

maxnum.toplot <- 100

##################################################################################
############################# start of TWAS loop #################################
##################################################################################
##################################################################################
for (tissue in tissue.list) {
  print(paste("##### starting tissue:",tissue))
  
  for (outcome.category in outcome.category.list) {
    print(paste("### starting outcome caterogy:", outcome.category))
    
    index <- 1
    ttlist <- list()
    ttlistdesign <- list()
  
    # read in filtered data from normalize_mono_data.R or other script
    # specify technical variables, other covariates, and outcome variables
    if (tissue =="dlpfc") {
      dge_filtered <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/normalization/dlpfc_filtered_only.rds")
      techvars <- techvars.dlpfc
      
      ### remove outlier subjects
      removesubs.dlpfc <- c("20634274","23690880","50103967","32697960","20886846","43596435","16322424","30544882","74522154","21246218","74284255","11390174","69924281","85578107","11475462","00402800","01797756","24141372")
      dge_filtered <- dge_filtered[,which(dge_filtered$samples$projid %nin% removesubs.dlpfc)]
      saveRDS(dge_filtered,file="dlpfc_dge_filtered_used.rds")
      
      if(outcome.category=="pathology"){
      indepvec <- indepvec.pathology
      covars <- covars.dlpfc.pathology
      } else if (outcome.category=="cognition") {
        indepvec <- indepvec.cognition
        covars <- covars.dlpfc.cognition
      } } else if (tissue=="monocyte") {
        dge_filtered <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/normalization/monocytes_filtered_only.rds")
        techvars <- techvars.mono
        
        ### remove outlier subjects
        removesubs.mono <- c("00246264","35072859","91804757","21135680","69866926","50303145")
        dge_filtered <- dge_filtered[,which(dge_filtered$samples$projid %nin% removesubs.mono)]
        saveRDS(dge_filtered,file="mono_dge_filtered_used.rds")
        
        if(outcome.category=="pathology"){
          indepvec <- indepvec.pathology
          covars <- covars.mono.pathology
        } else if (outcome.category=="cognition") {
          indepvec <- indepvec.cognition
          covars <- covars.mono.cognition
        } }
  
      # set parameters for plotting individual effects
      sigthres <- 0.05
      ## determine cont. index automatically
      if (length(indepvec)==1) { 
        tablist <- length(table(ROSmaster[,indepvec]))
      } else {
        tablist <- apply(ROSmaster[,indepvec],2,table) %>%
          lapply(length) %>%
          as.numeric()
      }
      cont.index <- tablist > 2 #TRUE if numeric indepvec
  
  for (indep in indepvec) {
    print(paste("## starting phenotype:",indep))
    
    dge_fortwas <- dge_filtered # preserve full tissue-specific QC'd DGEList object for plotting later
    
    designvars <- c(techvars,covars,indep) 
    
    # special line to accommodate microglial phenotype analysis
    # removes study variable, since it is only measured in MAP subjects
    if (indep %in% c("vm3123","pput3123","it3123","mf3123")){
      designvars <- designvars[-which(designvars=="study")]
    }
    
    toharvest <- designvars[which(designvars %nin% names(dge_fortwas$samples))]
    
    dge_fortwas$samples[,toharvest] <- ROSmaster[match(dge_fortwas$samples$projid,ROSmaster$projid),
                                                 match(toharvest,names(ROSmaster))]

    
    dge_fortwas <- dge_fortwas[,which(complete.cases(dge_fortwas$samples)==T)]
    
    # remove levels of batch that have been brought to zero
    dge_fortwas$samples$batch <- as.factor(as.character(dge_fortwas$samples$batch))
    
    if (cont.index[index]==T) {
      dge_fortwas$samples[,indep] <- as.numeric(dge_fortwas$samples[,indep]) 
      phenvalues <- mean(dge_fortwas$samples[,indep])
      phenprops <- sd(dge_fortwas$samples[,indep]) } else {
        dge_fortwas$samples[,indep] <- as.factor(dge_fortwas$samples[,indep])
        phenostats <- describe(dge_fortwas$samples[,indep])$values
        phenvalues <- paste(phenostats[[1]],collapse = "_")
        phenprops <- paste(phenostats[[2]],collapse = "_")
      }

    designtext <- paste("model.matrix( ~ ",paste(designvars,sep="",collapse=" + "),", data=dge_fortwas$samples)")
    design <- eval(parse(text=designtext))
    
    # calculate normalization factors, run voom
    print("## voom")
    dge_fortwas <- calcNormFactors(dge_fortwas)
    
    if (tissue=="monocyte") {
      v <- voomWithQualityWeights(dge_fortwas,design)
    } else {
      v <- voom(dge_fortwas,design)
    }
    
    n <- dim(v)[2] # track input sample size
    
    # run TWAS robust linear model
    print("## lmFit")
    lmod <- lmFit(v, design, method = "robust",maxit=20000)
    eb <- eBayes(fit=lmod,robust = T)
    
    coefnum <- grep(indep,colnames(eb))
    
    if (length(coefnum) > 1) {
    tt <- topTable(eb,coef=coefnum,20000,sort.by = "F")  
    } else {
    tt <- topTable(eb,coef=coefnum,20000,sort.by = "p")
    }
  
    tt$se <- eb$stdev.unscaled[,coefnum]*sqrt(eb$s2.post)
    tt$n <- n
    tt$pheno <- indep
    tt$gene <- rownames(tt)
    tt$hugo <- all_genes$external_gene_id[match(tt$gene,all_genes$ensembl_gene_id)]
    tt$location <- all_genes$mid[match(tt$gene,all_genes$ensembl_gene_id)]
    tt$chr <- all_genes$chromosome_name[match(tt$gene,all_genes$ensembl_gene_id)]
    tt$phenvalues <- phenvalues
    tt$phenprops <- phenprops
    
    ttlist[[indep]] <- tt
    ttlistdesign[[indep]] <- lmod$design

    print("## complete")
    index <- index+1
  }
  
      print(paste("### saving summary statistics for:",outcome.category,"/",tissue))
  saveRDS(ttlistdesign,file=paste0(outdir,"TWAS_designmatrices_",tissue,"_",outcome.category,".rds"))
  saveRDS(ttlist,file=paste0(outdir,"TWAS_results_",tissue,"_",outcome.category,".rds"))
  
  #################################################
  ######### START OF MANHATTAN PLOTTING ###########
  #################################################
  
  #### visualize results and save plots
  ttlist2 <- lapply(ttlist, FUN = function(x) {
    z <- subset(x, is.na(chr)==F)
    z$chr <- ifelse(z$chr=="X",23,z$chr)
    z$chr <- ifelse(z$chr=="MT",24,z$chr)
    z$transpos <- z$location + as.numeric(z$chr)*1e12
    z$transpos <- rank(z$transpos)
    z
  })
  
  index=1
  plotlist <- list()
  
  numtolabel <- 20
  pthres <- 1
  
  print(paste("### generating manhattan plots for:",outcome.category,"/",tissue))
  for (tt in seq(1,length(ttlist2)) ) {
    
    ttlist2[[tt]]$chr <- factor(ttlist2[[tt]]$chr, levels=naturalsort(unique(ttlist2[[tt]]$chr)))
    startloc <- 1
    xmins <- NULL
    xmaxs <- NULL
    
    for (chr in c(1:24)) {
      xmins[chr] <- startloc
      xmaxs[chr] <- startloc + table(ttlist2[[tt]]$chr)[chr]
      startloc <- xmaxs[chr]
    }
    
    xmids <- xmins + (xmaxs-xmins)/2
    ylabpos <- log10(min(c(ttlist2[[tt]]$adj.P.Val)))-1
    
    outcomevar <- names(ttlist2)[[tt]]
    modelvars <- colnames(ttlistdesign[[tt]])[-which(colnames(ttlistdesign[[tt]]) %in% c("(Intercept)",outcomevar))]
    
    if (outcomevar %in% c("vm3123","pput3123","it3123","mf3123")) {
      modelvars <- modelvars[-which(modelvars=="study")]
    }
    
    plotlist[[index]] <- ggplot(data=ttlist2[[tt]],
                                aes(y=-log(adj.P.Val,base=10)*sign(logFC),x=transpos))+
      geom_point()+
      annotate("rect",
               xmin=xmins, 
               xmax=xmaxs, 
               ymin=-Inf, 
               ymax=Inf,
               fill=rep(c("grey","white"),12),
               alpha=0.2)+
      annotate("text",x=xmids, y=ylabpos-0.2,label=c(1:22,"X","MT"),size=3)+
      theme_classic()+
      geom_hline(yintercept = -log(0.05,base=10),col="orange",lty=4)+
      geom_hline(yintercept = log(0.05,base=10),col="orange",lty=4)+
      geom_hline(yintercept =0,col="blue")+
      geom_text_repel(data=ttlist2[[tt]][1:numtolabel,],
                      aes(label=hugo))+
      theme(axis.line=element_blank(),
            axis.text.x=element_blank(),
            axis.ticks.x = element_blank())+
      labs(title=paste("Outcome = ",outcomevar,"; tissue = ",tissue," (n=",ttlist2[[tt]]$n[1],")",sep=""),
           subtitle=paste("covariates = ",paste(modelvars,collapse=", "),sep=""),
           y="-log10(FDR p-value)",
           x="Chromosome")
    
    index=index+1
  }
  
  ## save plots
  require(gridExtra)
  
  # manhattan
  print(paste("### saving manhattan plots for:",outcome.category,"/",tissue))
  pdf(paste0(outdir,"Manhattan_plot_",tissue,"_",outcome.category,".pdf"),width=14,height=9, onefile = T)
  for (a in seq(1,length(plotlist))) { print(plotlist[[a]]) }
  dev.off()
  
  #######################################################
  ########## START OF TOP RESULT SCATTERPLOTS ###########
  #######################################################
  
  print(paste("### generating individual significant result plots for:",outcome.category,"/",tissue))
  ##### get tables of significant results:
  allres <- do.call(rbind,lapply(ttlist2,function(x) { 
    y <- x[which(x[,"adj.P.Val"] < sigthres ),]
    ylength <- nrow(y)
    if (ylength > maxnum.toplot) {
      y[1:maxnum.toplot,]
    } else {
      y
    }
    }))
  
  if (dim(allres)[1]==0) { next }
  
  allres$trait <- gsub(rownames(allres),pattern="\\..*",replacement="")
  
  datplotlist <- list()
  for (trait in unique(allres$trait)) {
    print(paste("## getting normalized data for:",trait))
    
    # generate partial residuals for plotting
    covs <- ttlistdesign[[trait]][,-c(1,grep(trait,colnames(ttlistdesign[[trait]])))]
    design2 <- as.matrix(ttlistdesign[[trait]][,grep(trait,colnames(ttlistdesign[[trait]]))])
    design2 <- cbind(rep(1,length(design2)),design2)
    colnames(design2) <- c("(intercept)",trait)
    
    dge_forplot <- dge_filtered[,rownames(design2)]
    
    v2 <- voom(dge_forplot)
    normdata <- removeBatchEffect(v2$E, covariates=covs, design=design2)
    
    norm <- as.data.frame(t(normdata))
    norm$projid <- rownames(norm)
    m.norm <- merge(norm,ROSmaster,by="projid")
    
    genes <- allres$gene[which(allres$trait==trait & allres$gene %in% names(m.norm))]
    RM2 <- m.norm[,c(trait,genes)]
    
    RMmelt <- melt(RM2,id.vars = trait)
    RMmelt$hugo <- all_genes$external_gene_id[match(RMmelt$variable,all_genes$ensembl_gene_id)]
    
    if (cont.index[which(indepvec==trait)]==TRUE) {
      datplotlist[[trait]] <- ggplot(data=RMmelt, aes_string(y="value", x=trait))+
        geom_point()+
        geom_smooth(method="lm")+
        facet_wrap(~hugo,scales="free_y")+
        theme_bw()+
        labs(y="Expression (log2(CPM))")
      
    } else if (cont.index[which(indepvec==trait)]==FALSE) {
      datplotlist[[trait]] <- ggplot(data=RMmelt, aes_string(y="value", x=trait, group=trait))+
        geom_boxplot(outlier.colour=NA)+
        geom_jitter(col="darkgrey",width=0.1,alpha=0.75,pch=19)+
        facet_wrap(~hugo, scales = "free_y")+
        theme_bw()+
        labs(y="Expression (log2(CPM))")+
        theme(legend.position="none")
    }
  }
  
  print(paste("### saving individual significant result plots for:",outcome.category,"/",tissue))
  pdf(paste0(outdir,"Individual_effects_",tissue,"_",outcome.category,".pdf"),width=16,height=16, onefile = T)
  for (a in seq(1,length(datplotlist))) { print(datplotlist[[a]]) }
  dev.off()
  }
}



###################################################
################ FULL LOOP END ####################
###################################################
library(corrplot)
library(rms)
library(BRETIGEA)

## pathology
ttd <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/TWAS_results_dlpfc_pathology.rds")
ttdc_raw <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/TWAS_results_dlpfc_pathology_celltypes.rds")
ttdc <- lapply(ttdc_raw,function(tt){
  tt <- tt[,c(2:16)]
  tt <- tt[!duplicated(tt[,"gene"]),]
  rownames(tt) <- tt[,"gene"]
  tt
  })
ttm <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/TWAS_results_monocyte_pathology.rds")

ttdp <- ttd
ttdcp <- ttdc
ttmp <- ttm

nameset1 <- c("Neuritic plaques","Diffuse plaques","Total AB","PHF tau","NFT","Gross cerebral infarcts","Micro Cerebral infarcts","Arteriolosclerosis","Cerebral AA","Cerebral atherosclerosis","Lewy body stage","Hippocampal sclerosis","TDP-43","PD Dx","Patho AD","PAM VM Caudate","PAM post. putamen","PAM IT","PAM MF")


## cognition
ttd <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/newcog/TWAS_results_dlpfc_cognition.rds")
ttdc_raw <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/newcog/TWAS_results_dlpfc_cognition_celltypes.rds")
ttdc <- lapply(ttdc_raw,function(tt){
  tt <- tt[,c(2:16)]
  tt <- tt[!duplicated(tt[,"gene"]),]
  rownames(tt) <- tt[,"gene"]
  tt
})
ttm <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/newcog/TWAS_results_monocyte_cognition.rds")

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


siglistm <- lapply(ttm,function(x) { x$hugo[which(x$adj.P.Val < 0.05)] })

intersect(intersect(siglistm$plaq_n_sqrt,siglistm$plaq_d_sqrt),siglistm$amyloid_sqrt)
stest1 <- supertest(siglistm,degree=c(2,3))
plot.msets(stest1,Layout = "landscape",
           keep.empty.intersections = F,
           show.elements = T,
           margin = c(1,10,10,5),
           color.on = "black",
           sort.by = "size")


### overlapping genes for PAM phenotypes
allpam <- c(siglistm$mf3123,siglistm$it3123,siglistm$pput3123,siglistm$vm3123)
length(unique(allpam))

lapply(ttm[,])


### overlapping genes for amyloid phenotype
allamy <- c(siglistm$plaq_n_sqrt,siglistm$plaq_d_sqrt,siglistm$amyloid_sqrt)
length(unique(allamy))


ggplot(data=subset(sigres,tissue!="dlpfc"),aes(y=nsig,x=pheno,fill=sig))+
  geom_bar(stat = "identity",position = "identity",alpha=0.4)+
  facet_wrap(~tissue,nrow=1,scales = "free")+
  geom_text(aes(label=nsig),hjust=-0.1,size=2.5)+
  coord_flip()+
  scale_fill_aaas()+
  theme_classic()


########## cognition
siglistmc <- lapply(ttmc,function(x) { x$hugo[which(x$adj.P.Val < 0.05)] })
stest2 <- supertest(siglistmc,degree=c(1,2,3,4,5))
plot.msets(stest2,Layout = "landscape",
           keep.empty.intersections = F,
           show.elements = T,
           margin = c(1,10,10,5),
           color.on = "black",
           sort.by = "size")

siglistmc2 <- lapply(ttmc,function(x) { x$hugo[which(x$adj.P.Val < 0.1)] })
stest22 <- supertest(siglistmc2,degree=c(1,2,3,4,5))
plot.msets(stest22,Layout = "landscape",
           keep.empty.intersections = F,
           show.elements = T,
           margin = c(1,10,10,5),
           color.on = "black",
           sort.by = "size")


siglistmc2 <- lapply(ttmc,function(x) { x$hugo[which(x$adj.P.Val < 0.1)] })
stest22 <- supertest(siglistmc2,degree=c(1,2,3,4,5))
plot.msets(stest22,Layout = "landscape",
           keep.empty.intersections = F,
           show.elements = T,
           margin = c(1,10,10,5),
           color.on = "black",
           sort.by = "size")


#### plot volcano plots for PAM phenotypes and diffuse plaques and atherosclerosis

p1 <- ggplot(ttm$mf3123,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title = "A. PAM (MF)")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=ttm$mf3123[1:5,],aes(label=hugo))+
  theme_minimal()
p2 <- ggplot(ttm$it3123,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title = "B. PAM (IT)")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=ttm$it3123[1:5,],aes(label=hugo))+
  theme_minimal()
p3 <- ggplot(ttm$vm3123,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title="C. PAM (VM)")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=ttm$vm3123[1:5,],aes(label=hugo))+
  theme_minimal()
p4 <- ggplot(ttm$pput3123,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title="D. PAM (PPUT)")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=ttm$pput3123[1:5,],aes(label=hugo))+
  theme_minimal()
p5 <- ggplot(ttm$plaq_d_sqrt,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title="E. Diffuse plaques")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=ttm$plaq_d_sqrt[1:5,],aes(label=hugo))+
  theme_minimal()
p6 <- ggplot(ttm$cvda_4gp2,aes(y=-log10(adj.P.Val),x=logFC))+
  geom_point(aes(col=as.factor(ifelse(adj.P.Val<0.05,1,0))),show.legend = F)+
  scale_color_aaas()+
  labs(y="-log10(FDR p-value)",title="F. Cereb. atheroscl.")+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text_repel(data=ttm$cvda_4gp2[1:5,],aes(label=hugo))+
  theme_minimal()

plot_grid(p1,p2,p3,p4,p5,p6,ncol=3)


###############################################################
#                                                             #
#                 enrichment LOOP begins                      #
#                                                             #
###############################################################
library(cowplot)
library(readr)
library(dplyr)
library(ggplot2)
library(tidyr)
library(magrittr)
library(AnnotationDbi)
library(annotate)
library(GO.db)
library(tmod)

source("/Users/dfelsky/Documents/GT_RNAseq/NEW_analyses_Jan21/three_sample_data/enrichment/make_GOset_file.R")

make.geneGO.list(gene.ids=ttm$plaq_n_sqrt$hugo,
                 out.file = "/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/GOset_file_monocyte.rds",maxGOgroupSize = 500, minGOgroupSize = 25)
make.geneGO.list(gene.ids=ttd$plaq_n_sqrt$hugo,
                 out.file = "/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/GOset_file_dlpfc.rds",maxGOgroupSize = 500, minGOgroupSize = 25)

GOfiled <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/GOset_file_dlpfc.rds")
GOfilem <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/GOset_file_monocyte.rds")

GOl <- list(mono=GOfilem,dlpfc=GOfiled,dlpfc_cells=GOfiled)
TTl <- list(mono=ttm,dlpfc=ttd,dlpfc_cells=ttdc)

setwd(outdir)
minsiggenes <- 10

for (tissue in c("mono","dlpfc","dlpfc_cells")) {
  
  # retrieve summary stats and GO file
  ttlist <- TTl[[tissue]]
  Gofile <- GOl[[tissue]]
  
  for (qthreshold in c(0.05,0.1)) {
  # get significant genes for each trait
  labeledgenes <- lapply(ttlist, function(x) {
    hugo <- x$hugo[which(x$adj.P.Val < qthreshold)]
    pheno <- x$pheno[which(x$adj.P.Val < qthreshold)]
    return(cbind(hugo,pheno))
  })
  
  siggroups <- as.data.frame(do.call(rbind,labeledgenes))
  
  # remove sets with fewer than minsiggenes
  groupdefs <- subset(siggroups, pheno %nin% names(which(table(siggroups$pheno) < minsiggenes)))
  
  enrichres <- lapply(unique(groupdefs$pheno), function(pheno) {
    hitListGenes <- as.character(groupdefs$hugo[which(groupdefs$pheno==pheno)])
    sortedGenes <- groupdefs$hugo
    result <- tbl_df(tmodHGtest(fg = hitListGenes, 
                                bg = sortedGenes, 
                                mset=GOl[[tissue]],
                                qval = 1.01, 
                                filter = F))
    
    result %<>% rowwise() %>% mutate(aspect = Ontology(ID))
    result$rank <- 1:nrow(result)
    result %<>% dplyr::select(Title, 
                              overlap = b, 
                              setSize=B, 
                              hitListSize = n,
                              genesWithGOAnnotation = N, 
                              P.Value, 
                              adj.P.Val, 
                              everything()) %>% arrange(adj.P.Val)
    return(result)
  })
  
  names(enrichres) <- unique(groupdefs$pheno)
  enrichresdf <- lapply(enrichres,as.data.frame)
  
  allresenrich <- do.call(rbind,enrichresdf)
  allresenrich$pheno <- gsub(rownames(allresenrich),pattern="[[:punct:]].*",replacement="")
  
  saveRDS(allresenrich,file=paste0(outdir,"Enrichment_results_",tissue,"_pathology_FDR",qthreshold,".rds"))
  }
  
  
  ###### AUC based enrichment, sorted positive to negative t-value
  allauclist_tval <- lapply(ttlist, function(pheno,GOfile=Gofile) {
    orderedgenes <- pheno$hugo[order(pheno$t,decreasing = T)]
    res <- as.data.frame(tmodUtest(orderedgenes,
            mset=GOfile,
            qval=1.01,
            filter=F))
    res$pheno <- pheno$pheno[1]
    res
  })
  
  ###### AUC based enrichment, sorted only by p-value
  allauclist_pval <- lapply(ttlist, function(pheno,GOfile=Gofile) {
    orderedgenes <- pheno$hugo
    res <- as.data.frame(tmodUtest(orderedgenes,
                                   mset=GOfile,
                                   qval=1.01,
                                   filter=F))
    res$pheno <- pheno$pheno[1]
    res
  })
  
  allaucres_pval <- do.call(rbind,allauclist_pval)
  rownames(allaucres) <- NULL
  
  allaucres_tval <- do.call(rbind,allauclist_tval)
  rownames(allaucres) <- NULL
  
  saveRDS(allaucres_pval,file=paste0(outdir,"AUCEnrichment_results_",tissue,"_pathology_Psorted.rds"))
  saveRDS(allaucres_tval,file=paste0(outdir,"AUCEnrichment_results_",tissue,"_pathology_Tsorted.rds"))
}


####################################################################
###### plot number of sig effects in both tissues per pathology
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
for (pheno in names(ttd)) {
  if (index==1) {
    y <- ttd[[pheno]]
    names(y)[c(1:9,14:15)] <- paste0("DLPFC_",names(y)[c(1:9,14:15)],"_",pheno)
  } else {
    z <- ttd[[pheno]]
    names(z)[c(1:9,14:15)] <- paste0("DLPFC_",names(z)[c(1:9,14:15)],"_",pheno)
    y <- merge(y,z,by=c("gene","hugo","location","chr"))
  }
  index <- index+1
}
ttdwide <- y

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

################## actual between-tissue gene expression correlations
mono <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/normalization/monocytes_techandagesex_residuals.rds")
removesubs.mono <- c("00246264","35072859","91804757","21135680","69866926","50303145")
mono <- mono[,which(colnames(mono) %nin% removesubs.mono)]
dlpfc <- readRDS("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/normalization/dlpfc_techandagesex_residuals.rds")
removesubs.dlpfc <- c("20634274","23690880","50103967","32697960","20886846","43596435","16322424","30544882","74522154","21246218","74284255","11390174","69924281","85578107","11475462","00402800","01797756","24141372")
dlpfc <- dlpfc[,which(colnames(dlpfc) %nin% removesubs.dlpfc)]
dlpfc.nocells <- dlpfc

### correct for cell types
genesymbols <- all_genes$external_gene_name[match(rownames(dlpfc),all_genes$ensembl_gene_id)]
rownames(dlpfc) <- genesymbols
dlpfc <- dlpfc[!is.na(genesymbols),]
dlpfc.adj <- adjustBrainCells(dlpfc,species = "human")
dlpfc <- dlpfc.adj$expression
rownames(dlpfc) <- all_genes$ensembl_gene_id[match(rownames(dlpfc),all_genes$external_gene_name)]
dlpfc.celltypes <- dlpfc

# commonlists no cell type correction
commonsubs <- intersect(colnames(mono),colnames(dlpfc.nocells))
commongenes <- intersect(rownames(mono),rownames(dlpfc.nocells))
# commonlists with cell type correction
commonsubs.cells <- intersect(colnames(mono),colnames(dlpfc.celltypes))
commongenes.cells <- intersect(rownames(mono),rownames(dlpfc.celltypes))

# new matched dfs no cells
m2 <- t(mono[commongenes,commonsubs])
d2 <- t(dlpfc.nocells[commongenes,commonsubs])
# new matched dfs with cells
m2.cells <- t(mono[commongenes.cells,commonsubs.cells])
d2.cells <- t(dlpfc.celltypes[commongenes.cells,commonsubs.cells])

# correlations no cells
allcors <- lapply(commongenes, FUN=function(x) { cor.test(m2[,x],d2[,x],method = "spearman") })
allcor_r <- as.numeric(unlist(lapply(allcors,FUN=function(x) { x$estimate })))
allcor_p <- as.numeric(unlist(lapply(allcors,FUN=function(x) { x$p.value })))
# correlations with cells
allcors.cells <- lapply(commongenes.cells, FUN=function(x) { cor.test(m2.cells[,x],d2.cells[,x],method = "spearman") })
allcor_r.cells <- as.numeric(unlist(lapply(allcors.cells,FUN=function(x) { x$estimate })))
allcor_p.cells <- as.numeric(unlist(lapply(allcors.cells,FUN=function(x) { x$p.value })))

# allinfo no cells
allinfo <- data.frame(gene=commongenes,
                      cor_r=allcor_r,
                      cor_p=allcor_p)
allinfo$cor_fdr <- p.adjust(allinfo$cor_p,method="BH")
allinfo$hugo <- all_genes$external_gene_id[match(allinfo$gene,all_genes$ensembl_gene_id)]
allinfo$cor_sig <- ifelse(allinfo$cor_p<0.05,1,0)
allinfo$cor_sig <- as.factor(ifelse(allinfo$cor_fdr<0.05,2,allinfo$cor_sig))
# allinfo with cells
allinfo.cells <- data.frame(gene=commongenes.cells,
                            cor_r_celltype=allcor_r.cells,
                            cor_p_celltype=allcor_p.cells)
allinfo.cells$cor_fdr_celltype <- p.adjust(allinfo.cells$cor_p,method="BH")
allinfo.cells$hugo <- all_genes$external_gene_id[match(allinfo.cells$gene,all_genes$ensembl_gene_id)]
allinfo.cells$cor_sig_celltype <- ifelse(allinfo.cells$cor_p_celltype<0.05,1,0)
allinfo.cells$cor_sig_celltype <- as.factor(ifelse(allinfo.cells$cor_fdr_celltype<0.05,2,allinfo.cells$cor_sig_celltype))

#merge
allinfo.both <- merge(allinfo,allinfo.cells,by=c("gene","hugo"))

library(ggsci)
library(ggrepel)
library(ggpubr)
library(gridExtra)

# no cell type correction
ggplot(data=allinfo.both, aes(x=cor_r,fill=cor_sig))+
  geom_histogram(bins=500)+
  scale_fill_manual(values = colors()[c(12,30,54)])+
  geom_vline(xintercept = c(0,median(allinfo$cor_r)),col=c("black","yellow"),lty=c(2,4))+
  annotate(x=c(-0.25,-0.19,0.01,0.2,0.34),y=c(13,25,50,25,13),label=as.numeric(table(as.numeric(as.character(allinfo$cor_sig))*sign(allinfo$cor_r))),geom = "text",col=colors()[c(54,30,24,30,54)])+
  theme_minimal()

# with cell type correction
ggplot(data=allinfo.both, aes(x=cor_r_celltype,fill=cor_sig_celltype))+
  geom_histogram(bins=500)+
  scale_fill_manual(values = colors()[c(12,30,54)])+
  geom_vline(xintercept = c(0,median(allinfo.both$cor_r_celltype)),col=c("black","yellow"),lty=c(2,4))+
  annotate(x=c(-0.25,-0.19,0.01,0.2,0.34),y=c(13,25,50,25,13),label=as.numeric(table(as.numeric(as.character(allinfo.both$cor_sig_celltype))*sign(allinfo.both$cor_r_celltype))),geom = "text",col=colors()[c(54,30,24,30,54)])+
  theme_minimal()


##### merge info - THIS WILL CUT OUT GENES NOT TRANSLATED FROM ENSMBL IDs to HGNC SYMBOLS
allwide1 <- merge(ttmwide,ttdwide,by=c("gene","hugo","location","chr"))
allwide2 <- merge(allwide1,ttdcwide,by=c("gene","hugo","location","chr"))
allinfo2 <- merge(allinfo.both,allwide2,by=c("gene","hugo"))
saveRDS(allinfo2,file="/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/newcog/Merged_summary_and_correlation_data_spearman_cognition.rds")

####### cross-sample overlap per phenotype
library(SuperExactTest)
#ttdorig <- ttd
#ttd <- ttdc # run this swap to test the overlap with corrected effects - there is more overlap

commongenes <- ad$gene

index <- 1
stestlist <- list()
for (sigthres in c(0.05,0.1,0.2)) {
  for (pheno in names(ttm)) { ##### this needs to be adjusted for cognitive outcomes because monocyte and dlpfc have different outcome names (i.e. at_last_visit vs. at_draw)
    ttmp <- ttm[[pheno]][commongenes,]
    ttdp <- ttdc[[pheno]][commongenes,]
    
    monosig <- ttmp$gene[which(ttmp$adj.P.Val < sigthres)]
    dlpfcsig <- ttdp$gene[which(ttdp$adj.P.Val < sigthres)]
  
    if(length(monosig)==0 & length(dlpfcsig)==0) { 
      stestlist[[index]] <- data.frame(pheno=pheno,
                                       FDR=as.character(sigthres),
                                       dlpfcsig=0,
                                       monosig=0,
                                       overlap=NA,
                                       expected.overlap=NA,
                                       p=NA)
      } else {
      
    sig_sets <- list(mono=monosig,
                   dlpfc=dlpfcsig)
  
    sumobj <- summary(supertest(sig_sets,n=length(commongenes)))
    
    stestlist[[index]] <- data.frame(pheno=pheno,
                            FDR=as.character(sigthres),
                            dlpfcsig=as.numeric(sumobj$set.sizes["dlpfc"]),
                            monosig=as.numeric(sumobj$set.sizes["mono"]),
                            overlap=as.numeric(sumobj$otab[3]),
                            expected.overlap=as.numeric(sumobj$etab[3]),
                            p=as.numeric(sumobj$P.value[3]))
      }
    index <- index + 1
  }
}

setdata <- as.data.frame(do.call(rbind,stestlist))
setdata

setdata$pheno <- factor(setdata$pheno,levels=names(ttm)[c(1:5,15,9,6:8,10,12:13,11,14,16:19)])
ggplot(data=setdata,aes(y=-log10(p),x=pheno,col=FDR))+
  geom_point(aes(size=overlap))+
  scale_color_aaas()+
  coord_flip()+
  geom_hline(yintercept = -log10(0.05),lty=2,col="red")+
  geom_text(aes(label=overlap),col="black",nudge_y = 0.05,nudge_x = 0.6,size=2)+
  theme_hc()
  #theme(axis.text.x=element_text(angle = -90, hjust = 0))


# SIGNIFICANT GENE OVERLAP IS NOT SIGNIFICANT FOR ANY PHENOTYPE AT MULTIPLE THRESHOLDS (except TDP43 nominal)

####### cross-sample correlations: DETERMINED THAT logFC is not useful as a metric because some genes will have larger magnitude effects when not standardized by their standard deviation/error. ALSO, cortype does not affect results, so spearman is best
allcordats <- list()
index <- 1

for (pheno in names(ttm)) {
  for (abso in c(TRUE,FALSE)) {
    mvar <- paste0("MONO_t_",pheno)
    dvar <- paste0("DLPFC_cells_t_",pheno)
    
    if (abso==TRUE) {
      mdcor <- cor.test(abs(allinfo2[,mvar]),abs(allinfo2[,dvar]), method="spearman")
    } else {
      mdcor <- cor.test(allinfo2[,mvar],allinfo2[,dvar],method="spearman")
    }

    cplotdat <- c(mdcor$estimate,p=mdcor$p.value,pheno=pheno,abs=abso)
    allcordats[[index]] <- cplotdat
    
    index <- index + 1
  }
}



aplotdat <- as.data.frame(do.call(rbind,allcordats))
aplotdat <- within(aplotdat,{
  rho <- as.numeric(as.character(rho))
  p <- as.numeric(as.character(p))
  pheno <- as.factor(pheno)
  abs <- as.character(abs)
})

library(ggrepel)
library(ggsci)

ggplot(data=aplotdat,aes(y=rho,x=pheno,col=abs))+
  geom_point(aes(size=-log10(p),shape=abs))+
  geom_text_repel(data=subset(aplotdat,p < 0.05/nrow(aplotdat) ), aes(label=signif(p,digits=3)),nudge_x = 0.2)+
  geom_hline(yintercept = 0)+
  scale_color_aaas()+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))

#### scatterplots of individual correlations
phenosofinterest <- as.character(unique(aplotdat$pheno[which(aplotdat$p< 0.05/nrow(aplotdat))]))

plotlist1 <- list()
for (pheno in phenosofinterest) {
    mvar <- paste0("MONO_t_",pheno)
    dvar <- paste0("DLPFC_cells_t_",pheno)

aestext <- paste("aes(x=",mvar,",y=",dvar,")")
aeseval <- eval(parse(text=aestext))

plotlist1[[pheno]] <- ggplot(data=allinfo2, aeseval)+
  geom_point()+
  geom_smooth(method="lm")+
  theme_minimal()

}

plotlist1$parkdx

#### correlation of inter-tissue correlation with monocyte gene effects

#### phrased: "to test if effects of monocyte genes on brain phenotypes are due to mutual regulation or separate processes"
# read in eQTL data from matrixeqtl in ROS/MAP
load("/Users/dfelsky/Documents/monocyte_twas/eqtl/v2/output/alleqtls.Rdata")

dlpfceqtl <- subset(alleqtl, tissue=="dlpfc")
monoeqtl <- subset(alleqtl, tissue=="mono")

monotopeqtls <- lapply(unique(monoeqtl$gene), function(gene) {
  datsub <- monoeqtl[which(monoeqtl$gene==gene),]
  datsub[which(datsub$p==min(datsub$p)),][1,]
})
monotopeqtls <- do.call(rbind,monotopeqtls)

dlpfctopeqtls <- lapply(unique(dlpfceqtl$gene), function(gene) {
  datsub <- dlpfceqtl[which(dlpfceqtl$gene==gene),]
  datsub[which(datsub$p==min(datsub$p)),][1,]
})
dlpfctopeqtls <- do.call(rbind,dlpfctopeqtls)

dqtl <- dlpfctopeqtls[,c(4,6,7,8)]
names(dqtl)[2:4] <- c("t_dlpfc_eqtl","p_dlpfc_eqtl","fdr_dlpfc_eqtl")
mqtl <- monotopeqtls[,c(4,6,7,8)]
names(mqtl)[2:4] <- c("t_mono_eqtl","p_mono_eqtl","fdr_mono_eqtl")

alltopeqtl <- merge(dqtl,mqtl,by="gene",all=T)
alltopeqtl$fishers_p_eqtl <- apply(alltopeqtl[,c(3,6)], 1, function(x) {
  if (is.na(x[1]) || is.na(x[2])) { NA
    } else {
      metp <- metap::sumlog(x)
      metp$p
      }
  }
  )


### merge info again to include eqtl data
allinfo3 <- merge(allinfo2,alltopeqtl,by="gene",all.x=T)
allinfo3$mono_eqtl_yes <- ifelse(allinfo3$fdr_mono_eqtl < 0.05 ,1,0)
allinfo3$dlpfc_eqtl_yes <- ifelse(allinfo3$fdr_dlpfc_eqtl < 0.05 ,1,0)
allinfo3$both_eqtl_yes <- ifelse(allinfo3$fdr_mono_eqtl & allinfo3$fdr_dlpfc_eqtl < 0.05,1,0)

saveRDS(allinfo3,file="/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/Merged_summary_and_correlation_data_witheQTLs_anddlpfccells.rds")

####### read in GTEX data for whole blood, cells, and brain
gtex_allfiles <- list.files("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/eQTL/GTExv8/GTEx_Analysis_v8_eQTL/",full.names = T)

tomatch <- c("Whole_Blood","Brain","Cells")
gtexfiles <- unique(grep(paste(tomatch,collapse="|"),gtex_allfiles, value=TRUE))

allegeneslist <- lapply(gtexfiles, function(x){
  y <- read.table(x,sep="\t",header=T)
  y$tissue <- gsub("\\..*","",basename(x))
  y
})

allegenes <- do.call(rbind,allegeneslist)
rownames(allegenes) <- NULL
allegenes$gene <- gsub("\\..*","",allegenes$gene_id)
allegenes2 <- allegenes[,c("gene","maf","slope","slope_se","pval_perm","pval_beta","qval","pval_nominal_threshold","log2_aFC","tissue")]

library(reshape2)

allegenescast <- reshape(allegenes2,
                         v.names = names(allegenes2)[2:9],
                         timevar = "tissue",
                         idvar = "gene",
                         direction="wide")

allinfo4 <- merge(allinfo3,allegenescast,by="gene",all.x=T)
saveRDS(allinfo4,file="/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/Merged_summary_and_correlation_data_witheQTLs_anddlpfccells_andGTExv8.rds")

######## check correlations of mono effects on each phenotype against eQTL data AND cross-tissue correlations

dset <- subset(ad,cor_p < 0.05)
allcordats2 <- list()
index <- 1

for (pheno in names(ttm)) {
  for (abso in c(TRUE,FALSE)) {
    mvar <- paste0("MONO_t_",pheno)
    
    if (abso==TRUE) {
      mdcor_corr <- cor.test(abs(dset[,mvar]),abs(dset[,"cor_r"]), method="spearman")
      mdcor_corp <- cor.test(abs(dset[,mvar]),abs(-log10(dset[,"cor_p"])), method="spearman")
      mdcor_monoeqtlt <- cor.test(abs(dset[,mvar]),abs(dset[,"t_mono_eqtl"]), method="spearman")
      mdcor_dlpfceqtlt <- cor.test(abs(dset[,mvar]),abs(dset[,"t_dlpfc_eqtl"]), method="spearman")
      mdcor_fisheqtlp <- cor.test(abs(dset[,mvar]),abs(-log10(dset[,"fishers_p_eqtl"])), method="spearman")
    } else {
      mdcor_corr <- cor.test(dset[,mvar],dset[,"cor_r"], method="spearman")
      mdcor_corp <- cor.test(dset[,mvar],-log10(dset[,"cor_p"]), method="spearman")
      mdcor_monoeqtlt <- cor.test(dset[,mvar],dset[,"t_mono_eqtl"], method="spearman")
      mdcor_dlpfceqtlt <- cor.test(dset[,mvar],dset[,"t_dlpfc_eqtl"], method="spearman")
      mdcor_fisheqtlp <- cor.test(dset[,mvar],-log10(dset[,"fishers_p_eqtl"]), method="spearman")
    }
    
    cplotdat <- as.data.frame(rbind(c(mdcor_corr$estimate,p=mdcor_corr$p.value,pheno=pheno,abs=abso),
                                    c(mdcor_corp$estimate,p=mdcor_corp$p.value,pheno=pheno,abs=abso),
                                    c(mdcor_monoeqtlt$estimate,p=mdcor_monoeqtlt$p.value,pheno=pheno,abs=abso),
                                    c(mdcor_dlpfceqtlt$estimate,p=mdcor_dlpfceqtlt$p.value,pheno=pheno,abs=abso),
                                   c(mdcor_fisheqtlp$estimate,p=mdcor_fisheqtlp$p.value,pheno=pheno,abs=abso)))
    
    
    cplotdat$contrast <- factor(c("Cross-tissue cor. rho","Cross-tissue cor. p-value","monocyte eQTL t","DLPFC eQTL t","eQTL Fisher's p-value"))
    
    allcordats2[[index]] <- cplotdat
    index <- index + 1
  }
}



aplotdat2 <- as.data.frame(do.call(rbind,allcordats2))
aplotdat2 <- within(aplotdat2,{
  rho <- as.numeric(rho)
  p <- as.numeric(p)
  pheno <- as.factor(pheno)
  abs <- as.character(abs)
})

library(ggrepel)
library(ggsci)

ggplot(data=aplotdat2,aes(y=rho,x=pheno,col=abs))+
  geom_point(aes(size=-log10(p),shape=abs))+
  geom_text_repel(data=subset(aplotdat2,p < 0.05/nrow(aplotdat2) ), aes(label=signif(p,digits=3)),nudge_x = 0.2)+
  geom_hline(yintercept = 0)+
  facet_wrap(~ contrast)+
  scale_color_aaas()+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))


ggplot(data=subset(aplotdat2,contrast=="Cross-tissue cor. rho"),aes(y=rho,x=pheno,col=abs))+
  geom_point(aes(size=-log10(p),shape=abs))+
  geom_text_repel(data=subset(aplotdat2,p < 0.05/nrow(aplotdat2) & contrast=="Cross-tissue cor. rho" ), aes(label=signif(p,digits=3)),nudge_x = 0.2)+
  geom_hline(yintercept = 0)+
  scale_color_aaas()+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))



######## test significant mono genes for pathology against cor_r_cells
# 

mfgenes <- ad$gene[which(ad$MONO_adj.P.Val_mf3123<0.05)]

boxplot(ad$cor_r_celltype ~ ifelse(ad$MONO_adj.P.Val_plaq_d_sqrt<0.1,1,0))

ggplot(data=ad, aes(y=MONO_t_tdp_stage4,x=DLPFC_cells_t_tdp_stage4))+
  geom_point()+
  geom_smooth(method="lm")+
  theme_minimal()


#### Plot cross-tissue correlation with eQTL data

eqtl_cor <- ggplot(data=allinfo4, aes(x=cor_r,y=abs(slope.Brain_Frontal_Cortex_BA9/slope_se.Brain_Frontal_Cortex_BA9),fill=cor_sig))+
  geom_point(aes(col=cor_sig))+
  scale_color_manual(values = colors()[c(12,30,54)])+
  geom_text_repel(data=subset(allinfo4,cor_fdr < 0.05),aes(label=hugo),size=2.5)+
  geom_vline(xintercept = c(0,median(allinfo4$cor_r)),col=c("black","orange"),lty=c(2,4))+
  theme_minimal()

eqtl_cor_celltype <- ggplot(data=allinfo4, aes(x=cor_r_celltype,y=abs(slope.Brain_Frontal_Cortex_BA9/slope_se.Brain_Frontal_Cortex_BA9),fill=cor_sig_celltype))+
  geom_point(aes(col=cor_sig_celltype))+
  scale_color_manual(values = colors()[c(12,30,54)])+
  geom_text_repel(data=subset(allinfo4,cor_fdr_celltype < 0.05),aes(label=hugo),size=2.5)+
  geom_vline(xintercept = c(0,median(allinfo4$cor_r_celltype)),col=c("black","orange"),lty=c(2,4))+
  theme_minimal()

cowplot::plot_grid(eqtl_cor,eqtl_cor_celltype,ncol=2,labels = c("No cell type correction","With cell type correction"))

### cross-tabs

qtlvars <- grep("slope",names(ad),value=T)
qtlvars1 <- grep("slope_se",qtlvars,value=T)
qtlvars2 <- grep("slope_se",qtlvars,value=T,invert = T)

allres <- list()
for (tnum in seq(1,length(qtlvars1))) {
  tstat <- abs(ad[,qtlvars2[tnum]]/ad[,qtlvars1[tnum]])
   ctest <- cor.test(ad$cor_r_celltype,tstat)
   allres[[tnum]] <- c(ctest$estimate,ctest$p.value)
}
 
pdat1 <- as.data.frame(do.call(rbind,allres))
names(pdat1)[2] <- "p"
pdat1$tissue <- gsub(".*\\.","",qtlvars1)

pdat1$tissue <- factor(pdat1$tissue,levels=pdat1$tissue[order(pdat1$p)])
ggplot(data=pdat1, aes(y=cor,x=tissue))+
  geom_bar(stat="identity")+
  coord_flip()+
  theme_minimal()
  

boxplot(data=ad, abs(slope.Brain_Frontal_Cortex_BA9/slope_se.Brain_Frontal_Cortex_BA9) ~ cor_sig_celltype)
summary(lm(data=ad, abs(slope.Brain_Frontal_Cortex_BA9/slope_se.Brain_Frontal_Cortex_BA9) ~ cor_sig_celltype))


ad$gtex_brainQTL <- ifelse(ad$pval_perm.Brain_Frontal_Cortex_BA9<0.05,1,0)

fisher.test(table(ad$gtex_brainQTL,ad$cor_sig_celltype))

ggplot(data=ad,aes(y=abs(slope.Brain_Frontal_Cortex_BA9/slope_se.Brain_Frontal_Cortex_BA9),x=cor_sig_celltype))+
  geom_violin(draw_quantiles = T)+
  geom_jitter(width=0.1,alpha=0.2,size=0.5)+
  theme_minimal()


summary(lm(data=ad, abs(slope.Brain_Frontal_Cortex_BA9/slope_se.Brain_Frontal_Cortex_BA9) ~ cor_sig_celltype))
