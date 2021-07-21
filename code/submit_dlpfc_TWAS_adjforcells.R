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


roslong <- read.csv("/Users/dfelsky/Documents/data/ROSMAP_phenotype_data/update_10162020/dataset_978_long_10-13-2020.csv",colClasses = c(projid="character"))

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

##################################
#### set parameters for TWAS  ####
##################################
outdir <- "/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/TWAS/Jan2021/newcog/"
setwd(outdir)

outcome.category.list <- c("cognition")
tissue.list <- c("dlpfc")

techvars.mono <- c("batch","PCT_USABLE_BASES","PERCENT_DUPLICATION","MEDIAN_3PRIME_BIAS","study","PCT_PF_READS_ALIGNED","ESTIMATED_LIBRARY_SIZE")
techvars.dlpfc <- c("batch", "MEDIAN_CV_COVERAGE", "PCT_RIBOSOMAL_BASES", "PCT_CODING_BASES", "PCT_UTR_BASES", "LOG_ESTIMATED_LIBRARY_SIZE", "LOG_PF_READS_ALIGNED", "MEDIAN_5PRIME_TO_3PRIME_BIAS", "PCT_PF_READS_ALIGNED", "study", "PERCENT_DUPLICATION", "MEDIAN_3PRIME_BIAS", "PCT_INTERGENIC_BASES")

covars.mono.pathology <- c("msex","age_death","pmi","age_draw")
covars.mono.cognition <- c("msex","age_draw","educ")

covars.dlpfc.pathology <- c("msex","pmi","age_death")
covars.dlpfc.cognition <- c("educ","msex","age_death","age_at_visit_at_lastvisit")

indepvec.pathology <- c("plaq_n_sqrt","plaq_d_sqrt","amyloid_sqrt","tangles_sqrt","nft_sqrt","ci_num2_gct","ci_num2_mct","arteriol_scler","caa_4gp","cvda_4gp2","dlbdx","hspath_any","tdp_stage4","parkdx","pathoAD","vm3123","pput3123","it3123","mf3123") 

indepvec.cognition <- names(newdatsub2)[-c(1,9)] #c("cogn_global_random_slope","cogn_ep_random_slope","cogn_po_random_slope","cogn_ps_random_slope","cogn_se_random_slope","cogn_wo_random_slope")

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
  
      saveRDS(dge_filtered,file="dlpfc_dge_filtered_used_celltypes.rds")
      
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
        removesubs.mono <- c("00246264","35072859","91804757","21135680","50303145","50107907","10100574","51520126")
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
      
      ### adjust for cell types
      genesymbols <- all_genes$external_gene_name[match(rownames(v$E),all_genes$ensembl_gene_id)]
      rownames(v$E) <- genesymbols
      v$E <- v$E[!is.na(genesymbols),]
      v$weights <- v$weights[!is.na(genesymbols),]
      dlpfc_celladj <- adjustBrainCells(v$E,species = "human")
      v$E <- dlpfc_celladj$expression
      rownames(v$E) <- all_genes$ensembl_gene_id[match(rownames(v$E),all_genes$external_gene_name)]
      
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
      tt$gene <- tt$ID
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
    saveRDS(ttlistdesign,file=paste0(outdir,"TWAS_designmatrices_",tissue,"_",outcome.category,"_celltypes.rds"))
    saveRDS(ttlist,file=paste0(outdir,"TWAS_results_",tissue,"_",outcome.category,"_celltypes.rds"))
    
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
    pdf(paste0(outdir,"Manhattan_plot_",tissue,"_",outcome.category,"_celltypes.pdf"),width=14,height=9, onefile = T)
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
    pdf(paste0(outdir,"Individual_effects_",tissue,"_",outcome.category,"_celltypes.pdf"),width=16,height=16, onefile = T)
    for (a in seq(1,length(datplotlist))) { print(datplotlist[[a]]) }
    dev.off()
  }
}

