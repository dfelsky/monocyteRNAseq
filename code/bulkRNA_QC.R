### bulk RNAseq QC function
bulk_qc <- function(logoffset=0.5,
                    pc.bygroup=FALSE,
                    min.count=15,
                    IQR_sub_threshold=3,
                    IQR_obs_threshold=3,
                    IQR_highprune_threshold=1.5,
                    IQR_varprune_threshold=1.5,
                    variance.prune=FALSE,
                    high.prune=FALSE,
                    coerce.obs=FALSE,
                    ydge) {
  
  require(edgeR)
  require(rms)
  require(limma)
  
  # save input params
  print("Saving input parameters")
  returnparams <- data.frame(paramter=c("logoffset","pc.bygroup","min.count","IQR_sub_threshold","IQR_obs_threshold","IQR_highprune_threshold","IQR_varprune_threshold","variance.prune","high.prune","coerce.obs"),
                             value=as.character(c(logoffset,pc.bygroup,min.count,IQR_sub_threshold,IQR_obs_threshold,IQR_highprune_threshold,IQR_varprune_threshold,variance.prune,high.prune,coerce.obs)),
                             description=c("count offset applied before voom transformation - counts are transformed before identification of individual outlying observations",
                                           "logical indicator for whether or not to perform PCA on all subjects together (FALSE) or separately (TRUE; default) according to group variable in DGElist object",
                                           "threshold at which all genes with lower median expected count (before transformation) will be pruned, not inclusive",
                                           "multiplier of IQR beyond which any observed PC value will be tagged as an outlier during subject pruning",
                                           "multiplier of IQR at which any observed voom-transformed count value will be deemed a single outlying observation and coerced up to minimum, or down to maximum, of remaining observed values",
                                           "multiplier of IQR at which median voom-transformed count for a gene, out of the entire set of genes, will be deemed an outlier and removed, if high.prune=TRUE",
                                           "multiplier of IQR at which voom-transformed count variance for a gene, out of the entire set of genes, will be deemed an outlier and removed, if variance.prune=TRUE",
                                           "prune genes based on excessive variance","prune genes based on excessively high counts",
                                           "coerce indiviual expression values toward centrality, per gene, based on outlier status defined by IQR_obs_threshold"))
  
  
  # remove lowly expressed genes using median_threshold before voom transform
  print("Removing genes with low absolute counts")
  y <- ydge 
  med <- apply(ydge$counts,1,median)
  keep <- (med >= min.count)
  y.genes <- y[keep,]
  
  # save stats for pruning
  sk <- describe(keep)
  returnprune <- paste0("Kept ",sk$values$frequency[2]," out of ",length(keep)," genes after pruning for median raw count.")
  print(returnprune)
  
  # if high.prune=TRUE then prune for overly expressed genes
  if (high.prune) {
    print("Removing genes with high median log2CPM")
    ygene.voom <- voom(y.genes)
    medvoom <- apply(ygene.voom$E,1,median)
    prune <- unique(boxplot(medvoom,range=IQR_highprune_threshold)$out)
    prune <- prune[which(prune > median(medvoom))]
    prunegenes <- rownames(ygene.voom$E)[which(medvoom %in% prune)]
    keep2 <- rownames(y.genes$counts) %nin% prunegenes
    y.genes <- y.genes[keep2,]
    
    hp <- describe(keep2)
    returnprunehigh <- paste0("Kept ",hp$values$frequency[2]," out of ",length(keep2)," genes after pruning for high median log2 expression.")
    print(returnprunehigh)
  } else { returnprunehigh <- "pruning for high expression not performed"
  print(returnprunehigh) }
  
  # if variance.prune=TRUE then prune based on outlying variance
  if (variance.prune) {
    print("Removing genes with high log2CPM variance")
    ygene.voom <- voom(y.genes)
    varvoom <- apply(ygene.voom$E,1,var)
    prune <- unique(boxplot(varvoom,range=IQR_varprune_threshold)$out)
    prune <- prune[which(prune > median(varvoom))]
    prunegenes <- rownames(ygene.voom$E)[which(varvoom %in% prune)]
    keep3 <- rownames(y.genes$counts) %nin% prunegenes
    y.genes <- y.genes[keep3,]
    
    vp <- describe(keep3)
    returnprunevar <- paste0("Kept ",vp$values$frequency[2]," out of ",length(keep3)," genes after pruning for high median log2 expression variance.")
    print(returnprunevar)
  } else { returnprunevar <- paste0("pruning for high variance not performed")
  print(returnprunevar) }
  
  # prune subjects based on large scale expression patterns
  if (pc.bygroup==TRUE) {
    print("pruning subjects by PCA according to batch")
    groupindex <- 1
    outsubs.bypc <- list()
    
    for (group in unique(y.genes$samples$group)) {
      ysub <- y.genes[,which(y.genes$samples$group==group)]
      mds <- plotMDS(ysub,top = 5000,plot=T,ndim = 2)
      
      pc1out <- unique(boxplot(mds$x,range=IQR_sub_threshold)$out)
      outsubs1 <- rownames(y.genes$samples)[which(mds$x %in% pc1out)]
      pc2out <- unique(boxplot(mds$y,range=IQR_sub_threshold)$out)
      outsubs2 <- rownames(y.genes$samples)[which(mds$y %in% pc2out)]
      
      outsubs.bypc[[groupindex]] <- unique(c(outsubs1,outsubs2))
      groupindex <- groupindex+1
    }
    outsubs <- unlist(outsubs.bypc)
  } else {
    print("pruning subjects by PCA NOT according to batch")
    mds <- plotMDS(y.genes,top=5000,plot=F,ndim=2)
    
    pc1out <- unique(boxplot(mds$x,range=IQR_sub_threshold)$out)
    outsubs1 <- rownames(y.genes$samples)[which(mds$x %in% pc1out)]
    pc2out <- unique(boxplot(mds$y,range=IQR_sub_threshold)$out)
    outsubs2 <- rownames(y.genes$samples)[which(mds$y %in% pc2out)]
    
    outsubs.bypc[[groupindex]] <- unique(c(outsubs1,outsubs2))
    
  }
  
  y.genes.subs <- y.genes[,colnames(y.genes$counts) %nin% outsubs]
  ygs <- y.genes.subs
  
  # save pruned subjects for reporting
  returnsub <- list(paste0("Removed ",length(outsubs)," subjects: ",paste0(outsubs,collapse = ", ")), paste0("Pre-QC n=",dim(y.genes)[2],". Post-QC n=",dim(y.genes.subs)[2]))
  print(returnsub)
  
  
  if (coerce.obs==T) { 
  # now,  coerce observations for outliers for each gene
  print("coercing extreme values to max/min indiviually per gene")
  
  # first, get numbers of outliers per gene
  outlier.obs.counts <- apply(log2(ygs$counts+logoffset),1,function(x) {
    outliers <- c(boxplot(x, plot=FALSE,range=IQR_obs_threshold)$out)
    outliers.high <- outliers[which(outliers > median(x))]
    outliers.low <- outliers[which(outliers < median(x))]
    return(c(length(outliers.low),length(outliers.high)))
  })
  
  low.outlier.range <- range(t(outlier.obs.counts)[,1])[2]
  low.outlier.quant <- quantile(t(outlier.obs.counts)[,1][which(t(outlier.obs.counts)[,1]>0)])
  high.outlier.range <- range(t(outlier.obs.counts)[,2])[2]
  high.outlier.quant <- quantile(t(outlier.obs.counts)[,2][which(t(outlier.obs.counts)[,2]>0)])
  n.low.outliers <- summary(t(outlier.obs.counts)[,1]>0)[3]
  n.high.outliers <- summary(t(outlier.obs.counts)[,2]>0)[3]
  n.total.genes <- ncol(outlier.obs.counts)
  
  print(paste0("number of genes detected with >0 low outliers = ",n.low.outliers,"/",n.total.genes))
  print(paste0("number of genes detected with >0 high outliers = ",n.high.outliers,"/",n.total.genes))
  print(paste0("number of observations adjusted for low outliers = 1-",low.outlier.range,". quantiles=",paste0(low.outlier.quant,collapse = ",")))
  print(paste0("number of observations adjusted for high outliers = 1-",high.outlier.range,". quantiles=",paste0(high.outlier.quant,collapse=",")))
  
  
  # perform coercion for individual observations across all genes
  newcounts <- t(apply(log2(ygs$counts+logoffset),1,function(x) {
    outliers <- c(boxplot(x, plot=FALSE,range=IQR_obs_threshold)$out)
    outliers.high <- outliers[which(outliers > median(x))]
    outliers.low <- outliers[which(outliers < median(x))]
    
    x2 <- x
    x2[x %in% outliers] <- NA
    lowval <- min(x2,na.rm=T)
    highval <- max(x2,na.rm=T)
    
    x[x %in% outliers.high] <- highval
    x[x %in% outliers.low] <- lowval
    (2^x)-logoffset
  }))
  
  ygs2 <- ygs
  ygs2$counts <- newcounts
  } else {
    print("coercion of individual observations not performed")
    ygs2 <- ygs
  }
  
  
  ts <- gsub(timestamp(quiet = T),pattern="[[:punct:]]|[[:space:]]",replacement="")
  print(paste0("completed ",timestamp(quiet=T)))
  
  ############ quick check to see if there are missing values for phenotypes
  
  
  y2 <- ygs2[,which(complete.cases(ygs2$samples)==T)]
  if (dim(y2)[2] < dim(ygs2)[2]) { print("WARNING: DGE object phenotype data contains missing values")}
  
  return(list(YDGE=ygs2,returnprune=returnprune,returnprunehigh=returnprunehigh,returnprunevar=returnprunevar,returnsub=returnsub,timestamp=ts,returnparams=returnparams))
}