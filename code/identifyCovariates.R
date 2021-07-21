# This file provides function to identify covariates
# that should be used to adjust for technical variables when
# running a TWAS. The functions are using voom/edgeR and require
# an  dgeList object with outlier samples and weakly
# transcribed genes already removed. Normalization factors
# should be calculated before calling these functions.

library(edgeR)
library(ggplot2)

# Returns matrix with correlations and p-values
getCorrelations <- function(data, x, y) {
  
  stopifnot(all(c(x, y) %in% colnames(data)))
  
  corMat <- matrix(NA, nrow=length(x), ncol=length(y)) # (symmetric)
  pMat   <- matrix(NA, nrow=length(x), ncol=length(y)) # (symmetric)
  rownames(corMat) <- x
  colnames(corMat) <- y
  rownames(pMat) <- x
  colnames(pMat) <- y
  
  for (i in 1:length(x)) {
    for (j in 1:length(y)) {
      
      # Both variables are identical
      if (x[i] == y[j]) {
        # do nothing - keep NA 
        
        # If both are numerical: Pearson correlation
      } else if (is.numeric(data[, x[i]]) & is.numeric(data[, y[j]])) {
        corMat[i, j] <- cor(data[, x[i]], data[, y[j]])
        pMat[i, j] <- cor.test(data[, x[i]], data[, y[j]])$p.value
        
        # One is numeric (x), one is factor (y): One-way anova
      } else if (is.numeric(data[, x[i]]) & is.factor(data[, y[j]])) {
        if ( shapiro.test(data[, x[i]])$p.value > 0.05){
          a <- summary(aov(data[, x[i]] ~ data[, y[j]]))
          corMat[i, j] <- sqrt(a[[1]][1, "Sum Sq"] / sum(a[[1]][, "Sum Sq"])) 
          pMat[i, j]   <- a[[1]][1, "Pr(>F)"]  
        } else {
          a <- kruskal.test(data[, x[i]] ~ data[, y[j]], data = data)
          corMat[i, j] <- sqrt(a$statistic / length(data[, x[i]]) )
          pMat[i, j]   <- a$p.value 
        }
        
        # One is numeric (y), one is factor (x): One-way anova
      } else if (is.factor(data[, x[i]]) & is.numeric(data[, y[j]])) {
         if ( shapiro.test(data[, y[j]])$p.value > 0.05){
          a <- summary(aov(data[, y[j]] ~ data[, x[i]]))
            corMat[i, j] <- sqrt(a[[1]][1, "Sum Sq"] / sum(a[[1]][, "Sum Sq"])) 
          pMat[i, j]   <- a[[1]][1, "Pr(>F)"]  
        } else {
          a <- kruskal.test(data[, y[j]] ~ data[, x[i]], data = data)
          corMat[i, j] <- sqrt(a$statistic / length(data[, x[i]]) )
          pMat[i, j]   <- a$p.value 
        }
        
        # Both are categorical: Fisher's exact test and Cramer's V
      } else if (is.factor(data[, x[i]]) & is.factor(data[, y[j]])) {
           suppressWarnings(corMat[i, j] <- sqrt(chisq.test(data[, x[i]], data[, y[j]], correct=FALSE)$statistic /
                                                (length(data[, x[i]]) * (min(length(unique(data[, x[i]])),length(unique(data[, y[j]]))) - 1))))
        if ( all( suppressWarnings(chisq.test( table(data[, x[i]], data[, y[j]]) )$expected ) >= 5) ){
          pMat[i, j] <- chisq.test(data[, x[i]], data[, y[j]], correct=FALSE)$p.value
        } else {
          pMat[i, j] <- fisher.test(data[, x[i]], data[, y[j]], simulate.p.value=TRUE)$p.value
        }
      }
    }
  }
  
  return(list(c=corMat, p=pMat))
}


# Return principal components
getPCs.dge <- function(dge, n=10) {
  exprs <- cpm(dge, log=TRUE, prior.count=1)
  pc <- prcomp(x=t(exprs), retx=TRUE, center=TRUE, scale.=FALSE)
  return(pc$x[, 1:n])
}
getPCs.matrix <- function(mat, n=10) {
  pc <- prcomp(x=t(mat), retx=TRUE, center=TRUE, scale.=FALSE)
  return(pc$x[, 1:n])
}

# Run voom and calculate residuals (and number of sign. genes)
getResiduals <- function(dge, vars, normalize="none", values=values) {
  model <- formula(paste("~ ", paste(vars, collapse=" + ")))
  design <- model.matrix(model, data=dge$samples)
  
  if("tpm" %in% values) {
    v <- t(log2(t(dge$counts)/dge$samples$norm.factors+0.01))
  }
  if("counts" %in% values) {
    v <- voom(dge, design, plot=FALSE, normalize.method=normalize)
  }
  fit <- lmFit(v, design)
  fit <- eBayes(fit)
  resids <- residuals(fit, y=v)
  
  numSig <- rep(NA, times=length(vars))
  names(numSig) <- vars
  
  getCoefNames <- function (var) { # return names of all coefficients of given variable
    if (is.numeric(dge$samples[, var])) {
      return(var)
    } else if (is.factor(dge$samples[, var])) {
      return(paste(var, levels(dge$samples[, var]), sep="")[-1])
    }
  }
  
  for (i in 1:length(vars)) {
    if (grepl(":", vars[i])) { # Interaction term
      iVars <- strsplit(vars[i], ":")[[1]]
      iVars1 <- getCoefNames(iVars[1])
      iVars2 <- getCoefNames(iVars[2])
      coefs <- levels(interaction(iVars1, iVars2, sep=":"))
    } else {
      coefs <- getCoefNames(vars[i])
    }
    numSig[i] <- sum(topTable(fit, coef=coefs, number=nrow(dge))$adj.P.Val <= 0.05)
  }
  
  return(list(resids=resids, numSig=numSig))
}


# Adds covariates associated with residual PCs iteratively
# Interaction effects can be given as "varsFix", but not as optional variables as
# it is unclear how to calculate their correlation with PCs.
detectSignCovariates <- function (dge, varsFix, varsOpt, normalize="none", npc=10, maxIter=10, values=values) {
  
  varsSelected <- varsFix
  varsOut <- varsOpt
  iter <- 0
  update <- TRUE
  
  summary <- data.frame()
  
  while (update == TRUE & iter <= maxIter) {
    indInteraction <- grepl(":", varsSelected)
    resids <- getResiduals(dge, vars=c(varsSelected, varsFix), normalize=normalize, values=values)
    pcs <- getPCs.matrix(resids$resids, n=npc)
    cor <- getCorrelations(data=cbind(dge$samples, pcs) , x=c(varsSelected[!indInteraction], varsOut) , y=colnames(pcs))
    minP <- apply(cor$p, 1, min)
    
    numSig <- rep(NA, times=length(c(varsSelected, varsOut)))
    names(numSig) <- c(varsSelected, varsOut)
    numSig[names(resids$numSig)] <- resids$numSig
    pcminp=minP[c(varsSelected, varsOut)]
    names(pcminp) <- c(varsSelected, varsOut)
    summary <- rbind(summary, data.frame(variable=c(varsSelected, varsOut),
                                         iteration=iter,
                                         selected=c(varsSelected, varsOut) %in% varsSelected,
                                         pc=pcminp,
                                         numSig=numSig))
    
    if (length(varsOut) > 0) {
      minP <- minP[names(minP) %in% varsOut]
      nextVar <- names(minP)[which.min(minP)]
      if (minP[nextVar] > 0.05) {
        update <- FALSE
      } else {
        varsOut <- setdiff(varsOut, nextVar)
        varsSelected <- c(varsSelected, nextVar)
        iter <- iter + 1
      }
    } else {
      update <- FALSE
    }
  }
  
  row.names(summary) <- NULL
  return(summary)
}
