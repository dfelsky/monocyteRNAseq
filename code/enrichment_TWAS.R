#### perform enrichment on all TWAS results
#### module enrichment #####
library(org.Hs.eg.db)
library(gprofiler2)
library(rrvgo)

### read TWAS results
alltt <- readRDS("output/all_TWAS_results.rds")

### get summary numbers of significant DE genes
enrichment_results <- lapply(alltt, function(tissue) {
  lapply(tissue, function(pheno) {
    
    siggenesup <- pheno$gene[which(pheno$adj.P.Val < 0.05 & pheno$t > 0)]
    if (length(siggenesup)==1) { gup <- siggenes 
    } else if (length(siggenesup)==0) { gup <- NA 
    } else {
      gup <- gost(query=siggenesup,custom_bg = pheno$gene, significant = T)
      }
    
    siggenesdown <- pheno$gene[which(pheno$adj.P.Val < 0.05 & pheno$t < 0)]
    if (length(siggenesdown)==1) { gdown <- siggenesdown
    } else if (length(siggenesdown)==0) { gdown <- NA
    } else {
      gdown <- gost(query=siggenesdown,custom_bg = pheno$gene, significant = T)
    }
    return(list(up=gup,down=gdown))
    })
  })

saveRDS(enrichment_results,"output/all_TWAS_enrichment_results.rds")


#### get similarity in GO terms, BP
sim_method <- "Rel"
sim_threshold <- 0.7

enrichment_results <- readRDS("output/all_TWAS_enrichment_results.rds")

BPresults <- list()
for (tissue in names(enrichment_results)[1:2]) {  ## only perform full analyis for monocyte and monocyte_blood. dlpfc too many sig terms
  print(paste0("tissue: ",tissue))
  for (pheno in names(enrichment_results[[tissue]])) {
    print(paste0("pheno: ",pheno))
    phenores <- enrichment_results[[tissue]][[pheno]]    
    
    BPresults[[tissue]][[pheno]] <- lapply(phenores, function(updown) {
          
          if (is.list(updown)==T) {
            ressub <- subset(updown$result, updown$result$source %in% c("GO:BP"))
            
            if (nrow(ressub) == 1) {
            return(list(raw_results=BP_res_sub))
              } else {
                simMatrix <- calculateSimMatrix(ressub$term_id,
                                            orgdb="org.Hs.eg.db",
                                            ont="BP",
                                            method=sim_method)
                scores <- setNames(-log10(ressub$p_value), ressub$term_id)
            
                reducedTerms <- reduceSimMatrix(simMatrix,
                                            scores,
                                            threshold=sim_threshold,
                                            orgdb="org.Hs.eg.db")
            
                return(list(raw_results=ressub,
                            simMatrix=simMatrix,
                            scores=scores,
                            reducedTerms=reducedTerms,
                            meta=list(method=sim_method,
                                      threshold=sim_threshold)))
              }
          } else {
            print("skipping because enrichment result is not a list")
          }
        })
  }
}

saveRDS(BPresults,"output/TWAS_enrichment_results_bothmonocyte_BP.rds")
