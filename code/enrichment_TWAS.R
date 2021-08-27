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
    if (length(siggenesup)==1) { gup <- siggenesup 
    } else if (length(siggenesup)==0) { gup <- NA 
    } else {
      gup <- gost(query=siggenesup,custom_bg = pheno$gene, significant = T,sources = "GO:BP")
      }
    
    siggenesdown <- pheno$gene[which(pheno$adj.P.Val < 0.05 & pheno$t < 0)]
    if (length(siggenesdown)==1) { gdown <- siggenesdown
    } else if (length(siggenesdown)==0) { gdown <- NA
    } else {
      gdown <- gost(query=siggenesdown,custom_bg = pheno$gene, significant = T, sources="GO:BP")
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
for (tissue in names(enrichment_results)[2]) {  ## only perform full analysis for monocyte_blood
  print(paste0("tissue: ",tissue))
  for (pheno in names(enrichment_results[[tissue]])) {
    print(paste0("pheno: ",pheno))
    phenores <- enrichment_results[[tissue]][[pheno]]    
    
    BPresults[[tissue]][[pheno]] <- lapply(phenores, function(updown) {
          
          if (is.list(updown)==T) {
            ressub <- subset(updown$result, updown$result$source %in% c("GO:BP"))
            ressub <- subset(ressub,significant=="TRUE")
            
            if (nrow(ressub) < 1) { return(NA)
            } else if (nrow(ressub) == 1)  {
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

#saveRDS(BPresults,"output/TWAS_enrichment_results_monocyte_blood_BP.rds")

#################################
#################################
#################################

BPresults <- readRDS("output/TWAS_enrichment_results_monocyte_blood_BP.rds")

ressub <- eres$vm3123$up$result

simMatrix <- calculateSimMatrix(ressub$term_id,
                                orgdb="org.Hs.eg.db",
                                ont="BP",
                                method="Rel")

scores <- setNames(-log10(ressub$p_value), ressub$term_id)

reducedTerms <- reduceSimMatrix(simMatrix,
                                scores,
                                threshold=sim_threshold,
                                orgdb="org.Hs.eg.db")


# heatmapPlot(simMatrix,
#             reducedTerms,
#              annotateParent=TRUE,
#             annotationLabel="parentTerm",
#             fontsize=6)

treemapPlot(reducedTerms)


####################################
#### DOSE and clusterProfiler version ### below is underdeveloped - used for rank-based GSEA, which is inot implemented in gprofiler. This analysis is likely unnecessary since the enrichment of significant FDR<0.05 gene sets with gprofiler should suffice. Also gprofiler interface is better for looking at a multitude of reference ontologies/gene sets.
#### module enrichment #####
library(clusterProfiler)
library(DOSE)
library(AnnotationDbi)
library(org.Hs.eg.db)
library(enrichplot)
library(gridExtra)
library(rms)
library(data.table)
library(ggraph)

geneList <- -log10(ttall_unlabelled$monocyte$mf3123$P.Value)
names(geneList) <- ttall_unlabelled$monocyte$mf3123$gene
forgo <- mapIds(org.Hs.eg.db,
                         keys=names(geneList), 
                         column="ENTREZID",
                         keytype="ENSEMBL",
                         multiVals="first")
names(geneList) <- forgo


ego <- gseGO(geneList     = geneList,
              OrgDb        = org.Hs.eg.db,
              ont          = "BP",
              minGSSize    = 20,
              maxGSSize    = 300,
              pvalueCutoff = 0.05,
             eps=0,
              verbose      = TRUE)
