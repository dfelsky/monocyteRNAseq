#### module enrichment #####
library(org.Hs.eg.db)
library(gprofiler2)
library(rrvgo)

#### Monocytes
#### load module definitions
net <- readRDS("output/WGCNA/mono/net_noPAM.rds")

####
modlist <- lapply(unique(net$colors), function(module) {
  names(net$colors)[which(net$colors==module)]
})
names(modlist) <- unique(net$colors)
modlist$grey <- NULL

multi_gostres <- gost(query=modlist,
                      multi_query = F,
                      custom_bg = names(net$colors))

#### get similarity in GO terms, BP
enrich_results <- list()
sim_method <- "Rel"
sim_threshold <- 0.7

BP_res <- subset(multi_gostres$result, multi_gostres$result$source %in% c("GO:BP"))

for (module in unique(BP_res$query)) {
  
  BP_res_sub <- subset(BP_res, query==module)
  
  if (nrow(BP_res_sub) == 1) {
    enrich_results[[module]] <- list(raw_results=BP_res_sub)
  } else {
  
  simMatrix <- calculateSimMatrix(BP_res_sub$term_id,
                                  orgdb="org.Hs.eg.db",
                                  ont="BP",
                                  method=sim_method)
  
  scores <- setNames(-log10(BP_res_sub$p_value), BP_res_sub$term_id)
  
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=sim_threshold,
                                  orgdb="org.Hs.eg.db")
  
  enrich_results[[module]] <- list(raw_results=BP_res_sub,
                                   simMatrix=simMatrix,
                                   scores=scores,
                                   reducedTerms=reducedTerms,
                                   meta=list(method=sim_method,
                                             threshold=sim_threshold))
  }
}

saveRDS(enrich_results,"output/WGCNA/mono/enrichment_results.rds")

# heatmapPlot(simMatrix,
#             reducedTerms,
#             annotateParent=TRUE,
#             annotationLabel="parentTerm",
#             fontsize=6)
# 
# treemapPlot(reducedTerms)


#### DLPFC
#### load module definitions
net_dlpfc <- readRDS("output/WGCNA/dlpfc/net_noPAM_celltypes.rds")

####
modlist <- lapply(unique(net_dlpfc$colors), function(module) {
  names(net_dlpfc$colors)[which(net_dlpfc$colors==module)]
})
names(modlist) <- unique(net_dlpfc$colors)
modlist$grey <- NULL

multi_gostres <- gost(query=modlist,
                      multi_query = F,
                      custom_bg = names(net_dlpfc$colors))

#### get similarity in GO terms, BP
enrich_results <- list()
sim_method <- "Rel"
sim_threshold <- 0.7

BP_res <- subset(multi_gostres$result, multi_gostres$result$source %in% c("GO:BP"))

for (module in unique(BP_res$query)) {
  
  BP_res_sub <- subset(BP_res, query==module)
  
  if (nrow(BP_res_sub) == 1) {
    enrich_results[[module]] <- list(raw_results=BP_res_sub)
  } else {
    
    simMatrix <- calculateSimMatrix(BP_res_sub$term_id,
                                    orgdb="org.Hs.eg.db",
                                    ont="BP",
                                    method=sim_method)
    
    scores <- setNames(-log10(BP_res_sub$p_value), BP_res_sub$term_id)
    
    reducedTerms <- reduceSimMatrix(simMatrix,
                                    scores,
                                    threshold=sim_threshold,
                                    orgdb="org.Hs.eg.db")
    
    enrich_results[[module]] <- list(raw_results=BP_res_sub,
                                     simMatrix=simMatrix,
                                     scores=scores,
                                     reducedTerms=reducedTerms,
                                     meta=list(method=sim_method,
                                               threshold=sim_threshold))
  }
}

saveRDS(enrich_results,"output/WGCNA/dlpfc/enrichment_results.rds")

######### summarize and plot enrichment results

eres_mono <- readRDS("output/WGCNA/mono/enrichment_results.rds")
eres_dlpfc <- readRDS("output/WGCNA/dlpfc/enrichment_results.rds")




