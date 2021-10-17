#### perform enrichment on all TWAS results
#### module enrichment #####
library(rms)
library(org.Hs.eg.db)
library(gprofiler2)
library(rrvgo)
library(GGally)
library(network)
library(sna)
library(ggnetwork)
library(reshape2)
library(stringr)
library(ggsci)
library(gplots)
library(ggrepel)
library(ggthemes)
library(ggpubr)
library(ggsci)
library(ggraph)
library(clusterProfiler)
library(DOSE)
library(AnnotationDbi)
library(enrichplot)
library(gridExtra)
library(data.table)

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


########################################################################
########################################################################
#### DOSE and clusterProfiler version 
### below is underdeveloped - used for rank-based GSEA, which is inot implemented in gprofiler. This analysis is likely unnecessary since the enrichment of significant FDR<0.05 gene sets with gprofiler should suffice. Also gprofiler interface is better for looking at a multitude of reference ontologies/gene sets.
#### module enrichment #####
###### overrepresentation
####
universe <- unique(as.vector(na.omit(mapIds(org.Hs.eg.db,
                                        keys=alltt$monocyte_blood$plaq_n_sqrt$gene, 
                                        column="ENTREZID",
                                        keytype="ENSEMBL",
                                        multiVals="first"))))

phenolist <- lapply(unique(names(alltt$monocyte_blood)), function(pheno) {
  ats <- alltt$monocyte_blood[[pheno]]
  
  upensg <- ats$gene[which(ats$adj.P.Val < 0.05 & ats$t > 0)]
  if(length(upensg)<1) { upentrez <- NA 
  } else {
  upentrez <- as.vector(na.omit(mapIds(org.Hs.eg.db,
                           keys=upensg, 
                           column="ENTREZID",
                           keytype="ENSEMBL",
                           multiVals="first")))
  }
  
  downensg <- ats$gene[which(ats$adj.P.Val < 0.05 & ats$t < 0)]
  if(length(downensg)<1) { downentrez <- NA 
  } else {
  downentrez <- as.vector(na.omit(mapIds(org.Hs.eg.db,
                                      keys=downensg, 
                                      column="ENTREZID",
                                      keytype="ENSEMBL",
                                      multiVals="first")))
  }
  return(list(up=upentrez,
         down=downentrez))
})
names(phenolist) <- unique(names(alltt$monocyte_blood))
phenolist1 <- unlist(phenolist,recursive=F)

ck1 <- compareCluster(geneCluster = phenolist1,
                      fun = "enrichGO",
                      OrgDb='org.Hs.eg.db',
                      ont="BP",
                      minGSSize = 20,
                      maxGSSize = 300,
                      universe=universe)


ckr <- setReadable(ck1, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
#saveRDS(ckr,"output/TWAS_enrichment_overrepClusProfiler_monocyteblood.rds")
ckr <- readRDS("output/TWAS_enrichment_overrepClusProfiler_monocyteblood.rds")

### RRVGO and network plot
ck1df <- as.data.frame(ckr)
ck1df$Cluster <- as.character(ck1df$Cluster)
ck1df_sub <- subset(ck1df, Cluster %in% names(which(table(ck1df$Cluster)>0)))

ck1rvgo <- lapply(unique(ck1df_sub$Cluster), function(x) {
  y <- subset(ck1df_sub, Cluster==x)
  if (nrow(y) < 2) {
    return(y)
  } else {
  simMatrix <- calculateSimMatrix(y$ID,
                                  orgdb="org.Hs.eg.db",
                                  ont="BP",
                                  method="Rel")
  scores <- setNames(-log10(y$pvalue), y$ID)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.8,
                                  orgdb="org.Hs.eg.db")
  keepIDs <- unique(reducedTerms$parentTerm)
  return(subset(y, Description %in% keepIDs))
  }
})
ck1all <- do.call(rbind, ck1rvgo)

ck1all$pheno <- gsub("\\..*","",ck1all$Cluster)
ck1all$updown <- gsub(".*\\.","",ck1all$Cluster)

saveRDS(ck1all,"output/TWAS_enrichment_overrepClusProfiler_monocyteblood_rrvgo_simplified_08.rds")
ck1all <- readRDS("output/TWAS_enrichment_overrepClusProfiler_monocyteblood_rrvgo_simplified_08.rds")

edges_forplot <- network(ck1all[,c("Description","pheno")],directed=F,loops=FALSE)

x <- data.frame(vertices=network.vertex.names(edges_forplot))
set.edge.attribute(edges_forplot,"updown", ck1all$updown)
set.edge.attribute(edges_forplot,"p.adjust", ck1all$p.adjust)

n <- ggnetwork(edges_forplot, 
               layout="fruchtermanreingold",
               niter=5000,
               cool.exp=3,
               cell.jitter=0.5)
n$node_category <- ifelse(n$vertex.names %in% c(ck1all$pheno), "pheno","GO")
n$node_size<- ifelse(n$vertex.names %in% c(ck1all$pheno), 5,1)

GOtolabel <- names(table(ck1all$Description))[which(table(ck1all$Description)>1)]
nameref <- readRDS("output/variable_name_reference.rds")
n$pheno_names <- nameref$varnames[match(n$vertex.names,nameref$mono.variable)]
n$GO_names <- n$vertex.names
n$GO_names[duplicated(n$GO_names)] <- NA

#pdf("paper/figures/TWAS_results_enrichment_hypergeometric_network.pdf",w=12,h=12)
ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(size=-log(p.adjust),col=updown),
             curvature=0,
             show.legend = T) +
  geom_nodes(aes(col=node_category, size=node_size*10)) +
  geom_label(data=subset(n, node_category=="pheno"),aes(label=pheno_names))+
  geom_text_repel(data=subset(n, vertex.names %in% GOtolabel),aes(label=GO_names),size=3)+
  scale_color_jco()+
  scale_size_continuous(range=c(0.1,3))+
  theme_blank()+
  labs(title = "overrepresentation FDR < 0.05 hypergeometric enrichment, BP, rrvgo simplified, 0.8")+
  theme(legend.position = "none")
#dev.off()



######### built-in simplify function...
cks <- clusterProfiler::simplify(ckr)
cksdf <- as.data.frame(cks)

enrichplot::cnetplot(cks,node_labels="category")
enrichplot::dotplot(cks)

#########################################
#########################################
#########################################
###### rank-based GSEA
rankedlist <- lapply(unique(names(alltt$monocyte_blood)), function(pheno) {
  ats <- alltt$monocyte_blood[[pheno]]
  
  ats <- ats[order(ats$t,decreasing=T),]
  
  geneList <- ats$t
  names(geneList) <- ats$gene
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
  
  egor <- setReadable(ego, OrgDb = org.Hs.eg.db, keyType="ENTREZID")
  egors <- simplify(egor)
  
  return(list(raw=egor,
              simplified=egors))

})
names(rankedlist) <- unique(names(alltt$monocyte_blood))
saveRDS(rankedlist, "output/TWAS_enrichment_rankGSEA_monocyteblood.rds")

#### example GSEA plot
toplot <- rankedlist$mf3123$simplified

egodf <- as.data.frame(toplot)
top5cat <- egodf$ID[1:5]
gseaplot2(toplot,geneSetID = top5cat)

##############
#### custom network plot
rankedlist <- readRDS("output/TWAS_enrichment_rankGSEA_monocyteblood.rds")

alldf <- lapply(rankedlist,function(x){
  as.data.frame(x$simplified)
})

alldf <- do.call(rbind,alldf)
alldf$pheno <- gsub("\\..*","",rownames(alldf))

alldfrvgo <- lapply(unique(alldf$pheno), function(x) {
  lapply(list("up","down"), function(updown) {
    if (updown=="up") {
      y <- subset(alldf, pheno==x & NES > 0)
    } else {
      y <- subset(alldf, pheno==x & NES < 0)
    }
    
    if (nrow(y) < 2) {
      return(y)
    } else {
  
  simMatrix <- calculateSimMatrix(y$ID,
                                  orgdb="org.Hs.eg.db",
                                  ont="BP",
                                  method="Rel")
  scores <- setNames(-log10(y$pvalue), y$ID)
  reducedTerms <- reduceSimMatrix(simMatrix,
                                  scores,
                                  threshold=0.8,
                                  orgdb="org.Hs.eg.db")
  
  keepIDs <- unique(reducedTerms$parentTerm)
  return(subset(y, Description %in% keepIDs))
    }
  })
})

dfr1 <- unlist(alldfrvgo,recursive = F)
dfr <- do.call(rbind,dfr1)

#saveRDS(dfr,"output/TWAS_enrichment_rankGSEA_monocyteblood_rrvgo_simplified_08.rds")
dfr <- readRDS("output/TWAS_enrichment_rankGSEA_monocyteblood_rrvgo_simplified_08.rds")

edges_forplot <- network(dfr[,c("Description","pheno")],directed=F,loops=FALSE)

x <- data.frame(vertices=network.vertex.names(edges_forplot))
set.edge.attribute(edges_forplot,"NES", dfr$NES)
set.edge.attribute(edges_forplot,"updown", ifelse(dfr$NES>0,"up","down"))
set.edge.attribute(edges_forplot,"p.adjust", dfr$p.adjust)

n <- ggnetwork(edges_forplot,
               layout="fruchtermanreingold",
               niter=5000,
               cool.exp=3,
               cell.jitter=0.5)

n$node_category <- ifelse(n$vertex.names %in% c(dfr$pheno), "pheno","GO")
n$node_size<- ifelse(n$vertex.names %in% c(dfr$pheno), 5,1)

GOtolabel <- names(table(dfr$Description))[which(table(dfr$Description)>1)]
nameref <- readRDS("output/variable_name_reference.rds")
n$pheno_names <- nameref$varnames[match(n$vertex.names,nameref$mono.variable)]
n$GO_names <- n$vertex.names
n$GO_names[duplicated(n$GO_names)] <- NA

#pdf("paper/figures/TWAS_results_enrichment_GSEA_network.pdf",w=12,h=12)
ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(aes(size=-log(p.adjust),col=updown),
             curvature=0.05,
             show.legend = T) +
  geom_nodes(aes(col=node_category, size=node_size*100)) +
  geom_label(data=subset(n, node_category=="pheno"),aes(label=pheno_names))+
  geom_text_repel(data=subset(n, vertex.names %in% GOtolabel),aes(label=GO_names),size=3)+
  scale_color_jco()+
  scale_size_continuous(range=c(0.1,3))+
  theme_blank()+
  labs(title = "GSEA rank enrichment, BP, rrvgo simplified, 0.8")+
  theme(legend.position = "none")
#dev.off()




testcast <- reshape2::dcast(data=dfr,formula = pheno ~ Description,value.var = "NES")
rownames(testcast) <- testcast$pheno
testcast$pheno <- NULL
testcast[is.na(testcast)] <- 0
testhm <- heatmap(as.matrix(testcast))


dfr_plot <- dfr
dfr_plot$phenof <- factor(dfr_plot$pheno,levels=rownames(testcast)[testhm$rowInd])
dfr_plot$Descf <- factor(dfr_plot$Description,levels=colnames(testcast)[testhm$colInd])

pdf("paper/figures/TWAS_results_enrichment_GSEA_onlytopNES2_heatmap.pdf",w=12,h=12)
subset(dfr_plot, NES > 2 | NES < -2) %>%
ggplot(aes(x=phenof,y=Descf,fill=NES))+
  geom_tile()+
  scale_fill_gradient2_tableau()+
  theme_minimal()+
  theme(axis.text.x=element_text(angle = -45, hjust = 0))
dev.off()
