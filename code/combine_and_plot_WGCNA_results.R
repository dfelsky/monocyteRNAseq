### analyze WGCNA output
### assess which set of results has most enrichments and use that in downstream analyses:

enrichfiles <- list.files(pattern = "module_enrichment_results",recursive = T)

allenrich <- lapply(enrichfiles, function(x) {
  y <- readRDS(x)
  y <- subset(y, adj.P.Val < 0.05)
  return(list(file=x,
           uniquemodswithsig=length(unique(y$module)),
           totalenrich=nrow(y),
           pdist=quantile(y$adj.P.Val,seq(0,0.1,0.01))
           ))
})

aer <- do.call(rbind,allenrich)





setwd("/Users/dfelsky/Documents/monocyte_twas/integrated_analyses_v4/WGCNA/")

datalist_mono <- readRDS("mono/datalist.rds")
dat_mono <- t(datalist_mono$all)
datalist_dlpfc <- readRDS("dlpfc/datalist.rds")
dat_dlpfc <- t(datalist_dlpfc$all)

net_mono <- readRDS("mono/net_all_signed_DeepSplit3_power12.rds")
net_dlpfc <- readRDS("dlpfc/net_all_signed_DeepSplit3_power20.rds")

#### module overlaps
mod_overlap <- overlapTableUsingKME(dat1 = dat_mono,
                                    dat2 = dat_dlpfc,
                                    MEs1 = net_mono$MEs,
                                    MEs2 = net_dlpfc$MEs,
                                    name1 = "mono",
                                    name2 = "dlpfc",
                                    colorh1 = net_mono$colors,
                                    colorh2 = net_dlpfc$colors)

overlapYN <- apply(mod_overlap$PvaluesHypergeo,2,function(x) {
  ifelse(x < 0.0005,1,0)
})

library(corrplot)

corrplot(overlapYN,is.corr = F,method="color")
heatmap(overlapYN)
quantile(mono_TOM,seq(0.95,1,by=0.001))


############
library(igraph)

mono_TOM <- readRDS("mono/TOM_power12_withoutsuboutliers_withgenetreeoutliers.rds")
colnames(mono_TOM) <- colnames(dat_mono)
rownames(mono_TOM) <- colnames(dat_mono)
adj <- mono_TOM
adj[adj > 0.08] = 1
adj[adj != 1] = 0
g <- graph.adjacency(adj,mode="undirected",weighted=TRUE,diag=FALSE)
g <- simplify(g,remove.multiple=TRUE,remove.loops=TRUE) 

#E(g)[which(E(g)$weight<0)]$color <- "darkblue"
#E(g)[which(E(g)$weight>0)]$color <- "darkred"

#E(g)$weight <- abs(E(g)$weight)

V(g)$color <- net_mono$colors
# remove unconnected nodes
g <- delete.vertices(g, degree(g)==0)
g <- delete.vertices(g, V(g)$color=="grey")
# Amplify or decrease the width of the edges
# edgeweights <- E(g)$weight * 2.0

# Convert the graph adjacency object into a minimum spanning tree based on Prim's algorithm
#mst <- mst(g, algorithm="prim")

# Plot the tree object
plot(
  g,
  layout=layout.fruchterman.reingold,
  edge.curved=TRUE,
  asp=FALSE,
  edge.width=0.1,
  edge.arrow.mode=0,vertex.label=NA,vertex.size=3)



###################################################
##### module - trait associations

allnetfiles <- list.files(pattern = "net_")

for (netfile in allnetfiles) {
  net <- readRDS(netfile)
  dset <- strsplit(netfile,"_")[[1]][2]
  is.signed <- strsplit(netfile,"_")[[1]][3]
  deepsplit <- gsub(strsplit(netfile,"_")[[1]][4],pattern="DeepSplit",replacement="")
  
  ### trait-ME associations
  ME <- net$MEs
  ME$projid <- rownames(net$MEs)
  
  ### vilas cell types
  library(stringr)
  vcells <- read.csv("/Users/dfelsky/Documents/data/vilas_deconv_2020/deconv_estimates_20200723.csv",as.is=T,header = T)
  rownames(vcells) <- vcells$X
  vcells$X <- NULL
  vcells <- as.data.frame(t(vcells))
  rownames(vcells) <- as.character(gsub(rownames(vcells),pattern="X",replacement=""))
  rownames(vcells) <- str_pad(rownames(vcells),8,pad="0",side="left")
  vcells$projid <- rownames(vcells)
  
  md1 <- merge(ME,ROSmaster,by="projid")
  md <- merge(md1,vcells,by="projid")
  #md <- merge(md,v$targets,by="projid",all.x=T)
  
  library(broom)
  egenelist <- names(net$MEs) ### might want to subset for example to top 20
  cont.phenos <- c("plaq_n_sqrt","plaq_d_sqrt","amyloid_sqrt","tangles_sqrt","nft_sqrt","arteriol_scler","caa_4gp","cvda_4gp2","dlbdx","hspath_any","tdp_stage4","vm3123","pput3123","it3123","mf3123") #names(vcells)[-48] #
  bin.phenos <- c("smoking","ad_reagan","ci_num2_gct","ci_num2_mct") # binary phenotypes
  covars <- c("age_death","pmi","msex")
  
  
  index <- 1
  pvalues <- NULL
  bvalues <- NULL
  nvalues <- NULL
  phenovalues <- NULL
  egenevalues <- NULL
  for (pheno in cont.phenos) {
    for (egene in egenelist) {
      form <- formula(paste(pheno,"~",egene,"+",paste0(covars,collapse = " + ")))
      mod <- lm(data=md, form)
      pvalues[index] <- tidy(mod)$p.value[2]
      bvalues[index] <- coef(mod)[2]
      nvalues[index] <- dim(mod$model)[1]
      egenevalues[index] <- egene
      phenovalues[index] <- pheno
      index <- index + 1
    }
  }
  
  results <- data.frame(pheno=phenovalues,
                        snp=egenevalues,
                        b=bvalues,
                        p=pvalues,
                        n=nvalues)
  
  bonfT <- 0.05/nrow(results)
  
  pdf(file=paste0("ME_pathology_associations_",dset,"_",is.signed,"_deepSplit",deepsplit,".pdf"))
  print(ggplot(data=results, aes(x=pheno,y=-log10(p)*sign(b),group=snp))+
          geom_bar(stat="identity",position="dodge")+
          geom_hline(yintercept = c(-log10(bonfT),log10(bonfT)),col="red",lty=2)+ 
          geom_hline(yintercept = c(-log10(0.05),log10(0.05)),col="blue",lty=3)+ 
          geom_text(data=subset(results,p<0.05),aes(label=snp))+
          theme_minimal()+
          theme(axis.text.x=element_text(angle = -45, hjust = 0)))
  dev.off()
  
}


