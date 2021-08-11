### Network plot of monocyte modules
library(GGally)
library(network)
library(sna)
library(ggplot2)
library(ggnetwork)
library(reshape2)
library(stringr)
library(ggsci)
library(gplots)
library(Hmisc)

load("output/WGCNA/dlpfc/dlpfc_signed_signed_minKMEtoStay03_noPAM_celltypes-block.1.RData")
net <- readRDS("output/WGCNA/dlpfc/net_noPAM_celltypes.rds")

TOM <- as.matrix(TOM)
dimnames(TOM) <- list(names(net$colors),names(net$colors))

netdf <- melt(lower.tri.remove(TOM)) # matrix is symmetrical, we designate half of it at NA with lower.tri.remove
names(netdf)[1:2] <- c("gene1","gene2")
netdf <- netdf[complete.cases(netdf),]

netdf_sub <- subset(netdf, value > 0 & gene1!=gene2)
netdf_sub <- subset(netdf_sub, value > quantile(netdf_sub$value,0.99))

## remove grey
netdf_sub <- subset(netdf_sub, gene1 %nin% names(net$colors)[which(net$colors=="grey")] & gene2 %nin% names(net$colors)[which(net$colors=="grey")])

edges <- netdf_sub[,c("gene1","gene2")]
edges$gene1 <- as.character(edges$gene1)
edges$gene2 <- as.character(edges$gene2)

edges_forplot <- network(edges,directed=F,loops=FALSE)
x <- data.frame(vertices=network.vertex.names(edges_forplot))
set.edge.attribute(edges_forplot,"similarity", netdf_sub$value)

n <- ggnetwork(edges_forplot)
n$module <- net$colors[match(n$vertex.names,names(net$colors))]

ucols <- unique(n$module)
uhex <-  col2hex(ucols)
colpal <- uhex
names(colpal) <- ucols

ggplot(n, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(size=0.05,curvature=0,col="grey85",show.legend = F) +
  geom_nodes(aes(colour=module)) +
  scale_size_continuous(range=c(0.1,3))+
  scale_colour_manual(values=colpal)+
  theme_blank()


saveRDS(list(n=n,colpal=colpal),file="output/WGCNA/dlpfc/networkplot_objects.rds")