# Figure 4
net_mono <- readRDS("output/WGCNA/mono/net.rds")
net_dlpfc <- readRDS("output/WGCNA/dlpfc/net_noPAM_celltypes.rds")

tiff("paper/figures/Figure4_dendrogram_mono_noPAM.tif",w=4,h=2,units="in",res=600)
plotDendroAndColors(net_mono$dendrograms[[1]],
                    net_mono$colors,
                    dendroLabels = FALSE,
                    hang = 0.02,
                    addGuide = F,
                    guideHang = 0.05,
                    main = "Monocytes")
dev.off()


net_mono <- readRDS("output/WGCNA/mono/net_noPAM.rds")

tiff("paper/figures/Figure4_dendrogram_dlpfc_noPAM.tif",w=4,h=2,units="in",res=600)
plotDendroAndColors(net_dlpfc$dendrograms[[1]],
                    net_dlpfc$colors,
                    dendroLabels = FALSE,
                    hang = 0.02,
                    addGuide = F,
                    guideHang = 0.05,
                    main = "DLPFC")
dev.off()


#### plot for significant gene effects distribution in modules
moddist_siggenes <- lapply(ttall_unlabelled$monocyte_blood[c(2,8,9,16:19)],function(pheno) {
  sapply(c("up","down"), function(updown) {
    if (updown=="up") {
      sgenes <- pheno$hugo[which(pheno$adj.P.Val < 0.05 & pheno$t > 0 )]
      table(net_mono$colors[match(sgenes,names(net_mono$colors))])
    } else {
  sgenes <- pheno$hugo[which(pheno$adj.P.Val < 0.05 & pheno$t < 0)]
  table(net_mono$colors[match(sgenes,names(net_mono$colors))])
    }
  })
})






