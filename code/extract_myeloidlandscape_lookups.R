# extract genes for myeloidlandscape
library(SuperExactTest)

### read TWAS results
ttall_unlabelled <- readRDS("output/all_TWAS_results.rds")

### get overlapping sets for different phenotypes:
## cognition
cogsig <- list("MMSE"=ttall_unlabelled$monocyte_blood$cts_mmse30_at_draw$gene[which(ttall_unlabelled$monocyte_blood$cts_mmse30_at_draw$adj.P.Val<0.05)],
               "Global"=ttall_unlabelled$monocyte_blood$cogn_global_at_draw$gene[which(ttall_unlabelled$monocyte_blood$cogn_global_at_draw$adj.P.Val<0.05)],
               "Working memory"=ttall_unlabelled$monocyte_blood$cogn_wo_at_draw$gene[which(ttall_unlabelled$monocyte_blood$cogn_wo_at_draw$adj.P.Val<0.05)],
               "Semantic memory"=ttall_unlabelled$monocyte_blood$cogn_se_at_draw$gene[which(ttall_unlabelled$monocyte_blood$cogn_se_at_draw$adj.P.Val<0.05)],
               "Perceptual speed"=ttall_unlabelled$monocyte_blood$cogn_ps_at_draw$gene[which(ttall_unlabelled$monocyte_blood$cogn_ps_at_draw$adj.P.Val<0.05)],
               "Perceptual orientation"=ttall_unlabelled$monocyte_blood$cogn_po_at_draw$gene[which(ttall_unlabelled$monocyte_blood$cogn_po_at_draw$adj.P.Val<0.05)],
               "Episodic memory"=ttall_unlabelled$monocyte_blood$cogn_ep_at_draw$gene[which(ttall_unlabelled$monocyte_blood$cogn_ep_at_draw$adj.P.Val<0.05)])

cogstest <- supertest(cogsig,n=nrow(ttall_unlabelled$monocyte_blood$mf3123))

## microglia
pamsig <- list(MF=ttall_unlabelled$monocyte_blood$mf3123$hugo[which(ttall_unlabelled$monocyte_blood$mf3123$adj.P.Val<0.05)],
               IT=ttall_unlabelled$monocyte_blood$it3123$hugo[which(ttall_unlabelled$monocyte_blood$it3123$adj.P.Val<0.05)],
               "Post. putamen"=ttall_unlabelled$monocyte_blood$pput3123$hugo[which(ttall_unlabelled$monocyte_blood$pput3123$adj.P.Val<0.05)],
               "VM caudate"=ttall_unlabelled$monocyte_blood$vm3123$hugo[which(ttall_unlabelled$monocyte_blood$vm3123$adj.P.Val<0.05)])

pamstest <- supertest(pamsig,n=nrow(ttall_unlabelled$monocyte_blood$mf3123))
pst <- summary(pamstest)
overlapgenespam <- pst$Table$Elements[which(pst$Table$Degree==max(pst$Table$Degree))]
finalpamoverlap <- data.frame(pheno="PAM",hugo=unlist(strsplit(overlapgenespam,split=", ")))

## amyloid
amyloidsig <- list("Total AB"=ttall_unlabelled$monocyte_blood$amyloid_sqrt$gene[which(ttall_unlabelled$monocyte_blood$amyloid_sqrt$adj.P.Val<0.05)],
                   "Neuritic plaques"=ttall_unlabelled$monocyte_blood$plaq_n_sqrt$gene[which(ttall_unlabelled$monocyte_blood$plaq_n_sqrt$adj.P.Val<0.05)],
                   "Diffuse plaques"=ttall_unlabelled$monocyte_blood$plaq_d_sqrt$gene[which(ttall_unlabelled$monocyte_blood$plaq_d_sqrt$adj.P.Val<0.05)],
                   "Cerebral AA"=ttall_unlabelled$monocyte_blood$caa_4gp$gene[which(ttall_unlabelled$monocyte_blood$caa_4gp$adj.P.Val<0.05)],
                   "Patho AD"=ttall_unlabelled$monocyte_blood$pathoAD$gene[which(ttall_unlabelled$monocyte_blood$pathoAD$adj.P.Val<0.05)])

amyloidstest <- supertest(amyloidsig,n=nrow(ttall_unlabelled$monocyte_blood$mf3123))

ast <- summary(amyloidstest)
overlapgenesamy <- ast$Table$Elements[which(ast$Table$Degree==max(ast$Table$Degree))]
finalamyoverlap <- data.frame(pheno="amyloid",hugo=unlist(strsplit(overlapgenesamy,split=", ")))
                          
## arteriolosclerosis
arteriolgenes =ttall_unlabelled$monocyte_blood$arteriol_scler$gene[which(ttall_unlabelled$monocyte_blood$arteriol_scler$adj.P.Val<0.05)]

############ ALL GENES
allsig <- lapply(ttall_unlabelled$monocyte_blood, function(x) {
  xsub <- subset(x, adj.P.Val < 0.05)
  xsub[,c("pheno","gene","hugo","t","adj.P.Val")]
})
allsig <- do.call(rbind,allsig)
rownames(allsig) <- NULL

write.csv(allsig,file="output/mono_blood_for_myeloidlandscape.csv",row.names = F,quote=F)

