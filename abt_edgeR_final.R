library(edgeR)
library(ggplot2)
library(gplots)
library(reshape)
library(eulerr)
library(Vennerable)
library(venn)
library(Glimma)
library(AnnotationDbi)
library(org.Mm.eg.db)
library(annotate)
library(lattice)
library(ggpubr)
library(ggsignif)
library(FactoMineR)
library(fifer)
library(GOSemSim)

####Prepare data########################################################
tissue = "Liver"

workdir = paste0("/home/labs/amit/tomerlan/ABT_project/", tissue)

if(tissue == "Skin"){
  rawdata = read.delim(paste0("/home/labs/amit/tomerlan/ABT_project/data_files/raw_data_skin_Symbol.txt"), check.names = F, stringsAsFactors = F)[, -c(2:9)]
  colnames(rawdata) = c("Symbol", "ABT1", "ABT3", "ABT4", "ABT5", "DM1", "DM3", "DM4", "DM5", "Y1", "Y2", "Y3", "Y4")
  rawdata = rawdata[, -13] #remove Y4 which is off
  rawdata = rawdata[, c(1, 10, 11, 12, 6, 7, 8, 9, 2, 3, 4, 5)] #reorder conditions
  rawdata = rawdata[, -c(5, 9)] #remove possibly swaped samples
  groups = factor(c(rep("Y", 3), rep("DM", 3), rep("ABT", 3)), levels = c("Y", "DM", "ABT"))
  y = DGEList(counts = rawdata[, -1], genes = rawdata[, 1], group = groups) #make DGE object
  y$samples$group = relevel(y$samples$group, ref = "Y")
  idfound = y$genes$genes %in% mappedRkeys(org.Mm.egSYMBOL)
  y = y[idfound,]
  egSYMBOL = toTable(org.Mm.egSYMBOL)
  m = match(y$genes$genes, egSYMBOL$symbol)
  y$genes$EntrezGene = egSYMBOL$gene_id[m] ##add ENTREZ annotation
} else {
  rawdata = read.delim(paste0("/home/labs/amit/tomerlan/ABT_project/data_files/raw_data_RefSeq.txt"), check.names = F, stringsAsFactors = F)[, -c(2:9)]
  rawdata = rawdata[, c(1, grep(tissue, colnames(rawdata)))]
  colnames(rawdata) = c("RefSeq", "ABT1", "ABT3", "ABT4", "ABT5", "DM1", "DM3", "DM4", "DM5", "Y1", "Y2", "Y3", "Y4")
  rawdata = rawdata[, c(1, 10, 11, 12, 13, 6, 7, 8, 9, 2, 3, 4, 5)] #reorder conditions
  rawdata = rawdata[, -c(6, 10)] #remove possibly swaped samples
  groups = factor(c(rep("Y", 4), rep("DM", 3), rep("ABT", 3)), levels = c("Y", "DM", "ABT"))
  y = DGEList(counts = rawdata[, -1], genes = rawdata[, 1], group = groups)
  y$samples$group = relevel(y$samples$group, ref = "Y")
  idfound = y$genes$genes %in% mappedRkeys(org.Mm.egREFSEQ)
  y = y[idfound,]
  egREFSEQ = toTable(org.Mm.egREFSEQ)
  m = match(y$genes$genes, egREFSEQ$accession)
  y$genes$EntrezGene = egREFSEQ$gene_id[m] ##add ENTREZ annotation
  egSYMBOL = toTable(org.Mm.egSYMBOL)
  m = match(y$genes$EntrezGene, egSYMBOL$gene_id)
  y$genes$genes = egSYMBOL$symbol[m]  ##replace RefSeq with SYMBOL annotation
}


####Filtering and normalization################################################
#Different RefSeq transcripts for the same gene symbol count predominantly the same reads. So we keep one transcript for each gene symbol. We choose the transcript with highest overall count:
o = order(rowSums(y$counts), decreasing = TRUE) #order in desending 
y = y[o,]
d = duplicated(y$genes$genes)
y = y[!d,]
nrow(y)

hist(log10(cpm(y)), breaks = 1000)
sample_cutoff = 6
abline(v = log10(sample_cutoff), lty = 2, col = 2)
keep = rowSums(cpm(y) > sample_cutoff) >= 2 #find good genes
table(keep)
y = y[keep, keep.lib.sizes = FALSE] #keep only good genes

y$samples$lib.size = colSums(y$counts) #Recompute the library sizes
y = calcNormFactors(y) #TMM normalization is applied to this dataset to account for compositional difference between the libraries

rownames(y$counts) <- rownames(y$genes) <- y$genes$EntrezGene #Use Entrez Gene IDs as row names
# 
# y$genes$EntrezGene <- NULL #remove Entrez Gene IDs as column

plotMD(cpm(y, log = TRUE), column = 1)
abline(h = 0, col = "red", lty = 2, lwd = 2)

y$samples


####Data exploration######################################################
colors = rep(c("blue", "red", "green"), 2)
plotMDS(y, col = colors[groups])
dev.copy(pdf, file = paste0(workdir, "/MDS_",tissue, ".pdf"))
dev.off()


####DE analysis classic approach###########################################################
alpha = 0.05

y = estimateDisp(y, robust = F)
y$common.dispersion
plotBCV(y)

YvsDM = exactTest(y, pair = c("Y", "DM"))
summary(decideTests(YvsDM))

YvsABT = exactTest(y, pair = c("Y", "ABT"))
summary(decideTests(YvsABT))

DMvsABT = exactTest(y, pair = c("DM", "ABT"))
summary(decideTests(DMvsABT))

save(YvsDM, file = paste0(workdir, "/YvsDM.RData"))
save(YvsABT, file = paste0(workdir, "/YvsABT.RData"))
save(DMvsABT, file = paste0(workdir, "/DMvsABT.RData"))

cor = cor.test(x = YvsDM$table[,"logFC"], y = DMvsABT$table[,"logFC"], method = "spearman")$estimate


####Senescence signature#####################################
sasp2 = c("Mmp1a", "Mmp3", "Mmp10", "Mmp12", "Mmp13", "Mmp14",
         "Timp1", "Serpine1", "Serpinb2", "Plat", "Plau", "Ctsb",
         "Icam1", "Icam5", "Tnfrsf11b", "Tnfrsf1a", "Tnfrsf23",
         "Tnfrsf22", "Fas", "Tnfrsf1b", "Plaur", "Il6st", "Egfr", "Fn1",
         "Areg", "Ereg", "Egf", "Fgf2", "Hgf","Fgf7", "Vegfa", "Ang",
         "Kitlg", "Pigf", "Igfbp2", "Igfbp3", "Igfbp4", "Igfbp5", "Csf2", "Mif", 
         "Igfbp7", "Il1", "Il6", "Il7", "Il13", "Il15", "Cxcl1", "Cxcl2", "Cxcl12", "Ccl8", "Ccl3", "Ccl20",
         "Ccl16-ps", "Ccl26")

sasp = c("Mmp1a", "Mmp3", "Mmp10", "Mmp12", "Mmp13", "Mmp14",
          "Timp1", "Serpine1", "Serpinb2", "Plat", "Plau", "Ctsb",
          "Icam1", "Icam5", "Tnfrsf11b", "Tnfrsf1a", "Tnfrsf23",
          "Tnfrsf22", "Fas", "Tnfrsf1b", "Plaur", "Il6st", "Egfr", "Fn1",
          "Areg", "Ereg", "Egf", "Fgf2", "Hgf","Fgf7", "Vegfa", "Ang",
          "Kitlg", "Pigf", "Igfbp2", "Igfbp3", "Igfbp4", "Igfbp5", "Csf2", "Mif", 
          "Igfbp7",  "Il1a", "Il1b", "Il6", "Il7", "Il13", "Il15", 
          "Cxcl1", "Cxcl2", "Cxcl12", "Cxcl15", "Ccl8", "Ccl3", "Ccl20", "Ccl16-ps", "Ccl26")

sasp = rownames(y[which(y$genes[, "genes"] %in% sasp),]$genes)

sasp_FC = data.frame(gene = y$gene[sasp, "genes"], UT = YvsDM$table[sasp, "logFC"], UT_pvalue = YvsDM$table[sasp, "PValue"], T = YvsABT$table[sasp, "logFC"], T_pvalue = YvsABT$table[sasp, "PValue"], UTT_pvalue = DMvsABT$table[sasp, "PValue"])
# [DMvsABT$table[sasp, "PValue"] <= alpha,]
sasp_score = c(mean(abs(YvsDM$table[sasp, "logFC"])), mean(abs(YvsABT$table[sasp, "logFC"])))
sasp_t = t.test(sasp_FC[,"UT"], sasp_FC[, "T"], paired = T)$p.value
# 
# ggplot(melt(sasp_FC[, -c(3, 5, 6)]), aes(x = factor(gene, levels = sasp_FC$gene[order(YvsDM$table[sasp, "logFC"])]), y = value, fill = factor(variable))) + 
#   geom_bar(position = "dodge", stat = "identity") + 
#   xlab("FC over young") + 
#   ylab("SASP") + 
#   ggtitle(tissue) +
#   guides(fill = guide_legend("Condition")) +
#   coord_flip() +
#   ggsave(file = paste0(workdir, "/sasp_FC_", tissue,".pdf"), width = 20, height = 30, units = "cm")

write.csv(sasp_FC, file = paste0(workdir, "/sasp_FC_", tissue, ".csv"))
write.csv(c(sasp_score, sasp_t, cor), file = paste0(workdir, "/sasp_score_", tissue,".csv"))


#########################
sen = c("Trp53bp1", "Trp53", "Cdkn2a", "Cdkn2b", "Cdkn1c", "Cdkn1a", "Hmga1", "Hmgb1", 
        "Lmnb1", "Tnf", "Il1a", "Il1b", "Il10", "Ccl5", 
        "Cxcl10",  "Cxcl9", "Bhlhe40", "Notch3", "Dpp4")
sen = rownames(y[which(y$genes[, "genes"] %in% sen),]$genes)

sen_FC = data.frame(gene = y$gene[sen, "genes"], UT = YvsDM$table[sen, "logFC"], UT_pvalue = YvsDM$table[sen, "PValue"], T = YvsABT$table[sen, "logFC"], T_pvalue = YvsABT$table[sen, "PValue"], UTT_pvalue = DMvsABT$table[sen, "PValue"])
# [DMvsABT$table[sen, "PValue"] <= alpha,]
sen_score = c(mean(abs(YvsDM$table[sen, "logFC"])), mean(abs(YvsABT$table[sen, "logFC"])))
sen_t = t.test(sen_FC[,"UT"], sen_FC[, "T"], paired = T)$p.value

# ggplot(melt(sen_FC[, -c(3, 5, 6)]), aes(x = factor(gene, levels = sen_FC$gene[order(YvsDM$table[sen, "logFC"])]), y = value, fill = factor(variable))) + 
#   geom_bar(position = "dodge", stat = "identity") + 
#   xlab("FC over young") + 
#   ylab("SASP") + 
#   ggtitle(tissue) +
#   guides(fill = guide_legend("Condition")) +
#   coord_flip() +
#   ggsave(file = paste0(workdir, "/sen_FC_", tissue,".pdf"), width = 20, height = 30, units = "cm")

write.csv(sen_FC, file = paste0(workdir, "/sen_FC_", tissue,".csv"))
write.csv(c(sen_score, sen_t), file = paste0(workdir, "/sen_score_", tissue,".csv"))


######GO and KEGG analysis
go = goana(DMvsABT, species = "Mm", FDR = alpha)

down = topGO(go, ontology = NULL, number = 2000)[which(topGO(go, ontology = NULL, number = 2000)[, 7] < alpha),]
write.csv(down, file = paste0(workdir, "/go_down_", tissue,".csv"))
up = topGO(go, ontology = NULL, number = 2000)[which(topGO(go, ontology = NULL, number = 2000)[, 6] < alpha),]
write.csv(up, file = paste0(workdir, "/go_up_", tissue,".csv"))

kegg = kegga(DMvsABT, species = "Mm", plot = T); keggtogene = getGeneKEGGLinks("mmu")

KEGG_down = topKEGG(kegg, number = 2000)[which(topKEGG(kegg, number = 2000)[, 6] < alpha),]
write.csv(KEGG_down, file = paste0(workdir, "/kegg_down_", tissue,".csv"))
KEGG_up = topKEGG(kegg, number = 2000)[which(topKEGG(kegg, number = 2000)[, 5] < alpha),]
write.csv(KEGG_up, file = paste0(workdir, "/kegg_up_", tissue,".csv"))


####All tissues#####################################
workdir2 = "/home/labs/amit/tomerlan/ABT_project/all_tissues"
workdir3 = "/home/labs/amit/tomerlan/ABT_project"

tissues = c("Kidney", "Liver", "Lung", "Skin")

alpha = 0.05

##FC heatmap
all_genes = table(unlist(apply(as.matrix(paste0("/home/labs/amit/tomerlan/ABT_project/", tissues, "/DMvsABT.RData")), 1, 
                               function(n) get(load(n))$genes$EntrezGene)))

common_genes = names(all_genes[which(all_genes >= (length(tissues) - 3))])

all_FC = get(load(paste0("/home/labs/amit/tomerlan/ABT_project/", tissues[1], "/DMvsABT.RData")))$table[common_genes, "logFC"]
apply(as.matrix(paste0("/home/labs/amit/tomerlan/ABT_project/", tissues[-1], "/DMvsABT.RData")), 1, 
      function(n) all_FC <<- cbind(all_FC, get(load(n))$table[common_genes, "logFC"]))
colnames(all_FC) = tissues
rownames(all_FC) = common_genes

all_affect_genes = unique(unlist(apply(as.matrix(paste0("/home/labs/amit/tomerlan/ABT_project/", tissues, "/DMvsABT.RData")), 1, 
                                       function(n) rownames(topTags(get(load(n)), sort.by = "logFC", p.value = alpha, n = 3000)))))

all_FC = all_FC[intersect(all_affect_genes, common_genes),]
all_FC[which(is.na(all_FC))] = 0

clustering_method = "kmeans"
FC_clusts = 2
# all_FC_clusters = cutree(hclust(dist(all_FC), method = "ward.D2"), FC_clusts)
all_FC_clusters = kmeans(x = all_FC, centers = FC_clusts)$cluster
all_FC = all_FC[names(sort(all_FC_clusters)),]

all_FC_melted = melt(as.matrix(all_FC)); all_FC_melted$X1 = factor(all_FC_melted$X1, levels = rownames(all_FC))

cluster_lines = data.frame(y1 = cumsum(table(sort(all_FC_clusters))), y2 = cumsum(table(sort(all_FC_clusters))), x1 = 0.5, x2 = 4.5)

ggplot(all_FC_melted, aes(x = X2, y = X1, fill = value)) +
  geom_tile() +
  scale_fill_gradientn(colours = colorRampPalette(c("royalblue3", "white", "firebrick3"))(1000), values = scales::rescale(c(min(as.matrix(all_FC), na.rm = T), min(as.matrix(all_FC), na.rm = T)/3 ,  0,  max(as.matrix(all_FC)/3, na.rm = T), max(as.matrix(all_FC), na.rm = T)))) +
  labs(x = NULL, y = NULL, fill = "log(Vehicle/ABT-737)") +
  # guides(color = guide_legend("my title")) +
  theme(legend.title = element_text(size = 4),
        legend.text = element_text(size = 4),
        axis.line.x = element_blank(),
        axis.line.y = element_blank(),
        axis.text.x = element_blank(),
        axis.text.y = element_blank(),
        axis.ticks.x = element_blank(),
        axis.ticks.y = element_blank(),
        plot.title = element_text(size = 4),
        axis.title = element_text(size = 4, face = "bold"),
        legend.position = "bottom",
        legend.key.size = unit(x = 0.2, unit = "cm")) +
  geom_segment(data = cluster_lines, aes(y = y1, yend = y2, x = x1, xend = x2), color = "black", size = 0.4, inherit.aes = F) +
  ggsave(file = paste0(workdir2, "/heatmap_", clustering_method, ".pdf"), width = 6, height = 10, units = "cm")

# annot_clust = data.frame(group = factor(all_FC_clusters, labels = unique(all_FC_clusters)))
# annot_clust_colors = list(group = brewer.pal(FC_clusts, "Set3"))
# names(annot_clust_colors$group) = unique(all_FC_clusters)
# pheatmap(all_FC, border_color = NA, col = colorRampPalette(c("royalblue3", "white", "firebrick3"))(1000), cluster_cols = F, cluster_rows = F, clustering_method = "ward.D2", cutree_rows = FC_clusts, treeheight_row = 20, show_rownames = F, annotation_row = annot_clust, annotation_colors = annot_clust_colors)
# dev.copy(pdf, file = paste0(workdir2, "/FC_heatmap3.pdf"))
# dev.off()


##GO analysis for FC heatmap
hsGO = godata('org.Mm.eg.db', ont = "BP")

FC_go = list()
for(n in 1:FC_clusts) {FC_go[[n]] = goana(common_genes[which(common_genes %in% names(which(all_FC_clusters == n)))], universe = common_genes, species = "Mm", FDR = alpha)}

FC_go_sig = list()
# for(n in 1:FC_clusts) {FC_go_sig[[n]] = rownames(FC_go[[n]])[which(FC_go[[n]][, "P.DE"] < alpha)]}
for(n in 1:FC_clusts) {FC_go_sig[[n]] = topGO(FC_go[[n]], ontology = NULL, number = 2000)[which(topGO(FC_go[[n]], ontology = NULL, number = 2000)[, "P.DE"] < alpha),]}

go_table = list()
for(n in 1:2) {go_table[[n]] = cbind(rep.int(n, times = c(length(FC_go_sig[[n]][, 2]))), FC_go_sig[[n]])}

go_table = do.call(rbind, go_table)

write.csv(go_table, file = paste0(workdir2, "/go_", clustering_method, ".csv"))

# FC_go_sim = list()
# for(n in 1:FC_clusts) {FC_go_sim[[n]] = mgoSim(FC_go_sig[[n]], FC_go_sig[[n]], semData = hsGO, measure = "Wang", combine = NULL)}
# 
# GO_clusts = 50
# FC_go_sim_clust = list()
# for(n in c(1:FC_clusts)) {if(length(colnames(FC_go_sim[[n]])) >= GO_clusts) {FC_go_sim_clust[[n]] = cutree(hclust(dist(FC_go_sim[[n]]), method = "ward.D2"), k = GO_clusts)} else {FC_go_sim_clust[[n]] = colnames(FC_go_sim[[n]])}}
# # FC_go_sim_clust[[n]] = cutree(hclust(dist(FC_go_sim[[n]]), method = "ward.D2"), k = GO_clusts)
# 
# FC_go_sim_top = list()
# for(n in 1:FC_clusts) {FC_go_sim_top[[n]] = apply(as.matrix(1:GO_clusts), 1, function(m) FC_go[[n]][names(which(FC_go_sim_clust[[n]] == m)),][which.min(FC_go[[n]][names(which(FC_go_sim_clust[[n]] == m)), "P.DE"]), c("Term", "P.DE")])}
# for(n in 1:FC_clusts) {FC_go_sim_top[[n]] = do.call(rbind, FC_go_sim_top[[n]])}
# go_table = do.call(rbind, FC_go_sim_top)
# g= c(); for(n in 1:FC_clusts) {g = c(g, length(FC_go_sim_top[[n]][, 1]))}
# go_table = cbind(rep.int(1:FC_clusts, times = g), go_table)
# colnames(go_table) = c("Cluster", "GO Tern", "P-Value")
# go_table[,"P-Value"] = signif(go_table[,"P-Value"], 2)
# write.csv(go_table, file = paste0(workdir2, "/go_", clustering_method, "_" , FC_clusts, ".csv"))


##KEGG analysis for FC heatmap
FC_kegg = list()
for(n in 1:FC_clusts) {FC_kegg[[n]] = kegga(common_genes[which(common_genes %in% names(which(all_FC_clusters == n)))], universe = common_genes, species = "Mm", FDR = alpha)}

FC_kegg_sig = list()
# for(n in 1:FC_clusts) {FC_kegg_sig[[n]] = rownames(FC_kegg[[n]])[which(FC_kegg[[n]][, "P.DE"] < alpha)]}
for(n in 1:FC_clusts) {FC_kegg_sig[[n]] = topKEGG(FC_kegg[[n]], number = 2000)[which(topKEGG(FC_kegg[[n]], number = 2000)[, "P.DE"] < alpha),]}

kegg_table = list()
for(n in 1:2) {kegg_table[[n]] = cbind(rep.int(n, times = c(length(FC_kegg_sig[[n]][,2]))), FC_kegg_sig[[n]])}

kegg_table = do.call(rbind, kegg_table)
write.csv(kegg_table, file = paste0(workdir2, "/kegg_", clustering_method, ".csv"))


###########SASP score
###main figure
scores = data.frame(t(apply(as.matrix(tissues), 1, function(n) read.csv(paste0("/home/labs/amit/tomerlan/ABT_project/",n,"/sasp_score_",n,".csv"))[,"x"])))
scores = cbind(tissues, scores)
colnames(scores) = c("Tissue", "Vehicle", "ABT-737", "p-value", "cor")

write.csv(scores, file = paste0(workdir2, "/scores.csv"))

scores[,"Vehicle"] = scores[,"Vehicle"] - scores[,"ABT-737"]

ggplot(melt(scores, id.vars = c("Tissue", "p-value")), aes(x = factor(Tissue, levels = tissues), y = value, fill = factor(variable))) +
  geom_bar(stat = "identity") +
  xlab("Tissue") +
  ylab("Geometric mean of EG/Young") +
  ggtitle("SASP score") +
  scale_fill_manual(values = c("cornflowerblue", "orange")) +
  guides(fill = guide_legend("Experimental group")) +
  annotate('text', x = c(1,2,3,4), y = scores[,"Vehicle"] + scores[,"ABT-737"] + 0.1, label = c("**", "*", "NS", "***" )) +
  ggsave(file = paste0(workdir2, "/sasp_score.pdf"))

sasp_gg = apply(as.matrix(tissues), 1, function(n) data.frame(data.frame(read.csv(paste0("/home/labs/amit/tomerlan/ABT_project/", n, "/sasp_FC_", n, ".csv"), sep = ",", header = T))[,-1], tissue = rep(n, nrow(data.frame(read.csv(paste0("/home/labs/amit/tomerlan/ABT_project/", n, "/sasp_FC_", n, ".csv"), sep = ",", header = T))))))
# good_sasp = as.character(unique(unlist(apply(as.matrix(1:4), 1, function(n) sasp_gg[[n]][which(sasp_gg[[n]][,"UTT_pvalue"] < alpha), "gene"]))))
all_sasp = as.character(unique(unlist(apply(as.matrix(1:4), 1, function(n) sasp_gg[[n]][, "gene"]))))
for(n in 1:4) {
  rownames(sasp_gg[[n]]) = sasp_gg[[n]][, "gene"]
}

sasp_gg_all = list()
for(m in 1:4){
  sasp_gg_all[[m]] = data.frame(matrix(ncol = 3, nrow = length(all_sasp)))
  rownames(sasp_gg_all[[m]]) = all_sasp
  colnames(sasp_gg_all[[m]]) = c("Vehicle", "ABT-737", "p-value")
  for(n in all_sasp){
    if(n %in% rownames(sasp_gg[[m]])){
      sasp_gg_all[[m]][n,] = sasp_gg[[m]][n, c(2, 4, 6)]
    }
  }
}
sasp_all = do.call(cbind, sasp_gg_all)
write.csv(sasp_all, file = paste0(workdir2, "/sasp_all.csv"))


sasp_gg = do.call(rbind, sasp_gg_all)
sasp_gg = cbind(sasp_gg, c(rep(tissues[[1]], 40), rep(tissues[[2]], 40), rep(tissues[[3]], 40), rep(tissues[[4]], 40)), rep(rownames(sasp_gg_all[[1]]), 4))
sasp_gg = sasp_gg[,-3]
colnames(sasp_gg) = c("Vehicle", "ABT-737", "Tissue", "Gene")
sasp_gg = melt(sasp_gg, id.vars = c("Tissue","Gene"))

ggplot(sasp_gg, aes(x = Gene, y = value, fill = variable)) +
  geom_bar(position = "dodge", stat = "identity") +
  theme_light() +
  scale_fill_manual(values = c("cornflowerblue", "orange")) +
  xlab("Gene") +
  ylab("logFC(EG/Young)") +
  theme(axis.line.x = element_line(size = 0.1),
        axis.line.y = element_line(size = 0.1),
        axis.ticks.x = element_line(size = 0.1),
        axis.ticks.y = element_line(size = 0.1),
        text = element_text(size = 4),
        axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
        axis.text.y = element_text(size = 4),
        strip.background = element_blank(),
        strip.text.x = element_blank(),
        legend.key.size = unit(1, "mm")) +
    coord_flip() +
    # geom_signif(comparisons = list(c("ABT-737", "Vehicle")),manual = T, inherit.aes = F, aes(group = variable), map_signif_level = F, annotation = "1", textsize = 2, margin_top = 0.4) +
    # geom_signif(annotation = c("1","**"), y_position=c(1.8, 6.8), xmin=c(0.875, 1.875), xmax=c(1.125, 2.125), tip_length = 0.01) +
  # geom_signif(data=data.frame(x=c(0.875, 1.875), xend=c(1.125, 2.125),y=c(5.8, 8.5), annotation=c("**", "NS")),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
  facet_wrap(~Tissue, nrow = 1, ncol = 4) +
  guides(fill = guide_legend(title = "Experimental group")) +
  ggsave(file = paste0(workdir2, "/sasp_sen.pdf"))


# for(n in 1:4) {
#   rownames(sasp_gg[[n]]) = sasp_gg[[n]][, "gene"]
#   sasp_gg[[n]] = sasp_gg[[n]][, c(2, 4, 6)]
#   sasp_gg[[n]] = sasp_gg[[n]][all_sasp,]
#   rownames(sasp_gg[[n]]) = all_sasp
#   colnames(sasp_gg[[n]]) = c("Vehicle", "ABT-737", "p-value")
#   sasp_gg[[n]][is.na(sasp_gg[[n]][,"p-value"]), "p-value"] = 1
# }

# sasp_all = do.call(cbind, sasp_gg)
# write.csv(sasp_all, file = paste0(workdir2, "/sasp_all.csv"))

#####supp figures
# k = 4
# sasp_heat = sasp_gg[[k]]
# for(n in 1:nrow(sasp_gg[[k]])){
#   if(sasp_gg[[k]][n, "p-value"] < 0.05 & sasp_gg[[k]][n,"p-value"] > 0.005){sasp_heat[n,"p-value"] = "*"} 
#   else if(sasp_gg[[k]][n, "p-value"] < 0.005 & sasp_gg[[k]][n,"p-value"] > 0.0005){sasp_heat[n,"p-value"] = "**"}
#   else if(sasp_gg[[k]][n, "p-value"] < 0.0005){sasp_heat[n,"p-value"] = "***"}
#   else {sasp_heat[n,"p-value"] = "NS"}
#  }
# 
# zlim = c(-max(unlist(sasp_gg), na.rm = T),  max(unlist(sasp_gg), na.rm = T))
# 
# pdf(paste0(workdir2, "/heatmap_Skin.pdf"), width = 3, height = 7.5)
# par(mar = c(5, 5, 5, 5))
# image(1:ncol(sasp_heat[,c("Vehicle", "ABT-737")]), 1:nrow(sasp_heat[,c("Vehicle", "ABT-737")]), t(sasp_heat[,c("Vehicle", "ABT-737")]), col =  colorRampPalette(c("royalblue3", "white", "firebrick3"))(1000), zlim = zlim, axes = F, xlab = NA, ylab = NA)
# axis(1, 1:ncol(sasp_heat[,c("Vehicle", "ABT-737")]), colnames(sasp_heat[,c("Vehicle", "ABT-737")]), par(las = 2), lwd = 0, lwd.ticks = 1)
# # axis(2, 1:nrow(sasp_heat[,c("Vehicle", "ABT-737")]), rownames(sasp_heat[,c("Vehicle", "ABT-737")]), par(las = 2), lwd = 0, lwd.ticks = 1)
# axis(4, 1:nrow(sasp_heat[,c("Vehicle", "ABT-737")]), sasp_heat[,"p-value"], lwd = 0, lwd.ticks = 1)
# for (z in 1:ncol(sasp_heat[,c("Vehicle", "ABT-737")])){for(y in 1:nrow(sasp_heat[,c("Vehicle", "ABT-737")])){text(z, y, signif(sasp_heat[, c("Vehicle", "ABT-737")][y, z], 2))}}
# box(which = "plot", lty = "solid")
# dev.off()
# 
# 
# ##
# sen_gg = apply(as.matrix(tissues), 1, function(n) data.frame(data.frame(read.csv(paste0("/home/labs/amit/tomerlan/ABT_project/", n, "/sen_FC_", n, ".csv"), sep = ",", header = T))[,-1], tissue = rep(n, nrow(data.frame(read.csv(paste0("/home/labs/amit/tomerlan/ABT_project/", n, "/sen_FC_", n, ".csv"), sep = ",", header = T))))))
# good_sen = as.character(unique(unlist(apply(as.matrix(1:4), 1, function(n) sen_gg[[n]][which(sen_gg[[n]][,"UTT_pvalue"] < alpha), "gene"]))))
# all_sen = as.character(unique(unlist(apply(as.matrix(1:4), 1, function(n) sen_gg[[n]][, "gene"]))))
# 
# for(n in 1:4) {rownames(sen_gg[[n]]) = sen_gg[[n]][, "gene"]
#   sen_gg[[n]] = sen_gg[[n]][, c(2, 4, 6)]
#   sen_gg[[n]] = sen_gg[[n]][all_sen,]
#   rownames(sen_gg[[n]]) = all_sen
#   colnames(sen_gg[[n]]) = c("Vehicle", "ABT-737", "p-value")
#   sen_gg[[n]][is.na(sen_gg[[n]][,"p-value"]), "p-value"] = 1
# }
# 
# sen_all = do.call(cbind, sen_gg)
# write.csv(sen_all, file = paste0(workdir2, "/sen_all.csv"))
# 
# k = 4
# sen_heat = sen_gg[[k]]
# for(n in 1:nrow(sen_gg[[k]])){
#   if(sen_gg[[k]][n, "p-value"] < 0.05 & sen_gg[[k]][n,"p-value"] > 0.005){sen_heat[n,"p-value"] = "*"} 
#   else if(sen_gg[[k]][n, "p-value"] < 0.005 & sen_gg[[k]][n,"p-value"] > 0.0005){sen_heat[n,"p-value"] = "**"}
#   else if(sen_gg[[k]][n, "p-value"] < 0.0005){sen_heat[n,"p-value"] = "***"}
#   else {sen_heat[n,"p-value"] = "NS"}
# }
# 
# zlim = c(-max(unlist(sen_gg), na.rm = T),  max(unlist(sen_gg), na.rm = T))
# 
# pdf(paste0(workdir2, "/sen_Skin.pdf"), width = 3, height = 7.5)
# par(mar = c(5, 5, 5, 5))
# image(1:ncol(sen_heat[,c("Vehicle", "ABT-737")]), 1:nrow(sen_heat[,c("Vehicle", "ABT-737")]), t(sen_heat[,c("Vehicle", "ABT-737")]), col =  colorRampPalette(c("royalblue3", "white", "firebrick3"))(1000), zlim = zlim, axes = F, xlab = NA, ylab = NA)
# axis(1, 1:ncol(sen_heat[,c("Vehicle", "ABT-737")]), colnames(sen_heat[,c("Vehicle", "ABT-737")]), par(las = 2), lwd = 0, lwd.ticks = 1)
# # axis(2, 1:nrow(sen_heat[,c("Vehicle", "ABT-737")]), rownames(sen_heat[,c("Vehicle", "ABT-737")]), par(las = 2), lwd = 0, lwd.ticks = 1)
# axis(4, 1:nrow(sen_heat[,c("Vehicle", "ABT-737")]), sen_heat[,"p-value"], lwd = 0, lwd.ticks = 1)
# for (z in 1:ncol(sen_heat[,c("Vehicle", "ABT-737")])){for(y in 1:nrow(sen_heat[,c("Vehicle", "ABT-737")])){text(z, y, signif(sen_heat[, c("Vehicle", "ABT-737")][y, z], 2))}}
# box(which = "plot", lty = "solid")
# dev.off()
# 
# 
# # ####################################
# sen_gg = apply(as.matrix(tissues),1, function(n) data.frame(data.frame(read.csv(paste0("/home/labs/amit/tomerlan/ABT_project/", n, "/sen_FC_", n, ".csv"), sep = ",", header = T))[,-1], tissue = rep(n, nrow(data.frame(read.csv(paste0("/home/labs/amit/tomerlan/ABT_project/", n, "/sen_FC_", n, ".csv"), sep = ",", header = T))))))
# good_sen = as.character(unique(unlist(apply(as.matrix(1:4), 1, function(n) sen_gg[[n]][which(sen_gg[[n]][,"UTT_pvalue"]<alpha),"gene"]))))
# sen_gg = do.call(rbind, sen_gg)
# sen_gg = sen_gg[, -c(3, 5, 6)]

# sasp_gg = apply(as.matrix(tissues),1, function(n) data.frame(data.frame(read.csv(paste0("/home/labs/amit/tomerlan/ABT_project/", n, "/sasp_FC_", n, ".csv"), sep = ",", header = T))[,-1], tissue = rep(n, nrow(data.frame(read.csv(paste0("/home/labs/amit/tomerlan/ABT_project/", n, "/sasp_FC_", n, ".csv"), sep = ",", header = T))))))
# good_sasp = as.character(unique(unlist(apply(as.matrix(1:4), 1, function(n) sasp_gg[[n]][which(sasp_gg[[n]][,"UTT_pvalue"]<alpha),"gene"]))))
# sasp_gg = do.call(rbind, sasp_gg)
# sasp_gg = sasp_gg[, -c(3, 5, 6)]

# sasp_sen = rbind(sasp_gg)
# colnames(sasp_sen) = c("Gene","Vehicle", "ABT-737", "Tissue")

# levels = c("Ccl3", "Ccl5", "Ccl8", "Cxcl9", "Cxcl10", "Cxcl12", "Il1b", "Icam1", "Tnf", "Il6st", "Mmp3", "Mmp14", "Timp1", "Igfbp2", "Igfbp5", "Igfbp7", "Tnfrsf1a", "Tnfrsf1b", "Tnfrsf23", "Ang", "Egfr", "Cdkn1a", "Cdkn1c", "Hmgb1")
# 
# ggplot(melt(sasp_sen, id.vars = c("Tissue", "Gene")), aes(x = factor(Gene), y = value, fill = variable)) +
#   geom_bar(position = "dodge", stat = "identity") +
#   theme_light() +
#   scale_fill_manual(values = c("cornflowerblue", "orange")) +
#   xlab("Gene") +
#   ylab("logFC(EG/Young)") +
#   theme(axis.line.x = element_line(size = 0.1),
#         axis.line.y = element_line(size = 0.1),
#         axis.ticks.x = element_line(size = 0.1),
#         axis.ticks.y = element_line(size = 0.1),
#         text = element_text(size = 4),
#         axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
#         axis.text.y = element_text(size = 4),
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         legend.key.size = unit(1, "mm")) +
#     coord_flip() +
#     # geom_signif(comparisons = list(c("ABT-737", "Vehicle")),manual = T, inherit.aes = F, aes(group = variable), map_signif_level = F, annotation = "1", textsize = 2, margin_top = 0.4) +
#     # geom_signif(annotation = c("1","**"), y_position=c(1.8, 6.8), xmin=c(0.875, 1.875), xmax=c(1.125, 2.125), tip_length = 0.01) +
#   # geom_signif(data=data.frame(x=c(0.875, 1.875), xend=c(1.125, 2.125),y=c(5.8, 8.5), annotation=c("**", "NS")),aes(x=x,xend=xend, y=y, yend=y, annotation=annotation)) +
#   facet_wrap(~Tissue, nrow = 1, ncol = 4) +
#   guides(fill = guide_legend(title = "Experimental group")) +
#   ggsave(file = paste0(workdir2, "/sasp_sen.pdf"))
# 
# y_position = 1.05*max(mean(norm[gene, which(colnames(norm) %in% "Y")]), mean(norm[gene, which(colnames(norm) %in% "DM")]))
# 
# ##############################################
# sasp_gg = apply(as.matrix(tissues),1, function(n) data.frame(data.frame(read.csv(paste0("/home/labs/amit/tomerlan/ABT_project/", n, "/sasp_FC_", n, ".csv"), sep = ",", header = T))[,-1], tissue = rep(n, nrow(data.frame(read.csv(paste0("/home/labs/amit/tomerlan/ABT_project/", n, "/sasp_FC_", n, ".csv"), sep = ",", header = T))))))
# good_sasp = as.character(unique(unlist(apply(as.matrix(1:4), 1, function(n) sasp_gg[[n]][which(sasp_gg[[n]][,"UTT_pvalue"]<alpha),"gene"]))))
# sasp_gg = do.call(rbind, sasp_gg)
# sasp_gg = sasp_gg[, -c(3, 5)]
# sen_gg = apply(as.matrix(tissues),1, function(n) data.frame(data.frame(read.csv(paste0("/home/labs/amit/tomerlan/ABT_project/", n, "/sen_FC_", n, ".csv"), sep = ",", header = T))[,-1], tissue = rep(n, nrow(data.frame(read.csv(paste0("/home/labs/amit/tomerlan/ABT_project/", n, "/sen_FC_", n, ".csv"), sep = ",", header = T))))))
# sen_gg = do.call(rbind, sen_gg)
# sen_gg = sen_gg[, -c(3, 5)]
# 
# sasp_sen = rbind(sasp_gg, sen_gg)
# colnames(sasp_sen) = c("Gene", "Vehicle", "ABT-737", "p-value", "Tissue")
#   
# ggplot(melt(sasp_sen, id.vars = c("Tissue", "Gene")), aes(y = factor(Gene), x = variable, fill = value)) +
#   geom_tile(color = "white") +
#   theme_bw() +
#   theme(axis.line.x = element_line(size = 0.1),
#         axis.line.y = element_line(size = 0.1),
#         axis.ticks.x = element_line(size = 0.1),
#         axis.ticks.y = element_line(size = 0.1),
#         text = element_text(size = 4),
#         axis.text.x = element_text(angle = 90, hjust = 1, size = 4),
#         axis.text.y = element_text(size = 4),
#         strip.background = element_blank(),
#         strip.text.x = element_blank(),
#         legend.key.size = unit(1, "mm")) +
#   scale_fill_gradientn(colours = colorRampPalette(c("royalblue3", "white", "firebrick3"))(1000), values = scales::rescale(c(min(as.matrix(sasp_sen[,2:3]), na.rm = T), min(as.matrix(sasp_sen[,2:3]), na.rm = T)/3 ,  0,  max(as.matrix(sasp_sen[,2:3])/3, na.rm = T), max(as.matrix(sasp_sen[,2:3]), na.rm = T))), limits=c(min(as.matrix(sasp_sen[,2:3])),max(as.matrix(sasp_sen[,2:3])))) +
#   geom_text(aes(label = signif(value, 2)), size=1) +
#   facet_wrap(~Tissue, nrow = 1, ncol = 4) +
#   # scale_y_discrete(breaks = seq(1, 40, by = 2), labels = as.character(1:20))
#   ggsave(file = paste0(workdir2, "/sasp_sen_heatmap.pdf"), width = 100, height = 150, units = "mm")

  
###########Euler
universe = unique(unlist(c(rownames(get(load(paste0(workdir3, "/Kidney/YvsABT.RData")))$genes), 
                           rownames(get(load(paste0(workdir3, "/Liver/YvsABT.RData")))$genes),
                           rownames(get(load(paste0(workdir3, "/Lung/YvsABT.RData")))$genes),
                           rownames(get(load(paste0(workdir3, "/Skin/YvsABT.RData")))$genes))))

YvsABT_sig = unique(c(rownames(topTags(get(load(paste0(workdir3, "/Kidney/YvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)),
                      rownames(topTags(get(load(paste0(workdir3, "/Liver/YvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)),
                      rownames(topTags(get(load(paste0(workdir3, "/Lung/YvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)),
                      rownames(topTags(get(load(paste0(workdir3, "/Skin/YvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000))))

YvsDM_sig = unique(c(rownames(topTags(get(load(paste0(workdir3, "/Kidney/YvsDM.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)),
                     rownames(topTags(get(load(paste0(workdir3, "/Liver/YvsDM.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)),
                     rownames(topTags(get(load(paste0(workdir3, "/Lung/YvsDM.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)),
                     rownames(topTags(get(load(paste0(workdir3, "/Skin/YvsDM.RData"))), sort.by = "logFC", p.value = alpha, n = 3000))))

DMvsABT_sig = unique(c(rownames(topTags(get(load(paste0(workdir3, "/Kidney/DMvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)),
                       rownames(topTags(get(load(paste0(workdir3, "/Liver/DMvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)),
                       rownames(topTags(get(load(paste0(workdir3, "/Lung/DMvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)),
                       rownames(topTags(get(load(paste0(workdir3, "/Skin/DMvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000))))

# euler = euler(list("Y-DM" = YvsDM_sig, "Y-ABT" = YvsABT_sig))
# pdf(file = paste0(workdir2, "/euler.pdf"))
# plot(euler, fill = c("dodgerblue4", "darkgoldenrod1", "cornsilk4", "orange"), border = "transparent", fontsize = 10, counts = list(fontsize = 10), main = "Differential genes")
# dev.off()

euler = euler(list("Y-DM" = YvsDM_sig, "DM-ABT" = DMvsABT_sig, "Detected genes" = universe))
hyp = signif(phyper(length(intersect(DMvsABT_sig, YvsDM_sig)), length(YvsDM_sig), length(universe) - length(YvsDM_sig), length(DMvsABT_sig), lower.tail = FALSE, log.p = FALSE), 2)
pdf(file = paste0(workdir2, "/euler.pdf"))
plot(euler, fill = c("dodgerblue4", "darkgoldenrod1", "cornsilk4", "orange"), border = "transparent", fontsize = 10, counts = list(fontsize = 10), main = paste("Differential genes\nP-value for hypergeometric test =", hyp))
dev.off()

DMvsABT_sig_list = list(rownames(topTags(get(load(paste0(workdir3, "/Kidney/DMvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)),
                        rownames(topTags(get(load(paste0(workdir3, "/Liver/DMvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)),
                        rownames(topTags(get(load(paste0(workdir3, "/Lung/DMvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)),
                        rownames(topTags(get(load(paste0(workdir3, "/Skin/DMvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)))

venn(DMvsABT_sig_list, ilab = TRUE, zcolor = "style", snames = tissues)
dev.copy(pdf, file = paste0(workdir2, "/venn.pdf"))
dev.off()


############FC scatter
FC_scatter_ABT = apply(as.matrix(1:length(tissues)), 1, function(n) cbind(get(load(paste0(workdir3, "/", tissues[n], "/YvsABT.RData")))$table[intersect(rownames(topTags(get(load(paste0(workdir3, "/", tissues[n], "/YvsDM.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)), rownames(topTags(get(load(paste0(workdir3, "/", tissues[n], "/DMvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000))), "logFC"], get(load(paste0(workdir3, "/", tissues[n], "/YvsDM.RData")))$table[intersect(rownames(topTags(get(load(paste0(workdir3, "/", tissues[n], "/YvsDM.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)), rownames(topTags(get(load(paste0(workdir3, "/", tissues[n], "/DMvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000))), "logFC"]))
FC_scatter_ABT_all = data.frame()
for(n in 1:length(tissues)){FC_scatter_ABT_all  = rbind(FC_scatter_ABT_all, FC_scatter_ABT[[n]])}

# FC_scatter_ABT = apply(as.matrix(1:length(tissues)), 1, function(n) cbind(get(load(paste0(workdir3, "/", tissues[n], "/YvsABT.RData")))$table[rownames(topTags(get(load(paste0(workdir3, "/", tissues[n], "/DMvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)), "logFC"], get(load(paste0(workdir3, "/", tissues[n], "/YvsDM.RData")))$table[rownames(topTags(get(load(paste0(workdir3, "/", tissues[n], "/DMvsABT.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)), "logFC"]))
# FC_scatter_ABT_all = data.frame()
# for(n in 1:length(tissues)){FC_scatter_ABT_all  = rbind(FC_scatter_ABT_all, FC_scatter_ABT[[n]])}

FC_scatter = apply(as.matrix(1:length(tissues)), 1, function(n) cbind(get(load(paste0(workdir3, "/", tissues[n], "/YvsABT.RData")))$table[rownames(topTags(get(load(paste0(workdir3, "/", tissues[n], "/YvsDM.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)), "logFC"], get(load(paste0(workdir3, "/", tissues[n], "/YvsDM.RData")))$table[rownames(topTags(get(load(paste0(workdir3, "/", tissues[n], "/YvsDM.RData"))), sort.by = "logFC", p.value = alpha, n = 3000)), "logFC"]))
FC_scatter_all = data.frame()
for(n in 1:length(tissues)){FC_scatter_all  = rbind(FC_scatter_all, FC_scatter[[n]])}

FC_scatter_total = apply(as.matrix(1:length(tissues)), 1, function(n) cbind(get(load(paste0(workdir3, "/", tissues[n], "/YvsABT.RData")))$table[rownames(topTags(get(load(paste0(workdir3, "/", tissues[n], "/YvsDM.RData"))), sort.by = "logFC", p.value = 1, n = 20000)), "logFC"], get(load(paste0(workdir3, "/", tissues[n], "/YvsDM.RData")))$table[rownames(topTags(get(load(paste0(workdir3, "/", tissues[n], "/YvsDM.RData"))), sort.by = "logFC", p.value = 1, n = 20000)), "logFC"]))
FC_scatter_total_all = data.frame()
for(n in 1:length(tissues)){FC_scatter_total_all  = rbind(FC_scatter_total_all, FC_scatter_total[[n]])}

rot = matrix(c(cos(-45*pi/180), sin(-45*pi/180), -sin(-45*pi/180),cos(-45*pi/180)), 2, 2)
ff = as.data.frame(t(apply(as.matrix(FC_scatter_all), 1, function(x) x%*%rot)))

ggplot(data = FC_scatter_all, aes(x = V2, y = V1)) + 
  # geom_point(data = FC_scatter_total_all, aes(x = V2, y = V1), color = "grey", alpha = 0.4, size = 0.0005, shape = 20) +
  geom_point(data = FC_scatter_all, aes(x = V2, y = V1), color = "cornflowerblue", alpha = 0.8, size = 0.2, shape = 20) +
  geom_point(data = FC_scatter_ABT_all, aes(x = V2, y = V1), color = "lightpink", alpha = 0.4, size = 0.2, shape = 20) +
  # stat_smooth(method = "lm", aes(x = V2, y = V1), color = "grey", level = 0.95, size = 0.2) +
  # ggtitle(paste0("R^2 = ", signif(summary(lm(FC_scatter_all[, 1] ~ FC_scatter_all[, 2]))$r.squared, 2), "\nTilt = ", signif(summary(lm(FC_scatter_all[, 1] ~ FC_scatter_all[, 2]))$coefficients[2,1], 2), "\nP-value ~ ", signif(summary(lm(ff[, 1] ~ ff[, 2]))$coefficients[2, 4], 2))) +
  # ylab("Log(ABT-737/Young)") + 
  # xlab("Log(Vehicle/Young)") +
  theme_classic() +
  theme(axis.line.x = element_line(size = 0.1),
        axis.line.y = element_line(size = 0.1),
        axis.ticks.x = element_line(size = 0.1),
        axis.ticks.y = element_line(size = 0.1),
        text = element_text(size = 4),
        axis.text.x = element_text(size = 4),
        axis.text.y = element_text(size = 4),
        plot.title = element_text(size = 4),
        axis.title = element_blank()) +
  geom_abline(intercept = 0, slope = 1, size = 0.1) +
  geom_abline(intercept = 0, slope = 0, size = 0.1) +
  scale_x_continuous(limits = c(-12, 12), breaks = c(-12, -6, 0 , 6, 12)) +
  scale_y_continuous(limits = c(-12, 12), breaks = c(-12, -6, 0, 6, 12)) +
  # geom_abline(intercept = summary(lm(FC_scatter_all[, 1] ~ FC_scatter_all[, 2]))$coefficients[1, 1], slope = summary(lm(FC_scatter_all[, 1] ~ FC_scatter_all[, 2]))$coefficients[2, 1], size = 0.1, color = "black", linetype = 2) +
  # geom_rangeframe(data = data.frame(x = c(-15, 15), y = c(-15, 15)), aes(x, y), size = 0.1)
  ggsave(file = paste0(workdir2, "/FC_scatter.pdf"), width = 41, height = 39, units = "mm")

summary(lm(FC_scatter_all[, 1] ~ FC_scatter_all[, 2]))
summary(lm(ff[, 1] ~ ff[, 2]))
