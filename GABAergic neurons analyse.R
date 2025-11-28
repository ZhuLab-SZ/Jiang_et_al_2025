###################################################################################################################GABA neuron
GABA_neurons <- GABAergic_neurons
table(GABA_neurons$group)
GABA_neurons <- NormalizeData(GABA_neurons, normalization.method = "LogNormalize", scale.factor = 10000)
GABA_neurons <- FindVariableFeatures(GABA_neurons, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GABA_neurons)
GABA_neurons <- ScaleData(GABA_neurons, features = all.genes)
GABA_neurons <- RunPCA(GABA_neurons, features = VariableFeatures(object = GABA_neurons))
print(GABA_neurons[["pca"]], dims = 1:5, nfeatures = 5)  #查看前五个高变基因
ElbowPlot(GABA_neurons)
VizDimLoadings(GABA_neurons, dims = 1:15, reduction = "pca")
DimHeatmap(GABA_neurons, dims = 1:10, cells = 200, balanced = TRUE)
GABA_neurons <- FindNeighbors(GABA_neurons, dims = 1:15)
GABA_neurons <- FindClusters(GABA_neurons, resolution = 0.07)
GABA_neurons <- RunUMAP(GABA_neurons, dims = 1:10)
DimPlot(GABA_neurons, reduction = "umap", label = TRUE, pt.size = 0.2)
VlnPlot(GABA_neurons, features = fea, group.by = "group") + stat_compare_means()
Idents(GABA_neurons) <- GABA_neurons$group
Idents(object = GABA_neurons) <- "seurat_clusters"

GABA_neurons <- RunTSNE(GABA_neurons, dims = 1:10 )
DimPlot(GABA_neurons,label = T,reduction = 'tsne',pt.size =0.2)

VlnPlot(GABA_neurons, features = "Hcn1", group.by = "group") + stat_compare_means()
expression_data <- FetchData(
  GABA_neurons,
  vars = c("Hcn1", "group"),
  slot = "data"  # 
)
head(expression_data)
fea1= c("Gad2")
FeaturePlot(GABA_neurons,reduction = "tsne", label = TRUE, ncol = 1,features =fea1,min.cutoff = 0.2,max.cutoff = 0.5)
############################################################################################################################
Idents(object = GABA_neurons) <- "seurat_clusters"
GABA_neurons_0 <- subset(x = GABA_neurons, idents = "0")
Idents(GABA_neurons_0) <- GABA_neurons_0$group
markers.0 <- FindMarkers(GABA_neurons_0, 
                         ident.1 = "STIM", ident.2 = "CTRL", 
                         min.pct = 0.05, logfc.threshold = 0.01,
)
markers.0$p_val_adj<- p.adjust(markers.0$p_val, method = "BH")
head(markers.0)
summary(markers.0, n = 10)
total_genes_0 <- sum(rowSums(GetAssayData(GABA_neurons_0, layer = "counts") > 0) >= 0.1 * ncol(GABA_neurons_0))
table(markers.0$p_val_adj < 0.05) 
prop.table(table(markers.0$p_val_adj < 0.05)) 
table(sign(markers.0$avg_log2FC)) # 
prop.table(table(sign(markers.0$avg_log2FC))) # 
sig_markers.0 <- subset(markers.0, p_val_adj < 0.05)
write.csv(sig_markers.0, file = "E:/2025-05/sig_markers.0_BH.csv")
up_genes_0 <- sum(sig_markers.0$avg_log2FC > 0)
down_genes_0<- sum(sig_markers.0$avg_log2FC < 0)
print(up_genes_0)
print(down_genes_0)
up_ratio_0 <- up_genes_0 / nrow(sig_markers.0)
down_ratio_0 <- down_genes_0 / nrow(sig_markers.0)
cat("The proportion of significantly up-regulated_0 genes is", round(up_ratio_0, 2), "\n")
cat("The proportion of significantly down-regulated_0 genes is", round(down_ratio_0, 2), "\n")
###########################################################################################################################
GABA_neurons_1 <- subset(x = GABA_neurons, idents = "1")
Idents(GABA_neurons_1) <- GABA_neurons_1$group
markers.1 <- FindMarkers(GABA_neurons_1, 
                         ident.1 = "STIM", ident.2 = "CTRL", 
                         min.pct = 0.05, logfc.threshold = 0.01,
)
markers.1$p_val_adj<- p.adjust(markers.1$p_val, method = "BH")
head(markers.1)
summary(markers.1, n = 10)
total_genes_1 <- sum(rowSums(GetAssayData(GABA_neurons_1, layer = "counts") > 0) >= 0.1 * ncol(GABA_neurons_1))
table(markers.1$p_val_adj < 0.05) # 
prop.table(table(markers.1$p_val_adj < 0.05)) #
table(sign(markers.1$avg_log2FC)) # 
prop.table(table(sign(markers.1$avg_log2FC))) # 
sig_markers.1 <- subset(markers.1, p_val_adj < 0.05)
write.csv(sig_markers.1, file = "E:/2025-05/sig_markers.1_BH.csv")
up_genes_1 <- sum(sig_markers.1$avg_log2FC > 0)
down_genes_1 <- sum(sig_markers.1$avg_log2FC < 0)
print(up_genes_1)
print(down_genes_1)
up_ratio_1 <- up_genes_1 / nrow(sig_markers.1)
down_ratio_1 <- down_genes_1 / nrow(sig_markers.1)
cat("The proportion of significantly up-regulated_1 genes is", round(up_ratio_1, 2), "\n")
cat("The proportion of significantly down-regulated_1 genes is", round(down_ratio_1, 2), "\n")
#########################################################################################################################
GABA_neurons_2 <- subset(x = GABA_neurons, idents = "2")
Idents(GABA_neurons_2) <- GABA_neurons_2$group
markers.2 <- FindMarkers(GABA_neurons_2, 
                         ident.1 = "STIM", ident.2 = "CTRL", 
                         min.pct = 0.05, logfc.threshold = 0.01,
)
markers.2$p_val_adj<- p.adjust(markers.2$p_val, method = "BH")
head(markers.2)
summary(markers.2, n = 10)
total_genes_2 <- sum(rowSums(GetAssayData(GABA_neurons_2, layer = "counts") > 0) >= 0.1 * ncol(GABA_neurons_2))
table(markers.2$p_val_adj < 0.05) # 
prop.table(table(markers.2$p_val_adj < 0.05)) # 
table(sign(markers.2$avg_log2FC)) # 
prop.table(table(sign(markers.2$avg_log2FC))) # 
sig_markers.2 <- subset(markers.2, p_val_adj < 0.05)
up_genes_2 <- sum(sig_markers.2$avg_log2FC > 0)
down_genes_2 <- sum(sig_markers.2$avg_log2FC < 0)
print(up_genes_2)
print(down_genes_2)
up_ratio_2 <- up_genes_2 / nrow(sig_markers.2)
down_ratio_2 <- down_genes_2 / nrow(sig_markers.2)
cat("The proportion of significantly up-regulated_2 genes is", round(up_ratio_2, 2), "\n")
cat("The proportion of significantly down-regulated_2 genes is", round(down_ratio_2, 2), "\n")
#######################################################################################################################
GABA_neurons_3 <- subset(x = GABA_neurons, idents = "3")
Idents(GABA_neurons_3) <- GABA_neurons_3$group
markers.3 <- FindMarkers(GABA_neurons_3, 
                         ident.1 = "STIM", ident.2 = "CTRL", 
                         min.pct = 0.05, logfc.threshold = 0.01,
)
markers.3$p_val_adj<- p.adjust(markers.3$p_val, method = "BH")
head(markers.3)
summary(markers.3, n = 10)
total_genes_3 <- sum(rowSums(GetAssayData(GABA_neurons_3, layer = "counts") > 0) >= 0.1 * ncol(GABA_neurons_3))
table(markers.3$p_val_adj < 0.05) # 
prop.table(table(markers.3$p_val_adj < 0.05)) # 
table(sign(markers.3$avg_log2FC)) # 
prop.table(table(sign(markers.3$avg_log2FC))) # 
sig_markers.3 <- subset(markers.3, p_val_adj < 0.05)
write.csv(sig_markers.3, file = "E:/V3/2025-05/sig_markers.3_BH.csv")
up_genes_3 <- sum(sig_markers.3$avg_log2FC > 0)
down_genes_3<- sum(sig_markers.3$avg_log2FC < 0)
print(up_genes_3)
print(down_genes_3)
up_ratio_3 <- up_genes_3 / nrow(sig_markers.3)
down_ratio_3 <- down_genes_3 / nrow(sig_markers.3)
cat("The proportion of significantly up-regulated_3 genes is", round(up_ratio_3, 2), "\n")
cat("The proportion of significantly down-regulated_3 genes is", round(down_ratio_3, 2), "\n")
######################################################################################################################
GABA_neurons_4 <- subset(x = GABA_neurons, idents = "4")
Idents(GABA_neurons_4) <- GABA_neurons_4$group
markers.4 <- FindMarkers(GABA_neurons_4, 
                         ident.1 = "STIM", ident.2 = "CTRL", 
                         min.pct = 0.05, logfc.threshold = 0.01,
)
markers.4$p_val_adj<- p.adjust(markers.4$p_val, method = "BH")
head(markers.4)
summary(markers.4, n = 10)
total_genes_4 <- sum(rowSums(GetAssayData(GABA_neurons_4, layer = "counts") > 0) >= 0.1 * ncol(GABA_neurons_4))
table(markers.4$p_val_adj < 0.05) 
prop.table(table(markers.4$p_val_adj < 0.05)) 
table(sign(markers.4$avg_log2FC)) 
prop.table(table(sign(markers.4$avg_log2FC))) # 
sig_markers.4 <- subset(markers.4, p_val_adj < 0.05)
write.csv(sig_markers.4, file = "E:V3/2025-05/sig_markers.4_BH.csv")
up_genes_4 <- sum(sig_markers.4$avg_log2FC > 0)
down_genes_4<- sum(sig_markers.4$avg_log2FC < 0)
print(up_genes_4)
print(down_genes_4)
up_ratio_4 <- up_genes_4 / nrow(sig_markers.4)
down_ratio_4 <- down_genes_4 / nrow(sig_markers.4)
cat("The proportion of significantly up-regulated_4 genes is", round(up_ratio_4, 2), "\n")
cat("The proportion of significantly down-regulated_4 genes is", round(down_ratio_4, 2), "\n")
########################################################################################################################
GABA_neurons_5 <- subset(x = GABA_neurons, idents = "5")
Idents(GABA_neurons_5) <- GABA_neurons_5$group
markers.5 <- FindMarkers(GABA_neurons_5, 
                         ident.1 = "STIM", ident.2 = "CTRL", 
                         min.pct = 0.05, logfc.threshold = 0.01,
)
markers.5$p_val_adj<- p.adjust(markers.5$p_val, method = "BH")
head(markers.5)
summary(markers.5, n = 10)
total_genes_5 <- sum(rowSums(GetAssayData(GABA_neurons_5, layer = "counts") > 0) >= 0.1 * ncol(GABA_neurons_5))
table(markers.5$p_val_adj < 0.05) # 
prop.table(table(markers.5$p_val_adj < 0.05)) # 
table(sign(markers.5$avg_log2FC)) # 
prop.table(table(sign(markers.5$avg_log2FC))) # 
sig_markers.5 <- subset(markers.5, p_val_adj < 0.05)
write.csv(sig_markers.5, file = "E:/2025-05/sig_markers.5_BH.csv")
up_genes_5 <- sum(sig_markers.5$avg_log2FC > 0)
down_genes_5<- sum(sig_markers.5$avg_log2FC < 0)
print(up_genes_5)
print(down_genes_5)
up_ratio_5 <- up_genes_5 / nrow(sig_markers.5)
down_ratio_5 <- down_genes_5 / nrow(sig_markers.5)
cat("The proportion of significantly up-regulated_5 genes is", round(up_ratio_5, 2), "\n")
cat("The proportion of significantly down-regulated_5 genes is", round(down_ratio_5, 2), "\n")
########################################################################################################################
GABA_neurons_6 <- subset(x = GABA_neurons, idents = "6")
Idents(GABA_neurons_6) <- GABA_neurons_6$group
markers.6 <- FindMarkers(GABA_neurons_6, 
                         ident.1 = "STIM", ident.2 = "CTRL", 
                         min.pct = 0.2, logfc.threshold = 0.01,
                         test.use = "negbinom"
)
markers.6$p_val_adj<- p.adjust(markers.6$p_val, method = "BH")
head(markers.6)
summary(markers.6, n = 10)
total_genes_6 <- sum(rowSums(GetAssayData(GABA_neurons_6, layer = "counts") > 0) >= 0.1 * ncol(GABA_neurons_6))
table(markers.6$p_val_adj < 0.05) # 
prop.table(table(markers.6$p_val_adj < 0.05)) # 
table(sign(markers.6$avg_log2FC)) # 
prop.table(table(sign(markers.6$avg_log2FC))) # 
sig_markers.6 <- subset(markers.6, p_val_adj < 0.05)
#write.csv(sig_markers.6, file = "E:/V3/2025-05/sig_markers.6_BH.csv")
up_genes_6 <- sum(sig_markers.6$avg_log2FC > 0)
down_genes_6<- sum(sig_markers.6$avg_log2FC < 0)
print(up_genes_6)
print(down_genes_6)
up_ratio_6 <- up_genes_6 / nrow(sig_markers.6)
down_ratio_6 <- down_genes_6 / nrow(sig_markers.6)
cat("The proportion of significantly up-regulated_6 genes is", round(up_ratio_6, 2), "\n")
cat("The proportion of significantly down-regulated_6 genes is", round(down_ratio_6, 2), "\n")
##########################################################################################################################
GABA_neurons_7 <- subset(x = GABA_neurons, idents = "7")
Idents(GABA_neurons_7) <- GABA_neurons_7$group
markers.7 <- FindMarkers(GABA_neurons_7, 
                         ident.1 = "STIM", ident.2 = "CTRL", 
                         min.pct = 0.05, logfc.threshold = 0.01,
                         test.use = "negbinom"
)
markers.7$p_val_adj<- p.adjust(markers.7$p_val, method = "BH")
head(markers.7)
summary(markers.7, n = 10)
total_genes_7 <- sum(rowSums(GetAssayData(GABA_neurons_7, layer = "counts") > 0) >= 0.1 * ncol(GABA_neurons_7))
table(markers.7$p_val_adj < 0.05) 
prop.table(table(markers.7$p_val_adj < 0.05)) # 
table(sign(markers.7$avg_log2FC)) # 
prop.table(table(sign(markers.7$avg_log2FC))) # 
sig_markers.7 <- subset(markers.7, p_val_adj < 0.05)
write.csv(sig_markers.7, file = "E:/2025-05/sig_markers.7_BH.csv")
up_genes_7 <- sum(sig_markers.7$avg_log2FC > 0)
down_genes_7<- sum(sig_markers.7$avg_log2FC < 0)
print(up_genes_7)
print(down_genes_7)
up_ratio_7 <- up_genes_7 / nrow(sig_markers.7)
down_ratio_7 <- down_genes_7 / nrow(sig_markers.7)
cat("The proportion of significantly up-regulated_7 genes is", round(up_ratio_7, 2), "\n")
cat("The proportion of significantly down-regulated_7 genes is", round(down_ratio_7, 2), "\n")