library(sp)
library(SeuratObject)
library(Seurat)
library(tidyverse)
library(patchwork)
library(dplyr)
library(ggpubr)
library(purrr)
library(broom)
library(ggplot2)
library(clusterProfiler)
library(org.Mm.eg.db)
library(ggrepel)
library(EnhancedVolcano)
library(withr)
library(pheatmap)
library(gtable)
library(tidyr)
library(rstatix)
library(AUCell)
library(readxl)   
library(stringr) 
library(aPEAR)
library(cols4all)
library(DOSE)
library(Rtsne)
library(vcfR)
install.packages("broom")
BiocManager::install("biomaRt")
library(biomaRt)
install.packages("parallelly")
install.packages("BiocManager")
BiocManager::install("GOSemSim")
install.packages("promises")
###################################################################
data_dirs <- c("H:/R/RROW/CD1", "H:/R/RROW/CD2", 
               "H:/R/RROW/HFD2", "H:/R/RROW/HFD2")
group_labels <- c("CD", "CD", "HFD", "HFD")
seurat_list <- list()
for (i in 1:length(data_dirs)){
  data <- Read10X(data.dir = data_dirs[i])
  seurat_object <- CreateSeuratObject(counts = data, project = paste0(i), min.cells = 100)
  seurat_object$group <- group_labels[i] 
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^mt-")
  seurat_object <- subset(seurat_object, subset = nFeature_RNA > 200  & percent.mt < 5 & nCount_RNA >= 500)
  seurat_list[[paste0("Sample_", i)]] <- seurat_object
}
combined.dataset <- merge(seurat_list[[1]], y = seurat_list[2:length(seurat_list)])
DefaultAssay(combined.dataset) = "RNA"
combined.dataset <- NormalizeData(combined.dataset, normalization.method = "LogNormalize", scale.factor = 10000)
# 
head(combined.dataset[[]])
combined.dataset <- FindVariableFeatures(combined.dataset, selection.method = "vst", nfeatures = 2000)
# Identify the 10 most highly variable genes
top10_combined.dataset <- head(VariableFeatures(combined.dataset), 10)
# plot variable features with and without labels
plot1 <- VariableFeaturePlot(combined.dataset)
plot2 <- LabelPoints(plot = plot1, points = top10_combined.dataset, repel = TRUE)
all.genes <- rownames(combined.dataset)
combined.dataset <- ScaleData(combined.dataset, features = all.genes)

combined.dataset <- RunPCA(combined.dataset, features = VariableFeatures(object = combined.dataset))
print(combined.dataset[["pca"]], dims = 1:15, nfeatures = 5)  
ElbowPlot(combined.dataset)
VizDimLoadings(combined.dataset, dims = 1:15, reduction = "pca")
DimHeatmap(combined.dataset, dims = 1:15, cells = 200, balanced = TRUE)

combined.dataset <- FindNeighbors(combined.dataset, dims = 1:15)
combined.dataset <- FindClusters(combined.dataset, resolution = 0.1)
combined.dataset <- RunTSNE(combined.dataset, dims = 1:10 )
DimPlot(combined.dataset, reduction = "tsne", 
        group.by = "seurat_clusters", label = TRUE, pt.size = 0.8)
DimPlot(combined.dataset,label = T,reduction = 'tsne',pt.size =0.2)

combined.dataset <- JoinLayers(combined.dataset)
combined.dataset_3<- subset(combined.dataset, idents = c("0","6","8","1", "2", "3","5","7","9","10"))
current_idents <- Idents(combined.dataset_3)
new_idents <- as.character(current_idents)
new_idents[new_idents %in% c("0", "6")] <- "0_6"  
Idents(combined.dataset_3) <- new_idents
combined.dataset_3.markers <- FindAllMarkers(combined.dataset_3, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
head(combined.dataset_3.markers)
combined.dataset_3.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
top10_combined.dataset_3 <- combined.dataset_3.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(combined.dataset_3, features = top10_combined.dataset_3$gene) + NoLegend()+
  scale_fill_gradient2(low = "lightgray", high = "red1")
DoHeatmap(combined.dataset_3, features = top10_combined.dataset_3$gene) + 
  NoLegend() +
  scale_fill_gradient2(low = "lightgray", mid = "white", high = "red1", midpoint = 0)
top5_combined.dataset_3 <- combined.dataset_3.markers %>% group_by(cluster) %>% top_n(n = 5, wt = avg_log2FC)
DotPlot(
  object = combined.dataset_3,                  
  features = gene2  ,  # deg list
  group.by = "seurat_clusters",                # cluster
  cols = c("lightgray", "red"),                #
  dot.scale = 6,                               # 
  scale = TRUE,                                # 
  cluster.idents = FALSE                       # 
) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1)  # 
  )

gene1 <- c("Gad1", "Gja1", "Mog", "Pdgfra", "Fam216b", "Selplg", "Slc17a6", "Flt1", "Rgs5")  
VlnPlot(object = combined.dataset_3,   
        features = gene1,  # 
        split.plot = TRUE,  # 
        split.by = "group",   
        assay = "RNA",   
        pt.size = 0.5)  
tsne_params <- list(perplexity = 10, learning_rate = 200)
combined.dataset_3 <- RunTSNE(combined.dataset_3, dims = 1:10 )
DimPlot(combined.dataset_3, reduction = "tsne",  
        group.by = "seurat_clusters", label = TRUE, pt.size = 0.8)
DimPlot(combined.dataset_3,label = T,reduction = 'tsne',pt.size =0.2)
gene2 <- c("Gad1", "Gad2", "Slc32a1", 
           "Ppp1r1b","Stard5",
           "Gja1", "Mog", "Pdgfra", "Fam216b", "Selplg", "Slc17a6", "Flt1", "Rgs5") 
DotPlot(
  object = combined.dataset_3,
  features = gene2,
  group.by = "seurat_clusters",
  cols = c("lightgray", "red"),
  col.min = 0,
  col.max = 1.7,
  dot.scale = 6,
  scale = TRUE,
  scale.max = 50,
  cluster.idents = FALSE
) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
  )
tsne_params <- list(perplexity = 10, learning_rate = 200)
combined.dataset_3 <- RunTSNE(combined.dataset_3, dims = 1:10, do.fast = TRUE, tsne.params = tsne_params)
DimPlot(combined.dataset_3, reduction = "tsne",  
        group.by = "seurat_clusters", label = TRUE, pt.size = 1)
DimPlot(combined.dataset_3, reduction = "tsne",
        group.by = "group", label = TRUE, pt.size = 0.8)

colnames(combined.dataset_3@meta.data)
head(combined.dataset_3@meta.data)
DimPlot(combined.dataset_3, reduction = "tsne", group.by = "orig.ident")
colors <- c("1" = "green", "2" = "blue", "3" = "cadetblue4", "4" = "coral")
p11 <- DimPlot(combined.dataset_3, reduction = "tsne", group.by = "orig.ident",pt.size = 0.8)
p11 + scale_color_manual(values = colors)
Idents(object = combined.dataset_3) <- "seurat_clusters"
levels(x = combined.dataset_3)
DoHeatmap(combined.dataset_3, features = top10_combined.dataset_3$gene) + NoLegend()+
  scale_fill_gradient(low = "white", high = "red1")
DoHeatmap(combined.dataset_3, features = top10_combined.dataset_3$gene) + NoLegend()
selected_data <- subset(combined.dataset_3, idents = c("0_6","1", "2", "3","5","7","9","8"))
DoHeatmap(selected_data, features = top10_combined.dataset_3$gene) + NoLegend() + 
  theme(axis.text.y = element_text(size = 2))
DefaultAssay(combined.dataset_3) = "RNA"
selected_genes <- head(top10_combined.dataset_3$gene, 5)
##########################################################################################################    GABAergic_neurons
#GABAergic_neurons
Idents(object = combined.dataset) <- "seurat_clusters"
GABAergic_neurons <- subset(combined.dataset_3, idents = c("0,6"))
Idents(GABAergic_neurons) <- GABAergic_neurons$group
markers.GABAergic_neurons <- FindMarkers(GABAergic_neurons, 
                                  ident.1 = "STIM", ident.2 = "CTRL", 
                                  min.pct = 0.05, logfc.threshold = 0.01,
)
markers.GABAergic_neurons$p_val_adj<- p.adjust(markers.GABAergic_neurons$p_val, method = "BH")
head(markers.GABAergic_neurons )
summary(markers.GABAergic_neurons, n = 10)
total_genes_GABAergic_neurons <- sum(rowSums(GetAssayData(GABAergic_neurons, 
                                                          layer = "counts") > 0) >= 0.1 * ncol(GABAergic_neurons))
table(markers.GABAergic_neurons$p_val_adj < 0.05) 
prop.table(table(markers.GABAergic_neurons$p_val_adj < 0.05)) 
table(sign(markers.GABAergic_neurons$avg_log2FC)) 
prop.table(table(sign(markers.GABAergic_neurons$avg_log2FC))) # 
sig_markers.GABAergic_neurons <- subset(markers.GABAergic_neurons, p_val_adj < 0.05)
write.csv(sig_markers.GABAergic_neurons, 
          file = "E:/sig_markers.GABAergic_neurons.csv")
up_genes_GABAergic_neurons <- sum(sig_markers.GABAergic_neurons$avg_log2FC > 0)
down_genes_GABAergic_neurons <- sum(sig_markers.GABAergic_neurons$avg_log2FC < 0)
print(up_genes_GABAergic_neurons)
print(down_genes_GABAergic_neurons)
up_ratio_GABAergic_neurons <- up_genes_GABAergic_neurons / nrow(sig_markers.GABAergic_neurons)
down_ratio_GABAergic_neurons <- down_genes_GABAergic_neurons / nrow(sig_markers.GABAergic_neurons)
cat("The proportion of significantly up-regulated genes is", round(up_ratio_GABAergic_neurons, 2), "\n")
cat("The proportion of significantly down-regulated genes is", round(down_ratio_GABAergic_neurons, 2), "\n")
total_genes_GABAergic_neurons <- sum(rowSums
                                     (GetAssayData(GABAergic_neurons, layer = "data") > 0) >= 0.1 * ncol(GABAergic_neurons))
ratio_GABAergic_neurons_changed <- nrow(sig_markers.GABAergic_neurons) / total_genes_GABAergic_neurons
print(ratio_GABAergic_neurons_changed)
num_cells_STIM_GABAergic_neurons <- sum(Idents(GABAergic_neurons) == "STIM")
num_cells_CTRL_GABAergic_neurons <- sum(Idents(GABAergic_neurons) == "CTRL")
cat("STIM cell：", num_cells_STIM_GABAergic_neurons, "\n")
cat("CTRL cell：", num_cells_CTRL_GABAergic_neurons, "\n")
cell_counts_sample_GABAergic_neurons <- table(GABAergic_neurons$orig.ident)
print(cell_counts_sample_GABAergic_neurons)
genes_GABAergic_neurons <- rownames(sig_markers.GABAergic_neurons)
geneList_GABAergic_neurons <- bitr(genes_GABAergic_neurons, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# GO
ego_GABAergic_neurons <- enrichGO(gene         = geneList_GABAergic_neurons$ENTREZID,
                           OrgDb        = org.Mm.eg.db,
                           keyType      = "ENTREZID",
                           ont          = "ALL", # Biological Process
                           pAdjustMethod = "BH", # Benjamini & Hochberg (1995) method for p-value adjustment
                           qvalueCutoff = 0.05, # q-value cutoff
                           readable     = TRUE)
head(ego_GABAergic_neurons@result[, c("pvalue")])
print(ego_GABAergic_neurons)
dotplot(ego_GABAergic_neurons, showCategory = 10)
write.csv(ego_GABAergic_neurons@result, file = "E:/enrichment_results-GABAergic_neurons.csv")
###################################################################################################### Glutamatergic_neuron
#Glutamatergic neurons
table(Glutamatergic_neurons$group)
Glutamatergic_neurons <- subset(x = combined.dataset_3, idents = "8")
Idents(Glutamatergic_neurons) <- Glutamatergic_neurons$group
head(Glutamatergic_neurons@meta.data)
cell_counts <- table(Glutamatergic_neurons@meta.data$group)
print(cell_counts)
markers.Glutamatergic_neurons <- FindMarkers(Glutamatergic_neurons, 
                                  ident.1 = "STIM", ident.2 = "CTRL", 
                                  min.pct = 0.2, logfc.threshold = 0.01,
)
write.csv(markers.Glutamatergic_neurons, file = "E:/markers.Glutamatergic_neurons.csv")
markers.Glutamatergic_neurons$p_val_adj<- p.adjust(markers.Glutamatergic_neurons$p_val, method = "BH")
head(markers.Glutamatergic_neurons )
summary(markers.Glutamatergic_neurons, n = 10)
total_genes_Glutamatergic_neurons <- sum(rowSums(GetAssayData(Glutamatergic_neurons, 
                                    layer = "counts") > 0) >= 0.1 * ncol(Glutamatergic_neurons))
table(markers.Glutamatergic_neurons$p_val_adj < 0.05) 
prop.table(table(markers.Glutamatergic_neurons$p_val_adj < 0.05)) 
table(sign(markers.Glutamatergic_neurons$avg_log2FC)) 
prop.table(table(sign(markers.Glutamatergic_neurons$avg_log2FC))) # 
sig_markers.Glutamatergic_neurons <- subset(markers.Glutamatergic_neurons, p_val_adj < 0.05)
write.csv(sig_markers.Glutamatergic_neurons, file = "E:/sig_markers.Glutamatergic_neurons_BH_0.2.csv")
up_genes_Glutamatergic_neurons <- sum(sig_markers.Glutamatergic_neurons$avg_log2FC > 0)
down_genes_Glutamatergic_neurons <- sum(sig_markers.Glutamatergic_neurons$avg_log2FC < 0)
up_ratio_Glutamatergic_neurons <- up_genes_Glutamatergic_neurons / nrow(sig_markers.Glutamatergic_neurons)
down_ratio_Glutamatergic_neurons <- down_genes_Glutamatergic_neurons / nrow(sig_markers.Glutamatergic_neurons)
print(up_genes_Glutamatergic_neurons)
print(down_genes_Glutamatergic_neurons)
cat("The proportion of significantly up-regulated genes is", round(up_ratio_Glutamatergic_neurons, 2), "\n")
cat("The proportion of significantly down-regulated genes is", round(down_ratio_Glutamatergic_neurons, 2), "\n")
total_genes_Glutamatergic_neurons <- 
  sum(rowSums(GetAssayData(Glutamatergic_neurons, layer = "data") > 0) >= 0.1 * ncol(Glutamatergic_neurons))
ratio_Glutamatergic_neurons_changed <- nrow(sig_markers.Glutamatergic_neurons) / total_genes_Glutamatergic_neurons
print(ratio_Glutamatergic_neurons_changed)
num_cells_STIM_Glutamatergic_neurons <- sum(Idents(Glutamatergic_neurons) == "STIM")
num_cells_CTRL_Glutamatergic_neurons <- sum(Idents(Glutamatergic_neurons) == "CTRL")
cat("STIM cell：", num_cells_STIM_Glutamatergic_neurons, "\n")
cat("CTRL cell：", num_cells_CTRL_Glutamatergic_neurons, "\n")
cell_counts_sample_Glutamatergic_neurons <- table(Glutamatergic_neurons$orig.ident)
print(cell_counts_sample_Glutamatergic_neurons)
fea = c("Hcn1")
VlnPlot(Glutamatergic_neurons, features = fea, group.by = "group") + stat_compare_means()
expression_data_Glu <- FetchData(Glutamatergic_neurons, vars = c("Hcn1", "group"))  
head(expression_data_Glu)  
write.csv(expression_data_Glu, file = "E:/Hcn1-expression_Glu.csv")
############################################################################################################Astrocytes
table(Astrocytes$group)
Astrocytes <- subset(x = combined.dataset_3, idents = "1")
Idents(Astrocytes) <- Astrocytes$group
markers.Astrocytes <- FindMarkers(Astrocytes, 
                              ident.1 = "STIM", ident.2 = "CTRL", 
                              min.pct = 0.05, logfc.threshold = 0.01,
)
markers.Astrocytes$p_val_adj<- p.adjust(markers.Astrocytes$p_val, method = "BH")
head(markers.Astrocytes )
summary(markers.Astrocytes, n = 10)
total_genes_Astrocytes <- sum(rowSums(GetAssayData(Astrocytes, layer = "counts") > 0) >= 0.1 * ncol(Astrocytes))
table(markers.Astrocytes$p_val_adj < 0.05)
prop.table(table(markers.Astrocytes$p_val_adj < 0.05)) 
table(sign(markers.Astrocytes$avg_log2FC)) 
prop.table(table(sign(markers.Astrocytes$avg_log2FC))) # 
sig_markers.Astrocytes <- subset(markers.Astrocytes, p_val_adj < 0.05)
write.csv(sig_markers.Astrocytes, file = "数/sig_markers.Astrocytes_BH.csv")
up_Astrocytes <- subset(sig_markers.Astrocytes,p_val_adj < 0.05&avg_log2FC > 0)
down_Astrocytes <- subset(sig_markers.Astrocytes,p_val_adj < 0.05&avg_log2FC < 0)
sig_markers.Astrocytes <- subset(markers.Astrocytes, p_val_adj < 0.05)
write.csv(down_Astrocytes, file = "E:/down_Astrocytes_BH.csv")
print(up_genes_Astrocytes)
print(down_genes_Astrocytes)
up_ratio_Astrocytes <- up_genes_Astrocytes / nrow(sig_markers.Astrocytes)
down_ratio_Astrocytes <- down_genes_Astrocytes / nrow(sig_markers.Astrocytes)
cat("The proportion of significantly up-regulated genes is", round(up_ratio_Astrocytes, 2), "\n")
cat("The proportion of significantly down-regulated genes is", round(down_ratio_Astrocytes, 2), "\n")
total_genes_Astrocytes <- sum(rowSums(GetAssayData(Astrocytes, layer = "data") > 0) >= 0.1 * ncol(Astrocytes))
ratio_Astrocytes_changed <- nrow(sig_markers.Astrocytes) / total_genes_Astrocytes
print(ratio_Astrocytes_changed)
num_cells_STIM_Astrocytes <- sum(Idents(Astrocytes) == "STIM")
num_cells_CTRL_Astrocytes <- sum(Idents(Astrocytes) == "CTRL")
cat("STIM cell：", num_cells_STIM_Astrocytes, "\n")
cat("CTRL cell：", num_cells_CTRL_Astrocytes, "\n")
cell_counts_sample_Astrocytes <- table(Astrocytes$orig.ident)
print(cell_counts_sample_Astrocytes)
genes_Astrocytes <- rownames(sig_markers.Astrocytes)
geneList_Astrocytes <- bitr(genes_Astrocytes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# GO
ego_Astrocytes <- enrichGO(gene         = geneList_Astrocytes$ENTREZID,
                              OrgDb        = org.Mm.eg.db,
                              keyType      = "ENTREZID",
                              ont          = "ALL", # Biological Process
                              pAdjustMethod = "BH", # Benjamini & Hochberg (1995) method for p-value adjustment
                              qvalueCutoff = 0.05, # q-value cutoff
                              readable     = TRUE)
head(ego_Astrocytes@result[, c("pvalue")])
print(ego_Astrocytes)
dotplot(ego_Astrocytes, showCategory = 10)
write.csv(ego_Astrocytes@result, file = "E:/enrichment_results-Astrocytes.csv")
###########################################################################################################oligodendroglioma 
table(oligodendroglioma$group)
oligodendroglioma <- subset(x = combined.dataset_3, idents = "2")
Idents(oligodendroglioma) <- oligodendroglioma$group
markers.oligodendroglioma  <- FindMarkers(oligodendroglioma, 
                                  ident.1 = "STIM", ident.2 = "CTRL", 
                                  min.pct = 0.05, logfc.threshold = 0.01,
)
markers.oligodendroglioma$p_val_adj<- p.adjust(markers.oligodendroglioma$p_val, method = "BH")
head(markers.oligodendroglioma)
summary(markers.oligodendroglioma, n = 10)
total_genes_oligodendroglioma <- 
  sum(rowSums(GetAssayData(oligodendroglioma, layer = "counts") > 0) >= 0.1 * ncol(oligodendroglioma))
table(markers.oligodendroglioma$p_val_adj < 0.05) 
prop.table(table(markers.oligodendroglioma$p_val_adj < 0.05)) 
table(sign(markers.oligodendroglioma$avg_log2FC)) 
prop.table(table(sign(markers.oligodendroglioma$avg_log2FC))) # 
sig_markers.oligodendroglioma <- subset(markers.oligodendroglioma, p_val_adj < 0.05)
write.csv(sig_markers.oligodendroglioma, file = "E:/sig_markers.oligodendroglioma_BH.csv")
up_genes_oligodendroglioma <- sum(sig_markers.oligodendroglioma$avg_log2FC > 0)
down_genes_oligodendroglioma <- sum(sig_markers.oligodendroglioma$avg_log2FC < 0)
print(up_genes_oligodendroglioma)
print(down_genes_oligodendroglioma)
up_ratio_oligodendroglioma <- up_genes_Astrocytes / nrow(sig_markers.oligodendroglioma)
down_ratio_oligodendroglioma <- down_genes_Astrocytes / nrow(sig_markers.oligodendroglioma)
cat("The proportion of significantly up-regulated_oligodendroglioma genes is", round(up_ratio_oligodendroglioma, 2), "\n")
cat("The proportion of significantly down-regulated_oligodendroglioma genes is", round(down_ratio_oligodendroglioma, 2), "\n")
total_genes_oligodendroglioma <- 
  sum(rowSums(GetAssayData(oligodendroglioma, layer = "data") > 0) >= 0.1 * ncol(oligodendroglioma))
ratio_oligodendroglioma_changed <- nrow(sig_markers.oligodendroglioma) / total_genes_oligodendroglioma
print(ratio_oligodendroglioma_changed)
num_cells_STIM_oligodendroglioma <- sum(Idents(oligodendroglioma) == "STIM")
num_cells_CTRL_oligodendroglioma <- sum(Idents(oligodendroglioma) == "CTRL")
cat("STIM组的细胞数目：", num_cells_STIM_oligodendroglioma, "\n")
cat("CTRL组的细胞数目：", num_cells_CTRL_oligodendroglioma, "\n")
cell_counts_sample_oligodendroglioma <- table(oligodendroglioma$orig.ident)
print(cell_counts_sample_oligodendroglioma)
genes_oligodendroglioma <- rownames(sig_markers.oligodendroglioma)
geneList_oligodendroglioma <- bitr(genes_oligodendroglioma, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# GO
ego_oligodendroglioma <- enrichGO(gene         = geneList_oligodendroglioma$ENTREZID,
                           OrgDb        = org.Mm.eg.db,
                           keyType      = "ENTREZID",
                           ont          = "ALL", # Biological Process
                           pAdjustMethod = "BH", # Benjamini & Hochberg (1995) method for p-value adjustment
                           qvalueCutoff = 0.05, # q-value cutoff
                           readable     = TRUE)
head(ego_oligodendroglioma@result[, c("pvalue")])
print(ego_oligodendroglioma)
dotplot(ego_oligodendroglioma, showCategory = 10)
write.csv(ego_oligodendroglioma@result, file = "E:/enrichment_results-oligodendroglioma.csv")
######################################################################################################### Oligodendrocyte progenitor cells
table(OPC$group)
OPC <- subset(x = combined.dataset_3, idents = "3")
Idents(OPC) <- OPC$group
markers.OPC  <- FindMarkers(OPC, 
                                          ident.1 = "STIM", ident.2 = "CTRL", 
                                          min.pct = 0.05, logfc.threshold = 0.01,
)

markers.OPC$p_val_adj<- p.adjust(markers.OPC$p_val, method = "BH")
head(markers.OPC)
summary(markers.OPC, n = 10)
total_genes_OPC <- sum(rowSums(GetAssayData(OPC, layer = "counts") > 0) >= 0.1 * ncol(OPC))
table(markers.OPC$p_val_adj < 0.05)
prop.table(table(markers.OPC$p_val_adj < 0.05))
table(sign(markers.OPC$avg_log2FC)) 
prop.table(table(sign(markers.OPC$avg_log2FC))) # 
sig_markers.OPC <- subset(markers.OPC, p_val_adj < 0.05)
write.csv(sig_markers.OPC, file = "E:/sig_markers.OPC_BH.csv")
up_genes_OPC <- sum(sig_markers.OPC$avg_log2FC > 0)
down_genes_OPC <- sum(sig_markers.OPC$avg_log2FC < 0)
print(up_genes_OPC)
print(down_genes_OPC)
up_ratio_OPC <- up_genes_OPC / nrow(sig_markers.OPC)
down_ratio_OPC <- down_genes_OPC / nrow(sig_markers.OPC)
cat("The proportion of significantly up-regulated_OPC genes is", round(up_ratio_OPC, 2), "\n")
cat("The proportion of significantly down-regulated_OPC genes is", round(down_ratio_OPC, 2), "\n")
total_genes_OPC <- sum(rowSums(GetAssayData(OPC, layer = "data") > 0) >= 0.1 * ncol(OPC))
ratio_OPC_changed <- nrow(sig_markers.OPC) / total_genes_OPC
print(ratio_OPC_changed)
num_cells_STIM_OPC <- sum(Idents(OPC) == "STIM")
num_cells_CTRL_OPC <- sum(Idents(OPC) == "CTRL")
cat("STIM cell：", num_cells_STIM_OPC, "\n")
cat("CTRL cell：", num_cells_CTRL_OPC, "\n")
cell_counts_sample_OPC <- table(OPC$orig.ident)
print(cell_counts_sample_OPC)
genes_OPC <- rownames(sig_markers.OPC)
geneList_OPC <- bitr(genes_OPC, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# GO
ego_OPC <- enrichGO(gene         = geneList_OPC$ENTREZID,
                                  OrgDb        = org.Mm.eg.db,
                                  keyType      = "ENTREZID",
                                  ont          = "ALL", # Biological Process
                                  pAdjustMethod = "BH", # Benjamini & Hochberg (1995) method for p-value adjustment
                                  qvalueCutoff = 0.05, # q-value cutoff
                                  readable     = TRUE)
head(ego_OPC@result[, c("pvalue")])
print(ego_OPC)
dotplot(ego_OPC, showCategory = 10)
write.csv(ego_OPC@result, file = "E:/enrichment_results-OPC.csv")
##########################################################################################################         Epen cell
#Epen cell
table(Epen$group)
Epen <- subset(x = combined.dataset_3, idents = "5")
Idents(Epen) <- Epen$group
markers.Epen  <- FindMarkers(Epen, 
                            ident.1 = "STIM", ident.2 = "CTRL", 
                            min.pct = 0.05, logfc.threshold = 0.01,
)
markers.Epen$p_val_adj<- p.adjust(markers.Epen$p_val, method = "BH")
head(markers.Epen)
summary(markers.Epen, n = 10)
total_genes_Epen <- sum(rowSums(GetAssayData(Epen, layer = "counts") > 0) >= 0.1 * ncol(Epen))
table(markers.Epen$p_val_adj < 0.05) 
prop.table(table(markers.Epen$p_val_adj < 0.05))
table(sign(markers.Epen$avg_log2FC)) 
prop.table(table(sign(markers.Epen$avg_log2FC))) # 
sig_markers.Epen <- subset(markers.Epen, p_val_adj < 0.05)
write.csv(sig_markers.Epen, file = "E:/sig_markers.Epen_BH.csv")
up_genes_Epen <- sum(sig_markers.Epen$avg_log2FC > 0)
down_genes_Epen<- sum(sig_markers.Epen$avg_log2FC < 0)
print(up_genes_Epen)
print(down_genes_Epen)
up_ratio_Epen <- up_genes_Epen / nrow(sig_markers.Epen)
down_ratio_Epen <- down_genes_Epen / nrow(sig_markers.Epen)
cat("The proportion of significantly up-regulated_Epen genes is", round(up_ratio_Epen, 2), "\n")
cat("The proportion of significantly down-regulated_Epen genes is", round(down_ratio_Epen, 2), "\n")
total_genes_Epen <- sum(rowSums(GetAssayData(Epen, layer = "data") > 0) >= 0.1 * ncol(Epen))
ratio_Epen_changed <- nrow(sig_markers.Epen) / total_genes_Epen
print(ratio_Epen_changed)
num_cells_STIM_Epen <- sum(Idents(Epen) == "STIM")
num_cells_CTRL_Epen <- sum(Idents(Epen) == "CTRL")
cat("STIM cell：", num_cells_STIM_Epen, "\n")
cat("CTRL cell：", num_cells_CTRL_Epen, "\n")
cell_counts_sample_Epen <- table(Epen$orig.ident)
print(cell_counts_sample_Epen)
genes_Epen <- rownames(sig_markers.Epen)
geneList_Epen <- bitr(genes_Epen, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# GO
ego_Epen <- enrichGO(gene         = geneList_Epen$ENTREZID,
                    OrgDb        = org.Mm.eg.db,
                    keyType      = "ENTREZID",
                    ont          = "ALL", # Biological Process
                    pAdjustMethod = "BH", # Benjamini & Hochberg (1995) method for p-value adjustment
                    qvalueCutoff = 0.05, # q-value cutoff
                    readable     = TRUE)
head(ego_Epen@result[, c("pvalue")])
print(ego_Epen)
dotplot(ego_Epen, showCategory = 10)
write.csv(ego_Epen@result, file = "E:/enrichment_results-Epen.csv")
###########################################################################################################        microglia
table(microglia$group)
microglia <- subset(x = combined.dataset_3, idents = "7")
Idents(microglia) <- microglia$group
markers.microglia <- FindMarkers(microglia, 
                             ident.1 = "STIM", ident.2 = "CTRL", 
                             min.pct = 0.05, logfc.threshold = 0.01,
)
markers.microglia$p_val_adj<- p.adjust(markers.microglia$p_val, method = "BH")
head(markers.microglia)
summary(markers.microglia, n = 10)
total_genes_microglia <- sum(rowSums(GetAssayData(microglia, layer = "counts") > 0) >= 0.1 * ncol(microglia))
table(markers.microglia$p_val_adj < 0.05) 
prop.table(table(markers.microglia$p_val_adj < 0.05)) 
table(sign(markers.microglia$avg_log2FC)) 
prop.table(table(sign(markers.microglia$avg_log2FC))) # 
sig_markers.microglia <- subset(markers.microglia, p_val_adj < 0.05)
up_microglia <- subset(sig_markers.microglia,p_val_adj < 0.05&avg_log2FC > 0)
down_microglia<- subset(sig_markers.microglia,p_val_adj < 0.05&avg_log2FC < 0)
write.csv(down_microglia, file = "E:/down_microglia.csv")
up_genes_microglia <- sum(sig_markers.microglia$avg_log2FC > 0)
down_genes_microglia<- sum(sig_markers.microglia$avg_log2FC < 0)
print(up_genes_microglia)
print(down_genes_microglia)
up_ratio_microglia <- up_genes_microglia / nrow(sig_markers.microglia)
down_ratio_microglia <- down_genes_microglia / nrow(sig_markers.microglia)
cat("The proportion of significantly up-regulated_microglia genes is", round(up_ratio_microglia, 2), "\n")
cat("The proportion of significantly down-regulated_microglia genes is", round(down_ratio_microglia, 2), "\n")
total_genes_microglia <- sum(rowSums(GetAssayData(microglia, layer = "data") > 0) >= 0.1 * ncol(microglia))
ratio_microglia_changed <- nrow(sig_markers.microglia) / total_genes_microglia
print(ratio_microglia_changed)
num_cells_STIM_microglia <- sum(Idents(microglia) == "STIM")
num_cells_CTRL_microglia <- sum(Idents(microglia) == "CTRL")
cat("STIM cell：", num_cells_STIM_microglia, "\n")
cat("CTRL cell：", num_cells_CTRL_microglia, "\n")
cell_counts_sample_microglia <- table(microglia$orig.ident)
print(cell_counts_sample_microglia)
genes_microglia <- rownames(sig_markers.microglia)
geneList_microglia <- bitr(genes_microglia, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# GO
ego_microglia <- enrichGO(gene         = geneList_microglia$ENTREZID,
                     OrgDb        = org.Mm.eg.db,
                     keyType      = "ENTREZID",
                     ont          = "ALL", # Biological Process
                     pAdjustMethod = "BH", # Benjamini & Hochberg (1995) method for p-value adjustment
                     qvalueCutoff = 0.05, # q-value cutoff
                     readable     = TRUE)
head(ego_microglia@result[, c("pvalue")])
print(ego_microglia)
dotplot(ego_microglia, showCategory = 10)
write.csv(ego_microglia@result, file = "E:/enrichment_results-microglia.csv")
###########################################################################################################     endothelial cells
table(endothelial$group)
endothelial <- subset(x = combined.dataset_3, idents = "9")
table(endothelial$group)
Idents(endothelial) <- endothelial$group
markers.endothelial <- FindMarkers(endothelial, 
                                 ident.1 = "STIM", ident.2 = "CTRL", 
                                 min.pct = 0.05, logfc.threshold = 0.01,
)
markers.endothelial$p_val_adj<- p.adjust(markers.endothelial$p_val, method = "BH")
head(markers.endothelial)
summary(markers.endothelial, n = 10)
total_genes_endothelial <- sum(rowSums(GetAssayData(endothelial, layer = "counts") > 0) >= 0.1 * ncol(endothelial))
table(markers.endothelial$p_val_adj < 0.05) 
prop.table(table(markers.endothelial$p_val_adj < 0.05)) 
table(sign(markers.endothelial$avg_log2FC)) 
prop.table(table(sign(markers.endothelial$avg_log2FC))) # 
sig_markers.endothelial <- subset(markers.endothelial, p_val_adj < 0.05)
write.csv(sig_markers.endothelial, file = "E:/sig_markers.endothelial_BH.csv")
up_genes_endothelial <- sum(sig_markers.endothelial$avg_log2FC > 0)
down_genes_endothelial<- sum(sig_markers.endothelial$avg_log2FC < 0)
print(up_genes_endothelial)
print(down_genes_endothelial)
up_ratio_endothelial <- up_genes_endothelial / nrow(sig_markers.endothelial)
down_ratio_endothelial <- down_genes_endothelial / nrow(sig_markers.endothelial)
cat("The proportion of significantly up-regulated_endothelial genes is", round(up_ratio_endothelial, 2), "\n")
cat("The proportion of significantly down-regulated_endothelial genes is", round(down_ratio_endothelial, 2), "\n")
total_genes_endothelial <- sum(rowSums(GetAssayData(endothelial, layer = "data") > 0) >= 0.1 * ncol(endothelial))
ratio_endothelial_changed <- nrow(sig_markers.endothelial) / total_genes_endothelial
print(ratio_endothelial_changed)
num_cells_STIM_endothelial <- sum(Idents(endothelial) == "STIM")
num_cells_CTRL_endothelial <- sum(Idents(endothelial) == "CTRL")
cat("STIM cell：", num_cells_STIM_endothelial, "\n")
cat("CTRL cell：", num_cells_CTRL_endothelial, "\n")
cell_counts_sample_endothelial <- table(endothelial$orig.ident)
print(cell_counts_sample_endothelial)
genes_endothelial <- rownames(sig_markers.endothelial)
geneList_endothelial <- bitr(genes_endothelial, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# GO
ego_endothelial <- enrichGO(gene         = geneList_endothelial$ENTREZID,
                          OrgDb        = org.Mm.eg.db,
                          keyType      = "ENTREZID",
                          ont          = "ALL", # Biological Process
                          pAdjustMethod = "BH", # Benjamini & Hochberg (1995) method for p-value adjustment
                          qvalueCutoff = 0.05, # q-value cutoff
                          readable     = TRUE)
head(ego_endothelial@result[, c("pvalue")])
print(ego_endothelial)
dotplot(ego_endothelial, showCategory = 10)
write.csv(ego_endothelial@result, file = "E:/enrichment_results-endothelial.csv")
##########################################################################################################	mural
#mural
table(mural$group)
mural <- subset(x = combined.dataset_3, idents = "10")
Idents(mural) <- mural$group
markers.mural <- FindMarkers(mural, 
                                   ident.1 = "STIM", ident.2 = "CTRL", 
                                   min.pct = 0.05, logfc.threshold = 0.01,
)

markers.mural$p_val_adj<- p.adjust(markers.mural$p_val, method = "BH")
head(markers.mural)
summary(markers.mural, n = 10)
total_genes_mural <- sum(rowSums(GetAssayData(mural, layer = "counts") > 0) >= 0.1 * ncol(mural))
table(markers.mural$p_val_adj < 0.05) 
prop.table(table(markers.mural$p_val_adj < 0.05)) 
table(sign(markers.mural$avg_log2FC)) 
prop.table(table(sign(markers.mural$avg_log2FC))) # 
sig_markers.mural <- subset(markers.mural, p_val_adj < 0.05)
write.csv(sig_markers.mural, file = "E:/sig_markers.mural_BH.csv")
up_genes_mural <- sum(sig_markers.mural$avg_log2FC > 0)
down_genes_mural<- sum(sig_markers.mural$avg_log2FC < 0)
print(up_genes_mural)
print(down_genes_mural)
up_ratio_mural <- up_genes_mural / nrow(sig_markers.mural)
down_ratio_mural <- down_genes_mural / nrow(sig_markers.mural)
cat("The proportion of significantly up-regulated_mural genes is", round(up_ratio_mural, 2), "\n")
cat("The proportion of significantly down-regulated_mural genes is", round(down_ratio_mural, 2), "\n")
total_genes_mural <- sum(rowSums(GetAssayData(mural, layer = "data") > 0) >= 0.1 * ncol(mural))
ratio_mural_changed <- nrow(sig_markers.mural) / total_genes_mural
print(ratio_mural_changed)
num_cells_STIM_mural <- sum(Idents(mural) == "STIM")
num_cells_CTRL_mural <- sum(Idents(mural) == "CTRL")
cat("STIM cell：", num_cells_STIM_mural, "\n")
cat("CTRL cell：", num_cells_CTRL_mural, "\n")
cell_counts_sample_mural <- table(mural$orig.ident)
print(cell_counts_sample_mural)
genes_mural <- rownames(sig_markers.mural)
geneList_mural <- bitr(genes_mural, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# GO
ego_mural <- enrichGO(gene         = geneList_mural$ENTREZID,
                            OrgDb        = org.Mm.eg.db,
                            keyType      = "ENTREZID",
                            ont          = "ALL", # Biological Process
                            pAdjustMethod = "BH", # Benjamini & Hochberg (1995) method for p-value adjustment
                            qvalueCutoff = 0.05, # q-value cutoff
                            readable     = TRUE)
head(ego_mural@result[, c("pvalue")])
print(ego_mural)
dotplot(ego_mural, showCategory = 10)
write.csv(ego_mural@result, file = "E:/enrichment_results-mural.csv")
################################################################################################################HCN1
neuron_data <- subset(x = combined.dataset_3, idents = c("0,6,8"))
Idents(object = neuron_data) <- "seurat_clusters"
Gad2.positive.cells <- subset(neuron_data, subset =  Gad2 > 0)
Gad2.negative.cells <- subset(neuron_data, subset = Gad2 == 0)
Gad2.positive.cells <- subset(neuron_data, subset =  Gad2 > 0)
Hcn1.positive.cells <- subset(neuron_data, subset = Hcn1 > 0)
Hcn1.Gad2 <- subset(Hcn1.positive.cells, subset = Hcn1 > 0)
fea = c("Gad2")
fea = c("Hcn1","Hcn2","Hcn3","Hcn4")
DotPlot(Gad2.positive.cells, features = fea) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
Hcn2.negative.cells <- subset(neuron_data, subset = Hcn2 > 0)
Hcn3.positive.cells <- subset(neuron_data, subset = Hcn3 > 0)
Hcn4.negative.cells <- subset(neuron_data, subset = Hcn4 > 0)
Gad2_expression <- FetchData(neuron_data, vars = "Gad2")  
neuron_data$Gad2_group <- ifelse(Gad2_expression > 0, "Expressing", "Not Expressing")  
DotPlot(object = neuron_data, features = features, group.by = "Gad2_group")
DotPlot(object = neuron_data,features = fea,group.by = "Gad2_group", col.min=-1,col.max = 2,scale = FALSE,,
        scale.min = 0,scale.max=80,cols =  c("lightgrey", "red"))+ 
  theme(axis.text.x = element_text(angle=0, hjust=.5, vjust=.5))+coord_flip()
dot_plot <- DotPlot(object = neuron_data, features = fea, group.by = "Gad2_group",
                    col.min = -1, col.max = 2, scale = FALSE,
                    scale.min = 0, scale.max = 0, cols = c("lightgrey", "red")) +  # 设置标尺最大值为5
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  coord_flip()
dot_plot
dot_plot_data <- dot_plot$data
percent_expressed <- dot_plot_data$pct.exp
average_expression <- dot_plot_data$avg.exp
#######################################################################################################################     AUCcell
inflammatory_gene <- readxl::read_xlsx("E:/C0028754_disease_gda_summary_ID.xlsx")
gene_human_obsety <- as.list(inflammatory_gene)
combined.dataset_3 <- AddModuleScore(object = combined.dataset_3, features = gene_human_obsety, ctrl = 100, name = 'CD_Features')
colnames(combined.dataset_3@meta.data)
colnames(combined.dataset_3@meta.data)[10] <- 'inflammatory_score'
VlnPlot(combined.dataset_3,features= 'inflammatory_score')
inflammatory_scores <- combined.dataset_3@meta.data$inflammatory_score
head(inflammatory_scores)
scores_df <- data.frame(Cell = rownames(combined.dataset_3@meta.data), Inflammatory_Score = inflammatory_scores)
head(scores_df)
write.csv(scores_df, "E:/inflammatory_scores_2.csv", row.names = FALSE)
cluster_info <- combined.dataset_3@meta.data$seurat_clusters
inflammatory_scores <- combined.dataset_3@meta.data$inflammatory_score
scores_df <- data.frame(Cluster = cluster_info, Inflammatory_Score = inflammatory_scores)
cluster_summary <- scores_df %>%
  group_by(Cluster) %>%
  summarise(
    Mean_Score = mean(Inflammatory_Score, na.rm = TRUE),
    Median_Score = median(Inflammatory_Score, na.rm = TRUE),
    SD_Score = sd(Inflammatory_Score, na.rm = TRUE)
  )
print(cluster_summary)
cluster_scores <- scores_df %>%
  group_by(Cluster) %>%
  summarise(Inflammatory_Scores = list(Inflammatory_Score))
print(cluster_scores)
write.csv(cluster_summary, "E:/cluster_summary_2.csv", row.names = FALSE)
expanded_cluster_scores <- cluster_scores %>%
  unnest(cols = Inflammatory_Scores)
head(expanded_cluster_scores)
write.csv(expanded_cluster_scores, "E:/expanded_cluster_scores_2.csv", row.names = FALSE)
####################################################################################################################
install.packages("BiocManager")
BiocManager::install("Seurat")
inflammatory_gene <- readxl::read_xlsx("E:/Oebsity_gene_Summary_list.xlsx")  
colnames(inflammatory_gene) <- trimws(colnames(inflammatory_gene)) # 去除可能的空格  
gene_column <- colnames(inflammatory_gene)[1]  
capitalize_first_letter <- function(x) {  
  return(sapply(x, function(y) {  
    paste0(toupper(substring(y, 1, 1)), tolower(substring(y, 2)))  
  }))  
}  

inflammatory_gene[[gene_column]] <- capitalize_first_letter(inflammatory_gene[[gene_column]])  
write.csv(inflammatory_gene, "E:/Oebsity_gene_Summary_list_2.csv", row.names = FALSE)
head(inflammatory_gene$gene_symbol)
germ <- JoinLayers(germ)
exprMatrix <- as.matrix(GetAssayData(GABAergic_neurons, assay = "RNA", layer = "data"))
dim(exprMatrix)
cells_rankings <- AUCell_buildRankings(exprMatrix, plotStats = TRUE)
geneSets <- list(inflammatory_gene = inflammatory_gene$gene_symbol)
cells_AUC <- AUCell_calcAUC(geneSets, cells_rankings, aucMaxRank=3851)
cells_assignment <- AUCell_exploreThresholds(cells_AUC, thrP = 0.01, plotHist = TRUE, assign=TRUE)
print(cells_assignment)
auc_matrix <- getAUC(cells_AUC)
write.csv(auc_matrix, "E:/auc_GABAergic_neurons_50_2.csv", row.names = FALSE)
################################################################################################################
numeric_cells_AUC <- as.numeric(auc_matrix)
global_k1_threshold <- cells_assignment$inflammatory_gene$aucThr$thresholds["Global_k1", "threshold"]
df <- data.frame(AUC = numeric_cells_AUC)
ggplot(df, aes(x = AUC)) +
  geom_histogram(aes(y = ..density..), bins = 100, fill = "lightblue", color = "black") +
  geom_density(color = "red", size = 1) +
  geom_vline(xintercept = global_k1_threshold, linetype = "dashed", color = "blue", size = 1) +
  labs(title = "Fibroblasts",
       x = "AUC",
       y = "Density") +
  theme_minimal() +
  annotate("text", x = global_k1_threshold, y = 1, label = paste("Global_k1 Threshold =", round(global_k1_threshold, 4)), color = "blue", angle = 90, vjust = -0.5)
################################################################################################################
cluster_scores <- sapply(levels(Idents(Astrocytes)), function(cluster) {
  cells_in_cluster <- WhichCells(Astrocytes, idents = cluster)
  mean(getAUC(cells_AUC)[, cells_in_cluster])
})
print(cluster_scores)
library(readxl)
genes1 <- read_excel("E:/sig_markers.oligodendroglioma-45.xlsx", col_names = FALSE)
genes2 <- read_excel("E:/Oebsity_gene_Summary_list_2.xlsx", col_names = FALSE)
common_genes_sig_markers.OPC <- intersect(genes1[[1]], genes2[[1]])
print(common_genes_sig_markers.OPC)
write.csv(common_genes_sig_markers.OPC, file = "E:/sig_markers.OPC.csv")
intersection_size <- length(intersect(genes1[[1]], genes2[[1]]))
union_size <- length(union(genes1[[1]], genes2[[1]]))
print(union_size)
jaccard_similarity <- intersection_size / union_size
print(intersection_size)
print(paste("Jaccard:", jaccard_similarity))















