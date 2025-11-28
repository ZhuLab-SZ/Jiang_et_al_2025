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
remove.packages("SeuratObject")
install.packages("SeuratObject")

GABAergic_neurons <- NormalizeData(GABAergic_neurons, normalization.method = "LogNormalize", scale.factor = 10000)
GABAergic_neurons <- FindVariableFeatures(GABAergic_neurons, selection.method = "vst", nfeatures = 2000)
all.genes <- rownames(GABAergic_neurons)
GABAergic_neurons <- ScaleData(GABAergic_neurons, features = all.genes)
GABAergic_neurons <- RunPCA(GABAergic_neurons, features = VariableFeatures(object = GABAergic_neurons))
print(GABAergic_neurons[["pca"]], dims = 1:5, nfeatures = 5)  
ElbowPlot(GABAergic_neurons)
VizDimLoadings(GABAergic_neurons, dims = 1:15, reduction = "pca")
DimHeatmap(GABAergic_neurons, dims = 1:10, cells = 200, balanced = TRUE)
GABAergic_neurons <- FindNeighbors(GABAergic_neurons, dims = 1:15)
GABAergic_neurons <- FindClusters(GABAergic_neurons, resolution = 0.07)
GABAergic_neurons <- RunTSNE(GABAergic_neurons, dims = 1:10)
DimPlot(GABAergic_neurons, reduction = "tsne", label = TRUE, pt.size = 0.2)

fea = c("Gad2")
VlnPlot(GABAergic_neurons, features = "Gad2", group.by = "group") + stat_compare_means()
GABAergic_neurons.markers <- FindAllMarkers(GABAergic_neurons, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.1)
head(GABAergic_neurons.markers)
GABAergic_neurons.markers %>% group_by(cluster) %>% top_n(n = 10, wt = avg_log2FC)
top10_GABAergic_neurons <- GABAergic_neurons.markers %>% group_by(cluster) %>% top_n(n = 20, wt = avg_log2FC)
DoHeatmap(GABAergic_neurons, features = top10_GABAergic_neurons$gene) + NoLegend()+
  scale_fill_gradient(low = "white", high = "red1")

total_genes_GABAergic_neurons <- sum(rowSums(GetAssayData(GABAergic_neurons, 
                                                          layer = "counts") > 0) >= 0.1 * ncol(GABAergic_neurons))
Idents(object = GABAergic_neurons) <- "seurat_clusters"
Idents(object = GABAergic_neurons) <- "group"
markers.GABAergic_neurons <- FindMarkers(GABAergic_neurons, 
                                         ident.1 = "STIM", ident.2 = "CTRL", 
                                         min.pct = 0.05, logfc.threshold = 0
)

markers.GABAergic_neurons$p_val_adj<- p.adjust(markers.GABAergic_neurons$p_val, method = "BH")
table(markers.GABAergic_neurons$p_val_adj < 0.05) 
prop.table(table(markers.GABAergic_neurons$p_val_adj < 0.05)) 
table(sign(markers.GABAergic_neurons$avg_log2FC)) 
prop.table(table(sign(markers.GABAergic_neurons$avg_log2FC))) # 
sig_markers.GABAergic_neurons <- subset(markers.GABAergic_neurons, p_val_adj < 0.05)
write.csv(sig_markers.GABAergic_neurons, file = "E:/姜少磊/R/V3/2025-04/sig_markers.GABAergic_neurons_BH校正.csv")

fea2 = c("Ina","Matn2","Syndig1","Rnf220","Sel1l3","Tshz1","Drd3","Unc5cl","Calb2", "Eomes","Slc17a6",
         "Hcn1","Gfra1","Gad1","Gad2","Slc32a1"
         ,"Drd1","Drd2","Ppp1r1b","Stard5")

plot_data <- GABAergic_neurons
plot_data$seurat_clusters <- factor(
  as.character(plot_data$seurat_clusters),
  levels = c("0", "1", "2", "3", "4", "5", "6", "7"))
Idents(plot_data) <- plot_data$seurat_clusters
DotPlot(plot_data, features = fea2,cluster.idents = FALSE, 
        col.min = -1, col.max = 3, scale = TRUE,scale.max = 75, cols = c("white", "red1")
) + 
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    axis.text.y = element_text(face = "italic")
  )

interested_genes <- c("Slc1a2","Slc1a3","Apoe","Sox6","Gpc5",
                      "Grm3","Mertk","Rora","Usp50","Has1","Prkd1",
                      "Syt4","Atp6v0c","Gphn","Cbarp","Selenok","Mdga2",
                      "Slitrk1","Nr4a1","Bex1","Arf5","Lpcat4",
                      "Tiparp","Prlr","Lrrtm2","Bax","Rab11a","Leprotl1",
                      "Bcl2l2","Tmem158","Scg2","Gpc3","Arrdc3",
                      "Gdpd2","Bag6","Scg2",
                      "Grid2","Gabrb1","Snap25",
                      "Grm5","Gria4","Egr1","Kcnn2","Gabrg2","Grik2",
                      "Atad1","Nalcn","Gabra5",
                      "Hcn1","Gad2"
                      ,"Gabra2","Gria2","Rims4","Gabrb3","Scn1a"
                      ,"Gabrg3","Gabbr1")
interested_genes <- c("Hcn1")
markers.GABA_neuron_2 <- subset(markers.GABAergic_neurons, p_val_adj < 1)
p <- EnhancedVolcano(
  markers.GABA_neuron_2,
  x = "avg_log2FC",
  y = "p_val_adj",
  lab = rownames(markers.GABA_neuron_2),
  selectLab = interested_genes,  
  pCutoff = 0.05,     
  FCcutoff = 0.01,         
  cutoffLineWidth = 0.7,      
  cutoffLineType = "twodash", 
  xlim = c(-0.5, 0.5), 
  ylim = c(0, -log10(10e-20)), 
  pointSize = 1,        
  labSize = 6,      
  xlab = bquote(~Log[2] ~ "fold change"),    
  ylab = bquote(~-Log[10] ~ italic(p_value_adj)), 
  axisLabSize = 12,     
  title = "Gene differential expression",    
  titleLabSize = 16,    
  subtitle = bquote(italic("Volcano plot")),  
  subtitleLabSize = 14, 
  legendLabSize = 11,   # legend大小
  col = c( "darkgreen","grey60","orange", "blue" ), 
  colAlpha = 1,       
  gridlines.major = FALSE,     
  gridlines.minor = FALSE
)
p

markers.GABA_neuron$names <- rownames(markers.GABA_neuron)
sig_dge.all <- subset(markers.GABA_neuron,p_val_adj < 0.05 & abs(avg_log2FC) > 0.15)
markers.GABA_neuron <- markers.GABA_neuron %>% mutate(Difference = pct.1-pct.2)
markers.GABA_neuron$group = ""
for (i in 1:nrow(markers.GABA_neuron)) {
  if (markers.GABA_neuron$avg_log2FC[i] >= 1 & markers.GABA_neuron$Difference[i] >= 0.2 & markers.GABA_neuron$pct.2[i] <= 0.05){
    markers.GABA_neuron$group[i] = "up"
  }
  else if(markers.GABA_neuron$avg_log2FC[i] <= (-1) & markers.GABA_neuron$Difference[i] <= -0.2 & markers.GABA_neuron$pct.1[i] <= 0.05){
    markers.GABA_neuron$group[i] = "down"
  }
  else{
    markers.GABA_neuron$group[i] = "no"
  }
}

markers.GABA_neuron$label <- ifelse(markers.GABA_neuron$group %in% c("up","down"), markers.GABA_neuron$names,"")
label <- c("Slc1a2","Slc1a3","Apoe","Sox6","Gpc5",
           "Grm3","Mertk","Rora","Usp50","Has1","Prkd1",
           "Syt4","Atp6v0c","Gphn","Cbarp","Selenok","Mdga2",
           "Slitrk1","Nr4a1","Bex1","Arf5","Lpcat4",
           "Tiparp","Prlr","Lrrtm2","Bax","Rab11a","Leprotl1",
           "Bcl2l2","Tmem158","Scg2","Gpc3","Arrdc3",
           "Gdpd2","Bag6","Scg2",
           "Grid2","Gabrb1","Snap25",
           "Grm5","Gria4","Egr1","Kcnn2","Gabrg2","Grik2"
           ,"Atad1","Nalcn","Gabra5","Hcn1","Gad2"
           ,"Gabra2","Gria2","Rims4","Gabrb3","Scn1a"
           ,"Gabrg3","Gabbr1")
label <- c("Hcn1","Gad2"
)
markers.GABA_neuron$label <- ifelse(markers.GABA_neuron$names %in% label, markers.GABA_neuron$names,"")
ggplot(markers.GABA_neuron,aes(x=Difference,y=avg_log2FC))+
  geom_point(size = 0.9,aes(color = group))+
  scale_color_manual(values = c("#cb5a54","#1176bd","#7b7f80"))+
  geom_label_repel(aes(label = label),fill = "#cb7935",
                   segment.size = 0.1,size = 1,max.overlaps = 10,
                   fontface="bold", label.size = NA )+
  geom_vline(xintercept = 0, linetype = 2)+
  geom_hline(yintercept = 0, linetype = 2)+
  theme_classic()+
  xlab("∆ Percentage Difference")+
  ylab("Log-Fold Change")

geom_text_repel(
  mapping = NULL,
  data = NULL,
  stat = "identity",
  position = "identity",
  parse = FALSE,
  ...,
  box.padding = 0.25,
  point.padding = 1e-06,
  min.segment.length = 0.5,
  arrow = NULL,
  force = 1,
  force_pull = 1,
  max.time = 0.5,
  max.iter = 10000,
  max.overlaps = getOption("ggrepel.max.overlaps", default = 10),
  nudge_x = 0,
  nudge_y = 0,
  xlim = c(NA, NA),
  ylim = c(NA, NA),
  na.rm = FALSE,
  show.legend = NA,
  direction = c("both", "y", "x"),
  seed = NA,
  verbose = FALSE,
  inherit.aes = TRUE
)
up_genes_GABA_neuron  <-subset(sig_markers.GABAergic_neurons, avg_log2FC > 0)
down_genes_GABA_neuron <-subset(sig_markers.GABAergic_neurons, avg_log2FC < 0)
up_genes_GABAergic_neurons <- sum(sig_markers.GABAergic_neurons$avg_log2FC > 0)
down_genes_GABAergic_neurons <- sum(sig_markers.GABAergic_neurons$avg_log2FC < 0)
print(up_genes_GABAergic_neurons)
print(down_genes_GABAergic_neurons)
################################################################################################### GABA neuron GO- down
downregulated_genes_GABA <- rownames(sig_markers.GABAergic_neurons)[sig_markers.GABAergic_neurons$avg_log2FC < 0]
geneList_downregulated_genes_GABA <- bitr(downregulated_genes_GABA, 
                                          fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# GO
ego_GABA_down <- enrichGO(gene         = geneList_downregulated_genes_GABA$ENTREZID,
                          OrgDb        = org.Mm.eg.db,
                          keyType      = "ENTREZID",
                          ont          = "all", 
                          pAdjustMethod = "BH", 
                          qvalueCutoff = 0.05, 
                          readable     = TRUE)
head(ego_GABA_down@result[, c("pvalue")])
print(ego_GABA_down)
dotplot(ego_GABA_down, showCategory = 10)
#write.table(ego,file="F:/RNA_Seq/N3-N4-H2-H4-MUTAN/enrichment_results-down-GABA-6.csv")
write.csv(ego_GABA_down@result, file = "E:/enrichment_results-down-GABA_BH.csv")

selected_ego <- ego_GABA_down@result[c(177,185,236,405), ]
selected_ego <- new("enrichResult", result = selected_ego)
#write.csv(selected_ego, file = "E:/selected_ego.csv")
color_params <- list(foldChange = geneList_named)
cnetplot(selected_ego, 
         showCategory = 4,
         foldChange = color_params,
         layout = "fr",
         colorEdge = TRUE,
         circular = FALSE,
         node_label = "all",
         color_gene = "grey",
         cex_category = 0.2,
         cex_gene = 0.2,
         node_label_size = NULL,
         cex_label_category = 0.2,
         cex_label_gene = 0.2,
         max.overlaps = 100)
ggrepel::geom_text_repel(aes(label = Hcn1), size = 3, box.padding = 0.5)
p1 <- enrichmentNetwork(selected_ego@result,
                        colorBy = 'pvalue',
                        colorType = c("pval"), 
                        nodeSize = 'setSize',
                        fontSize = 1.5,
                        drawEllipses = TRUE,
                        verbose = TRUE)
p1
p1 + scale_color_continuous_c4a_div('benedictus', reverse = T)
################################################################################################ GABA neuron KEGG-down
kegg_GABA_down <- enrichKEGG(gene = geneList_downregulated_genes_GABA$ENTREZID,
                             organism = 'mmu', 
                             keyType = 'kegg',
                             pAdjustMethod = "BH",
                             qvalueCutoff = 0.05,
                             use_internal_data = FALSE)
eKEGG_down <- setReadable(kegg_GABA_down, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
kegg_terms_down <- eKEGG_down@result
write.csv(kegg_terms_down, file = "E:/kegg_terms_down.csv")
metabolism_terms_down <- kegg_terms_down[grepl("Metabolism", kegg_terms_down$category), ]
target_terms_down <- metabolism_terms[
  grepl("Carbohydrate metabolism|Lipid metabolism", metabolism_terms$subcategory, ignore.case = TRUE),
]
all_genes_down <- unlist(strsplit(target_terms_down$geneID, "/"))
unique_genes_down <- unique(all_genes_down)
print(unique_genes_down)
num_unique_genes_down <- length(unique_genes_down)
#print(paste("Unique genes in Carbohydrate & Lipid metabolism pathways:", num_unique_genes_down))
################################################################################################## GABA neuron GO-up
upregulated_genes_GABA <- rownames(sig_markers.GABAergic_neurons)[sig_markers.GABAergic_neurons$avg_log2FC > 0]
geneList_upregulated_genes_GABA <- bitr(upregulated_genes_GABA, 
                                        fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Mm.eg.db)

# GO
ego_GABA_up <- enrichGO(gene         = geneList_upregulated_genes_GABA$ENTREZID,
                        OrgDb        = org.Mm.eg.db,
                        keyType      = "ENTREZID",
                        ont          = "all", 
                        pAdjustMethod = "BH", 
                        qvalueCutoff = 0.05, 
                        readable     = TRUE)


head(ego_GABA_up@result[, c("pvalue")])
print(ego_GABA_up)
dotplot(ego_GABA_up, showCategory = 10)
write.csv(ego_GABA_up@result, file = "E:/enrichment_results-up-GABA_BH.csv")
##################################################################################################      GABA KEGG-up
kegg_GABA_up <- enrichKEGG(gene = geneList_upregulated_genes_GABA$ENTREZID,
                           organism = 'mmu', 
                           keyType = 'kegg',
                           pAdjustMethod = "BH",
                           qvalueCutoff = 0.05,
                           
                           use_internal_data = FALSE)

eKEGG_up <- setReadable(kegg_GABA_up, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
kegg_terms_up <- eKEGG_up@result
write.csv(kegg_terms_up, file = "E:/GABA-UP/kegg_terms_up.csv")
metabolism_terms_up <- kegg_terms_up[grepl("Metabolism", kegg_terms_up$category), ]
target_terms_up <- metabolism_terms_up[
  grepl("Carbohydrate metabolism|Lipid metabolism", metabolism_terms_up$subcategory, ignore.case = TRUE),
]
all_genes_up <- unlist(strsplit(target_terms_up$geneID, "/"))
unique_genes_up <- unique(all_genes_up)
print(unique_genes_up)
num_unique_genes_up <- length(unique_genes_up)
print(paste("Unique genes in Carbohydrate & Lipid metabolism pathways:", num_unique_genes_up))
head(kegg_GABA_up@result[, c("pvalue")])
print(kegg_GABA_up)
dotplot(kegg_GABA_up, showCategory = 10)
write.csv(kegg_GABA_up@result, file = "E:/kegg_results-kegg_GABA_up-5.csv")
eKEGG_up <- setReadable(kegg_GABA_up, OrgDb = org.Mm.eg.db, keyType="ENTREZID")
write.csv(eKEGG_up@result,file = 'E:/DEGs_KEGG_up.CSV')
##################################################################################################################HCN1
Idents(object = GABAergic_neurons) <- "seurat_clusters"
Gad2.positive.cells <- subset(GABAergic_neurons, subset =  Gad2 > 0)
Gad2.negative.cells <- subset(GABAergic_neurons, subset = Gad2 == 0)
Hcn1.positive.cells <- subset(GABAergic_neurons, subset = Hcn1 > 0)
Hcn1.Gad2 <- subset(Hcn1.positive.cells, subset = Hcn1 > 0)
fea = c("Gad2")
fea = c("Hcn1","Hcn2","Hcn3","Hcn4")
dot_plot_co_expression_data <- DotPlot(Gad2.positive.cells, features = fea) + 
  theme_minimal() + 
  theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
dot_plot_co_expression_data
co_expression_data <- dot_plot_co_expression_data$data
Hcn2.negative.cells <- subset(GABAergic_neurons, subset = Hcn2 > 0)
Hcn3.positive.cells <- subset(GABAergic_neurons, subset = Hcn3 > 0)
Hcn4.negative.cells <- subset(GABAergic_neurons, subset = Hcn4 > 0)

Gad2_expression <- FetchData(GABAergic_neurons, vars = "Gad2")  
GABAergic_neurons$Gad2_group <- ifelse(Gad2_expression > 0, "Expressing", "Not Expressing")   
DotPlot(object = GABAergic_neurons, features = "Gad2", group.by = "Gad2_group")
DotPlot(object = GABAergic_neurons,features = fea,group.by = "Gad2_group", col.min=-1,col.max = 5,scale = F,,
        scale.min = 0,scale.max=100,cols =  c("lightgrey", "red"))+ 
  theme(axis.text.x = element_text(angle=0, hjust=.5, vjust=.5))+coord_flip()
dot_plot <- DotPlot(object = GABAergic_neurons, features = fea, group.by = "Gad2_group",
                    col.min = -1, col.max = 3, scale = FALSE,
                    scale.min = 0, scale.max = 100, cols = c("lightgrey", "red")) +  # 设置标尺最大值为5
  theme(axis.text.x = element_text(angle = 0, hjust = 0.5, vjust = 0.5)) +
  coord_flip()
dot_plot
dot_plot_data <- dot_plot$data
percent_expressed <- dot_plot_data$pct.exp
average_expression <- dot_plot_data$avg.exp