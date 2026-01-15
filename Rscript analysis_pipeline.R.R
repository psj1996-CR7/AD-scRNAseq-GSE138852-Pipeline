############################################################
# scRNA-seq analysis pipeline for GSE138852
# Author: PRASHANTH S JAVALI
# Requirements: R >= 4.1, Seurat v4, Monocle3
############################################################

############################
# 1. SETUP & REQUISITES
############################
set.seed(1234)

# Load libraries
suppressPackageStartupMessages({
  library(GEOquery)
  library(Seurat)
  library(dplyr)
  library(ggplot2)
  library(Matrix)
  library(patchwork)
  library(clusterProfiler)
  library(msigdbr)
  library(SingleR)
  library(celldex)
  library(ComplexHeatmap)
  library(tidyverse)
  library(circlize)
  library(enrichplot)
  library(monocle3)
  library(SeuratWrappers)
})

# Create directory structure
outdir <- "./outputs"
dirs <- c(file.path(outdir, "figures"), 
          file.path(outdir, "tables"), 
          file.path(outdir, "objects"))
lapply(dirs, dir.create, showWarnings = FALSE, recursive = TRUE)

figdir <- dirs[1]
tabdir <- dirs[2]
objdir <- dirs[3]

############################
# 2. DATA LOADING & OBJECT CREATION
############################

# Load raw counts and metadata
counts <- read.csv("GSE138852_counts.csv.gz", row.names = 1, check.names = FALSE)
meta   <- read.csv("GSE138852_covariates.csv.gz", row.names = 1, check.names = FALSE)

# Ensure cell names match between matrices
common_cells <- intersect(colnames(counts), rownames(meta))
counts <- as.matrix(counts[, common_cells])
meta   <- meta[common_cells, , drop = FALSE]

# Initialize Seurat Object
seu <- CreateSeuratObject(
  counts = counts,
  meta.data = meta,
  project = "GSE138852",
  min.cells = 3,
  min.features = 200
)

############################
# 3. QUALITY CONTROL & FILTERING
############################

# Calculate Mitochondrial percentage
seu[["percent.mt"]] <- PercentageFeatureSet(seu, pattern = "^MT-")

# Visualize QC metrics
qc_plot <- VlnPlot(seu, features = c("nFeature_RNA", "nCount_RNA", "percent.mt"), 
                   ncol = 3, pt.size = 0.1)
ggsave(file.path(figdir, "QC_violin.jpeg"), qc_plot, width = 12, height = 4, dpi = 600)

# Apply filtering thresholds
seu <- subset(
  seu,
  subset = nFeature_RNA > 200 &
    nFeature_RNA < 1600 &
    nCount_RNA < 3000 &
    percent.mt < 10
)

############################
# 4. NORMALIZATION & DIM REDUCTION
############################

# SCTransform (includes normalization and scaling)
seu <- SCTransform(seu, vars.to.regress = "percent.mt", verbose = FALSE)

# Run PCA, UMAP and Clustering
seu <- RunPCA(seu, verbose = FALSE)
seu <- RunUMAP(seu, dims = 1:30)
seu <- FindNeighbors(seu, dims = 1:30)
seu <- FindClusters(seu, resolution = 0.5)

# Save Cluster UMAP
umap_clusters <- DimPlot(seu, reduction = "umap", label = TRUE) +
  theme_classic() + ggtitle("GSE138852 – Clusters")
ggsave(file.path(figdir, "UMAP_clusters.jpeg"), umap_clusters, width = 6, height = 5, dpi = 600)

############################
# 5. AUTOMATED CELL TYPE ANNOTATION
############################

ref <- celldex::HumanPrimaryCellAtlasData()
pred <- SingleR(test = GetAssayData(seu, slot = "data"), ref = ref, labels = ref$label.main)
seu$celltype <- pred$labels

# Save CellType UMAP
umap_celltypes <- DimPlot(seu, group.by = "celltype", label = TRUE, repel = TRUE) +
  theme_classic()
ggsave(file.path(figdir, "UMAP_celltypes.jpeg"), umap_celltypes, width = 7, height = 6, dpi = 600)

############################
# 6. DIFFERENTIAL EXPRESSION (DE)
############################

# A. Find Markers for all clusters
markers <- FindAllMarkers(seu, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
write.csv(markers, file.path(tabdir, "cluster_markers.csv"), row.names = FALSE)

# B. AD vs Control (Global)
Idents(seu) <- "oupSample.batchCond"
deg_AD_vs_Control <- FindMarkers(seu, ident.1 = "AD", ident.2 = "ct", 
                                 min.pct = 0.25, logfc.threshold = 0.25)
deg_AD_vs_Control$gene <- rownames(deg_AD_vs_Control)
deg_AD_vs_ct <- deg_AD_vs_Control %>% arrange(desc(avg_log2FC))
write.csv(deg_AD_vs_ct, file.path(tabdir, "DEG_AD_vs_ct.csv"), row.names = FALSE)

# C. AD vs Control (Per Cell Type)
deg_list <- list()
for (ct in unique(seu$celltype)) {
  cells_ct <- WhichCells(seu, expression = celltype == ct)
  if (length(cells_ct) < 50) next
  seu_ct <- subset(seu, cells = cells_ct)
  if (length(unique(seu_ct$oupSample.batchCond)) < 2) next
  
  deg <- FindMarkers(seu_ct, ident.1 = "AD", ident.2 = "ct", min.pct = 0.2)
  deg$gene <- rownames(deg)
  deg$celltype <- ct
  deg_list[[ct]] <- deg
}
deg_celltype <- bind_rows(deg_list)
write.csv(deg_celltype, file.path(tabdir, "DEG_AD_vs_ct_by_celltype.csv"), row.names = FALSE)

############################
# 7. EXPRESSION VISUALIZATION
############################

# Top Genes Violin Plots
top_genes <- deg_AD_vs_ct %>% filter(p_val_adj < 0.05) %>% slice_head(n = 6) %>% pull(gene)

vln <- VlnPlot(seu, features = top_genes, group.by = "oupSample.batchCond", pt.size = 0, ncol = 3) &
  theme_classic() & theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave(file.path(figdir, "Violin_AD_vs_ct.jpeg"), vln, width = 10, height = 6, dpi = 600)

# Heatmap of Top Cluster Markers
top_cluster_genes <- markers %>% group_by(cluster) %>% slice_max(order_by = avg_log2FC, n = 5) %>% pull(gene) %>% unique()
avg_expr <- AverageExpression(seu, features = top_cluster_genes, group.by = "celltype", assays = "SCT")$SCT
avg_expr_scaled <- t(scale(t(avg_expr)))
avg_expr_scaled[is.na(avg_expr_scaled)] <- 0

celltype_colors <- structure(circlize::rand_color(ncol(avg_expr_scaled)), names = colnames(avg_expr_scaled))
top_anno <- HeatmapAnnotation(CellType = colnames(avg_expr_scaled), col = list(CellType = celltype_colors))

jpeg(file.path(figdir, "Heatmap_celltype_markers.jpeg"), width = 14, height = 10, units = "in", res = 600)
Heatmap(avg_expr_scaled, name = "Z-score", top_annotation = top_anno, cluster_columns = FALSE)
dev.off()

############################
# 8. GENE SET ENRICHMENT ANALYSIS (GSEA)
############################

# Prepare ranked gene list
gene_list <- deg_AD_vs_ct$avg_log2FC
names(gene_list) <- deg_AD_vs_ct$gene
gene_list <- sort(gene_list, decreasing = TRUE)

# GSEA with GO Biological Processes
m_t2g_go <- msigdbr(species = "Homo sapiens", category = "C5", subcategory = "BP") %>% 
  dplyr::select(gs_name, gene_symbol)

gsea_res <- GSEA(gene_list, TERM2GENE = m_t2g_go, pvalueCutoff = 0.05, verbose = FALSE)
write.csv(as.data.frame(gsea_res), file.path(tabdir, "GSEA_AD_vs_ct_Hallmark.csv"), row.names = FALSE)

# GSEA Visualizations
gsea_dot <- dotplot(gsea_res, showCategory = 15, split = ".sign") + facet_grid(.~.sign) + theme_bw()
ggsave(file.path(figdir, "GSEA_dotplot.jpeg"), gsea_dot, width = 10, height = 8, dpi = 600)

ridge_p <- ridgeplot(gsea_res, showCategory = 10) + theme_classic()
ggsave(file.path(figdir, "GSEA_RidgePlot.jpeg"), ridge_p, width = 10, height = 8, dpi = 600)

cnet_p <- cnetplot(gsea_res, showCategory = 5, foldChange = gene_list)
ggsave(file.path(figdir, "GSEA_CNET_Corrected.jpeg"), cnet_p, width = 12, height = 10, dpi = 600)

# Heatmap for Protein Folding pathway genes
protein_folding_genes <- as.data.frame(gsea_res) %>% 
  filter(ID == "GOBP_PROTEIN_FOLDING") %>% pull(core_enrichment) %>% strsplit("/") %>% unlist()

avg_expr_folding <- AverageExpression(seu, features = protein_folding_genes, group.by = "celltype", assays = "SCT")$SCT
avg_expr_folding_scaled <- t(scale(t(avg_expr_folding)))

jpeg(file.path(figdir, "Heatmap_Protein_Folding_Genes.jpeg"), width = 10, height = 8, units = "in", res = 600)
Heatmap(avg_expr_folding_scaled, name = "Z-score", column_title = "Protein Folding Genes", cluster_columns = TRUE)
dev.off()



########################## Optional but effective Downstream analysis ############################

############################
# 9. TRAJECTORY ANALYSIS (MONOCLE3)
############################

# Conversion and Trajectory Graph Learning
cds <- as.cell_data_set(seu)
cds <- cluster_cells(cds)
cds <- learn_graph(cds)

# To visualize this, you would typically use:
plot_cells(cds, color_cells_by = "celltype", label_groups_by_cluster = FALSE)

############################
# 10. CELL-CELL INTERACTION (CellChat)
############################
# We compare the signaling "landscape" between AD and Control

# 1. Prepare Data for CellChat
# Ensure celltype labels are set as Idents
Idents(seu) <- "celltype"
data.input  <- GetAssayData(seu, assay = "SCT", slot = "data")
labels      <- Idents(seu)
meta        <- seu@meta.data

# 2. Create CellChat object
cellchat <- createCellChat(object = data.input, meta = meta, group.by = "celltype")

# 3. Set the database (Human)
CellChatDB <- CellChatDB.human 
# Use 'Secreted Signaling' for AD-related cytokine/chemokine analysis
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling") 
cellchat@DB <- CellChatDB.use

# 4. Preprocessing
cellchat <- subsetData(cellchat)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat)

# 5. Compute Communication Probability
cellchat <- computeCommunProb(cellchat, raw.use = TRUE)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# 6. Visualization: Interaction Overview

jpeg(file.path(figdir, "CellChat_Interaction_Network.jpeg"), width = 8, height = 8, units = "in", res = 600)
groupSize <- as.numeric(table(cellchat@idents))
netVisual_circle(cellchat@net$count, vertex.weight = groupSize, weight.scale = T, label.edge= F, title.name = "Number of interactions")
dev.off()

# 7. Signaling Role Analysis (identifying dominant senders/receivers)
cellchat <- netAnalysis_computeCentrality(cellchat, slot.name = "netP")
role_plot <- netAnalysis_signalingRole_scatter(cellchat)
ggsave(file.path(figdir, "CellChat_Signaling_Roles.jpeg"), role_plot, width = 10, height = 8)

############################
# 11. GSEA INTERPRETATION LOGIC
############################
# This part filters the GSEA results to identify actionable biological themes

gsea_results_df <- as.data.frame(gsea_res)

# A. Identify Neurodegeneration-related Pathways
neuro_pathways <- gsea_results_df %>%
  filter(grepl("NEURO|ALZHEIMER|SYNAPSE|AXON", ID, ignore.case = TRUE)) %>%
  arrange(p.adjust)

write.csv(neuro_pathways, file.path(tabdir, "GSEA_Neuro_Specific.csv"))

# B. Identify Inflammatory/Immune Response
immune_pathways <- gsea_results_df %>%
  filter(grepl("CYTOKINE|IMMUNE|INFLAMMATORY|INTERFERON", ID, ignore.case = TRUE)) %>%
  arrange(p.adjust)

write.csv(immune_pathways, file.path(tabdir, "GSEA_Immune_Specific.csv"))

# C. Plotting GSEA Running Score for a key Neuro pathway

if(nrow(neuro_pathways) > 0){
  top_neuro <- neuro_pathways$ID[1]
  p_gsea <- gseaplot2(gsea_res, geneSetID = top_neuro, title = top_neuro)
  ggsave(file.path(figdir, paste0("GSEA_Top_Neuro_", top_neuro, ".jpeg")), p_gsea, width = 8, height = 6)
}



############################
# FINAL SAVE
############################
saveRDS(seu, file = file.path(objdir, "GSE138852_final_processed.rds"))


