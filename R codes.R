# Load necessary libraries
library(Seurat)
library(loupeR)
library(ggplot2)
library(patchwork)
library(SeuratDisk)
library(SeuratData)
library(RColorBrewer)
library(CellChat)
library(ggprism)
library(cowplot)
library(openxlsx)
library(RDAVIDWebService)
library(biomaRt)
library(SpatialExperiment)
library(SingleCellExperiment)
library(scater)
library(scran)
library(Matrix)
library(CARD)
library(monocle3)
library(dplyr)
library(ggplotify)
library(plotly)
library(beepr)
library(LSD)
library(ggimage)
library(ggprism)

# Set working directory
setwd("~/Dropbox/01_visium/data_and_code/data")

# Load spatial transcriptomics data
data_dir <- getwd()
seurat_obj <- Load10X_Spatial(
  data.dir = data_dir,
  filename = "filtered_feature_bc_matrix.h5",
  assay = "Spatial",
  slice = "tissue_hires_image",
  filter.matrix = TRUE,
  to.upper = FALSE,
  image = NULL
)

# Preprocessing
seurat_obj <- NormalizeData(seurat_obj, normalization.method = "LogNormalize", scale.factor = 10000)
seurat_obj <- FindVariableFeatures(seurat_obj, selection.method = "vst", nfeatures = 2000)

# Scaling data
all_genes <- rownames(seurat_obj)
seurat_obj <- ScaleData(seurat_obj, features = all_genes)

# Adjusting scale factors for image resolution
current_scale_factor <- seurat_obj@images[["tissue_hires_image"]]@scale.factors[["lowres"]]
old_dpi <- 72  # Original image DPI
new_dpi <- 300  # Adjusted image DPI

# Calculate new scale factor
new_scale_factor <- current_scale_factor * (new_dpi / old_dpi)
print(new_scale_factor)

# Update scale factors in Seurat object
seurat_obj@images[["tissue_hires_image"]]@scale.factors[["lowres"]] <- new_scale_factor

# Plot variable features
top10 <- head(VariableFeatures(seurat_obj), 10)
plot1 <- VariableFeaturePlot(seurat_obj, pt.size = 0.5)
plot2 <- LabelPoints(plot = plot1, points = top10, repel = TRUE)
plot2

# PCA and clustering
seurat_obj <- RunPCA(seurat_obj, features = VariableFeatures(object = seurat_obj))
seurat_obj <- FindNeighbors(seurat_obj, dims = 1:10)
seurat_obj <- FindClusters(seurat_obj, resolution = c(1.0, 0.8, 0.6, 0.4, 0.2, 0.1))
seurat_obj <- RunUMAP(seurat_obj, dims = 1:10)

# UMAP visualization
umap1 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "Spatial_snn_res.1")
umap2 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "Spatial_snn_res.0.6")
umap3 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "Spatial_snn_res.0.2")
umap4 <- DimPlot(seurat_obj, reduction = "umap", pt.size = 0.2, group.by = "Spatial_snn_res.0.1")
(umap1 + umap2) / (umap3 + umap4)
umap2

# Adjust cluster labels
seurat_obj$Spatial_snn_res.0.6[seurat_obj$Spatial_snn_res.0.6 == 7] <- 5
seurat_obj$Spatial_snn_res.1 <- as.character(seurat_obj$Spatial_snn_res.1)
seurat_obj$Spatial_snn_res.1[seurat_obj$Spatial_snn_res.0.6 == 0] <- "ZG"
seurat_obj$Spatial_snn_res.1[seurat_obj$Spatial_snn_res.0.6 == 1] <- "CT"
seurat_obj$Spatial_snn_res.1[seurat_obj$Spatial_snn_res.0.6 == 2] <- "ZF"
seurat_obj$Spatial_snn_res.1[seurat_obj$Spatial_snn_res.0.6 == 4] <- "ZX"
seurat_obj$Spatial_snn_res.1[seurat_obj$Spatial_snn_res.0.6 == 5] <- "Medulla"
seurat_obj$Spatial_snn_res.1[seurat_obj$Spatial_snn_res.0.6 == 3] <- "WAT"
seurat_obj$Spatial_snn_res.1[seurat_obj$Spatial_snn_res.0.6 == 6] <- "BAT"
seurat_obj$Spatial_snn_res.1 <- factor(seurat_obj$Spatial_snn_res.1, levels = c("ZG", "ZF", "ZX", "Medulla", "CT", "WAT", "BAT"))

# Identify markers
Idents(seurat_obj) <- seurat_obj$Spatial_snn_res.1
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_filtered <- markers[markers$avg_log2FC > 1 & markers$p_val_adj < 0.05, ]

# Select top 10 markers per cluster
markers_top10 <- markers_filtered %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Generate heatmap
top_genes <- unique(markers_top10$gene)
DoHeatmap(seurat_obj, features = top_genes, size = 4)

# Save Seurat object
SaveH5Seurat(seurat_obj, overwrite = TRUE)

# Visualization
library(RColorBrewer)
num_clusters <- length(unique(seurat_obj$Spatial_snn_res.0.6))
color_codes <- brewer.pal(num_clusters, "Set1")
p1 <- DimPlot(seurat_obj, reduction = "umap", label = TRUE)
p2 <- SpatialDimPlot(seurat_obj, label = TRUE, label.size = 3, pt.size.factor = 1.3, images = "tissue_hires_image")
p1 + p2

# Re-identify markers with adjusted clusters
Idents(seurat_obj) <- seurat_obj$Spatial_snn_res.0.6
markers <- FindAllMarkers(seurat_obj, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25)
markers_top10 <- markers %>%
  group_by(cluster) %>%
  top_n(n = 10, wt = avg_log2FC)

# Generate heatmap with top markers
top_genes <- unique(markers_top10$gene)
pdf("heatmap_seurat5.pdf", width = 9, height = 9)
DoHeatmap(seurat_obj, features = top_genes, size = 4, raster = FALSE)
dev.off()

# Save marker genes to Excel
wb <- createWorkbook()
clusters <- unique(markers_filtered$cluster)
for (i in seq_along(clusters)) {
  cluster_id <- clusters[i]
  addWorksheet(wb, paste0("Cluster_", cluster_id))
  cluster_markers <- markers_filtered[markers_filtered$cluster == cluster_id, ]
  cluster_markers <- cluster_markers[order(cluster_markers$p_val_adj), ]
  writeData(wb, sheet = i, cluster_markers, rowNames = FALSE)
}
saveWorkbook(wb, "Marker_Genes.xlsx", overwrite = TRUE)

# Monocle3 Trajectory Analysis ------------------------------------------

# Subset Seurat object if necessary (e.g., specific clusters)
# Here, we're using the existing 'seurat_obj' or you can subset it as per your analysis
# For example, to focus on certain clusters:
# subset_seurat <- subset(seurat_obj, idents = c("Cluster1", "Cluster2"))

# Convert Seurat object to CellDataSet for Monocle3
library(monocle3)
cds <- as.cell_data_set(seurat_obj)

# Add gene short names if missing
fData(cds)$gene_short_name <- rownames(fData(cds))

# Assign partitions (here assuming one partition)
recreate_partition <- rep(1, length(cds@colData@rownames))
names(recreate_partition) <- cds@colData@rownames
recreate_partition <- as.factor(recreate_partition)
cds@clusters@listData[["UMAP"]][["partitions"]] <- recreate_partition

# Use Seurat clusters
list_cluster <- Idents(seurat_obj)
cds@clusters@listData[["UMAP"]][["clusters"]] <- list_cluster

# Assign UMAP embeddings
cds@int_colData@listData$reducedDims$UMAP <- seurat_obj@reductions$umap@cell.embeddings

# Plot before trajectory inference
cluster_before_trajectory <- plot_cells(cds,
                                        label_groups_by_cluster = FALSE,
                                        group_label_size = 5) +
  theme(legend.position = "right")

cluster_names <- plot_cells(cds,
                            color_cells_by = "Spatial_snn_res.1",
                            label_groups_by_cluster = FALSE,
                            group_label_size = 5) +
  theme(legend.position = "right")

cluster_before_trajectory | cluster_names

# Learn the trajectory graph
cds <- learn_graph(cds, use_partition = FALSE)

plot_cells(cds,
           color_cells_by = 'Spatial_snn_res.1',
           label_groups_by_cluster = TRUE,
           label_branch_points = FALSE,
           label_roots = TRUE,
           label_leaves = FALSE,
           group_label_size = 5)

# Order cells in pseudotime
root_cells <- colnames(cds[, clusters(cds) == 1])  # Adjust cluster number as needed
cds <- order_cells(cds, reduction_method = 'UMAP', root_cells = root_cells)

# Plot pseudotime
pseudotime_plot <- plot_cells(cds,
                              reduction_method = "UMAP",
                              color_cells_by = 'pseudotime',
                              label_groups_by_cluster = FALSE,
                              label_branch_points = FALSE,
                              label_roots = FALSE,
                              label_leaves = FALSE,
                              cell_size = 0.85) +
  theme_prism(base_size = 9)
pseudotime_plot

# Add pseudotime to Seurat object
seurat_obj$pseudotime <- pseudotime(cds)

# Visualize pseudotime in Seurat
FeaturePlot(seurat_obj, features = "pseudotime", label = TRUE) +
  theme_prism(base_size = 9)

# Identify genes that change over pseudotime
deg_pseudotime <- graph_test(cds, neighbor_graph = 'principal_graph', cores = 8)
deg_pseudotime <- deg_pseudotime %>%
  arrange(q_value) %>%
  filter(status == 'OK')

# Select top genes (adjust as needed)
top_genes_pseudotime <- head(deg_pseudotime$gene_short_name, 10)

# Plot expression of top genes over pseudotime
FeaturePlot(seurat_obj, features = top_genes_pseudotime) +
  theme_prism(base_size = 9)

# Combine plots
library(cowplot)
combined_plot <- cowplot::plot_grid(
  cluster_before_trajectory + theme_prism(base_size = 9),
  pseudotime_plot,
  ncol = 1,
  labels = c("A", "B")
)
ggsave(combined_plot, file = "trajectory_analysis.pdf", width = 10, height = 8)

# End of Monocle3 Trajectory Analysis ------------------------------------

# DAVID analysis
library(RDAVIDWebService)
library(biomaRt)
library(openxlsx)

# Convert gene symbols to Entrez IDs
mart <- useMart("ensembl", dataset = "mmusculus_gene_ensembl")
genes <- markers_filtered$gene
genes <- getLDS(attributes = c("mgi_symbol"), filters = "mgi_symbol", values = genes, mart = mart,
                attributesL = c("mgi_symbol", "entrezgene_id"), martL = mart, uniqueRows = TRUE)
markers_annotated <- merge(markers_filtered, genes, by.x = "gene", by.y = "MGI.symbol", all = TRUE)
markers_annotated$geneid <- markers_annotated$NCBI.gene..formerly.Entrezgene..ID

# Perform DAVID functional annotation
email <- "your_email@example.com"  # Replace with your email
david <- DAVIDWebService(email = email)
annotation_categories <- c("GOTERM_BP_DIRECT", "GOTERM_MF_DIRECT", "GOTERM_CC_DIRECT", "KEGG_PATHWAY")

results_list <- list()
for (category in annotation_categories) {
  setAnnotationCategories(david, category)
  clusters <- unique(markers_annotated$cluster)
  cluster_results <- list()
  for (cluster_id in clusters) {
    genes_in_cluster <- markers_annotated$geneid[markers_annotated$cluster == cluster_id]
    genes_in_cluster <- na.omit(genes_in_cluster)
    if (length(genes_in_cluster) == 0) {
      genes_in_cluster <- c(1)
    }
    addList(david, genes_in_cluster, idType = "ENTREZ_GENE_ID", listName = paste0("Cluster_", cluster_id), listType = "Gene")
    result <- getFunctionalAnnotationChart(david)
    if (nrow(result) == 0) {
      result <- data.frame(Term = "No GO term", Count = 1, Genes = NA)
    }
    cluster_results[[as.character(cluster_id)]] <- result
  }
  results_list[[category]] <- cluster_results
}

# Save DAVID results to Excel
for (category in names(results_list)) {
  wb <- createWorkbook()
  for (cluster_id in names(results_list[[category]])) {
    addWorksheet(wb, paste0("Cluster_", cluster_id))
    writeData(wb, sheet = paste0("Cluster_", cluster_id), results_list[[category]][[cluster_id]], rowNames = FALSE)
  }
  saveWorkbook(wb, paste0("DAVID_", category, ".xlsx"), overwrite = TRUE)
}

# Plot DAVID results
library(ggplot2)
library(ggprism)
for (category in names(results_list)) {
  combined_results <- do.call(rbind, lapply(results_list[[category]], function(df) {
    df$Cluster <- names(df)
    return(df)
  }))
  ggplot(combined_results, aes(x = Cluster, y = Term, size = Count)) +
    geom_point(color = "#377EB8", alpha = -log2(combined_results$PValue)) +
    labs(x = "", y = "GO terms", title = category, size = "Gene Count") +
    theme_prism(base_size = 10) +
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

# Save marker genes to Excel
wb <- createWorkbook()
clusters <- unique(markers_filtered$cluster)
for (i in seq_along(clusters)) {
  cluster_id <- clusters[i]
  addWorksheet(wb, paste0("Cluster_", cluster_id))
  cluster_markers <- markers_filtered[markers_filtered$cluster == cluster_id, ]
  cluster_markers <- cluster_markers[order(cluster_markers$p_val_adj), ]
  writeData(wb, sheet = i, cluster_markers, rowNames = FALSE)
}
saveWorkbook(wb, "Marker_Genes.xlsx", overwrite = TRUE)

# CellChat analysis
library(CellChat)
cellchat <- createCellChat(object = seurat_obj, group.by = "Spatial_snn_res.0.6", datatype = "spatial")
CellChatDB <- CellChatDB.mouse
CellChatDB.use <- subsetDB(CellChatDB, search = "Secreted Signaling", key = "annotation")
cellchat@DB <- CellChatDB.use
cellchat <- subsetData(cellchat)
future::plan("multisession", workers = 4)
cellchat <- identifyOverExpressedGenes(cellchat)
cellchat <- identifyOverExpressedInteractions(cellchat, variable.both = TRUE)
cellchat <- computeCommunProb(cellchat)
cellchat <- filterCommunication(cellchat, min.cells = 10)
cellchat <- computeCommunProbPathway(cellchat)
cellchat <- aggregateNet(cellchat)

# Save CellChat object
save(cellchat, file = "cellchat.RData")

# Visualize interaction networks
group_size <- as.numeric(table(cellchat@idents))
par(mfrow = c(1, 2), xpd = TRUE)
netVisual_circle(cellchat@net$count, vertex.weight = rowSums(cellchat@net$count), weight.scale = TRUE, label.edge = FALSE, title.name = "Number of interactions")
netVisual_circle(cellchat@net$weight, vertex.weight = rowSums(cellchat@net$weight), weight.scale = TRUE, label.edge = FALSE, title.name = "Interaction weights/strength")

# Generate plots for selected pathways
pathways <- unique(cellchat@netP$pathway_name)
for (pathway in pathways) {
  netVisual_aggregate(cellchat, signaling = pathway, layout = "circle", title.space = 3, vertex.label.cex = 0.8, pt.title = 2)
  spatialFeaturePlot(cellchat, signaling = pathway, point.size = 0.5)
}

# Save trajectory analysis plots
# (Assuming that trajectory analysis code is included elsewhere)

# Deconvolution using CARD
library(CARD)
# Load marker list and spatial data
marker_list <- read.xlsx("marker_list.xlsx", sheet = 1)
spatial_count <- assay(spe)
spatial_location <- data.frame(x = spe@colData$imagecol, y = spe@colData$imagerow)

# Create CARD object
CARD_obj <- createCARDfreeObject(
  markerList = marker_list,
  spatial_count = spatial_count,
  spatial_location = spatial_location,
  minCountGene = 200,
  minCountSpot = 20
)

# Run CARD
CARD_obj <- CARD_refFree(CARD_obj)

# Visualization
colors <- brewer.pal(9, "Set1")
pie_plot <- CARD.visualize.pie(CARD_obj@Proportion_CARD, CARD_obj@spatial_location, colors = colors, radius = 180)

# Save plots
ggsave(pie_plot, filename = "CARD_pie_chart.pdf", width = 10, height = 8)
