saveRDS(thesis_seu, "GBM_single_cell_analysis/outputs/objects/thesis_seu.RDS")
saveRDS(Lseu, "GBM_single_cell_analysis/outputs/objects/Lseu.RDS")
saveRDS(Mseu, "GBM_single_cell_analysis/outputs/objects/Mseu.RDS")
saveRDS(cd4s, "GBM_single_cell_analysis/outputs/objects/cd4s.RDS")
saveRDS(cd8s, "GBM_single_cell_analysis/outputs/objects/cd8s.RDS")
# Load libraries ####

library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(cowplot)
library(SingleR)
library(scRepertoire)
library(openxlsx)
library(CelliD)
library(ggpubr)
library(HGNChelper)
library(clustree)
library(patchwork)
library(tibble)
library(EnhancedVolcano)
library(ComplexHeatmap)
library(pheatmap)
library(DESeq2)
library(scCustomize)
library(lme4)
library(lmerTest)
library(SingleCellExperiment)
library(slingshot)
library(RColorBrewer)
library(viridis)
library(clusterProfiler)
library(org.Hs.eg.db)

set.seed(2025)

# Read in data

thesis_seu <- readRDS("GBM_single_cell_analysis/outputs/objects/thesis_seu.RDS")


elbow <- ElbowPlot(thesis_seu, reduction = "pca", ndims = 50) +
  ggtitle("Percentage variance by principal component") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12,face = "bold")
  )  +
  geom_vline(xintercept = 25.5, color = "red", linetype = "dashed", linewidth = 0.5)

ggsave("GBM_single_cell_analysis/outputs/clustering/elbow_plot.png",
       elbow,
       height = 6, width = 6)

# Further automatic annotation ####
data("HgProteinCodingGenes")

thesis_annot <- JoinLayers(thesis_seu)
fvar <- head(VariableFeatures(thesis_annot), 2000)
thesis_annot <- subset(thesis_annot, features = fvar)
thesis_annot <- subset(thesis_annot, features = intersect(HgProteinCodingGenes, rownames(thesis_annot)))
rm(HgProteinCodingGenes, fvar)
gc()
options(future.globals.maxSize = 7 * 1024^3)
thesis_annot <- RunMCA(thesis_annot)

panglao <- read_tsv("GBM_single_cell_analysis/outputs/clustering/PanglaoDB_markers.tsv")

panglao <- panglao %>%
  filter(str_detect(species, "Hs")) %>%
  group_by(`cell type`) %>%  
  summarise(geneset = list(`official gene symbol`))

gs <- setNames(panglao$geneset, panglao$`cell type`)
gs <- gs[sapply(gs, length) >= 10]
HGT_gs <- RunCellHGT(thesis_annot, pathways = gs, dims = 1:50, n.features = 200)
gs_prediction <- rownames(HGT_gs)[apply(HGT_gs, 2, which.max)]
gs_prediction_signif <- ifelse(apply(HGT_gs, 2, max)>2, yes = gs_prediction, "unassigned")

table(gs_prediction_signif)
thesis_seu$CelliD.labels <- gs_prediction_signif

Idents(thesis_seu) <- thesis_seu$harm_clusters_pc35_res3

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R"); 
thesis_seu <- run_sctype(thesis_seu, assay = "RNA", scaled = TRUE, known_tissue_type=c("Brain", "Immune system"), custom_marker_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",name="ScType.labels")

rm(gs, HGT_gs, gs_prediction, gs_prediction_signif, panglao, thesis_annot)

GBMap_3 <- read.xlsx("GBM_single_cell_analysis/outputs/clustering/GBMap_3.xlsx", colNames=F)
colnames(GBMap_3) <- c("cellName", "marker")

GBMap_3 <- GBMap_3 %>%
  group_by(cellName) %>%
  summarise(geneSymbolmore1 = paste(marker, collapse = ","), .groups = "drop")

GBMap_3$tissueType <- "Glioblastoma"
GBMap_3$geneSymbolmore2 <- NA

GBMap_3 <- GBMap_3 %>%
  select(tissueType, cellName, geneSymbolmore1, geneSymbolmore2)

write.xlsx(GBMap_3, "GBM_single_cell_analysis/outputs/clustering/GBMap_3DB.xlsx")

GBMap_4 <- read.xlsx("GBM_single_cell_analysis/outputs/clustering/GBMap_4.xlsx", colNames=F)
colnames(GBMap_4) <- c("cellName", "marker")

GBMap_4 <- GBMap_4 %>%
  group_by(cellName) %>%
  summarise(geneSymbolmore1 = paste(marker, collapse = ","), .groups = "drop")

GBMap_4$tissueType <- "Glioblastoma"
GBMap_4$geneSymbolmore2 <- NA

GBMap_4 <- GBMap_4 %>%
  select(tissueType, cellName, geneSymbolmore1, geneSymbolmore2)

write.xlsx(GBMap_4, "GBM_single_cell_analysis/outputs/clustering/GBMap_4DB.xlsx")

rm(GBMap_3, GBMap_4)

thesis_seu <- run_sctype(thesis_seu, assay = "RNA", scaled = TRUE, known_tissue_type = "Glioblastoma", custom_marker_file = "GBM_single_cell_analysis/outputs/clustering/GBMap_3DB.xlsx", name = "GBMap_3.labels")

thesis_seu <- run_sctype(thesis_seu, assay = "RNA", scaled = TRUE, known_tissue_type = "Glioblastoma", custom_marker_file = "GBM_single_cell_analysis/outputs/clustering/GBMap_4DB.xlsx", name = "GBMap_4.labels")

thesis_seu@meta.data <- thesis_seu@meta.data %>%
  select(1:which(names(thesis_seu@meta.data) == "SingleR.labels"), CelliD.labels, ScType.labels, GBMap_3.labels, GBMap_4.labels, everything())

thesis_seu@meta.data <- thesis_seu@meta.data %>%
  select(-c(harm_clusters_pc30_res1, harm_clusters_pc30_res2, harm_clusters_pc30_res3, harm_clusters_pc40_res2,
            harm_clusters_pc40_res2.5, harm_clusters_pc35_res2, harm_clusters_pc35_res2.5))
  

# Primary clustering ####

thesis_seu <- JoinLayers(thesis_seu)

pcs <- seq(15, 40, by = 5)
resolutions <- seq(0.2, 1.2, by = 0.2)

for (pc in pcs){
  cat("Constructing SNN graph,", pc, "dimensions", "\n")
  thesis_seu <- FindNeighbors(thesis_seu,
                              reduction = "harmony",
                              dims = 1:pc,
                              graph.name = paste0("harmony_snn_pc", pc))
  cat("Generating UMAP,", pc, "dimensions", "\n")
  thesis_seu <- RunUMAP(thesis_seu,
                        reduction = "harmony",
                        dims = 1:pc,
                        reduction.name = paste0("umap.harm_pc", pc),
                        graph.name = paste0("harmony_snn_pc", pc))
  
  for (res in resolutions){
    cat("Clustering,", pc, "dimensions, resolution", res, "\n")
    
    plots_pdf <- list()
    
    output_dir <- paste0("GBM_single_cell_analysis/outputs/clustering/pc", pc, "/pc", pc, "_res", res, "/")
    dir.create(output_dir, recursive=T, showWarnings=F)
    
    thesis_seu <- FindClusters(thesis_seu,
                               resolution = res,
                               graph.name = paste0("harmony_snn_pc", pc),
                               cluster.name = paste0("harm_clusters_pc", pc, "_res", res))
    
    p <- DimPlot(thesis_seu,
                 reduction = paste0("umap.harm_pc", pc),
                 group.by = paste0("harm_clusters_pc", pc, "_res", res),
                 label = TRUE, label.size = 4.5) +
      ggtitle(paste0("Clusters (pc", pc, "_res", res, ")")) +
      labs(x = "umap_1", y = "umap_2") +
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
    
    plots_pdf[[1]] <- p
    
    ggsave(paste0(output_dir, "harm_clusters_pc", pc, "_res", res, ".png"),
           p,
           height = 8.5, width = 10)
    
    cat("Plotting annotations...", "\n")
    
    annotations <- c("SingleR.labels", "CelliD.labels", "ScType.labels", "GBMap_3.labels", "GBMap_4.labels")
    
    for (annot in seq_along(annotations)){
      p <- DimPlot(thesis_seu,
                   reduction = paste0("umap.harm_pc", pc),
                   group.by = paste0(annotations[[annot]]),
                   label = TRUE, label.size = 3) +
        ggtitle(paste0(annotations[[annot]], " (pc", pc, "_res", res, ")")) +
        labs(x = "umap_1", y = "umap_2") +
        guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
      
      plots_pdf[[(annot + 1)]] <- p
      
      ggsave(paste0(output_dir, annotations[[annot]], "_pc", pc, "_res", res, ".png"),
             p,
             height = 8.5, width = 11)
      
    }
    
    cat("Running differential gene expression analysis for generated clusters,", pc, "dimensions, resolution", res, "\n")
    
    markers <- FindAllMarkers(thesis_seu, 
                              logfc.threshold = 0.5, 
                              min.pct = 0.1, 
                              only.pos = TRUE)
    
    clusterMarkers <- markers %>%
      mutate(
        diff_pct = (pct.1 - pct.2) * 100,
        signif = -log10(p_val_adj),
        signif = ifelse(is.infinite(signif) | signif > 300, 300, signif)) %>%
      group_by(cluster) %>%
      arrange(cluster, desc(diff_pct)) %>%
      slice_head(n = 20) %>%
      ungroup()
    
    m.data <- thesis_seu@meta.data %>%
      arrange(seurat_clusters)
    
    clusters <- unique(m.data$seurat_clusters)
    cluster_plotlist <- list()
    
    cat("Plotting clusters...", "\n")
    
    for (i in seq_along(clusters)){
      cluster_id <- clusters[[i]]
      cluster_cells <- length(which(thesis_seu$seurat_clusters == cluster_id))
      cluster_size <- round(cluster_cells / nrow(thesis_seu@meta.data) * 100, 2)
      
      topmarkers <- clusterMarkers %>%
        filter(cluster == cluster_id) %>%
        arrange(desc(diff_pct))
      
      p <- ggplot(topmarkers, aes(x = diff_pct, y = reorder(gene, diff_pct), size = signif)) +
        # Trailing lines from y-axis to dots
        geom_segment(aes(x = 0, xend = diff_pct, y = reorder(gene, diff_pct), yend = reorder(gene, diff_pct)),
                     color = "gray", size = 0.5, alpha = 0.5) +  # Line settings
        
        # Dot plot
        geom_point(aes(color = avg_log2FC, x = diff_pct, y = gene), alpha = 1) +  
        labs(title = paste0(cluster_id, "; n = ", cluster_cells, "; ", cluster_size, "%")) +
        scale_color_gradient(
          low = "#132132", 
          high = "steelblue1", 
          limits = range(clusterMarkers$avg_log2FC, na.rm = TRUE),
          name = "Fold change"
        ) +
        scale_size_continuous(name = "Significance",
                              range = c(0.05, 2),
                              limits = range(clusterMarkers$signif, na.rm = TRUE)) +  # Adjust dot size range
        xlim(0, 100) +
        theme_cowplot() +
        theme(axis.text.y = element_text(size = 4.5),
              axis.text.x = element_text(size = 5.5),
              axis.title.x = element_blank(),
              axis.title.y = element_blank(),
              plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
              legend.title = element_text(size = 5, face = "bold"),
              legend.text = element_text(size = 5),
              legend.position = "none",
              plot.background = element_rect(fill = "white", color = NA))
      
      cluster_plotlist[[i]] <- p
      
    }
    
    legend_plot <- ggplot(topmarkers, aes(x = diff_pct, y = reorder(gene, diff_pct), size = signif)) +
      geom_segment(aes(x = 0, xend = diff_pct, y = reorder(gene, diff_pct), yend = reorder(gene, diff_pct)),
                   color = "gray", size = 0.5, alpha = 0.5) +
      geom_point(aes(color = avg_log2FC, x = diff_pct, y = gene), alpha = 1) +  
      labs(title = paste0(cluster_id, "; n = ", cluster_cells, "; ", cluster_size, "%")) +
      scale_color_gradient(
        low = "#132132", 
        high = "steelblue1", 
        limits = range(clusterMarkers$avg_log2FC, na.rm = TRUE),
        name = "Fold change"
      ) +
      scale_size_continuous(name = "Significance",
                            range = c(0.2, 3),
                            limits = range(clusterMarkers$signif, na.rm = TRUE)) +  # Adjust dot size range
      xlim(0, 100) +
      theme_cowplot() +
      theme(axis.text = element_text(size = 5), 
            axis.title = element_text(size = 7, face = "bold"),
            plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
            legend.title = element_text(size = 10, face = "bold"),
            legend.text = element_text(size = 10),
            legend.position = "right",
            plot.background = element_rect(fill = "white", color = NA))
    
    legend <- get_legend(legend_plot)
    
    cluster_grid <- plot_grid(plotlist = cluster_plotlist, ncol = 5)
    plot <- plot_grid(cluster_grid, legend,
                      rel_widths = c(5,1), nrow = 1, align = "vh")
    final_plot <- ggdraw(plot) +
      draw_label("%(cluster) - %(rest)", x = 0.5, y = 0.01, vjust = 0.5, fontface = "bold", size = 14) +
      draw_label("Gene", x = 0.01, y = 0.5, angle = 90, vjust = 0.5, fontface = "bold", size = 14)
    
    plots_pdf[[(annot + 1 + 1)]] <- final_plot
    
    pdf_path <- paste0(output_dir, "outputs_pc", pc, "_res", res, ".pdf")
    pdf(pdf_path,
        height = 12, width = 14)
    for (p in plots_pdf) {
      print(p)
    }
    dev.off()
  }
}

# Clustree ####

output_dir <- paste("GBM_single_cell_analysis/outputs/clustering/clustree/")
dir.create(output_dir, recursive=T, showWarnings=F)

dims <- seq(15, 40, by = 5)
fixed_dim <- list()
for (dim in dims){
  p <- clustree(thesis_seu, prefix = paste0("harm_clusters_pc", dim, "_res"))
  fixed_dim[[dim]] <- p
}

pdf_path <- paste0(output_dir, "fixed_dim.pdf")
pdf(pdf_path,
    height = 10, width = 14)
for (p in fixed_dim) {
  print(p)
}
dev.off()

resolutions <- seq(0.2, 1.2, by = 0.2)
fixed_res <- list()
for (res in resolutions){
  for (dim in dims){
    copy_colname <- paste0("harm_clusters_pc", dim, "_res", res)
    new_col <- paste0("harm_clusters_res", res, "_pc", dim)
    thesis_seu@meta.data[[new_col]] <- thesis_seu@meta.data[[copy_colname]]
  }
  
  p <- clustree(thesis_seu, prefix = paste0("harm_clusters_res", res, "_pc"))
  fixed_res[[as.character(res)]] <- p
  
  for (dim in dims){
    new_col <- paste0("pc", dim)
    thesis_seu@meta.data[[new_col]] <- NULL
  }
}

pdf_path <- paste0(output_dir, "fixed_res.pdf")
pdf(pdf_path,
    height = 10, width = 14)
for (p in fixed_res) {
  print(p)
}
dev.off()

thesis_seu@meta.data <- thesis_seu@meta.data[, !grepl("^harm_clusters_res", colnames(thesis_seu@meta.data))]
Idents(thesis_seu) <- thesis_seu$harm_clusters_pc25_res0.6
thesis_seu$seurat_clusters <- thesis_seu$harm_clusters_pc25_res0.6

# Assign primary labels (merged for Lymphoid, Myeloid and Tumour/Brain)

thesis_seu$Primary.labels <- NA
lymphoid_clusters <- c(0, 1, 11, 12, 14)
myeloid_clusters <- c(2, 3, 4, 5, 6, 8, 13)
brain_clusters <- c(7, 9, 10, 15)

thesis_seu$Primary.labels[thesis_seu$harm_clusters_pc25_res0.6 %in% lymphoid_clusters] <- "Lymphoid"
thesis_seu$Primary.labels[thesis_seu$harm_clusters_pc25_res0.6 %in% myeloid_clusters] <- "Myeloid"
thesis_seu$Primary.labels[thesis_seu$harm_clusters_pc25_res0.6 %in% brain_clusters] <- "Tumour/Brain"

thesis_seu$Primary.labels[thesis_seu$Primary.harm_clusters %in% lymphoid_clusters] <- "Lymphoid"
thesis_seu$Primary.labels[thesis_seu$Primary.harm_clusters %in% myeloid_clusters] <- "Myeloid"
thesis_seu$Primary.labels[thesis_seu$Primary.harm_clusters %in% brain_clusters] <- "Tumour/Brain"


# Remove unused clusters/reductions

thesis_seu@meta.data <- thesis_seu@meta.data %>%
  rename("Primary.harm_clusters" = "harm_clusters_pc25_res0.6") %>%
  select(-starts_with("harm_clusters"))

thesis_seu@reductions$Primary.harm_umap <- thesis_seu@reductions$umap.harm_pc25
to_remove <- grep("^umap\\.harm", names(thesis_seu@reductions), value = TRUE)
thesis_seu@reductions[to_remove] <- NULL

thesis_seu@graphs$Primary.harm_snn <- thesis_seu@graphs$harmony_snn_pc25
to_remove <- grep("^harmony_snn", names(thesis_seu@graphs), value = TRUE)
thesis_seu@graphs[to_remove] <- NULL


# Temporary re-integration and clustering by lineage to move outliers ####
thesis_seu <- readRDS("GBM_single_cell_analysis/outputs/thesis_seu.RDS")

load("GBM_single_cell_analysis/BatchProcess.RData")

subseu <- list(
  Lymphoid = subset(thesis_seu, subset = Primary.labels == "Lymphoid"),
  Myeloid = subset(thesis_seu, subset = Primary.labels == "Myeloid")
)

nrow(subseu$Lymphoid@meta.data)
nrow(subseu$Myeloid@meta.data)

seuNames <- names(subseu)

gc()
options(future.globals.maxSize = 7 * 1024^3)

for (obj in seuNames){
  
  cat("Processing", obj, "\n")
  
  seu <- subseu[[obj]]
  seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)
  seu <- BatchProcess(seu)
  
  cat("Re-integrating", "\n")
  
  seu <- IntegrateLayers(
    object = seu,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony_temp",
    verbose = FALSE)
  
  seu <- JoinLayers(seu)
  
  cat("Constructing SNN graph", "\n")
  seu <- FindNeighbors(seu,
                       reduction = "harmony_temp",
                       dims = 1:25,
                       graph.name = "harmony_temp_snn_25")
  
  cat("Generating UMAP", "\n")
  seu <- RunUMAP(seu,
                 reduction = "harmony_temp",
                 dims = 1:25,
                 reduction.name = "umap.harm_temp",
                 graph.name = "harmony_temp_snn_25")
  
  cat("Clustering", "\n")
  seu <- FindClusters(seu,
                      resolution = 1,
                      graph.name = "harmony_temp_snn_25",
                      cluster.name = "harmony_temp_clusters")
  
  output_dir <- paste0("GBM_single_cell_analysis/outputs/clustering/", obj, "/temp/")
  dir.create(output_dir, recursive=T, showWarnings=F)
  
  plots_pdf <- list()
  
  p <- DimPlot(seu,
               reduction = "umap.harm_temp",
               group.by = "harmony_temp_clusters",
               label = TRUE, label.size = 4.5) +
    ggtitle(paste0("Temporary clusters (", obj, ", pc25, res1)")) +
    labs(x = "umap_1", y = "umap_2") +
    guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
  
  plots_pdf[[1]] <- p
  
  ggsave(paste0(output_dir, "harmony_temp_clusters.png"),
         p,
         height = 8.5, width = 10)
  
  cat("Plotting annotations...", "\n")
  
  annotations <- c("SingleR.labels", "CelliD.labels", "ScType.labels", "GBMap_3.labels", "GBMap_4.labels")
  
  for (annot in seq_along(annotations)){
    p <- DimPlot(seu,
                 reduction = "umap.harm_temp",
                 group.by = paste0(annotations[[annot]]),
                 label = TRUE, label.size = 3) +
      ggtitle(paste0(obj, " ", annotations[[annot]], " (pc25_res1")) +
      labs(x = "umap_1", y = "umap_2") +
      guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
    
    plots_pdf[[(annot + 1)]] <- p
    
    ggsave(paste0(output_dir, obj, "_", annotations[[annot]], "_pc25_res1.png"),
           p,
           height = 8.5, width = 11)
    
  }
  
  cat("Running differential gene expression analysis for generated clusters", "\n")
  
  markers <- FindAllMarkers(seu, 
                            logfc.threshold = 0.5, 
                            min.pct = 0.1, 
                            only.pos = TRUE)
  
  clusterMarkers <- markers %>%
    mutate(
      diff_pct = (pct.1 - pct.2) * 100,
      signif = -log10(p_val_adj),
      signif = ifelse(is.infinite(signif) | signif > 300, 300, signif)) %>%
    group_by(cluster) %>%
    arrange(cluster, desc(diff_pct)) %>%
    slice_head(n = 20) %>%
    ungroup()
  
  m.data <- seu@meta.data %>%
    arrange(seurat_clusters)
  
  clusters <- unique(m.data$seurat_clusters)
  cluster_plotlist <- list()
  
  cat("Plotting clusters...", "\n")
  
  for (i in seq_along(clusters)){
    cluster_id <- clusters[[i]]
    cluster_cells <- length(which(seu$seurat_clusters == cluster_id))
    cluster_size <- round(cluster_cells / nrow(seu@meta.data) * 100, 2)
    
    topmarkers <- clusterMarkers %>%
      filter(cluster == cluster_id) %>%
      arrange(desc(diff_pct))
    
    p <- ggplot(topmarkers, aes(x = diff_pct, y = reorder(gene, diff_pct), size = signif)) +
      # Trailing lines from y-axis to dots
      geom_segment(aes(x = 0, xend = diff_pct, y = reorder(gene, diff_pct), yend = reorder(gene, diff_pct)),
                   color = "gray", size = 0.5, alpha = 0.5) +  # Line settings
      
      # Dot plot
      geom_point(aes(color = avg_log2FC, x = diff_pct, y = gene), alpha = 1) +  
      labs(title = paste0(cluster_id, "; n = ", cluster_cells, "; ", cluster_size, "%")) +
      scale_color_gradient(
        low = "#132132", 
        high = "steelblue1", 
        limits = range(clusterMarkers$avg_log2FC, na.rm = TRUE),
        name = "Fold change"
      ) +
      scale_size_continuous(name = "Significance",
                            range = c(0.05, 2),
                            limits = range(clusterMarkers$signif, na.rm = TRUE)) +  # Adjust dot size range
      xlim(0, 100) +
      theme_cowplot() +
      theme(axis.text.y = element_text(size = 4.5),
            axis.text.x = element_text(size = 5.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
            legend.title = element_text(size = 5, face = "bold"),
            legend.text = element_text(size = 5),
            legend.position = "none",
            plot.background = element_rect(fill = "white", color = NA))
    
    cluster_plotlist[[i]] <- p
    
  }
  
  legend_plot <- ggplot(topmarkers, aes(x = diff_pct, y = reorder(gene, diff_pct), size = signif)) +
    geom_segment(aes(x = 0, xend = diff_pct, y = reorder(gene, diff_pct), yend = reorder(gene, diff_pct)),
                 color = "gray", size = 0.5, alpha = 0.5) +
    geom_point(aes(color = avg_log2FC, x = diff_pct, y = gene), alpha = 1) +  
    labs(title = paste0(cluster_id, "; n = ", cluster_cells, "; ", cluster_size, "%")) +
    scale_color_gradient(
      low = "#132132", 
      high = "steelblue1", 
      limits = range(clusterMarkers$avg_log2FC, na.rm = TRUE),
      name = "Fold change"
    ) +
    scale_size_continuous(name = "Significance",
                          range = c(0.2, 3),
                          limits = range(clusterMarkers$signif, na.rm = TRUE)) +  # Adjust dot size range
    xlim(0, 100) +
    theme_cowplot() +
    theme(axis.text = element_text(size = 5), 
          axis.title = element_text(size = 7, face = "bold"),
          plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 10, face = "bold"),
          legend.text = element_text(size = 10),
          legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA))
  
  legend <- get_legend(legend_plot)
  
  cluster_grid <- plot_grid(plotlist = cluster_plotlist, ncol = 5)
  plot <- plot_grid(cluster_grid, legend,
                    rel_widths = c(5,1), nrow = 1, align = "vh")
  final_plot <- ggdraw(plot) +
    draw_label("%(cluster) - %(rest)", x = 0.5, y = 0.01, vjust = 0.5, fontface = "bold", size = 14) +
    draw_label("Gene", x = 0.01, y = 0.5, angle = 90, vjust = 0.5, fontface = "bold", size = 14)
  
  plots_pdf[[(annot + 1 + 1)]] <- final_plot
  
  pdf_path <- paste0(output_dir, "temp_outputs.pdf")
  pdf(pdf_path,
      height = 12, width = 14)
  for (p in plots_pdf) {
    print(p)
  }
  dev.off()
  
  subseu[[obj]] <- seu
  
  rm(cluster_grid, cluster_plotlist, clusterMarkers, final_plot, legend,
     legend_plot, m.data, markers, topmarkers, p, plot, plots_pdf, seu)
  
  cat(obj, "temp clustering complete", "\n")
}

# Exploratory analysis ####

misplaced_macs <- subset(subseu[["Lymphoid"]],
                         subset = harmony_temp_clusters == 13 & !is.na(cloneSize))

thesis_seu$Primary.labels[Cells(thesis_seu) %in% Cells(misplaced_tcrs)] <- "Lymphoid"

misplaced_tils <- subset(subseu[["Myeloid"]],
                         subset = harmony_temp_clusters == 13)


thesis_seu$Primary.labels[Cells(thesis_seu) %in% Cells(misplaced_tcrs)] <- "Lymphoid"


misplaced_tcrs <- subset(thesis_seu,
                         subset = Primary.labels != "Lymphoid" & !is.na(cloneSize))

cd3e_expr <- GetAssayData(misplaced_tcrs, slot = "data")["CD3E", ]
nrow(cd3e_expr)
sum(cd3e_expr > 1)
View(cd3e_expr)

clusterTen$CD3_status <- ifelse(
  cd3e_expr > 0, 
  "Positive",
  ifelse(cd3e_expr == 0, 
         "Negative", 
         NA)
)
rm(misplaced_macs, misplaced_tils, misplaced_tcrs)

# Move back to temp clustering as necessary


# Explore c9
genes <- c("TYROBP", "CCL4", "CD3E", "C1QC", "APOE")

c10 <- subset(subseu$Lymphoid,
              subset = Secondary.harm_clusters == "10")

tilOrMac <- lapply(genes, function(gene){
  df <- FetchData(c10, vars = c("Secondary.harm_clusters", gene)) %>%
    mutate(gene = gene,
           expressed = .data[[gene]] > 0) %>%
    group_by(Secondary.harm_clusters, gene) %>%
    summarise(proportion = mean(expressed)) %>%
    ungroup()
}) %>% bind_rows()

unclearCells <- subset(subseu$Lymphoid,
                       subset = seurat_clusters == "9")

Idents(unclearCells) <- unclearCells$Subclusters_9

degs <- FindAllMarkers(unclearCells, 
                       logfc.threshold = 0.5, 
                       min.pct = 0.1,
                       only.pos = TRUE)

degs_grouped <- degs %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  ungroup()

degs_grouped <- degs_grouped %>%
  mutate(diff_pct = pct.1 - pct.2) %>%
  arrange(diff_pct)

clusterTen <- subset(subseu$Myeloid,
                     subset = harmony_temp_clusters == "10")

immune <- subset(thesis_seu,
                 subset = Primary.labels %in% c("Lymphoid", "Myeloid"))

Idents(immune) <- immune$Primary.labels

degs <- FindAllMarkers(immune,
                       logfc.threshold = 0.5, 
                       min.pct = 0.1, 
                       only.pos = TRUE)

degs_grouped <- degs %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  ungroup()

rm(immune)

clusterTen <- ScaleData(clusterTen, features = degs_grouped$gene, verbose = FALSE, scale = TRUE)
clusterTen <- RunPCA(clusterTen, npcs = 30, features = degs_grouped$gene, verbose = FALSE, maxit = 1000)
elbow <- ElbowPlot(clusterTen, ndims = 30)
clusterTen <- FindNeighbors(clusterTen, dims = 1:7, graph.name = "tilOrMac_snn_pc7")
clusterTen <- FindClusters(clusterTen, resolution = 0.1, graph.name = "tilOrMac_snn_pc7", cluster.name = "tilOrMac_clusters")
clusterTen <- RunUMAP(clusterTen, dims = 1:7, reduction.name = "umap.tilOrMac_pc7", graph.name = "tilOrMac_snn_pc7")

DimPlot(clusterTen, group.by = "tilOrMac_clusters")

clusterTen <- subset(Lseu,
                     subset = Secondary.harm_clusters == "10")

genes <- c("CD3E", "CD4" , "CD8B", "IL7R")

cd4cd8 <- lapply(genes, function(gene){
  df <- FetchData(clusterTen, vars = c("Secondary.harm_clusters", gene)) %>%
    mutate(gene = gene,
           expressed = .data[[gene]] > 0) %>%
    group_by(Secondary.harm_clusters, gene) %>%
    summarise(proportion = mean(expressed)) %>%
    ungroup()
}) %>% bind_rows() %>%
  arrange(Secondary.harm_clusters)

RidgePlot(clusterTen, features = "CD4-hashtag", group.by = "Secondary.harm_clusters")

sum(clusterTen$tilOrMac_clusters == 2)

clusterTen$nk_cd8_clusters <- as.character(clusterTen$nk_cd8_clusters)
clusterTen$nk_cd8_clusters[clusterTen$nk_cd8_clusters == "0"] <- "10_0"
clusterTen$nk_cd8_clusters <- factor(clusterTen$nk_cd8_clusters)

clusterTen$nk_cd8_clusters <- as.character(clusterTen$nk_cd8_clusters)
clusterTen$nk_cd8_clusters[clusterTen$nk_cd8_clusters == "1"] <- "10_1"
clusterTen$nk_cd8_clusters <- factor(clusterTen$nk_cd8_clusters)

# Re-integrate for secondary clustering ####

subseu <- list(
  Lymphoid = subset(thesis_seu, subset = Primary.labels == "Lymphoid"),
  Myeloid = subset(thesis_seu, subset = Primary.labels == "Myeloid")
)

table(subseu$Lymphoid$Primary.labels)
table(subseu$Myeloid$Primary.labels)
table(thesis_seu$Primary.labels[Cells(subseu$Lymphoid)])
table(thesis_seu$Primary.labels[Cells(subseu$Myeloid)])

seuNames <- names(subseu)

for (obj in seuNames){
  
  cat("Processing", obj, "\n")
  
  seu <- subseu[[obj]]
  seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)
  seu <- BatchProcess(seu)
  
  cat("Re-integrating", "\n")
  
  seu <- IntegrateLayers(
    object = seu,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony_subset",
    verbose = FALSE)
  
  elbow <- ElbowPlot(seu, reduction = "pca", ndims = 50) +
    ggtitle(paste0("Percentage variance by principal component (", obj, ")")) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 12,face = "bold")
    )
  
  output_dir <- paste0("GBM_single_cell_analysis/outputs/clustering/", obj, "/")
  dir.create(output_dir, recursive=T, showWarnings=F)
  
  ggsave(paste0(output_dir, obj, "_elbow_plot.png"),
         elbow,
         height = 6, width = 6)
  
  cat("Saving...", "\n")
  
  saveRDS(seu, paste0("GBM_single_cell_analysis/outputs/", obj, ".RDS"))
  
  subseu[[obj]] <- seu
  
  rm(elbow, seu)
  
  cat("Complete", "\n")
}

# Secondary clustering ####

for (obj in seuNames){
  
  cat("Processing", obj, "cells", "\n")
  
  seu <- subseu[[obj]]
  seu <- JoinLayers(seu)
  
  pcs <- seq(10, 30, by = 5)
  resolutions <- seq(0.2, 1.2, by = 0.2)
  
  for (pc in pcs){
    cat("Constructing SNN graph,", pc, "dimensions", "\n")
    seu <- FindNeighbors(seu,
                         reduction = "harmony_subset",
                         dims = 1:pc,
                         graph.name = paste0("harmony_subset_snn_pc", pc))
    
    cat("Generating UMAP,", pc, "dimensions", "\n")
    seu <- RunUMAP(seu,
                   reduction = "harmony_subset",
                   dims = 1:pc,
                   reduction.name = paste0("umap.harm_subset_pc", pc),
                   graph.name = paste0("harmony_subset_snn_pc", pc))
    
    for (res in resolutions){
      cat("Clustering,", pc, "dimensions, resolution", res, "\n")
      
      plots_pdf <- list()
      
      output_dir <- paste0("GBM_single_cell_analysis/outputs/clustering/", obj, "/pc", pc, "/pc", pc, "_res", res, "/")
      dir.create(output_dir, recursive=T, showWarnings=F)
      
     seu <- FindClusters(seu,
                         resolution = res,
                         graph.name = paste0("harmony_subset_snn_pc", pc),
                         cluster.name = paste0("harm_subset_clusters_pc", pc, "_res", res))
     
     p <- DimPlot(seu,
                  reduction = paste0("umap.harm_subset_pc", pc),
                  group.by = paste0("harm_subset_clusters_pc", pc, "_res", res),
                  label = TRUE, label.size = 4.5) +
       ggtitle(paste0(obj, " clusters (pc", pc, "_res", res, ")")) +
       labs(x = "umap_1", y = "umap_2") +
       guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
     
     plots_pdf[[1]] <- p
     
     ggsave(paste0(output_dir, obj, "_clusters_pc", pc, "_res", res, ".png"),
            p,
            height = 8.5, width = 10)
     
     cat("Plotting annotations...", "\n")
     
     annotations <- c("SingleR.labels", "CelliD.labels", "ScType.labels", "GBMap_3.labels", "GBMap_4.labels")
     
     for (annot in seq_along(annotations)){
       p <- DimPlot(seu,
                    reduction = paste0("umap.harm_subset_pc", pc),
                    group.by = paste0(annotations[[annot]]),
                    label = TRUE, label.size = 3) +
         ggtitle(paste0(obj, " ", annotations[[annot]], " (pc", pc, "_res", res, ")")) +
         labs(x = "umap_1", y = "umap_2") +
         guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
       
       plots_pdf[[(annot + 1)]] <- p
       
       ggsave(paste0(output_dir, obj, "_", annotations[[annot]], "_pc", pc, "_res", res, ".png"),
              p,
              height = 8.5, width = 11)
       
     }
     
     cat("Running differential gene expression analysis for generated clusters,", pc, "dimensions, resolution", res, "\n")
     
     markers <- FindAllMarkers(seu, 
                               logfc.threshold = 0.5, 
                               min.pct = 0.1, 
                               only.pos = TRUE)
     
     clusterMarkers <- markers %>%
       mutate(
         diff_pct = (pct.1 - pct.2) * 100,
         signif = -log10(p_val_adj),
         signif = ifelse(is.infinite(signif) | signif > 300, 300, signif)) %>%
       group_by(cluster) %>%
       arrange(cluster, desc(diff_pct)) %>%
       slice_head(n = 20) %>%
       ungroup()
     
     m.data <- seu@meta.data %>%
       arrange(seurat_clusters)
     
     clusters <- unique(m.data$seurat_clusters)
     cluster_plotlist <- list()
     
     cat("Plotting clusters...", "\n")
     
     for (i in seq_along(clusters)){
       cluster_id <- clusters[[i]]
       cluster_cells <- length(which(seu$seurat_clusters == cluster_id))
       cluster_size <- round(cluster_cells / nrow(seu@meta.data) * 100, 2)
       
       topmarkers <- clusterMarkers %>%
         filter(cluster == cluster_id) %>%
         arrange(desc(diff_pct))
       
       p <- ggplot(topmarkers, aes(x = diff_pct, y = reorder(gene, diff_pct), size = signif)) +
         # Trailing lines from y-axis to dots
         geom_segment(aes(x = 0, xend = diff_pct, y = reorder(gene, diff_pct), yend = reorder(gene, diff_pct)),
                      color = "gray", size = 0.5, alpha = 0.5) +  # Line settings
         
         # Dot plot
         geom_point(aes(color = avg_log2FC, x = diff_pct, y = gene), alpha = 1) +  
         labs(title = paste0(cluster_id, "; n = ", cluster_cells, "; ", cluster_size, "%")) +
         scale_color_gradient(
           low = "#132132", 
           high = "steelblue1", 
           limits = range(clusterMarkers$avg_log2FC, na.rm = TRUE),
           name = "Fold change"
         ) +
         scale_size_continuous(name = "Significance",
                               range = c(0.05, 2),
                               limits = range(clusterMarkers$signif, na.rm = TRUE)) +  # Adjust dot size range
         xlim(0, 100) +
         theme_cowplot() +
         theme(axis.text.y = element_text(size = 4.5),
               axis.text.x = element_text(size = 5.5),
               axis.title.x = element_blank(),
               axis.title.y = element_blank(),
               plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
               legend.title = element_text(size = 5, face = "bold"),
               legend.text = element_text(size = 5),
               legend.position = "none",
               plot.background = element_rect(fill = "white", color = NA))
       
       cluster_plotlist[[i]] <- p
       
     }
     
     legend_plot <- ggplot(topmarkers, aes(x = diff_pct, y = reorder(gene, diff_pct), size = signif)) +
       geom_segment(aes(x = 0, xend = diff_pct, y = reorder(gene, diff_pct), yend = reorder(gene, diff_pct)),
                    color = "gray", size = 0.5, alpha = 0.5) +
       geom_point(aes(color = avg_log2FC, x = diff_pct, y = gene), alpha = 1) +  
       labs(title = paste0(cluster_id, "; n = ", cluster_cells, "; ", cluster_size, "%")) +
       scale_color_gradient(
         low = "#132132", 
         high = "steelblue1", 
         limits = range(clusterMarkers$avg_log2FC, na.rm = TRUE),
         name = "Fold change"
       ) +
       scale_size_continuous(name = "Significance",
                             range = c(0.2, 3),
                             limits = range(clusterMarkers$signif, na.rm = TRUE)) +  # Adjust dot size range
       xlim(0, 100) +
       theme_cowplot() +
       theme(axis.text = element_text(size = 5), 
             axis.title = element_text(size = 7, face = "bold"),
             plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
             legend.title = element_text(size = 10, face = "bold"),
             legend.text = element_text(size = 10),
             legend.position = "right",
             plot.background = element_rect(fill = "white", color = NA))
     
     legend <- get_legend(legend_plot)
     
     cluster_grid <- plot_grid(plotlist = cluster_plotlist, ncol = 5)
     plot <- plot_grid(cluster_grid, legend,
                       rel_widths = c(5,1), nrow = 1, align = "vh")
     final_plot <- ggdraw(plot) +
       draw_label("%(cluster) - %(rest)", x = 0.5, y = 0.01, vjust = 0.5, fontface = "bold", size = 14) +
       draw_label("Gene", x = 0.01, y = 0.5, angle = 90, vjust = 0.5, fontface = "bold", size = 14)
     
     plots_pdf[[(annot + 1 + 1)]] <- final_plot
     
     pdf_path <- paste0(output_dir, "outputs_pc", pc, "_res", res, ".pdf")
     pdf(pdf_path,
         height = 12, width = 14)
     for (p in plots_pdf) {
       print(p)
     }
     dev.off()
    }
  }
  
  subseu[[obj]] <- seu
  
  rm(cluster_grid, cluster_plotlist, clusterMarkers, final_plot, legend,
     legend_plot, m.data, markers, topmarkers, p, plot, plots_pdf, seu)
  
  cat(obj, "clustering complete", "\n")
}

# Secondary clustree ####

fixed_dim <- as.list(seuNames)
fixed_res <- as.list(seuNames)

for (obj in seuNames){
  dims <- seq(10, 30, by = 5)
  fixed_dim[[obj]] <- list()
  
  output_dir <- paste0("GBM_single_cell_analysis/outputs/clustering/clustree/", obj, "/")
  dir.create(output_dir, recursive=T, showWarnings=F)
  
  for (dim in dims){
    p <- clustree(subseu[[obj]], prefix = paste0("harm_subset_clusters_pc", dim, "_res")) 
    fixed_dim[[obj]][[dim]] <- p
  }
  
  pdf_path <- paste0(output_dir, obj,"_fixed_dim.pdf")
  pdf(pdf_path,
    height = 10, width = 14)
  for (p in fixed_dim) {
  print(p)
    }
  dev.off()
  
  resolutions <- seq(0.2, 1.2, by = 0.2)
  fixed_res[[obj]] <- list()
  
  for (res in resolutions){
    for (dim in dims){
      copy_colname <- paste0("harm_subset_clusters_pc", dim, "_res", res)
      new_col <- paste0("harm_subset_clusters_res", res, "_pc", dim)
      subseu[[obj]]@meta.data[[new_col]] <-  subseu[[obj]]@meta.data[[copy_colname]]
      }
    
    p <- clustree(subseu[[obj]], prefix = paste0("harm_subset_clusters_res", res, "_pc"))
    fixed_res[[obj]][[as.character(res)]] <- p
    
    for (dim in dims){
    new_col <- paste0("pc", dim)
    subseu[[obj]]@meta.data[[new_col]] <- NULL
    }
    
    }
  
  pdf_path <- paste0(output_dir, obj, "_fixed_res.pdf")
  pdf(pdf_path,
      height = 10, width = 14)
  for (p in fixed_res[[obj]]) {
    print(p)
    }
  
  dev.off()
  
}

rm(fixed_dim, fixed_res, p)

Lseu <- subseu$Lymphoid
Mseu <- subseu$Myeloid
################################### Return to here to manage myeloid object ####
# Lymphoid annotation ####

# Remove unused clusters/reductions ####

Lseu@meta.data <- Lseu@meta.data[, !grepl("^harm_subset_clusters_res", colnames(Lseu@meta.data))]
Idents(Lseu) <- Lseu$harm_subset_clusters_pc25_res0.8
Lseu$seurat_clusters <- Lseu$harm_subset_clusters_pc25_res0.8

Lseu$Secondary.harm_clusters <- Lseu$harm_subset_clusters_pc25_res0.8
Lseu@meta.data <- Lseu@meta.data %>%
  select(-starts_with("harm_subset_clusters"))

Lseu@reductions$Secondary.harm_umap <- Lseu@reductions$umap.harm_subset_pc25
to_remove <- grep("^umap\\.harm_subset", names(Lseu@reductions), value = TRUE)
Lseu@reductions[to_remove] <- NULL

Lseu@graphs$Secondary.harm_snn <- Lseu@graphs$harmony_subset_snn_pc25
to_remove <- grep("^harmony_subset", names(Lseu@graphs), value = TRUE)
Lseu@graphs[to_remove] <- NULL

# Overwrite elbow plot

elbow <- ElbowPlot(Lseu, reduction = "pca", ndims = 50) +
  ggtitle("Percentage variance by principal component (Lymphoid)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12,face = "bold")
  )   +
  geom_vline(xintercept = 25.5, color = "red", linetype = "dashed", linewidth = 0.5)

output_dir <- paste0("GBM_single_cell_analysis/outputs/clustering/Lymphoid/")
dir.create(output_dir, recursive=T, showWarnings=F)

ggsave(paste0(output_dir, "Lymphoid_elbow_plot.png"),
       elbow,
       height = 6, width = 6)

# Check proportions by cluster

genes <- c("CD3E", "CD4", "CD8B")

cd4cd8 <- lapply(genes, function(gene){
  df <- FetchData(Lseu, vars = c("Secondary.harm_clusters", gene)) %>%
    mutate(gene = gene,
           expressed = .data[[gene]] > 0) %>%
    group_by(Secondary.harm_clusters, gene) %>%
    summarise(proportion = mean(expressed)) %>%
    ungroup()
}) %>% bind_rows() %>%
  arrange(Secondary.harm_clusters)

View(cd4cd8)



# Lymphoid secondary UMAP and greying out clusters ####
Lseu$Secondary.harm_clusters <- factor(
  Lseu$Secondary.harm_clusters,
  levels = sort(as.numeric(as.character(unique(Lseu$Secondary.harm_clusters))))
)

p <- DimPlot(Lseu,
             reduction = "Secondary.harm_umap",
             group.by = "Secondary.harm_clusters",
             label = TRUE, label.size = 4) +
  ggtitle("Clustering of lymphoid cells") +
  labs(x = "umap_1", y = "umap_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/lymphoid_clusters.png",
       p,
       height = 6.2, width = 7.75)



col_umap_data <- Lseu

orig_colours <- scales::hue_pal()(length(unique(col_umap_data$Secondary.harm_clusters)))
names(orig_colours) <- levels(col_umap_data$Secondary.harm_clusters)

col_umap_data$Secondary.harm_clusters <- as.character(col_umap_data$Secondary.harm_clusters)
col_umap_data$Secondary.harm_clusters[col_umap_data$Secondary.final_clusters == "12_1"] <- "12_1"
col_umap_data$Secondary.harm_clusters <- factor(col_umap_data$Secondary.harm_clusters)
orig_colours["12_1"] <- orig_colours["12"]



keep_coloured <- c("0", "1", "2", "3", "4", "5", "6", "7", "8", "9", "11", "12_1", "14")
make_grey <- setdiff(levels(col_umap_data$Secondary.harm_clusters), keep_coloured)
orig_colours[make_grey] <- "grey80"

col_umap <- DimPlot(col_umap_data,
                    reduction = "Secondary.harm_umap",
                    group.by = "Secondary.harm_clusters",
                    cols = orig_colours,
                    label = F) +
  labs(x = "umap_1", y = "umap_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/lymphoid_clusters_greyandcolour.png",
       col_umap,
       height = 6.2, width = 7.75)

rm(col_umap, col_umap_data)

# DGEA by cluster for final assignments to present marker genes ####

Lseu$Secondary.final_clusters <- Lseu$Secondary.harm_clusters
Lseu$Secondary.harm_clusters <- as.character(Lseu$Secondary.harm_clusters)
Lseu$Secondary.harm_clusters[Lseu$Secondary.harm_clusters %in% c("12_0", "12_1")] <- "12"
Lseu$Secondary.harm_clusters <- factor(Lseu$Secondary.harm_clusters)

Idents(Lseu) <- Lseu$Secondary.harm_clusters
markers <- FindAllMarkers(Lseu, 
                          logfc.threshold = 0.5, 
                          min.pct = 0.1, 
                          only.pos = TRUE)

clusterMarkers <- markers %>%
  mutate(
    diff_pct = (pct.1 - pct.2) * 100,
    signif = -log10(p_val_adj),
    signif = ifelse(is.infinite(signif) | signif > 300, 300, signif)) %>%
  group_by(cluster) %>%
  arrange(cluster, desc(diff_pct)) %>%
  slice_head(n = 20) %>%
  ungroup()

m.data <- Lseu@meta.data %>%
  arrange(Secondary.harm_clusters)

clusters <- unique(m.data$Secondary.harm_clusters)
cluster_plotlist <- list()

cat("Plotting clusters...", "\n")

for (i in seq_along(clusters)){
  cluster_id <- clusters[[i]]
  cluster_cells <- length(which(Lseu$Secondary.harm_clusters == cluster_id))
  cluster_size <- round(cluster_cells / nrow(Lseu@meta.data) * 100, 2)
  
  topmarkers <- clusterMarkers %>%
    filter(cluster == cluster_id) %>%
    arrange(desc(diff_pct))
  
  p <- ggplot(topmarkers, aes(x = diff_pct, y = reorder(gene, diff_pct), size = signif)) +
    # Trailing lines from y-axis to dots
    geom_segment(aes(x = 0, xend = diff_pct, y = reorder(gene, diff_pct), yend = reorder(gene, diff_pct)),
                 color = "gray", size = 0.5, alpha = 0.5) +  # Line settings
    
    # Dot plot
    geom_point(aes(color = avg_log2FC, x = diff_pct, y = gene), alpha = 1) +  
    labs(title = paste0(cluster_id, "; n = ", cluster_cells, "; ", cluster_size, "%")) +
    scale_color_gradient(
      low = "#132132", 
      high = "steelblue1", 
      limits = range(clusterMarkers$avg_log2FC, na.rm = TRUE),
      name = "Fold change"
    ) +
    scale_size_continuous(name = "Significance",
                          range = c(0.05, 2),
                          limits = range(clusterMarkers$signif, na.rm = TRUE)) +  # Adjust dot size range
    xlim(0, 100) +
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 4.5),
          axis.text.x = element_text(size = 5.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 5, face = "bold"),
          legend.text = element_text(size = 5),
          legend.position = "none",
          plot.background = element_rect(fill = "white", color = NA))
  
  cluster_plotlist[[i]] <- p
  
}

legend_plot <- ggplot(topmarkers, aes(x = diff_pct, y = reorder(gene, diff_pct), size = signif)) +
  geom_segment(aes(x = 0, xend = diff_pct, y = reorder(gene, diff_pct), yend = reorder(gene, diff_pct)),
               color = "gray", size = 0.5, alpha = 0.5) +
  geom_point(aes(color = avg_log2FC, x = diff_pct, y = gene), alpha = 1) +  
  labs(title = paste0(cluster_id, "; n = ", cluster_cells, "; ", cluster_size, "%")) +
  scale_color_gradient(
    low = "#132132", 
    high = "steelblue1", 
    limits = range(clusterMarkers$avg_log2FC, na.rm = TRUE),
    name = "Fold change"
  ) +
  scale_size_continuous(name = "Significance",
                        range = c(0.2, 3),
                        limits = range(clusterMarkers$signif, na.rm = TRUE)) +  # Adjust dot size range
  xlim(0, 100) +
  theme_cowplot() +
  theme(axis.text = element_text(size = 5), 
        axis.title = element_text(size = 7, face = "bold"),
        plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "right",
        plot.background = element_rect(fill = "white", color = NA))

legend <- get_legend(legend_plot)

unique(Lseu$Secondary.harm_clusters) # to get positions in plotlist - note, factors

initial_plots <- list(cluster_plotlist[[7]],
                      cluster_plotlist[[13]],
                      cluster_plotlist[[8]],
                      cluster_plotlist[[14]],
                      cluster_plotlist[[11]]
)

cluster_grid <- plot_grid(plotlist = initial_plots, ncol = 2)

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/initial_assignment_plots.png",
       cluster_grid,
       height = 7, width = 6)

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/initial_plots_legend.png",
       legend,
       height = 3, width = 1.5)

# NK subclustering using Ruiz-Moreno DEGs ####

cluster12_RM <- subset(Lseu,
                       subset = Secondary.harm_clusters == "12")

featurelist <- read.xlsx("GBM_single_cell_analysis/outputs/clustering/NK_vs_CD8_genes.xlsx", colNames = T)
cluster12_RM <- ScaleData(cluster12_RM, features = featurelist$Gene, verbose = FALSE, scale = TRUE)
cluster12_RM <- RunPCA(cluster12_RM, npcs = 25, features = featurelist$Gene, verbose = FALSE, maxit = 1000)
elbow <- ElbowPlot(cluster12_RM, ndims = 25)
cluster12_RM <- FindNeighbors(cluster12_RM, dims = 1:10, graph.name = "c12_snn")
cluster12_RM <- FindClusters(cluster12_RM, resolution = 0.3, graph.name = "c12_snn", cluster.name = "c12_clusters")
cluster12_RM <- RunUMAP(cluster12_RM, dims = 1:10, reduction.name = "umap.c12", graph.name = "c12_snn")

DimPlot(cluster12_RM, reduction = "umap.c12", group.by = "c12_clusters", label = T)
genes <- c("CD3E", "CD8B", "IL7R", "FCER1G", "KLRF1", "KIR2DL1")

c12_exp <- lapply(genes, function(gene){
  df <- FetchData(cluster12_RM, vars = c("c12_clusters", gene)) %>%
    mutate(gene = gene,
           expressed = .data[[gene]] > 0) %>%
    group_by(c12_clusters, gene) %>%
    summarise(proportion = mean(expressed)) %>%
    ungroup()
}) %>% bind_rows() %>%
  arrange(c12_clusters)

View(c12_exp)

FeaturePlot(cluster12_RM, reduction = "umap.c12", features = genes, ncol = 3)

cluster12_RM$c12_clusters <- as.character(cluster12_RM$c12_clusters)
cluster12_RM$c12_clusters[cluster12_RM$c12_clusters == "0"] <- "12_0"
cluster12_RM$c12_clusters <- factor(cluster12_RM$c12_clusters)

cluster12_RM$c12_clusters <- as.character(cluster12_RM$c12_clusters)
cluster12_RM$c12_clusters[cluster12_RM$c12_clusters == "1"] <- "12_1"
cluster12_RM$c12_clusters <- factor(cluster12_RM$c12_clusters)

# Plot subclusters:

elbow_edited <- elbow +
  ggtitle("Percentage variance by principal component") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12,face = "bold")
  ) +
  geom_vline(xintercept = 10.5, color = "red", linetype = "dashed", linewidth = 0.5)

output_dir <- paste0("GBM_single_cell_analysis/outputs/clustering/Lymphoid/NK_subclustering/")
dir.create(output_dir, showWarnings = F, recursive = T)

ggsave(paste0(output_dir, "c12_elbow.png"),
       elbow_edited,
       height = 5.5, width = 6)

nk_umap <- DimPlot(cluster12_RM, reduction = "umap.c12", group.by = "c12_clusters", label = F, label.size = 5) +
  labs(x = "umap_1", y = "umap_2") +
  ggtitle("NK cell/NK-like T cell subclustering") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 10))

cluster12_RM$clonalFrequency_scaled <- log1p(cluster12_RM$clonalFrequency)
p2 <- FeaturePlot(cluster12_RM, reduction = "umap.c12", features = "clonalFrequency_scaled", pt.size = 0.7) +
  labs(x = "umap_1", y = "umap_2") +
  ggtitle("TCR clonal frequency") +
  theme(plot.title = element_text(size = 12, face = "bold"),
        axis.title = element_text(size = 10))

ggsave(paste0(output_dir, "c12_cloneSize.png"),
       p2,
       height = 2.2, width = 3)

genes <- c("FCER1G", "KLRF1", "KIR2DL1", "CD3E", "CD8B", "IL7R")

fp_list <- FeaturePlot(cluster12_RM, reduction = "umap.c12", features = genes)

fp_labeled <- lapply(fp_list, function(p) {
  p + labs(x = "umap_1", y = "umap_2")
})

nk_featureplots <- wrap_plots(fp_labeled, ncol = 3)

ggsave(paste0(output_dir, "nk_featureplots.png"),
       nk_featureplots,
       height = 5, width = 10)

# Join subclusters back to main lymphoid object

c12_clusters_df <- data.frame(c12_clusters = cluster12_RM$c12_clusters, cell_id = names(cluster12_RM$c12_clusters))
lymphoid_meta <- Lseu@meta.data %>%
  tibble::rownames_to_column(var = "cell_id")

lymphoid_meta <- left_join(lymphoid_meta, c12_clusters_df, by = "cell_id")
rownames(lymphoid_meta) <- lymphoid_meta$cell_id
lymphoid_meta$cell_id <- NULL
Lseu@meta.data <- lymphoid_meta

Lseu$c12_clusters <- as.character(Lseu$c12_clusters)
idx <- Lseu$Secondary.harm_clusters == "12"
Lseu$Secondary.final_clusters[idx] <- Lseu$c12_clusters[idx]
Lseu$Secondary.final_clusters <- as.factor(Lseu$Secondary.harm_clusters)
Lseu$c12_clusters <- NULL
Lseu$seurat_clusters <- Lseu$Secondary.final_clusters

rm(c12_clusters_df, c12_exp, cluster12_RM, elbow, elbow_edited, featurelist,
   fp_labeled, fp_list, lymphoid_meta, nk_featureplots, nk_umap, idx)


# Plot broad automatic assignments of remaining cells ####

tils <- Lseu

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R"); 
tils <- run_sctype(tils, assay = "RNA", scaled = TRUE, known_tissue_type="Immune system", custom_marker_file="https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/ScTypeDB_short.xlsx",name="ScType.labels")

tils <- run_sctype(tils, assay = "RNA", scaled = TRUE, known_tissue_type = "Glioblastoma", custom_marker_file = "GBM_single_cell_analysis/outputs/clustering/GBMap_4DB.xlsx", name = "GBMap_4.labels")

tils$SingleR.simple <- NA
tils$ScType.simple <- NA
tils$GBMap.simple <- NA

SingleR.cd4 <- c("CD4+ Tcm", "Tregs", "CD4+ Tem", "CD4+ T-cells")
SingleR.cd8 <- c("CD8+ Tcm", "CD8+ T-cells", "CD8+ Tem")
ScType.cd4 <- c("Naive CD4+ T cells", "Memory CD4+ T cells")
ScType.cd8 <- "CD8+ NKT-like cells"
GBMap.cd4 <- c("CD4 rest", "CD4 INF", "Reg T")
GBMap.cd8 <- c("CD8 cytotoxic",  "CD8 NK sig")

tils$Secondary.harm_clusters <- as.character(tils$Secondary.harm_clusters)
tils$Secondary.harm_clusters[tils$Secondary.final_clusters == "12_0"] <- "12_0"
non_tils <- c("10", "12_0", "13")
tils$SingleR.simple[tils$SingleR.labels %in% SingleR.cd4] <- "CD4+ T cells"
tils$SingleR.simple[tils$SingleR.labels %in% SingleR.cd8] <- "CD8+ T cells"
tils$SingleR.simple[is.na(tils$SingleR.simple)] <- "Other/indeterminate"
tils$SingleR.simple[tils$Secondary.final_clusters %in% non_tils] <- "Excluded"

tils$ScType.simple[tils$ScType.labels %in% ScType.cd4] <- "CD4+ T cells"
tils$ScType.simple[tils$ScType.labels %in% ScType.cd8] <- "CD8+ T cells"
tils$ScType.simple[is.na(tils$ScType.simple)] <- "Other/indeterminate"
tils$ScType.simple[tils$Secondary.final_clusters %in% non_tils] <- "Excluded"

tils$GBMap.simple[tils$GBMap_4.labels %in% GBMap.cd8] <- "CD8+ T cells"
tils$GBMap.simple[tils$GBMap_4.labels %in% GBMap.cd4] <- "CD4+ T cells"
tils$GBMap.simple[is.na(tils$GBMap.simple)] <- "Other/indeterminate"
tils$GBMap.simple[tils$Secondary.final_clusters %in% non_tils] <- "Excluded"

# Grey out umap based on any column (example for ScTypeDB plotting) ####

tils$ScType.simple <- factor(tils$ScType.simple,
                              levels = c("CD8+ T cells",
                                         "CD4+ T cells",
                                         "Other/indeterminate",
                                         "Excluded")
)

group_levels <- levels(tils$ScType.simple)
orig_colours <- setNames(scales::hue_pal()(length(group_levels)), group_levels)

orig_colours <- c("#F8766D", "#00C0AF", "#B983FF", "grey80") # to match later figure

keep_coloured <- setdiff(levels(tils$ScType.simple), "Excluded")
orig_colours[!names(orig_colours) %in% keep_coloured] <- "grey80"

col_tils <- DimPlot(tils,
                    reduction = "Secondary.harm_umap",
                    group.by = "ScType.simple",
                    cols = orig_colours,
                    label = F) +
  ggtitle("ScTypeDB") +
  labs(x = "umap_1", y = "umap_2") +
  theme(plot.title = element_text(size = 24, face = "bold"),
        axis.title = element_text(size = 20)) +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/ScType_greyandcolour.png",
       col_tils,
       height = 6, width = 8.5)

# Check CD4/8 expression in clusters ####

genes <- c("CD3E", "CD4", "CD8B")

byCluster <- lapply(genes, function(gene){
  df <- FetchData(Lseu, vars = c("Secondary.harm_clusters", gene)) %>%
    mutate(gene = gene,
           expressed = .data[[gene]] > 0) %>%
    group_by(Secondary.harm_clusters, gene) %>%
    summarise(proportion = mean(expressed)) %>%
    ungroup()
}) %>% bind_rows() %>%
  arrange(Secondary.harm_clusters)
View(byCluster)

# Get DEGs for CD4/CD8 subclustering ####

cd3e_expr <- GetAssayData(Lseu, slot = "data")["CD3E", ]
cd4_expr <- GetAssayData(Lseu, slot = "data")["CD4", ]
cd8a_expr <- GetAssayData(Lseu, slot = "data")["CD8A", ]
cd8b_expr <- GetAssayData(Lseu, slot = "data")["CD8B", ]

Lseu$cd4cd8_status <- ifelse(
  cd3e_expr > 0 & cd4_expr > 0 & cd8b_expr == 0,
  "CD4",
  ifelse(
    cd3e_expr > 0 & cd4_expr == 0 & cd8b_expr > 0,
    "CD8",
    NA)
)

Idents(Lseu) <- "cd4cd8_status"

deg_cells <- subset(Lseu,
                    idents = c("CD4", "CD8"))

degs <- FindAllMarkers(deg_cells, 
                       logfc.threshold = 0.5, 
                       min.pct = 0.1,
                       only.pos = TRUE)

degs_grouped <- degs %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(p_val_adj) %>%
  ungroup()

rm(deg_cells)

degs_grouped$avg_log2FC[degs_grouped$cluster == "CD4"] <- degs_grouped$avg_log2FC[degs_grouped$cluster == "CD4"] * -1

cd4cd8_volcano <- EnhancedVolcano(degs_grouped,
                                  lab = degs_grouped$gene,
                                  x = 'avg_log2FC',
                                  y = 'p_val_adj',
                                  xlab = "Log2FC",
                                  ylab = "-Log10 Adjusted p-value",
                                  pCutoff = 0.05,
                                  FCcutoff = 0.5,
                                  pointSize = 2,
                                  labSize = 3.5,
                                  max.overlaps = 35,
                                  title = "          CD4+ CD8B- ←→ CD4- CD8B+",
                                  subtitle = NULL,
                                  legendPosition = "none",
                                  drawConnectors = TRUE,
                                  colConnectors = "gray"
)

output_dir <- paste0("GBM_single_cell_analysis/outputs/clustering/Lymphoid/cd4cd8_subclustering/")
dir.create(output_dir, showWarnings=F, recursive=T)
ggsave(paste0(output_dir, "cd4cd8_volcano.png"),
       cd4cd8_volcano,
       height = 6, width = 6)

# Check CD4 expression in CD3+CD8A-CD8B- ####

Lseu$cd8_status <- ifelse(
  cd3e_expr > 0 & cd8a_expr == 0 & cd8b_expr == 0,
  "Negative",
  "Other"
) 

cd8_neg <- subset(Lseu,
                  subset = cd8_status == "Negative")


Idents(cd8_neg) <- cd8_neg$cd8_status
VlnPlot(cd8_neg, features = "CD4")

pct.cd4_in_cd8negs <- (sum(cd4_expr > 0) / nrow(cd8_neg@meta.data)) * 100

# Get CD4-hashtag info ####

Lseu$CD4_hashtag <- FetchData(Lseu, vars = "CD4-hashtag", assay = "CITE")[,1]

pctHT_in_cd8s <- (sum(Lseu$cd4cd8_status == "CD8" & Lseu$CD4_hashtag > 0.8, na.rm=T) / sum(Lseu$cd4cd8_status == "CD8", na.rm=T)) * 100

meanHT_in_cd8s <- mean(Lseu$CD4_hashtag[Lseu$cd4cd8_status == "CD8" & Lseu$CD4_hashtag < 0.8], na.rm = T)

cells <- subset(Lseu,
                subset = !is.na(cd4cd8_status))
cells$cd4cd8_status[cells$cd4cd8_status == "CD4"] <- "CD4+ CD8B-"
cells$cd4cd8_status[cells$cd4cd8_status == "CD8"] <- "CD4- CD8B+"

# Plot HT expression in CD4+ CD8B- and CD4- CD8B+ cells
cd4HT_ridge <- RidgePlot(
  cells,
  features = "CD4-hashtag",
  sort = TRUE,
  ncol = 2,
  group.by = "cd4cd8_status",
  cols = c("#00C0AF", "#F8766D")
) +
  ggtitle("CD4-hashtag expression") +
  theme(
    plot.title = element_text(size = 14,
                              hjust = 0.5,
                              face = "bold"),
    axis.title.x = element_text(size = 12,
                                hjust = 0.5,
                                face = "bold"),
    axis.title.y = element_text(size = 12,
                                hjust = 0.5,
                                face = "bold")) +
  NoLegend() +
  geom_vline(xintercept = 0.8, color = "red", linetype = "dashed", linewidth = 0.75)

ggsave(paste0(output_dir, "cd4HT_cutoff.png"),
       cd4HT_ridge,
       height = 3, width = 10)

# Plot CD4 expression HT+ cells
ht_cells <- subset(Lseu,
                   subset = CD4_hashtag > 0.8)

cd4_expr_ht <- GetAssayData(ht_cells, slot = "data")["CD4",]

pct.cd4_in_HTs <- (sum(cd4_expr_ht > 0) / nrow(ht_cells@meta.data)) * 100

# Subcluster by CD4 and CD8 #### 

# Note - interactive loop, prompting required throughout

Lseu$cd4cd8_clusters <- NA
pc_res_used <- list()

mixed_clusters <- c("1", "3", "5", "6", "7", "8", "9", "11", "14") # note, unable to split or scale cluster 14 due to ?small size - ran new pca on previous scaling

for (cluster in mixed_clusters){
  cat("Processing cluster", cluster, "\n")
  
  cells <- subset(Lseu,
                  subset = Secondary.harm_clusters == cluster)
  
  cat("Re-integrating", "\n")
  
  cells[["RNA"]] <- split(cells[["RNA"]], f = cells$orig.ident)
  
  cells <- ScaleData(cells, features = degs_grouped$gene, verbose = FALSE, scale = TRUE)
  
  cells <- RunPCA(cells, npcs = 30, features = degs_grouped$gene, verbose = FALSE, maxit = 1000)
  
  cells <- IntegrateLayers(
    object = cells,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony_cd4cd8",
    verbose = FALSE)
  
  repeat { # To change number of pcs, accept a resolution below and repeat when prompted
    
    # Inspect Elbow plot and define number of PCs to take forward
    print(ElbowPlot(cells, ndims = 30) +
            ggtitle(paste("Variance by PC - cluster", cluster)))
    
    pcs <- as.integer(readline(prompt = paste0("Enter number of PCs to use for cluster ", cluster, ": ")))
    
    elbow <- ElbowPlot(cells, ndims = 30) +
      ggtitle(paste("Variance by PC - cluster", cluster)) +
      theme_classic() +
      theme(
        plot.title = element_text(size = 18, hjust = 0.5, face = "bold"),
        axis.title = element_text(size = 12,face = "bold")
      ) +
      geom_vline(xintercept = (pcs + 0.5), color = "red", linetype = "dashed", linewidth = 0.5)
    
    cat("Saving elbow plot", "\n")
    output_dir <- paste0("GBM_single_cell_analysis/outputs/clustering/Lymphoid/cd4cd8_subclustering/cluster_", cluster, "/")
    dir.create(output_dir, recursive=T, showWarnings=F)
    
    ggsave(paste0(output_dir, "c", cluster, "_elbow.png"),
           elbow,
           height = 5.5, width = 6)
    
    cells <- JoinLayers(cells)
    cells <- FindNeighbors(cells, dims = 1:pcs, reduction = "harmony_cd4cd8", graph.name = paste0("cd4cd8_snn"))
    
    #Add CD4-hashtag info
    cells$CD4_hashtag <- FetchData(cells, vars = "CD4-hashtag", assay = "CITE")[,1]
    
    repeat { # Loop to adjust resolution
      res <- as.numeric(readline(prompt = "Enter clustering resolution: "))
      
      cat("Clustering...", "\n")
      cells <- FindClusters(cells, resolution = res, graph.name = "cd4cd8_snn", cluster.name = "cd4cd8_clusters")
      
      #Rename clusters
      cells$cd4cd8_clusters <- as.character(cells$cd4cd8_clusters)
      cells$cd4cd8_clusters <- paste0(cluster, "_", cells$cd4cd8_clusters)
      cells$cd4cd8_clusters <- factor(cells$cd4cd8_clusters)
      
      cells <- RunUMAP(cells, dims = 1:pcs, reduction = "harmony_cd4cd8", reduction.name = "umap.cd4cd8", graph.name = "cd4cd8_snn")
      
      #Inspect UMAP 
      genes <- c("CD4", "IL7R", "CD40LG", "GPR183", "CD8B", "CD8A", "GZMH", "NKG7")
      
      print(DimPlot(cells, reduction = "umap.cd4cd8", group.by = "cd4cd8_clusters", label = T, label.size = 5) +
              labs(x = "umap_1", y = "umap_2") +
              ggtitle(paste0("Cluster ", cluster, ", res = ", res)))
      
      readline(prompt = "To move on, press Enter ")
      
      #Inspect key genes, CD4-hashtag and CD4/CD8 expression in subclusters
      fp_list <- FeaturePlot(cells, reduction = "umap.cd4cd8", features = genes)
      
      fp_labeled <- lapply(fp_list, function(p) {
        p + labs(x = "umap_1", y = "umap_2")
      })
      
      featureplots <- wrap_plots(fp_labeled, ncol = 3)
      
      ridgeplot <- RidgePlot(
        cells,
        features = "CD4-hashtag",
        sort = TRUE,
        ncol = 2,
        group.by = "cd4cd8_clusters"
      ) +
        ggtitle(paste0("Cluster ", cluster, " - CD4 hashtag by subcluster")) +
        theme(
          plot.title = element_text(size = 18,
                                    hjust = 0.5,
                                    face = "bold"),
          axis.title.x = element_text(hjust = 0.5,
                                      face = "bold"),
          axis.title.y = element_text(hjust = 0.5,
                                      face = "bold")) +
        NoLegend()
      
      print(ridgeplot + featureplots & NoLegend())
      
      genes <- c("CD3E", "CD4", "CD8B")
      
      cd4cd8 <- lapply(genes, function(gene){
        df <- FetchData(cells, vars = c("cd4cd8_clusters", gene)) %>%
          mutate(gene = gene,
                 expressed = .data[[gene]] > 0) %>%
          group_by(cd4cd8_clusters, gene) %>%
          summarise(proportion = mean(expressed)) %>%
          ungroup()
      }) %>% bind_rows() %>%
        arrange(cd4cd8_clusters)
      
      cd4cd8_ratio <- cd4cd8 %>%
        pivot_wider(names_from = gene, values_from = proportion) %>%
        rowwise() %>%
        mutate(
          cd4_over_cd8 = CD4 / CD8B,
          cd8_over_cd4 = CD8B / CD4,
          num_cells = sum(cells$cd4cd8_clusters == cd4cd8_clusters),
          diff_pctHT_cells_CD8s = ((sum(cells$cd4cd8_clusters == cd4cd8_clusters & cells$CD4_hashtag > 0.8) / num_cells) * 100) - pctHT_in_cd8s,
          diff_meanHT_cells_CD8s = mean(cells$CD4_hashtag[cells$cd4cd8_clusters == cd4cd8_clusters]) - meanHT_in_cd8s
        ) %>%
        ungroup()
      
      View(cd4cd8_ratio)
      
      write.xlsx(cd4cd8_ratio, paste0(output_dir, "c", cluster, "_subclusters_cd4cd8.xlsx"))
      
      decision <- readline(prompt = "Type 'r' to repeat with new resolution, or press Enter to continue: ")
      
      if (tolower(decision) != "r") {
        break  # exit the resolution repeat loop and move on
      }
    }
    
    pc_res_used[[cluster]] <- paste0("pc", pcs, "_res", res)
    
    decision <- readline(prompt = "Type 'p' to repeat with new pcs, or press Enter to continue: ")
    
    if (tolower(decision) != "p") {
      break  # exit the pcs repeat loop and move on
    }
  }
  
  cat("Plotting and returning subclusters", "\n")
  
  #Plot subclusters:
  umap <- DimPlot(cells, reduction = "umap.cd4cd8", group.by = "cd4cd8_clusters", label = T, label.size = 5) +
    labs(x = "umap_1", y = "umap_2") +
    ggtitle(paste0("Cluster ", cluster, ", res = ", res))
  
  ggsave(paste0(output_dir, "c", cluster, "_subclusters_umap.png"),
         umap,
         height = 5, width = 6.7)
  
  #Plot key genes:
  genes <- c("CD4", "IL7R", "CD40LG", "GPR183", "CD8B", "CD8A", "GZMH", "NKG7")
  
  fp_list <- FeaturePlot(cells, reduction = "umap.cd4cd8", features = genes)
  
  fp_labeled <- lapply(fp_list, function(p) {
    p + labs(x = "umap_1", y = "umap_2")
  })
  
  featureplots <- wrap_plots(fp_labeled, ncol = 4)
  
  ggsave(paste0(output_dir, "c", cluster, "_featureplots.png"),
         featureplots,
         height = 5, width = 10)
  
  #Plot CD4-hashtag
  ridgeplot <- RidgePlot(
    cells,
    features = "CD4-hashtag",
    sort = TRUE,
    ncol = 2,
    group.by = "cd4cd8_clusters"
  ) +
    ggtitle(paste0("Cluster ", cluster, " - CD4 hashtag by subcluster")) +
    theme(
      plot.title = element_text(size = 18,
                                hjust = 0.5,
                                face = "bold"),
      axis.title.x = element_text(hjust = 0.5,
                                  face = "bold"),
      axis.title.y = element_text(hjust = 0.5,
                                  face = "bold")) +
    NoLegend()
  
  ggsave(paste0(output_dir, "c", cluster, "_cd4_hashtag.png"),
         ridgeplot,
         height = 5, width = 10)
  
  # Join subclusters back to main lymphoid object
  cd4cd8_clusters_vec <- cells$cd4cd8_clusters
  names(cd4cd8_clusters_vec) <- names(cells$cd4cd8_clusters)
  Lseu$cd4cd8_clusters[names(cd4cd8_clusters_vec)] <- cd4cd8_clusters_vec
  
}

View(Lseu@meta.data)

Lseu$temp <- as.character(Lseu$Secondary.final_clusters)
Lseu$Secondary.final_clusters <- as.character(Lseu$Secondary.harm_clusters)
inds <- Lseu$temp %in% c("12_0", "12_1")
Lseu$Secondary.final_clusters[inds] <- Lseu$temp[inds]
Lseu$cd4cd8_clusters <- Lseu$cd4cd8_clusters - 1 # to re-assign first subclusters as _0 and so on

idx <- !is.na(Lseu$cd4cd8_clusters)
Lseu$Secondary.final_clusters[idx] <- paste0(Lseu$Secondary.final_clusters[idx], "_", Lseu$cd4cd8_clusters[idx])
Lseu$temp <- NULL

rm(cells, elbow, featureplots, fp_list, fp_labeled, ridgeplot, umap, lymphoid_meta,
   cd4cd8_clusters_df, cd4cd8_ratio, cd4cd8_volcano, cd4cd8, degs_grouped, idx)

# Gather and plot data on subclusters ####

genes <- c("CD3E", "CD4", "CD8B")

cd4cd8 <- lapply(genes, function(gene){
  df <- FetchData(Lseu, vars = c("Secondary.final_clusters", gene)) %>%
    mutate(gene = gene,
           expressed = .data[[gene]] > 0) %>%
    group_by(Secondary.final_clusters, gene) %>%
    summarise(proportion = mean(expressed)) %>%
    ungroup()
}) %>% bind_rows() %>%
  arrange(Secondary.final_clusters)

subclusterData <- cd4cd8 %>%
  pivot_wider(names_from = gene, values_from = proportion) %>%
  rowwise() %>%
  mutate(
    cd4_over_cd8 = CD4 / CD8B,
    cd8_over_cd4 = CD8B / CD4,
    num_cells = sum(Lseu$Secondary.final_clusters == Secondary.final_clusters),
    diff_pctHT_cells_CD8s = ((sum(Lseu$Secondary.final_clusters == Secondary.final_clusters & Lseu$CD4_hashtag > 0.8) / num_cells) * 100) - pctHT_in_cd8s,
    diff_meanHT_cells_CD8s = mean(Lseu$CD4_hashtag[Lseu$Secondary.final_clusters == Secondary.final_clusters]) - meanHT_in_cd8s
  ) %>%
  ungroup()

exclusions <- c("10", "12_0", "13")
plotData <- subclusterData %>%
  filter(!Secondary.final_clusters %in% exclusions) %>%
  mutate(
    pct.cd4 = CD4 * 100,
    pct.cd8 = CD8B * 100,
    cd4_ratio = cd4_over_cd8 * -1,
    cd8_ratio = cd8_over_cd4,
    final_ratio = ifelse(abs(cd4_ratio) > abs(cd8_ratio), cd4_ratio, cd8_ratio),
    proportion = pmax(CD4, CD8B) * 100,
  ) %>%
  select(Secondary.final_clusters, final_ratio, proportion, num_cells, diff_pctHT_cells_CD8s, diff_meanHT_cells_CD8s)

assignments <- c("Assigned on fold change", "Assigned on fold change", "Assigned on CD4-hashtag", "Assigned on fold change", "Assigned on fold change", "Assigned on CD4-hashtag",
                 "Assigned on CD4-hashtag", "Assigned on CD4-hashtag", "Assigned on CD4-hashtag", "Assigned on fold change", "Assigned on fold change", "Assigned on CD4-hashtag",
                 "Assigned on CD4-hashtag", "Indeterminate", "Indeterminate", "Assigned on fold change", "Assigned on CD4-hashtag", "Assigned on CD4-hashtag",
                 "Assigned on fold change", "Assigned on CD4-hashtag", "Indeterminate", "Assigned on fold change", "Indeterminate", "Assigned on CD4-hashtag",
                 "Assigned on fold change", "Assigned on CD4-hashtag", "Assigned on fold change", "Assigned on CD4-hashtag")

plotData$assignments <- assignments
plotData$assignments <- factor(plotData$assignments, levels = c("Assigned on fold change", "Assigned on CD4-hashtag", "Indeterminate"))

subclusterPlot <- ggplot(plotData, aes(x = final_ratio, y = diff_pctHT_cells_CD8s, 
                                       size = num_cells, color = assignments)) +
  annotate("rect",
           xmin = -10, xmax = 0,
           ymin = -Inf, ymax = 70,
           fill = "lightgrey", alpha = 0.4) +
  annotate("rect",
           xmin = 0, xmax = 10,
           ymin = -Inf, ymax = Inf,
           fill = "lightgrey", alpha = 0.4) +
  geom_vline(xintercept = 0, color = "black", linetype = "dashed", linewidth = 0.5) +  # Vertical line at x = 0
  geom_vline(xintercept = 10, color = "#2ca02c", linetype = "dashed") +
  geom_vline(xintercept = -10, color = "#2ca02c", linetype = "dashed") +
  geom_point(alpha = 0.5) +
  ggrepel::geom_label_repel(aes(label = Secondary.final_clusters),
                            size = 3.5,
                            fill = "white",
                            color = "black",
                            max.overlaps = Inf,
                            box.padding = 0.5,
                            segment.color = "grey30",
                            segment.size = 0.3,
                            min.segment.length = 0,
                            nudge_y = 0.1,
                            nudge_x = 0,
                            inherit.aes = TRUE) +
  scale_size_continuous(range = c(1, 10)) +
  annotate("segment", x = 0, xend = -10,
           y = 70, yend = 70, color = "#ff9900", linetype = "dashed") +
  scale_color_manual(values = c(
    "Assigned on fold change" = "#2ca02c",
    "Assigned on CD4-hashtag" = "#ff9900",
    "Indeterminate" = "#e31a1c")) +
  guides(color = guide_legend(override.aes = list(size = 5))) +  # Increase legend dot size
  theme_bw() +
  theme(axis.line.x = element_line(linewidth = 0.5, color = "black"),
        axis.title.x = element_text(hjust = 0.21, size = 14),
        axis.title.y = element_text(size = 14),
        axis.text = element_text(size = 12),
        legend.title = element_text(size = 14),
        legend.text = element_text(size = 12),
        plot.title = element_text(face = "bold", size = 18, hjust = 0.035)) +
  labs(
    x = "CD4-CD8B bias (± fold change)",
    y = "% CD4-hashtag expression (background-corrected)",
    size = "Number of cells",
    color = "Assignment"
  ) +
  ggtitle("CD4+ clusters ←→ CD8B+ clusters") +
  scale_x_continuous(breaks = seq(-40, 60, by = 10))


ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/cd4cd8_subclustering/subclusterPlot.png",
       subclusterPlot,
       height = 6.5, width = 10.5)

plotData <- plotData %>%
  arrange(assignments)

green <- as.character(plotData[[1]][1:11])
amber <- as.character(plotData[[1]][12:24])
red <- as.character(plotData[[1]][25:28])

Lseu$assignments <- NA
Lseu$assignments[Lseu$Secondary.final_clusters %in% green] <- "Assigned on fold change"
Lseu$assignments[Lseu$Secondary.final_clusters %in% amber] <- "Assigned on CD4-hashtag"
Lseu$assignments[Lseu$Secondary.final_clusters %in% red] <- "Excluded"
Lseu$assignments[is.na(Lseu$assignments)] <- "NA"



p <- DimPlot(Lseu, reduction = "Secondary.harm_umap", group.by = "assignments", alpha = 0.3, label = F) +
  scale_color_manual(values = c(
  "Assigned on fold change" = "#2ca02c",
  "Assigned on CD4-hashtag" = "#ff9900",
  "Excluded" = "#e31a1c",
  "NA" = "grey80")) & NoLegend()

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/cd4cd8_subclustering/assignments_umap.png",
       p,
       height = 5.5, width = 5.8)

# CD4-hashtag by cluster and extracting colours from umap ####

umap_plot <- DimPlot(Lseu, reduction = "Secondary.harm_umap", group.by = "Secondary.harm_clusters")
cluster_ids <- levels(factor(Lseu$Secondary.harm_clusters))

# Get colours from DimPlot manually
umap_plot <- DimPlot(Lseu, reduction = "Secondary.harm_umap", group.by = "Secondary.harm_clusters")
umap_colours <- ggplot_build(umap_plot)$data[[1]]$colour
names(umap_colours) <- umap_plot$data$Secondary.harm_clusters

# Deduplicate to get one colour per cluster ID
umap_colours_named <- umap_colours[!duplicated(names(umap_colours))]

col_list <- c("7" = "#00C0AF", "2" = "#C99800", "9" = "#00B0F6",
              "0" = "#F8766D", "3" = "#A3A500", "1" = "#E58700",
              "8" = "#00BCD8", "5" = "#00BA38", "4" = "#6BB100",
              "6" = "#00BF7D", "10" = "#619CFF", "12" = "#E76BF3",
              "11" = "#B983FF", "14" = "#FF67A4", "13" = "#FD61D1")

order <- rev(c("0", "1_0", "1_1", "1_2", "1_3", "2",
               "3_0", "3_1", "3_2", "3_3", "4", "5_0",
               "5_1", "6_0", "6_1", "6_2", "7_0", "7_1", "8_0", "8_1", "8_2",
               "9_0", "9_1", "11_0", "11_1", "12_1", "14_0", "14_1"))

# Extract base cluster numbers
base_numbers <- sub("_.*", "", order)

# Map colors
mapped_colours <- rev(setNames(col_list[base_numbers], order))

tils <- subset(Lseu,
               subset = Secondary.final_clusters %in% order)

tils$Secondary.final_clusters <- factor(tils$Secondary.final_clusters,
                                        levels = order)
cd4HT_subcluster_ridge <- RidgePlot(
  tils,
  features = "CD4-hashtag",
  sort = FALSE,
  group.by = "Secondary.final_clusters",
  cols = mapped_colours,
) +
  theme(
    plot.title = element_text(size = 14,
                              hjust = 0.5,
                              face = "bold"),
    axis.title.x = element_text(size = 12,
                                hjust = 0.5,
                                face = "bold"),
    axis.title.y = element_text(size = 12,
                                hjust = 0.5,
                                face = "bold")) +
  NoLegend() +
  ggtitle("CD4-hashtag expression by cluster") +
  geom_vline(xintercept = 0.8, color = "red", linetype = "dashed", linewidth = 0.5)

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/cd4cd8_subclustering/cd4HT_subcluster_ridge.png",
       cd4HT_subcluster_ridge,
       height = 9, width = 4)

# Make secondary label assignments ####

cd8s <- c("0", "2", "4","6_0", "7_0", "11_0", "12_1", "14_0")
cd4s <- c("1_0", "1_1", "1_2", "1_3", "3_0", "3_1", "5_0", "5_1", "6_1", "8_0", "8_1", "8_2", "9_0", "9_1", "11_1", "14_1")
ind <- c("3_2", "3_3", "6_2", "7_1", "10")

Lseu$Secondary.labels <- NA
Lseu$Secondary.labels[Lseu$Secondary.final_clusters == "13"] <- "B cells"
Lseu$Secondary.labels[Lseu$Secondary.final_clusters == "12_0"] <- "NK cells"
Lseu$Secondary.labels[Lseu$Secondary.final_clusters %in% cd8s] <- "CD8+ T cells"
Lseu$Secondary.labels[Lseu$Secondary.final_clusters %in% cd4s] <- "CD4+ T cells"
Lseu$Secondary.labels[Lseu$Secondary.final_clusters %in% ind] <- "Other/indeterminate"
Lseu$Secondary.labels <- as.character(Lseu$Secondary.labels)
Lseu$Secondary.labels[Lseu$Secondary.labels == "Indeterminate"] <- "Other/indeterminate"
# Plot by Secondary.labels

Lseu$Secondary.labels <- factor(
  Lseu$Secondary.labels,
  levels = c(
    "CD8+ T cells",
    "CD4+ T cells",
    "NK cells",
    "B cells",
    "Other/indeterminate"
  )
)

cols <- c("#F8766D", "#00C0AF", "#FD61D1", "#E58700",  "#B983FF")

p <- DimPlot(Lseu,
             reduction = "Secondary.harm_umap",
             group.by = "Secondary.labels",
             label = F, label.size = 4, cols = cols) +
  ggtitle("Final lymphoid cell type assignments") +
  labs(x = "umap_1", y = "umap_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/cd4cd8_subclustering/cd4cd8_final.png",
       p,
       height = 5.5, width = 7.5)

# Explore gene expression in excluded subclusters ####

clusters <- c("3_2", "3_3", "6_2", "7_1")
genes <- c("CD3E", "TRBC1", "TRAC", "CD8B", "KLRC4", "NKG7", "CD4", "TNFRSF18", "CD40LG")
excluded <- subset(Lseu, Secondary.final_clusters %in% clusters)

p <- DimPlot(excluded,
             reduction = "Secondary.harm_umap",
             group.by = "Secondary.final_clusters") +
  ggtitle("Indeterminate subclusters") +
  labs(x = "umap_1", y = "umap_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))


ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/cd4cd8_subclustering/indeterminate_umap.png",
       p,
       height = 4, width = 5.5)

features <- FeaturePlot(excluded,
                        reduction = "Secondary.harm_umap",
                        features = genes) +
  labs(x = "umap_1", y = "umap_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
features <- lapply(features, function(p) {
  p +
    labs(x = "umap_1", y = "umap_2") +
    theme(
      axis.text = element_blank()
      )
  })
features <- wrap_plots(features)

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/cd4cd8_subclustering/indeterminate_features.png",
       features,
       height = 6, width = 8)

# Excluded T/myeloid indeterminate cluster:

indeterminate <- subset(Lseu, Secondary.harm_clusters == "10")

p <- DimPlot(indeterminate,
             reduction = "Secondary.harm_umap",
             group.by = "Secondary.harm_clusters", cols = "#619CFF", label=F) +
  ggtitle("Cluster 10") +
  labs(x = "umap_1", y = "umap_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) & NoLegend()


ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/clusterTen_umap.png",
       p,
       height = 4, width = 5)

genes <- c("CD3E","TRBC1", "IL7R", "APOE", "C1QC", "TYROBP")

features <- FeaturePlot(indeterminate,
                        reduction = "Secondary.harm_umap",
                        features = genes) +
  labs(x = "umap_1", y = "umap_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
features <- lapply(features, function(p) {
  p +
    labs(x = "umap_1", y = "umap_2") +
    theme(
      axis.text = element_blank()
    )
})
features <- wrap_plots(features)

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/clusterTen_features.png",
       features,
       height = 5, width = 9)

# Explore proportions of secondary labels ####

bySample <- Lseu@meta.data %>%
  select(MULTI_ID, Secondary.labels) %>%
  group_by(MULTI_ID, Secondary.labels) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = Secondary.labels, values_from = count, values_fill = 0) 

global_lymphoid <- bySample %>%
  mutate(All = `B cells` + `CD4+ T cells` + `CD8+ T cells` + `Other/indeterminate` + `NK cells`,
         p.b = `B cells` / All * 100,
         p.cd4 = `CD4+ T cells` / All * 100,
         p.cd8 = `CD8+ T cells` / All * 100,
         p.ind = `Other/indeterminate` / All * 100,
         p.nk = `NK cells` / All * 100,
         Donor = str_extract(MULTI_ID, "^.{3}"),
         Timepoint = ifelse(grepl("Primary", MULTI_ID), "Primary", "Recurrence1"),
         Tissue = ifelse(grepl("Tumour", MULTI_ID), "Tumour", "PBZ")) %>%
  select(MULTI_ID, Donor, Timepoint, Tissue, p.b, p.cd4, p.cd8, p.ind, p.nk, All) %>%
  filter(Donor != "N01",
         Timepoint == "Primary")

t_props <- bySample %>%
  mutate(
    All = `CD4+ T cells` + `CD8+ T cells`,
    p.cd4 = `CD4+ T cells` / All * 100,
    p.cd8 = `CD8+ T cells` / All * 100,
    Donor = str_extract(MULTI_ID, "^.{3}"),
    Timepoint = ifelse(grepl("Primary", MULTI_ID), "Primary", "Recurrence1"),
    Tissue = ifelse(grepl("Tumour", MULTI_ID), "Tumour", "PBZ")) %>%
  select(MULTI_ID, Donor, Timepoint, Tissue, p.cd4, p.cd8, All) %>%
  filter(Donor != "N01",
         Timepoint == "Primary")

primary <- global_lymphoid %>%
  group_by(Timepoint) %>%
  summarise(medb = median(p.b),
            lowerb = quantile(p.b, 0.25),
            higherb = quantile(p.b, 0.75),
            medcd4 = median(p.cd4),
            lowercd4 = quantile(p.cd4, 0.25),
            highercd4 = quantile(p.cd4, 0.75),
            medcd8 = median(p.cd8),
            lowercd8 = quantile(p.cd8, 0.25),
            highercd8 = quantile(p.cd8, 0.75),
            mednk = median(p.nk),
            lowernk = quantile(p.nk, 0.25),
            highernk = quantile(p.nk, 0.75))

primary_t <- t_props %>%
  group_by(Timepoint) %>%
  summarise(medcd4 = median(p.cd4),
            lowercd4 = quantile(p.cd4, 0.25),
            highercd4 = quantile(p.cd4, 0.75),
            medcd8 = median(p.cd8),
            lowercd8 = quantile(p.cd8, 0.25),
            highercd8 = quantile(p.cd8, 0.75)
            )

model <- lmer(p.nk ~ Tissue + (1 | Donor), data = global_lymphoid)
summary(model)

# Plot proportions of secondary labels by region ####

byRegion <- Lseu@meta.data %>%
  select(MULTI_Region, Secondary.labels) %>%
  group_by(MULTI_Region, Secondary.labels) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = Secondary.labels, values_from = count, values_fill = 0) %>%
  mutate(All = `B cells` + `CD4+ T cells` + `CD8+ T cells` + `Other/indeterminate` + `NK cells`,
         p.b = `B cells` / All * 100,
         p.cd4 = `CD4+ T cells` / All * 100,
         p.cd8 = `CD8+ T cells` / All * 100,
         p.ind = `Other/indeterminate` / All * 100,
         p.nk = `NK cells` / All * 100) %>%
  select(MULTI_Region, p.b, p.cd4, p.cd8, p.ind, p.nk) %>%
  pivot_longer(cols = c(p.b, p.cd4, p.cd8, p.ind, p.nk), names_to = "Type", values_to = "Proportion") %>%
  mutate(Type = factor(Type, levels = rev(c("p.cd8", "p.cd4", "p.nk", "p.b", "p.ind"))),
         MULTI_Region = case_when(
           MULTI_Region == "Primary-Tumour" ~ "Primary Tumour",
           MULTI_Region == "Primary-PBZ" ~ "Primary PBZ",
           MULTI_Region == "Recurrence1-Tumour" ~ "Recurrent Tumour",
           MULTI_Region == "Recurrence1-PBZ" ~ "Recurrent PBZ",
           TRUE ~ MULTI_Region
         ),
         MULTI_Region = factor(MULTI_Region, levels = c("Primary Tumour", "Primary PBZ", "Recurrent Tumour", "Recurrent PBZ"))
         )

# Plot
cols <- rev(c("#F8766D", "#00C0AF",  "#FD61D1", "#E58700",  "#B983FF"))
regionPlot <- ggplot(
  byRegion, aes(x = MULTI_Region, y = Proportion, fill = Type)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = cols,
                    labels = rev(c("CD8+ T cells", "CD4+ T cells", "NK cells", "B cells", "Other/indeterminate"))) +
  labs(x = "Sample", y = "Proportions by region", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold"),
        strip.text = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none")

# Legend - note plot incorrect
cols <- c("#F8766D", "#00C0AF",  "#FD61D1", "#E58700",  "#B983FF")
dummyRegion <- ggplot(
  byRegion, aes(x = MULTI_Region, y = Proportion, fill = Type)
) +
  geom_bar(stat = "identity") +
  scale_fill_manual(values = cols,
                    labels = c("CD8+ T cells", "CD4+ T cells", "NK cells", "B cells", "Other/indeterminate")) +
  labs(x = "Sample", y = "Lymphoid cell proportions", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold"),
        strip.text = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.title = element_text(hjust = 0.5, face = "bold"))

dummyRegion
region_legend <- get_legend(dummyRegion)

grid.newpage()
grid.draw(region_legend)

output_dir <- "GBM_single_cell_analysis/outputs/clustering/Lymphoid/"
ggsave(paste0(output_dir, "region_lymphoid_proportions.png"), regionPlot,
       height = 4, width = 1.5)
ggsave(paste0(output_dir, "region_lymphoid_proportions_legend.png"), region_legend,
       height = 4, width = 4)

props <- data.frame(table(Lseu$MULTI_ID, Lseu$Secondary.labels)) %>%
  pivot_wider(names_from = Var2, values_from = Freq) %>%
  select(Var1, `CD4+ T cells`, `CD8+ T cells`) %>%
  mutate(T = `CD4+ T cells` + `CD8+ T cells`,
         Donor = str_extract(Var1, "^.{3}"),
         Timepoint = ifelse(grepl("Primary", Var1), "Primary", "Recurrence1"),
         Tissue = ifelse(grepl("Tumour", Var1), "Tumour", "PBZ")) %>%
  filter(T > 50,
         Timepoint == "Primary")

props$Tissue <- factor(props$Tissue, levels = c("PBZ", "Tumour"))
model <- glmer(cbind(`CD8+ T cells`, T - `CD8+ T cells`) ~ Tissue + (1 | Donor),
               family = binomial,
               data = props)
summary(model)

# Plot proportions by donor ####

byDonor <- Lseu@meta.data %>%
  select(MULTI_Donor, Secondary.labels) %>%
  group_by(MULTI_Donor, Secondary.labels) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = Secondary.labels, values_from = count, values_fill = 0) %>%
  mutate(All = `B cells` + `CD4+ T cells` + `CD8+ T cells` + `Other/indeterminate` + `NK cells`,
         p.b = `B cells` / All * 100,
         p.cd4 = `CD4+ T cells` / All * 100,
         p.cd8 = `CD8+ T cells` / All * 100,
         p.ind = `Other/indeterminate` / All * 100,
         p.nk = `NK cells` / All * 100) %>%
  select(MULTI_Donor, p.b, p.cd4, p.cd8, p.ind, p.nk) %>%
  pivot_longer(cols = c(p.b, p.cd4, p.cd8, p.ind, p.nk), names_to = "Type", values_to = "Proportion") %>%
  mutate(Type = factor(Type, levels = rev(c("p.cd8", "p.cd4", "p.nk", "p.b", "p.ind"))),
         MULTI_Donor = factor(MULTI_Donor, levels = c("N03", "N02", "N06", "N07", "N08", "N05"))) %>%
  filter(!is.na(MULTI_Donor))

# Plot
cols <- rev(c("#F8766D", "#00C0AF",  "#FD61D1", "#E58700", "#B983FF"))
donorPlot <- ggplot(
  byDonor, aes(x = MULTI_Donor, y = Proportion, fill = Type)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = cols,
                    labels = rev(c("CD8+ T cells", "CD4+ T cells", "NK cells", "B cells", "Other/indeterminate"))) +
  labs(x = "Sample", y = "Proportions by donor", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold", size = 14),
        strip.text = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none")

ggsave(paste0(output_dir, "donor_lymphoid_proportions.png"), donorPlot,
       height = 4, width = 2.5)

# Re-integrate for tertiary clustering ####

cd8s <- subset(Lseu,
               subset = Lseu$Secondary.labels == "CD8+ T cells")

cd4s <- subset(Lseu,
               subset = Lseu$Secondary.labels == "CD4+ T cells")

subseu <- list(
  CD8 = cd8s,
  CD4 = cd4s
)

seuNames <- names(subseu)

for (obj in seuNames){
  
  cat("Processing", obj, "\n")
  
  seu <- subseu[[obj]]
  seu[["RNA"]] <- split(seu[["RNA"]], f = seu$orig.ident)
  seu <- BatchProcess(seu)
  
  cat("Re-integrating", "\n")
  
  seu <- IntegrateLayers(
    object = seu,
    method = HarmonyIntegration,
    orig.reduction = "pca",
    new.reduction = "harmony_tertiary",
    verbose = FALSE)
  
  elbow <- ElbowPlot(seu, reduction = "pca", ndims = 50) +
    ggtitle(paste0("Percentage variance by principal component (", obj, ")")) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
      axis.title = element_text(size = 12,face = "bold")
    )
  
  output_dir <- paste0("GBM_single_cell_analysis/outputs/clustering/Lymphoid/", obj, "/")
  dir.create(output_dir, recursive=T, showWarnings=F)
  
  ggsave(paste0(output_dir, obj, "_elbow_plot.png"),
         elbow,
         height = 6, width = 6)
  
  cat("Saving...", "\n")
  
  saveRDS(seu, paste0("GBM_single_cell_analysis/outputs/objects/", obj, ".RDS"))
  
  subseu[[obj]] <- seu
  
  rm(elbow, seu)
  
  cat("Complete", "\n")
}

# Tertiary clustering ####

for (obj in seuNames){
  
  cat("Processing", obj, "cells", "\n")
  
  seu <- subseu[[obj]]
  seu <- JoinLayers(seu)
  
  #pcs <- seq(10, 30, by = 5)
  #resolutions <- seq(0.2, 1.2, by = 0.2)
  pcs <- 15
  res <- 0.6
  for (pc in pcs){
    cat("Constructing SNN graph,", pc, "dimensions", "\n")
    seu <- FindNeighbors(seu,
                         reduction = "harmony_tertiary",
                         dims = 1:pc,
                         graph.name = paste0("harmony_tertiary_snn_pc", pc))
    
    cat("Generating UMAP,", pc, "dimensions", "\n")
    seu <- RunUMAP(seu,
                   reduction = "harmony_tertiary",
                   dims = 1:pc,
                   reduction.name = paste0("umap.harm_tertiary_pc", pc),
                   graph.name = paste0("harmony_tertiary_snn_pc", pc))
    
    for (res in resolutions){
      cat("Clustering,", pc, "dimensions, resolution", res, "\n")
      
      plots_pdf <- list()
      
      output_dir <- paste0("GBM_single_cell_analysis/outputs/clustering/Lymphoid/", obj, "/pc", pc, "/pc", pc, "_res", res, "/")
      dir.create(output_dir, recursive=T, showWarnings=F)
      
      seu <- FindClusters(seu,
                          resolution = res,
                          graph.name = paste0("harmony_tertiary_snn_pc", pc),
                          cluster.name = paste0("harm_tertiary_clusters_pc", pc, "_res", res))
      
      p <- DimPlot(seu,
                   reduction = paste0("umap.harm_tertiary_pc", pc),
                   group.by = paste0("harm_tertiary_clusters_pc", pc, "_res", res),
                   label = TRUE, label.size = 4.5) +
        ggtitle(paste0(obj, " clusters (pc", pc, "_res", res, ")")) +
        labs(x = "umap_1", y = "umap_2") +
        guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
      
      plots_pdf[[1]] <- p
      
      ggsave(paste0(output_dir, obj, "_clusters_pc", pc, "_res", res, ".png"),
             p,
             height = 8.5, width = 10)
      
      cat("Plotting annotations...", "\n")
      
      annotations <- c("SingleR.labels", "CelliD.labels", "ScType.labels", "GBMap_3.labels", "GBMap_4.labels")
      
      for (annot in seq_along(annotations)){
        p <- DimPlot(seu,
                     reduction = paste0("umap.harm_tertiary_pc", pc),
                     group.by = paste0(annotations[[annot]]),
                     label = TRUE, label.size = 3) +
          ggtitle(paste0(obj, " ", annotations[[annot]], " (pc", pc, "_res", res, ")")) +
          labs(x = "umap_1", y = "umap_2") +
          guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
        
        plots_pdf[[(annot + 1)]] <- p
        
        ggsave(paste0(output_dir, obj, "_", annotations[[annot]], "_pc", pc, "_res", res, ".png"),
               p,
               height = 8.5, width = 11)
        
      }
      
      cat("Running differential gene expression analysis for generated clusters,", pc, "dimensions, resolution", res, "\n")
      
      markers <- FindAllMarkers(seu, 
                                logfc.threshold = 0.5, 
                                min.pct = 0.1, 
                                only.pos = TRUE)
      
      clusterMarkers <- markers %>%
        mutate(
          diff_pct = (pct.1 - pct.2) * 100,
          signif = -log10(p_val_adj),
          signif = ifelse(is.infinite(signif) | signif > 300, 300, signif)) %>%
        group_by(cluster) %>%
        arrange(cluster, desc(diff_pct)) %>%
        slice_head(n = 20) %>%
        ungroup()
      
      m.data <- seu@meta.data %>%
        arrange(seurat_clusters)
      
      clusters <- unique(m.data$seurat_clusters)
      cluster_plotlist <- list()
      
      cat("Plotting clusters...", "\n")
      
      for (i in seq_along(clusters)){
        cluster_id <- clusters[[i]]
        cluster_cells <- length(which(seu$seurat_clusters == cluster_id))
        cluster_size <- round(cluster_cells / nrow(seu@meta.data) * 100, 2)
        
        topmarkers <- clusterMarkers %>%
          filter(cluster == cluster_id) %>%
          arrange(desc(diff_pct))
        
        p <- ggplot(topmarkers, aes(x = diff_pct, y = reorder(gene, diff_pct), size = signif)) +
          # Trailing lines from y-axis to dots
          geom_segment(aes(x = 0, xend = diff_pct, y = reorder(gene, diff_pct), yend = reorder(gene, diff_pct)),
                       color = "gray", size = 0.5, alpha = 0.5) +  # Line settings
          
          # Dot plot
          geom_point(aes(color = avg_log2FC, x = diff_pct, y = gene), alpha = 1) +  
          labs(title = paste0(cluster_id, "; n = ", cluster_cells, "; ", cluster_size, "%")) +
          scale_color_gradient(
            low = "#132132", 
            high = "steelblue1", 
            limits = range(clusterMarkers$avg_log2FC, na.rm = TRUE),
            name = "Fold change"
          ) +
          scale_size_continuous(name = "Significance",
                                range = c(0.05, 2),
                                limits = range(clusterMarkers$signif, na.rm = TRUE)) +  # Adjust dot size range
          xlim(0, 100) +
          theme_cowplot() +
          theme(axis.text.y = element_text(size = 4.5),
                axis.text.x = element_text(size = 5.5),
                axis.title.x = element_blank(),
                axis.title.y = element_blank(),
                plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
                legend.title = element_text(size = 5, face = "bold"),
                legend.text = element_text(size = 5),
                legend.position = "none",
                plot.background = element_rect(fill = "white", color = NA))
        
        cluster_plotlist[[i]] <- p
        
      }
      
      legend_plot <- ggplot(topmarkers, aes(x = diff_pct, y = reorder(gene, diff_pct), size = signif)) +
        geom_segment(aes(x = 0, xend = diff_pct, y = reorder(gene, diff_pct), yend = reorder(gene, diff_pct)),
                     color = "gray", size = 0.5, alpha = 0.5) +
        geom_point(aes(color = avg_log2FC, x = diff_pct, y = gene), alpha = 1) +  
        labs(title = paste0(cluster_id, "; n = ", cluster_cells, "; ", cluster_size, "%")) +
        scale_color_gradient(
          low = "#132132", 
          high = "steelblue1", 
          limits = range(clusterMarkers$avg_log2FC, na.rm = TRUE),
          name = "Fold change"
        ) +
        scale_size_continuous(name = "Significance",
                              range = c(0.2, 3),
                              limits = range(clusterMarkers$signif, na.rm = TRUE)) +  # Adjust dot size range
        xlim(0, 100) +
        theme_cowplot() +
        theme(axis.text = element_text(size = 5), 
              axis.title = element_text(size = 7, face = "bold"),
              plot.title = element_text(size = 8, face = "bold", hjust = 0.5),
              legend.title = element_text(size = 10, face = "bold"),
              legend.text = element_text(size = 10),
              legend.position = "right",
              plot.background = element_rect(fill = "white", color = NA))
      
      legend <- get_legend(legend_plot)
      
      cluster_grid <- plot_grid(plotlist = cluster_plotlist, ncol = 5)
      plot <- plot_grid(cluster_grid, legend,
                        rel_widths = c(5,1), nrow = 1, align = "vh")
      final_plot <- ggdraw(plot) +
        draw_label("%(cluster) - %(rest)", x = 0.5, y = 0.01, vjust = 0.5, fontface = "bold", size = 14) +
        draw_label("Gene", x = 0.01, y = 0.5, angle = 90, vjust = 0.5, fontface = "bold", size = 14)
      
      plots_pdf[[(annot + 1 + 1)]] <- final_plot
      
      pdf_path <- paste0(output_dir, "outputs_pc", pc, "_res", res, ".pdf")
      pdf(pdf_path,
          height = 12, width = 14)
      for (p in plots_pdf) {
        print(p)
      }
      dev.off()
    }
  }
  
  subseu[[obj]] <- seu
  
  rm(cluster_grid, cluster_plotlist, clusterMarkers, final_plot, legend,
     legend_plot, m.data, markers, topmarkers, p, plot, plots_pdf, seu)
  
  cat(obj, "clustering complete", "\n")
}

# Tertiary clustree ####

fixed_dim <- as.list(seuNames)
fixed_res <- as.list(seuNames)

for (obj in seuNames){
  dims <- seq(10, 30, by = 5)
  fixed_dim[[obj]] <- list()
  
  output_dir <- paste0("GBM_single_cell_analysis/outputs/clustering/clustree/Lymphoid/", obj, "/")
  dir.create(output_dir, recursive=T, showWarnings=F)
  
  for (dim in dims){
    p <- clustree(subseu[[obj]], prefix = paste0("harm_tertiary_clusters_pc", dim, "_res")) 
    fixed_dim[[obj]][[dim]] <- p
  }
  
  pdf_path <- paste0(output_dir, obj,"_fixed_dim.pdf")
  pdf(pdf_path,
      height = 10, width = 14)
  for (p in fixed_dim) {
    print(p)
  }
  dev.off()
  
  resolutions <- seq(0.2, 1.2, by = 0.2)
  fixed_res[[obj]] <- list()
  
  for (res in resolutions){
    for (dim in dims){
      copy_colname <- paste0("harm_tertiary_clusters_pc", dim, "_res", res)
      new_col <- paste0("harm_tertiary_clusters_res", res, "_pc", dim)
      subseu[[obj]]@meta.data[[new_col]] <-  subseu[[obj]]@meta.data[[copy_colname]]
    }
    
    p <- clustree(subseu[[obj]], prefix = paste0("harm_tertiary_clusters_res", res, "_pc"))
    fixed_res[[obj]][[as.character(res)]] <- p
    
    for (dim in dims){
      new_col <- paste0("pc", dim)
      subseu[[obj]]@meta.data[[new_col]] <- NULL
    }
    
  }
  
  pdf_path <- paste0(output_dir, obj, "_fixed_res.pdf")
  pdf(pdf_path,
      height = 10, width = 14)
  for (p in fixed_res[[obj]]) {
    print(p)
  }
  
  dev.off()
  
}

rm(fixed_dim, fixed_res, p)

# Remove unused clusters/reductions ####

for (obj in seuNames){
  seu <- subseu[[obj]]
  
  seu@meta.data <- seu@meta.data[, !grepl("^harm_tertiary_clusters_res", colnames(seu@meta.data))]
  Idents(seu) <- seu$harm_tertiary_clusters_pc15_res0.6
  seu$seurat_clusters <- seu$harm_tertiary_clusters_pc15_res0.6
  
  seu$Tertiary.harm_clusters <- seu$harm_tertiary_clusters_pc15_res0.6
  seu@meta.data <- seu@meta.data %>%
    select(-starts_with("harm_tertiary_clusters"))
  
  seu@reductions$Tertiary.harm_umap <- seu@reductions$umap.harm_tertiary_pc15
  to_remove <- grep("^umap\\.harm_tertiary", names(seu@reductions), value = TRUE)
  seu@reductions[to_remove] <- NULL
  
  seu@graphs$Tertiary.harm_snn <- seu@graphs$harmony_tertiary_snn_pc15
  to_remove <- grep("^harmony_tertiary", names(seu@graphs), value = TRUE)
  seu@graphs[to_remove] <- NULL
  
  subseu[[obj]] <- seu
}

cd8s <- subseu$CD8
cd4s <- subseu$CD4

# Annotation ####

GBMap_T <- read.xlsx("GBM_single_cell_analysis/outputs/clustering/GBMap_TcellModules.xlsx", colNames=F)
colnames(GBMap_T) <- c("cellName", "marker")

GBMap_T <- GBMap_T %>%
  group_by(cellName) %>%
  summarise(geneSymbolmore1 = paste(marker, collapse = ","), .groups = "drop")

GBMap_T$tissueType <- "Glioblastoma"
GBMap_T$geneSymbolmore2 <- NA

GBMap_T <- GBMap_T %>%
  select(tissueType, cellName, geneSymbolmore1, geneSymbolmore2)

write.xlsx(GBMap_T, "GBM_single_cell_analysis/outputs/clustering/GBMap_TcellModulesDB.xlsx")

cd8s <- run_sctype(cd8s, assay = "RNA", scaled = TRUE, known_tissue_type = "Glioblastoma", custom_marker_file = "GBM_single_cell_analysis/outputs/clustering/GBMap_TcellModulesDB.xlsx", name = "GBMap_TcellModules")
cd4s <- run_sctype(cd4s, assay = "RNA", scaled = TRUE, known_tissue_type = "Glioblastoma", custom_marker_file = "GBM_single_cell_analysis/outputs/clustering/GBMap_TcellModulesDB.xlsx", name = "GBMap_TcellModules")

PBT <- read.xlsx("GBM_single_cell_analysis/outputs/clustering/PBT_sigs.xlsx", colNames=F)
colnames(PBT) <- c("cellName", "marker")

PBT <- PBT %>%
  group_by(cellName) %>%
  summarise(geneSymbolmore1 = paste(marker, collapse = ","), .groups = "drop")

PBT$tissueType <- "Glioblastoma"
PBT$geneSymbolmore2 <- NA

PBT <- PBT %>%
  select(tissueType, cellName, geneSymbolmore1, geneSymbolmore2)

write.xlsx(PBT, "GBM_single_cell_analysis/outputs/clustering/PBT_sigsDB.xlsx")

cd8s <- run_sctype(cd8s, assay = "RNA", scaled = TRUE, known_tissue_type = "Glioblastoma", custom_marker_file = "GBM_single_cell_analysis/outputs/clustering/PBT_sigsDB.xlsx", name = "PBT_sigs")
cd4s <- run_sctype(cd4s, assay = "RNA", scaled = TRUE, known_tissue_type = "Glioblastoma", custom_marker_file = "GBM_single_cell_analysis/outputs/clustering/PBT_sigsDB.xlsx", name = "PBT_sigs")

cd8s <- run_sctype(cd8s, assay = "RNA", scaled = TRUE, known_tissue_type = "Glioblastoma", custom_marker_file = "GBM_single_cell_analysis/outputs/clustering/GBMap_4DB.xlsx", name = "GBMap_4.labels")
cd4s <- run_sctype(cd4s, assay = "RNA", scaled = TRUE, known_tissue_type = "Glioblastoma", custom_marker_file = "GBM_single_cell_analysis/outputs/clustering/GBMap_4DB.xlsx", name = "GBMap_4.labels")

# Explore ####

p <- DimPlot(cd8s,
             reduction = "Tertiary.harm_umap",
             group.by = "MULTI_Region",
             label = F, label.size = 4) +
  ggtitle("CD8s grouped by tumour region") +
  labs(x = "umap_1", y = "umap_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/cd8_region.png",
       p,
       height = 6.2, width = 7.75)

features <- FeaturePlot(cd8s, reduction = "Tertiary.harm_umap",
            features = c("CD3E", "CD8A", "NKG7", "CCL4", "GZMK", "PDCD1", "HLA-DRB1", "ITGAE", "CD69", "IFIT3",
                         "FCGR3A", "FGFBP2","KLRB1", "IL7R", "HSPA6", "TOX")) +
  labs(x = "umap_1", y = "umap_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
features <- lapply(features, function(p) {
  p + theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
})
features <- wrap_plots(features)

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/features.png",
       features,
       height = 6, width = 12)


clones <- DimPlot(cd8s, reduction = "Tertiary.harm_umap", group.by = "cloneSize")
ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clones.png",
       clones,
       height = 3, width = 5.5)

# DGEA ####

subseu <- list(CD8 = cd8s, CD4 = cd4s)

for (obj in seuNames){
  seu <- subseu[[obj]]
  
  markers <- FindAllMarkers(seu, 
                            logfc.threshold = 0.5, 
                            min.pct = 0.1, 
                            only.pos = TRUE)
  
  clusterMarkers <- markers %>%
    mutate(
      diff_pct = (pct.1 - pct.2) * 100,
      signif = -log10(p_val_adj),
      signif = ifelse(is.infinite(signif) | signif > 300, 300, signif)) %>%
    group_by(cluster) %>%
    arrange(cluster, desc(diff_pct)) %>%
    slice_head(n = 20) %>%
    ungroup()
  
  m.data <- seu@meta.data %>%
    arrange(seurat_clusters)
  
  clusters <- unique(m.data$seurat_clusters)
  cluster_plotlist <- list()
  
  cat("Plotting clusters...", "\n")
  
  for (i in seq_along(clusters)){
    cluster_id <- clusters[[i]]
    cluster_cells <- length(which(seu$seurat_clusters == cluster_id))
    cluster_size <- round(cluster_cells / nrow(seu@meta.data) * 100, 2)
    
    topmarkers <- clusterMarkers %>%
      filter(cluster == cluster_id) %>%
      arrange(desc(avg_log2FC))
    
    p <- ggplot(topmarkers, aes(x = avg_log2FC, y = reorder(gene, avg_log2FC), size = signif)) +
      # Trailing lines from y-axis to dots
      geom_segment(aes(x = 0, xend = avg_log2FC, y = reorder(gene, avg_log2FC), yend = reorder(gene, avg_log2FC)),
                   color = "gray", size = 0.5, alpha = 0.5) +  # Line settings
      
      # Dot plot
      geom_point(aes(color = diff_pct, x = avg_log2FC, y = gene), alpha = 1) +  
      labs(title = paste0(cluster_id, "; n = ", cluster_cells, "; ", cluster_size, "%")) +
      scale_color_gradient(
        low = "#132132", 
        high = "steelblue1", 
        limits = range(clusterMarkers$diff_pct, na.rm = TRUE),
        name = "Fold change"
      ) +
      scale_size_continuous(name = "Significance",
                            range = c(0.05, 2),
                            limits = range(clusterMarkers$signif, na.rm = TRUE)) +  # Adjust dot size range
      xlim(0, 10) +
      theme_cowplot() +
      theme(axis.text.y = element_text(size = 4.5),
            axis.text.x = element_text(size = 5.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
            legend.title = element_text(size = 5, face = "bold"),
            legend.text = element_text(size = 5),
            legend.position = "none",
            plot.background = element_rect(fill = "white", color = NA))
    
    cluster_plotlist[[i]] <- p
    
  }
  
  legend_plot <- ggplot(topmarkers, aes(x = avg_log2FC, y = reorder(gene, avg_log2FC), size = signif)) +
    # Trailing lines from y-axis to dots
    geom_segment(aes(x = 0, xend = avg_log2FC, y = reorder(gene, avg_log2FC), yend = reorder(gene, avg_log2FC)),
                 color = "gray", size = 0.5, alpha = 0.5) +  # Line settings
    
    # Dot plot
    geom_point(aes(color = diff_pct, x = avg_log2FC, y = gene), alpha = 1) +  
    labs(title = paste0(cluster_id, "; n = ", cluster_cells, "; ", cluster_size, "%")) +
    scale_color_gradient(
      low = "#132132", 
      high = "steelblue1", 
      limits = range(clusterMarkers$diff_pct, na.rm = TRUE),
      name = "Differential percentage expression"
    ) +
    scale_size_continuous(name = "Significance",
                          range = c(0.05, 2),
                          limits = range(clusterMarkers$signif, na.rm = TRUE)) +  # Adjust dot size range
    xlim(0, 10) +
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 4.5),
          axis.text.x = element_text(size = 5.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 5, face = "bold"),
          legend.text = element_text(size = 5),
          legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA))
  
  legend <- get_legend(legend_plot)
  
  cluster_grid <- plot_grid(plotlist = cluster_plotlist, ncol = 4)
  plot <- plot_grid(cluster_grid, legend,
                    rel_widths = c(5,1), nrow = 1, align = "vh")
  final_plot <- ggdraw(plot) +
    draw_label("log2FC", x = 0.5, y = 0.01, vjust = 0.5, fontface = "bold", size = 14) +
    draw_label("Gene", x = 0.01, y = 0.5, angle = 90, vjust = 0.5, fontface = "bold", size = 14)
  

  
  pdf(paste0("GBM_single_cell_analysis/outputs/clustering/Lymphoid/", obj, "/", obj, "_logFC_degs.pdf"),
      height = 4, width = 8)
  print(final_plot)
  dev.off()

  subseu[[obj]] <- seu
}
# CD8 tertiary labels ####

cd8s$Tertiary.labels <- NA
cd8s$Tertiary.labels[cd8s$Tertiary.harm_clusters == "0"] <- "GZMK-high Tem"
cd8s$Tertiary.labels[cd8s$Tertiary.harm_clusters == "1"] <- "Activated Tcm"
cd8s$Tertiary.labels[cd8s$Tertiary.harm_clusters == "2"] <- "Resting memory"
cd8s$Tertiary.labels[cd8s$Tertiary.harm_clusters == "3"] <- "IFN response"
cd8s$Tertiary.labels[cd8s$Tertiary.harm_clusters == "4"] <- "TOX-high"
cd8s$Tertiary.labels[cd8s$Tertiary.harm_clusters == "5"] <- "HSP-high"
cd8s$Tertiary.labels[cd8s$Tertiary.harm_clusters == "6"] <- "NK-like effector"
cd8s$Tertiary.labels[cd8s$Tertiary.harm_clusters == "7"] <- "Cycling"

# Plot CD8s by Tertiary.labels
cd8s$Tertiary.labels <- factor(
  cd8s$Tertiary.labels,
  levels = c(
    "GZMK-high Tem",
    "Activated Tcm",
    "Resting memory",
    "IFN response",
    "TOX-high",
    "HSP-high",
    "NK-like effector",
    "Cycling"
  )
)
p <- DimPlot(cd8s,
             reduction = "Tertiary.harm_umap",
             group.by = "Tertiary.labels",
             label = T, label.size = 5) +
  ggtitle("Cluster analysis of CD8+ T cells") +
  labs(x = "umap_1", y = "umap_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) & NoLegend()


ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/cd8s_umap.png",
       p,
       height = 5.5, width = 7)

# CD8 clusters - visualise gene expression ####

# Heatmap

all_markers_pct10 <- FindAllMarkers(object = cd8s, logfc.threshold = 0.5, min.pct = 0.1, only.pos = TRUE)
all_markers_pct10 <- all_markers_pct10 %>%
  filter(p_val_adj < 0.05)
top_markers_pct10 <- all_markers_pct10 %>%
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n=50)

all_markers_pct50 <- FindAllMarkers(object = cd8s, logfc.threshold = 0.5, min.pct = 0.5, only.pos = TRUE)
top_markers_pct50 <- all_markers_pct50 %>%
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n=50)
top_markers_pct50 <- all_markers_pct50 %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n=50)

plot_markers <- top_markers_pct10 %>%
  pull(gene) %>%
  unique()

palette <- colorRampPalette(
  c("blue", 
    "black",
    "yellow"
  )
)(100)

heatmap <- DoHeatmap(cd8s,
          features = plot_markers,
          group.by = "Tertiary.harm_clusters") +
  scale_fill_gradientn(colors = palette) +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "left")

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clusters_heatmap.png",
       heatmap,
       height = 6, width = 6)

genes <- c("CD3E", "CD3D", "CD3G", "CD8A", "CD8B", "TRBC1", "TRAC", # CD8s
           "GZMK", "GZMA", "GZMB", "GZMH", "PRF1", "IFNG", "TNF", # Effector
           "CD69", "HLA-DRB1", "JUN", "FOS", "JUNB", "FOSB", "DUSP2", # early activation
           "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "TOX", "BATF", # activation/exhaustion
            "ITGAE", "MYO7A", "GPR25", "ZNF683", "KRT86", "RBPJ", # residency
           "ATP10D", "ENTPD1", "KIR2DL4", "LAYN", "HTRA1", "CD70", # CD8 neoantigen reactivity (Lowery)
            "NKG7", "FGFBP2", "CX3CR1", "FCGR3A", "ADGRG1", "PTGDS", "KLRB1", # NK-like
           "IL7R", "GPR183", "LMNA", "NR4A3", "TCF7", "MGAT4A", "CD55", "EGR1") #Tcm

palette2 <- colorRampPalette(
  c("blue",    
    "#E6E0EB",
    "firebrick3" 
  )
)(100)

orig_colours <- scales::hue_pal()(length(unique(cd8s$Tertiary.harm_clusters)))
names(orig_colours) <- levels(cd8s$Tertiary.harm_clusters)

dotplot <- Clustered_DotPlot(seurat_object = cd8s,
                  features = genes, 
                  cluster_feature = FALSE,
                  cluster_ident = FALSE,
                  colors_use_ident = orig_colours,
                  plot_km_elbow = FALSE,
                  colors_use_exp = palette2)

png("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/cd8_clusters_dotplot.png",
    width = 1500, height = 2500, res = 300)
ComplexHeatmap::draw(dotplot)
dev.off()


# CD8 GO term analysis ####

# DGEA

Idents(cd8s) <- cd8s$Tertiary.labels

markers <- FindAllMarkers(cd8s, 
                          logfc.threshold = 0.5, 
                          min.pct = 0.1, 
                          only.pos = TRUE)

# Create a named list of gene vectors
gene_list <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  arrange(desc(avg_log2FC)) %>%
  summarise(genes = list(gene)) %>%
  deframe()

# Convert gene symbols to Entrez IDs
gene_list_entrez <- lapply(gene_list, function(genes) {
  bitr(genes, fromType = "SYMBOL", toType = "ENTREZID", OrgDb = org.Hs.eg.db)$ENTREZID
})

# Background genes
universe_ids <- bitr(unique(markers$gene), fromType = "SYMBOL",
                     toType = "ENTREZID", OrgDb = org.Hs.eg.db)

comp <- compareCluster(geneCluster = gene_list_entrez,
                       fun = "enrichGO",
                       universe = universe_ids$ENTREZID,
                       OrgDb = org.Hs.eg.db,
                       keyType = "ENTREZID",
                       ont = "BP",
                       pAdjustMethod = "BH",
                       pvalueCutoff = 0.05,
                       qvalueCutoff = 0.05)

go_results_df <- as.data.frame(comp) %>%
  arrange(Cluster, desc(RichFactor))

write.xlsx(go_results_df, "GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/cd8_go_terms.xlsx")

# Plot
go_dot <- dotplot(comp, showCategory = 3) +
  theme(axis.text.x = element_text(size = 8.5, angle = 45, hjust = 1, vjust = 1),
        axis.text.y = element_text(size = 7.5),
        axis.title.x = element_blank())
go_dot

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/cd8_go_dot.png",
       go_dot,
       height = 10, width = 7)

# Plot proportions of CD8 clusters by region ####

byRegion <- cd8s@meta.data %>%
  dplyr::select(MULTI_Region, Tertiary.labels) %>%
  group_by(MULTI_Region, Tertiary.labels) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = Tertiary.labels, values_from = count, values_fill = 0) %>%
  mutate(All = `GZMK-high Tem` +
           `Activated Tcm` +
           `Resting memory` +
           `IFN response` +
           `TOX-high` +
           `HSP-high` +
           `NK-like effector` +
           `Cycling`,
         p.0 = `GZMK-high Tem` / All * 100,
         p.1 = `Activated Tcm` / All * 100,
         p.2 = `Resting memory` / All * 100,
         p.3 = `IFN response` / All * 100,
         p.4 = `TOX-high` / All * 100,
         p.5 = `HSP-high` / All * 100,
         p.6 = `NK-like effector` / All * 100,
         p.7 = `Cycling` / All * 100) %>%
  dplyr::select(MULTI_Region, p.0, p.1, p.2, p.3, p.4, p.5, p.6, p.7) %>%
  pivot_longer(cols = c(p.0, p.1, p.2, p.3, p.4, p.5, p.6, p.7), names_to = "Type", values_to = "Proportion") %>%
  mutate(Type = factor(Type, levels = rev(c("p.0", "p.1", "p.2", "p.3", "p.4", "p.5", "p.6", "p.7"))),
         MULTI_Region = case_when(
           MULTI_Region == "Primary-Tumour" ~ "Primary Tumour",
           MULTI_Region == "Primary-PBZ" ~ "Primary PBZ",
           MULTI_Region == "Recurrence1-Tumour" ~ "Recurrent Tumour",
           MULTI_Region == "Recurrence1-PBZ" ~ "Recurrent PBZ",
           TRUE ~ MULTI_Region
         ),
         MULTI_Region = factor(MULTI_Region, levels = c("Primary Tumour", "Primary PBZ", "Recurrent Tumour", "Recurrent PBZ"))
  )

# Plot

cols <- rev(scales::hue_pal()(length(unique(cd8s$Tertiary.labels))))

regionPlot <- ggplot(
  byRegion, aes(x = MULTI_Region, y = Proportion, fill = Type)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = cols,
                    labels = rev(c("GZMK-high Tem", "Activated Tcm", "Resting memory", "IFN response",
                                   "TOX-high", "HSP-high", "NK-like effector", "Cycling"))) +
  labs(x = "Sample", y = "Proportions by region", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold", size = 12),
        strip.text = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none")

# Legend - note plot incorrect
cols <- scales::hue_pal()(length(unique(cd8s$Tertiary.labels)))
dummyRegion <- ggplot(
  byRegion, aes(x = MULTI_Region, y = Proportion, fill = Type)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = cols,
                    labels = c("GZMK-high Tem", "Activated Tcm", "Resting memory", "IFN response",
                                   "TOX-high", "HSP-high", "NK-like effector", "Cycling")) +
  labs(x = "Sample", y = "Proportions by region", fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold"),
        strip.text = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.title = element_text(hjust = 0.5, face = "bold"))


dummyRegion
region_legend <- get_legend(dummyRegion)

grid.newpage()
grid.draw(region_legend)

output_dir <- "GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/"
ggsave(paste0(output_dir, "region_cd8_cluster_proportions.png"), regionPlot,
       height = 4, width = 1.5)
ggsave(paste0(output_dir, "region_legend.png"), region_legend,
       height = 4, width = 4)

# Plot proportions of CD8 clusters by donor ####

byDonor <- cd8s@meta.data %>%
  select(MULTI_Donor, Tertiary.labels) %>%
  group_by(MULTI_Donor, Tertiary.labels) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = Tertiary.labels, values_from = count, values_fill = 0) %>%
  mutate(All = `GZMK-high Tem` +
           `Activated Tcm` +
           `Resting memory` +
           `IFN response` +
           `TOX-high` +
           `HSP-high` +
           `NK-like effector` +
           `Cycling`,
         p.0 = `GZMK-high Tem` / All * 100,
         p.1 = `Activated Tcm` / All * 100,
         p.2 = `Resting memory` / All * 100,
         p.3 = `IFN response` / All * 100,
         p.4 = `TOX-high` / All * 100,
         p.5 = `HSP-high` / All * 100,
         p.6 = `NK-like effector` / All * 100,
         p.7 = `Cycling` / All * 100) %>%
  select(MULTI_Donor, p.0, p.1, p.2, p.3, p.4, p.5, p.6, p.7) %>%
  pivot_longer(cols = c(p.0, p.1, p.2, p.3, p.4, p.5, p.6, p.7), names_to = "Type", values_to = "Proportion") %>%
  mutate(Type = factor(Type, levels = rev(c("p.0", "p.1", "p.2", "p.3", "p.4", "p.5", "p.6", "p.7"))),
         MULTI_Donor = factor(MULTI_Donor, levels = c("N03", "N02", "N06", "N07", "N08", "N05"))) %>%
  filter(!is.na(MULTI_Donor))

# Plot

cols <- rev(scales::hue_pal()(length(unique(cd8s$Tertiary.labels))))

donorPlot <- ggplot(
  byDonor, aes(x = MULTI_Donor, y = Proportion, fill = Type)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = cols,
                    labels = rev(c("GZMK-high Tem", "Activated Tcm", "Resting memory", "IFN response",
                                   "TOX-high", "HSP-high", "NK-like effector", "Cycling"))) +
  labs(x = "Sample", y = "Proportions by donor", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold", size = 14),
        strip.text = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none")

ggsave(paste0(output_dir, "donor_cd8_cluster_proportions.png"), donorPlot,
       height = 4, width = 2.5)


# Cluster forest plot by region ####

primary <- cd8s@meta.data %>%
  filter(grepl("Primary", MULTI_Region))

table(primary$Tertiary.labels)
table(primary$MULTI_Region, primary$Tertiary.labels)

cd8data <- cd8s@meta.data %>%
  select(MULTI_ID, Tertiary.labels) %>%
  group_by(MULTI_ID, Tertiary.labels) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = Tertiary.labels, values_from = count, values_fill = 0) %>%
  mutate(
    All = `GZMK-high Tem` +
      `Activated Tcm` +
      `Resting memory` +
      `IFN response` +
      `TOX-high` +
      `HSP-high` +
      `NK-like effector` +
      `Cycling`,
    Donor = str_extract(MULTI_ID, "^.{3}"),
    Timepoint = ifelse(grepl("Primary", MULTI_ID), "Primary", "Recurrence1"),
    Tissue = ifelse(grepl("Tumour", MULTI_ID), "Tumour", "PBZ")) %>%
  select(MULTI_ID, Donor, Timepoint, Tissue, everything()) 

primary <- cd8data %>%
  filter(Timepoint == "Primary")

celltypes <- c("GZMK-high Tem",
               "Activated Tcm",
               "Resting memory",
               "IFN response",
               "TOX-high",
               "HSP-high",
               "NK-like effector",
               "Cycling")

min_nonzero <- min(primary %>% select(any_of(celltypes)) %>% unlist() %>% .[. > 0])
primary_clr_input <- primary[, celltypes] /100 # switch back to proportions
min_nonzero <- min(primary_clr_input[primary_clr_input > 0])
pseudocount <- min_nonzero / 100
primary_clr_adj <- primary_clr_input
primary_clr_adj[primary_clr_adj == 0] <- pseudocount
primary_clr_adj <- primary_clr_adj / rowSums(primary_clr_adj)

clr_data <- clr(primary_clr_adj)
clr_df <- cbind(primary[, c("Donor", "Tissue")], as.data.frame(clr_data))

sub_df <- data.frame(
  CellType = character(length(celltypes)),
  Estimate = numeric(length(celltypes)),
  SE = numeric(length(celltypes)),
  p_value = numeric(length(celltypes)),
  Upper = numeric(length(celltypes)),
  Lower = numeric(length(celltypes)),
  Donor = numeric(length(celltypes)),
  stringsAsFactors = FALSE
)

for (i in seq_along(celltypes)) {
  var <- celltypes[i]
  model_formula <- as.formula(paste0("`", var, "` ~ Tissue + (1 | Donor)"))
  model <- lmer(model_formula, data = clr_df)
  
  est <- fixef(model)["TissueTumour"]
  se <- summary(model)$coefficients["TissueTumour", "Std. Error"]
  donor <- attr(summary(model)$varcor$Donor, "stddev")
  
  sub_df[i, "CellType"] <- var
  sub_df[i, "Estimate"] <- est
  sub_df[i, "SE"] <- se
  sub_df[i, "p_value"] <- ifelse(summary(model)$coefficients["TissueTumour", "Pr(>|t|)"] < 0.0001, "p < 0.0001", signif(summary(model)$coefficients["TissueTumour", "Pr(>|t|)"], 2))
  sub_df[i, "Upper"] <- est + 1.96 * se
  sub_df[i, "Lower"] <- est - 1.96 * se
  sub_df[i, "Donor"] <- donor
}

sub_df$FoldChange <- ifelse(sub_df$Estimate < 0, -exp(abs(sub_df$Estimate)), exp(sub_df$Estimate))
sub_df$Upper_FC <- ifelse(sub_df$Upper < 0, -exp(abs(sub_df$Upper)), exp(sub_df$Upper))
sub_df$Lower_FC <- ifelse(sub_df$Lower < 0, -exp(abs(sub_df$Lower)), exp(sub_df$Lower))
sub_df$Donor_FC <- exp(sub_df$Donor)

sub_df <- sub_df %>%
  arrange(FoldChange) %>%
  mutate(CellType = factor(CellType, levels = unique(CellType)))

write.xlsx(sub_df, "GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/cd8s_region_forest_models.xlsx")

orig_colours <- scales::hue_pal()(length(unique(cd8s$Tertiary.labels)))
cols <- c(orig_colours[5+1], # after inspecting result; 
          orig_colours[7+1], # each cluster number + 1 for position in orig_colours
          orig_colours[2+1],
          orig_colours[3+1],
          orig_colours[4+1],
          orig_colours[6+1],
          orig_colours[1+1],
          orig_colours[0+1]
)
          
subDiff <- ggplot(sub_df, aes(x = CellType, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(color = rev(cols), size = 5) +
  geom_linerange(aes(ymin = Lower, ymax = Upper),
                 color = rev(cols), size = 1) +
  geom_text(
    aes(label = paste0("p = ", p_value), y = max(Upper) + 0.3), hjust = 0,
    color = "black",
    size = 5
  ) +
  coord_flip() +
  labs(y = "Centered log-ratio change", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.35),
    plot.background = element_rect(fill = "white")) +
  scale_y_continuous(limits = c(min(sub_df$Lower), max(sub_df$Upper) + 1),
                     breaks = seq(-4, 4, 1),
                     labels = seq(-4, 4, 1)
  )

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/cd8s_region_forest.png",
       subDiff,
       height = 3.5, width = 8.5)

# Cluster forest plot by donor ####

min_nonzero <- min(cd8data %>% select(any_of(celltypes)) %>% unlist() %>% .[. > 0])
primary_clr_input <- cd8data[, celltypes] /100 # switch back to proportions
min_nonzero <- min(primary_clr_input[primary_clr_input > 0])
pseudocount <- min_nonzero / 100
primary_clr_adj <- primary_clr_input
primary_clr_adj[primary_clr_adj == 0] <- pseudocount
primary_clr_adj <- primary_clr_adj / rowSums(primary_clr_adj)

clr_data <- clr(primary_clr_adj)
clr_df <- cbind(cd8data[, c("Donor", "Tissue", "Timepoint")], as.data.frame(clr_data))

sub_df <- data.frame(
  CellType = character(length(celltypes)),
  Estimate = numeric(length(celltypes)),
  SE = numeric(length(celltypes)),
  p_value = character(length(celltypes)),
  Upper = numeric(length(celltypes)),
  Lower = numeric(length(celltypes)),
  Tissue_estimate = numeric(length(celltypes)),
  stringsAsFactors = FALSE
)

for (i in seq_along(celltypes)) {
  var <- celltypes[i]
  model_formula <- as.formula(paste0("`", var, "` ~ Timepoint + Tissue"))
  model <- lm(model_formula, data = clr_df)
  
  coef_summary <- summary(model)$coefficients
  
  est <- coef_summary["TimepointRecurrence1", "Estimate"]
  se <- coef_summary["TimepointRecurrence1", "Std. Error"]
  p_val <- coef_summary["TimepointRecurrence1", "Pr(>|t|)"]
  tissue <- coef_summary["TissueTumour", "Estimate"]
  
  sub_df[i, "CellType"] <- var
  sub_df[i, "Estimate"] <- est
  sub_df[i, "SE"] <- se
  sub_df[i, "p_value"] <- ifelse(p_val < 0.0001, "p < 0.0001", signif(p_val, 2))
  sub_df[i, "Upper"] <- est + 1.96 * se
  sub_df[i, "Lower"] <- est - 1.96 * se
  sub_df[i, "Tissue_estimate"] <- tissue
}

sub_df$FoldChange <- ifelse(sub_df$Estimate < 0, -exp(abs(sub_df$Estimate)), exp(sub_df$Estimate))
sub_df$Upper_FC <- ifelse(sub_df$Upper < 0, -exp(abs(sub_df$Upper)), exp(sub_df$Upper))
sub_df$Lower_FC <- ifelse(sub_df$Lower < 0, -exp(abs(sub_df$Lower)), exp(sub_df$Lower))
sub_df$Tissue_FC <- exp(sub_df$Tissue_estimate)

sub_df <- sub_df %>%
  arrange(FoldChange) %>%
  mutate(CellType = factor(CellType, levels = unique(CellType)))

write.xlsx(sub_df, "GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/cd8s_timepoint_forest_models.xlsx")

orig_colours <- scales::hue_pal()(length(unique(cd8s$Tertiary.labels)))
cols <- c(orig_colours[3+1], # after inspecting result; 
          orig_colours[0+1], # each cluster number + 1 for position in orig_colours
          orig_colours[7+1],
          orig_colours[6+1],
          orig_colours[4+1],
          orig_colours[5+1],
          orig_colours[2+1],
          orig_colours[1+1]
)

subDiff <- ggplot(sub_df, aes(x = CellType, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(color = rev(cols), size = 5) +
  coord_flip() +
  labs(y = "Centered log-ratio change", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.4),
    plot.background = element_rect(fill = "white")) +
  scale_y_continuous(limits = c(-3, 4),
                     breaks = seq(-3, 4, 1),
                     labels = seq(-3, 4, 1)
  )

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/cd8s_timepoint_forest.png",
       subDiff,
       height = 3.5, width = 8.5)

# Compare CD8 clonality by cluster ####
cd8s$clonalfreq_nt_cd8s <- table(cd8s$CTnt)[cd8s$CTnt]

cd8s$clonalFrequency_2 <- as.numeric(as.character(cd8s$clonalFrequency))

cd8s$clonalFrequency_2[is.na(cd8s$clonalFrequency_2)] <- 0

str(cd8s$clonalFrequency_2)

embeddings <- Embeddings(cd8s, reduction = "Tertiary.harm_umap")
df <- data.frame(
  x = embeddings[, 1],
  y = embeddings[, 2],
  clonalFrequency = cd8s$clonalFrequency_2,
  cluster = cd8s$Tertiary.labels
)

str(df$clonalFrequency)

p <- ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(size = clonalFrequency, color = cluster), alpha = 0.7) +
  scale_size_continuous(range = c(0.1, 5),
                        breaks = c(80, 50, 20, 5, 1, 0),
                        labels = c("80", "50", "20", "5", "1", "N/A")) +
  labs(x = "umap_1", y = "umap_2", size = "Clonal frequency", color = "Cluster") +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 4)))

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/clonal_umap.png",
       p,
       height = 5, width = 7)

orig_colours <- scales::hue_pal()(length(unique(cd8s$Tertiary.labels)))

p <- clonalQuant(cd8s,
                 cloneCall = "nt",
                 chain = "both",
                 group.by = "Tertiary.labels",
                 scale = T) +
  labs(y = "Unique clones (%)", fill = "Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = orig_colours)

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/clonalquant_clusters.png",
       p,
       height = 3, width = 3)

p <- clonalProportion(cd8s,
                      cloneCall = "nt",
                      group.by = "Tertiary.labels",
                      clonalSplit = c(1, 10, 50, 100, 500, 1000)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  scale_fill_brewer(
    name = "Clonal Indices",
    palette = "Set2",
    labels = c("Largest clone", "Clones 2–10", "Clones 11–50", 
               "Clones 51–100", "Clones 101–500", "Clones 501–1000")
  )

prop <- clonalProportion(cd8s,
                         cloneCall = "nt",
                         group.by = "Tertiary.labels",
                         clonalSplit = c(1, 10, 50, 100, 500, 1000),
                         exportTable=T) %>%
  data.frame() %>%
  mutate(total = rowSums(across(everything()))) %>%
  mutate(across(1:6, ~ .x / total * 100),
         `4` = sum(c_across(1:4)))
prop

isna <- cd8s@meta.data %>%
  filter(!is.na(cloneSize))
table(isna$Tertiary.labels)

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/clonalproportion_clusters.png",
       p,
       height = 3, width = 4.6)


cols <- c(orig_colours[1+1], # after inspecting result; 
          orig_colours[7+1], # each cluster number + 1 for position in orig_colours
          orig_colours[0+1],
          orig_colours[5+1],
          orig_colours[3+1],
          orig_colours[6+1],
          orig_colours[2+1],
          orig_colours[4+1])

metric_labels <- c(
  shannon = "Shannon Index",
  inv.simpson = "Inverse Simpson Index"
)

# Subset clusters with n > 100 for diversity

cluster_counts <- cd8s@meta.data %>%
  filter(!is.na(cloneSize)) %>%
  group_by(Tertiary.labels) %>%
  tally() %>%
  filter(n >= 100)

tcrs <- subset(cd8s, subset = Tertiary.labels %in% cluster_counts$Tertiary.labels)

cols <- c(orig_colours[1+1], # after inspecting result; 
          orig_colours[0+1], # each cluster number + 1 for position in orig_colours
          orig_colours[5+1],
          orig_colours[3+1],
          orig_colours[6+1],
          orig_colours[2+1])

set.seed(2025)
p <- clonalDiversity(tcrs,
                     cloneCall = "nt",
                     group.by = "Tertiary.labels",
                     metrics = c("shannon", "inv.simpson")) +
  scale_fill_manual(values = cols) +
  facet_wrap(~variable, scales = "free", ncol = length(c("shannon", "inv.simpson")), 
             labeller = labeller(variable = metric_labels)) +
  theme(legend.position = "none")
set.seed(2025)
p
ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/clonaldiversity_clusters.png",
       p,
       height = 3.6, width = 3.5)

set.seed(2025)
clonalDiversity(tcrs,
                cloneCall = "nt",
                group.by = "Tertiary.labels",
                metrics = c("shannon", "inv.simpson"),
                exportTable = T)

p <- clonalOverlap(cd8s, cloneCall = "nt", chain = "both", method = "morisita", group.by = "Tertiary.labels") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/clonaloverlap_clusters.png",
       p,
       height = 3.5, width = 5.5)

p <- clonalCompare(cd8s,
                   samples = c("GZMK-high Tem", "IFN response"),
                   top.clones = 10,
                   cloneCall="nt", 
                   graph = "alluvial",
                   relabel.clones = T) & NoLegend()
cols <- setNames(palette36.colors(16), levels(p$data$clones))
p <- p + scale_fill_manual(values = cols)


q <- clonalCompare(cd8s,
                   samples = c("Activated Tcm", "Resting memory"),
                   top.clones = 10,
                   cloneCall="nt", 
                   graph = "alluvial",
                   relabel.clones = T) & NoLegend()
cols <- setNames(palette36.colors(16), levels(q$data$clones))
q <- q + scale_fill_manual(values = cols)

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/compare_gzmk_ifn.png",
       p,
       height = 3, width = 4)
ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/compare_tcm_resting.png",
       q,
       height = 3, width = 4)

# DGEA between expanded and unexpanded cells ####

backup <- cd8s

expanded <- sort(unique(cd8s$clonalFrequency))
expanded <- expanded[-1]
expanded <- rev(expanded)
degs <- list()

for (exp in seq_along(expanded)){

  cd8s@meta.data <- cd8s@meta.data %>%
  mutate(clonalExp = case_when(
      clonalFrequency >= expanded[[exp]] ~ "Expanded",
      clonalFrequency == 1 ~ "Unexpanded",
      TRUE ~ NA_character_
    )
  )
  
  Idents(cd8s) <- cd8s$clonalExp
  
  clones <- length(unique(cd8s$CTnt[!is.na(cd8s$clonalFrequency) & cd8s$clonalFrequency >= expanded[[exp]]]))
  
  markers <- FindMarkers(cd8s,
                         ident.1 = "Expanded",
                         ident.2 = "Unexpanded",
                         logfc.threshold = 0.5,
                         min.pct = 0.1,
                         only.pos = T)
  markers <- markers %>%
    filter(p_val_adj <= 0.05) %>%
    mutate(clonotypes = clones,
           clonalFrequency = expanded[[exp]],
           cells = sum(cd8s$clonalFrequency == expanded[[exp]], na.rm = T),
           diff_pct = pct.1 - pct.2,
           signif = -log10(p_val_adj),
           signif = ifelse(is.infinite(signif) | signif > 300, 300, signif),
           gene = rownames(.)) %>%
    arrange(desc(avg_log2FC))
  
  degs <- rbind(degs, markers)
}

summary <- as.data.frame(table(degs$gene, degs$clonalFrequency)) %>%
  pivot_wider(names_from = Var2, values_from = Freq)

freq_by_clonotype <- degs %>%
  select(clonotypes, clonalFrequency, cells) %>%
  unique()

clones <- length(unique(cd8s$CTnt[!is.na(cd8s$clonalFrequency) & cd8s$clonalFrequency == 1]))
cells <- sum(cd8s$clonalFrequency == 1, na.rm = T)
add <- c(1244, 1, 1268)
freq_by_clonotype <- rbind(freq_by_clonotype, add) %>%
  mutate(cumul_cells = cumsum(cells))

freq_by_clonotype <- freq_by_clonotype[-26, ]

freq_by_clonotype <- freq_by_clonotype %>%
  select(clonalFrequency, cumul_cells)

degs <- degs %>%
  group_by(gene) %>%
  slice_min(order_by = clonalFrequency, n = 1) %>%
  ungroup() %>%
  arrange(clonotypes, desc(avg_log2FC))

degs <- left_join(degs, freq_by_clonotype, by = "clonalFrequency")

freq_by_clonotype <- degs %>%
  select(clonotypes, clonalFrequency, cells, cumul_cells) %>%
  unique()

palette_15 <- c(
  "#E6194B", # Bright Red
  "#3CB44B", # Medium Green
  "gold", # Bright Yellow
  "#4363D8", # Blue
  "#F58231", # Orange
  "#800000", # Purple
  "#46F0F0", # Cyan/Turquoise
  "#F032E6", # Magenta
  "#BCF60C", # Lime Green
  "#FABEBE", # Light Pink
  "#008080", # Teal
  "#E6BEFF", # Lavender
  "#9A6324", # Brown
  "black", # 
 "#911EB4"  # Maroon
)

p <- ggplot(freq_by_clonotype, aes(x = clonalFrequency, y = cumul_cells)) +
  geom_point(aes(color = factor(clonalFrequency))) +  # treat clonalFrequency as a factor
  labs(
    y = "Number of expanded cells analysed",
    x = "Clonal frequency cut-off",
    color = "Clonal Frequency"  # legend title
  ) +
  theme_bw() +
  scale_x_continuous(breaks = seq(0, max(freq_by_clonotype$clonalFrequency), by = 10)) +
  scale_y_continuous(breaks = seq(0, max(freq_by_clonotype$cumul_cells) + 100, by = 100)) +
  scale_color_manual(values = palette_15) & NoLegend()
p

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/cellnum_expanded_degs.png",
       p,
       height = 3.5, width = 5)


p <- ggplot(degs, aes(x = avg_log2FC, y = reorder(gene, -clonalFrequency), size = diff_pct, color = factor(clonalFrequency))) +
  geom_point(alpha = 1) +
  scale_color_manual(values = palette_15, name = "Clonal frequency cut-off") +
  # Trailing lines from y-axis to dots
  geom_segment(aes(x = 0, xend = avg_log2FC, y = reorder(gene, -clonalFrequency), yend = reorder(gene, -clonalFrequency)),
               color = "gray", size = 0.5, alpha = 0.5) +  # Line settings
  
  # Dot plot
  scale_size_continuous(name = "Differential prevalence",
                        range = c(1, 5),
                        limits = range(degs$diff_pct, na.rm = TRUE)) +  # Adjust dot size range
  xlim(0, 4) +
  labs(title = "DEGs between expanded and unexpanded CD8+ T cells", x = "Log2(fold change)", y = "Gene", color = "Clonal frequency cut-off") +
  theme_cowplot() +
  theme(axis.text.y = element_text(size = 8),
        axis.text.x = element_text(size = 10),
        axis.title.x = element_text(),
        axis.title.y = element_text(),
        plot.title = element_text(size = 12, face = "bold", hjust = 0.5),
        legend.title = element_text(size = 10, face = "bold"),
        legend.text = element_text(size = 10),
        legend.position = "right",
        plot.background = element_rect(fill = "white", color = NA))

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/expanded_degs.png",
       p,
       height = 5, width = 8)




# DGEA across regions by cluster ####

primary <- subset(cd8s,
                  grepl("Primary", MULTI_Region))
Idents(primary) <- "Tertiary.labels"
deseq2_markers <- list()
wilcoxon_markers <- list()
clusters <- levels(primary$Tertiary.labels)

for (cl in clusters) {
  cells <- WhichCells(primary, idents = cl)
  cluster <- subset(primary, cells = cells)
  
  cluster$MULTI_Region <- factor(cluster$MULTI_Region, levels = c("Primary-Tumour", "Primary-PBZ"))
  Idents(cluster) <- "MULTI_Region"
  
  # Pseudobulk by MULTI_Donor within region
  counts <- GetAssayData(cluster, slot = "counts", assay = "RNA")
  meta <- cluster@meta.data
  
  # Aggregate counts per donor-region pair
  meta$group_id <- paste(meta$MULTI_Donor, meta$MULTI_Region, sep = "_")
  
  pseudobulk <- as.data.frame(as.matrix(counts)) %>%
    t() %>%
    as.data.frame() %>%
    mutate(group_id = meta$group_id) %>%
    group_by(group_id) %>%
    summarise(across(everything(), sum)) %>%
    column_to_rownames("group_id") %>%
    t()
  
  # Create colData
  group_split <- strsplit(colnames(pseudobulk), "_")
  coldata <- data.frame(
    donor = sapply(group_split, `[`, 1),
    region = sapply(group_split, `[`, 2)
  )
  rownames(coldata) <- colnames(pseudobulk)
  coldata$region <- factor(coldata$region, levels = c("Primary-Tumour", "Primary-PBZ"))
  
  # Only proceed if both conditions are present and ≥10 cells
  region_cell_counts <- table(cluster$MULTI_Region)
  if (all(c("Primary-Tumour", "Primary-PBZ") %in% names(region_cell_counts)) && all(region_cell_counts >= 10)) {
    
    dds <- DESeqDataSetFromMatrix(countData = pseudobulk,
                                  colData = coldata,
                                  design = ~ region)
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("region", "Primary-Tumour", "Primary-PBZ"))
    
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df$cluster <- cl
    
   deseq2_markers[[cl]] <- res_df
    
    # Wilcoxon to compare
    
    markers <- FindMarkers(cluster,
                           ident.1 = "Primary-Tumour",
                           ident.2 = "Primary-PBZ",
                           logfc.threshold = 0.5,
                           min.pct = 0.1,
                           only.pos = FALSE)
    
    markers_filtered <- markers %>%
      mutate(gene = rownames(.)) %>%
      slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
      bind_rows(
        markers %>%
          mutate(gene = rownames(.)) %>%
          slice_min(order_by = avg_log2FC, n = 20, with_ties = FALSE)
      ) %>%
      mutate(cluster = cl)
    
    wilcoxon_markers[[cl]] <- markers_filtered
  }
}

deseq2_markers_df <- bind_rows(deseq2_markers) %>%
  group_by(cluster) %>%
  slice_max(order_by = log2FoldChange, n = 20) %>%
  bind_rows(
    bind_rows(deseq2_markers) %>%
      group_by(cluster) %>%
      slice_min(order_by = log2FoldChange, n = 20)
  )
  
wilcoxon_markers_df <- bind_rows(wilcoxon_markers) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  bind_rows(
    bind_rows(wilcoxon_markers) %>%
      group_by(cluster) %>%
      slice_min(order_by = avg_log2FC, n = 20)
  )

deseq2_bind <- deseq2_markers_df %>%
  mutate(log2FC = log2FoldChange,
         method = "DESeq2") %>%
  select(cluster, gene, log2FC, method) %>%
  arrange(cluster, desc(log2FC))

wilcoxon_bind <- wilcoxon_markers_df %>%
  mutate(log2FC = avg_log2FC,
         method = "Wilcoxon") %>%
  select(cluster, gene, log2FC, method) %>%
  arrange(cluster, desc(log2FC))

degs <- bind_rows(deseq2_bind, wilcoxon_bind) %>%
  group_by(cluster, gene) %>%
  mutate(max_log2FC = max(log2FC, na.rm = T)) %>%
  ungroup() %>%
  arrange(cluster, desc(max_log2FC), gene)

degs <- degs %>%
  filter(gene %in% table$gene) %>%
  arrange(cluster, desc(max_log2FC)) %>%
  group_by(cluster, gene) %>%
  filter(n() == 2,
         !(n() == 2 & method == "Wilcoxon")) %>%
  ungroup() %>%
  mutate(order = case_when(
    cluster == "GZMK-high Tem" ~ 1,
    cluster == "Activated Tcm" ~ 2,
    cluster == "Resting memory" ~ 3,
    cluster == "TOX-high" ~ 4,
    cluster == "NK-like effector" ~ 5
  )) %>%
  select(order, cluster, gene, log2FC) %>%
  arrange(order, desc(log2FC))

mat <- pivot_wider(degs, names_from = cluster, values_from = log2FC)

mat$`TOX-high`[32] <- 3.514376
mat$`TOX-high`[8] <- 2.785413

mat <- mat[-c(44, 51), ]

mat <- mat %>%
  column_to_rownames("gene") %>%
  select(-order) %>%
  as.matrix()
  
palette <- colorRampPalette(
  c("blue",    
    "#E6E0EB",
    "firebrick3" 
  )
)(100)

orig_colours <- scales::hue_pal()(length(unique(cd8s$Tertiary.harm_clusters)))

annotation_col <- data.frame(
  Cluster = colnames(mat)
)
rownames(annotation_col) <- colnames(mat)

cluster_colors <- setNames(orig_colours[c(1,2,3,5,7)], unique(annotation_col$Cluster))

hm <- pheatmap(mat,
         na_col = "grey90",
         color = palette,
         cluster_rows = FALSE,
         cluster_cols = FALSE,
         annotation_col = annotation_col,
         annotation_legend = F,
         annotation_colors = list(Cluster = cluster_colors))

pdf("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/cd8_region_degs.pdf",
    width = 3, height = 9)
draw(hm)
dev.off()

# CD8 trajectory analysis ####

# Global

Idents(cd8s) <- cd8s$Tertiary.labels

umap <- Embeddings(cd8s, reduction = "Tertiary.harm_umap")
clusters <- Idents(cd8s)
counts_mat <- cd8s[["RNA"]]@layers$counts

colnames(counts_mat) <- colnames(cd8s)

sce_global <- SingleCellExperiment(assays = list(counts = counts_mat))

# Add UMAP coordinates from Seurat
reducedDims(sce_global)$UMAP <- Embeddings(cd8s, reduction = "Tertiary.harm_umap")

# Add cluster labels
sce_global$cluster <- Idents(cd8s)

# Add MULTI_Region

multi_region <- cd8s$MULTI_Region
multi_region <- multi_region[colnames(sce_global)]
sce_global$MULTI_Region <- multi_region

# Run Slingshot
sce_global <- slingshot(sce_global, clusterLabels = 'cluster', reducedDim = 'UMAP')

# Pseudotime plot coloured by pseudotime

png("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/CD8_global_pseudotime.png",
    width = 1000, height = 700)

par(mar = c(5.1, 5.5, 4.1, 2.1))  # more space on left (2nd value)

ptime_matrix <- slingPseudotime(sce_global)

# Combine pseudotime across all lineages (minimum per cell)
ptime_combined <- apply(ptime_matrix, 1, function(x) min(x, na.rm = TRUE))

# Set pseudotime to NA where all lineages were NA
ptime_combined[!is.finite(ptime_combined)] <- NA

# Initialize color vector
colors <- rep("lightgray", length(ptime_combined))

# Only apply viridis scale to non-NA values
non_na_idx <- which(!is.na(ptime_combined))
colors[non_na_idx] <- viridis::viridis(100)[as.numeric(cut(ptime_combined[non_na_idx], breaks = 100))]

#plot
plot(reducedDims(sce_global)$UMAP,
     col = colors,
     pch = 16,
     asp = 1,
     xlab = "umap_1", 
     ylab = "umap_2",
     cex.lab = 2,
     cex.axis = 1.5)

lines(SlingshotDataSet(sce_global), lwd = 4, col = 'black')

dev.off()

# Coloured by cluster

png("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/CD8_global_slingshot.png",
    width = 1000, height = 700)

par(mar = c(5.1, 5.5, 4.1, 2.1))  # more space on left (2nd value)

colors <- orig_colours[as.numeric(factor(clusters))]
plot(reducedDims(sce_global)$UMAP, 
     col = colors, 
     pch = 16, 
     asp = 1, 
     xlab = "umap_1", 
     ylab = "umap_2",
     cex.lab = 2,
     cex.axis = 1.5)

lines(SlingshotDataSet(sce_global), lwd = 4, col = 'black')

dev.off()

# Coloured by region

png("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/CD8_global_region.png",
    width = 1000, height = 700)

par(mar = c(5.1, 5.5, 4.1, 2.1))  # more space on left (2nd value)

region_cols <- c("darkgreen", "darkseagreen1","orangered3", "peachpuff")
colors <- region_cols[as.numeric(factor(cd8s$MULTI_Region))]
plot(reducedDims(sce_global)$UMAP, 
     col = colors, 
     pch = 16, 
     asp = 1, 
     xlab = "umap_1", 
     ylab = "umap_2",
     cex.lab = 2,
     cex.axis = 1.5,
     cex = 1.8)

lines(SlingshotDataSet(sce_global), lwd = 4, col = 'black')

dev.off()

# Retrieve pseudotime data to Seurat


pseudotime <- slingPseudotime(sce_global)

pseudotime_min <- apply(pseudotime, 1, function(x) {
  if (all(is.na(x))) NA else min(x, na.rm = TRUE)
})

cd8s$pseudotime <- pseudotime_min
View(cd8s@meta.data)

# Focused

subset_cd8s <- subset(cd8s, idents = c("Activated Tcm", "Resting memory"))
umap <- Embeddings(subset_cd8s, reduction = "Tertiary.harm_umap")
clusters <- Idents(subset_cd8s)

counts_mat <- subset_cd8s[["RNA"]]@layers$counts
colnames(counts_mat) <- colnames(subset_cd8s)

sce <- SingleCellExperiment(assays = list(counts = counts_mat))

# Add UMAP coordinates from Seurat
reducedDims(sce)$UMAP <- Embeddings(subset_cd8s, reduction = "Tertiary.harm_umap")

# Add cluster labels
sce$cluster <- Idents(subset_cd8s)

# Add MULTI_Region

multi_region <- subset_cd8s$MULTI_Region
multi_region <- multi_region[colnames(sce)]
sce$MULTI_Region <- multi_region

# Run Slingshot
sce <- slingshot(sce, clusterLabels = 'cluster', reducedDim = 'UMAP')

# Pseudotime plot coloured by pseudotime

png("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/CD8_focused_pseudotime.png",
    width = 1000, height = 700)

par(mar = c(5.1, 5.5, 4.1, 2.1))  # more space on left (2nd value)

ptime <- slingPseudotime(sce)[,1]
colors_ptime <- viridis::viridis(n = length(ptime))[rank(ptime)]
plot(reducedDims(sce)$UMAP,
     col = colors_ptime,
     pch = 16,
     asp = 1,
     xlab = "umap_1", 
     ylab = "umap_2",
     cex.lab = 2,
     cex.axis = 1.5)

lines(SlingshotDataSet(sce), lwd = 4, col = 'black')

dev.off()

# Coloured by cluster

png("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/CD8_focused_slingshot.png",
    width = 1000, height = 700)

par(mar = c(5.1, 5.5, 4.1, 2.1))  # more space on left (2nd value)

colors <- orig_colours[2:3]
colors <- colors[as.numeric(factor(clusters))]
plot(reducedDims(sce)$UMAP, 
     col = colors, 
     pch = 16, 
     asp = 1, 
     xlab = "umap_1", 
     ylab = "umap_2",
     cex.lab = 2,
     cex.axis = 1.5)

lines(SlingshotDataSet(sce), lwd = 4, col = 'black')

dev.off()

# Coloured by region

png("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/CD8_focused_region.png",
    width = 1000, height = 700)

par(mar = c(5.1, 5.5, 4.1, 2.1))  # more space on left (2nd value)

region_cols <- c("darkgreen", "darkseagreen1","orangered3", "peachpuff")
colors <- region_cols[as.numeric(factor(subset_cd8s$MULTI_Region))]
plot(reducedDims(sce)$UMAP, 
     col = colors, 
     pch = 16, 
     asp = 1, 
     xlab = "umap_1", 
     ylab = "umap_2",
     cex.lab = 2,
     cex.axis = 1.5,
     cex = 1.8)

lines(SlingshotDataSet(sce), lwd = 4, col = 'black')

dev.off()

# Legend
ptime <- slingPseudotime(sce)[, 1]

png("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/pseudotime_legend.png",
    width = 200, height = 400)

par(mar = c(5, 4, 4, 6))  # margins: bottom, left, top, right

# Generate 100 pseudotime values from 0 to 1
legend_pts <- seq(0, 1, length.out = 100)
legend_cols <- viridis(100)

# Blank plot for legend space
plot(1, type = 'n', axes = FALSE, xlab = '', ylab = '',
     xlim = c(0, 1), ylim = c(0, 1))

# Draw color bar
rect(
  xleft = -0.5,
  ybottom = legend_pts[-length(legend_pts)],
  xright = 0.7,
  ytop = legend_pts[-1],
  col = legend_cols,
  border = NA
)

# Add axis and labels
axis(4, las = 1, at = seq(0, 1, by = 0.25), labels = seq(0, 1, by = 0.25), cex.axis = 1.4)

# Title for legend
mtext("Pseudotime", side = 3, line = 3, cex = 1.8)

dev.off()

# Retrieve pseudotime data to Seurat

pseudotime <- slingPseudotime(sce)

all(rownames(subset_cd8s@meta.data) == rownames(pseudotime)) # Check TRUE before proceeding
subset_cd8s$pseudotime <- pseudotime[,1]
View(subset_cd8s@meta.data)

# Vln Plot

cd8s$pseudotime <- (cd8s$pseudotime - min(cd8s$pseudotime, na.rm = TRUE)) / 
  (max(cd8s$pseudotime, na.rm = TRUE) - min(cd8s$pseudotime, na.rm = TRUE))

subset_cd8s$pseudotime <- (subset_cd8s$pseudotime - min(subset_cd8s$pseudotime, na.rm = TRUE)) / 
  (max(subset_cd8s$pseudotime, na.rm = TRUE) - min(subset_cd8s$pseudotime, na.rm = TRUE))

vln <- VlnPlot(cd8s,
               features = "pseudotime",
               group.by = "Tertiary.labels",
               pt.size = 0,
               cols = orig_colours) +
  labs(title = "Pseudotime across clusters") +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
    )

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/vln_pseudotime_clusters.png",
       vln,
       height = 5, width = 7)
    
subset_cd8s$region_cluster <- paste(subset_cd8s$MULTI_Region, subset_cd8s$Tertiary.labels, sep = "_")

cluster_order <- c(
  "Activated Tcm", 
  "Resting memory"
)

region_order <- c("Primary-PBZ", "Primary-Tumour", "Recurrence1-PBZ", "Recurrence1-Tumour")

ordered_levels <- as.vector(outer(region_order, cluster_order, paste, sep = "_"))

subset_cd8s$region_cluster <- factor(subset_cd8s$region_cluster, levels = ordered_levels)

cells_of_interest <- subset(cd8s, Tertiary.labels %in% c("Activated Tcm", "Resting memory"))

interest_cols <- c(orig_colours[2], orig_colours[2], orig_colours[2], orig_colours[2],
                   orig_colours[3], orig_colours[3], orig_colours[3], orig_colours[3])

vln <- VlnPlot(subset_cd8s,
               features = "pseudotime",
               group.by = "region_cluster",
               pt.size = 0,
               cols = interest_cols) +
  labs(title = "Pseudotime by region (activated Tcm/resting memory)") +
  scale_x_discrete(labels = c(
    "Primary-PBZ_Activated Tcm" = "Primary PBZ",
    "Primary-Tumour_Activated Tcm" = "Primary Tumour",
    "Recurrence1-PBZ_Activated Tcm" = "Recurrent PBZ",
    "Recurrence1-Tumour_Activated Tcm" = "Recurrent Tumour",
    "Primary-PBZ_Resting memory" = "Primary PBZ",
    "Primary-Tumour_Resting memory" = "Primary Tumour",
    "Recurrence1-PBZ_Resting memory" = "Recurrent PBZ",
    "Recurrence1-Tumour_Resting memory" = "Recurrent Tumour"
  )) +
  theme(
    legend.position = "none",
    axis.title.x = element_blank()
  )

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/vln_pseudotime_focused_regions.png",
       vln,
       height = 5, width = 7)

# Compare CD8 clonality across regions ####

region_cols <- c("darkgreen", "darkseagreen1","orangered3", "peachpuff")

p <- DimPlot(cd8s,
             reduction = "Tertiary.harm_umap",
             group.by = "MULTI_Region",
             label = F, repel = F, label.size = 5, pt.size = 2) +
  labs(x = "umap_1", y = "umap_2") +
  theme(
    plot.title = element_blank()
  ) +
  scale_color_manual(values = region_cols,
                    labels = c("Primary Tumour", "Primary PBZ", "Recurrent Tumour", "Recurrent PBZ"))
p
ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/regions_umap.png",
       p,
       height = 5, width = 7)

cd8s$MULTI_Region <- factor(cd8s$MULTI_Region,
                            levels = c("Primary-Tumour", "Primary-PBZ", "Recurrence1-Tumour", "Recurrence1-PBZ")
)

# clonalQuant data for significance testing:

table <- as.data.frame(
  clonalQuant(cd8s,
              cloneCall = "nt",
              chain = "both",
              group.by = "MULTI_ID",
              scale = T,
              exportTable = T)
) %>%
  mutate(Donor = str_extract(MULTI_ID, "^.{3}"),
         Timepoint = ifelse(grepl("Primary", MULTI_ID), "Primary", "Recurrence1"),
         Tissue = ifelse(grepl("Tumour", MULTI_ID), "Tumour", "PBZ"))

model <- glmer(cbind(contigs, total - contigs) ~ Tissue + (1 | Donor), family = binomial, data = table)
summary(model)

cd8s$MULTI_Region <- factor(cd8s$MULTI_Region,
                            levels = c("Primary-Tumour", "Primary-PBZ", "Recurrence1-Tumour", "Recurrence1-PBZ"),
                      )
# For clonalQuant plot

p <- clonalQuant(cd8s,
                 cloneCall = "nt",
                 chain = "both",
                 group.by = "MULTI_Region",
                 scale = T) +
  labs(y = "Unique clones (%)", fill = "Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_x_discrete(labels = c("Primary Tumour", "Primary PBZ", "Recurrent Tumour", "Recurrent PBZ")) +
  scale_fill_manual(values = region_cols,
                    labels = c("Primary Tumour", "Primary PBZ", "Recurrent Tumour", "Recurrent PBZ"))
p

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/clonalquant_regions.png",
       p,
       height = 3, width = 3)

legend_plot <- clonalQuant(cd8s,
                 cloneCall = "nt",
                 chain = "both",
                 group.by = "MULTI_Region",
                 scale = T) +
  labs(y = "Unique clones (%)", fill = "Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "right") +
  scale_x_discrete(labels = c("Primary Tumour", "Primary PBZ", "Recurrent Tumour", "Recurrent PBZ")) +
  scale_fill_manual(values = region_cols,
                    labels = c("Primary Tumour", "Primary PBZ", "Recurrent Tumour", "Recurrent PBZ"))

legend <- get_legend(legend_plot)
grid.newpage()
grid.draw(legend)

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/clonalquant_regions_legend.png",
       legend,
       height = 3, width = 3)
backup <- cd8s
cd8s$region_cluster <- paste(cd8s$MULTI_Region, cd8s$Tertiary.labels, sep = "_")
primarycd8s <- subset(cd8s, grepl("Primary", cd8s$MULTI_Region))
primarycd8s <- subset(primarycd8s, !is.na(cd8s$cloneSize))
unique(primarycd8s$region_cluster)

clonalDiversity(primarycd8s,
                cloneCall = "nt",
                group.by = "region_cluster",
                metrics = c("shannon", "inv.simpson"))
# Diversity

metric_labels <- c(
  shannon = "Shannon Index",
  inv.simpson = "Inverse Simpson Index"
)

set.seed(2025)
p <- clonalDiversity(cd8s,
                     cloneCall = "nt",
                     group.by = "MULTI_Donor",
                     metrics = c("shannon", "inv.simpson")) +
  scale_fill_manual(values = region_cols) +
  facet_wrap(~variable, scales = "free", ncol = length(c("shannon", "inv.simpson")), 
             labeller = labeller(variable = metric_labels)) +
  theme(legend.position = "none")

set.seed(2025)
p

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/clonaldiversity_regions.png",
       p,
       height = 3.5, width = 3.5)

# Overlap

p <- clonalOverlap(cd8s, cloneCall = "nt", chain = "both", method = "morisita", group.by = "MULTI_ID") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
p
ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/clonaloverlap_regions.png",
       p,
       height = 5, width = 7.5)

Idents(cd8s) <- cd8s$MULTI_ID
clonalCompare(cd8s,
              samples = c("N05-Recurrence1-Tumour1", "N05-Recurrence1-PBZ1"),
              top.clones = 10,
              cloneCall="nt", 
              graph = "alluvial",
              relabel.clones = T) & NoLegend()

p <- clonalCompare(cd8s,
                   samples = c("N08-Primary-Tumour1", "N08-Primary-Tumour4"),
                   top.clones = 10,
                   cloneCall="nt", 
                   graph = "alluvial",
                   relabel.clones = T) & NoLegend()
cols <- setNames(palette36.colors(18), levels(p$data$clones))
p <- p + scale_fill_manual(values = cols)

# Example for raw data:

table <- clonalCompare(cd8s,
                        samples = c("N08-Primary-Tumour1", "N08-Primary-Tumour4"),
                        top.clones = 10,
                        cloneCall="nt", 
                        graph = "alluvial",
                        relabel.clones = T,
              exportTable = T)

sum(table$Proportion[table$Sample == "N08-Primary-Tumour4"])

q <- clonalCompare(cd8s,
                   samples = c("N05-Recurrence1-Tumour1", "N05-Recurrence1-PBZ1"),
                   top.clones = 10,
                   cloneCall="nt", 
                   graph = "alluvial",
                   relabel.clones = T) & NoLegend()
cols <- setNames(palette36.colors(18), levels(q$data$clones))
q <- q + scale_fill_manual(values = cols)

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/compare_n08.png",
       p,
       height = 3, width = 4)
ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/compare_n05.png",
       q,
       height = 3, width = 4)


# TCR aa positions

Idents(cd8s) <- cd8s$MULTI_Region

p <- positionalEntropy(cd8s, 
          chain = "TRB", 
          aa.length = 20,
          method = "shannon") +
  scale_color_manual(values = region_cols,
                    labels = c("Primary Tumour", "Primary PBZ", "Recurrent Tumour", "Recurrent PBZ")) +
  NoLegend()

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/posEnt_region.png",
       p,
       height = 3, width = 4)

Idents(cd8s) <- cd8s$MULTI_Donor

q <- percentKmer(cd8s, 
            cloneCall = "aa",
            chain = "TRB", 
            motif.length = 4, 
            top.motifs = 10) +
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust = 1))
q
ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/pKmer_region.png",
       q,
       height = 3, width = 4)



viz <- vizGenes(cd8s, 
         x.axis = "TRBV",
         y.axis = NULL,
         plot = "barplot",  
         scale = TRUE) +
  theme(axis.text.x = element_text(size = 6, angle = 45, hjust = 1, vjust = 1))

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD8/clonal/cd8_vizgenes.png",
       viz,
       height = 5.6, width = 5)

# CD4 tertiary labels ####

cd4s$Tertiary.labels <- NA
cd4s$Tertiary.labels[cd4s$Tertiary.harm_clusters == "0"] <- "CD4-CTL"
cd4s$Tertiary.labels[cd4s$Tertiary.harm_clusters == "1"] <- "Resting memory"
cd4s$Tertiary.labels[cd4s$Tertiary.harm_clusters == "2"] <- "Early activated"
cd4s$Tertiary.labels[cd4s$Tertiary.harm_clusters == "3"] <- "Tcm"
cd4s$Tertiary.labels[cd4s$Tertiary.harm_clusters == "4"] <- "Th17"
cd4s$Tertiary.labels[cd4s$Tertiary.harm_clusters == "5"] <- "Treg"
cd4s$Tertiary.labels[cd4s$Tertiary.harm_clusters == "6"] <- "HSP-high"
cd4s$Tertiary.labels[cd4s$Tertiary.harm_clusters == "7"] <- "Cycling"

# Plot CD4s by Tertiary.labels
cd4s$Tertiary.labels <- factor(
  cd4s$Tertiary.labels,
  levels = c(
    "CD4-CTL",
    "Resting memory",
    "Early activated",
    "Tcm",
    "Th17",
    "Treg",
    "HSP-high",
    "Cycling"
  )
)

p <- DimPlot(cd4s,
             reduction = "Tertiary.harm_umap",
             group.by = "Tertiary.labels",
             label = T, repel = F, label.size = 5) +
  ggtitle("Cluster analysis of CD4+ T cells") +
  labs(x = "umap_1", y = "umap_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3))) &
  NoLegend()


ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD4/cd4s_umap.png",
       p,
       height = 5.5, width = 7)

# CD4 clusters - visualise gene expression ####

all_markers_pct10 <- FindAllMarkers(object = cd4s, logfc.threshold = 0.5, min.pct = 0.1, only.pos = TRUE)
all_markers_pct10 <- all_markers_pct10 %>%
  filter(p_val_adj < 0.05)
top_markers_pct10 <- all_markers_pct10 %>%
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n=50)

all_markers_pct50 <- FindAllMarkers(object = cd4s, logfc.threshold = 0.5, min.pct = 0.5, only.pos = TRUE)
top_markers_pct50 <- all_markers_pct50 %>%
  group_by(cluster) %>% 
  slice_max(avg_log2FC, n=50)
top_markers_pct50 <- all_markers_pct50 %>%
  group_by(cluster) %>%
  slice_max(avg_log2FC, n=50)

plot_markers <- top_markers_pct10 %>%
  pull(gene) %>%
  unique()

palette <- colorRampPalette(
  c("blue", 
    "black",
    "yellow"
  )
)(100)

heatmap <- DoHeatmap(cd4s,
                     features = plot_markers,
                     group.by = "Tertiary.harm_clusters") +
  scale_fill_gradientn(colors = palette) +
  theme(axis.text.y = element_blank(),
        axis.title.x = element_blank(),
        legend.position = "left")

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD4/clusters_heatmap.png",
       heatmap,
       height = 6, width = 6)

# Compare CD4 clonality by cluster ####

cd4s$clonalfreq_nt_cd4s <- table(cd4s$CTnt)[cd4s$CTnt]

cd4s$clonalFrequency_2 <- as.numeric(as.character(cd4s$clonalFrequency))

cd4s$clonalFrequency_2[is.na(cd4s$clonalFrequency_2)] <- 0

str(cd4s$clonalFrequency_2)

embeddings <- Embeddings(cd4s, reduction = "Tertiary.harm_umap")
df <- data.frame(
  x = embeddings[, 1],
  y = embeddings[, 2],
  clonalFrequency = cd4s$clonalFrequency_2,
  cluster = cd4s$Tertiary.labels
)

str(df$clonalFrequency)

p <- ggplot(df, aes(x = x, y = y)) +
  geom_point(aes(size = clonalFrequency, color = cluster), alpha = 0.7) +
  scale_size_continuous(range = c(0.1, 5),
                        breaks = c(80, 50, 20, 5, 1, 0),
                        labels = c("80", "50", "20", "5", "1", "N/A")) +
  labs(x = "umap_1", y = "umap_2", size = "Clonal frequency", color = "Cluster") +
  theme_classic() +
  guides(color = guide_legend(override.aes = list(size = 4)))
p

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD4/clonal/clonal_umap.png",
       p,
       height = 5, width = 7)

orig_colours <- scales::hue_pal()(length(unique(cd4s$Tertiary.labels)))

p <- clonalQuant(cd4s,
                 cloneCall = "nt",
                 chain = "both",
                 group.by = "Tertiary.labels",
                 scale = T) +
  labs(y = "Unique clones (%)", fill = "Region") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank(),
        legend.position = "none") +
  scale_fill_manual(values = orig_colours)

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD4/clonal/clonalquant_clusters.png",
       p,
       height = 3, width = 3)

p <- clonalProportion(cd4s,
                      cloneCall = "nt",
                      group.by = "Tertiary.labels",
                      clonalSplit = c(1, 10, 50, 100, 500, 1000)) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1),
        axis.title.x = element_blank()) +
  scale_fill_brewer(
    name = "Clonal Indices",
    palette = "Set2",
    labels = c("Largest clone", "Clones 2–10", "Clones 11–50", 
               "Clones 51–100", "Clones 101–500", "Clones 501–1000")
  )

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD4/clonal/clonalproportion_clusters.png",
       p,
       height = 3, width = 4.6)

cols <- c(orig_colours[0+1], # after inspecting result; 
          orig_colours[7+1], # each cluster number + 1 for position in orig_colours
          orig_colours[2+1],
          orig_colours[6+1],
          orig_colours[1+1],
          orig_colours[3+1],
          orig_colours[4+1],
          orig_colours[5+1])

metric_labels <- c(
  shannon = "Shannon Index",
  inv.simpson = "Inverse Simpson Index"
)

# Subset clusters with n > 100 for diversity

cluster_counts <- cd4s@meta.data %>%
  filter(!is.na(cloneSize)) %>%
  group_by(Tertiary.labels) %>%
  tally() %>%
  filter(n >= 100)

tcrs <- subset(cd4s, subset = Tertiary.labels %in% cluster_counts$Tertiary.labels)

cols <- c(orig_colours[0+1], # after inspecting result; 
          orig_colours[2+1], # each cluster number + 1 for position in orig_colours
          orig_colours[1+1],
          orig_colours[3+1],
          orig_colours[4+1],
          orig_colours[5+1])

set.seed(2025)
p <- clonalDiversity(tcrs,
                     cloneCall = "nt",
                     group.by = "Tertiary.labels",
                     metrics = c("shannon", "inv.simpson")) +
  scale_fill_manual(values = cols) +
  facet_wrap(~variable, scales = "free", ncol = length(c("shannon", "inv.simpson")), 
             labeller = labeller(variable = metric_labels)) +
  theme(legend.position = "none")
set.seed(2025)
ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD4/clonal/clonaldiversity_clusters.png",
       p,
       height = 3.5, width = 3.6)


p <- clonalOverlap(cd4s, cloneCall = "nt", chain = "both", method = "morisita", group.by = "Tertiary.labels") +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD4/clonal/clonaloverlap_clusters.png",
       p,
       height = 3.5, width = 5.5)

Idents(cd4s) <- cd4s$Tertiary.labels
p <- clonalCompare(cd4s,
                   samples = c("CD4-CTL", "Resting memory"),
                   top.clones = 10,
                   cloneCall="nt", 
                   graph = "alluvial",
                   relabel.clones = T) & NoLegend()
cols <- setNames(palette36.colors(18), levels(p$data$clones))
p <- p + scale_fill_manual(values = cols)
p


ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD4/clonal/compare_ctl_resting.png",
       p,
       height = 3, width = 4)

# Plot proportions of CD4 clusters by region ####

byRegion <- cd4s@meta.data %>%
  select(MULTI_Region, Tertiary.labels) %>%
  group_by(MULTI_Region, Tertiary.labels) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = Tertiary.labels, values_from = count, values_fill = 0) %>%
  mutate(All = `CD4-CTL` +
           `Resting memory` +
           `Early activated` +
           `Tcm` +
           `Th17` +
           `Treg` +
           `HSP-high` +
           `Cycling`,
         p.0 = `CD4-CTL` / All * 100,
         p.1 = `Resting memory` / All * 100,
         p.2 = `Early activated` / All * 100,
         p.3 = `Tcm` / All * 100,
         p.4 = `Th17` / All * 100,
         p.5 = `Treg` / All * 100,
         p.6 = `HSP-high` / All * 100,
         p.7 = `Cycling` / All * 100) %>%
  select(MULTI_Region, p.0, p.1, p.2, p.3, p.4, p.5, p.6, p.7) %>%
  pivot_longer(cols = c(p.0, p.1, p.2, p.3, p.4, p.5, p.6, p.7), names_to = "Type", values_to = "Proportion") %>%
  mutate(Type = factor(Type, levels = rev(c("p.0", "p.1", "p.2", "p.3", "p.4", "p.5", "p.6", "p.7"))),
         MULTI_Region = case_when(
           MULTI_Region == "Primary-Tumour" ~ "Primary Tumour",
           MULTI_Region == "Primary-PBZ" ~ "Primary PBZ",
           MULTI_Region == "Recurrence1-Tumour" ~ "Recurrent Tumour",
           MULTI_Region == "Recurrence1-PBZ" ~ "Recurrent PBZ",
           TRUE ~ MULTI_Region
         ),
         MULTI_Region = factor(MULTI_Region, levels = c("Primary Tumour", "Primary PBZ", "Recurrent Tumour", "Recurrent PBZ"))
  )

# Plot

cols <- rev(scales::hue_pal()(length(unique(cd4s$Tertiary.labels))))

regionPlot <- ggplot(
  byRegion, aes(x = MULTI_Region, y = Proportion, fill = Type)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = cols,
                    labels = rev(c("CD4-CTL", "Resting memory", "Early activated",
                                   "Tcm", "Th17", "Treg", "HSP-high", "Cycling"))) +
  labs(x = "Sample", y = "Proportions by region", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold", size = 12),
        strip.text = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none")

# Legend - note plot incorrect
cols <- scales::hue_pal()(length(unique(cd4s$Tertiary.labels)))
dummyRegion <- ggplot(
  byRegion, aes(x = MULTI_Region, y = Proportion, fill = Type)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = cols,
                    labels = c("CD4-CTL", "Resting memory", "Early activated",
                               "Tcm", "Th17", "Treg", "HSP-high", "Cycling")) +
  labs(x = "Sample", y = "Proportions by region", fill = "Cluster") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold"),
        strip.text = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.title = element_text(hjust = 0.5, face = "bold"))

dummyRegion
region_legend <- get_legend(dummyRegion)

grid.newpage()
grid.draw(region_legend)

output_dir <- "GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD4/"
ggsave(paste0(output_dir, "region_cd4_cluster_proportions.png"), regionPlot,
       height = 4, width = 1.5)
ggsave(paste0(output_dir, "region_legend.png"), region_legend,
       height = 4, width = 4)

# Plot proportions of CD4 clusters by donor ####

byDonor <- cd4s@meta.data %>%
  select(MULTI_Donor, Tertiary.labels) %>%
  group_by(MULTI_Donor, Tertiary.labels) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = Tertiary.labels, values_from = count, values_fill = 0) %>%
  mutate(All = `CD4-CTL` +
           `Resting memory` +
           `Early activated` +
           `Tcm` +
           `Th17` +
           `Treg` +
           `HSP-high` +
           `Cycling`,
         p.0 = `CD4-CTL` / All * 100,
         p.1 = `Resting memory` / All * 100,
         p.2 = `Early activated` / All * 100,
         p.3 = `Tcm` / All * 100,
         p.4 = `Th17` / All * 100,
         p.5 = `Treg` / All * 100,
         p.6 = `HSP-high` / All * 100,
         p.7 = `Cycling` / All * 100) %>%
  select(MULTI_Donor, p.0, p.1, p.2, p.3, p.4, p.5, p.6, p.7) %>%
  pivot_longer(cols = c(p.0, p.1, p.2, p.3, p.4, p.5, p.6, p.7), names_to = "Type", values_to = "Proportion") %>%
  mutate(Type = factor(Type, levels = rev(c("p.0", "p.1", "p.2", "p.3", "p.4", "p.5", "p.6", "p.7"))),
         MULTI_Donor = factor(MULTI_Donor, levels = c("N03", "N02", "N06", "N07", "N08", "N05"))) %>%
  filter(!is.na(MULTI_Donor))

# Plot

cols <- rev(scales::hue_pal()(length(unique(cd4s$Tertiary.labels))))

donorPlot <- ggplot(
  byDonor, aes(x = MULTI_Donor, y = Proportion, fill = Type)
) +
  geom_bar(stat = "identity", color = "black", linewidth = 0.3) +
  scale_fill_manual(values = cols,
                    labels = rev(c("CD4-CTL", "Resting memory", "Early activated",
                                   "Tcm", "Th17", "Treg", "HSP-high", "Cycling"))) +
  labs(x = "Sample", y = "Proportions by donor", fill = "Cell Type") +
  theme_minimal() +
  theme(axis.text.x = element_text(angle = 90, color = "black", vjust = 0.5, hjust = 0.5),
        axis.title.x = element_blank(),
        axis.text.y = element_blank(),
        axis.title.y = element_text(color = "black", face = "bold", size = 14),
        strip.text = element_blank(),
        plot.background = element_rect(fill = "white", color = NA),
        legend.position = "none")

ggsave(paste0(output_dir, "donor_cd4_cluster_proportions.png"), donorPlot,
       height = 4, width = 2.5)


# Cluster forest plot by region ####

primary <- cd4s@meta.data %>%
  filter(grepl("Primary", MULTI_Region))

table(primary$Tertiary.labels)
table(primary$MULTI_Region, primary$Tertiary.labels)

cd4data <- cd4s@meta.data %>%
  select(MULTI_ID, Tertiary.labels) %>%
  group_by(MULTI_ID, Tertiary.labels) %>%
  summarise(count = n()) %>%
  ungroup() %>%
  pivot_wider(names_from = Tertiary.labels, values_from = count, values_fill = 0) %>%
  mutate(
    All = `CD4-CTL` +
      `Resting memory` +
      `Early activated` +
      `Tcm` +
      `Th17` +
      `Treg` +
      `HSP-high` +
      `Cycling`,
    Donor = str_extract(MULTI_ID, "^.{3}"),
    Timepoint = ifelse(grepl("Primary", MULTI_ID), "Primary", "Recurrence1"),
    Tissue = ifelse(grepl("Tumour", MULTI_ID), "Tumour", "PBZ")) %>%
  select(MULTI_ID, Donor, Timepoint, Tissue, everything()) 

primary <- cd4data %>%
  filter(Timepoint == "Primary")

celltypes <- c("CD4-CTL",
                 "Resting memory",
                 "Early activated",
                 "Tcm",
                 "Th17",
                 "Treg",
                 "HSP-high",
                 "Cycling")

primary_clr_input <- primary[, celltypes] /100 # switch back to proportions
min_nonzero <- min(primary_clr_input[primary_clr_input > 0])
pseudocount <- min_nonzero / 100
primary_clr_adj <- primary_clr_input
primary_clr_adj[primary_clr_adj == 0] <- pseudocount
primary_clr_adj <- primary_clr_adj / rowSums(primary_clr_adj)

clr_data <- clr(primary_clr_adj)
clr_df <- cbind(primary[, c("Donor", "Tissue")], as.data.frame(clr_data))

sub_df <- data.frame(
  CellType = character(length(celltypes)),
  Estimate = numeric(length(celltypes)),
  SE = numeric(length(celltypes)),
  p_value = numeric(length(celltypes)),
  Upper = numeric(length(celltypes)),
  Lower = numeric(length(celltypes)),
  Donor = numeric(length(celltypes)),
  stringsAsFactors = FALSE
)

for (i in seq_along(celltypes)) {
  var <- celltypes[i]
  model_formula <- as.formula(paste0("`", var, "` ~ Tissue + (1 | Donor)"))
  model <- lmer(model_formula, data = clr_df)
  
  est <- fixef(model)["TissueTumour"]
  se <- summary(model)$coefficients["TissueTumour", "Std. Error"]
  donor <- attr(summary(model)$varcor$Donor, "stddev")
  
  sub_df[i, "CellType"] <- var
  sub_df[i, "Estimate"] <- est
  sub_df[i, "SE"] <- se
  sub_df[i, "p_value"] <- ifelse(summary(model)$coefficients["TissueTumour", "Pr(>|t|)"] < 0.0001, "p < 0.0001", signif(summary(model)$coefficients["TissueTumour", "Pr(>|t|)"], 2))
  sub_df[i, "Upper"] <- est + 1.96 * se
  sub_df[i, "Lower"] <- est - 1.96 * se
  sub_df[i, "Donor"] <- donor
}

sub_df$FoldChange <- ifelse(sub_df$Estimate < 0, -exp(abs(sub_df$Estimate)), exp(sub_df$Estimate))
sub_df$Upper_FC <- ifelse(sub_df$Upper < 0, -exp(abs(sub_df$Upper)), exp(sub_df$Upper))
sub_df$Lower_FC <- ifelse(sub_df$Lower < 0, -exp(abs(sub_df$Lower)), exp(sub_df$Lower))
sub_df$Donor_FC <- exp(sub_df$Donor)

sub_df <- sub_df %>%
  arrange(FoldChange) %>%
  mutate(CellType = factor(CellType, levels = unique(CellType)))

write.xlsx(sub_df, "GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD4/cd4s_region_forest_models.xlsx")

orig_colours <- scales::hue_pal()(length(unique(cd4s$Tertiary.labels)))
cols <- c(orig_colours[6+1], # after inspecting result; 
          orig_colours[1+1], # each cluster number + 1 for position in orig_colours
          orig_colours[5+1],
          orig_colours[7+1],
          orig_colours[0+1],
          orig_colours[3+1],
          orig_colours[4+1],
          orig_colours[2+1]
)

subDiff <- ggplot(sub_df, aes(x = CellType, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(color = rev(cols), size = 5) +
  geom_linerange(aes(ymin = Lower, ymax = Upper),
                 color = rev(cols), size = 1) +
  geom_text(
    aes(label = paste0("p = ", p_value), y = max(Upper) + 0.3), hjust = 0,
    color = "black",
    size = 5
  ) +
  coord_flip() +
  labs(y = "Centered log-ratio change", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.35),
    plot.background = element_rect(fill = "white")) +
  scale_y_continuous(limits = c(min(sub_df$Lower), max(sub_df$Upper) + 1),
                     breaks = seq(-4, 4, 1),
                     labels = seq(-4, 4, 1)
  )

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD4/cd4s_region_forest.png",
       subDiff,
       height = 3.5, width = 8.5)

# Cluster forest plot by donor ####

primary_clr_input <- cd4data[, celltypes] /100 # switch back to proportions
min_nonzero <- min(primary_clr_input[primary_clr_input > 0])
pseudocount <- min_nonzero / 100
primary_clr_adj <- primary_clr_input
primary_clr_adj[primary_clr_adj == 0] <- pseudocount
primary_clr_adj <- primary_clr_adj / rowSums(primary_clr_adj)

clr_data <- clr(primary_clr_adj)
clr_df <- cbind(cd4data[, c("Donor", "Tissue", "Timepoint")], as.data.frame(clr_data))

sub_df <- data.frame(
  CellType = character(length(celltypes)),
  Estimate = numeric(length(celltypes)),
  SE = numeric(length(celltypes)),
  p_value = character(length(celltypes)),
  Upper = numeric(length(celltypes)),
  Lower = numeric(length(celltypes)),
  Tissue_estimate = numeric(length(celltypes)),
  stringsAsFactors = FALSE
)

for (i in seq_along(celltypes)) {
  var <- celltypes[i]
  model_formula <- as.formula(paste0("`", var, "` ~ Timepoint + Tissue"))
  model <- lm(model_formula, data = clr_df)
  
  coef_summary <- summary(model)$coefficients
  
  est <- coef_summary["TimepointRecurrence1", "Estimate"]
  se <- coef_summary["TimepointRecurrence1", "Std. Error"]
  p_val <- coef_summary["TimepointRecurrence1", "Pr(>|t|)"]
  tissue <- coef_summary["TissueTumour", "Estimate"]
  
  sub_df[i, "CellType"] <- var
  sub_df[i, "Estimate"] <- est
  sub_df[i, "SE"] <- se
  sub_df[i, "p_value"] <- ifelse(p_val < 0.0001, "p < 0.0001", signif(p_val, 2))
  sub_df[i, "Upper"] <- est + 1.96 * se
  sub_df[i, "Lower"] <- est - 1.96 * se
  sub_df[i, "Tissue_estimate"] <- tissue
}

sub_df$FoldChange <- ifelse(sub_df$Estimate < 0, -exp(abs(sub_df$Estimate)), exp(sub_df$Estimate))
sub_df$Upper_FC <- ifelse(sub_df$Upper < 0, -exp(abs(sub_df$Upper)), exp(sub_df$Upper))
sub_df$Lower_FC <- ifelse(sub_df$Lower < 0, -exp(abs(sub_df$Lower)), exp(sub_df$Lower))
sub_df$Tissue_FC <- exp(sub_df$Tissue_estimate)

sub_df <- sub_df %>%
  arrange(FoldChange) %>%
  mutate(CellType = factor(CellType, levels = unique(CellType)))

write.xlsx(sub_df, "GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD4/cd4s_timepoint_forest_models.xlsx")

orig_colours <- scales::hue_pal()(length(unique(cd8s$Tertiary.labels)))
cols <- c(orig_colours[0+1], # after inspecting result; 
          orig_colours[5+1], # each cluster number + 1 for position in orig_colours
          orig_colours[2+1],
          orig_colours[4+1],
          orig_colours[1+1],
          orig_colours[7+1],
          orig_colours[3+1],
          orig_colours[6+1]
)

subDiff <- ggplot(sub_df, aes(x = CellType, y = Estimate)) +
  geom_hline(yintercept = 0, linetype = "dashed", size = 1) +
  geom_point(color = rev(cols), size = 5) +
  coord_flip() +
  labs(y = "Centered log-ratio change", x = NULL) +
  theme_minimal() +
  theme(
    axis.text.y = element_blank(),
    axis.text.x = element_text(size = 15),
    axis.title.x = element_text(size = 19, hjust = 0.4),
    plot.background = element_rect(fill = "white")) +
  scale_y_continuous(limits = c(-3, 4),
                     breaks = seq(-3, 4, 1),
                     labels = seq(-3, 4, 1)
  )

ggsave("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD4/cd4s_timepoint_forest.png",
       subDiff,
       height = 3.5, width = 8.5)

# DGEA across regions by cluster ####

primary <- subset(cd4s,
                  grepl("Primary", MULTI_Region))
Idents(primary) <- "Tertiary.labels"
deseq2_markers <- list()
wilcoxon_markers <- list()
clusters <- levels(primary$Tertiary.labels)

for (cl in clusters) {
  cells <- WhichCells(primary, idents = cl)
  cluster <- subset(primary, cells = cells)
  
  cluster$MULTI_Region <- factor(cluster$MULTI_Region, levels = c("Primary-Tumour", "Primary-PBZ"))
  Idents(cluster) <- "MULTI_Region"
  
  # Pseudobulk by MULTI_Donor within region
  counts <- GetAssayData(cluster, slot = "counts", assay = "RNA")
  meta <- cluster@meta.data
  
  # Aggregate counts per donor-region pair
  meta$group_id <- paste(meta$MULTI_Donor, meta$MULTI_Region, sep = "_")
  
  pseudobulk <- as.data.frame(as.matrix(counts)) %>%
    t() %>%
    as.data.frame() %>%
    mutate(group_id = meta$group_id) %>%
    group_by(group_id) %>%
    summarise(across(everything(), sum)) %>%
    column_to_rownames("group_id") %>%
    t()
  
  # Create colData
  group_split <- strsplit(colnames(pseudobulk), "_")
  coldata <- data.frame(
    donor = sapply(group_split, `[`, 1),
    region = sapply(group_split, `[`, 2)
  )
  rownames(coldata) <- colnames(pseudobulk)
  coldata$region <- factor(coldata$region, levels = c("Primary-Tumour", "Primary-PBZ"))
  
  # Only proceed if both conditions are present and ≥10 cells
  region_cell_counts <- table(cluster$MULTI_Region)
  if (all(c("Primary-Tumour", "Primary-PBZ") %in% names(region_cell_counts)) && all(region_cell_counts >= 10)) {
    
    dds <- DESeqDataSetFromMatrix(countData = pseudobulk,
                                  colData = coldata,
                                  design = ~ region)
    dds <- DESeq(dds)
    res <- results(dds, contrast = c("region", "Primary-Tumour", "Primary-PBZ"))
    
    res_df <- as.data.frame(res)
    res_df$gene <- rownames(res_df)
    res_df$cluster <- cl
    
    deseq2_markers[[cl]] <- res_df
    
    # Wilcoxon to compare
    
    markers <- FindMarkers(cluster,
                           ident.1 = "Primary-Tumour",
                           ident.2 = "Primary-PBZ",
                           logfc.threshold = 0.5,
                           min.pct = 0.1,
                           only.pos = FALSE)
    
    markers_filtered <- markers %>%
      mutate(gene = rownames(.)) %>%
      slice_max(order_by = avg_log2FC, n = 20, with_ties = FALSE) %>%
      bind_rows(
        markers %>%
          mutate(gene = rownames(.)) %>%
          slice_min(order_by = avg_log2FC, n = 20, with_ties = FALSE)
      ) %>%
      mutate(cluster = cl)
    
    wilcoxon_markers[[cl]] <- markers_filtered
  }
}

# Combine all results
deseq2_markers_df <- bind_rows(deseq2_markers) %>%
  group_by(cluster) %>%
  slice_max(order_by = log2FoldChange, n = 20) %>%
  bind_rows(
    bind_rows(deseq2_markers) %>%
      group_by(cluster) %>%
      slice_min(order_by = log2FoldChange, n = 20)
  )


wilcoxon_markers_df <- bind_rows(wilcoxon_markers) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 20) %>%
  bind_rows(
    bind_rows(wilcoxon_markers) %>%
      group_by(cluster) %>%
      slice_min(order_by = avg_log2FC, n = 20)
  )

deseq2_bind <- deseq2_markers_df %>%
  mutate(log2FC = log2FoldChange,
         method = "DESeq2") %>%
  select(cluster, gene, log2FC, method) %>%
  arrange(cluster, desc(log2FC))

wilcoxon_bind <- wilcoxon_markers_df %>%
  mutate(log2FC = avg_log2FC,
         method = "Wilcoxon") %>%
  select(cluster, gene, log2FC, method) %>%
  arrange(cluster, desc(log2FC))

degs <- bind_rows(deseq2_bind, wilcoxon_bind) %>%
  group_by(cluster, gene) %>%
  mutate(max_log2FC = max(log2FC, na.rm = T)) %>%
  ungroup() %>%
  arrange(cluster, desc(max_log2FC), gene)

degs <- degs %>%
  arrange(cluster, desc(max_log2FC)) %>%
  group_by(cluster, gene) %>%
  filter(n() == 2,
         !(n() == 2 & method == "Wilcoxon")) %>%
  ungroup() %>%
  mutate(order = case_when(
    cluster == "CD4-CTL" ~ 1,
    cluster == "Early activated" ~ 2,
    cluster == "Tcm" ~ 3,
    cluster == "Th17" ~ 4,
    cluster == "Treg" ~ 5
  )) %>%
  select(order, cluster, gene, log2FC) %>%
  arrange(order, desc(log2FC))

any(duplicated(degs$gene))

mat <- pivot_wider(degs, names_from = cluster, values_from = log2FC)

mat <- mat %>%
  column_to_rownames("gene") %>%
  select(-order) %>%
  as.matrix()

palette <- colorRampPalette(
  c("blue",    
    "#E6E0EB",
    "firebrick3" 
  )
)(100)

orig_colours <- scales::hue_pal()(length(unique(cd8s$Tertiary.harm_clusters)))

annotation_col <- data.frame(
  Cluster = colnames(mat)
)
rownames(annotation_col) <- colnames(mat)

# Create a named vector of colors for your clusters
cluster_colors <- setNames(orig_colours[c(1,2,3,5,7)], unique(annotation_col$Cluster))

hm <- pheatmap(mat,
               na_col = "grey90",
               color = palette,
               cluster_rows = FALSE,
               cluster_cols = FALSE,
               annotation_col = annotation_col,
               annotation_legend = F,
               annotation_colors = list(Cluster = cluster_colors))

pdf("GBM_single_cell_analysis/outputs/clustering/Lymphoid/CD4/cd4_region_degs.pdf",
    width = 3, height = 9)
draw(hm)
dev.off()

# Myeloid analysis ####
# Remove unused clusters/reductions ####

Mseu@meta.data <- Mseu@meta.data[, !grepl("^harm_subset_clusters_res", colnames(Mseu@meta.data))]
Idents(Mseu) <- Mseu$harm_subset_clusters_pc20_res1
Mseu$seurat_clusters <- Mseu$harm_subset_clusters_pc20_res1

Mseu$Secondary.harm_clusters <- Mseu$harm_subset_clusters_pc20_res1
Mseu@meta.data <- Mseu@meta.data %>%
  select(-starts_with("harm_subset_clusters"))

Mseu@reductions$Secondary.harm_umap <- Mseu@reductions$umap.harm_subset_pc20
to_remove <- grep("^umap\\.harm_subset", names(Mseu@reductions), value = TRUE)
Mseu@reductions[to_remove] <- NULL

Mseu@graphs$Secondary.harm_snn <- Mseu@graphs$harmony_subset_snn_pc25
to_remove <- grep("^harmony_subset", names(Mseu@graphs), value = TRUE)
Mseu@graphs[to_remove] <- NULL

# Overwrite elbow plot

elbow <- ElbowPlot(Mseu, reduction = "pca", ndims = 50) +
  ggtitle("Percentage variance by principal component (Myeloid)") +
  theme_classic() +
  theme(
    plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12,face = "bold")
  )   +
  geom_vline(xintercept = 25.5, color = "red", linetype = "dashed", linewidth = 0.5)

output_dir <- paste0("GBM_single_cell_analysis/outputs/clustering/Myeloid/")
dir.create(output_dir, recursive=T, showWarnings=F)

ggsave(paste0(output_dir, "Myeloid_elbow_plot.png"),
       elbow,
       height = 6, width = 6)

# Repeat automatic annotation ####

Mseu$GBMap_3_Myeloid <- NA
Mseu$GBMap_4_Myeloid <- NA
Mseu$GBMap_MyeloidModules <- NA

GBMap_MyeloidModules <- read.xlsx("GBM_single_cell_analysis/outputs/clustering/GBMap_MyeloidModules.xlsx", colNames=F)
colnames(GBMap_MyeloidModules) <- c("cellName", "marker")

GBMap_MyeloidModules <- GBMap_MyeloidModules %>%
  group_by(cellName) %>%
  summarise(geneSymbolmore1 = paste(marker, collapse = ","), .groups = "drop")

GBMap_MyeloidModules$tissueType <- "Glioblastoma"
GBMap_MyeloidModules$geneSymbolmore2 <- NA

GBMap_MyeloidModules <- GBMap_MyeloidModules %>%
  select(tissueType, cellName, geneSymbolmore1, geneSymbolmore2)

write.xlsx(GBMap_MyeloidModules, "GBM_single_cell_analysis/outputs/clustering/GBMap_MyeloidModulesDB.xlsx")

source("https://raw.githubusercontent.com/IanevskiAleksandr/sc-type/master/R/sctype_wrapper.R"); 
Mseu <- run_sctype(Mseu, assay = "RNA", scaled = TRUE, known_tissue_type = "Glioblastoma", custom_marker_file = "GBM_single_cell_analysis/outputs/clustering/GBMap_3Myeloid.xlsx", name = "GBMap_3_Myeloid")
Mseu <- run_sctype(Mseu, assay = "RNA", scaled = TRUE, known_tissue_type = "Glioblastoma", custom_marker_file = "GBM_single_cell_analysis/outputs/clustering/GBMap_4Myeloid.xlsx", name = "GBMap_4_Myeloid")
Mseu <- run_sctype(Mseu, assay = "RNA", scaled = TRUE, known_tissue_type = "Glioblastoma", custom_marker_file = "GBM_single_cell_analysis/outputs/clustering/GBMap_MyeloidModulesDB.xlsx", name = "GBMap_MyeloidModules")


# UMAP ####

Mseu$Secondary.harm_clusters <- factor(
  Mseu$Secondary.harm_clusters,
  levels = sort(as.numeric(as.character(unique(Mseu$Secondary.harm_clusters))))
)

# Tidy up three tiny clusters 14, 15, 16

Mseu$Secondary.harm_clusters[Mseu$Secondary.harm_clusters == "15"] <- "2"
Mseu$Secondary.harm_clusters[Mseu$Secondary.harm_clusters == "16"] <- "4"
Mseu$Secondary.harm_clusters[Mseu$Secondary.harm_clusters == "14"] <- "3"


p <- DimPlot(Mseu,
             reduction = "Secondary.harm_umap",
             group.by = "MULTI_Region",
             label = F, label.size = 4) +
  ggtitle("Myeloid cells grouped by tumour region") +
  labs(x = "umap_1", y = "umap_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))

ggsave("GBM_single_cell_analysis/outputs/clustering/Myeloid/myeloid_region.png",
       p,
       height = 6.2, width = 9)


# DGEA ####


markers <- FindAllMarkers(Mseu, 
                            logfc.threshold = 0.5, 
                            min.pct = 0.1, 
                            only.pos = TRUE)
  
clusterMarkers <- markers %>%
    mutate(
      diff_pct = (pct.1 - pct.2) * 100,
      signif = -log10(p_val_adj),
      signif = ifelse(is.infinite(signif) | signif > 300, 300, signif)) %>%
    group_by(cluster) %>%
    arrange(cluster, desc(diff_pct)) %>%
    slice_head(n = 20) %>%
    ungroup()
  
m.data <- Mseu@meta.data %>%
    arrange(seurat_clusters)
  
clusters <- unique(m.data$seurat_clusters)
cluster_plotlist <- list()
  
cat("Plotting clusters...", "\n")
  
for (i in seq_along(clusters)){
    cluster_id <- clusters[[i]]
    cluster_cells <- length(which(Mseu$seurat_clusters == cluster_id))
    cluster_size <- round(cluster_cells / nrow(Mseu@meta.data) * 100, 2)
    
    topmarkers <- clusterMarkers %>%
      filter(cluster == cluster_id) %>%
      arrange(desc(avg_log2FC))
    
    p <- ggplot(topmarkers, aes(x = avg_log2FC, y = reorder(gene, avg_log2FC), size = signif)) +
      # Trailing lines from y-axis to dots
      geom_segment(aes(x = 0, xend = avg_log2FC, y = reorder(gene, avg_log2FC), yend = reorder(gene, avg_log2FC)),
                   color = "gray", size = 0.5, alpha = 0.5) +  # Line settings
      
      # Dot plot
      geom_point(aes(color = diff_pct, x = avg_log2FC, y = gene), alpha = 1) +  
      labs(title = paste0(cluster_id, "; n = ", cluster_cells, "; ", cluster_size, "%")) +
      scale_color_gradient(
        low = "#132132", 
        high = "steelblue1", 
        limits = range(clusterMarkers$diff_pct, na.rm = TRUE),
        name = "Fold change"
      ) +
      scale_size_continuous(name = "Significance",
                            range = c(0.05, 2),
                            limits = range(clusterMarkers$signif, na.rm = TRUE)) +  # Adjust dot size range
      xlim(0, 10) +
      theme_cowplot() +
      theme(axis.text.y = element_text(size = 4.5),
            axis.text.x = element_text(size = 5.5),
            axis.title.x = element_blank(),
            axis.title.y = element_blank(),
            plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
            legend.title = element_text(size = 5, face = "bold"),
            legend.text = element_text(size = 5),
            legend.position = "none",
            plot.background = element_rect(fill = "white", color = NA))
    
    cluster_plotlist[[i]] <- p
    
  }
  
legend_plot <- ggplot(topmarkers, aes(x = avg_log2FC, y = reorder(gene, avg_log2FC), size = signif)) +
    # Trailing lines from y-axis to dots
    geom_segment(aes(x = 0, xend = avg_log2FC, y = reorder(gene, avg_log2FC), yend = reorder(gene, avg_log2FC)),
                 color = "gray", size = 0.5, alpha = 0.5) +  # Line settings
    
    # Dot plot
    geom_point(aes(color = diff_pct, x = avg_log2FC, y = gene), alpha = 1) +  
    labs(title = paste0(cluster_id, "; n = ", cluster_cells, "; ", cluster_size, "%")) +
    scale_color_gradient(
      low = "#132132", 
      high = "steelblue1", 
      limits = range(clusterMarkers$diff_pct, na.rm = TRUE),
      name = "Differential percentage expression"
    ) +
    scale_size_continuous(name = "Significance",
                          range = c(0.05, 2),
                          limits = range(clusterMarkers$signif, na.rm = TRUE)) +  # Adjust dot size range
    xlim(0, 10) +
    theme_cowplot() +
    theme(axis.text.y = element_text(size = 4.5),
          axis.text.x = element_text(size = 5.5),
          axis.title.x = element_blank(),
          axis.title.y = element_blank(),
          plot.title = element_text(size = 10, face = "bold", hjust = 0.5),
          legend.title = element_text(size = 5, face = "bold"),
          legend.text = element_text(size = 5),
          legend.position = "right",
          plot.background = element_rect(fill = "white", color = NA))
  
legend <- get_legend(legend_plot)
  
cluster_grid <- plot_grid(plotlist = cluster_plotlist, ncol = 5)
plot <- plot_grid(cluster_grid, legend,
                    rel_widths = c(5,1), nrow = 1, align = "vh")
final_plot <- ggdraw(plot) +
    draw_label("log2FC", x = 0.5, y = 0.01, vjust = 0.5, fontface = "bold", size = 14) +
    draw_label("Gene", x = 0.01, y = 0.5, angle = 90, vjust = 0.5, fontface = "bold", size = 14)
  
  
  
pdf(paste0("GBM_single_cell_analysis/outputs/clustering/Myeloid/myeloid_logFC_degs.pdf"),
      height = 10, width = 12)
  print(final_plot)
  dev.off()



#Exploratory analysis ####

c13 <- subset(Mseu,
              subset = Secondary.harm_clusters == "13")

genes <- c()
cd4cd8 <- lapply(genes, function(gene){
  df <- FetchData(Mseu, vars = c("Secondary.harm_clusters", gene)) %>%
    mutate(gene = gene,
           expressed = .data[[gene]] > 0) %>%
    group_by(Secondary.harm_clusters, gene) %>%
    summarise(proportion = mean(expressed)) %>%
    ungroup()
}) %>% bind_rows() %>%
  arrange(Secondary.harm_clusters)

FeaturePlot(Mseu, reduction = "Secondary.harm_umap",
            features = c("HLA-DQA1", "MS4A6A", "TGFBI", "HLA-DRB5", "CD74"))

DimPlot(Mseu, reduction = "Secondary.harm_umap", group.by = "Jackson.labels")



features <- FeaturePlot(Mseu, reduction = "Secondary.harm_umap",
                        features = c("FCER1A", "S100A4", "VCAN", "S100A6", "FN1", "LGALS3", "THBS1", "CD55", "MT2A", "LGALS1", "P2RY12", "CCL4L2", "CD69", "CD68", "CD163", "HLA-DRB1", "HLA-DQA1", "C1QA", "IFI6", "HSPH1")) +
  labs(x = "umap_1", y = "umap_2") +
  guides(color = guide_legend(ncol = 1, override.aes = list(size = 3)))
features <- lapply(features, function(p) {
  p + theme(
    axis.title = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank()
  )
})
features <- wrap_plots(features)

ggsave("GBM_single_cell_analysis/outputs/clustering/Myeloid/features.png",
       features,
       height = 6, width = 12)

JacksonDegs <- read.xlsx("GBM_single_cell_analysis/outputs/clustering/Jackson_DEGs.xlsx", colNames=F)
colnames(JacksonDegs) <- c("cellName", "posmarker", "negmarker")

JacksonDegs <- JacksonDegs %>%
  group_by(cellName) %>%
  summarise(geneSymbolmore1 = paste(na.omit(posmarker), collapse = ","),
            geneSymbolmore2 = paste(rev(na.omit(negmarker)), collapse = ","), .groups = "drop") %>%
  mutate(tissueType = "Glioblastoma") %>%
  select(tissueType, cellName, geneSymbolmore1, geneSymbolmore2)

write.xlsx(JacksonDegs, "GBM_single_cell_analysis/outputs/clustering/Jackson_sigs.xlsx")

rm(JacksonDegs)

Mseu <- run_sctype(Mseu, assay = "RNA", scaled = TRUE, known_tissue_type = "Glioblastoma", custom_marker_file = "GBM_single_cell_analysis/outputs/clustering/Jackson_sigs.xlsx", name = "Jackson.labels")




