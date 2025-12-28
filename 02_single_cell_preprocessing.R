# Load libraries ####

library(Seurat)
library(ggplot2)
library(tidyverse)
library(dplyr)
library(purrr)
library(cowplot)
library(patchwork)
library(reshape2)
library(sctransform)
library(scRepertoire)
library(scCustomize)
library(openxlsx)

set.seed(2025)

# Read in data ####

# Find files and create sample information dataframe
samples <- list.files("GBM_single_cell_analysis/data/GEX_CITE")

sampleInfo <- data.frame(Name=samples)
sampleInfo$Batch <- c("GBM01", "GBM01", "GBM02", "GBM02")
sampleInfo$Cell_types <- c("CD45+ CD3-", "CD45+ CD3+", "CD45+ CD3-", "CD45+ CD3+")
sampleInfo$Cells_sorted <- c(62830, 28225, 66541, 6622)

rm(samples)

objectList <- list()
for (sample in seq_along(sampleInfo)){
  
  df <- Read10X(
    file.path("GBM_single_cell_analysis/data/GEX_CITE/",
              sampleInfo$Name[sample],
              "filtered_feature_barcode_matrix"))
  
  objectList[[sample]] <- CreateSeuratObject(
    counts = df[["Gene Expression"]],
    project = sampleInfo$Name[sample],
    min.cells = 3,
    min.features = 0)
  
  objectList[[sample]][["percent.mt"]] <- PercentageFeatureSet(
    objectList[[sample]],
    pattern = "^MT-")
  
  rm(df)
  } # read in files and create Seurat objects

for (sample in seq_along(sampleInfo)){
  
  sampleInfo$Before_filtering[sample] <- dim(objectList[[sample]])[2]
  
} # get cell numbers before filtering

# QC plots ####

nCount_PlotList <- list()
for (sample in seq_len(nrow(sampleInfo))){
  
  p <- VlnPlot(
    objectList[[sample]],
    features = "nCount_RNA",
    pt.size = 0
    ) +
    ggtitle(sampleInfo$Name[[sample]]) +
    xlab("Number of cells") +
    ylab("Number of UMI counts") +
    ylim(0, 43000) +
    theme(legend.position = "none") +
    geom_hline(yintercept = 500,
               colour="red",
               linewidth = 1,
               linetype = "dotted") +
    geom_hline(yintercept = 7000,
               colour="red",
               linewidth = 1,
               linetype = "dotted")
  
  nCount_PlotList[[sample]] <- p
} # Plot UMI counts

nFeature_PlotList <- list()
for (sample in seq_len(nrow(sampleInfo))){
  
  p <- VlnPlot(
    objectList[[sample]],
    features = "nFeature_RNA",
    pt.size = 0
    ) +
    ggtitle(sampleInfo$Name[[sample]]) +
    xlab("Number of cells") +
    ylab("Number of genes") +
    ylim(0, 8000) +
    theme(legend.position = "none") +
    geom_hline(yintercept = 500,
               colour="red",
               linewidth = 1,
               linetype = "dotted") +
    geom_hline(yintercept = 2800,
               colour="red",
               linewidth = 1,
               linetype = "dotted")
  
  nFeature_PlotList[[sample]] <- p
} # Plot nFeatures

percent.mt_PlotList <- list()
for (sample in seq_len(nrow(sampleInfo))){
  
  p <- VlnPlot(
    objectList[[sample]],
    features = "percent.mt",
    pt.size = 0
    ) +
    ggtitle(sampleInfo$Name[[sample]]) +
    xlab("Number of cells") +
    ylab("Mitochondrial gene counts (% total)") +
    ylim(0, 100) +
    theme(legend.position = "none") +
    geom_hline(yintercept = 10,
               colour="red",
               linewidth = 1,
               linetype = "dotted")
  
  percent.mt_PlotList[[sample]] <- p
} # Plot percent.mt

# Create QC violin plot grid
Violin_Grid <- plot_grid(
  plot_grid(plotlist = nCount_PlotList,
            ncol = 4),
  plot_grid(plotlist = nFeature_PlotList,
            ncol = 4),
  plot_grid(plotlist = percent.mt_PlotList,
            ncol = 4),
  ncol = 1
)

ggsave(
  filename = "GBM_single_cell_analysis/outputs/QC_Violin_Grid.png",
  plot = Violin_Grid,
  width = 10,
  height = 12
  )

rm(p, nCount_PlotList, nFeature_PlotList, percent.mt_PlotList, Violin_Grid)

nCount_nFeature_scatter_PlotList <- list()
for (sample in seq_len(nrow(sampleInfo))){
  p <- FeatureScatter(
    objectList[[sample]],
    feature1 = "nCount_RNA",
    feature2 = "nFeature_RNA"
    ) +
    ggtitle(sampleInfo$Name[[sample]]) +
    ylab("Number of genes") +
    xlab("Number of UMI counts") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 22,
                                    hjust = 0.5,
                                    face = "bold"),
      axis.title = element_text(size = 15),
      legend.position = "none")
  
  nCount_nFeature_scatter_PlotList[[sample]] <- p
  
} # Plot count/feature scatter
nCount_nFeature_scatter_Grid <- plot_grid(
  plotlist = nCount_nFeature_scatter_PlotList,
  ncol = 4
  )

nCount_mt_scatter_PlotList <- list()

for (sample in seq_len(nrow(sampleInfo))){
  
  p <- FeatureScatter(
    objectList[[sample]],
    feature1 = "nCount_RNA",
    feature2 = "percent.mt"
    ) +
    ggtitle(sampleInfo$Name[[sample]]) +
    ylab("Mitochondrial gene counts (% total)") +
    xlab("Number of UMI counts") +
    theme_classic() +
    theme(
      plot.title = element_text(size = 22,
                                hjust = 0.5,
                                face = "bold"),
      axis.title = element_text(size = 15),
      legend.position = "none"
      )
  
  nCount_mt_scatter_PlotList[[sample]] <- p
  
} # Plot count/mt scatter
nCount_mt_scatter_Grid <- plot_grid(
  plotlist = nCount_mt_scatter_PlotList,
  ncol = 4
  )

# Create QC scatter plot grid
Scatter_Grid <- plot_grid(
  nCount_nFeature_scatter_Grid, nCount_mt_scatter_Grid,
  ncol = 1)

ggsave(
  filename = "GBM_single_cell_analysis/outputs/QC_Scatter_Grid.png",
  plot = Scatter_Grid,
  width = 12, height = 8
  )
rm(p, nCount_nFeature_scatter_PlotList, nCount_mt_scatter_PlotList, nCount_nFeature_scatter_Grid, nCount_mt_scatter_Grid, Scatter_Grid)

# Filter and batch process data ####

objectsFiltered <- list()
for (sample in seq_len(nrow(sampleInfo))){
  
  objectsFiltered[[sample]] <- subset(
    objectList[[sample]],
    subset =
      nCount_RNA > 500 &
      nCount_RNA < 10000 &
      nFeature_RNA > 500 &
      nFeature_RNA < 3500 &
      percent.mt < 10
    )
} # filter Seurat objects, triplicate to test regression
objects_nCount <- objectsFiltered
objects_mt <- objectsFiltered
rm(objectList)

for (sample in seq_len(nrow(sampleInfo))){
  sampleInfo$After_filtering[sample] <- dim(objectsFiltered[[sample]])[2]
} # get cell numbers after filtering

BatchProcess <- function(
    seu,
    regress = NULL,
    n.pcs = 30,
    resolution = 0.5,
    umap.dims = 20
    ) {  
  
  seu %>%
    NormalizeData() %>%
    FindVariableFeatures(verbose = FALSE) %>%
    ScaleData(vars.to.regress = regress) %>%
    RunPCA(features = VariableFeatures(seu)) %>%
    FindNeighbors(dims = 1:n.pcs) %>%
    FindClusters(resolution = resolution) %>%
    RunUMAP(dims = 1:umap.dims)
} # wrap batch process functions

objectsFiltered <- lapply(objectsFiltered, BatchProcess)

objects_nCount <- lapply(objects_nCount, BatchProcess, regress = "nCount_RNA")

objects_mt <- lapply(objects_mt, BatchProcess, regress = "percent.mt")

Elbow_PlotList <- list()
for (object in seq_along(objectsFiltered)){
  p <- ElbowPlot(objectsFiltered[[object]], ndims = 50) +
    ggtitle(sampleInfo$Name[[object]]) +
    theme_classic() +
    theme(
      plot.title = element_text(size = 22,
                                hjust = 0.5,
                                face = "bold"),
      axis.title = element_text(size = 15))
  
  Elbow_PlotList[[object]] <- p
} # Batch elbow plots

Elbow_Grid <- plot_grid(plotlist = Elbow_PlotList, ncol = 2)
ggsave("GBM_single_cell_analysis/outputs/QC_Elbow_Grid.png",
       Elbow_Grid,
       height = 10, width = 10)
rm(p, Elbow_PlotList, Elbow_Grid)

ref.se <- celldex::BlueprintEncodeData()
AnnotateSingleR <- function(seu, ref.se){
  pred_seu <- SingleR(
    test = as.SingleCellExperiment(seu),
    ref = ref.se,
    assay.type.test = 1,
    labels = ref.se$label.fine
  )
  
  prediction_table <- data.frame(
    table(
      pred_seu$labels))
  colnames(prediction_table) <- c("CellType", "Number")
  
  seu[["SingleR.labels"]] <- pred_seu$labels
  rm(pred_seu, prediction_table)
  return(seu)
} # wrap SingleR function

all_objects <- list(objectsFiltered, objects_nCount, objects_mt)
for (list in seq_along(all_objects)){
  
  for (object in seq_along(all_objects[[list]])){
    
    seu <- all_objects[[list]][[object]]
    seu <- AnnotateSingleR(seu, ref.se = ref.se)
    all_objects[[list]][[object]] <- seu
  }
} # nested loop to annotate all objects

objectsFiltered <- all_objects[[1]]
objects_nCount <- all_objects[[2]]
objects_mt <- all_objects[[3]]
rm(seu, ref.se, all_objects)

# Look for unwanted variance - nCount ####

nCount_unregressed_scatter_PlotList <- vector("list", length(objectsFiltered))
for (object in seq_along(objectsFiltered)) {
  
  nCount_unregressed_scatter_PlotList[[object]] <- vector("list", 5)  # Ensure sublist exists
  
  for (pc in seq_len(5)){
    
    seu <- objectsFiltered[[object]]
    Idents(seu) <- seu$seurat_clusters
    
    p <- FeatureScatter(
      seu,
      feature1 = paste0("PC_", pc),
      feature2 = "nCount_RNA"
      ) +
      ggtitle(paste0("PC_",
                     pc,
                     " before regression | R = ",
                     round(cor(seu$nCount_RNA,
                               seu@reductions$pca@cell.embeddings[, pc],
                               method = "pearson"), 2))) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(y = "Number of UMI counts") +
      NoLegend()
    
    nCount_unregressed_scatter_PlotList[[object]][[pc]] <- p
    objectsFiltered[[object]] <- seu
  }
} # plot first 5 PCs without nCount regression

nCount_regressed_scatter_PlotList <- vector("list", length(objects_nCount))
for (object in seq_along(objects_nCount)){
  
  nCount_regressed_scatter_PlotList[[object]] <- vector("list", 5)
  
  for (pc in seq_len(5)){
    
    seu <- objects_nCount[[object]]
    Idents(seu) <- seu$seurat_clusters
    
    p <- FeatureScatter(
      seu,
      feature1 = paste0("PC_", pc),
      feature2 = "nCount_RNA") +
      ggtitle(paste0("PC_",
                     pc, 
                     " after regression | R = ",
                    round(cor(seu$nCount_RNA,
                              seu@reductions$pca@cell.embeddings[, pc],
                              method = "pearson"), 2))) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(y = "Number of UMI counts") +
      NoLegend()
    
    nCount_regressed_scatter_PlotList[[object]][[pc]] <- p
    objects_nCount[[object]] <- seu
  }
} # Plot with nCount regression

nCount_regression_PCs_Grid <- list()
for (plot in seq_along(nCount_unregressed_scatter_PlotList)) {
  
  title_plot <- ggdraw() + 
    draw_label(paste0(sampleInfo$Name[plot], " - effect of UMI count regression by top 5 principal components"), fontface = "bold", size = 24, hjust = 0.5) +
    theme(plot.background = element_rect(fill = "white", color = NA))
  
  combined_plot <- plot_grid(
    plot_grid(plotlist = nCount_unregressed_scatter_PlotList[[plot]],
              ncol = 1),
    plot_grid(plotlist = nCount_regressed_scatter_PlotList[[plot]],
              ncol = 1),
    ncol = 2
  )
  
  nCount_regression_PCs_Grid[[plot]] <- plot_grid(
    title_plot, combined_plot,
    ncol = 1,
    rel_heights = c(0.1, 1)
    ) 
  
  ggsave(
    filename = paste0("GBM_single_cell_analysis/outputs/nCount_regression_PCs_Grid_", 
                      sampleInfo$Name[plot],
                      ".png"),
    plot = nCount_regression_PCs_Grid[[plot]],  
    height = 16, 
    width = 12
  )
} # Plot and save nCount regression PC grids
rm(p, seu, title_plot, combined_plot, nCount_unregressed_scatter_PlotList, nCount_regressed_scatter_PlotList, nCount_regression_PCs_Grid)

nCountByCluster_unregressed_PlotList <- list()
for (object in seq_along(objectsFiltered)){
  
  seu <- objectsFiltered[[object]]
  Idents(seu) <- seu$seurat_clusters
  
  p <- VlnPlot(seu, "nCount_RNA", pt.size = 0) +
    ggtitle(paste0(sampleInfo$Name[[object]],
                   " before regression")) +
    xlab("Cluster") +
    ylab("UMI Count") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    NoLegend()
  
  nCountByCluster_unregressed_PlotList[[object]] <- p
  
  objectsFiltered[[object]] <- seu
} # Plot UMIs by cluster before regression

nCountByCluster_regressed_PlotList <- list()
for (object in seq_along(objects_nCount)){
  
  seu <- objects_nCount[[object]]
  Idents(seu) <- seu$seurat_clusters
  
  p <- VlnPlot(seu, "nCount_RNA", pt.size = 0) +
    ggtitle(paste0(sampleInfo$Name[[object]], " after regression")) +
    xlab("Cluster") +
    ylab("UMI Count") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5, face = "bold")) +
    NoLegend()
  
  nCountByCluster_regressed_PlotList[[object]] <- p
  
  objects_nCount[[object]] <- seu
} # After regression

# Plot and save nCount by cluster grid:

title_plot <- ggdraw() + 
  draw_label("Effect of UMI count regression on clustering",
             fontface = "bold",
             size = 24,
             hjust = 0.5) +
  theme(plot.background = element_rect(fill = "white",
                                       color = NA)
        )
  
combined_plot <- plot_grid(
  plot_grid(plotlist = nCountByCluster_unregressed_PlotList,
            ncol = 1),
  plot_grid(plotlist = nCountByCluster_regressed_PlotList,
            ncol = 1),
  ncol = 2
  )
  
nCountByCluster_Grid <- plot_grid(
  title_plot, combined_plot,
  ncol = 1,
  rel_heights = c(0.1, 1)
  )
  
ggsave(filename = "GBM_single_cell_analysis/outputs/nCountByCluster_Grid.png",
       plot = nCountByCluster_Grid,
       height = 16, width = 12)

rm(p, seu, title_plot, combined_plot, nCountByCluster_unregressed_PlotList, nCountByCluster_regressed_PlotList, nCountByCluster_Grid)

df_objectsFiltered <- list()
for (object in seq_along(objectsFiltered)){
  
  seu <- objectsFiltered[[object]]
  Idents(seu) <- seu$SingleR.labels
  prop_seu <- prop.table(table(seu$seurat_clusters,
                               seu$SingleR.labels), margin = 1)
  df_seu <- as.data.frame(prop_seu)
  colnames(df_seu) <- c("Cluster", "Cell_Type", "Proportion")
  df_seu$Condition <- "Before regression"
  df_objectsFiltered[[object]] <- df_seu
  objectsFiltered[[object]] <- seu
  
  rm(seu, prop_seu, df_seu)
} # Create df for plotting celltypes by cluster before nCount regression

df_objects_nCount <- list()
for (object in seq_along(objects_nCount)){
  
  seu <- objects_nCount[[object]]
  Idents(seu) <- seu$SingleR.labels
  prop_seu <- prop.table(table(seu$seurat_clusters,
                               seu$SingleR.labels), margin = 1)
  df_seu <- as.data.frame(prop_seu)
  colnames(df_seu) <- c("Cluster", "Cell_Type", "Proportion")
  df_seu$Condition <- "After regression"
  df_objects_nCount[[object]] <- df_seu
  objects_nCount[[object]] <- seu
  
  rm(seu, prop_seu, df_seu)
} # After regression

df_nCount <- list()
for (object in seq_along(df_objectsFiltered)) {
  
  df_nCount[[object]] <- rbind(df_objectsFiltered[[object]], df_objects_nCount[[object]])
  df_nCount[[object]]$Condition <- factor(df_nCount[[object]]$Condition, levels = c("Before regression", "After regression"))
} # Combine for plotting

for (data in seq_along(df_nCount)){
  p <- ggplot(
    df_nCount[[data]],
    aes(x = Cluster,
        y = Proportion,
        fill = Cell_Type)
    ) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~Condition, scales = "free_x") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white",
                                      color = NA),  # white background
      plot.background = element_rect(fill = "white",
                                     color = NA),  # white background for the plot
      text = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      strip.text = element_text(size = 16)
    ) +
    labs(title = paste0(sampleInfo$Name[data], " - effect of UMI count regression on cell type proportions"),
         x = "Cluster",
         y = "Proportion",
         fill = "Cell Type")
  
  ggsave(
    paste0("GBM_single_cell_analysis/outputs/CellTypes_nCountRegression_",
                sampleInfo$Name[data],
                ".png"),
    p,
    height = 12, width = 12
    )
} # Plot and save cell types by cluster before and after nCount regression
rm(p, df_objects_nCount, df_nCount, objects_nCount)

# Look for unwanted variance - mt ####

mt_unregressed_scatter_PlotList <- vector("list", length(objectsFiltered))
for (object in seq_along(objectsFiltered)){
  
  mt_unregressed_scatter_PlotList[[object]] <- vector("list", 5)
  
  for (pc in seq_len(5)){
    seu <- objectsFiltered[[object]]
    Idents(seu) <- seu$seurat_clusters
    p <- FeatureScatter(
      seu,
      feature1 = paste0("PC_", pc),
      feature2 = "percent.mt"
      ) +
      ggtitle(paste0("PC_", pc, 
                    " before regression | R = ",
                    round(cor(seu$percent.mt,
                              seu@reductions$pca@cell.embeddings[, pc],
                              method = "pearson"), 2))) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(y = "Mitochondrial gene counts (% total)") +
      NoLegend()
    
    mt_unregressed_scatter_PlotList[[object]][[pc]] <- p
    objectsFiltered[[object]] <- seu
  }
} # plot first 5 PCs without mt regression

mt_regressed_scatter_PlotList <- vector("list", length(objects_mt))
for (object in seq_along(objects_mt)){
  
  mt_regressed_scatter_PlotList[[object]] <- vector("list", 5)
  for (pc in seq_len(5)){
    
    seu <- objects_mt[[object]]
    Idents(seu) <- seu$seurat_clusters
    p <- FeatureScatter(
      seu,
      feature1 = paste0("PC_", pc),
      feature2 = "percent.mt"
      ) +
      ggtitle(paste0("PC_", pc,
                    " after regression | R = ",
                    round(cor(seu$percent.mt,
                              seu@reductions$pca@cell.embeddings[, pc],
                              method = "pearson"), 2))) +
      theme_classic() +
      theme(plot.title = element_text(hjust = 0.5)) +
      labs(y = "Mitochondrial gene counts (% total)") +
      NoLegend()
    
    mt_regressed_scatter_PlotList[[object]][[pc]] <- p
    objects_mt[[object]] <- seu
  }
} # plot with mt regression

mt_regression_PCs_Grid <- list()
for (plot in seq_along(mt_unregressed_scatter_PlotList)) {
  
  title_plot <- ggdraw() +
    draw_label(
      paste0(sampleInfo$Name[plot],
             " - effect of mitochondrial regression by top 5 principal components"),
      fontface = "bold",
      size = 24,
      hjust = 0.5
      ) +
    theme(plot.background = element_rect(fill = "white",
                                         color = NA))
  
  combined_plot <- plot_grid(
    plot_grid(plotlist = mt_unregressed_scatter_PlotList[[plot]],
              ncol = 1),
    plot_grid(plotlist = mt_regressed_scatter_PlotList[[plot]],
              ncol = 1),
    ncol = 2
  )
  
  mt_regression_PCs_Grid[[plot]] <- plot_grid(
    title_plot, combined_plot,
    ncol = 1,
    rel_heights = c(0.1, 1)
    ) 
  
  ggsave(
    filename = paste0("GBM_single_cell_analysis/outputs/mt_regression_PCs_Grid_", 
                      sampleInfo$Name[plot],
                      ".png"),
    plot = mt_regression_PCs_Grid[[plot]],  
    height = 16, 
    width = 12
  )
} # plot and save mt regression grids
rm(p, seu, mt_unregressed_scatter_PlotList, mt_regressed_scatter_PlotList, title_plot, combined_plot, mt_regression_PCs_Grid)

mtByCluster_unregressed_PlotList <- list()
for (object in seq_along(objectsFiltered)){
  
  seu <- objectsFiltered[[object]]
  Idents(seu) <- seu$seurat_clusters
  p <- VlnPlot(seu, "percent.mt", pt.size = 0) +
    ggtitle(paste0(sampleInfo$Name[[object]], " before regression")) +
    xlab("Cluster") +
    ylab("Mitochondrial gene counts (% total)") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,
                                    face = "bold")) +
    NoLegend()
  
  mtByCluster_unregressed_PlotList[[object]] <- p
  objectsFiltered[[object]] <- seu
} # Plot mt by cluster before regression

mtByCluster_regressed_PlotList <- list()
for (object in seq_along(objects_mt)){
  
  seu <- objects_mt[[object]]
  Idents(seu) <- seu$seurat_clusters
  p <- VlnPlot(seu, "percent.mt", pt.size = 0) +
    ggtitle(paste0(sampleInfo$Name[[object]], " after regression")) +
    xlab("Cluster") +
    ylab("Mitochondrial gene counts (% total)") +
    theme_classic() +
    theme(plot.title = element_text(hjust = 0.5,
                                    face = "bold")) +
    NoLegend()
  
  mtByCluster_regressed_PlotList[[object]] <- p
  objects_mt[[object]] <- seu
} # After regression

title_plot <- ggdraw() +
  draw_label(
    "Effect of mitochondrial regression on clustering",
    fontface = "bold",
    size = 24,
    hjust = 0.5
    ) +
  theme(plot.background = element_rect(fill = "white",
                                       color = NA))

combined_plot <- plot_grid(
  plot_grid(plotlist = mtByCluster_unregressed_PlotList,
            ncol = 1),
  plot_grid(plotlist = mtByCluster_regressed_PlotList,
            ncol = 1),
  ncol = 2  
)

mtByCluster_Grid <- plot_grid(
  title_plot,
  combined_plot,
  ncol = 1,
  rel_heights = c(0.1, 1)
  )
  
ggsave("GBM_single_cell_analysis/outputs/mtByCluster_Grid.png",
       plot = mtByCluster_Grid,
       height = 16, width = 12) # plot and save mt by cluster

rm(p, seu, mtByCluster_unregressed_PlotList, mtByCluster_regressed_PlotList, title_plot, combined_plot, mtByCluster_Grid)

df_objects_mt <- list()
for (object in seq_along(objects_mt)){
  
  seu <- objects_mt[[object]]
  Idents(seu) <- seu$SingleR.labels
  prop_seu <- prop.table(table(seu$seurat_clusters,
                               seu$SingleR.labels), margin = 1)
  df_seu <- as.data.frame(prop_seu)
  colnames(df_seu) <- c("Cluster", "Cell_Type", "Proportion")
  df_seu$Condition <- "After regression"
  df_objects_mt[[object]] <- df_seu
  objects_mt[[object]] <- seu
  rm(seu, prop_seu, df_seu)
} # Create df for plotting celltypes after mt regression

df_mt <- list()
for (object in seq_along(df_objectsFiltered)) {
  
  df_mt[[object]] <- rbind(df_objectsFiltered[[object]], df_objects_mt[[object]])
  df_mt[[object]]$Condition <- factor(df_mt[[object]]$Condition,
                                      levels = c("Before regression", "After regression"))
} # Combine for plotting

for (data in seq_along(df_mt)){
  p <- ggplot(
    df_mt[[data]],
    aes(x = Cluster,
        y = Proportion,
        fill = Cell_Type)
    ) +
    geom_bar(stat = "identity", position = "fill") +
    facet_wrap(~Condition, scales = "free_x") +
    theme_minimal() +
    theme(
      panel.background = element_rect(fill = "white",
                                      color = NA),  # white background
      plot.background = element_rect(fill = "white",
                                     color = NA),  # white background for the plot
      text = element_text(size = 14),
      plot.title = element_text(hjust = 0.5, size = 18, face = "bold"),
      strip.text = element_text(size = 16)
    ) +
    labs(title = paste0(sampleInfo$Name[data],
                        " - effect of mitochondrial regression on cell type proportions"),
         x = "Cluster",
         y = "Proportion",
         fill = "Cell Type")
  
  ggsave(paste0("GBM_single_cell_analysis/outputs/CellTypes_mtRegression_", 
                sampleInfo$Name[data],
                ".png"),
         p,
         height = 12, width = 12)
} # Plot and save cell types by cluster before and after mt regression
rm(p, df_objectsFiltered, df_objects_mt, df_mt, objects_mt)

# Add in CITE data ####

for (sample in seq_len(nrow(sampleInfo))){
  
  seu <- objectsFiltered[[sample]]
  df <- Read10X(file.path("GBM_single_cell_analysis/data/GEX_CITE",
                          sampleInfo$Name[sample],
                          "filtered_feature_barcode_matrix"))
  CITE_df <- df[["Antibody Capture"]]
  CITE_df <- CITE_df[, Cells(seu)]
  
  if (sample == 1 | sample == 2) {
    rownames(CITE_df) <- c("N01-Primary-PBZ1", "N02-Primary-Tumour1", "N03-Primary-Tumour1", "N05-Recurrence1_Tumour1", "N06-Primary-Tumour1", "N07-Primary-Tumour1", "N08-Primary-Tumour1", "CD4-hashtag", "CD8-hashtag")
  } else if (sample == 3 | sample == 4) {
    rownames(CITE_df) <- c("N01-Primary-PBZ2", "N02-Primary-PBZ1", "N03-Primary-PBZ1", "N05-Recurrence1-PBZ1", "N06-Primary-PBZ1", "N07-Primary-PBZ2", "N08-Primary-Tumour4", "CD45-int-hashtag", "CD4-hashtag", "CD8-hashtag")
  }
  
  seu[["CITE"]] <- CreateAssayObject(CITE_df)
  seu <- NormalizeData(seu, assay = "CITE", normalization.method = "CLR")
  objectsFiltered[[sample]] <- seu
  
  rm(seu, df, CITE_df)
} # Conditional loop to name hashtags by histology

for (object in seq_along(objectsFiltered)) {
  
  seu <- objectsFiltered[[object]]
  CITE_data <- seu[["CITE"]]@counts
  donor_idx <- grep("N0", rownames(CITE_data))
  donorData <- CITE_data[donor_idx, ]
  CITE_data <- CITE_data[-donor_idx, ]
  seu[["Donor"]] <- CreateAssayObject(donorData)
  seu <- NormalizeData(seu, assay = "Donor", normalization.method = "CLR")
  seu[["CITE"]] <- CreateAssayObject(CITE_data)
  seu <- NormalizeData(seu, assay = "CITE", normalization.method = "CLR")
  objectsFiltered[[object]] <- seu
  
  rm(seu, CITE_data, donor_idx, donorData)
} # Move donor hashtags to Donor assay (CITE hashtags remain in CITE assay)

hashtagCelltype_PlotList <- list()
index <- 1
for (sample in which(sampleInfo$Cell_types == "CD45+ CD3+")) {

  seu <- objectsFiltered[[sample]]
  
  Idents(seu) <- seu$SingleR.labels
  
  p <- RidgePlot(
    seu,
    features = "CD4-hashtag",
    sort = TRUE,
    ncol = 2
    ) +
    ggtitle(paste0(sampleInfo$Name[sample], " - CD4 hashtag by cell type")) +
    theme(
      plot.title = element_text(size = 18,
                                hjust = 0.5,
                                face = "bold"),
      axis.title.x = element_text(hjust = 0.5,
                                  face = "bold"),
      axis.title.y = element_text(hjust = 0.5,
                                  face = "bold")) +
    NoLegend()
  
  hashtagCelltype_PlotList[[index]] <- p
  index <- index + 1
  
  p <- RidgePlot(
    seu,
    features = "CD8-hashtag",
    ncol = 2
    ) +
    ggtitle(paste0(sampleInfo$Name[sample], " - CD8 hashtag by cell type")) +
    theme(
      plot.title = element_text(size = 18,
                                hjust = 0.5,
                                face = "bold"),
      axis.title.x = element_text(hjust = 0.5,
                                  face = "bold"),
      axis.title.y = element_text(hjust = 0.5, 
                                  face = "bold")) +
    NoLegend()
  
  hashtagCelltype_PlotList[[index]] <- p
  index <- index + 1
  objectsFiltered[[sample]] <- seu
  
  rm(p, seu)
} # Plot and save CD4/CD8 hashtags

hashtagCelltype_Grid <- plot_grid(
  plotlist = hashtagCelltype_PlotList,
  ncol = 2
)

ggsave("GBM_single_cell_analysis/outputs/hashtagCelltype_Grid.png",
       hashtagCelltype_Grid,
       height = 14, width = 18
       )
rm(index, hashtagCelltype_PlotList, hashtagCelltype_Grid)

# Demultiplexing ####

# Settings from LJI Github:

objectsDemux <- list()
for (object in seq_along(objectsFiltered)){
  
  seu <- objectsFiltered[[object]]
    
  seu <- MULTIseqDemux(
    seu,
    assay = "Donor",
    autoThresh = TRUE,
    maxiter = 10
    )
  
  tvar <- as.character(seu$MULTI_ID) %in% c('Negative', 'Doublet')
  
  seu$MULTI_classification.global <- ifelse(tvar, as.character(seu$MULTI_ID), 'Singlet')
  
  objectsDemux[[object]] <- seu
  
  rm(seu, tvar)
}

# Fine-tune donor ID calls as per LJI Github ####

for (object in seq_along(objectsDemux)){
  
  print(paste("Fine-tuning demultiplexing for sample", sampleInfo$Name[[object]]))
  
  seu <- objectsDemux[[object]]
  
  annot <- seu@meta.data[, grepl("MULTI|HTO|hash|origlib", colnames(seu[[]]))]
  
  ids_umis <- data.frame(data.table::rbindlist(lapply(1:nrow(annot), function(thiscell){
    assigned_id <- as.character(annot[thiscell, "MULTI_ID"])
    all_ids <- seu@assays$Donor@counts[, rownames(annot[thiscell, ])]
    
    # Exclude the assigned ID to get the rest
    rest_ids <- all_ids[!names(all_ids) %in% assigned_id]
    
    assigned_umi <- all_ids[assigned_id] # take hashtag ID (NA if Negative or Doublet)
    
    y <- data.frame(
      hash.ID = names(which.max(all_ids)), # this is the highest regardless of the assigned
      # first ID's UMI, take assigned; if not assigned, taking the first in rest
      firstUMI = unname(ifelse(is.na(assigned_umi), max(rest_ids), assigned_umi)),
      # second ID, take the top in rest; if not assigned, taking the second in rest
      secondID = unname(ifelse(is.na(assigned_umi), names(sort(rest_ids))[length(rest_ids) - 1], names(which.max(rest_ids)))),
      secondUMI = unname(ifelse(is.na(assigned_umi), sort(rest_ids)[length(rest_ids) - 1], max(rest_ids)))
    ); colnames(y) <- paste0("MULTI_ID", "_", colnames(y))
    y
  })), row.names = rownames(annot))

  annot <- cbind(annot,ids_umis)
  
  cnames <- grep("MULTI_ID", colnames(annot), value = TRUE)
  
  annot$HT_FoldChange <- annot[, grep("firstUMI$",cnames,v=T)] / annot[, grep("secondUMI$",cnames,v=T)]
  
  annot$HT_FoldChange[is.infinite(annot$HT_FoldChange)] <- max(annot$HT_FoldChange[is.finite(annot$HT_FoldChange)])
  
  annot$HT_FoldChange[is.nan(annot$HT_FoldChange)] <- min(annot$HT_FoldChange, na.rm = TRUE)
  
  # Assigned id does not have the highest UMI count
  annot$HT_FoldChange_conflict <- annot[, grep("firstUMI$",cnames,v=T)] < annot[, grep("secondUMI$",cnames,v=T)]
  
  # Swapping the donor identity of conflict cells that are not doublets
  tvar <- annot$HT_FoldChange_conflict & annot$MULTI_ID != "Doublet"
  if(any(tvar)){ 
    cat("No. of false highest assigned:", sum(tvar), "\n")
    print(table(annot[tvar, ]$MULTI_ID))
    
    # Correct MULTI_ID and MULTI_classification for the conflict cells
    annot[tvar, ]$MULTI_ID <- annot[tvar, ]$MULTI_ID_secondID
    annot[tvar, ]$MULTI_classification <- annot[tvar, ]$MULTI_ID_secondID
    
    # Adjust HT_FoldChange for swapped cells - because the assigned ID was the second-highest by UMI
    annot[tvar, ]$HT_FoldChange <- annot[tvar, grep("secondUMI$",cnames, v=T)] / annot[tvar, grep("firstUMI$",cnames, v=T)]
  }
  
  # Assign singlets with a fold change < 3 as doublets
  tvar <- annot$HT_FoldChange < 3 & annot$MULTI_ID != "Negative"
  if(any(tvar)){
    cat("Low-confidence singlets re-assigned:", sum(tvar), "\n")
    annot$MULTI_ID <- ifelse(annot$HT_FoldChange < 3 & annot$MULTI_ID != "Negative", "Doublet", as.character(annot$MULTI_ID))
  }
  
  # Rescue any negative calls that have good fold change
  tvar <- annot$HT_FoldChange >= 3 & annot$MULTI_classification.global == "Negative"
  if(any(tvar)){
    cat("Negative declassification:", sum(tvar), "\n")
    annot$MULTI_ID <- ifelse(tvar, annot$MULTI_ID_hash.ID, as.character(annot$MULTI_ID))
    annot$MULTI_classification <- ifelse(tvar, annot$MULTI_ID_hash.ID, as.character(annot$MULTI_classification))
  }

  # Rescue any doublet calls that have good fold change
  tvar <- annot$HT_FoldChange >= 3 & annot$MULTI_classification.global == "Doublet"
  if(any(tvar)){
    cat("Doublet declassification:", sum(tvar), "\n")
    annot$MULTI_ID <- ifelse(tvar, annot$MULTI_ID_hash.ID, as.character(annot$MULTI_ID))
    annot$MULTI_classification <- ifelse(tvar, annot$MULTI_ID_hash.ID, as.character(annot$MULTI_classification))
  }
  
  # Re-label MULTI_classification.global
  tvar <- as.character(annot$MULTI_ID) %in% c("Negative", "Doublet") # mark the doublets and negatives
  annot$MULTI_classification.global <- ifelse(tvar, as.character(annot$MULTI_ID), 'Singlet') # newly classified re-named as singlets
  
  #Distribution of the fold changes
  print(summary(annot$HT_FoldChange))
  
  #summary table of the donor 
  print(paste("Sample", sampleInfo$Name[[object]], "donor summary table:"))
        
  print(table(annot[, "MULTI_ID"]))
  
  # Assign back to seu
  seu$MULTI_classification.global <- annot$MULTI_classification.global
  seu$MULTI_classification <- annot$MULTI_classification
  seu$HT_FoldChange <- annot$HT_FoldChange
  seu$MULTI_ID <- annot$MULTI_ID
  
  # Capture numbers
  after_doublet <-
    length(seu$MULTI_classification.global) -
    sum(seu$MULTI_classification.global == "Doublet")
  
  after_demux <-
    after_doublet -
    sum(seu$MULTI_classification.global == "Negative")
  
  sampleInfo[object, "After_doublet_removal"] <- after_doublet
  sampleInfo[object, "After_demultiplexing"] <- after_demux
  
  objectsDemux[[object]] <- seu
}

rm(annot, ids_umis, tvar, cnames, after_doublet, after_demux, objectsFiltered)

# Check final doublet calls in the GEX

Doublet_Vln_PlotList <- list()

index <- 1

for (i in seq_along(objectsDemux)){
  seu <- objectsDemux[[i]]
  p <- VlnPlot(
    seu,
    features = "nFeature_RNA",
    group.by = "MULTI_classification.global",
    pt.size = 0) +
    labs(title = paste0(sampleInfo$Name[[i]]) , y = "Number of genes") +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      legend.position = "none")
  
  Doublet_Vln_PlotList[[index]] <- p
  
  index <- index + 1
  
  p <- VlnPlot(
    seu,
    features = "nCount_RNA",
    group.by = "MULTI_classification.global",
    pt.size = 0) +
    labs(title = paste0(sampleInfo$Name[[i]]) , y = "Number of UMI counts") +
    theme(
      axis.title.x = element_blank(),
      axis.title.y = element_text(face = "bold"),
      legend.position = "none")
  
  Doublet_Vln_PlotList[[index]] <- p
  
  index <- index + 1
}

title <- ggdraw() +
  draw_label("Final doublet, negative and singlet calls by gene and UMI count",
             fontface = "bold",
             size = 20, hjust = 0.5) +
  theme(plot.background = element_rect(fill = "white", color = NA))

Doublet_Vln_Grid <- plot_grid(
  title,
  plot_grid(plotlist = Doublet_Vln_PlotList,
            ncol = 2),
  ncol = 1,
  rel_heights = c(0.1, 1)
)

ggsave("GBM_single_cell_analysis/outputs/Doublet_Vln_Grid.png",
       Doublet_Vln_Grid,
       height = 14, width = 10)

rm(p, seu, Doublet_Vln_Grid, Doublet_Vln_PlotList, title)

# Filter demultiplexed singlets, capture cell numbers and divide into myeloid and T cell objects ####

recovered <- list()
for (object in seq_along(objectsDemux)){
  seu <- objectsDemux[[object]]
  seu <- subset(
    seu,
    subset = MULTI_classification.global == "Singlet")
  seu$MULTI_classification <- NULL
  objectsDemux[[object]] <- seu
  recovered[[object]] <- as.data.frame(table(seu$orig.ident, seu$MULTI_ID))
}

recovered <- bind_rows(recovered)

write.xlsx(recovered, "GBM_single_cell_analysis/outputs/recovered_numbers.xlsx")

objectsTCell <- objectsDemux[which(sampleInfo$Cell_types == "CD45+ CD3+")]
objectsMyeloid <- objectsDemux[which(sampleInfo$Cell_types == "CD45+ CD3-")]

rm(objectsDemux, seu)

# Relationship between cells loaded and cells applied to donor

hashed_numbers <- data.frame(read.xlsx("GBM_single_cell_analysis/outputs/hashed_numbers.xlsx"))

plot <- ggplot(hashed_numbers, aes(x = sorted, y = recovered, color = batch)) + 
  geom_smooth(method = "lm", color = "black", aes(group = 1), se = TRUE, alpha = 0.2) +
  geom_point() + 
  scale_x_log10(n.breaks = 6) +
  scale_y_log10(n.breaks = 5) +
  theme_cowplot() + 
  ggtitle("Cell recovery by sample") +
  labs(y = "Number of cells recovered", x = "Number of cells sorted") +
  theme(
    plot.title = element_text(size = 16, hjust = 0.5, face = "bold"),
    axis.title = element_text(size = 12, face = "bold"),
    legend.title = element_text(size = 12, face = "bold"),
    plot.background = element_rect(fill = "white", color = NA),
    panel.background = element_rect(fill = "white", color = NA) # Ensure background is white
  ) +
  scale_color_discrete(name = "Batch")

ggsave("GBM_single_cell_analysis/outputs/QC/cell_recovery_by_sample.png",
       plot = plot, 
       width = 5, height = 4.5, dpi = 300, units = "in", bg = "white")

rm(recovered, hashed_numbers)

# Add in TCR data ####

TCR_samplenames <- sampleInfo$Name[which(sampleInfo$Cell_types == "CD45+ CD3+")]
sampleInfo$TCR_recovered <- NA

for (object in seq_along(objectsTCell)){
  seu <- objectsTCell[[object]]
  tcrData <- read.csv(
    paste0("GBM_single_cell_analysis/data/TCR/",
           TCR_samplenames[object],
           ".csv"))
  combined.TCR <- combineTCR(tcrData, 
                             samples = TCR_samplenames[object],
                             removeNA = FALSE, 
                             removeMulti = FALSE, 
                             filterMulti = FALSE)
  
  seu <- RenameCells(seu, add.cell.id = TCR_samplenames[object])
  
  seu <- combineExpression(combined.TCR, seu, cloneCall = "nt", chain = "both")
  
  num_TCR <- sum(!is.na(seu$CTgene))
  index <- which(sampleInfo$Name == as.character(unique(seu$orig.ident)))
  sampleInfo$TCR_recovered[index] <- num_TCR
  objectsTCell[[object]] <- seu
  
}

rm(num_TCR, index, tcrData, combined.TCR)


TCR_QC_PlotList <- list()

for(object in seq_along(objectsTCell)){
  
  seu <- objectsTCell[[object]]
  seu$clonalFrequency_scaled <- log1p(seu$clonalFrequency)

  p1 <- FeaturePlot(seu, features = "TRBC1", pt.size = 0.7) +
    labs(title = paste0(TCR_samplenames[[object]], " - TRBC1 gene expression")) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  
  p2 <- FeaturePlot(seu, features = "clonalFrequency_scaled", pt.size = 0.7) +
    labs(title = paste0(TCR_samplenames[[object]], " - clonal frequency")) +
    theme(
      plot.title = element_text(size = 20, face = "bold", hjust = 0.5)
    )
  
  p3 <- plot_grid(
    p1, p2,
    ncol = 2
  )
  
  TCR_QC_PlotList[[object]] <- p3
  
  seu$clonalFrequency_scaled <- NULL
  rm(p1, p2, p3)

}

TCR_QC_plot <- plot_grid(
  plotlist = TCR_QC_PlotList,
  ncol = 1
)
  
ggsave("GBM_single_cell_analysis/outputs/QC/TCR_QC_plot2.png",
       TCR_QC_plot,
       height = 10, width = 12)

rm(seu, TCR_QC_PlotList, TCR_QC_plot)

objectsMyeloid_ready_to_integrate <- list()
objectsTCell_ready_to_integrate <- list()

for (object in seq_along(objectsMyeloid)){
  
  myeloid_object <- objectsMyeloid[[object]] # Relies on equal TCell and Myeloid batches
  rownames(myeloid_object@meta.data) <- make.unique(rownames(myeloid_object@meta.data))
  
  extra_columns <- c("CTgene", "CTnt", "CTaa", "CTstrict", "clonalProportion", "clonalFrequency", "cloneSize")
  for (col in extra_columns){
    myeloid_object@meta.data[[col]] <- NA # make columns suitable for merging
  }
  
  tcell_object <- objectsTCell[[object]]
  rownames(tcell_object@meta.data) <- make.unique(rownames(tcell_object@meta.data))
  
  objectsMyeloid_ready_to_integrate[[object]] <- myeloid_object
  objectsTCell_ready_to_integrate[[object]] <- tcell_object
}

rm(objectsTCell, objectsMyeloid, tcell_object,
   myeloid_object, extra_columns)

# Integrate and merge ####

gc()
options(future.globals.maxSize = 6.5 * 1024^3)

# global - integ and merge

global_integ <- merge(x = objectsTCell_ready_to_integrate[[1]], y = c(objectsTCell_ready_to_integrate[[2]],
                                                                     objectsMyeloid_ready_to_integrate[[1]],
                                                                     objectsMyeloid_ready_to_integrate[[2]]))

global_integ <- BatchProcess(global_integ)
global_integ <- IntegrateLayers(
  object = global_integ,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = FALSE)
global_integ <- IntegrateLayers(
  object = global_integ,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE)
global_integ <- IntegrateLayers(
  object = global_integ,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE)

global_merge <- merge(x = objectsTCell_ready_to_integrate[[1]], y = c(objectsTCell_ready_to_integrate[[2]],
                                                                      objectsMyeloid_ready_to_integrate[[1]],
                                                                      objectsMyeloid_ready_to_integrate[[2]]))
global_merge <- JoinLayers(global_merge)
global_merge <- BatchProcess(global_merge, resolution = 2, umap.dims = 30)
  
# T cells - integ and merge

tcell_integ <- Reduce(merge, objectsTCell_ready_to_integrate)
tcell_integ <- BatchProcess(tcell_integ)
tcell_integ <- IntegrateLayers(
  object = tcell_integ,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE)
tcell_integ <- IntegrateLayers(
  method = RPCAIntegration,
  orig.reduction = "pca",
  object = tcell_integ,
  new.reduction = "integrated.rpca",
  verbose = FALSE)
tcell_integ <- IntegrateLayers(
  object = tcell_integ,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE)

tcell_integ <- JoinLayers(tcell_integ)

tcell_merge <- Reduce(merge, objectsTCell_ready_to_integrate)
tcell_merge <- JoinLayers(tcell_merge)
tcell_merge <- BatchProcess(tcell_merge, resolution = 2, umap.dims = 30)

# Myeloid cells - integ and merge

myeloid_integ <- Reduce(merge, objectsMyeloid_ready_to_integrate)
myeloid_integ <- BatchProcess(myeloid_integ)
myeloid_integ <- IntegrateLayers(
  object = myeloid_integ,
  method = CCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.cca",
  verbose = FALSE)
myeloid_integ <- IntegrateLayers(
  object = myeloid_integ,
  method = RPCAIntegration,
  orig.reduction = "pca",
  new.reduction = "integrated.rpca",
  verbose = FALSE)
myeloid_integ <- IntegrateLayers(
  object = myeloid_integ,
  method = HarmonyIntegration,
  orig.reduction = "pca",
  new.reduction = "harmony",
  verbose = FALSE)

myeloid_integ <- JoinLayers(myeloid_integ)

myeloid_merge <- Reduce(merge, objectsMyeloid_ready_to_integrate)
myeloid_merge <- JoinLayers(myeloid_merge)
myeloid_merge <- BatchProcess(myeloid_merge, resolution = 2, umap.dims = 30)

rm(objectsTCell_ready_to_integrate, objectsMyeloid_ready_to_integrate)

# Tidy up metadata:

wrangle <- list(global_integ)

for (object in seq_along(wrangle)){
  m.data <- wrangle[[object]]@meta.data
  m.data <- m.data %>%
    select(
      -c("MULTI_classification.global")
    ) %>%
    mutate(
      MULTI_Donor = str_extract(MULTI_ID, "^[^-]+"),
      MULTI_Region = str_match(MULTI_ID, "^[^-]+-([^-]+-[A-Za-z]+)")[,2]
    ) %>%
    select(
      1:which(names(m.data) == "MULTI_ID"),
      "MULTI_ID",
      "MULTI_Donor",
      "MULTI_Region",
      "HT_FoldChange",
      everything()
    )
  wrangle[[object]]@meta.data <- m.data
}

thesis_seu <- wrangle[[1]]

rm(wrangle, m.data)

####### NB if reinstating other integrations, adjust cluster and reduction names
 global_integ - all batches integrated x3 (use reductions)

global_integ <- FindNeighbors(global_integ, reduction = "integrated.cca", dims = 1:30)
global_integ <- FindClusters(global_integ, resolution = 2, cluster.name = "cca_clusters")
global_integ <- RunUMAP(global_integ, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
umap_global_cca <- DimPlot(global_integ, reduction = "umap.cca", group.by = "orig.ident") +
  ggtitle("CCA") +
  xlab("umap_1") +
  ylab("umap_2") &
  NoLegend()

global_integ <- FindNeighbors(global_integ, reduction = "integrated.rpca", dims = 1:30)
global_integ <- FindClusters(global_integ, resolution = 2, cluster.name = "rpca_clusters")
global_integ <- RunUMAP(global_integ, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
umap_global_rpca <- DimPlot(global_integ, reduction = "umap.rpca", group.by = "orig.ident") +
  ggtitle("RPCA") +
  xlab("umap_1") +
  ylab("umap_2") &
  NoLegend()

thesis_seu <- FindNeighbors(thesis_seu, reduction = "harmony", dims = 1:25)
thesis_seu <- FindClusters(thesis_seu, resolution = 0.5, graph.name = "harmony_snn_pc30", cluster.name = "harm_clusters_pc30_res0.5")
thesis_seu <- RunUMAP(thesis_seu, reduction = "harmony", dims = 1:25, reduction.name = "umap.harm_pc25")
umap_global_harmony <- DimPlot(thesis_seu, reduction = "umap.harm_pc30", group.by = "harm_clusters_pc30_res0.5", label = T) +
  ggtitle("Harmony") +
  xlab("umap_1") +
  ylab("umap_2") &
  NoLegend()

# global_merge

umap_global_merge <- DimPlot(global_merge, reduction = "umap", group.by = "orig.ident") +
  ggtitle("No integration") +
  xlab("umap_1") +
  ylab("umap_2")

# compare:

title <- ggdraw() + 
  draw_label(paste0("Integration methods across all batches"), fontface = "bold", size = 24, hjust = 0.5) +
  theme(plot.background = element_rect(fill = "white", color = NA))
  
plot <- umap_global_cca + umap_global_rpca + umap_global_harmony + umap_global_merge

final_plot <- plot_grid(
  title, plot,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

ggsave("GBM_single_cell_analysis/outputs/QC/integration_plots_1.png",
       final_plot,
       height = 9.6, width = 12)

rm(final_plot, title)

tcell_integ <- FindNeighbors(tcell_integ, reduction = "integrated.cca", dims = 1:30)
tcell_integ <- FindClusters(tcell_integ, resolution = 2, cluster.name = "cca_clusters")
tcell_integ <- RunUMAP(tcell_integ, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
umap_tcell_cca <- DimPlot(tcell_integ, reduction = "umap.cca", group.by = "orig.ident")

tcell_integ <- FindNeighbors(tcell_integ, reduction = "integrated.rpca", dims = 1:30)
tcell_integ <- FindClusters(tcell_integ, resolution = 2, cluster.name = "rpca_clusters")
tcell_integ <- RunUMAP(tcell_integ, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
umap_tcell_rpca <- DimPlot(tcell_integ, reduction = "umap.rpca", group.by = "orig.ident")

tcell_integ <- FindNeighbors(tcell_integ, reduction = "harmony", dims = 1:30)
tcell_integ <- FindClusters(tcell_integ, resolution = 2, cluster.name = "harmony_clusters")
tcell_integ <- RunUMAP(tcell_integ, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
umap_tcell_harmony <- DimPlot(tcell_integ,
                              reduction = "umap.harmony",
                              group.by = "orig.ident",
                              pt.size = 0.1,
                              cols = c("1P" = "#7CAE00", "2P" = "#C77CFF")) +
  ggtitle("Batch 1P/2P") +
  xlab("umap_1") +
  ylab("umap_2") &
  NoLegend()

 tcell_merge

umap_tcell_merge <- DimPlot(tcell_merge, reduction = "umap", group.by = "orig.ident")

 myeloidIntegrated (3 methods) - use reductions
myeloid_integ <- FindNeighbors(myeloid_integ, reduction = "integrated.cca", dims = 1:30)
myeloid_integ <- FindClusters(myeloid_integ, resolution = 2, cluster.name = "cca_clusters")
myeloid_integ <- RunUMAP(myeloid_integ, reduction = "integrated.cca", dims = 1:30, reduction.name = "umap.cca")
umap_myeloid_cca <- DimPlot(myeloid_integ, reduction = "umap.cca", group.by = "orig.ident")

myeloid_integ <- FindNeighbors(myeloid_integ, reduction = "integrated.rpca", dims = 1:30)
myeloid_integ <- FindClusters(myeloid_integ, resolution = 2, cluster.name = "rpca_clusters")
myeloid_integ <- RunUMAP(myeloid_integ, reduction = "integrated.rpca", dims = 1:30, reduction.name = "umap.rpca")
umap_myeloid_rpca <- DimPlot(myeloid_integ, reduction = "umap.rpca", group.by = "orig.ident")

myeloid_integ <- FindNeighbors(myeloid_integ, reduction = "harmony", dims = 1:30)
myeloid_integ <- FindClusters(myeloid_integ, resolution = 2, cluster.name = "harmony_clusters")
myeloid_integ <- RunUMAP(myeloid_integ, reduction = "harmony", dims = 1:30, reduction.name = "umap.harmony")
umap_myeloid_harmony <- DimPlot(myeloid_integ,
                                reduction = "umap.harmony",
                                group.by = "orig.ident",
                                pt.size = 0.1) +
  ggtitle("Batch 1N/2N") +
  xlab("umap_1") +
  ylab("umap_2") &
  NoLegend()

 myeloid_merge

umap_myeloid_merge <- DimPlot(myeloid_merge, reduction = "umap", group.by = "orig.ident")

# Compare:

umap2_global_harmony <- DimPlot(global_integ,
                               reduction = "umap.harmony",
                               group.by = "orig.ident",
                               pt.size = 0.0025) +
  ggtitle("All batches") +
  xlab("umap_1") +
  ylab("umap_2")

title <- ggdraw() + 
  draw_label(paste0("Harmony integration by cell type versus all batches"), fontface = "bold", size = 24, hjust = 0.5) +
  theme(plot.background = element_rect(fill = "white", color = NA))

plot <- umap_myeloid_harmony + umap_tcell_harmony + umap2_global_harmony

final_plot <- plot_grid(
  title, plot,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

ggsave("GBM_single_cell_analysis/outputs/QC/integration_plots_2.png",
       final_plot,
       height = 4.8, width = 12)

rm(umap_global_cca, umap_global_merge, umap_global_rpca, umap_myeloid_cca,
   umap_myeloid_harmony, umap_myeloid_merge, umap_myeloid_rpca, umap_tcell_cca,
   umap_tcell_harmony, umap_tcell_merge, umap_tcell_rpca)

# Save ####

write.xlsx(sampleInfo, "GBM_single_cell_analysis/outputs/sampleInfo.xlsx")
saveRDS(objectsMyeloid_ready_to_integrate[[1]], "GBM_single_cell_analysis/outputs/objectsMyeloid_ready_to_integrate[[1]].RDS")
saveRDS(objectsMyeloid_ready_to_integrate[[2]], "GBM_single_cell_analysis/outputs/objectsMyeloid_ready_to_integrate[[2]].RDS")
saveRDS(objectsTCell_ready_to_integrate[[1]], "GBM_single_cell_analysis/outputs/objectsTCell_ready_to_integrate[[1]].RDS")
saveRDS(objectsTCell_ready_to_integrate[[2]], "GBM_single_cell_analysis/outputs/objectsTCell_ready_to_integrate[[2]].RDS")
saveRDS(thesis_seu, "GBM_single_cell_analysis/outputs/thesis_seu.RDS")

save.image("GBM_single_cell_analysis/outputs/thesis_single_cell_environment.RData")
