#Code for glioma single-cell RNA-seq analysis:

# Glioma Single-cell RNA-seq Analysis Pipeline
# Using Seurat for preprocessing, clustering, and annotation

# Load required packages
library(Seurat)
library(tidyverse)
library(patchwork)
library(harmony)
library(cowplot)
library(clustree)
library(ComplexHeatmap)
library(circlize)
library(dittoSeq)
library(openxlsx)

# Set directories
results_dir <- "Results"
figures_dir <- "Figures"
dir.create(results_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(figures_dir, showWarnings = FALSE, recursive = TRUE)

# Define color palette for visualizations
color_palette <- c(
  '#E5D2DD', '#53A85F', '#F1BB72', '#F3B1A0', '#D6E7A3',
'#57C3F3', '#E95C59', '#E59CC4', '#AB3282', '#BD956A',
'#6778AE', '#91D0BE', '#B53E2B', '#712820', '#DCC1DD',
'#CCE0F5', '#CCC9E6', '#625D9E', '#68A180', '#3A6963',
'#968175', '#F0E68C', '#FFFFE0', '#EE82EE', '#FF6347',
'#8B0000', '#FF4500', '#FFD700', '#2E8B57', '#FF69B4',
'#6A5ACD', '#8FBC8F', '#FF8C00', '#FF1493', '#3CB371',
'#FFB6C1', '#B8860B', '#FF00FF', '#FFA07A', '#20B2AA',
'#87CEFA', '#BA55D3', '#CD5C5C', '#FA8072', '#48D1CC'
)

# Glioma-specific marker genes
glioma_markers <- list(
  Astrocytes = c("GFAP", "AQP4"),
  T_cells = c("CD3D", "CD3E", "CD4", "CD8A"),
  Microglia = c("TMEM119", "P2RY12", "CX3CR1"),
  Macrophages = c("CD68", "CD163", "CSF1R", "SIGLEC9"),
  Monocytes = c("CD14", "FCGR3A"),
  NK_cells = c("NCAM1", "KLRD1", "NKG7"),
  B_cells = c("CD19", "CD79A", "MS4A1"),
  Oligodendrocytes = c("MBP", "MOG", "PLP1"),
  OPCs = c("PDGFRA", "CSPG4"),
  Dendritic_cells = c("ITGAX", "CLEC9A", "CD1C", "LILRA4"),
  Endothelial = c("PECAM1", "VWF", "CLDN5"),
  Mural_cells = c("ACTA2", "PDGFRB", "RGS5"),
  Mast_cells = c("KIT", "TPSAB1", "CPA3"),
  Cycling_cells = c("MKI67", "TOP2A", "PCNA"),
  #Malignant = c("EPCAM", "KRT18", "KRT8") # CNV-based annotation preferred for malignant cells
)

# Quality control and preprocessing
preprocess_sc_data <- function(seurat_object) {
  cat("Starting quality control...\n")
  
  # Calculate quality metrics
  seurat_object[["percent.mt"]] <- PercentageFeatureSet(seurat_object, pattern = "^MT-")
  seurat_object[["percent.rp"]] <- PercentageFeatureSet(seurat_object, pattern = "^RP[SL]")
  seurat_object[["percent.hb"]] <- PercentageFeatureSet(seurat_object, pattern = "^HB[^(P)]")
  
  # Generate QC plot
  qc_plot <- VlnPlot(
    object = seurat_object,
    features = c("nCount_RNA", "nFeature_RNA", "percent.mt"),
    cols = color_palette,
    ncol = 3,
    group.by = "orig.ident",
    pt.size = 0
  )
  
  ggsave(
    file.path(figures_dir, "1_pre_QC.pdf"),
    plot = qc_plot,
    width = 16,
    height = 4
  )
  
  # Filter cells
  cat("Filtering cells...\n")
  filtered_object <- subset(seurat_object,
                            subset = nFeature_RNA >= 300 & nFeature_RNA <= 6000 &
                              percent.mt >= 0 & percent.mt <= 15 &
                              percent.hb >= 0 & percent.hb <= 0.1 &
                              percent.rp >= 1 & percent.rp <= 100
  )
  
  cat("Initial cells:", ncol(seurat_object), "\n")
  cat("After filtering:", ncol(filtered_object), "\n")
  
  return(filtered_object)
}

# Data integration with Harmony
integrate_data <- function(seurat_object) {
  cat("Normalizing and integrating data...\n")
  
  # Standard preprocessing
  seurat_object <- NormalizeData(seurat_object,
                                 normalization.method = "LogNormalize",
                                 scale.factor = 10000
  )
  seurat_object <- FindVariableFeatures(seurat_object,
                                        selection.method = "vst",
                                        nfeatures = 2000
  )
  seurat_object <- ScaleData(seurat_object, features = rownames(seurat_object))
  seurat_object <- RunPCA(seurat_object,
                          features = VariableFeatures(object = seurat_object)
  )
  
  # Harmony integration
  cat("Running Harmony integration...\n")
  integrated_object <- RunHarmony(seurat_object,
                                  group.by.vars = "orig.ident",
                                  plot_convergence = TRUE,
                                  max.iter = 20,
                                  reduction.save = "harmony"
  )
  
  return(integrated_object)
}

# Clustering and UMAP
cluster_cells <- function(seurat_object, dims = 30) {
  cat("Finding neighbors and clustering...\n")
  
  # Find neighbors
  seurat_object <- FindNeighbors(seurat_object,
                                 reduction = "harmony",
                                 dims = 1:dims
  )
  
  # Test multiple resolutions
  resolutions <- c(0.01, 0.05, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1)
  
  for (res in resolutions) {
    seurat_object <- FindClusters(seurat_object,
                                  resolution = res,
                                  algorithm = 1
    )
  }
  
  # Generate resolution comparison plots
  pdf(file.path(figures_dir, "Clustering_Resolution_Comparison.pdf"),
      width = 14, height = 10
  )
  
  # DimPlots for different resolutions
  p1 <- plot_grid(
    ncol = 3,
    DimPlot(seurat_object, reduction = "umap", group.by = "RNA_snn_res.0.1") +
      ggtitle("Resolution 0.1"),
    DimPlot(seurat_object, reduction = "umap", group.by = "RNA_snn_res.0.3") +
      ggtitle("Resolution 0.3"),
    DimPlot(seurat_object, reduction = "umap", group.by = "RNA_snn_res.0.6") +
      ggtitle("Resolution 0.6")
  )
  
  print(p1)
  
  # Clustree plot
  p2 <- clustree(seurat_object@meta.data, prefix = "RNA_snn_res.")
  print(p2)
  
  dev.off()
  
  # Run UMAP and t-SNE
  cat("Running UMAP and t-SNE...\n")
  seurat_object <- RunUMAP(seurat_object,
                           reduction = "harmony",
                           dims = 1:dims
  )
  seurat_object <- RunTSNE(seurat_object,
                           reduction = "harmony",
                           dims = 1:dims
  )
  
  # Save intermediate results
  saveRDS(seurat_object,
          file = file.path(results_dir, "integrated_seurat_object.rds")
  )
  
  return(seurat_object)
}

# Find marker genes
find_marker_genes <- function(seurat_object, ident_column = "seurat_clusters") {
  cat("Finding marker genes...\n")
  
  Idents(seurat_object) <- seurat_object@meta.data[[ident_column]]
  
  all_markers <- FindAllMarkers(
    object = seurat_object,
    only.pos = TRUE,
    min.pct = 0.25,
    logfc.threshold = 0.25,
    test.use = "wilcox"
  )
  
  # Save markers
  write.xlsx(all_markers,
             file = file.path(results_dir, "all_markers.xlsx")
  )
  
  save(all_markers,
       file = file.path(results_dir, "all_markers.Rdata")
  )
  
  return(all_markers)
}

# Generate marker heatmap
create_marker_heatmap <- function(seurat_object, marker_genes, n_top = 5) {
  cat("Creating marker heatmap...\n")
  
  # Select top markers per cluster
  markers_df <- find_marker_genes(seurat_object)
  top_markers <- markers_df %>%
    group_by(cluster) %>%
    top_n(n = n_top, wt = avg_log2FC) %>%
    as.data.frame()
  
  # Calculate average expression
  avg_expr <- AverageExpression(seurat_object,
                                features = top_markers$gene,
                                group.by = "seurat_clusters",
                                slot = "data"
  )$RNA
  
  # Scale and prepare for heatmap
  scaled_expr <- t(scale(t(avg_expr)))
  
  # Color setup
  heatmap_colors <- colorRamp2(c(-2, 0, 2), c("#0da9ce", "white", "#e74a32"))
  
  # Create heatmap
  pdf(file.path(figures_dir, "Marker_Heatmap.pdf"),
      width = 20, height = 18
  )
  
  Heatmap(scaled_expr,
          name = "Expression",
          col = heatmap_colors,
          cluster_columns = FALSE,
          cluster_rows = FALSE,
          column_names_side = "top",
          column_names_rot = 60,
          rect_gp = gpar(col = "white", lwd = 1.5)
  )
  
  dev.off()
}

# Annotate cell types
annotate_cells <- function(seurat_object, resolution = "RNA_snn_res.0.6") {
  cat("Annotating cell types...\n")
  
  # Set cluster identity
  seurat_object <- SetIdent(seurat_object, value = resolution)
  
  # Define cluster to cell type mapping
  cluster_mapping <- c(
    "0" = "T_cells",
    "1" = "Malignant",
    "2" = "Microglia",
    "3" = "Microglia",
    "4" = "NK_cells",
    "5" = "T_cells",
    "6" = "Mural_cells",
    "7" = "Cycling_cells",
    "8" = "Astrocytes",
    "9" = "Mural_cells",
    "10" = "Oligodendrocytes",
    "11" = "Endothelial",
    "12" = "Monocytes",
    "13" = "Astrocytes",
    "14" = "Dendritic_cells",
    "15" = "B_cells",
    "16" = "Macrophages",
    "17" = "T_cells",
    "18" = "Cycling_cells",
    "19" = "Astrocytes",
    "20" = "Mural_cells",
    "21" = "T_cells",
    "22" = "Mast_cells",
    "23" = "Oligodendrocytes",
    "24" = "Cycling_cells",
    "25" = "Cycling_cells",
    "26" = "Mural_cells",
    "27" = "T_cells",
    "28" = "Astrocytes"
  )
  
  # Apply mapping
  seurat_object@meta.data$cell_type <- cluster_mapping[as.character(seurat_object@meta.data[[resolution]])]
  Idents(seurat_object) <- seurat_object@meta.data$cell_type
  
  return(seurat_object)
}

# Generate visualizations
generate_visualizations <- function(seurat_object) {
  cat("Generating visualizations...\n")
  
  # UMAP by sample
  p1 <- DimPlot(seurat_object,
                reduction = "umap",
                group.by = "orig.ident",
                cols = color_palette,
                label = FALSE
  ) +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
  
  ggsave(
    file.path(figures_dir, "UMAP_by_Sample.pdf"),
    plot = p1,
    width = 7,
    height = 4
  )
  
  # UMAP by cell type
  p2 <- DimPlot(seurat_object,
                reduction = "umap",
                group.by = "cell_type",
                cols = color_palette,
                label = TRUE
  ) +
    theme_bw() +
    theme(
      legend.position = "right",
      axis.title.x = element_text(size = 12),
      axis.title.y = element_text(size = 12)
    )
  
  ggsave(
    file.path(figures_dir, "UMAP_by_CellType.pdf"),
    plot = p2,
    width = 6.5,
    height = 4
  )
  
  # Dot plot of marker genes
  p3 <- DotPlot(seurat_object,
                features = glioma_markers,
                assay = "RNA"
  ) +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1),
      axis.text = element_text(size = 10, color = "black")
    )
  
  ggsave(
    file.path(figures_dir, "Marker_DotPlot.pdf"),
    plot = p3,
    width = 12,
    height = 5
  )
  
  # Cell composition bar plots
  p4 <- dittoBarPlot(seurat_object,
                     "cell_type",
                     group.by = "orig.ident"
  ) +
    guides(col = guide_legend(ncol = 1)) +
    theme(legend.position = "right")
  
  p5 <- dittoBarPlot(seurat_object,
                     "orig.ident",
                     group.by = "cell_type"
  ) +
    guides(col = guide_legend(ncol = 1)) +
    theme(legend.position = "right")
  
  combined_plot <- p4 | p5
  ggsave(
    file.path(figures_dir, "Cell_Composition.pdf"),
    plot = combined_plot,
    width = 14,
    height = 6
  )
}

# Main analysis pipeline
run_glioma_scrnaseq_analysis <- function(input_seurat_object) {
  cat("Starting Glioma Single-cell Analysis Pipeline\n")
  cat("============================================\n\n")
  
  # Step 1: Quality control
  cat("Step 1: Quality Control\n")
  filtered_data <- preprocess_sc_data(input_seurat_object)
  
  # Step 2: Integration
  cat("\nStep 2: Data Integration\n")
  integrated_data <- integrate_data(filtered_data)
  
  # Step 3: Clustering
  cat("\nStep 3: Clustering\n")
  clustered_data <- cluster_cells(integrated_data, dims = 30)
  
  # Step 4: Annotation
  cat("\nStep 4: Cell Type Annotation\n")
  annotated_data <- annotate_cells(clustered_data)
  
  # Step 5: Marker gene analysis
  cat("\nStep 5: Marker Gene Analysis\n")
  create_marker_heatmap(annotated_data, glioma_markers, n_top = 5)
  
  # Step 6: Visualization
  cat("\nStep 6: Visualization\n")
  generate_visualizations(annotated_data)
  
  # Save final object
  cat("\nSaving final results...\n")
  saveRDS(annotated_data,
          file = file.path(results_dir, "final_annotated_seurat.rds")
  )
  
  # Save metadata
  write.csv(annotated_data@meta.data,
            file = file.path(results_dir, "cell_metadata.csv"),
            row.names = TRUE
  )
  
  cat("\nAnalysis complete!\n")
  cat("Results saved to:", results_dir, "\n")
  cat("Figures saved to:", figures_dir, "\n")
  
  return(annotated_data)
}

# Example usage (uncomment and modify as needed)
if (FALSE) {
  # Load your single-cell data
  # input_data <- readRDS("your_glioma_scRNA_data.rds")
  
  # Run the analysis pipeline
  # results <- run_glioma_scrnaseq_analysis(input_data)
}


## English Documentation

# Glioma Single-cell RNA-seq Analysis

## Overview
Complete pipeline for processing, analyzing, and visualizing glioma single-cell RNA-seq data using Seurat package.

## Key Features
- **Quality Control**: Filtering based on gene counts and mitochondrial percentage
- **Data Integration**: Harmony batch correction for multi-sample datasets
- **Clustering**: Multi-resolution clustering with UMAP visualization
- **Cell Annotation**: Glioma-specific cell type identification
- **Marker Analysis**: Identification of cell type-specific genes
- **Visualization**: UMAP plots, heatmaps, and composition bar plots

## Input Requirements
Seurat object containing:
  - Raw count matrix
- Metadata with sample IDs (`orig.ident`)
- Cell barcodes and gene symbols

## Processing Steps
1. **Quality Control**: Filter low-quality cells
2. **Normalization**: Log normalization and scaling
3. **Integration**: Harmony batch correction
4. **Clustering**: Multi-resolution Louvain clustering
5. **Annotation**: Assign glioma cell types
6. **Marker Identification**: Find cluster-specific genes
7. **Visualization**: Generate publication-ready figures

## Output Files
### Results Directory:
- `integrated_seurat_object.rds` - Intermediate Seurat object
- `final_annotated_seurat.rds` - Final annotated object
- `all_markers.xlsx` - Marker genes per cluster
- `cell_metadata.csv` - Complete cell metadata

### Figures Directory:
- QC violin plots
- UMAP visualizations (by sample and cell type)
- Marker gene dot plots
- Heatmap of top marker genes
- Cell composition bar plots

## Cell Types Annotated
- **Malignant cells**: Glioma tumor cells
- **T cells**: CD4+ and CD8+ T lymphocytes
- **Microglia**: Brain-resident immune cells
- **Macrophages**: Infiltrating immune cells
- **Astrocytes**: Glial support cells
- **Oligodendrocytes**: Myelin-producing cells
- **Endothelial cells**: Blood vessel lining
- **Mural cells**: Pericytes and smooth muscle
- **Dendritic cells**: Antigen-presenting cells
- **Cycling cells**: Proliferating cells

## Usage
r
# Load your Seurat object
# your_data <- Read10X("path/to/data") %>% CreateSeuratObject()

# Run analysis pipeline
source("glioma_scrnaseq_analysis.R")
results <- run_glioma_scrnaseq_analysis(your_data)


## Dependencies
- Seurat (single-cell analysis)
- Harmony (batch correction)
- tidyverse (data manipulation)
- ComplexHeatmap (heatmap visualization)
- dittoSeq (composition plots)
- openxlsx (Excel export)

## Notes
- Pipeline optimized for glioma datasets
- Adjust clustering resolution based on dataset size
- Marker gene lists can be customized for specific glioma subtypes
- All visualizations use colorblind-friendly palettes

## References
- Seurat: https://satijalab.org/seurat/
  - Harmony: https://github.com/immunogenomics/harmony
- Glioma cell atlas references integrated into annotation scheme