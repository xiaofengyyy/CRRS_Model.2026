
# Single-cell CNV Inference using CopyKAT Package
# IMPORTANT: Uses published CopyKAT package for CNV inference

# Load required packages
library(Seurat)
library(copykat)
library(dplyr)
library(ggplot2)
library(patchwork)
library(AnnoProbe)

# Prepare single-cell data
prepare_sc_data <- function(sce) {
  # Create Seurat object
  sce_obj <- CreateSeuratObject(
    counts = sce@assays$RNA$counts,
    meta.data = sce@meta.data
  )
  
  # Standard preprocessing
  sce_obj <- NormalizeData(sce_obj, 
                           normalization.method = "LogNormalize", 
                           scale.factor = 1e4)
  sce_obj <- FindVariableFeatures(sce_obj, selection.method = "vst", nfeatures = 2000)
  sce_obj <- ScaleData(sce_obj)
  sce_obj <- RunPCA(sce_obj, npcs = 50)
  
  return(sce_obj)
}

# Prepare CopyKAT input files
prepare_copykat_input <- function(sce_obj, celltype_column = "celltype") {
  # Set cell identities
  Idents(sce_obj) <- sce_obj@meta.data[[celltype_column]]
  
  # Get expression matrix
  exp_matrix <- GetAssayData(sce_obj, layer = 'counts', assay = 'RNA')
  
  # Prepare group information
  group_info <- data.frame(
    cell_id = colnames(exp_matrix),
    cell_type = Idents(sce_obj)
  )
  
  # Prepare gene annotation (remove sex chromosomes optionally)
  gene_info <- annoGene(rownames(exp_matrix), "SYMBOL", 'human')
  gene_info <- gene_info[, c("SYMBOL", "chr", "start", "end")]
  gene_info <- gene_info[!duplicated(gene_info$SYMBOL), ]
  
  # Filter and order expression matrix
  exp_matrix <- exp_matrix[rownames(exp_matrix) %in% gene_info$SYMBOL, ]
  exp_matrix <- exp_matrix[match(gene_info$SYMBOL, rownames(exp_matrix)), ]
  
  # Identify normal cells (reference)
  normal_cells <- colnames(exp_matrix)[grepl("T cell|normal|stroma", 
                                             group_info$cell_type, 
                                             ignore.case = TRUE)]
  
  if (length(normal_cells) == 0) {
    warning("No normal cells identified for reference. Using all cells.")
    normal_cells <- NULL
  }
  
  return(list(
    exp_matrix = exp_matrix,
    group_info = group_info,
    gene_info = gene_info,
    normal_cells = normal_cells
  ))
}

# Run CopyKAT analysis
run_copykat <- function(exp_matrix, sample_name, normal_cells = NULL, n_cores = 4) {
  
  copykat_result <- copykat(
    rawmat = as.data.frame(as.matrix(exp_matrix)),
    id.type = "S",
    ngene.chr = 5,
    win.size = 25,
    KS.cut = 0.1,
    sam.name = sample_name,
    distance = "euclidean",
    norm.cell.names = normal_cells,
    n.cores = n_cores,
    output.seg = FALSE
  )
  
  return(copykat_result)
}

# Process and visualize results
process_copykat_results <- function(sce_obj, copykat_result) {
  # Extract predictions
  predictions <- data.frame(copykat_result$prediction)
  
  # Add predictions to metadata
  sce_obj$copykat_pred <- predictions$copykat.pred[match(
    colnames(sce_obj), 
    predictions$cell.names
  )]
  
  # Remove NA predictions
  sce_obj <- sce_obj[, !is.na(sce_obj$copykat_pred)]
  
  # Calculate statistics
  stats <- table(sce_obj$copykat_pred)
  malignant_ratio <- sum(sce_obj$copykat_pred == "aneuploid") / ncol(sce_obj)
  
  cat("CopyKAT Prediction Summary:\n")
  cat("Total cells:", ncol(sce_obj), "\n")
  cat("Aneuploid (malignant):", sum(sce_obj$copykat_pred == "aneuploid"), 
      sprintf("(%.1f%%)", malignant_ratio*100), "\n")
  cat("Diploid (normal):", sum(sce_obj$copykat_pred == "diploid"), "\n")
  
  return(sce_obj)
}

# Generate visualizations
visualize_copykat <- function(sce_obj, output_dir = "copykat_results") {
  
  dir.create(output_dir, showWarnings = FALSE, recursive = TRUE)
  
  # UMAP colored by CopyKAT prediction
  p1 <- DimPlot(sce_obj, 
                group.by = "copykat_pred",
                pt.size = 1,
                label = FALSE,
                raster = FALSE) +
    labs(title = "CopyKAT Predictions",
         subtitle = "Aneuploid (malignant) vs Diploid (normal)") +
    theme_minimal()
  
  # UMAP colored by original cell type
  if ("celltype" %in% colnames(sce_obj@meta.data)) {
    p2 <- DimPlot(sce_obj,
                  group.by = "celltype",
                  pt.size = 1,
                  label = FALSE,
                  raster = FALSE) +
      labs(title = "Original Cell Type Annotation") +
      theme_minimal()
    
    # Combined plot
    combined_plot <- p2 | p1
    ggsave(file.path(output_dir, "copykat_vs_celltype.pdf"), 
           combined_plot, width = 12, height = 5)
  }
  
  # Save individual plots
  ggsave(file.path(output_dir, "copykat_predictions_umap.pdf"), 
         p1, width = 6, height = 5)
  
  # Bar plot of predictions
  pred_counts <- as.data.frame(table(sce_obj$copykat_pred))
  colnames(pred_counts) <- c("Prediction", "Count")
  
  p3 <- ggplot(pred_counts, aes(x = Prediction, y = Count, fill = Prediction)) +
    geom_bar(stat = "identity") +
    geom_text(aes(label = Count), vjust = -0.5) +
    scale_fill_manual(values = c("aneuploid" = "#E41A1C", "diploid" = "#377EB8")) +
    labs(title = "CopyKAT Prediction Distribution",
         x = "Prediction", y = "Number of Cells") +
    theme_minimal() +
    theme(legend.position = "none")
  
  ggsave(file.path(output_dir, "copykat_prediction_distribution.pdf"), 
         p3, width = 5, height = 5)
  
  # Save results
  write.csv(as.data.frame(table(sce_obj$copykat_pred)),
            file.path(output_dir, "copykat_prediction_summary.csv"),
            row.names = FALSE)
}

# Main analysis pipeline
analyze_scn_cnv <- function(sce, sample_name = "sample", 
                            celltype_column = "celltype",
                            n_cores = 4) {
  
  cat("Starting CopyKAT CNV analysis for:", sample_name, "\n")
  cat("===============================================\n")
  
  # Step 1: Prepare single-cell data
  cat("1. Preparing single-cell data...\n")
  sce_obj <- prepare_sc_data(sce)
  
  # Step 2: Prepare CopyKAT input
  cat("2. Preparing CopyKAT input files...\n")
  input_data <- prepare_copykat_input(sce_obj, celltype_column)
  
  cat("   Expression matrix dimensions:", dim(input_data$exp_matrix), "\n")
  cat("   Reference (normal) cells identified:", length(input_data$normal_cells), "\n")
  
  # Step 3: Run CopyKAT
  cat("3. Running CopyKAT analysis...\n")
  copykat_result <- run_copykat(
    exp_matrix = input_data$exp_matrix,
    sample_name = sample_name,
    normal_cells = input_data$normal_cells,
    n_cores = n_cores
  )
  
  # Save raw results
  save(copykat_result, 
       file = file.path("copykat_results", paste0(sample_name, "_copykat_raw.Rdata")))
  
  # Step 4: Process results
  cat("4. Processing CopyKAT results...\n")
  sce_obj <- process_copykat_results(sce_obj, copykat_result)
  
  # Step 5: Visualize
  cat("5. Generating visualizations...\n")
  visualize_copykat(sce_obj, output_dir = "copykat_results")
  
  # Step 6: Save final annotated object
  cat("6. Saving final results...\n")
  saveRDS(sce_obj, file.path("copykat_results", paste0(sample_name, "_annotated.rds")))
  
  cat("\nAnalysis complete!\n")
  cat("Results saved to: copykat_results/\n")
  
  return(list(
    seurat_obj = sce_obj,
    copykat_result = copykat_result
  ))
}

# Example usage (uncomment and modify as needed)
if (FALSE) {
  # Load your single-cell data
  # sce <- readRDS("your_single_cell_data.rds")
  
  # Run analysis
  results <- analyze_scn_cnv(
    sce = sce,
    sample_name = "Glioma_GSE270109",
    celltype_column = "celltype",
    n_cores = 12
  )
  
  # Access results
  seurat_obj <- results$seurat_obj
  copykat_result <- results$copykat_result
}
```

## English Documentation

# Single-cell CNV Analysis with CopyKAT

## Overview
Uses **published CopyKAT package** to infer copy number variations from scRNA-seq data. Identifies malignant vs. normal cells based on aneuploidy patterns.

## Key Features
- **CNV Inference**: Genome-wide copy number at 5MB resolution
- **Cell Classification**: Distinguishes aneuploid (malignant) vs. diploid (normal) cells
- **Quality Control**: Automatic reference cell selection
- **Visualization**: UMAP plots with CNV predictions

## Input Requirements
Single-cell data (Seurat object) with:
  - RNA count matrix
- Cell type annotations (column name specified)

## Processing Steps
1. **Data Preparation**: Filter genes, prepare expression matrix
2. **Reference Selection**: Identify normal cells (T cells/stroma)
3. **CopyKAT Analysis**: Run CNV inference
4. **Classification**: Label cells as aneuploid/diploid
5. **Visualization**: UMAP plots and statistics

## Output Files
- `copykat_results/sample_copykat_raw.Rdata`: Raw CopyKAT results
- `copykat_results/sample_annotated.rds`: Annotated Seurat object
- `copykat_predictions_umap.pdf`: UMAP with CNV predictions
- `copykat_prediction_summary.csv`: Cell count statistics

## Usage
```r
# Run analysis
results <- analyze_scn_cnv(
  sce = your_seurat_object,
  sample_name = "SampleName",
  celltype_column = "celltype_column_name",
  n_cores = 12
)
```

## Important Notes
- Uses **published CopyKAT algorithms**, not custom methods
- Normal cells needed as reference (auto-detected from T/stromal cells)
- Requires sufficient cell numbers (>100 recommended)
- Computation intensive - adjust `n_cores` based on system

## References
- CopyKAT: Delineating copy number and clonal substructure in human tumors from single-cell transcriptomes
- Package: https://github.com/navinlabcode/copykat