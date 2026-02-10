
# Tumor Microenvironment Analysis with IOBR Package
# IMPORTANT: Uses published IOBR package for immune cell deconvolution

# Install and load required packages
install_and_load_packages <- function() {
  # Install IOBR from GitHub if needed
  if (!requireNamespace("IOBR", quietly = TRUE)) {
    devtools::install_github("IOBR/IOBR")
  }
  
  # Required packages
  required_packages <- c(
    "tidyverse", "IOBR", "openxlsx", "beepr",
    "survival", "survminer"
  )
  
  # Install missing packages
  for (pkg in required_packages) {
    if (!requireNamespace(pkg, quietly = TRUE)) {
      install.packages(pkg)
    }
    library(pkg, character.only = TRUE)
  }
}

# Main TME analysis function
analyze_tumor_microenvironment <- function(expr_matrix, sample_name = "TCGA_LGG") {
  
  cat("Starting Tumor Microenvironment Analysis\n")
  cat("=======================================\n")
  cat("Sample:", sample_name, "\n")
  cat("Expression matrix dimensions:", dim(expr_matrix), "\n\n")
  
  # Ensure expression matrix is in correct format
  if (min(expr_matrix) < 0) {
    warning("Negative values detected. Converting to positive...")
    expr_matrix <- expr_matrix - min(expr_matrix) + 0.1
  }
  
  # Method 1: MCPcounter
  cat("1. Running MCPcounter...\n")
  im_mcpcounter <- deconvo_tme(
    eset = expr_matrix,
    method = "mcpcounter"
  )
  
  # Method 2: EPIC
  cat("2. Running EPIC...\n")
  im_epic <- deconvo_tme(
    eset = expr_matrix,
    method = "epic",
    arrays = FALSE
  )
  
  # Method 3: xCell
  cat("3. Running xCell...\n")
  im_xcell <- deconvo_tme(
    eset = expr_matrix,
    method = "xcell",
    arrays = FALSE
  )
  
  # Method 4: CIBERSORT
  cat("4. Running CIBERSORT...\n")
  im_cibersort <- deconvo_tme(
    eset = expr_matrix,
    method = "cibersort",
    arrays = FALSE,
    perm = 1000
  )
  
  # Method 5: Immunophenoscore (IPS)
  cat("5. Running Immunophenoscore...\n")
  im_ips <- deconvo_tme(
    eset = expr_matrix,
    method = "ips",
    plot = FALSE
  )
  
  # Method 6: quanTIseq
  cat("6. Running quanTIseq...\n")
  im_quantiseq <- deconvo_tme(
    eset = expr_matrix,
    method = "quantiseq",
    scale_mrna = TRUE
  )
  
  # Method 7: ESTIMATE
  cat("7. Running ESTIMATE...\n")
  im_estimate <- deconvo_tme(
    eset = expr_matrix,
    method = "estimate"
  )
  
  # Method 8: TIMER
  cat("8. Running TIMER...\n")
  im_timer <- deconvo_tme(
    eset = expr_matrix,
    method = "timer",
    group_list = rep("lgg", ncol(expr_matrix))
  )
  
  # Method 9: ssGSEA (if signature provided)
  cat("9. Running ssGSEA...\n")
  
  # Load ssGSEA signature if available
  ssgsea_file <- "InputData/ssGSEA28.Rdata"
  if (file.exists(ssgsea_file)) {
    load(ssgsea_file)
    im_ssgsea <- calculate_sig_score(
      eset = expr_matrix,
      signature = cellMarker,  # Must be loaded from ssGSEA28.Rdata
      method = "ssgsea"
    )
  } else {
    warning("ssGSEA signature file not found. Skipping ssGSEA analysis.")
    im_ssgsea <- NULL
  }
  
  # Save all results
  results <- list(
    mcpcounter = im_mcpcounter,
    epic = im_epic,
    xcell = im_xcell,
    cibersort = im_cibersort,
    ips = im_ips,
    quantiseq = im_quantiseq,
    estimate = im_estimate,
    timer = im_timer,
    ssgsea = im_ssgsea
  )
  
  # Save as RData
  save_file <- paste0("Results/TME_", sample_name, "_Results.Rdata")
  save(results, file = save_file)
  cat("\nResults saved to:", save_file, "\n")
  
  # Export to Excel
  export_to_excel(results, sample_name)
  
  # Completion signal
  beep("fanfare")
  
  return(results)
}

# Export results to Excel
export_to_excel <- function(results, sample_name) {
  cat("\nExporting results to Excel...\n")
  
  # Filter out NULL results
  valid_results <- results[!sapply(results, is.null)]
  
  # Create workbook
  wb <- createWorkbook()
  
  # Add each result as a sheet
  for (method_name in names(valid_results)) {
    addWorksheet(wb, method_name)
    writeData(wb, sheet = method_name, valid_results[[method_name]])
  }
  
  # Save workbook
  excel_file <- paste0("Results/TME_", sample_name, "_Results.xlsx")
  saveWorkbook(wb, excel_file, overwrite = TRUE)
  cat("Excel file saved to:", excel_file, "\n")
}

# Combine all TME results
combine_tme_results <- function(results) {
  cat("Combining all TME results...\n")
  
  # Start with first non-null result
  valid_results <- results[!sapply(results, is.null)]
  combined <- valid_results[[1]]
  
  # Join all results
  for (i in 2:length(valid_results)) {
    combined <- combined %>% 
      inner_join(valid_results[[i]], by = "ID")
  }
  
  cat("Combined matrix dimensions:", dim(combined), "\n")
  return(combined)
}

# Generate summary statistics
summarize_tme_results <- function(results) {
  cat("\nGenerating TME analysis summary...\n")
  
  summary_stats <- data.frame(
    Method = names(results),
    Samples = sapply(results, function(x) if(!is.null(x)) nrow(x) else NA),
    Cell_Types = sapply(results, function(x) if(!is.null(x)) ncol(x)-1 else NA),
    Status = sapply(results, function(x) if(!is.null(x)) "Complete" else "Skipped")
  )
  
  print(summary_stats)
  write.csv(summary_stats, "Results/TME_Analysis_Summary.csv", row.names = FALSE)
  
  return(summary_stats)
}

# Generate basic visualization
visualize_tme <- function(results, method = "cibersort", top_n = 12) {
  
  if (!method %in% names(results) || is.null(results[[method]])) {
    warning(paste("Method", method, "not available."))
    return(NULL)
  }
  
  # Create bar plot for top samples
  p <- cell_bar_plot(
    input = results[[method]][1:min(top_n, nrow(results[[method]])),],
    title = paste(method, "Cell Fractions"),
    pattern = "SYMBOL"
  )
  
  # Save plot
  ggsave(
    paste0("Figures/", method, "_Cell_Fractions.pdf"),
    p, width = 10, height = 6
  )
  
  return(p)
}

# Example usage function
example_tme_analysis <- function() {
  cat("Example Tumor Microenvironment Analysis Pipeline\n")
  cat("================================================\n\n")
  
  # Step 1: Load expression data
  cat("Step 1: Loading expression data...\n")
  
  # Example: Load TCGA-LGG data (modify as needed)
  load("Results/Data.pre.101.ML.Rdata")
  expr_matrix <- list_train_vali_Data[["TCGA.LGG"]][, -c(1:3)] %>% 
    t() %>% 
    as.data.frame()
  
  # Convert to appropriate format if needed
  expr_matrix <- 2^expr_matrix - 1
  
  cat("Expression matrix loaded. Dimensions:", dim(expr_matrix), "\n\n")
  
  # Step 2: Run TME analysis
  cat("Step 2: Running TME analysis...\n")
  tme_results <- analyze_tumor_microenvironment(
    expr_matrix = expr_matrix,
    sample_name = "TCGA_LGG"
  )
  
  # Step 3: Combine results
  cat("\nStep 3: Combining results...\n")
  combined_results <- combine_tme_results(tme_results)
  
  # Step 4: Generate summary
  cat("\nStep 4: Generating summary...\n")
  summary_stats <- summarize_tme_results(tme_results)
  
  # Step 5: Visualize
  cat("\nStep 5: Creating visualizations...\n")
  dir.create("Figures", showWarnings = FALSE, recursive = TRUE)
  
  # Visualize CIBERSORT results
  if (!is.null(tme_results$cibersort)) {
    visualize_tme(tme_results, method = "cibersort", top_n = 12)
  }
  
  cat("\nAnalysis complete!\n")
  cat("Check 'Results/' and 'Figures/' directories for outputs.\n")
  
  return(list(
    results = tme_results,
    combined = combined_results,
    summary = summary_stats
  ))
}

# Run analysis (uncomment to use)
if (FALSE) {
  # Install and load packages
  install_and_load_packages()
  
  # Run example analysis
  analysis_output <- example_tme_analysis()
}
```

## English Documentation

# Tumor Microenvironment Analysis with IOBR

## Overview
Uses **published IOBR package** to perform 9 immune deconvolution methods from gene expression data.

## Methods Included
1. **MCPcounter** - Cell type quantification
2. **EPIC** - Estimating Proportions of Immune Cells
3. **xCell** - Cell type enrichment
4. **CIBERSORT** - Relative subset recognition
5. **Immunophenoscore (IPS)** - Immunogenicity score
6. **quanTIseq** - Quantifying immune context
7. **ESTIMATE** - Stromal/immune scores
8. **TIMER** - Tumor Immune Estimation Resource
9. **ssGSEA** - Single-sample Gene Set Enrichment (requires signature)

## Input Requirements
Gene expression matrix (samples × genes):
  - Raw counts or normalized expression
- Gene symbols as row names
- Sample IDs as column names

## Usage
```r
# Run complete analysis
source("tme_analysis.R")
analysis_results <- example_tme_analysis()

# Or run custom analysis:
tme_results <- analyze_tumor_microenvironment(
  expr_matrix = your_expression_matrix,
  sample_name = "Your_Sample_Name"
)
```

## Output Files
### Results Directory:
- `TME_TCGA_LGG_Results.Rdata` - All results (R format)
- `TME_TCGA_LGG_Results.xlsx` - Excel with separate sheets
- `TME_Analysis_Summary.csv` - Method statistics

### Figures Directory:
- `cibersort_Cell_Fractions.pdf` - Sample visualization

## Key Features
- **Standardized Output**: Uniform format across all methods
- **Easy Integration**: Results can be merged seamlessly
- **Comprehensive**: Covers major published algorithms
- **Citation Ready**: Includes reference information

## Note on Algorithms
All methods are **implemented by IOBR package authors**. This script provides a wrapper for consistent analysis workflow.

## References
Zeng D, et al. IOBR: Multi-Omics Immuno-Oncology Biological Research to Decode Tumor Microenvironment and Signatures. Frontiers in Immunology. 2021. doi: 10.3389/fimmu.2021.687975

*Note: We are users of the IOBR package, not developers*