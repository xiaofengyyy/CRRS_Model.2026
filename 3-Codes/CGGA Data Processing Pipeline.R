# CGGA Data Processing Pipeline
# Chinese Glioma Genome Atlas datasets (CGGA301 and CGGA325)
library(survival)
library(tidyverse)
library(data.table)

# Function to process CGGA dataset
process_cgga_dataset <- function(dataset_name) {
  # Dataset-specific file paths and parameters
  if (dataset_name == "CGGA301") {
    clinical_file <- "CGGA.mRNA_array_301_clinical.20200506.txt"
    expression_file <- "CGGA.mRNA_array_301_gene_level.20200506.txt"
    time_divisor <- 30  # Convert days to months
    
  } else if (dataset_name == "CGGA325") {
    clinical_file <- "CGGA.mRNAseq_325_clinical.20200506.txt"
    expression_file <- "CGGA.mRNAseq_325.RSEM-genes.20200506.txt"
    time_divisor <- 30  # Convert days to months
    
  } else {
    stop("Unknown dataset. Use 'CGGA301' or 'CGGA325'")
  }
  
  # Load and process clinical data
  cat("Processing", dataset_name, "dataset...\n")
  
  clinical_data <- fread(clinical_file, data.table = FALSE)
  expression_data <- fread(expression_file, data.table = FALSE) %>%
    column_to_rownames("Gene_Name")
  
  # Prepare expression matrix
  expr_matrix <- t(expression_data) %>% as.data.frame()
  
  # Prepare clinical metadata
  if (dataset_name == "CGGA301") {
    meta <- clinical_data %>%
      dplyr::select(CGGA_ID, OS.time = OS, OS = Censor..alive.0..dead.1.)
  } else if (dataset_name == "CGGA325") {
    meta <- clinical_data %>%
      dplyr::select(CGGA_ID, OS.time = OS, OS = Censor..alive.0..dead.1.)
  }
  
  rownames(meta) <- meta$CGGA_ID
  
  # Find common samples
  common_samples <- intersect(rownames(meta), rownames(expr_matrix))
  meta <- meta[common_samples, ]
  expr_matrix <- expr_matrix[common_samples, ]
  
  # Verify alignment
  if (!identical(rownames(meta), rownames(expr_matrix))) {
    stop("Sample alignment failed for", dataset_name)
  }
  
  # Prepare final survival-expression matrix
  survival_info <- meta %>%
    dplyr::select(ID = CGGA_ID, OS.time, OS)
  
  sur_exp <- cbind(survival_info, expr_matrix) %>%
    dplyr::select(ID, OS.time, OS, everything())
  
  # Process survival time
  sur_exp$OS.time <- as.numeric(sur_exp$OS.time) / time_divisor
  
  # Filter valid samples
  sur_exp <- sur_exp %>%
    na.omit() %>%
    filter(OS.time > 0)
  
  # Apply log2 transformation for RNA-seq data (CGGA325)
  if (dataset_name == "CGGA325") {
    sur_exp[, 4:ncol(sur_exp)] <- log2(sur_exp[, 4:ncol(sur_exp)] + 1)
  }
  
  # Save results
  output_file <- paste0(dataset_name, "_sur_exp.Rdata")
  save(sur_exp, file = output_file)
  
  cat("Successfully processed", nrow(sur_exp), "samples from", dataset_name, "\n")
  cat("Output saved to:", output_file, "\n\n")
  
  return(sur_exp)
}

# Function to filter GBM samples from processed data
filter_gbm_samples <- function(dataset_name) {
  cat("Filtering GBM samples from", dataset_name, "...\n")
  
  # Load clinical data to identify GBM samples
  if (dataset_name == "CGGA301") {
    clinical_file <- "CGGA.mRNA_array_301_clinical.20200506.txt"
    expression_file <- "CGGA.mRNA_array_301_gene_level.20200506.txt"
  } else if (dataset_name == "CGGA325") {
    clinical_file <- "CGGA.mRNAseq_325_clinical.20200506.txt"
    expression_file <- "CGGA.mRNAseq_325.RSEM-genes.20200506.txt"
    count_file <- "CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt"
  } else {
    stop("Unknown dataset. Use 'CGGA301' or 'CGGA325'")
  }
  
  # Load clinical data
  clinical_data <- fread(clinical_file, data.table = FALSE)
  
  # Identify GBM samples
  gbm_samples <- clinical_data %>%
    filter(str_detect(toupper(Histology), "GBM|GLIOBLASTOMA"))
  
  gbm_ids <- gbm_samples$CGGA_ID
  
  # Load processed data
  load(paste0(dataset_name, "_sur_exp.Rdata"))
  
  # Filter for GBM samples
  sur_exp_gbm <- sur_exp %>%
    filter(ID %in% gbm_ids) %>%
    na.omit()
  
  # Load and filter expression matrices
  if (dataset_name == "CGGA301") {
    expression_data <- fread(expression_file, data.table = FALSE) %>%
      column_to_rownames("Gene_Name")
    expression_gbm <- expression_data[, gbm_ids]
    
    # Save GBM subset
    save(sur_exp_gbm, clinical_data, gbm_samples, expression_gbm,
         file = paste0(dataset_name, "_sur_exp_all_GBM.Rdata"))
    
  } else if (dataset_name == "CGGA325") {
    expression_rsem <- fread(expression_file, data.table = FALSE) %>%
      column_to_rownames("Gene_Name")
    expression_count <- fread(count_file, data.table = FALSE) %>%
      column_to_rownames("gene_name")
    
    expression_rsem_gbm <- expression_rsem[, gbm_ids]
    expression_count_gbm <- expression_count[, gbm_ids]
    
    # Save GBM subset
    save(sur_exp_gbm, clinical_data, gbm_samples, 
         expression_rsem_gbm, expression_count_gbm,
         file = paste0(dataset_name, "_sur_exp_all_GBM.Rdata"))
  }
  
  cat("GBM filtering completed:", nrow(sur_exp_gbm), "GBM samples identified\n")
  cat("Output saved to:", paste0(dataset_name, "_sur_exp_all_GBM.Rdata"), "\n\n")
  
  return(sur_exp_gbm)
}

# Main processing pipeline
process_cgga_pipeline <- function() {
  cat("Starting CGGA Data Processing Pipeline\n")
  cat("=====================================\n\n")
  
  # Process CGGA301 dataset
  cgga301_data <- process_cgga_dataset("CGGA301")
  cgga301_gbm <- filter_gbm_samples("CGGA301")
  
  # Process CGGA325 dataset
  cgga325_data <- process_cgga_dataset("CGGA325")
  cgga325_gbm <- filter_gbm_samples("CGGA325")
  
  # Generate summary report
  cat("Pipeline Summary Report:\n")
  cat("-----------------------\n")
  cat("CGGA301 - Total samples:", nrow(cgga301_data), "\n")
  cat("CGGA301 - GBM samples:", nrow(cgga301_gbm), "\n")
  cat("CGGA325 - Total samples:", nrow(cgga325_data), "\n")
  cat("CGGA325 - GBM samples:", nrow(cgga325_gbm), "\n")
  cat("Total GBM samples:", nrow(cgga301_gbm) + nrow(cgga325_gbm), "\n\n")
  
  # Save combined GBM dataset (optional)
  combined_gbm <- bind_rows(
    cgga301_gbm %>% mutate(Dataset = "CGGA301"),
    cgga325_gbm %>% mutate(Dataset = "CGGA325")
  )
  
  save(combined_gbm, file = "CGGA_combined_GBM.Rdata")
  cat("Combined GBM dataset saved to: CGGA_combined_GBM.Rdata\n")
  
  return(list(
    CGGA301 = cgga301_data,
    CGGA301_GBM = cgga301_gbm,
    CGGA325 = cgga325_data,
    CGGA325_GBM = cgga325_gbm,
    Combined_GBM = combined_gbm
  ))
}

# Execute pipeline
results <- process_cgga_pipeline()

# Save session information
writeLines(capture.output(sessionInfo()), "CGGA_processing_session_info.txt")
cat("\nSession information saved to: CGGA_processing_session_info.txt\n")
```

## English Documentation

# CGGA Data Processing Pipeline

## Overview
Processes Chinese Glioma Genome Atlas (CGGA) datasets (CGGA301 and CGGA325), integrating gene expression data with clinical survival information and filtering for glioblastoma (GBM) samples.

## Input Files Required
**CGGA301:**
- `CGGA.mRNA_array_301_clinical.20200506.txt` - Clinical data
- `CGGA.mRNA_array_301_gene_level.20200506.txt` - Gene expression data

**CGGA325:**
- `CGGA.mRNAseq_325_clinical.20200506.txt` - Clinical data
- `CGGA.mRNAseq_325.RSEM-genes.20200506.txt` - RSEM normalized counts
- `CGGA.mRNAseq_325.Read_Counts-genes.20220620.txt` - Raw read counts

## Key Features
1. **Automated Processing**: Unified functions for both CGGA datasets
2. **Data Integration**: Matches expression and clinical data by sample IDs
3. **Time Conversion**: Converts survival time from days to months
4. **GBM Filtering**: Identifies glioblastoma samples using histology information
5. **Normalization**: Applies log2 transformation to RNA-seq data (CGGA325)
6. **Quality Control**: Verifies sample alignment and filters invalid data

## Output Files
- `CGGA301_sur_exp.Rdata` - Processed CGGA301 dataset
- `CGGA301_sur_exp_all_GBM.Rdata` - GBM subset of CGGA301
- `CGGA325_sur_exp.Rdata` - Processed CGGA325 dataset
- `CGGA325_sur_exp_all_GBM.Rdata` - GBM subset of CGGA325
- `CGGA_combined_GBM.Rdata` - Combined GBM samples from both datasets

## Data Structure
Each processed RData file contains a dataframe with:
- First 3 columns: `ID`, `OS.time` (months), `OS` (0=alive, 1=dead)
- Remaining columns: Gene expression values

## Processing Steps
1. **Data Loading**: Load clinical and expression data
2. **Sample Matching**: Identify common samples between datasets
3. **Time Processing**: Convert survival time units
4. **Data Filtering**: Remove samples with missing or invalid data
5. **GBM Identification**: Filter samples based on histology
6. **Data Transformation**: Apply log2 normalization to RNA-seq data
7. **Output Generation**: Save processed datasets

## Usage
1. Ensure all required input files are in the working directory
2. Run the script: `source("CGGA_processing.R")`
3. The pipeline will automatically process both datasets and generate outputs

## Notes
- CGGA301 uses microarray data, CGGA325 uses RNA-seq data
- Survival time is converted from days to months (divided by 30)
- GBM samples are identified using histology field containing "GBM" or "glioblastoma"
- The script generates a comprehensive summary report of processing results

## Dependencies
- survival
- tidyverse
- data.table

## Quality Control
- Verifies sample ID alignment between datasets
- Reports processing statistics and sample counts
- Saves session information for reproducibility

#################################
# CGGA693 Data Processing Pipeline
# Chinese Glioma Genome Atlas dataset 693 (RNA-seq)
#################################
library(tidyverse)
library(data.table)

process_cgga693 <- function() {
  cat("Starting CGGA693 Data Processing\n")
  cat("================================\n")
  
  # Set working directory (adjust as needed)
  # setwd("./Download")  # Uncomment if data is in specific directory
  
  # Load clinical data
  cat("1. Loading clinical data...\n")
  clinical_data <- fread("CGGA.mRNAseq_693_clinical.20200506.txt", data.table = FALSE)
  
  # Load expression data - RSEM normalized counts
  cat("2. Loading expression data...\n")
  rsem_data <- fread("CGGA.mRNAseq_693.RSEM-genes.20200506.txt", data.table = FALSE) %>%
    column_to_rownames("Gene_Name")
  
  # Load raw count data (optional, for differential expression analysis)
  cat("3. Loading raw count data...\n")
  count_data <- fread("CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt", data.table = FALSE) %>%
    column_to_rownames("gene_name")
  
  # Verify sample consistency between datasets
  cat("4. Verifying dataset consistency...\n")
  clinical_samples <- clinical_data$CGGA_ID
  rsem_samples <- colnames(rsem_data)
  count_samples <- colnames(count_data)
  
  cat("   Clinical data samples:", length(clinical_samples), "\n")
  cat("   RSEM data samples:", length(rsem_samples), "\n")
  cat("   Count data samples:", length(count_samples), "\n")
  
  # Prepare expression matrix
  expr_matrix <- t(rsem_data) %>% as.data.frame()
  
  # Prepare survival metadata
  cat("5. Preparing survival data...\n")
  meta <- clinical_data %>%
    dplyr::select(CGGA_ID, OS.time = OS, OS = Censor..alive.0..dead.1.) %>%
    column_to_rownames("CGGA_ID")
  
  # Find common samples between expression and clinical data
  common_samples <- intersect(rownames(meta), rownames(expr_matrix))
  
  cat("   Common samples found:", length(common_samples), "\n")
  
  if (length(common_samples) == 0) {
    stop("No common samples found between clinical and expression data")
  }
  
  # Subset data to common samples
  meta <- meta[common_samples, , drop = FALSE]
  expr_matrix <- expr_matrix[common_samples, ]
  
  # Verify alignment
  if (!identical(rownames(meta), rownames(expr_matrix))) {
    stop("Sample alignment failed")
  }
  
  cat("   Sample alignment verified: TRUE\n")
  
  # Create survival-expression dataframe
  survival_info <- data.frame(
    ID = rownames(meta),
    OS.time = as.numeric(meta$OS.time),
    OS = meta$OS,
    row.names = rownames(meta)
  )
  
  # Combine survival info with expression data
  sur_exp <- cbind(survival_info, expr_matrix) %>%
    dplyr::select(ID, OS.time, OS, everything())
  
  # Process survival time (convert days to months)
  sur_exp$OS.time <- sur_exp$OS.time / 30
  
  # Filter samples
  cat("6. Filtering samples...\n")
  
  # Remove samples with missing data
  sur_exp <- sur_exp %>% na.omit()
  
  # Remove samples with non-positive survival time
  initial_samples <- nrow(sur_exp)
  sur_exp <- sur_exp %>% filter(OS.time > 0)
  filtered_samples <- initial_samples - nrow(sur_exp)
  
  cat("   Samples removed (OS.time <= 0):", filtered_samples, "\n")
  cat("   Remaining samples:", nrow(sur_exp), "\n")
  
  # Apply log2 transformation to expression data
  cat("7. Applying log2 transformation...\n")
  gene_start_col <- 4  # Columns after ID, OS.time, OS
  sur_exp[, gene_start_col:ncol(sur_exp)] <- log2(sur_exp[, gene_start_col:ncol(sur_exp)] + 1)
  
  # Check expression range
  expr_range <- range(sur_exp[, gene_start_col:ncol(sur_exp)], na.rm = TRUE)
  cat("   Expression range after log2:", round(expr_range[1], 2), "to", round(expr_range[2], 2), "\n")
  
  # Save processed data
  cat("8. Saving results...\n")
  
  # Save basic survival-expression matrix
  save(sur_exp, file = "CGGA693_sur_exp.Rdata")
  cat("   Saved: CGGA693_sur_exp.Rdata\n")
  
  # Save comprehensive dataset
  save(sur_exp, clinical_data, count_data, rsem_data, 
       file = "CGGA693_complete_dataset.Rdata")
  cat("   Saved: CGGA693_complete_dataset.Rdata\n")
  
  # Generate summary statistics
  cat("\nProcessing Summary:\n")
  cat("-------------------\n")
  cat("Total samples processed:", nrow(sur_exp), "\n")
  cat("Total genes:", ncol(sur_exp) - 3, "\n")
  cat("Survival status distribution:\n")
  print(table(sur_exp$OS))
  cat("Median survival time (months):", median(sur_exp$OS.time), "\n")
  cat("Mean survival time (months):", round(mean(sur_exp$OS.time), 1), "\n")
  cat("Survival time range (months):", round(range(sur_exp$OS.time), 1), "\n")
  
  # Optional: Save sample metadata
  sample_info <- data.frame(
    Sample_ID = sur_exp$ID,
    OS_status = sur_exp$OS,
    OS_time_months = round(sur_exp$OS.time, 1),
    stringsAsFactors = FALSE
  )
  
  write.csv(sample_info, "CGGA693_sample_info.csv", row.names = FALSE)
  cat("   Sample metadata saved: CGGA693_sample_info.csv\n")
  
  return(list(
    sur_exp = sur_exp,
    clinical = clinical_data,
    counts = count_data,
    rsem = rsem_data
  ))
}

# Function to filter GBM samples from CGGA693
filter_cgga693_gbm <- function() {
  cat("\nFiltering GBM samples from CGGA693...\n")
  
  # Load clinical data
  clinical_data <- fread("CGGA.mRNAseq_693_clinical.20200506.txt", data.table = FALSE)
  
  # Identify GBM samples
  gbm_samples <- clinical_data %>%
    filter(str_detect(toupper(Histology), "GBM|GLIOBLASTOMA"))
  
  gbm_ids <- gbm_samples$CGGA_ID
  
  cat("   GBM samples identified:", length(gbm_ids), "\n")
  
  # Load processed data
  load("CGGA693_sur_exp.Rdata")
  
  # Filter for GBM samples
  sur_exp_gbm <- sur_exp %>%
    filter(ID %in% gbm_ids) %>%
    na.omit()
  
  cat("   GBM samples with complete data:", nrow(sur_exp_gbm), "\n")
  
  # Filter expression matrices
  rsem_data <- fread("CGGA.mRNAseq_693.RSEM-genes.20200506.txt", data.table = FALSE) %>%
    column_to_rownames("Gene_Name")
  
  count_data <- fread("CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt", data.table = FALSE) %>%
    column_to_rownames("gene_name")
  
  rsem_gbm <- rsem_data[, gbm_ids]
  count_gbm <- count_data[, gbm_ids]
  
  # Save GBM subset
  save(sur_exp_gbm, clinical_data, gbm_samples, rsem_gbm, count_gbm,
       file = "CGGA693_GBM_dataset.Rdata")
  
  cat("   GBM dataset saved: CGGA693_GBM_dataset.Rdata\n")
  
  return(sur_exp_gbm)
}

# Main execution
if (interactive()) {
  # Process CGGA693 dataset
  cgga693_results <- process_cgga693()
  
  # Filter GBM samples
  cgga693_gbm <- filter_cgga693_gbm()
  
  cat("\nCGGA693 Processing Complete!\n")
  cat("=============================\n")
  cat("Output files generated:\n")
  cat("1. CGGA693_sur_exp.Rdata - Processed survival-expression matrix\n")
  cat("2. CGGA693_complete_dataset.Rdata - Comprehensive dataset\n")
  cat("3. CGGA693_sample_info.csv - Sample metadata\n")
  cat("4. CGGA693_GBM_dataset.Rdata - GBM subset\n")
}
```

## English Documentation

# CGGA693 Data Processing

## Overview
Processes the CGGA693 RNA-seq dataset from the Chinese Glioma Genome Atlas, integrating gene expression data with clinical survival information.

## Input Files Required
1. `CGGA.mRNAseq_693_clinical.20200506.txt` - Clinical and survival data
2. `CGGA.mRNAseq_693.RSEM-genes.20200506.txt` - RSEM normalized gene expression
3. `CGGA.mRNAseq_693.Read_Counts-genes.20220620.txt` - Raw read counts

## Key Processing Steps
1. **Data Loading**: Load clinical and expression datasets
2. **Sample Matching**: Identify common samples across datasets
3. **Data Integration**: Merge survival information with expression profiles
4. **Time Conversion**: Convert survival time from days to months
5. **Quality Filtering**: Remove samples with missing or invalid data
6. **Normalization**: Apply log2 transformation to expression values
7. **GBM Filtering**: Extract glioblastoma samples (optional)
8. **Output Generation**: Save processed datasets

## Output Files
- `CGGA693_sur_exp.Rdata` - Integrated survival-expression matrix
- `CGGA693_complete_dataset.Rdata` - All raw and processed data
- `CGGA693_sample_info.csv` - Sample metadata
- `CGGA693_GBM_dataset.Rdata` - GBM-specific subset

## Data Structure
The main output (`CGGA693_sur_exp.Rdata`) contains:
- `ID`: Sample identifier (CGGA_ID)
- `OS.time`: Overall survival time (months)
- `OS`: Survival status (0=alive, 1=dead)
- Columns 4+: Log2-transformed gene expression values

## Key Features
- **Comprehensive Processing**: Handles both RSEM normalized and raw count data
- **Quality Control**: Verifies sample alignment and filters invalid data
- **Data Transformation**: Applies log2(x+1) transformation to expression data
- **GBM Subsetting**: Optional extraction of glioblastoma samples
- **Detailed Logging**: Provides step-by-step processing information
- **Summary Statistics**: Reports key metrics about the processed dataset

## Usage
1. Place all required input files in the working directory
2. Run the script: `source("CGGA693_processing.R")`
3. The pipeline will automatically process the data and generate outputs
4. For interactive use, run: `process_cgga693()` and `filter_cgga693_gbm()`

## Data Specifications
- **Dataset**: CGGA693 mRNA sequencing
- **Technology**: RNA-seq
- **Normalization**: RSEM + log2 transformation
- **Survival Time**: Converted to months (days/30)
- **Sample Filtering**: Removes samples with OS.time ≤ 0

## Dependencies
- tidyverse
- data.table

## Notes
- The script includes both full dataset processing and GBM-specific filtering
- All expression values are log2-transformed (log2(x+1))
- Sample alignment is verified at multiple stages
- The pipeline generates comprehensive summary statistics