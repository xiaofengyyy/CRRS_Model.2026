# TCGA-LGG Data Processing Pipeline
## Overview
This script downloads and processes TCGA-LGG (Lower Grade Glioma) data from the Xena platform, performing essential preprocessing steps for downstream analysis.

## Key Features
- **Data Download**: Automatically retrieves count data, TPM data, clinical information, and survival data
- **ID Conversion**: Converts Ensembl gene IDs to standardized gene symbols
- **Quality Filtering**: Removes lowly-expressed genes (expressed in <50% of samples)
- **Sample Annotation**: Classifies samples as tumor or normal based on TCGA barcodes
- **Data Output**: Saves processed matrices in RData format for easy loading

## Requirements

### R Packages
- **CRAN**: `tidyverse`, `survival`, `survminer`, `glmnet`, `timeROC`
- **Bioconductor**: `limma`, `DESeq2`, `edgeR`, `SummarizedExperiment`, `clusterProfiler`, `org.Hs.eg.db`
- **GitHub**: `AnnoProbe` (for gene annotation)

## Files Generated
The script creates the following files in `TCGA_LGG_Data/` directory:
  1. `TCGA-LGG_phenotype.tsv.gz` - Clinical metadata
2. `TCGA-LGG_survival.tsv` - Survival data
3. `TCGA-LGG_counts.tsv.gz` - Raw count matrix
4. `TCGA-LGG_tpm.tsv.gz` - TPM normalized matrix
5. `TCGA-LGG_processed.Rdata` - Final processed data (main output)

## Output Data Structure
The RData file contains:
  - `exp_counts_filtered`: Filtered count matrix (genes × samples)
- `exp_tpm_filtered`: Filtered TPM matrix (genes × samples)
- `Group`: Factor indicating tumor/normal status
- `clinical`: Clinical metadata dataframe
- `surv`: Survival data dataframe
- `proj`: Project identifier ("TCGA-LGG")

## Usage
1. Run the script: `source("TCGA_LGG_processing.R")`
2. Load processed data: `load("TCGA-LGG_processed.Rdata")`
3. The environment will contain all processed objects ready for analysis

## Data Filtering Criteria
- **Gene filtering**: Retain genes expressed (count > 0) in >50% of samples
- **Sample classification**: 
  - Tumor: TCGA barcode positions 14-15 = 01-09
- Normal: TCGA barcode positions 14-15 = 10-19

## Note
This pipeline focuses on protein-coding genes only. Non-coding RNAs are excluded from the final matrices.
###########################################
# Load required packages
options("repos" = c(CRAN = "http://mirrors.tuna.tsinghua.edu.cn/CRAN/"))

if (!require("BiocManager")) install.packages("BiocManager", update = FALSE, ask = FALSE)
options(BioC_mirror = "http://mirrors.tuna.tsinghua.edu.cn/bioconductor/")

# Define packages
cran_packages <- c('tidyverse', 'survival', 'survminer', 'glmnet', 'timeROC')
bioconductor_packages <- c("limma", "DESeq2", "edgeR", "SummarizedExperiment", "clusterProfiler", "org.Hs.eg.db")

# Install and load CRAN packages
for (pkg in cran_packages) {
  if (!require(pkg, character.only = TRUE)) {
    install.packages(pkg, ask = FALSE, update = FALSE)
    require(pkg, character.only = TRUE)
  }
}

# Install and load Bioconductor packages
for (pkg in bioconductor_packages) {
  if (!require(pkg, character.only = TRUE)) {
    BiocManager::install(pkg, ask = FALSE, update = FALSE)
    require(pkg, character.only = TRUE)
  }
}

# Install AnnoProbe from GitHub if needed
if (!require("AnnoProbe")) {
  devtools::install_github("jmzeng1314/AnnoProbe")
}
library(AnnoProbe)

# Set project and create directory
proj <- "TCGA-LGG"
dir.create("TCGA_LGG_Data", showWarnings = FALSE)
setwd("TCGA_LGG_Data")

# Load clinical and survival data
clinical <- read.delim(paste0(proj, "_phenotype.tsv.gz"), fill = TRUE, header = TRUE, sep = "\t")
surv <- read.delim(paste0(proj, "_survival.tsv"), header = TRUE)

# Process count data
count_data <- read.table(paste0(proj, "_counts.tsv.gz"), check.names = FALSE, row.names = 1, header = TRUE)
exp_counts <- as.matrix(2^count_data - 1)
exp_counts <- round(exp_counts)

# Process TPM data
tpm_data <- read.table(paste0(proj, "_tpm.tsv.gz"), check.names = FALSE, row.names = 1, header = TRUE)
exp_tpm <- as.matrix(2^tpm_data - 1)

# Extract Ensembl IDs (remove version numbers)
rownames(exp_counts) <- str_split(rownames(exp_counts), "\\.", simplify = TRUE)[, 1]
rownames(exp_tpm) <- str_split(rownames(exp_tpm), "\\.", simplify = TRUE)[, 1]

# Annotate genes and convert Ensembl IDs to gene symbols
gene_anno <- annoGene(rownames(exp_counts), ID_type = "ENSEMBL")
protein_coding_genes <- gene_anno[gene_anno$biotypes == "protein_coding", ]

# Convert count matrix to gene symbols
exp_counts_symbol <- trans_array(exp_counts, ids = protein_coding_genes, from = "ENSEMBL", to = "SYMBOL")

# Convert TPM matrix to gene symbols
exp_tpm_symbol <- trans_array(exp_tpm, ids = protein_coding_genes, from = "ENSEMBL", to = "SYMBOL")

# Filter genes: keep genes expressed in more than 50% of samples
filter_genes <- function(exp_matrix) {
  exp_matrix[apply(exp_matrix, 1, function(x) sum(x > 0) > 0.5 * ncol(exp_matrix)), ]
}

exp_counts_filtered <- filter_genes(exp_counts_symbol)
exp_tpm_filtered <- filter_genes(exp_tpm_symbol)

cat("Count data dimensions after filtering:", dim(exp_counts_filtered), "\n")
cat("TPM data dimensions after filtering:", dim(exp_tpm_filtered), "\n")

# Create sample groups based on TCGA barcode (01-09: tumor, 10-19: normal)
create_sample_groups <- function(sample_names) {
  sample_type <- as.numeric(str_sub(sample_names, 14, 15))
  group <- ifelse(sample_type < 10, 'tumor', 'normal')
  factor(group, levels = c("normal", "tumor"))
}

Group <- create_sample_groups(colnames(exp_counts_filtered))
table(Group)

# Save processed data
save(
  exp_counts_filtered,
  exp_tpm_filtered,
  Group,
  proj,
  clinical,
  surv,
  file = paste0(proj, "_processed.Rdata")
)

cat("Data processing completed. Files saved to:", getwd(), "\n")

