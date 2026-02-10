#Code for GSE16011 data processing:

# GSE16011 Data Processing Pipeline
# GEO Accession: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE16011

library(tidyverse)
library(data.table)
library(openxlsx)

# Load GSE16011 expression data
load("GSE16011_data.RData")  # Contains expression matrix with gene symbols
exp <- GSE16011[[2]]  # Expression matrix

# Load survival data
os.sf <- fread("brain.cancer.glioma.gse16011.hgu133plus2_entrezcdf.tsv", data.table = FALSE)
os.sf2 <- fread("PRECOG.GSE16011.HGU133Plus2_EntrezCDF.OS_base", data.table = FALSE)

# Load clinical annotations for GBM filtering
cli <- read.xlsx("GSE16011_sample_anno.xlsx")

# Filter for GBM (glioblastoma) samples only
gbm_samples <- cli %>%
  filter(str_detect(`histology:ch1`, "GBM|glioblastoma")) %>%
  pull(Array)  # Assuming 'Array' column contains sample IDs

# Prepare survival data
pd <- os.sf
surdata <- pd %>%
  dplyr::select(Array, OS.time = `overall survival (months):ch1`, OS = `survival status:ch1`) %>%
  dplyr::rename(ID = Array)

# Filter samples with valid survival time
surdata <- surdata %>%
  filter(OS.time > 0)

# Convert OS status to binary (1 = Dead, 0 = Alive)
surdata$OS <- ifelse(surdata$OS == "D", 1, 0)

# Convert survival time from months to days (if needed, multiplied by 30.44)
# surdata$OS.time <- surdata$OS.time * 30.44  # For days
surdata$OS.time <- surdata$OS.time * 12  # For months to years (12 months/year)

# Prepare expression matrix
exp_df <- as.data.frame(t(exp))
exp_df$ID <- rownames(exp_df)

# Filter for GBM samples only
exp_gbm <- exp_df %>%
  filter(ID %in% gbm_samples)

surdata_gbm <- surdata %>%
  filter(ID %in% gbm_samples)

# Find common samples between expression and survival data
common_samples <- intersect(exp_gbm$ID, surdata_gbm$ID)

# Subset data to common samples
exp_common <- exp_gbm %>%
  filter(ID %in% common_samples) %>%
  column_to_rownames("ID")

surdata_common <- surdata_gbm %>%
  filter(ID %in% common_samples) %>%
  column_to_rownames("ID")

# Verify sample alignment
identical(rownames(exp_common), rownames(surdata_common))

# Combine survival and expression data
sur_exp <- cbind(surdata_common, exp_common)

# Save processed data
save(sur_exp, file = "GSE16011_GBM_processed.Rdata")

# Summary statistics
cat("GSE16011 Processing Summary:\n")
cat("Total samples in expression data:", nrow(exp_df), "\n")
cat("GBM samples identified:", length(gbm_samples), "\n")
cat("Samples with survival data:", nrow(surdata), "\n")
cat("Final processed GBM samples:", nrow(sur_exp), "\n")
cat("Number of genes:", ncol(exp_common), "\n")

# Optional: Save sample metadata
sample_info <- data.frame(
  Sample_ID = rownames(sur_exp),
  OS_status = sur_exp$OS,
  OS_time = sur_exp$OS.time,
  Histology = "GBM"
)

write.csv(sample_info, file = "GSE16011_GBM_sample_info.csv", row.names = FALSE)
cat("\nSample metadata saved to: GSE16011_GBM_sample_info.csv\n")
```

## Documentation for GSE16011 Glioma Data Processing

## Overview
Processes GSE16011 microarray data from GEO for glioblastoma (GBM) samples, integrating expression profiles with clinical survival information.

## Input Files Required
1. `GSE16011_data.RData` - Expression matrix with gene symbols
2. `brain.cancer.glioma.gse16011.hgu133plus2_entrezcdf.tsv` - Survival data
3. `PRECOG.GSE16011.HGU133Plus2_EntrezCDF.OS_base` - Additional survival data (optional)
4. `GSE16011_sample_anno.xlsx` - Clinical annotations

## Processing Steps
1. **Data Loading**: Load expression matrix and clinical data
2. **Sample Filtering**: Select only GBM (glioblastoma) samples using histology annotations
3. **Survival Data Preparation**:
   - Convert OS status to binary (1=Dead, 0=Alive)
   - Filter samples with valid survival time
   - Convert time units as needed
4. **Data Integration**: Match expression and survival data by sample IDs
5. **Quality Control**: Verify sample alignment between datasets
6. **Output Generation**: Save integrated dataset for downstream analysis

## Output Files
- `GSE16011_GBM_processed.Rdata` - Integrated survival + expression matrix
- `GSE16011_GBM_sample_info.csv` - Sample metadata

## Data Structure
The processed RData file contains a dataframe `sur_exp` with:
- First 3 columns: `ID`, `OS.time`, `OS`
- Remaining columns: Gene expression values (gene symbols as column names)

## Key Features
- Focuses specifically on GBM samples
- Handles survival time unit conversions
- Ensures proper sample matching between datasets
- Provides summary statistics of processing results

## Dependencies
- tidyverse
- data.table
- openxlsx

## Usage
1. Place all input files in working directory
2. Run script: `source("GSE16011_processing.R")`
3. Load processed data: `load("GSE16011_GBM_processed.Rdata")`

## Notes
- Survival time is converted to years (multiplied by 12)
- Only samples with positive survival time are retained
- Sample IDs must match between expression and clinical data

#################################
# GSE43378 Data Processing Pipeline
# GEO Accession: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE43378
#################################
library(tidyverse)

# Load GSE43378 expression and clinical data
load("GSE43378_data.RData")  # Contains expression matrix with gene symbols
exp <- GSE43378[[2]]  # Expression matrix
pd <- GSE43378[[3]]   # Phenotype/clinical data

# Explore clinical variables
cat("Available clinical variables:\n")
print(colnames(pd))

cat("\nSample sources:\n")
print(table(pd$source_name_ch1))

cat("\nHistology distribution:\n")
print(table(pd$`histology:ch1`))

cat("\nPathological diagnosis:\n")
print(table(pd$`pathological diagnosis:ch1`))

# Prepare survival data
surdata <- pd %>%
  dplyr::select(
    ID = geo_accession,
    OS.time = `overall survival (months):ch1`,
    diagnosis = `pathological diagnosis:ch1`,
    histology = `histology:ch1`,
    grade = `who grade:ch1`,
    outcome = `outcome:ch1`
  )

# Convert OS status to binary (1 = Dead, 0 = Alive)
surdata$OS <- ifelse(surdata$outcome == "DEAD", 1, 0)

# Filter for glioma samples only (adjust criteria as needed)
glioma_samples <- surdata %>%
  filter(
    str_detect(tolower(diagnosis), "glioma|glioblastoma|astrocytoma|oligodendroglioma") |
    str_detect(tolower(histology), "glioma|glioblastoma|astrocytoma|oligodendroglioma")
  )

cat("\nGlioma samples identified:", nrow(glioma_samples), "\n")
cat("Diagnosis distribution in glioma samples:\n")
print(table(glioma_samples$diagnosis))

# Select relevant columns for final dataset
surdata_clean <- glioma_samples %>%
  dplyr::select(ID, OS.time, OS, diagnosis, histology, grade)

# Prepare expression matrix
exp_df <- as.data.frame(t(exp))
exp_df$ID <- rownames(exp_df)

# Filter for glioma samples only
exp_glioma <- exp_df %>%
  filter(ID %in% surdata_clean$ID)

# Find common samples between expression and survival data
common_samples <- intersect(exp_glioma$ID, surdata_clean$ID)

# Subset data to common samples
exp_common <- exp_glioma %>%
  filter(ID %in% common_samples) %>%
  column_to_rownames("ID")

surdata_common <- surdata_clean %>%
  filter(ID %in% common_samples) %>%
  column_to_rownames("ID")

# Order samples consistently
exp_common <- exp_common[rownames(surdata_common), ]

# Verify sample alignment
if (identical(rownames(exp_common), rownames(surdata_common))) {
  cat("\nSample alignment verified: TRUE\n")
} else {
  cat("\nWarning: Sample IDs do not match perfectly\n")
  # Reorder to ensure alignment
  exp_common <- exp_common[rownames(surdata_common), ]
}

# Process survival time (convert days to months)
surdata_common$OS.time <- as.numeric(surdata_common$OS.time)
surdata_common$OS.time <- surdata_common$OS.time / 30  # Convert days to months

# Filter samples with valid survival time
surdata_common <- surdata_common %>%
  filter(OS.time > 0)

# Update expression matrix to match filtered samples
exp_common <- exp_common[rownames(surdata_common), ]

# Combine survival and expression data
sur_exp <- cbind(surdata_common, exp_common)

# Save processed data
save(sur_exp, file = "GSE43378_glioma_processed.Rdata")

# Summary statistics
cat("\nGSE43378 Processing Summary:\n")
cat("Total samples in expression data:", nrow(exp_df), "\n")
cat("Glioma samples identified:", nrow(glioma_samples), "\n")
cat("Final processed samples:", nrow(sur_exp), "\n")
cat("Number of genes:", ncol(exp_common), "\n")
cat("Survival status distribution:\n")
print(table(sur_exp$OS))
cat("Median survival time (months):", median(sur_exp$OS.time, na.rm = TRUE), "\n")

# Optional: Save sample metadata with glioma subtypes
sample_info <- data.frame(
  Sample_ID = rownames(sur_exp),
  OS_status = sur_exp$OS,
  OS_time_months = round(sur_exp$OS.time, 2),
  Diagnosis = sur_exp$diagnosis,
  Histology = sur_exp$histology,
  Grade = sur_exp$grade,
  stringsAsFactors = FALSE
)

write.csv(sample_info, file = "GSE43378_glioma_sample_info.csv", row.names = FALSE)
cat("\nSample metadata saved to: GSE43378_glioma_sample_info.csv\n")

# Optional: Create glioma subtype-specific datasets
if (FALSE) {  # Set to TRUE to generate subtype files
  # Glioblastoma samples
  gbm_samples <- sur_exp %>%
    filter(str_detect(tolower(diagnosis), "glioblastoma|gbm")) %>%
    rownames()
  
  if (length(gbm_samples) > 0) {
    sur_exp_gbm <- sur_exp[gbm_samples, ]
    save(sur_exp_gbm, file = "GSE43378_glioblastoma_processed.Rdata")
    cat("Glioblastoma subset saved:", length(gbm_samples), "samples\n")
  }
  
  # Low-grade glioma samples
  lgg_samples <- sur_exp %>%
    filter(str_detect(tolower(grade), "ii|grade 2")) %>%
    rownames()
  
  if (length(lgg_samples) > 0) {
    sur_exp_lgg <- sur_exp[lgg_samples, ]
    save(sur_exp_lgg, file = "GSE43378_lgg_processed.Rdata")
    cat("LGG subset saved:", length(lgg_samples), "samples\n")
  }
}


## English Documentation
# GSE43378 Glioma Data Processing

## Overview
Processes GSE43378 microarray data from GEO for glioma samples, integrating gene expression profiles with detailed clinical and survival information.

## Input Files Required
1. `GSE43378_data.RData` - Contains expression matrix and phenotype data

## Processing Steps
1. **Data Loading**: Load expression matrix and clinical annotations
2. **Clinical Data Exploration**: Examine available clinical variables and distributions
3. **Survival Data Preparation**:
   - Select relevant clinical variables
   - Convert OS status to binary (1=Dead, 0=Alive)
   - Convert survival time from days to months
4. **Sample Filtering**: Select glioma samples based on diagnosis/histology
5. **Data Integration**: Match expression and clinical data by sample IDs
6. **Quality Control**: Verify sample alignment and filter invalid data
7. **Output Generation**: Save integrated dataset and sample metadata

## Output Files
- `GSE43378_glioma_processed.Rdata` - Integrated clinical + expression matrix
- `GSE43378_glioma_sample_info.csv` - Detailed sample metadata

## Data Structure
The processed RData file contains a dataframe `sur_exp` with:
- First 6 columns: Clinical variables (`ID`, `OS.time`, `OS`, `diagnosis`, `histology`, `grade`)
- Remaining columns: Gene expression values

## Clinical Variables Included
- `OS.time`: Overall survival time (months)
- `OS`: Survival status (0=Alive, 1=Dead)
- `diagnosis`: Pathological diagnosis
- `histology`: Histological classification
- `grade`: WHO tumor grade

## Key Features
- Automatic glioma sample identification based on diagnosis/histology
- Survival time unit conversion (days to months)
- Comprehensive quality checks and alignment verification
- Detailed processing summary with statistics
- Optional subtype-specific dataset generation

## Dependencies
- tidyverse

## Usage
1. Ensure `GSE43378_data.RData` is in working directory
2. Run script: `source("GSE43378_processing.R")`
3. Load processed data: `load("GSE43378_glioma_processed.Rdata")`

## Notes
- Survival time is converted from days to months (divided by 30)
- Only samples with positive survival time are retained
- Sample filtering criteria can be adjusted based on specific research needs
- The script includes optional code for generating glioma subtype-specific datasets

## Quality Control
- Verifies sample ID alignment between expression and clinical data
- Provides distribution statistics for clinical variables
- Reports processing summary with sample counts and survival statistics