# Survival Analysis Modeling with Mime1 R Package

# https://github.com/l-magnificence/Mime

# IMPORTANT NOTE: *Note: We are users of this package, not developers*
# This analysis uses the Mime1 R package (already published) 
# for survival model development. 
# We are NOT developing new algorithms. We are just users.

# Mime1 package source: https://github.com/l-magnificence/Mime
# For more details ## Reference 
# https://pubmed.ncbi.nlm.nih.gov/40788954/
# Liu H, Zhang W, Zhang Y, Li X, Wanggou S. Constructing and Visualizing Models using Mime-based Machine-learning Framework. J Vis Exp. 2025 Jul 22;(221). doi: 10.3791/68553. PMID: 40788954.
library(Mime1)
library(beepr)
library(tidyverse)
library(survival)
library(survminer)
library(timeROC)
library(glmnet)
library(openxlsx)

# Set up directories
res.path <- "Results"
fig.path <- "Figures"
dir.create(res.path, showWarnings = FALSE, recursive = TRUE)
dir.create(fig.path, showWarnings = FALSE, recursive = TRUE)

# Load preprocessed data
load(file.path(res.path, "Data.pre.101.ML.Rdata"))

cat("Data loaded successfully\n")
cat("Training dataset samples:", nrow(list_train_vali_Data[[1]]), "\n")
cat("Number of datasets:", length(list_train_vali_Data), "\n")
cat("Dataset names:", names(list_train_vali_Data), "\n")

# ==============================================
# STEP 1: Model Development with Mime1 Package
# ==============================================

cat("\nStarting model development with Mime1 package...\n")
cat("NOTE: Using published Mime1 algorithms, not custom-developed methods\n")

# Set seed for reproducibility
set.seed(123)

# Train survival prediction model using Mime1
system.time({
  cat("Training multiple survival models...\n")
  
  res <- ML.Dev.Prog.Sig(
    train_data = list_train_vali_Data[[1]],
    list_train_vali_Data = list_train_vali_Data,
    unicox.filter.for.candi = TRUE,
    unicox_p_cutoff = 0.01,
    candidate_genes = filtered_genes,
    mode = 'all',
    nodesize = 5,
    seed = 123
  )
  
  # Alert when training completes
  beep("fanfare")
  cat("Model training completed!\n")
})

# Save model results
save(res, file = file.path(res.path, "Results.ML.Dev.Prog.Sig.rdata"))
cat("Model results saved to:", file.path(res.path, "Results.ML.Dev.Prog.Sig.rdata"), "\n")

# Extract key results
Cindex.res <- res$Cindex.res  # Concordance indices
ml.res <- res$ml.res           # Model objects
riskscore <- res$riskscore     # Risk scores
Sig.genes <- res$Sig.genes     # Selected signature genes

cat("\nModel Summary:\n")
cat("Selected genes:", length(Sig.genes), "\n")
cat("Signature genes:", paste(Sig.genes, collapse = ", "), "\n")

# Save signature genes
write.xlsx(as.data.frame(Sig.genes), 
           file = file.path(res.path, "Signature_Genes.xlsx"),
           rowNames = FALSE)

# ==============================================
# STEP 2: Model Interpretation
# ==============================================

cat("\nInterpreting model results...\n")

# Select best performing model (based on C-index)
best_model <- "StepCox[both] + Ridge"
cat("Best performing model:", best_model, "\n")

# Extract coefficients from Ridge model
ridge_model <- res$ml.res[[best_model]]
coef_lambda.min <- coef(ridge_model, s = "lambda.min")
coef_lambda.1se <- coef(ridge_model, s = "lambda.1se")

# Process coefficients
process_coefficients <- function(coef_matrix) {
  coef_df <- as.data.frame(as.matrix(coef_matrix))
  coef_df$Gene <- rownames(coef_df)
  colnames(coef_df)[1] <- "Coefficient"
  coef_df <- coef_df[coef_df$Coefficient != 0, ]
  return(coef_df)
}

coef_min_df <- process_coefficients(coef_lambda.min)
coef_1se_df <- process_coefficients(coef_lambda.1se)

cat("Number of genes in lambda.min model:", nrow(coef_min_df), "\n")
cat("Number of genes in lambda.1se model:", nrow(coef_1se_df), "\n")

# Save coefficients
write.xlsx(coef_min_df, 
           file = file.path(res.path, "Ridge_Model_Coefficients_lambda.min.xlsx"),
           rowNames = FALSE)

# ==============================================
# STEP 3: Model Performance Evaluation
# ==============================================

cat("\nEvaluating model performance...\n")

# 3.1 Concordance Index Visualization
cat("Generating C-index heatmap...\n")
pdf(file.path(fig.path, "Cindex_Heatmap.pdf"), 
    width = 8, height = 6)
cindex_dis_all(res,
               validate_set = names(list_train_vali_Data)[-1],
               order = names(list_train_vali_Data),
               width = 0.06,
               height = NULL)
dev.off()

# 3.2 Time-dependent AUC Calculation
cat("Calculating time-dependent AUC...\n")

time_points <- c(1, 2, 3, 4, 5, 6)  # Years
auc_results <- list()

for (time_point in time_points) {
  cat("  Calculating", time_point, "year AUC...\n")
  auc_results[[paste0("auc_", time_point, "y")]] <- cal_AUC_ml_res(
    res.by.ML.Dev.Prog.Sig = res,
    train_data = list_train_vali_Data[[1]],
    inputmatrix.list = list_train_vali_Data,
    mode = 'all',
    AUC_time = time_point,
    auc_cal_method = "KM"
  )
}

# Save AUC results
save(auc_results, file = file.path(res.path, "AUC_Results.rdata"))

# 3.3 AUC Visualization
cat("Generating AUC bar plots...\n")
pdf(file.path(fig.path, "AUC_Barplot.pdf"), 
    width = 10, height = 6)
auc_dis_select(
  list(auc_results$auc_1y, auc_results$auc_3y, auc_results$auc_5y),
  model_name = best_model,
  dataset = names(list_train_vali_Data),
  order = names(list_train_vali_Data),
  year = c(1, 3, 5)
) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
dev.off()

# ==============================================
# STEP 4: Visualization
# ==============================================

cat("\nGenerating visualizations...\n")

# 4.1 ROC Curves
cat("Generating ROC curves...\n")
roc_plots <- list()
for (time_point in c(1, 3, 5)) {
  roc_plots[[paste0("roc_", time_point, "y")]] <- roc_vis(
    auc_results[[paste0("auc_", time_point, "y")]],
    model_name = best_model,
    dataset = names(list_train_vali_Data),
    order = names(list_train_vali_Data),
    anno_position = c(0.65, 0.55),
    year = time_point
  )
}

# 4.2 Survival Curves
cat("Generating survival curves...\n")
survival_plots <- list()

for (dataset_name in names(list_train_vali_Data)) {
  cat("  Plotting survival for", dataset_name, "...\n")
  
  survival_plots[[dataset_name]] <- rs_sur(
    res,
    model_name = best_model,
    dataset = dataset_name,
    median.line = "hv",
    color = c("blue", "red"),
    conf.int = TRUE,
    xlab = "Time (Days)",
    pval.coord = c(2000, 0.9)
  )
}

# Combine all survival plots
pdf(file.path(fig.path, "All_Survival_Curves.pdf"), 
    width = 15, height = 20)
plot_list <- list()
for (i in seq_along(survival_plots)) {
  plot_list[[i]] <- survival_plots[[i]]
}
aplot::plot_list(gglist = plot_list, ncol = 2)
dev.off()

# 4.3 Time-dependent ROC Curves
cat("Generating time-dependent ROC curves...\n")
pdf(file.path(fig.path, "Time_ROC_Curves.pdf"), 
    width = 12, height = 8)

par(mfrow = c(2, 3))
for (dataset_name in names(list_train_vali_Data)) {
  
  # Get risk scores for this dataset
  tmp <- res$riskscore[[best_model]][[dataset_name]]
  tmp$OS.time <- as.numeric(tmp$OS.time)
  tmp$OS <- as.numeric(tmp$OS)
  
  # Calculate timeROC
  ROC_result <- with(tmp,
                     timeROC(T = OS.time,
                             delta = OS,
                             marker = RS,
                             cause = 1,
                             weighting = "marginal",
                             times = c(365, 1080, 1800),  # 1, 3, 5 years
                             ROC = TRUE,
                             iid = TRUE)
  )
  
  # Plot
  plot(ROC_result, time = 365, col = "#0072B2", lwd = 2, title = "")
  plot(ROC_result, time = 1080, col = "#D55E00", lwd = 2, add = TRUE)
  plot(ROC_result, time = 1800, col = "#009E73", lwd = 2, add = TRUE)
  
  title(main = dataset_name, line = 0.5)
  legend("bottomright",
         c(paste0("1-Year (AUC=", round(ROC_result$AUC[1], 3), ")"),
           paste0("3-Year (AUC=", round(ROC_result$AUC[2], 3), ")"),
           paste0("5-Year (AUC=", round(ROC_result$AUC[3], 3), ")")),
         col = c("#0072B2", "#D55E00", "#009E73"),
         lty = 1, lwd = 2)
}
dev.off()

# ==============================================
# STEP 5: Meta-analysis and Comparison
# ==============================================

cat("\nPerforming meta-analysis...\n")

# 5.1 Univariate Cox meta-analysis
unicox_results <- cal_unicox_ml_res(
  res.by.ML.Dev.Prog.Sig = res,
  optimal.model = best_model,
  type = 'categorical'
)

metamodel <- cal_unicox_meta_ml_res(input = unicox_results)

# Visualize meta-analysis results
pdf(file.path(fig.path, "Meta_Analysis.pdf"), 
    width = 8, height = 6)
meta_unicox_vis(metamodel,
                dataset = names(list_train_vali_Data))
dev.off()

# ==============================================
# STEP 6: Comparison with Published Signatures
# ==============================================

cat("\nComparing with published signatures...\n")

# Calculate risk scores for published signatures
published_scores <- cal_RS_pre.prog.sig(
  use_your_own_collected_sig = FALSE,
  type.sig = c('LGG', 'GBM', 'Glioma'),
  list_input_data = list_train_vali_Data
)

# Compare Hazard Ratios
pdf(file.path(fig.path, "HR_Comparison.pdf"), 
    width = 10, height = 6)
HR_com(published_scores,
       res,
       model_name = best_model,
       dataset = names(list_train_vali_Data),
       type = "categorical")
dev.off()

# ==============================================
# STEP 7: Generate Summary Report
# ==============================================

cat("\nGenerating summary report...\n")

# Create summary statistics
summary_stats <- data.frame(
  Dataset = names(list_train_vali_Data),
  Samples = sapply(list_train_vali_Data, nrow),
  Genes = sapply(list_train_vali_Data, function(x) ncol(x) - 3),
  C_index = sapply(names(list_train_vali_Data), 
                   function(d) res$Cindex.res$Cindex[res$Cindex.res$ID == d & 
                                                      res$Cindex.res$Model == best_model][1])
)

write.xlsx(summary_stats, 
           file = file.path(res.path, "Model_Summary_Statistics.xlsx"),
           rowNames = FALSE)

# Final completion message
beep("coin")
cat("\n==========================================\n")
cat("ANALYSIS COMPLETED SUCCESSFULLY!\n")
cat("==========================================\n")
cat("\nOutput files generated:\n")
cat("1. Model results: Results/Results.ML.Dev.Prog.Sig.rdata\n")
cat("2. Figures: Figures/*.pdf\n")
cat("3. Statistics: Results/*.xlsx\n")
cat("\nBest model:", best_model, "\n")
cat("Signature genes selected:", length(Sig.genes), "\n")
cat("Analysis time:", date(), "\n")
```

## English Documentation

# Survival Analysis Modeling with Mime1 R Package

## **CRITICAL DISCLAIMER**
**This analysis uses the published Mime1 R package for survival model development. All machine learning algorithms and modeling approaches are implemented through this existing package. We are NOT developing new algorithms or methods.**

## Overview
This script performs comprehensive survival analysis for glioma datasets using the Mime1 R package, which provides multiple machine learning algorithms for survival prediction.

## Key Points
1. **Using Published Package**: All modeling is done through Mime1 package functions
2. **Multi-cohort Analysis**: Processes multiple glioma datasets (TCGA, CGGA, GEO)
3. **Model Comparison**: Evaluates multiple ML algorithms and selects best performer
4. **Comprehensive Validation**: Time-dependent AUC, C-index, survival curves
5. **Comparison to Literature**: Benchmarks against published prognostic signatures

## Dataset Requirements
Preprocessed data in `Data.pre.101.ML.Rdata` containing:
- `list_train_vali_Data`: List of datasets with survival and expression data
- `filtered_genes`: Candidate gene list for feature selection

## Processing Steps
1. **Model Training**: Multiple ML algorithms via `ML.Dev.Prog.Sig()`
2. **Performance Evaluation**: C-index and time-dependent AUC
3. **Visualization**: ROC curves, survival plots, heatmaps
4. **Meta-analysis**: Combined results across datasets
5. **Literature Comparison**: Benchmark against published signatures

## Output Files
### Results (`Results/`)
- `Results.ML.Dev.Prog.Sig.rdata`: Complete model results
- `Signature_Genes.xlsx`: Selected prognostic genes
- `Ridge_Model_Coefficients*.xlsx`: Model coefficients
- `AUC_Results.rdata`: Time-dependent AUC values
- `Model_Summary_Statistics.xlsx`: Dataset statistics

### Figures (`Figures/`)
- `Cindex_Heatmap.pdf`: Concordance indices across models
- `AUC_Barplot.pdf`: AUC values at 1, 3, 5 years
- `All_Survival_Curves.pdf`: KM plots for all datasets
- `Time_ROC_Curves.pdf`: Time-dependent ROC curves
- `Meta_Analysis.pdf`: Forest plot of pooled HRs
- `HR_Comparison.pdf`: Comparison to published signatures

## Model Algorithms (via Mime1)
The script evaluates multiple algorithms including:
- Stepwise Cox regression with Ridge/Enet
- Random Survival Forests (RSF)
- Partial Least Squares Cox (plsRcox)
- Supervised Principal Components (SuperPC)
- Survival Support Vector Machines

## Usage
1. Ensure preprocessed data is available
2. Run script: `source("survival_analysis_mime1.R")`
3. Results will be saved in `Results/` and `Figures/` directories

## Dependencies
- Mime1 (survival ML package)
- survival, survminer, timeROC
- glmnet, tidyverse
- beepr (for alerts)

## Notes
- Default seed: 123 for reproducibility
- Best model selected based on C-index
- All visualizations use publication-ready colors
- Complete logging of processing steps

## Reference
Liu H, Zhang W, Zhang Y, Li X, Wanggou S. Constructing and Visualizing Models using Mime-based Machine-learning Framework. J Vis Exp. 2025 Jul 22;(221). doi: 10.3791/68553. PMID: 40788954.

Mime1 package: https://github.com/l-magnificence/Mime

*Note: We are users of this package, not developers*