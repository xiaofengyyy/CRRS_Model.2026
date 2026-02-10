# Project Directory Setup
# Sets up standardized directory structure for analysis projects
setup_project_directories <- function(work_path = ".") {
  # Set working directory
  setwd(work_path)
  
  # Define directory paths
  code_path <- file.path(work_path, "Codes")
  data_path <- file.path(work_path, "InputData")
  results_path <- file.path(work_path, "Results")
  figures_path <- file.path(work_path, "Figures")
  
  # Create directories if they don't exist
  create_directory <- function(path) {
    if (!dir.exists(path)) {
      dir.create(path, recursive = TRUE, showWarnings = FALSE)
      cat("Created directory:", path, "\n")
    } else {
      cat("Directory already exists:", path, "\n")
    }
  }
  
  # Create all directories
  cat("Setting up project directory structure...\n")
  cat("=========================================\n")
  
  create_directory(code_path)
  create_directory(data_path)
  create_directory(results_path)
  create_directory(figures_path)
  
  # Set global options
  options(stringsAsFactors = FALSE)
  
  # Load utility functions if available
  lib_file <- "scRNA_scripts/lib.R"
  if (file.exists(lib_file)) {
    cat("Loading utility functions from:", lib_file, "\n")
    source(lib_file)
  } else {
    cat("Utility file not found:", lib_file, "\n")
  }
  
  # Return directory paths as list
  directory_list <- list(
    work_path = work_path,
    code_path = code_path,
    data_path = data_path,
    results_path = results_path,
    figures_path = figures_path
  )
  
  cat("\nProject setup complete!\n")
  cat("Directory structure:\n")
  cat("-", work_path, "/ (working directory)\n")
  cat("  |- Codes/       (analysis scripts)\n")
  cat("  |- InputData/   (raw and processed data)\n")
  cat("  |- Results/     (analysis outputs)\n")
  cat("  |- Figures/     (generated plots)\n")
  
  return(directory_list)
}

# Example usage (uncomment to use):
# directories <- setup_project_directories(".")
```

## English Documentation

# Project Directory Setup

## Overview
Sets up standardized directory structure for reproducible analysis projects.

## Directory Structure
Creates the following directory structure:
  ```
Working Directory/
  ├── Codes/         # Analysis scripts and source code
  ├── InputData/     # Raw and processed input data
  ├── Results/       # Analysis results and outputs
  ├── Figures/       # Generated plots and visualizations
  ```

## Usage
```r
# Basic setup in current directory
directories <- setup_project_directories(".")

# Setup in specific directory
# directories <- setup_project_directories("/path/to/your/project")
```

## Features
- **Automatic Creation**: Creates directories if they don't exist
- **Path Management**: Returns organized path list for reference
- **Utility Loading**: Attempts to load shared utility functions
- **Configuration**: Sets default R options for consistency

## Returns
List containing paths to all created directories for easy reference in scripts.

## Notes
1. Run this script at the beginning of any analysis project
2. Keep raw data in `InputData/`
3. Store scripts in `Codes/`
4. Save outputs to `Results/` and `Figures/`