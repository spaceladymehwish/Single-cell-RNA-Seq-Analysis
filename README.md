Seurat Single-Cell RNA-seq Analysis Pipeline
This repository provides a complete, modular pipeline for analyzing single-cell RNA sequencing (scRNA-seq) data using the Seurat package in R. The pipeline includes quality control, normalization, dimensionality reduction, clustering, differential expression, and visualization, tailored for exploratory and publication-ready analysis.

🔧 Requirements
R (≥ 4.0)
Seurat (v4+)
Tidyverse
patchwork
Matrix
cowplot
Install dependencies in R:
install.packages("Seurat")
install.packages("tidyverse")
install.packages("patchwork")
install.packages("cowplot")
Or use an R environment manager like renv or Conda (environment.yml provided if applicable).

🚀 Usage
Clone the repository
cd 
Load your data
Place your raw/filtered matrix (e.g., 10x filtered_feature_bc_matrix/) in the data/ folder.

Run the pipeline step-by-step
Execute each R script in order or combine them in an RMarkdown notebook.

Example in R console:
source("scripts/01_qc.R")
source("scripts/02_normalization.R")
...
Outputs
Plots (UMAP, feature plots, violin plots)
Tables (marker genes, QC stats)
Final annotated Seurat object for downstream use
