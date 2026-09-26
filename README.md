# Single-Cell RNA-seq Analysis of Hepatocellular Carcinoma

## Overview

This project analyzes publicly available single-cell RNA-seq data from hepatocellular carcinoma (HCC), based on the study:

> Lu et al. (2022). *A single-cell atlas of the multicellular ecosystem of primary and metastatic hepatocellular carcinoma*. Nature Communications, 13, 4594.

**Dataset:** GSE149614

## Objective

The aim of this project is to analyze HCC single-cell RNA-seq data and identify the major cellular populations present in the tumor microenvironment.

## Analysis workflow

The analysis was performed in R using Seurat and includes:

1. Loading the dataset
2. Quality control
3. Normalization
4. Identification of highly variable genes
5. Scaling and PCA
6. Clustering
7. UMAP dimensionality reduction
8. Marker-gene identification
9. Cell-type annotation

## Dataset

The original study generated single-cell RNA-seq data from HCC patients and profiled 71,915 cells across multiple tissue sites.

For this project, a subset of the dataset was analyzed because of available computational resources.

The processed count matrix and metadata were obtained from the Gene Expression Omnibus (GEO).

**GEO accession:** GSE149614

## Cell populations

The analysis identified populations corresponding to:

* T cells and regulatory T cells
* NK and cytotoxic lymphocytes
* Macrophages
* Endothelial cells
* Fibroblasts / CAF-like cells
* Smooth muscle / pericyte-like cells
* Hepatocytes
* Cycling / proliferating cells

## Repository structure

```text
HCC-scRNAseq-GSE149614/
│
├── README.md
├── .gitignore
├── HCC-scRNAseq-GSE149614.Rproj
│
├── R/
│   ├── 01_full_pipeline.R
│   ├── loading data.R
│   ├── QC.R
│   └── normalization to annotation.R
│
├── data/
│   └── raw/
│
└── results/
```

## Software

* R
* Seurat
* ggplot2
* dplyr

## Note

This project is a reproducible learning and portfolio analysis based on the publicly available processed data from the original study. It does not attempt to reproduce every computational and experimental result reported in the publication.

Cell-type annotations were assigned primarily using canonical marker genes and should be considered working annotations.

## Reference

Lu, Y. et al. (2022). A single-cell atlas of the multicellular ecosystem of primary and metastatic hepatocellular carcinoma. *Nature Communications*, 13, 4594.

https://doi.org/10.1038/s41467-022-32283-3

## Data source

GEO: https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE149614
