
---

# Spatial Transcriptomics Analysis of Human Breast Cancer (10x Visium)

**Author:** Varun Umesh Gowda   
**Platform:** 10x Genomics Visium   
**Species:** Human  
**Tissue:** Invasive Ductal Carcinoma (IDC), Fresh Frozen   
**Status:** Actively developed (core analyses validated)    

---

## Overview

This project implements a **reproducible spatial transcriptomics pipeline** for analyzing **10x Genomics Visium breast cancer data**, with a focus on:

* Spatial clustering and localization
* Data-driven and literature-supported cell type annotation
* Identification and validation of **TLS-like (Tertiary Lymphoid Structure–like) lymphoid aggregates**
* Spatial autocorrelation and neighborhood enrichment analyses

The pipeline emphasizes **biological validation**, **transparent assumptions**, and **interpretable spatial results**, rather than black-box annotation.

---

## Dataset

* **Source:** 10x Genomics / BioIVT
* **Sample:** Human Invasive Ductal Carcinoma (IDC)
* **Preparation:** Fresh frozen tissue, cryosectioned following Visium Spatial Protocols
* **Data Type:** Spot-level spatial gene expression with histology-aligned coordinates

**Data access:**
- 10x Genomics Visium Spatial Gene Expression datasets  
- https://www.10xgenomics.com/datasets/human-breast-cancer-visium-fresh-frozen-whole-transcriptome-1-standard

Due to file size constraints, raw Visium output files are not included in this repository.
Users should download the dataset directly from 10x Genomics and provide the path to the
expression matrix when running the pipeline.

---

## Project Goals

1. Learn and implement a full **spatial transcriptomics analysis workflow**
2. Go beyond clustering to **validate biological meaning using spatial evidence**
3. Identify and characterize **tumor-associated immune niches**, including TLS-like structures
4. Build a **portfolio-ready, modular, and extensible pipeline**

---

## Pipeline Structure (Current Version)

### **Part 1 – Spatial Validation & Annotation (This Repository)** ✅

Implemented and validated in this version:

* Spatial clustering (Leiden)
* Marker-based and mean-expression–based annotation
* Module scoring for immune and epithelial programs
* Spatial localization plots (tissue vs dots)
* Moran’s I spatial autocorrelation (targeted genes)
* Spatial neighborhood enrichment (cluster adjacency)
* Quantitative summaries for cluster-level interpretation

### **Part 2 – Downstream Analyses (Planned / In Progress)** ⏳

Code exists but results are under review:

* Pathway enrichment (MSigDB-based)
* Spatial domain composition
* Reference-based deconvolution

These components will be added in a future release after biological verification.

---

## Repository Contents

```
.
├── SP_claude.py
│   └─ Core preprocessing, clustering, annotation, and spatial analysis logic
│
├── spatial_plots.py
│   └─ Helper utilities for spatial visualization, module scoring,
│      Moran’s I, neighborhood enrichment, and summary exports
│
├── environment.yml
│   └─ Conda environment for reproducibility
│
├── README.md
│   └─ Project documentation (this file)
```

---

## Key Analyses and Outputs

### 1. Spatial Clustering & Localization

* Leiden clustering on spatial transcriptomics data
* Tissue-aligned visualization of clusters
* Cluster-specific spatial localization plots

**Outputs**

* `spatial_clusters_overview.png`
* `spatial_cluster_0_localization.png`

---

### 2. Cell Type & Niche Annotation

* Canonical immune and epithelial marker panels
* Mean-expression–based evaluation across clusters
* Manual curation supported by literature

**Key annotated programs**

* B cells
* T cells
* Plasma cells
* Epithelial cells
* TLS-associated chemokines

---

### 3. TLS-like Structure Identification

TLS-like lymphoid aggregates were identified based on:

* Co-expression of **CCL19, CCL21, CXCL13, LTB**
* Enrichment of B, T, and plasma cell markers
* Spatial co-localization on tissue
* Significant spatial autocorrelation (Moran’s I)
* Neighborhood enrichment patterns

> **Note:** Clusters are labeled as **“TLS-like”** unless canonical TLS histological organization is confirmed.

**Outputs**

* `spatial_module_TLS_enrichment.png`
* `moranI_tls_genes.csv`
* `cluster_0_score_summary.csv`

---

### 4. Spatial Autocorrelation (Moran’s I)

* Targeted Moran’s I analysis for immune and TLS-related genes
* Identifies genes with non-random spatial expression patterns

---

### 5. Neighborhood Enrichment Analysis

* Evaluates spatial adjacency between clusters
* Reveals immune–tumor and immune–immune interactions

**Outputs**

* `nhood_enrichment.png`
* `nhood_enrichment_zscores.csv` (if supported by Squidpy version)

---

## How to Run

### 1. Create Conda Environment

```bash
conda env create -f environment.yml
conda activate spatial
```

### 2. Spatial Preprocessing & Clustering

**Purpose:**
Preprocess raw 10x Visium data, perform clustering, and generate a processed AnnData object for downstream spatial validation.

**Key output:**
`adata_spatial_processed.h5ad`

```bash
python spatial_preprocessing_and_annotation.py \
  --data_path path/to/filtered_feature_bc_matrix.h5 \
  --out_dir results/
```



### **3. Spatial Validation & Visualization**

**Purpose:**
Validate cluster annotations using spatial localization, module enrichment, Moran’s I, and neighborhood enrichment.

**Key outputs:**

* Spatial gene & module localization plots
* Moran’s I tables
* Neighborhood enrichment heatmaps

```bash
python spatial_validation_and_visualization.py \
  --adata_path results/adata_spatial_processed.h5ad \
  --outdir results/spatial_validation \
  --cluster_key leiden \
  --target_cluster 0
```
---

## Software & Tools

* **Scanpy**
* **Squidpy**
* **NumPy / Pandas**
* **Matplotlib**
* **Python 3.10+**
* **Conda**

---

## Design Philosophy

* Modular scripts (core logic vs plotting helpers)
* Explicit biological assumptions
* Clear separation between **validated** and **experimental** analyses
* Outputs designed for **manual inspection and interpretation**

---

## Current Status & Roadmap

✔ Spatial clustering and annotation validated   
✔ TLS-like niche identification supported by multiple spatial metrics   
⏳ Pathway enrichment under biological review    
⏳ Deconvolution pending reference validation

Future updates will be versioned and documented.

---

## Citation & Use

This repository is intended for:

* Educational purposes
* Methodological demonstration
* Portfolio and reproducibility examples

If you use or adapt this pipeline, please cite appropriately.

---

## Contact

**Varun Umesh Gowda**
MS Bioinformatics, Northeastern University  
Bioinformatics Research Assistant, Brigham & Women’s Hospital   
📧 [gowda.var@northeastern.edu](mailto:gowda.var@northeastern.edu)
🔗 LinkedIn: *www.linkedin.com/in/varun-u-gowda*

---

If you want, next I can:

* Create a **GitHub project description + tags**
* Write a **LinkedIn post** announcing this project
* Help version this as **v2.0 (Validated Spatial Annotation)**
