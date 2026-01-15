# scRNA-seq Analysis Pipeline for GSE138852 (Alzheimer's Disease)

This repository contains a comprehensive single-cell RNA sequencing (scRNA-seq) analysis pipeline for the dataset **GSE138852**. The workflow covers everything from quality control to trajectory analysis and cell-cell communication.

## Analysis Workflow
1.  **QC & Filtering**: Removal of low-quality cells based on MT-percentage and feature counts.
2.  **Normalization**: Using `SCTransform` to regress out mitochondrial variation.
3.  **Annotation**: Automated cell type identification via `SingleR`.
4.  **Differential Expression**: Identification of markers for clusters and AD vs. Control comparisons.
5.  **Enrichment**: GSEA using MSigDB (GO Biological Processes).
6.  **Advanced Dynamics**: Trajectory inference with `Monocle3` and signaling analysis with `CellChat`.

## Requirements
- R >= 4.1
- Seurat v4
- Monocle3
- CellChat
- clusterProfiler

## Usage
1. Clone the repository.
2. Place `GSE138852_counts.csv.gz` and `GSE138852_covariates.csv.gz` in the root directory.
3. Run `Rscript analysis_pipeline.R`.

## Outputs
Results are automatically organized into:
- `/outputs/figures`: UMAPs, Heatmaps, Dotplots, and Networks.
- `/outputs/tables`: DEGs, GSEA results, and Pathway lists.
- `/outputs/objects`: Final processed Seurat object.
