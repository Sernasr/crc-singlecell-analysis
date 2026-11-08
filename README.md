# Colorectal Cancer Single-Cell RNA-seq Analysis

## Overview

This project explores single-cell RNA sequencing (scRNA-seq) data from colorectal cancer (CRC), with a focus on cellular heterogeneity, microsatellite instability (MSI) status, and differences between MSI-high (MSI-H) and microsatellite-stable (MSS) tumors.

Two published datasets are analyzed separately and compared to investigate immune cell populations, tissue-specific distributions, and consistency of cell annotations.

## Datasets

- **Zhang et al. (2018) — GSE108989:** Analysis of T-cell populations, including CD4+ and CD8+ T cells and their subclusters.
- **Zhang et al. (2020) — GSE146771:** Analysis of major cell populations, including T cells, B cells, myeloid cells, innate lymphoid cells, epithelial cells, and fibroblasts.

Both datasets are publicly available through the NCBI Gene Expression Omnibus (GEO).

## Analysis

The project includes:

- Exploratory analysis of cell populations and patient characteristics.
- Visualization of cellular heterogeneity using t-SNE.
- Comparison of cell-type proportions across tissues.
- Patient-level comparisons between MSI-H and MSS tumors.
- Statistical testing using Mann–Whitney U and Welch's t-test, with Benjamini–Hochberg FDR correction.
- Cross-dataset comparison of shared patients and cell annotations.

## Repository Structure

- `notebooks/` — Jupyter notebooks for the 2018 dataset, 2020 dataset, and cross-dataset comparison.
- `results/GSE108989/` — Figures and statistical results from Zhang 2018.
- `results/GSE146771/` — Figures and statistical results from Zhang 2020.
- `results/dataset_comparison/` — Cross-dataset comparison figures and tables.

## Interpretation and Limitations

The analyses are exploratory. Statistical comparisons are based on patient-level cell proportions, and the limited number of MSI-H patients restricts statistical power.

The datasets contain overlapping patients and cells and therefore should not be treated as independent validation cohorts. Cell proportions also reflect the original cell-sorting and sampling strategies rather than unbiased whole-tissue composition.

## Tools

Python, Jupyter Notebook, pandas, NumPy, Matplotlib, Seaborn, SciPy, and scikit-learn.

## Data Availability

The original datasets are publicly accessible through GEO:

- [GSE108989](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE108989)
- [GSE146771](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE146771)

Large raw data and expression matrices are not included in this repository.
