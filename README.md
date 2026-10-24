---
title: "Pipeline: non-additive genetics of tree architecture in columnar apples"
author: "Nuri Güvencli"
date: "24.10.2026"
---

## Overview
This repository contains the R analysis workflow used in the study: "Epistatic interactions underlie phenotypic diversity of tree architecture in columnar apple." 
The workflow processes tree architectural data and SNP genotype data to partition genetic variance and map epistatic interactions for five tree architecture traits: 
Internode Length (INL), Long Shoot Angle (LSA), Trunk Branching Zone Proportion (TBZP), Trunk Conicity (TC), and Trunk Length (TL).

## Pipeline Workflow:

* Outlier removal from raw phenotypic data.
* Estimation of breeding values (BLUPs) and narrow-sense heritability (h2) using GBLUP model with spatial field corrections.
* Principal Component Analysis (PCA) of estimated breeding values.
* Comparative assessment of h2 models (e.g., impact of including the Co locus as a fixed effect).
* GWAS (BLINK) performed on deregressed estimated breeding values (dGEBVs), with genomic control before the Benjamini-Hochberg FDR.
* SNP × Co-locus and SNP × SNP epistatic interactions using Wald tests (at most the 30 most significant GWAS SNPs per trait-year).
* Generation of Circos plots for visual QTL mapping.
* Estimation of genetic correlations (rg) across traits and years.
* Visualization of data, (phenotype density distribution, correlation heatmaps, Variance component barplots, Genome-wide interaction map).

## Project Structure
To execute the pipeline, ensure your local directory is organized as follows:
* `/SNP Chip`: Contains genomic input files (e.g., `GM_v11.csv`, `GT_filtered_numeric_transposed.csv`).
* `/Lidar/Single_trait_walkthrough`: Root directory for phenotypic data (`/outlier_removed_raw`) and subdirectories for specific model outputs (`/deregression_output`, `/MV5_Rg_output`, etc.).

## Dependencies
* **Core R Packages:** `asreml`, `ASRgenomics`, `GAPIT`, `BLINK`, `tidyverse`, `circlize`, `ComplexHeatmap`, `svglite`.

## Correction (October 2026)
An error found at the proof stage was corrected in `6_GBLUP.R`, `13_SNPxCo_GBLUP.R` and `14_Circosplot.R`.

* `6_GBLUP.R`: GEBVs are now the prediction of the genetic term only (`predict(..., only = "vm(Genotype, Ginv)")`). Before, `predict()` also returned the intercept and the averaged population effect, which the deregression (`7_Deregression.R`) then divided by each genotype's accuracy.
* `13_SNPxCo_GBLUP.R` and `14_Circosplot.R`: genomic control (lambda = median chi-square / 0.455) is applied to the GWAS p-values before the Benjamini-Hochberg FDR. The joint interaction model is fitted with at most the 30 most significant SNPs per trait-year (this affects only TBZP 2018) and without SNP × SNP terms when only one SNP is significant.

Affected results are Figures 1 and 5, Figures S3 and S7 and Tables S3 to S6. Heritabilities, variance components and genetic correlations are unchanged.

---
*For detailed information regarding statistical methodology, experimental design, and biological interpretation, please refer to the project manuscript.*
