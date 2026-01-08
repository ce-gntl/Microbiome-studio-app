# Microbiome Analysis Shiny Application

An interactive R Shiny application for end-to-end microbiome data analysis, including compositional analysis, diversity metrics, PERMANOVA testing, and differential abundance analysis using ANCOM-BC.

## Table of Contents
  
- [Overview](#overview)
- [Features](#features)
- [Application Modules](#application-modules)
- [Input Data Format](#input-data-format)
- [Outputs](#outputs)
- [Technologies Used](#technologies-used)


## Overview

This application provides an interactive framework for analyzing microbiome sequencing data. Users can upload OTU abundance tables and associated metadata, explore compositional patterns, calculate diversity metrics, and perform differential abundance testing through a unified web interface.
The application is built using R Shiny and leverages established microbiome analysis packages such as phyloseq, vegan, and ANCOMBC.

## Features
-	Upload and validate microbiome count and metadata files
-	Interactive compositional visualizations
-	Alpha and beta diversity analysis with statistical testing
-	Differential abundance testing with compositional bias correction
-	Downloadable tables and publication-ready plots
-	Modular and scalable codebase

## Application Modules
### Data Upload

Upload .csv files containing:
-	OTU abundance data
-	Sample metadata
-	Data are combined into a phyloseq object for downstream analysis.
  
### Compositional Analysis

Interactive visualizations:

-	Venn diagrams
-	Heatmaps
-	Bar plots
-	Visualizations update dynamically based on:
-	Detection thresholds
-	Prevalence thresholds
  
### Alpha Diversity

- Calculates common alpha diversity indices
- Enables group comparisons using the Wilcoxon–Mann–Whitney test
  
### Beta Diversity and PERMANOVA

Ordination via PCoA using:
-	Bray–Curtis distances
-	Jaccard distances
	-	Includes:
	  -	Dispersion analysis (PERMDISP, betadisper)
	  -	PERMANOVA testing (adonis2)
### ANCOM-BC Analysis
Performs:
-	Normalization
-	Compositional bias correction
-	Differential abundance testing
Outputs:
-	Log fold change plots with confidence intervals
-	Exportable result tables

## Input Data Format
### OTU Table
- **Rows**: OTUs
- **Columns**: Samples
- **Values**: Raw or normalized counts
### Taxonomy Table 
- **Rows**: OTUs
- **Columns**: Taxonomic rank
### Metadata Table
- **Rows**: Samples
- **Columns**: Experimental variables (e.g., group, treatment, timepoint)

> Files must be provided in .csv format and contain matching sample IDs.

## Outputs
- Interactive tables (sortable, searchable, downloadable)
- Alpha diversity plots and statistics
- Beta diversity ordination plots
- PERMANOVA and dispersion results
- Differential abundance plots and tables from ANCOM-BC
  
> All tables are rendered using R’s DT package.

## Technologies Used
- R Shiny
- phyloseq
- vegan
- ANCOMBC
- ggplot2
- DT
