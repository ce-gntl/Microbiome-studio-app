# Microbiome-studio-app
The application is organized into five main modules, each designed to support a complete microbiome data analysis workflow.

1. Data Upload

Allows users to upload .csv files containing:

OTU abundance data

Patient metadata

Uploaded data are used to construct a phyloseq object, which serves as the foundation for all downstream analyses.

2. Compositional Analysis

Provides interactive visualizations, including:

Venn diagrams

Heatmaps

Bar plots

Visualizations dynamically update based on user-defined:

Detection thresholds

Prevalence thresholds

3. Alpha Diversity

Computes standard alpha diversity indices

Enables group comparisons using the Wilcoxon–Mann–Whitney test

Results are presented through interactive tables and plots

4. Beta Diversity and PERMANOVA

Generates ordination plots using:

PCoA on Bray–Curtis or Jaccard distances

Includes:

Dispersion analysis via PERMDISP (betadisper)

Statistical testing using PERMANOVA (adonis2)

Outputs are visual and exportable

5. ANCOM-BC Analysis

Performs:

Data normalization

Compositional bias correction

Differential abundance testing

Produces:

Log fold change plots with confidence intervals

Downloadable result tables

Interactivity and Data Export

All tabular outputs are rendered using R’s DT package

Features include:

Sorting

Searching

Direct downloading of tables

Implementation Details

The application was tested locally

Designed for future server deployment

Codebase is modularized into separate scripts:

ui.R

server.R

Supporting function scripts

This structure ensures maintainability, scalability, and ease of future development.
