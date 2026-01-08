# Microbiome-studio-app

The application is organized into five main modules:
•	Data Upload, which allows the user to upload .csv files containing information on the abundances of each OTU and patient metadata. These are used to create a phyloseq object that will then be used in subsequent analyses.
•	Compositional Analysis, with interactive visualizations such as Venn diagrams, heatmaps, and barplots that dynamically update with user-selected detection and prevalence thresholds.
•	Alpha Diversity, which produces diversity indices and allows comparisons between groups using the Wilcoxon-Mann-Whitney test.
•	Beta Diversity and PERMANOVA, which includes the generation of an ordination plot (PCoA on Bray-Curtis or Jaccard distances), dispersion analysis (PERMDISP with the "betadisper" function), and significance testing using "adonis2".
•	ANCOM-BC analysis, which performs data normalization, compositional bias correction and differential abundance testing, generating graphical outputs (log fold change with confidence intervals) and exportable tables.
The displayed data is made interactive using R's DT package, allowing sorting, searching, and direct downloading. The application was tested locally and designed for future deployment on a server. The entire code was organized into separate scripts (UI, server, and functions) to ensure the project's maintainability and scalability.
