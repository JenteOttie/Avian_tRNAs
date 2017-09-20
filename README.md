# Avian_tRNAs
Exploring the evolution of tRNAs in bird genomes

This project contains python and R-scripts to analyze the tRNA content of bird genomes. The workflow is described below.

# 1. GETTING tRNA INFORMATION

*1a. GENERAL PATTERNS*

tRNAscan.py: Run tRNAscan-SE on all available bird genomes (combination of three scripts)

Heatmap.R: Visualize tRNA content of bird genomes in a heatmap

*1b. FILTERING OUT TE-DERIVED tRNAs*

The general overview of tRNA content in bird genomes revealed that certain species contain extremely high numbers of particular tRNAs (e.g., Isoleucine in Trogon). These tRNAs are derived from transposable elements (TEs). The following scripts were used to filter these TE-derived tRNAs out.

*WORK IN PROGRESS*

# 2. GENOMIC DISTRIBUTION

*2a. GENERAL LOCATIONS ON CHROMOSOMES*

Plot_tRNAs_on_Chromosomes.R: Creates figure with tRNAs mapped to chicken chromosomes (other species have specific scripts)

*2b. CLUSTERING ON CHROMOSOMES*

First distances between consecutive tRNAs are calculated (tRNA_Distances.py). The output of this analysis is then used to calculate the median distance, both genome-wide and chromosome-specific (Median_Distance.R). The genome-wide median distance is used as a threshold to check for clustering of tRNAs (tRNA_Clusters.py). The nature of these clusters - size and content - is visualized in a histogram (tRNA_Cluster_Contents.R)

tRNA_Distances.py: Extract genomic locations of tRNAs from specific text-file and calculate distance between consecutive tRNAs

tRNA_Cluster_Contents.R: Creates Histogram that shows number of tRNAs per cluster and whether the clusters are homogeneous (i.e. coding for same aminoacids) or heterogeneous (i.e. coding for different aminoacids)
