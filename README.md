# Avian_tRNAs
Exploring the evolution of tRNAs in bird genomes

This project contains python and R-scripts to analyze the tRNA content of bird genomes. The workflow is described below

1. GETTING tRNA INFORMATION

tRNAscan.py: Run tRNAscan-SE on all available bird genomes (combination of three scripts)

2. GENOMIC DISTRIBUTION

2a. GENERAL LOCATIONS ON CHROMOSOMES

tRNA_Clusters.py: Extract genomic locations of tRNAs from specific text-file and check if they cluster together

2b. CLUSTERING ON CHROMOSOMES

First distances between consecutive tRNAs are calculated (tRNA_Clusters.py). The output of this analysis is then used to calculate the median distance, genome-wide and chromosome-specific (

Plot_tRNAs_on_Chromosomes.R: Creates plot with tRNAs mapped to chicken chromosomes (other species have specific scripts)

tRNA_Cluster_Contents.R: Creates Histogram that shows number of tRNAs per cluster and whether the clusters are homogeneous (i.e. coding for same aminoacids) or heterogeneous (i.e. coding for different aminoacids)
