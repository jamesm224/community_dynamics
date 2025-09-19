# Community Dynamics Working Title

## Introduction

Here is the code associated with this project. It includes all of the code to reproduce all bioinformatic analysis performed in this paper. It contains code for genome assembly, RNASeq analysis, metagenomic, metatranscriptomic, and additional analysis.

## Publication

Citation [Blank]

## Packages used:
Pipelines: 
   - Snakemake
   - bbtools
   - bowtie2
   - HTSeq
   - Python

Data Analysis: There were a few others that weren't used much but here are most of them!
   - ggplot2 v3.5.2
   - dplyr v1.1.4
   - scales v1.4.0
   - pheatmap v1.0.13
   - tidyplots v0.3.1
   - tidyr v1.3.1
   - edgeR v4.2.2
   - MicrobiotaProcess v1.16.1
   - clusterProfiler v4.12.6
   - KEGGREST v1.44.1
   - ComplexUpset v1.3.3

## File Structure
   - pipelines
      - HIFI-genome-closing-improved - Contains genome closing workflow
      - rnaseq-absolute-pipeline - Contains scripts for RNASeq analysis
   - scripts
      - analyze_metaG_data - Contains metagenomic data analysis
      - analyze_metaT_data - Contains metatranscriptomic data analysis
      - replicate analysis - Contains Thalassospira strain-level analysis

## Reproducing Data Analysis
1. Download raw data from NCBI (BioProject: )
   
2. Process RNA Seq data through a series of scripts to process RNASeq (Nhi could you fill this out when you get the chance please?)
  - The following Github link was used: https://github.com/nhinvo/rnaseq-absolute-pipeline
  - The pipeline utilizes the following packages:
      
3. Process metagenomic data using bowtie and read map to the known reference genomes
  - Genomes used in this study
    - Marinobacter ()
    - Thalassospira ()
    - Alteromonas ()
    - Pseudohoeflea ()
    - MED4 () 
    - Thermus Thermophilus ()
      
  - Scripts used in data analysis (scripts/analyze_metaG_data/metagenomic_analysis_code.Rmd)
    - Contain QC plots, comparison between FCM and metagenomic data
    - Visualize abundance of reads using bowtie2
    - Visualize abundance of reads corrected with and without estimated extraction efficiency
      
4. Process the RNASeq data. The remainder of Data Analysis was performed in R.
  - Process raw data to obtain edgeR data.
    - Relative read counts were used for Differential Expression Analysis
    - Includes ClusterProfiler Pathway Enrichement Analysis
      
  - Analyze Differential Expression Analysis results for RNASeq data
    - Includes heatmaps, pathway plots, and LFC analysis
      
  - Analyze the absolute count RNASeq data
    - Includes the conversion of relative to absolute count data
    - Incorporates a customizable script for generating gene diagrams for all desired KEGG pathways
