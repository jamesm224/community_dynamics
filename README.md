# Community Dynamics Working Title

## Introduction

Here is the code associated with this project. It includes all of the code to reproduce all bioinformatic analysis performed in this paper. It contains code for genome assembly, RNASeq analysis, metagenomic, metatranscriptomic, and additional analysis.

## Publication

Citation [Blank]

## File Structure
* bin
* other folders

## Reproducing Data Analysis
1. Download raw data from NCBI (BioProject: )
2. Process RNA Seq data through a series of scripts to process RNASeq
  - The following Github link was used: https://github.com/nhinvo/rnaseq-absolute-pipeline
  - The pipeline utilizes the following packages:
    - Snakemake
    - bbtools
    - bowtie2
    - HTSeq
    - Python
3. Process metagenomic data using bowtie and read map to the known reference genomes
  - Genomes used in this study
    - Marinobacter ()
    - Thalassospira ()
    - Alteromonas ()
    - Pseudohoeflea ()
    - MED4 () 
    - Thermus Thermophilus ()
  - Scripts used in data analysis ()
    - Contain QC plots, comparison between FCM and metagenomic data
    - Visualize abundance of reads using bowtie2
    - Visualize abundance of reads corrected by estimated extraction efficiency
4. Process the RNASeq data. The remainder of Data Analysis was performed in R.
  - Process raw data to obtain edgeR data.
    - Relative read counts were used for Differential Expression Analysis
    - Includes ClusterProfiler Pathway Enrichement Analysis
  - Analyze Differential Expression Analysis results for RNASeq data
    - Includes heatmaps, pathway plots, and LFC analysis
  - Analyze the absolute count RNASeq data
    - Includes the conversion of relative to absolute count data
    - Incorporates a customizable script for generating gene diagrams for all desired KEGG pathways
