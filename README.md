# Community Dynamics Working Title

## Introduction

Here is the code associated with this project. It includes all of the code to reproduce all bioinformatic analysis performed in this paper. It contains code for RNASeq Analysis, metagenomic, metatranscriptomic, and additional analysis.

## Publication

Citation [Blank]

## Table of Contents
* Overall Structure
* 

## Reproducing Data Analysis
1. Download raw data from NCBI (BioProject: )
2. Process RNA Seq data through a series of scripts to process RNASeq
  - The following Github link was used: https://github.com/nhinvo/rnaseq-absolute-pipeline
  - The pipeline utilizes the following packages:
      a. Snakemake
      b. bbtools
      c. bowtie2
      d. HTSeq
      e. Python
3. Process metagenomic data using bowtie and read map to the known reference genomes
  - Marinobacter ()
  - Thalassospira ()
  - Alteromonas ()
  - Pseudohoeflea ()
  - MED4 () 
  - Thermus Thermophilus ()
4. Process the RNASeq data using edgeR. The remainder of Data Analysis was performed in R.
  - All code is located here

## Packages and Programs Used


