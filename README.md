# The shared and distinct roles of Prochlorococcus and co-occurring heterotrophic bacteria in regulating community dynamics

## Introduction
Prochlorococcus is one of the most important photosynthetic organisms in Earth's oceans. While this microbe is abundant, there still remains a limited understanding of how Prochlorococcus interacts with heterotrophs. Here we isolated four of the most abundant strains associated with MED4 and through the use of absolute quantification of RNA, DNA, and cell counts we examine how these heterotrophs impact the growth curve of Prochlorococcus. Here is the code associated with this project. It includes all of the code to reproduce all bioinformatic analysis performed in this paper. It contains code for genome assembly, RNASeq processing, metagenomic, metatranscriptomic, and additional analysis.

## Publication
1. Ziegler, C.A., Mullet, J.I., Coe, A., Vo, N.N., Salcedo, E., Arrigan, D.M., Parker, S.M., Chisholm, S.W. (2025). The shared and distinct roles of Prochlorococcus and co-occurring heterotrophs in regulating community fitness, as revealed by synthetic communities. (In Preparation).

## Packages 
RNAseq Pre-processing Pipeline: 
   - Snakemake v7.32.4  
   - bbtools v39.18
   - bowtie2 v2.5.4
   - HTSeq v2.0.9
   - Python v3.12.2

Genome Assembly Pipeline: 
  - samtools v1.22.1
  - metaFlye v2.9.5
  - MMseqs2 v17.b804f
  - GTDB v207

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
### Install FASTQ files from SRA 
1. Download raw data from NCBI (BioProject: )

### Genome Assembly  
Process PacBio long read sequencing samples through a series of scripts located in `pipelines/HIFI-genome-closing-improved` to assemble reference genomes.  

**Installation**  
  - Install: samtools v1.22.1, metaFlye v2.9.5, MMseqs2 v17.b804f, GTDB v207

**Setup and Run Scripts**
  - Set up paths to output directories and input files in each script
  - Edit resource specifications and conda package names 
  - Run scripts to process the samples


### RNAseq Pre-processing
Process meta-transcriptomics samples through Snakemake pipeline located in `pipelines/rnaseq-absolute-pipeline` to obtain internal standard normalized read counts.  

**Installation**  
  - Install Snakemake v7.32.4 and Conda (or Mamba)

**Set up Snakemake Pipeline**  
  - Create samples.tsv file: 
    - Required columns: 
      1. "sample" - Unique sample name 
      2. "forward read" - Absolute path to forward reads .fastq file 
      3. "reverse read" - Absolute path to reverse reads .fastq file
    - Optional: any additional columns with sample metadata 
  - Create internal_standard_concentration.tsv file:
    - Required columns: 
      1. "standard_name" - Unique name of internal standard added. 
          - **Important**: Column value should match with name of sequence in provided FASTA and GFF file
      2. "standard_group" - Group that standard is in 
      3. "concentration (ng/ul)" - Concentration of standard added 
      4. "volume_added (ul)" - Volume of standard added 
  - Create cell_count.tsv file:
    - Required columns: 
      1. "sample" - Unique sample name. Should match with "sample" column from samples.tsv	
      2. "total_cell_count" - Count of cells in sample 
  - Edit config.yaml file:
    1. Edit names (yaml keys) of reference genomes and their paths
    2. Edit path to folder to store intermediate files: "scratch directory"
    3. Edit length of sequenced read (from fastq file): "read length"
    4. Edit minimum number of samples a standard has to be in: "minimum standard sample"
  - Edit profile/config.yaml file:
    1. Edit partition name in "default-resources" - "partition"
    2. Edit any other resources as needed 
  - Edit run_RNAseq_SM.sbatch file: 
    1. Edit slurm SBATCH resource specifications as needed (e.g. time, partition)

**Running Snakemake pipeline**  
  1. Run pipeline by: `sbatch run_RNAseq_SM.sbatch`
      - Note: create logs/ folder before submitting job 
  2. Check log files in logs/ folder 

### RNAseq Post-processing
Process the RNASeq data. The remainder of Data Analysis was performed in R.
This data analysis was developed to be easily reproducible and understand.
All corresponding .Rmd files also include a visual .md file that displays the results from the R scripts in a user friendly way.
The code is located here: scripts/analyze_metaT_data

   - 01_differential_expression_analysis: Process raw data to obtain differentially expressed genes
      - Process raw data to obtain edgeR data.
      - Relative read counts were used for Differential Expression Analysis
        
   - 02_processing_relative_RNASeq_data: Analyze differentially expressed genes
      - Contains overview heatmaps, dotplots, COG overview, LFC analysis
      - Includes ClusterProfiler Pathway Enrichement Analysis
        
   - 03_analyze_KEGG_pathways: Pathway analysis
      - Contains KEGG and pathway analysis for heterotrophs and Prochlorococcus
      - Has overview heatmaps and specific pathway plot breakdown
        
  - 04_absolute_RNA_analysis: Examining cell count corrected data
     - Converts data to transcript per cell counts
     - Incorporates a customizable script for generating gene diagrams for all desired KEGG pathways
     - Contains a overview KEGG heatmap and heatmap for gene presence/absence analysis
       
### Metagenomics Analysis   
Process metagenomic data using bowtie and read map to the known reference genomes  
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

### Replicate Analysis   
We had two experimental trials using identical Thalassospira strains. 
To ensure that between multiple experimental trials we were getting the same results, we compared the differential expression results.
The data includes the outputted edgeR results for Prochlorococcus and Thalassospira in both trials.
Scripts are located here: scripts/replicate_analysis/thalassospira_replicate_analysis
