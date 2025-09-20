---
title: "01_differential_expression_analysis"
date: "2025-09-19"
---

##### 1. Load Packages
```{r}
library(ggplot2)
library(dplyr)
library(tidyr)
library(tidyverse)
library(cowplot)
library(pheatmap)
library(tidyplots)
library(edgeR)
library(MicrobiotaProcess)
library(clusterProfiler)
library(enrichplot)
library(Cairo)
library(ComplexUpset)
library(clusterProfiler)
library(KEGGREST)
library(RColorBrewer)
library(readxl)
```

##### 2. Load Files
```{r}
# Input RNASeq File - post internal standard correction
input_file <- "updated_comm_dyn_wide_counts.tsv"

# Combined metadata key #
metadata_table <- "pooled_key.csv"

# Annotations file of all genomes #
annotation_file <- 'updated_comm_dyn_annotations_v2.xlsx'
annotations <- read_excel(annotation_file)
annotations

# Subsetted annotations df for easier parsing #
small_annotations <- subset(annotations, select = c('ID','em_Preferred_name','em_desc'))


```



##### 3. QC Analysis
```{r,fig.width=16, fig.height=10}
##### 1. Preprocess samples and metadata for edgeR input #####
process_input_file_QC <- function(input_file,metadata_table,control,treatment,experiment_day,genome) {
  # Load Input File #
  raw_counts <- read.csv(file = input_file,sep = '\t')
  raw_counts<- raw_counts%>% filter(organism == genome)

  # Remove outliers #
  raw_counts <- subset(raw_counts, select = -c(`MED4ax_all_B_DAY2`,`MED4ax_1907_C_DAY2`,organism))
  
  # Clean df # 
  rownames(raw_counts) <- raw_counts$ID
  #raw_counts <- subset(raw_counts, select = -ID)

  # Load metadata and filter to comparison groups # 
  metadata_table <- read.csv(file = metadata_table)

  # Convert df to long format #
  df_long <- pivot_longer(raw_counts, 
                        cols = starts_with("MED4ax_"),
                        names_to = c("experiment_type", "rep", "day"),
                        names_pattern = "MED4ax_([\\w\\d]+)_(\\w+)_(DAY\\d+)",
                        values_to = "count")
  # Convert rep and day to factors for easier manipulation #
  df_long<- df_long%>% filter(day != experiment_day)
  df_long$rep <- as.factor(df_long$rep)
  df_long$day <- as.factor(df_long$day)
  df_long$experiment_type <- as.factor(df_long$experiment_type)

  # Return Cleaned Table #
  return(df_long)
}

##### 2. Plot Boxplot for each gene in Analysis #####
plot_boxplot<- function(data_processing_output_QC) {
  # Graph Data #
  overview_boxplot <- ggplot(data_processing_output_QC, aes(x = rep, y = count,color=rep)) +
    geom_boxplot() +
    facet_grid(day ~ experiment_type, scales = "free_x") +  # Split by day and experiment type, free scales for better visualization
    labs(x = "Replicate", y = "Gene Count") +
    theme_classic() +
    theme(
      axis.text.x = element_text(family = "Helvetica",color='black',size=16),
      #axis.text.x = element_text(family = "Helvetica",color='black',size=16,angle = 90,hjust=0.9,vjust=0.5),
      #axis.text.x = element_blank(),
      axis.text.y = element_text(family = "Helvetica",color='black',size=16),
      #axis.title.x = element_text(family = "Helvetica",color='black',size=22),
      axis.title.x = element_blank(),
      axis.title.y = element_text(family = "Helvetica",color='black',size=16),
      strip.text = element_text(family = "Helvetica",color='black',size=16),
      legend.text = element_text(family = "Helvetica",size = 16),
      legend.title = element_text(family = "Helvetica",size = 16),
      panel.border=element_blank(),
      panel.background = element_rect(colour = "black", size=1),
      legend.position='top')+
    scale_color_brewer(palette = "Set2")+
    scale_y_log10()

  return(overview_boxplot)
}

plot_correlation_plot<- function(data_processing_output_QC) {
  
  data_processing_output_QC$experiment_type <- gsub("1904", "alteromonas", data_processing_output_QC$experiment_type)
  data_processing_output_QC$experiment_type <- gsub("1905", "thalassospira_1", data_processing_output_QC$experiment_type)
  data_processing_output_QC$experiment_type <- gsub("1907", "thalassospira_2", data_processing_output_QC$experiment_type)
  data_processing_output_QC$experiment_type <- gsub("1940", "pseudohoeflea", data_processing_output_QC$experiment_type)
  data_processing_output_QC$experiment_type <- gsub("1943", "marinobacter", data_processing_output_QC$experiment_type)
  data_processing_output_QC$experiment_type <- gsub("all", "community", data_processing_output_QC$experiment_type)
  
  data_processing_output_QC$experiment_type <- factor(data_processing_output_QC$experiment_type, levels=c('control','marinobacter','thalassospira_1','thalassospira_2','pseudohoeflea','alteromonas','community'))
  
  
  # List Replicate Comparisons #
  comparisons <- list(c("A", "B"),c("B", "C"),c("A", "C"))
  
  # Create scatter plots for A vs B, A vs C, B vs C #
  # Note: Ensure all columns are present in the wide format df #
  # Change out A, B, or C depending on which two you want to compare #
  df_scatter <- do.call(rbind, lapply(comparisons, function(comp_group) 
    {df_temp <- data_processing_output_QC %>%
      filter(rep %in% comp_group) %>%
      pivot_wider(names_from = rep, values_from = count) %>%
      drop_na() %>%
      mutate(comparison_group = paste(comp_group[1], "vs", comp_group[2])) %>%
      rename(rep1 = !!comp_group[1], rep2 = !!comp_group[2]) %>%
      mutate(log10_rep1 = log10(rep1),log10_rep2 = log10(rep2))
    }))
  
  # Make comparison_group a character #
  df_scatter$comparison_group <- as.character(df_scatter$comparison_group)
  
  # Calculate Spearman correlation for each combination of day and experiment_type #
  correlation_df <- df_scatter %>%
    complete(day, experiment_type, comparison_group = comparison_group, fill = list(spearman_coefficient = NA)) %>%
    group_by(comparison_group,day,experiment_type) %>%
    summarise(spearman_coefficient = cor(log10_rep1, log10_rep2, method = "spearman"))
  
  # Add spearman results to df #
  df_scatter <- df_scatter %>%
    left_join(correlation_df, by = c("comparison_group","day","experiment_type"))
  
  # Summarize Spearman correlation coefficients by comparison group
  spearman_summary_table <- correlation_df %>%
    group_by(day,experiment_type) %>%
    summarize(spearman1 = spearman_coefficient[comparison_group == "A vs B"],  # Ensure correct order for A vs B
              spearman2 = spearman_coefficient[comparison_group == "A vs C"],  # Ensure correct order for A vs C
              spearman3 = spearman_coefficient[comparison_group == "B vs C"]) %>%
    mutate(spearman_text = paste0(
      "A vs B: ", round(spearman1, 2), "\n",
      "A vs C: ", round(spearman2, 2), "\n",
      "B vs C: ", round(spearman3, 2)))
  
  
  # Plot scatterplots for each comparison type
  correlation_plot <- ggplot(df_scatter, aes(x = log10_rep1, y = log10_rep2,color=comparison_group)) +
    geom_point(size=1,alpha=0.5) +  
    facet_grid(day ~ experiment_type) +
    labs(x = "log10 Replicate 1", y = "log10 Replicate 2") +
    theme_classic()+
    theme(
      axis.text.x = element_text(family = "Helvetica",color='black',size=16),
      axis.text.y = element_text(family = "Helvetica",color='black',size=16),
      axis.title.x = element_text(family = "Helvetica",color='black',size=16),
      plot.title = element_text(family = "Helvetica",color='black',size=16,hjust = 0.5),
      axis.title.y = element_text(family = "Helvetica",color='black',size=16),
      strip.text = element_text(family = "Helvetica",color='black',size=16),
      legend.text = element_text(family = "Helvetica",size = 16),
      legend.title = element_text(family = "Helvetica",size = 16),
      #legend.key.size = unit(10, "cm"),  # Increase legend key size
      panel.border=element_blank(),
      panel.background = element_rect(colour = "black", size=1),
      legend.spacing = unit(1.5, "cm"),
      legend.position='top')+
    scale_color_brewer(palette = "Set2")+
    guides(color = guide_legend(override.aes = list(size = 5),title = "Replicate Comparison: "))+
    geom_text(data = spearman_summary_table,aes(x = Inf, y = Inf, label = spearman_text),
              #hjust = 2, vjust = 1,family = "Helvetica",color='black',inherit.aes = FALSE)+
              hjust = 1.1, vjust = 4.9,family = "Helvetica",color='black',inherit.aes = FALSE)+
    scale_y_continuous(limits = c(2,14),breaks = c(4,8,12))+
    scale_x_continuous(limits = c(2,14),breaks = c(4,8,12))
  
  return(correlation_plot)
}

##### Define Inputs #####
# Subset to one genome #
genome<-"MED4"

# What day do you want to remove? - we removed day 6 because cells were entering stationary #
experiment_day <- "DAY6"

##### Run Preprocessing Step #####
data_processing_output_QC <-process_input_file_QC(input_file,metadata_table,control,treatment,experiment_day,genome)
head(data_processing_output_QC)

##### Generate Replicate Boxplot #####
boxplot <- plot_boxplot(data_processing_output_QC)
boxplot

##### Generate Replicate Correlation Plots #####
# For this plot you have to go to the function and comment out the different replicates for different comparisons #
correlation_plot <- plot_correlation_plot(data_processing_output_QC)
correlation_plot

```

##### 4. edgeR analysis
```{r}
##### 1. Preprocess samples and metadata for edgeR input #####
process_input_file <- function(input_file,metadata_table,control,treatment,day,genome) {
  # Load Input File #
  raw_counts <- read.csv(file = input_file,sep = '\t')
  raw_counts<- raw_counts%>% filter(organism == genome)

  # Remove outliers #
  raw_counts <- subset(raw_counts, select = -`MED4ax_all_B_DAY2`)
  raw_counts <- subset(raw_counts, select = -`MED4ax_1907_C_DAY2`)
  raw_counts <- subset(raw_counts, select = -organism)
  
  # Clean df # 
  rownames(raw_counts) <- raw_counts$ID
  raw_counts <- subset(raw_counts, select = -ID)

  # Load metadata and filter to comparison groups # 
  metadata_table <- read.csv(file = metadata_table)
  metadata_table<- metadata_table%>% filter((experimental_trial == control)|(experimental_trial == treatment))
  metadata_table<- metadata_table%>% filter(sampling_day == day)
  rownames(metadata_table) <- metadata_table$sample
  
  # Ensure metadata and gene count tables have the same rows/columns as each other #
  metadata_table <- metadata_table[rownames(metadata_table) %in% colnames(raw_counts), ]
  matching_columns <- colnames(raw_counts) %in% rownames(metadata_table)
  raw_counts_filtered <- raw_counts[, matching_columns]

  # Replace NA values with 0 #
  raw_counts_filtered[is.na(raw_counts_filtered)] <- 0
  metadata_table <- metadata_table[match(colnames(raw_counts_filtered), rownames(metadata_table)), ]
  
  # Return Cleaned Tables #
  return(list(cleaned_raw_counts = raw_counts_filtered, cleaned_metadata = metadata_table))
}

##### 2. Run EdgeR on Input Data #####
run_edgeR <- function(cleaned_metadata,cleaned_raw_counts,small_annotations,genome) {
  
  # Extract experimental groups from metadata - ensure control is the reference group #
  group <- cleaned_metadata$experimental_trial
  group <- factor(group)
  group <- relevel(group, ref = control)

  # Create edgeR object #
  edgeR_output <- DGEList(counts = cleaned_raw_counts, group = group)
  
  # Calculate the normalization factors for edgeR - using TMM #
  edgeR_output <- calcNormFactors(edgeR_output,method='TMM')

  # Create a design matrix for analysis #
  design <- model.matrix(~ group, data = cleaned_metadata)
  
  # Estimate dispersion for samples #
  edgeR_output <- estimateDisp(edgeR_output, design)
  
  # Fit the model using glmQLFit and Robust=True #
  fit <- glmQLFit(edgeR_output, design,robust=TRUE)
  
  # Obtain statistical significance #
  results <- glmQLFTest(fit)
  
  # Obtain results #
  results_df <- topTags(results, n = Inf)$table
  results_df$ID <- rownames(results_df)

  # Align Output with annotations #
  results_w_annot <- merge(results_df, small_annotations, by = "ID")
  
  # Optional Filter out non significant genes #
  #results_w_annot<- results_w_annot%>% filter(FDR < 0.05 )
  results_w_annot
  return(results_w_annot)
}

##### Define Appropriate Input Files for EdgeR Analysis #####
##### Make sure to define control vs treatment - control results in negative LFC and can be defined by user #####
##### Define day for comparison and genome/organism that you want to perform edgeR on #####

# Define Reference Genome #
genome<-"MED4"
#genome<-"Alteromonas_1904"
#genome<-"Pseudohoeflea_1940"
#genome<-"Thalassospira_1907"
#genome<-"Marinobacter_1943"

# Define control group (negative LFC) #
#control <-"control"
control <- "marinobacter"
#control <- "alteromonas"
#control <- "thalassospira_1"
#control <- "thalassospira_2"
#control <- "pseudohoeflea"

# Define treatment group (positive LFC) #
#treatment <- "marinobacter"
#treatment <- "alteromonas"
#treatment <- "thalassospira_1"
#treatment <- "thalassospira_2"
#treatment <- "pseudohoeflea"
treatment <- "community"

# Choose the data that you want to analyze #
#day <- "2"
#day <- "4"
day <- "5"

##### Run Preprocessing Step #####
data_processing_output <-process_input_file(input_file,metadata_table,control,treatment,day,genome)

# Output metadata #
cleaned_metadata <- data_processing_output$cleaned_metadata
head(cleaned_metadata)

# Output the raw count data - in appropriate format for edgeR #
cleaned_raw_counts <- data_processing_output$cleaned_raw_counts
head(cleaned_raw_counts)

##### Run EdgeR Step #####
edgeR_output <- run_edgeR(cleaned_metadata,cleaned_raw_counts,small_annotations,control)
head(edgeR_output)

##### Save output files #####
# Finalized gene tables were created using excel after running all edgeR data #
#write.csv(edgeR_output, "comm_dyn_edgeR_exp1_exp_analysis_Pro_day2.csv", row.names = FALSE)
```
