02_processing_relative_RNASeq_data
================
2025-09-19

##### 1. Load Packages

``` r
library(ggplot2)
library(dplyr)
```

    ## 
    ## Attaching package: 'dplyr'

    ## The following objects are masked from 'package:stats':
    ## 
    ##     filter, lag

    ## The following objects are masked from 'package:base':
    ## 
    ##     intersect, setdiff, setequal, union

``` r
library(tidyr)
library(tidyverse)
```

    ## ── Attaching core tidyverse packages ──────────────────────── tidyverse 2.0.0 ──
    ## ✔ forcats   1.0.0     ✔ readr     2.1.5
    ## ✔ lubridate 1.9.4     ✔ stringr   1.5.1
    ## ✔ purrr     1.1.0     ✔ tibble    3.3.0

    ## ── Conflicts ────────────────────────────────────────── tidyverse_conflicts() ──
    ## ✖ dplyr::filter() masks stats::filter()
    ## ✖ dplyr::lag()    masks stats::lag()
    ## ℹ Use the conflicted package (<http://conflicted.r-lib.org/>) to force all conflicts to become errors

``` r
library(cowplot)
```

    ## 
    ## Attaching package: 'cowplot'
    ## 
    ## The following object is masked from 'package:lubridate':
    ## 
    ##     stamp

``` r
library(pheatmap)
library(tidyplots)
```

    ## 
    ## Attaching package: 'tidyplots'
    ## 
    ## The following object is masked from 'package:cowplot':
    ## 
    ##     save_plot

``` r
library(edgeR)
```

    ## Loading required package: limma

``` r
library(MicrobiotaProcess)
```

    ## MicrobiotaProcess v1.16.1 For help:
    ## https://github.com/YuLab-SMU/MicrobiotaProcess/issues
    ## 
    ## If you use MicrobiotaProcess in published research, please cite the
    ## paper:
    ## 
    ## Shuangbin Xu, Li Zhan, Wenli Tang, Qianwen Wang, Zehan Dai, Lang Zhou,
    ## Tingze Feng, Meijun Chen, Tianzhi Wu, Erqiang Hu, Guangchuang Yu.
    ## MicrobiotaProcess: A comprehensive R package for deep mining
    ## microbiome. The Innovation. 2023, 4(2):100388. doi:
    ## 10.1016/j.xinn.2023.100388
    ## 
    ## Export the citation to BibTex by citation('MicrobiotaProcess')
    ## 
    ## This message can be suppressed by:
    ## suppressPackageStartupMessages(library(MicrobiotaProcess))
    ## 
    ## Attaching package: 'MicrobiotaProcess'
    ## 
    ## The following object is masked from 'package:tidyr':
    ## 
    ##     extract
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(clusterProfiler)
```

    ## 
    ## clusterProfiler v4.12.6 Learn more at https://yulab-smu.top/contribution-knowledge-mining/
    ## 
    ## Please cite:
    ## 
    ## G Yu. Thirteen years of clusterProfiler. The Innovation. 2024,
    ## 5(6):100722
    ## 
    ## Attaching package: 'clusterProfiler'
    ## 
    ## The following object is masked from 'package:purrr':
    ## 
    ##     simplify
    ## 
    ## The following object is masked from 'package:stats':
    ## 
    ##     filter

``` r
library(enrichplot)
library(Cairo)
library(ComplexUpset)
library(clusterProfiler)
library(KEGGREST)
library(readxl)
library(RColorBrewer)
```

##### 2. Load Files

``` r
# Annotations file of all genomes #
annotation_file <- 'updated_comm_dyn_annotations_v2.xlsx'
annotations <- read_excel(annotation_file)
annotations
```

    ## # A tibble: 19,501 × 48
    ##    em_ID     organism ID    seq_id source type   start    end score strand phase
    ##    <chr>     <chr>    <chr> <chr>  <chr>  <chr>  <dbl>  <dbl> <chr> <chr>  <dbl>
    ##  1 batch10_… Alterom… 1_18… batch… Prodi… CDS   2.21e6 2.21e6 69.00 -          0
    ##  2 batch10_… Alterom… 1_18… batch… Prodi… CDS   2.21e6 2.22e6 123.… -          0
    ##  3 batch10_… Alterom… 1_36… batch… Prodi… CDS   4.24e6 4.25e6 151.… +          0
    ##  4 batch10_… Marinob… 1_15… batch… Prodi… CDS   1.71e6 1.71e6 89.99 +          0
    ##  5 batch10_… Thalass… 1_39… batch… Prodi… CDS   4.39e5 4.40e5 127.… -          0
    ##  6 batch10.… Pseudoh… 1_15… batch… Prodi… CDS   1.64e6 1.64e6 115.… -          0
    ##  7 batch10_… Alterom… 1_11… batch… Prodi… CDS   1.32e6 1.33e6 376.… +          0
    ##  8 batch10_… Thalass… 1_35… batch… Prodi… CDS   3.85e6 3.85e6 422.… -          0
    ##  9 batch10_… Marinob… 1_29… batch… Prodi… CDS   3.23e6 3.23e6 458.… -          0
    ## 10 batch10.… Pseudoh… 1_46… batch… Prodi… CDS   4.97e5 5.00e5 555.… -          0
    ## # ℹ 19,491 more rows
    ## # ℹ 37 more variables: attributes <chr>, conf <dbl>, cscore <dbl>,
    ## #   em_BRITE <chr>, em_BiGG_Reaction <chr>, em_CAZy <chr>, em_COG_cat <chr>,
    ## #   em_EC <chr>, em_GOs <chr>, em_KEGG_Module <chr>, em_KEGG_Pathway <chr>,
    ## #   em_KEGG_Reaction <chr>, em_KEGG_TC <chr>, em_KEGG_ko <chr>,
    ## #   em_KEGG_rclass <chr>, em_OGs <chr>, em_PFAMs <chr>,
    ## #   em_Preferred_name <chr>, em_desc <chr>, em_evalue <dbl>, …

``` r
cog_file <- 'cog_key.txt'

small_annotations <- subset(annotations, select = c('ID','em_Preferred_name','em_desc'))

# Load edgeR cleaned up output files #
# Pro Gene Expression: Control vs Monoxenic # 
input_gene_table_file_pro <- "comm_dyn_edgeR_Pro_v2_overview_gene_table.csv"

# Pro Gene Expression: Monoxenic vs Community # 
input_gene_table_file_pro_comm <- "comm_dyn_edgeR_pro_comm_normalized_overview_gene_table_v2.csv"

# Het Gene Expression: Community vs Monoxenic # 
input_gene_table_file_hets <- "comm_dyn_het_overview_gene_table_v2.csv"

# KEGG Gene Key #
kegg_gene_key <- 'comm_dyn_kegg_gene_key.csv'
```

##### 3. Heatmap Plot Functions

``` r
# Since there are several functions - I grouped the functions in one chunk so it is easier to run the remaining chunks #

##### 1. Generate Heatmaps on Input Data for Pro#####
prepare_heatmap_dfs_pro <- function(input_gene_table_file,experiment_day,new_organism) {
  
  # Load Gene Table #
  gene_table <- read.csv(file = input_gene_table_file,header = TRUE, check.names = FALSE)

  # Convert table to a Long Format #
  df_long <- gene_table %>%
    pivot_longer(
      cols = contains("_LFC") | contains("_pvalue")| contains("_padj"),
      names_to = c("treatment", ".value"),
      names_pattern = "(.*)_(.*)")
  
  # Obtain the genes that have are significant - as defined by LFC > 1.5 or LFC < -1.5 AND padj < 0.05 #
  df_long <- df_long %>%
    mutate(significance = ifelse(padj < 0.05 & abs(LFC) > 1.2, "significant", "not significant"))
  
  # Remove additional technical replicate #
  df_long<- df_long%>% filter(treatment != 'Thalassospira_1')

  # Rename samples for interpretability #
  df_long$treatment <- gsub("Thalassospira_2", "Thalassospira", df_long$treatment)
  df_long$organism <- gsub("Alteromonas_1904", "Alteromonas", df_long$organism)
  df_long$organism <- gsub("Marinobacter_1943", "Marinobacter", df_long$organism)
  df_long$organism <- gsub("Thalassospira_1907", "Thalassospira", df_long$organism)
  df_long$organism <- gsub("Pseudohoeflea_1940", "pseudohoeflea", df_long$organism)
  df_long$day <- gsub("2", "Day2", df_long$day)
  df_long$day <- gsub("4", "Day4", df_long$day)
  df_long$day <- gsub("5", "Day5", df_long$day)

  df_long <- df_long %>% mutate(updated_treatment_name = paste(treatment, day, sep = "_"))
  df_long<- df_long%>% filter(organism == new_organism)
  df_long
  
  # Define the significant genes accross each treatment
  significant_genes <- df_long %>%
    filter(significance == "significant") %>%  # Keep only significant genes
    select(ID,em_Preferred_name,em_desc) %>%                        # Select only the gene identifier column
    distinct()
  
  custom_order <- c("Marinobacter", "Thalassospira", "Pseudohoeflea", "Alteromonas", "Community")
  # Reorder the dataframe based on the custom order
  df_long_reordered <- df_long[order(match(df_long$treatment, custom_order)), ]
  df_long_reordered <- df_long_reordered[order(df_long_reordered$day), ]

  # Clean DF #
  df_long_reordered <- subset(df_long_reordered, select = -c(em_desc,day,organism,treatment,em_Preferred_name,pvalue,padj,significance))
  df_long_reordered

  # Convert DF to a wide format for Heatmap # 
  df_wide <- df_long_reordered %>%
    pivot_wider(names_from = updated_treatment_name,values_from = LFC)
  
  # Only keep genes with at least one significant gene across all treatments #
  df_wide<- df_wide%>%filter(ID%in%significant_genes$ID)%>%
    as.data.frame()
  rownames(df_wide) <- df_wide$ID
  
  # Remove Sample column #
  df_wide <- subset(df_wide, select = -c(ID))
  
  # Transpose DF for Heatmap #
  transposed_df <- t(df_wide)%>%
    as.data.frame()
  
  # Return Cleaned Tables #
  return(list(transposed_df = transposed_df, long_df = df_long))
}

##### 2. Generate Heatmaps on Input Data for Hets #####
prepare_heatmap_dfs_hets <- function(input_gene_table_file_hets,new_organism) {
  
  # Load Gene Table #
  het_gene_table <- read.csv(file = input_gene_table_file_hets,header = TRUE, check.names = FALSE)
  het_gene_table<- het_gene_table%>% filter(treatment == new_organism)

  # Convert table to a Long Format #
  df_long <- het_gene_table %>%
    pivot_longer(
      cols = contains("_LFC") | contains("_pvalue")| contains("_padj"),
      names_to = c("day", ".value"),
      names_pattern = "(.*)_(.*)")
  
  # Obtain the genes that have are significant - as defined by LFC > 1.5 or LFC < -1.5 AND padj < 0.05 #
  df_long <- df_long %>%
    mutate(significance = ifelse(padj < 0.05 & abs(LFC) > 1.5, "significant", "not significant"))
  
  # Define the significant genes accross each treatment
  significant_genes <- df_long %>%
    filter(significance == "significant") %>%  # Keep only significant genes
    select(ID,em_Preferred_name) %>%                        # Select only the gene identifier column
    distinct()

  # Clean DF #
  df_long_reordered <- df_long[order(df_long$day), ]
  df_long_reordered <- subset(df_long_reordered, select = -c(em_desc,treatment,em_Preferred_name,pvalue,padj,significance))

  # Convert DF to a wide format for Heatmap # 
  df_wide <- df_long_reordered %>%
    pivot_wider(names_from = day,values_from = LFC)
  
  # Only keep genes with at least one significant gene accross all treatments #
  df_wide<- df_wide%>%filter(ID%in%significant_genes$ID)%>%
    as.data.frame()
  rownames(df_wide) <- df_wide$ID
  
  # Remove Sample column #
  df_wide <- subset(df_wide, select = -c(ID))
  
  # Transpose DF for Heatmap #
  transposed_df <- t(df_wide)%>%
    as.data.frame()
  
  # Return Cleaned Tables #
  return(list(transposed_df = transposed_df, long_df = df_long))
}

##### Generate Heatmaps and Gene Clusters #####
generate_heatmaps <- function(transposed_df,num_clusters) {
  
  # Define color scheme #
  neg_colors <- colorRampPalette(c("blue", "white"))(25)  # From blue to white for negative values
  pos_colors <- colorRampPalette(c("white", "red"))(25)   # From white to red for positive values
  color_palette <- c(neg_colors, pos_colors[-1])
  
  # Can adjust colors #
  #breaks <- c(seq(-15, 0, length.out = 25), seq(0, 15, length.out = 24)[-1])
  #breaks <- c(seq(-20, 0, length.out = 25), seq(0, 20, length.out = 24)[-1])
  breaks <- c(seq(-5, 0, length.out = 25), seq(0, 5, length.out = 24)[-1])
  
  
  # Plot Heatmap using pheatmap package #
  heatmap <- pheatmap(transposed_df,cluster_rows = FALSE,show_colnames = FALSE,color = color_palette,
                    clustering_distance_rows = "euclidean", clustering_distance_cols = "euclidean",clustering_method = "ward.D2",
                    #clustering_distance_rows = "correlation", clustering_distance_cols = "correlation",clustering_method = "ward.D2",
                    breaks=breaks)
  
  # Extract row clustering information (dendrogram)
  gene_clusters <- cutree(heatmap$tree_col, k = num_clusters)  # k = 2 for two main clusters
  #swapped_gene_clusters<- ifelse(gene_clusters == 1, 2, 1)
  #gene_clusters

  # Convert col_clusters to a data frame for annotation
  col_annotation <- data.frame(Cluster = factor(gene_clusters))
  
  # Step 2: Adjust the factor levels to ensure Cluster 1 comes first, then 2, etc.
  #col_annotation$Cluster <- factor(col_annotation$Cluster, levels = sort(unique(gene_clusters)))

  # Set row names of annotation to match column names of the original matrix
  rownames(col_annotation) <- colnames(transposed_df)

  # Step 3: Plot the heatmap with column cluster annotations
  updated_heatmap <- pheatmap(transposed_df,cluster_rows = FALSE,show_colnames = FALSE,breaks = breaks,
                              annotation_col = col_annotation,cluster_cols = heatmap$tree_col,
                              annotation_names_col= FALSE,color = color_palette)

  return(list(updated_heatmap = updated_heatmap, gene_clusters = gene_clusters))
}


##### 3. Extract Gene Cluster Information #####
extract_cluster_information <- function(gene_clusters,num_clusters,annotation_file,cog_file,new_organism) {
  
  # Load Input Files #
  annotations <- read_excel(annotation_file)
  annotations$organism <- gsub("Alteromonas_1904", "Alteromonas", annotations$organism)
  annotations$organism <- gsub("Marinobacter_1943", "Marinobacter", annotations$organism)
  annotations$organism <- gsub("Thalassospira_1907", "Thalassospira", annotations$organism)
  annotations$organism <- gsub("Pseudohoeflea_1940", "pseudohoeflea", annotations$organism)
  annotations<- annotations%>% filter(organism == new_organism)
  cog_key <- read.csv(file = cog_file,sep = '\t')

  annotations$em_COG_cat <- str_replace_all(annotations$em_COG_cat, "None", "-")
  annotations <- annotations %>%
    mutate(em_COG_cat = strsplit(as.character(em_COG_cat), "")) %>%
    unnest(em_COG_cat) %>%
    as.data.frame()
  
  # Obtain Merged Annotation File #
  merged_annotations <- merge(annotations, cog_key, by = "em_COG_cat", all.x = TRUE)

  # Create empty list #
  filtered_genes <- list()

  # Loop through gene cluster values #
  for (i in 1:num_clusters) {
    # Obtain genes in cluster #
    cluster_genes <- names(gene_clusters[gene_clusters == i])
  
    # Filter df for only genes in that cluster #
    filtered_genes[[i]] <- merged_annotations %>% filter(ID %in% cluster_genes)
  
    # Add cluster name to df #
    filtered_genes[[i]]$cluster <- paste0("Cluster ", i)
  }

  # Combine all of the gene clusters together #
  combined_df <- bind_rows(filtered_genes)

  return(combined_df)
}

##### 4. Plot DEG COG Distribution #####
graph_overview_bar_chart <- function(cluster_output,cluster_colors) {
  # Overview Counts #
  cluster_COG_counts <- cluster_output %>%
    dplyr::group_by(updated_COG_name,cluster) %>%
    dplyr::summarise(counts = n()) %>% 
    as.data.frame()
  cluster_COG_counts <- cluster_COG_counts %>% mutate(counts = ifelse(cluster == "Cluster 1", counts * -1, counts))
  
  cluster_COG_counts$updated_COG_name <- fct_rev(cluster_COG_counts$updated_COG_name)
  # Plot #
  stacked_bar_cog_graph<-
    ggplot(cluster_COG_counts, aes(x=updated_COG_name,y=counts,fill=cluster))+
    geom_bar(stat = "identity", position = "identity", color = "black", width = 0.7) +
    geom_hline(yintercept = 0,color='black')+
    theme_classic()+
    theme(
      axis.text.x = element_text(family = "Helvetica",color='black',size=16),
      axis.text.y = element_text(family = "Helvetica",color='black',size=16),
      axis.title.x = element_text(family = "Helvetica",color='black',size=16),
      axis.title.y = element_blank(),
      strip.text = element_text(family = "Helvetica",color='black',size=16),
      legend.text = element_text(family = "Helvetica",size = 16),
      legend.title = element_text(family = "Helvetica",size = 16),
      panel.border=element_blank(),
      panel.background = element_rect(colour = "black", size=1),
      legend.position='right')+
    #scale_y_continuous(breaks = seq(-300, 200, 50),labels = abs,limits=c(-300,200))+
    scale_y_continuous(breaks = seq(-350, 300, 50),labels = abs,limits=c(-350,300))+
    scale_fill_manual(values = cluster_colors) +
    ylab('number of differentially expressed genes')+
    coord_flip()

  return(stacked_bar_cog_graph)
}


##### 5. Visualize Dot Plot for Pro #####
graph_overview_dotplot_pro <- function(long_df,new_organism,annotation_file,cog_file) {
  # Load Input Files #
  annotations<- annotations%>% filter(organism == new_organism)
  cog_key <- read.csv(file = cog_file,sep = '\t')
  annotations$em_COG_cat <- str_replace_all(annotations$em_COG_cat, "None", "-")
  annotations <- annotations %>%
    mutate(em_COG_cat = strsplit(as.character(em_COG_cat), "")) %>%
    unnest(em_COG_cat) %>%
    as.data.frame()
  
  # Obtain Merged Annotation File #
  merged_annotations <- merge(annotations, cog_key, by = "em_COG_cat", all.x = TRUE)
  merged_sig_annotations <- merge(merged_annotations, long_df, by = "ID")
  merged_sig_annotations<- merged_sig_annotations%>% filter(significance == 'significant')
  merged_sig_annotations

  # Overview Counts #
  grouped_sig_genes <- merged_sig_annotations %>%
    dplyr::group_by(treatment,updated_COG_name,day) %>%
    dplyr::summarise(counts = n(),
              .groups = 'drop') %>%
    as.data.frame()
  grouped_sig_genes
  
  grouped_sig_genes$updated_COG_name <- fct_rev(grouped_sig_genes$updated_COG_name)
  grouped_sig_genes$treatment = factor(grouped_sig_genes$treatment, levels=c('Marinobacter','Alteromonas','Thalassospira','Pseudohoeflea','Community'))
  dot_plot<- ggplot(grouped_sig_genes, aes(x = treatment, y = updated_COG_name, size=counts,color = updated_COG_name)) +
    geom_point(show.legend = TRUE) +
    scale_size_continuous(range = c(2, 15)) +
    theme_classic()+
    theme(
      axis.text.x = element_text(family = "Helvetica",color='black',size=16,angle = 90,hjust=0.9,vjust=0.5),
      axis.text.y = element_text(family = "Helvetica",color='black',size=16),
      axis.title.x = element_blank(),
      axis.title.y = element_blank(),
      strip.text = element_text(family = "Helvetica",color='black',size=16),
      legend.text = element_text(family = "Helvetica",size = 16),
      legend.title = element_text(family = "Helvetica",size = 16),
      panel.border=element_blank(),
      panel.background = element_rect(colour = "black", size=1),
      legend.position='right')+
    guides(size = guide_legend(title='number of significant genes')) +
    facet_grid(cols=vars(day),space="free")+
    scale_color_manual(values = cog_colors,guide='none')
  
  dot_plot

  return(dot_plot)
  
}
##### 6. Upset Plot for Pro #####
graph_upset_plot <- function(input_df) {
  input_df<- input_df%>% filter(treatment != 'community')
  input_df<- input_df%>% filter(significance == 'significant')
  input_df

  subsetted_df <- subset(input_df, select = c(ID,treatment,day))
  subsetted_df<- subsetted_df%>% filter(day == 'Day4')
  subsetted_df
  
  wide_df <- subsetted_df %>%
    mutate(ID = as.character(ID)) %>%
    mutate(presence = 1L) %>%
    pivot_wider(names_from = treatment, values_from = presence, values_fill = list(presence = 0L)) %>%
    # Ensure that all the presence/absence columns are integers
    mutate(across(where(is.numeric), as.integer)) %>% 
    as.data.frame()
  
  # Some categories have zero significant genes - if that is the case you can set these values here #
  wide_df$Community <- 0
  wide_df$Alteromonas <- 0
  
  wide_df$Community <- as.integer(wide_df$Community)
  wide_df$Alteromonas <- as.integer(wide_df$Alteromonas)
  wide_df
  
  upset_plot <-upset(wide_df, 
                    intersect = c("Thalassospira", "Pseudohoeflea", "Alteromonas", "Community","Marinobacter"),
                    width_ratio=0.2,
                    themes=upset_modify_themes(
                      list(
                        'intersections_matrix'=theme(axis.text.x = element_blank(),
                                                    axis.text.y = element_text(family = "Helvetica",color='black',size=12),
                                                    axis.title.x = element_blank(),
                                                    axis.title.y = element_blank()),
                        'overall_sizes'=theme(axis.text.y = element_blank(),
                                              axis.text.x = element_text(family = "Helvetica",color='black',size=8),
                                              axis.title.x = element_text(family = "Helvetica",color='black',size=12),
                                              axis.title.y = element_blank()),
                        'Intersection size'=theme(axis.text.y = element_text(family = "Helvetica",color='black',size=12),
                                              axis.text.x = element_blank(),
                                              axis.title.x = element_blank(),
                                              axis.title.y = element_text(family = "Helvetica",color='black',size=12)))))
  return(upset_plot)
  
}
```

##### 3a. Analyzing Signficant Genes for Pro in Control vs Monoxenic

``` r
##### 1. Run Preprocessing Step #####
# Define Input Files and Variables #
# This can be easily modified to focus data analysis on 1 day - just filter for experiment day = # of interest #
#experiment_day <- "2"
#experiment_day <- "4"
experiment_day <- "5"
new_organism <- "MED4"

# Run function command #
data_processing_output_pro <-prepare_heatmap_dfs_pro(input_gene_table_file_pro,experiment_day,new_organism)
heatmap_input_pro <- data_processing_output_pro$transposed_df
long_df_pro <- data_processing_output_pro$long_df
#heatmap_input_pro

##### 2. Generate Heatmaps and Gene Clusters #####
# Define number of clusters #
num_clusters <- 2

# Run function command #
heatmap_output_pro <-generate_heatmaps(heatmap_input_pro,num_clusters)
```

![](additional_files/unnamed-chunk-4-1.png)<!-- -->![](additional_files/unnamed-chunk-4-2.png)<!-- -->

``` r
updated_heatmap_pro <- heatmap_output_pro$updated_heatmap
gene_clusters_pro <- heatmap_output_pro$gene_clusters
#gene_clusters_pro
#updated_heatmap_pro

##### 3. Extract Gene Cluster Information #####
# Extract Cluster Outputs #
cluster_output_pro <-extract_cluster_information(gene_clusters_pro,num_clusters,annotation_file,cog_file,new_organism)
head(cluster_output_pro)
```

    ##   em_COG_cat            em_ID organism                 ID      seq_id
    ## 1          - MED4_genome_1283     MED4 1_1283_MED4_genome MED4_genome
    ## 2          -  MED4_genome_807     MED4  1_807_MED4_genome MED4_genome
    ## 3          - MED4_genome_1140     MED4 1_1140_MED4_genome MED4_genome
    ## 4          - MED4_genome_1598     MED4 1_1598_MED4_genome MED4_genome
    ## 5          - MED4_genome_1143     MED4 1_1143_MED4_genome MED4_genome
    ## 6          -  MED4_genome_364     MED4  1_364_MED4_genome MED4_genome
    ##            source type   start     end score strand phase
    ## 1 Prodigal_v2.6.3  CDS 1084749 1084931 16.58      -     0
    ## 2 Prodigal_v2.6.3  CDS  699665  699847 22.66      +     0
    ## 3 Prodigal_v2.6.3  CDS  963752  964045 11.52      +     0
    ## 4 Prodigal_v2.6.3  CDS 1348681 1348830  9.25      -     0
    ## 5 Prodigal_v2.6.3  CDS  965124  965303  7.77      +     0
    ## 6 Prodigal_v2.6.3  CDS  335406  335666 25.01      +     0
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                               attributes
    ## 1 ID=1_1283_MED4_genome;partial=00;start_type=ATG;rbs_motif=TAA;rbs_spacer=3bp;gc_cont=0.339;conf=97.84;score=16.58;cscore=10.14;sscore=6.44;rscore=3.17;uscore=2.04;tscore=1.89;em_ID=MED4_genome_1283;em_target=1471459.JFLJ01000130_gene1527;em_score=83.6;em_evalue=8.3e-21;em_tcov=100.0;em_OGs=2B2B7@1|root,31UVB@2|Bacteria,1GIQ3@1117|Cyanobacteria,1MMY7@1212|Prochloraceae;em_COG_cat=None;em_desc=None;em_max_annot_lvl=1117|Cyanobacteria;em_Preferred_name=
    ## 2           ID=1_807_MED4_genome;partial=00;start_type=ATG;rbs_motif=TAAAA;rbs_spacer=9bp;gc_cont=0.240;conf=99.46;score=22.66;cscore=17.27;sscore=5.40;rscore=3.41;uscore=0.10;tscore=1.89;em_ID=MED4_genome_807;em_target=146891.A9601_07931;em_score=84.3;em_evalue=4.11e-21;em_tcov=100.0;em_OGs=2A6WV@1|root,30VRR@2|Bacteria,1GI71@1117|Cyanobacteria,1MP3T@1212|Prochloraceae;em_COG_cat=None;em_desc=None;em_max_annot_lvl=1117|Cyanobacteria;em_Preferred_name=
    ## 3         ID=1_1140_MED4_genome;partial=00;start_type=ATG;rbs_motif=TAAAA;rbs_spacer=3bp;gc_cont=0.228;conf=93.39;score=11.52;cscore=1.14;sscore=10.38;rscore=5.04;uscore=2.72;tscore=2.63;em_ID=MED4_genome_1140;em_target=1501268.EW14_1142;em_score=163.0;em_evalue=2.03e-51;em_tcov=100.0;em_OGs=2B8YE@1|root,3228Y@2|Bacteria,1GN2T@1117|Cyanobacteria,1MNRJ@1212|Prochloraceae;em_COG_cat=None;em_desc=None;em_max_annot_lvl=1117|Cyanobacteria;em_Preferred_name=
    ## 4 ID=1_1598_MED4_genome;partial=00;start_type=ATG;rbs_motif=TAAA;rbs_spacer=13bp;gc_cont=0.267;conf=89.34;score=9.25;cscore=4.73;sscore=4.51;rscore=1.36;uscore=2.26;tscore=1.54;em_ID=MED4_genome_1598;em_target=1471459.JFLJ01000195_gene1335;em_score=65.9;em_evalue=4.32e-14;em_tcov=78.0;em_OGs=2A38G@1|root,30RQ8@2|Bacteria,1GQ22@1117|Cyanobacteria,1MPBY@1212|Prochloraceae;em_COG_cat=None;em_desc=None;em_max_annot_lvl=1117|Cyanobacteria;em_Preferred_name=
    ## 5             ID=1_1143_MED4_genome;partial=00;start_type=ATG;rbs_motif=AAA;rbs_spacer=4bp;gc_cont=0.278;conf=85.64;score=7.77;cscore=3.79;sscore=3.98;rscore=0.61;uscore=1.51;tscore=1.86;em_ID=MED4_genome_1143;em_target=167546.P9301_11111;em_score=94.7;em_evalue=2.91e-25;em_tcov=100.0;em_OGs=2A6MC@1|root,30VF6@2|Bacteria,1GI6B@1117|Cyanobacteria,1MMVF@1212|Prochloraceae;em_COG_cat=None;em_desc=None;em_max_annot_lvl=1117|Cyanobacteria;em_Preferred_name=
    ## 6           ID=1_364_MED4_genome;partial=00;start_type=TTG;rbs_motif=TAAAA;rbs_spacer=10bp;gc_cont=0.299;conf=99.68;score=25.01;cscore=25.92;sscore=-0.92;rscore=4.73;uscore=-0.64;tscore=-5.01;em_ID=MED4_genome_364;em_target=64471.sync_0955;em_score=54.3;em_evalue=2.26e-08;em_tcov=75.5;em_OGs=2B9HR@1|root,322VW@2|Bacteria,1GNHX@1117|Cyanobacteria,1H0ZF@1129|Synechococcus;em_COG_cat=None;em_desc=None;em_max_annot_lvl=1117|Cyanobacteria;em_Preferred_name=
    ##    conf cscore em_BRITE em_BiGG_Reaction em_CAZy em_EC em_GOs em_KEGG_Module
    ## 1 97.84  10.14     <NA>             <NA>    <NA>  <NA>   <NA>           <NA>
    ## 2 99.46  17.27     <NA>             <NA>    <NA>  <NA>   <NA>           <NA>
    ## 3 93.39   1.14     <NA>             <NA>    <NA>  <NA>   <NA>           <NA>
    ## 4 89.34   4.73     <NA>             <NA>    <NA>  <NA>   <NA>           <NA>
    ## 5 85.64   3.79     <NA>             <NA>    <NA>  <NA>   <NA>           <NA>
    ## 6 99.68  25.92     <NA>             <NA>    <NA>  <NA>   <NA>           <NA>
    ##   em_KEGG_Pathway em_KEGG_Reaction em_KEGG_TC em_KEGG_ko em_KEGG_rclass
    ## 1            <NA>             <NA>       <NA>       <NA>           <NA>
    ## 2            <NA>             <NA>       <NA>       <NA>           <NA>
    ## 3            <NA>             <NA>       <NA>       <NA>           <NA>
    ## 4            <NA>             <NA>       <NA>       <NA>           <NA>
    ## 5            <NA>             <NA>       <NA>       <NA>           <NA>
    ## 6            <NA>             <NA>       <NA>       <NA>           <NA>
    ##                                                                            em_OGs
    ## 1 2B2B7@1|root,31UVB@2|Bacteria,1GIQ3@1117|Cyanobacteria,1MMY7@1212|Prochloraceae
    ## 2 2A6WV@1|root,30VRR@2|Bacteria,1GI71@1117|Cyanobacteria,1MP3T@1212|Prochloraceae
    ## 3 2B8YE@1|root,3228Y@2|Bacteria,1GN2T@1117|Cyanobacteria,1MNRJ@1212|Prochloraceae
    ## 4 2A38G@1|root,30RQ8@2|Bacteria,1GQ22@1117|Cyanobacteria,1MPBY@1212|Prochloraceae
    ## 5 2A6MC@1|root,30VF6@2|Bacteria,1GI6B@1117|Cyanobacteria,1MMVF@1212|Prochloraceae
    ## 6 2B9HR@1|root,322VW@2|Bacteria,1GNHX@1117|Cyanobacteria,1H0ZF@1129|Synechococcus
    ##   em_PFAMs em_Preferred_name em_desc em_evalue   em_max_annot_lvl em_score
    ## 1     <NA>              <NA>    None  8.30e-21 1117|Cyanobacteria     83.6
    ## 2     <NA>              <NA>    None  4.11e-21 1117|Cyanobacteria     84.3
    ## 3     <NA>              <NA>    None  2.03e-51 1117|Cyanobacteria    163.0
    ## 4     <NA>              <NA>    None  4.32e-14 1117|Cyanobacteria     65.9
    ## 5     <NA>              <NA>    None  2.91e-25 1117|Cyanobacteria     94.7
    ## 6     <NA>              <NA>    None  2.26e-08 1117|Cyanobacteria     54.3
    ##                       em_target em_tcov gc_cont partial rbs_motif rbs_spacer
    ## 1 1471459.JFLJ01000130_gene1527   100.0   0.339       0       TAA        3bp
    ## 2            146891.A9601_07931   100.0   0.240       0     TAAAA        9bp
    ## 3             1501268.EW14_1142   100.0   0.228       0     TAAAA        3bp
    ## 4 1471459.JFLJ01000195_gene1335    78.0   0.267       0      TAAA       13bp
    ## 5            167546.P9301_11111   100.0   0.278       0       AAA        4bp
    ## 6               64471.sync_0955    75.5   0.299       0     TAAAA       10bp
    ##   rscore sscore start_type tscore uscore product kegg updated_gene_name contig
    ## 1   3.17   6.44        ATG   1.89   2.04      NA <NA>              <NA>     NA
    ## 2   3.41   5.40        ATG   1.89   0.10      NA <NA>              <NA>      1
    ## 3   5.04  10.38        ATG   2.63   2.72      NA <NA>              <NA>     NA
    ## 4   1.36   4.51        ATG   1.54   2.26      NA <NA>              <NA>     NA
    ## 5   0.61   3.98        ATG   1.86   1.51      NA <NA>              <NA>     NA
    ## 6   4.73  -0.92        TTG  -5.01  -0.64      NA <NA>              <NA>      1
    ##   updated_COG_name   cluster
    ## 1           No COG Cluster 1
    ## 2           No COG Cluster 1
    ## 3           No COG Cluster 1
    ## 4           No COG Cluster 1
    ## 5           No COG Cluster 1
    ## 6           No COG Cluster 1

``` r
##### 4. Graph Overview #####

cluster_colors <- c(`Cluster 1`= '#bf87fd',`Cluster 2`= '#900072',`Cluster 3`= '#66a557',
                    `Cluster 4`= '#ce595f',`Cluster 5`= '#7485c8',`Cluster 6`= '#bc8b3c',
                    `Cluster 7`= '#cb65a0',`Cluster 8`= '#47b19b',`Cluster 9`= '#725ac2')

# Extract Cluster Outputs #
bar_chart_pro <-graph_overview_bar_chart(cluster_output_pro,cluster_colors)
```

    ## `summarise()` has grouped output by 'updated_COG_name'. You can override using
    ## the `.groups` argument.

    ## Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
bar_chart_pro
```

![](additional_files/unnamed-chunk-4-3.png)<!-- -->

``` r
##### 5. Graph Overview #####

cog_colors <- c(`C (Energy production and conversion)`= '#A6CEE3',
                `D (Cell cycle control, cell division, chromosome partitioning)`= '#1F78B4',
                `E (Amino acid transport and metabolism)`='#B2DF8A',
                `F (Nucleotide transport and metabolism)`= '#33A02C',
                `G (Carbohydrate transport and metabolism)`= '#FB9A99',
                `H (Coenzyme transport and metabolism)`= '#E31A1C',
                `I (Lipid transport and metabolism)`= '#FDBF6F',
                `J (Translation,  ribosomal structure and biogenesis)`= '#FF7F00',
                `K (Transcription)`= '#CAB2D6',
                `L (Replication, recombination and repair)`= '#6A3D9A',
                `M (Cell wall/membrane/envelope biogenesis)`= '#FFFF99',
                `N (Cell motility)`= '#E6AB02',
                `O (Posttranslational modification, protein turnover, chaperones)`= '#8bc5c5',
                `P (Inorganic ion transport and metabolism)`= '#1B9E77',
                `Q (Secondary metabolites biosynthesis, transport and catabolism)`= '#dca7c2',
                `S (Function unknown)`='grey',
                `T (Signal transduction mechanisms)`= '#E7298A',
                `U (Intracellular trafficking, secretion, and vesicular transport)`= '#f2c36e',
                `V (Defense mechanisms)`= '#A6761D')

# Extract Cluster Outputs #
dotplot_chart_pro <-graph_overview_dotplot_pro(long_df_pro,new_organism,annotations,cog_file)
dotplot_chart_pro
```

![](additional_files/unnamed-chunk-4-4.png)<!-- -->

``` r
##### Upset Plot #####
upset_plot_pro <- graph_upset_plot(long_df_pro)
upset_plot_pro
```

![](additional_files/unnamed-chunk-4-5.png)<!-- -->

``` r
##### Save Files as SVGs #####
# Save Heatmap #
#ggsave("abundance_plot.png", plot = updated_heatmap_pro, width = 14, height = 8)

#Cairo(file = "abundance_plot.svg", type = "svg", width = 80, height = 40, units = "in", dpi = 350)
#print(updated_heatmap_pro)
#dev.off()
```

##### 3b. Analyzing Signficant Genes for Pro in Monoxenic vs Community

``` r
##### 1. Run Preprocessing Step #####
# This can be easily modified to focus data analysis on 1 day - just filter for experiment day = # of interest #
experiment_day <- "5"
new_organism <- "MED4"

# Run function command #
data_processing_output_pro_comm <-prepare_heatmap_dfs_pro(input_gene_table_file_pro_comm,experiment_day,new_organism)
heatmap_input_pro_comm <- data_processing_output_pro_comm$transposed_df
long_df_pro_comm <- data_processing_output_pro_comm$long_df

head(long_df_pro_comm)
```

    ## # A tibble: 6 × 11
    ##   ID    day   organism em_Preferred_name em_desc treatment    LFC  pvalue   padj
    ##   <chr> <chr> <chr>    <chr>             <chr>   <chr>      <dbl>   <dbl>  <dbl>
    ## 1 1_1_… Day5  MED4     dnaN              Confer… Marinoba…  0.247 0.167   0.532 
    ## 2 1_1_… Day5  MED4     dnaN              Confer… Alteromo…  0.391 0.00568 0.0612
    ## 3 1_1_… Day5  MED4     dnaN              Confer… Pseudoho…  0.178 0.215   0.583 
    ## 4 1_1_… Day5  MED4     dnaN              Confer… Thalasso…  0.177 0.311   0.730 
    ## 5 1_10… Day5  MED4     ftsY              Involv… Marinoba… -0.189 0.220   0.548 
    ## 6 1_10… Day5  MED4     ftsY              Involv… Alteromo… -0.199 0.201   0.390 
    ## # ℹ 2 more variables: significance <chr>, updated_treatment_name <chr>

``` r
#heatmap_input_pro_comm

##### 2. Generate Heatmaps and Gene Clusters #####
# Define number of clusters #
num_clusters <- 2

# Run function command #
heatmap_output_pro_comm <-generate_heatmaps(heatmap_input_pro_comm,num_clusters)
```

![](additional_files/unnamed-chunk-5-1.png)<!-- -->![](additional_files/unnamed-chunk-5-2.png)<!-- -->

``` r
updated_heatmap_pro_comm <- heatmap_output_pro_comm$updated_heatmap
gene_clusters_pro_comm <- heatmap_output_pro_comm$gene_clusters
#gene_clusters_pro
#updated_heatmap_pro_comm

##### 3. Extract Gene Cluster Information #####
# Extract Cluster Outputs #
cluster_output_pro_comm <-extract_cluster_information(gene_clusters_pro_comm,num_clusters,annotation_file,cog_file,new_organism)
head(cluster_output_pro_comm)
```

    ##   em_COG_cat            em_ID organism                 ID      seq_id
    ## 1          - MED4_genome_1598     MED4 1_1598_MED4_genome MED4_genome
    ## 2          - MED4_genome_1148     MED4 1_1148_MED4_genome MED4_genome
    ## 3          - MED4_genome_1607     MED4 1_1607_MED4_genome MED4_genome
    ## 4          -  MED4_genome_899     MED4  1_899_MED4_genome MED4_genome
    ## 5          -  MED4_genome_907     MED4  1_907_MED4_genome MED4_genome
    ## 6          - MED4_genome_1175     MED4 1_1175_MED4_genome MED4_genome
    ##            source type   start     end score strand phase
    ## 1 Prodigal_v2.6.3  CDS 1348681 1348830  9.25      -     0
    ## 2 Prodigal_v2.6.3  CDS  966541  966648  4.93      +     0
    ## 3 Prodigal_v2.6.3  CDS 1352535 1353047 70.20      -     0
    ## 4 Prodigal_v2.6.3  CDS  772787  772966 16.42      -     0
    ## 5 Prodigal_v2.6.3  CDS  774815  775090 23.94      -     0
    ## 6 Prodigal_v2.6.3  CDS  980848  981120 13.38      -     0
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                               attributes
    ## 1 ID=1_1598_MED4_genome;partial=00;start_type=ATG;rbs_motif=TAAA;rbs_spacer=13bp;gc_cont=0.267;conf=89.34;score=9.25;cscore=4.73;sscore=4.51;rscore=1.36;uscore=2.26;tscore=1.54;em_ID=MED4_genome_1598;em_target=1471459.JFLJ01000195_gene1335;em_score=65.9;em_evalue=4.32e-14;em_tcov=78.0;em_OGs=2A38G@1|root,30RQ8@2|Bacteria,1GQ22@1117|Cyanobacteria,1MPBY@1212|Prochloraceae;em_COG_cat=None;em_desc=None;em_max_annot_lvl=1117|Cyanobacteria;em_Preferred_name=
    ## 2                   ID=1_1148_MED4_genome;partial=00;start_type=ATG;rbs_motif=TAAAA;rbs_spacer=5bp;gc_cont=0.296;conf=75.66;score=4.93;cscore=0.11;sscore=4.82;rscore=1.99;uscore=1.73;tscore=1.10;em_ID=MED4_genome_1148;em_target=1471522.JFNU01000141_gene1689;em_score=63.2;em_evalue=2.05e-13;em_tcov=100.0;em_OGs=2B3VN@1|root,31WJI@2|Bacteria,1GJ4C@1117|Cyanobacteria,1MPE2@1212|Prochloraceae;em_COG_cat=None;em_desc=None;em_max_annot_lvl=1117|Cyanobacteria
    ## 3          ID=1_1607_MED4_genome;partial=00;start_type=ATG;rbs_motif=None;rbs_spacer=None;gc_cont=0.257;conf=100.00;score=70.20;cscore=68.94;sscore=1.26;rscore=-4.59;uscore=1.95;tscore=2.63;em_ID=MED4_genome_1607;em_target=59919.PMM1410;em_score=337.0;em_evalue=7.55e-118;em_tcov=100.0;em_OGs=2B8ZK@1|root,322A8@2|Bacteria,1GN3R@1117|Cyanobacteria,1MNTC@1212|Prochloraceae;em_COG_cat=None;em_desc=None;em_max_annot_lvl=1117|Cyanobacteria;em_Preferred_name=
    ## 4          ID=1_899_MED4_genome;partial=00;start_type=ATG;rbs_motif=TAA;rbs_spacer=11bp;gc_cont=0.317;conf=97.76;score=16.42;cscore=15.35;sscore=1.07;rscore=0.98;uscore=-1.77;tscore=1.86;em_ID=MED4_genome_899;em_target=74546.PMT9312_1507;em_score=112.0;em_evalue=1.95e-32;em_tcov=100.0;em_OGs=291A4@1|root,322F6@2|Bacteria,1GN72@1117|Cyanobacteria,1MNYW@1212|Prochloraceae;em_COG_cat=None;em_desc=None;em_max_annot_lvl=1117|Cyanobacteria;em_Preferred_name=
    ## 5         ID=1_907_MED4_genome;partial=00;start_type=TTG;rbs_motif=None;rbs_spacer=None;gc_cont=0.290;conf=99.59;score=23.94;cscore=32.03;sscore=-8.09;rscore=-4.59;uscore=2.16;tscore=-5.01;em_ID=MED4_genome_907;em_target=167546.P9301_15891;em_score=94.4;em_evalue=1.86e-24;em_tcov=97.0;em_OGs=2DJRR@1|root,30738@2|Bacteria,1GHIA@1117|Cyanobacteria,1MN68@1212|Prochloraceae;em_COG_cat=None;em_desc=None;em_max_annot_lvl=1117|Cyanobacteria;em_Preferred_name=
    ## 6            ID=1_1175_MED4_genome;partial=00;start_type=TTG;rbs_motif=TAA;rbs_spacer=6bp;gc_cont=0.227;conf=95.59;score=13.38;cscore=18.76;sscore=-5.38;rscore=0.58;uscore=-2.17;tscore=-5.01;em_ID=MED4_genome_1175;em_target=74545.EU96_0042;em_score=80.5;em_evalue=7.04e-19;em_tcov=98.7;em_OGs=2B8FP@1|root,321QS@2|Bacteria,1GNS6@1117|Cyanobacteria,1MMXH@1212|Prochloraceae;em_COG_cat=None;em_desc=None;em_max_annot_lvl=1117|Cyanobacteria;em_Preferred_name=
    ##     conf cscore em_BRITE em_BiGG_Reaction em_CAZy em_EC em_GOs em_KEGG_Module
    ## 1  89.34   4.73     <NA>             <NA>    <NA>  <NA>   <NA>           <NA>
    ## 2  75.66   0.11     <NA>             <NA>    <NA>  <NA>   <NA>           <NA>
    ## 3 100.00  68.94     <NA>             <NA>    <NA>  <NA>   <NA>           <NA>
    ## 4  97.76  15.35     <NA>             <NA>    <NA>  <NA>   <NA>           <NA>
    ## 5  99.59  32.03     <NA>             <NA>    <NA>  <NA>   <NA>           <NA>
    ## 6  95.59  18.76     <NA>             <NA>    <NA>  <NA>   <NA>           <NA>
    ##   em_KEGG_Pathway em_KEGG_Reaction em_KEGG_TC em_KEGG_ko em_KEGG_rclass
    ## 1            <NA>             <NA>       <NA>       <NA>           <NA>
    ## 2            <NA>             <NA>       <NA>       <NA>           <NA>
    ## 3            <NA>             <NA>       <NA>       <NA>           <NA>
    ## 4            <NA>             <NA>       <NA>       <NA>           <NA>
    ## 5            <NA>             <NA>       <NA>       <NA>           <NA>
    ## 6            <NA>             <NA>       <NA>       <NA>           <NA>
    ##                                                                            em_OGs
    ## 1 2A38G@1|root,30RQ8@2|Bacteria,1GQ22@1117|Cyanobacteria,1MPBY@1212|Prochloraceae
    ## 2 2B3VN@1|root,31WJI@2|Bacteria,1GJ4C@1117|Cyanobacteria,1MPE2@1212|Prochloraceae
    ## 3 2B8ZK@1|root,322A8@2|Bacteria,1GN3R@1117|Cyanobacteria,1MNTC@1212|Prochloraceae
    ## 4 291A4@1|root,322F6@2|Bacteria,1GN72@1117|Cyanobacteria,1MNYW@1212|Prochloraceae
    ## 5 2DJRR@1|root,30738@2|Bacteria,1GHIA@1117|Cyanobacteria,1MN68@1212|Prochloraceae
    ## 6 2B8FP@1|root,321QS@2|Bacteria,1GNS6@1117|Cyanobacteria,1MMXH@1212|Prochloraceae
    ##   em_PFAMs em_Preferred_name em_desc em_evalue   em_max_annot_lvl em_score
    ## 1     <NA>              <NA>    None  4.32e-14 1117|Cyanobacteria     65.9
    ## 2     <NA>              <NA>    None  2.05e-13 1117|Cyanobacteria     63.2
    ## 3     <NA>              <NA>    None 7.55e-118 1117|Cyanobacteria    337.0
    ## 4     <NA>              <NA>    None  1.95e-32 1117|Cyanobacteria    112.0
    ## 5     <NA>              <NA>    None  1.86e-24 1117|Cyanobacteria     94.4
    ## 6     <NA>              <NA>    None  7.04e-19 1117|Cyanobacteria     80.5
    ##                       em_target em_tcov gc_cont partial rbs_motif rbs_spacer
    ## 1 1471459.JFLJ01000195_gene1335    78.0   0.267       0      TAAA       13bp
    ## 2 1471522.JFNU01000141_gene1689   100.0   0.296       0     TAAAA        5bp
    ## 3                 59919.PMM1410   100.0   0.257       0      None       None
    ## 4            74546.PMT9312_1507   100.0   0.317       0       TAA       11bp
    ## 5            167546.P9301_15891    97.0   0.290       0      None       None
    ## 6               74545.EU96_0042    98.7   0.227       0       TAA        6bp
    ##   rscore sscore start_type tscore uscore product kegg updated_gene_name contig
    ## 1   1.36   4.51        ATG   1.54   2.26      NA <NA>              <NA>     NA
    ## 2   1.99   4.82        ATG   1.10   1.73      NA <NA>              <NA>     NA
    ## 3  -4.59   1.26        ATG   2.63   1.95      NA <NA>              <NA>     NA
    ## 4   0.98   1.07        ATG   1.86  -1.77      NA <NA>              <NA>      1
    ## 5  -4.59  -8.09        TTG  -5.01   2.16      NA <NA>              <NA>      1
    ## 6   0.58  -5.38        TTG  -5.01  -2.17      NA <NA>              <NA>     NA
    ##   updated_COG_name   cluster
    ## 1           No COG Cluster 1
    ## 2           No COG Cluster 1
    ## 3           No COG Cluster 1
    ## 4           No COG Cluster 1
    ## 5           No COG Cluster 1
    ## 6           No COG Cluster 1

##### 3c. Analyzing Signficant Genes for Hets

``` r
##### 1. Run Preprocessing Step #####

# Subset for heterotroph of interest #
het_organism <- "Marinobacter"
#het_organism <- "Alteromonas"
#het_organism <- "Thalassospira"
#het_organism <- "Pseudohoeflea"


# Run function command #
data_processing_output_hets <-prepare_heatmap_dfs_hets(input_gene_table_file_hets,het_organism)
#data_processing_output_hets

heatmap_input_hets <- data_processing_output_hets$transposed_df
long_df_hets <- data_processing_output_hets$long_df

#long_df_hets
#heatmap_input_hets

##### 2. Generate Heatmaps and Gene Clusters #####
# Define number of clusters - this can change as desired #
num_clusters <- 2

# Run function command #
heatmap_output_hets <-generate_heatmaps(heatmap_input_hets,num_clusters)
```

![](additional_files/unnamed-chunk-6-1.png)<!-- -->![](additional_files/unnamed-chunk-6-2.png)<!-- -->

``` r
updated_heatmap_hets <- heatmap_output_hets$updated_heatmap
gene_clusters_hets <- heatmap_output_hets$gene_clusters
head(gene_clusters_hets)
```

    ## 1_1062_batch10_bc1022_metaflye_ccs_contig_1_Marinobacter 
    ##                                                        1 
    ## 1_1063_batch10_bc1022_metaflye_ccs_contig_1_Marinobacter 
    ##                                                        1 
    ## 1_1102_batch10_bc1022_metaflye_ccs_contig_1_Marinobacter 
    ##                                                        1 
    ## 1_1132_batch10_bc1022_metaflye_ccs_contig_1_Marinobacter 
    ##                                                        1 
    ## 1_1133_batch10_bc1022_metaflye_ccs_contig_1_Marinobacter 
    ##                                                        1 
    ## 1_1160_batch10_bc1022_metaflye_ccs_contig_1_Marinobacter 
    ##                                                        1

``` r
#updated_heatmap_hets

#heatmap_output_hets

##### 3. Extract Gene Cluster Information #####

# Extract Cluster Outputs #
cluster_output_hets <-extract_cluster_information(gene_clusters_hets,num_clusters,annotation_file,cog_file,het_organism)
#cluster_output_hets

##### 4. Graph Overview #####

cluster_colors <- c(`Cluster 1`= '#bf87fd',`Cluster 2`= '#900072',`Cluster 3`= '#66a557',
                    `Cluster 4`= '#ce595f',`Cluster 5`= '#7485c8',`Cluster 6`= '#bc8b3c',
                    `Cluster 7`= '#cb65a0',`Cluster 8`= '#47b19b',`Cluster 9`= '#725ac2')

# Extract Cluster Outputs #
bar_chart_hets <-graph_overview_bar_chart(cluster_output_hets,cluster_colors)
```

    ## `summarise()` has grouped output by 'updated_COG_name'. You can override using
    ## the `.groups` argument.

``` r
bar_chart_hets
```

![](additional_files/unnamed-chunk-6-3.png)<!-- -->

##### 4. Pro Pathway Enrichment Analysis

``` r
##### 1. Extract Gene Cluster Information #####
cluster_profiler <- function(cluster_output,new_organism,category_show,cog_colors,kegg_gene_key_file,cluster_name) {
  kegg_gene_key <- read.csv(file = kegg_gene_key_file,header = TRUE, check.names = FALSE)
  
  names(kegg_gene_key)[names(kegg_gene_key) == "pathway_name"] <- "Description"
  
  subsetted_KEGG_pathway <- subset(cluster_output, select = c(kegg,cluster))
  names(subsetted_KEGG_pathway)[names(subsetted_KEGG_pathway) == "kegg"] <- "em_KEGG_ko"
  
  # Keep only first KEGG #
  #subsetted_KEGG_pathway <- subsetted_KEGG_pathway %>%
  #  separate(em_KEGG_ko, into = c("em_KEGG_ko"), sep = ",", remove = TRUE)
  cluster_ko_list <- split(subsetted_KEGG_pathway$em_KEGG_ko, subsetted_KEGG_pathway$cluster)
  cluster_ko_list
  
  # You can now use this list in compareCluster
  gene_enrich_result <- compareCluster(geneClusters = cluster_ko_list,fun = "enrichKEGG",organism = "ko")
  
  output_df <- as.data.frame(gene_enrich_result)
  #output_df<- output_df%>% filter(Cluster == 'Cluster 1')
  output_df<- output_df%>% filter(Cluster == cluster_name)
  output_df<- output_df%>% filter(Description != 'Carbon metabolism')
  output_df<- output_df%>% filter(Description != 'Biosynthesis of cofactors')
  output_df<- output_df%>% filter(Description != 'Biosynthesis of amino acids')
  #output_df<- output_df%>% filter(RichFactor >= 0.0)
  top_20_data <- output_df %>%arrange(desc(RichFactor)) %>%slice_head(n = 20)
  
  top_20_data$subcategory <- gsub("Global and overview maps", "Other", top_20_data$subcategory)
  top_20_data$subcategory <- gsub("Folding, sorting and degradation", "Other", top_20_data$subcategory)
  top_20_data$subcategory <- gsub("Metabolism of other amino acids", "Amino acid metabolism", top_20_data$subcategory)
  
  #merged_key_counts <- merge(output_df, kegg_gene_key, by = "Description")
  merged_key_counts <- left_join(top_20_data, kegg_gene_key, by = "Description")
  
  # Visualize the result
  gene_cluster_enrichment_analysis_plot <- ggplot(merged_key_counts, aes(x = RichFactor, y = reorder(as.factor(Description), RichFactor), color = functional_pathway)) +
    geom_point(size=6) +
    theme_classic()+
    ggtitle(new_organism) +
    labs(color = "KEGG Function")+
    scale_color_manual(values = cog_colors)+
    facet_grid(rows=vars(functional_pathway),scales = 'free_y',space = "free_y")+
    theme(
      axis.text.x = element_text(family = "Helvetica",color='black',size=16),
      axis.text.y = element_text(family = "Helvetica",color='black',size=16),
      axis.title.x = element_text(family = "Helvetica",color='black',size=16),
      axis.title.y = element_blank(),
      plot.title = element_text(hjust = 0.5,family = "Helvetica",color='black',size=16),
      strip.text = element_blank(),
      strip.background = element_blank(),
      legend.text = element_text(family = "Helvetica",size = 16),
      legend.title = element_text(family = "Helvetica",size = 16),
      panel.border=element_blank(),
      panel.background = element_rect(colour = "black", size=1),
      legend.position='right')+
    #scale_x_continuous(limits = c(0, 0.2))
    scale_x_continuous(limits = c(0, 0.7))
    
  return(gene_cluster_enrichment_analysis_plot)
}

cog_colors <- c(`Amino Acid Metabolism`= '#688e26',
                `Central Metabolism`= '#cf4d6f',
                `Cofactors and Vitamins`= '#4da1a9',
                `Fatty Acids`= '#ff8200',
                `Movement and Cell regulation`= '#ffc100',
                `Replication and repair`= '#a36d90',
                `Translation`='#b1a0cf',
                `Other`='#2660a4')

##### 1. Cluster Profiler for Pro#####
#temp_title_pro <- "Functional enrichment of gene clusters for"
#plot_title_pro <- paste(temp_title_pro, new_organism)
plot_title_pro <- "Top 20 Pathways Enriched in Xenic Community 1.2 LFC"
#plot_title_pro <- "Top 20 Pathways Enriched in Axenic Community 1.2 LFC"

# Choose the gene cluster you want to run pathway enrichment #
#cluster_name <- 'Cluster 2'
cluster_name <- 'Cluster 1'

# Number of Categories to show #
category_show_pro = 25

cluster_profiler_output_pro <- cluster_profiler(cluster_output_pro,plot_title_pro,category_show_pro,cog_colors,kegg_gene_key,cluster_name)
cluster_profiler_output_pro
```

![](additional_files/unnamed-chunk-7-1.png)<!-- -->
