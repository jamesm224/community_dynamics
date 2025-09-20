---
title: "02_processing_relative_RNASeq_data"
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
library(readxl)
library(RColorBrewer)
```

##### 2. Load Files
```{r}
# Annotations file of all genomes #
annotation_file <- 'updated_comm_dyn_annotations_v2.xlsx'
annotations <- read_excel(annotation_file)
annotations

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
```{r}
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
```{r,fig.width=16, fig.height=8}

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
updated_heatmap_pro <- heatmap_output_pro$updated_heatmap
gene_clusters_pro <- heatmap_output_pro$gene_clusters
#gene_clusters_pro
#updated_heatmap_pro

##### 3. Extract Gene Cluster Information #####
# Extract Cluster Outputs #
cluster_output_pro <-extract_cluster_information(gene_clusters_pro,num_clusters,annotation_file,cog_file,new_organism)
head(cluster_output_pro)

##### 4. Graph Overview #####

cluster_colors <- c(`Cluster 1`= '#bf87fd',`Cluster 2`= '#900072',`Cluster 3`= '#66a557',
                    `Cluster 4`= '#ce595f',`Cluster 5`= '#7485c8',`Cluster 6`= '#bc8b3c',
                    `Cluster 7`= '#cb65a0',`Cluster 8`= '#47b19b',`Cluster 9`= '#725ac2')

# Extract Cluster Outputs #
bar_chart_pro <-graph_overview_bar_chart(cluster_output_pro,cluster_colors)
bar_chart_pro

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

##### Upset Plot #####
upset_plot_pro <- graph_upset_plot(long_df_pro)
upset_plot_pro

##### Save Files as SVGs #####
# Save Heatmap #
#ggsave("abundance_plot.png", plot = updated_heatmap_pro, width = 14, height = 8)

#Cairo(file = "abundance_plot.svg", type = "svg", width = 80, height = 40, units = "in", dpi = 350)
#print(updated_heatmap_pro)
#dev.off()

```


##### 3b. Analyzing Signficant Genes for Pro in Monoxenic vs Community
```{r,fig.width=16, fig.height=8}

##### 1. Run Preprocessing Step #####
# This can be easily modified to focus data analysis on 1 day - just filter for experiment day = # of interest #
experiment_day <- "5"
new_organism <- "MED4"

# Run function command #
data_processing_output_pro_comm <-prepare_heatmap_dfs_pro(input_gene_table_file_pro_comm,experiment_day,new_organism)
heatmap_input_pro_comm <- data_processing_output_pro_comm$transposed_df
long_df_pro_comm <- data_processing_output_pro_comm$long_df

head(long_df_pro_comm)
#heatmap_input_pro_comm

##### 2. Generate Heatmaps and Gene Clusters #####
# Define number of clusters #
num_clusters <- 2

# Run function command #
heatmap_output_pro_comm <-generate_heatmaps(heatmap_input_pro_comm,num_clusters)
updated_heatmap_pro_comm <- heatmap_output_pro_comm$updated_heatmap
gene_clusters_pro_comm <- heatmap_output_pro_comm$gene_clusters
#gene_clusters_pro
#updated_heatmap_pro_comm

##### 3. Extract Gene Cluster Information #####
# Extract Cluster Outputs #
cluster_output_pro_comm <-extract_cluster_information(gene_clusters_pro_comm,num_clusters,annotation_file,cog_file,new_organism)
head(cluster_output_pro_comm)

```


##### 3c. Analyzing Signficant Genes for Hets
```{r,fig.width=16, fig.height=8}

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
updated_heatmap_hets <- heatmap_output_hets$updated_heatmap
gene_clusters_hets <- heatmap_output_hets$gene_clusters
head(gene_clusters_hets)
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
bar_chart_hets

```

##### 4. Pro Pathway Enrichment Analysis
```{r,fig.width=12, fig.height=8}
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
