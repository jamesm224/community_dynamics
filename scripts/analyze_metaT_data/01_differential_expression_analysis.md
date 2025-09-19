01_differential_expression_analysis
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
library(RColorBrewer)
library(readxl)
```

##### 2. Load Files

``` r
# Input RNASeq File - post internal standard correction
input_file <- "/Users/jamesmullet/updated_comm_dyn_wide_counts.tsv"

# Combined metadata key #
metadata_table <- "/Users/jamesmullet/pooled_key.csv"

# Annotations file of all genomes #
annotation_file <- '/Users/jamesmullet/updated_comm_dyn_annotations_v2.xlsx'
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
# Subsetted annotations df for easier parsing #
small_annotations <- subset(annotations, select = c('ID','em_Preferred_name','em_desc'))
```

##### 3. QC Analysis

``` r
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
```

    ## # A tibble: 6 × 5
    ##   ID                 experiment_type rep   day       count
    ##   <chr>              <fct>           <fct> <fct>     <dbl>
    ## 1 1_1000_MED4_genome 1904            A     DAY4   4352283.
    ## 2 1_1000_MED4_genome 1904            A     DAY5  11950035.
    ## 3 1_1000_MED4_genome 1904            B     DAY2  45623604.
    ## 4 1_1000_MED4_genome 1904            B     DAY4   6630337.
    ## 5 1_1000_MED4_genome 1904            B     DAY5  12140950.
    ## 6 1_1000_MED4_genome 1904            C     DAY2   2062564.

``` r
##### Generate Replicate Boxplot #####
boxplot <- plot_boxplot(data_processing_output_QC)
```

    ## Warning: The `size` argument of `element_rect()` is deprecated as of ggplot2 3.4.0.
    ## ℹ Please use the `linewidth` argument instead.
    ## This warning is displayed once every 8 hours.
    ## Call `lifecycle::last_lifecycle_warnings()` to see where this warning was
    ## generated.

``` r
boxplot
```

    ## Warning in scale_y_log10(): log-10 transformation introduced infinite values.

    ## Warning: Removed 574 rows containing non-finite outside the scale range
    ## (`stat_boxplot()`).

![](01_differential_expression_analysis_files/figure-gfm/unnamed-chunk-3-1.png)<!-- -->

``` r
##### Generate Replicate Correlation Plots #####
# For this plot you have to go to the function and comment out the different replicates for different comparisons #
correlation_plot <- plot_correlation_plot(data_processing_output_QC)
```

    ## `summarise()` has grouped output by 'comparison_group', 'day'. You can override
    ## using the `.groups` argument.
    ## `summarise()` has grouped output by 'day'. You can override using the `.groups`
    ## argument.

``` r
correlation_plot
```

    ## Warning: Removed 3 rows containing missing values or values outside the scale range
    ## (`geom_point()`).

![](01_differential_expression_analysis_files/figure-gfm/unnamed-chunk-3-2.png)<!-- -->

##### 4. edgeR analysis

``` r
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
```

    ##                                sample experimental_trial sampling_day replicate
    ## MED4ax_1943_A_DAY5 MED4ax_1943_A_DAY5       marinobacter            5         A
    ## MED4ax_1943_B_DAY5 MED4ax_1943_B_DAY5       marinobacter            5         B
    ## MED4ax_1943_C_DAY5 MED4ax_1943_C_DAY5       marinobacter            5         C
    ## MED4ax_all_A_DAY5   MED4ax_all_A_DAY5          community            5         A
    ## MED4ax_all_B_DAY5   MED4ax_all_B_DAY5          community            5         B
    ## MED4ax_all_C_DAY5   MED4ax_all_C_DAY5          community            5         C

``` r
# Output the raw count data - in appropriate format for edgeR #
cleaned_raw_counts <- data_processing_output$cleaned_raw_counts
head(cleaned_raw_counts)
```

    ##                    MED4ax_1943_A_DAY5 MED4ax_1943_B_DAY5 MED4ax_1943_C_DAY5
    ## 1_1000_MED4_genome          9544766.8           13422834           12147962
    ## 1_1001_MED4_genome          1115336.0            1988190            2934148
    ## 1_1002_MED4_genome        181756782.9          274302991          315748735
    ## 1_1003_MED4_genome        132184761.3          233537741          164847642
    ## 1_1004_MED4_genome         32930850.9           70731217           51039369
    ## 1_1005_MED4_genome           687695.4            2031139            1770901
    ##                    MED4ax_all_A_DAY5 MED4ax_all_B_DAY5 MED4ax_all_C_DAY5
    ## 1_1000_MED4_genome          16517089          14258406          16420943
    ## 1_1001_MED4_genome           3235673           3947291           3512896
    ## 1_1002_MED4_genome         400169467         414210020         443813531
    ## 1_1003_MED4_genome         300018285         291966483         334853886
    ## 1_1004_MED4_genome          78528316          75071373          81605911
    ## 1_1005_MED4_genome           3124687           2528810           2352655

``` r
##### Run EdgeR Step #####
edgeR_output <- run_edgeR(cleaned_metadata,cleaned_raw_counts,small_annotations,control)
head(edgeR_output)
```

    ##                   ID      logFC    logCPM         F    PValue       FDR
    ## 1    1_1_MED4_genome  0.2466676 6.2165134 2.3689835 0.1671253 0.5308681
    ## 2   1_10_MED4_genome -0.1887718 2.7934442 1.8152799 0.2193479 0.5476295
    ## 3  1_100_MED4_genome -0.1388609 4.6247835 1.3253216 0.2869828 0.5752072
    ## 4 1_1000_MED4_genome -0.1111369 2.8069912 0.6143547 0.4585357 0.6979771
    ## 5 1_1001_MED4_genome  0.3438132 0.4490681 1.6931528 0.2338773 0.5543051
    ## 6 1_1002_MED4_genome  0.1892703 7.3957442 1.6362080 0.2411213 0.5579744
    ##   em_Preferred_name
    ## 1              dnaN
    ## 2              ftsY
    ## 3               rbn
    ## 4              suhB
    ## 5              hisZ
    ## 6              htpG
    ##                                                                                                                                                                                                                                                                                                                                                                                                                                                                 em_desc
    ## 1 Confers DNA tethering and processivity to DNA polymerases and other proteins. Acts as a clamp, forming a ring around DNA (a reaction catalyzed by the clamp-loading complex) which diffuses in an ATP-independent manner freely and bidirectionally along dsDNA. Initially characterized for its ability to contact the catalytic subunit of DNA polymerase III (Pol III), a complex, multichain enzyme responsible for most of the replicative synthesis in bacteria
    ## 2                                                                                                                                                                                                                                               Involved in targeting and insertion of nascent membrane proteins into the cytoplasmic membrane. Acts as a receptor for the complex formed by the signal recognition particle (SRP) and the ribosome-nascent chain (RNC)
    ## 3                                                                                                                                                                                                                                                                                                                                                                                                                              Serum resistance locus BrkB-like protein
    ## 4                                                                                                                                                                                                                                                                                                                                                                                                                   Belongs to the inositol monophosphatase superfamily
    ## 5                                                                                                                                                                                                                                                                                                                       Required for the first step of histidine biosynthesis. May allow the feedback regulation of ATP phosphoribosyltransferase activity by histidine
    ## 6                                                                                                                                                                                                                                                                                                                                                                                                                                                Chaperone protein HtpG

``` r
##### Save output files #####
# Finalized gene tables were created using excel after running all edgeR data #
#write.csv(edgeR_output, "comm_dyn_edgeR_exp1_exp_analysis_Pro_day2.csv", row.names = FALSE)
```
