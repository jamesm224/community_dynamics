##### 1. Load Packages
library(ggplot2)
library(ggsci)
library(dplyr)
library(scales)
library(pheatmap)
library(viridis)
library(tidyplots)
library(readxl)

##### 2. Define Functions for QC Checks
### a. Plot Total DNA Reads in each sequenced sample ###
plot_total_dna_reads<- function(input_dna_reads_file) {
  
  # Load File #
  fastqc <- read.csv(file = input_dna_reads_file)
  fastqc$sequences_millions <- as.numeric(fastqc$sequences_millions)
  
  # Plot Data #
  fastqc_check<- ggplot(fastqc, aes(x=Sample.Name,y=sequences_millions)) + 
    geom_bar(stat="identity",position = 'dodge',color='white',fill='#4682B4')+
    theme_classic()+
    theme(axis.text.x = element_blank(),
          axis.text.y = element_text(family = "Helvetica",color='black',size=28),
          axis.title.x = element_blank(),
          axis.title.y = element_text(family = "Helvetica",color='black',size=28),
          legend.text = element_text(family = "Helvetica",size = 28),
          legend.title = element_text(family = "Helvetica",size = 28),
          strip.text = element_text(color='Black', size = 28),
          panel.border=element_blank(),
          legend.position='top')+
    ylab("Number of reads (million)")+
    geom_hline(yintercept = 1, color = "black", linetype = "solid", linewidth = 1) +
    scale_y_continuous(expand = c(0, 0))
  return(fastqc_check)
}
##### b. Plot number of internal DNA standards recovered #####
plot_dna_standards<- function(dna_standards_file) {
  # Load File #
  dna_standards <- read.csv(file = dna_standards_file)
  dna_standards_graph<- ggplot(dna_standards, aes(x=sample_name,y=`genome_equivalents`)) + 
    geom_bar(stat="identity",position = 'dodge',color='black',fill='#CBC3E3')+
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90,hjust=0.9,vjust=0.5,family = "Helvetica",color='black',size=16),
          axis.text.y = element_text(family = "Helvetica",color='black',size=22),
          axis.title.x = element_text(family = "Helvetica",color='black',size=22),
          axis.title.y = element_text(family = "Helvetica",color='black',size=22),
          legend.text = element_text(family = "Helvetica",size = 22),
          legend.title = element_text(family = "Helvetica",size = 22),
          strip.text = element_text(color='Black', size = 22),
          panel.border=element_blank(),
          legend.position='top')+
    xlab("")+
    ylab("genome equivalents")+
    scale_y_continuous(expand = c(0, 0))
  return(dna_standards_graph)
}

##### c. Plot number of internal DNA standards recovered #####
plot_metaG_FCM_comparison<- function(fcm_counts_file) {
  # Load File #
  fcm_counts <- read.csv(file = fcm_counts_file)
  
  fcm_scatter<-ggplot(fcm_counts, aes(x=new_name,y=het_ratio,fill=experiment,color=experiment))+
    geom_line(aes(group=experiment),size=1)+
    geom_point(size=3)+
    #geom_errorbar(aes(ymin = mean_absolute_counts - se_y, ymax = mean_absolute_counts + se_y), width = 0.06,size=0.5)+  # Error bar layer
    theme_classic()+
    theme(axis.text.x = element_text(angle = 90,hjust=0.9,vjust=0.5,family = "Helvetica",color='black',size=16),
          #axis.text.x = element_text(family = "Helvetica",color='black',size=16),
          axis.text.y = element_text(family = "Helvetica",color='black',size=22),
          axis.title.x = element_text(family = "Helvetica",color='black',size=22),
          axis.title.y = element_text(family = "Helvetica",color='black',size=22),
          legend.text = element_text(family = "Helvetica",size = 18),
          legend.title = element_text(family = "Helvetica",size = 18),
          strip.text = element_text(color='Black', size = 22),
          panel.border=element_blank(),
          legend.position='top')+
    xlab("")+
    ylab("relative abundance of heterotrophs")
    #facet_grid(cols=vars(experiment_day),scales="free",space="free")
    #scale_color_manual(values = alpha(c("#ff3145","#eb91e3","#4f2c7f",
    #                                    "#007bc0","#00e674","#ff7f00","#cbbbc2")))+
    #scale_y_log10(limits = c(50000, 100000000))+
    #annotation_logticks(sides="l")
  return(fcm_scatter)
}

##### 3. Generate QC Plots
##### a. Define Input Files #####
input_dna_reads_file='/Users/jamesmullet/dna_community_dynamics_data.csv'
dna_standards_file='/Users/jamesmullet/dna_standard_reads_mapped.csv'
fcm_counts_file = '/Users/jamesmullet/fcm_metaG_comparison.csv'

##### b. Generate Total DNA Reads Plot #####
dna_counts_plot<- plot_total_dna_reads(input_dna_reads_file)
dna_counts_plot

##### c. Generate Number of Recovered DNA Standards #####
internal_standard_plot <- plot_dna_standards(dna_standards_file)
internal_standard_plot

##### d. FCM vs metaG Comparison #####
fcm_metaG_plot <- plot_metaG_FCM_comparison(fcm_counts_file)
fcm_metaG_plot
