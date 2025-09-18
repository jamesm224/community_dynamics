#!/usr/bin/env bash

# purpose: to extract assemblies from genome closing results 

eval "$(conda shell.bash hook)"  
conda activate seqkit

input=data/contig_data_filtered.tsv 

outdir=data/extracted_genomes
all_genome_outdir=data/all_extracted_genomes
mkdir -p ${outdir}
mkdir -p ${all_genome_outdir}

 {
    read;  # ignore/skip the first line (column names)
    while read line; do
        # obtain data from table 
        batch=`echo "$line" | cut -d$'\t' -f1 | tr -d '\r'`
        barcode=`echo "$line" | cut -d$'\t' -f2 | tr -d '\r'`
        assembly_method=`echo "$line" | cut -d$'\t' -f3 | tr -d '\r'`
        contig=`echo "$line" | cut -d$'\t' -f4 | xargs`
        assembly_path=`echo "$line" | cut -d$'\t' -f7 | xargs`
        classification=`echo "$line" | cut -d$'\t' -f10 | xargs`  # xargs = better method to remove leading + trailing white spaces 
        

        outfile=${outdir}/${batch}/${barcode}/${batch}_${barcode}_${assembly_method}_${contig}.fna

        echo ${batch}_${barcode}_${assembly_method}_${contig}.fna
        seqkit grep -p ${contig} ${assembly_path} -o ${outfile}  # obtain contig from assembly 
        cp ${outfile} ${all_genome_outdir}  # copy to dir with all files
    done
 } < ${input}
