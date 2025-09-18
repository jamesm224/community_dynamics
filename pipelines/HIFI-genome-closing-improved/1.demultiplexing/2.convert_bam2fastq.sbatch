#!/usr/bin/env bash
#SBATCH --job-name=bam2fastq
#SBATCH --time 1-0                         
#SBATCH -p sched_mit_chisholm              
#SBATCH -c 1                               
#SBATCH -N 1
#SBATCH --mem 10G
#SBATCH --array=1-8%8
#SBATCH -o logs/bam2fastq.%a.out
#SBATCH -e logs/bam2fastq.%a.err

out_dir="/nfs/chisholmlab001/kve/2023_HIFI_Genome_Closing_Project/demultiplexing/all"

batch=$(sed "${SLURM_ARRAY_TASK_ID}q;d" demultiplex.tsv | cut -f1 | tr -d '\r')
fpath=$(sed "${SLURM_ARRAY_TASK_ID}q;d" demultiplex.tsv | cut -f2 | tr -d '\r')

echo "convert hifi bams to fastqs"

eval "$(conda shell.bash hook)"
conda activate samtools

samtools --version

outfile="${out_dir}/batch10.1.${batch}.hifi_readset.fastq"
samtools fastq ${fpath} > ${outfile}
