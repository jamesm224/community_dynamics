#!/usr/bin/env bash
#SBATCH --job-name=assem
#SBATCH --time 1-0                         
#SBATCH -p sched_mit_chisholm              
#SBATCH -c 20                                
#SBATCH -n 1
#SBATCH -N 1
#SBATCH --mem 250G
#SBATCH --array=1-8%8
#SBATCH -o logs/assembly.%a.out
#SBATCH -e logs/assembly.%a.err

date

out_dir="/nfs/chisholmlab001/kve/2023_HIFI_Genome_Closing_Project/assembly"

batch='batch10'
barcode=$(sed "${SLURM_ARRAY_TASK_ID}q;d" assembly.tsv | cut -f1 | tr -d '\r')
ccs_fastq_path=$(sed "${SLURM_ARRAY_TASK_ID}q;d" assembly.tsv | cut -f2 | tr -d '\r')  # ccs is hifi 

# ccs_fastq_path=$(sed "${SLURM_ARRAY_TASK_ID}q;d" assembly.tsv | cut -f3 | tr -d '\r')
# subreadset_fastq_path=$(sed "${SLURM_ARRAY_TASK_ID}q;d" assembly.tsv | cut -f4 | tr -d '\r')

echo ${batch}
echo ${barcode}
echo ${hifi_fastq_path}

flye_ccs_outdir="${out_dir}/${batch}/${barcode}/flye_ccs"
metaflye_ccs_outdir="${out_dir}/${batch}/${barcode}/metaflye_ccs"
canu_ccs_outdir="${out_dir}/${batch}/${barcode}/canu_ccs"
spades_ccs_outdir="${out_dir}/${batch}/${barcode}/spades_ccs"
circlator_ccs_outdir="${out_dir}/${batch}/${barcode}/circlator_ccs"

mkdir -p ${flye_ccs_outdir}
mkdir -p ${metaflye_ccs_outdir}
mkdir -p ${canu_ccs_outdir}
mkdir -p ${spades_ccs_outdir}


### flye and metaflye ### 
echo "flye assembly on CCS"

eval "$(conda shell.bash hook)"
conda activate flye
date
flye --version
flye --pacbio-hifi ${ccs_fastq_path} --out-dir ${flye_ccs_outdir} -t ${SLURM_CPUS_PER_TASK}

echo "metaflye assembly on CCS"
date
flye --version
flye --meta --pacbio-hifi ${ccs_fastq_path} --out-dir ${metaflye_ccs_outdir} -t ${SLURM_CPUS_PER_TASK}

source deactivate


### canu ### 
echo "canu assembly on CCS"
date
/home/kve/bin/canu-2.2/bin/canu --version
/home/kve/bin/canu-2.2/bin/canu \
   -p assembly \
   -d ${canu_ccs_outdir} \
   genomeSize=1.5m \
   -pacbio-hifi ${ccs_fastq_path} \
   useGrid=false \
   batThreads="${SLURM_CPUS_PER_TASK}"

### spades and circlator ### 
echo "SPAdes assembly on CCS"

echo "circlator on spades output"
conda activate spades

date
spades.py --version
spades.py -t ${SLURM_CPUS_PER_TASK} -s ${ccs_fastq_path} -o ${spades_ccs_outdir}

source deactivate

echo "circlator on spades output"
conda activate circlator
date
circlator version
circlator all --threads ${SLURM_CPUS_PER_TASK} ${spades_ccs_outdir}/scaffolds.fasta ${ccs_fastq_path} ${circlator_ccs_outdir}

date