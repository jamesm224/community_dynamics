from pathlib import Path 

demultiplexed_outdir = "/nfs/chisholmlab001/kve/2023_HIFI_Genome_Closing_Project/demultiplexing/all"

demultiplexed_data = {}
for batch10_path in Path(demultiplexed_outdir).glob('batch10*'):
    barcode = batch10_path.stem.split('.')[2] 
    demultiplexed_data[barcode] = batch10_path

with open('assembly.tsv', 'w') as file:
    for batch, fpath in demultiplexed_data.items():
        file.write(f"{barcode}\t{fpath}\n")