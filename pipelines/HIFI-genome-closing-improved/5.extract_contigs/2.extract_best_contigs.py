"""
To extract the best contig for a barcode/sample. 
    "best" assemblies = metaflye assemblies
    metaflye > flye > circlator > canu 
    ccs > subreadset
"""

import shutil
import os
from pathlib import Path
import pandas as pd


def main():
    """
    """
    assembly_info = 'data/contig_data_filtered.tsv'
    extracted_contig_path = 'data/all_extracted_genomes'
    best_genome_outdir = 'data/best_asesmblies'
    Path(best_genome_outdir).mkdir(parents=True, exist_ok=True)

    df = pd.read_csv(assembly_info, sep='\t')
    # remove rows without classification 
    df = df.dropna(subset=['taxonomic label'])

    order_mapping = {'metaflye_ccs':0,
                'flye_ccs':1,
                'metaflye_subreadset':2, 
                'flye_subreadset':3, 
                'circlator_ccs':4, 
                'circlator_subreadset':5, 
                'canu_ccs':6, 
                'canu_subreadset':7}

    # assign ranking based on assembly method (for sorting and obtaining 'best' assembly)
    df['ranking']=df['assembly method'].map(order_mapping)
    df['classification'] = df['taxonomic label'].apply(lambda x: x.strip())  # remove trailing + leading white space, for better grouping 

    # obtain best assembly within each group
    dfs=[]
    groups=df.groupby(['batch', 'barcode', 'classification'])
    for group_name, gdf in groups:
        gdf=gdf.sort_values(by=['ranking'])
        dfs.append(gdf.head(1))  # top based on assembly method rating 
    df=pd.concat(dfs)

    for index, row in df.iterrows():
        # copy over best assembly 
        row=row.str.strip()
        fname=(f"{row['batch']}_{row['barcode']}_{row['assembly method']}_{row['contig']}.fna")
        source=Path(f"{extracted_contig_path}/{fname}")
        destination=best_genome_outdir
        try: shutil.copy(source, destination)
        except Exception as e: print(e)

        # rename file
        classification=row['classification'].replace(' ', '_')
        old_name=f"{best_genome_outdir}/{fname}"
        new_name=f"{best_genome_outdir}/{source.stem}_{classification}.fasta"
        try: os.rename(old_name, new_name)
        except Exception as e: print(e)



if __name__ == "__main__":
    main()