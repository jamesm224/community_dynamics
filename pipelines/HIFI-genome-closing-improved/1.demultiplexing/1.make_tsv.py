from pathlib import Path 

# path to raw sequencing dir from BMC 
raw_seq_dir = "/nfs/chisholmlab001/chisholmlab/experiment_repository/2023/230823ChiA/PacBio/1_A01"
 
bc_data = {}
for fpath in Path(raw_seq_dir).glob('bc*/m64408e_231129_202258*.bam'):
    batch = fpath.stem.split('--')[1] 
    bc_data[batch] = fpath

with open("demultiplex.tsv", "w") as file:
    for batch, fpath in bc_data.items():
        file.write(f"{batch}\t{fpath}\n")
