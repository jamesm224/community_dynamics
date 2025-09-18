from pathlib import Path

wd = Path("/nfs/chisholmlab001/kve/2023_HIFI_Genome_Closing_Project/assembly")

flye_paths = wd.glob("batch10*/*/*flye_*/assembly.fasta")
canu_paths = wd.glob("batch10*/*/canu_*/assembly.contigs.fasta")

assemblies = []
for group in [flye_paths, canu_paths]:
    assemblies.extend([
        {
            "batch": path.parent.parent.parent.name,
            "barcode": path.parent.parent.name,
            "assembly method": path.parent.name,
            "path": path,
        }
        for path in group
    ])

for assembly in assemblies:
    new_file = wd / 'all_assemblies' / f"{assembly['batch']}.{assembly['barcode']}.{assembly['assembly method']}.fasta"
    new_file.hardlink_to(assembly['path'])