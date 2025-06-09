from pathlib import Path
import csv

input_csv = Path("proteins_c1_cro_1st.csv")
output_dir = Path("phrog_fastas")
output_dir.mkdir(exist_ok=True)

phrog_groups = {}

with input_csv.open() as f: 
    reader = csv.DictReader(f) #Aquest commando llegeix directament el csv creant un diccionari i per tant, linea a linea busca el phrog_id i mira si esta ja ficat en el diccionary,
    for row in reader:
        phrog_id = row["phrog_id"]
        if phrog_id not in phrog_groups:
            phrog_groups[phrog_id] = []
        phrog_groups[phrog_id].append(row)

#Writting one fasta file per phrog
for phrog_id, proteins in phrog_groups.items():
    fasta_path = output_dir / f"{phrog_id}.faa"
    with fasta_path.open("w") as fasta_file:
        for protein in proteins:
            header = f">{protein['protein_id']}_{protein['genome']}"
            sequence = protein['aa_seq']
            fasta_file.write(f"{header}\n{sequence}\n")
print(f"✅ Generated {len(phrog_groups)} FASTA files in {output_dir}")
