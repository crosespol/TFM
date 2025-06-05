from Bio import SeqIO
from pathlib import Path
import re
import csv
ref_phrogs = Path("/home/carlesroses/PHROGs/Reference_Phrogs/Reference_PHROGs_c1_cro_0.csv")
Genomes = Path("/home/carlesroses/PHROGs/GenomesDB")
phrog_proteins = [] # list of dictionaries
output_file = Path("proteins_c1_cro_1st.csv")
Total_genomes = 0
processed_genomes = 0
#genomes_without_gbf = 0
#genomes_processed_from_gff = 0

with ref_phrogs.open() as f:
    reader = csv.DictReader(f)
    reference_phrogs = {row['phrog_ID'] for row in reader}

print(f"✅ Loaded {len(reference_phrogs)} reference PHROGs")

for folder in Genomes.iterdir():
    if folder.is_dir():
        genome_id = folder.name
        Total_genomes += 1
        # print(f"Processing genome: {genome_id}")
        
        '''
        gbf_path = folder / f"{genome_id}.gbf"
        if gbf_path.exists():
            processed_genomes += 1
            # print(f"Parsing {gbf_path}...")
            with gbf_path.open() as gbf_file:
                for record in SeqIO.parse(gbf_file, "genbank"):
                    # print(f"Genome length: {len(record.seq)}")
                    # print(f"Number of features: {len(record.features)}")
                    for feature in record.features:
                        if feature.type == "CDS":
                            qualifiers = feature.qualifiers
                            inference_text = " ".join(qualifiers.get("inference", []))
                            match = re.search(r"phrog_\d+", inference_text)
                            if match:
                                phrog_id = match.group() #returns the object of the match, which in this case is the phrog_xxx
                                if phrog_id not in reference_phrogs:
                                    continue
                                protein_id = qualifiers.get("protein_id", qualifiers.get("locus_tag", ["unknown"]))[0]
                                start = int(feature.location.start)
                                end = int(feature.location.end)
                                strand = "+" if feature.location.strand == 1 else "-"
                                aa_seq = qualifiers.get("translation", [""])[0]

                                phrog_proteins.append({
                                    'protein_id': protein_id,
                                    'phrog_id': phrog_id,
                                    'genome': genome_id,
                                    'start': start,
                                    'end': end,
                                    'strand': strand,
                                    'aa_seq': aa_seq
                                })
                                print(f"Total genomes: {Total_genomes} | Processed genomes: {processed_genomes} | Missing genomes: {Total_genomes - processed_genomes} | ⚠️ Genomes without .gbf: {genomes_without_gbf} | Genomes from .gff {genomes_processed_from_gff}", end="\r")
'''
        #else:
        # print(f"⚠️ No .gbf file found for genome {genome_id}")
        #genomes_without_gbf +=1
        gff_path = folder / f"{genome_id}.gff"
        faa_path = folder / f"{genome_id}.faa"
        if gff_path.exists() and faa_path.exists():
            processed_genomes += 1
            #Carreguem totes les sequencies de proteines per després buscarles desde les proteins ID que treiem del gff.
            faa_proteins = {
                record.id.split()[0]: str(record.seq) #el que fem aquí es guardar el protein ID junt amb la sequencia de la proteina. Els headers dels faa, norlament tenen més info de la proteina que volem evitar
                for record in SeqIO.parse(faa_path, "fasta") # amb "fasta" li diem a biopython que el archiu esta en format fasta
            }

            # Ara que tenim totes les proteines del archiu faa carregades, anem a buscar les que volem dins del archiu gff. on agafarm les propteines que volem dels nostres PHROGs de referencia
            with gff_path.open() as gff_file:
                for line in gff_file:
                    if line.startswith("#") or "\tCDS\t" not in line: # el format dels gff sempre comenÇa amb comentaris i despres la resta. Hem de buscar lineas que siguin CDS (protein coding genes) que es el que ens interessa, i després d'aquesta linea extreurem el locus_tag
                        continue
                    fields = line.strip().split("\t") #Ho fem servir per accedir a continuació a cada camp de la linea
                    start = int(fields[3])
                    end = int(fields[4])
                    strand = fields[6]
                    atributes = fields[8]
                    #Aquí hem accedit i extret la info de cada una de las columnas del gff file.
                    #Ara creem un diccionari temporal amb tuples, que es reiniciara cadda cop que accedim a una altre gff file.
                    attrs = dict(
                        item.split("=") for item in atributes.split(";") if "=" in item
                    )
                    locus_tag = attrs.get('locus_tag', 'unknown')
                    inference = attrs.get('inference', '')
                    match = re.search(r'phrog_\d+', inference) #busquem dintre del camp the inference si hi ha un phrog_XXXX
                    if not match:
                        continue #Els continues skipegen tot el codi de sota i tornen al look a fer check de la seguent linea
                    phrog_id = match.group()
                    if phrog_id not in reference_phrogs:
                        continue

                    aa_seq = faa_proteins.get(locus_tag, '')
                    if not aa_seq:
                        continue

                    phrog_proteins.append({
                        'protein_id': locus_tag,
                        'phrog_id': phrog_id,
                        'genome': genome_id,
                        'start': start,
                        'end': end,
                        'strand': strand,
                        'aa_seq': aa_seq
                    })
        #genomes_processed_from_gff += 1
        print(f"Total genomes: {Total_genomes} | Processed genomes: {processed_genomes} | Missing genomes: {Total_genomes - processed_genomes}", end="\r")

print()  # newline after the \r printing
with output_file.open("w", newline="") as csvfile:
    fieldnames = ['protein_id', 'phrog_id', 'genome', 'start', 'end', 'strand', 'aa_seq']
    writer = csv.DictWriter(csvfile, fieldnames=fieldnames)

    writer.writeheader()
    for protein in phrog_proteins:
        writer.writerow(protein)

print(f"✅ Saved {len(phrog_proteins)} entries to: {output_file}")

fasta_path = Path("reference_proteins.faa")
with fasta_path.open("w") as fasta_file:
    for entry in phrog_proteins:
        header = f">{entry['protein_id']}_{entry['phrog_id']}_{entry['genome']}"
        sequence = entry['aa_seq']
        fasta_file.write(f"{header}\n{sequence}\n")

print(f"✅ FASTA file written: {fasta_path}")

'''
def grouping_prot_phrog(proteins_csv):
    input_csv = Path(proteins_csv)
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

grouping_prot_phrog(output_file)
'''