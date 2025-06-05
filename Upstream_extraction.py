from pathlib import Path
import re
import csv
from Bio import SeqIO
from collections import defaultdict
clusters_path = Path("/home/carlesroses/PHROGs/usearch_results_byPHROG/clusters_0.8/phrog_4_clusters.uc")
#cluster_proteins = {}
genome_to_prots = defaultdict(list)
Genomes = Path("/home/carlesroses/PHROGs/GenomesDB")
total_proteins = 0
written_sequences = 0

with clusters_path.open() as f: #Canviem la manera per agilitzar el proces
    for line in f:
        if line.startswith(("S", "H")):
            fields = line.strip().split("\t")
            cluster_id = fields[1]
            prot_id = fields[8]

            # Extract genome from prot_id using regex (robustly handles underscores)
            match = re.match(r"^(.+?)_\d+_\1$", prot_id)
            if not match:
                print(f"⚠️ Skipping bad prot_id: {prot_id}")
                continue
            genome_id = match.group(1)

            genome_to_prots[genome_id].append((cluster_id, prot_id))


'''
with clusters_path.open() as f:
    for line in f:
        if line.startswith(("S", "H")):
            fields = line.strip().split("\t")
            cluster_id = fields[1]  # Cluster ID from column 2
            if line.startswith("S"):
                seed_id = fields[8]
                cluster_proteins[cluster_id] = {
                    "seed": seed_id,
                    "hits": []
                }
            else:  # line starts with "H"
                hit_id = fields[8]
                cluster_proteins[cluster_id]["hits"].append(hit_id)
#print(cluster_proteins)
'''
#Funcio per parsejar el gff file de cada un dels genomes
def parse_gff(gff_path):
    cds_by_locus = {} #Diccionari per mapejar cada locus_tag, amb la respectiva start, end i strand
    cds_positions = [] #List per posicionar els CDS que estan Upstream de la proteina del cluster en la que estem centrats ara
    with gff_path.open() as f:
        for line in f:
            if line.startswith("#") or "\tCDS\t" not in line: #Les lines que comencin amb almoadilla no les volem
                continue # els continues no fan el que esperem, si no que tornen al for
            fields = line.strip().split("\t") #Splitegem la linea per poder accedir als camps corresponents amb les dades que ens interessen
            start = int(fields[3])
            end = int(fields[4])
            strand = fields[6]
            attributes = fields[8]
            attrs = dict(item.split("=") for item in attributes.split(";") if "=" in item) # Splitegem la linea de atributes per poder capturar el locus tag
            locus_tag = attrs.get("locus_tag") # Capturem el locus tag
            if locus_tag:
                cds_by_locus[locus_tag] = (start, end, strand) #We create a dictionary that we will acces through the locus tag.
                cds_positions.append((start, end, strand)) #We need this list of tupples, because when we look for the upstreamd CDS, we have to make sure that we look on the same strand, as CDS may be on differents strands.
    cds_positions.sort()
    return cds_by_locus, cds_positions

# Function to load the full genome sequence from .fna
def load_genome_sequence(fna_path):
    return next(SeqIO.parse(fna_path, "fasta")).seq #In this case, we need to use the next command, otherwise we are just activating the SeqIO processor, which would retunr nothing readable.

# Function to get upstream sequence from genome
def get_upstream(start, strand, genome_seq, cds_positions): #We get this values when we call the function. Start and strand will come from the cds_by_locus, the genome_seq from the previous function that loads the genome, and cds position, again from the previous function where the list is deffined
    if strand == "+":
        upstream_end = start - 1
        upstream_start = max(0, upstream_end - 250) # This is to avoid to go backwards to negative when lookin for the upstream regions
        for s, e, _ in reversed(cds_positions): # This is the trikiest part of the whole code: We are ignoring the strand here
            if e < upstream_end: #S and E stands for start and end of the list cds_positions. We are checking if that position of the cds, ends before the current CDS starts.
                if e >= upstream_start: #And here we check if that same position starts after the upstream area. In case that this condition does not happens, means that there is an intergeni region of 250 without a CDS. NO OVERLAPPING
                    upstream_start = e + 1 #If there is a CDS, truncate it, add one base to the end of the UPstream CDS, and save that positon as the upstream region. NO OVERLAPPING. Upstream region so it begins after the previous CDS ends
                break # When we find a CDS that satisfies the condition we break. OUT OF the cds_positions
        return genome_seq[upstream_start:upstream_end].upper() #In case the for goes through all the cds_position list, the for returns the 250 bases upstream, because there is no CDS overlapping the UPSTREAM REGION.
    elif strand == "-":
        upstream_start = start #Here we don't put a +1, because of python indexing, due to the 0 in the first position, the start is already the correct index to begin with
        upstream_end = min(len(genome_seq), upstream_start + 250)
        for s, e, _ in cds_positions:
            if s > upstream_start:
                if s <= upstream_end:
                    upstream_end = s - 1 #Here, we replicate the same logic, but instead of resting, we ad bases, because we look the upstream region, adding bases to the CDS location. On the minus strand, upstream of a CDS is numerically after the CDS start coordinates.
                break
        return genome_seq[upstream_start:upstream_end].reverse_complement().upper()



cluster_dir = Path("/home/carlesroses/PHROGs/upstream_regions") #We define the path to create a direcoty in PHROGs directory
cluster_dir.mkdir(exist_ok=True)

for folder in Genomes.iterdir():
    if folder.is_dir():
        genome_id = folder.name
        gff_path = folder / f"{genome_id}.gff"
        fna_path = folder / f"{genome_id}.fna"
        
        if gff_path.exists() and fna_path.exists():
        
        
        # Parse GFF and genome
            cds_by_locus, cds_positions = parse_gff(gff_path) #In here we are assigning the returned values of the function, with the two variables. We are using the same name but they could be compleatly different.
            #print(f"📘 Example locus_tags from {genome_id}:")
            #for tag in list(cds_by_locus.keys())[:10]:
                #print(" -", tag)

            genome_seq = load_genome_sequence(fna_path) #Loading the genome
            short_genome_id = genome_id.split(".")[0]
            
            r'''
            # For each cluster, process seed and hit proteins from this genome
            for cluster_id, members in cluster_proteins.items(): #Assigning two variables to the dictionary definned at the beggining called cluster_proteins
                for prot_id in [members["seed"]] + members["hits"]:
                    total_proteins += 1
                    match = re.match(r"^(.+?)_(\d+)_\1$", prot_id) #Hem hagut de ferho aixi perque si no teniem un problema per captuyrar els items.
                    

                    if not match:
                        print(f"⚠️ Skipping: unexpected prot_id format: {prot_id}")
                        continue
                    prot_genome = match.group(1)
                    locus_tag = f"{prot_genome}_{match.group(2)}" #Here we are capturing the groups created on line 98: group 1 is de genome_id, group 2, the numbers, gropu 3 the repeated genome_id

                    #print(f"{cluster_id}: {prot_id} → {prot_genome}, {locus_tag}")

                    if prot_genome != short_genome_id:
                        continue

                    if locus_tag not in cds_by_locus:
                        continue
                    
                    #if genome_id not in prot_id: #Genome ID is defined through al the folders in GENOME_DB. We just to focus on the genomes of the defined proteins, not all the genomes. So we are checking all the genomes
                     #   continue  # skip if protein not from current genome. If it is, lets get the locus_tag
                    #locus_tag = prot_id.split("_")[1]  # extract the locus tag from prot_id. We separete this MF612072_00008, to get only the 00008, which is the CDS lcation
                    #if locus_tag not in cds_by_locus: #Now we check if the locus tag is in the genome we just found the proteins are from, in the previous comand
                     #   continue  # no CDS info found, try another locus tag
                    
                    
                    #print(f"🧬 Locus tag found: {locus_tag} → {cds_by_locus[locus_tag]}")
                    start, end, strand = cds_by_locus[locus_tag] #Gettin the values of the specific locus tag, and definning the variables
                    upstream_seq = get_upstream(start if strand == "+" else end, strand, genome_seq, cds_positions) #Calling the function to get the upstream region

                    # Create output directory and file per cluster
                    
                    output_path = cluster_dir / f"cluster_{cluster_id}.fna"

                    with output_path.open("a") as out_fasta:
                        out_fasta.write(f">{prot_id}\n{upstream_seq}\n")
                        '''
            if genome_id not in genome_to_prots:
                continue  # Skip genomes with no relevant proteins

            for cluster_id, prot_id in genome_to_prots[genome_id]:
                total_proteins += 1
                match = re.match(r"^(.+?)_(\d+)_\1$", prot_id)
                if not match:
                    continue
                prot_genome = match.group(1)
                locus_tag = f"{prot_genome}_{match.group(2)}"

                if locus_tag not in cds_by_locus:
                    continue

                # Extract and write upstream sequence
                start, end, strand = cds_by_locus[locus_tag]
                upstream_seq = get_upstream(start if strand == "+" else end, strand, genome_seq, cds_positions)

                output_path = cluster_dir / f"cluster_{cluster_id}.fna"
                with output_path.open("a") as out_fasta:
                    out_fasta.write(f">{prot_id}\n{upstream_seq}\n")

                    written_sequences += 1
                print(f"🧬 Written: {written_sequences} | Processed: {total_proteins}", end="\r", flush=True)

print()  # newline to clear last progress line
print(f"🎉 Done! Total processed: {total_proteins}")
print(f"✅ Total upstream sequences written: {written_sequences}")
