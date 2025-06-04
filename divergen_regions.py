import csv
import os
import re
reference_PHROGs = set()
Genomes = "/home/carlesroses/PHROGs/GenomesDB"


def read_PHROGs(direction, target_set, col_index):
    with open(direction, "r") as file:
        header = file.readline()
        for line in file:
            columns = line.strip().split(",")
            if len(columns) > col_index:
                target_set.add(columns[col_index].strip())


def extract_features(path):
    cds_list = []
    with open(path, "r") as genome:
        for line in genome:
            if line.startswith("#"):
                continue    
            columns = line.strip().split("\t")
            if len(columns) > 8 and columns[2] == "CDS":
                genome = columns[0] 
                start = int(columns[3])
                end = int(columns[4])
                strand = columns[6]
                atributes = columns[8] 
                phrog_id = None
                match = re.search(r"phrog_\d+", atributes)
                if match:
                    phrog_id = match.group()

                cds_list.append({
                    'genome': genome,
                    'start': start,
                    'end': end,
                    'strand': strand,
                    'phrog': phrog_id
                })                       
    return cds_list

def find_divergent(features, reference_PHROGs):
    global intergenic_regions
    global counter_new_PHROGs
    global reference_counts
    global new_phrog_counts
    for i in range(len(features) - 1):
        gene1 = features[i]
        gene2 = features[i + 1]     
        intergenic_size = gene2['start'] - gene1['end']
        if intergenic_size > 15 and (gene1['strand'] == "-" and gene2['strand'] == "+"):
            intergenic_regions += 1
            if gene1['phrog'] and gene1['phrog'] in reference_PHROGs and gene2['phrog'] and gene2['phrog'] not in reference_PHROGs:
                new_PHROGs.add((gene2['genome'], gene2['phrog']))
                counter_new_PHROGs += 1
                reference_counts[gene1['phrog']] = reference_counts.get(gene1['phrog'], 0) + 1
                new_phrog_counts[gene2['phrog']] = new_phrog_counts.get(gene2['phrog'], 0) + 1

            #Mirror
            elif gene2['phrog'] and gene2['phrog']in reference_PHROGs and gene1['phrog'] and gene1['phrog'] not in reference_PHROGs:
                new_PHROGs.add((gene1['genome'], gene1['phrog']))
                counter_new_PHROGs += 1
                reference_counts[gene2['phrog']] = reference_counts.get(gene2['phrog'], 0) + 1
                new_phrog_counts[gene1['phrog']] = new_phrog_counts.get(gene1['phrog'], 0) + 1
            
            if i + 2 < len(features):
                gene3 = features[i + 2]
                if gene3['phrog'] and gene3['strand'] == "+" and gene3['phrog'] in reference_PHROGs:
                    if gene1['phrog'] and gene1['phrog'] not in reference_PHROGs:
                        new_PHROGs.add((gene1['genome'], gene1['phrog']))
                        counter_new_PHROGs += 1
                        reference_counts[gene3['phrog']] = reference_counts.get(gene3['phrog'], 0) + 1
                        new_phrog_counts[gene1['phrog']] = new_phrog_counts.get(gene1['phrog'], 0) + 1

                    if gene2['phrog'] and gene2['phrog'] not in reference_PHROGs:
                        new_PHROGs.add((gene2['genome'], gene2['phrog']))
                        counter_new_PHROGs += 1
                        new_phrog_counts[gene2['phrog']] = new_phrog_counts.get(gene2['phrog'], 0) + 1
            #Mirror
            if i - 1 >= 0:
                gene0 = features[i - 1]
                if gene0['phrog'] and gene0['strand'] == "-" and gene0['phrog'] in reference_PHROGs:
                    if gene1['phrog'] and gene1['phrog'] not in reference_PHROGs:
                        new_PHROGs.add((gene1['genome'], gene1['phrog']))
                        counter_new_PHROGs += 1
                        reference_counts[gene0['phrog']] = reference_counts.get(gene0['phrog'], 0) + 1
                        new_phrog_counts[gene1['phrog']] = new_phrog_counts.get(gene1['phrog'], 0) + 1
                    if gene2['phrog'] and gene2['phrog'] not in reference_PHROGs:
                        new_PHROGs.add((gene2['genome'], gene2['phrog']))
                        counter_new_PHROGs += 1
                        new_phrog_counts[gene2['phrog']] = new_phrog_counts.get(gene2['phrog'], 0) + 1
i = 0
reference_PHROGs = set()
undiscovered = True
while undiscovered:    

    if i == 0:
        read_PHROGs("/home/carlesroses/PHROGs/366+handpick.csv", reference_PHROGs, 1) 

    processed_genomes = 0
    counter_new_PHROGs = 0
    reference_counts = {}
    new_phrog_counts = {}
    new_PHROGs = set()
    genome_summary = []
    total_intergenic_regions = 0
    for genome_folder in os.listdir(Genomes):
        genome_path = os.path.join(Genomes, genome_folder)
        if os.path.isdir(genome_path):  # skip the .fasta files
            for file in os.listdir(genome_path):
                if file.endswith(".gff"):
                    name_genome = file.replace(".gff", "")      
                    gff_path = os.path.join(genome_path, file)
                    CDS_list = extract_features(gff_path)
                    processed_genomes += 1
                    intergenic_regions = 0
                    find_divergent(CDS_list, reference_PHROGs)
                    total_intergenic_regions += intergenic_regions
                    genome_summary.append({
                        'genome': name_genome,
                        'intergenic_regions': intergenic_regions
                    })
                    print(f"Processed genomes {processed_genomes}", end="\r")
    

    print(f"Processed genomes {processed_genomes}")
    print(f"Total new PHROGs found: {len(new_PHROGs)}")
    print(f"New phrogs appear in {counter_new_PHROGs} intergenic regions")
    unique_phrogs = {phrog for _, phrog in new_PHROGs if phrog is not None}
    print(f"Total unique PHROG IDs: {len(unique_phrogs)}")
    

    genome_info = {g['genome']: g for g in genome_summary}

    if not new_PHROGs:
        print("⚠️ No new PHROGs detected — skipping file write.")
        
    else:
        phrog_records = []      
        for genome, phrog in new_PHROGs:
            if phrog:
                intergenic = genome_info.get(genome, {}).get('intergenic_regions', 0)
                times = new_phrog_counts.get(phrog, 0)
                phrog_records.append((genome, intergenic, phrog, times))

        phrog_records_sorted = sorted(phrog_records, key=lambda x: x[3],  reverse=True)
        already_written = set()
        with open(f"New_PHROGs_from_366+handpick_{i}.csv", "w") as output:
            output.write("genome,intergenic_regions,phrog_ID,times\n")
            for genome, intergenic, phrog, times in phrog_records_sorted:
                if phrog not in already_written:
                    output.write(f"{genome},{intergenic},{phrog},{times}\n")
                    already_written.add(phrog)

        print(f"📄 Results saved to 'New_PHROGs_from_366+handpick_{i}.csv'")

        # Save reference PHROG appearance counts
        
        # Merge counts from reference and new PHROGs
        reference_phrog_total_counts = {}

        # Add reference PHROGs that were used in this iteration
        for phrog, count in reference_counts.items():
            reference_phrog_total_counts[phrog] = reference_phrog_total_counts.get(phrog, 0) + count

        # Add new PHROGs that are now becoming references
        
        added_PHROGs = 0
        for phrog, count in new_phrog_counts.items():
            if count >= 10:
                reference_PHROGs.add(phrog)  # ensure it's also in the reference list
                added_PHROGs += 1
        if added_PHROGs >= 10:
            print(f"This round {added_PHROGs} PHROGs have been added to the reference set")
            # Save the updated reference set with counts
            sorted_phrogs = sorted(reference_PHROGs, key=lambda phrog: reference_phrog_total_counts.get(phrog, 0), reverse=True)
            with open(f"Reference_PHROGs_366+handpick_{i}.csv", "w") as ref_output:
                ref_output.write("phrog_ID,appearances\n")
                for phrog in sorted_phrogs:
                    count = reference_phrog_total_counts.get(phrog, 0)
                    ref_output.write(f"{phrog},{count}\n")

            print(f"📄 Reference PHROGs with counts saved to 'Reference_PHROGs_366+handpick_{i}.csv'")

            i += 1
        else:
            undiscovered = False  
            print(f"Total intergenic regions: {total_intergenic_regions}")
            print("⚠️ Less than 10 PHROGs detected — skipping file write.")