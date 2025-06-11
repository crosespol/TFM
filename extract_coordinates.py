from pathlib import Path

GFF_BASE_DIR = Path("/home/carlesroses/PHROGs/GenomesDB")
LOCUS_TAG_FILE = Path("/home/carlesroses/PHROGs/locus_tag_phrog4.txt")

def extract_end_position(gff_path, locus_tag):
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) != 9:
                continue
            if fields[2] != "CDS":
                continue
            if f"locus_tag={locus_tag}" in fields[8]:
                return int(fields[4])
    return None

with open(LOCUS_TAG_FILE) as f:
    for line in f:
        full = line.strip()
        if not full:
            continue

        # Exemple: OR896311_00008_OR896311
        parts = full.split("_")
        genome_id = parts[0]  # OR896311
        locus_tag = f"{parts[0]}_{parts[1]}"  # OR896311_00008

        gff_path = GFF_BASE_DIR / genome_id / f"{genome_id}.gff"

        if not gff_path.exists():
            print(f"[WARNING] GFF not found for {genome_id}")
            continue

        end = extract_end_position(gff_path, locus_tag)
        if end:
            print(f"{genome_id}\t{end}")
        else:
            print(f"[ERROR] locus_tag {genome_id} not found in {gff_path}")
