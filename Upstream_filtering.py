from Bio import SeqIO, Align
from pathlib import Path
from calculate_identity import calculate_identity

cluster_dir = Path("/home/carlesroses/PHROGs/upstream_regions_filtered") #We define the path to create a direcoty in PHROGs directory
cluster_dir.mkdir(exist_ok=True)

# Directory with cluster files
upstream_dir = Path("upstream_regions")
identity_threshold = 0.9

# Prepare aligner
aligner = Align.PairwiseAligner()
aligner.mode = "global"

# Prepare summary
summary = []

for fasta_file in sorted(upstream_dir.glob("cluster_*.fna")):
    sequences = list(SeqIO.parse(fasta_file, "fasta"))
    to_remove = set()

    for i, seq1 in enumerate(sequences):
        if i in to_remove or len(seq1.seq) == 0:
            continue
        for j in range(i + 1, len(sequences)):
            if j in to_remove:
                continue
            seq2 = sequences[j]
            if str(seq1.seq) == str(seq2.seq) or len(seq2.seq) == 0:
                to_remove.add(j)
                continue
            alignment = aligner.align(str(seq1.seq), str(seq2.seq))
            identity = calculate_identity(alignment)
            if identity >= identity_threshold:
                to_remove.add(j)

    filtered = [seq for i, seq in enumerate(sequences) if i not in to_remove]
    output_path = cluster_dir / f"{fasta_file.stem}_{identity_threshold}_filtered.fna"
    SeqIO.write(filtered, output_path, "fasta")

    # Print summary per file
    print(f"📦 {fasta_file.name}: {len(sequences)} → {len(filtered)}")

    summary.append((fasta_file.name, len(sequences), len(filtered)))

# Optional: save summary
summary_path = Path(f"/home/carlesroses/PHROGs/filtered_th{identity_threshold}_summary.csv")
with open(summary_path, "w") as out:
    out.write("file,original,filtered\n")
    for name, original, filtered in sorted(summary, key=lambda x: x[1], reverse=True):
        out.write(f"{name},{original},{filtered}\n")

print(f"\n✅ Done! Summary saved to: {summary_path}")
