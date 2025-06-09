from collections import Counter

uc_file = "/home/carlesroses/PHROGs/usearch_results_byPHROG/clusters_0.7/phrog_4_clusters.uc"

cluster_counts = Counter()

with open(uc_file) as f:
    for line in f:
        if line.startswith("S"):
            fields = line.strip().split("\t")
            rep_seq = fields[8]  # Seed's own ID
            cluster_counts[rep_seq] += 1  # Start cluster with seed
        elif line.startswith("H"):
            fields = line.strip().split("\t")
            rep_seq = fields[9]  # ID of the representative sequence
            cluster_counts[rep_seq] += 1  # Add hit to cluster

# Sort clusters by size in descending order
sorted_clusters = sorted(cluster_counts.items(), key=lambda x: x[1], reverse=True)

# Print result
print(f"{'Cluster ID':<35} {'Size'}")
print("-" * 45)
for cluster_id, size in sorted_clusters:
    print(f"{cluster_id:<35} {size}")
