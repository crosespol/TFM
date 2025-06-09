#!/bin/bash

# Input/output paths
input_faa=~/PHROGs/phrog_fastas/phrog_4.faa
sorted_faa=~/PHROGs/phrog_fastas_sorted/sorted_phrog_4.faa
centroids_out=~/PHROGs/usearch_results_byPHROG/centroids_0.7/phrog_4_centroids.faa
clusters_out=~/PHROGs/usearch_results_byPHROG/clusters_0.7/phrog_4_clusters.uc

# Create output folders
mkdir -p ~/PHROGs/phrog_fastas_sorted
mkdir -p ~/PHROGs/usearch_results_byPHROG/centroids_0.7
mkdir -p ~/PHROGs/usearch_results_byPHROG/clusters_0.7

# Sort by length
usearch -sortbylength "$input_faa" -fastaout "$sorted_faa"

# Run clustering at 70% identity
usearch --cluster_fast "$sorted_faa" \
    -id 0.7 \
    -centroids "$centroids_out" \
    -uc "$clusters_out"

echo "✅ USEARCH clustering complete for phrog_4 at 70% identity"
