#!/bin/bash

# Create output base folder
mkdir -p ~/PHROGs/usearch_results_byPHROG
mkdir -p ~/PHROGs/phrog_fastas/sorted
cd ~/PHROGs/phrog_fastas || exit 1

for file in *.faa; do
	usearch -sortbylength "$file" -fastaout ~/PHROGs/phrog_fastas/sorted/"sorted_$file"
done
cd ~/PHROGs/phrog_fastas/sorted || exit 1
# Loop over thresholds
for threshold in 0.6 0.7 0.8 0.9; do
	mkdir -p ~/PHROGs/usearch_results_byPHROG/centroids_$threshold 
	mkdir -p ~/PHROGs/usearch_results_byPHROG/clusters_$threshold

	for file in sorted_*.faa; do
		base=$(basename "$file" .faa)
		base=${base#sorted_}  # remove leading 'sorted_'

		centroids_out=~/PHROGs/usearch_results_byPHROG/centroids_$threshold/${base}_centroids.faa
		clusters_out=~/PHROGs/usearch_results_byPHROG/clusters_$threshold/${base}_clusters.uc

		usearch --cluster_fast "$file" -id $threshold -centroids "$centroids_out" -uc "$clusters_out"

		echo "Clustering done for PHROG: $base"
		
	done

	echo "✅ Clustering completed for threshold: $threshold"
done

echo "All clusterings complete"

# Count sequences in centroid files
for threshold in 0.6 0.7 0.8 0.9; do
	echo "Centroid counts for threshold $threshold:"
	for file in ~/PHROGs/usearch_results_byPHROG/centroids_$threshold/*.faa; do
		count=$(grep -c "^>" "$file")
		echo "$(basename "$file"): $count centroids"
	done
done
