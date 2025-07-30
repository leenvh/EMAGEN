#!/bin/bash

# Path to the folder containing FASTQ files
fastq_folder="/mnt/storage13/ahri/MARS/Raw_fastq"
# Path to the reference genome
ref="/mnt/storage13/ahri/plasmodium_falciparum/Pfalciparum.genome.fasta"
# Path to the main pipeline script
script="/mnt/storage13/ahri/MARS/vcf_updated/scripts/run_variant_calling.sh"

# Loop through all R1 files in the folder
for r1 in "$fastq_folder"/*_R1_001.fastq.gz; do
    # Determine the corresponding R2 file
    r2="${r1/_R1_001.fastq.gz/_R2_001.fastq.gz}"
    # Extract the sample name (remove suffixes after '_S' for clarity)
    sample_name=$(basename "$r1" | sed 's/_R1_001.fastq.gz//')

    # Debug: Print the current files being processed
    echo "Processing sample: $sample_name"
    echo "R1 file: $r1"
    echo "R2 file: $r2"

    # Check if the corresponding R2 file exists
    if [[ ! -f "$r2" ]]; then
        echo "R2 file missing for $sample_name. Skipping."
        continue
    fi

    # Run the main pipeline script
    bash "$script" "$ref" "$r1" "$r2" "$sample_name"
done
