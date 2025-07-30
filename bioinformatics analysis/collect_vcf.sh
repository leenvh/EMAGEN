#!/bin/bash

# Destination directory for the collected VCF files
DEST_DIR="../sample_vcf/"

# Create the destination directory if it doesn't exist
mkdir -p "$DEST_DIR"

# Find all vc files and process them
find ../per_sample_vcf -type f -name "variants.g.vcf.gz" | while read -r filepath; do
  # Get the directory name of the current BAM file
  dir_name=$(dirname "$filepath")
  base_name=$(basename "$dir_name")

  # Copy and rename the vcf file to the destination directory
  cp "$filepath" "$DEST_DIR/${base_name}.g.vcf.gz"
done

