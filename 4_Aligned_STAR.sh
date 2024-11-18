#!/bin/bash

# Define the output directory
path=$(pwd) 
INPUT_DIR="${path}/trimmed_filtered"  # Ensure there's no trailing slash
OUTPUT_DIR="${path}/Aligned_OUT"  # Replace with your actual output directory
GENOME_DIR="${path}/Genome_Data"  # Ensure this points to your genome directory

# Make sure the output directory exists
mkdir -p "$OUTPUT_DIR"  # Create the output directory if it doesn't exist

# Read sample names from fastQFileList file
filesList=$(cat fastQlist | sort | uniq)
#echo "Sample list: $filesList"

# Loop through each sample name in the list
echo "$filesList" | while read file; do
    # Construct input file names
    inputFileR1="${INPUT_DIR}/${file}_R1.trimmed.fastq.gz"
    inputFileR2="${INPUT_DIR}/${file}_R2.trimmed.fastq.gz"

#--runThreadN 40

    # Run the STAR alignment command
    STAR --runThreadN 20 \
         --genomeDir "$GENOME_DIR" \
         --readFilesIn "$inputFileR1" "$inputFileR2" \
         --outFileNamePrefix "${OUTPUT_DIR}/${file}_" \
         --outSAMunmapped Within \
         --outSAMattributes Standard
    
done

echo "All samples processed"

