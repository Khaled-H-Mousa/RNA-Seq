#!/bin/bash

# Define paths
path=$(pwd)                 # Current directory (or specify another path)
inputFolder="${path}/raw"   # Folder with raw reads
filteredDir="${path}/fastqc_output"  # Folder for FastQC output

# Create output directory if it doesn't exist
mkdir -p "$filteredDir"

# Read sample names from fastQFileList file
filesList=$(cat fastQList | sort | uniq)

# Loop through each sample name in the list
echo "$filesList" | while read file; do
    echo "Processing sample: $file"

    # Input files
    forwardFile="${inputFolder}/${file}_R1.fastq.gz"
    reverseFile="${inputFolder}/${file}_R2.fastq.gz"

    # Run FastQC for paired-end reads
    fastqc -o "$filteredDir" -t 8 "$forwardFile" "$reverseFile"
done

echo "All samples processed."

