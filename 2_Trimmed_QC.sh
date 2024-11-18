#!/bin/bash

# Define paths
path=$(pwd)                 # Current directory (or specify another path)
inputFolder="${path}/raw/"   # Folder with raw reads
filteredDir="${path}/trimmomatic_filtered/"  # Folder for trimmed output

# Create output directory if it doesn't exist
mkdir -p "$filteredDir"

# Read sample names from fastQFileList file
filesList=$(cat fastQList | sort | uniq)
echo "Samples to process"

# Loop through each sample name in the list
echo "$filesList" | while read file
do
    echo "Processing sample"

    # Run Trimmomatic for paired-end trimming
    trimmomatic PE -threads 8 -phred33 \
                   "${inputFolder}${file}_R1.fastq.gz" "${inputFolder}${file}_R2.fastq.gz" \
                   "${filteredDir}${file}_R1.trimmed.fastq.gz" "${filteredDir}${file}_R1_un.trimmed.fastq.gz" \
                   "${filteredDir}${file}_R2.trimmed.fastq.gz" "${filteredDir}${file}_R2_un.trimmed.fastq.gz" \
                   SLIDINGWINDOW:4:20 MINLEN:25 ILLUMINACLIP:NexteraPE-PE.fa:2:40:15
done

echo "All samples processed"

