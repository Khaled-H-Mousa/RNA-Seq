#!/bin/bash

# This script is used to count the number of reads mapped to each gene using featureCounts

# Define paths
path=$(pwd) 
ExperimentName="/${path}/raw"
genomeName="${path}/Genome_Data/reference_genome.fasta"
AnnotationFolder="${path}/Genome_Data/annotation.gtf"

inputFolder="${path}/Processed_BAM"  # Path to the directory containing sorted BAM files
outFolder="${path}/Count" # Output directory for featureCounts results

# Create the output folder if it doesn't exist
mkdir -p "$outFolder"

# Run featureCounts for each sorted BAM file in the input folder
for sample in "$inputFolder"/*_Aligned.sorted.bam  # Assumes BAM files are named with *_sorted.bam
do
    # Extract the sample name (removing the "_sorted.bam" suffix)
    sampleName=$(basename "$sample" "_Aligned.sorted.bam")

    # Run featureCounts
    featureCounts -p -O -T 8 -a "$AnnotationFolder" \
                  -o "$outFolder"/"$sampleName"_featureCounts.txt "$sample"

    echo "Processed $sampleName"
done

echo "All samples processed successfully."
