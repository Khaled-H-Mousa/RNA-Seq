#!/bin/bash


# Define paths to genome output directory, genome FASTA file, and annotation GTF file
path=$(pwd) 
genomeDir="${path}/Genome_Data"        # Output directory for the genome index
genomeFastaFile="${path}/Genome_Data/Sesamum_indicum.fasta"     # Input genome FASTA file
sjdbGTFfile="${path}/Genome_Data/Sesamum_indicum.gtf"          # Annotation GTF file


# Run STAR genomeGenerate command
STAR --runThreadN 8 \
     --runMode genomeGenerate \
     --genomeDir "$genomeDir" \
     --genomeFastaFiles "$genomeFastaFile" \
     --sjdbGTFfile "$sjdbGTFfile" \
     --genomeSAindexNbases 13 \
     --sjdbOverhang 99

echo "STAR genome index generation complete."

