# RNA-Seq Analysis Pipeline

RNA-Seq analysis is a powerful approach for studying transcriptomics, enabling the identification and quantification of RNA molecules in a sample. The process includes quality control of raw reads, trimming low-quality regions, aligning reads to a reference genome, and quantifying gene expression. This workflow provides insights into gene activity, alternative splicing, and differential expression, making it essential for understanding functional genomics and transcriptome dynamics.

## Steps in the Pipeline

1. **Quality Control (FastQC)**
2. **Read Trimming (Trimmomatic)**
3. **Genome Indexing (STAR)**
4. **Read Alignment (STAR)**
5. **Format Conversion (samtools)**
6. **Feature Counting (featureCounts)**

---

## Prerequisites

The following tools are required to run the pipeline:
- **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/):** For quality control of raw sequencing reads.
- **[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic):** For trimming adapters and low-quality bases.
- **[STAR](https://github.com/alexdobin/STAR):** For genome indexing and alignment.
- **[samtools](http://www.htslib.org/):** For Manipulation of SAM/BAM files.
- **[featureCounts](http://subread.sourceforge.net):** For read quantification.

### Input Requirements
- **Raw Reads**: Paired-end FASTQ files.
- **Reference Genome**: FASTA file.
- **Annotation File**: GTF file.

---

## Steps in the Pipeline

### 1. Quality Control (FastQC)
Evaluate the quality of raw RNA-Seq reads.
```bash
fastqc -o qc_reports/ sample_R1.fastq.gz sample_R2.fastq.gz
```
- **Input**: Raw FASTQ files.
- **Output**: HTML quality reports in the `qc_reports/` directory.

---

### 2. Read Trimming (Trimmomatic)
Remove adapters and low-quality bases.
```bash
trimmomatic PE -threads 4 -phred33 \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    trimmed_R1_paired.fastq.gz trimmed_R1_unpaired.fastq.gz \
    trimmed_R2_paired.fastq.gz trimmed_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:25
```
- **Input**: Raw FASTQ files.
- **Output**: Trimmed paired and unpaired FASTQ files.

---

### 3. Genome Indexing (STAR)
Generate a genome index for alignment.
```bash
STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --genomeDir genome_index/ \
     --genomeFastaFiles reference_genome.fasta \
     --sjdbGTFfile annotation.gtf \
     --sjdbOverhang 100
```
- **Input**: Reference genome (FASTA) and annotation (GTF).
- **Output**: Genome index files in the `genome_index/` directory.

---

### 4. Read Alignment (STAR)
Align trimmed reads to the indexed reference genome.
```bash
STAR --runThreadN 4 \
     --genomeDir genome_index/ \
     --readFilesIn trimmed_R1_paired.fastq.gz trimmed_R2_paired.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix aligned/sample_ \
     --outSAMtype BAM SortedByCoordinate
```
- **Input**: Trimmed reads and genome index.
- **Output**: Sorted BAM file in the `aligned/` directory.

---

### **5. Format Conversion**
Convert alignment files (SAM to BAM), sort, and index them.
```bash
samtools view -@ 4 -Sb input.sam | samtools sort -o sorted_output.bam
samtools index sorted_output.bam
```
- **Input**: SAM file.
- **Output**: Sorted BAM file in the `Processed_BAM/` directory.

---

### 6. Feature Counting (featureCounts)
Count aligned reads for each gene.
```bash
featureCounts -p -O -T 4 \
    -a annotation.gtf \
    -o counts/example_featureCounts_output.txt \
    aligned/sample_Aligned.sortedByCoord.out.bam
```
- **Input**: Sorted BAM file and annotation (GTF).
- **Output**: Gene counts in `counts/` directory.

---
## Folder Structure
```plaintext
├── raw_data/               # Raw FASTQ files
├── qc_reports/             # Quality control reports
├── trimmed_reads/          # Trimmed paired and unpaired reads
├── genome_index/           # Indexed genome (HISAT2/STAR)
├── aligned/                # Sorted BAM files from alignment
├── counts/                 # Gene-level read counts
```

---

## Conclusion

This pipeline provides a structured approach to RNA-Seq data analysis, integrating quality control, trimming, alignment, and quantification. It is versatile and can be customized for different organisms and experimental designs, offering researchers a robust framework for transcriptome analysis.
