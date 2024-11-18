# RNA-Seq Analysis Pipeline

This repository provides a comprehensive RNA-Seq analysis workflow, from raw data quality control to gene-level expression quantification. The pipeline ensures reproducibility and efficiency, using industry-standard tools and structured steps to process high-throughput sequencing data.

---

## Overview of the Workflow

1. **Quality Control**: Assess raw reads for quality using FastQC.
2. **Trimming**: Remove adapters and low-quality bases with Trimmomatic.
3. **Genome Indexing**: Prepare the reference genome for alignment using HISAT2 or STAR.
4. **Read Alignment**: Map reads to the reference genome using HISAT2 or STAR.
5. **Format Conversion**: Convert and process alignment files using samtools.
6. **Feature Counting**: Quantify gene-level expression using featureCounts.

---

## Required Tools

### Prerequisites:
- **[FastQC](https://www.bioinformatics.babraham.ac.uk/projects/fastqc/):** Quality control for raw reads.
- **[Trimmomatic](http://www.usadellab.org/cms/?page=trimmomatic):** Read trimming and adapter removal.
- **[gffread](https://github.com/gpertea/gffread):** Conversion of GFF files to GTF format.
- **[HISAT2](https://daehwankimlab.github.io/hisat2/):** Splice-aware alignment.
- **[STAR](https://github.com/alexdobin/STAR):** High-performance RNA-Seq alignment.
- **[samtools](http://www.htslib.org/):** Manipulation of SAM/BAM files.
- **[featureCounts](http://subread.sourceforge.net):** Gene-level read quantification.

---

## Detailed Workflow

### **1. Quality Control**
Script: `1-Trimmomatic_QC.sh`

- **Purpose**: Assess raw sequencing data quality and prepare for trimming.
- **Steps**:
  - Use FastQC to evaluate read quality.
  - Generate comprehensive multi-sample QC reports using multiQC.

Command Example:
```bash
fastqc -o qc_reports/ sample_R1.fastq.gz sample_R2.fastq.gz
```

---

### **2. Trimming**
Script: `1-Trimmomatic_QC.sh`

- **Purpose**: Remove adapters and low-quality regions from raw reads.
- **Tools**: Trimmomatic.
- **Steps**:
  - Set up paths for input/output directories.
  - Perform paired-end trimming to create paired and unpaired output reads.

Command Example:
```bash
trimmomatic PE -threads 4 -phred33 \
    sample_R1.fastq.gz sample_R2.fastq.gz \
    trimmed_R1_paired.fastq.gz trimmed_R1_unpaired.fastq.gz \
    trimmed_R2_paired.fastq.gz trimmed_R2_unpaired.fastq.gz \
    ILLUMINACLIP:adapters.fa:2:30:10 SLIDINGWINDOW:4:20 MINLEN:25
```

---

### **3. Genome Indexing**
Scripts: `3-GenomeIndexHISAT2.sh` and `3-GenomeIndexSTAR.sh`

- **Purpose**: Prepare the reference genome for alignment.
- **Tools**: HISAT2 or STAR.
- **Steps**:
  - Extract splice sites and exons for HISAT2.
  - Use STAR to generate a genome index.

Command Example (STAR):
```bash
STAR --runThreadN 4 \
     --runMode genomeGenerate \
     --genomeDir genome_index/ \
     --genomeFastaFiles reference_genome.fasta \
     --sjdbGTFfile annotation.gtf \
     --sjdbOverhang 100
```

---

### **4. Read Alignment**
Scripts: `4-HISAT2_Align.sh` and `4-STAR_Align.sh`

- **Purpose**: Align trimmed reads to the reference genome.
- **Tools**: HISAT2 or STAR.
- **Steps**:
  - Map paired-end reads to the indexed genome.
  - Output sorted BAM files for downstream analysis.

Command Example (STAR):
```bash
STAR --runThreadN 4 \
     --genomeDir genome_index/ \
     --readFilesIn trimmed_R1_paired.fastq.gz trimmed_R2_paired.fastq.gz \
     --readFilesCommand zcat \
     --outFileNamePrefix aligned/sample_ \
     --outSAMtype BAM SortedByCoordinate
```

---

### **5. Format Conversion**
Script: `5-sam-to-bam.sh`

- **Purpose**: Convert alignment files (SAM to BAM), sort, and index them.
- **Tools**: samtools.
- **Steps**:
  - Convert SAM to BAM format.
  - Sort and index BAM files for efficient access.

Command Example:
```bash
samtools view -@ 4 -Sb input.sam | samtools sort -o sorted_output.bam
samtools index sorted_output.bam
```

---

### **6. Feature Counting**
Script: `6-featureCounts.sh`

- **Purpose**: Quantify gene-level expression using sorted BAM files.
- **Tools**: featureCounts.
- **Steps**:
  - Use GTF annotation to count aligned reads for each feature.
  - Generate a tab-delimited file with gene counts.

Command Example:
```bash
featureCounts -p -O -T 4 \
    -a annotation.gtf \
    -o counts/example_featureCounts_output.txt \
    aligned/sample_Aligned.sortedByCoord.out.bam
```

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
