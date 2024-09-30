# Single-cell RNA-seq Analysis of Equine Endometrium Before and After Implantation

## Overview
This repository contains the workflow and scripts for analyzing single-cell RNA sequencing (scRNA-seq) data from equine endometrium at two key stages of gestation: Days 33 and 42. The analysis involves:
1. Downloading and preparing the reference genome and gene annotation files.
2. Filtering the gene annotations to include relevant biotypes.
3. Creating a custom reference genome using `Cell Ranger` for read alignment.
4. Processing raw sequencing files (FASTQ) to generate count matrices for downstream analysis.

## Workflow Summary
The process begins by downloading genome and annotation data from Ensembl, followed by filtering the annotations to retain specific gene biotypes. A custom reference genome is then built using `Cell Ranger`, which is subsequently used to process the raw sequencing data and create count matrices.

### Key Steps:
1. **Download Reference Genome and Annotation Files**  
   Obtain the reference genome (FASTA) and gene annotations (GTF) from Ensembl for *Equus caballus*. These files serve as the foundation for aligning the sequencing reads.

2. **Filter GTF Annotations**  
   Filter the GTF file to include relevant biotypes like protein-coding genes, long non-coding RNAs (lncRNAs), and immune-related genes. This ensures that only the relevant transcripts are used in the analysis.

3. **Create a Custom Reference Genome**  
   Use `Cell Ranger` to create a reference genome that can be used to align the single-cell RNA sequencing data. This includes both the FASTA file and the filtered GTF file.

4. **Process the Raw Data**  
   Align the raw FASTQ files to the custom reference genome and generate a count matrix for downstream analysis. Each sample is processed separately using `Cell Ranger`.

## Setup Instructions

### Requirements
To run this workflow, you will need the following:
- `Cell Ranger` (v7.1.0 or higher)
- Basic utilities like `wget` and `gunzip` for file download and decompression
- Access to *Equus caballus* genome and GTF annotation files from Ensembl (release 109)

### Step-by-Step Instructions

#### 1. Download the Reference Genome and GTF Files
Download the genome sequence (FASTA) and gene annotation (GTF) files from Ensembl.

```bash
# Download the reference genome (FASTA)
wget https://ftp.ensembl.org/pub/release-109/fasta/equus_caballus/dna/Equus_caballus.EquCab3.0.dna.toplevel.fa.gz
gunzip Equus_caballus.EquCab3.0.dna.toplevel.fa.gz

# Download the gene annotations (GTF)
wget https://ftp.ensembl.org/pub/release-109/gtf/equus_caballus/Equus_caballus.EquCab3.0.109.gtf.gz
gunzip Equus_caballus.EquCab3.0.109.gtf.gz
```

#### 2. Filter the GTF File
Filter the gene annotations to keep only the biotypes that are necessary for the analysis, such as protein-coding genes and long non-coding RNAs.

```bash
cellranger mkgtf   Equus_caballus.EquCab3.0.109.gtf Equus_caballus.EquCab3.0.109.filtered.gtf   --attribute=gene_biotype:protein_coding   --attribute=gene_biotype:lncRNA   --attribute=gene_biotype:antisense   --attribute=gene_biotype:IG_V_gene   --attribute=gene_biotype:TR_V_gene   --attribute=gene_biotype:TR_V_pseudogene   --attribute=gene_biotype:IG_J_gene   --attribute=gene_biotype:IG_J_pseudogene   --attribute=gene_biotype:IG_C_gene   --attribute=gene_biotype:IG_C_pseudogene
```

#### 3. Create the Custom Reference Genome
Use the filtered GTF and the downloaded FASTA file to build a reference genome that will be used in the analysis.

```bash
cellranger mkref   --genome=EquCab3.0   --fasta=Equus_caballus.EquCab3.0.dna.toplevel.fa   --genes=Equus_caballus.EquCab3.0.109.filtered.gtf   --ref-version=1.0.0
```

#### 4. Process the FASTQ Files
Process the raw sequencing data (FASTQ) using the custom reference genome to generate count matrices. For each sample, use the following command:

```bash
cellranger count --id=Sample_XYZ --fastqs=/path/to/FASTQ/ --sample=Sample_XYZ --transcriptome=/path/to/EquCab/EquCab3
```

# Additional Information
- **Custom Reference Genome**: The reference genome created here will be saved in the `EquCab3.0` directory and will be used for all subsequent alignment and counting steps.
- **Count Matrices**: The output from the `cellranger count` command will include a matrix file that contains the number of reads for each gene in each cell, which can be used in downstream analysis, including clustering and differential expression.
- **Data Anonymization**: All data used in this analysis have been anonymized to remove personal and sensitive information, ensuring compliance with data privacy regulations.

## Contact Information
For any questions or issues related to this workflow, please contact:  
**Dawid Tobolski** at [tobola28@gmail.com]

