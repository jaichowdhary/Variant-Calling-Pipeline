# Variant Calling Pipeline

This repository contains a reproducible, command-line variant calling pipeline implemented in Bash.  
The pipeline aligns short paired-end reads to a reference genome and produces variant calls in VCF format.

Tools used: 
1.BWA 
2.SAMtools
3.BCFtools


# Pipeline Overview

1. Index reference genome (`bwa index`, `samtools faidx`)
2. Align reads using `bwa mem`
3. Convert SAM to BAM (`samtools view`)
4. Sort and index BAM (`samtools sort`, `samtools index`)
5. Generate genotype likelihood estimates (`bcftools mpileup`)
6. Call variants (`bcftools call`)

# Requirements

The following tools must be installed and available in PATH:
- `bwa`
- `samtools`
- `bcftools`

# Example installation via conda: 

```bash
conda create -n vc_pipeline bwa samtools bcftools -c bioconda -c conda-forge
conda activate vc_pipeline
```

# Usage
bash scripts/run_pipeline.sh \
  -r reference.fna \
  -1 sample_R1.fastq \
  -2 sample_R2.fastq \
  -s sample_name \
  [-t threads]

# Notes
Pipeline was created as a learning tool to successfully align practice datasets and produce real world variant callign outputs.


