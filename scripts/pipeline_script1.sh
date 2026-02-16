#!/usr/bin/env bash
set -euo pipefail

#ECOLI VARIANT CALLING PIPELINE DRAFT 1

#defined reference genome
REF="ref_genome/ecoli_ref.fna"

#defined raw reads
READ1="raw_reads/SRR2584857_1.fastq"
READ2="raw_reads/SRR2584857_2.fastq"

#defined aligned read output files
SAM="aligned_reads/ecoli_aligned_reads.sam"
BAM="aligned_reads/ecoli_aligned_reads.bam"
SORTED_BAM="aligned_reads/ecoli_aligned_reads.sorted.bam"

#indexing of reference genome
bwa index "$REF" 
samtools faidx "$REF"

#alignment of raw reads to reference genome
bwa mem -t 4 "$REF" "$READ1" "$READ2" > "$SAM"

#formatting of aligned reads
samtools view -bS "$SAM" > "$BAM" 
samtools sort -@ 4 -o "$SORTED_BAM" "$BAM"
samtools index "$SORTED_BAM"

#Use bcftools to generate pileup output 
PILEUP="variants/ecoli.pileup.bcf"
bcftools mpileup -f "$REF" -Ob -o "$PILEUP" "$SORTED_BAM" 

#Use bcftools to call actual variants
VARIANTS="variants/ecoli_variants.vcf"
bcftools call -mv -Ov -o "$VARIANTS" "$PILEUP"


