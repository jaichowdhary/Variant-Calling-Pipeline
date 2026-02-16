#!/usr/bin/env bash
set -euo pipefail

#ECOLI VARIANT CALLING PIPELINE DRAFT 2

 
#Updates: Added CLI Arguments and error handling
# using getopts for reproducibility

#defined variable for sample name, and all file outpus
SAMPLE=""

#defined reference genome
REF=""

#defined variable for number of threads used
THREADS=4
#defined raw reads
READ1=""
READ2=""

#adding CLI arguments 

while getopts ":r:1:2:s:t:" opt; do
    case "$opt" in
        r) REF="$OPTARG" ;;
        1) READ1="$OPTARG" ;;
        2) READ2="$OPTARG" ;;
        s) SAMPLE="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        \?) echo "-$OPTARG is an invalid flag. Correct flags : -r -1 -2 -s -t"; exit 1;;
        :) echo "Option -$OPTARG requires an argument."; exit 1;;
    esac
done

#error handling for empty ref_genome and raw_read files

if [[ -z "$REF" || -z "$READ1" || -z "$READ2" || -z "$SAMPLE" ]]; then
  echo "Correct Usage: $0 -r reference.fna -1 Read1.fastq -2 Read2.fastq -s sample_name -t thread_count"
  exit 1
fi

#defined aligned read output files
SAM="aligned_reads/${SAMPLE}.sam"
BAM="aligned_reads/${SAMPLE}.bam"
SORTED_BAM="aligned_reads/${SAMPLE}.sorted.bam"

#indexing of reference genome
bwa index "$REF" 
samtools faidx "$REF"

#alignment of raw reads to reference genome
bwa mem -t "$THREADS" "$REF" "$READ1" "$READ2" > "$SAM"

#formatting of aligned reads
samtools view -bS "$SAM" > "$BAM" 
samtools sort -@ "$THREADS" -o "$SORTED_BAM" "$BAM"
samtools index "$SORTED_BAM"

#Use bcftools to generate pileup output 
PILEUP="variants/${SAMPLE}.pileup.bcf"
bcftools mpileup -f "$REF" -Ob -o "$PILEUP" "$SORTED_BAM" 

#Use bcftools to call actual variants
VARIANTS="variants/${SAMPLE}.vcf"
bcftools call -mv -Ov -o "$VARIANTS" "$PILEUP"

