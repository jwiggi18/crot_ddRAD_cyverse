#!/bin/bash

#variant calling workflow to produce alignment with anolis caroleninsis genome for female C. collaris

#catch errors
set -e 

echo "before running this script index the reference genome"

echo " before running this script cd ~/Library/Mobile\ Documents/com~apple~CloudDocs/2.research/crots/ddRAD/"

#reference genome
genome=vcf/data/ref_genome/acar.fasta

#folder name within each results folder type (sam, bam, bcf, vcf). file structure = ddRAD/vcf/results/sam/all_only_female
location=all_only_female

#use a loop to run the variant calling workflow on each FASTQ file
#assign the name of the FASTQ file we are working with to a variable called fq1 and tell the script to echo the file name
for fq1 in vcf/data/processed3_f_m/cc_f/*.fq
    do
    echo "working with file $fq1"
#extract the base name of the file (excluding the path and .fastq extension) and assign it to a new variable called base
    base=$(basename $fq1 .fq)
    echo "base name is $base"

    fq1=vcf/data/processed3_f_m/cc_f/${base}.fq

    sam=vcf/results/sam/${location}/${base}.aligned.sam

    bam=vcf/results/bam/${location}/${base}.aligned.bam

    sorted_bam=vcf/results/bam/${location}/${base}.aligned.sorted.bam

    raw_bcf=vcf/results/bcf/${location}/${base}_raw.bcf

    variants=vcf/results/bcf/${location}/${base}_variants.vcf

    final_variants=vcf/results/vcf/${location}/${base}_final_variants.vcf

#align the reads to the reference genome and output a sam file
    bwa mem $genome $fq1 > $sam

#convert the SAM file to BAM format
    samtools view -S -b $sam > $bam

#sort the BAM file
    samtools sort -o $sorted_bam $bam

#index the BAM file for display purpose
    samtools index $sorted_bam

#calculate the read coverage of positions in the genome
    bcftools mpileup -O b -o $raw_bcf -f $genome $sorted_bam

#call SNVs with bcftools
    bcftools call --ploidy 2 -m -o $variants $raw_bcf

#filter and report the SNVs in variant calling format (VCF)
    vcfutils.pl varFilter $variants > $final_variants

    done