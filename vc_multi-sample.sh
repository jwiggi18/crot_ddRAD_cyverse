

#produce a multipsample variant call taht contains all the information on all samples

genome=vcf/data/ref_genome/acar.fasta

input_location=vcf/data/results/vcf

results_location=vcf/data/results/bam/all_only_female

# output type (-O) uncompressed VCF (v)
bcftools mpileup -O v -f ${genome} ${input_location}/*.bam > ${results_location}/f_genotypes.vcf

# norm: output only first reference of record that is present multiple times (-d) all records are compatable (all)
bcftools call --ploidy 2 -vm -O v | bcftools norm -O v -f ${genome} -d all - > ${results_location}/f_variants.vcf
