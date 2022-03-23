#this is sooo sloppy... needs a lot of work

#produce a multipsample variant call taht contains all the information on all samples for females

genome=vcf/data/ref_genome/acar.fasta

f_input_location=vcf/results/bam/all_only_female

f_results_location=vcf/results/vcf/all_only_female

# output type (-O) uncompressed VCF (v)
bcftools mpileup -O v -f ${genome} ${f_input_location}/*.sorted.bam > ${f_results_location}/f_genotypes.vcf

# norm: output only first reference of record that is present multiple times (-d) all records are compatable (all)
bcftools call --ploidy 2 -vm -O v ${f_results_location}/f_genotypes.vcf | bcftools norm -O v -f ${genome} -d all - > ${f_results_location}/f_variants.vcf

#produce a multipsample variant call taht contains all the information on all samples for males

genome=vcf/data/ref_genome/acar.fasta

m_input_location=vcf/results/bam/all_only_male

m_results_location=vcf/results/vcf/all_only_male

# output type (-O) uncompressed VCF (v)
bcftools mpileup -O v -f ${genome} ${m_input_location}/*.sorted.bam > ${m_results_location}/m_genotypes.vcf

# norm: output only first reference of record that is present multiple times (-d) all records are compatable (all)
bcftools call --ploidy 2 -vm -O v ${m_results_location}/m_genotypes.vcf | bcftools norm -O v -f ${genome} -d all - > ${m_results_location}/m_variants.vcf

#explore the vcfs
#total number of SNPS
#males
bcftools view -v snps ${m_results_location}/m_variants.vcf | grep -v "^#" | wc -l
#output: 333809

# total number of unique positions, indicating that several sites have two or more alternate alleles
bcftools view -v snps ${m_results_location}/m_variants.vcf | grep -v "^#" | cut -f2 | sort -u | wc -l
#output 329686

#output 304491
#females
bcftools view -v snps ${f_results_location}/f_variants.vcf | grep -v "^#" | wc -l
#output 307677

# total number of unique positions, indicating that several sites have two or more alternate alleles
bcftools view -v snps ${f_results_location}/f_variants.vcf | grep -v "^#" | cut -f2 | sort -u | wc -l
#output 304491



#Can I create a merged bam file for all female genotypes? output type (-O) uncompressed BAM (b)
bcftools mpileup -O b -f ${genome} ${f_input_location}/*.sorted.bam > ${f_input_location}/f_alignment.bam

bcftools index ${f_input_location}/f_alignment.bam