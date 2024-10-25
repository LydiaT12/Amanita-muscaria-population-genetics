#/nesi/nobackup/uoo03761/whole_genome_resequencing/01.VariantCalling/10.MATphasing/phased_out/temp
bgzip HD.S1.phased.vcf #compress variants
tabix -p vcf HD.S1.phased.vcf.gz # index them
cat ../../Amanita_scaffold1.fasta | vcf-consensus -i -H 1 -s S1 HD.S1.phased.vcf.gz  > temp.fa #vcftools tool to create fasta from phased vcf H1 is first haplotype
samtools faidx temp.fa index fasta
samtools faidx temp.fa Scaffold1_Amanita_muscaria:10-10 #Get position 10 in fasta (1-10 would be one to ten)
