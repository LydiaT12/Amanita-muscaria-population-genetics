# Outgroups

I grabbed the assemblies of A. phalloides (2023 version), A. muscaria (old) and A. muscaria var guessowii from NCBI.

## Make vcf files with the outgroup genomes

I have three (potential) outgroup genomes in fasta format. I want to compare these with my genome assembly and create vcf files which can be combined with my main dataset. 

Download/compile GSAlign:

```
module load Miniconda3
conda create -n GSAenvt
source activate GSAenvt
conda install gsalign
conda install -c conda-forge -c bioconda gsalign
git clone https://github.com/hsinnan75/GSAlign.git
cd GSAlign/
make
```

I'm going to align to the unmasked genome and apply the mask before combining with wider dataset. First step is to index the reference. 

```
bin/bwt_index A_muscaria_genome_2024.fasta A_muscaria_index
```

Now I do the alignments:

```
bin/GSAlign -i A_muscaria_index -q ../A_muscaria_old_genome.fasta -o A_muscaria_old_aligned
# Summary output: GSAlign identifies 784053 SNVs, 52790 insertions, and 54108 deletions [A_muscaria_old_aligned.vcf].

bin/GSAlign -i A_muscaria_index -q ../A_muscaria_guessowii_genome.fasta -o A_muscaria_guessowii_aligned
# Summary output: GSAlign identifies 1844505 SNVs, 130523 insertions, and 135712 deletions [A_muscaria_guessowii_aligned.vcf]

bin/GSAlign -i A_muscaria_index -q ../A_phalloides_genome.fasta -o A_phalloides_aligned
# Summary output: GSAlign identifies 22229 SNVs, 699 insertions, and 732 deletions [A_phalloides_aligned.vcf].
```

The summary has fewer differences with phalloides than either of the old muscaria genomes? That must surely be an artefact of misassemblies. 





tail of a vcf:
```
Scaffold240_Amanita_muscaria    785     .       G       A       100     *       TYPE=SUBSTITUTE
Scaffold240_Amanita_muscaria    787     .       T       C       100     *       TYPE=SUBSTITUTE
Scaffold240_Amanita_muscaria    795     .       C       CA      100     *       TYPE=INSERT
```

I was getting the error message "FILTER '*' is not defined in the header"

```
sed -i "s/"*"/PASS/" A_muscaria_old_aligned.copy.vcf
```

```
module load BCFtools
module load SAMtools

bgzip A_muscaria_guessowii_aligned.copy.vcf 
bcftools index A_muscaria_guessowii_aligned.copy.vcf.gz

bgzip A_muscaria_old_aligned.copy.vcf
bgzip A_phalloides_aligned.copy.vcf 
bcftools index A_phalloides_aligned.copy.vcf.gz 
bcftools index A_muscaria_old_aligned.vcf.gz
```

<!--
```
bcftools norm --check-ref x --fasta-ref GSAlign/A_muscaria_genome_2024.fasta A_muscaria_old_aligned.copy.vcf.gz > A_muscaria_old.aligned.copy.norm.vcf.gz
# Output: Lines   total/split/joined/realigned/skipped:   890951/0/0/12861/200408
bcftools norm --check-ref x --fasta-ref GSAlign/A_muscaria_genome_2024.fasta A_phalloides_aligned.copy.vcf.gz > A_phalloides.aligned.copy.norm.vcf.gz
# Lines   total/split/joined/realigned/skipped:   23660/0/0/22/14450
bcftools norm --check-ref x --fasta-ref GSAlign/A_muscaria_genome_2024.fasta A_muscaria_guessowii_aligned.copy.vcf.gz  > A_muscaria_guessowii.aligned.copy.norm.vcf.gz
# Lines   total/split/joined/realigned/skipped:   2110740/0/0/24427/915876
```
-->

```
bcftools norm --check-ref x --fasta-ref GSAlign/A_muscaria_genome_2024.fasta A_muscaria_old_aligned.copy.vcf.gz -Oz > A_muscaria_old.2.aligned.copy.norm.vcf.gz
# Output: Lines   total/split/joined/realigned/skipped:   890951/0/0/12861/200408
bcftools norm --check-ref x --fasta-ref GSAlign/A_muscaria_genome_2024.fasta A_phalloides_aligned.copy.vcf.gz -Oz > A_phalloides.2.aligned.copy.norm.vcf.gz
# Lines   total/split/joined/realigned/skipped:   23660/0/0/22/14450
bcftools norm --check-ref x --fasta-ref GSAlign/A_muscaria_genome_2024.fasta A_muscaria_guessowii_aligned.copy.vcf.gz  -Oz > A_muscaria_guessowii.2.aligned.copy.norm.vcf.gz
# Lines   total/split/joined/realigned/skipped:   2110740/0/0/24427/915876
```

```
bcftools index A_phalloides.aligned.copy.norm.vcf.gz 
bcftools index A_muscaria_old.aligned.copy.norm.vcf.gz
bcftools index A_muscaria_guessowii.aligned.copy.norm.vcf.gz


bcftools index A_muscaria_guessowii.2.aligned.copy.norm.vcf.gz 


bcftools merge A_muscaria_old.2.aligned.copy.norm.vcf.gz A_phalloides.2.aligned.copy.norm.vcf.gz -Oz -o merged2.vcf.gz
bcftools index merged2.vcf.gz 
bcftools merge merged2.vcf.gz A_muscaria_guessowii.2.aligned.copy.norm.vcf.gz -Oz -o merged3.vcf.gz
```

<!--
```
bcftools merge A_muscaria_old.aligned.copy.norm.vcf.gz A_phalloides.aligned.copy.norm.vcf.gz -Oz -o merged.vcf.gz --no-index
```
It seemed to produce something


Tried adding -Oz to the norm function. 
bcftools merge A_muscaria_old.2.aligned.copy.norm.vcf.gz A_phalloides.2.aligned.copy.norm.vcf.gz -Oz -o merged2.vcf.gz
> Failed to open A_muscaria_old.2.aligned.copy.norm.vcf.gz: could not load index

tried indexing again
bcftools index A_muscaria_old.2.aligned.copy.norm.vcf.gz 
[turly826@wbn001 8.outgroups]$ bcftools index A_phalloides.2.aligned.copy.norm.vcf.gz 
[turly826@wbn001 8.outgroups]$ bcftools merge A_muscaria_old.2.aligned.copy.norm.vcf.gz A_phalloides.2.aligned.copy.norm.vcf.gz -Oz -o merged2.vcf.gz
-->


So the file merged3.vcf.gz contains the vcf's for all three outgroup genomes. 

Index it:
```
bcftools index merged3.vcf.gz
```


Moce folders. Subset just the first 28 contigs
```
bcftools view -Oz -r Scaffold1_Amanita_muscaria,Scaffold2_Amanita_muscaria,Scaffold3_Amanita_muscaria,Scaffold4_Amanita_muscaria,Scaffold5_Amanita_muscaria,Scaffold6_Amanita_muscaria,Scaffold7_Amanita_muscaria,Scaffold8_Amanita_muscaria,Scaffold9_Amanita_muscaria,Scaffold10_Amanita_muscaria,Scaffold11_Amanita_muscaria,Scaffold12_Amanita_muscaria,Scaffold13_Amanita_muscaria,Scaffold14_Amanita_muscaria,Scaffold15_Amanita_muscaria,Scaffold16_Amanita_muscaria,Scaffold17_Amanita_muscaria,Scaffold18_Amanita_muscaria,Scaffold19_Amanita_muscaria,Scaffold20_Amanita_muscaria,Scaffold21_Amanita_muscaria,Scaffold22_Amanita_muscaria,Scaffold23_Amanita_muscaria,Scaffold24_Amanita_muscaria,Scaffold25_Amanita_muscaria,Scaffold26_Amanita_muscaria,Scaffold27_Amanita_muscaria,Scaffold28_Amanita_muscaria ../../01.VariantCalling/8.outgroups/merged3.vcf.gz  > outgroups_1to28.vcf.gz
```
cp ../../02.vcftools/scaffolds_1to28_filtered.vcf.gz samples_1to28.vcf.gz
Now I index both 

bcftools index outgroups_1to28.vcf.gz 
bcftools index samples_1to28.vcf.gz


