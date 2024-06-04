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




module load BCFtools
module load SAMtools
bgzip A_muscaria_guessowii_aligned.vcf 
bcftools index A_muscaria_guessowii_aligned.vcf.gz

bgzip A_muscaria_old_aligned.vcf
bgzip A_phalloides_aligned.vcf 
bcftools index A_muscaria_old_aligned.vcf.gz 
bcftools index A_phalloides_aligned.vcf.gz
