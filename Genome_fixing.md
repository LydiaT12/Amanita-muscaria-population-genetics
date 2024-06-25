```
# Comparing these two, the masking isn't the (full) problem). Reads just aren't mapping properly between 524K and 545K
# This lines up perfectly with my missassembled region in SeqMonk
# There is a break in the long reads which looks like a misassembly ast ~523K. Another just under 530. Some sus behavior around 542K
# I gotta wonder what the SeqMonk stuff looks like without masking. 
# Regardless, I can't blame everything on the masking, and should look at assembling. 
#510-515K might be the MAT locus?
```



# Long reads in the bad "gap"

This takes an initial "bait" bit of sequence and assembles genome bits around it. 
What would be a good bait? Maybe a MAT sequence from another species. But maybe that isn't exactly what exists here. 
Maybe the more genuine-looking sequence on each side?
Maybe one or more long reads?

Potential long reads:
512748..529802 (17kbp) [left edge] *good choice, as it has a large left overlap
520300..534223 (13.9) [leftedge]
520734..533185 (12.4) [left]
518904..538467 (19.5) [left]

525994..546631 (20.6) [right]
535561..565318 (29.7) [right] *good choice as it has a large right overlap
534430..566984 (32.5) [right] *good choice as it has a large right overlap


Okay. How do I find these long reads?

sam format: column 4 has the left-mapping condition. 

```
awk '$4==512748' mapped_long_reads_2.sam > mitreads.txt
```
That got the entry where the left coordinate is 512748, as desired.
```
awk '$4==535561' mapped_long_reads_2.sam >> mitreads.txt
awk '$4==534430' mapped_long_reads_2.sam >> mitreads.txt
```
Yup. Deleted one entry not on Scaffold 1, where the same coordinate had turned up twise. 
I think that what I want are the three read IDs:
```
240e12e9-5bb4-44f3-a871-c6b1cf457472
6b8e6cad-c009-44ea-b2d1-14bce6a007e5
85d99cee-938b-432e-9788-26f6d971b78b
```

If I go into the fasta, each entry should start with @numbersofthereadid. I can grep the IDs that I want + the following 3 lines to get the full fastq read. 

```
grep 240e12e9-5bb4-44f3-a871-c6b1cf457472 longreads.fastq -A 4 >> gapreads.fastq
grep 6b8e6cad-c009-44ea-b2d1-14bce6a007e5 longreads.fastq -A 4 >> gapreads.fastq
grep 85d99cee-938b-432e-9788-26f6d971b78b longreads.fastq -A 4 >> gapreads.fastq
```

Now gapreads.fastq is a fastq file containing the three long reads which justify the filling of this gap. 




# MITObim
<!--
Docker/singularity commands:

Sooo... to make this docker run as a singularity, I first load singularity
```module load Singularity```
Set the directory to work in
```WORKING_DIR=/home/turly826/nobackup/nanopore_reference_genome/long_read_mapping/mitbim_files```
Download and convert the Docker container as a Singularity image:
```singularity pull mitobim.sif docker://chrishah/mitobim```
Give the container access to the directory:
```singularity run --bind $PWD mitobim.sif```


```
module load Singularity
WORKING_DIR=/home/turly826/nobackup/nanopore_reference_genome/long_read_mapping/mitbim_files
singularity pull mitobim.sif docker://chrishah/mitobim # this should be a set-up-only step (I believe)
singularity run --bind $PWD mitobim.sif
```

The image had a lot of empty files. The links to the docker images were broken. idk, maybe the image is broken?


## Singularity 2
-->
```
module load Singularity
singularity pull docker://chrishah/mitobim
```

```
singularity inspect mitobim_latest.sif
#Data:
org.label-schema.build-arch: amd64
org.label-schema.build-date: Wednesday_12_June_2024_14:32:30_NZST
org.label-schema.schema-version: 1.0
org.label-schema.usage.singularity.deffile.bootstrap: docker
org.label-schema.usage.singularity.deffile.from: chrishah/mitobim
org.label-schema.usage.singularity.version: 3.11.3
```

Do I ... just run it?
```
singularity run mitobim_latest.sif 
MITObim.pl
  MITObim.pl -sample MITOtest -ref Amamus_LR1 -readpool ../../trimmed_H12_1.fastq --quick gapread.fasta -end 10 --clean &> log_test1

# -readpool: are the pool of short reads to assemble with. These need to be interleaved (There should be a way to convert)
# --quick: the bait read (I'm using a single carefully chosen long read)
# -end: number of iterations to do
# -clean: discard intermediate steps (sure)
# -sample and -ref: ID names
```

Ran 5 iterations in <2 hours (not a huge amount less). 

Loads of outputs. final assembly written to
```
/nesi/nobackup/uoo03761/nanopore_reference_genome/long_read_mapping/mitbim_files/iteration5/MITOtest-Amamus_LR1-it5_noIUPAC.fasta
```
readpool contains 222904 reads
assembly contains 1 contig(s)
contig length: 17779

## What next?

MITObim assembled something. I don't know if it's any good or if it fills the gap. 

1. I want to compare it with the genome assembly - just map this little assembly agains scaffold 1. I can do this in MUMmer and plot it.
2. I want to assemble MITObim bits with the other starting baits
3. I want to compare the four assemblies; three from mitobim plus scaffold 1. I wanna align them all first. Maybe with MUMmer, and look at the overlaps in table format?



