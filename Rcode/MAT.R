## Read in sequences and tidy them up
setwd("")
library(Biostrings)
library(ape)
library(adegenet)

# Read in and clean the fasta file
raw_fasta <- readLines("HD_phasedset.a.named.fasta")
clean_fasta <- gsub("[*.]", "-", raw_fasta) # Replace '*' or '.' with '-' (these are all representing indel gaps)
writeLines(clean_fasta, "temp_HDa.fasta")
HDa <- readDNAStringSet("temp_HDa.fasta")



## Plot edit distances
# HDa contains all the sequences, without removing duplicates
seqs <- HDa

# Custom function to define distances where N matches any nucleotide.
fuzzy_distance <- function(seq1, seq2) {
  # Split each sequence into an array of characters
  s1 <- strsplit(seq1, "")[[1]]
  s2 <- strsplit(seq2, "")[[1]]
  # Sequences must be same length
  if (length(s1) != length(s2)) return(FALSE)
  
  # Compare each position allowing "N"
   x= (s1 == s2 | s1 == "N" | s2 == "N")
   sum(!x) # Return the number of sites at which the two strings differ
}


# Now we go through and calculate all pairwise distances
seq_dist <- matrix(NA, nrow = length(seqs), ncol = length(seqs))
seqchar <- as.character(seqs)


for (i in 1:length(seqs)) {
  if (i>1){
    for (j in 1:(i-1)) {
      seq_dist[i,j] <- fuzzy_distance(seqchar[i], seqchar[j])
      seq_dist[j,i] <- seq_dist[i,j]
    }}
}





# Plot the distribution of distances against the distribution of distances within a dikaryon


# Pull out the inside differences (both sequences from the same dikaryon)
inside_dists <- matrix(NA, nrow = length(seqs)/2, ncol = 1)
colnames(inside_dists) <- "Distance"
for (i in 1:(length(seqs)/2)){
  inside_dists[i] <- seq_dist[2*i-1,2*i]
}
inside_dists <- as.data.frame(inside_dists)


ggplot(dist_df, aes(x=Distance)) +
  geom_histogram(aes(y = ..density..), binwidth = 1)

ggplot(inside_dists, aes(x=Distance)) +
  geom_histogram(aes(y = ..density..), binwidth = 1)


df <- rbind(
  data.frame(value = dist_df$Distance, group = "All"),
  data.frame(value = inside_dists$Distance, group = "Dikaryon pair")
)

# Plot overlayed histograms
ggplot(df, aes(x = value, fill = group)) +
  geom_histogram(alpha = 0.5, position = "identity", aes(y = ..density..), binwidth = 1) +
  scale_fill_manual(values = c("dodgerblue", "firebrick")) +
  theme_minimal() +
  xlab("String distance") +
  ylab("Denisity")




## K-means clustering

kmers <- oligonucleotideFrequency(HDa, width = 3)

library(ggplot2)
library(reshape2)
library(adegenet)
nk <- 30
myMat <- matrix(nrow=10, ncol=nk)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(kmers, n.pca = 40, choose.n.clust = FALSE,  max.n.clust = nk)
  myMat[i,] <- grp$Kstat
}

my_df <- melt(myMat)
colnames(my_df)[1:3] <- c("Group", "K", "BIC")
my_df$K <- as.factor(my_df$K)
head(my_df)

p1 <- ggplot(my_df, aes(x = K, y = BIC))
p1 <- p1 + geom_boxplot()
p1 <- p1 + theme_bw()
p1 <- p1 + xlab("Number of groups (K)")
p1


