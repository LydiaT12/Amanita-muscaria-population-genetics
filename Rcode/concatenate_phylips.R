# Concatenate the phylip files of three gene into one phylip file with all (up to three) genes aligned and concatenated for each individual


setwd("C:/Users/Lydia/OneDrive - University of Otago/Documents/PhD/PhD - Botany/R/treefiles")

library(ape)

# Read sequences
gene1 <- read.dna("LSU.phy", format="interleaved")
gene2 <- read.dna("ITS_out.phy", format="interleaved")
gene3 <- read.dna("betatubulin_out.phy", format="interleaved")

new_names_LSU <- c("New_genome","GAL15330","GAL16775","GAL15776","GAL4247","GAL14284","GAL16735","GAL2814","GAL5895","GAL15453","GAL15461","GAL8950","GAL5900","GAL5946","GAL2810","GAL3169","GAL3688","GAL16654","GAL5505","GAL6027","GAL4302")
rownames(gene1) <- new_names_LSU

new_names_ITS <- c("A31452","A31445","A80048","A1539","A4220","A45863","A45843","A45785","A45820","A45840","A45883","A49100","A44761","A45060","A30961","A30962","A30963","A30964","A30965","A30976","A30977","A30978","A30985","A30981","A30982","A30987","New_genome","GAL15330","GAL16775","GAL15776","GAL4247","GAL14284","GAL16735","GAL2814","GAL5895","GAL15453","GAL15461","GAL8950","GAL5900","GAL5946","GAL2810","GAL3169","GAL3688","GAL16654","GAL5505","GAL6027","GAL4302","pantherina")  # New names in order
rownames(gene2) <- new_names_ITS

new_names_beta <- c("A506","A1539","A4220","A30978","A30981","A30982","A45863","A45843","A45785","A45820","A45840","A45883","A94100","A44761","A45060","A30961","A30962","A30963","A30976","A30977","A30985","A30964","A30965","A31452","A31445","A80048","A30987","New_genome","GAL15330","GAL16775","GAL15776","GAL4247","GAL14284","GAL16735","GAL2814","GAL5895","GAL15453","GAL15330","GAL16654","GAL5505","GAL4302","pantherina")
rownames(gene3) <- new_names_beta


# Find max length per gene
max_lens <- c(ncol(gene1), ncol(gene2), ncol(gene3))

# Pad missing taxa with gaps
all_names <- unique(c(rownames(gene1), rownames(gene2), rownames(gene3)))
pad_with_gaps <- function(data, max_len) {
  full_data <- matrix("-", nrow=length(all_names), ncol=max_len)
  rownames(full_data) <- all_names
  full_data[rownames(data), ] <- as.character(data)
  return(full_data)
}

gene1 <- pad_with_gaps(gene1, max_lens[1])
gene2 <- pad_with_gaps(gene2, max_lens[2])
gene3 <- pad_with_gaps(gene3, max_lens[3])

# Concatenate
concatenated <- cbind(gene1, gene2, gene3)

# Save as PHYLIP
write.dna(concatenated, file="concatenated_renamed.phy", format="interleaved")
