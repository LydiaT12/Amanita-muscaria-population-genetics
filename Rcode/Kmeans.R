setwd("C:/Users/Lydia/OneDrive - University of Otago/Documents/PhD/PhD - Botany/R")
library(poppr)
library(vcfR)
#library(Biostrings)
library(ape)
## Read in the appropriate dataset and forest groups:
vcf_file <- read.vcfR("filtered_dataset_noclone_thinned.vcf.gz") #nuclear genome, no clone duplicates
# vcf_file <- read.vcfR("mitogenome_filtered_a2.noclone.vcf.gz") #mitogenome, no clone duplicates
variantDataset <- vcfR2genlight(vcf_file) #convert the file to a genlight
# variantDataset <- fasta2genlight("HD_phasedset.c.fasta") #HD phased region c (all individuals)

#forest_cl <- c("F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W") #forest sites with all clones
Forest <- c("F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","T","T","T","T","T","T","T","T","T","T","T","T","T","T","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W") #forest sites for the clone-corrected dataset



# Try K-mean clustering to look for population structure. First K=3
# fileclust <- find.clusters(variantDataset, n.pca = 20) #adegenet - choose >40 PCs
#Chose 40 PCs and 3 clusters



## BIC plot - I skip this



## DAPC - do and visualise the groupings for k = 2,3,4,5

library(ggplot2)
library(reshape2)
my_k <- 2:5   #values that I try # for nuclear genome
# my_k <- seq(2, 8, by = 2) # values for mitogenome

grp_l <- vector(mode = "list", length = length(my_k))
dapc_l <- vector(mode = "list", length = length(my_k))

for(i in 1:length(dapc_l)){
  # set.seed(9)
  grp_l[[i]] <- find.clusters(variantDataset, n.pca = 40, n.clust = my_k[i])
  dapc_l[[i]] <- dapc(variantDataset, pop = grp_l[[i]]$grp, n.pca = 40, n.da = my_k[i])
}

#grp_l contains the group assignment info for each k
#dapc_l contains the info needed to scatterplot the groups

my_pal <- RColorBrewer::brewer.pal(n=8, name = "Dark2")


my_df <- as.data.frame(dapc_l[[2]]$ind.coord)
my_df$Group <- dapc_l[[2]]$grp
my_df$Forest <- Forest
head(my_df)

p2 <- ggplot(my_df, aes(x = LD1, y = LD2, color = Group, fill = Group, shape = Forest))
p2 <- p2 + geom_point(size = 4)
p2 <- p2 + theme_bw()
#p2 <- p2 + scale_color_manual(values=c(my_pal))
#p2 <- p2 + scale_fill_manual(values=c(paste(my_pal, "66", sep = "")))
p2




## Make dataframes

tmp <- as.data.frame(dapc_l[[1]]$posterior)
tmp$K <- my_k[1]
tmp$Isolate <- rownames(tmp)
tmp <- melt(tmp, id = c("Isolate", "K"))
names(tmp)[3:4] <- c("Group", "Posterior")
tmp$Region <- Forest
#head(tmp)
my_df <- tmp

for(i in 2:length(dapc_l)){
  tmp <- as.data.frame(dapc_l[[i]]$posterior)
  tmp$K <- my_k[i]
  tmp$Isolate <- rownames(tmp)
  tmp <- melt(tmp, id = c("Isolate", "K"))
  names(tmp)[3:4] <- c("Group", "Posterior")
  tmp$Region <- Forest
  
  my_df <- rbind(my_df, tmp)
}





## Plot group assignment

grp.labs <- paste("K =", my_k)
names(grp.labs) <- my_k

p3 <- ggplot(my_df, aes(x = Isolate, y = Posterior, fill = Group))
p3 <- p3 + geom_bar(stat = "identity")
p3 <- p3 + facet_grid(K ~ Region, scales = "free_x", space = "free", 
                      labeller = labeller(K = grp.labs))
p3 <- p3 + theme_bw()
p3 <- p3 + ylab("Posterior membership probability")
p3 <- p3 + theme(legend.position='none')
p3 <- p3 + scale_fill_manual(values=c(my_pal))
#p3 <- p3 + theme(axis.text.x = element_text(angle = 90, hjust = 1, size = 8))
#p3 <- p3 + theme(axis.text.x = element_blank())
p3




### Plot BIC

library(ggplot2)
library(reshape2)
library(adegenet)
nk <- 8 # Number of K to try. Values: 8 for nuclear, 8 for mitogenome, 40 for HDc
myMat <- matrix(nrow=10, ncol=nk)
colnames(myMat) <- 1:ncol(myMat)
for(i in 1:nrow(myMat)){
  grp <- find.clusters(variantDataset, n.pca = 40, choose.n.clust = FALSE,  max.n.clust = nk)
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



# Make pretty joint plot
library(ggpubr)
ggarrange(ggarrange(p1,
                    p2,
                    ncol = 2, labels = c("A", "B")),
          p3,
          nrow = 2,
          labels = c("", "C"),
          heights = c(2, 3)
)
