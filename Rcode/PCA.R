library(poppr)
library(vcfR)

vcf_file <- read.vcfR("filtered_dataset_noclone_thinned.vcf.gz")
#vcf_file <- read.vcfR("mitogenome_filtered_a2.noclone.vcf.gz")
#vcf_file <- read.vcfR("HDlocus.noclone.vcf.gz")
#vcf_file <- read.vcfR("PRlocus.noclone.vcf.gz")


variantDatasetlig <- vcfR2genlight(vcf_file)


groups  <- c("F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","T","T","T","T","T","T","T","T","T","T","T","T","T","T","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W","W")
# Groups records which population each sample belongs to. groups_cl below records the data with the inclusion of clone individuals. 

sites <- c("Flagstaff 1", "Flagstaff 2", "Sullivans Dan", "Fir Trader", "Wakari Creek")

pca_test1 <- glPca(variantDatasetlig, nf=40)
library(ggplot2)
library(ggpubr)
#scatter(pca_test1, xax = 2, yax = 3)
#loadingplot(pca_test1)
#plot(pca_test1$scores)
#p4 <- barplot(pca_test1$eig, main="eigenvalues", col=heat.colors(length(pca_test1$eig)))

mydf2  <- as.data.frame(pca_test1$eig)
mydf2$PC <- 1:43  #88 for nuclear, 25 for mito, 65 for HD, 43 for PR
head(mydf2)
colnames(mydf2) <- c("eigenvalue","PC")
p4 <- ggplot(mydf2, aes(x = PC, y = eigenvalue))+
  geom_col()
#p4

my_df1 <- as.data.frame(pca_test1$scores)
head(my_df1)
my_df1$Group <- groups
shapes <- c(21, 17, 18, 15, 25) 
p1 <- ggplot(my_df1, aes(x = PC1, y = PC2, color = Group, fill = Group, shape = Group)) 
p1 <- p1 + geom_point(size = 4)
p1 <- p1 + scale_fill_discrete(labels = sites)
p1 <- p1 + scale_color_discrete(labels = sites)
p1 <- p1 + scale_shape_manual(values = shapes, labels = sites)

p2 <- ggplot(my_df1, aes(x = PC3, y = PC2, color = Group, fill = Group, shape = Group)) 
p2 <- p2 + geom_point(size = 4)
p2 <- p2 + scale_fill_discrete(labels = sites)
p2 <- p2 + scale_color_discrete(labels = sites)
p2 <- p2 + scale_shape_manual(values = shapes, labels = sites)

p3 <- ggplot(my_df1, aes(x = PC1, y = PC3, color = Group, fill = Group, shape = Group)) 
p3 <- p3 + geom_point(size = 4)
p3 <- p3 + scale_fill_discrete(labels = sites)
p3 <- p3 + scale_color_discrete(labels = sites)
p3 <- p3 + scale_shape_manual(values = shapes, labels = sites)


PCplot <- ggarrange(p1 + rremove("legend"), p2 + rremove("legend"), p3 + rremove("legend"), p4, 
                    ncol = 2, nrow = 2,
                    labels = c("A", "B", "C", "D"),
                    common.legend = TRUE
)


PCplot
