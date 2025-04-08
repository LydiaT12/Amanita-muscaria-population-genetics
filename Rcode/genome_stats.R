## CALCULATE FSFS
# I want count data. These have already been arranged so that COUNT1 is always the smaller value

## Function to arrange data in count_data into a new dataframe with COUNT1 <= COUNT2
#sfscounts_normalised <- count_data
#for (ii in 1:300000){
#  if (count_data$COUNT1[ii] > count_data$COUNT2[ii]){
#    sfscounts_normalised$COUNT1[ii] <- count_data$COUNT2[ii]
#    sfscounts_normalised$COUNT2[ii] <- count_data$COUNT1[ii]
#  }
#}


count_data <- read.delim("sfscounts_normalised.tsv", header=TRUE, sep="")
head(count_data)
count_data <- count_data[count_data$COUNT1>0, ]


bins <- 0:89

# Get the number of observations of each allele count
hist_info <- hist(count_data$COUNT1, breaks = bins, plot = FALSE)
#hist_info

# Define the function that the observations are expected to follow, and its RSS
WFexpect_dist <- function(theta, count) {
  theta/count + theta/(89*2-count)
}
RSS <- function(theta) {
  bb = 0
  for (ii in 2:89){
  bb = bb+(WFexpect_dist(theta,ii)-hist_info$counts[ii])^2
  }
  return(bb)
}

# Find the optimal value of theta. This gives 39269 when including singletons and 49931 when singletons are excluded. 
optim(40000, RSS, method="CG")


library(ggplot2)

ggplot(count_data, aes(x = COUNT1)) +
  geom_histogram(bins=89) +#there are 89 individuals (having removed dup clones), so  the SFS should run from 1 to 89
   stat_function(fun = WFexpect_dist, args = list(theta = 49931), 
                color = "blue", linewidth = 1) +
  xlab("Minor allele count") + ylab("Frequency")







## Calculate distribution of f

# I have a dataset with columns named:
# Chromosome    Position    Observed frequencies of homozygous allele 1/heterozygous/homozygous allele 2    Expected homo/heterozygous freqs under HWE   Chi2 value for HWE    P for following a HWE distribution

library(ggplot2)

# To get the heterozygosity frequencies, I read in the data as / seperated. Column 2 $HET contains the observed heterozygote frequencies
het_data2 <- read.delim("Amanita_nucl.hwe", header=TRUE, sep="/")
head(het_data2)


## Now calculate f
#$HET is the observed heterozygote count, $HET.1 is the expected heterozygote count
nrow <- 307380
het_data2$f <- replicate(nrow, 0)

het_data2$f <- 1- (het_data2$HET / het_data2$HET.1)
head(het_data2)
 het_data2 <- het_data2[het_data2$HET.1 > 2, ] # This can filter out low-count alleles

ggplot(het_data2, aes(x = f)) +
  geom_histogram(binwidth=0.02)

  
  
