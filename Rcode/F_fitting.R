# Read in the data
# To get the heterozygosity frequencies, I read in the data as / seperated. Column 2 $HET contains the observed heterozygote frequencies
het_data <- read.delim("Amanita_nucl.hwe", header=TRUE, sep="")
head(het_data)
het_temp <- het_data$OBS.HOM1.HET.HOM2.
write.table(het_temp, file='hetobs.tsv', sep = '\t', quote = FALSE , row.names = FALSE)
new_het <- read.delim("hetobs.tsv", header=FALSE, sep="/")
head(new_het)

#Find the number of individuals sampled at each SNP
new_het$nind <- new_het$V1+new_het$V2+new_het$V3

# Estimate the (major/reported) allele frequency at each SNP
new_het$w_est <- (2*new_het$V1 + new_het$V2)/(2*new_het$nind)

#For each SNP, calculate the chi test values as a function of param

chisquareteststat <- function(data, param) {
  
  new_het$E1 <- new_het$nind*((1-param)*new_het$w_est*new_het$w_est+param*new_het$w_est)
  new_het$E2 <- 2*new_het$nind*(1-param)*new_het$w_est*(1-new_het$w_est)
  new_het$E3 <- new_het$nind*(1-param)*(1-new_het$w_est)*(1-new_het$w_est)+new_het$nind*param*(1-new_het$w_est)
 (new_het$V1 - new_het$E1)^2/(new_het$E1) + (new_het$V2 - new_het$E2)^2/(new_het$E2) + (new_het$V3 - new_het$E3)^2/(new_het$E3)
}

#Check that the chi-square stats look plausible
teststats <- chisquareteststat(new_het, 0.999)
head(teststats)

#I want to minimise the difference between the observed number of outlier SNPs and the expected number at a given level of confidence alpha and it's corresponding chi2 value. For a biallelic locus I have 1 degree of freedom

confidencefunc <- function(alpha, param, data){
  nSNP <- nrow(data)
  teststats <- chisquareteststat(data, param)
  if (alpha==0.975){chi2=5.024}
  else if (alpha==0.95){chi2=3.841}
  else if (alpha==0.9){chi2=2.706}
  else {chi2=0}
   nextreme <- sum(!is.na(teststats) & (teststats > chi2)) #Get the number of SNPs with significantly large test stats
     return(abs(nextreme - (1-alpha)*nSNP)) #I want to minimise the difference between the number of significant SNPs and the number expected at the level of significance testing
}


optimise(confidencefunc, c(0,1), alpha = 0.95, data = new_het, maximum = FALSE)

>$minimum
>[1] 0.03220299

>$objective
>[1] 3947
