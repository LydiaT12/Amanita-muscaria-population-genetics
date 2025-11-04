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
  
  
    (new_het$V1 - new_het$nind*((1-param)*new_het$w_est*new_het$w_est+param*new_het$w_est))^2/(new_het$nind*(1-param)*new_het$w_est*new_het$w_est+new_het$nind*param*new_het$w_est) + (new_het$V2 - 2*new_het$nind*(1-param)*new_het$w_est*(1-new_het$w_est))^2/(2*new_het$nind*(1-param)*new_het$w_est*(1-new_het$w_est)) + (new_het$V3 - new_het$nind*(1-param)*(1-new_het$w_est)*(1-new_het$w_est)+new_het$nind*param*(1-new_het$w_est))^2/(new_het$nind*(1-param)*(1-new_het$w_est)*(1-new_het$w_est)+new_het$nind*param*(1-new_het$w_est))

}

#Check that the chi-square stats look plausible
teststats <- chisquareteststat(new_het, 0.999)
head(teststats)

#I want to minimise the difference between the observed number of outlier SNPs and the expected number at a given level of confidence alpha and it's corresponding chi2 value. For a biallelic locus I have 1 degree of freedom

confidencefunc <- function(alpha, param, data){
#confidencefunc <- function(param){
 # alpha <- 0.975
  #data <- new_het
  nSNP <- nrow(data)
  teststats <- chisquareteststat(data, param)
  if (alpha==0.975){chi2=0.001}
  else if (alpha==0.95){chi2=0.004}
  else if (alpha==0.9){chi2=0.016}
  else {chi2=0}
   nextreme <- sum(!is.na(teststats) & teststats > chi2) #Get the number of SNPs with significantly large test stats
  return(abs(nextreme - (1-alpha)*nSNP)) #I want to minimise the difference between the number of significant SNPs and the number expected at the level of singificance testing
}


optimise(confidencefunc, c(-1,10), alpha = 0.95, data = new_het)

>$minimum
>[1] 1.142099
>$objective
>[1] 415
