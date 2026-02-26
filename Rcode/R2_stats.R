## Plot and fit with r^2 data

## Import data

LD_data1 <- read.delim("scaffold36.geno.ld", header=TRUE, sep="")
LD_data$dist <- abs(LD_data$POS1 - LD_data$POS2)
LD_data$r2 <- LD_data$R.2


## Plot values

#Plot binned averages over top
library(ggplot2)
(ggplot(LD_data, aes(x=dist,y=r2)) +
  geom_point(alpha = 0.4) +
  stat_summary_bin(fun.y='mean', bins=100,
                   color='orange', size=2, geom='point')) +
   xlab("Distance (bp)") +
  ylab(expression(r^2))

#I'll restrict the dataset to distance <3000 to just zoom in on the curvey bit. 

LD_data3 <- LD_data[LD_data$dist < 3000, ]
(ggplot(LD_data3, aes(x=dist,y=r2)) +
  geom_point(alpha = 0.4) +
  stat_summary_bin(fun.y='mean', bins=100,
                   color='orange', size=2, geom='point')) +
   xlab("Distance (bp)") +
  ylab(expression(r^2))




## Now fit to 1/(1+a*dist)
 



 df2 <- aggregate(LD_data, #the data frame
                 by=list(cut(LD_data$dist,breaks=seq(min(LD_data$dist), max(LD_data$dist), by = 100))), #the bins (see below)
                 mean)

head(df2)

# fit model
model <- nls(r2 ~ 1/(1+a * dist), data = df2, start = list(a = 10))
summary(model)
# This estimates a = 0.002



#I'll restrict the dataset to distance <10000 to just zoom in on the curvey bit. 

LD_data3 <- LD_data[LD_data$dist < 10000, ]

df3 <- aggregate(LD_data3, #the data frame
                 by=list(cut(LD_data3$dist,breaks=500)), #the bins (see below)
                 mean)

head(df3)

# fit model
model3 <- nls(r2 ~ 1/(1+a * dist), data = df3, start = list(a = 10))
summary(model3)
# This estimates a = 0.004
#These estimates differ by a factor of 2, which is arguably big - but also they're of the same order of magnitude. 
