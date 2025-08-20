outlierdf <- as.data.frame(outliertest2.het)
forests <- c("F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","F","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","H","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","S","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","T","W","W","W","W","W","W","W","W","W","W","W","W","W","W", "W","W", "W", "W", "W", "W", "W", "W","W","W", "W","W")
outlierdf$forests <- forests
aggregate(outlierdf, list(Forest = outlierdf$forests), mean)


#  Forest FID IID   O.HOM. E.HOM.    N.NM.           F forests
#1      F  NA  NA 18019.12  18140 22422.76 -0.02892918      NA
#2      H  NA  NA 18081.56  18140 22422.89 -0.01435539      NA
#3      S  NA  NA 18006.72  18140 22422.89 -0.03184000      NA
#4      T  NA  NA 18678.24  18140 22422.88  0.12505912      NA
#5      W  NA  NA 18057.08  18140 22422.69 -0.02002846      NA


modelAOV <- aov(F ~ forests,data = outlierdf)
summary(modelAOV)

#            Df Sum Sq Mean Sq F value   Pr(>F)    
#forests      4 0.3117 0.07792   7.106 5.03e-05 ***
#Residuals   91 0.9977 0.01096                     
#---
#Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1

TukeyHSD(modelAOV, "forests", ordered= FALSE, conf.level=.95)
#  Tukey multiple comparisons of means
#    95% family-wise confidence level
#
#Fit: aov(formula = F ~ forests, data = outlierdf)
#
#$forests
#            diff         lwr         upr     p adj
#H-F  0.014573788 -0.08398793  0.11313550 0.9938805
#S-F -0.002910824 -0.10147254  0.09565089 0.9999895
#T-F  0.153988294  0.05402847  0.25394812 0.0004237
#W-F  0.008900715 -0.08199806  0.09979949 0.9987680
#S-H -0.017484611 -0.11462810  0.07965887 0.9870794
#T-H  0.139414507  0.04085279  0.23797622 0.0014858
#W-H -0.005673073 -0.09503208  0.08368594 0.9997777
#T-S  0.156899118  0.05833741  0.25546083 0.0002489
#W-S  0.011811538 -0.07754747  0.10117055 0.9960278
#W-T -0.145087579 -0.23598635 -0.05418881 0.0002380

T has significant differences with all the other sites




alternatively, to handle the unequal variances:
oneway.test(F ~ forests, data = outlierdf, var.equal = FALSE)
#	One-way analysis of means (not assuming equal variances)
# data:  F and forests F = 2.621, num df = 4.000, denom df = 43.521, p-value = 0.04768
install.packages("rstatix")
library(rstatix)
games_howell_test(F ~ forests, data = outlierdf)
