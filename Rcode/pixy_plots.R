library(ggplot2)
library(tidyverse)
library(ggpubr)

# Provide path to input. Can be pi or Dxy.
inp<-read.table("1to28win5000_pi.txt",sep="\t",header=T)

# Find the chromosome names and order them: first numerical order, then any non-numerical chromosomes
#   e.g., chr1, chr2, chr22, chrX
chroms <- unique(inp$chromosome)
chrOrder <- sort(chroms)
inp$chrOrder <- factor(inp$chromosome,levels=chrOrder)



#####
# Plot pi for each population found in the input file
# Saves a copy of each plot in the working directory
if("avg_pi" %in% colnames(inp)){
  pops <- unique(inp$pop)
  for (p in pops){
    thisPop <- subset(inp, pop == p)
    # Plot stats along all chromosomes:
    popPlot <- ggplot(thisPop, aes(window_pos_1, avg_pi, color=chrOrder)) +
      geom_point()+
      facet_grid(. ~ chrOrder)+
      labs(title=paste("Pi for population", p))+
      labs(x="Position of window start", y="Pi")+
      scale_color_manual(values=rep(c("black","gray"),ceiling((length(chrOrder)/2))))+
      theme_classic()+
      theme(legend.position = "none")
    ggsave(paste("piplot_", p,".png", sep=""), plot = popPlot, device = "png", dpi = 300)
  }
} else {
  print("Pi not found in this file")
}

# Plot Dxy for each combination of populations found in the input file
# Saves a copy of each plot in the working directory
if("avg_dxy" %in% colnames(inp)){
  # Get each unique combination of populations
  pops <- unique(inp[c("pop1", "pop2")])
  for (p in 1:nrow(pops)){
    combo <- pops[p,]
    thisPop <- subset(inp, pop1 == combo$pop1[[1]] & pop2 == combo$pop2[[1]])
    # Plot stats along all chromosomes:
    popPlot <- ggplot(thisPop, aes(window_pos_1, avg_dxy, color=chrOrder)) +
      geom_point()+
      facet_grid(. ~ chrOrder)+
      labs(title=paste("Dxy for", combo$pop1[[1]], "&", combo$pop2[[1]]))+
      labs(x="Position of window start", y="Dxy")+
      theme(legend.position = "none")+
      scale_color_manual(values=rep(c("black","gray"),ceiling((length(chrOrder)/2))))+
      theme_classic()+
      theme(legend.position = "none")
    ggsave(paste("dxyplot_", combo$pop1[[1]], "_", combo$pop2[[1]],".png", sep=""), plot = popPlot, device = "png", dpi = 300)
  }
} else {
  print("Dxy not found in this file")
}











###
## Re-order the three pixy files into one dataframe pixy_df (for scaffold 1 only. Remove the .1 to read in other files)




pixy_to_long <- function(pixy_files){
  
  pixy_df <- list()
  
  for(i in 1:length(pixy_files)){
    
    stat_file_type <- gsub(".*_|.txt", "", pixy_files[i])
#    stat_file_type <- gsub(".*_|.1.txt", "", pixy_files[i])
    
    if(stat_file_type == "pi"){
      
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value") %>%
        rename(pop1 = pop) %>%
        mutate(pop2 = NA)
      
      pixy_df[[i]] <- df
      
      
    } else{
      
      df <- read_delim(pixy_files[i], delim = "\t")
      df <- df %>%
        gather(-pop1, -pop2, -window_pos_1, -window_pos_2, -chromosome,
               key = "statistic", value = "value")
      pixy_df[[i]] <- df
      
    }
    
  }
  
  bind_rows(pixy_df) %>%
    arrange(pop1, pop2, chromosome, window_pos_1, statistic)
  
}

#pixy_files <- c("pixy_dxy.1.txt", "pixy_pi.1.txt", "pixy_fst.1.txt")
#pixy_df <- pixy_to_long(pixy_files)

pixy_m1000 <- c("mitowin1000_dxy.txt", "mitowin1000_pi.txt", "mitowin1000_fst.txt")
pixy_df2 <- pixy_to_long(pixy_m1000)



###
## use the pixy_df
head(pixy_df)


# custom labeller for special characters in pi/dxy/fst
pixy_labeller <- as_labeller(c(avg_pi = "pi",
                               avg_dxy = "D[XY]",
                               avg_wc_fst = "F[ST]"),
                             default = label_parsed)


### plotting summary statistics along a single chromosome
p2 <- pixy_df2 %>%
#  filter(chromosome == 1) %>%
  filter(statistic %in% c("avg_pi", "avg_dxy", "avg_wc_fst")) %>%
  mutate(chr_position = ((window_pos_1 + window_pos_2)/2)/1000000) %>%
  ggplot(aes(x = chr_position, y = value, color = statistic))+
 # geom_line(size = 0.25)+
  geom_point(size = 0.5)+
  facet_grid(statistic ~ .,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value))+
  xlab("Position on Chromosome (Mb)")+
  ylab("Statistic Value")+
  theme_bw()+
  theme(panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position = "none",
        strip.text = element_text(size = 14))+ #this changes the size of the label text
  scale_x_continuous(expand = c(0, 0))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1))+
  scale_color_brewer(palette = "Set1")

p2

## Plotting side-by-side different windows
pixyplot <- ggarrange(p1 , p2, 
                    ncol = 1, nrow = 2,
                    labels = c("A", "B")
)

pixyplot



### 
## Plotting distribution of pi

#pi_in<-read.table("1to28win5000_pi.txt",sep="\t",header=T)
pi_in<-read.table("mitowin1000_pi.txt",sep="\t",header=T)
pops <- unique(pi_in$pop)
for (p in pops){
  thisPop <- subset(pi_in, pop == p)
#  thisPop_df <- as.df(thisPop)
  p_tmp <- ggplot(thisPop, aes(x = avg_pi)) +
    geom_histogram(binwidth=0.005)

}
head(pi_in)



piF<- subset(pi_in, pop == "F")
piF <- piF[is.finite(piF$avg_pi)==TRUE, ]
pf <- ggplot(piF, aes(x = avg_pi)) +
  geom_histogram(binwidth=0.005)
pi_estF <- mean(piF$avg_pi)
pi_estF #0.0906178 #0.0916449 (5000 1-28, 1000 1)

piH<- subset(pi_in, pop == "H")
piH <- piH[is.finite(piH$avg_pi)==TRUE, ]
ph <- ggplot(piH, aes(x = avg_pi)) +
  geom_histogram(binwidth=0.005)
pi_estH <- mean(piH$avg_pi)
pi_estH # 0.09011165 #0.0944988


piS<- subset(pi_in, pop == "S")
piS <- piS[is.finite(piS$avg_pi)==TRUE, ]
ps <- ggplot(piS, aes(x = avg_pi)) +
  geom_histogram(binwidth=0.005)
pi_estS <- mean(piS$avg_pi)
pi_estS #0.09200834 #0.09247049


piT<- subset(pi_in, pop == "T")
piT <- piT[is.finite(piT$avg_pi)==TRUE, ]
pt <- ggplot(piT, aes(x = avg_pi)) +
  geom_histogram(binwidth=0.005)
pi_estT <- mean(piT$avg_pi)
pi_estT #0.07696297 #0.08258985


piW<- subset(pi_in, pop == "W")
piW <- piW[is.finite(piW$avg_pi)==TRUE, ]
head(piW)
pw <- ggplot(piW, aes(x = avg_pi)) +
  geom_histogram(binwidth=0.005)
pi_estW <- mean(piW$avg_pi)
pi_estW #0.09050138 #0.09164038

pi <- pi_in[is.finite(pi_in$avg_pi)==TRUE, ]
p <- ggplot(pi, aes(x = avg_pi)) +
  geom_histogram(binwidth=0.005)
pi_est <- mean(pi$avg_pi)
pi_est # 0.08804043 #0.0905689

library(ggplot2)
library(ggpubr)


piplot <- ggarrange(pf, ph, ps, pt, pw, p,
                    ncol = 2, nrow = 3,
                    labels = c("Flagstaff 1", "Flagstaff 2", "Sullivans Dam", "Fir Trader", "Wakari Creek", "Combined")
)

piplot





res_anova <- aov(avg_pi ~ pop, data = pi_in)
summary(res_anova)

#Df Sum Sq Mean Sq F value Pr(>F)    
#pop             4  0.546 0.13640   82.96 <2e-16 ***
#  Residuals   17545 28.846 0.00164                   
#---
#  Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
#2280 observations deleted due to missingness


head(pixy_df2)

scaffoldnames <- c("Scaffold1_Amanita_muscaria","Scaffold2_Amanita_muscaria","Scaffold3_Amanita_muscaria","Scaffold4_Amanita_muscaria","Scaffold5_Amanita_muscaria","Scaffold6_Amanita_muscaria","Scaffold7_Amanita_muscaria","Scaffold8_Amanita_muscaria","Scaffold9_Amanita_muscaria","Scaffold10_Amanita_muscaria","Scaffold11_Amanita_muscaria","Scaffold12_Amanita_muscaria","Scaffold13_Amanita_muscaria","Scaffold14_Amanita_muscaria","Scaffold15_Amanita_muscaria","Scaffold16_Amanita_muscaria","Scaffold17_Amanita_muscaria","Scaffold18_Amanita_muscaria","Scaffold19_Amanita_muscaria","Scaffold20_Amanita_muscaria","Scaffold21_Amanita_muscaria","Scaffold22_Amanita_muscaria","Scaffold23_Amanita_muscaria","Scaffold24_Amanita_muscaria","Scaffold25_Amanita_muscaria","Scaffold26_Amanita_muscaria","Scaffold27_Amanita_muscaria","Scaffold28_Amanita_muscaria")

#### 
# plotting summary statistics across all chromosomes
p2 <- pixy_df2 %>%
#  mutate(chrom_color_group = case_when(as.numeric(chromosome) %% 2 != 0 ~ "even",
 #                                      TRUE ~ "odd" )) %>%
  mutate(chromosome = factor(chromosome, levels = scaffoldnames)) %>%
  filter(statistic %in% c("avg_pi", "avg_dxy", "avg_wc_fst")) %>%
 # ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value, color = chrom_color_group))+
  ggplot(aes(x = (window_pos_1 + window_pos_2)/2, y = value))+
  geom_point(size = 0.5, alpha = 0.5, stroke = 0)+
  facet_grid(statistic ~ chromosome,
             scales = "free_y", switch = "x", space = "free_x",
             labeller = labeller(statistic = pixy_labeller,
                                 value = label_value))+
  xlab("Supercontig")+
  ylab("Statistic Value")+
  scale_color_manual(values = c("grey50", "black"))+
  theme_classic()+
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank(),
        panel.spacing = unit(0.1, "cm"),
        strip.background = element_blank(),
        strip.placement = "outside",
        legend.position ="none")+
  scale_x_continuous(expand = c(0, 0)) +
  scale_y_continuous(expand = c(0, 0), limits = c(0,1))
p2









#### Significance testing for the random subsets which control for sample size
#Summary: in all cases, differences between T and any/all other sites are significant and there are no significant differences between any other sites. 

df <- as.data.frame(read.table("subA_pi.txt", header = TRUE))
modelAOV <- aov(avg_pi ~ pop, data = df)
summary(modelAOV)
             Df Sum Sq  Mean Sq F value   Pr(>F)    
pop           4 0.0204 0.005109   7.334 8.07e-06 ***
Residuals   970 0.6757 0.000697                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
5 observations deleted due to missingness
> TukeyHSD(modelAOV, "pop", ordered= FALSE, conf.level=.95)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = avg_pi ~ pop, data = df)

$pop
             diff          lwr          upr     p adj
H-F  0.0037336835 -0.003571129  0.011038496 0.6298331
S-F  0.0024936772 -0.004811136  0.009798490 0.8840873
T-F -0.0092643330 -0.016569146 -0.001959520 0.0049927
W-F  0.0004212924 -0.006883520  0.007726105 0.9998615
S-H -0.0012400063 -0.008544819  0.006064806 0.9904853
T-H -0.0129980165 -0.020302829 -0.005693204 0.0000133
W-H -0.0033123911 -0.010617204  0.003992422 0.7283253
T-S -0.0117580102 -0.019062823 -0.004453197 0.0001177
W-S -0.0020723847 -0.009377198  0.005232428 0.9377070
W-T  0.0096856254  0.002380813  0.016990438 0.0028197

> df <- as.data.frame(read.table("subB_pi.txt", header = TRUE))
> modelAOV <- aov(avg_pi ~ pop, data = df)
> summary(modelAOV)
             Df Sum Sq  Mean Sq F value   Pr(>F)    
pop           4 0.0179 0.004467    6.38 4.55e-05 ***
Residuals   970 0.6792 0.000700                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
5 observations deleted due to missingness
> TukeyHSD(modelAOV, "pop", ordered= FALSE, conf.level=.95)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = avg_pi ~ pop, data = df)

$pop
             diff          lwr          upr     p adj
H-F  0.0011538597 -0.006169891  0.008477610 0.9928419
S-F  0.0009218186 -0.006401932  0.008245569 0.9969921
T-F -0.0099823919 -0.017306142 -0.002658642 0.0019263
W-F  0.0006300528 -0.006693697  0.007953803 0.9993240
S-H -0.0002320411 -0.007555791  0.007091709 0.9999873
T-H -0.0111362516 -0.018460002 -0.003812501 0.0003401
W-H -0.0005238069 -0.007847557  0.006799943 0.9996745
T-S -0.0109042106 -0.018227961 -0.003580460 0.0004892
W-S -0.0002917658 -0.007615516  0.007031984 0.9999683
W-T  0.0106124447  0.003288694  0.017936195 0.0007646


> df <- as.data.frame(read.table("subC_pi.txt", header = TRUE))
> modelAOV <- aov(avg_pi ~ pop, data = df)
> summary(modelAOV)
             Df Sum Sq  Mean Sq F value   Pr(>F)    
pop           4 0.0134 0.003353   4.918 0.000628 ***
Residuals   970 0.6613 0.000682                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
5 observations deleted due to missingness
> TukeyHSD(modelAOV, "pop", ordered= FALSE, conf.level=.95)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = avg_pi ~ pop, data = df)

$pop
             diff           lwr           upr     p adj
H-F -0.0002661619 -0.0074930596  0.0069607357 0.9999768
S-F -0.0021441644 -0.0093710620  0.0050827333 0.9273397
T-F -0.0101684666 -0.0173953643 -0.0029415689 0.0012090
W-F -0.0025624371 -0.0097893348  0.0046644606 0.8690833
S-H -0.0018780024 -0.0091049001  0.0053488952 0.9542012
T-H -0.0099023047 -0.0171292023 -0.0026754070 0.0017870
W-H -0.0022962752 -0.0095231728  0.0049306225 0.9084331
T-S -0.0080243023 -0.0152511999 -0.0007974046 0.0208355
W-S -0.0004182727 -0.0076451704  0.0068086249 0.9998595
W-T  0.0076060295  0.0003791318  0.0148329272 0.0334079

> df <- as.data.frame(read.table("subD_pi.txt", header = TRUE))
> modelAOV <- aov(avg_pi ~ pop, data = df)
> summary(modelAOV)
             Df Sum Sq  Mean Sq F value   Pr(>F)    
pop           4 0.0223 0.005582   7.872 3.03e-06 ***
Residuals   970 0.6878 0.000709                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
5 observations deleted due to missingness
> TukeyHSD(modelAOV, "pop", ordered= FALSE, conf.level=.95)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = avg_pi ~ pop, data = df)

$pop
             diff          lwr          upr     p adj
H-F  0.0020541436 -0.005315763  0.009424050 0.9414213
S-F  0.0001275903 -0.007242316  0.007497497 0.9999989
T-F -0.0096493840 -0.017019290 -0.002279478 0.0033335
W-F  0.0043671840 -0.003002722  0.011737090 0.4852124
S-H -0.0019265533 -0.009296460  0.005443353 0.9532292
T-H -0.0117035276 -0.019073434 -0.004333621 0.0001531
W-H  0.0023130404 -0.005056866  0.009682947 0.9121255
T-S -0.0097769743 -0.017146881 -0.002407068 0.0028002
W-S  0.0042395937 -0.003130313  0.011609500 0.5157149
W-T  0.0140165680  0.006646662  0.021386474 0.0000024

> df <- as.data.frame(read.table("subE_pi.txt", header = TRUE))
> modelAOV <- aov(avg_pi ~ pop, data = df)
> summary(modelAOV)
             Df Sum Sq  Mean Sq F value   Pr(>F)    
pop           4 0.0243 0.006075   8.739 6.24e-07 ***
Residuals   970 0.6742 0.000695                     
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
5 observations deleted due to missingness
> TukeyHSD(modelAOV, "pop", ordered= FALSE, conf.level=.95)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = avg_pi ~ pop, data = df)

$pop
             diff          lwr          upr     p adj
H-F  0.0011557983 -0.006141256  0.008452852 0.9926944
S-F -0.0009903981 -0.008287452  0.006306656 0.9959693
T-F -0.0099823919 -0.017279446 -0.002685338 0.0018286
W-F  0.0051889119 -0.002108142  0.012485966 0.2950844
S-H -0.0021461964 -0.009443251  0.005150858 0.9294724
T-H -0.0111381902 -0.018435244 -0.003841136 0.0003178
W-H  0.0040331136 -0.003263941  0.011330168 0.5558526
T-S -0.0089919938 -0.016289048 -0.001694940 0.0070365
W-S  0.0061793100 -0.001117744  0.013476364 0.1411460
W-T  0.0151713038  0.007874250  0.022468358 0.0000002

> df <- as.data.frame(read.table("subF_pi.txt", header = TRUE))
> modelAOV <- aov(avg_pi ~ pop, data = df)
> summary(modelAOV)
             Df Sum Sq  Mean Sq F value  Pr(>F)    
pop           4 0.0237 0.005921   8.429 1.1e-06 ***
Residuals   970 0.6814 0.000702                    
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
5 observations deleted due to missingness
> TukeyHSD(modelAOV, "pop", ordered= FALSE, conf.level=.95)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = avg_pi ~ pop, data = df)

$pop
             diff          lwr          upr     p adj
H-F  0.0040444236 -0.003291024  0.011379871 0.5582642
S-F  0.0004806959 -0.006854751  0.007816143 0.9997699
T-F -0.0095882054 -0.016923653 -0.002252758 0.0034073
W-F  0.0036506985 -0.003684749  0.010986146 0.6534213
S-H -0.0035637277 -0.010899175  0.003771720 0.6739613
T-H -0.0136326290 -0.020968076 -0.006297182 0.0000045
W-H -0.0003937251 -0.007729172  0.006941722 0.9998959
T-S -0.0100689013 -0.017404349 -0.002733454 0.0017420
W-S  0.0031700026 -0.004165445  0.010505450 0.7623112
W-T  0.0132389039  0.005903457  0.020574351 0.0000095


## Stats on the main pi comparison
> df <- as.data.frame(read.table("../1to28win5000_pi.txt", header = TRUE))
> modelAOV <- aov(avg_pi ~ pop, data = df)
> summary(modelAOV)
               Df Sum Sq Mean Sq F value Pr(>F)    
pop             4  0.546 0.13640   82.96 <2e-16 ***
Residuals   17545 28.846 0.00164                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
2280 observations deleted due to missingness
> TukeyHSD(modelAOV, "pop", ordered= FALSE, conf.level=.95)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = avg_pi ~ pop, data = df)

$pop
             diff           lwr          upr     p adj
H-F -0.0005061426 -0.0031466080  0.002134323 0.9850672
S-F  0.0013905440 -0.0012499214  0.004031009 0.6038526
T-F -0.0136548291 -0.0162952945 -0.011014364 0.0000000
W-F -0.0001164204 -0.0027568858  0.002524045 0.9999529
S-H  0.0018966865 -0.0007437789  0.004537152 0.2860382
T-H -0.0131486865 -0.0157891519 -0.010508221 0.0000000
W-H  0.0003897222 -0.0022507432  0.003030188 0.9944757
T-S -0.0150453731 -0.0176858385 -0.012404908 0.0000000
W-S -0.0015069643 -0.0041474298  0.001133501 0.5252587
W-T  0.0135384087  0.0108979433  0.016178874 0.0000000


## And the main pi comparison with the "outliers" removed

> df <- as.data.frame(read.table("../nooutliers/1to28win5000_outliers_pi.txt", header = TRUE))
> modelAOV <- aov(avg_pi ~ pop, data = df)
> summary(modelAOV)
               Df Sum Sq Mean Sq F value Pr(>F)    
pop             5   5.23  1.0460     611 <2e-16 ***
Residuals   21054  36.04  0.0017                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
2736 observations deleted due to missingness
> TukeyHSD(modelAOV, "pop", ordered= FALSE, conf.level=.95)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = avg_pi ~ pop, data = df)

$pop
             diff           lwr          upr     p adj
H-F -0.0005061426 -0.0033208708  0.002308586 0.9957159
O-F -0.0429347342 -0.0457494624 -0.040120006 0.0000000
S-F  0.0013905440 -0.0014241842  0.004205272 0.7222706
T-F -0.0082855590 -0.0111002872 -0.005470831 0.0000000
W-F -0.0001164204 -0.0029311486  0.002698308 0.9999968
O-H -0.0424285916 -0.0452433198 -0.039613863 0.0000000
S-H  0.0018966865 -0.0009180417  0.004711415 0.3894912
T-H -0.0077794165 -0.0105941447 -0.004964688 0.0000000
W-H  0.0003897222 -0.0024250060  0.003204450 0.9987706
S-O  0.0443252782  0.0415105500  0.047140006 0.0000000
T-O  0.0346491752  0.0318344470  0.037463903 0.0000000
W-O  0.0428183138  0.0400035856  0.045633042 0.0000000
T-S -0.0096761030 -0.0124908312 -0.006861375 0.0000000
W-S -0.0015069643 -0.0043216925  0.001307764 0.6475238
W-T  0.0081691387  0.0053544105  0.010983867 0.0000000

> df <- df[df$pop != "O", ]
> modelAOV <- aov(avg_pi ~ pop, data = df)
> summary(modelAOV)
               Df Sum Sq Mean Sq F value Pr(>F)    
pop             4  0.209 0.05226   31.42 <2e-16 ***
Residuals   17545 29.178 0.00166                   
---
Signif. codes:  0 ‘***’ 0.001 ‘**’ 0.01 ‘*’ 0.05 ‘.’ 0.1 ‘ ’ 1
2280 observations deleted due to missingness
> TukeyHSD(modelAOV, "pop", ordered= FALSE, conf.level=.95)
  Tukey multiple comparisons of means
    95% family-wise confidence level

Fit: aov(formula = avg_pi ~ pop, data = df)

$pop
             diff           lwr          upr     p adj
H-F -0.0005061426 -0.0031617779  0.002149493 0.9853847
S-F  0.0013905440 -0.0012650913  0.004046179 0.6091966
T-F -0.0082855590 -0.0109411943 -0.005629924 0.0000000
W-F -0.0001164204 -0.0027720557  0.002539215 0.9999540
S-H  0.0018966865 -0.0007589488  0.004552322 0.2917627
T-H -0.0077794165 -0.0104350518 -0.005123781 0.0000000
W-H  0.0003897222 -0.0022659131  0.003045358 0.9945962
T-S -0.0096761030 -0.0123317383 -0.007020468 0.0000000
W-S -0.0015069643 -0.0041625996  0.001148671 0.5310540
W-T  0.0081691387  0.0055135034  0.010824774 0.0000000
