library(ggplot2)
library(tidyverse)
library(ggpubr)
setwd("C:/Users/Lydia/OneDrive - University of Otago/Documents/PhD/PhD - Botany/R/Data")

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
