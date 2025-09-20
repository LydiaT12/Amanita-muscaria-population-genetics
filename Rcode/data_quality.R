MQ <- as.numeric(c(20.00,30.0,40.0,50.0,55.0,60.0))
NSNPmq <- as.numeric(c(1269,1147,862,711,662,394))
kinmq <- as.numeric(c(0.333,0.343,0.376, 0.406, 0.419, 0.474))

plot(MQ, kinmq, type="l", col="blue")
par(new=TRUE)
plot(MQ, NSNPmq, type="l", col="purple")


udepth <- as.numeric(c(10000,5000,4000,3500,3000,2800,2600,2500, 2300))
NSNPdp <- as.numeric(c(1201,1166,985,771,686,597,549,500,442))
kindp <- as.numeric(c(0.330,0.328,0.371,0.405,0.415,0.427,0.433,0.440,0.449))

plot(udepth, kindp, type="l", col="blue")
par(new=TRUE)
plot(udepth, NSNPdp, type="l", col="purple")


filtSNP <- as.numeric(c(834,603,563,544,351,313,286,253))
filtkin <- as.numeric(c(0.376,0.439,0.438,0.442,0.486,0.492,0.494,0.496))

plot(filtSNP, filtkin, type="b", col="blue")
par(new=TRUE)
plot(NSNPdp,kindp, type="b", col="purple")
par(new=TRUE)
# plot(NSNPmq,kinmq, type="b", col="green", xlim xlab="Number SNPs retained", ylab="Average kinship of clones")

contigSNP <- as.numeric(c(280,245))
contigkin <- as.numeric(c(0.494,0.495))

SNPs <- c(filtSNP,NSNPdp,NSNPmq)
kins <- c(filtkin,kindp,kinmq)
type <- c("f","f","f","f","f","f","f","f","d","d","d","d","d","d","d","d","d","q","q","q","q","q","q")

qual_df <- as.data.frame(SNPs) 
qual_df$kins <- kins
qual_df$type <- as.factor(type)

library(ggplot2)
ggplot(qual_df, aes(x = SNPs, y = kins, colour = type, linetype = type)) +
  geom_line(size=1) +
  geom_point(size=3)+
  scale_x_reverse()+
  xlab("Number of SNPs in filtered dataset (1000s)")+
  ylab("Average kinship coefficient of hypothesised clone pairs")+
  labs(colour="Dataset filtered by:")+
  scale_colour_discrete(limits = c("d", "q", "f"), labels = c("Read depth", "Mapping quality", "Both"))+
  guides(linetype=FALSE)


