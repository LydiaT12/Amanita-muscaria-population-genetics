library(ggplot2)

FH_F <- read.delim("FandH_F.frq", header=TRUE)
FH_H <- read.delim("FandH_H.frq", header=TRUE)
FHfreq <- data.frame(f = FH_F$FREQ1, h= FH_H$FREQ1)

fh <- ggplot(FHfreq, aes(x=f, y=h)) + 
 geom_bin2d(bins=10)+
 scale_fill_gradient(trans = "log", limits = c(1, 3000))

FS_F <- read.delim("FandS_F.frq", header=TRUE)
FS_S <- read.delim("FandS_S.frq", header=TRUE)
FSfreq <- data.frame(f = FS_F$FREQ1, s= FS_S$FREQ1)

fs <- ggplot(FSfreq, aes(x=f, y=s)) + 
  geom_bin2d(bins=10)+
  scale_fill_gradient(trans = "log", limits = c(1, 3500))

FT_F <- read.delim("FandT_F.frq", header=TRUE)
FT_T <- read.delim("FandT_T.frq", header=TRUE)
FTfreq <- data.frame(f = FT_F$FREQ1, t= FT_T$FREQ1)

ft <- ggplot(FTfreq, aes(x=f, y=t)) + 
  geom_bin2d(bins=10)+
  scale_fill_gradient(trans = "log", limits = c(1, 3500))

FW_F <- read.delim("FandW_F.frq", header=TRUE)
FW_W <- read.delim("FandW_W.frq", header=TRUE)
FWfreq <- data.frame(f = FW_F$FREQ1, w= FW_W$FREQ1)

fw <- ggplot(FWfreq, aes(x=f, y=w)) + 
  geom_bin2d(bins=10)+
  scale_fill_gradient(trans = "log", limits = c(1, 3500))



HS_H <- read.delim("HandS_H.frq", header=TRUE)
HS_S <- read.delim("HandS_S.frq", header=TRUE)
HSfreq <- data.frame(h = HS_H$FREQ1, s= HS_S$FREQ1)

hs <- ggplot(HSfreq, aes(x=h, y=s)) + 
  geom_bin2d(bins=10)+
  scale_fill_gradient(trans = "log", limits = c(1, 3500))

HT_H <- read.delim("HandT_H.frq", header=TRUE)
HT_T <- read.delim("HandT_T.frq", header=TRUE)
HTfreq <- data.frame(h = HT_H$FREQ1, t= HT_T$FREQ1)

ht <- ggplot(HTfreq, aes(x=h, y=t)) + 
  geom_bin2d(bins=10)+
  scale_fill_gradient(trans = "log", limits = c(1, 3500))

HW_H <- read.delim("HandW_H.frq", header=TRUE)
HW_W <- read.delim("HandW_W.frq", header=TRUE)
HWfreq <- data.frame(h = HW_H$FREQ1, w= HW_W$FREQ1)

hw <- ggplot(HWfreq, aes(x=h, y=w)) + 
  geom_bin2d(bins=10)+
  scale_fill_gradient(trans = "log", limits = c(1, 3500))



ST_S <- read.delim("SandT_S.frq", header=TRUE)
ST_T <- read.delim("SandT_T.frq", header=TRUE)
STfreq <- data.frame(s = ST_S$FREQ1, t= ST_T$FREQ1)

st <- ggplot(STfreq, aes(x=s, y=t)) + 
  geom_bin2d(bins=10)+
  scale_fill_gradient(trans = "log", limits = c(1, 3500))


SW_S <- read.delim("SandW_S.frq", header=TRUE)
SW_W <- read.delim("SandW_W.frq", header=TRUE)
SWfreq <- data.frame(s = SW_S$FREQ1, w= SW_W$FREQ1)

sw <- ggplot(SWfreq, aes(x=s, y=w)) + 
  geom_bin2d(bins=10)+
  scale_fill_gradient(trans = "log", limits = c(1, 3500))

TW_T <- read.delim("TandW_T.frq", header=TRUE)
TW_W <- read.delim("TandW_W.frq", header=TRUE)
head(TW_T)

TWfreq <- data.frame(t = TW_T$FREQ1, w= TW_W$FREQ1)
# TWfreq <- TWfreq[(TWfreq$t<0.95 | TWfreq$w<0.97), ]

tw <- ggplot(TWfreq, aes(x=t, y=w)) + 
  geom_bin2d(bins=10)+
  scale_fill_gradient(trans = "log", limits = c(1, 3500))



library(ggpubr)

ggarrange(
  fs, fs, st, fw, hs, ht, hw, st, sw, tw, NULL, NULL, nrow = 3, ncol = 4, labels = c("F-H", "F-S", "F-T", "F-W", "H-S", "H-T", "H-W", "S-T", "S-W", "T-W", "", ""),
   common.legend = TRUE, legend = "bottom")


legend_plot <- ggplot() +
  geom_blank() + 
  theme_void() +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(title = "Legend Title"))

plot_grid / legend_plot

