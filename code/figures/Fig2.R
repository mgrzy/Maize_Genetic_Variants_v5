library(tidyverse)
library(patchwork)
library(data.table)
library(ggsci)
library(scales)
library(ggrepel)
library(ggformula)

# Set plot basic parameters 
theme_set(theme_classic(base_size = 16))
theme_update(axis.text = element_text(colour = "black"))

my_pal <- pal_d3("category20")(11)
my_pal[5] <- "#7F7F7FFF"
my_pal[8] <- "#9467BDFF"


dat <- fread("data/TableS1.csv")

### Map

### Plot PCA 
# Prepare data for PCA plot
datO <- dat[dat$Group_As_On_Fig2=="Other",]
datA <- dat[!dat$Group_As_On_Fig2=="Other",]

dat2 <- rbind(datO, datA)
dat2$Taxa <- ""
dat2$Taxa[838] <- "B73"

#Plot
g1 <- ggplot() + 
  geom_point(data=datO, aes(PC1, PC2, colour=Group_As_On_Fig2), size=1.5, alpha=0.5) +
  geom_point(data=datA, aes(PC1, PC2, colour=Group_As_On_Fig2), size=1.5) +
  geom_text_repel(data=dat2, aes(PC1, PC2, label=Taxa), box.padding = 0.5, max.overlaps = Inf, size=4) + 
  xlab("PC1 (2.6%)") + ylab("PC2 (1.9%)") + 
  scale_color_manual(values = my_pal, name="Group")
  
g2 <- ggplot() + 
  geom_point(data=datO, aes(PC2, PC3, colour=Group_As_On_Fig2), size=1.5, alpha=0.5) +
  geom_point(data=datA, aes(PC2, PC3, colour=Group_As_On_Fig2), size=1.5) +
  geom_text_repel(data=dat2, aes(PC2, PC3, label=Taxa), box.padding = 0.5, max.overlaps = Inf, size=4) + 
  xlab("PC2 (1.9%)") + ylab("PC3 (1.3%)") + 
  scale_color_manual(values = my_pal, name="Group")

# Combine two PCA plot
gAll <- g1 + g2 & theme(legend.position = "top") 
gAll <- gAll + plot_layout(guides = "collect")

### LD decay plot

LD <- fread("../BigResults/maizesnpV5/LD/LD.csv.gz", data.table = F)
  
gLD <- ggplot() + 
  geom_spline(data=LD[!LD$Group=="WiDiv",], aes(Dist, LD, colour=Group), size=1.5) + 
  geom_spline(data=LD[LD$Group=="WiDiv",], aes(Dist, LD, colour=Group), size=1.5, linetype=2) + 
  scale_x_continuous(labels = paste0(c(50, 100, 150, 200), " Kb"),
                     breaks = 10^3 * c(50, 100, 150, 200)) + 
  xlab("Distance (Kb)") + 
  scale_color_manual(values = pal_d3()(6)[c(1,2,4,3,5,6)]) + 
  theme(legend.position=c(.9, .8)) +
  scale_alpha_identity()

### Nucleotide diversity plot

div <- read.csv("data/Diversity/Diversity.csv")

gDiv <- ggplot(div, aes(Group, PI, fill=Group)) +
  geom_boxplot(width=0.5) +
  scale_fill_manual(values = pal_d3()(6)[c(1,2,4,3,5,6)]) +
  theme(legend.position = "none") +
  ylab(expression(pi)) 

#
g <-  gLD + gDiv

  gg <- gAll / g 
  
  gg <- gg + plot_annotation(tag_levels = list(c("a", "", "b", "c")))
  
  ggsave(plot = gg, "results/plots/Fig2.png", width = 12, height = 10)
  ggsave(plot = gg, "results/plots/Fig2.eps", width = 12, height = 10, device="eps")