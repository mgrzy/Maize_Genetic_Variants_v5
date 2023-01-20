library(tidyverse)
library(patchwork)
library(data.table)
library(ggsci)
library(scales)
library(ggrepel)
library(ggformula)
library(maps)
library(geosphere)

### Set plot basic parameters 
theme_set(theme_classic(base_size = 16))
theme_update(axis.text = element_text(colour = "black"))

my_pal <- pal_d3("category20")(11)
my_pal[5] <- "#7F7F7FFF"
my_pal[8] <- "#9467BDFF"

###
dat <- fread("data/TableS1.csv")

### Map

### Plot PCA 
# Prepare data for PCA plot
datO <- dat[dat$Group_As_On_Fig2=="Other",]
datA <- dat[!dat$Group_As_On_Fig2=="Other",]
datA$Group_As_On_Fig2[datA$Group_As_On_Fig2=="Wild relatives"] <- "Other wild relatives"

dat2 <- rbind(datO, datA)
dat2$Taxa <- ""
dat2$Taxa[838] <- "B73"

#Plot
g1 <- ggplot() + 
  geom_point(data=datO, aes(PC1, PC2, colour=Group_As_On_Fig2), size=1.5, alpha=0.5) +
  geom_point(data=datA, aes(PC1, PC2, colour=Group_As_On_Fig2), size=1.5) +
  geom_text_repel(data=dat2, aes(PC1, PC2, label=Taxa), box.padding = 0.5, max.overlaps = Inf, size=4) + 
  xlab("PC1 (2.6%)") + ylab("PC2 (1.9%)") + 
  scale_color_manual(values = my_pal, name="Group") +
  theme(legend.spacing.y = unit(0.1, 'cm'),legend.position = "top") + 
  guides(fill = guide_legend(byrow = TRUE))

g2 <- ggplot() + 
  geom_point(data=datO, aes(PC2, PC3, colour=Group_As_On_Fig2), size=1.5, alpha=0.5) +
  geom_point(data=datA, aes(PC2, PC3, colour=Group_As_On_Fig2), size=1.5) +
  geom_text_repel(data=dat2, aes(PC2, PC3, label=Taxa), box.padding = 0.5, max.overlaps = Inf, size=4) + 
  xlab("PC2 (1.9%)") + ylab("PC3 (1.3%)") + 
  scale_color_manual(values = my_pal, name="Group") + 
  theme(legend.spacing.y = unit(0.1, 'cm'),legend.position = "top") + 
  guides(fill = guide_legend(byrow = TRUE))

# Combine two PCA plot
gAll <- g1 + g2 & theme(legend.position = "top", legend.spacing.y = unit(0.1, 'cm')) 
gAll <- gAll + plot_layout(guides = "collect")

### LD decay plot

LD <- fread("../BigResults/maizesnpV5/LD/LD.csv.gz", data.table = F)

LD$Group[LD$Group=="Temperate"] <- "Northern temperate"
  
gLD <- ggplot() + 
  geom_spline(data=LD[!LD$Group=="WiDiv",], aes(Dist, LD, colour=Group), size=1.5) + 
  geom_spline(data=LD[LD$Group=="WiDiv",], aes(Dist, LD, colour=Group), size=1.5, linetype=2) + 
  scale_x_continuous(labels = paste0(c(50, 100, 150, 200), " Kb"),
                     breaks = 10^3 * c(50, 100, 150, 200)) + 
  xlab("Distance (Kb)") + 
  scale_color_manual(values = pal_d3()(6)[c(1,2,4,3,5,6)]) + 
  theme(legend.position=c(.8, .8), legend.spacing.y = unit(0.01, 'cm')) +
  scale_alpha_identity()

### Nucleotide diversity plot

div <- read.csv("../BigResults/maizesnpV5/NucleotideDiversity/Summary/ND.csv.gz")
colnames(div)[c(1,5)] <- c("Group", "PI")

div$Group[div$Group=="Bugeater"] <- "WiDiv"

gDiv <- ggplot(div, aes(Group, PI, fill=Group)) +
  geom_boxplot(width=0.5) +
  scale_fill_manual(values = pal_d3()(6)[c(1,2,4,3,5,6)]) +
  theme(legend.position = "none", axis.text.x = element_text(angle = 10)) +
  ylab(expression(pi)) 

### Map
world <- map_data("world")
dat2 <- read.csv("data/meta/all1515.csv")

b <- as.data.frame(table(dat2$CountryOfOrigin))

cn <- data.frame(region=unique(dat2$CountryOfOrigin))
cn <- world[world$region %in% cn$region,]

cn <- cn %>% 
  group_by(region) %>% 
  group_modify(~ data.frame(centroid(cbind(.x$long, .x$lat))))

colnames(b)[1] <- "region"

b <- merge(b, cn, by="region")
colnames(b)[2] <- "n"

b$lon[32] <- c(-100)
b$lat[32] <- c(40)
b$lon[27] <- 22
b$lat[27] <- -30
b$lon[6] <- -110
b$lat[6] <- 55

world <- world[world$lat>-55,]

gm <- ggplot(world, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group), fill="white", colour = "lightgray", size=0.1) + 
  theme(panel.background = element_rect(fill = 'lightblue'), legend.position = "top") +
  geom_point(data=b, aes(x=lon, y=lat, size=n), colour="green4") + 
  xlab("Longitude") + 
  ylab("Latitude")

### Merge all plots
gB <-  gLD + gDiv
gg <- gAll / gB
  
gAll <- gm / gg + plot_layout(heights = c(1,2))
gAll <- gAll + plot_annotation(tag_levels = list(c("a", "b", "", "c", "d")))  

ggsave(plot = gAll, "results/figures/Fig2.png", width = 15, height = 18)
ggsave(plot = gAll, "results/figures/Fig2.eps", width = 15, height = 18, device=cairo_ps)
ggsave(plot = gAll, "results/figures/Fig2.pdf", width = 15, height = 18)
