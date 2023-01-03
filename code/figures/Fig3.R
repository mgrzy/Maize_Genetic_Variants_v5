library(data.table)
library(tidyverse)
library(patchwork)
library(scales)
library(jcolors)
library(ggsci)

theme_set(theme_classic(base_size = 16))
theme_update(axis.text.x = element_text(colour = "black"), axis.text.y = element_text(colour = "black"))

###

gwas.dat <- fread("../BigResults/maizesnpV5/gwas/DTS_2020_RNA_v4.csv.gz", data.table = F)

nCHR <- length(unique(gwas.dat$CHROM))
gwas.dat$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(gwas.dat$CHROM))){
  nbp[i] <- max(gwas.dat[gwas.dat$CHROM == i,]$POS)
  gwas.dat[gwas.dat$CHROM == i,"BPcum"] <- gwas.dat[gwas.dat$CHROM == i,"POS"] + s
  s <- s + nbp[i]
}

axis.set <- gwas.dat %>% 
  group_by(CHROM) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

gwas.dat.RNA <- gwas.dat[order(gwas.dat$DaysToSilk.MLM)[1:200000],]

gGWAS1 <- ggplot() + 
  geom_point(data=gwas.dat.RNA, aes(BPcum, -log10(DaysToSilk.MLM), colour=factor(CHROM, levels = c(1:10)))) + 
  geom_hline(yintercept = -log10(0.01/428487), linetype=2) + 
  scale_color_manual(values = rep(pal_d3()(3), 5)) + 
  annotate("text", label="paste(italic(MADS69))", y=12, x=713290558, parse=T, size=5) +
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  ylab(expression(-log[10](p-value))) + 
  xlab("Chromosome") + 
  ylim(1,13)

###

gwas.dat <- fread("../BigResults/maizesnpV5/gwas/DaysToSilk_2020.csv.gz", data.table = F)
gwas.dat <- gwas.dat[,c(1,2,3,8,11)]

nCHR <- length(unique(gwas.dat$CHROM))
gwas.dat$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(gwas.dat$CHROM))){
  nbp[i] <- max(gwas.dat[gwas.dat$CHROM == i,]$POS)
  gwas.dat[gwas.dat$CHROM == i,"BPcum"] <- gwas.dat[gwas.dat$CHROM == i,"POS"] + s
  s <- s + nbp[i]
}

axis.set <- gwas.dat %>% 
  group_by(CHROM) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

gwas.dat2 <- gwas.dat[order(gwas.dat$DaysToSilk.MLM)[1:200000],]

gGWAS2 <- ggplot() + 
  geom_point(data=gwas.dat2, aes(BPcum, -log10(DaysToSilk.MLM), colour=factor(CHROM, levels = c(1:10)))) + 
  geom_hline(yintercept = -log10(0.01/16634049), linetype=2) + 
  scale_color_manual(values = rep(pal_d3()(3), 5)) + 
  annotate("text", label="paste(italic(MADS69))", y=11, x=713290558, parse=T, size=5) +
  annotate("text", label="paste(italic(ZCN8))", y=12, x=1758239801, parse=T, size=5) +
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  ylab(expression(-log[10](p-value))) + 
  xlab("Chromosome") + 
  ylim(3,13)

###

gGWAS <- gGWAS1 / gGWAS2 + plot_annotation(tag_levels = "a")

ggsave(plot = gGWAS, "results/figures/Fig3.png", width = 12, height = 10)
ggsave(plot = gGWAS, "results/figures/Fig3.eps", width = 12, height = 10, device=cairo_ps)
