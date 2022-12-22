library(data.table)
library(tidyverse)
library(ggsci)
library(patchwork)
library(scales)

theme_set(theme_classic(base_size = 16))
theme_update(axis.text = element_text(colour = "black"))

### SNP denisty across genome
snps100kb <- fread("../BigResults/maizesnpV5/diversity/SNPsPer100kb005.txt", data.table = F)
colnames(snps100kb) <- c("CHROM", "P1", "P2", "N")

snps100kb$POS <- apply(snps100kb[,2:3], 1, mean)

snps100kb$CHROM <- as.numeric(gsub("chr","", snps100kb$CHROM))

CHR <- length(unique(snps100kb$CHROM))
snps100kb$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(snps100kb$CHROM))){
  nbp[i] <- max(snps100kb[snps100kb$CHROM == i,]$POS)
  snps100kb[snps100kb$CHROM == i,"BPcum"] <- snps100kb[snps100kb$CHROM == i,"POS"] + s
  s <- s + nbp[i]
}

axis.set <- snps100kb %>% 
  group_by(CHROM) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

addKb <- function(x) {format(paste0(x, "/kb"))}

g1 <- ggplot() + 
  geom_point(data=snps100kb, aes(BPcum, N, colour=factor(CHROM, levels = c(1:10))), size=1) + 
  scale_x_continuous(label = axis.set$CHROM, breaks = axis.set$center) + 
  theme(legend.position = "none") + 
  scale_color_manual(values = rep(pal_d3()(3), 5)) + 
  #scale_y_continuous(labels = addKb) + 
  ylab("SNPs") + 
  xlab("Chromosome")

### SNP MAF
snpMAF <- fread("../BigResults/maizesnpV5/MAF/stat_1_stats_snps.txt.gz")

g2 <- ggplot(snpMAF, aes(x=V5)) + 
  geom_histogram() + 
  scale_y_continuous(labels = label_number(suffix = " M", scale = 1e-6)) + 
  xlab("Minor Allele Frequency") + 
  ylab("Frequency")

### Indels denisty across genome
indels100kb <- fread("../BigResults/maizesnpV5/diversity/IndelssPer100kb005.txt", data.table = F)
colnames(indels100kb) <- c("CHROM", "P1", "P2", "N")

indels100kb$POS <- apply(indels100kb[,2:3], 1, mean)

indels100kb$CHROM <- as.numeric(gsub("chr","", indels100kb$CHROM))

CHR <- length(unique(indels100kb$CHROM))
indels100kb$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(indels100kb$CHROM))){
  nbp[i] <- max(indels100kb[indels100kb$CHROM == i,]$POS)
  indels100kb[indels100kb$CHROM == i,"BPcum"] <- indels100kb[indels100kb$CHROM == i,"POS"] + s
  s <- s + nbp[i]
}

axis.set2 <- indels100kb %>% 
  group_by(CHROM) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

g3 <- ggplot() + 
  geom_point(data=indels100kb, aes(BPcum, N, colour=factor(CHROM, levels = c(1:10))), size=1) + 
  scale_x_continuous(label = axis.set2$CHROM, breaks = axis.set2$center) + 
  theme(legend.position = "none") + 
  scale_color_manual(values = rep(pal_d3()(3), 5)) + 
  #scale_y_continuous(labels = addKb) + 
  ylab("InDels") + 
  xlab("Chromosome")

### Inels MAF
indelMAF <- fread("../BigResults/maizesnpV5/MAF/stat_1_stats_indels.txt.gz")

g4 <- ggplot(indelMAF, aes(x=V5)) + 
  geom_histogram() + 
  scale_y_continuous(labels = label_number(suffix = " M", scale = 1e-6)) + 
  xlab("Minor Allele Frequency") +
  ylab("Frequency")

### LD score
res <- read.csv("../BigResults/maizesnpV5/LD/LDscore1MB.csv")
res <- na.omit(res)

res$CHROM <- as.numeric(gsub("chr","", res$CHROM))

CHR <- length(unique(res$CHROM))
res$BPcum <- NA
s <- 0
nbp <- c()

for (i in sort(unique(res$CHROM))){
  nbp[i] <- max(res[res$CHROM == i,]$POS)
  res[res$CHROM == i,"BPcum"] <- res[res$CHROM == i,"POS"] + s
  s <- s + nbp[i]
}

axis.set3 <- res %>% 
  group_by(CHROM) %>% 
  summarize(center = (max(BPcum) + min(BPcum)) / 2, maxBP=max(BPcum))

cent <- read.csv("data/meta/Centromers.csv", header = T)

z <- c(0)
zz <- c()
for (i in 1:9){
  z <- sum(z, nbp[i])
  zz <- c(zz, z)
}

for (i in 2:10){
  cent$Start[i] <- cent$Start[i] + zz[i-1]
  cent$End[i] <- cent$End[i] + zz[i-1]
}

g5 <- ggplot() + 
  geom_point(data=res, aes(BPcum, LD, colour=factor(CHROM, levels = c(1:10))), size=1) + 
  geom_point(data=cent, aes(x=Start, y=0),size=6, colour="black", shape=17) + 
  scale_x_continuous(label = axis.set3$CHROM, breaks = axis.set3$center) + 
  # geom_segment(aes(x=404105482, xend=405905334, y=0, yend=0.4)) + 
  theme(legend.position = "none") + 
  scale_color_manual(values = rep(pal_d3()(3), 5)) + 
  ylab("Linkage Disequilibrium") + 
  xlab("Chromosome")

### Variants in genomic regions

region <- read.csv("data/SNPregion.csv")
region$Region[3] <- "CDS"
region$Percent[3] <- 0.4

region$Region<- factor(region$Region, levels = c("Upstream", "5_UTR", "CDS", "Intron", "3_UTR", "Downstream", "Intragenic"))

g6 <- ggplot(region, aes(Region, Percent, fill=Region)) +
  geom_col() + 
  scale_fill_d3(label=c("Downstream", "5'UTR", "CDS", "Intron", "3'UTR", "Upstream", "Intragenic")) + 
  theme(legend.position=c(0.35, 0.72), 
        axis.text.x  = element_blank(), 
        axis.ticks.x = element_blank()) + 
  xlab("") + 
  ylab("Proportion of all variants") + 
  guides(fill=guide_legend(ncol=2))


gSNP <- g1 + g2 + plot_layout(widths = c(4,1))
gIndel <- g3 + g4 + plot_layout(widths = c(4,1))
gB <- g5 + g6 + plot_layout(widths = c(4,1))
gSNP / gIndel / gB + plot_annotation(tag_levels = "a")

###
ggsave("results/figures/Fig1.eps", width = 20, height = 10, device="eps")
ggsave("results/figures/Fig1.png", width = 20, height = 10, device="png")
