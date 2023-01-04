library(data.table)
library(tidyverse)
library(patchwork)
library(scales)
library(jcolors)
library(ggsci)
library(bigmemory)
library(ggformula)

### Plot setting
theme_set(theme_classic(base_size = 16))
theme_update(axis.text = element_text(colour = "black"))
my_pal <- pal_d3()(3)

### Gene models in GFF format
gff <- fread("../../BigData/MaizeV5/Zm-B73-REFERENCE-NAM-5.0_Zm00001eb.1.gff3.gz")
colnames(gff) <- c("chr", "assemly", "type", "start", "end", "nothing", "string", "nothing2", "name")
gff$pos <- apply(gff[,4:5],1, mean)

### Phenotypic data
pheno <- read.csv("data/Bugeater2020.csv")
pheno <- pheno %>%
  select(Taxa, DaysToSilk)
ind <- read.table("../../BigData/WiDivGeno/WGSv5/bugeater/rMVP/Bugeater.geno.ind")
colnames(ind) <- "Taxa"
pheno <- plyr::join(ind, pheno, by="Taxa")

### Genotypic data
geno <- attach.big.matrix("../../BigData/WiDivGeno/WGSv5/bugeater/rMVP/Bugeater.geno.desc")
map <- fread("../../BigData/WiDivGeno/WGSv5/bugeater/rMVP/Bugeater.geno.map")

### GWAS results 
gwas.dat <- fread("../BigResults/maizesnpV5/gwas/DaysToSilk_2020.csv.gz", data.table = F)
gwas.dat <- gwas.dat[,c(1,2,3,8,11)]

##############
### MADS69 ###
##############
mads69cord <- c(161172047, 161199909)
mads69pos <- mean(mads69cord)

mads69 <- gwas.dat %>%
  filter(CHROM==3, POS > mads69cord[1]-500000 &  POS < mads69cord[2]+1000000)
  
### Local LD data
mads69LD <- fread("data/Local_LD/mads69_LD_3-161177471.ld")
mads69LD <- mads69LD[,6:7]
colnames(mads69LD)[1] <- "SNP"
  
mads69 <- merge(mads69, mads69LD, by="SNP")
mads69 <- mads69[order(mads69$R2),]

### Local GWAS results with LD in color
gmads69 <- ggplot() + 
  geom_point(data=mads69, aes(POS, -log10(DaysToSilk.MLM), colour=R2)) + 
  geom_point(data=mads69[mads69$SNP=="chr3_161177471",], aes(POS, -log10(DaysToSilk.MLM)), shape=25, size=5,
             fill="red") + 
  geom_vline(xintercept = mads69pos, linetype=2, colour="red") + 
  annotate("text", x=mads69pos+120000, y=12, label="paste(italic(MADS69))", size=5, parse=T) + 
  geom_hline(yintercept = -log10(0.01/26631280), linetype=2) + 
  ylab(expression(-log[10](p-value))) + 
  xlab("Chromosome 3") + 
  scale_x_continuous(label = label_bytes("kB")) + 
  scale_colour_gradientn(colours=rev(jcolors('rainbow')), name="LD") + 
  theme(legend.position = "top", 
        axis.text.x = element_blank(), 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), 
        legend.text = element_text(angle = 45, hjust = 1, vjust = 1.25))

### Gene track
mads69locBig <- gff %>%
  filter(chr=="chr3", start > mads69cord[1]-500000 &  end < mads69cord[2]+1000000, type=="gene")

mads69locBig$nothing[grep("Zm00001eb143080", mads69locBig$name)] <- "red"
mads69locBig$nothing[!grep("Zm00001eb143080", mads69locBig$name)] <- "black"

gmads69trackBig <- ggplot() + 
  geom_point(data=mads69locBig[mads69locBig$string=="+",], aes(x=pos, y=1.05, colour=nothing), shape="\u25BA", size=6) +
  geom_point(data=mads69locBig[mads69locBig$string=="-",], aes(x=pos, y=0.95, colour=nothing), shape="\u25C4", size=6) +
  ylim(0.9, 1.1) + 
  scale_color_manual(values = c("black", "red")) +
  theme(legend.position = "none", 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank()) + 
  scale_x_continuous(labels = paste0(c("161,000", "161,500", "162,000"), " Kb"),
                     breaks = 10^3 * c(161000, 161500, 162000)) +
xlab("Chromosome 3")

### SNP effect plot as boxplot
pheno$chr3_161177471 <- geno[which(gwas.dat$SNP=="chr3_161177471"),]
pheno$chr3_161177471[pheno$chr3_161177471==0] <- "T"
pheno$chr3_161177471[pheno$chr3_161177471==2] <- "C"
pheno$chr3_161177471[pheno$chr3_161177471==1] <- NA
pheno$chr3_161177471 <- factor(pheno$chr3_161177471, levels = c("T", "C"))

gmads69Box <- ggplot(na.omit(pheno), aes(chr3_161177471, DaysToSilk, fill=chr3_161177471)) +
  geom_boxplot(width=0.5) +
  annotate("text", y=55, x=c(1, 2), label=c("n = 419", "n = 371")) +
  scale_fill_manual(values = my_pal[c(1,3)]) + 
  theme(legend.position = "none") + 
  xlab("chr3:161,177,471")
  
### Proportion of SNP in different groups
mads69SNPprop <- read.csv("data/SNP_proportion/mads69SNPprop.csv")
mads69SNPprop <- pivot_longer(mads69SNPprop, cols = 2:3)
mads69SNPprop$Group <- factor(mads69SNPprop$Group, 
                              levels = c("Z.mexicana","Z.parviglumis", "Tropical", 
                                         "China", "SS", "NSS", "IDT", "Europe"))
mads69SNPprop$name <- factor(mads69SNPprop$name, levels = c("T", "C"))
colnames(mads69SNPprop)[2] <- "SNP"

gmads69prop <- ggplot(mads69SNPprop, aes(Group, value, fill=SNP)) +
  geom_col(colour="black") + 
  theme(legend.position = "top", axis.text.x = element_text(angle = 15, vjust = 1.1, hjust = 1)) +
  scale_fill_manual(values = my_pal[c(1,3)]) + 
  ylab("Freqency")

### MADS gene model
mads69gen <- gff[grep("Zm00001eb143080", gff$name)]
mads69gen <- mads69gen %>%
  filter(type %in% c("gene", "CDS"))

gmads69track <- ggplot() + 
  geom_segment(data=mads69gen[mads69gen$type=="gene",], aes(y=0, yend=0, x=start-100, xend=end)) +
  # arrow = arrow(length = unit(0.5, "cm"), ends="first", type="closed")) +
  geom_segment(data=mads69gen[mads69gen$type=="CDS",], aes(y=0, yend=0, x=start, xend=end), size=10) +
  #ylim(0., 1.1) + 
  scale_color_manual(values = c("black", "red")) +
  theme(legend.position = "none", 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank()) + 
  scale_x_continuous(labels = paste0(c("161,160", "161,180", "161,200"), " Kb"),
                     breaks = 10^3 * c(161160, 161180, 161200), 
                     limits = c(mads69cord[1]-15000,  mads69cord[2]+15000)) + 
  xlab("Chromosome 3")

### Nucleotide diversity in MADS69 locus
divMads69 <- fread("../BigResults/maizesnpV5/NucleotideDiversity/MADS69/pixy_pi.txt")

divMads69 <- divMads69 %>%
  filter(pop %in% c("Temp", "Tropical", "Z.parviglumis") & 
           window_pos_1 > mads69cord[1]-15000 & window_pos_2 < mads69cord[2]+15000)

divMads69$avg_pi[is.na(divMads69$avg_pi)] <- 0
divMads69$pop[divMads69$pop=="Tropical"] <- "Tropical maize"
divMads69$pop[divMads69$pop=="Temp"] <- "Temperate maize"
divMads69$pop <- factor(divMads69$pop, levels = c("Z.parviglumis", "Tropical maize", "Temperate maize"))

gpimads69 <- ggplot(divMads69, aes(window_pos_1, avg_pi, colour=pop, group=pop)) + 
  geom_spline(spar = 0.25) + 
  scale_color_d3(name="Group") + 
  theme(legend.direction="horizontal", legend.position = c(0.4, 1.1),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  ylab(expression(pi))

### Assemble all MADS69 plots
gmads69All <- gmads69 / gmads69trackBig / gpimads69 /gmads69track / (gmads69Box + gmads69prop) + 
  plot_layout(heights = c(4, 0.5, 2, 0.5, 4))

############
### ZCN8 ###
############

### ZCN8 position
zcn8pos <- c(126678863, 126680377)
z8 <- mean(zcn8pos)

### ZNC8 gene model
zcn8loc <- gff[grep("Zm00001eb353250", gff$name)]
zcn8loc <- zcn8loc %>%
  filter(type %in% c("gene", "CDS"))

### ZCN8 locus GWAS results 
ZCN8 <- gwas.dat %>%
  filter(CHROM==8, POS > zcn8pos[1]-500000 &  POS < zcn8pos[2]+500000)

### Local LD for ZCN8
zcn8LD <- fread("data/Local_LD/zcn8_LD_chr8_126660665.ld")
zcn8LD <- zcn8LD[,6:7]
colnames(zcn8LD)[1] <- "SNP"

ZCN8 <- merge(ZCN8, zcn8LD, by="SNP")
ZCN8 <- ZCN8[order(ZCN8$R2),]

### Local GWAS result with LD
gZR <- ggplot() + 
  geom_point(data=ZCN8, aes(POS, -log10(DaysToSilk.MLM), colour=R2)) + 
  geom_point(data=ZCN8[ZCN8$SNP=="chr8_126660665",], aes(POS, -log10(DaysToSilk.MLM)), shape=25, size=5,
             fill="red") + 
  geom_hline(yintercept = -log10(0.01/16634049), linetype=2) + 
  geom_vline(xintercept = z8, linetype=2, colour="red") + 
  annotate("text", x=z8+80000, y=12, label="paste(italic(ZCN8))", size=5, parse=T) + 
  theme(legend.position = "none", 
        axis.title.x = element_blank(), 
        axis.ticks.x = element_blank(), axis.text.x = element_blank(), 
        legend.text = element_text(angle = 45, hjust = 1, vjust = 1.25)) + 
  ylab(expression(-log[10](p-value))) + 
  scale_colour_gradientn(colours=rev(jcolors('rainbow')), name="LD") + 
  theme(legend.position = "top")

### Local gene track
zcn8locBig <- gff %>%
  filter(chr=="chr8", start > zcn8pos[1]-500000 &  end < zcn8pos[2]+500000, type=="gene")

zcn8locBig$nothing <- "black"
zcn8locBig$nothing[grep(c("Zm00001eb353250"), zcn8locBig$name)] <- "red"
zcn8locBig <- rbind(zcn8locBig[-c(13)], zcn8locBig[c(13),])

gzcn8trackBig <- ggplot() + 
  geom_point(data=zcn8locBig[zcn8locBig$string=="+",], aes(x=pos, y=1.05, colour=nothing), shape="\u25BA", size=6) +
  geom_point(data=zcn8locBig[zcn8locBig$string=="-",], aes(x=pos, y=0.95, colour=nothing), shape="\u25C4", size=6) + 
  ylim(0.9, 1.1) + 
  scale_color_manual(values = c("black", "red")) +
  theme(legend.position = "none", 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank()) + 
  scale_x_continuous(labels = paste0(c("126,200", "126,400", "126,600", "126,800", "127,000"), " Kb"),
                     breaks = 10^3 * c(126200, 126400, 126600, 126800, 127000)) + 
  xlab("Chromosome 8")

### Nucleotide diversity in MADS69 locus
divZCN8 <- fread("../BigResults/maizesnpV5/NucleotideDiversity/ZCN8/ZCN8_PI_pi.txt")

divZCN8 <- divZCN8 %>%
  filter(pop %in% c("Temp", "Tropical", "Z.parviglumis") & 
           window_pos_1 > 126654363 & window_pos_2 < 126695370)

divZCN8$avg_pi[is.na(divZCN8$avg_pi)] <- 0
divZCN8$pop[divZCN8$pop=="Tropical"] <- "Tropical maize"
divZCN8$pop[divZCN8$pop=="Temp"] <- "Temperate maize"
divZCN8$pop <- factor(divZCN8$pop, levels = c("Z.parviglumis", "Tropical maize", "Temperate maize"))

gpizcn8 <- ggplot(divZCN8, aes(window_pos_1, avg_pi, colour=pop, group=pop)) + 
  geom_spline(spar = 0.25) + 
  scale_color_d3(name="Group") + 
  theme(legend.direction="horizontal", legend.position = c(0.4, 1.1),
        axis.title.x = element_blank(), 
        axis.text.x = element_blank(), 
        axis.ticks.x = element_blank()) + 
  ylab(expression(pi))

### ZCN8 gene track
gzcn8track <- ggplot() + 
  geom_segment(data=zcn8loc[zcn8loc$type=="gene",], aes(y=0, yend=0, x=start-100, xend=end)) + 
  #arrow = arrow(length = unit(0.5, "cm"), ends="first", type="closed")) +
  geom_segment(data=zcn8loc[zcn8loc$type=="CDS",], aes(y=0, yend=0, x=start, xend=end), size=10) +
  #ylim(0., 1.1) + 
  scale_color_manual(values = c("black", "red")) +
  theme(legend.position = "none", 
        axis.text.y = element_blank(), 
        axis.title.y = element_blank(), 
        axis.ticks.y = element_blank()) + 
  scale_x_continuous(labels = paste0(c("126,660", "126,670", "126,680", "126,690"), " Kb"),
                     breaks = 10^3 * c(126660, 126670, 126680, 126690), limits = c(zcn8pos[1]-25000,  zcn8pos[2]+15000)) + 
  xlab("Chromosome 8")

### SNP effect for ZCN8
pheno$chr8_126660665 <- geno[which(gwas.dat$SNP=="chr8_126660665"),]
pheno$chr8_126660665[pheno$chr8_126660665==0] <- "C"
pheno$chr8_126660665[pheno$chr8_126660665==2] <- "T"
pheno$chr8_126660665[pheno$chr8_126660665==1] <- NA
pheno$chr8_126660665 <- factor(pheno$chr8_126660665, levels = c("C", "T"))

gZCN8box <- ggplot(na.omit(pheno), aes(chr8_126660665, DaysToSilk, fill=chr8_126660665)) +
  geom_boxplot(width=0.5) +
  scale_fill_manual(values = my_pal[c(1,3)]) + 
  annotate("text", y=55, x=c(1, 2), label=c("n = 688", "n = 102")) + 
  theme(legend.position = "non") + 
  xlab("chr8:126,660,665") + 
  ylab("Days to silk")

### Proportion of allele in different groups 
zcn8SNPprop <- read.csv("data/SNP_proportion/zcn8SNPprop.csv")
zcn8SNPprop <- pivot_longer(zcn8SNPprop, cols = 2:3)
zcn8SNPprop$Group <- factor(zcn8SNPprop$Group, 
                            levels = c("Z.mexicana","Z.parviglumis", "Tropical", "China", "SS", "NSS", "IDT", "Europe"))
zcn8SNPprop$name <- factor(zcn8SNPprop$name, levels = c("C", "T"))
colnames(zcn8SNPprop)[2] <- "SNP"

gzcn8SNPprop <- ggplot(zcn8SNPprop, aes(Group, value, fill=SNP)) +
  geom_col(colour="black") + 
  theme(legend.position = "top", axis.text.x = element_text(angle = 15, vjust = 1.1, hjust = 1)) +
  scale_fill_manual(values = my_pal[c(1,3)]) + 
  ylab("Freqency")

###
gZCN8All <-gZR / gzcn8trackBig /gpizcn8 /gzcn8track / (gZCN8box + gzcn8SNPprop) + 
  plot_layout(heights = c(4,0.5,2,0.5, 4))

gAll <- gmads69All | gZCN8All
gAll <- gAll + plot_annotation(tag_levels = list(c("a", "", "", "", "", "", "b")))
ggsave(plot = gAll, "results/figures/Fig4.png", width = 20, height = 15)
ggsave(plot = gAll, "results/figures/Fig4.eps", width = 20, height = 15, device=cairo_ps)
