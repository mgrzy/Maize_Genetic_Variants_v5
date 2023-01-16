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

### GWAS results 
gwas.dat <- fread("../BigResults/maizesnpV5/gwas/DaysToSilk_2020.csv.gz", data.table = F)
gwas.dat <- gwas.dat[,c(1,2,3,8,11)]

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
ggplot() + 
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

zcn8locBig <- gff %>%
  filter(chr=="chr8", start > zcn8pos[1]-500000 &  end < zcn8pos[2]+500000, type=="gene")

z <- c()
for (i in 1:nrow(zcn8locBig)) {
  q <- ZCN8 %>% filter(POS>zcn8locBig$start[i] & POS<zcn8locBig$end[i])
  z <- rbind.data.frame(z, q)
  print(q)
  rm(q)
}

zz <- ZCN8[!ZCN8$SNP %in% z$SNP, ]

ggplot() + 
  geom_point(data=zz, aes(POS, -log10(DaysToSilk.MLM)), colour="grey50", alpha=0.55) + 
  geom_point(data=z, aes(POS, -log10(DaysToSilk.MLM), colour=R2), size=3) + 
  geom_point(data=ZCN8[ZCN8$SNP=="chr8_126660665",], aes(POS, -log10(DaysToSilk.MLM)), shape=25, size=5,
             fill="grey50") + 
  geom_hline(yintercept = -log10(0.01/16634049), linetype=2) + 
  geom_vline(xintercept = z8, linetype=2, colour="red") + 
  annotate("text", x=z8+50000, y=12, label="paste(italic(ZCN8))", size=5, parse=T) + 
  theme(legend.text = element_text(angle = 45, hjust = 1, vjust = 1.25)) + 
  ylab(expression(-log[10](p-value))) + 
  xlab("Chromosome 8") + 
  scale_colour_gradientn(colours=rev(jcolors('rainbow')), name="LD") + 
  scale_x_continuous(labels = paste0(c("126,200", "126,400", "126,600", "126,800", "127,000"), " Kb"),
                     breaks = 10^3 * c(126200, 126400, 126600, 126800, 127000)) + 
  theme(legend.position = "top")

ggsave(filename = "results/figures/FigS10.png", width = 8, height = 6)
