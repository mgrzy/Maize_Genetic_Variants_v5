library(rMVP)

args <- commandArgs()
k <- as.numeric(args[6])
print(k)

ph <- read.csv("../data/DTS.csv")

geno <- attach.big.matrix("/scratch/Bugeater25.geno.desc")
map <- data.table::fread("/scratch/Bugeater25.geno.map", data.table = F)
Kin <- attach.big.matrix("/scratch/Bugeater25.kin.desc")
pc <-  attach.big.matrix("/scratch/Bugeater25.pc.desc")[]

trait <- ph[,c(1, k)]

colnames(trait)

imMVP <- MVP(
  phe=trait,
  geno=geno,
  map=map,
  CV.FarmCPU=pc, 
  CV.MLM = pc, 
  K = Kin,
  priority="speed",
  ncpus=16,
  maxLoop=10,
  method.bin="FaST-LMM",
  method=c("MLM", "FarmCPU"), 
  file.output = F, 
  p.threshold = 0.05/nrow(map)
)
imMVP <- cbind(imMVP$map, imMVP$mlm.results, imMVP$farmcpu.results)
data.table::fwrite(imMVP, paste0("../results/", colnames(trait)[2], "_All.csv.gz"))

