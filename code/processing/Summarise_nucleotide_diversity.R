library(data.table)

lf <- list.files("../BigResults/maizesnpV5/NucleotideDiversity/Whole_Genome/", pattern = "pi.txt.gz")
lf <- paste0("../BigResults/maizesnpV5/NucleotideDiversity/Whole_Genome/", lf)

lst <- lapply(lf, fread)
dt <- rbindlist(lst)
dt <- na.omit(dt)

dt <- dt[!dt$pop %in% c("China", "Other", "Rest", "tripsacum", "Temp")]
dt$pop[dt$pop=="Ntemp"] <- "Northern\ntemperate"
dt <- dt[dt$avg_pi<0.1,]

fwrite(dt, "../BigResults/maizesnpV5/NucleotideDiversity/Summary/ND.csv.gz")
