library(tidyverse)
library(maps)
library(geosphere)
library(mxmaps)
library(scatterpie)
library(ggsci)
library(patchwork)

theme_set(theme_classic(base_size = 16))
theme_update(axis.text = element_text(colour = "black"))

USA <- map_data("state")
USA <- USA[,c(1,2,3,5)]
data(mxstate.map)
mxstate.map <- mxstate.map[,c(1,2,7,8)]
USMex <- rbind(USA, mxstate.map)

######

dat <- read.csv("data/meta/Country.csv")
colnames(dat)[1] <- "Taxa"

dat2 <- read.csv("data/TableS1.csv")
colnames(dat2)[1] <- "Taxa"

dat.Other <- dat %>%
  filter(State_Country_Origin %in% c("Advanta Technology Limited", "Asgrow seed company", "Cargill, Inc., ExPVP",
                                     "Claeys Semences", "DEKALB Plant Genetics", "ExPVP", "Garst Seed Company",
                                     "Holden's Foundation Seeds, Inc.", "Northrup King Company", "Novartis Seeds, Inc.",
                                     "Pioneer Pioneer Hi-Bred International, Inc.", 
                                     "United AgriSeeds, Inc.", "United States",
                                     "USDA"))

dat.Other <- merge(dat.Other, dat2[,c(1,4)], by="Taxa")

b2 <- as.data.frame(table(dat.Other$State_Country_Origin, dat.Other$Group_As_On_Fig2))
colnames(b2)[1] <- "region"
b2$region <- "none"

b2$lon <- -70
b2$lat <- 30

b2 <- b2 %>% group_by(region, Var2, lon, lat) %>% summarise(Freq=sum(Freq))

dat.State <- dat %>% 
  filter(State_Country_Origin %in% c("Alabama", "Arkansas", "Connecticut", "Delaware", "Florida", "Georgia", 
                                     "Illinois", "Indiana", "Iowa", "Kansas", "Kentucky", "Louisiana", "Michigan", 
                                     "Minnesota", "Mississippi", "Missouri", "Nebraska", "New Jersey", "North Carolina",
                                     "North Dakota", "Ohio", "Pennsylvania", "South Carolina", "South Dakota", "Tennessee",
                                     "Texas", "Virginia", "Wisconsin"))
dat.State <- merge(dat.State, dat2[,c(1,4)])

b <- as.data.frame(table(dat.State$State_Country_Origin, dat.State$Group_As_On_Fig2))
colnames(b)[1] <- "region"
b$region <- tolower(b$region)

z <- USA %>% 
  group_by(region) %>% 
  group_modify(~ data.frame(centroid(cbind(.x$long, .x$lat))))

b <- merge(b, z, by="region")

b <- rbind.data.frame(b, b2)

b.n <- b %>% group_by(region, lon, lat) %>% summarise(n=sum(Freq))

b <- b %>% pivot_wider(id_cols = c(region, lon, lat), names_from = Var2, values_from = Freq)
b$r <- abs(rnorm(apply(b[4:8], 1, sum)))

g1 <- ggplot(USA, aes(x = long, y = lat)) +
  geom_polygon(aes(group = group), fill="white", colour = "lightgray", size=0.1) +
  geom_scatterpie(aes(x=lon, y=lat, r=1), data=b,
                  cols=colnames(b)[4:8]) + coord_equal() + 
  theme(panel.background = element_rect(fill = 'lightblue'), legend.position = "top") + 
  scale_fill_d3(name="Group")

g2 <- ggplot() +
  geom_polygon(data=USA, aes(group = group, x = long, y = lat), fill="white", colour = "lightgray", size=0.1) +
  geom_point(data=b.n, aes(x=lon, y=lat, size=n), colour="green4") + 
  theme(panel.background = element_rect(fill = 'lightblue'), legend.position = "top") + 
  scale_size(range=c(1,10), breaks = c(10,50,100,150,200))

gg <- g1 | g2 
gg + plot_annotation(tag_levels = "a")
