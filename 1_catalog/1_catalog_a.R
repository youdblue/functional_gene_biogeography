############
######
library(sf)
library(ggspatial)
library(ggplot2)
library(patchwork)
######
map_china <- read_sf("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json")
map_china$sample <- "Not sampled"
map_china[c(12,14,18,19,20,24), "sample"] <- "Sampled"
map_sampled <- read_sf("https://geo.datav.aliyun.com/areas_v3/bound/100000_full.json")[c(12,14,18,19,20,24),]
metadata1 <- read.table("data/metadata.txt", header=T)[,c("Long","Lat")]
samples <- st_as_sf(metadata1, coords=c("Long","Lat"), crs=st_crs(4326))
metadata2 <- read.table("data/metadata.txt", header=T)[,c("Long","Lat","Location")]
long <- aggregate(metadata2$Long, by=list(type=metadata2$Location), mean)
lat <- aggregate(metadata2$Lat, by=list(type=metadata2$Location), mean)
longlat <- cbind(long, lat)
colnames(longlat) <- c("sites","Long","sites","Lat")
sites <- st_as_sf(longlat, coords=c("Long","Lat"), crs=st_crs(4326))[,3]
######
p1 <- ggplot(map_china) +
  geom_sf(aes(fill=sample), color="black") +
  scale_fill_manual(values=c("white","gray")) +
  theme_bw() +
  coord_sf(label_axes=list(bottom=NA, right=NA, top=NA, left=NA)) +
  theme(panel.grid=element_blank(), legend.position=c(0.3,0.2), legend.title=element_blank(), 
        legend.key.size=unit(8, "pt"), legend.text=element_text(size=8), 
        plot.margin=margin(t=0, r=0, b=0, l=0), aspect.ratio=1, 
        legend.background=element_rect(fill="transparent"))
p11 <- ggplot() +
  geom_sf(data=map_sampled, fill="gray98", color="black") +
  geom_sf(data=samples, color="orange", shape=16, size=1) +
  geom_sf(data=sites, color="black", shape=0, size=2) +
  annotate("text", x=116, y=38, label="No. of samples: 90", size=2.5) +
  annotate("text", x=116, y=39, label="No. of sites: 18", size=2.5) +
  geom_point(aes(x=114, y=36), color="orange", shape=16, size=2) +
  geom_point(aes(x=114, y=37), color="black", shape=0, size=2) +
  annotate("text", x=116, y=36, label="Samples", size=2.5) +
  annotate("text", x=116, y=37, label="Sites", size=2.5) +
  annotation_scale(location="br", style="ticks", width_hint=0.2, 
                   height=unit(0.12, "cm"), pad_y=unit(0.3, "cm")) +
  annotate("text", x=117, y=32.1, label="Anhui\n3 sites\n20 samples", size=2) +
  annotate("text", x=116, y=28.5, label="Jiangxi\n3 sites\n13 samples", size=2) +
  annotate("text", x=113.5, y=23.5, label="Guangdong\n3 sites\n16 samples", size=2) +
  annotate("text", x=112, y=28, label="Hunan\n2 sites\n6 samples", size=2) +
  annotate("text", x=109.2, y=24.3, label="Guangxi\n3 sites\n12 samples", size=2) +
  annotate("text", x=107.5, y=27.5, label="Guizhou\n4 sites\n23 samples", size=2) +
  theme_bw() +
  coord_sf(label_axes=list(bottom="E", right="N", top=NA, left=NA)) +
  scale_x_continuous(limits=c(103,120), breaks=seq(105,120,5)) +
  scale_y_continuous(limits=c(20,39), breaks=seq(20,39,5)) +
  theme(panel.grid=element_blank(), axis.text=element_text(color="black"), 
        axis.title=element_blank(), text=element_text(size=10))
p111 <- p11 + inset_element(p1, -0.0863, 0.6, 0.6, 1)
ggsave("pdf/1_catalog_a.pdf", p111, width=4, height=6)


