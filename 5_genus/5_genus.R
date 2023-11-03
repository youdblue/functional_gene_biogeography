############
######
library(dplyr)
library(ggplot2)
library(ggrepel)
library(ggh4x)
######
s_geuns <- read.table("data/s.abun.gene.genus.txt", sep="\t")
n_geuns <- read.table("data/n.abun.gene.genus.txt", sep="\t")
f_geuns <- read.table("data/fe.abun.gene.genus.txt", sep="\t")
genus <- rbind(s_geuns, n_geuns, f_geuns)
colnames(genus) <- c("abundance","gene","tax")
genus$gene[genus$gene=="amoA_A"] <- "amoA"
genus$gene[genus$gene=="amoA_B"] <- "amoA"
genus$gene[genus$gene=="MtoA"] <- "mtoA"
genus$gene[genus$gene=="MtrC_TIGR03507"] <- "mtrC"
genus$gene[genus$gene=="Cyc2_repCluster1"] <- "cyc2_1"
genus$gene[genus$gene=="Cyc2_repCluster2"] <- "cyc2_2"
genus$gene[genus$gene=="Cyc2_repCluster3"] <- "cyc2_3"
genus_sum <- aggregate(genus$abundance, list(genus$tax, genus$gene), sum)
colnames(genus_sum) <- c("tax","gene","sum")
######
gene_genus <- filter(genus_sum, gene %in% c(
               "DFE_0461","mtrC", 
               "mtoA","cyc2_3","cyc2_2","cyc2_1", 
               "nifH","nosZ","norB","nirB","nirS","nirK","narG", 
               "nxrB","hao","amoA", 
               "asrB","dsrB","dsrA","aprA","sat", 
               "soxB","fccB","sqr"))
order <- rev(c("DFE_0461","mtrC", 
               "mtoA","cyc2_3","cyc2_2","cyc2_1", 
               "nifH","nosZ","norB","nirB","nirS","nirK","narG", 
               "nxrB","hao","amoA", 
               "asrB","dsrB","dsrA","aprA","sat", 
               "soxB","fccB","sqr"))
gene_genus <- gene_genus %>% group_by(gene) %>% mutate(abundance=sum/sum(sum))
gene_genus <- filter(gene_genus, !tax %in% "zsglj")
gene_genus <- gene_genus %>% group_by(gene) %>% slice_max(order_by=abundance, n=6)
gene_genus <- gene_genus %>% group_by(gene) %>% mutate(max=max(abundance))
######
zgcts <- c("#303030","#4e1892","#1f3696","#276893","#569597","#c65306","#25386b","#4e5f45","#dbce54","#757570",
           "#5a5c5b","#af5e53","#7ba1a8","#e3efd1","#006e5f","#43454a","#6d7358","#304758","#d7c16b","#aec4b7",
           "#363532","#1b54f2","#c4473d","#c35655","#e4cf8e","#6a6834","#eadcd6","#6493af","#88aea3","#17507d",
           "#4f5355","#b0b7ac","#546b83","#ebe8db","#5d828a","#5c8987","#b6b196","#b49436","#6d614a","#704d4e",
           "#e7693f","#e8853b","#c77a3a","#cad4ba","#0041a5","#85794f","#b7b278","#e7e5d0","#3d6e53","#d54b44",
           "#a9b08f","#973444","#793d56","#e1bda2","#c5bfad","#f5f5dc","#afc8ba","#c1a299","#e9db39","#a71368",
           "#3c5e91","#dea87a","#da9558","#a22076","#ab96c5","#c4c3cb","#c9ae8c","#eacdd1","#e1dbcd","#d5b884",
           "#647370","#674950","#455667","#31678d","#90caaf","#005b5a","#2b5e7d","#5a4c4c","#643441","#2578b5",
           "#fcb1aa","#949c97","#bed2b6","#f2de76","#2ec3e7","#37444b","#ce9335","#625c52","#a03e28","#c43739",
           "#d0853d","#e47542","#4d1919","#b8c8b7","#796f54","#fffafa","#79485a","#d1e3db","#9c6680","#dc143c",
           "#cc3536","#008e59","#dd3b44","#45554a","#3f3f3c","#3e3c3d","#c03f3c","#585a57","#bb1c33","#507883",
           "#89303f","#eb652d","#93a2a9","#dbc7a6","#748a8d","#bca590","#a9987c","#a54358","#c3a6cb","#857e95",
           "#eea5d1","#b8844f")
#colors <- sample(zgcts, size=122)
zgs <- c("#d6a01d","#d6a01d","#d6a01d","#d6a01d", 
         "#d6a01d","#d6a01d","#d6a01d","#d6a01d", 
         "#bf3553","#bf3553","#bf3553","#bf3553","#bf3553", 
         "#bf3553","#bf3553","#bf3553","#bf3553","#bf3553", 
         "#229453","#229453","#229453","#229453","#229453","#229453")
strip <- strip_themed(background_x=elem_list_rect(fill=zgs))
######
cp <- coord_polar(clip="off")
cp$is_free <- function() T
p1 <- ggplot(gene_genus) +
  geom_bar(aes(x=tax, y=abundance, fill=tax), stat="identity") +
  geom_text_repel(aes(x=tax, y=max, label=stringr::str_wrap(tax, 10)), 
                  size=3.5, color="black", fontface="italic", max.time=100, 
                  direction="y", segment.color=NA) +
  cp +
  facet_wrap2(~factor(gene, levels=order), scales="free", ncol=4, strip=strip) +
  scale_fill_manual(values=colors) +
  xlab(NULL) + ylab("Relative abundance") +
  theme_bw() +
  theme(legend.position="none", aspect.ratio=1, 
        panel.spacing=unit(0.1, "cm"), 
        strip.text=element_text(size=15, color="white", face="bold.italic"), 
        axis.title.y=element_text(size=20, vjust=3), 
        axis.text.y=element_text(size=10, color="black"), 
        axis.text.x=element_blank())
ggsave("pdf/5_genus.pdf", p1, width=13, height=18)


