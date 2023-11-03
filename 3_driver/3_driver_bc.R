############
######
library(vegan)
library(ape)
library(ggplot2)
library(cowplot)
######
fun <- read.table("data/kegg.abundance.txt", header=T, row.names=1)
tax <- read.table("data/species.abundance.tsv", header=T, row.names=1)
metadata <- read.table("data/metadata.txt", header=T)
fun_dist <- vegdist(t(fun), method="bray")
tax_dist <- vegdist(t(tax), method="bray")
fun_pcoa <- pcoa(fun_dist)
tax_pcoa <- pcoa(tax_dist)
df2 <- data.frame(fun_pcoa$vectors)
df3 <- data.frame(tax_pcoa$vectors)
df2$group <- metadata$Location
df3$group <- metadata$Location
df2$lat <- metadata$Lat
df3$lat <- metadata$Lat
ggg <- colorRampPalette(c("#00bfc4","#f8766d"))(18)
x_label_fun <- round(fun_pcoa$values$Rel_corr_eig[1]*100, 2)
x_label_tax <- round(tax_pcoa$values$Rel_corr_eig[1]*100, 2)
y_label_fun <- round(fun_pcoa$values$Rel_corr_eig[2]*100, 2)
y_label_tax <- round(tax_pcoa$values$Rel_corr_eig[2]*100, 2)
######
fun_ano <- anosim(fun_dist, metadata$Location, permutations=999, distance="bray")
tax_ano <- anosim(tax_dist, metadata$Location, permutations=999, distance="bray")
summary(fun_ano)
summary(tax_ano)
fun_test <- adonis2(t(fun) ~ Location, metadata, permutations=999, method="bray")
tax_test <- adonis2(t(tax) ~ Location, metadata, permutations=999, method="bray")
fun_r2 <- round(fun_test$R2[1], 2)
tax_r2 <- round(tax_test$R2[1], 2)
fun_p <- fun_test$`Pr(>F)`[1]
tax_p <- tax_test$`Pr(>F)`[1]
######
p2 <- ggplot(df2) +
  geom_vline(xintercept=0, color="gray", linewidth=0.5) +
  geom_hline(yintercept=0, color="gray", linewidth=0.5) +
  geom_point(aes(x=Axis.1, y=Axis.2, color=lat, shape=group)) +
  theme_bw() +
  scale_color_gradientn(colors=ggg) +
  scale_shape_manual(values=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)) +
  theme(panel.grid=element_blank(), text=element_text(size=15), 
        legend.position="none", axis.text=element_text(color="black")) +
  labs(x=paste0("PCoA1 (", x_label_fun, "%)"), y=paste0("PCoA2 (", y_label_fun, "%)"), 
       title=bquote(paste(R^2, " = ", .(fun_r2), ", P-value = ", .(fun_p))))
p3 <- ggplot(df3) +
  geom_vline(xintercept=0, color="gray", linewidth=0.5) +
  geom_hline(yintercept=0, color="gray", linewidth=0.5) +
  geom_point(aes(x=Axis.1, y=Axis.2, color=lat, shape=group)) +
  theme_bw() +
  scale_color_gradientn(colors=ggg) +
  scale_shape_manual(values=c(0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18)) +
  scale_x_continuous(limits=c(-0.5,0.52), breaks=seq(-0.5,0.5,0.25)) +
  theme(panel.grid=element_blank(), text=element_text(size=15), 
        legend.position="none", axis.text=element_text(color="black")) +
  labs(x=paste0("PCoA1 (", x_label_tax, "%)"), y=paste0("PCoA2 (", y_label_tax, "%)"), 
       title=bquote(paste(R^2, " = ", .(tax_r2), ", P-value = ", .(tax_p))))
p23 <- plot_grid(p2, p3, ncol=1, align="h", labels=c("B","C"))
ggsave("pdf/3_driver_bc.pdf", p23, width=4, height=8)


