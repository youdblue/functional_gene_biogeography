############
######
library(vegan)
library(ggplot2)
library(ggbreak)
######
dt1 <- read.table("data/catalog.abundance.txt", header=T, row.names=1)
catalog2 <- read.table("data/kegg.abundance.txt", header=T, row.names=1)
catalog3 <- read.table("data/species.abundance.tsv", header=T, row.names=1)
sp2 <- specaccum(t(catalog2), method="random", permutations=100)
sp3 <- specaccum(t(catalog3), method="random", permutations=100)
dt2 <- data.frame(sp2$richness, sp2$sites, sp2$sd)
dt3 <- data.frame(sp3$richness, sp3$sites, sp3$sd)
dt1$lable <- "Nr"
dt2$lable <- "Fun"
dt3$lable <- "Tax"
colnames(dt1) <- c("richness","sites","sd","lable")
colnames(dt2) <- c("richness","sites","sd","lable")
colnames(dt3) <- c("richness","sites","sd","lable")
df2 <- rbind(dt1, dt2, dt3)
df2$lable <- factor(df2$lable, levels=c("Nr","Fun","Tax"))
######
p2 <- ggplot(df2) +
  geom_point(aes(x=sites, y=richness, color=lable), size=0.1, alpha=0.6) +
  geom_errorbar(aes(x=sites, ymin=richness-sd, ymax=richness+sd, color=lable), alpha=0.8, linewidth=0.3) +
  scale_colour_manual(values=c("black","purple","blue")) +
  scale_x_continuous(limits=c(0,100), breaks=seq(0,100,20)) +
  scale_y_continuous(limits=c(7000,90000000), breaks=c(7000,12000)) +
  annotate("rect", xmin=1.5, xmax=100, ymin=50000, ymax=3000000, fill="white") +
  theme_classic() +
  theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"), 
        axis.text.y.right=element_blank(), axis.ticks.y.right=element_blank(), 
        axis.line.y.right=element_blank(), legend.title=element_blank(), 
        text=element_text(size=10), legend.text=element_text(size=8), 
        legend.key.size=unit(12, "pt"), legend.background=element_rect(fill="transparent")) +
  labs(x="Number of samples", y="Cumulative number")
lg <- cowplot::get_legend(p2)
p22 <- p2 + ggbreak::scale_y_break(c(12300,17000), space=0.01, scales=1, ticklabels=c(20000,50000), expand=F) +
  scale_y_break(c(52000,1000000), space=0.01, scales=1, ticklabels=c(4000000,80000000), expand=F) +
  theme(legend.position="none")
p222 <- ggplotify::as.ggplot(print(p22))
p2222 <- p222 + ggimage::geom_subview(x=0.85, y=0.3, subview=lg)
ggsave("pdf/1_catalog_b.pdf", p2222, width=4, height=3)


