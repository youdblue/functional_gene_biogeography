############
######
library(NST)
library(reshape2)
library(ggplot2)
library(cowplot)
######
fun <- read.table("data/kegg.abundance.txt", header=T, row.names=1, check.names=F)
tax <- read.table("data/species.abundance.tsv", header=T, row.names=1, check.names=F)
metadata <- read.table("data/metadata.txt", header=T, row.names=1)
group <- data.frame(metadata[,"Location"])
rownames(group) <- rownames(metadata)
colnames(group) <- "group"
######
tnst_fun <- tNST(comm=t(fun), group=group, dist.method="bray", abundance.weighted=T, 
                 rand=1000, output.rand=F, nworker=12, LB=F, 
                 null.model="PF")
tnst_tax <- tNST(comm=t(tax), group=group, dist.method="bray", abundance.weighted=T, 
                 rand=1000, output.rand=F, nworker=12, LB=F, 
                 null.model="PF")
######
df5 <- data.frame(tnst_fun$index.pair$ST.ij.bray, tnst_tax$index.pair$ST.ij.bray)
colnames(df5) <- c("Function","Taxonomy")
df55 <- melt(df5, measure.vars=c("Function","Taxonomy"), variable.name="Var", value.name="ST")
######
p5 <- ggplot(df55, aes(x=Var, y=ST)) +
  geom_violin(aes(fill=Var, alpha=0.5), color="white", width=0.8) +
  geom_boxplot(width=0.15, linewidth=0.5) +
  theme_classic() +
  scale_fill_manual(values=c("purple","blue")) +
  theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"), 
        text=element_text(size=15), legend.position="none", 
        axis.text.x=element_text(size=15)) +
  scale_y_continuous(limits=c(0,1), breaks=seq(0,1,0.25)) +
  ylab("Stochasticity Ratio") + xlab(NULL)
p55 <- plot_grid(p5, labels=c("E"))
ggsave("pdf/3_driver_e.pdf", p55, width=4, height=5)


