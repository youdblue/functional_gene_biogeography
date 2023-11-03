############
######
library(vegan)
library(reshape2)
library(ggplot2)
######
fun <- read.table("data/kegg.abundance.txt", header=T, row.names=1)
tax <- read.table("data/species.abundance.tsv", header=T, row.names=1)
fun_dist <- as.matrix(vegdist(t(fun), method="bray"))
tax_dist <- as.matrix(vegdist(t(tax), method="bray"))
fun_dist[upper.tri(fun_dist)] <- NA
tax_dist[upper.tri(tax_dist)] <- NA
diag(fun_dist) <- NA
diag(tax_dist) <- NA
fun_melt <- melt(fun_dist)
tax_melt <- melt(tax_dist)
colnames(fun_melt) <- c("Sample1","Sample2","Dis_fun")
colnames(tax_melt) <- c("Sample3","Sample4","Dis_tax")
fun_Dis <- na.omit(fun_melt)
tax_Dis <- na.omit(tax_melt)
Dis <- cbind(fun_Dis, tax_Dis)
df5 <- Dis[,c("Dis_fun","Dis_tax")]
colnames(df5) <- c("Fun","Tax")
######
lm <- lm(formula=Fun ~ Tax, data=df5)
p <- signif(coef(summary(lm))["Tax","Pr(>|t|)"], 2)
######
p5 <- ggplot(df5) +
  geom_point(aes(x=Tax, y=Fun), color="gray", size=0.1, alpha=0.8) +
  geom_smooth(aes(x=Tax, y=Fun), method="lm", formula=y ~ x, color="black", linewidth=0.5) +
  scale_x_continuous(limits=c(0,1), breaks=seq(0,1,0.3)) +
  scale_y_continuous(limits=c(0,0.7), breaks=seq(0,0.7,0.2)) +
  theme_classic() +
  theme(legend.position="none", text=element_text(size=10), 
        axis.text=element_text(color="black"), axis.ticks=element_line(color="black")) +
  labs(x="Taxonomic dissimilarity", y="Functional dissimilarity")
ggsave("pdf/1_catalog_e.pdf", p5, width=4, height=3)


