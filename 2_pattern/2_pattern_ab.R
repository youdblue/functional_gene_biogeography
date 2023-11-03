############
######
library(vegan)
library(reshape2)
library(ggplot2)
library(cowplot)
######
fun <- read.table("data/kegg.abundance.txt", header=T, row.names=1)
fun_s <- read.table("data/sgene.abundance.txt", header=T, row.names=1)
fun_n <- read.table("data/ngene.abundance.txt", header=T, row.names=1)
fun_f <- read.table("data/fegene.abundance.txt", header=T, row.names=1)
tax <- read.table("data/species.abundance.tsv", header=T, row.names=1)
tax_s <- read.table("data/sspecies.abundance.tsv", header=T, row.names=1)
tax_n <- read.table("data/nspecies.abundance.tsv", header=T, row.names=1)
tax_f <- read.table("data/fespecies.abundance.tsv", header=T, row.names=1)
fun_shannon <- diversity(t(fun), "shannon")
fun_shannon_s <- diversity(t(fun_s), "shannon")
fun_shannon_n <- diversity(t(fun_n), "shannon")
fun_shannon_f <- diversity(t(fun_f), "shannon")
tax_shannon <- diversity(t(tax), "shannon")
tax_shannon_s <- diversity(t(tax_s), "shannon")
tax_shannon_n <- diversity(t(tax_n), "shannon")
tax_shannon_f <- diversity(t(tax_f), "shannon")
metadata <- read.table("data/metadata.txt", header=T)
lat <- data.frame(metadata$Lat)
df_fun <- cbind(fun_shannon, fun_shannon_s, fun_shannon_n, fun_shannon_f, lat)
df_tax <- cbind(tax_shannon, tax_shannon_s, tax_shannon_n, tax_shannon_f, lat)
colnames(df_fun) <- c("All","S","N","Fe","Lat")
colnames(df_tax) <- c("All","S","N","Fe","Lat")
df1 <- melt(df_fun, id="Lat", variable.name="Var", value.name="Diversity")
df2 <- melt(df_tax, id="Lat", variable.name="Var", value.name="Diversity")
######
lm_fun <- lm(formula=Diversity ~ Lat, data=df1, subset=(Var=="All"))
lm_fun_s <- lm(formula=Diversity ~ Lat, data=df1, subset=(Var=="S"))
lm_fun_n <- lm(formula=Diversity ~ Lat, data=df1, subset=(Var=="N"))
lm_fun_f <- lm(formula=Diversity ~ Lat, data=df1, subset=(Var=="Fe"))
lm_tax <- lm(formula=Diversity ~ Lat, data=df2, subset=(Var=="All"))
lm_tax_s <- lm(formula=Diversity ~ Lat, data=df2, subset=(Var=="S"))
lm_tax_n <- lm(formula=Diversity ~ Lat, data=df2, subset=(Var=="N"))
lm_tax_f <- lm(formula=Diversity ~ Lat, data=df2, subset=(Var=="Fe"))
summary(lm_fun)
summary(lm_fun_s)
summary(lm_fun_n)
summary(lm_fun_f)
summary(lm_tax)
summary(lm_tax_s)
summary(lm_tax_n)
summary(lm_tax_f)
######
p1 <- ggplot(df1) +
  geom_point(aes(x=Lat, y=Diversity, color=Var), size=1, alpha=0.8) +
  geom_smooth(aes(x=Lat, y=Diversity, color=Var), method="lm", formula=y ~ x, linewidth=1) +
  scale_colour_manual(values=c("#4E4D4D","#d6a01d","#bf3553","#229453")) +
  scale_x_continuous(limits=c(22,34), breaks=seq(22,34,4)) +
  scale_y_continuous(limits=c(1,8), breaks=seq(1,8,2)) +
  theme_classic() +
  theme(legend.position=c(0.9,0.2), legend.title=element_blank(), text=element_text(size=15), 
        axis.text=element_text(color="black"), axis.ticks=element_line(color="black")) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  labs(x="Latitude", y="Functional diversity")
p2 <- ggplot(df2) +
  geom_point(aes(x=Lat, y=Diversity, color=Var), size=1, alpha=0.8) +
  geom_smooth(aes(x=Lat, y=Diversity, color=Var), method="lm", formula=y ~ x, linewidth=1) +
  scale_colour_manual(values=c("#4E4D4D","#d6a01d","#bf3553","#229453")) +
  scale_x_continuous(limits=c(22,34), breaks=seq(22,34,4)) +
  scale_y_continuous(limits=c(1,7), breaks=seq(1,7,2)) +
  theme_classic() +
  theme(legend.position=c(0.9,0.2), legend.title=element_blank(), text=element_text(size=15), 
        axis.text=element_text(color="black"), axis.ticks=element_line(color="black")) +
  guides(color=guide_legend(override.aes=list(fill=NA))) +
  labs(x="Latitude", y="Taxonomic diversity")
p12 <- plot_grid(p1, p2, ncol=1, align="h", labels=c("A","B"))
ggsave("pdf/2_pattern_ab.pdf", p12, width=4, height=8)


