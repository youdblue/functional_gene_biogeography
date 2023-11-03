############
######
library(vegan)
library(reshape2)
library(geosphere)
library(ggplot2)
library(cowplot)
######
fun <- read.table("data/kegg.abundance.txt", header=T, row.names=1)
tax <- read.table("data/species.abundance.tsv", header=T, row.names=1)
fun_mean <- apply(fun, 2, function(x){x/sum(x)})
tax_mean <- apply(tax, 2, function(x){x/sum(x)})
fun$MRA <- apply(fun_mean, 1, function(x){mean(x)})
tax$MRA <- apply(tax_mean, 1, function(x){mean(x)})
fun_abun <- fun[fun$MRA>0.001, c(1:90)]
fun_rare <- fun[fun$MRA<0.00001, c(1:90)]
tax_abun <- tax[tax$MRA>0.001, c(1:90)]
tax_rare <- tax[tax$MRA<0.00001, c(1:90)]
fun_dis_abun <- as.matrix(vegdist(t(fun_abun), method="bray"))
fun_dis_rare <- as.matrix(vegdist(t(fun_rare), method="bray"))
tax_dis_abun <- as.matrix(vegdist(t(tax_abun), method="bray"))
tax_dis_rare <- as.matrix(vegdist(t(tax_rare), method="bray"))
fun_dis_abun[upper.tri(fun_dis_abun)] <- NA
fun_dis_rare[upper.tri(fun_dis_rare)] <- NA
tax_dis_abun[upper.tri(tax_dis_abun)] <- NA
tax_dis_rare[upper.tri(tax_dis_rare)] <- NA
diag(fun_dis_abun) <- NA
diag(fun_dis_rare) <- NA
diag(tax_dis_abun) <- NA
diag(tax_dis_rare) <- NA
fun_melt_abun <- melt(fun_dis_abun)
fun_melt_rare <- melt(fun_dis_rare)
tax_melt_abun <- melt(tax_dis_abun)
tax_melt_rare <- melt(tax_dis_rare)
fun_Dis_abun <- na.omit(fun_melt_abun)
fun_Dis_rare <- na.omit(fun_melt_rare)
tax_Dis_abun <- na.omit(tax_melt_abun)
tax_Dis_rare <- na.omit(tax_melt_rare)
fun_Dis_abun$Sim_abun <- log(1-fun_Dis_abun$value, 10)
fun_Dis_rare$Sim_rare <- log(1-fun_Dis_rare$value, 10)
tax_Dis_abun$Sim_abun <- log(1-tax_Dis_abun$value, 10)
tax_Dis_rare$Sim_rare <- log(1-tax_Dis_rare$value, 10)
######
metadata <- read.table("data/metadata.txt", header=T)
lonlat <- cbind(metadata$Long, metadata$Lat)
Distm <- distm(lonlat, fun=distVincentyEllipsoid)
rownames(Distm) <- metadata$Sample
colnames(Distm) <- metadata$Sample
Distm[upper.tri(Distm)] <- NA
diag(Distm) <- NA
Dist_melt <- melt(Distm)
colnames(Dist_melt) <- c("Sample1","Sample2","Dist")
Dist <- na.omit(Dist_melt)
Dist$Dist_log <- log(Dist$Dist, 10)
######
df5 <- cbind(Dist, fun_Dis_abun, fun_Dis_rare)[,c("Dist_log","Sim_abun","Sim_rare")]
df6 <- cbind(Dist, tax_Dis_abun, tax_Dis_rare)[,c("Dist_log","Sim_abun","Sim_rare")]
colnames(df5) <- c("Dist","abun","rare")
colnames(df6) <- c("Dist","abun","rare")
df55 <- melt(df5, id="Dist", variable.name="Var", value.name="Sim")
df66 <- melt(df6, id="Dist", variable.name="Var", value.name="Sim")
######
lm_fun_abun <- lm(formula=Sim ~ Dist, data=df55, subset=(Var=="abun"))
lm_fun_rare <- lm(formula=Sim ~ Dist, data=df55, subset=(Var=="rare"))
lm_tax_abun <- lm(formula=Sim ~ Dist, data=df66, subset=(Var=="abun"))
lm_tax_rare <- lm(formula=Sim ~ Dist, data=df66, subset=(Var=="rare"))
summary(lm_fun_abun)
summary(lm_fun_rare)
summary(lm_tax_abun)
summary(lm_tax_rare)
######
p5 <- ggplot(df55) +
  geom_point(aes(x=Dist, y=Sim, color=Var), alpha=0.3, size=0.5) +
  geom_smooth(aes(x=Dist, y=Sim, color=Var), method="lm", formula=y ~ x, linewidth=1) +
  scale_colour_manual(values=c("#00bfc4","#f8766d")) +
  annotate("text", x=1, y=-1.6, label="S = -0.019 ***", color="#00bfc4", size=4) +
  annotate("text", x=1, y=-1.8, label="S = -0.089 ***", color="#f8766d", size=4) +
  scale_x_continuous(limits=c(0,6.2), breaks=seq(0,6.2,2)) +
  scale_y_continuous(limits=c(-2,0), breaks=seq(-2,0,0.5)) +
  theme_classic() +
  theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"), 
        text=element_text(size=15), legend.position="none") +
  labs(x=expression(paste(Log[10], " (distance) (m)")), y=expression(paste(Log[10], " (similarity)")))
p6 <- ggplot(df66) +
  geom_point(aes(x=Dist, y=Sim, color=Var), alpha=0.3, size=0.5) +
  geom_smooth(aes(x=Dist, y=Sim, color=Var), method="lm", formula=y ~ x, linewidth=1) +
  scale_colour_manual(values=c("#00bfc4","#f8766d")) +
  annotate("text", x=1, y=-1.6, label="S = -0.122 ***", color="#00bfc4", size=4) +
  annotate("text", x=1, y=-1.8, label="S = -0.089 ***", color="#f8766d", size=4) +
  scale_x_continuous(limits=c(0,6.2), breaks=seq(0,6.2,2)) +
  scale_y_continuous(limits=c(-2,0), breaks=seq(-2,0,0.5)) +
  theme_classic() +
  theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"), 
        text=element_text(size=15), legend.position="none") +
  labs(x=expression(paste(Log[10], " (distance) (m)")), y=expression(paste(Log[10], " (similarity)")))
p56 <- plot_grid(p5, p6, ncol=1, align="h", labels=c("E","F"))
ggsave("pdf/2_pattern_ef.pdf", p56, width=4, height=8)


