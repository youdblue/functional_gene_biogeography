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
fun_dis <- as.matrix(vegdist(t(fun), method="bray"))
tax_dis <- as.matrix(vegdist(t(tax), method="bray"))
fun_dis[upper.tri(fun_dis)] <- NA
tax_dis[upper.tri(tax_dis)] <- NA
diag(fun_dis) <- NA
diag(tax_dis) <- NA
fun_melt <- melt(fun_dis)
tax_melt <- melt(tax_dis)
fun_Dis <- na.omit(fun_melt)
tax_Dis <- na.omit(tax_melt)
fun_Dis$Sim <- log(1-fun_Dis$value, 10)
tax_Dis$Sim <- log(1-tax_Dis$value, 10)
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
df3 <- cbind(Dist, fun_Dis)[,c("Dist_log","Sim")]
df4 <- cbind(Dist, tax_Dis)[,c("Dist_log","Sim")]
colnames(df3) <- c("Dist","Sim")
colnames(df4) <- c("Dist","Sim")
df33 <- melt(df3, id="Dist", variable.name="Var", value.name="Sim")
df44 <- melt(df4, id="Dist", variable.name="Var", value.name="Sim")
######
lm_fun <- lm(formula=Sim ~ Dist, data=df33)
lm_fun_1 <- lm(formula=Sim ~ Dist, data=df33, subset=(Dist<3))
lm_fun_10 <- lm(formula=Sim ~ Dist, data=df33, subset=(Dist>4))
lm_tax <- lm(formula=Sim ~ Dist, data=df44)
lm_tax_1 <- lm(formula=Sim ~ Dist, data=df44, subset=(Dist<3))
lm_tax_10 <- lm(formula=Sim ~ Dist, data=df44, subset=(Dist>4))
summary(lm_fun)
summary(lm_fun_1)
summary(lm_fun_10)
summary(lm_tax)
summary(lm_tax_1)
summary(lm_tax_10)
######
p3 <- ggplot(df33) +
  geom_point(aes(x=Dist, y=Sim), color="gray", alpha=0.5, size=0.5) +
  geom_smooth(aes(x=Dist, y=Sim), color="black", method="lm", formula=y ~ x, linewidth=0.8) +
  geom_smooth(data=subset(df33, Dist<3), aes(x=Dist, y=Sim), color="blue", method="lm", formula=y ~ x, linewidth=0.8) +
  geom_smooth(data=subset(df33, Dist>4), aes(x=Dist, y=Sim), color="purple", method="lm", formula=y ~ x, linewidth=0.8) +
  annotate("text", x=1, y=-0.45, label="S = -0.056 ***", color="blue", size=4) +
  annotate("text", x=1, y=-0.5, label="S = -0.030 ***", color="black", size=4) +
  annotate("text", x=1, y=-0.55, label="S = -0.056 ***", color="purple", size=4) +
  scale_x_continuous(limits=c(0,6.2), breaks=seq(0,6.2,2)) +
  scale_y_continuous(limits=c(-0.6,0), breaks=seq(-0.6,0,0.2)) +
  theme_classic() +
  theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"), 
        text=element_text(size=15)) +
  labs(x=expression(paste(Log[10], " (distance) (m)")), y=expression(paste(Log[10], " (similarity)")))
p4 <- ggplot(df44) +
  geom_point(aes(x=Dist, y=Sim), color="gray", alpha=0.5, size=0.5) +
  geom_smooth(aes(x=Dist, y=Sim), color="black", method="lm", formula=y ~ x, linewidth=0.8) +
  geom_smooth(data=subset(df44, Dist<3), aes(x=Dist, y=Sim), color="blue", method="lm", formula=y ~ x, linewidth=0.8) +
  geom_smooth(data=subset(df44, Dist>4), aes(x=Dist, y=Sim), color="purple", method="lm", formula=y ~ x, linewidth=0.8) +
  annotate("text", x=1, y=-1.2, label="S = -0.11 ***", color="blue", size=4) +
  annotate("text", x=1, y=-1.3, label="S = -0.12 ***", color="black", size=4) +
  annotate("text", x=1, y=-1.4, label="S = -0.20 ***", color="purple", size=4) +
  scale_x_continuous(limits=c(0,6.2), breaks=seq(0,6.2,2)) +
  scale_y_continuous(limits=c(-1.5,0), breaks=seq(-1.5,0,0.5)) +
  theme_classic() +
  theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"), 
        text=element_text(size=15)) +
  labs(x=expression(paste(Log[10], " (distance) (m)")), y=expression(paste(Log[10], " (similarity)")))
p34 <- plot_grid(p3, p4, ncol=1, align="h", labels=c("C","D"))
ggsave("pdf/2_pattern_cd.pdf", p34, width=4, height=8)


