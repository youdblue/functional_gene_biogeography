############
######
#install.packages("devtools")
#devtools::install_github("hannet91/ggcor")
library(vegan)
library(dplyr)
library(ggcor)
library(ggplot2)
library(patchwork)
library(cowplot)
######
fun <- read.table("data/kegg.abundance.txt", header=T, row.names=1)
tax <- read.table("data/species.abundance.tsv", header=T, row.names=1)
fun_div <- diversity(t(fun), "shannon")
tax_div <- diversity(t(tax), "shannon")
geo <- read.table("data/metadata.txt", header=T)[,20:21]
env <- read.table("data/metadata.txt", header=T)[,c(2:18,20,21)]
env <- cbind(fun_div, tax_div, env)
env <- decostand(env, "standardize")
df1 <- data.frame(t(fun), t(tax))
######
geo_dis <- as.matrix(vegdist(geo, method="euclidean"))
fun_dis <- as.matrix(vegdist(t(fun), method="bray"))
tax_dis <- as.matrix(vegdist(t(tax), method="bray"))
summary(bioenv(t(fun), env[,3:17], method="spearman", index="bray"))
summary(bioenv(t(tax), env[,3:17], method="spearman", index="bray"))
env_dis_fun <- as.matrix(vegdist(env[,c(4,5,12,17)], method="euclidean"))
env_dis_tax <- as.matrix(vegdist(env[,c(1,4,5)], method="euclidean"))
mantel_fun_geo <- mantel.partial(fun_dis, geo_dis, env_dis_fun, method="spearman")
mantel_fun_env <- mantel.partial(fun_dis, env_dis_fun, geo_dis, method="spearman")
mantel_tax_geo <- mantel.partial(tax_dis, geo_dis, env_dis_tax, method="spearman")
mantel_tax_env <- mantel.partial(tax_dis, env_dis_tax, geo_dis, method="spearman")
######
mantel <- mantel_test(df1, env, spec.select=list(Taxonomy=10262:59402, Function=1:10261))
#mantel <- mantel_test(df1, env, env.ctrl=env[18:19], mantel.fun="mantel.partial", spec.select=list(Taxonomy=10262:59402, Function=1:10261))
mantel <- mutate(mantel, 
  rd=cut(r, breaks=c(-Inf, 0.1, 0.2, Inf), labels=c("< 0.1", "0.1 - 0.2", "> 0.2")), 
  pd=cut(p.value, breaks=c(-Inf, 0.001, 0.01, 0.05, Inf), 
         labels=c("< 0.001", "0.001 - 0.01", "0.01 - 0.05", "> 0.05")))
extra.params <- extra_params(
  spec.label=text_params(size=5, color="black"), 
  spec.point=point_params(size=2, shape=23, fill="red", color="red"), 
  env.point=point_params(size=2, fill="gray", color="gray"), 
  link.params=link_params(env.point.hjust=-1, env.point.vjust=-0.4, 
                          spec.point.hjust=0, spec.point.vjust=0))
ggg <- colorRampPalette(c("#00bfc4","white","#f8766d"))(400)
######
p1 <- quickcor(env, cor.test=T, type="upper", grid.size=0.6, grid.colour="gray") +
  geom_square(color="white") +
  geom_mark(r=NA, size=4, color="white", nudge_y=-0.1) +
  geom_diag_label(nudge_x=-0.55) +
  remove_axis("y") +
  expand_axis(x=-3.5) +
  add_link(mantel, aes(color=pd, size=rd), extra.params=extra.params, curvature=0.1) +
  scale_fill_gradientn(colors=ggg, limits=c(-1,1), breaks=seq(-1,1,0.5)) +
  scale_size_manual(values=c(0.5,1,1.5)) +
  scale_color_manual(values=c("#4F9546","#615E89","#ACAA43","gray")) +
  guides(size=guide_legend(title="Mantel's r", order=2), 
         color=guide_legend(title="Mantel's p", order=3, 
                            override.aes=list(linewidth=1.5), ncol=4), 
         fill=guide_colorbar(title="Spearman's r", order=1, ticks=F)) +
  theme(legend.key=element_blank(), axis.ticks=element_blank(), 
        legend.direction="vertical", legend.position="bottom", 
        legend.box="horizontal", legend.box.just="bottom", 
        legend.background=element_rect(fill="transparent"))
p11 <- cowplot::get_legend(p1)
p111 <- p1 + theme(legend.position="none", plot.margin=margin(0,10,30,5))
p1111 <- p111 + inset_element(p11, -0.1, 0, 1, 0.1, on_top=F)
ggsave("pdf/3_driver_a.pdf", p1111, width=8, height=8)


