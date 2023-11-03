############
######
library(vegan)
library(ggvenn)
library(cowplot)
######
fun <- read.table("data/kegg.abundance.txt", header=T, row.names=1)
fun <- decostand(fun, "hellinger", MARGIN=2)
tax <- read.table("data/species.abundance.tsv", header=T, row.names=1)
tax <- decostand(tax, "hellinger", MARGIN=2)
env <- read.table("data/metadata.txt", header=T)[,c(2:18)]
env <- decostand(env, "standardize")
geo <- read.table("data/metadata.txt", header=T)[,c(20:21)]
pcnm <- data.frame(pcnm(dist(geo))$vectors)
######
decorana(t(fun))
decorana(t(tax))
fun_env_0 <- rda(t(fun) ~ 1, data=env)
fun_env_1 <- rda(t(fun) ~ ., data=env)
fun_geo_0 <- rda(t(fun) ~ 1, data=pcnm)
fun_geo_1 <- rda(t(fun) ~ ., data=pcnm)
tax_env_0 <- rda(t(tax) ~ 1, data=env)
tax_env_1 <- rda(t(tax) ~ ., data=env)
tax_geo_0 <- rda(t(tax) ~ 1, data=pcnm)
tax_geo_1 <- rda(t(tax) ~ ., data=pcnm)
fun_env_sel <- ordiR2step(fun_env_0, scope=formula(fun_env_1), direction="forward")
fun_geo_sel <- ordiR2step(fun_geo_0, scope=formula(fun_geo_1), direction="forward")
tax_env_sel <- ordiR2step(tax_env_0, scope=formula(tax_env_1), direction="forward")
tax_geo_sel <- ordiR2step(tax_geo_0, scope=formula(tax_geo_1), direction="forward")
fun_env_sel$call
fun_geo_sel$call
tax_env_sel$call
tax_geo_sel$call
fun_env_sub <- subset(env, select=names(fun_env_sel$terminfo$ordered))
fun_geo_sub <- subset(pcnm, select=names(fun_geo_sel$terminfo$ordered))
tax_env_sub <- subset(env, select=names(tax_env_sel$terminfo$ordered))
tax_geo_sub <- subset(pcnm, select=names(tax_geo_sel$terminfo$ordered))
######
fun_vpa <- varpart(t(fun), fun_env_sub, fun_geo_sub)
tax_vpa <- varpart(t(tax), tax_env_sub, tax_geo_sub)
anova.cca(rda(t(fun), fun_env_sub, fun_geo_sub))
anova.cca(rda(t(fun), fun_geo_sub, fun_env_sub))
anova.cca(rda(t(tax), tax_env_sub, tax_geo_sub))
anova.cca(rda(t(tax), tax_geo_sub, tax_env_sub))
######
fun_vpa$part$indfract
tax_vpa$part$indfract
df4_1 <-list(Env=c("10 %", "12 %"), Geo=c("31 %", "12 %"))
df4_2 <-list(Env=c("12 %", "16 %"), Geo=c("29 %", "16 %"))
df4_11 <- list_to_data_frame(df4_1)
df4_22 <- list_to_data_frame(df4_2)
######
p4_1 <- ggplot(df4_11) +
  geom_venn(aes(A=Env, B=Geo, label=key), fill_color=c("#00bfc4","#f8766d"), 
            fill_alpha=0.6, stroke_color="white", 
            set_name_size=3, stroke_size=1.5, text_size=3) +
  annotate("text", x=0, y=-1.3, label="Unexplained = 47 %", size=3) +
  annotate("text", x=0, y=1.9, label="Function", size=4.5, color="purple", fontface="bold") +
  annotate("text", x=0, y=2.4, label="Function", size=1, color="white") +
  coord_fixed() +
  theme_void()
p4_2 <- ggplot(df4_22) +
  geom_venn(aes(A=Env, B=Geo, label=key), fill_color=c("#00bfc4","#f8766d"), 
            fill_alpha=0.6, stroke_color="white", 
            set_name_size=3, stroke_size=1.5, text_size=3) +
  annotate("text", x=0, y=-1.3, label="Unexplained = 43 %", size=3) +
  annotate("text", x=0, y=1.9, label="Taxonomy", size=4.5, color="blue", fontface="bold") +
  annotate("text", x=0, y=2.4, label="Taxonomy", size=1, color="white") +
  coord_fixed() +
  theme_void()
p4 <- plot_grid(p4_1, p4_2, nrow=1, align="h", labels=c("D",""), scale=1.3, label_y=0.95)
ggsave("pdf/3_driver_d.pdf", p4, width=4, height=3)


