############
######
library(vegan)
library(dplyr)
library(cowplot)
library(ggplot2)
library(rfPermute)
library(randomForest)
######
fun <- read.table("data/kegg.abundance.txt", header=T, row.names=1)
tax <- read.table("data/species.abundance.tsv", header=T, row.names=1)
fun_shannon <- diversity(t(fun), "shannon")
tax_shannon <- diversity(t(tax), "shannon")
env <- read.table("data/metadata.txt", header=T)[,c(2:18,20,21)]
env$fun <- fun_shannon
env$tax <- tax_shannon
######
set.seed(1234)
model_fun <- rfPermute(env$fun~., data=env[,1:20], ntree=1000)
summary(model_fun)
pred_fun <- predict(model_fun, env[,1:20])
plot(env$fun, pred_fun)
abline(a=0, b=1)
set.seed(1234)
model_tax <- rfPermute(env$tax~., data=env[,c(1:19,21)], ntree=1000)
summary(model_tax)
pred_tax <- predict(model_tax, env[,c(1:19,21)])
plot(env$tax, pred_tax)
abline(a=0, b=1)
######
plotImportance(model_fun, scale=T, sig.only=T)
plotImportance(model_tax, scale=T, sig.only=T)
set.seed(1234)
randomForest(env$fun~., data=env[,1:20], importance=T, ntree=1000)
set.seed(1234)
randomForest(env$tax~., data=env[,c(1:19,21)], importance=T, ntree=1000)
s1_1 <- data.frame(importance(model_fun, scale=T))[,1:2]
s1_2 <- data.frame(importance(model_tax, scale=T))[,1:2]
s1_1 <- cbind(rownames(s1_1), s1_1)
s1_2 <- cbind(rownames(s1_2), s1_2)
colnames(s1_1) <- c("env","imp","p")
colnames(s1_2) <- c("env","imp","p")
s1_1$lable <- "Env"
s1_1[1:2, "lable"] <- "Geo"
s1_2$lable <- "Env"
s1_2[1:2, "lable"] <- "Geo"
s1_1 <- mutate(s1_1, pl=cut(p, breaks=c(0, 0.01, 0.05, Inf), labels=c("**", "*", "n.s.")))
s1_2 <- mutate(s1_2, pl=cut(p, breaks=c(0, 0.01, 0.05, Inf), labels=c("**", "*", "n.s.")))
######
ps1_1 <- ggplot(s1_1[1:10,]) +
  geom_col(aes(x=reorder(env, -imp), y=imp, fill=lable)) +
  geom_text(aes(x=env, y=imp+1, label=pl), size=5) +
  scale_fill_manual(values=c("#00bfc4","#f8766d")) +
  scale_y_continuous(limits=c(0,30), breaks=seq(0,30,5), expand=c(0,0)) +
  labs(x=NULL, y="Increased in MSE (%)", fill=expression(paste(R^2, " = 0.67"))) +
  theme_classic() +
  theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"), 
        legend.position=c(0.8,0.8), legend.key.size=unit(20, "pt"), 
        text=element_text(size=15), axis.text.x=element_text(angle=45, hjust=0.8))
ps1_2 <- ggplot(s1_2[1:14,]) +
  geom_col(aes(x=reorder(env, -imp), y=imp, fill=lable)) +
  geom_text(aes(x=env, y=imp+1, label=pl), size=5) +
  scale_fill_manual(values=c("#00bfc4","#f8766d")) +
  scale_y_continuous(limits=c(0,30), breaks=seq(0,30,5), expand=c(0,0)) +
  labs(x=NULL, y="Increased in MSE (%)", fill=expression(paste(R^2, " = 0.75"))) +
  theme_classic() +
  theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"), 
        legend.position=c(0.8,0.8), legend.key.size=unit(20, "pt"), 
        text=element_text(size=15), axis.text.x=element_text(angle=45, hjust=0.8))
ps1 <- plot_grid(ps1_1, ps1_2, nrow=1, align="h", labels=c("A","B"))
ggsave("pdf/s3_1.pdf", ps1, width=8, height=4)


