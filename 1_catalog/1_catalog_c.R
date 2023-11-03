############
######
library(ggplot2)
library(ggbreak)
######
df3 <- data.frame(x=c("EggNOG","KEGG","SCycDB","NCycDB","FeGenie"), 
                  y=c(54282385/81602336,8850519/81602336,659561/81602336,275441/81602336,57105/81602336), 
                  z=c(74,18,0.9,0.42,0.14))
df3$x <- factor(df3$x, levels=c("EggNOG","KEGG","SCycDB","NCycDB","FeGenie"))
######
p3 <- ggplot(df3) +
  geom_bar(aes(x=x, y=y*100), stat="identity", fill=c("#4E4D4D","#4E4D4D","#d6a01d","#bf3553","#229453")) +
  geom_text(aes(x=x, y=z, label=signif(y*100, 2))) +
  theme_classic() +
  theme(axis.text=element_text(color="black"), axis.ticks=element_line(color="black"), 
        axis.text.y.right=element_blank(), axis.ticks.y.right=element_blank(), axis.line.y.right=element_blank()) +
  labs(x="Database", y="Percentage (%)") +
  scale_y_continuous(limits=c(0,80), breaks=c(0,0.5,1)) +
  scale_y_break(c(1,10), space=0.2, scales=0.8, ticklabels=c(10,40,70), expand=F)
ggsave("pdf/1_catalog_c.pdf", p3, width=4, height=3, onefile=F)


