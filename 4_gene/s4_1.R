############
######
library(linkET)
library(reshape2)
library(ggplot2)
library(patchwork)
library(cowplot)
library(openxlsx)
######
env <- read.table("data/metadata.txt", header=T)[,c(2:18,20,21)]
s_gene_d <- read.table("data/sgene.diversity.tsv", header=T, row.names=1)
n_gene_d <- read.table("data/ngene.diversity.tsv", header=T, row.names=1)
f_gene_d <- read.table("data/fegene.diversity.tsv", header=T, row.names=1)
gene_d <- data.frame(cbind(s_gene_d, n_gene_d, f_gene_d))
gene_d$amoA <- rowSums(gene_d[,c("amoA_A","amoA_B")])
key_gene_d <- gene_d[,c("DFE_0461","MtrC_TIGR03507", 
                    "MtoA","Cyc2_repCluster3","Cyc2_repCluster2","Cyc2_repCluster1", 
                    "nifH","nosZ","norB","nirB","nirS","nirK","narG", 
                    "nxrB","hao","amoA", 
                    "asrB","dsrB","dsrA","aprA","sat", 
                    "soxB","fccB","sqr")]
colnames(key_gene_d) <- c("DFE_0461","mtrC", 
                        "mtoA","cyc2_3","cyc2_2","cyc2_1", 
                        "nifH","nosZ","norB","nirB","nirS","nirK","narG", 
                        "nxrB","hao","amoA", 
                        "asrB","dsrB","dsrA","aprA","sat", 
                        "soxB","fccB","sqr")
write.xlsx(key_gene_d, "table/t4_1.xlsx", rowNames=T)
######
rf_d <- random_forest(key_gene_d, env, seed=1234)
gene_mean_d <- colMeans(key_gene_d)/4
bar_data_d <- cbind(rf_d$explained, gene_mean_d)
bar_data_d <- melt(bar_data_d, measure.vars=c("explained","gene_mean_d"), variable.name="Var", value.name="value")
ggg <- colorRampPalette(c("#00bfc4","white","#f8766d"))(400)
zgs <- c("#d6a01d","#d6a01d","#d6a01d","#d6a01d","#d6a01d","#d6a01d","#d6a01d","#d6a01d", 
         "#bf3553","#bf3553","#bf3553","#bf3553","#bf3553", 
         "#bf3553","#bf3553","#bf3553","#bf3553","#bf3553", 
         "#229453","#229453","#229453","#229453","#229453","#229453")
orderc <- c(rep("gray", 24), rev(zgs))
######
p2_1 <- correlate(key_gene_d, env, method="spearman", adjust=F) %>%
  qcorrplot(extra_mat=list(importance=rf_d$importance, pvalue=rf_d$p), 
            grid_col="gray80", grid_size=0.2) +
  geom_tile(color="gray80", linewidth=0.2, 
            data=function(data){data[data$p<0.05, , drop=F]}) +
  scale_fill_gradientn(colors=ggg, limits=c(-1,1), breaks=seq(-1,1,0.5)) +
  geom_point(aes(size=importance), fill=NA, shape=21, stroke=0.8, 
             data=function(data){data[data$pvalue<0.05, , drop=F]}) +
  theme_void() +
  coord_flip() +
  xlab(NULL) + ylab(NULL) +
  guides(size=guide_legend(title="Importance", order=1), 
         fill=guide_colorbar(title="Spearman's r", order=2, ticks=F, title.vjust=2, 
                             barwidth=1, barheight=5)) +
  theme(legend.position="right", legend.key=element_rect(color="white", fill="white"), 
        axis.text=element_text(color="black"), axis.ticks=element_blank(), 
        axis.text.x=element_text(color=zgs, angle=30, face="italic"), 
        plot.margin=margin(0,0,0,0), aspect.ratio=19/24)
p2_2 <- ggplot(bar_data_d, aes(x=name, y=value, fill=forcats::fct_inorder(interaction(name,Var)))) +
  geom_hline(yintercept=c(20,40,60,80), linetype="dashed", color="gray", linewidth=0.5) +
  geom_col(position="dodge") +
  scale_x_discrete(limits=rev(names(key_gene))) +
  scale_fill_manual(values=orderc) +
  theme_grey() +
  xlab(NULL) +
  labs(y="Explained (%)") +
  scale_y_continuous(limits=c(0,90), breaks=c(0,20,40,60,80), expand=c(0,0), 
                     sec.axis=sec_axis(~.*4, breaks=c(0,80,160,240,320), name="Diversity")) +
  theme(axis.text.x=element_blank(), axis.ticks=element_blank(), 
        axis.text=element_text(color="black"), panel.background=element_blank(), 
        axis.title.y=element_text(size=12, vjust=-5), axis.text.y=element_text(size=10), 
        axis.title.y.right=element_text(size=12, vjust=4), 
        plot.margin=margin(0,0,0,0), legend.position="none")
p2 <- p2_2 / plot_spacer() / p2_1 + plot_layout(heights=c(1,-0.2,4))
p22 <- cowplot::get_legend(p2)
p222 <- p2 + theme(legend.position="none", plot.margin=margin(0,50,0,0))
p2222 <- p222 + inset_element(p22, 1.08, 0.2, 1.1, 0.8)
ggsave("pdf/s4_1.pdf", p2222, width=10, height=8.5)


