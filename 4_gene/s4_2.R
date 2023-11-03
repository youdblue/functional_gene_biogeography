############
######
library(psych)
library(ggcorrplot)
library(patchwork)
library(openxlsx)
######
s_gene <- read.table("data/sgene.abundance.txt", header=T, row.names=1, check.names=F)
n_gene <- read.table("data/ngene.abundance.txt", header=T, row.names=1, check.names=F)
f_gene <- read.table("data/fegene.abundance.txt", header=T, row.names=1, check.names=F)
gene <- data.frame(t(rbind(s_gene, n_gene, f_gene)))
gene$amoA <- rowSums(gene[,c("amoA_A","amoA_B")])
key_gene <- gene[,c("DFE_0461","MtrC_TIGR03507", 
                    "MtoA","Cyc2_repCluster3","Cyc2_repCluster2","Cyc2_repCluster1", 
                    "nifH","nosZ","norB","nirB","nirS","nirK","narG", 
                    "nxrB","hao","amoA", 
                    "asrB","dsrB","dsrA","aprA","sat", 
                    "soxB","fccB","sqr")]
colnames(key_gene) <- c("DFE_0461","mtrC", 
                        "mtoA","cyc2_3","cyc2_2","cyc2_1", 
                        "nifH","nosZ","norB","nirB","nirS","nirK","narG", 
                        "nxrB","hao","amoA", 
                        "asrB","dsrB","dsrA","aprA","sat", 
                        "soxB","fccB","sqr")
write.xlsx(round(key_gene, 2), "table/t4_2.xlsx", rowNames=T)
key_gene_f <- key_gene[,1:6]
key_gene_n <- key_gene[,7:16]
key_gene_s <- key_gene[,17:24]
######
gene_cor_fn <- corr.test(key_gene_f, key_gene_n, method="spearman", adjust="BH")
gene_r_fn <- gene_cor_fn$r
gene_p_fn <- gene_cor_fn$p
gene_cor_sn <- corr.test(key_gene_s, key_gene_n, method="spearman", adjust="BH")
gene_r_sn <- gene_cor_sn$r
gene_p_sn <- gene_cor_sn$p
gene_cor_fs <- corr.test(key_gene_f, key_gene_s, method="spearman", adjust="BH")
gene_r_fs <- gene_cor_fs$r
gene_p_fs <- gene_cor_fs$p
gene_cor_nn <- corr.test(key_gene_n, key_gene_n, method="spearman", adjust="BH")
gene_r_nn <- gene_cor_nn$r
gene_p_nn <- gene_cor_nn$p
######
ps2_1 <- ggcorrplot(gene_r_fn, method="square", type="full", 
           outline.col="white", colors=c("#00bfc4","white","#f8766d"), 
           show.legend=T, legend.title="Spearman's r", 
           p.mat=gene_p_fn, sig.level=0.05, insig="pch", pch=4, pch.cex=3, pch.col="black") +
           theme(axis.text.x=element_text(face="italic", color="black", size=10, angle=30), 
                 axis.text.y=element_text(face="italic", color="black", size=10), 
                 panel.grid.major=element_blank()) +
           guides(fill=guide_colorbar(ticks=F, title.vjust=2, barwidth=1.2, barheight=8))
ps2_2 <- ggcorrplot(gene_r_sn, method="square", type="full", 
           outline.col="white", colors=c("#00bfc4","white","#f8766d"), 
           show.legend=T, legend.title="Spearman's r", 
           p.mat=gene_p_sn, sig.level=0.05, insig="pch", pch=4, pch.cex=3, pch.col="black") +
           theme(axis.text.x=element_text(face="italic", color="black", size=10, angle=30), 
                 axis.text.y=element_text(face="italic", color="black", size=10), 
                 panel.grid.major=element_blank()) +
           guides(fill=guide_colorbar(ticks=F, title.vjust=2, barwidth=1.2, barheight=8))
ps2_3 <- ggcorrplot(gene_r_fs, method="square", type="full", 
           outline.col="white", colors=c("#00bfc4","white","#f8766d"), 
           show.legend=T, legend.title="Spearman's r", 
           p.mat=gene_p_fs, sig.level=0.05, insig="pch", pch=4, pch.cex=3, pch.col="black") +
           theme(axis.text.x=element_text(face="italic", color="black", size=10, angle=30), 
                 axis.text.y=element_text(face="italic", color="black", size=10), 
                 panel.grid.major=element_blank()) +
           guides(fill=guide_colorbar(ticks=F, title.vjust=2, barwidth=1.2, barheight=8))
ps2_4 <- ggcorrplot(gene_r_nn, method="square", type="full", 
           outline.col="white", colors=c("#00bfc4","white","#f8766d"), 
           show.legend=T, legend.title="Spearman's r", 
           p.mat=gene_p_nn, sig.level=0.05, insig="pch", pch=4, pch.cex=3, pch.col="black") +
           theme(axis.text.x=element_text(face="italic", color="black", size=10, angle=30), 
                 axis.text.y=element_text(face="italic", color="black", size=10), 
                 panel.grid.major=element_blank()) +
           guides(fill=guide_colorbar(ticks=F, title.vjust=2, barwidth=1.2, barheight=8))
ps2 <- ps2_1 + ps2_2 + ps2_3 + ps2_4 + 
       plot_layout(ncol=2, guides="collect") + plot_annotation(tag_levels="A")
ggsave("pdf/s4_2.pdf", ps2, width=8, height=8)


