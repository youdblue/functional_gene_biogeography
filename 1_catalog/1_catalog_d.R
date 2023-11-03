############
######
library(ggplot2)
######
phylum <- read.table("data/phylum.abundance.txt", header=T, row.names=1, sep="\t")
phylum_mean <- data.frame(rowMeans(phylum))
phylum_mean <- cbind(rownames(phylum_mean), phylum_mean)
colnames(phylum_mean) <- c("tax","phylum")
rownames(phylum_mean) <- NULL
phylum_order <- phylum_mean[order(phylum_mean$phylum, decreasing=T),]
phylum_order1 <- phylum_order[1:10,]
phylum_order2 <- phylum_order[-1:-10,]
phylum_order3 <- cbind("Others", data.frame(sum(phylum_order2$phylum)))
colnames(phylum_order3) <- c("tax","phylum")
df4 <- rbind(phylum_order1, phylum_order3)
df4$tax <- factor(df4$tax, 
  levels=c("Alphaproteobacteria","Betaproteobacteria","Gammaproteobacteria",
  "Deltaproteobacteria","Pseudomonadota_others","Candidatus Thermoplasmatota",
  "Actinomycetota","Acidobacteriota","Bacillota","Nitrospirota","Others"))
######
p4 <- ggplot(df4) +
  geom_bar(aes(x="Content", y=phylum, fill=tax), stat="identity") +
  coord_polar(theta="y") +
  scale_fill_brewer(palette="Spectral") +
  theme_void() +
  labs(fill="Taxonomy") +
  theme(legend.title=element_text(size=10), legend.text=element_text(size=8), 
        legend.key.size=unit(10, "pt"))
ggsave("pdf/1_catalog_d.pdf", p4, width=4, height=3)


