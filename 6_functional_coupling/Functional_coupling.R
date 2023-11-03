############
######
library(ggsci)
library(eulerr)
library(ggplot2)
library(ggimage)
library(ggplotify)
library(ComplexUpset)
######
bins <- read.table("data/raw_data.txt", header=T, row.names=1, sep="\t")
colnames(bins)[1:7] <- c("NAR","NIR","NOR","NOS","SO","FO","DN")
genes <- rev(colnames(bins)[1:6])
######
bins$Phylum[bins$Phylum=="Nitrospirota_A"] <- "Nitrospirota"
bins$Phylum[bins$Phylum=="Firmicutes_A"] <- "Firmicutes"
bins$Phylum[bins$Phylum=="Firmicutes_B"] <- "Firmicutes"
bins$Phylum[bins$Phylum=="Firmicutes_D"] <- "Firmicutes"
bins$Phylum[bins$Phylum=="Firmicutes_E"] <- "Firmicutes"
bins$Phylum[bins$Phylum=="Firmicutes_G"] <- "Firmicutes"
bins$Phylum[bins$Phylum=="Desulfobacterota_B"] <- "Desulfobacterota"
bins$Phylum[bins$Phylum=="Desulfobacterota_C"] <- "Desulfobacterota"
bins$Phylum[bins$Phylum=="Desulfobacterota_D"] <- "Desulfobacterota"
bins$Phylum[bins$Phylum=="Desulfobacterota_E"] <- "Desulfobacterota"
bins$Phylum[bins$Phylum=="Desulfobacterota_F"] <- "Desulfobacterota"
bins$Phylum[bins$Phylum=="Cyanobacteria"] <- "Others"
bins$Phylum[bins$Phylum=="Armatimonadota"] <- "Others"
bins$Phylum[bins$Phylum=="Thermodesulfobiota"] <- "Others"
bins$Phylum[bins$Phylum=="Patescibacteria"] <- "Others"
bins$Phylum[bins$Phylum=="Nanoarchaeota"] <- "Others"
bins$Phylum[bins$Phylum=="Micrarchaeota"] <- "Others"
bins$Phylum[bins$Phylum=="Nitrospinota"] <- "Others"
bins$Phylum[bins$Phylum=="Bdellovibrionota"] <- "Others"
bins$Phylum[bins$Phylum=="Campylobacterota"] <- "Others"
bins$Phylum[bins$Phylum=="Chlamydiota"] <- "Others"
bins$Phylum[bins$Phylum=="Dependentiae"] <- "Others"
bins$Phylum[bins$Phylum=="FCPU426"] <- "Others"
bins$Phylum[bins$Phylum=="Halobacteriota"] <- "Others"
bins$Phylum[bins$Phylum=="Omnitrophota"] <- "Others"
bins$Phylum[bins$Phylum=="Ratteibacteria"] <- "Others"
bins$Phylum[bins$Phylum=="Thermotogota"] <- "Others"
bins$Phylum[bins$Phylum=="UBA10199"] <- "Others"
bins$Phylum[bins$Phylum=="SpSt-1190"] <- "Others"
######
NAR <- sum(bins[bins$NAR==1 & bins$NIR==0 & bins$NOR==0 & bins$NOS==0,][,"coverm"])
NIR <- sum(bins[bins$NAR==0 & bins$NIR==1 & bins$NOR==0 & bins$NOS==0,][,"coverm"])
NOR <- sum(bins[bins$NAR==0 & bins$NIR==0 & bins$NOR==1 & bins$NOS==0,][,"coverm"])
NOS <- sum(bins[bins$NAR==0 & bins$NIR==0 & bins$NOR==0 & bins$NOS==1,][,"coverm"])
NAR_NIR <- sum(bins[bins$NAR==1 & bins$NIR==1 & bins$NOR==0 & bins$NOS==0,][,"coverm"])
NAR_NOR <- sum(bins[bins$NAR==1 & bins$NIR==0 & bins$NOR==1 & bins$NOS==0,][,"coverm"])
NAR_NOS <- sum(bins[bins$NAR==1 & bins$NIR==0 & bins$NOR==0 & bins$NOS==1,][,"coverm"])
NIR_NOR <- sum(bins[bins$NAR==0 & bins$NIR==1 & bins$NOR==1 & bins$NOS==0,][,"coverm"])
NIR_NOS <- sum(bins[bins$NAR==0 & bins$NIR==1 & bins$NOR==0 & bins$NOS==1,][,"coverm"])
NOR_NOS <- sum(bins[bins$NAR==0 & bins$NIR==0 & bins$NOR==1 & bins$NOS==1,][,"coverm"])
NAR_NIR_NOR <- sum(bins[bins$NAR==1 & bins$NIR==1 & bins$NOR==1 & bins$NOS==0,][,"coverm"])
NAR_NIR_NOS <- sum(bins[bins$NAR==1 & bins$NIR==1 & bins$NOR==0 & bins$NOS==1,][,"coverm"])
NAR_NOR_NOS <- sum(bins[bins$NAR==1 & bins$NIR==0 & bins$NOR==1 & bins$NOS==1,][,"coverm"])
NIR_NOR_NOS <- sum(bins[bins$NAR==0 & bins$NIR==1 & bins$NOR==1 & bins$NOS==1,][,"coverm"])
NAR_NIR_NOR_NOS <- sum(bins[bins$NAR==1 & bins$NIR==1 
                            & bins$NOR==1 & bins$NOS==1,][,"coverm"])
DN <- sum(bins[bins$DN==1 & bins$SO==0 & bins$FO==0,][,"coverm"])
SO <- sum(bins[bins$DN==0 & bins$SO==1 & bins$FO==0,][,"coverm"])
FO <- sum(bins[bins$DN==0 & bins$SO==0 & bins$FO==1,][,"coverm"])
DN_SO <- sum(bins[bins$DN==1 & bins$SO==1 & bins$FO==0,][,"coverm"])
DN_FO <- sum(bins[bins$DN==1 & bins$SO==0 & bins$FO==1,][,"coverm"])
SO_FO <- sum(bins[bins$DN==0 & bins$SO==1 & bins$FO==1,][,"coverm"])
DN_SO_FO <- sum(bins[bins$DN==1 & bins$SO==1 & bins$FO==1,][,"coverm"])
######
sum1 <- sum(NAR, NIR, NOR, NOS, NAR_NIR, NAR_NOR, NAR_NOS, NIR_NOR, NIR_NOS, NOR_NOS, 
            NAR_NIR_NOR, NAR_NIR_NOS, NAR_NOR_NOS, NIR_NOR_NOS, NAR_NIR_NOR_NOS)*0.01
sum2 <- sum(DN, SO, FO, DN_SO, DN_FO, SO_FO, DN_SO_FO)*0.01
venn1 <- round(c("NAR"=NAR/sum1, "NIR"=NIR/sum1, "NOR"=NOR/sum1, "NOS"=NOS/sum1, 
           "NAR&NIR"=NAR_NIR/sum1, "NAR&NOR"=NAR_NOR/sum1, "NAR&NOS"=NAR_NOS/sum1, 
           "NIR&NOR"=NIR_NOR/sum1, "NIR&NOS"=NIR_NOS/sum1, "NOR&NOS"=NOR_NOS/sum1, 
           "NAR&NIR&NOR"=NAR_NIR_NOR/sum1, "NAR&NIR&NOS"=NAR_NIR_NOS/sum1, 
           "NAR&NOR&NOS"=NAR_NOR_NOS/sum1, "NIR&NOR&NOS"=NIR_NOR_NOS/sum1, 
           "NAR&NIR&NOR&NOS"=NAR_NIR_NOR_NOS/sum1), 1)
venn2 <- round(c("DN"=DN/sum2, "SO"=SO/sum2, "FO"=FO/sum2, "DN&SO"=DN_SO/sum2, 
           "DN&FO"=DN_FO/sum2, "SO&FO"=SO_FO/sum2, "DN&SO&FO"=DN_SO_FO/sum2), 1)
######
p1 <- upset(bins, genes, name="Gene sets", sort_sets=F, set_sizes=F, guides="over", 
           sort_intersections_by=c("degree"), sort_intersections="ascending", 
           n_intersections=200, min_size=1, max_size=3000, height_ratio=0.5, 
           mode="inclusive_intersection", 
       stripes=upset_stripes(colors=c("white")), 
       matrix=(intersection_matrix(
              geom=geom_point(shape="circle filled", size=3)) +
              scale_color_manual(
              values=c("NAR"="#d6a01d", "NIR"="#d6a01d", "NOR"="#d6a01d", 
                       "NOS"="#d6a01d", "SO"="#bf3553", "FO"="#229453"), 
              guide=guide_legend(override.aes=list(shape="circle")))), 
       queries=list(upset_query(set="NAR", fill="#d6a01d"), 
                    upset_query(set="NIR", fill="#d6a01d"), 
                    upset_query(set="NOR", fill="#d6a01d"), 
                    upset_query(set="NOS", fill="#d6a01d"), 
                    upset_query(set="SO", fill="#bf3553"), 
                    upset_query(set="FO", fill="#229453")), 
       annotations=list(
        "Completeness"=ggplot(mapping=aes(y=Completeness)) +
                       geom_jitter(size=0.1, alpha=0.3, height=0) +
                       geom_violin(linewidth=0.4, color="#6799CD", alpha=0.8) +
                       theme(axis.text=element_text(size=10, color="black"), 
                             axis.title=element_text(size=14, face="bold")) +
                       scale_y_continuous(limits=c(50,100))), 
       base_annotations=list(
        "intersection"=intersection_size(mode="inclusive_intersection", 
         counts=T, mapping=aes(fill=Phylum), text=list(size=2.5), 
         text_colors=c(on_background="black", on_bar="black")) +
         theme(axis.text=element_text(size=10, color="black"), 
               axis.title=element_text(size=14, face="bold"), 
               legend.text=element_text(size=12), 
               legend.key.size=unit(20, "pt"), 
               legend.title=element_text(size=14, face="bold")) +
         scale_y_continuous(limits=c(0,1500), breaks=seq(0,1500,250)) +
         guides(fill=guide_legend(ncol=1)) +
         ylab("MAGs in gene set") +
         scale_fill_ucscgb())) +
      patchwork::plot_layout(heights=c(0.25,1,0.5)) +
      theme(axis.text=element_text(size=12, color="black"), 
            axis.title=element_text(size=14, face="bold"))
p2 <- plot(euler(venn1, shape="ellipse"), 
           labels=list(font=2, cex=0.8), 
           edges=list(col="black", lex=1), 
           quantities=list(type="counts", cex=0.7), 
           fills=list(fill=c("white","#E1E0E1","#BDDEEC","#ED8183"), alpha=0.7))
p3 <- plot(euler(venn2, shape="ellipse"), 
           labels=list(font=2, cex=0.8), 
           edges=list(col="black", lex=1), 
           quantities=list(type="counts", cex=0.7), 
           fills=list(fill=c("#bf3553","#d6a01d","#229453"), alpha=0.7))
######
p22 <- as.ggplot(p2)
p33 <- as.ggplot(p3)
p11 <- p1 + geom_subview(subview=p22, x=24, y=15, w=15, h=15) +
            geom_subview(subview=p33, x=42, y=15, w=15, h=15)
ggsave("pdf/Functional_coupling.pdf", p11, width=12, height=10)
######
p4 <- upset(bins, genes, name="Gene sets", sort_sets=F, set_sizes=F, guides="over", 
           sort_intersections_by=c("degree"), sort_intersections="ascending", 
           n_intersections=200, min_size=1, max_size=3000, height_ratio=0.5, 
           mode="exclusive_intersection", 
       base_annotations=list(
        "intersection"=intersection_size(mode="exclusive_intersection", 
         counts=T, mapping=aes(fill=Phylum), text=list(size=4), 
         text_colors=c(on_background="black", on_bar="black")) +
         guides(fill=guide_legend(ncol=1)) +
         ylab("MAGs in gene set") +
         scale_fill_ucscgb()))
p4 <- upset(bins, genes, name="Gene sets", sort_sets=F, set_sizes=F, guides="over", 
           sort_intersections_by=c("degree"), sort_intersections=F, 
           mode="exclusive_intersection", height_ratio=0.5, 
       stripes=upset_stripes(colors=c("white")), 
       intersections=list("NAR","NIR","NOR","NOS", 
                     c("NAR","NIR"),c("NAR","NOR"),c("NAR","NOS"), 
                     c("NIR","NOR"),c("NIR","NOS"),c("NOR","NOS"), 
                     c("NAR","NIR","NOR"),c("NAR","NOR","NOS"),c("NIR","NOR","NOS")), 
       annotations=list(
        "Completeness"=ggplot(mapping=aes(y=Completeness)) +
                       geom_jitter(size=0.2, height=0) +
                       geom_violin(linewidth=0.5, color="#6799CD", alpha=0.8) +
                       theme(axis.text=element_text(size=10, color="black"), 
                             axis.title=element_text(size=14, face="bold")) +
                       scale_y_continuous(limits=c(50,100), breaks=seq(50,100,25))), 
       base_annotations=list(
        "intersection"=intersection_size(mode="exclusive_intersection", 
         counts=T, text=list(size=4), mapping=aes(fill=Phylum), 
         text_mapping=aes(label=paste0(!!upset_text_percentage(digits=1, sep=" "), "\n", 
                         "(", !!get_size_mode("exclusive_intersection"), ")")), 
         text_colors=c(on_background="black", on_bar="black")) +
         theme(axis.text=element_text(size=10, color="black"), 
               axis.title=element_text(size=14, face="bold"), 
               legend.text=element_text(size=12), 
               legend.key.size=unit(20, "pt"), 
               legend.title=element_text(size=14, face="bold")) +
         scale_y_continuous(limits=c(0,70), breaks=seq(0,69,20)) +
         guides(fill=guide_legend(ncol=1)) +
         ylab("MAGs in gene set") +
         scale_fill_ucscgb())) +
      patchwork::plot_layout(heights=c(0.25,1,0.5)) +
      theme(axis.text=element_text(size=12, color="black"), 
            axis.title=element_text(size=14, face="bold"))
ggsave("pdf/Functional_n_coupling.pdf", p4, width=12, height=10)


