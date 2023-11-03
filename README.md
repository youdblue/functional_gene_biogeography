# functional_gene_biogeography
This repository contains R codes used to generate results and figures for the unpublished manuscript "Biogeography and functional coupling of iron, sulfur and nitrogen cycling genes in acid mine drainage sediment microbiomes".

1_catalog  
1_catalog_a.R: R script used to generate the map of sampling sites.  
1_catalog_b.R: R script used to generate the accumulation curves of gene catalog.  
1_catalog_c.R: R script used to generate the relative abundance of gene catalog annotated by different databases.  
1_catalog_d.R: R script used to generate the taxonomic composition of gene catalog at phylum level.  
1_catalog_e.R: R script used to generate the relationship between taxonomic dissimilarity and functional dissimilarity.  

2_pattern  
2_pattern_ab.R: R script used to generate the relationship between functional and taxonomic diversity and latitude.  
2_pattern_cd.R: R script used to generate the distance decay relationship of functional and taxonomic traits at different geographic distance level.  
2_pattern_ef.R: R script used to generate the distance decay relationship of functional and taxonomic traits of generalists and specialists.  

3_driver  
3_driver_a.R: R script used to generate the relationship between biotic and abiotic variables by correlation analysis and mantel tests.  
3_driver_bc.R: R script used to perform principal coordinate analysis (PCoA) of functional and taxonomic traits colored by latitude.  
3_driver_d.R: R script used to perform variance partitioning analysis (VPA).  
3_driver_e.R: R script used to perform null model analysis.  
s3_1.R: R script used to perform random forest analysis.  

4_gene  
4_gene.R: R script used to generate the relationship between the abundance of key functional genes and geo-environmental factors.  
s4_1.R: R script used to generate the relationship between the diversity of key functional genes and geo-environmental factors.  
s4_2.R: R script used to perform correlation analysis of the abundance of iron, sulfur and nitrogen cycling genes.  

5_genus  
5_genus.R: R script used to generate the taxonomic composition of key functional genes at genus level.  

6_functional_coupling  
Functional_coupling.R: R script used to show the coupling mechanisms of denitrification genes and iron or sulfur oxidation genes in metagenome-assembled genomes.  
