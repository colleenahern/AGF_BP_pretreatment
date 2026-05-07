# 12/11/2024
# Colleen Ahern
# Using qiime2 output from Magda_BPs_experiment to create NMDS plot

if (!requireNamespace("devtools", quietly = TRUE)){install.packages("devtools")}
devtools::install_github("jbisanz/qiime2R")

library(qiime2R)
metadata_a_grouped <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast/metadata_a_grouped.txt", header=T, stringsAsFactors=F, sep="\t")

setwd("/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast")

pswd <- qza_to_phyloseq(features="arc2analysis/16S-table-nocmb-rarefied-10000_filtered.qza",
                        tree="arc2analysis/rooted-16S-tree-filteredSVs.qza", 
                        taxonomy = "output_arc2/qiime2/input/taxonomy.qza", 
                        metadata= "metadata_a_grouped_nmds.txt")

pswd.prop <- transform_sample_counts(pswd, function(otu) otu/sum(otu))

ord.nmds.bray_wd <- ordinate(pswd.prop, method="NMDS", distance="bray")

library(ggforce)
# deg_plot <- plot_ordination(pswd.prop, ord.nmds.bray_wd, color="Substrate", shape = "Day", title="Bray NMDS") + scale_shape_binned(name = "Day", limits=c(0,90), breaks=c(0,45,90)) + geom_point(size = 4) 
deg_plot <- plot_ordination(pswd.prop, ord.nmds.bray_wd, color = "Group", title="Bray NMDS") + geom_point(size = 3) 

deg_plot <- plot_ordination(pswd.prop, ord.nmds.bray_wd, color = "Testcond", shape = "Time", title="Bray NMDS") + geom_point(size = 4) 

deg_plot

