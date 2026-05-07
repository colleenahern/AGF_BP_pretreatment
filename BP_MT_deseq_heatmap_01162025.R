# 01/16/2025
# Making heatmaps from my BPs_MT_deseq.R file
# DESeq analysis excel file in /Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_deseq folder

library(readr)
library(tidyverse)

# Plot heatmaps
# Malte is against a TPM filter but our lab always uses one?
# Do this without a TPM filter for now -- 01/16/2025 maybe incorporate one later

## BP: Import metadata file
asm <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/bp_metadata_01132025.txt")
rr <- c("I13-7_MT_PLANC_BSE_S13", "I174-7_MT_PLANC_BSE_S14")
asm <- subset(asm, !grepl("I13-7_MT_PLANC_BSE_S13|I174-7_MT_PLANC_BSE_S14", asm$name))
asm$name <- gsub("-", ".", asm$name)

## BP: Import list of differentially regulated genes for heatmap
vary_43 <- read.csv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_deseq/BP_DESeq_filt_01162025.csv", header = TRUE)

## BP: Import list of gene annotations
my_gene_col_44 <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/gene_annotation_metassembly.tsv")

nn <- colnames(vary_43)
mm <- nn[grepl("padj", nn)]
vary_43$padjmin <- apply(vary_43[,mm],1,min, na.rm = TRUE) 
vary_43 <- merge(vary_43, my_gene_col_44, by = "GeneID", all = T)
vary_43 <- vary_43[order(vary_43$padjmin,decreasing=FALSE),]

lala <- vary_43[vary_43$padjmin<0.05,]

newnames <- lapply(
  colnames(vary_43[,c(3,10,17,24,31,38,45,52)]),
  function(x) bquote(.(x)))
newnames
newnames <- gsub("_log2FoldChange", "", newnames)
newnames
newnames <- gsub("_v_", "\nvs.\n", newnames)
head(vary_43)
nrow(vary_43)
#quantile(rowSums(vary_43))
rownum_from <- 1
rownum_to <- 100

vary_44 <- vary_43[c(rownum_from:rownum_to),]
data_matrix_44 <- vary_44[,c(1,3,10,17,24,31,38,45,52,58,65)]
data_matrix_44
data_matrix_44$Protein <- ifelse(is.na(data_matrix_44$Description) == TRUE | data_matrix_44$Description == '-', 'hypothetical protein', data_matrix_44$Description)
data_matrix_44$Protein_short <- gsub("\\..*","",data_matrix_44$Protein)

lala <- grepl('GO', my_gene_col_44$Description)

# ann_cols <- data_matrix_44$GeneID
ann_cols <- data_matrix_44$Protein
#install.packages("RColorBrewer")
library("RColorBrewer")
display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n=12, name="Set3")
display.brewer.pal(n = 8, name = 'Set2')
brewer.pal(n=8, name="Set2")
display.brewer.pal(n = 9, name = 'Pastel1')
brewer.pal(n=9, name="Pastel1")
display.brewer.pal(n = 8, name = 'Pastel2')
brewer.pal(n=8, name="Pastel2")
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n=8, name="Dark2")
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n=12, name="Paired")

#ann_colors = list("GH Family"=c("GH1"="#00CC99","GH10"="#984EA3" ,"GH11"="#FF7F00","GH114"="#FFFF33","GH115"="#4DAF4A","GH13_28"="#DECBE4","GH16"="#996633","GH18"="#999999","GH2"="#E5D8BD","GH24"="#80B1D3","GH26"="#FDB462","GH3"="#B3DE69","GH30_5"="#FCCDE5","GH31"="#E41A1C","GH32"="#CCEBC5","GH39"="#F781BF","GH43"="#FB8072","GH45"="#FFCC00","GH48"="#7FC97F","GH5_1"="#666666","GH5_4"="#999999","GH5_5"="#F0027F","GH5_7"="#FF9933","GH6"="#FB8072","GH8"="#386CB0","GH9"="#FDCDAC","GH95"="#984EA3"),"Multiple Annotations"=c("Y"="#1F78B4","N"="#A6CEE3"))

#install.packages("pheatmap")
library("pheatmap")
library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.8, hjust = .5, rot = 360, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)

bp_regulation_heatmap<-pheatmap(data_matrix_44[,-c(1,10:13)], col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=TRUE, cluster_cols=F, 
                                fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_col = as.character(newnames), labels_row = as.expression(data_matrix_44$Protein_short), border_color=NA,breaks=seq(-4, 6, length.out=300), angle_col = 0)


save_pheatmap_png <- function(x, filename, width=1200, height=2000, res = 1500) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(bp_regulation_heatmap, "bp_regulation_heatmap.png")

########################################################################################################################
#CAZymes only

## BP: Import list of differentially regulated genes for heatmap
vary_43 <- read.csv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_deseq/BP_DESeq_filt_01162025.csv", header = TRUE)

## BP: Import list of gene annotations
my_gene_col_44 <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/gene_annotation_metassembly.tsv")
# filter for CAZymes
dbcananno <- my_gene_col_44[is.na(my_gene_col_44$dbcan_annotations) == FALSE,] # 7215 genes, selects the genes with dbcan hits (aka cazymes)
dbcananno_ce <- dbcananno[grepl("CE",dbcananno$dbcan_annotations),] #562 genes

nn <- colnames(vary_43)
mm <- nn[grepl("padj", nn)]
vary_43$padjmin <- apply(vary_43[,mm],1,min, na.rm = TRUE) 
vary_43 <- merge(vary_43, dbcananno, by = "GeneID") # 7154 genes
#vary_43 <- merge(vary_43, dbcananno_ce, by = "GeneID") #558 genes
vary_43 <- vary_43[order(vary_43$padjmin,decreasing=FALSE),] 

lala <- vary_43[vary_43$padjmin>0.05,] # 6449 cazyme genes, 521 CE genes

newnames <- lapply(
  colnames(vary_43[,c(3,10,17,24,31,38,45,52)]),
  function(x) bquote(.(x)))
newnames
newnames <- gsub("_log2FoldChange", "", newnames)
newnames
newnames <- gsub("_v_", "\nvs.\n", newnames)
head(vary_43)
nrow(vary_43)
#quantile(rowSums(vary_43))
rownum_from <- 1
rownum_to <- 100

vary_44 <- vary_43[c(rownum_from:rownum_to),]
data_matrix_44 <- vary_44[,c(1,3,10,17,24,31,38,45,52,58,65,80)]
data_matrix_44
#data_matrix_44$Protein <- ifelse(is.na(data_matrix_44$Description) == TRUE | data_matrix_44$Description == '-', 'hypothetical protein', data_matrix_44$Description)

ann_cols <- data_matrix_44$dbcan_annotations
#install.packages("RColorBrewer")
library("RColorBrewer")
display.brewer.pal(n = 12, name = 'Set3')
brewer.pal(n=12, name="Set3")
display.brewer.pal(n = 8, name = 'Set2')
brewer.pal(n=8, name="Set2")
display.brewer.pal(n = 9, name = 'Pastel1')
brewer.pal(n=9, name="Pastel1")
display.brewer.pal(n = 8, name = 'Pastel2')
brewer.pal(n=8, name="Pastel2")
display.brewer.pal(n = 8, name = 'Dark2')
brewer.pal(n=8, name="Dark2")
display.brewer.pal(n = 12, name = 'Paired')
brewer.pal(n=12, name="Paired")

#install.packages("pheatmap")
library("pheatmap")
library(grid)
draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 0.8, hjust = .5, rot = 360, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)

bp_regulation_heatmap<-pheatmap(data_matrix_44[,-c(1,10:12)], col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=TRUE, cluster_cols=F, 
                                fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_col = as.character(newnames), labels_row = as.expression(data_matrix_44$dbcan_annotations), border_color=NA,breaks=seq(-4, 6, length.out=300), angle_col = 0)


save_pheatmap_png <- function(x, filename, width=1200, height=2000, res = 1500) {
  png(filename, width = width, height = height, res = res)
  grid::grid.newpage()
  grid::grid.draw(x$gtable)
  dev.off()
}

save_pheatmap_png(bp_regulation_heatmap, "bp_regulation_heatmap.png")