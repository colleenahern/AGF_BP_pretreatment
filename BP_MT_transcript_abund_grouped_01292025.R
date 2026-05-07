# Colleen Ahern
# 01/29/2025
# I want to normalize my raw gene counts and add together gene groups (COGs, CAZymes, CAZyme type, etc.)
# not sure what kind of normalisation to use... I see several metatranscriptomics papers using FPKM or RPKM
# but I have no idea why they would choose those over TPM
# Let's start with TPM

install.packages("remotes")
remotes::install_github("davidrequena/drfun", force = TRUE)
library('DRnaSeq')
library(readr)

## Count matrix input
bpcoldata <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/bp_metadata_01132025.txt")
bpcoldata$name <- gsub("-", ".", bpcoldata$name)

bpcts <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/read_counts_gene_metassembly.tsv")

gl <- bpcts$Length

rownames(bpcts) <- bpcts[,1]
bpcts <- bpcts[,-(1:6)]
bpcoldata
rownames(bpcoldata) <- bpcoldata[,1]
bpcoldata <- bpcoldata[,-(1)]
head(bpcts)
bpcoldata

## Examine the count matrix and column data to see if they are consistent in terms of sample order
head(bpcts, 2)
bpcoldata
bpcoldata2 <- bpcoldata[rownames(bpcoldata) %in% colnames(bpcts), ]

## Rearrange
all(rownames(bpcoldata2) %in% colnames(bpcts))
all(rownames(bpcoldata2) == colnames(bpcts))

bpcts <- bpcts[, rownames(bpcoldata2)]
all(rownames(bpcoldata2) == colnames(bpcts))

# Generate tpm matrix
bptpm <- as.data.frame(tpm(bpcts, gl))
bptpm$GeneID <- rownames(bptpm)

# Input gene annotation info
gene_anno <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/gene_annotation_metassembly.tsv")

bptpm_geneanno <- merge(bptpm, gene_anno, by = "GeneID", all = TRUE)

library(dplyr)
library(stringr)
library(tidyr)
aa <- bptpm_geneanno %>%
  select(GeneID, everything())%>%
  mutate(dbcan_annotations=str_split(dbcan_annotations, "\\s*\\|\\s*")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
rownames(aa) <- aa$GeneID
aa <- aa[,-1]

unique(bptpm_geneanno$dbcan_annotations)
unique(aa$dbcan_annotations1)

caztable_gen <- data.frame(matrix(0, nrow = 1, ncol = 20))
colnames(caztable_gen) <- colnames(aa[,2:21])

aasub <- as.data.frame(aa[1:300,])
rownames(aasub) <- aasub$GeneID
aasub <- aasub[,-1]
aasub_caz <- aasub[is.na(aasub$dbcan_annotations1) == FALSE,]


# Trial run on subset of the data
for (j in 1:ncol(caztable_gen)) {
    for (i in 1:nrow(aasub)) {
      if(is.na(aasub$dbcan_annotations1[i]) == FALSE) {
        caztable_gen[j] = caztable_gen[j] + aasub[i,j]
      }
      else {
        caztable_gen[j] = caztable_gen[j] + 0
      }
    }
  }


# I THINK Malte gave me the raw counts so use my tpm normalized counts
# So now run this on the entire aa matrix
caztab_gen <- data.frame(matrix(0, nrow = 1, ncol = 20))
colnames(caztab_gen) <- colnames(aa[,1:20])

for (j in 1:ncol(caztab_gen)) {
    for (i in 1:nrow(aa)) {
      if(is.na(aa$dbcan_annotations1[i]) == FALSE) {
        caztab_gen[j] = caztab_gen[j] + aa[i,j]
      }
      else {
        caztab_gen[j] = caztab_gen[j] + 0
      }
    }
  }
    
# colnames(caztab_gen) <- bpcoldata2$condition_cba   no don't do this add another row with the conditions
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
caztab_gen_filt <- caztab_gen[,!(names(caztab_gen) %in% drops)]

bpcoldata3 <- bpcoldata2[,c("condition_cba","group")]
bpcoldata3 <- bpcoldata3[!(rownames(bpcoldata3) %in% drops),]

caztab_gen_avg <- data.frame(matrix(0, nrow = 1, ncol = 8))
colnames(caztab_gen_avg) <- (unique(bpcoldata2$condition_cba))

caztab_gen_avg$`CA + AS t1` <- mean(caztab_gen_filt$I5_MT_PLANC_BSE_S5,caztab_gen_filt$I7_CA15_1_MT_PLANC)
caztab_gen_avg$`CA + S3 + AS t1` <- mean(caztab_gen_filt$I9_CA61_1_MT_PLANC,caztab_gen_filt$I10_CA62_1_MT_PLANC)
caztab_gen_avg$`CA + AS t2` <- mean(caztab_gen_filt$I11_CA12_2_MT_PLANC,caztab_gen_filt$I12_CA14_2_MT_PLANC)
caztab_gen_avg$`CA + S3 + AS t2` <- mean(caztab_gen_filt$I1_CA61_2_MT_PLANC,caztab_gen_filt$I15.8_MT_PLANC_BSE_S15)

caztab_gen_avg$`PHA + AS t1` <- mean(caztab_gen_filt$I17_PH11_1_MT_PLANC,caztab_gen_filt$I16_PH8_1_MT_PLANC)
caztab_gen_avg$`PHA + G1 + AS t1` <- mean(caztab_gen_filt$I18_PH46_1_MT_PLANC,caztab_gen_filt$I2_PH47_1_MT_PLANC,caztab_gen_filt$I19_PH48_1_MT_PLANC)
caztab_gen_avg$`PHA + AS t2` <- mean(caztab_gen_filt$I21_PH11_2_MT_PLANC,caztab_gen_filt$I20_PH8_2_MT_PLANC)
caztab_gen_avg$`PHA + G1 + AS t2` <- mean(caztab_gen_filt$I3_PH46_2_MT_PLANC,caztab_gen_filt$I22_PH47_2_MT_PLANC,caztab_gen_filt$I4_PH48_2_MT_PLANC)

# wait redo this to make sure I did it correct
## AHH this actually does it correctly, use this method instead of above
caztab_gen_avg2 <- data.frame(matrix(0, nrow = 1, ncol = 8))
colnames(caztab_gen_avg2) <- (unique(bpcoldata2$condition_cba))

caztab_gen_avg2[,'CA + AS t1'] <- mean(as.matrix(caztab_gen_filt[, names(caztab_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + AS t1',])]))
caztab_gen_avg2[,'CA + S3 + AS t1'] <- mean(as.matrix(caztab_gen_filt[, names(caztab_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + S3 + AS t1',])]))
caztab_gen_avg2[,'CA + AS t2'] <- mean(as.matrix(caztab_gen_filt[, names(caztab_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + AS t2',])]))
caztab_gen_avg2[,'CA + S3 + AS t2'] <- mean(as.matrix(caztab_gen_filt[, names(caztab_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + S3 + AS t2',])]))

caztab_gen_avg2[,'PHA + AS t1'] <- mean(as.matrix(caztab_gen_filt[, names(caztab_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + AS t1',])]))
caztab_gen_avg2[,'PHA + G1 + AS t1'] <- mean(as.matrix(caztab_gen_filt[, names(caztab_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + G1 + AS t1',])]))
caztab_gen_avg2[,'PHA + AS t2'] <- mean(as.matrix(caztab_gen_filt[, names(caztab_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + AS t2',])]))
caztab_gen_avg2[,'PHA + G1 + AS t2'] <- mean(as.matrix(caztab_gen_filt[, names(caztab_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + G1 + AS t2',])]))

########################################################################################################################
## Make heatmap
## BP: Import metadata file

## BP: TPM-normalized ranscripts -- CAZymes added together
vary_43 <- caztab_gen_avg2


newnames <- lapply(
  colnames(vary_43),
  function(x) bquote(.(x)))
newnames
head(vary_43)
nrow(vary_43)
#quantile(rowSums(vary_43))

data_matrix_44 <- vary_43
data_matrix_44

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

my_palette <- colorRampPalette(c("black", "red"))(n = 299)

bp_caztpm_heatmap<-pheatmap(data_matrix_44, col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=F, cluster_cols=F, 
                                fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_row = NULL, labels_col = as.character(newnames), border_color=NA,breaks=seq(3000, 6000, length.out=300), angle_col = 45)





                             