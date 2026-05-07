# Colleen Ahern
# 05/19/2025

# Continuation of BP_MT_transabund_grouped_CAZyme_03182025.R script
# I want to normalizing (or not normalize, depending if I need to do one-way anova/wilcoxin test or ancomBC) 
# my raw gene counts and add together CAZyme SUBtype
# before I just did broad CAZyme type (GH, CE, etc.)
# Continue using TPM normalization if I normalize

# SKIP TO LINE 790 FOR MOST UPDATED

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

# split multiple dbcan annotations into individual columns for group counting
aa <- bptpm_geneanno %>%
  select(GeneID, everything())%>%
  mutate(dbcan_annotations=str_split(dbcan_annotations, "\\s*\\|\\s*")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:46,1)]

unique(bptpm_geneanno$dbcan_annotations)
ab<-as.vector(as.matrix(aa[,grepl( "dbcan_annotations" , names(aa))]))
unique(ab)
length(unique(ab))

# make subset of aa data for a test grouping run
bb <- aa[is.na(aa$dbcan_annotations3) == FALSE,]
aasub <- bb[1:5,]
aasub$dbcan_annotations1[is.na(aasub$dbcan_annotations_b1)] <- "-"
aasubccs <- aasub[,c(42:45)]

uab <-unique(ab)
uab
uab <- uab[!is.na(uab)]
rn <- append(uab, "-")
rn
dbcantable_gen <- data.frame(matrix(0, nrow = length(unique(rn)), ncol = 20))
colnames(dbcantable_gen) <- colnames(aa[,1:20])
rownames(dbcantable_gen) <- rn

# Trial run on subset of the data
for (j in 1:ncol(dbcantable_gen)) {
  for (k in 1:length(aasub[,grepl( "dbcan_annotations" , names(aasub))])) { # Number of separate cazyme types broken up
    for (i in 1:nrow(aasub)) {
      if(is.na(aasub[i,paste("dbcan_annotations",k, sep = "")]) == FALSE) {
        dbcantable_gen[rn[rownames(dbcantable_gen) %in% aasub[i,paste("dbcan_annotations",k, sep = "")]],j] = dbcantable_gen[rn[rownames(dbcantable_gen) %in% aasub[i,paste("dbcan_annotations",k, sep = "")]],j] + aasub[i,j]
      }
      else {
        dbcantable_gen = dbcantable_gen
      }
    }
  }
}


## Now do it on the entire data
aa$dbcan_annotations1[is.na(aa$dbcan_annotations1)] <- "-"

uab <-unique(ab)
uab
uab <- uab[!is.na(uab)]
rn <- append(uab, "-")
rn
dbcantable_gen <- data.frame(matrix(0, nrow = length(unique(rn)), ncol = 20))
colnames(dbcantable_gen) <- colnames(aa[,1:20])
rownames(dbcantable_gen) <- rn

for (j in 1:ncol(dbcantable_gen)) {
  for (k in 1:length(aa[,grepl( "dbcan_annotations" , names(aa))])) { # Number of separate cazyme types broken up
    for (i in 1:nrow(aa)) {
      if(is.na(aa[i,paste("dbcan_annotations",k, sep = "")]) == FALSE) {
        dbcantable_gen[rn[rownames(dbcantable_gen) %in% aa[i,paste("dbcan_annotations",k, sep = "")]],j] = dbcantable_gen[rn[rownames(dbcantable_gen) %in% aa[i,paste("dbcan_annotations",k, sep = "")]],j] + aa[i,j]
      }
      else {
        dbcantable_gen = dbcantable_gen
      }
    }
  }
}
dbcantable_gen$CAZID <- rownames(dbcantable_gen)
write_tsv(dbcantable_gen,'/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_tpm_05192025.tsv')
dbcantable_genr <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_tpm_05192025.tsv')
dbcantable_genr <- as.data.frame(dbcantable_genr)
rownames(dbcantable_genr) <- dbcantable_genr$CAZID
drop <- "CAZID"
dbcantable_genr <- dbcantable_genr[,!names(dbcantable_genr) %in% drop]

# Taking an average is useful for heatmaps, but I wouldn't use the average grouped counts for statistical testing
# Take average of replicates
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
dbcantable_gen_filt <- dbcantable_genr[,!names(dbcantable_genr) %in% drops]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(dbcantable_gen_filt),]
dbcantable_gen_avg <- data.frame(matrix(0, nrow = length(rn), ncol = 8))
colnames(dbcantable_gen_avg) <- (unique(bpcoldata3$condition_cba))
rownames(dbcantable_gen_avg) <- rownames(dbcantable_gen_filt)

dbcantable_gen_avg[,'CA + AS t1'] <- rowMeans(as.matrix(dbcantable_gen_filt[, names(dbcantable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + AS t1',])]))
dbcantable_gen_avg[,'CA + S3 + AS t1'] <- rowMeans(as.matrix(dbcantable_gen_filt[, names(dbcantable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + S3 + AS t1',])]))
dbcantable_gen_avg[,'CA + AS t2'] <- rowMeans(as.matrix(dbcantable_gen_filt[, names(dbcantable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + AS t2',])]))
dbcantable_gen_avg[,'CA + S3 + AS t2'] <- rowMeans(as.matrix(dbcantable_gen_filt[, names(dbcantable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + S3 + AS t2',])]))

dbcantable_gen_avg[,'PHA + AS t1'] <- rowMeans(as.matrix(dbcantable_gen_filt[, names(dbcantable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + AS t1',])]))
dbcantable_gen_avg[,'PHA + G1 + AS t1'] <- rowMeans(as.matrix(dbcantable_gen_filt[, names(dbcantable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + G1 + AS t1',])]))
dbcantable_gen_avg[,'PHA + AS t2'] <- rowMeans(as.matrix(dbcantable_gen_filt[, names(dbcantable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + AS t2',])]))
dbcantable_gen_avg[,'PHA + G1 + AS t2'] <- rowMeans(as.matrix(dbcantable_gen_filt[, names(dbcantable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + G1 + AS t2',])]))

########################################################################################################################
## Make heatmap
## BP: Import metadata file

## BP: TPM-normalized ranscripts -- CAZymes added together
vary_43 <- dbcantable_gen_avg

newnames <- lapply(
  colnames(vary_43),
  function(x) bquote(.(x)))
newnames
head(vary_43)
nrow(vary_43)
#quantile(rowSums(vary_43))

data_matrix_44 <- vary_43
drop <- c("-")
data_matrix_44 <- data_matrix_44[!(row.names(data_matrix_44) %in% drop),]

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

dbcantpm_heatmap<-pheatmap(data_matrix_44, col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=F, cluster_cols=F, 
                              fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_row = rownames(data_matrix_44), labels_col = as.character(newnames), border_color=NA,breaks=seq(0, 200, length.out=300), angle_col = 45)

########################################################################################################################################
########################################################################################################################################
# If you want raw counts instead of TPM normalized

library('DRnaSeq')
library(readr)

## Count matrix input
bpcoldata <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/bp_metadata_01132025.txt")
bpcts <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/read_counts_gene_metassembly.tsv")

bpcts <- bpcts %>%
  column_to_rownames("GeneID")
bpcts <- bpcts[,-(1:5)]
head(bpcts)

bpcoldata <- bpcoldata %>%
  column_to_rownames("name")
head(bpcoldata)

## Examine the count matrix and column data to see if they are consistent in terms of sample order
head(bpcts, 2)
bpcoldata
bpcoldata2 <- bpcoldata[rownames(bpcoldata) %in% colnames(bpcts), ]

## Rearrange
all(rownames(bpcoldata2) %in% colnames(bpcts))
all(rownames(bpcoldata2) == colnames(bpcts))

bpcts <- bpcts[, rownames(bpcoldata2)]
all(rownames(bpcoldata2) == colnames(bpcts))
bpcts$GeneID <- rownames(bpcts)

# Input gene annotation info
gene_anno <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/gene_annotation_metassembly.tsv")

bpcts_geneanno <- merge(bpcts, gene_anno, by = "GeneID", all = TRUE)
bpcts_geneanno2 <- bpcts_geneanno
bpcts_geneanno2$eggNOG_OGs2 <- bpcts_geneanno2$eggNOG_OGs
bpcts_geneanno2 <- bpcts_geneanno2[,c(1:25,44,26:43)]
bpcts_geneanno2$eggNOG_OGs2 <- gsub(",.*?\\|", "", bpcts_geneanno2$eggNOG_OGs)
bpcts_geneanno2$eggNOG_OGs2 <- gsub(".*root", "", bpcts_geneanno2$eggNOG_OGs2)
bpcts_geneanno2$eggNOG_OGs2 <- gsub("([[:lower:]])([[:upper:]])", "\\1 \\2", bpcts_geneanno2$eggNOG_OGs2)
bpcts_geneanno2$eggNOG_OGs2 <- gsub("unclassified", " unclassified", bpcts_geneanno2$eggNOG_OGs2)
# bpcts_geneanno2$KEGG_ko <- gsub("ko:","",bpcts_geneanno2$KEGG_ko)

library(dplyr)
library(stringr)
library(tidyr)

# dbcan annotations
aa <- bpcts_geneanno2 %>%
  dplyr::select(GeneID, everything())%>%
  mutate(dbcan_annotations=str_split(dbcan_annotations, "\\s*\\|\\s*")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:47,1)]

unique(bpcts_geneanno2$dbcan_annotations)
length(unique(bpcts_geneanno2$dbcan_annotations))
ab<-as.vector(as.matrix(aa[,grepl( "dbcan_annotations" , names(aa))]))
unique(ab)
length(unique(ab))

# make subset of aa data for a test grouping run
bb <- aa[is.na(aa$dbcan_annotations3) == FALSE,]
aasub <- bb[1:5,]
aasub$dbcan_annotations1[is.na(aasub$dbcan_annotations_b1)] <- "-"
aasubccs <- aasub[,c(25,43:46)]

uab <-unique(ab)
uab
uab <- uab[!is.na(uab)]
rn <- append(uab, "-")
rn
dbcantable_gen <- data.frame(matrix(0, nrow = length(rn), ncol = 20))
colnames(dbcantable_gen) <- colnames(aa[,1:20])
rownames(dbcantable_gen) <- rn

# Trial run on subset of the data
for (j in 1:ncol(dbcantable_gen)) {
  for (k in 1:length(names(aasub)[grepl( "dbcan_annotations" , names(aasub))])) { # Number of separate cazyme types broken up
    for (i in 1:nrow(aasub)) {
      if(is.na(aasub[i,paste("dbcan_annotations",k, sep = "")]) == FALSE) {
        dbcantable_gen[rn[rn %in% aasub[i,paste("dbcan_annotations",k, sep = "")]],j] = dbcantable_gen[rn[rn %in% aasub[i,paste("dbcan_annotations",k, sep = "")]],j] + aasub[i,j]
      }
      else {
        dbcantable_gen = dbcantable_gen
      }
    }
  }
}

for (k in 1:length(names(aasub)[grepl( "dbcan_annotations" , names(aasub))])) { # Number of separate CAZyme types broken up
  for (i in 1:nrow(aasub)) {
    if(aasub[i,paste("dbcan_annotations",k, sep = "")] %in% rownames(dbcantable_gen) == TRUE ) {
      dbcantable_gen[rn[rn %in% aasub[i,paste("dbcan_annotations",k, sep = "")]],"taxa"] <- paste(dbcantable_gen[rn[rn %in% aasub[i,paste("dbcan_annotations",k, sep = "")]],"taxa"],"; ", aasub$eggNOG_OGs2[i], sep = "")
    }
    else {
      dbcantable_gen = dbcantable_gen
    }
  }
}


## Now do it on the entire data
aa$dbcan_annotations1[is.na(aa$dbcan_annotations1)] <- "-"

ab<-as.vector(as.matrix(aa[,grepl( "dbcan_annotations" , names(aa))]))
uab <-unique(ab)
uab
uab <- uab[!is.na(uab)]
uab
rn <- uab
dbcantable_gen_raw <- data.frame(matrix(0, nrow = length(rn), ncol = 20))
colnames(dbcantable_gen_raw) <- colnames(aa[,1:20])
rownames(dbcantable_gen_raw) <- rn

for (j in 1:ncol(dbcantable_gen_raw)) {
  for (k in 1:length(names(aa)[grepl( "dbcan_annotations" , names(aa))])) { # Number of separate cazyme types broken up
    for (i in 1:nrow(aa)) {
      if(is.na(aa[i,paste("dbcan_annotations",k, sep = "")]) == FALSE) {
        dbcantable_gen_raw[rn[rn %in% aa[i,paste("dbcan_annotations",k, sep = "")]],j] = dbcantable_gen_raw[rn[rn %in% aa[i,paste("dbcan_annotations",k, sep = "")]],j] + aa[i,j]
      }
      else {
        dbcantable_gen_raw = dbcantable_gen_raw
      }
    }
  }
}

dbcantable_gen_raw$CAZID <- rownames(dbcantable_gen_raw)
write_tsv(dbcantable_gen_raw, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_raw_05192025.tsv') # I have checked this multiple times, even on 10/03/2025 and I trust that this is accurate
dbcantable_gen_rawr <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_raw_05192025.tsv')
dbcantable_gen_rawr <- as.data.frame(dbcantable_gen_rawr)
rownames(dbcantable_gen_rawr) <- dbcantable_gen_rawr$CAZID
drop <- "CAZID"
dbcantable_gen_rawr <- dbcantable_gen_rawr[,!(names(dbcantable_gen_rawr) %in% drop)]

for (k in 1:length(aa[,grepl("dbcan_annotations" , names(aa))])) { # Number of separate CAZyme types broken up
  for (i in 1:nrow(aa)) {
    if(aa[i,paste("dbcan_annotations",k, sep = "")] %in% rownames(dbcantable_gen_rawr) == TRUE ) {
      dbcantable_gen_rawr[rn[rownames(dbcantable_gen_rawr) %in% aa[i,paste("dbcan_annotations",k, sep = "")]],"taxa"] <- paste(dbcantable_gen_rawr[rn[rownames(dbcantable_gen_rawr) %in% aa[i,paste("dbcan_annotations",k, sep = "")]],"taxa"],"; ", aa$eggNOG_OGs2[i], sep = "")
    }
    else {
      dbcantable_gen_rawr = dbcantable_gen_rawr
    }
  }
}

dbcantable_gen_rawr$CAZID <- rownames(dbcantable_gen_rawr)
write_tsv(dbcantable_gen_rawr, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_raw_taxa_05202025.tsv')
dbcantable_gen_rawrt <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_raw_taxa_05202025.tsv')
dbcantable_gen_rawrt <- dbcantable_gen_rawrt %>%
  column_to_rownames("CAZID")

# drop any outlier samples (based on PCA analysis)
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
dbcantable_gen_rawrt_filt <- dbcantable_gen_rawrt[,!(names(dbcantable_gen_rawrt) %in% drops)]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(dbcantable_gen_rawrt_filt),]

#######################################################################################################################################
#### CA + S3 + AS vs. CA + AS comparison

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "CA + AS" | bpcoldata3$group == "CA + S3 + AS",]
drop <- "-"

dbcantable_raw_filt_sub <- dbcantable_gen_rawrt_filt[,colnames(dbcantable_gen_rawrt_filt) %in% rownames(bpcoldata3_sub)]
dbcantable_raw_filt_sub_noNA <- dbcantable_raw_filt_sub[!(row.names(dbcantable_raw_filt_sub) %in% drop), ]
dbcantable_raw_taxa <- dbcantable_gen_rawrt_filt[!(row.names(dbcantable_gen_rawrt_filt) %in% drop),]

colnames(dbcantable_raw_filt_sub) == rownames(bpcoldata3_sub)
colnames(dbcantable_raw_filt_sub_noNA) == rownames(bpcoldata3_sub)

##############################################################################################################################
#### Now do for PHA + G1 + AS vs. PHA + AS comparison

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "PHA + G1 + AS" | bpcoldata3$group == "PHA + AS",]
drop <- "-"

dbcantable_raw_filt_sub <- dbcantable_gen_rawrt_filt[,colnames(dbcantable_gen_rawrt_filt) %in% rownames(bpcoldata3_sub)]
dbcantable_raw_filt_sub_noNA <- dbcantable_raw_filt_sub[!(row.names(dbcantable_raw_filt_sub) %in% drop), ]
dbcantable_raw_taxa <- dbcantable_gen_rawrt_filt[!(row.names(dbcantable_gen_rawrt_filt) %in% drop),]

colnames(dbcantable_raw_filt_sub) == rownames(bpcoldata3_sub)
colnames(dbcantable_raw_filt_sub_noNA) == rownames(bpcoldata3_sub)

#######################################################################################################################################
# ANCOMBC

library(phyloseq)
library(ANCOMBC)
library(tidyverse)
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))
library(dplyr)
library(ggplot2)

otumat = as.matrix(dbcantable_raw_filt_sub_noNA)

taxmat = matrix(sample(letters, 50, replace = TRUE), nrow = nrow(otumat), ncol = 8)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Taxa")
taxmat[,"Genus"] <- rownames(taxmat)
all(rownames(dbcantable_raw_taxa) == rownames(taxmat))
taxmat[,"Taxa"] <- dbcantable_raw_taxa[, "taxa"]
taxmat[,"Taxa"] <- gsub("NA; ","", taxmat[,"Taxa"])
taxmat

class(otumat)
class(taxmat)

library("phyloseq")
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU
TAX

sampledata <- sample_data(bpcoldata3_sub)
sampledata
physeq = phyloseq(OTU, TAX, sampledata)
physeq
physeq@sam_data

library(ANCOMBC)
out = ancombc(data = physeq, tax_level = "Genus", 
              formula = "group", 
              p_adj_method = "BH",
              conserve = TRUE,
              alpha = 0.01,
              verbose = TRUE)

res = out$res
res_global = out$res_global

## ANCOMBC primary result
# LFC
tab_lfc = res$lfc
col_name = c("CAZyme", "Intercept", "CA + S3 + AS")
col_name = c("CAZyme", "Intercept", "PHA + G1 + AS")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# SE
tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Test statistic
tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# P-values
tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Adjusted p-values
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Differentially abundant taxa
tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant CAZymes from the Primary Result")

cadiff <- tab_diff[tab_diff$`CA + S3 + AS` == TRUE,]
phadiff <- tab_diff[tab_diff$`PHA + G1 + AS` == TRUE,]

tab_lfc_filt <- tab_lfc[tab_lfc$CAZyme %in% cadiff$CAZyme,]
colnames(tab_lfc_filt) <- c("CAZyme", "Intercept LFC", "CA + S3 + AS LFC")
tab_lfc_filt <- tab_lfc[tab_lfc$CAZyme %in% phadiff$CAZyme,]
colnames(tab_lfc_filt) <- c("CAZyme", "Intercept LFC", "PHA + G1 + AS LFC")

tab_se_filt <- tab_se[tab_se$CAZyme %in% cadiff$CAZyme,]
colnames(tab_se_filt) <- c("CAZyme", "Intercept SE", "CA + S3 + AS SE")
tab_se_filt <- tab_se[tab_se$CAZyme %in% phadiff$CAZyme,]
colnames(tab_se_filt) <- c("CAZyme", "Intercept SE", "PHA + G1 + AS SE")

tab_w_filt <- tab_w[tab_w$CAZyme %in% cadiff$CAZyme,]
colnames(tab_w_filt) <- c("CAZyme", "Intercept W", "CA + S3 + AS W")
tab_w_filt <- tab_w[tab_w$CAZyme %in% phadiff$CAZyme,]
colnames(tab_w_filt) <- c("CAZyme", "Intercept W", "PHA + G1 + AS W")

tab_p_filt <- tab_p[tab_p$CAZyme %in% cadiff$CAZyme,]
colnames(tab_p_filt) <- c("CAZyme", "Intercept P", "CA + S3 + AS P")
tab_p_filt <- tab_p[tab_p$CAZyme %in% phadiff$CAZyme,]
colnames(tab_p_filt) <- c("CAZyme", "Intercept P", "PHA + G1 + AS P")

tab_q_filt <- tab_q[tab_q$CAZyme %in% cadiff$CAZyme,]
colnames(tab_q_filt) <- c("CAZyme", "Intercept Q", "CA + S3 + AS Q")
tab_q_filt <- tab_q[tab_q$CAZyme %in% phadiff$CAZyme,]
colnames(tab_q_filt) <- c("CAZyme", "Intercept Q", "PHA + G1 + AS Q")

tab_diff_filt <- tab_diff[tab_diff$CAZyme %in% cadiff$CAZyme,]
colnames(tab_diff_filt) <- c("CAZyme", "Intercept diff", "CA + S3 + AS diff")
tab_diff_filt <- tab_diff[tab_diff$CAZyme %in% phadiff$CAZyme,]
colnames(tab_diff_filt) <- c("CAZyme", "Intercept diff", "PHA + G1 + AS diff")

#put all data frames into list
df_list <- list(tab_lfc_filt, tab_se_filt, tab_w_filt, tab_p_filt, tab_q_filt, tab_diff_filt)

#merge all data frames in list
res_filt <- df_list %>% purrr::reduce(full_join, by='CAZyme')

# save list
#writexl::write_xlsx(res_filt, "")

### Visualization of differentially abundant taxa heatmap

sig_taxa = tab_diff %>%
  dplyr::filter(`CA + S3 + AS` == TRUE) %>%
  .$CAZyme
sig_taxa = tab_diff %>%
  dplyr::filter(`PHA + G1 + AS` == TRUE) %>%
  .$CAZyme

tab_lfc = res$lfc
col_name = c("CAZyme", "Intercept", "CA + S3 + AS vs. CA + AS")
col_name = c("CAZyme", "Intercept", "PHA + G1 + AS vs. PHA + AS")
colnames(tab_lfc) = col_name
tab_lfc <- tab_lfc[,-2]

df_CAS3AS = tab_lfc %>%
  filter(CAZyme %in% sig_taxa)

dbcantable_raw_taxa$taxa <- gsub("NA; ","",dbcantable_raw_taxa$taxa)
dbcantable_raw_taxa$CAZyme <- rownames(dbcantable_raw_taxa)
df_CAS3AS_fin <- merge(df_CAS3AS, dbcantable_raw_taxa[, c("taxa","CAZyme")], by = "CAZyme")
write_tsv(df_CAS3AS_fin, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/CAS3AS_CAZyme_subs_ancombc_05212025.tsv")
bbb <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/CAS3AS_CAZyme_subs_ancombc_05212025.tsv")

df_PHAG1AS = tab_lfc %>%
  filter(CAZyme %in% sig_taxa) 

dbcantable_raw_taxa$taxa <- gsub("NA; ","",dbcantable_raw_taxa$taxa)
dbcantable_raw_taxa$CAZyme <- rownames(dbcantable_raw_taxa)
df_PHAG1AS_fin <- merge(df_PHAG1AS, dbcantable_raw_taxa[, c("taxa","CAZyme")], by = "CAZyme")
write_tsv(df_PHAG1AS_fin, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/PHAG1AS_CAZyme_subs_ancombc_05212025.tsv")
bbb <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/PHAG1AS_CAZyme_subs_ancombc_05212025.tsv")

#####################################################################################################################################
# Make heatmap

df_heat = df_CAS3AS %>%
  pivot_longer(cols = -one_of("CAZyme"),
               names_to = "region", values_to = "value") %>%
  mutate(value = round(value, 2))
df_heat$CAZyme = factor(df_heat$CAZyme, levels = sort(sig_taxa))
# df_heatmod <- merge(df_heat, modtab,  by = "Module", all.x = TRUE)

df_heat = df_PHAG1AS %>%
  pivot_longer(cols = -one_of("CAZyme"),
               names_to = "region", values_to = "value") %>%
  mutate(value = round(value, 2))
df_heat$CAZyme = factor(df_heat$CAZyme, levels = sort(sig_taxa))
# df_heatmod <- merge(df_heat, modtab,  by = "Module", all.x = TRUE)

lo = floor(min(df_heat$value))
up = ceiling(max(df_heat$value))
mid = (lo + up)/2
# df_heat_filt <- df_heat[abs(df_heat$value) > 1,]
# df_heatmod_filt <- df_heatmod[abs(df_heatmod$value) > 1,]
# df_heatmod_filt$Label <- paste(df_heatmod_filt$Module, ": ", df_heatmod_filt$Function, sep = "")

# df_heatmod_filt2 <- df_heatmod_filt[,c(1:3,6)]
# write_xlsx(df_heatmod_filt2, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAS3AS_KEGGancombc_05052025.xlsx")
# 
# df_heatmod_filt2 <- df_heatmod_filt[,c(1:3,6)]
# write_xlsx(df_heatmod_filt2, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/PHAG1AS_KEGGancombc_05052025.xlsx")

p_heat = df_heat %>%
  ggplot(aes(x = region, y = CAZyme, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(region, CAZyme, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, size = 5,
       title = "Log fold changes for globally significant CAZymes") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 14))
p_heat





####################################################################
# 10/05/2025
# making volcano plots instead of heatmmaps to show differentially expressed CAZymes
# taken from https://biostatsquid.com/volcano-plots-r-tutorial/

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

# for CAS3AS vs CAAS comparison
# Import Ancombc results

tab_lfc2 <- tab_lfc %>%
  rename("CA + S3 + AS vs. CA + AS" = "LFC")
tab_lfc2 <- tab_lfc %>%
  rename("CA + S3 + AS" = "LFC")

tab_q2 <- tab_q %>%
  rename("CA + S3 + AS" = "padj")

df <- merge(tab_lfc2, tab_q2, by = "CAZyme")

ggplot(data = df, aes(x = LFC, y = -log10(padj))) +
  geom_point()

# Add threshold lines
ggplot(data = df, aes(x = LFC, y = -log10(padj))) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point() 

# Biostatsquid theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

ggplot(data = df, aes(x = LFC, y = -log10(padj))) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point() 

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
df$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$LFC > 1 & df$padj < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$LFC < -1 & df$padj < 0.05] <- "DOWN"
 
# Explore a bit
head(df[order(df$padj) & df$diffexpressed == 'DOWN', ])

ggplot(data = df, aes(x = LFC, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  labs(color='Differential Expression') 


# Edit axis labels and limits
ggplot(data = df, aes(x = LFC, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 50), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) # to customise the breaks in the x axis

# Note. with coord_cartesian() even if we have genes with p-values or log2FC ourside our limits, they will still be plotted.

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
dlup <- df[df$diffexpressed == "UP",]
df$delabel <- ifelse(df$CAZyme %in% head(dlup[order(dlup$padj), "CAZyme"], 10), df$CAZyme, NA)

dldown <- df[df$diffexpressed == "DOWN",]
df$delabel <- ifelse(df$CAZyme %in% head(dldown[order(dldown$padj), "CAZyme"], 10), df$CAZyme, df$delabel)

library(ggrepel)

ggplot(data = df, aes(x = LFC, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 50), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) + # to customise the breaks in the x axis
  ggtitle('CAZyme expression in S3 + CA + AS vs. \nCA + AS comparison') + # Plot title 
  geom_text_repel(max.overlaps = Inf) # To show all labels 

dff <- df[is.na(df$delabel) == FALSE,]
dbcantable_raw_taxa <- dbcantable_raw_taxa %>%
  rownames_to_column("CAZyme")
dfff <- merge(dff,dbcantable_raw_taxa, by =)

## 10/05/2025 ^ ok so for the graph above I am labeling the 10 most significantly expressed (lowest padj values) upregulated CAZyme subgroup and 
# the 10 most significantly expressed (lowest padj values) downregulated CAZyme subgroup and 
# This also includes an abs(1) LFC cutoff, which excludes some of the lowest p value CAZ subgroups, so maybe change that? Talk to Magda
# see if we want to keep the LFC cutoff or get rid of it

# for G1PHAAS vs PHAAS comparison

tab_lfc2 <- tab_lfc %>%
  rename("PHA + G1 + AS vs. PHA + AS" = "LFC")

tab_q2 <- tab_q %>%
  rename("PHA + G1 + AS" = "padj")

df <- merge(tab_lfc2, tab_q2, by = "CAZyme")

theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
df$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$LFC > 1 & df$padj < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$LFC < -1 & df$padj < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
dlup <- df[df$diffexpressed == "UP",]
df$delabel <- ifelse(df$CAZyme %in% head(dlup[order(dlup$padj), "CAZyme"], 10), df$CAZyme, NA)

dldown <- df[df$diffexpressed == "DOWN",]
df$delabel <- ifelse(df$CAZyme %in% head(dldown[order(dldown$padj), "CAZyme"], 10), df$CAZyme, df$delabel)

ggplot(data = df, aes(x = LFC, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 15), xlim = c(-6, 6)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-6, 6, 2)) + # to customise the breaks in the x axis
  ggtitle('CAZyme expression in G1 + PHA + AS vs. \nPHA + AS comparison') + # Plot title 
  geom_text_repel(max.overlaps = Inf) # To show all labels 




























############################################################################################################################################
############################################################################################################################################
############################################################################################################################################
# 10/07/2025
# Fixing my taxa associated with each group of interest (in this case CAZyme subgroups) because Malte said I was using the wrong file in our meeting
# Need to use tax_gtdb_metassembly.tsv that links contigs to taxa
# Need to use bowtie2_metaassembly.stranded2.counts.txt to link genes to contigs
library(readr)
library(tidyr)
library(dplyr)
library(tidyverse)

taxa_gtdb <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/tax_gtdb_metassembly.tsv")
contig_gene <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/input/bowtie2_metassembly.stranded2.counts.txt", header = FALSE)
contig_gene <- contig_gene[-c(1),]
names(contig_gene) <- contig_gene[1,]
contig_gene <- contig_gene[-c(1),]
names(contig_gene) <- gsub("mapping/bowtie.meta_t.","",names(contig_gene))
names(contig_gene) <- gsub(".sorted.bam","",names(contig_gene))


taxa_gtdb2 <- taxa_gtdb
taxa_gtdb2 <- taxa_gtdb2 %>%
  column_to_rownames("contig")
taxa_gtdb3 <- taxa_gtdb2 %>%
  unite("taxa", domain, phylum, class, order, family, genus, species, sep = " ", remove = FALSE)

taxa_gtdb3 <- taxa_gtdb3 %>%
  rownames_to_column("Contig")

# contig_gene <- contig_gene %>%
#   rename("Contig" = "Chr")
contig_gene <- contig_gene %>%
  rename("Chr" = "Contig")

contig_gene2 <- merge(contig_gene, taxa_gtdb3, all.x = TRUE)

bpcts <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/read_counts_gene_metassembly.tsv")
bpcts <- bpcts[order(match(bpcts$GeneID, contig_gene2$Geneid)), ]
all(bpcts[,7:26] == contig_gene2[,7:26])

# ok this confirmed that the bowtie2_metassembly.stranded2.counts.txt counts are the same as the read_counts_gene_metassembly.tsv, phew
# ok now repeat analysis using contig_gene2 instead of bpcts cus they have the same info, but contig_gene2 has the taxa info added on

library('DRnaSeq')
library(readr)

## Count matrix input
bpcoldata <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/bp_metadata_01132025.txt")
# bpcts <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/read_counts_gene_metassembly.tsv")
bpcts <- contig_gene2
tax_anno <- bpcts[,c(2,27:34)]
tax_anno <- tax_anno %>%
  rename("GeneID" = "Geneid") %>%
  as.data.frame()

bpcts <- bpcts %>%
  column_to_rownames("Geneid")
bpcts <- bpcts[,-c(1:5,26:33)]
head(bpcts)

bpcoldata <- bpcoldata %>%
  column_to_rownames("name")
head(bpcoldata)

## Examine the count matrix and column data to see if they are consistent in terms of sample order
head(bpcts, 2)
bpcoldata
bpcoldata2 <- bpcoldata[rownames(bpcoldata) %in% colnames(bpcts), ]

## Rearrange
all(rownames(bpcoldata2) %in% colnames(bpcts))
all(rownames(bpcoldata2) == colnames(bpcts))

bpcts <- bpcts[, rownames(bpcoldata2)]
all(rownames(bpcoldata2) == colnames(bpcts))
bpcts$GeneID <- rownames(bpcts)

# Input gene annotation info
gene_anno <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/gene_annotation_metassembly.tsv")

bpcts_geneanno <- merge(bpcts, gene_anno, by = "GeneID", all = TRUE)
bpcts_geneanno2 <- merge(bpcts_geneanno, tax_anno, by = "GeneID", all = TRUE)
# bpcts_geneanno2$KEGG_ko <- gsub("ko:","",bpcts_geneanno2$KEGG_ko)

library(dplyr)
library(stringr)
library(tidyr)

# dbcan annotations
aa <- bpcts_geneanno2 %>%
  dplyr::select(GeneID, everything())%>%
  mutate(dbcan_annotations=str_split(dbcan_annotations, "\\s*\\|\\s*")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:54,1)]

unique(bpcts_geneanno2$dbcan_annotations)
length(unique(bpcts_geneanno2$dbcan_annotations))
ab<-as.vector(as.matrix(aa[,grepl( "dbcan_annotations" , names(aa))]))
unique(ab)
length(unique(ab))

# make subset of aa data for a test grouping run
bb <- aa[is.na(aa$dbcan_annotations3) == FALSE,]
aasub <- bb[1:5,]
aasub$dbcan_annotations1[is.na(aasub$dbcan_annotations_b1)] <- "-"
aasubccs <- aasub[,c(25,43:46)]

uab <-unique(ab)
uab
uab <- uab[!is.na(uab)]
rn <- append(uab, "-")
rn
dbcantable_gen <- data.frame(matrix(0, nrow = length(rn), ncol = 20))
colnames(dbcantable_gen) <- colnames(aa[,1:20])
rownames(dbcantable_gen) <- rn

# Trial run on subset of the data
for (j in 1:ncol(dbcantable_gen)) {
  for (k in 1:length(names(aasub)[grepl( "dbcan_annotations" , names(aasub))])) { # Number of separate cazyme types broken up
    for (i in 1:nrow(aasub)) {
      if(is.na(aasub[i,paste("dbcan_annotations",k, sep = "")]) == FALSE) {
        dbcantable_gen[rn[rn %in% aasub[i,paste("dbcan_annotations",k, sep = "")]],j] = dbcantable_gen[rn[rn %in% aasub[i,paste("dbcan_annotations",k, sep = "")]],j] + aasub[i,j]
      }
      else {
        dbcantable_gen = dbcantable_gen
      }
    }
  }
}

# i made dbcantable_gen_raw_ttax to test on a subset of dbcantable_gen_raw since I didnt want to wipe the current data. It's the same
dbcantable_gen_raw_ttax <- dbcantable_gen_raw
dbcantable_gen_raw_ttax$taxaf <- NA
for (k in 1:length(aasub[,grepl("dbcan_annotations" , names(aasub))])) { # Number of separate CAZyme types broken up
  for (i in 1:nrow(aasub)) {
    if(aasub[i,paste("dbcan_annotations",k, sep = "")] %in% rn == TRUE ) {
      dbcantable_gen_raw_ttax[rn[rn %in% aasub[i,paste("dbcan_annotations",k, sep = "")]],"taxaf"] <- paste(dbcantable_gen_raw_ttax[rn[rn %in% aasub[i,paste("dbcan_annotations",k, sep = "")]],"taxaf"],"; ", aasub$taxa[i], sep = "")
    }
    else {
      dbcantable_gen_raw_ttax = dbcantable_gen_raw_ttax
    }
  }
}


## Now do it on the entire data
aa$dbcan_annotations1[is.na(aa$dbcan_annotations1)] <- "-"

ab<-as.vector(as.matrix(aa[,grepl( "dbcan_annotations" , names(aa))]))
uab <-unique(ab)
uab
uab <- uab[!is.na(uab)]
uab
rn <- uab
dbcantable_gen_raw <- data.frame(matrix(0, nrow = length(rn), ncol = 20))
colnames(dbcantable_gen_raw) <- colnames(aa[,1:20])
rownames(dbcantable_gen_raw) <- rn
length(rn)

for (j in 1:ncol(dbcantable_gen_raw)) {
  for (k in 1:length(names(aa)[grepl( "dbcan_annotations" , names(aa))])) { # Number of separate cazyme types broken up
    for (i in 1:nrow(aa)) {
      if(is.na(aa[i,paste("dbcan_annotations",k, sep = "")]) == FALSE) {
        dbcantable_gen_raw[rn[rn %in% aa[i,paste("dbcan_annotations",k, sep = "")]],j] <- as.numeric(dbcantable_gen_raw[rn[rn %in% aa[i,paste("dbcan_annotations",k, sep = "")]],j]) + as.numeric(aa[i,j])
      }
      else {
        dbcantable_gen_raw = dbcantable_gen_raw
      }
    }
  }
}

dbcantable_gen_raw$CAZID <- rownames(dbcantable_gen_raw)
write_tsv(dbcantable_gen_raw, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_raw_05192025.tsv') # I have checked this multiple times, even on 10/03/2025 and I trust that this is accurate
dbcantable_gen_rawr <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_raw_05192025.tsv')
dbcantable_gen_rawr <- as.data.frame(dbcantable_gen_rawr)
rownames(dbcantable_gen_rawr) <- dbcantable_gen_rawr$CAZID
drop <- "CAZID"
dbcantable_gen_rawr <- dbcantable_gen_rawr[,!(names(dbcantable_gen_rawr) %in% drop)]
dbcantable_gen_rawr <- dbcantable_gen_rawr[rownames(dbcantable_gen_raw),]
all(dbcantable_gen_rawr == dbcantable_gen_raw)

# for (k in 1:length(aa[,grepl("dbcan_annotations" , names(aa))])) { # Number of separate CAZyme types broken up
#   for (i in 1:nrow(aa)) {
#     if(aa[i,paste("dbcan_annotations",k, sep = "")] %in% rownames(dbcantable_gen_rawr) == TRUE ) {
#       dbcantable_gen_rawr[rn[rownames(dbcantable_gen_rawr) %in% aa[i,paste("dbcan_annotations",k, sep = "")]],"taxa"] <- paste(dbcantable_gen_rawr[rn[rownames(dbcantable_gen_rawr) %in% aa[i,paste("dbcan_annotations",k, sep = "")]],"taxa"],"; ", aa$eggNOG_OGs2[i], sep = "")
#     }
#     else {
#       dbcantable_gen_rawr = dbcantable_gen_rawr
#     }
#   }
# }

dbcantable_gen_raw$taxaf <- NA
for (k in 1:length(aa[,grepl("dbcan_annotations" , names(aa))])) { # Number of separate CAZyme types broken up
  for (i in 1:nrow(aa)) {
    if(aa[i,paste("dbcan_annotations",k, sep = "")] %in% rn == TRUE ) {
      dbcantable_gen_raw[rn[rn %in% aa[i,paste("dbcan_annotations",k, sep = "")]],"taxaf"] <- paste(dbcantable_gen_raw[rn[rn %in% aa[i,paste("dbcan_annotations",k, sep = "")]],"taxaf"],"; ", aa$taxa[i], sep = "")
    }
    else {
      dbcantable_gen_raw = dbcantable_gen_raw
    }
  }
}

# for (k in 1:length(aasub[,grepl("dbcan_annotations" , names(aasub))])) { # Number of separate CAZyme types broken up
#   for (i in 1:nrow(aasub)) {
#     if(aasub[i,paste("dbcan_annotations",k, sep = "")] %in% rn == TRUE ) {
#       dbcantable_gen_raw_ttax[rn[rn %in% aasub[i,paste("dbcan_annotations",k, sep = "")]],"taxaf"] <- paste(dbcantable_gen_raw_ttax[rn[rn %in% aasub[i,paste("dbcan_annotations",k, sep = "")]],"taxaf"],"; ", aasub$taxa[i], sep = "")
#     }
#     else {
#       dbcantable_gen_raw_ttax = dbcantable_gen_raw_ttax
#     }
#   }
# }
# 10/10/2025 add-in to save my updated dbcantable_gen_raw with the correct taxa
dbcantable_gen_raw$CAZID <- rownames(dbcantable_gen_raw)
write_tsv(dbcantable_gen_raw, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_raw_CORRECT_taxa_10102025.tsv')
dbcantable_gen_rawr <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_raw_CORRECT_taxa_10102025.tsv")
dbcantable_gen_rawr <- dbcantable_gen_rawr %>%
  column_to_rownames("CAZID")

# # removed "NA; " at beginning of taxaf
# dbcantable_gen_raw$taxaf <- gsub("NA; ","",dbcantable_gen_raw$taxaf)
# write_tsv(dbcantable_gen_raw, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_raw_CORRECT_taxa_noNAstart_10102025.tsv')
# 01/30/2026 restarting from here bc the gsub command is wrong
dbcantable_gen_rawr$taxaf <- substr(dbcantable_gen_rawr$taxaf, 5, nchar(dbcantable_gen_rawr$taxaf))
dbcantable_gen_rawrt <- dbcantable_gen_rawr

# dbcantable_gen_rawr$CAZID <- rownames(dbcantable_gen_rawr)
# write_tsv(dbcantable_gen_rawr, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_raw_taxa_05202025.tsv')
# dbcantable_gen_rawrt <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_raw_taxa_05202025.tsv')
# dbcantable_gen_rawrt <- dbcantable_gen_rawrt %>%
#   column_to_rownames("CAZID")

dbcantable_gen_rawrt <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_raw_CORRECT_taxa_noNAstart_NOnonCAZymes_10102025.tsv') # big nonCAZ category removed cus it slows it down
dbcantable_gen_rawrt <- dbcantable_gen_rawrt %>%
  column_to_rownames("CAZID")

# drop any outlier samples (based on PCA analysis)
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
dbcantable_gen_rawrt_filt <- dbcantable_gen_rawrt[,!(names(dbcantable_gen_rawrt) %in% drops)]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(dbcantable_gen_rawrt_filt),]

#######################################################################################################################################
#### CA + S3 + AS vs. CA + AS comparison

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "CA + AS" | bpcoldata3$group == "CA + S3 + AS",]
drop <- "-"

dbcantable_raw_filt_sub <- dbcantable_gen_rawrt_filt[,colnames(dbcantable_gen_rawrt_filt) %in% rownames(bpcoldata3_sub)]
dbcantable_raw_filt_sub_noNA <- dbcantable_raw_filt_sub[!(row.names(dbcantable_raw_filt_sub) %in% drop), ]
dbcantable_raw_taxa <- dbcantable_gen_rawrt_filt[!(row.names(dbcantable_gen_rawrt_filt) %in% drop),]

colnames(dbcantable_raw_filt_sub) == rownames(bpcoldata3_sub)
colnames(dbcantable_raw_filt_sub_noNA) == rownames(bpcoldata3_sub)

##############################################################################################################################
#### Now do for PHA + G1 + AS vs. PHA + AS comparison

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "PHA + G1 + AS" | bpcoldata3$group == "PHA + AS",]
drop <- "-"

dbcantable_raw_filt_sub <- dbcantable_gen_rawrt_filt[,colnames(dbcantable_gen_rawrt_filt) %in% rownames(bpcoldata3_sub)]
dbcantable_raw_filt_sub_noNA <- dbcantable_raw_filt_sub[!(row.names(dbcantable_raw_filt_sub) %in% drop), ]
dbcantable_raw_taxa <- dbcantable_gen_rawrt_filt[!(row.names(dbcantable_gen_rawrt_filt) %in% drop),]

colnames(dbcantable_raw_filt_sub) == rownames(bpcoldata3_sub)
colnames(dbcantable_raw_filt_sub_noNA) == rownames(bpcoldata3_sub)

#######################################################################################################################################
# ANCOMBC

library(phyloseq)
library(ANCOMBC)
library(tidyverse)
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))
library(dplyr)
library(ggplot2)

otumat = as.matrix(dbcantable_raw_filt_sub_noNA)

taxmat = matrix(sample(letters, 50, replace = TRUE), nrow = nrow(otumat), ncol = 8)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Taxa")
taxmat[,"Genus"] <- rownames(taxmat)
all(rownames(dbcantable_raw_taxa) == rownames(taxmat))
taxmat[,"Taxa"] <- dbcantable_raw_taxa[, "taxaf"]
taxmat

class(otumat)
class(taxmat)

library("phyloseq")
OTU = otu_table(otumat, taxa_are_rows = TRUE)
TAX = tax_table(taxmat)
OTU
TAX

sampledata <- sample_data(bpcoldata3_sub)
sampledata
physeq = phyloseq(OTU, TAX, sampledata)
physeq
physeq@sam_data

library(ANCOMBC)
out = ancombc(data = physeq, tax_level = "Genus", 
              formula = "group", 
              p_adj_method = "BH",
              conserve = TRUE,
              alpha = 0.01,
              verbose = TRUE)

res = out$res
res_global = out$res_global

## ANCOMBC primary result
# LFC
tab_lfc = res$lfc
col_name = c("CAZyme", "Intercept", "CAS3AS")
col_name = c("CAZyme", "Intercept", "PHAG1AS")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# SE
tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Test statistic
tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# P-values
tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Adjusted p-values
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Differentially abundant taxa
tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant CAZymes from the Primary Result")

cadiff <- tab_diff[tab_diff$`CA + S3 + AS` == TRUE,]
phadiff <- tab_diff[tab_diff$`PHA + G1 + AS` == TRUE,]

tab_lfc_filt <- tab_lfc[tab_lfc$CAZyme %in% cadiff$CAZyme,]
colnames(tab_lfc_filt) <- c("CAZyme", "Intercept LFC", "CA + S3 + AS LFC")
tab_lfc_filt <- tab_lfc[tab_lfc$CAZyme %in% phadiff$CAZyme,]
colnames(tab_lfc_filt) <- c("CAZyme", "Intercept LFC", "PHA + G1 + AS LFC")

tab_se_filt <- tab_se[tab_se$CAZyme %in% cadiff$CAZyme,]
colnames(tab_se_filt) <- c("CAZyme", "Intercept SE", "CA + S3 + AS SE")
tab_se_filt <- tab_se[tab_se$CAZyme %in% phadiff$CAZyme,]
colnames(tab_se_filt) <- c("CAZyme", "Intercept SE", "PHA + G1 + AS SE")

tab_w_filt <- tab_w[tab_w$CAZyme %in% cadiff$CAZyme,]
colnames(tab_w_filt) <- c("CAZyme", "Intercept W", "CA + S3 + AS W")
tab_w_filt <- tab_w[tab_w$CAZyme %in% phadiff$CAZyme,]
colnames(tab_w_filt) <- c("CAZyme", "Intercept W", "PHA + G1 + AS W")

tab_p_filt <- tab_p[tab_p$CAZyme %in% cadiff$CAZyme,]
colnames(tab_p_filt) <- c("CAZyme", "Intercept P", "CA + S3 + AS P")
tab_p_filt <- tab_p[tab_p$CAZyme %in% phadiff$CAZyme,]
colnames(tab_p_filt) <- c("CAZyme", "Intercept P", "PHA + G1 + AS P")

tab_q_filt <- tab_q[tab_q$CAZyme %in% cadiff$CAZyme,]
colnames(tab_q_filt) <- c("CAZyme", "Intercept Q", "CA + S3 + AS Q")
tab_q_filt <- tab_q[tab_q$CAZyme %in% phadiff$CAZyme,]
colnames(tab_q_filt) <- c("CAZyme", "Intercept Q", "PHA + G1 + AS Q")

tab_diff_filt <- tab_diff[tab_diff$CAZyme %in% cadiff$CAZyme,]
colnames(tab_diff_filt) <- c("CAZyme", "Intercept diff", "CA + S3 + AS diff")
tab_diff_filt <- tab_diff[tab_diff$CAZyme %in% phadiff$CAZyme,]
colnames(tab_diff_filt) <- c("CAZyme", "Intercept diff", "PHA + G1 + AS diff")

#put all data frames into list
df_list <- list(tab_lfc_filt, tab_se_filt, tab_w_filt, tab_p_filt, tab_q_filt, tab_diff_filt)

#merge all data frames in list
res_filt <- df_list %>% purrr::reduce(full_join, by='CAZyme')

# save list
#writexl::write_xlsx(res_filt, "")

### Visualization of differentially abundant taxa heatmap

sig_taxa = tab_diff %>%
  dplyr::filter(`CA + S3 + AS` == TRUE) %>%
  .$CAZyme
sig_taxa = tab_diff %>%
  dplyr::filter(`PHA + G1 + AS` == TRUE) %>%
  .$CAZyme

tab_lfc = res$lfc
col_name = c("CAZyme", "Intercept", "CA + S3 + AS vs. CA + AS")
col_name = c("CAZyme", "Intercept", "PHA + G1 + AS vs. PHA + AS")
colnames(tab_lfc) = col_name
tab_lfc <- tab_lfc[,-2]

df_CAS3AS = tab_lfc %>%
  filter(CAZyme %in% sig_taxa)

dbcantable_raw_taxa$taxa <- gsub("NA; ","",dbcantable_raw_taxa$taxa)
dbcantable_raw_taxa$CAZyme <- rownames(dbcantable_raw_taxa)
df_CAS3AS_fin <- merge(df_CAS3AS, dbcantable_raw_taxa[, c("taxa","CAZyme")], by = "CAZyme")
write_tsv(df_CAS3AS_fin, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/CAS3AS_CAZyme_subs_ancombc_05212025.tsv")
bbb <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/CAS3AS_CAZyme_subs_ancombc_05212025.tsv")

df_PHAG1AS = tab_lfc %>%
  filter(CAZyme %in% sig_taxa) 

dbcantable_raw_taxa$taxa <- gsub("NA; ","",dbcantable_raw_taxa$taxa)
dbcantable_raw_taxa$CAZyme <- rownames(dbcantable_raw_taxa)
df_PHAG1AS_fin <- merge(df_PHAG1AS, dbcantable_raw_taxa[, c("taxa","CAZyme")], by = "CAZyme")
write_tsv(df_PHAG1AS_fin, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/PHAG1AS_CAZyme_subs_ancombc_05212025.tsv")
bbb <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/PHAG1AS_CAZyme_subs_ancombc_05212025.tsv")

# 01/30/2026 updated to save padj lfc taxa in one data frame
tab_lfc <- tab_lfc[,c(1,3)]
tab_lfc <- tab_lfc %>%
  rename("CAS3AS" = "CAS3AS_lfc")
tab_q <- tab_q[,c(1,3)]
tab_q <- tab_q %>%
  rename("CAS3AS" = "CAS3AS_padj")
df_CAS3AS <- merge(tab_lfc, tab_q, by = "CAZyme")
dbcantable_raw_taxa <- dbcantable_raw_taxa %>%
  rownames_to_column("CAZyme")
taxtab <- dbcantable_raw_taxa[,c(1,20)]
df_CAS3AS2 <- merge(df_CAS3AS, taxtab, by = "CAZyme")
write_tsv(df_CAS3AS2, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/df_CAS3AS2_cazymelfcpadjtaxacorrect_01302026.tsv")

tab_lfc <- tab_lfc[,c(1,3)]
tab_lfc <- tab_lfc %>%
  rename("PHAG1AS" = "PHAG1AS_lfc")
tab_q <- tab_q[,c(1,3)]
tab_q <- tab_q %>%
  rename("PHAG1AS" = "PHAG1AS_padj")
df_PHAG1AS <- merge(tab_lfc, tab_q, by = "CAZyme")
dbcantable_raw_taxa <- dbcantable_raw_taxa %>%
  rownames_to_column("CAZyme")
taxtab <- dbPHAntable_raw_taxa[,c(1,20)]
df_PHAG1AS2 <- merge(df_PHAG1AS, taxtab, by = "CAZyme")
write_tsv(df_PHAG1AS2, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/df_PHAG1AS2_cazymelfcpadjtaxacorrect_01302026.tsv")


#####################################################################################################################################
# Make heatmap

df_heat = df_CAS3AS %>%
  pivot_longer(cols = -one_of("CAZyme"),
               names_to = "region", values_to = "value") %>%
  mutate(value = round(value, 2))
df_heat$CAZyme = factor(df_heat$CAZyme, levels = sort(sig_taxa))
# df_heatmod <- merge(df_heat, modtab,  by = "Module", all.x = TRUE)

df_heat = df_PHAG1AS %>%
  pivot_longer(cols = -one_of("CAZyme"),
               names_to = "region", values_to = "value") %>%
  mutate(value = round(value, 2))
df_heat$CAZyme = factor(df_heat$CAZyme, levels = sort(sig_taxa))
# df_heatmod <- merge(df_heat, modtab,  by = "Module", all.x = TRUE)

lo = floor(min(df_heat$value))
up = ceiling(max(df_heat$value))
mid = (lo + up)/2
# df_heat_filt <- df_heat[abs(df_heat$value) > 1,]
# df_heatmod_filt <- df_heatmod[abs(df_heatmod$value) > 1,]
# df_heatmod_filt$Label <- paste(df_heatmod_filt$Module, ": ", df_heatmod_filt$Function, sep = "")

# df_heatmod_filt2 <- df_heatmod_filt[,c(1:3,6)]
# write_xlsx(df_heatmod_filt2, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAS3AS_KEGGancombc_05052025.xlsx")
# 
# df_heatmod_filt2 <- df_heatmod_filt[,c(1:3,6)]
# write_xlsx(df_heatmod_filt2, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/PHAG1AS_KEGGancombc_05052025.xlsx")

p_heat = df_heat %>%
  ggplot(aes(x = region, y = CAZyme, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(region, CAZyme, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, size = 5,
       title = "Log fold changes for globally significant CAZymes") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 14))
p_heat





####################################################################
# 10/05/2025
# making volcano plots instead of heatmmaps to show differentially expressed CAZymes
# taken from https://biostatsquid.com/volcano-plots-r-tutorial/

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

# for CAS3AS vs CAAS comparison
# Import Ancombc results
df_CAS3AS2r <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/df_CAS3AS2_cazymelfcpadjtaxacorrect_01302026.tsv")
# df1 <- df_CAS3AS2r %>%
#   rename("CAS3AS_lfc" = "LFC", "CAS3AS_padj" = "padj")
df1 <- df_CAS3AS2r %>%
  rename("LFC" = "CAS3AS_lfc", "padj" = "CAS3AS_padj")

# tab_lfc2 <- tab_lfc %>%
#   rename("CA + S3 + AS vs. CA + AS" = "LFC")
# 
# tab_q2 <- tab_q %>%
#   rename("CA + S3 + AS" = "padj")
# 
# df <- merge(tab_lfc2, tab_q2, by = "CAZyme")

ggplot(data = df1, aes(x = LFC, y = -log10(padj))) +
  geom_point()

# Add threshold lines
ggplot(data = df1, aes(x = LFC, y = -log10(padj))) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point() 

# Biostatsquid theme
theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

ggplot(data = df1, aes(x = LFC, y = -log10(padj))) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point() 

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
df1$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df1$diffexpressed[df1$LFC > 1 & df1$padj < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df1$diffexpressed[df1$LFC < -1 & df1$padj < 0.05] <- "DOWN"

# Explore a bit
head(df1[order(df1$padj) & df1$diffexpressed == 'DOWN', ])

ggplot(data = df1, aes(x = LFC, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  labs(color='Differential Expression') 


# Edit axis labels and limits
ggplot(data = df1, aes(x = LFC, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 50), xlim = c(-10, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-10, 10, 2)) # to customise the breaks in the x axis

# Note. with coord_cartesian() even if we have genes with p-values or log2FC ourside our limits, they will still be plotted.

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
dlup <- df1[df1$diffexpressed == "UP",]
df1$delabel <- ifelse(df1$CAZyme %in% as.matrix(head(dlup[order(dlup$padj), "CAZyme"], 10)), df1$CAZyme, NA)

dldown <- df1[df1$diffexpressed == "DOWN",]
df1$delabel <- ifelse(df1$CAZyme %in% as.matrix(head(dldown[order(dldown$padj), "CAZyme"], 10)), df1$CAZyme, df1$delabel)

ggplot(data = df1, aes(x = LFC, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 35), xlim = c(-6, 7)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-6, 7, 2)) + # to customise the breaks in the x axis
  geom_label_repel(size = 4, show.legend = FALSE, min.segment.length = unit(0, 'lines')) +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), legend.title = element_blank(), legend.position = "bottom")

## 10/05/2025 ^ ok so for the graph above I am labeling the 10 most significantly expressed (lowest padj values) upregulated CAZyme subgroup and 
# the 10 most significantly expressed (lowest padj values) downregulated CAZyme subgroup and 
# This also includes an abs(1) LFC cutoff, which excludes some of the lowest p value CAZ subgroups, so maybe change that? Talk to Magda
# see if we want to keep the LFC cutoff or get rid of it

# for G1PHAAS vs PHAAS comparison
# Import ancombc results
df_PHAG1AS2r <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/df_PHAG1AS2_cazymelfcpadjtaxacorrect_01302026.tsv")
# df <- df_PHAG1AS2r %>%
#   rename("PHAG1AS_lfc" = "LFC", "PHAG1AS_padj" = "padj")
df <- df_PHAG1AS2r %>%
  rename("LFC" = "PHAG1AS_lfc", "padj" = "PHAG1AS_padj")
  
# tab_lfc2 <- tab_lfc %>%
#   rename("PHA + G1 + AS vs. PHA + AS" = "LFC")
# 
# tab_q2 <- tab_q %>%
#   rename("PHA + G1 + AS" = "padj")
# 
# df <- merge(tab_lfc2, tab_q2, by = "CAZyme")

theme_set(theme_classic(base_size = 20) +
            theme(
              axis.title.y = element_text(face = "bold", margin = margin(0,20,0,0), size = rel(1.1), color = 'black'),
              axis.title.x = element_text(hjust = 0.5, face = "bold", margin = margin(20,0,0,0), size = rel(1.1), color = 'black'),
              plot.title = element_text(hjust = 0.5)
            ))

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
df$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$LFC > 1 & df$padj < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$LFC < -1 & df$padj < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
dlup <- df[df$diffexpressed == "UP",]
df$delabel <- ifelse(df$CAZyme %in% as.matrix(head(dlup[order(dlup$padj), "CAZyme"], 10)), df$CAZyme, NA)

dldown <- df[df$diffexpressed == "DOWN",]
df$delabel <- ifelse(df$CAZyme %in% as.matrix(head(dldown[order(dldown$padj), "CAZyme"], 10)), df$CAZyme, df$delabel)

ggplot(data = df, aes(x = LFC, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 10), xlim = c(-5, 4)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-5, 4, 2)) + # to customise the breaks in the x axis
  geom_label_repel(size = 4, show.legend = FALSE, min.segment.length = unit(0, 'lines')) +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), legend.title = element_blank(), legend.position = "bottom")


## 02/02/2026 making smaller volcano plot and heatmap instead of giant volcano plot - maybe better for plot?
# for S3 CA AS vs CA AS
library(ggrepel)

np1 <- ggplot(data = df1, aes(x = LFC, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 40), xlim = c(-6, 7)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(title = "CA + S3 + AS vs. CA + AS", color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-6, 7, 2)) + # to customise the breaks in the x axis
  geom_label_repel(size = 3, show.legend = FALSE, min.segment.length = unit(0, 'lines'), max.overlaps = Inf) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 12), 
        legend.title = element_blank(), legend.text = element_text(size = 12), legend.position = "bottom", plot.margin = margin(5, 5, 5, 5))
np1

# for G1 PHA AS vs PHA AS
np2 <- ggplot(data = df, aes(x = LFC, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 10), xlim = c(-5, 5)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(title = "PHA + G1 + AS vs. PHA + AS", color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-5, 5, 2)) + # to customise the breaks in the x axis 
  scale_y_continuous(breaks = seq(-0, 10, 2)) + # to customise the breaks in the x axis
  geom_label_repel(size = 3, show.legend = FALSE, min.segment.length = unit(0, 'lines'), max.overlaps = Inf) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 12), 
        legend.title = element_blank(), legend.text = element_text(size = 12), legend.position = "bottom", plot.margin = margin(5, 5, 5, 5))
np2

library(ggpubr)
np3 <- ggarrange(np1, np2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
np3



# 05/06/2026 make volcano plots for Maaslin2 output instead of ancombc1

taxa_gtdb <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/tax_gtdb_metassembly.tsv")
contig_gene <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/input/bowtie2_metassembly.stranded2.counts.txt", header = FALSE)
contig_gene <- contig_gene[-c(1),]
names(contig_gene) <- contig_gene[1,]
contig_gene <- contig_gene[-c(1),]
names(contig_gene) <- gsub("mapping/bowtie.meta_t.","",names(contig_gene))
names(contig_gene) <- gsub(".sorted.bam","",names(contig_gene))


taxa_gtdb2 <- taxa_gtdb
taxa_gtdb2 <- taxa_gtdb2 %>%
  column_to_rownames("contig")
taxa_gtdb3 <- taxa_gtdb2 %>%
  unite("taxa", domain, phylum, class, order, family, genus, species, sep = " ", remove = FALSE)

taxa_gtdb3 <- taxa_gtdb3 %>%
  rownames_to_column("Contig")

contig_gene <- contig_gene %>%
  rename("Contig" = "Chr")

contig_gene2 <- merge(contig_gene, taxa_gtdb3, all.x = TRUE)

bpcts <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/read_counts_gene_metassembly.tsv")
bpcts <- bpcts[order(match(bpcts$GeneID, contig_gene2$Geneid)), ]
all(bpcts[,7:26] == contig_gene2[,7:26])

# ok this confirmed that the bowtie2_metassembly.stranded2.counts.txt counts are the same as the read_counts_gene_metassembly.tsv, phew
# ok now repeat analysis using contig_gene2 instead of bpcts cus they have the same info, but contig_gene2 has the taxa info added on

library('DRnaSeq')
library(readr)

## Count matrix input
bpcoldata <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/bp_metadata_01132025.txt")
# bpcts <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/read_counts_gene_metassembly.tsv")
bpcts <- contig_gene2
tax_anno <- bpcts[,c(2,27:34)]
tax_anno <- tax_anno %>%
  rename("GeneID" = "Geneid") %>%
  as.data.frame()

bpcts <- bpcts %>%
  column_to_rownames("Geneid")
bpcts <- bpcts[,-c(1:5,26:33)]
head(bpcts)

bpcoldata <- bpcoldata %>%
  column_to_rownames("name")
head(bpcoldata)

## Examine the count matrix and column data to see if they are consistent in terms of sample order
head(bpcts, 2)
bpcoldata
bpcoldata2 <- bpcoldata[rownames(bpcoldata) %in% colnames(bpcts), ]

## Rearrange
all(rownames(bpcoldata2) %in% colnames(bpcts))
all(rownames(bpcoldata2) == colnames(bpcts))

bpcts <- bpcts[, rownames(bpcoldata2)]
all(rownames(bpcoldata2) == colnames(bpcts))
bpcts$GeneID <- rownames(bpcts)

# Input gene annotation info
gene_anno <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/gene_annotation_metassembly.tsv")

bpcts_geneanno <- merge(bpcts, gene_anno, by = "GeneID", all = TRUE)
bpcts_geneanno2 <- merge(bpcts_geneanno, tax_anno, by = "GeneID", all = TRUE)

dbcantable_gen_rawrt <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_raw_CORRECT_taxa_noNAstart_NOnonCAZymes_10102025.tsv') # big nonCAZ category removed cus it slows it down
dbcantable_gen_rawrt <- dbcantable_gen_rawrt %>%
  column_to_rownames("CAZID")

# drop any outlier samples (based on PCA analysis)
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC", "taxaf")
dbcantable_gen_rawrt_filt <- dbcantable_gen_rawrt[,!(names(dbcantable_gen_rawrt) %in% drops)]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(dbcantable_gen_rawrt_filt),]

# Run Maaslin2 instead of ANCOMBC
library(Maaslin2)
library(dplyr)
library(readr)

bpcoldata3$group <- gsub(" + ", "_", bpcoldata3$group, fixed = TRUE)
bpcoldata3$name <- row.names(bpcoldata3)

# Transpose so rows = samples, columns = features
dbcantable_gen_rawrt_filt_t <- as.data.frame(t(dbcantable_gen_rawrt_filt))
all(rownames(dbcantable_gen_rawrt_filt_t) == rownames(bpcoldata3))

# --- Comparison 1: CA_S3_AS vs CA_AS ---
meta1 <- bpcoldata3 |> filter(group %in% c("CA_S3_AS", "CA_AS"))
data1 <- dbcantable_gen_rawrt_filt_t[rownames(meta1), ]

fit1 <- Maaslin2(
  input_data      = data1,
  input_metadata  = meta1,
  output          = "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/CAZ_module_maaslin2_CA_05062026",
  fixed_effects   = "group",
  reference       = "group,CA_AS",
  normalization   = "TSS",
  transform       = "LOG",
  analysis_method = "LM",
  min_prevalence  = 0.1,
  min_abundance   = 0.0,
  cores           = 1
)

# --- Comparison 2: PHA_G1_AS vs PHA_AS ---
meta2 <- bpcoldata3 |> filter(group %in% c("PHA_G1_AS", "PHA_AS"))
data2 <- dbcantable_gen_rawrt_filt_t[rownames(meta2), ]

fit2 <- Maaslin2(
  input_data      = data2,
  input_metadata  = meta2,
  output          = "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/CAZ_module_maaslin2_PHA_05062026",
  fixed_effects   = "group",
  reference       = "group,PHA_AS",
  normalization   = "TSS",
  transform       = "LOG",
  analysis_method = "LM",
  min_prevalence  = 0.1,
  min_abundance   = 0.0,
  cores           = 1
)

# # --- Comparison 3: PHA_G1_AS vs PHA_AS ---
# bpcoldata3$condition_cba <- gsub(" + ", "_", bpcoldata3$condition_cba, fixed = TRUE)
# bpcoldata3$condition_cba <- gsub(" ", "_", bpcoldata3$condition_cba)
# 
# meta3 <- bpcoldata3 |> filter(condition_cba %in% c("PHA_G1_AS_t2", "PHA_G1_AS_t1"))
# data3 <- dbcantable_gen_rawrt_filt_t[rownames(meta3), ]
# 
# fit3 <- Maaslin2(
#   input_data      = data3,
#   input_metadata  = meta3,
#   output          = "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/KEGG_module_maaslin2_PHAG1ASt2t1_05052026",
#   fixed_effects   = "condition_cba",
#   reference       = "condition_cba,PHA_G1_AS_t1",
#   normalization   = "TSS",
#   transform       = "LOG",
#   analysis_method = "LM",
#   min_prevalence  = 0.1,
#   min_abundance   = 0.0,
#   cores           = 1
# )

# make heatmap
library(dplyr)
library(pheatmap)

# ── Extract MaAsLin2 results ──────────────────────────────────────────────────
res_CA  <- fit1$results
res_PHA <- fit2$results
# res_PHAb <- fit3$results

# ── Build LFC and q-value matrices ───────────────────────────────────────────
all_cogs <- union(res_CA$feature, res_PHA$feature)
# all_cogs <- union(res_CA$feature, union(res_PHA$feature, res_PHAb$feature))

lfc_CA  <- setNames(res_CA$coef,  res_CA$feature)
lfc_PHA <- setNames(res_PHA$coef, res_PHA$feature)
# lfc_PHAb <- setNames(res_PHAb$coef, res_PHAb$feature)
q_CA    <- setNames(res_CA$qval,  res_CA$feature)
q_PHA   <- setNames(res_PHA$qval, res_PHA$feature)
# q_PHAb   <- setNames(res_PHAb$qval, res_PHAb$feature)

plot_mat <- data.frame(
  CA  = lfc_CA[all_cogs],
  PHA = lfc_PHA[all_cogs],
  # PHAb = lfc_PHAb[all_cogs],
  row.names = all_cogs
)
plot_mat[is.na(plot_mat)] <- 0

padj_mat <- data.frame(
  CA  = q_CA[all_cogs],
  PHA = q_PHA[all_cogs],
  # PHAb = q_PHAb[all_cogs],
  row.names = all_cogs
)
padj_mat[is.na(padj_mat)] <- 1

combined_mat <- data.frame(
  lfc_CA  = lfc_CA[all_cogs],
  q_CA    = q_CA[all_cogs],
  lfc_PHA = lfc_PHA[all_cogs],
  q_PHA   = q_PHA[all_cogs],
  # lfc_PHAb = lfc_PHAb[all_cogs],
  # q_PHAb   = q_PHAb[all_cogs],
  row.names = all_cogs
)

combined_mat$lfc_CA[is.na(combined_mat$lfc_CA)]   <- 0
combined_mat$lfc_PHA[is.na(combined_mat$lfc_PHA)] <- 0
# combined_mat$lfc_PHAb[is.na(combined_mat$lfc_PHAb)] <- 0
combined_mat$q_CA[is.na(combined_mat$q_CA)]   <- 1
combined_mat$q_PHA[is.na(combined_mat$q_PHA)] <- 1
# combined_mat$q_PHAb[is.na(combined_mat$q_PHAb)] <- 1
combined_mat$CAZyme <- rownames(combined_mat)

write_tsv(combined_mat, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/CAZ_maaslin2_CAPHA_coeffpadj_05062026.tsv")
# write_tsv(combined_mat, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/KEGG_module_maaslin2_CAPHAPHAb_coeffpadj_05052026.tsv")
combined_matr <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/CAZ_maaslin2_CAPHA_coeffpadj_05062026.tsv")
# combined_matrb <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/KEGG_module_maaslin2_CAPHAPHAb_coeffpadj_05052026.tsv")


library(ggrepel)
df <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/CAZ_maaslin2_CAPHA_coeffpadj_05062026.tsv")

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
df$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$lfc_CA > 0 & df$q_CA < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$lfc_CA < 0 & df$q_CA < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
dlup <- df[df$diffexpressed == "UP",]
df$delabel <- ifelse(df$CAZyme %in% as.matrix(head(dlup[order(dlup$q_CA), "CAZyme"], 10)), df$CAZyme, NA)

dldown <- df[df$diffexpressed == "DOWN",]
df$delabel <- ifelse(df$CAZyme %in% as.matrix(head(dldown[order(dldown$q_CA), "CAZyme"], 10)), df$CAZyme, df$delabel)

np1 <- ggplot(data = df, aes(x = lfc_CA, y = -log10(q_CA), col = diffexpressed, label = delabel)) +
  # geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 4), xlim = c(-6, 6)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(title = "CA + S3 + AS vs. CA + AS", color = 'Differential Expression', #legend_title, 
       x = expression("Coefficient"), y = expression("-log"[10]*"q-value")) + 
  scale_x_continuous(breaks = seq(-6, 6, 2)) + # to customise the breaks in the x axis
  geom_label_repel(size = 3, show.legend = FALSE, min.segment.length = unit(0, 'lines'), max.overlaps = Inf) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 12), 
        legend.title = element_blank(), legend.text = element_text(size = 12), legend.position = "bottom", plot.margin = margin(5, 5, 5, 5))
np1

# for G1 PHA AS vs PHA AS

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
df$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$lfc_PHA > 0 & df$q_PHA < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$lfc_PHA < 0 & df$q_PHA < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in PHAse they are not)
dlup <- df[df$diffexpressed == "UP",]
df$delabel <- ifelse(df$CAZyme %in% as.matrix(head(dlup[order(dlup$q_PHA), "CAZyme"], 20)), df$CAZyme, NA)

dldown <- df[df$diffexpressed == "DOWN",]
df$delabel <- ifelse(df$CAZyme %in% as.matrix(head(dldown[order(dldown$q_PHA), "CAZyme"], 10)), df$CAZyme, df$delabel)

np2 <- ggplot(data = df, aes(x = lfc_PHA, y = -log10(q_PHA), col = diffexpressed, label = delabel)) +
  # geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 4), xlim = c(-6, 6)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(title = "PHA + G1 + AS vs. PHA + AS", color = 'Differential Expression', #legend_title, 
       x = expression("Coefficient"), y = expression("-log"[10]*"q-value")) + 
  scale_x_continuous(breaks = seq(-6, 6, 2)) + # to customise the breaks in the x axis 
  geom_label_repel(size = 3, show.legend = FALSE, min.segment.length = unit(0, 'lines'), max.overlaps = Inf) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 12), 
        legend.title = element_blank(), legend.text = element_text(size = 12), legend.position = "bottom", plot.margin = margin(5, 5, 5, 5))
np2

library(ggpubr)
np3 <- ggarrange(np1, np2, ncol=2, nrow=1, common.legend = TRUE, legend="bottom", labels = c("A)", "B)"))
np3
ggsave(plot = np3, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Chord_plots/CAZ_volcfig_05062026.pdf", width = 12, height = 6, units = "in", dpi = 300, bg = "white")















# just quick check rerunning ancombc1 to make sure i didn;t change anything except for the maaslin2 run
library(ANCOMBC)
library(phyloseq)
library(dplyr)
library(readr)

counts_t <- dbcantable_gen_rawrt_filt_t

# ── Helper: build phyloseq object ─────────────────────────────────────────────
make_phyloseq <- function(counts_t, metadata) {
  otu  <- otu_table(t(counts_t), taxa_are_rows = TRUE)
  samp <- sample_data(metadata)
  phyloseq(otu, samp)
}

# ── Comparison 1: CA_S3_AS vs CA_AS ──────────────────────────────────────────
meta1 <- bpcoldata3 |> filter(group %in% c("CA_S3_AS", "CA_AS"))
meta1$group <- factor(meta1$group, levels = c("CA_AS", "CA_S3_AS"))
data1 <- dbcantable_gen_rawrt_filt_t[rownames(meta1), ]
ps1 <- make_phyloseq(data1, meta1)

fit1 <- ancombc(
  data             = ps1,
  formula          = "group",
  p_adj_method     = "BH",
  prv_cut          = 0.10,
  lib_cut          = 0,
  group            = "group",
  struc_zero       = FALSE,
  neg_lb           = FALSE,
  tol              = 1e-05,
  max_iter         = 100,
  conserve         = TRUE,
  alpha            = 0.05,
  global           = FALSE
)

# ── Comparison 2: PHA_G1_AS vs PHA_AS ────────────────────────────────────────
meta2 <- bpcoldata3 |> filter(group %in% c("PHA_G1_AS", "PHA_AS"))
meta2$group <- factor(meta2$group, levels = c("PHA_AS", "PHA_G1_AS"))
data2 <- dbcantable_gen_rawrt_filt_t[rownames(meta2), ]
ps2 <- make_phyloseq(data2, meta2)

fit2 <- ancombc(
  data             = ps2,
  formula          = "group",
  p_adj_method     = "BH",
  prv_cut          = 0.10,
  lib_cut          = 0,
  group            = "group",
  struc_zero       = FALSE,
  neg_lb           = FALSE,
  tol              = 1e-05,
  max_iter         = 100,
  conserve         = TRUE,
  alpha            = 0.05,
  global           = FALSE
)

# ── Extract results ───────────────────────────────────────────────────────────
res_CA <- data.frame(
  feature = fit1$res$lfc$taxon,
  coef    = fit1$res$lfc$groupCA_S3_AS,
  qval    = fit1$res$q_val$groupCA_S3_AS
)

res_PHA <- data.frame(
  feature = fit2$res$lfc$taxon,
  coef    = fit2$res$lfc$groupPHA_G1_AS,
  qval    = fit2$res$q_val$groupPHA_G1_AS
)

# ── Build LFC and q-value matrices ────────────────────────────────────────────
all_cogs <- union(res_CA$feature, res_PHA$feature)

lfc_CA  <- setNames(res_CA$coef, res_CA$feature)
lfc_PHA <- setNames(res_PHA$coef, res_PHA$feature)
q_CA    <- setNames(res_CA$qval, res_CA$feature)
q_PHA   <- setNames(res_PHA$qval, res_PHA$feature)

combined_mat <- data.frame(
  lfc_CA  = lfc_CA[all_cogs],
  q_CA    = q_CA[all_cogs],
  lfc_PHA = lfc_PHA[all_cogs],
  q_PHA   = q_PHA[all_cogs],
  row.names = all_cogs
)

combined_mat$lfc_CA[is.na(combined_mat$lfc_CA)]   <- 0
combined_mat$lfc_PHA[is.na(combined_mat$lfc_PHA)] <- 0
combined_mat$q_CA[is.na(combined_mat$q_CA)]       <- 1
combined_mat$q_PHA[is.na(combined_mat$q_PHA)]     <- 1
combined_mat$CAZyme <- rownames(combined_mat)

library(ggrepel)
df <- combined_mat

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
df$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$lfc_CA > 1 & df$q_CA < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$lfc_CA < -1 & df$q_CA < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in case they are not)
dlup <- df[df$diffexpressed == "UP",]
df$delabel <- ifelse(df$CAZyme %in% as.matrix(head(dlup[order(dlup$q_CA), "CAZyme"], 10)), df$CAZyme, NA)

dldown <- df[df$diffexpressed == "DOWN",]
df$delabel <- ifelse(df$CAZyme %in% as.matrix(head(dldown[order(dldown$q_CA), "CAZyme"], 10)), df$CAZyme, df$delabel)

np1 <- ggplot(data = df, aes(x = lfc_CA, y = -log10(q_CA), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 40), xlim = c(-6, 7)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(title = "CA + S3 + AS vs. CA + AS", color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-6, 7, 2)) + # to customise the breaks in the x axis
  geom_label_repel(size = 3, show.legend = FALSE, min.segment.length = unit(0, 'lines'), max.overlaps = Inf) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 12), 
        legend.title = element_blank(), legend.text = element_text(size = 12), legend.position = "bottom", plot.margin = margin(5, 5, 5, 5))
np1

# for G1 PHA AS vs PHA AS

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
df$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$lfc_PHA > 1 & df$q_PHA < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
df$diffexpressed[df$lfc_PHA < -1 & df$q_PHA < 0.05] <- "DOWN"

# Create a new column "delabel" to de, that will contain the name of the top 30 differentially expressed genes (NA in PHAse they are not)
dlup <- df[df$diffexpressed == "UP",]
df$delabel <- ifelse(df$CAZyme %in% as.matrix(head(dlup[order(dlup$q_PHA), "CAZyme"], 22)), df$CAZyme, NA)

dldown <- df[df$diffexpressed == "DOWN",]
df$delabel <- ifelse(df$CAZyme %in% as.matrix(head(dldown[order(dldown$q_PHA), "CAZyme"], 22)), df$CAZyme, df$delabel)

np2 <- ggplot(data = df, aes(x = lfc_PHA, y = -log10(q_PHA), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 10), xlim = c(-6, 6)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(title = "PHA + G1 + AS vs. PHA + AS", color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-6, 6, 2)) + # to customise the breaks in the x axis 
  geom_label_repel(size = 3, show.legend = FALSE, min.segment.length = unit(0, 'lines'), max.overlaps = Inf) +
  theme_classic() +
  theme(plot.title = element_text(hjust = 0.5, face = "bold", size = 14), axis.text = element_text(size = 12), axis.title = element_text(size = 12), 
        legend.title = element_blank(), legend.text = element_text(size = 12), legend.position = "bottom", plot.margin = margin(5, 5, 5, 5))
np2













# make heatmap of CAZymes
df_CAS3AS2r <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/df_CAS3AS2_cazymelfcpadjtaxacorrect_01302026.tsv") %>%
  rename("taxaf" = "taxaCAS3") %>%
  arrange(CAS3AS_padj)
df_CAS3AS2rf <- df_CAS3AS2r[abs(df_CAS3AS2r$CAS3AS_lfc) > 1 & df_CAS3AS2r$CAS3AS_padj < 0.05,] %>%
  arrange(CAS3AS_padj)

df_PHAG1AS2r <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/df_PHAG1AS2_cazymelfcpadjtaxacorrect_01302026.tsv") %>%
  rename("taxaf" = "taxaPHAG1") %>%
  arrange(PHAG1AS_padj)
df_PHAG1AS2rf <- df_PHAG1AS2r[abs(df_PHAG1AS2r$PHAG1AS_lfc) > 1 & df_PHAG1AS2r$PHAG1AS_padj < 0.05,] %>%
  arrange(PHAG1AS_padj)

# Finding the CAZsubs where at least either CAS3AS or PHAG1AS has lfc > 1 and padj < 0.05
# Make the padj < 0 ones have lfc = 0 so they don't show expression in heatmap
df_figg <- merge(df_CAS3AS2r, df_PHAG1AS2r, by = "CAZyme", all = TRUE)
# lala <- df_figg[which(df_figg$CAS3AS_padj < 0.05 | df_figg$PHAG1AS_padj < 0.05),]
# lalaf <- lala[which(abs(lala$CAS3AS_lfc) > 1 | abs(lala$PHAG1AS_lfc) > 1),]
df_figgf <- df_figg[which((abs(df_figg$CAS3AS_lfc) > 2 & df_figg$CAS3AS_padj < 0.05) | (abs(df_figg$PHAG1AS_lfc) > 2 & df_figg$PHAG1AS_padj < 0.05)),]
df_figgf <- df_figg[which(abs(df_figg$CAS3AS_lfc) > 2 & df_figg$CAS3AS_padj < 0.05),]

df_figgf2 <- df_figgf
df_figgf2[is.na(df_figgf2$PHAG1AS_padj) == TRUE, "PHAG1AS_padj"] <- 100

for (j in c(2,5)) {
  for (i in 1:nrow(df_figgf2)) {
    if (df_figgf2[i,j+1] < 0.05) {
      df_figgf2[i,j] = df_figgf2[i,j]
    } else {
      df_figgf2[i,j] = 0
    }
  }
}


library(tidyverse)

newnames <- lapply(
  c("CA + S3 + AS\nvs. CA + AS", "PHA + G1 + AS\nvs. PHA + AS"),
  function(x) bquote(.(x)))
newnames

library(pheatmap)
library(grid)

draw_colnames_45 <- function (coln, gaps, ...) {
  coord = pheatmap:::find_coordinates(length(coln), gaps)
  x = coord$coord - 0.5 * coord$size
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 3, hjust = .5, rot = 360, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))
# rm(draw_colnames_45)

my_palette <- colorRampPalette(c("blue", "black", "red"))(n = 299)

# pheatmap(as.matrix(df_figb), cellheight=10, cluster_cols=FALSE,
#          fontsize = 10, show_rownames=TRUE, fontsize_row=8, fontsize_col=8,
#          labels_col = as.expression(collab), labels_row = as.expression(rownames(df_figb)), border_color=NA, 
#          annotation_row = annotation_row, annotation_colors = ann_colors, breaks=seq(-10, 10, length.out=300), angle_col = 0)

ph <- pheatmap(as.matrix(df_figgf2[,c(2,5)]), col=my_palette, cellwidth = 80, cellheight=10, cluster_cols=FALSE,
               fontsize = 10, show_rownames=TRUE,
               labels_col = as.expression(newnames), 
               labels_row = df_figgf2$CAZyme,
               annotation_names_row = FALSE )
ph