# Colleen Ahern
# 02/10/2025

# Continuation of BP_MT_transcript_abund_grouped_02062025 script
# I want to normalize my raw gene counts and add together KEGG groups instead of COGs
# Continue using TPM normalization

# do cazyme type later

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

# let's now do the COG annotations
aa <- bptpm_geneanno %>%
  select(GeneID, everything())%>%
  mutate(KEGG_Module=str_split(KEGG_Module, ",")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:60,1)]
#rownames(aa) <- 1:nrow(aa)

unique(bptpm_geneanno$KEGG_Module)
ab<-as.vector(as.matrix(aa[,grepl( "KEGG_Module" , names(aa))]))
unique(ab)
length(unique(ab))

# make subset of aa data
aasub <- as.data.frame(aa[1:11,])
aasub$KEGG_Module1[is.na(aasub$KEGG_Module1)] <- "-"
aasubccs <- aasub[,c(33:50)]

uab <-unique(ab)
uab <- uab[!is.na(uab)]
rn <- uab
keggmodtable_gen <- data.frame(matrix(0, nrow = length(unique(uab)), ncol = 20))
colnames(keggmodtable_gen) <- colnames(aa[,1:20])
rownames(keggmodtable_gen) <- uab

# Trial run on subset of the data
for (j in 1:ncol(keggmodtable_gen)) {
  for (k in 1:length(aasub[,grepl( "KEGG_Module" , names(aasub))])) { # Number of separate KEGG modules broken up
    for (i in 1:nrow(aasub)) {
      if(is.na(aasub[i,paste("KEGG_Module",k, sep = "")]) == FALSE) {
        keggmodtable_gen[rn[rownames(keggmodtable_gen) %in% aasub[i,paste("KEGG_Module",k, sep = "")]],j] = keggmodtable_gen[rn[rownames(keggmodtable_gen) %in% aasub[i,paste("KEGG_Module",k, sep = "")]],j] + aasub[i,j]
      }
      else {
        keggmodtable_gen = keggmodtable_gen
      }
    }
  }
}


## Now do it on the entire data
aa$KEGG_Module1[is.na(aa$KEGG_Module1)] <- "-"

uab <-unique(ab)
uab <- uab[!is.na(uab)]
rn <- uab
rn
keggmodtable_gen <- data.frame(matrix(0, nrow = length(unique(uab)), ncol = 20))
colnames(keggmodtable_gen) <- colnames(aa[,1:20])
rownames(keggmodtable_gen) <- uab

for (j in 1:ncol(keggmodtable_gen)) {
  for (k in 1:length(aa[,grepl( "KEGG_Module" , names(aa))])) { # Number of separate KEGG modules broken up
    for (i in 1:nrow(aa)) {
      if(is.na(aa[i,paste("KEGG_Module",k, sep = "")]) == FALSE) {
        keggmodtable_gen[rn[rownames(keggmodtable_gen) %in% aa[i,paste("KEGG_Module",k, sep = "")]],j] = keggmodtable_gen[rn[rownames(keggmodtable_gen) %in% aa[i,paste("KEGG_Module",k, sep = "")]],j] + aa[i,j]
      }
      else {
        keggmodtable_gen = keggmodtable_gen
      }
    }
  }
}

write_tsv(keggmodtable_gen, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggmodtable_gen_02102025.tsv')
keggmodtable_gen <- read.csv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggmodtable_gen_02102025.tsv', header = TRUE, sep = "")

# Take average of replicates
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
keggmodtable_gen_filt <- keggmodtable_gen[,!(names(keggmodtable_gen) %in% drops)]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(keggmodtable_gen_filt),]
keggmodtable_gen_avg <- data.frame(matrix(0, nrow = length(unique(uab)), ncol = 8))
colnames(keggmodtable_gen_avg) <- (unique(bpcoldata3$condition_cba))
rownames(keggmodtable_gen_avg) <- rownames(keggmodtable_gen_filt)

keggmodtable_gen_avg[,'CA + AS t1'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + AS t1',])]))
keggmodtable_gen_avg[,'CA + S3 + AS t1'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + S3 + AS t1',])]))
keggmodtable_gen_avg[,'CA + AS t2'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + AS t2',])]))
keggmodtable_gen_avg[,'CA + S3 + AS t2'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + S3 + AS t2',])]))

keggmodtable_gen_avg[,'PHA + AS t1'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + AS t1',])]))
keggmodtable_gen_avg[,'PHA + G1 + AS t1'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + G1 + AS t1',])]))
keggmodtable_gen_avg[,'PHA + AS t2'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + AS t2',])]))
keggmodtable_gen_avg[,'PHA + G1 + AS t2'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + G1 + AS t2',])]))

########################################################################################################################
## Make heatmap
## BP: Import metadata file

## BP: TPM-normalized ranscripts -- CAZymes added together
vary_43 <- keggmodtable_gen_avg


newnames <- lapply(
  colnames(vary_43),
  function(x) bquote(.(x)))
newnames
head(vary_43)
nrow(vary_43)
#quantile(rowSums(vary_43))

data_matrix_44 <- vary_43
drop <- c("-", "S")
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

bp_keggmodtpm_heatmap<-pheatmap(data_matrix_44, col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=F, cluster_cols=F, 
                            fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_row = rownames(data_matrix_44), labels_col = as.character(newnames), border_color=NA,breaks=seq(0, 100000, length.out=300), angle_col = 45)


########################################################################################################################################################
# aldex on kegg groups
# do kegg groups on raw data

## Count matrix input
bpcoldata <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/bp_metadata_01132025.txt")
bpcoldata$name <- gsub("-", ".", bpcoldata$name)

bpcts <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/read_counts_gene_metassembly.tsv")

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

# Input gene annotation info
gene_anno <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/gene_annotation_metassembly.tsv")

bpcts$GeneID <- rownames(bpcts)
bpcts_geneanno <- merge(bpcts, gene_anno, by = "GeneID", all = TRUE)

library(dplyr)
library(stringr)
library(tidyr)

# let's now do the COG annotations
aa <- bpcts_geneanno %>%
  select(GeneID, everything())%>%
  mutate(KEGG_Module=str_split(KEGG_Module, ",")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:60,1)]
#rownames(aa) <- 1:nrow(aa)

unique(bpcts_geneanno$KEGG_Module)
ab<-as.vector(as.matrix(aa[,grepl( "KEGG_Module" , names(aa))]))
unique(ab)
length(unique(ab))

# make subset of aa data
aasub <- as.data.frame(aa[1:11,])
aasub$KEGG_Module1[is.na(aasub$KEGG_Module1)] <- "-"
aasubccs <- aasub[,c(33:50)]

uab <-unique(ab)
uab <- uab[!is.na(uab)]
rn <- uab
rn
keggmodtable_gen <- data.frame(matrix(0, nrow = length(unique(uab)), ncol = 20))
colnames(keggmodtable_gen) <- colnames(aa[,1:20])
rownames(keggmodtable_gen) <- uab

# Trial run on subset of the data
for (j in 1:ncol(keggmodtable_gen)) {
  for (k in 1:length(aasub[,grepl( "KEGG_Module" , names(aasub))])) { # Number of separate KEGG modules broken up
    for (i in 1:nrow(aasub)) {
      if(is.na(aasub[i,paste("KEGG_Module",k, sep = "")]) == FALSE) {
        keggmodtable_gen[rn[rownames(keggmodtable_gen) %in% aasub[i,paste("KEGG_Module",k, sep = "")]],j] = keggmodtable_gen[rn[rownames(keggmodtable_gen) %in% aasub[i,paste("KEGG_Module",k, sep = "")]],j] + aasub[i,j]
      }
      else {
        keggmodtable_gen = keggmodtable_gen
      }
    }
  }
}


## Now do it on the entire data
aa$KEGG_Module1[is.na(aa$KEGG_Module1)] <- "-"

uab <-unique(ab)
uab <- uab[!is.na(uab)]
rn <- uab
rn
keggmodtable_gen <- data.frame(matrix(0, nrow = length(unique(uab)), ncol = 20))
colnames(keggmodtable_gen) <- colnames(aa[,1:20])
rownames(keggmodtable_gen) <- uab

for (j in 1:ncol(keggmodtable_gen)) {
  for (k in 1:length(aa[,grepl( "KEGG_Module" , names(aa))])) { # Number of separate KEGG modules broken up
    for (i in 1:nrow(aa)) {
      if(is.na(aa[i,paste("KEGG_Module",k, sep = "")]) == FALSE) {
        keggmodtable_gen[rn[rownames(keggmodtable_gen) %in% aa[i,paste("KEGG_Module",k, sep = "")]],j] = keggmodtable_gen[rn[rownames(keggmodtable_gen) %in% aa[i,paste("KEGG_Module",k, sep = "")]],j] + aa[i,j]
      }
      else {
        keggmodtable_gen = keggmodtable_gen
      }
    }
  }
}

write_tsv(keggmodtable_gen, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggmodtable_gen_raw_03062025.tsv')
keggmodtable_gen <- read.csv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggmodtable_gen_raw_03062025.tsv', header = TRUE, sep = "")

# Take average of replicates
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
keggmodtable_gen_filt <- keggmodtable_gen[,!(names(keggmodtable_gen) %in% drops)]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(keggmodtable_gen_filt),]
# bpcoldata3 <- bpcoldata[bpcoldata$name %in% names(keggmodtable_gen_filt),]
# rownames(bpcoldata3) <- bpcoldata3$name
# bpcoldata3 <- bpcoldata3[,-1]

# keggmodtable_gen_avg <- data.frame(matrix(0, nrow = length(unique(uab)), ncol = 8))
# colnames(keggmodtable_gen_avg) <- (unique(bpcoldata3$condition_cba))
# rownames(keggmodtable_gen_avg) <- rownames(keggmodtable_gen_filt)
# 
# keggmodtable_gen_avg[,'CA + AS t1'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + AS t1',])]))
# keggmodtable_gen_avg[,'CA + S3 + AS t1'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + S3 + AS t1',])]))
# keggmodtable_gen_avg[,'CA + AS t2'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + AS t2',])]))
# keggmodtable_gen_avg[,'CA + S3 + AS t2'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + S3 + AS t2',])]))
# 
# keggmodtable_gen_avg[,'PHA + AS t1'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + AS t1',])]))
# keggmodtable_gen_avg[,'PHA + G1 + AS t1'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + G1 + AS t1',])]))
# keggmodtable_gen_avg[,'PHA + AS t2'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + AS t2',])]))
# keggmodtable_gen_avg[,'PHA + G1 + AS t2'] <- rowMeans(as.matrix(keggmodtable_gen_filt[, names(keggmodtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + G1 + AS t2',])]))

#### now do aldex# Now do aldex2 on the grouped raw counts 


library(ALDEx2)
bpselex <- keggmodtable_gen_filt
drop <- "-"
bpselex <- bpselex[!(row.names(bpselex) %in% drop), ]
names(keggmodtable_gen_filt) == rownames(bpcoldata3)

#bpx <- aldex(bpselex, bpconds, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, denom="iqlr", verbose=TRUE)

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "CA + AS" | bpcoldata3$group == "CA + S3 + AS",]

keggmodtable_gen_filt_sub <- keggmodtable_gen_filt[,colnames(keggmodtable_gen_filt) %in% rownames(bpcoldata3_sub)]
keggmodtable_gen_filt_sub_noNA <- keggmodtable_gen_filt_sub[!(row.names(keggmodtable_gen_filt_sub) %in% drop), ]

colnames(keggmodtable_gen_filt_sub) == rownames(bpcoldata3_sub)
colnames(keggmodtable_gen_filt_sub_noNA) == rownames(bpcoldata3_sub)

condsbp <- bpcoldata3_sub$group
condsbp

#bp <- aldex(dbcantable_gen_raw_filt_sub, condsbp, mc.samples=16, test="t", effect=TRUE, include.sample.summary=FALSE, denom="iqlr", verbose=FALSE)
bp <- aldex(keggmodtable_gen_filt_sub_noNA, condsbp, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, denom="iqlr", verbose=FALSE)

par(mfrow=c(1,3))
aldex.plot(bp, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference", main='Bland-Altman plot')
aldex.plot(bp, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference", main='Effect plot')
aldex.plot(bp, type="volcano", test="welch", xlab="Difference",
           ylab="-1(log10(q))", main='Volcano plot') 

# modular approach ??
x <- aldex.clr(keggmodtable_gen_filt_sub_noNA, condsbp, mc.samples=128, denom="all", verbose=F, gamma=1e-3)
x.tt <- aldex.ttest(x, hist.plot=F, paired.test=FALSE, verbose=FALSE)
x.effect <- aldex.effect(x, CI=T, verbose=F, include.sample.summary=F, 
                         paired.test=FALSE, glm.conds=NULL, useMC=F)
x.all <- data.frame(x.tt,x.effect)
par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="welch")
aldex.plot(x.all, type="MW", test="welch")

head(bp)

efbp1 <- bp[bp$effect >= 1,]
head(efbp1)
efbp1$keggmod <- row.names(efbp1)
efbp1_sig <- efbp1[efbp1$wi.eBH < 0.05,]

##############################################################################################################################
# Now do for PHA + G1 + AS vs. PHA + AS comparison
bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "PHA + G1 + AS" | bpcoldata3$group == "PHA + AS",]

keggmodtable_gen_filt_sub <- keggmodtable_gen_filt[,colnames(keggmodtable_gen_filt) %in% rownames(bpcoldata3_sub)]
keggmodtable_gen_filt_sub_noNA <- keggmodtable_gen_filt_sub[!(row.names(keggmodtable_gen_filt_sub) %in% drop), ]

colnames(keggmodtable_gen_filt_sub) == rownames(bpcoldata3_sub)
colnames(keggmodtable_gen_filt_sub_noNA) == rownames(bpcoldata3_sub)

condsbp <- bpcoldata3_sub$group
condsbp

#bp <- aldex(dbcantable_gen_raw_filt_sub, condsbp, mc.samples=16, test="t", effect=TRUE, include.sample.summary=FALSE, denom="iqlr", verbose=FALSE)
bp <- aldex(keggmodtable_gen_filt_sub_noNA, condsbp, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, denom="iqlr", verbose=FALSE)

par(mfrow=c(1,3))
aldex.plot(bp, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference", main='Bland-Altman plot')
aldex.plot(bp, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference", main='Effect plot')
aldex.plot(bp, type="volcano", test="welch", xlab="Difference",
           ylab="-1(log10(q))", main='Volcano plot') 

# modular approach ??
x <- aldex.clr(keggmodtable_gen_filt_sub_noNA, condsbp, mc.samples=128, denom="all", verbose=F, gamma=1e-3)
x.tt <- aldex.ttest(x, hist.plot=F, paired.test=FALSE, verbose=FALSE)
x.effect <- aldex.effect(x, CI=T, verbose=F, include.sample.summary=F, 
                         paired.test=FALSE, glm.conds=NULL, useMC=F)
x.all <- data.frame(x.tt,x.effect)
par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="welch")
aldex.plot(x.all, type="MW", test="welch")

head(bp)

efbp1 <- bp[bp$effect >= 1,]
head(efbp1)
efbp1$keggmod <- row.names(efbp1)
efbp1_sig <- efbp1[efbp1$wi.eBH < 0.05,]



