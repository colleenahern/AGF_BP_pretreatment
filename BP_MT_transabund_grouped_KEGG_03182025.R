# Colleen Ahern
# 03/18/2025

# Finalized script for normalizing my raw gene counts and add together KEGG modules
# Continue using TPM normalization
library(tidyverse)
library(DRnaSeq)
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

# split multiple KEGG annotations into individual columns for group counting
aa <- bptpm_geneanno %>%
  select(GeneID, everything())%>%
  mutate(KEGG_Module=str_split(KEGG_Module, ",")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:60,1)]

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


########################################################################################################################################################
# Aldex2 on kegg groups
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
  dplyr::select(GeneID, everything())%>%
  mutate(KEGG_Module=str_split(KEGG_Module, ",")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:60,1)]

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
keggmodtable_gen$modID <- rownames(keggmodtable_gen)
write_tsv(keggmodtable_gen, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggmodtable_gen_raw_03182025.tsv')
keggmodtable_gen <- read.csv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggmodtable_gen_raw_03182025.tsv', header = TRUE, sep = "")
rownames(keggmodtable_gen) <- keggmodtable_gen$modID


# drop any outlier samples (based on PCA analysis)
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
keggmodtable_gen_filt <- keggmodtable_gen[,!(names(keggmodtable_gen) %in% drops)]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(keggmodtable_gen_filt),]

#######################################################################################################################################
#### CA + S3 + AS vs. CA + AS comparison

library(ALDEx2)

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "CA + AS" | bpcoldata3$group == "CA + S3 + AS",]
drop <- "-"

keggmodtable_gen_filt_sub <- keggmodtable_gen_filt[,colnames(keggmodtable_gen_filt) %in% rownames(bpcoldata3_sub)]
keggmodtable_gen_filt_sub_noNA <- keggmodtable_gen_filt_sub[!(row.names(keggmodtable_gen_filt_sub) %in% drop), ]

colnames(keggmodtable_gen_filt_sub) == rownames(bpcoldata3_sub)
colnames(keggmodtable_gen_filt_sub_noNA) == rownames(bpcoldata3_sub)

## aldex2
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
drop <- "-"

keggmodtable_gen_filt_sub <- keggmodtable_gen_filt[,colnames(keggmodtable_gen_filt) %in% rownames(bpcoldata3_sub)]
keggmodtable_gen_filt_sub_noNA <- keggmodtable_gen_filt_sub[!(row.names(keggmodtable_gen_filt_sub) %in% drop), ]

colnames(keggmodtable_gen_filt_sub) == rownames(bpcoldata3_sub)
colnames(keggmodtable_gen_filt_sub_noNA) == rownames(bpcoldata3_sub)

# aldex2
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

#######################################################################################################################################
# 03/18/2025
# Also trying out ANCOM-BC because I dont know how to get a heatmap with aldex
# need to convert data into a phyloseq object
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

otumat = as.matrix(keggmodtable_gen_filt_sub_noNA[,!(names(keggmodtable_gen_filt_sub_noNA) %in% drop)])

taxmat = matrix(sample(letters, 86, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species")
taxmat[,"Genus"] <- rownames(taxmat)
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
col_name = c("Taxon", "Intercept", "CA film + S3 + AS")
col_name = c("Taxon", "Intercept", "PHA + G1 + AS")
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
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")

cadiff <- tab_diff[tab_diff$`CA film + S3 + AS` == TRUE,]
phadiff <- tab_diff[tab_diff$`PHA + G1 + AS` == TRUE,]
# ok results??

tab_lfc_filt <- tab_lfc[tab_lfc$Taxon %in% cadiff$Taxon,]
colnames(tab_lfc_filt) <- c("Taxon", "Intercept LFC", "CA film + S3 + AS LFC")
tab_lfc_filt <- tab_lfc[tab_lfc$Taxon %in% phadiff$Taxon,]
colnames(tab_lfc_filt) <- c("Taxon", "Intercept LFC", "PHA + G1 + AS LFC")

tab_se_filt <- tab_se[tab_se$Taxon %in% cadiff$Taxon,]
colnames(tab_se_filt) <- c("Taxon", "Intercept SE", "CA film + S3 + AS SE")
tab_se_filt <- tab_se[tab_se$Taxon %in% phadiff$Taxon,]
colnames(tab_se_filt) <- c("Taxon", "Intercept SE", "PHA + G1 + AS SE")

tab_w_filt <- tab_w[tab_w$Taxon %in% cadiff$Taxon,]
colnames(tab_w_filt) <- c("Taxon", "Intercept W", "CA film + S3 + AS W")
tab_w_filt <- tab_w[tab_w$Taxon %in% phadiff$Taxon,]
colnames(tab_w_filt) <- c("Taxon", "Intercept W", "PHA + G1 + AS W")

tab_p_filt <- tab_p[tab_p$Taxon %in% cadiff$Taxon,]
colnames(tab_p_filt) <- c("Taxon", "Intercept P", "CA film + S3 + AS P")
tab_p_filt <- tab_p[tab_p$Taxon %in% phadiff$Taxon,]
colnames(tab_p_filt) <- c("Taxon", "Intercept P", "PHA + G1 + AS P")

tab_q_filt <- tab_q[tab_q$Taxon %in% cadiff$Taxon,]
colnames(tab_q_filt) <- c("Taxon", "Intercept Q", "CA film + S3 + AS Q")
tab_q_filt <- tab_q[tab_q$Taxon %in% phadiff$Taxon,]
colnames(tab_q_filt) <- c("Taxon", "Intercept Q", "PHA + G1 + AS Q")

tab_diff_filt <- tab_diff[tab_diff$Taxon %in% cadiff$Taxon,]
colnames(tab_diff_filt) <- c("Taxon", "Intercept diff", "CA film + S3 + AS diff")
tab_diff_filt <- tab_diff[tab_diff$Taxon %in% phadiff$Taxon,]
colnames(tab_diff_filt) <- c("Taxon", "Intercept diff", "PHA + G1 + AS diff")

#put all data frames into list
df_list <- list(tab_lfc_filt, tab_se_filt, tab_w_filt, tab_p_filt, tab_q_filt, tab_diff_filt)

#merge all data frames in list
res_filt <- df_list %>% purrr::reduce(full_join, by='Taxon')

# save list
#writexl::write_xlsx(res_filt, "")

### Visualization of differentially abundant taxa heatmap

sig_taxa = tab_diff %>%
  dplyr::filter(`CA film + S3 + AS` == TRUE) %>%
  .$Taxon
sig_taxa = tab_diff %>%
  dplyr::filter(`PHA + G1 + AS` == TRUE) %>%
  .$Taxon

tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "CA film + S3 + AS vs. CA film + AS")
col_name = c("Taxon", "Intercept", "PHA + G1 + AS vs. PHA + AS")
colnames(tab_lfc) = col_name
tab_lfc <- tab_lfc[,-2]

df_CAS3AS = tab_lfc %>%
  filter(Taxon %in% sig_taxa) %>%
  dplyr::rename(Module = Taxon)
df_PHAG1AS = tab_lfc %>%
  filter(Taxon %in% sig_taxa) %>%
  dplyr::rename(Module = Taxon)

modtab <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/KEGG_mod_funcs_05052025.txt", header = TRUE)

df_CAS3ASmod <- merge(df_CAS3AS, modtab,  by = "Module", all.x = TRUE)
df_PHAG1ASmod <- merge(df_PHAG1AS, modtab,  by = "Module", all.x = TRUE)

df_heat = df_CAS3AS %>%
  pivot_longer(cols = -one_of("Module"),
               names_to = "region", values_to = "value") %>%
  mutate(value = round(value, 2))
df_heat$Module = factor(df_heat$Module, levels = sort(sig_taxa))
df_heatmod <- merge(df_heat, modtab,  by = "Module", all.x = TRUE)

df_heat = df_PHAG1AS %>%
  pivot_longer(cols = -one_of("Module"),
               names_to = "region", values_to = "value") %>%
  mutate(value = round(value, 2))
df_heat$Module = factor(df_heat$Module, levels = sort(sig_taxa))
df_heatmod <- merge(df_heat, modtab,  by = "Module", all.x = TRUE)

lo = floor(min(df_heat$value))
up = ceiling(max(df_heat$value))
mid = (lo + up)/2
df_heat_filt <- df_heat[abs(df_heat$value) > 1,]
df_heatmod_filt <- df_heatmod[abs(df_heatmod$value) > 1,]
df_heatmod_filt$Label <- paste(df_heatmod_filt$Module, ": ", df_heatmod_filt$Function, sep = "")

# Manually add annotations for unannotated modules
# CA S3 AS
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00772"] <- gsub("NA", "HupT-HupR (hydrogenase synthesis regulation) two-component regulatory system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00772"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00686"] <- gsub("NA", "Toll-like receptor signaling", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00686"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00685"] <- gsub("NA", "Apoptotic machinery", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00685"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00683"] <- gsub("NA", "Hippo signaling", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00683"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00678"] <- gsub("NA", "Hedgehog signaling", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00678"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00670"] <- gsub("NA", "Mce transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00670"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00669"] <- gsub("NA", "gamma-Hexachlorocyclohexane transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00669"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00662"] <- gsub("NA", "Hk1-Rrp1 (glycerol uptake and utilization) two-component regulatory system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00662"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00604"] <- gsub("NA", "Trehalose transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00604"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00599"] <- gsub("NA", "Inositol-phosphate transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00599"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00587"] <- gsub("NA", "Arginine/lysine/histidine/glutamine transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00587"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00584"] <- gsub("NA", "Acetoin utilization transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00584"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00518"] <- gsub("NA", "GlnK-GlnL (glutamine utilization) two-component regulatory system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00518"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00508"] <- gsub("NA", "PixL-PixGH (positive phototaxis) two-component regulatory system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00508"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00502"] <- gsub("NA", "GlrK-GlrR (amino sugar metabolism) two-component regulatory system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00502"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00490"] <- gsub("NA", "MalK-MalR (malate transport) two-component regulatory system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00490"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00475"] <- gsub("NA", "BarA-UvrY (central carbon metabolism) two-component regulatory system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00475"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00461"] <- gsub("NA", "MtrB-MtrA (osmotic stress response) two-component regulatory system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00461"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00430"] <- gsub("NA", "Exon junction complex (EJC)", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00430"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00429"] <- gsub("NA", "Competence-related DNA transformation transporter", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00429"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00412"] <- gsub("NA", "ESCRT-III complex", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00412"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00403"] <- gsub("NA", "HRD1/SEL1 ERAD complex", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00403"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00400"] <- gsub("NA", "p97-Ufd1-Npl4 complex", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00400"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00395"] <- gsub("NA", "Decapping complex", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00395"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00322"] <- gsub("NA", "Neutral amino acid transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00322"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00298"] <- gsub("NA", "Multidrug/hemolysin transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00298"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00297"] <- gsub("NA", "DNA-PK complex", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00297"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00296"] <- gsub("NA", "BER complex", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00296"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00280"] <- gsub("NA", "Phosphotransferase system, glucitol/sorbitol-specific II component", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00280"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00274"] <- gsub("NA", "Phosphotransferase system, mannitol-specific II component", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00274"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00259"] <- gsub("NA", "Heme transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00259"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00251"] <- gsub("NA", "Teichoic acid transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00251"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00232"] <- gsub("NA", "General L-amino acid transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00232"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00205"] <- gsub("NA", "N-Acetylglucosamine transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00205"]) # from https://www.nature.com/articles/s41522-025-00679-w

df_heatmod_filt2 <- df_heatmod_filt[,c(1:3,6)]
write_xlsx(df_heatmod_filt2, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAS3AS_KEGGancombc_05052025.xlsx")


# PHA G1 AS
# natkeggmod <- read.csv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/KeggModuleStepHitInfo.csv")
# natkeggmod$Module.step <- gsub("\\+.*","",natkeggmod$Module.step)
# 
# for (i in nrow(df_heatmod_filt)) {
#   if(is.na(df_heatmod_filt$Function[i] == TRUE)) {
#     df_heatmod_filt$Label[df_heatmod_filt$Module == df_heatmod_filt$Module[i]] <- gsub("NA",natkeggmod$Module , df_heatmod_filt$Label[df_heatmod_filt$Module =="M00205"]) # from https://www.nature.com/articles/s41522-025-00679-w
#     
#   }
# }
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00809"] <- gsub("NA", "Phosphotransferase system, glucose-specific II component", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00809"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00765"] <- gsub("NA", "Multidrug resistance, efflux pump Bmr", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00765"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00678"] <- gsub("NA", "Hedgehog signaling", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00678"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00610"] <- gsub("NA", "Phosphotransferase system, D-glucosaminate-specific II component", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00610"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00566"] <- gsub("NA", "Dipeptide transport system, Firmicutes", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00566"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00502"] <- gsub("NA", "GlrK-GlrR (amino sugar metabolism) two-component regulatory system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00502"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00495"] <- gsub("NA", "AgrC-AgrA (exoprotein synthesis) two-component regulatory system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00495"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00492"] <- gsub("NA", "LytS-LytR two-component regulatory system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00492"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00403"] <- gsub("NA", "HRD1/SEL1 ERAD complex", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00403"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00400"] <- gsub("NA", "p97-Ufd1-Npl4 complex", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00400"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00391"] <- gsub("NA", "Exosome, eukaryotes", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00391"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00390"] <- gsub("NA", "Exosome, archaea", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00390"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00343"] <- gsub("NA", "Archaeal proteasome", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00343"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00342"] <- gsub("NA", "Bacterial proteasome", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00342"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00322"] <- gsub("NA", "Neutral amino acid transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00322"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00295"] <- gsub("NA", "BRCA1-associated genome surveillance complex (BASC)", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00295"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00288"] <- gsub("NA", "RPA complex", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00288"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00280"] <- gsub("NA", "Phosphotransferase system, glucitol/sorbitol-specific II component", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00280"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00246"] <- gsub("NA", "Nickel transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00246"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00245"] <- gsub("NA", "Cobalt/nickel transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00245"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00244"] <- gsub("NA", "Putative zinc/manganese transport system", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00244"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00184"] <- gsub("NA", "RNA polymerase, archaea", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00184"]) # from https://www.nature.com/articles/s41522-025-00679-w
df_heatmod_filt$Label[df_heatmod_filt$Module =="M00177"] <- gsub("NA", "Ribosome, eukaryotes", df_heatmod_filt$Label[df_heatmod_filt$Module =="M00177"]) # from https://www.nature.com/articles/s41522-025-00679-w

df_heatmod_filt2 <- df_heatmod_filt[,c(1:3,6)]
write_xlsx(df_heatmod_filt2, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/PHAG1AS_KEGGancombc_05052025.xlsx")

p_heat = df_heat_filt %>%
  ggplot(aes(x = region, y = Module, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(region, Module, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, size = 5,
       title = "Log fold changes for globally significant modules") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 14))
p_heat

p_heatmod = df_heatmod_filt %>%
  ggplot(aes(x = region, y = Label, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(region, Label, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, size = 5,
       title = "Log fold changes for globally significant modules") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 14))
p_heatmod
