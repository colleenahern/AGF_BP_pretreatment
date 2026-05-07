# Colleen Ahern
# 10/30/2025

# Magda asked that I provide her log fold change values for KEGG ko's instead of modules
# and to include taxonomic information
# Fixing my taxa associated with each group of interest (in this case KEGG kos) because Malte said I was using the wrong file in our meeting
# Need to use tax_gtdb_metassembly.tsv that links contigs to taxa
# Need to use bowtie2_metaassembly.stranded2.counts.txt to link genes to contigs
# The same way I fixed this for the CAZyme subgroups
# Also make volcano plots instead of heatmaps to visualize the data

# We will use raw counts, not TPM normalized counts since I will be using ANCOM-BC
library(tidyverse)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

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
bpcts_geneanno2$KEGG_ko <- gsub("ko:","",bpcts_geneanno2$KEGG_ko)

library(dplyr)
library(stringr)
library(tidyr)

# let's now do the KEGG ko annotations
aa <- bpcts_geneanno2 %>%
  dplyr::select(GeneID, everything())%>%
  mutate(KEGG_ko=str_split(KEGG_ko, ",")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:63,1)]

unique(bpcts_geneanno2$KEGG_ko)
length(unique(bpcts_geneanno2$KEGG_ko))
ab<-as.vector(as.matrix(aa[,grepl("KEGG_ko",names(aa))]))
unique(ab)
length(unique(ab))

# make subset of aa data
bb <- aa[is.na(aa$KEGG_ko3) == FALSE,]
aasub <- as.data.frame(bb[1:5,])
aasub$KEGG_ko1[is.na(aasub$KEGG_ko1)] <- "-"
aasubccs <- aasub[,c(31:43)]
# this above subset didn't have overlapping KO's to make sure that when there is one KO represented by multile genes, those gene counts are added together
aasub <- as.data.frame(aa[c(162484,179646),])
aasub$KEGG_ko1[is.na(aasub$KEGG_ko1)] <- "-"
aasubccs <- aasub[,c(31:43)]

uab <-unique(ab)
uab
rn <- uab[!is.na(uab)]
rn

keggkotable_raw <- data.frame(matrix(0, nrow = length(rn), ncol = 20))
colnames(keggkotable_raw) <- colnames(aasub[,1:20])
rownames(keggkotable_raw) <- rn
keggkotable_raw$taxa <- NA

# Trial run on subset of the data
for (j in 1:(ncol(keggkotable_raw)-1)) {
  for (k in 1:length(names(aasub)[grepl("KEGG_ko" , names(aasub))])) { # Number of separate KEGG kos broken up
    for (i in 1:nrow(aasub)) {
      if(is.na(aasub[i,paste("KEGG_ko",k, sep = "")]) == FALSE) {
        keggkotable_raw[rn[rn %in% aasub[i,paste("KEGG_ko",k, sep = "")]],j] = keggkotable_raw[rn[rn %in% aasub[i,paste("KEGG_ko",k, sep = "")]],j] + as.numeric(aasub[i,j])
      }
      else {
        keggkotable_raw = keggkotable_raw
      }
    }
  }
}

for (k in 1:length(aasub[,grepl("KEGG_ko" , names(aasub))])) { # Number of separate KEGG kos broken up
  for (i in 1:nrow(aasub)) {
    if(aasub[i,paste("KEGG_ko",k, sep = "")] %in% rownames(keggkotable_raw) == TRUE ) {
      keggkotable_raw[rn[rn %in% aasub[i,paste("KEGG_ko",k, sep = "")]],"taxa"] <- paste(keggkotable_raw[rn[rn %in% aasub[i,paste("KEGG_ko",k, sep = "")]],"taxa"],"; ", aasub$taxa[i], sep = "")
    }
    else {
      keggkotable_raw = keggkotable_raw
    }
  }
}


## Now do it on the entire data
aa$KEGG_ko1[is.na(aa$KEGG_ko1)] <- "-"

uab <-unique(ab)
uab
rn <- uab[!is.na(uab)]
rn
length(rn)

keggkotable_raw <- data.frame(matrix(0, nrow = length(unique(rn)), ncol = 20))
colnames(keggkotable_raw) <- colnames(aa[,1:20])
rownames(keggkotable_raw) <- rn
keggkotable_raw$taxa <- NA

for (j in 1:(ncol(keggkotable_raw)-1)) {
  for (k in 1:length(names(aa)[grepl("KEGG_ko" , names(aa))])) { # Number of separate KEGG kos broken up
    for (i in 1:nrow(aa)) {
      if(is.na(aa[i,paste("KEGG_ko",k, sep = "")]) == FALSE) {
        keggkotable_raw[rn[rn %in% aa[i,paste("KEGG_ko",k, sep = "")]],j] = keggkotable_raw[rn[rn %in% aa[i,paste("KEGG_ko",k, sep = "")]],j] + as.numeric(aa[i,j])
      }
      else {
        keggkotable_raw = keggkotable_raw
      }
    }
  }
}

for (k in 1:length(aa[,grepl("KEGG_ko" , names(aa))])) { # Number of separate KEGG kos broken up
  for (i in 1:nrow(aa)) {
    if(aa[i,paste("KEGG_ko",k, sep = "")] %in% rownames(keggkotable_raw) == TRUE ) {
      keggkotable_raw[rn[rn %in% aa[i,paste("KEGG_ko",k, sep = "")]],"taxa"] <- paste(keggkotable_raw[rn[rn %in% aa[i,paste("KEGG_ko",k, sep = "")]],"taxa"],"; ", aa$taxa[i], sep = "")
    }
    else {
      keggkotable_raw = keggkotable_raw
    }
  }
}

keggkotable_raw$koID <- rownames(keggkotable_raw) 
write_tsv(keggkotable_raw, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggkotable_raw_CORRECT_taxa_11042025.tsv')
keggkotable_raw_r <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggkotable_raw_CORRECT_taxa_11042025.tsv')
keggkotable_raw_r <- as.data.frame(keggkotable_raw_r)
rownames(keggkotable_raw_r) <- keggkotable_raw_r$koID

# Remove the results for "-" since I think this is crashing my computer, plus I don't care about the taxa associated with non-CAZymes
drop <- c("-")
keggkotable_rawrrr <- keggkotable_raw_r[!(keggkotable_raw_r$koID %in% drop),]
write_tsv(keggkotable_rawrrr, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggkotable_raw_CORRECT_taxa_NOunknown_11042025.tsv')

keggkotable_rawrrrr <- keggkotable_rawrrr
keggkotable_rawrrrr$taxa <- substr(keggkotable_rawrrrr$taxa, 4, nchar(keggkotable_rawrrrr$taxa))
write_tsv(keggkotable_rawrrrr, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggkotable_raw_CORRECT_taxa_NOunknown_noNAstart_11042025.tsv')
keggkotable_rawrrrr <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggkotable_raw_CORRECT_taxa_NOunknown_noNAstart_11042025.tsv')
keggkotable_rawrrrr <- as.data.frame(keggkotable_rawrrrr)
row.names(keggkotable_rawrrrr) <- NULL

taxaf <- keggkotable_rawrrrr[,c("koID","taxa")]
taxaf <- taxaf %>%
  column_to_rownames("koID")

keggkotable_rawrt <- keggkotable_rawrrrr %>%
  column_to_rownames("koID")

# drop any outlier samples (based on PCA analysis)
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC","taxa")
keggkotable_rawrt_filt <- keggkotable_rawrt[,!(names(keggkotable_rawrt) %in% drops)]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(keggkotable_rawrt_filt),]

#######################################################################################################################################
#### CA + S3 + AS vs. CA + AS comparison

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "CA + AS" | bpcoldata3$group == "CA + S3 + AS",]

keggkotable_raw_filt_sub <- keggkotable_rawrt_filt[,colnames(keggkotable_rawrt_filt) %in% rownames(bpcoldata3_sub)]

all(colnames(keggkotable_raw_filt_sub) == rownames(bpcoldata3_sub))

##############################################################################################################################
# Now do for PHA + G1 + AS vs. PHA + AS comparison

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "PHA + G1 + AS" | bpcoldata3$group == "PHA + AS",]

keggkotable_raw_filt_sub <- keggkotable_rawrt_filt[,colnames(keggkotable_rawrt_filt) %in% rownames(bpcoldata3_sub)]

all(colnames(keggkotable_raw_filt_sub) == rownames(bpcoldata3_sub))


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

otumat = as.matrix(keggkotable_raw_filt_sub)

taxmat = matrix(sample(letters, 62, replace = TRUE), nrow = nrow(otumat), ncol = 7)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species") # add on taxa later after analysis so it doesn't crash
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
col_name = c("KO", "Intercept", "CA + S3 + AS")
col_name = c("KO", "Intercept", "PHA + G1 + AS")
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
  datatable(caption = "Differentially Abundant KOs from the Primary Result")

cadiff <- tab_diff[tab_diff$`CA + S3 + AS` == TRUE,]
phadiff <- tab_diff[tab_diff$`PHA + G1 + AS` == TRUE,]

tab_lfc_filt <- tab_lfc[tab_lfc$KO %in% cadiff$KO,]
colnames(tab_lfc_filt) <- c("KO", "Intercept LFC", "CA + S3 + AS LFC")
tab_lfc_filt <- tab_lfc[tab_lfc$KO %in% phadiff$KO,]
colnames(tab_lfc_filt) <- c("KO", "Intercept LFC", "PHA + G1 + AS LFC")

tab_se_filt <- tab_se[tab_se$KO %in% cadiff$KO,]
colnames(tab_se_filt) <- c("KO", "Intercept SE", "CA + S3 + AS SE")
tab_se_filt <- tab_se[tab_se$KO %in% phadiff$KO,]
colnames(tab_se_filt) <- c("KO", "Intercept SE", "PHA + G1 + AS SE")

tab_w_filt <- tab_w[tab_w$KO %in% cadiff$KO,]
colnames(tab_w_filt) <- c("KO", "Intercept W", "CA + S3 + AS W")
tab_w_filt <- tab_w[tab_w$KO %in% phadiff$KO,]
colnames(tab_w_filt) <- c("KO", "Intercept W", "PHA + G1 + AS W")

tab_p_filt <- tab_p[tab_p$KO %in% cadiff$KO,]
colnames(tab_p_filt) <- c("KO", "Intercept P", "CA + S3 + AS P")
tab_p_filt <- tab_p[tab_p$KO %in% phadiff$KO,]
colnames(tab_p_filt) <- c("KO", "Intercept P", "PHA + G1 + AS P")

tab_q_filt <- tab_q[tab_q$KO %in% cadiff$KO,]
colnames(tab_q_filt) <- c("KO", "Intercept Q", "CA + S3 + AS Q")
tab_q_filt <- tab_q[tab_q$KO %in% phadiff$KO,]
colnames(tab_q_filt) <- c("KO", "Intercept Q", "PHA + G1 + AS Q")

tab_diff_filt <- tab_diff[tab_diff$KO %in% cadiff$KO,]
colnames(tab_diff_filt) <- c("KO", "Intercept diff", "CA + S3 + AS diff")
tab_diff_filt <- tab_diff[tab_diff$KO %in% phadiff$KO,]
colnames(tab_diff_filt) <- c("KO", "Intercept diff", "PHA + G1 + AS diff")

#put all data frames into list
# df_list <- list(tab_lfc_filt, tab_se_filt, tab_w_filt, tab_p_filt, tab_q_filt, tab_diff_filt)

#merge all data frames in list
# res_filt <- df_list %>% purrr::reduce(full_join, by='KO')

# save list
#writexl::write_xlsx(res_filt, "")

### Visualization of differentially abundant taxa heatmap

# sig_taxa = tab_diff %>%
#   dplyr::filter(`CA + S3 + AS` == TRUE) %>%
#   .$KO
# sig_taxa = tab_diff %>%
#   dplyr::filter(`PHA + G1 + AS` == TRUE) %>%
#   .$KO

tab_lfc = res$lfc
col_name = c("KO", "Intercept", "CA + S3 + AS vs. CA + AS")
col_name = c("KO", "Intercept", "PHA + G1 + AS vs. PHA + AS")
colnames(tab_lfc) = col_name
tab_lfc <- tab_lfc[,-2]

taxaf$KO <- rownames(taxaf)
# df_CAS3AS = tab_lfc %>%
#   filter(Taxon %in% sig_taxa) %>%
#   dplyr::rename(koID = Taxon)

tab_lfc_CAS3AS <- merge(tab_lfc, taxaf, by = "KO")
write_tsv(tab_lfc_CAS3AS, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/tab_lfc_CAS3AS_KEGGko_ancombc_11052025.tsv")
tab_lfc_CAS3AS <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/tab_lfc_CAS3AS_KEGGko_ancombc_11052025.tsv")

tab_lfc_CAS3AS2 <- tab_lfc_CAS3AS %>%
  rename("LFC" = "CA + S3 + AS vs. CA + AS")

tab_q2 <- tab_q %>%
  rename( "padj" = "CA + S3 + AS")
tab_q2 <- tab_q2[,-2]

df_volc_CAS3AS2 <- merge(tab_lfc_CAS3AS2, tab_q2, by = "KO")
write_tsv(df_volc_CAS3AS2, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/df_volc_CAS3AS2_KEGGko_ancombc_11052025.tsv") # use this for volcano plots
df_volc_CAS3AS2 <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/df_volc_CAS3AS2_KEGGko_ancombc_11052025.tsv")
kegganno <- read_csv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/KEGG_ko_annolist_12152025.csv")
df_volc_CAS3AS2_ko <- merge(df_volc_CAS3AS2, kegganno, by = "KO", all.x = TRUE)
df_volc_CAS3AS2_ko <- df_volc_CAS3AS2_ko[order(df_volc_CAS3AS2_ko$padj),]
write_tsv(df_volc_CAS3AS2_ko, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/df_volc_CAS3AS2_koanno_12162025.tsv") # but update this to do analysis not excluding those two samples !! and then redo this
df_volc_CAS3AS2_ko <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/df_volc_CAS3AS2_koanno_12162025.tsv")
df_volc_CAS3AS2_ko_up <- df_volc_CAS3AS2_ko[df_volc_CAS3AS2_ko$LFC>0,]
df_volc_CAS3AS2_ko_down <- df_volc_CAS3AS2_ko[df_volc_CAS3AS2_ko$LFC<0,]
write_tsv(df_volc_CAS3AS2_ko_up, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/df_volc_CAS3AS2_koanno_UP_12162025.tsv")

# PHA
tab_lfc_PHAG1AS <- merge(tab_lfc, taxaf, by = "KO")
write_tsv(tab_lfc_PHAG1AS, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/tab_lfc_PHAG1AS_KEGGko_ancombc_11052025.tsv")
tab_lfc_PHAG1AS <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/tab_lfc_PHAG1AS_KEGGko_ancombc_11052025.tsv")

tab_lfc_PHAG1AS2 <- tab_lfc_PHAG1AS %>%
  rename("LFC" = "PHA + G1 + AS vs. PHA + AS")

tab_q2 <- tab_q %>%
  rename("padj" = "PHA + G1 + AS")
tab_q2 <- tab_q2[,-2]

df_volc_PHAG1AS2 <- merge(tab_lfc_PHAG1AS2, tab_q2, by = "KO")
write_tsv(df_volc_PHAG1AS2, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/df_volc_PHAG1AS2_KEGGko_ancombc_11052025.tsv") # use this for volcano plots
df_volc_PHAG1AS2 <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/df_volc_PHAG1AS2_KEGGko_ancombc_11052025.tsv")
kegganno <- read_csv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/KEGG_ko_annolist_12152025.csv")
df_volc_PHAG1AS2_ko <- merge(df_volc_PHAG1AS2, kegganno, by = "KO", all.x = TRUE)
df_volc_PHAG1AS2_ko <- df_volc_PHAG1AS2_ko[order(df_volc_PHAG1AS2_ko$padj),]
write_tsv(df_volc_PHAG1AS2_ko, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/df_volc_PHAG1AS2_koanno_12162025.tsv") # but update this to do analysis not excluding those two samples !! and then redo this
df_volc_PHAG1AS2_ko <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/df_volc_PHAG1AS2_koanno_12162025.tsv")
df_volc_PHAG1AS2_ko_up <- df_volc_PHAG1AS2_ko[df_volc_PHAG1AS2_ko$LFC > 0,]
df_volc_PHAG1AS2_ko_down <- df_volc_PHAG1AS2_ko[df_volc_PHAG1AS2_ko$LFC < 0,]
write_tsv(df_volc_PHAG1AS2_ko_up,  "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/df_volc_PHAG1AS2_koanno_UP_12162025.tsv")

########################################################################################################################################
# making volcano plots instead of heatmmaps to show differentially expressed CAZymes
# taken from https://biostatsquid.com/volcano-plots-r-tutorial/

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

# Import Ancombc results

# for CAS3AS vs CAAS comparison

dfkocas3 <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/df_volc_CAS3AS2_KEGGko_ancombc_11052025.tsv")

ggplot(data = dfkocas3, aes(x = LFC, y = -log10(padj))) +
  geom_point()

# Add threshold lines
ggplot(data = dfkocas3, aes(x = LFC, y = -log10(padj))) +
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

ggplot(data = dfkocas3, aes(x = LFC, y = -log10(padj))) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point() 

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
dfkocas3$diffexpressed <- "NO"

# if log2Foldchange > 1 and pvalue < 0.05, set as "UP"
dfkocas3$diffexpressed[dfkocas3$LFC > 1 & dfkocas3$padj < 0.05] <- "UP"

# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
dfkocas3$diffexpressed[dfkocas3$LFC < -1 & dfkocas3$padj < 0.05] <- "DOWN"

# Explore a bit
head(dfkocas3[order(dfkocas3$padj) & dfkocas3$diffexpressed == 'DOWN', ])

ggplot(data = dfkocas3, aes(x = LFC, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  labs(color='Differential Expression') 


# Edit axis labels and limits
ggplot(data = dfkocas3, aes(x = LFC, y = -log10(padj), col = diffexpressed)) +
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
dlup <- dfkocas3[dfkocas3$diffexpressed == "UP",]
dfkocas3$delabel <- ifelse(dfkocas3$KO %in% as.matrix(head(dlup[order(dlup$padj), "KO"], 10)), dfkocas3$KO, NA)

dldown <- dfkocas3[dfkocas3$diffexpressed == "DOWN",]
dfkocas3$delabel <- ifelse(dfkocas3$KO %in% as.matrix(head(dldown[order(dldown$padj), "KO"], 10)), dfkocas3$KO, dfkocas3$delabel)

ggplot(data = dfkocas3, aes(x = LFC, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 300), xlim = c(-6, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-6, 10, 2)) + # to customise the breaks in the x axis
  #geom_point() + geom_text(show.legend = FALSE, nudge_x = .2, nudge_y = 1)
  geom_label_repel(size = 4, show.legend = FALSE, min.segment.length = unit(0, 'lines')) +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), legend.title = element_blank(), legend.position = "bottom")


## 10/05/2025 ^ ok so for the graph above I am labeling the 10 most significantly expressed (lowest padj values) upregulated CAZyme subgroup and 
# the 10 most significantly expressed (lowest padj values) downregulated CAZyme subgroup and 
# This also includes an abs(1) LFC cutoff, which excludes some of the lowest p value CAZ subgroups, so maybe change that? Talk to Magda
# see if we want to keep the LFC cutoff or get rid of it

# for G1PHAAS vs PHAAS comparison

dfkophag1 <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/df_volc_PHAG1AS2_KEGGko_ancombc_11052025.tsv")

ggplot(data = dfkophag1, aes(x = LFC, y = -log10(padj))) +
  geom_point()

# Add threshold lines
ggplot(data = dfkophag1, aes(x = LFC, y = -log10(padj))) +
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

ggplot(data = dfkophag1, aes(x = LFC, y = -log10(padj))) +
  geom_vline(xintercept = c(-0.6, 0.6), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') + 
  geom_point() 

# Add a column to the data frame to specify if they are UP- or DOWN- regulated (log2fc respectively positive or negative)<br /><br /><br />
dfkophag1$diffexpressed <- "NO"

# if log2Foldchange > 0.6 and pvalue < 0.05, set as "UP"
dfkophag1$diffexpressed[dfkophag1$LFC > 1 & dfkophag1$padj < 0.05] <- "UP"

# if log2Foldchange < -0.6 and pvalue < 0.05, set as "DOWN"
dfkophag1$diffexpressed[dfkophag1$LFC < -1 & dfkophag1$padj < 0.05] <- "DOWN"

# Explore a bit
head(dfkophag1[order(dfkophag1$padj) & dfkophag1$diffexpressed == 'DOWN', ])

ggplot(data = dfkophag1, aes(x = LFC, y = -log10(padj), col = diffexpressed)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  labs(color='Differential Expression') 


# Edit axis labels and limits
ggplot(data = dfkophag1, aes(x = LFC, y = -log10(padj), col = diffexpressed)) +
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
dlup <- dfkophag1[dfkophag1$diffexpressed == "UP",]
dfkophag1$delabel <- ifelse(dfkophag1$KO %in% as.matrix(head(dlup[order(dlup$padj), "KO"], 10)), dfkophag1$KO, NA)

dldown <- dfkophag1[dfkophag1$diffexpressed == "DOWN",]
dfkophag1$delabel <- ifelse(dfkophag1$KO %in% as.matrix(head(dldown[order(dldown$padj), "KO"], 10)), dfkophag1$KO, dfkophag1$delabel)

ggplot(data = dfkophag1, aes(x = LFC, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 100), xlim = c(-5, 5)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-5, 5, 2)) + # to customise the breaks in the x axis
  #geom_point() + geom_text(show.legend = FALSE, nudge_x = .2, nudge_y = 1)
  geom_label_repel(size = 4, show.legend = FALSE, min.segment.length = unit(0, 'lines')) +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), legend.title = element_blank(), legend.position = "bottom")


# 02/05/2026 make final figure of these plots then combine with the CAZsub volcano plots
# for S3 CA AS vs CA AS
kocas3kop <- ggplot(data = dfkocas3, aes(x = LFC, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 300), xlim = c(-6, 10)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(title = "CA + S3 + AS vs. CA + AS", color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-6, 10, 2)) + # to customise the breaks in the x axis
  geom_label_repel(size = 3, show.legend = FALSE, min.segment.length = unit(0, 'lines'), max.overlaps = Inf) +
  theme_classic() +
  theme(plot.title = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 12), legend.title = element_blank(), 
        legend.text = element_text(size = 12), legend.position = "bottom", plot.margin = margin(23, 5, 5, 5))
kocas3kop

# for G1 PHA AS vs PHA AS
kophag1kop <- ggplot(data = dfkophag1, aes(x = LFC, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 1) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 100), xlim = c(-5, 5)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(title = "PHA + G1 + AS vs. PHA + AS", color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-5, 5, 2)) + # to customise the breaks in the x axis
  geom_label_repel(size = 3, show.legend = FALSE, min.segment.length = unit(0, 'lines'), max.overlaps = Inf) +
  theme_classic() +
  theme(plot.title = element_blank(), axis.text = element_text(size = 12), axis.title = element_text(size = 12), legend.title = element_blank(), 
        legend.text = element_text(size = 12), legend.position = "bottom", plot.margin = margin(23, 5, 5, 5))
kophag1kop

library(ggpubr)
kop <- ggarrange(kocas3kop, kophag1kop, ncol=2, nrow=1, common.legend = TRUE, legend="bottom")
kop 

fin_kop <- ggarrange(np1, np2, kocas3kop, kophag1kop, ncol=2, nrow=2, common.legend = TRUE, legend="bottom", labels = c("A)", "B)", "C)", "D)"))
fin_kop

ggsave(plot = fin_kop, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/CAZsub_KO_volcplotbothcond_02062026.pdf", width = 10, height = 10, units = "in", dpi = 300)
ggsave(plot = fin_kop, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/CAZsub_KO_volcplotbothcond_02062026.tiff", width = 10, height = 10, units = "in", dpi = 200, bg="white")
















##
abc <- rownames(keggkotable_rawrt)
abc <- as.data.frame(abc)
writexl::write_xlsx(abc, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/KEGG_ko_list_11062025.xlsx")
abc <- abc %>%
  rename("KOid" = "abc")

def <- readxl::read_xlsx("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/KEGG_kolist_ALL_11062025.xlsx")
gh <- merge(abc, def, by = "KOid", all.x = TRUE)

lala <- df[df$diffexpressed == "UP" | df$diffexpressed == "DOWN",]
lala <- lala %>%
  rename("KOid" = "KO")
lala2 <- merge(lala, gh, by = "KOid")

lalatop <- df[is.na(df$delabel) == "FALSE",]
lalatop <- lalatop %>%
  rename("KOid" = "KO")
lalatop2 <- merge(lalatop, gh, by = "KOid")

write_tsv(lala2, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/PHAG1AS2_KEGGko_sig_labeled_11102025.tsv")
write_tsv(lala2, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAS3AS2_KEGGko_sig_labeled_11102025.tsv")

dfkophag1_anno <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/PHAG1AS2_KEGGko_sig_labeled_11102025.tsv")
dfkocas3_anno <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAS3AS2_KEGGko_sig_labeled_11102025.tsv")

dfkocas3_filt <- dfkocas3[dfkocas3$KO %in% dfkocas3_anno$KOid,]
all(dfkocas3_anno$KOid %in% dfkocas3$KO)

dfkocas3_anno_upsig <- dfkocas3_anno[dfkocas3_anno$LFC > 0 & dfkocas3_anno$padj < 0.05,]
dfkocas3_anno_upsig <- dfkocas3_anno_upsig %>%
                        arrange(padj)
dfkophag1_anno_upsig <- dfkophag1_anno[dfkophag1_anno$LFC > 0 & dfkophag1_anno$padj < 0.05,]
dfkophag1_anno_upsig <- dfkophag1_anno_upsig %>%
                        arrange(padj)









## 06/02/2025 merging my CAS3AS_KEGGko and PHAG1AS_KEGGko dataframes with the KEGG pathway-ko
# to link my ko's to their pathways
kopath <- read.delim('/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/keggmods_paths_06022025.txt')
kopath_filt <- kopath[!grepl("path:ko",kopath$Pathway),]
kopath_filt$KO <- gsub("ko:","",kopath_filt$KO)
kopath_filt$Pathway <- gsub("path:","",kopath_filt$Pathway)
colnames(kopath_filt) <- c("koID", "Pathway")
CAS3AS_KEGGkopath <- merge(kopath_filt,CAS3AS_KEGGko, by = "koID", all.y = TRUE)

library(dplyr)
CAS3AS_KEGGkopath2 <- CAS3AS_KEGGkopath %>% 
  group_by(koID) %>% 
  mutate(Pathway = paste(Pathway, collapse=",")) %>%
  distinct(koID,`CA + S3 + AS vs. CA + AS`,taxa,Pathway)

df_PHAG1AS = tab_lfc %>%
  filter(Taxon %in% sig_taxa) %>%
  dplyr::rename(koID = Taxon)

keggkotable_raw_taxa$taxa <- gsub("NA; ","",keggkotable_raw_taxa$taxa)
df_PHAG1AS_fin <- merge(df_PHAG1AS, keggkotable_raw_taxa[, c("taxa","koID")])
write_tsv(df_PHAG1AS_fin, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/PHAG1AS_KEGGko_ancombc_05082025.tsv")
PHAG1AS_KEGGko <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/PHAG1AS_KEGGko_ancombc_05082025.tsv")


CAS3AS_KEGGko_unique <- CAS3AS_KEGGko[!CAS3AS_KEGGko$koID %in% PHAG1AS_KEGGko$koID,]
write_tsv(CAS3AS_KEGGko_unique, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAS3AS_KEGGko_unique_ancombc_05212025.tsv")
PHAG1AS_KEGGko_unique <- PHAG1AS_KEGGko[!PHAG1AS_KEGGko$koID %in% CAS3AS_KEGGko$koID,]
write_tsv(PHAG1AS_KEGGko_unique, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/PHAG1AS_KEGGko_unique_ancombc_05212025.tsv")


CAPHAunique <- CAS3AS_KEGGko[CAS3AS_KEGGko$koID %in% PHAG1AS_KEGGko$koID,]
PHACAunique <- PHAG1AS_KEGGko[PHAG1AS_KEGGko$koID %in% CAS3AS_KEGGko$koID,]
write_tsv(CAPHAunique, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/PHAG1AS_CAS3AS_shared_KEGGkosshared_ancombc_05212025.tsv")


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
