# Colleen Ahern
# 01/06/2025

# Continuation of BP_MT_transabund_grouped_CAZyme_05192025.R script
# I want to try normalizing using Natalia Ivanova's method from JGI:
# dividing gene counts by the total amount of gene counts of that group (CAZymes)

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
# dbcantable_gen_rawrt <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_raw_taxa_05202025.tsv')  # don't use this table, has the right counts but wrong taxa
dbcantable_gen_rawrt <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_raw_CORRECT_taxa_10102025.tsv') # use this, has right counts and right taxa
dbcantable_gen_rawrt <- dbcantable_gen_rawrt %>%
  column_to_rownames("CAZID")

# drop any outlier samples (based on PCA analysis)
# drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC") 12/11/2025 let's try not dropping them...I don't think we have enough reason to but check with Magda
# dbcantable_gen_rawrt_filt <- dbcantable_gen_rawrt[,!(names(dbcantable_gen_rawrt) %in% drops)]
dbcantable_gen_rawrt_filt <- dbcantable_gen_rawrt
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(dbcantable_gen_rawrt_filt),]

#######################################################################################################################################
#### CA + S3 + AS vs. CA + AS comparison

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "CA + AS" | bpcoldata3$group == "CA + S3 + AS",]
drop <- "-"

dbcantable_raw_filt_sub <- dbcantable_gen_rawrt_filt[,colnames(dbcantable_gen_rawrt_filt) %in% rownames(bpcoldata3_sub)]
dbcantable_raw_filt_sub_noNA <- dbcantable_raw_filt_sub[!(row.names(dbcantable_raw_filt_sub) %in% drop), ]
dbcantable_raw_taxa <- dbcantable_gen_rawrt_filt[!(row.names(dbcantable_gen_rawrt_filt) %in% drop),]

all(colnames(dbcantable_raw_filt_sub) == rownames(bpcoldata3_sub))
all(colnames(dbcantable_raw_filt_sub_noNA) == rownames(bpcoldata3_sub))

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
# ANCOM-BC2

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
dbcantable_raw_taxa$taxaf <- gsub("^.{0,3}", "", dbcantable_raw_taxa$taxaf)
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

out2 = ancombc2(data = physeq, tax_level = "Genus", 
                fix_formula = "group", 
                p_adj_method = "holm",
                alpha = 0.05,
                verbose = TRUE)

res2 = out2$res

## ANCOMBC primary result - skip if doing ancombc2 analysis 
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

# Make and save dataframe of results
tab_lfc = res$lfc
col_name = c("CAZyme", "Intercept", "CA + S3 + AS vs. CA + AS")
col_name = c("CAZyme", "Intercept", "PHA + G1 + AS vs. PHA + AS")
colnames(tab_lfc) = col_name
tab_lfc <- tab_lfc[,-2]

df_CAS3AS = tab_lfc %>%
  filter(CAZyme %in% sig_taxa)

df_CAS3AS_fin <- merge(df_CAS3AS, dbcantable_raw_taxa[, c("taxa","CAZyme")], by = "CAZyme")
write_tsv(df_CAS3AS_fin, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/CAS3AS_CAZyme_subs_ancombc_05212025.tsv")
bbb <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAZ_subs/CAS3AS_CAZyme_subs_ancombc_05212025.tsv")

df_PHAG1AS = tab_lfc %>%
  filter(CAZyme %in% sig_taxa) 

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
# Making volcano plots instead of heatmmaps to show differentially expressed CAZymes
# taken from https://biostatsquid.com/volcano-plots-r-tutorial/

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

# Import Ancombc results

# for CAS3AS vs CAAS comparison

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

## ^ ok so for the graph above I am labeling the 10 most significantly expressed (lowest padj values) upregulated CAZyme subgroup and 
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







