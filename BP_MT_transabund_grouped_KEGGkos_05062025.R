# Colleen Ahern
# 05/06/2025

# Magda asked that I provide her log fold change values for KEGG ko's instead of modules
# and to include taxonomic information

# We will use raw counts, not TPM normalized counts since I will be using ANCOM-BC
library(tidyverse)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)

## Count matrix input
bpcoldata <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/bp_metadata_01132025.txt")
bpcoldata$name <- gsub("-", ".", bpcoldata$name)

bpcts <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/read_counts_gene_metassembly.tsv")

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
bpcts_geneanno2 <- bpcts_geneanno
bpcts_geneanno2$eggNOG_OGs2 <- bpcts_geneanno2$eggNOG_OGs
bpcts_geneanno2 <- bpcts_geneanno2[,c(1:25,44,26:43)]
bpcts_geneanno2$eggNOG_OGs2 <- gsub(",.*?\\|", "", bpcts_geneanno2$eggNOG_OGs)
bpcts_geneanno2$eggNOG_OGs2 <- gsub(".*root", "", bpcts_geneanno2$eggNOG_OGs2)
bpcts_geneanno2$eggNOG_OGs2 <- gsub("([[:lower:]])([[:upper:]])", "\\1 \\2", bpcts_geneanno2$eggNOG_OGs2)
bpcts_geneanno2$eggNOG_OGs2 <- gsub("unclassified", " unclassified", bpcts_geneanno2$eggNOG_OGs2)
bpcts_geneanno2$KEGG_ko <- gsub("ko:","",bpcts_geneanno2$KEGG_ko)

# let's now do the KEGG ko annotations
aa <- bpcts_geneanno2 %>%
  dplyr::select(GeneID, everything())%>%
  mutate(KEGG_ko=str_split(KEGG_ko, ",")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:56,1)]

unique(bpcts_geneanno2$KEGG_ko)
ab<-as.vector(as.matrix(aa[,grepl("KEGG_ko",names(aa))]))
unique(ab)
length(unique(ab))

# make subset of aa data
aasub <- as.data.frame(aa[1:10,])
aasub$KEGG_ko1[is.na(aasub$KEGG_ko1)] <- "-"
aasubccs <- aasub[,c(32:44)]

uab <-unique(ab)
uab
rn <- uab[!is.na(uab)]
rn

keggkotable_raw <- data.frame(matrix(0, nrow = length(unique(rn)), ncol = 20))
colnames(keggkotable_raw) <- colnames(aa[,1:20])
rownames(keggkotable_raw) <- rn
keggkotable_raw$taxa <- NA

# Trial run on subset of the data
for (j in 1:(ncol(keggkotable_raw)-1)) {
  for (k in 1:length(aasub[,grepl("KEGG_ko" , names(aasub))])) { # Number of separate KEGG kos broken up
    for (i in 1:nrow(aasub)) {
      if(is.na(aasub[i,paste("KEGG_ko",k, sep = "")]) == FALSE) {
        keggkotable_raw[rn[rownames(keggkotable_raw) %in% aasub[i,paste("KEGG_ko",k, sep = "")]],j] = keggkotable_raw[rn[rownames(keggkotable_raw) %in% aasub[i,paste("KEGG_ko",k, sep = "")]],j] + aasub[i,j]
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
      keggkotable_raw[rn[rownames(keggkotable_raw) %in% aasub[i,paste("KEGG_ko",k, sep = "")]],"taxa"] <- paste(keggkotable_raw[rn[rownames(keggkotable_raw) %in% aasub[i,paste("KEGG_ko",k, sep = "")]],"taxa"],"; ", aasub$eggNOG_OGs2[i], sep = "")
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

keggkotable_raw <- data.frame(matrix(0, nrow = length(unique(rn)), ncol = 20))
colnames(keggkotable_raw) <- colnames(aa[,1:20])
rownames(keggkotable_raw) <- rn
keggkotable_raw$taxa <- NA

for (j in 1:(ncol(keggkotable_raw)-1)) {
  for (k in 1:length(aa[,grepl("KEGG_ko" , names(aa))])) { # Number of separate KEGG kos broken up
    for (i in 1:nrow(aa)) {
      if(is.na(aa[i,paste("KEGG_ko",k, sep = "")]) == FALSE) {
        keggkotable_raw[rn[rownames(keggkotable_raw) %in% aa[i,paste("KEGG_ko",k, sep = "")]],j] = keggkotable_raw[rn[rownames(keggkotable_raw) %in% aa[i,paste("KEGG_ko",k, sep = "")]],j] + aa[i,j]
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
      keggkotable_raw[rn[rownames(keggkotable_raw) %in% aa[i,paste("KEGG_ko",k, sep = "")]],"taxa"] <- paste(keggkotable_raw[rn[rownames(keggkotable_raw) %in% aa[i,paste("KEGG_ko",k, sep = "")]],"taxa"],"; ", aa$eggNOG_OGs2[i], sep = "")
    }
    else {
      keggkotable_raw = keggkotable_raw
    }
  }
}

keggkotable_raw$koID <- rownames(keggkotable_raw) 
write_tsv(keggkotable_raw, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggkotable_raw_05072025.tsv')
keggkotable_raw_r <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggkotable_raw_05072025.tsv')
rownames(keggkotable_raw_r) <- keggkotable_raw_r$koID
keggkotable_raw_r <- as.data.frame(keggkotable_raw_r)

# drop any outlier samples (based on PCA analysis)
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
keggkotable_raw_filt <- keggkotable_raw_r[,!(names(keggkotable_raw_r) %in% drops)]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(keggkotable_raw_filt),]

#######################################################################################################################################
#### CA + S3 + AS vs. CA + AS comparison

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "CA + AS" | bpcoldata3$group == "CA + S3 + AS",]
drop <- "-"

keggkotable_raw_filt_sub <- keggkotable_raw_r[,colnames(keggkotable_raw_r) %in% rownames(bpcoldata3_sub)]
rownames(keggkotable_raw_filt_sub) <- rownames(keggkotable_raw_r)
keggkotable_raw_filt_sub_noNA <- keggkotable_raw_filt_sub[!(row.names(keggkotable_raw_filt_sub) %in% drop), ]
keggkotable_raw_taxa <- keggkotable_raw_r[!(row.names(keggkotable_raw_filt_sub) %in% drop),]

colnames(keggkotable_raw_filt_sub) == rownames(bpcoldata3_sub)
colnames(keggkotable_raw_filt_sub_noNA) == rownames(bpcoldata3_sub)


##############################################################################################################################
# Now do for PHA + G1 + AS vs. PHA + AS comparison

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "PHA + G1 + AS" | bpcoldata3$group == "PHA + AS",]
drop <- "-"

keggkotable_raw_filt_sub <- keggkotable_raw_r[,colnames(keggkotable_raw_r) %in% rownames(bpcoldata3_sub)]
rownames(keggkotable_raw_filt_sub) <- rownames(keggkotable_raw_r)
keggkotable_raw_filt_sub_noNA <- keggkotable_raw_filt_sub[!(rownames(keggkotable_raw_filt_sub) %in% drop), ]
keggkotable_raw_taxa <- keggkotable_raw_r[!(row.names(keggkotable_raw_filt_sub) %in% drop),]

colnames(keggkotable_raw_filt_sub) == rownames(bpcoldata3_sub)
colnames(keggkotable_raw_filt_sub_noNA) == rownames(bpcoldata3_sub)


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

otumat = as.matrix(keggkotable_raw_filt_sub_noNA[,!(names(keggkotable_raw_filt_sub_noNA) %in% drop)])

taxmat = matrix(sample(letters, 62, replace = TRUE), nrow = nrow(otumat), ncol = 8)
rownames(taxmat) <- rownames(otumat)
colnames(taxmat) <- c("Domain", "Phylum", "Class", "Order", "Family", "Genus", "Species", "Taxa")
taxmat[,"Genus"] <- rownames(taxmat)
all(keggkotable_raw_taxa$koID == rownames(taxmat))
all(rownames(keggkotable_raw_taxa) == rownames(taxmat))
taxmat[,"Taxa"] <- keggkotable_raw_taxa[, "taxa"]
taxmat[,"Taxa"] <- gsub("NA; ", "", taxmat[,"Taxa"])
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
col_name = c("Taxon", "Intercept", "CA + S3 + AS")
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
  datatable(caption = "Differentially Abundant KOs from the Primary Result")

cadiff <- tab_diff[tab_diff$`CA + S3 + AS` == TRUE,]
phadiff <- tab_diff[tab_diff$`PHA + G1 + AS` == TRUE,]


tab_lfc_filt <- tab_lfc[tab_lfc$Taxon %in% cadiff$Taxon,]
colnames(tab_lfc_filt) <- c("Taxon", "Intercept LFC", "CA + S3 + AS LFC")
tab_lfc_filt <- tab_lfc[tab_lfc$Taxon %in% phadiff$Taxon,]
colnames(tab_lfc_filt) <- c("Taxon", "Intercept LFC", "PHA + G1 + AS LFC")

tab_se_filt <- tab_se[tab_se$Taxon %in% cadiff$Taxon,]
colnames(tab_se_filt) <- c("Taxon", "Intercept SE", "CA + S3 + AS SE")
tab_se_filt <- tab_se[tab_se$Taxon %in% phadiff$Taxon,]
colnames(tab_se_filt) <- c("Taxon", "Intercept SE", "PHA + G1 + AS SE")

tab_w_filt <- tab_w[tab_w$Taxon %in% cadiff$Taxon,]
colnames(tab_w_filt) <- c("Taxon", "Intercept W", "CA + S3 + AS W")
tab_w_filt <- tab_w[tab_w$Taxon %in% phadiff$Taxon,]
colnames(tab_w_filt) <- c("Taxon", "Intercept W", "PHA + G1 + AS W")

tab_p_filt <- tab_p[tab_p$Taxon %in% cadiff$Taxon,]
colnames(tab_p_filt) <- c("Taxon", "Intercept P", "CA + S3 + AS P")
tab_p_filt <- tab_p[tab_p$Taxon %in% phadiff$Taxon,]
colnames(tab_p_filt) <- c("Taxon", "Intercept P", "PHA + G1 + AS P")

tab_q_filt <- tab_q[tab_q$Taxon %in% cadiff$Taxon,]
colnames(tab_q_filt) <- c("Taxon", "Intercept Q", "CA + S3 + AS Q")
tab_q_filt <- tab_q[tab_q$Taxon %in% phadiff$Taxon,]
colnames(tab_q_filt) <- c("Taxon", "Intercept Q", "PHA + G1 + AS Q")

tab_diff_filt <- tab_diff[tab_diff$Taxon %in% cadiff$Taxon,]
colnames(tab_diff_filt) <- c("Taxon", "Intercept diff", "CA + S3 + AS diff")
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
  dplyr::filter(`CA + S3 + AS` == TRUE) %>%
  .$Taxon
sig_taxa = tab_diff %>%
  dplyr::filter(`PHA + G1 + AS` == TRUE) %>%
  .$Taxon

tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "CA + S3 + AS vs. CA + AS")
col_name = c("Taxon", "Intercept", "PHA + G1 + AS vs. PHA + AS")
colnames(tab_lfc) = col_name
tab_lfc <- tab_lfc[,-2]

df_CAS3AS = tab_lfc %>%
  filter(Taxon %in% sig_taxa) %>%
  dplyr::rename(koID = Taxon)

keggkotable_raw_taxa$taxa <- gsub("NA; ","",keggkotable_raw_taxa$taxa)
df_CAS3AS_fin <- merge(df_CAS3AS, keggkotable_raw_taxa[, c("taxa","koID")])
write_tsv(df_CAS3AS_fin, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAS3AS_KEGGko_ancombc_05082025.tsv")
CAS3AS_KEGGko <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAS3AS_KEGGko_ancombc_05082025.tsv")

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
