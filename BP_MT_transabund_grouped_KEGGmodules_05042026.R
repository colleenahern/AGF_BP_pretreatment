# Colleen Ahern
# 05042026

# Get differential abundance changes for KEGG modules instead of pathways
# use Maaslin2 for differential analysis
# and to include taxonomic information and make heatmaps
# Based on BP_MT_transabund_grouped_KEGGpathways_12162025.R
# Using correct taxa
# Need to use tax_gtdb_metassembly.tsv that links contigs to taxa
# Need to use bowtie2_metaassembly.stranded2.counts.txt to link genes to contigs
# The same way I fixed this for the CAZyme subgroups
# We will use raw counts

library(tidyverse)
library(readr)
library(dplyr)
library(stringr)
library(tidyr)
library(Maaslin2)
library(pheatmap)

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
all(bpcts$GeneID %in% contig_gene2$Geneid)
all(contig_gene2$Geneid %in% bpcts$GeneID)
bpcts <- bpcts[order(match(bpcts$GeneID, contig_gene2$Geneid)), ]
all(bpcts[,7:26] == contig_gene2[,7:26])

# ok this confirmed that the bowtie2_metassembly.stranded2.counts.txt counts are the same as the read_counts_gene_metassembly.tsv, phew
# ok now repeat analysis using contig_gene2 instead of bpcts cus they have the same info, but contig_gene2 has the taxa info added on

## Count matrix input
bpcoldata <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/bp_metadata_01132025.txt")
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
all(colnames(bpcts) %in% rownames(bpcoldata2))
all(rownames(bpcoldata2) == colnames(bpcts))

bpcts <- bpcts[, rownames(bpcoldata2)]
all(rownames(bpcoldata2) == colnames(bpcts))
bpcts$GeneID <- rownames(bpcts)

# Input gene annotation info
gene_anno <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/gene_annotation_metassembly.tsv")

bpcts_geneanno <- merge(bpcts, gene_anno, by = "GeneID", all = TRUE)
bpcts_geneanno2 <- merge(bpcts_geneanno, tax_anno, by = "GeneID", all = TRUE)

# let's now do the KEGG module annotations
aa <- bpcts_geneanno2 %>%
  dplyr::select(GeneID, everything())%>%
  mutate(KEGG_Module=str_split(KEGG_Module, ",")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:68,1)]

unique(bpcts_geneanno2$KEGG_Module)
length(unique(bpcts_geneanno2$KEGG_Module))
ab<-as.vector(as.matrix(aa[,grepl("KEGG_Module",names(aa))]))
unique(ab)
length(unique(ab))

# make subset of aa data
bb <- aa[is.na(aa$KEGG_Module2) == FALSE,]
aasub <- as.data.frame(bb[1:5,])
aasub$KEGG_Module1[is.na(aasub$KEGG_Module1) == TRUE] <- "-"
aasubccs <- aasub[,c(33:39,60)]

uab <-unique(ab)
uab
length(uab)
rn <- uab[!is.na(uab)]
rn
length(rn)

keggmod_raw <- data.frame(matrix(0, nrow = length(rn), ncol = 20))
colnames(keggmod_raw) <- colnames(aasub[,1:20])
rownames(keggmod_raw) <- rn
keggmod_raw$taxa <- NA

# Trial run on subset of the data
for (j in 1:(ncol(keggmod_raw)-1)) {
  for (k in 1:length(names(aasub)[grepl("KEGG_Module" , names(aasub))])) { # Number of separate KEGG pathways broken up
    for (i in 1:nrow(aasub)) {
      if(is.na(aasub[i,paste("KEGG_Module",k, sep = "")]) == FALSE) {
        keggmod_raw[rn[rn %in% aasub[i,paste("KEGG_Module",k, sep = "")]],j] = keggmod_raw[rn[rn %in% aasub[i,paste("KEGG_Module",k, sep = "")]],j] + as.numeric(aasub[i,j])
      }
      else {
        keggmod_raw = keggmod_raw
      }
    }
  }
}

for (k in 1:length(aasub[,grepl("KEGG_Module" , names(aasub))])) { # Number of separate KEGG pathways broken up
  for (i in 1:nrow(aasub)) {
    if(aasub[i,paste("KEGG_Module",k, sep = "")] %in% rn == TRUE ) {
      keggmod_raw[rn[rn %in% aasub[i,paste("KEGG_Module",k, sep = "")]],"taxa"] <- paste(keggmod_raw[rn[rn %in% aasub[i,paste("KEGG_Module",k, sep = "")]],"taxa"],"; ", aasub$taxa[i], sep = "")
    }
    else {
      keggmod_raw = keggmod_raw
    }
  }
}


## Now do it on the entire data
aa$KEGG_Module1[is.na(aa$KEGG_Module1) == TRUE] <- "-"

uab <-unique(ab)
uab
length(uab)
rn <- uab[!is.na(uab)]
rn
length(rn)

keggmod_raw <- data.frame(matrix(0, nrow = length(rn), ncol = 20))
colnames(keggmod_raw) <- colnames(aa[,1:20])
rownames(keggmod_raw) <- rn
keggmod_raw$taxa <- NA

for (j in 1:(ncol(keggmod_raw)-1)) {
  for (k in 1:length(names(aa)[grepl("KEGG_Module" ,names(aa))])) { # Number of separate KEGG modules broken up
    for (i in 1:nrow(aa)) {
      if(is.na(aa[i,paste("KEGG_Module",k, sep = "")]) == FALSE) {
        keggmod_raw[rn[rn %in% aa[i,paste("KEGG_Module",k, sep = "")]],j] = keggmod_raw[rn[rn %in% aa[i,paste("KEGG_Module",k, sep = "")]],j] + as.numeric(aa[i,j])
      }
      else {
        keggmod_raw = keggmod_raw
      }
    }
  }
}

keggmod_raw2 <- keggmod_raw %>%
  rownames_to_column("KEGG_Module")
write_tsv(keggmod_raw2, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/KEGGmod_raw_notaxaanno_05042026.tsv")
keggmod_rawr <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/KEGGmod_raw_notaxaanno_05042026.tsv")

for (k in 1:length(aa[,grepl("KEGG_Module" ,names(aa))])) { # Number of separate KEGG modules broken up
  for (i in 1:nrow(aa)) {
    if(aa[i,paste("KEGG_Module",k, sep = "")] %in% rn == TRUE ) {
      keggmod_raw[rn[rn %in% aa[i,paste("KEGG_Module",k, sep = "")]],"taxa"] <- paste(keggmod_raw[rn[rn %in% aa[i,paste("KEGG_Module",k, sep = "")]],"taxa"],"; ", aa$taxa[i], sep = "")
    }
    else {
      keggmod_raw = keggmod_raw
    }
  }
}

# 05/05/2026 didn;t complete took too long on mac, do on pod later
keggmod_raw2 <- keggmod_raw %>%
  rownames_to_column("KEGG_Pathway")
write_tsv(keggmod_raw2, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggmodway_raw_taxaanno_0502026.tsv") # redid to make sure it's right - pretty sure its the exact same as 1219
keggmod_rawr <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggmodway_raw_taxaanno_05022026.tsv")

# keggmod_rawr <- keggmod_rawr %>%
#   column_to_rownames("KEGG_Module")

# Remove the results for "-" since I think this is crashing my computer, plus I don't care about the taxa associated with non-CAZymes
drop <- c("-")
keggmod_rawrr <- keggmod_rawr[!(keggmod_rawr$KEGG_Module %in% drop),]
keggmod_rawr <- keggmod_rawrr

# keggmod_rawrrr <- keggmod_rawrr
# keggmod_rawrrr$taxa <- substr(keggmod_rawrrr$taxa, 4, nchar(keggmod_rawrrr$taxa))
# write_tsv(keggmod_rawrrr, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggmodway_raw_taxaanno_NOunknown_noNAstart_05022026.tsv')
# keggmod_rawr <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggmodway_raw_taxaanno_NOunknown_noNAstart_05022026.tsv') # confirmed all(keggmod_rawrt_filt == keggmod_rawrorig_filt)
# keggmod_rawr <- as.data.frame(keggmod_rawrrr)

# taxaf <- keggmod_rawr[,c("KEGG_Pathway","taxa")]
# taxaf <- taxaf %>%
#   column_to_rownames("KEGG_Pathway")

keggmod_rawrt <- keggmod_rawr %>%
  column_to_rownames("KEGG_Module")

# drop any outlier samples (based on PCA analysis)
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC","taxa")
keggmod_rawrt_filt <- keggmod_rawrt[,!(names(keggmod_rawrt) %in% drops)]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(keggmod_rawrt_filt),]

# Run Maaslin2 instead of ANCOMBC
library(Maaslin2)
library(dplyr)
library(readr)

bpcoldata3$group <- gsub(" + ", "_", bpcoldata3$group, fixed = TRUE)
bpcoldata3$name <- row.names(bpcoldata3)

# Transpose so rows = samples, columns = features
keggmod_rawrt_filt_t <- as.data.frame(t(keggmod_rawrt_filt))
all(rownames(keggmod_rawrt_filt_t) == rownames(bpcoldata3))

# --- Comparison 1: CA_S3_AS vs CA_AS ---
meta1 <- bpcoldata3 |> filter(group %in% c("CA_S3_AS", "CA_AS"))
data1 <- keggmod_rawrt_filt_t[rownames(meta1), ]

fit1 <- Maaslin2(
  input_data      = data1,
  input_metadata  = meta1,
  output          = "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/KEGG_module_maaslin2_CA_05052026",
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
data2 <- keggmod_rawrt_filt_t[rownames(meta2), ]

fit2 <- Maaslin2(
  input_data      = data2,
  input_metadata  = meta2,
  output          = "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/KEGG_module_maaslin2_PHA_05022026",
  fixed_effects   = "group",
  reference       = "group,PHA_AS",
  normalization   = "TSS",
  transform       = "LOG",
  analysis_method = "LM",
  min_prevalence  = 0.1,
  min_abundance   = 0.0,
  cores           = 1
)

# --- Comparison 3: PHA_G1_AS vs PHA_AS ---
bpcoldata3$condition_cba <- gsub(" + ", "_", bpcoldata3$condition_cba, fixed = TRUE)
bpcoldata3$condition_cba <- gsub(" ", "_", bpcoldata3$condition_cba)

meta3 <- bpcoldata3 |> filter(condition_cba %in% c("PHA_G1_AS_t2", "PHA_G1_AS_t1"))
data3 <- keggmod_rawrt_filt_t[rownames(meta3), ]

fit3 <- Maaslin2(
  input_data      = data3,
  input_metadata  = meta3,
  output          = "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/KEGG_module_maaslin2_PHAG1ASt2t1_05052026",
  fixed_effects   = "condition_cba",
  reference       = "condition_cba,PHA_G1_AS_t1",
  normalization   = "TSS",
  transform       = "LOG",
  analysis_method = "LM",
  min_prevalence  = 0.1,
  min_abundance   = 0.0,
  cores           = 1
)

# make heatmap
library(dplyr)
library(pheatmap)

# ── Extract MaAsLin2 results ──────────────────────────────────────────────────
res_CA  <- fit1$results
res_PHA <- fit2$results
res_PHAb <- fit3$results

# ── Build LFC and q-value matrices ───────────────────────────────────────────
all_cogs <- union(res_CA$feature, union(res_PHA$feature, res_PHAb$feature))

lfc_CA  <- setNames(res_CA$coef,  res_CA$feature)
lfc_PHA <- setNames(res_PHA$coef, res_PHA$feature)
lfc_PHAb <- setNames(res_PHAb$coef, res_PHAb$feature)
q_CA    <- setNames(res_CA$qval,  res_CA$feature)
q_PHA   <- setNames(res_PHA$qval, res_PHA$feature)
q_PHAb   <- setNames(res_PHAb$qval, res_PHAb$feature)

plot_mat <- data.frame(
  CA  = lfc_CA[all_cogs],
  PHA = lfc_PHA[all_cogs],
  PHAb = lfc_PHAb[all_cogs],
  row.names = all_cogs
)
plot_mat[is.na(plot_mat)] <- 0

padj_mat <- data.frame(
  CA  = q_CA[all_cogs],
  PHA = q_PHA[all_cogs],
  PHAb = q_PHAb[all_cogs],
  row.names = all_cogs
)
padj_mat[is.na(padj_mat)] <- 1

combined_mat <- data.frame(
  lfc_CA  = lfc_CA[all_cogs],
  q_CA    = q_CA[all_cogs],
  lfc_PHA = lfc_PHA[all_cogs],
  q_PHA   = q_PHA[all_cogs],
  lfc_PHAb = lfc_PHAb[all_cogs],
  q_PHAb   = q_PHAb[all_cogs],
  row.names = all_cogs
)

combined_mat$lfc_CA[is.na(combined_mat$lfc_CA)]   <- 0
combined_mat$lfc_PHA[is.na(combined_mat$lfc_PHA)] <- 0
combined_mat$lfc_PHAb[is.na(combined_mat$lfc_PHAb)] <- 0
combined_mat$q_CA[is.na(combined_mat$q_CA)]   <- 1
combined_mat$q_PHA[is.na(combined_mat$q_PHA)] <- 1
combined_mat$q_PHAb[is.na(combined_mat$q_PHAb)] <- 1
combined_mat$KEGG_Module <- rownames(combined_mat)
write_tsv(combined_mat, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/KEGG_module_maaslin2_CAPHA_coeffpadj_05052026.tsv")
write_tsv(combined_mat, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/KEGG_module_maaslin2_CAPHAPHAb_coeffpadj_05052026.tsv")
combined_matr <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/KEGG_module_maaslin2_CAPHA_coeffpadj_05052026.tsv")
combined_matrb <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/KEGG_module_maaslin2_CAPHAPHAb_coeffpadj_05052026.tsv")

# make KEGG module definition table
library(KEGGREST)

module_ids <- names(keggList("module"))

# KEGGREST limits to 10 per query — loop in chunks
results <- list()
chunks <- split(module_ids, ceiling(seq_along(module_ids)/10))

for (i in seq_along(chunks)) {
  cat("Fetching chunk", i, "of", length(chunks), "\n")
  res <- keggGet(chunks[[i]])
  results <- c(results, res)
  Sys.sleep(0.3)  # be polite to the API
}

# Parse into table
module_table <- do.call(rbind, lapply(results, function(x) {
  class_split <- strsplit(x$CLASS, "; ")[[1]]
  data.frame(
    module_id   = x$ENTRY,
    module_name = x$NAME,
    category1   = ifelse(length(class_split) >= 1, class_split[1], NA),
    category2   = ifelse(length(class_split) >= 2, class_split[2], NA),
    category3   = ifelse(length(class_split) >= 3, class_split[3], NA),
    stringsAsFactors = FALSE
  )
}))

module_table <- module_table %>%
  rename(KEGG_Module = module_id)
write.table(module_table, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/KEGG_modules_list_with_categories_05052026.tsv", sep="\t", row.names=FALSE, quote=FALSE)
modtab <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/KEGG_modules_list_with_categories_05052026.tsv")

combmodmat <- merge(combined_matr, modtab, all.x = TRUE)
combmodmat_filt <- combmodmat[combmodmat$q_CA < 0.05 | combmodmat$q_PHA < 0.05,]
combmodmat_filt2 <- combmodmat_filt[is.na(combmodmat_filt$module_name) == FALSE,]
combmodmat_filt2$Finalname <- paste(combmodmat_filt2$KEGG_Module, " - ", combmodmat_filt2$module_name, sep = "")

combmodmatb <- merge(combined_matrb, modtab, all.x = TRUE)
combmodmatb_filt <- combmodmatb[combmodmatb$q_CA < 0.05 | combmodmatb$q_PHA < 0.05 | combmodmatb$q_PHAb < 0.05,]
combmodmatb_filt2 <- combmodmatb_filt[is.na(combmodmatb_filt$module_name) == FALSE,]

df <- combmodmat_filt2
df <- combmodmatb_filt2

# Build LFC matrix (rows = modules, cols = comparisons)
plot_mat <- data.frame(
  CA  = df$lfc_CA,
  PHA = df$lfc_PHA,
  PHAb = df$lfc_PHAb,
  row.names = df$module_name
)

# Build q-value matrix
padj_mat <- data.frame(
  CA  = df$q_CA,
  PHA = df$q_PHA,
  PHAb = df$q_PHAb,
  row.names = df$module_name
)

# Mask non-significant values
plot_mat_masked <- plot_mat
plot_mat_masked[padj_mat > 0.05] <- 0

# ── Column names ──────────────────────────────────────────────────────────────
colnames(plot_mat_masked) <- c("CA + S3 + AS vs.\nCA + AS",
                               "PHA + G1 + AS vs.\nPHA + AS",
                               "PHA + G1 + AS t2 vs.\nPHA + G1 + AS t1")

# ── Row clustering on unmasked values ────────────────────────────────────────
# row_order <- hclust(dist(plot_mat), method = "ward.D2")

# ── Row annotations ───────────────────────────────────────────────────────────
row_annot <- data.frame(
  Category    = df$category2,
  Subcategory = df$category3,
  row.names   = df$module_name
)

# ── Annotation colors ─────────────────────────────────────────────────────────
# Category2 colors
cat2_vals <- unique(df$category2)
cat2_palette <- c(
  "Carbohydrate metabolism"                    = "#5BA75B",
  "Amino acid metabolism"                      = "#4E9BB5",
  "Nucleotide metabolism"                      = "#E07B39",
  "Lipid metabolism"                           = "#9B6BB5",
  "Glycan metabolism"                          = "#D4A84B",
  "Energy metabolism"                          = "#C0392B",
  "Metabolism of cofactors and vitamins"       = "#1ABC9C",
  "Xenobiotics biodegradation"                 = "#7F8C8D",
  "Biosynthesis of other secondary metabolites"= "#E91E8C",
  "Biosynthesis of terpenoids and polyketides" = "#FF6B35",
  "Gene set"                                   = "#BDC3C7"
)
# Keep only colors for categories present in data
cat2_colors <- cat2_palette[names(cat2_palette) %in% cat2_vals]
# Auto-assign any missing categories
missing_cat2 <- setdiff(cat2_vals, names(cat2_palette))
if (length(missing_cat2) > 0) {
  extra <- setNames(rainbow(length(missing_cat2)), missing_cat2)
  cat2_colors <- c(cat2_colors, extra)
}

# Category3 colors — auto-generate from unique values
cat3_vals <- unique(df$category3)
cat3_colors <- setNames(
  colorRampPalette(c("#A8D8EA","#AA96DA","#FCBAD3","#FFFFD2",
                     "#B5EAD7","#FFD700","#FF9AA2","#C7CEEA"))(length(cat3_vals)),
  cat3_vals
)

annotation_colors <- list(
  Category    = cat2_colors,
  Subcategory = cat3_colors
)

# ── Heatmap ───────────────────────────────────────────────────────────────────
# Build colors directly from what's in the data
cat2_vals <- unique(as.character(df$category2))
cat2_vals <- cat2_vals[!is.na(cat2_vals)]

cat2_colors <- setNames(
  colorRampPalette(c("#5BA75B","#4E9BB5","#E07B39","#9B6BB5",
                     "#D4A84B","#C0392B","#1ABC9C","#7F8C8D",
                     "#E91E8C","#FF6B35","#BDC3C7"))(length(cat2_vals)),
  cat2_vals
)

row_annot <- data.frame(
  Category = as.character(df$category2),
  row.names = df$module_name
)

# Rename column to whatever you want the legend title to be
row_annot <- data.frame(
  "KEGG Module Category" = as.character(df$category2),  
  row.names = df$module_name,
  check.names = FALSE  # prevents R from replacing spaces with dots
)

# Update annotation colors list name to match
annotation_colors <- list("KEGG Category" = cat2_colors)

modph <- pheatmap(plot_mat_masked,
                  cluster_cols          = FALSE,
                  cluster_rows          = TRUE,
                  color                 = colorRampPalette(c("blue", "black", "red"))(100),
                  breaks                = seq(-2, 2, length.out = 101),
                  border_color          = "grey80",
                  annotation_row        = row_annot,
                  annotation_colors     = list(Category = cat2_colors),
                  annotation_names_row  = FALSE,
                  fontsize              = 12,
                  fontsize_row          = 10,
                  angle_col             = 0,
                  cellwidth             = 120,
                  cellheight            = 15,
                  treeheight_row        = 50)
modph

library(ggplot2)
library(ggpubr)
library(ggplotify)
library(cowplot)
library(grid)

# Convert pheatmaps to ggplot objects
cog_gg <- as.ggplot(cogph) +
  theme(plot.margin = margin(t = 1, b = 0, l = -534, r = 0))

mod_gg <- as.ggplot(modph) +
  theme(plot.margin = margin(t = 15, b = 10, l = 10, r = 0))

# Stack and label
combined <- plot_grid(cog_gg, mod_gg,
                      labels = c("A)", "B)"),
                      ncol = 1,
                      rel_heights = c(1, 2.6))  # adjust ratio based on row counts

combined

ggsave(plot = combined, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/KEGGmodcog_maaslin2_heatmap_05062026.pdf", width = 20, height = 20.5, units = "in", dpi = 300, bg = "white")
ggsave(plot = combined, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/KEGGmodcog_maaslin22_heatmap_05062026.tiff", width = 20, height = 20.5, units = "in", dpi = 200, bg="white")



                               
