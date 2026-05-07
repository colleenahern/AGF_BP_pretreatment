# Colleen Ahern
# 12/16/2025

# Get log fold change values for KEGG pathways instead of KO's
# and to include taxonomic information and make volcano plots or heatmaps or barplots whatever
# Fixing my taxa associated with each group of interest (in this case KEGG kos) because Malte said I was using the wrong file in our meeting
# Need to use tax_gtdb_metassembly.tsv that links contigs to taxa
# Need to use bowtie2_metaassembly.stranded2.counts.txt to link genes to contigs
# The same way I fixed this for the CAZyme subgroups

# We will use raw counts if there are enough pathways to use ANCOM-BC. If not, use TPM counts
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
all(bpcts$GeneID %in% contig_gene2$Geneid)
all(contig_gene2$Geneid %in% bpcts$GeneID)
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
all(colnames(bpcts) %in% rownames(bpcoldata2))
all(rownames(bpcoldata2) == colnames(bpcts))

bpcts <- bpcts[, rownames(bpcoldata2)]
all(rownames(bpcoldata2) == colnames(bpcts))
bpcts$GeneID <- rownames(bpcts)

# Input gene annotation info
gene_anno <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/gene_annotation_metassembly.tsv")

bpcts_geneanno <- merge(bpcts, gene_anno, by = "GeneID", all = TRUE)
bpcts_geneanno2 <- merge(bpcts_geneanno, tax_anno, by = "GeneID", all = TRUE)
bpcts_geneanno2$KEGG_Pathway <- str_remove(bpcts_geneanno2$KEGG_Pathway, ",map.*")
bpcts_geneanno2$KEGG_Pathway <- gsub("ko", "map", bpcts_geneanno2$KEGG_Pathway)

library(dplyr)
library(stringr)
library(tidyr)

# let's now do the KEGG ko annotations
aa <- bpcts_geneanno2 %>%
  dplyr::select(GeneID, everything())%>%
  mutate(KEGG_Pathway=str_split(KEGG_Pathway, ",")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:136,1)]

unique(bpcts_geneanno2$KEGG_Pathway)
length(unique(bpcts_geneanno2$KEGG_Pathway))
ab<-as.vector(as.matrix(aa[,grepl("KEGG_Pathway",names(aa))]))
unique(ab)
length(unique(ab))

# make subset of aa data
bb <- aa[is.na(aa$KEGG_Pathway3) == FALSE,]
aasub <- as.data.frame(bb[1:5,])
aasub$KEGG_Pathway1[is.na(aasub$KEGG_Pathway1) == TRUE] <- "-"
aasubccs <- aasub[,c(32:41)]

uab <-unique(ab)
uab
length(uab)
rn <- uab[!is.na(uab)]
rn
length(rn)

keggpath_raw <- data.frame(matrix(0, nrow = length(rn), ncol = 20))
colnames(keggpath_raw) <- colnames(aasub[,1:20])
rownames(keggpath_raw) <- rn
keggpath_raw$taxa <- NA

# Trial run on subset of the data
for (j in 1:(ncol(keggpath_raw)-1)) {
  for (k in 1:length(names(aasub)[grepl("KEGG_Pathway" , names(aasub))])) { # Number of separate KEGG pathways broken up
    for (i in 1:nrow(aasub)) {
      if(is.na(aasub[i,paste("KEGG_Pathway",k, sep = "")]) == FALSE) {
        keggpath_raw[rn[rn %in% aasub[i,paste("KEGG_Pathway",k, sep = "")]],j] = keggpath_raw[rn[rn %in% aasub[i,paste("KEGG_Pathway",k, sep = "")]],j] + as.numeric(aasub[i,j])
      }
      else {
        keggpath_raw = keggpath_raw
      }
    }
  }
}

for (k in 1:length(aasub[,grepl("KEGG_Pathway" , names(aasub))])) { # Number of separate KEGG pathways broken up
  for (i in 1:nrow(aasub)) {
    if(aasub[i,paste("KEGG_Pathway",k, sep = "")] %in% rn == TRUE ) {
      keggpath_raw[rn[rn %in% aasub[i,paste("KEGG_Pathway",k, sep = "")]],"taxa"] <- paste(keggpath_raw[rn[rn %in% aasub[i,paste("KEGG_Pathway",k, sep = "")]],"taxa"],"; ", aasub$taxa[i], sep = "")
    }
    else {
      keggpath_raw = keggpath_raw
    }
  }
}


## Now do it on the entire data
aa$KEGG_Pathway1[is.na(aa$KEGG_Pathway1) == TRUE] <- "-"

uab <-unique(ab)
uab
length(uab)
rn <- uab[!is.na(uab)]
rn
length(rn)

keggpath_raw <- data.frame(matrix(0, nrow = length(rn), ncol = 20))
colnames(keggpath_raw) <- colnames(aa[,1:20])
rownames(keggpath_raw) <- rn
keggpath_raw$taxa <- NA

for (j in 1:(ncol(keggpath_raw)-1)) {
  for (k in 1:length(names(aa)[grepl("KEGG_Pathway" , names(aa))])) { # Number of separate KEGG pathways broken up
    for (i in 1:nrow(aa)) {
      if(is.na(aa[i,paste("KEGG_Pathway",k, sep = "")]) == FALSE) {
        keggpath_raw[rn[rn %in% aa[i,paste("KEGG_Pathway",k, sep = "")]],j] = keggpath_raw[rn[rn %in% aa[i,paste("KEGG_Pathway",k, sep = "")]],j] + as.numeric(aa[i,j])
      }
      else {
        keggpath_raw = keggpath_raw
      }
    }
  }
}

keggpath_raw2 <- keggpath_raw %>%
  rownames_to_column("KEGG_Pathway")
write_tsv(keggpath_raw2, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_notaxaanno_12182025.tsv")

for (k in 1:length(aa[,grepl("KEGG_Pathway" , names(aa))])) { # Number of separate KEGG pathways broken up
  for (i in 1:nrow(aa)) {
    if(aa[i,paste("KEGG_Pathway",k, sep = "")] %in% rn == TRUE ) {
      keggpath_raw[rn[rn %in% aa[i,paste("KEGG_Pathway",k, sep = "")]],"taxa"] <- paste(keggpath_raw[rn[rn %in% aa[i,paste("KEGG_Pathway",k, sep = "")]],"taxa"],"; ", aa$taxa[i], sep = "")
    }
    else {
      keggpath_raw = keggpath_raw
    }
  }
}

keggpath_raw2 <- keggpath_raw %>%
  rownames_to_column("KEGG_Pathway")
write_tsv(keggpath_raw2, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_12192025.tsv")
write_tsv(keggpath_raw2, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_0502026.tsv") # redid to make sure it's right - pretty sure its the exact same as 1219
keggpath_rawr <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_12192025.tsv")
keggpath_rawr <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_05022026.tsv")

keggpath_rawr <- keggpath_rawr %>%
  column_to_rownames("KEGG_Pathway")
keggpath_rawr <- keggpath_raw2

# Remove the results for "-" since I think this is crashing my computer, plus I don't care about the taxa associated with non-CAZymes
drop <- c("-")
keggpath_rawrr <- keggpath_rawr[!(keggpath_rawr$KEGG_Pathway %in% drop),]
write_tsv(keggpath_rawrr, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_NOunknown_12192025.tsv')
write_tsv(keggpath_rawrr, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_NOunknown_05022026.tsv')




keggpath_rawrrr <- keggpath_rawrr
keggpath_rawrrr$taxa <- substr(keggpath_rawrrr$taxa, 4, nchar(keggpath_rawrrr$taxa))
write_tsv(keggpath_rawrrr, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_NOunknown_noNAstart_12192025.tsv')
write_tsv(keggpath_rawrrr, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_NOunknown_noNAstart_05022026.tsv')
keggpath_rawr <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_NOunknown_noNAstart_12192025.tsv')
keggpath_rawr <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_NOunknown_noNAstart_05022026.tsv') # confirmed all(keggpath_rawrt_filt == keggpath_rawrorig_filt)
keggpath_rawr <- as.data.frame(keggpath_rawr)
keggpath_rawr <- as.data.frame(keggpath_rawrrr)

taxaf <- keggpath_rawr[,c("KEGG_Pathway","taxa")]
taxaf <- taxaf %>%
  column_to_rownames("KEGG_Pathway")

keggpath_rawrt <- keggpath_rawr %>%
  column_to_rownames("KEGG_Pathway")

# drop any outlier samples (based on PCA analysis)
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC","taxa")
keggpath_rawrt_filt <- keggpath_rawrt[,!(names(keggpath_rawrt) %in% drops)]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(keggpath_rawrt_filt),]

# do ancombc1 - 05/02/2026 checking to make sure this matches my original 1219 results

library(ANCOMBC)
library(phyloseq)

# ── Helper function ───────────────────────────────────────────────────────────
run_ancombc <- function(mat, metadata, cond1_labels, cond2_labels) {
  
  # Subset to relevant samples
  keep_samples <- metadata$name[metadata$group %in% c(cond1_labels, cond2_labels)]
  meta_sub     <- metadata[metadata$name %in% keep_samples, ]
  mat_sub      <- mat[, keep_samples]
  
  # Collapse to binary condition label for the contrast
  meta_sub$contrast <- ifelse(meta_sub$group %in% cond1_labels, "treatment", "control")
  meta_sub$contrast <- factor(meta_sub$contrast, levels = c("control", "treatment"))
  
  # Build phyloseq object
  otu  <- otu_table(mat_sub, taxa_are_rows = TRUE)
  samp <- sample_data(meta_sub)
  rownames(samp) <- meta_sub$name
  ps   <- phyloseq(otu, samp)
  
  # Run ANCOM-BC (original)
  out <- ancombc(
    data         = ps,
    formula      = "contrast",
    p_adj_method = "BH",
    prv_cut      = 0,
    verbose      = TRUE
  )
  
  return(out$res)
}

bpcoldata3$group <- gsub(" + ", "_", bpcoldata3$group, fixed = TRUE)
bpcoldata3$name <- row.names(bpcoldata3)

# ── Run two contrasts ─────────────────────────────────────────────────────────
res_S3 <- run_ancombc(keggpath_rawrt_filt, bpcoldata3,
                      cond1_labels = "CA_S3_AS",
                      cond2_labels = "CA_AS")

res_G1 <- run_ancombc(keggpath_rawrt_filt, bpcoldata3,
                      cond1_labels = "PHA_G1_AS",
                      cond2_labels = "PHA_AS")

# ── Extract results ───────────────────────────────────────────────────────────
# ANCOM-BC stores results differently from ANCOM-BC2
plot_mat <- cbind(
  S3 = res_S3$lfc,
  G1 = res_G1$lfc
)
rownames(plot_mat) <- plot_mat$S3.taxon
plot_mat <- plot_mat[,-c(1,2,4,5)]

padj_mat <- cbind(
  S3 = res_S3$q_val,
  G1 = res_G1$q_val
)
rownames(padj_mat) <- padj_mat$S3.taxon
padj_mat <- padj_mat[,-c(1,2,4,5)]

volc_mat <- cbind(
    S3_lfc = res_S3$lfc,
    S3_padj = res_S3$q_val,
    G1_lfc = res_G1$lfc,
    G1_padj = res_G1$q_val
  )

volc_mat <- volc_mat[,-c(2,4,5,7,8,10,11)]
names(volc_mat) <- gsub(".contrasttreatment", "", names(volc_mat))
volc_mat <- volc_mat %>%
  rename(KEGG_pathway = S3_lfc.taxon)
write_tsv(volc_mat, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/KEGGpathway_ancombc1_lfcpadj_05022026.tsv")

# ── Significance stars ────────────────────────────────────────────────────────
sig_mat <- ifelse(padj_mat < 0.05, "*", "")
sig_mat[padj_mat < 0.01]  <- "**"
sig_mat[padj_mat < 0.001] <- "***"

# Set non-significant to NA for black cells
plot_mat_na <- plot_mat
plot_mat_na[padj_mat > 0.05] <- NA

# ── Heatmap ───────────────────────────────────────────────────────────────────
library(pheatmap)

cog_definitions <- c(
  J = "J - Translation, ribosomal structure and biogenesis",
  A = "A - RNA processing and modification",
  K = "K - Transcription",
  L = "L - Replication, recombination and repair",
  B = "B - Chromatin structure and dynamics",
  D = "D - Cell cycle control, cell division, chromosome partitioning",
  Y = "Y - Nuclear structure",
  V = "V - Defense mechanisms",
  T = "T - Signal transduction mechanisms",
  M = "M - Cell wall/membrane/envelope biogenesis",
  N = "N - Cell motility",
  Z = "Z - Cytoskeleton",
  W = "W - Extracellular structures",
  U = "U - Intracellular trafficking, secretion, and vesicular transport",
  O = "O - Posttranslational modification, protein turnover, chaperones",
  X = "X - Mobilome: prophages, transposons",
  C = "C - Energy production and conversion",
  G = "G - Carbohydrate transport and metabolism",
  E = "E - Amino acid transport and metabolism",
  F = "F - Nucleotide transport and metabolism",
  H = "H - Coenzyme transport and metabolism",
  I = "I - Lipid transport and metabolism",
  P = "P - Inorganic ion transport and metabolism",
  Q = "Q - Secondary metabolites biosynthesis, transport and catabolism",
  R = "R - General function prediction only",
  S = "S - Function unknown"
)

rownames(plot_mat_na) <- cog_definitions[rownames(plot_mat_na)]
rownames(plot_mat)    <- cog_definitions[rownames(plot_mat)]

row_order <- hclust(dist(plot_mat), method = "ward.D2")

colnames(plot_mat_na) <- c("CA + S3 + AS vs.\nCA + AS",
                           "PHA + G1 + AS vs.\nPHA + AS")

cog_domains <- data.frame(
  Category = c(
    "J" = "Information storage & processing",
    "A" = "Information storage & processing",
    "K" = "Information storage & processing",
    "L" = "Information storage & processing",
    "B" = "Information storage & processing",
    "D" = "Cellular processes & signaling",
    "Y" = "Cellular processes & signaling",
    "V" = "Cellular processes & signaling",
    "T" = "Cellular processes & signaling",
    "M" = "Cellular processes & signaling",
    "N" = "Cellular processes & signaling",
    "Z" = "Cellular processes & signaling",
    "W" = "Cellular processes & signaling",
    "U" = "Cellular processes & signaling",
    "O" = "Cellular processes & signaling",
    "X" = "Cellular processes & signaling",
    "C" = "Metabolism",
    "G" = "Metabolism",
    "E" = "Metabolism",
    "F" = "Metabolism",
    "H" = "Metabolism",
    "I" = "Metabolism",
    "P" = "Metabolism",
    "Q" = "Metabolism",
    "R" = "Poorly characterized",
    "S" = "Poorly characterized"
  )
)

cog_letters <- substr(rownames(plot_mat_na), 1, 1)
row_annot <- data.frame(
  Category  = cog_domains[cog_letters, "Category"],
  row.names = rownames(plot_mat_na)
)

domain_colors <- list(
  Category = c(
    "Information storage & processing" = "#4E9BB5",
    "Cellular processes & signaling"   = "#E07B39",
    "Metabolism"                        = "#5BA75B",
    "Poorly characterized"              = "#9B6BB5"
  )
)

cogph1 <- pheatmap(plot_mat_na,
                   cluster_cols         = FALSE,
                   cluster_rows         = row_order,
                   color                = colorRampPalette(c("blue", "black", "red"))(100),
                   breaks               = seq(-1, 1, length.out = 101),
                   na_col               = "black",
                   border_color         = "grey",
                   annotation_row       = row_annot,
                   annotation_colors    = domain_colors,
                   annotation_names_row = FALSE,
                   fontsize             = 12,
                   fontsize_row         = 10,
                   angle_col            = 0,
                   cellwidth            = 120,
                   cellheight           = 15)


cogph1

library(ggplot2)
ggsave(plot = cogph1, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/COG_ancombc1_heatmap_04292026.pdf", width = 12.5, height = 6, units = "in", dpi = 300)
ggsave(plot = cogph1, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/COG_ancombc1_heatmap_04292026.tiff", width = 12.5, height = 6, units = "in", dpi = 200, bg="white")



# 05/04/2026 run Maaslin2 instead of ANCOMBC
library(Maaslin2)
library(dplyr)
library(readr)

keggpath_rawr <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_NOunknown_noNAstart_05022026.tsv') # confirmed all(keggpath_rawrt_filt == keggpath_rawrorig_filt)
keggpath_rawr <- as.data.frame(keggpath_rawr)

taxaf <- keggpath_rawr[,c("KEGG_Pathway","taxa")]
taxaf <- taxaf %>%
  column_to_rownames("KEGG_Pathway")

keggpath_rawrt <- keggpath_rawr %>%
  column_to_rownames("KEGG_Pathway")

# drop any outlier samples (based on PCA analysis)
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC","taxa")
keggpath_rawrt_filt <- keggpath_rawrt[,!(names(keggpath_rawrt) %in% drops)]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(keggpath_rawrt_filt),]
bpcoldata3$group <- gsub(" + ", "_", bpcoldata3$group, fixed = TRUE)
bpcoldata3$name <- row.names(bpcoldata3)

# Transpose so rows = samples, columns = features
keggpath_rawrt_filt_t <- as.data.frame(t(keggpath_rawrt_filt))

# --- Comparison 1: CA_S3_AS vs CA_AS ---
meta1 <- bpcoldata3 |> filter(group %in% c("CA_S3_AS", "CA_AS"))
data1 <- keggpath_rawrt_filt_t[rownames(meta1), ]

fit1 <- Maaslin2(
  input_data      = data1,
  input_metadata  = meta1,
  output          = "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/KEGG_pathway_maaslin2_CA_05022026",
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
data2 <- keggpath_rawrt_filt_t[rownames(meta2), ]

fit2 <- Maaslin2(
  input_data      = data2,
  input_metadata  = meta2,
  output          = "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/KEGG_pathway_maaslin2_PHA_05022026",
  fixed_effects   = "group",
  reference       = "group,PHA_AS",
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

# ── Build LFC and q-value matrices ───────────────────────────────────────────
all_cogs <- union(res_CA$feature, res_PHA$feature)

lfc_CA  <- setNames(res_CA$coef,  res_CA$feature)
lfc_PHA <- setNames(res_PHA$coef, res_PHA$feature)
q_CA    <- setNames(res_CA$qval,  res_CA$feature)
q_PHA   <- setNames(res_PHA$qval, res_PHA$feature)

plot_mat <- data.frame(
  CA  = lfc_CA[all_cogs],
  PHA = lfc_PHA[all_cogs],
  row.names = all_cogs
)
plot_mat[is.na(plot_mat)] <- 0

padj_mat <- data.frame(
  CA  = q_CA[all_cogs],
  PHA = q_PHA[all_cogs],
  row.names = all_cogs
)
padj_mat[is.na(padj_mat)] <- 1

combined_mat <- data.frame(
  lfc_CA  = lfc_CA[all_cogs],
  q_CA    = q_CA[all_cogs],
  lfc_PHA = lfc_PHA[all_cogs],
  q_PHA   = q_PHA[all_cogs],
  row.names = all_cogs
)

combined_mat$lfc_CA[is.na(combined_mat$lfc_CA)]   <- 0
combined_mat$lfc_PHA[is.na(combined_mat$lfc_PHA)] <- 0
combined_mat$q_CA[is.na(combined_mat$q_CA)]   <- 1
combined_mat$q_PHA[is.na(combined_mat$q_PHA)] <- 1
combined_mat$KEGG_Pathway <- rownames(combined_mat)
write_tsv(combined_mat, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/KEGG_pathway_maaslin2_CAPHA_coeffpadj_05022026.tsv")


# ── Mask non-significant LFCs ─────────────────────────────────────────────────
plot_mat_masked <- plot_mat
plot_mat_masked[padj_mat > 0.05] <- 0

# ── Replace COG letters with full definitions ─────────────────────────────────
cog_definitions <- c(
  J = "J - Translation, ribosomal structure and biogenesis",
  A = "A - RNA processing and modification",
  K = "K - Transcription",
  L = "L - Replication, recombination and repair",
  B = "B - Chromatin structure and dynamics",
  D = "D - Cell cycle control, cell division, chromosome partitioning",
  Y = "Y - Nuclear structure",
  V = "V - Defense mechanisms",
  T = "T - Signal transduction mechanisms",
  M = "M - Cell wall/membrane/envelope biogenesis",
  N = "N - Cell motility",
  Z = "Z - Cytoskeleton",
  W = "W - Extracellular structures",
  U = "U - Intracellular trafficking, secretion, and vesicular transport",
  O = "O - Posttranslational modification, protein turnover, chaperones",
  X = "X - Mobilome: prophages, transposons",
  C = "C - Energy production and conversion",
  G = "G - Carbohydrate transport and metabolism",
  E = "E - Amino acid transport and metabolism",
  F = "F - Nucleotide transport and metabolism",
  H = "H - Coenzyme transport and metabolism",
  I = "I - Lipid transport and metabolism",
  P = "P - Inorganic ion transport and metabolism",
  Q = "Q - Secondary metabolites biosynthesis, transport and catabolism",
  R = "R - General function prediction only",
  S = "S - Function unknown"
)

rownames(plot_mat_masked) <- cog_definitions[rownames(plot_mat_masked)]
rownames(plot_mat)        <- cog_definitions[rownames(plot_mat)]

# ── Clustering on original (unmasked) LFC ──────────────────────────────────── actually nvm cluster on masked data makes more sense
row_order <- hclust(dist(plot_mat), method = "ward.D2")
row_order <- hclust(dist(plot_mat_masked), method = "ward.D2")


colnames(plot_mat_masked) <- c("CA + S3 + AS vs.\nCA + AS",
                               "PHA + G1 + AS vs.\nPHA + AS")

# ── Row annotations ───────────────────────────────────────────────────────────
cog_domains <- data.frame(
  Category = c(
    "J" = "Information storage & processing",
    "A" = "Information storage & processing",
    "K" = "Information storage & processing",
    "L" = "Information storage & processing",
    "B" = "Information storage & processing",
    "D" = "Cellular processes & signaling",
    "Y" = "Cellular processes & signaling",
    "V" = "Cellular processes & signaling",
    "T" = "Cellular processes & signaling",
    "M" = "Cellular processes & signaling",
    "N" = "Cellular processes & signaling",
    "Z" = "Cellular processes & signaling",
    "W" = "Cellular processes & signaling",
    "U" = "Cellular processes & signaling",
    "O" = "Cellular processes & signaling",
    "X" = "Cellular processes & signaling",
    "C" = "Metabolism",
    "G" = "Metabolism",
    "E" = "Metabolism",
    "F" = "Metabolism",
    "H" = "Metabolism",
    "I" = "Metabolism",
    "P" = "Metabolism",
    "Q" = "Metabolism",
    "R" = "Poorly characterized",
    "S" = "Poorly characterized"
  )
)

cog_letters <- substr(rownames(plot_mat_masked), 1, 1)
row_annot <- data.frame(
  Category = cog_domains[cog_letters, "Category"],
  row.names = rownames(plot_mat_masked)
)

domain_colors <- list(
  Category = c(
    "Information storage & processing" = "#4E9BB5",
    "Cellular processes & signaling"   = "#E07B39",
    "Metabolism"                        = "#5BA75B",
    "Poorly characterized"              = "#9B6BB5"
  )
)

# ── Heatmap ───────────────────────────────────────────────────────────────────
cogph <- pheatmap(plot_mat_masked,
                  cluster_cols             = FALSE,
                  # cluster_rows             = row_order,
                  cluster_rows             = TRUE,
                  color                    = colorRampPalette(c("blue", "black", "red"))(100),
                  breaks                   = seq(-1, 1, length.out = 101),
                  border_color             = "grey",
                  annotation_row           = row_annot,
                  annotation_colors        = domain_colors,
                  annotation_names_row     = FALSE,
                  fontsize                 = 12,
                  fontsize_row             = 10,
                  angle_col                = 0,
                  cellwidth                = 120,
                  cellheight               = 15,
                  legend_labels = "Coefficient")

cogph

library(ggplot2)
ggsave(plot = cogph, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/COG_maaslin2_heatmap_05042026.pdf", width = 12.5, height = 6, units = "in", dpi = 300)
ggsave(plot = cogph, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/COG_maaslin22_heatmap_05042026.tiff", width = 12.5, height = 6, units = "in", dpi = 200, bg="white")





















#######################################################################################################################################
#### CA + S3 + AS vs. CA + AS comparison

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "CA + AS" | bpcoldata3$group == "CA + S3 + AS",]
bpcoldata3_sub

keggpath_rawrt_filt_sub <- keggpath_rawrt_filt[,colnames(keggpath_rawrt_filt) %in% rownames(bpcoldata3_sub)]

all(colnames(keggpath_rawrt_filt_sub) == rownames(bpcoldata3_sub))

##############################################################################################################################
# Now do for PHA + G1 + AS vs. PHA + AS comparison

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "PHA + G1 + AS" | bpcoldata3$group == "PHA + AS",]

keggpath_rawrt_filt_sub <- keggpath_rawrt_filt[,colnames(keggpath_rawrt_filt) %in% rownames(bpcoldata3_sub)]

all(colnames(keggpath_rawrt_filt_sub) == rownames(bpcoldata3_sub))


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

otumat = as.matrix(keggpath_rawrt_filt_sub)

taxmat = matrix(sample(letters, 25, replace = TRUE), nrow = nrow(otumat), ncol = 7)
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
# res_global = out$res_global

## ANCOMBC primary result
# LFC
tab_lfc = res$lfc
col_name = c("KEGG_Pathway", "Intercept", "CA + S3 + AS")
col_name = c("KEGG_Pathway", "Intercept", "PHA + G1 + AS")
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
  datatable(caption = "Differentially Abundant KEGG Pathways from the Primary Result")

tab_lfc = res$lfc
col_name = c("KEGG_Pathway", "Intercept", "CAS3AS_v_CAAS")
col_name = c("KEGG_Pathway", "Intercept", "PHAG1AS_v_PHAAS")
colnames(tab_lfc) = col_name
tab_lfc <- tab_lfc[,-2]

taxaf$KEGG_Pathway <- rownames(taxaf)

# CA
tab_lfc_CAS3AS <- merge(tab_lfc, taxaf, by = "KEGG_Pathway")
write_tsv(tab_lfc_CAS3AS, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/tab_lfc_CAS3AS_KEGGpathway_ancombc_12192025.tsv")
tab_lfc_CAS3AS <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/tab_lfc_CAS3AS_KEGGpathway_ancombc_12192025.tsv")

tab_lfc_CAS3AS2 <- tab_lfc_CAS3AS %>%
  rename("LFC" = "CAS3AS_v_CAAS")

tab_q2 <- tab_q %>%
  rename( "padj" = "CA + S3 + AS")
tab_q2 <- tab_q2[,-2]

df_volc_CAS3AS2 <- merge(tab_lfc_CAS3AS2, tab_q2, by = "KEGG_Pathway")
write_tsv(df_volc_CAS3AS2, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/df_volc_CAS3AS2_KEGGpathway_ancombc_12192025.tsv") # use this for volcano plots
df_volc_CAS3AS2 <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/df_volc_CAS3AS2_KEGGpathway_ancombc_12192025.tsv")
keggpathanno <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/KEGG_pathway_anno_12192025.txt")
df_volc_CAS3AS2_path <- merge(df_volc_CAS3AS2, keggpathanno, by = "KEGG_Pathway", all.x = TRUE)
df_volc_CAS3AS2_path <- df_volc_CAS3AS2_path[order(df_volc_CAS3AS2_path$padj),]
write_tsv(df_volc_CAS3AS2_path, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/df_volc_CAS3AS2_KEGGpathway_anno_ancombc_12192025.tsv") 
df_volc_CAS3AS2_path <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/df_volc_CAS3AS2_KEGGpathway_anno_ancombc_12192025.tsv")
df_volc_CAS3AS2_path_up <- df_volc_CAS3AS2_path[df_volc_CAS3AS2_path$LFC>0,]

# PHA
tab_lfc_PHAG1AS <- merge(tab_lfc, taxaf, by = "KEGG_Pathway")
write_tsv(tab_lfc_PHAG1AS, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/tab_lfc_PHAG1AS_KEGGpathway_ancombc_12192025.tsv")
tab_lfc_PHAG1AS <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/tab_lfc_PHAG1AS_KEGGpathway_ancombc_12192025.tsv")

tab_lfc_PHAG1AS2 <- tab_lfc_PHAG1AS %>%
  rename("LFC" = "PHAG1AS_v_PHAAS")

tab_q2 <- tab_q %>%
  rename("padj" = "PHA + G1 + AS")
tab_q2 <- tab_q2[,-2]

df_volc_PHAG1AS2 <- merge(tab_lfc_PHAG1AS2, tab_q2, by = "KEGG_Pathway")
write_tsv(df_volc_PHAG1AS2, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/df_volc_PHAG1AS2_KEGGpathway_ancombc_12192025.tsv") # use this for volcano plots
df_volc_PHAG1AS2 <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/df_volc_PHAG1AS2_KEGGpathway_ancombc_12192025.tsv")
keggpathanno <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/KEGG_pathway_anno_12192025.txt")
df_volc_PHAG1AS2_path <- merge(df_volc_PHAG1AS2, keggpathanno, by = "KEGG_Pathway", all.x = TRUE)
df_volc_PHAG1AS2_path <- df_volc_PHAG1AS2_path[order(df_volc_PHAG1AS2_path$padj),]
write_tsv(df_volc_PHAG1AS2_path, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/df_volc_PHAG1AS2_KEGGpathway_anno_ancombc_12192025.tsv") # but update this to do analysis not excluding those two samples !! and then redo this
df_volc_PHAG1AS2_path <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/df_volc_PHAG1AS2_KEGGpathway_anno_ancombc_12192025.tsv")
df_volc_PHAG1AS2_path_up <- df_volc_PHAG1AS2_path[df_volc_PHAG1AS2_path$LFC > 0,]

########################################################################################################################################
# making volcano plots instead of heatmmaps to show differentially expressed CAZymes
# taken from https://biostatsquid.com/volcano-plots-r-tutorial/

# Loading relevant libraries 
library(tidyverse) # includes ggplot2, for data visualisation. dplyr, for data manipulation.
library(RColorBrewer) # for a colourful plot
library(ggrepel) # for nice annotations

# Import Ancombc results

# for CAS3AS vs CAAS comparison

df <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/df_volc_CAS3AS2_KEGGpathway_anno_ancombc_12192025.tsv")

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

# if log2Foldchange > 1 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$LFC > 1 & df$padj < 0.05] <- "UP"

# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
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
df$delabel <- ifelse(df$KEGG_Pathway %in% as.matrix(head(dlup[order(dlup$padj), "KEGG_Pathway"], 10)), df$KEGG_Pathway, NA)

dldown <- df[df$diffexpressed == "DOWN",]
df$delabel <- ifelse(df$KEGG_Pathway %in% as.matrix(head(dldown[order(dldown$padj), "KEGG_Pathway"], 10)), df$KEGG_Pathway, df$delabel)

ggplot(data = df, aes(x = LFC, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 60), xlim = c(-4, 4)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-4, 4, 2)) + # to customise the breaks in the x axis
  #geom_point() + geom_text(show.legend = FALSE, nudge_x = .2, nudge_y = 1)
  geom_label_repel(size = 4, show.legend = FALSE, min.segment.length = unit(0, 'lines')) +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), legend.title = element_blank(), legend.position = "bottom")

# Probably adjust the LFC and padj cutoff here - the KEGG pathways with higher LFCs but padj values close to
# 0.05 are likely not as important/relevant as the KEGG pathways with lower LFCs but padj valuesmuch less than 
# 0.05
# 12/19/2025  didn't save these graphs cus I want to adjust the LFC and padj cutoffs before I do so


## for G1PHAAS vs PHAAS comparison

df <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/df_volc_PHAG1AS2_KEGGpathway_anno_ancombc_12192025.tsv")

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

# if log2Foldchange > 1 and pvalue < 0.05, set as "UP"
df$diffexpressed[df$LFC > 1 & df$padj < 0.05] <- "UP"

# if log2Foldchange < -1 and pvalue < 0.05, set as "DOWN"
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
df$delabel <- ifelse(df$KEGG_Pathway %in% as.matrix(head(dlup[order(dlup$padj), "KEGG_Pathway"], 10)), df$KEGG_Pathway, NA)

dldown <- df[df$diffexpressed == "DOWN",]
df$delabel <- ifelse(df$KEGG_Pathway %in% as.matrix(head(dldown[order(dldown$padj), "KEGG_Pathway"], 10)), df$KEGG_Pathway, df$delabel)

ggplot(data = df, aes(x = LFC, y = -log10(padj), col = diffexpressed, label = delabel)) +
  geom_vline(xintercept = c(-1, 1), col = "gray", linetype = 'dashed') +
  geom_hline(yintercept = -log10(0.05), col = "gray", linetype = 'dashed') +
  geom_point(size = 2) +
  scale_color_manual(values = c("blue", "grey", "red"), # to set the colours of our variable
                     labels = c("Downregulated", "Not significant", "Upregulated")) +  # to set the labels in case we want to overwrite the categories from the dataframe (UP, DOWN, NO)
  coord_cartesian(ylim = c(0, 10), xlim = c(-2, 2)) + # since some genes can have minuslog10padj of inf, we set these limits
  labs(color = 'Differential Expression', #legend_title, 
       x = expression("log"[2]*"FC"), y = expression("-log"[10]*"p-value")) + 
  scale_x_continuous(breaks = seq(-2, 2, 2)) + # to customise the breaks in the x axis
  #geom_point() + geom_text(show.legend = FALSE, nudge_x = .2, nudge_y = 1)
  geom_label_repel(size = 4, show.legend = FALSE, min.segment.length = unit(0, 'lines')) +
  theme(axis.title.x = element_text(size = 18), axis.title.y = element_text(size = 18), legend.title = element_blank(), legend.position = "bottom")

# Probably adjust the LFC and padj cutoff here - the KEGG pathways with higher LFCs but padj values close to
# 0.05 are likely not as important/relevant as the KEGG pathways with lower LFCs but padj values much less than 
# 0.05
# 12/19/2025  didn't save these graphs cus I want to adjust the LFC and padj cutoffs before I do so
# lol the pathways for this comparison are super weird...measles?



###########################################################################################################
###########################################################################################################
###########################################################################################################
# If using ANCOM-BC is not appropriate, it's probably best to normalize the raw counts, aggregate them by 
# group, and then perform the appropriate statistical test (Wilcoxin, KW, whatever JGI recommends?)
# Let's do TPM normalization for now (see what JGI says for normalization method?)

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
all(bpcts$GeneID %in% contig_gene2$Geneid)
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
gl1 <- bpcts[,c("Geneid", "Length")]

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

# Generate a tpm matrix
all(rownames(bpcts) == gl1$Geneid)
rownames(bpcts) <- NULL
gl <- as.numeric(gl1$Length)

bptpm <- as.data.frame(tpm(bpctsm, gl))
bpcts$GeneID <- rownames(bpcts)

# Input gene annotation info
gene_anno <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/gene_annotation_metassembly.tsv")

bpcts_geneanno <- merge(bpcts, gene_anno, by = "GeneID", all = TRUE)
bpcts_geneanno2 <- merge(bpcts_geneanno, tax_anno, by = "GeneID", all = TRUE)
bpcts_geneanno2$KEGG_Pathway <- str_remove(bpcts_geneanno2$KEGG_Pathway, ",map.*")
bpcts_geneanno2$KEGG_Pathway <- gsub("ko", "map", bpcts_geneanno2$KEGG_Pathway)

library(dplyr)
library(stringr)
library(tidyr)

# let's now do the KEGG ko annotations
aa <- bpcts_geneanno2 %>%
  dplyr::select(GeneID, everything())%>%
  mutate(KEGG_Pathway=str_split(KEGG_Pathway, ",")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:136,1)]

unique(bpcts_geneanno2$KEGG_Pathway)
length(unique(bpcts_geneanno2$KEGG_Pathway))
ab<-as.vector(as.matrix(aa[,grepl("KEGG_Pathway",names(aa))]))
unique(ab)
length(unique(ab))

# make subset of aa data
bb <- aa[is.na(aa$KEGG_Pathway3) == FALSE,]
aasub <- as.data.frame(bb[1:5,])
aasub$KEGG_Pathway1[is.na(aasub$KEGG_Pathway1) == TRUE] <- "-"
aasubccs <- aasub[,c(32:41)]

uab <-unique(ab)
uab
length(uab)
rn <- uab[!is.na(uab)]
rn
length(rn)

keggpath_raw <- data.frame(matrix(0, nrow = length(rn), ncol = 20))
colnames(keggpath_raw) <- colnames(aasub[,1:20])
rownames(keggpath_raw) <- rn
keggpath_raw$taxa <- NA

# Trial run on subset of the data
for (j in 1:(ncol(keggpath_raw)-1)) {
  for (k in 1:length(names(aasub)[grepl("KEGG_Pathway" , names(aasub))])) { # Number of separate KEGG pathways broken up
    for (i in 1:nrow(aasub)) {
      if(is.na(aasub[i,paste("KEGG_Pathway",k, sep = "")]) == FALSE) {
        keggpath_raw[rn[rn %in% aasub[i,paste("KEGG_Pathway",k, sep = "")]],j] = keggpath_raw[rn[rn %in% aasub[i,paste("KEGG_Pathway",k, sep = "")]],j] + as.numeric(aasub[i,j])
      }
      else {
        keggpath_raw = keggpath_raw
      }
    }
  }
}

for (k in 1:length(aasub[,grepl("KEGG_Pathway" , names(aasub))])) { # Number of separate KEGG pathways broken up
  for (i in 1:nrow(aasub)) {
    if(aasub[i,paste("KEGG_Pathway",k, sep = "")] %in% rn == TRUE ) {
      keggpath_raw[rn[rn %in% aasub[i,paste("KEGG_Pathway",k, sep = "")]],"taxa"] <- paste(keggpath_raw[rn[rn %in% aasub[i,paste("KEGG_Pathway",k, sep = "")]],"taxa"],"; ", aasub$taxa[i], sep = "")
    }
    else {
      keggpath_raw = keggpath_raw
    }
  }
}


## Now do it on the entire data
aa$KEGG_Pathway1[is.na(aa$KEGG_Pathway1) == TRUE] <- "-"

uab <-unique(ab)
uab
length(uab)
rn <- uab[!is.na(uab)]
rn
length(rn)

keggpath_raw <- data.frame(matrix(0, nrow = length(rn), ncol = 20))
colnames(keggpath_raw) <- colnames(aa[,1:20])
rownames(keggpath_raw) <- rn
keggpath_raw$taxa <- NA

for (j in 1:(ncol(keggpath_raw)-1)) {
  for (k in 1:length(names(aa)[grepl("KEGG_Pathway" , names(aa))])) { # Number of separate KEGG pathways broken up
    for (i in 1:nrow(aa)) {
      if(is.na(aa[i,paste("KEGG_Pathway",k, sep = "")]) == FALSE) {
        keggpath_raw[rn[rn %in% aa[i,paste("KEGG_Pathway",k, sep = "")]],j] = keggpath_raw[rn[rn %in% aa[i,paste("KEGG_Pathway",k, sep = "")]],j] + as.numeric(aa[i,j])
      }
      else {
        keggpath_raw = keggpath_raw
      }
    }
  }
}

keggpath_raw2 <- keggpath_raw %>%
  rownames_to_column("KEGG_Pathway")
write_tsv(keggpath_raw2, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_notaxaanno_12182025.tsv")

for (k in 1:length(aa[,grepl("KEGG_Pathway" , names(aa))])) { # Number of separate KEGG pathways broken up
  for (i in 1:nrow(aa)) {
    if(aa[i,paste("KEGG_Pathway",k, sep = "")] %in% rn == TRUE ) {
      keggpath_raw[rn[rn %in% aa[i,paste("KEGG_Pathway",k, sep = "")]],"taxa"] <- paste(keggpath_raw[rn[rn %in% aa[i,paste("KEGG_Pathway",k, sep = "")]],"taxa"],"; ", aa$taxa[i], sep = "")
    }
    else {
      keggpath_raw = keggpath_raw
    }
  }
}

keggpath_raw2 <- keggpath_raw %>%
  rownames_to_column("KEGG_Pathway")
write_tsv(keggpath_raw2, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_12192025.tsv")
keggpath_rawr <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_12192025.tsv")
keggpath_rawr <- keggpath_rawr %>%
  column_to_rownames("KEGG_Pathway")

# Remove the results for "-" since I think this is crashing my computer, plus I don't care about the taxa associated with non-CAZymes
drop <- c("-")
keggpath_rawrr <- keggpath_rawr[!(keggpath_rawr$KEGG_Pathway %in% drop),]
write_tsv(keggpath_rawrr, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_NOunknown_12192025.tsv')

keggpath_rawrrr <- keggpath_rawrr
keggpath_rawrrr$taxa <- substr(keggpath_rawrrr$taxa, 4, nchar(keggpath_rawrrr$taxa))
write_tsv(keggpath_rawrrr, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_NOunknown_noNAstart_12192025.tsv')
keggpath_rawr <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/keggpathway_raw_taxaanno_NOunknown_noNAstart_12192025.tsv')
keggpath_rawr <- as.data.frame(keggpath_rawr)

taxaf <- keggpath_rawr[,c("KEGG_Pathway","taxa")]
taxaf <- taxaf %>%
  column_to_rownames("KEGG_Pathway")

keggpath_rawrt <- keggpath_rawr %>%
  column_to_rownames("KEGG_Pathway")

# drop any outlier samples (based on PCA analysis)
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC","taxa")
keggpath_rawrt_filt <- keggpath_rawrt[,!(names(keggpath_rawrt) %in% drops)]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(keggpath_rawrt_filt),]

#######################################################################################################################################
#### CA + S3 + AS vs. CA + AS comparison

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "CA + AS" | bpcoldata3$group == "CA + S3 + AS",]
bpcoldata3_sub

keggpath_rawrt_filt_sub <- keggpath_rawrt_filt[,colnames(keggpath_rawrt_filt) %in% rownames(bpcoldata3_sub)]

all(colnames(keggpath_rawrt_filt_sub) == rownames(bpcoldata3_sub))

##############################################################################################################################
# Now do for PHA + G1 + AS vs. PHA + AS comparison

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "PHA + G1 + AS" | bpcoldata3$group == "PHA + AS",]

keggpath_rawrt_filt_sub <- keggpath_rawrt_filt[,colnames(keggpath_rawrt_filt) %in% rownames(bpcoldata3_sub)]

all(colnames(keggpath_rawrt_filt_sub) == rownames(bpcoldata3_sub))

#######################################################################################################################################


