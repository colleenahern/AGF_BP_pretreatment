# Colleen Ahern
# 03/02/2026

# Continuation of BP_MT_transcript_abund_grouped_02062025.R script
# Normalizing my raw gene counts using JGI's method instead of tpm and add together COG groups specifically
# Divide gene counts by numer of gene counts assigned as COGs and then add those together
# and make a heatmap for PHA + G1 + AS vs. PHA + AS and CA + S3 + AS vs. CA + AS
# before I was separating by timepoint but I don't want to do that right now
# use these heatmaps for Magda's slides

library('DRnaSeq')
library(readr)

## Count matrix input
bpcoldata <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/bp_metadata_01132025.txt")

bpcts <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/read_counts_gene_metassembly.tsv")

bpcts <- bpcts %>%
  column_to_rownames("GeneID")
bpcts <- bpcts[,-(1:5)]
bpcoldata
bpcoldata <- bpcoldata %>%
  column_to_rownames("name")
head(bpcts)
bpcoldata

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

# Input gene annotation info
gene_anno <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/gene_annotation_metassembly.tsv")

bpcts$GeneID <- rownames(bpcts)
bpcts_geneanno <- merge(bpcts, gene_anno, by = "GeneID", all = TRUE)
bpcts_geneanno <- bpcts_geneanno %>%
  column_to_rownames("GeneID")

# Make vector that contains the sum of all COG-annotated gene counts for each sample
cogvec <- matrix(0, nrow = 1, ncol = length(names(bpcts_geneanno)[grepl("MT_PLANC", names(bpcts_geneanno))]))
cogvec <- as.data.frame(cogvec)
cogvec
names(cogvec) <- as.list(names(bpcts_geneanno)[grepl("MT_PLANC", names(bpcts_geneanno))])
all(names(cogvec) == names(bpcts_geneanno)[grepl("MT_PLANC", names(bpcts_geneanno))])
cogvec

# try on a subset of data
ctssub <- bpcts_geneanno[1:8,]
cogvectry <- cogvec
cogvectry

for (j in 1:length(names(ctssub)[grepl("MT_PLANC", names(ctssub))])) {
    for (i in 1:nrow(ctssub)) {
      if(is.na(ctssub[i,"COG_category"]) == FALSE) {
        cogvectry[j] = as.numeric(cogvectry[j]) + as.numeric(ctssub[i,j])
      }
      else {
        cogvectry = cogvectry
      }
    }
}

cogvec
# # Do on entire data - didn't end up using because it is computationally intensive
# for (j in 1:length(names(bpcts_geneanno)[grepl("MT_PLANC", names(bpcts_geneanno))])) {
#   for (i in 1:nrow(bpcts_geneanno)) {
#     if(is.na(bpcts_geneanno[i,"COG_category"]) == FALSE) {
#       cogvec[j] = as.numeric(cogvec[j]) + as.numeric(bpcts_geneanno[i,j])
#     }
#     else {
#       cogvec = cogvec
#     }
#   }
# }

# Actually above loop is way more computationally intensive than it needs to be. Start by removing all genes with NA values in COG
ctssub_cog <- ctssub[is.na(ctssub$COG_category) == FALSE,]
bpcts_geneanno_cog <- bpcts_geneanno[is.na(bpcts_geneanno$COG_category) == FALSE,]
cogvecb <- matrix(0, nrow = 1, ncol = length(names(bpcts_geneanno_cog)[grepl("MT_PLANC", names(bpcts_geneanno_cog))]))
cogvecb <- as.data.frame(cogvecb)
cogvecb
names(cogvecb) <- as.list(names(bpcts_geneanno_cog)[grepl("MT_PLANC", names(bpcts_geneanno_cog))])
all(names(cogvecb) == names(bpcts_geneanno_cog)[grepl("MT_PLANC", names(bpcts_geneanno_cog))])
cogvecb

# try on subset first
cogvecbtry <- cogvecb
for (j in 1:length(names(ctssub_cog)[grepl("MT_PLANC", names(ctssub_cog))])) {
  for (i in 1:nrow(ctssub_cog)) {
      cogvecbtry[j] = as.numeric(cogvecbtry[j]) + as.numeric(ctssub_cog[i,j])
  }
}

# Now do on entire data
for (j in 1:length(names(bpcts_geneanno_cog)[grepl("MT_PLANC", names(bpcts_geneanno_cog))])) {
  for (i in 1:nrow(bpcts_geneanno_cog)) {
    cogvecb[j] = as.numeric(cogvecb[j]) + as.numeric(bpcts_geneanno_cog[i,j])
  }
}

# Also try easier way - just sum all in each column and should get the same result - yup this does the same thing but takes way less time
cogvecb_2 <- cogvecb
all(names(cogvecb_2) == names(bpcts_geneanno_cog)[grepl("MT_PLANC", names(bpcts_geneanno_cog))])
for (j in 1:length(names(bpcts_geneanno_cog)[grepl("MT_PLANC", names(bpcts_geneanno_cog))])) {
  cogvecb_2[j] = sum(as.numeric(bpcts_geneanno_cog[,j]))
}

# Divide each count by the total number of COG counts in that sample
bpcts_geneanno_cognorm <- bpcts_geneanno_cog
all(names(bpcts_geneanno_cognorm)[grepl("MT_PLANC", names(bpcts_geneanno_cognorm))] == names(cogvecb))
for (j in 1:length(names(bpcts_geneanno_cognorm)[grepl("MT_PLANC", names(bpcts_geneanno_cognorm))])) {
  bpcts_geneanno_cognorm[,j] = as.numeric(bpcts_geneanno_cognorm[,j])/as.numeric(cogvecb[j])
}

# choose a pesudocount value for log2 transformation
sub <- data.matrix(bpcts_geneanno_cognorm[, grepl("MT_PLANC", colnames(bpcts_geneanno_cognorm)), drop = FALSE])
pseudct <- 0.5*min(sub[is.finite(sub) & sub > 0])

# log transform with pseudo count
bpcts_geneanno_cognormlogt <- bpcts_geneanno_cognorm
for (j in 1:length(names(bpcts_geneanno_cognormlogt)[grepl("MT_PLANC", names(bpcts_geneanno_cognormlogt))])) {
  bpcts_geneanno_cognormlogt[,j] = log2(as.numeric(bpcts_geneanno_cognorm[,j]) + pseudct)
}





library(dplyr)
library(stringr)
library(tidyr)
bpcts_geneanno_cognorm$GeneID <- row.names(bpcts_geneanno_cognorm)

# let's now do the COG annotations
aa <- bpcts_geneanno_cognorm %>%
  select(GeneID, everything())%>%
  mutate(COG_category=str_split(COG_category, "")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:48,1)]
#rownames(aa) <- 1:nrow(aa)

unique(bpcts_geneanno_cognorm$COG_category)
ab<-as.vector(as.matrix(aa[,grepl("COG_category",names(aa))]))
unique(ab)
length(unique(ab))

# make subset of aa data
aasub <- as.data.frame(aa[c(1:6,15),])
aasubccs <- aasub[,c(26:31)]

uab <-unique(ab)
uab
length(uab)
rn <- uab[!is.na(uab)]
rn
length(rn)

cogtable_gen <- data.frame(matrix(0, nrow = length(rn), ncol = 20))
colnames(cogtable_gen) <- colnames(aasub[,1:20])
rownames(cogtable_gen) <- rn
cogtable_gen
all(names(cogtable_gen) == names(bpcts_geneanno_cognorm)[1:20])

# Trial run on subset of the data
for (j in 1:ncol(cogtable_gen)) {
  for (k in 1:length(names(aasub)[grepl("COG_category" , names(aasub))])) { # Number of separate COG categories broken up
    for (i in 1:nrow(aasub)) {
      if(is.na(aasub[i,paste("COG_category",k, sep = "")]) == FALSE) {
        cogtable_gen[rn[rn %in% aasub[i,paste("COG_category",k, sep = "")]],j] = cogtable_gen[rn[rn %in% aasub[i,paste("COG_category",k, sep = "")]],j] + aasub[i,j]
      }
      else {
        cogtable_gen = cogtable_gen
      }
    }
  }
}


## Now do it on the entire data
uab <-unique(ab)
uab
length(uab)
rn <- uab[!is.na(uab)]
rn
length(rn)

cogtable_gen <- data.frame(matrix(0, nrow = length(rn), ncol = 20))
colnames(cogtable_gen) <- colnames(aa[,1:20])
rownames(cogtable_gen) <- rn
cogtable_gen
all(names(cogtable_gen) == names(bpcts_geneanno_cognorm)[1:20])

for (j in 1:ncol(cogtable_gen)) {
  for (k in 1:length(names(aa)[grepl("COG_category" , names(aa))])) { # Number of separate COG categories broken up
    for (i in 1:nrow(aa)) {
      if(is.na(aa[i,paste("COG_category",k, sep = "")]) == FALSE) {
        cogtable_gen[rn[rn %in% aa[i,paste("COG_category",k, sep = "")]],j] = cogtable_gen[rn[rn %in% aa[i,paste("COG_category",k, sep = "")]],j] + aa[i,j]
      }
      else {
        cogtable_gen = cogtable_gen
      }
    }
  }
}

cogtable_gen$COGid <- rownames(cogtable_gen)
write_tsv(cogtable_gen, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/cogtable_jginorm_03022026.tsv")
cogtable_genr <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/cogtable_jginorm_03022026.tsv")

cogtable_genr <- as.data.frame(cogtable_genr)
cogtable_genr <- cogtable_genr %>%
  column_to_rownames("COGid")

# Take average of replicates
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
cogtable_gen_filt <- cogtable_genr[,!(names(cogtable_genr) %in% drops)]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(cogtable_gen_filt),]

# cogtable_gen_avg <- data.frame(matrix(0, nrow = 25, ncol = 8))
# colnames(cogtable_gen_avg) <- unique(bpcoldata3$condition_cba)
# cogtable_gen_avg <- data.frame(matrix(0, nrow = 25, ncol = 4))
# colnames(cogtable_gen_avg) <- unique(bpcoldata3$group)
cogtable_gen_avg <- data.frame(matrix(0, nrow = 25, ncol = 12))
colnames(cogtable_gen_avg) <- c(unique(bpcoldata3$condition_cba),unique(bpcoldata3$group))
rownames(cogtable_gen_avg) <- rownames(cogtable_gen_filt)

cogtable_gen_avg[,'CA + AS t1'] <- rowMeans(as.matrix(cogtable_gen_filt[, names(cogtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + AS t1',])]))
cogtable_gen_avg[,'CA + S3 + AS t1'] <- rowMeans(as.matrix(cogtable_gen_filt[, names(cogtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + S3 + AS t1',])]))
cogtable_gen_avg[,'CA + AS t2'] <- rowMeans(as.matrix(cogtable_gen_filt[, names(cogtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + AS t2',])]))
cogtable_gen_avg[,'CA + S3 + AS t2'] <- rowMeans(as.matrix(cogtable_gen_filt[, names(cogtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + S3 + AS t2',])]))

cogtable_gen_avg[,'PHA + AS t1'] <- rowMeans(as.matrix(cogtable_gen_filt[, names(cogtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + AS t1',])]))
cogtable_gen_avg[,'PHA + G1 + AS t1'] <- rowMeans(as.matrix(cogtable_gen_filt[, names(cogtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + G1 + AS t1',])]))
cogtable_gen_avg[,'PHA + AS t2'] <- rowMeans(as.matrix(cogtable_gen_filt[, names(cogtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + AS t2',])]))
cogtable_gen_avg[,'PHA + G1 + AS t2'] <- rowMeans(as.matrix(cogtable_gen_filt[, names(cogtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + G1 + AS t2',])]))

cogtable_gen_avg[,'CA + AS'] <- rowMeans(as.matrix(cogtable_gen_filt[, names(cogtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$group == 'CA + AS',])]))
cogtable_gen_avg[,'CA + S3 + AS'] <- rowMeans(as.matrix(cogtable_gen_filt[, names(cogtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$group == 'CA + S3 + AS',])]))
cogtable_gen_avg[,'PHA + AS'] <- rowMeans(as.matrix(cogtable_gen_filt[, names(cogtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$group == 'PHA + AS',])]))
cogtable_gen_avg[,'PHA + G1 + AS'] <- rowMeans(as.matrix(cogtable_gen_filt[, names(cogtable_gen_filt) %in% rownames(bpcoldata3[bpcoldata3$group == 'PHA + G1 + AS',])]))

cogtable_gen_avg$COGid <- rownames(cogtable_gen_avg)
write_tsv(cogtable_gen_avg, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/COGtable_jginorm_avg_03022026.tsv")
########################################################################################################################
## Make heatmap
## BP: Import metadata file

## BP: JGI-normalized transcripts -- COGs added together
vary_43 <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/COGtable_jginorm_avg_03022026.tsv")
vary_43 <- vary_43 %>%
  column_to_rownames("COGid")

newnames <- lapply(
  colnames(vary_43),
  function(x) bquote(.(x)))
newnames
head(vary_43)
nrow(vary_43)

data_matrix_44 <- vary_43
drop <- c("-", "S")
data_matrix_44 <- data_matrix_44[!(row.names(data_matrix_44) %in% drop),]
data_matrix_44

data_matrix_44_log <- log(data_matrix_44 + 1)

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
  res = textGrob(coln, x = x, y = unit(1, "npc") - unit(3,"bigpts"), vjust = 2, hjust = 1, rot = 45, gp = gpar(...))
  return(res)}
assignInNamespace(x="draw_colnames", value="draw_colnames_45",
                  ns=asNamespace("pheatmap"))

my_palette <- colorRampPalette(c("black", "red"))(n = 299)

bp_cogtpm_heatmap1<-pheatmap(data_matrix_44[,1:8], col=my_palette, cellwidth=80,cellheight=15, cluster_rows=F, cluster_cols=F, 
                            fontsize=12,show_rownames=T, border_color = "black", labels_row = rownames(data_matrix_44), labels_col = as.character(newnames[1:8]), breaks=seq(0, 0.16, length.out=300))

bp_cogtpm_heatmap2<-pheatmap(data_matrix_44[,9:12], col=my_palette, cellwidth=80,cellheight=15, cluster_rows=F, cluster_cols=F, 
                             fontsize=12,show_rownames=T, labels_row = rownames(data_matrix_44), labels_col = as.character(newnames[9:12]), border_color="black",breaks=seq(0, 0.16, length.out=300), angle_col = 45)



ggsave(plot = bp_cogtpm_heatmap1, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/COG_t1t2_heatmap_03022026.pdf", width = 10, height = 7, units = "in", dpi = 300)
ggsave(plot = bp_cogtpm_heatmap1, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/COG_t1t2_heatmap_03022026.tiff", width = 10, height = 7, units = "in", dpi = 200, bg="white")

ggsave(plot = bp_cogtpm_heatmap2, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/COG_tcomb_heatmap_03022026.pdf", width = 6, height = 6.5, units = "in", dpi = 300)
ggsave(plot = bp_cogtpm_heatmap2, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/COG_tcomb_heatmap_03022026.tiff", width = 6, height = 6.5, units = "in", dpi = 200, bg="white")



####################################################################################################################################
# do raw counts - claude says aggregate them normalize/transform

library('DRnaSeq')
library(readr)
library(tidyverse)

## Count matrix input
bpcoldata <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/bp_metadata_01132025.txt")
bpcts <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/read_counts_gene_metassembly.tsv")


bpcts <- bpcts %>%
  column_to_rownames("GeneID")
bpcts <- bpcts[,-(1:5)]
bpcoldata
bpcoldata <- bpcoldata %>%
  column_to_rownames("name")
head(bpcts)
bpcoldata

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
  mutate(COG_category=str_split(COG_category, "")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:48,1)]
#rownames(aa) <- 1:nrow(aa)

unique(bpcts_geneanno$COG_category)
ab<-as.vector(as.matrix(aa[,names(aa)[grepl("COG_category", names(aa))]]))
unique(ab)
length(unique(ab))

# make subset of aa data
aasub <- aa[c(1,2,3,21,18393,198),]
aasub$COG_category1[is.na(aasub$COG_category1)] <- "-"
aasubccs <- aasub[,c(26:31)]

cogtable_gen <- data.frame(matrix(0, nrow = 25, ncol = 20))
colnames(cogtable_gen) <- colnames(aa[,1:20])
uab <-unique(ab)
uab <- uab[!is.na(uab)]
rn <- uab
rownames(cogtable_gen) <- rn

# Trial run on subset of the data
for (j in 1:ncol(cogtable_gen)) {
  for (k in 1:6) {
    for (i in 1:nrow(aasub)) {
      if(is.na(aasub[i,paste("COG_category",k, sep = "")]) == FALSE) {
        cogtable_gen[rn[rn %in% aasub[i,paste("COG_category",k, sep = "")]],j] = cogtable_gen[rn[rn %in% aasub[i,paste("COG_category",k, sep = "")]],j] + aasub[i,j]
      }
      else {
        cogtable_gen = cogtable_gen
      }
    }
  }
}


## Now do it on the entire data
aa$COG_category1[is.na(aa$COG_category1)] <- "-"

cogtable_gen <- data.frame(matrix(0, nrow = 25, ncol = 20))
colnames(cogtable_gen) <- colnames(aa[,1:20])
uab <-unique(ab)
uab <- uab[!is.na(uab)]
rn <- uab
rownames(cogtable_gen) <- rn

for (j in 1:ncol(cogtable_gen)) {
  for (k in 1:6) {
    for (i in 1:nrow(aa)) {
      if(is.na(aa[i,paste("COG_category",k, sep = "")]) == FALSE) {
        cogtable_gen[rn[rn %in% aa[i,paste("COG_category",k, sep = "")]],j] = cogtable_gen[rn[rn %in% aa[i,paste("COG_category",k, sep = "")]],j] + aa[i,j]
      }
      else {
        cogtable_gen = cogtable_gen
      }
    }
  }
}

cogtable_gen$COGid <- rownames(cogtable_gen)
# write_tsv(cogtable_gen, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/cogtable_gen_raw_03062025.tsv')
write_tsv(cogtable_gen, '/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/cogtable_gen_raw_04292026.tsv')
cogtable_gen <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/cogtable_gen_raw_04292026.tsv')

cogtable_gen <- cogtable_gen %>%
  column_to_rownames("COGid")

colSums(cogtable_gen)

# filter out samples/COGs I dont want
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC", "-", "Y")
cogtable_gen_filt <- cogtable_gen[!rownames(cogtable_gen) %in% drops,]
cogtable_gen_filt <- cogtable_gen_filt[,!names(cogtable_gen_filt) %in% drops]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(cogtable_gen_filt),]
bpcoldata3$name <- rownames(bpcoldata3)

bpcoldata3$condition_cba <- gsub("\\s+\\+\\s+", "_", bpcoldata3$condition_cba)
bpcoldata3$condition_cba <- gsub(" ", "_", bpcoldata3$condition_cba)
bpcoldata3$group <- gsub("\\s+\\+\\s+", "_", bpcoldata3$group)
all(names(cogtable_gen_filt) == rownames(bpcoldata3))

# check if any cell has a zero - need pseduocount if yes
any(cogtable_gen_filt == 0)

# Yes, so dont need to add pseudocount to handle any zeros
# mat_pseudo <- cogtable_gen + 0.5

# CLR per sample (apply across columns)
clr_mat <- apply(cogtable_gen_filt, 2, function(x) {
  log2(x / exp(mean(log(x))))
})

clr_mat


# Fit limma model
library(limma)

condition <- factor(bpcoldata3$condition_cba)
design <- model.matrix(~ 0 + condition)
colnames(design) <- levels(condition)
design

contrast_mat <- makeContrasts(
  S3 = (`CA_S3_AS_t1` + `CA_S3_AS_t2`) - (`CA_AS_t1` + `CA_AS_t2`),
  G1 = (`PHA_G1_AS_t1` + `PHA_G1_AS_t2`) - (`PHA_AS_t1` + `PHA_AS_t2`),
  levels = design
)

fit  <- lmFit(clr_mat, design)
fit2 <- contrasts.fit(fit, contrast_mat)
fit2 <- eBayes(fit2)

# 4. Extract results
res_S3 <- topTable(fit2, coef = "S3", number = Inf, sort.by = "none")
res_G1 <- topTable(fit2, coef = "G1", number = Inf, sort.by = "none")

# 5. Build plot matrix
plot_mat <- cbind(S3 = res_S3$logFC, G1 = res_G1$logFC)
rownames(plot_mat) <- rownames(clr_mat)

padj_mat <- cbind(S3 = res_S3$adj.P.Val, G1 = res_G1$adj.P.Val)
rownames(padj_mat) <- rownames(clr_mat)

sig_mat <- ifelse(padj_mat < 0.05, "*", "")
sig_mat[padj_mat < 0.01]  <- "**"
sig_mat[padj_mat < 0.001] <- "***"

# 6. Heatmap
pheatmap(plot_mat,
         cluster_cols     = FALSE,
         cluster_rows     = TRUE,
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D2",
         color  = colorRampPalette(c("blue", "black", "red"))(100),
         breaks = seq(-2, 2, length.out = 101),
         border_color     = "black",
         display_numbers  = sig_mat,
         number_color     = "black",
         fontsize_row     = 10,
         fontsize_col     = 10,
         angle_col        = 45,
         cellwidth = 100,
         main = "COG category log2FC with limma significance")






library(ANCOMBC)

# ANCOM-BC2 expects a phyloseq or TreeSummarizedExperiment object
library(phyloseq)

# Build phyloseq object from raw counts
otu  <- otu_table(mat_filt, taxa_are_rows = TRUE)
samp <- sample_data(metadata_filt)
rownames(samp) <- metadata_filt$sample
ps   <- phyloseq(otu, samp)

# Run ANCOM-BC2
out <- ancombc2(
  data          = ps,
  fix_formula   = "condition",
  p_adj_method  = "BH",
  prv_cut       = 0,        # no prevalence filter, already filtered
  group         = "condition",
  struc_zero    = FALSE,
  verbose       = TRUE
)

# Extract results
res <- out$res














# assuming your metadata has columns "sample" and "condition"
cond_mat <- sapply(unique(bpcoldata2$condition_cba), function(cond) {
  cols <- rownames(bpcoldata2)[bpcoldata2$condition_cba == cond]
  rowMeans(clr_mat[, cols, drop = FALSE])
})

dim(cond_mat)      
head(cond_mat)

# Condition means already computed in cond_mat
# 4 contrasts: S3 effect at t1, S3 effect at t2, G1 effect at t1, G1 effect at t2

plot_mat <- cbind(
  S3_t1 = cond_mat[, "CA + S3 + AS t1"] - cond_mat[, "CA + AS t1"],
  S3_t2 = cond_mat[, "CA + S3 + AS t2"] - cond_mat[, "CA + AS t2"],
  G1_t1 = cond_mat[, "PHA + G1 + AS t1"] - cond_mat[, "PHA + AS t1"],
  G1_t2 = cond_mat[, "PHA + G1 + AS t2"] - cond_mat[, "PHA + AS t2"]
)

dim(plot_mat) 
plot_mat


drop <- c("-", "Y")
plot_mat <- plot_mat[!rownames(plot_mat) %in% drop,]
plot_mat

library(pheatmap)
pheatmap(plot_mat,
         cluster_cols   = FALSE,
         cluster_rows   = TRUE,
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D2",
         color  = colorRampPalette(c("blue", "black", "red"))(100),
         breaks = seq(-2, 2, length.out = 101),
         border_color = 'black',
         fontsize_row  = 12,
         fontsize_col  = 12,
         angle_col     = 0)



drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
cogtable_gen_filt <- cogtable_gen[,!(names(cogtable_gen) %in% drops)]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(cogtable_gen_filt),]
bpcoldata3$name <- rownames(bpcoldata3)

#### Now do aldex2 on the grouped raw counts 

library(ALDEx2)

# ── Input ───────────────────────────────────────────────────────────────────
# mat: raw integer count matrix, 12 rows (COG categories) × 20 columns (samples)
# metadata: data frame with columns "sample" and "condition"
# ALDEx2 works on RAW counts — do not use your CLR matrix here

# ── 1. Run each contrast separately ─────────────────────────────────────────
# ALDEx2 takes a two-group vector for each contrast

run_aldex2 <- function(mat, metadata, cond1, cond2) {
  # Subset to relevant samples
  cols <- metadata$name[metadata$condition_cba %in% c(cond1, cond2)]
  conds <- metadata$condition_cba[metadata$condition_cba %in% c(cond1, cond2)]
  
  mat_sub <- mat[, cols]
  
  # ALDEx2 expects conditions as a character vector in same order as columns
  aldex_out <- aldex(mat_sub,
                     conditions = as.character(conds),
                     mc.samples = 1000,   # increase to 1000 for publication
                     test = "t",         # uses Welch t-test on CLR internally
                     effect = TRUE,
                     denom = "all")      # uses all features as reference (CLR)
  return(aldex_out)
}

# ── 2. Run all four contrasts ────────────────────────────────────────────────
res_S3_t1 <- run_aldex2(cogtable_gen_filt, bpcoldata3, "CA + S3 + AS t1", "CA + AS t1")
res_S3_t2 <- run_aldex2(cogtable_gen_filt, bpcoldata3, "CA + S3 + AS t2", "CA + AS t2")
res_G1_t1 <- run_aldex2(cogtable_gen_filt, bpcoldata3, "PHA + G1 + AS t1", "PHA + AS t1")
res_G1_t2 <- run_aldex2(cogtable_gen_filt, bpcoldata3, "PHA + G1 + AS t2", "PHA + AS t2")

# ── 3. Extract adjusted p-values and effect sizes ───────────────────────────
# wi.eBH = Wilcoxon BH-adjusted p-value
# we.eBH = Welch t-test BH-adjusted p-value  
# effect = effect size (like Cohen's d on CLR scale)

padj_mat <- cbind(
  S3_t1 = res_S3_t1$we.eBH,
  S3_t2 = res_S3_t2$we.eBH,
  G1_t1 = res_G1_t1$we.eBH,
  G1_t2 = res_G1_t2$we.eBH
)
rownames(padj_mat) <- rownames(cogtable_gen_filt)

effect_mat <- cbind(
  S3_t1 = res_S3_t1$effect,
  S3_t2 = res_S3_t2$effect,
  G1_t1 = res_G1_t1$effect,
  G1_t2 = res_G1_t2$effect
)
rownames(effect_mat) <- rownames(cogtable_gen_filt)

# ── 4. Significance stars ────────────────────────────────────────────────────
sig_mat <- ifelse(padj_mat < 0.05, "*", "")
sig_mat[padj_mat < 0.01] <- "**"
sig_mat[padj_mat < 0.001] <- "***"

# ── 5. Heatmap using CLR fold changes + ALDEx2 significance ─────────────────
# plot_mat is your log2FC matrix computed earlier from CLR condition means
library(pheatmap)
pheatmap(plot_mat,
         cluster_cols     = FALSE,
         cluster_rows     = TRUE,
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D2",
         color  = colorRampPalette(c("blue", "black", "red"))(100),
         breaks = seq(-2, 2, length.out = 101),
         border_color     = NA,
         display_numbers  = sig_mat,
         number_color     = "black",
         fontsize_row     = 10,
         fontsize_col     = 10,
         angle_col        = 45,
         main = "COG category log2FC with ALDEx2 significance")


# nothing was significant = try without timepoints

# ── Input ───────────────────────────────────────────────────────────────────
# mat: raw integer count matrix, 12 rows (COG categories) × 20 columns (samples)
# metadata: data frame with columns "sample" and "condition"
# ALDEx2 works on RAW counts — do not use your CLR matrix here

# ── 1. Run each contrast separately ─────────────────────────────────────────
# ALDEx2 takes a two-group vector for each contrast

run_aldex2 <- function(mat, metadata, cond1, cond2) {
  # Subset to relevant samples
  cols <- metadata$name[metadata$group %in% c(cond1, cond2)]
  conds <- metadata$group[metadata$group %in% c(cond1, cond2)]
  
  mat_sub <- mat[, cols]
  
  # ALDEx2 expects conditions as a character vector in same order as columns
  aldex_out <- aldex(mat_sub,
                     conditions = as.character(conds),
                     mc.samples = 1000,   # increase to 1000 for publication
                     test = "t",         # uses Welch t-test on CLR internally
                     effect = TRUE,
                     denom = "all")      # uses all features as reference (CLR)
  return(aldex_out)
}

# ── 2. Run all four contrasts ────────────────────────────────────────────────
res_S3 <- run_aldex2(cogtable_gen_filt, bpcoldata3, "CA + S3 + AS", "CA + AS")
res_G1 <- run_aldex2(cogtable_gen_filt, bpcoldata3, "PHA + G1 + AS", "PHA + AS")

cogtable_gen_filt2 <- cogtable_gen_filt[-3,]
res_S3b <- run_aldex2(cogtable_gen_filt2, bpcoldata3, "CA + S3 + AS", "CA + AS")
res_G1b <- run_aldex2(cogtable_gen_filt2, bpcoldata3, "PHA + G1 + AS", "PHA + AS")

# drop1 <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
# cogtable_gen_filt <- cogtable_gen[,!colnames(cogtable_gen) %in% drop1]
# res_S3 <- run_aldex2(cogtable_gen_filt, bpcoldata3, "CA + S3 + AS", "CA + AS")
# res_G1 <- run_aldex2(cogtable_gen_filt, bpcoldata3, "PHA + G1 + AS", "PHA + AS")

# ── 3. Extract adjusted p-values and effect sizes ───────────────────────────
# wi.eBH = Wilcoxon BH-adjusted p-value
# we.eBH = Welch t-test BH-adjusted p-value  
# effect = effect size (like Cohen's d on CLR scale)

padj_mat <- cbind(
  S3 = res_S3$we.eBH,
  G1 = res_G1$we.eBH
)
rownames(padj_mat) <- rownames(cogtable_gen_filt)

effect_mat <- cbind(
  S3 = res_S3$effect,
  G1 = res_G1$effect
)
rownames(effect_mat) <- rownames(cogtable_gen_filt)

plot_mat <- cbind(
  S3 = res_S3$diff.btw,
  G1 = res_G1$diff.btw
)
rownames(plot_mat) <- rownames(cogtable_gen_filt)

# ── 4. Significance stars ────────────────────────────────────────────────────
sig_mat <- ifelse(padj_mat < 0.05, "*", "")
sig_mat[padj_mat < 0.01] <- "**"
sig_mat[padj_mat < 0.001] <- "***"
sig_mat <- sig_mat[rownames(plot_mat), colnames(plot_mat)]

sig_mat <- matrix("", nrow = nrow(plot_mat), ncol = ncol(plot_mat),
                  dimnames = dimnames(plot_mat))

# Large effect, not significant — dagger symbol
sig_mat[abs(effect_mat) > 1 & padj_mat > 0.05] <- "†"

# Large effect AND significant — star
sig_mat[abs(effect_mat) > 1 & padj_mat < 0.05] <- "*"

# ── 5. Heatmap using CLR fold changes + ALDEx2 significance ─────────────────
# plot_mat is your log2FC matrix computed earlier from CLR condition means
library(pheatmap)
pheatmap(plot_mat,
         cluster_cols     = FALSE,
         cluster_rows     = TRUE,
         clustering_distance_rows = "correlation",
         clustering_method = "ward.D2",
         color  = colorRampPalette(c("blue", "black", "red"))(100),
         breaks = seq(-1, 1, length.out = 101),
         border_color     = "black",
         display_numbers  = sig_mat,
         number_color     = "white",
         fontsize_row     = 10,
         fontsize_col     = 10,
         angle_col        = 45,
         main = "COG category log2FC with ALDEx2 significance")

# aldex very conservative and not great for low sample size - almost nothing is significant

# let's try limma test on my clr transformed data





# try ancombc
# filter out samples/COGs I dont want
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC", "-", "Y")
cogtable_gen_filt <- cogtable_gen[!rownames(cogtable_gen) %in% drops,]
cogtable_gen_filt <- cogtable_gen_filt[,!names(cogtable_gen_filt) %in% drops]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(cogtable_gen_filt),]
bpcoldata3$name <- rownames(bpcoldata3)

bpcoldata3$condition_cba <- gsub("\\s+\\+\\s+", "_", bpcoldata3$condition_cba)
bpcoldata3$condition_cba <- gsub(" ", "_", bpcoldata3$condition_cba)
bpcoldata3$group <- gsub("\\s+\\+\\s+", "_", bpcoldata3$group)

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
  
  # Run ANCOM-BC2
  out <- ancombc2(
    data         = ps,
    fix_formula  = "contrast",
    group        = "contrast",
    p_adj_method = "BH",
    prv_cut      = 0,
    verbose      = TRUE
  )
  
  return(out$res)
}

# ── Run two contrasts ─────────────────────────────────────────────────────────
res_S3 <- run_ancombc(cogtable_gen_filt, bpcoldata3,
                      cond1_labels = "CA_S3_AS",
                      cond2_labels = "CA_AS")

res_G1 <- run_ancombc(cogtable_gen_filt, bpcoldata3,
                      cond1_labels = "PHA_G1_AS",
                      cond2_labels = "PHA_AS")

# ── Extract results ───────────────────────────────────────────────────────────
# lfc = log fold change (treatment vs control)
# q_contrasttreatment = BH adjusted p-value

plot_mat <- cbind(
  S3 = res_S3$lfc_contrasttreatment,
  G1 = res_G1$lfc_contrasttreatment
)
rownames(plot_mat) <- res_S3$taxon

padj_mat <- cbind(
  S3 = res_S3$q_contrasttreatment,
  G1 = res_G1$q_contrasttreatment
)
rownames(padj_mat) <- res_S3$taxon

# ── Significance stars ────────────────────────────────────────────────────────
sig_mat <- ifelse(padj_mat < 0.05, "*", "")
sig_mat[padj_mat < 0.01]  <- "**"
sig_mat[padj_mat < 0.001] <- "***"

# Set logFC to 0 where padj > 0.05
plot_mat_masked <- plot_mat
plot_mat_masked[padj_mat > 0.05] <- 0

# ── Heatmap ───────────────────────────────────────────────────────────────────
library(pheatmap)
# COG category definitions
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

# Replace rownames with full definitions
rownames(plot_mat_masked) <- cog_definitions[rownames(plot_mat_masked)]

# Compute clustering on original logFC values
row_order <- hclust(dist(plot_mat), method = "ward.D2")
colnames(plot_mat_masked) <- c("CA + S3 + AS vs.\nCA + AS", 
                               "PHA + G1 + AS vs.\nPHA + AS")

# ── 1. Define domain groupings ─────────────────────────────────────────────
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

# ── 2. Subset to only COGs present in your matrix ──────────────────────────
# rownames of plot_mat_masked are now full definitions e.g. "J - Translation..."
# extract the letter from the first character
cog_letters <- substr(rownames(plot_mat_masked), 1, 1)
row_annot <- data.frame(
  Category = cog_domains[cog_letters, "Category"],
  row.names = rownames(plot_mat_masked)
)

# ── 3. Define domain colors ────────────────────────────────────────────────
domain_colors <- list(
  Category = c(
    "Information storage & processing" = "#4E9BB5",
    "Cellular processes & signaling"   = "#E07B39",
    "Metabolism"                        = "#5BA75B",
    "Poorly characterized"              = "#9B6BB5"
  )
)


cogph <- pheatmap(plot_mat_masked,
         cluster_cols     = FALSE,
         cluster_rows     = row_order,
         clustering_distance_rows = "correlation",
         clustering_method        = "ward.D2",
         color        = colorRampPalette(c("blue", "black", "red"))(100),
         breaks       = seq(-1, 1, length.out = 101),
         border_color = "grey",
         annotation_row    = row_annot,
         annotation_colors = domain_colors,
         annotation_names_row = FALSE,
         fontsize          = 12,
         fontsize_row      = 10,
         angle_col       = 0,
         cellwidth         = 120,
         cellheight        = 15)

cogph

library(ggplot2)
ggsave(plot = cogph, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/COG_ancombc2_heatmap_04292026.pdf", width = 12.5, height = 6, units = "in", dpi = 300)
ggsave(plot = cogph, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/COG_ancombc2_heatmap_04292026.tiff", width = 12.5, height = 6, units = "in", dpi = 200, bg="white")


# do ancombc1

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

# ── Run two contrasts ─────────────────────────────────────────────────────────
res_S3 <- run_ancombc(cogtable_gen_filt, bpcoldata3,
                      cond1_labels = "CA_S3_AS",
                      cond2_labels = "CA_AS")

res_G1 <- run_ancombc(cogtable_gen_filt, bpcoldata3,
                      cond1_labels = "PHA_G1_AS",
                      cond2_labels = "PHA_AS")

# ── Extract results ───────────────────────────────────────────────────────────
# ANCOM-BC stores results differently from ANCOM-BC2
plot_mat <- cbind(
  S3 = res_S3$lfc,
  G1 = res_G1$lfc
)
rownames(plot_mat) <- rownames(res_S3$lfc$taxon)
rownames(plot_mat) <- plot_mat$S3.taxon
plot_mat <- plot_mat[,-c(1,2,4,5)]

padj_mat <- cbind(
  S3 = res_S3$q_val,
  G1 = res_G1$q_val
)
rownames(padj_mat) <- rownames(res_S3$q_val$taxon)
rownames(padj_mat) <- padj_mat$S3.taxon
padj_mat <- padj_mat[,-c(1,2,4,5)]

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

# ANCOMBC1 works well but malte said to also try maaslin2
# 05/04/2026

library(Maaslin2)
library(dplyr)
cogtable_gen <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/cogtable_gen_raw_04292026.tsv')

cogtable_gen <- cogtable_gen %>%
  column_to_rownames("COGid")

colSums(cogtable_gen)

# filter out samples/COGs I dont want
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC", "-", "Y")
cogtable_gen_filt <- cogtable_gen[!rownames(cogtable_gen) %in% drops,]
cogtable_gen_filt <- cogtable_gen_filt[,!names(cogtable_gen_filt) %in% drops]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(cogtable_gen_filt),]
bpcoldata3$name <- rownames(bpcoldata3)

bpcoldata3$condition_cba <- gsub("\\s+\\+\\s+", "_", bpcoldata3$condition_cba)
bpcoldata3$condition_cba <- gsub(" ", "_", bpcoldata3$condition_cba)
bpcoldata3$group <- gsub("\\s+\\+\\s+", "_", bpcoldata3$group)
all(names(cogtable_gen_filt) == rownames(bpcoldata3))

# Transpose so rows = samples, columns = features
cogtable_gen_filt_t <- as.data.frame(t(cogtable_gen_filt))

# --- Comparison 1: CA_S3_AS vs CA_AS ---
meta1 <- bpcoldata3 |> filter(group %in% c("CA_S3_AS", "CA_AS"))
data1 <- cogtable_gen_filt_t[rownames(meta1), ]

fit1 <- Maaslin2(
  input_data      = data1,
  input_metadata  = meta1,
  output          = "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/maaslin2_CA_05022026",
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
data2 <- cogtable_gen_filt_t[rownames(meta2), ]

fit2 <- Maaslin2(
  input_data      = data2,
  input_metadata  = meta2,
  output          = "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/maaslin2/maaslin2_PHA_05022026",
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

row_annot <- data.frame(
  "COG Category" = cog_domains[cog_letters, "Category"],
  row.names = rownames(plot_mat_masked),
  check.names = FALSE
)

domain_colors <- list("COG Category" = c(
  "Information storage & processing" = "#4E9BB5",
  "Cellular processes & signaling"   = "#E07B39",
  "Metabolism"                        = "#5BA75B",
  "Poorly characterized"              = "#9B6BB5"
))


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























bpselex <- cogtable_gen_filt
drop <- "-"
bpselex <- bpselex[!(row.names(bpselex) %in% drop), ]
names(cogtable_gen_filt) == rownames(bpcoldata3)

#bpx <- aldex(bpselex, bpconds, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, denom="iqlr", verbose=TRUE)

bpcoldata3_sub <- bpcoldata3[bpcoldata3$group == "CA + AS" | bpcoldata3$group == "CA + S3 + AS",]

cogtable_gen_filt_sub <- cogtable_gen_filt[,colnames(cogtable_gen_filt) %in% rownames(bpcoldata3_sub)]
cogtable_gen_filt_sub_noNA <- cogtable_gen_filt_sub[!(row.names(cogtable_gen_filt_sub) %in% drop), ]

colnames(cogtable_gen_filt_sub) == rownames(bpcoldata3_sub)
colnames(cogtable_gen_filt_sub_noNA) == rownames(bpcoldata3_sub)

condsbp <- bpcoldata3_sub$group
condsbp

bp <- aldex(cogtable_gen_filt_sub, condsbp, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, denom="iqlr", verbose=FALSE)
bp <- aldex(cogtable_gen_filt_sub_noNA, condsbp, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, denom="iqlr", verbose=FALSE)

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
