# Colleen Ahern
# 02/11/2025

# Continuation of BP_MT_transcript_abund_grouped_02102025 script
# I want to normalize my raw gene counts and add together CAZyme type instead of COGs/KEGG modules
# Continue using TPM normalization

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
bptpm_geneanno$dbcan_annotations_a <- gsub('[[:digit:]]+', '', bptpm_geneanno$dbcan_annotations)

bptpm_geneanno$dbcan_annotations_b <- gsub("GT[[:punct:]]GT","GT",
                                     gsub("CBM[[:punct:]]CBM","CBM",
                                     gsub("GH[[:punct:]]GH","GH", bptpm_geneanno$dbcan_annotations_a)))

library(dplyr)
library(stringr)
library(tidyr)

# debcan annotations
aa <- bptpm_geneanno %>%
  select(GeneID, everything())%>%
  mutate(dbcan_annotations_b=str_split(dbcan_annotations_b, "\\s*\\|\\s*")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:47,1)]

unique(bptpm_geneanno$dbcan_annotations_a)
ab<-as.vector(as.matrix(aa[,grepl( "dbcan_annotations_b" , names(aa))]))
unique(ab)
length(unique(ab))

# make subset of aa data
aasub <- as.data.frame(aa[3253:3258,])
aasub$dbcan_annotations_b1[is.na(aasub$dbcan_annotations_b1)] <- "-"
aasubccs <- aasub[,c(44:46)]

uab <-unique(ab)
uab <- uab[!is.na(uab)]
rn <- append(uab, "-")
rn
dbcantable_gen <- data.frame(matrix(0, nrow = length(unique(rn)), ncol = 20))
colnames(dbcantable_gen) <- colnames(aa[,1:20])
rownames(dbcantable_gen) <- rn

# Trial run on subset of the data
for (j in 1:ncol(dbcantable_gen)) {
  for (k in 1:length(aasub[,grepl( "dbcan_annotations_b" , names(aasub))])) { # Number of separate cazyme types broken up
    for (i in 1:nrow(aasub)) {
      if(is.na(aasub[i,paste("dbcan_annotations_b",k, sep = "")]) == FALSE) {
        dbcantable_gen[rn[rownames(dbcantable_gen) %in% aasub[i,paste("dbcan_annotations_b",k, sep = "")]],j] = dbcantable_gen[rn[rownames(dbcantable_gen) %in% aasub[i,paste("dbcan_annotations_b",k, sep = "")]],j] + aasub[i,j]
      }
      else {
        dbcantable_gen = dbcantable_gen
      }
    }
  }
}


## Now do it on the entire data
aa$dbcan_annotations_b1[is.na(aa$dbcan_annotations_b1)] <- "-"

uab <-unique(ab)
uab <- uab[!is.na(uab)]
rn <- append(uab, "-")
rn
dbcantable_gen <- data.frame(matrix(0, nrow = length(unique(rn)), ncol = 20))
colnames(dbcantable_gen) <- colnames(aa[,1:20])
rownames(dbcantable_gen) <- rn

for (j in 1:ncol(dbcantable_gen)) {
  for (k in 1:length(aa[,grepl( "dbcan_annotations_b" , names(aa))])) { # Number of separate KEGG modules broken up
    for (i in 1:nrow(aa)) {
      if(is.na(aa[i,paste("dbcan_annotations_b",k, sep = "")]) == FALSE) {
        dbcantable_gen[rn[rownames(dbcantable_gen) %in% aa[i,paste("dbcan_annotations_b",k, sep = "")]],j] = dbcantable_gen[rn[rownames(dbcantable_gen) %in% aa[i,paste("dbcan_annotations_b",k, sep = "")]],j] + aa[i,j]
      }
      else {
        dbcantable_gen = dbcantable_gen
      }
    }
  }
}
library(openxlsx)
write.xlsx(dbcantable_gen, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_gen.xlsx", rowNames = TRUE)
library(readxl)
dbcantable_gen <- as.data.frame(read_excel('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_gen.xlsx'))
rownames(dbcantable_gen) <- dbcantable_gen$...1
dbcantable_gen <- dbcantable_gen[, -1]


# Take average of replicates
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
dbcantable_gen_filt <- dbcantable_gen[,!(names(dbcantable_gen) %in% drops)]
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

bp_dbcantpm_heatmap<-pheatmap(data_matrix_44, col=my_palette, cellwidth=55,cellheight=6.6, cluster_rows=F, cluster_cols=F, 
                                fontsize=8,show_rownames=T, fontsize_row=8, fontsize_col=8, labels_row = rownames(data_matrix_44), labels_col = as.character(newnames), border_color=NA,breaks=seq(0, 800, length.out=300), angle_col = 45)


# Now let's do aldex2 or ancom-bc on the dbcantable_gen_filt 

# This is just following a tutorial from online
# https://bioconductor.statistik.tu-dortmund.de/packages/3.6/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.pdf

# tutorial
library(ALDEx2)
data(selex)
#subset for efficiency
selex <- selex[1201:1600,]
conds <- c(rep("NS", 7), rep("S", 7))

x <- aldex(selex, conds, mc.samples=16, test="t", effect=TRUE, include.sample.summary=FALSE, denom="iqlr", verbose=TRUE)
aldex.plot(x, type="MA", test="welch")
aldex.plot(x, type="MW", test="welch")

head(x)

ef1 <- x[x$effect >= 1,]
head(ef1)

##########################################################################################
# Test significance of differential abundance of different groups
# My MT data

library(ALDEx2)
#bpselex <- dbcantable_gen_filt
bpselex <- floor(dbcantable_gen_filt)
drop <- "-"
bpselex <- bpselex[!(row.names(bpselex) %in% drop), ]
names(dbcantable_gen_filt) == rownames(bpcoldata3)
bpconds <- bpcoldata3$condition_cba

bpx <- aldex(bpselex, bpconds, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, denom="iqlr", verbose=TRUE)

#######################################################################################################################################
################ Ok I think differential abundance doesn't work with TPM values, try on the cumulative raw values

library('DRnaSeq')
library(readr)

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
bpcts$GeneID <- rownames(bpcts)

# Input gene annotation info
gene_anno <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_analysis/MT_analysis/processed_tables/gene_annotation_metassembly.tsv")

bpcts_geneanno <- merge(bpcts, gene_anno, by = "GeneID", all = TRUE)
bpcts_geneanno$dbcan_annotations_a <- gsub('[[:digit:]]+', '', bpcts_geneanno$dbcan_annotations)

bpcts_geneanno$dbcan_annotations_b <- gsub("GT[[:punct:]]GT","GT",
                                           gsub("CBM[[:punct:]]CBM","CBM",
                                                gsub("GH[[:punct:]]GH","GH", bpcts_geneanno$dbcan_annotations_a)))

library(dplyr)
library(stringr)
library(tidyr)

# dbcan annotations
aa <- bpcts_geneanno %>%
  select(GeneID, everything())%>%
  mutate(dbcan_annotations_b=str_split(dbcan_annotations_b, "\\s*\\|\\s*")) %>% 
  unnest_wider(where(is.list), names_sep = "")
aa <- as.data.frame(aa)
aa <- aa[,c(2:47,1)]

unique(bpcts_geneanno$dbcan_annotations_a)
unique(bpcts_geneanno$dbcan_annotations_b)
ab<-as.vector(as.matrix(aa[,grepl( "dbcan_annotations_b" , names(aa))]))
unique(ab)
length(unique(ab))


## add raw counts for each cazyme type
aa$dbcan_annotations_b1[is.na(aa$dbcan_annotations_b1)] <- "-"

uab <-unique(ab)
uab <- uab[!is.na(uab)]
rn <- append(uab, "-")
rn
dbcantable_gen_raw <- data.frame(matrix(0, nrow = length(unique(rn)), ncol = 20))
colnames(dbcantable_gen_raw) <- colnames(aa[,1:20])
rownames(dbcantable_gen_raw) <- rn

for (j in 1:ncol(dbcantable_gen_raw)) {
  for (k in 1:length(aa[,grepl( "dbcan_annotations_b" , names(aa))])) { # Number of separate cazyme types broken up
    for (i in 1:nrow(aa)) {
      if(is.na(aa[i,paste("dbcan_annotations_b",k, sep = "")]) == FALSE) {
        dbcantable_gen_raw[rn[rownames(dbcantable_gen_raw) %in% aa[i,paste("dbcan_annotations_b",k, sep = "")]],j] = dbcantable_gen_raw[rn[rownames(dbcantable_gen_raw) %in% aa[i,paste("dbcan_annotations_b",k, sep = "")]],j] + aa[i,j]
      }
      else {
        dbcantable_gen_raw = dbcantable_gen_raw
      }
    }
  }
}
library(openxlsx)
write.xlsx(dbcantable_gen_raw, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_gen_raw.xlsx", rowNames = TRUE)
library(readxl)
dbcantable_gen_raw <- as.data.frame(read_excel('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_gen_raw.xlsx'))
rownames(dbcantable_gen_raw) <- dbcantable_gen_raw$...1
dbcantable_gen_raw <- dbcantable_gen_raw[, -1]

# Take average of replicates
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
dbcantable_gen_raw_filt <- dbcantable_gen_raw[,!(names(dbcantable_gen_raw) %in% drops)]
bpcoldata3 <- bpcoldata2[rownames(bpcoldata2) %in% names(dbcantable_gen_raw_filt),]
rn <- rownames(dbcantable_gen_raw)
rn
dbcantable_gen_raw_avg <- data.frame(matrix(0, nrow = length(rn), ncol = 8))
colnames(dbcantable_gen_raw_avg) <- (unique(bpcoldata3$condition_cba))
rownames(dbcantable_gen_raw_avg) <- rownames(dbcantable_gen_raw_filt)

dbcantable_gen_raw_avg[,'CA + AS t1'] <- rowMeans(as.matrix(dbcantable_gen_raw_filt[, names(dbcantable_gen_raw_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + AS t1',])]))
dbcantable_gen_raw_avg[,'CA + S3 + AS t1'] <- rowMeans(as.matrix(dbcantable_gen_raw_filt[, names(dbcantable_gen_raw_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + S3 + AS t1',])]))
dbcantable_gen_raw_avg[,'CA + AS t2'] <- rowMeans(as.matrix(dbcantable_gen_raw_filt[, names(dbcantable_gen_raw_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + AS t2',])]))
dbcantable_gen_raw_avg[,'CA + S3 + AS t2'] <- rowMeans(as.matrix(dbcantable_gen_raw_filt[, names(dbcantable_gen_raw_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'CA + S3 + AS t2',])]))

dbcantable_gen_raw_avg[,'PHA + AS t1'] <- rowMeans(as.matrix(dbcantable_gen_raw_filt[, names(dbcantable_gen_raw_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + AS t1',])]))
dbcantable_gen_raw_avg[,'PHA + G1 + AS t1'] <- rowMeans(as.matrix(dbcantable_gen_raw_filt[, names(dbcantable_gen_raw_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + G1 + AS t1',])]))
dbcantable_gen_raw_avg[,'PHA + AS t2'] <- rowMeans(as.matrix(dbcantable_gen_raw_filt[, names(dbcantable_gen_raw_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + AS t2',])]))
dbcantable_gen_raw_avg[,'PHA + G1 + AS t2'] <- rowMeans(as.matrix(dbcantable_gen_raw_filt[, names(dbcantable_gen_raw_filt) %in% rownames(bpcoldata3[bpcoldata3$condition_cba == 'PHA + G1 + AS t2',])]))

########################################################################################################################################
# Now do aldex2 on the grouped raw counts 


library(ALDEx2)
bpselex <- dbcantable_gen_raw_filt
drop <- "-"
bpselex <- bpselex[!(row.names(bpselex) %in% drop), ]
names(dbcantable_gen_raw_filt) == rownames(bpcoldata3)

#bpx <- aldex(bpselex, bpconds, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, denom="iqlr", verbose=TRUE)

bpcoldata3raw_sub <- bpcoldata3raw[bpcoldata3raw$group == "CA + AS" | bpcoldata3raw$group == "CA + S3 + AS",]

dbcantable_gen_raw_filt_sub <- dbcantable_gen_raw_filt[,colnames(dbcantable_gen_raw_filt) %in% rownames(bpcoldata3raw_sub)]
dbcantable_gen_raw_filt_sub_noNA <- dbcantable_gen_raw_filt_sub[!(row.names(dbcantable_gen_raw_filt_sub) %in% drop), ]

colnames(dbcantable_gen_raw_filt_sub) == rownames(bpcoldata3raw_sub)
colnames(dbcantable_gen_raw_filt_sub_noNA) == rownames(bpcoldata3raw_sub)

condsbp <- bpcoldata3raw_sub$group

bp <- aldex(dbcantable_gen_raw_filt_sub, condsbp, mc.samples=16, test="t", effect=TRUE, include.sample.summary=FALSE, denom="iqlr", verbose=FALSE)
bp <- aldex(dbcantable_gen_raw_filt_sub_noNA, condsbp, mc.samples=128, test="t", effect=TRUE, include.sample.summary=FALSE, denom="iqlr", verbose=FALSE)
aldex.plot(y, type="MA", test="welch")
aldex.plot(y, type="MW", test="welch")

head(y)

efy1 <- y[y$effect >= 1,]
head(efy1)
efy1$asvid <- row.names(efy1)
tt <- as.data.frame(pswd_film@tax_table)
tt$asvid <- row.names(tt)

efy1tax <- merge(efy1, tt, by = 'asvid')

