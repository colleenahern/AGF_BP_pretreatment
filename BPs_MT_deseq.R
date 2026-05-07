# 01/08/2024
# Complete DESeq analysis on gene count matrix from Malte for BPs experiment metatranscriptomic data

library("readr")
library(tidyverse)

## 01/13/2025 picking up here
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

bpcoldata2$time <- paste("t", bpcoldata2$timepoint, sep = "")
bpcoldata2$group2 <- bpcoldata2$group

## With the count matrix and the sample info, construct a DESeqDataSet:
library("DESeq2")
bpdds <- DESeqDataSetFromMatrix(countData = bpcts,
                                colData = bpcoldata2,
                                design = ~ condition_cba)
bpdds

vsdata <- vst(bpdds, blind=FALSE)
z <- plotPCA(vsdata, intgroup="condition_cba", returnData = TRUE)
# z + labs(color = "Condition") + theme_classic()
# z + geom_label(aes(label = name))


# Custom plot using ggplot2
library(ggplot2)
ggplot(z, aes(x = PC1, y = PC2, color = group2, shape = time)) +
  geom_point(size = 3) +
  labs(color = "Condition", shape = "Timepoint", x = "PC1: 52% variance", y = "PC2: 24% variance") +
  theme_classic()

#########
## do i remove my outliers for publication or keep??
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
bpctsf <- bpcts[,!(names(bpcts) %in% drops)]
bpcoldata2f <- bpcoldata2[!(row.names(bpcoldata2) %in% drops),]

all(rownames(bpcoldata2f) == colnames(bpctsf))

bpcoldata2f$time <- paste("t", bpcoldata2f$timepoint, sep = "")
bpcoldata2f$group2 <- bpcoldata2f$group

## With the count matrix and the sample info, construct a DESeqDataSet:
library("DESeq2")
bpdds <- DESeqDataSetFromMatrix(countData = bpctsf,
                                colData = bpcoldata2f,
                                design = ~ condition_cba)
bpdds

vsdata <- vst(bpdds, blind=FALSE)
z <- plotPCA(vsdata, intgroup="condition_cba", returnData = TRUE)
z <- plotPCA(vsdata, intgroup="condition_cba")
z
# z + labs(color = "Condition") + theme_classic()
# z + geom_label(aes(label = name))

# Custom plot using ggplot2
library(ggplot2)
ggplot(z, aes(x = PC1, y = PC2, color = group2, shape = time)) +
  geom_point(size = 3) +
  labs(color = "Condition", shape = "Timepoint", x = "PC1: 58% variance", y = "PC2: 24% variance") +
  theme_classic()






## Pre-filter
keepbp <- rowSums(counts(bpdds)) > 0
bpdds <- bpdds[keepbp,]   # reduced number of genes from 845606 to 791656

## Set the factor levels 
#bpdds$groupName <- factor(bpdds$groupName, levels = c("CA + AS t1","CA + S3 + AS t1","CA + AS t2", "CA + S3 + AS t2", "PHA + AS t1", "PHA + G1 + AS t1", "PHA + AS t2", "PHA + G1 + AS t2"))

## Differential expression analysis: compare each substrate to glucose
bpdds <- DESeq(bpdds)
bpres <- results(bpdds)
bpres
bpres <- results(bpdds, contrast=c("condition_cba","CA + S3 + AS t1",	"CA + AS t1"))
bpres

########################################################################################################################
## Do the S3/CA comparisons: 4 in total
## CA + S3 + AS t1 vs CA + AS t1
bpres_CAS3ASt1_CAASt1 <- results(bpdds, contrast=c("condition_cba","CA + S3 + AS t1","CA + AS t1"))
bpres_CAS3ASt1_CAASt1
bpresmat_CAS3ASt1_CAASt1 <- as.matrix(bpres_CAS3ASt1_CAASt1)
bpresmat_CAS3ASt1_CAASt1
colnames(bpresmat_CAS3ASt1_CAASt1) <- paste("CA.S3.AS.t1_v_CA.AS.t1", colnames(bpresmat_CAS3ASt1_CAASt1), sep = "_")
bpresmat_CAS3ASt1_CAASt1 <- as.data.frame(bpresmat_CAS3ASt1_CAASt1)
sum(bpresmat_CAS3ASt1_CAASt1$CA.S3.AS.t1_v_CA.AS.t1_padj < 0.05, na.rm=TRUE)

## CA + S3 + AS t2 vs 	CA + AS t2
bpres_CAS3ASt2_CAASt2 <- results(bpdds, contrast=c("condition_cba","CA + S3 + AS t2", "CA + AS t2"))
bpres_CAS3ASt2_CAASt2
bpresmat_CAS3ASt2_CAASt2 <- as.matrix(bpres_CAS3ASt2_CAASt2)
bpresmat_CAS3ASt2_CAASt2
colnames(bpresmat_CAS3ASt2_CAASt2) <- paste("CA.S3.AS.t2_v_CA.AS.t2", colnames(bpresmat_CAS3ASt2_CAASt2), sep = "_")
bpresmat_CAS3ASt2_CAASt2 <- as.data.frame(bpresmat_CAS3ASt2_CAASt2)
sum(bpresmat_CAS3ASt2_CAASt2$CA.S3.AS.t2_v_CA.AS.t2_padj < 0.05, na.rm=TRUE)

## CA + AS t2 vs CA + AS t1
bpres_CAASt2_CAASt1 <- results(bpdds, contrast=c("condition_cba","CA + AS t2","CA + AS t1"))
bpres_CAASt2_CAASt1
bpresmat_CAASt2_CAASt1 <- as.matrix(bpres_CAASt2_CAASt1)
bpresmat_CAASt2_CAASt1
colnames(bpresmat_CAASt2_CAASt1) <- paste("CA.AS.t2_v_CA.AS.t1", colnames(bpresmat_CAASt2_CAASt1), sep = "_")
bpresmat_CAASt2_CAASt1 <- as.data.frame(bpresmat_CAASt2_CAASt1)
sum(bpresmat_CAASt2_CAASt1$CA.AS.t2_v_CA.AS.t1_padj < 0.05, na.rm=TRUE)

## CA + S3 + AS t2 vs 	CA + S3 + AS t1
bpres_CAS3ASt2_CAS3ASt1 <- results(bpdds, contrast=c("condition_cba","CA + S3 + AS t2", "CA + S3 + AS t1"))
bpres_CAS3ASt2_CAS3ASt1
bpresmat_CAS3ASt2_CAS3ASt1 <- as.matrix(bpres_CAS3ASt2_CAS3ASt1)
bpresmat_CAS3ASt2_CAS3ASt1
colnames(bpresmat_CAS3ASt2_CAS3ASt1) <- paste("CA.S3.AS.t2_v_CA.S3.AS.t1", colnames(bpresmat_CAS3ASt2_CAS3ASt1), sep = "_")
bpresmat_CAS3ASt2_CAS3ASt1 <- as.data.frame(bpresmat_CAS3ASt2_CAS3ASt1)
sum(bpresmat_CAS3ASt2_CAS3ASt1$CA.S3.AS.t2_v_CA.S3.AS.t1_padj < 0.05, na.rm=TRUE)

########################################################################################################################
## Now do the G1/PHA comparisons: 4 in total

## PHA + G1 + AS t1 vs PHA + AS t1
bpres_PHAG1ASt1_PHAASt1 <- results(bpdds, contrast=c("condition_cba","PHA + G1 + AS t1","PHA + AS t1"))
bpres_PHAG1ASt1_PHAASt1
bpresmat_PHAG1ASt1_PHAASt1 <- as.matrix(bpres_PHAG1ASt1_PHAASt1)
bpresmat_PHAG1ASt1_PHAASt1
colnames(bpresmat_PHAG1ASt1_PHAASt1) <- paste("PHA.G1.AS.t1_v_PHA.AS.t1", colnames(bpresmat_PHAG1ASt1_PHAASt1), sep = "_")
bpresmat_PHAG1ASt1_PHAASt1 <- as.data.frame(bpresmat_PHAG1ASt1_PHAASt1)
sum(bpresmat_PHAG1ASt1_PHAASt1$PHA.G1.AS.t1_v_PHA.AS.t1_padj < 0.05, na.rm=TRUE)

## PHA + G1 + AS t2 vs 	PHA + AS t2
bpres_PHAG1ASt2_PHAASt2 <- results(bpdds, contrast=c("condition_cba","PHA + G1 + AS t2", "PHA + AS t2"))
bpres_PHAG1ASt2_PHAASt2
bpresmat_PHAG1ASt2_PHAASt2 <- as.matrix(bpres_PHAG1ASt2_PHAASt2)
bpresmat_PHAG1ASt2_PHAASt2
colnames(bpresmat_PHAG1ASt2_PHAASt2) <- paste("PHA.G1.AS.t2_v_PHA.AS.t2", colnames(bpresmat_PHAG1ASt2_PHAASt2), sep = "_")
bpresmat_PHAG1ASt2_PHAASt2 <- as.data.frame(bpresmat_PHAG1ASt2_PHAASt2)
sum(bpresmat_PHAG1ASt2_PHAASt2$PHA.G1.AS.t2_v_PHA.AS.t2_padj < 0.05, na.rm=TRUE)

## PHA + AS t2 vs PHA + AS t1
bpres_PHAASt2_PHAASt1 <- results(bpdds, contrast=c("condition_cba","PHA + AS t2","PHA + AS t1"))
bpres_PHAASt2_PHAASt1
bpresmat_PHAASt2_PHAASt1 <- as.matrix(bpres_PHAASt2_PHAASt1)
bpresmat_PHAASt2_PHAASt1
colnames(bpresmat_PHAASt2_PHAASt1) <- paste("PHA.AS.t2_v_PHA.AS.t1", colnames(bpresmat_PHAASt2_PHAASt1), sep = "_")
bpresmat_PHAASt2_PHAASt1 <- as.data.frame(bpresmat_PHAASt2_PHAASt1)
sum(bpresmat_PHAASt2_PHAASt1$PHA.AS.t2_v_PHA.AS.t1_padj < 0.05, na.rm=TRUE)

## PHA + G1 + AS t2 vs 	PHA + G1 + AS t1
bpres_PHAG1ASt2_PHAG1ASt1 <- results(bpdds, contrast=c("condition_cba","PHA + G1 + AS t2", "PHA + G1 + AS t1"))
bpres_PHAG1ASt2_PHAG1ASt1
bpresmat_PHAG1ASt2_PHAG1ASt1 <- as.matrix(bpres_PHAG1ASt2_PHAG1ASt1)
bpresmat_PHAG1ASt2_PHAG1ASt1
colnames(bpresmat_PHAG1ASt2_PHAG1ASt1) <- paste("PHA.G1.AS.t2_v_PHA.G1.AS.t1", colnames(bpresmat_PHAG1ASt2_PHAG1ASt1), sep = "_")
bpresmat_PHAG1ASt2_PHAG1ASt1 <- as.data.frame(bpresmat_PHAG1ASt2_PHAG1ASt1)
sum(bpresmat_PHAG1ASt2_PHAG1ASt1$PHA.G1.AS.t2_v_PHA.G1.AS.t1_padj < 0.05, na.rm=TRUE)

#######################################################################################################################

## Confirm that the order of genes in all 6 matrices is the same
all(rownames(bpresmat_CAS3ASt1_CAASt1) == rownames(bpresmat_CAS3ASt2_CAASt2))
all(rownames(bpresmat_CAS3ASt2_CAASt2) == rownames(bpresmat_CAASt2_CAASt1))
all(rownames(bpresmat_CAASt2_CAASt1) == rownames(bpresmat_CAASt2_CAASt1))
all(rownames(bpresmat_CAASt2_CAASt1) == rownames(bpresmat_CAS3ASt2_CAS3ASt1))
all(rownames(bpresmat_CAS3ASt2_CAS3ASt1) == rownames(bpresmat_PHAG1ASt1_PHAASt1))
all(rownames(bpresmat_PHAG1ASt1_PHAASt1) == rownames(bpresmat_PHAG1ASt2_PHAASt2))
all(rownames(bpresmat_PHAG1ASt2_PHAASt2) == rownames(bpresmat_PHAASt2_PHAASt1))
all(rownames(bpresmat_PHAASt2_PHAASt1) == rownames(bpresmat_PHAG1ASt2_PHAG1ASt1))

#####################################################################################################################
## Combine the matrices into one dataframe 
## BPs
BP_DESeq <- cbind(bpresmat_CAS3ASt1_CAASt1, bpresmat_CAS3ASt2_CAASt2)
BP_DESeq$CA.S3.AS.t1_v_CA.AS.t1_Sig <- ifelse(BP_DESeq$CA.S3.AS.t1_v_CA.AS.t1_padj <= 0.05,"TRUE","FALSE")
BP_DESeq <- BP_DESeq[,c(1:6,13,7:12)]
BP_DESeq$CA.S3.AS.t2_v_CA.AS.t2_Sig <- ifelse(BP_DESeq$CA.S3.AS.t2_v_CA.AS.t2_padj <= 0.05,"TRUE","FALSE")

BP_DESeq <- cbind(BP_DESeq, bpresmat_CAASt2_CAASt1)
BP_DESeq$CA.AS.t2_v_CA.AS.t1_Sig <- ifelse(BP_DESeq$CA.AS.t2_v_CA.AS.t1_padj <= 0.05,"TRUE","FALSE")

BP_DESeq <- cbind(BP_DESeq, bpresmat_CAS3ASt2_CAS3ASt1)
BP_DESeq$CA.S3.AS.t2_v_CA.S3.AS.t1_Sig <- ifelse(BP_DESeq$CA.S3.AS.t2_v_CA.S3.AS.t1_padj <= 0.05,"TRUE","FALSE")

BP_DESeq <- cbind(BP_DESeq, bpresmat_PHAG1ASt1_PHAASt1)
BP_DESeq$PHA.G1.AS.t1_v_PHA.AS.t1_Sig <- ifelse(BP_DESeq$PHA.G1.AS.t1_v_PHA.AS.t1_padj <= 0.05,"TRUE","FALSE")

BP_DESeq <- cbind(BP_DESeq, bpresmat_PHAG1ASt2_PHAASt2)
BP_DESeq$PHA.G1.AS.t2_v_PHA.AS.t2_Sig <- ifelse(BP_DESeq$PHA.G1.AS.t2_v_PHA.AS.t2_padj <= 0.05,"TRUE","FALSE")

BP_DESeq <- cbind(BP_DESeq, bpresmat_PHAASt2_PHAASt1)
BP_DESeq$PHA.AS.t2_v_PHA.AS.t1_Sig <- ifelse(BP_DESeq$PHA.AS.t2_v_PHA.AS.t1_padj <= 0.05,"TRUE","FALSE")

BP_DESeq <- cbind(BP_DESeq, bpresmat_PHAG1ASt2_PHAG1ASt1)
BP_DESeq$PHA.G1.AS.t2_v_PHA.G1.AS.t1_Sig <- ifelse(BP_DESeq$PHA.G1.AS.t2_v_PHA.G1.AS.t1_padj <= 0.05,"TRUE","FALSE")

BP_DESeq$GeneID <- rownames(BP_DESeq)
BP_DESeq <- BP_DESeq[,c(57,1:56)]

write_csv(BP_DESeq, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_deseq/BP_DESeq_01162025.csv")

a <- colnames(BP_DESeq)
BP_DESeq_filt <- BP_DESeq %>% select(matches("GeneID|FoldChange|padj|Sig", a))
write_csv(BP_DESeq, "/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_deseq/BP_DESeq_filt_01162025.csv")
