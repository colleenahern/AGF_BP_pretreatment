# 03/02/2026
# C. Ahern
# Using one-way ANOVA, Kruskall-Wallis, and Wilcoxin tests to determine significance of my low-feature (COG, CAZyme) condition_cbaed, 
# JGI-normalized counts
# Test both the tcombined and t1t2 results (break up PHA + G1 + AS t2 and PHA + G1 + AS t1) for my analysis because they seen quite different on the PCA

# Import my condition_cbaed TPM-normalized counts
library(tidyverse)
library(readxl)

# Import cog table from BP_MT_transcript_abund_grouped_COGs_03012026.R
# cogtable_genr from line 196

cogtable_genr <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/cogtable_jginorm_03022026.tsv")

cogtable_genr <- cogtable_genr %>%
  column_to_rownames("COGid")

# Drop any outlier samples you have (based on PCA analysis)

# bb = dbcantable_gen
# drop <- "-"
# bb <- bb[!(row.names(bb) %in% drop), ]
# drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
# bb_filt <- bb[,!(names(bb) %in% drops)]
# bb_filt$cat <- rownames(bb_filt)

bb = cogtable_genr
drop <- "-"
bb <- bb[!(row.names(bb) %in% drop), ]
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
bb_filt <- bb[,!(names(bb) %in% drops)]
bb_filt$cat <- rownames(bb_filt)

# import and filter metadata
bpcoldata <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/bp_metadata_01132025.txt")
# bpcoldata$name <- gsub("-", ".", bpcoldata$name)
bpcoldata3 <- bpcoldata[bpcoldata$name %in% names(bb_filt),]
# bpcoldata3$group2 <- bpcoldata3$condition_cba
# bpcoldata3$group2[c(1:2,5:6)] <- "CA + AS"
# bpcoldata3$group2[c(3,4,7,8)] <- "CA + S3 + AS"
# bpcoldata3$group2[c(9,10,14,15)] <- "PHA + AS"

meta = bpcoldata3
ff = bb_filt %>%
  pivot_longer(-cat) %>% 
  left_join(meta)

############################################################################################################################################
# One way anova test

# COG groups test
result <- aov(value ~ group, data = subset(ff, cat == "N")) # PHA + G1 + AS vs. PHA + AS padj = 0.0000490 (lower than)
summary(result)                                             # CA + S3 + AS vs. CA + AS padj = 0.2190922 (not different)
TukeyHSD(result)

result <- aov(value ~ condition_cba, data = subset(ff, cat == "N")) # PHA + G1 + AS t1-PHA + AS t1 padj = 0.0110926 (lower than); PHA + G1 + AS t2-PHA + AS t2 padj = 0.0421331
summary(result)                                             # CA + S3 + AS t1-CA + AS t1 padj = 0.9115042 (not different); CA + S3 + AS t2-CA + AS t2 padj = 0.8349689
TukeyHSD(result)

result <- aov(value ~ group, data = subset(ff, cat == "C")) # PHA + G1 + AS t1 vs. PHA + AS t1 padj = 0.0002871
summary(result)                                                     # PHA + G1 + AS t2 vs. PHA + AS t2 padj = 0.9999999
TukeyHSD(result)                                                    # CA + S3 + AS t2 vs. CA + AS t2 padj = 0.2025365
# CA + S3 + AS t1 vs. CA + AS t1 padj = 0.0032633

result <- aov(value ~ condition_cba, data = subset(ff, cat == "C")) # PHA + G1 + AS t1-PHA + AS t1  padj = 0.0009425
summary(result)                                                     # CA + S3 + AS t1-CA + AS t1 padj = 0.0415266
TukeyHSD(result)  

result <- aov(value ~ condition_cba, data = subset(ff, cat == "E")) 
summary(result)                                                     
TukeyHSD(result)  

result <- aov(value ~ group, data = subset(ff, cat == "E")) # PHA + G1 + AS t1-PHA + AS t1  padj = 0.0009425
summary(result)                                                     # CA + S3 + AS t1-CA + AS t1 padj = 0.0415266
TukeyHSD(result)  


############################################################################################################################################
# Kruskall/dunn test

kruskal.test(value ~ condition_cba, data = subset(ff, cat == "C"))

kruskal.test(value ~ condition_cba, data = subset(ff, cat == "GH"))
library(FSA)
install.packages('FSA')
dunnTest(value ~ condition_cba, data = subset(ff, cat == "GH"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance
dunnTest(value ~ condition_cba, data = subset(ff, cat == "CE"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance
dunnTest(value ~ condition_cba, data = subset(ff, cat == "GT"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance
dunnTest(value ~ condition_cba, data = subset(ff, cat == "CBM"), method = "bh") # no CA + S3 + AS vs. CA + AS significance BUT PHA + G1 + AS vs. PHA + AS significance
dunnTest(value ~ condition_cba, data = subset(ff, cat == "PL"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance
dunnTest(value ~ condition_cba, data = subset(ff, cat == "AA"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance

dunnTest(value ~ condition_cba, data = subset(ff, cat == "C"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance

#############################################################################################################################################
# Wilcoxin test

ffsub <- subset(ff, cat == "C")
pairwise.wilcox.test(ffsub$value, ffsub$group,     
                     p.adjust.method = "BH")   

ffsub <- subset(ff, cat == "C")
pairwise.wilcox.test(ffsub$value, ffsub$condition_cba,     
                     p.adjust.method = "BH")  

ffsub <- subset(ff, cat == "N")
pairwise.wilcox.test(ffsub$value, ffsub$group,    # CA + S3 + AS vs CA + AS significance!
                     p.adjust.method = "BH")      # no PHA + G1 + AS vs. PHA + AS significance

ffsub <- subset(ff, cat == "GT")
pairwise.wilcox.test(ffsub$value, ffsub$condition_cba,    # no CA + S3 + AS vs CA + AS significance
                     p.adjust.method = "BH")      # no PHA + G1 + AS vs. PHA + AS significance

ffsub <- subset(ff, cat == "CBM")
pairwise.wilcox.test(ffsub$value, ffsub$condition_cba,    # no CA + S3 + AS vs CA + AS significance
                     p.adjust.method = "BH")      # but PHA + G1 + AS vs. PHA + AS significance

ffsub <- subset(ff, cat == "PL")
pairwise.wilcox.test(ffsub$value, ffsub$condition_cba,   # no CA + S3 + AS vs CA + AS significance
                     p.adjust.method = "BH")    # no PHA + G1 + AS vs. PHA + AS significance

ffsub <- subset(ff, cat == "AA") 
pairwise.wilcox.test(ffsub$value, ffsub$condition_cba,    # no CA + S3 + AS vs CA + AS significance
                     p.adjust.method = "BH")      # no PHA + G1 + AS vs. PHA + AS significance





t.test(value ~ condition_cba, subset(ff, cat == "C" & 
         (condition_cba == "CA + AS t1" | 
            condition_cba == "CA + S3 + AS t1")))

t.test(value ~ condition_cba, subset(ff, cat == "C" & 
                                       (condition_cba == "PHA + AS t1" | 
                                          condition_cba == "PHA + G1 + AS t1")))

t.test(value ~ group, subset(ff, cat == "C" & 
                                       (group == "PHA + AS" | 
                                          group == "PHA + G1 + AS")))

t.test(value ~ group, subset(ff, cat == "C" & 
                               (group == "CA + AS" | 
                                  group == "CA + S3 + AS")))
