# 03/12/2025
# C. Ahern
# Using one-way ANOVA, Kruskall-Wallis, and Wilcoxin tests to determine significance of my low-feature (COG, CAZyme) grouped, 
# TPM-normalized counts

# Import my grouped TPM-normalized counts
library(tidyverse)
library(readxl)

dbcantable_gen <- as.data.frame(read_excel('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_grouped_reads/dbcantable_gen.xlsx'))
rownames(dbcantable_gen) <- dbcantable_gen$...1
dbcantable_gen <- dbcantable_gen[, -1]
bpcoldata <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/bp_metadata_01132025.txt")
bpcoldata$name <- gsub("-", ".", bpcoldata$name)
rownames(bpcoldata) <- bpcoldata[,1]
bpcoldata <- bpcoldata[,-(1)]

# Drop any outlier samples you have (based on PCA analysis)
bb = dbcantable_gen
drop <- "-"
bb <- bb[!(row.names(bb) %in% drop), ]
drops <- c("I6_CA14_1_MT_PLANC","I8_CA58_1_MT_PLANC")
bb_filt <- bb[,!(names(bb) %in% drops)]
bb_filt$cat <- rownames(bb_filt)

# import and filter metadata
bpcoldata <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/bp_metadata_01132025.txt")
bpcoldata$name <- gsub("-", ".", bpcoldata$name)
bpcoldata3 <- bpcoldata[bpcoldata$name %in% names(bb_filt),]
meta = bpcoldata3
ff = bb_filt %>%
  pivot_longer(-cat) %>% 
  left_join(meta)

############################################################################################################################################
# One way anova test

result <- aov(value ~ group, data = subset(ff, cat == "GH")) # no CA + S3 + AS vs. CA + AS significance?
summary(result)                                             # but PHA + G1 + AS vs. PHA + AS significance
TukeyHSD(result)

result <- aov(value ~ group, data = subset(ff, cat == "CE")) # CA + S3 + AS vs. CA + AS significance!
summary(result)                                             # no PHA + G1 + AS vs. PHA + AS significance
TukeyHSD(result)

result <- aov(value ~ group, data = subset(ff, cat == "GT")) # CA + S3 + AS vs. CA + AS significance! (? only test I see this?)
summary(result)                                             # no PHA + G1 + AS vs. PHA + AS significance
TukeyHSD(result)

result <- aov(value ~ group, data = subset(ff, cat == "CBM")) # no CA + S3 + AS vs. CA + AS significance
summary(result)                                             # but PHA + G1 + AS vs. PHA + AS significance
TukeyHSD(result)

result <- aov(value ~ group, data = subset(ff, cat == "PL")) # no CA + S3 + AS vs. CA + AS significance
summary(result)                                             # no PHA + G1 + AS vs. PHA + AS significance
TukeyHSD(result)

result <- aov(value ~ group, data = subset(ff, cat == "AA")) # no CA + S3 + AS vs. CA + AS significance
summary(result)                                             # no PHA + G1 + AS vs. PHA + AS significance
TukeyHSD(result)

############################################################################################################################################
# Kruskall/dunn test

kruskal.test(value ~ group, data = subset(ff, cat == "GH"))
library(FSA)
install.packages('FS')
dunnTest(value ~ group, data = subset(ff, cat == "GH"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance
dunnTest(value ~ group, data = subset(ff, cat == "CE"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance
dunnTest(value ~ group, data = subset(ff, cat == "GT"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance
dunnTest(value ~ group, data = subset(ff, cat == "CBM"), method = "bh") # no CA + S3 + AS vs. CA + AS significance BUT PHA + G1 + AS vs. PHA + AS significance
dunnTest(value ~ group, data = subset(ff, cat == "PL"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance
dunnTest(value ~ group, data = subset(ff, cat == "AA"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance

#############################################################################################################################################
# Wilcoxin test

ffsub <- subset(ff, cat == "GH")
pairwise.wilcox.test(ffsub$value, ffsub$group,     # no CA + S3 + AS vs. CA + AS significance?
                     p.adjust.method = "BH")      # but PHA + G1 + AS vs. CA + AS or CA + S3 + AS or PHA + AS is significant (0.019)

ffsub <- subset(ff, cat == "CE")
pairwise.wilcox.test(ffsub$value, ffsub$group,    # CA + S3 + AS vs CA + AS significance!
                     p.adjust.method = "BH")      # no PHA + G1 + AS vs. PHA + AS significance

ffsub <- subset(ff, cat == "GT")
pairwise.wilcox.test(ffsub$value, ffsub$group,    # no CA + S3 + AS vs CA + AS significance
                     p.adjust.method = "BH")      # no PHA + G1 + AS vs. PHA + AS significance

ffsub <- subset(ff, cat == "CBM")
pairwise.wilcox.test(ffsub$value, ffsub$group,    # no CA + S3 + AS vs CA + AS significance
                     p.adjust.method = "BH")      # but PHA + G1 + AS vs. PHA + AS significance

ffsub <- subset(ff, cat == "PL")
pairwise.wilcox.test(ffsub$value, ffsub$group,   # no CA + S3 + AS vs CA + AS significance
                     p.adjust.method = "BH")    # no PHA + G1 + AS vs. PHA + AS significance

ffsub <- subset(ff, cat == "AA") 
pairwise.wilcox.test(ffsub$value, ffsub$group,    # no CA + S3 + AS vs CA + AS significance
                     p.adjust.method = "BH")      # no PHA + G1 + AS vs. PHA + AS significance


