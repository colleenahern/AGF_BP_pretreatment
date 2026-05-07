# 03/20/2025
# C. Ahern
# Using one-way ANOVA, Kruskall-Wallis, and Wilcoxin tests to determine significance of my low-feature (COG, CAZyme) condition_cbaed, 
# TPM-normalized counts
# Break up PHA + G1 + AS t2 and PHA + G1 + AS t1 for my analysis because they seen quite different on the PCA
# and that might be why I'm not seeing any significant CAZYME transcriptional abundance

# Import my condition_cbaed TPM-normalized counts
library(tidyverse)
library(readxl)

dbcantable_gen <- as.data.frame(read_excel('/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/MT_condition_cbaed_reads/dbcantable_gen.xlsx'))
dbcantable_gen <- dbcantable_gen[-1,]
rownames(dbcantable_gen) <- dbcantable_gen$...1
dbcantable_gen <- dbcantable_gen[, -1]
bpcoldata <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/metatranscriptomics/bp_metadata_01132025.txt")

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
bpcoldata3$group2 <- bpcoldata3$condition_cba
bpcoldata3$group2[c(1:2,5:6)] <- "CA + AS"
bpcoldata3$group2[c(3,4,7,8)] <- "CA + S3 + AS"
bpcoldata3$group2[c(9,10,14,15)] <- "PHA + AS"





meta = bpcoldata3
ff = bb_filt %>%
  pivot_longer(-cat) %>% 
  left_join(meta)

############################################################################################################################################
# One way anova test

result <- aov(value ~ condition_cba, data = subset(ff, cat == "GH")) # no CA + S3 + AS vs. CA + AS significance?
summary(result)                                             # no PHA + G1 + AS vs. PHA + AS significance
TukeyHSD(result)

result <- aov(value ~ condition_cba, data = subset(ff, cat == "CE")) # CA + S3 + AS vs. CA + AS significance??? t2t1 tho
summary(result)                                             # no PHA + G1 + AS vs. PHA + AS significance
TukeyHSD(result)

result <- aov(value ~ condition_cba, data = subset(ff, cat == "GT")) # CA + S3 + AS vs. CA + AS significance?? t2t1 tho
summary(result)                                             # no PHA + G1 + AS vs. PHA + AS significance
TukeyHSD(result)

result <- aov(value ~ condition_cba, data = subset(ff, cat == "CBM")) # no CA + S3 + AS vs. CA + AS significance
summary(result)                                             # PHA + G1 + AS vs. PHA + AS significance (but it's lower) PHA + G1 + AS t2-PHA + AS t2 and PHA + G1 + AS t1-PHA + AS t1
TukeyHSD(result)

result <- aov(value ~ condition_cba, data = subset(ff, cat == "PL")) # no CA + S3 + AS vs. CA + AS significance
summary(result)                                             # no PHA + G1 + AS vs. PHA + AS significance
TukeyHSD(result)

result <- aov(value ~ condition_cba, data = subset(ff, cat == "AA")) # no CA + S3 + AS vs. CA + AS significance
summary(result)                                             # no PHA + G1 + AS vs. PHA + AS significance
TukeyHSD(result)


# One way anova test
meta = bpcoldata3
ff = bb_filt %>%
  pivot_longer(-cat) %>% 
  left_join(meta)

result <- aov(value ~ group2, data = subset(ff, cat == "GH")) # no CA + S3 + AS-CA + AS significance
summary(result)                                             # PHA + G1 + AS t1-PHA + AS significance !! no PHA + G1 + AS t2-PHA + AS significance
TukeyHSD(result)

result <- aov(value ~ group2, data = subset(ff, cat == "CE")) # CA + S3 + AS-CA + AS significance
summary(result)                                             # no PHA + G1 + AS t1-PHA + AS or PHA + G1 + AS t2-PHA + AS significance
TukeyHSD(result)

result <- aov(value ~ group2, data = subset(ff, cat == "GT")) # CA + S3 + AS-CA + AS significance 
summary(result)                                             # no PHA + G1 + AS t1-PHA + AS or PHA + G1 + AS t2-PHA + AS significance
TukeyHSD(result)

result <- aov(value ~ group2, data = subset(ff, cat == "CBM")) # no CA + S3 + AS-CA + AS significance
summary(result)                                             # PHA + G1 + AS t1-PHA + AS significance (lower) and PHA + G1 + AS t2-PHA + AS significance 

result <- aov(value ~ group2, data = subset(ff, cat == "PL")) # no CA + S3 + AS vs. CA + AS significance
summary(result)                                             # no PHA + G1 + AS vs. PHA + AS significance
TukeyHSD(result)

result <- aov(value ~ group2, data = subset(ff, cat == "AA")) # no CA + S3 + AS vs. CA + AS significance
summary(result)                                             # no PHA + G1 + AS vs. PHA + AS significance
TukeyHSD(result)



############################################################################################################################################
# Kruskall/dunn test

kruskal.test(value ~ condition_cba, data = subset(ff, cat == "GH"))
library(FSA)
install.packages('FS')
dunnTest(value ~ condition_cba, data = subset(ff, cat == "GH"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance
dunnTest(value ~ condition_cba, data = subset(ff, cat == "CE"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance
dunnTest(value ~ condition_cba, data = subset(ff, cat == "GT"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance
dunnTest(value ~ condition_cba, data = subset(ff, cat == "CBM"), method = "bh") # no CA + S3 + AS vs. CA + AS significance BUT PHA + G1 + AS vs. PHA + AS significance
dunnTest(value ~ condition_cba, data = subset(ff, cat == "PL"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance
dunnTest(value ~ condition_cba, data = subset(ff, cat == "AA"), method = "bh") # no CA + S3 + AS vs. CA + AS significance and no PHA + G1 + AS vs. PHA + AS significance

#############################################################################################################################################
# Wilcoxin test

ffsub <- subset(ff, cat == "GH")
pairwise.wilcox.test(ffsub$value, ffsub$condition_cba,     # no CA + S3 + AS vs. CA + AS significance?
                     p.adjust.method = "BH")      # but PHA + G1 + AS vs. CA + AS or CA + S3 + AS or PHA + AS is significant (0.019)

ffsub <- subset(ff, cat == "CE")
pairwise.wilcox.test(ffsub$value, ffsub$condition_cba,    # CA + S3 + AS vs CA + AS significance!
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


