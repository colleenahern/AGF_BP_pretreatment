# C. Ahern
# 12/13/2024
# Using ALDEX on the BP_film data

install.packages("devtools")
devtools::install_github("ggloor/ALDEx_bioc")

if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ALDEx2")

library(ALDEx2)

if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

# I downloaded the java 64x package onto my mac and then rJava downloaded correctly
install.packages('rJava')
library(rJava)
BiocManager::install("SELEX")

options(java.parameters="-Xmx1500M")
library(SELEX)

workDir = "./cache/"
selex.config(workingDir=workDir, maxThreadNumber=4)

# Extract example data from package, including XML annotation
exampleFiles = selex.exampledata(workDir)
# Load all sample files using XML database
selex.loadAnnotation(exampleFiles[3])

selexdat <- selex.sampleSummary()
r0train = selex.sample(seqName="R0.libraries",
                       + sampleName="R0.barcodeGC", round=0)
r0test = selex.sample(seqName="R0.libraries",
                      + sampleName="R0.barcodeCG", round=0)
r2 = selex.sample(seqName="R2.libraries",
                  + sampleName="ExdHox.R2", round=2)

data(selex)
#subset only the last 400 features for efficiency
selex.sub <- selex[1:400,]

conds <- c(rep("NS", 7), rep("S", 7))
x.all <- aldex(selex.sub, conds, mc.samples=16, test="t", effect=TRUE,
               include.sample.summary=FALSE, denom="all", verbose=FALSE, paired.test=FALSE)

par(mfrow=c(1,2))
aldex.plot(x.all, type="MA", test="welch", xlab="Log-ratio abundance",
           ylab="Difference")
aldex.plot(x.all, type="MW", test="welch", xlab="Dispersion",
           ylab="Difference")

---------------------------------------------------------------------------------------------------
  
  #  https://microbiome.github.io/OMA/containers.html#package-data  tutorial
  
  BiocManager::install("mia")
library(mia)
data(package="mia")
data("GlobalPatterns", package="mia")
GlobalPatterns

install.packages("BiocManager")
install.packages("tidyverse")
install.packages("ggplot2")
BiocManager::install("tidySummarizedExperiment")
library(tidySummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(BiocManager)
library(Biostrings)

# Import dataset
library(mia)
data(mia)
data("Tengeler2020", package = "mia")
tse <- Tengeler2020

# Show patient status by cohort
table(tse$patient_status, tse$cohort) %>%
  knitr::kable()

# get tse from tengelerimport.r
tse <- tse

# Show patient status by cohort
table(tse$patient_status, tse$cohort) %>%
  knitr::kable()

#Preparing the data for DAA: Before starting the analysis, it is recommended to reduce the size and complexity 
# of the data to make the results more reproducible. For this purpose, we agglomerate the features by genus and 
#filter them by a prevalence threshold of 10%.

# Agglomerate by genus and subset by prevalence
tse <- subsetByPrevalentTaxa(tse,
                             rank = "Genus",
                             prevalence = 10 / 100)

install.packages("remotes")
remotes::install_github("tnaake/MatrixQCvis")
library(mia)
# Transform count assay to relative abundances
tse <- transformSamples(tse,
                        assay_name = "counts",
                        #assay.type = "counts",
                        method = "relabundance")  #ignore the warning message about useNames

# Load package
library(ALDEx2)

# Generate Monte Carlo samples of the Dirichlet distribution for each sample.
# Convert each instance using the centered log-ratio transform.
# This is the input for all further analyses.
set.seed(123)
x <- aldex.clr(assay(tse), tse$patient_status)     

# calculates expected values of the Welch's t-test and Wilcoxon rank
# test on the data returned by aldex.clr
x_tt <- aldex.ttest(x, paired.test = FALSE, verbose = FALSE)

# Determines the median clr abundance of the feature in all samples and in
# groups, the median difference between the two groups, the median variation
# within each group and the effect size, which is the median of the ratio
# of the between group difference and the larger of the variance within groups
x_effect <- aldex.effect(x, CI = TRUE, verbose = FALSE)

# combine all outputs 
aldex_out <- data.frame(x_tt, x_effect)

par(mfrow = c(1, 2))

aldex.plot(aldex_out,
           type = "MA",
           test = "welch",
           xlab = "Log-ratio abundance",
           ylab = "Difference",
           cutoff = 0.05)

aldex.plot(aldex_out,
           type = "MW",
           test = "welch",
           xlab = "Dispersion",
           ylab = "Difference",
           cutoff = 0.05)

library(tidyverse)
aldex_out %>%
  rownames_to_column(var = "Genus") %>%
  # here we choose the wilcoxon output rather than t-test output
  filter(wi.eBH <= 0.05)  %>%
  dplyr::select(Genus, we.eBH, wi.eBH, effect, overlap) %>%
  knitr::kable()


# my data
library(qiime2R)
setwd("/Users/colleenahern/Documents/BASF/qiime2downloads")
pswd <- qza_to_phyloseq(features="Analysis2/16S-table-noplant-rarefied-10000_filtered.qza",
                        tree="Analysis2/rooted-16S-tree-filteredSVs.qza", 
                        taxonomy = "Analysis2/16S-rep-seqs-taxonomy.qza", 
                        metadata= "Analysis2/metadata2.txt")

tseme <- mia::makeTreeSummarizedExperimentFromPhyloseq(pswd)

library(ALDEx2)

set.seed(123)
z <- aldex.clr(assay(tseme), tseme$Day)     

# calculates expected values of the Welch's t-test and Wilcoxon rank
# test on the data returned by aldex.clr
z_tt <- aldex.ttest(z, paired.test = FALSE, verbose = FALSE)

# aldex.glm -- tests for multiple conditions while aldex.clr is only good for 2 conditions??
#conds <- tseme$Group
covariates <- as.data.frame(unique(tseme$Group))
colnames(covariates) <- "Group"
mm <- model.matrix(~Group, tseme$Group)
z <- aldex.clr(assay(tseme), mm)     
zglm <- aldex.glm(z)

# Determines the median clr abundance of the feature in all samples and in
# groups, the median difference between the two groups, the median variation
# within each group and the effect size, which is the median of the ratio
# of the between group difference and the larger of the variance within groups
z_effect <- aldex.effect(z, CI = TRUE, verbose = FALSE)

# combine all outputs 
aldex_outme <- data.frame(z_tt, z_effect)

par(mfrow = c(1, 2))

aldex.plot(aldex_outme,
           type = "MA",
           test = "welch",
           xlab = "Log-ratio abundance",
           ylab = "Difference",
           cutoff = 0.05)

aldex.plot(aldex_outme,
           type = "MW",
           test = "welch",
           xlab = "Dispersion",
           ylab = "Difference",
           cutoff = 0.05)

library(tidyverse)
aldex_outme %>%
  rownames_to_column(var = "Genus") %>%
  # here we choose the wilcoxon output rather than t-test output
  filter(wi.eBH <= 0.05)  %>%
  dplyr::select(Genus, we.eBH, wi.eBH, effect, overlap) %>%
  knitr::kable()


molten_tseme <- mia::meltAssay(tseme,
                               add_row_data = TRUE,
                               add_col_data = TRUE,
)
molten_tseme
metadata2 <- read.delim("/Users/colleenahern/Documents/BASF/qiime2downloads/Analysis2/metadata2.txt")



selex <- selex[1201:1600,]
covariates <- data.frame("A" = sample(0:1, 14, replace = TRUE),
                         "B" = c(rep(0, 7), rep(1, 7)))
mm <- model.matrix(~ A + B, covariates)
x <- aldex.clr(selex, mm, mc.samples=4, denom="all")
glm.test <- aldex.glm(x)
glm.eff <- aldex.glm.effect(x)
aldex.plot(glm.test, eff=glm.eff, contrast='B', type='MW', post.hoc='holm')
aldex


-------------------------------------------------------------------------------------------------------------------
  
  # Rodney
  df_env = read_xlsx('simple_study_metadata.xlsx') %>%
  arrange(SampleID) %>% column_to_rownames('SampleID')
df = fread('norm_arg_counts_e4.csv') %>% column_to_rownames('argtype')

df_env = df_env[colnames(df),]

condit_mm = data.frame("StudyID" = as_factor(df_env$StudyID),
                       "major_env" = as_factor(df_env$major_env),
                       "is_plastic" = df_env$is_plastic,
                       "degradation_observed" = df_env$degradation_observed)

rownames(condit_mm) = colnames(df) 
mm = model.matrix(~ StudyID + major_env + is_plastic, condit_mm)


r_df = as.matrix(round(df)) # Rounding data for clr
drx = aldex(r_df,mm,test='glm',verbose = T)
drx.effect = aldex.glm.effect(drx)

# -----------------------------------------------------------------------------------------------------------------------

# Me
# Create phyloseq object from bacteria data
# Okok I forgot that aldex2 is binary so I do need to subset my data.
# Let's subset by Testcond and create 2 subsets: one that has CA film + AS and CA film + S3 + AS, 
# And one that has PHA + AS and PHA + G1 + AS


library(tidySummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(BiocManager)
library(Biostrings)
library(mia)
library(ALDEx2)
library(qiime2R)
library(readxl)
library(phyloseq)

setwd("/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast")

pswd <- qza_to_phyloseq(features="caqiime2b/bac2analysis/16S-table-nocma-rarefied-10000_filtered.qza",
                        tree="caqiime2b/bac2analysis/rooted-16S-tree-filteredSVs.qza", 
                        taxonomy = "output_bac2/qiime2/input/taxonomy.qza", 
                        metadata= "metadata_b_12122024.txt")

# For film
pswd_film <- subset_samples(pswd, Experiment == "BP_Film") ### IS SUBSETTING BEFORE TRANSFORMING FOR NMDS OKAY?????

# df_env = read_xlsx('simple_study_metadata.xlsx') %>%
#   arrange(SampleID) %>% column_to_rownames('SampleID')

df_env <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast/metadata_b_12122024.txt") %>%
  arrange(sample.id) %>% column_to_rownames('sample.id')

dfenvsub <- df_env[df_env$Testcond == "CA film + S3 + AS" | df_env$Testcond == "CA film + AS",]

df <- as.data.frame(pswd_film@otu_table) # %>% column_to_rownames('argtype')

dfsub <- df[,colnames(df) %in% rownames(dfenvsub)]

condit_me <- metadata_b[metadata_b$sample.id %in% colnames(dfsub),]

#rownames(condit_mm) = colnames(df) 
rownames(condit_me) <- condit_me$sample.id 
mm = model.matrix(~Testcond, condit_mm)


#r_df = as.matrix(round(df)) # Rounding data for clr
r_df = as.matrix(round(dfsub)) # Rounding data for clr

#drx = aldex(r_df,mm,test='glm',verbose = T)
#drx = aldex.clr(r_df, mm)   # this one worked
drx = aldex.clr(dfsub, condit_me$Testcond)    # this one worked for subsetting for ttest.  uhh only this worked
# no idea if I'm doing this right


z_tt <- aldex.ttest(drx, paired.test = FALSE, verbose = TRUE)
z_effect <- aldex.effect(drx, CI = TRUE, verbose = TRUE)
aldex_out <- data.frame(z_tt, z_effect)
par(mfrow = c(1, 2))

aldex.plot(aldex_out,
           type = "MA",
           test = "welch",
           xlab = "Log-ratio abundance",
           ylab = "Difference",
           cutoff = 0.05)

aldex.plot(aldex_out,
           type = "MW",
           test = "welch",
           xlab = "Dispersion",
           ylab = "Difference",
           cutoff = 0.05)

aldex_out %>%
  rownames_to_column(var = "Genus") %>%
  # here we choose the wilcoxon output rather than t-test output
  filter(wi.eBH <= 0.05)  %>%
  dplyr::select(Genus, we.eBH, wi.eBH, effect, overlap) %>%
  knitr::kable()

aldex_out1 <- aldex_out %>%
  rownames_to_column(var = "Genus") %>%
  # here we choose the wilcoxon output rather than t-test output
  filter(wi.eBH <= 0.05)  %>%
  dplyr::select(Genus, we.eBH, wi.eBH, effect, overlap)

#-----------------------------------------------------------------------------------------------------------------------

library(tidySummarizedExperiment)
library(tidyverse)
library(ggplot2)
library(BiocManager)
library(Biostrings)
library(mia)
library(ALDEx2)
library(qiime2R)
library(readxl)

pswd <- qza_to_phyloseq(features="Analysis2/16S-table-noplant-rarefied-10000_filtered.qza",
                        tree="Analysis2/rooted-16S-tree-filteredSVs.qza", 
                        taxonomy = "Analysis2/16S-rep-seqs-taxonomy.qza", 
                        metadata= "Analysis2/metadata2.txt")

# df_env = read_xlsx('simple_study_metadata.xlsx') %>%
#   arrange(SampleID) %>% column_to_rownames('SampleID')

df_env <- read.delim("/Users/colleenahern/Documents/BASF/qiime2downloads/Analysis2/metadata2.txt") %>%
  arrange(sample.id) %>% column_to_rownames('sample.id')

dfenvsub1 <- df_env[df_env$Day == 90,]
dfenvsub2 <- dfenvsub1[dfenvsub1$Substrate == "Cellulose" | dfenvsub1$Substrate == "Big Cellulose" | dfenvsub1$Substrate == "Copolymer" | dfenvsub1$Substrate == "HDPE" | dfenvsub1$Substrate == "Blank",]
dfenvsub2$Category <- 0
dfenvsub2$Category[c(1:3,7:9)] <- "Biodegradation"
dfenvsub2$Category[c(4:6,10:15)] <- "No Biodegradation"

df <- as.data.frame(pswd@otu_table) # %>% column_to_rownames('argtype')
dfsub2 = as.data.frame(df[,rownames(dfenvsub2)])

condit_mm = data.frame("Substrate" = as_factor(dfenvsub2$Substrate),
                       "Day" = as_factor(dfenvsub2$Day),
                       "Replicate" = dfenvsub2$Replicate,
                       "Group" = dfenvsub2$Group,
                       "Category" = as_factor(dfenvsub2$Category))

rownames(condit_mm) = colnames(dfsub2) 
mm = model.matrix(~Category, condit_mm)
r_df = as.matrix(round(dfsub2)) # Rounding data for clr

drx = aldex(r_df,mm,test='glm',verbose = T)
drx = aldex.clr(r_df, mm)   # this one worked
# drx = aldex.clr(dfsub2, dfenvsub2$Category)    # this one worked for subsetting for ttest

glm.test <- aldex.glm(drx, mm)
drx.effect = aldex.glm.effect(drx)
drx.effect = as.matrix(drx.effect)
sig <- glm.test[,20]<0.05

glm.effect <- aldex.effect
aldex.plot(glm.effect[["B"]], test="effect", cutoff=2)
sig <- glm.test[,20]<0.05
points(glm.effect[["B"]]$diff.win[sig],
       glm.effect[["B"]]$diff.btw[sig], col="blue")
sig <- glm.test[,20]<0.2
points(glm.effect[["B"]]$diff.win[sig],
       glm.effect[["B"]]$diff.btw[sig], col="blue")

x_tt <- aldex.ttest(glm.test, paired.test = FALSE, verbose = FALSE)




z_tt <- aldex.ttest(drx, paired.test = FALSE, verbose = TRUE)
z_effect <- aldex.effect(drx, CI = TRUE, verbose = TRUE)
aldex_out <- data.frame(z_tt, z_effect)
par(mfrow = c(1, 2))

aldex.plot(aldex_out,
           type = "MA",
           test = "welch",
           xlab = "Log-ratio abundance",
           ylab = "Difference",
           cutoff = 0.05)

aldex.plot(aldex_out,
           type = "MW",
           test = "welch",
           xlab = "Dispersion",
           ylab = "Difference",
           cutoff = 0.05)

aldex_out %>%
  rownames_to_column(var = "Genus") %>%
  # here we choose the wilcoxon output rather than t-test output
  filter(wi.eBH <= 0.05)  %>%
  dplyr::select(Genus, we.eBH, wi.eBH, effect, overlap) %>%
  knitr::kable()








# diversity
pswd_film <- subset_samples(pswd, Experiment == "BP_Film"| Test == "AS original") ### IS SUBSETTING BEFORE TRANSFORMING FOR NMDS OKAY?????

plot_richness(pswd_film, x="Time", measures=c("Shannon", "Simpson"), color="Testcond") + geom_point(size = 4) + theme_bw() + theme(text = element_text(size = 30)) 
