# C. Ahern
# 12/13/2024
# Using ALDEX on the BP_film data
# v1 follows code Rodney gave me but idk if I'm doing it right so this is just following a tutorial from online
# https://bioconductor.statistik.tu-dortmund.de/packages/3.6/bioc/vignettes/ALDEx2/inst/doc/ALDEx2_vignette.pdf

# tutorial
library(ALDEx2)
data(selex)
#subset for efficiency
selex <- selex[1201:1600,]
conds <- c(rep("NS", 7), rep("S", 7))

x <- aldex(selex, conds, mc.samples=16, test="t", effect=TRUE, include.sample.summary=FALSE, denom="iqlr", verbose=FALSE)
aldex.plot(x, type="MA", test="welch")
aldex.plot(x, type="MW", test="welch")

head(x)

ef1 <- x[x$effect >= 1,]
head(ef1)


## my data
setwd("/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast")

pswd <- qza_to_phyloseq(features="caqiime2b/bac2analysis/16S-table-nocma-rarefied-10000_filtered.qza",
                        tree="caqiime2b/bac2analysis/rooted-16S-tree-filteredSVs.qza", 
                        taxonomy = "output_bac2/qiime2/input/taxonomy.qza", 
                        metadata= "metadata/metadata_b_12122024.txt")
metadata_b <- read.delim("metadata/metadata_b_12122024.txt")

# For film
pswd_film <- subset_samples(pswd, Experiment == "BP_Film")

df_env <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast/metadata_b_12122024.txt") %>%
  arrange(sample.id) %>% column_to_rownames('sample.id')

dfenvsub <- df_env[df_env$Testcond == "CA film + S3 + AS" | df_env$Testcond == "CA film + AS",]

df <- as.data.frame(pswd_film@otu_table) 

dfsub <- df[,colnames(df) %in% rownames(dfenvsub)]

mdsub <- metadata_b[metadata_b$sample.id %in% colnames(dfsub),]
colnames(dfsub) == mdsub$sample.id

dfsubf <- dfsub[,mdsub$sample.id]
colnames(dfsubf) == mdsub$sample.id

condsy <- mdsub$Testcond




y <- aldex(dfsubf, condsy, mc.samples=200, test="t", effect=TRUE, include.sample.summary=FALSE, denom="iqlr", verbose=FALSE)
aldex.plot(y, type="MA", test="welch")
aldex.plot(y, type="MW", test="welch")

head(y)

efy1 <- y[y$effect >= 1,]
head(efy1)
efy1$asvid <- row.names(efy1)
tt <- as.data.frame(pswd_film@tax_table)
tt$asvid <- row.names(tt)

efy1tax <- merge(efy1, tt, by = 'asvid')

