# C. B. Ahern
# 07/17/2025
# Using the results from /Users/colleenahern/Documents/Magda_BPs_experiment/Rcodes/BP_MT_transabund_grouped_KEGGkos_05062025.R
# Starting with a new script just to make organization cleaner

# What I want to do: merge the results from BP_MT_transabund_grouped_KEGGkos_05062025.R with the pathways correlated to the KEGG kos
# and add the KO and pathway annotations so I can look at what's being upregulated
library(readr)

CAS3AS_KEGGko_unique <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAS3AS_KEGGko_unique_ancombc_05212025.tsv")
PHAG1AS_KEGGko_unique <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/PHAG1AS_KEGGko_unique_ancombc_05212025.tsv")
kegmod_paths <- read.delim("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/keggmods_paths_06022025.txt")
kegmod_paths2 <- kegmod_paths[!grepl("path:ko", kegmod_paths$Pathway),]
