# C. Ahern
# 12/05/2025
# getting annotations for the KOs in my volcano plots

cas3ko <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/CAS3AS2_KEGGko_sig_labeled_11102025.tsv')
phag1ko <- read_tsv('/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/PHAG1AS2_KEGGko_sig_labeled_11102025.tsv')

cas3kovolc <- cas3ko[cas3ko$KOid %in% c("K14069","K14068","K22162","K19793","K18140","K06079","K18129","K21835","K06974","K16011"),]
phag1kovolc <- phag1ko[phag1ko$KOid %in% c("K04072","K06399","K02486","K07480","K07217","K10954","K19824","K12880","K06012","K07705"),]

cas3kovolc <- cas3kovolc[, c(1,2,4,7,3)]
phag1kovolc <- phag1kovolc[, c(1,2,4,7,3)]
