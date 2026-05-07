# C. B. Ahern
# 04/27/2026
# Making a KEGG KO heatmap instead of volcano plot

# import ancombc results
dfkocas3 <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/df_volc_CAS3AS2_KEGGko_ancombc_11052025.tsv")
dfkophag1 <- read_tsv("/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/df_volc_PHAG1AS2_KEGGko_ancombc_11052025.tsv")

dfkocas3sig <- dfkocas3[dfkocas3$padj < 0.05,]
dfkophag1sig <- dfkophag1[dfkophag1$padj < 0.05,]

write_tsv(dfkocas3sig, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/df_volc_CAS3AS2_KEGGko_ancombc_sig_04272026.tsv")
write_tsv(dfkophag1sig, "/Users/colleenahern/Documents/Magda_BPs_experiment/KEGG_mods/df_volc_PHAG1AS2_KEGGko_ancombc_sig_04272026.tsv")


dfkocas3siguniq <- dfkocas3sig[!dfkocas3sig$KO %in% dfkophag1sig$KO,]
dfkophag1siguniq <- dfkophag1sig[!dfkophag1sig$KO %in% dfkocas3sig$KO,]
