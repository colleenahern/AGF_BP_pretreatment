# 12/16/2024
# Colleen Ahern
# Trying ANCOM-BC on my data for differential abundance testing
# There are a few ASVs with effect size >1 for aldex2 but the p values are all >0.05 so I am not sure 
# what to do
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")
BiocManager::install("ANCOMBC")

knitr::opts_chunk$set(message = FALSE, warning = FALSE, comment = NA, 
                      fig.width = 6.25, fig.height = 5)
library(ANCOMBC)
library(tidyverse)
library(DT)
options(DT.options = list(
  initComplete = JS("function(settings, json) {",
                    "$(this.api().table().header()).css({'background-color': 
  '#000', 'color': '#fff'});","}")))
library(dplyr)
library(ggplot2)
library(phyloseq)
library(qiime2R)

################################################################################################################
### This is for BACTERIA data (I rarefied 10,000 cutoff)

# Bacteria: CA film + S3 + AS and CA film + AS

## load data
setwd("/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast")

pswd <- qza_to_phyloseq(features="caqiime2b/bac2analysis/16S-table-nocma-rarefied-10000_filtered.qza",
                        tree="caqiime2b/bac2analysis/rooted-16S-tree-filteredSVs.qza", 
                        taxonomy = "output_bac2/qiime2/input/taxonomy.qza", 
                        metadata= "metadata/metadata_b_12122024.txt")

sample_variables(pswd)
metadata_b <- read.delim("metadata/metadata_b_12122024.txt")

# Subset for film
pswd_film <- subset_samples(pswd, Experiment == "BP_Film") # subset only film experiment
print(pswd_film)

# Subset for CA film + S3 + AS and CA film + AS conditions only (all timepoints because replicates are limited)
pswd_filmsub <- subset_samples(pswd, Testcond == "CA film + S3 + AS" | Testcond == "CA film + AS")
print(pswd_filmsub)
unique(pswd_filmsub@sam_data$Testcond)

# Run ancombc function using the phyloseq object
# Since I have small sample sizes, I set "conserve = TRUE" and lowered the p value cutoff for significance 
# from 0.05 to 0.01
# But I can loosen these restraints if we want
out = ancombc(data = pswd_filmsub, tax_level = "Genus", 
              formula = "Testcond", 
              p_adj_method = "holm",
              conserve = TRUE,
              alpha = 0.01,
              verbose = TRUE)

res = out$res
res_global = out$res_global

## ANCOMBC primary result
# LFC
tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "CA film + S3 + AS")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# SE
tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Test statistic
tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# P-values
tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Adjusted p-values
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Differentially abundant taxa
tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")

cadiff <- tab_diff[tab_diff$`CA film + S3 + AS` == TRUE,]
# ok results??

tab_lfc_filt <- tab_lfc[tab_lfc$Taxon %in% cadiff$Taxon,]
colnames(tab_lfc_filt) <- c("Taxon", "Intercept LFC", "CA film + S3 + AS LFC")

tab_se_filt <- tab_se[tab_se$Taxon %in% cadiff$Taxon,]
colnames(tab_se_filt) <- c("Taxon", "Intercept SE", "CA film + S3 + AS SE")

tab_w_filt <- tab_w[tab_w$Taxon %in% cadiff$Taxon,]
colnames(tab_w_filt) <- c("Taxon", "Intercept W", "CA film + S3 + AS W")

tab_p_filt <- tab_p[tab_p$Taxon %in% cadiff$Taxon,]
colnames(tab_p_filt) <- c("Taxon", "Intercept P", "CA film + S3 + AS P")

tab_q_filt <- tab_q[tab_q$Taxon %in% cadiff$Taxon,]
colnames(tab_q_filt) <- c("Taxon", "Intercept Q", "CA film + S3 + AS Q")

tab_diff_filt <- tab_diff[tab_diff$Taxon %in% cadiff$Taxon,]
colnames(tab_diff_filt) <- c("Taxon", "Intercept diff", "CA film + S3 + AS diff")

#put all data frames into list
df_list <- list(tab_lfc_filt, tab_se_filt, tab_w_filt, tab_p_filt, tab_q_filt, tab_diff_filt)

#merge all data frames in list
res_filt <- df_list %>% reduce(full_join, by='Taxon')

# save list
writexl::write_xlsx(res_filt, "/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast/ANCOMBC/CAfilm_S3_AS_bacteria_ancombc.xlsx")




### Visualization of differentially abundant taxa heatmap
tab_diff = res$diff_abn
tab_diff %>% datatable(caption = "Differentially Abundant Taxa 
                       from the Global Test Result")

sig_taxa = tab_diff %>%
  dplyr::filter(`TestcondCA film + S3 + AS` == TRUE) %>%
  .$taxon

tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "CA film + S3 + AS vs. CA film + AS")
colnames(tab_lfc) = col_name
tab_lfc <- tab_lfc[,-2]

df_CAS3AS = tab_lfc %>%
  filter(Taxon %in% sig_taxa)

df_heat = df_CAS3AS %>%
  pivot_longer(cols = -one_of("Taxon"),
               names_to = "region", values_to = "value") %>%
  mutate(value = round(value, 2))
df_heat$Taxon = factor(df_heat$Taxon, levels = sort(sig_taxa))

lo = floor(min(df_heat$value))
up = ceiling(max(df_heat$value))
mid = (lo + up)/2
p_heat = df_heat %>%
  ggplot(aes(x = region, y = Taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(region, Taxon, label = value), color = "black", size = 7) +
  labs(x = NULL, y = NULL, size = 7,
       title = "Log fold changes for globally significant taxa") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))
p_heat

##################################################################################################################
# 12/17/2024
# Now do this for Bacteria: PHA film + AS and PHA film + G1 + AS

## load data
setwd("/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast")

pswd <- qza_to_phyloseq(features="caqiime2b/bac2analysis/16S-table-nocma-rarefied-10000_filtered.qza",
                        tree="caqiime2b/bac2analysis/rooted-16S-tree-filteredSVs.qza", 
                        taxonomy = "output_bac2/qiime2/input/taxonomy.qza", 
                        metadata= "metadata/metadata_b_12122024.txt")

sample_variables(pswd)
metadata_b <- read.delim("metadata/metadata_b_12122024.txt")

# Subset for film
pswd_film <- subset_samples(pswd, Experiment == "BP_Film") # subset only film experiment
print(pswd_filmsub)

# Subset for PHA film + G1 + AS and PHA film + AS conditions only (all timepoints because replicates are limited)
pswd_filmsub <- subset_samples(pswd, Testcond == "PHA film + G1 + AS" | Testcond == "PHA film + AS")
print(pswd_filmsub)
unique(pswd_filmsub@sam_data$Testcond)

# Run ancombc function using the phyloseq object
# Since I have small sample sizes, I set "conserve = TRUE" and lowered the p value cutoff for significance 
# from 0.05 to 0.01
# But I can loosen these restraints if we want
out = ancombc(data = pswd_filmsub, tax_level = "Genus", 
              formula = "Testcond", 
              p_adj_method = "holm",
              conserve = TRUE,
              alpha = 0.01,
              verbose = TRUE)

res = out$res
res_global = out$res_global

## ANCOMBC primary result
# LFC
tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "PHA film + G1 + AS")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# SE
tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Test statistic
tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# P-values
tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Adjusted p-values
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Differentially abundant taxa
tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")

cadiff <- tab_diff[tab_diff$`PHA film + G1 + AS` == TRUE,]

tab_lfc_filt <- tab_lfc[tab_lfc$Taxon %in% cadiff$Taxon,]
colnames(tab_lfc_filt) <- c("Taxon", "Intercept LFC", "PHA film + G1 + AS LFC")

tab_se_filt <- tab_se[tab_se$Taxon %in% cadiff$Taxon,]
colnames(tab_se_filt) <- c("Taxon", "Intercept SE", "PHA film + G1 + AS SE")

tab_w_filt <- tab_w[tab_w$Taxon %in% cadiff$Taxon,]
colnames(tab_w_filt) <- c("Taxon", "Intercept W", "PHA film + G1 + AS W")

tab_p_filt <- tab_p[tab_p$Taxon %in% cadiff$Taxon,]
colnames(tab_p_filt) <- c("Taxon", "Intercept P", "PHA film + G1 + AS P")

tab_q_filt <- tab_q[tab_q$Taxon %in% cadiff$Taxon,]
colnames(tab_q_filt) <- c("Taxon", "Intercept Q", "PHA film + G1 + AS Q")

tab_diff_filt <- tab_diff[tab_diff$Taxon %in% cadiff$Taxon,]
colnames(tab_diff_filt) <- c("Taxon", "Intercept diff", "PHA film + G1 + AS diff")

#put all data frames into list
df_list <- list(tab_lfc_filt, tab_se_filt, tab_w_filt, tab_p_filt, tab_q_filt, tab_diff_filt)

#merge all data frames in list
res_filt <- df_list %>% reduce(full_join, by='Taxon')

# save list
writexl::write_xlsx(res_filt, "/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast/ANCOMBC/PHAfilm_G1_AS_bacteria_ancombc.xlsx")




### Visualization of differentially abundant taxa heatmap
tab_diff = res$diff_abn
tab_diff %>% datatable(caption = "Differentially Abundant Taxa 
                       from the Global Test Result")

sig_taxa = tab_diff %>%
  dplyr::filter(`TestcondPHA film + G1 + AS` == TRUE) %>%
  .$taxon

tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "PHA film + G1 + AS vs. PHA film + AS")
colnames(tab_lfc) = col_name
tab_lfc <- tab_lfc[,-2]

df_PHAG1AS = tab_lfc %>%
  filter(Taxon %in% sig_taxa)

df_heat = df_PHAG1AS %>%
  pivot_longer(cols = -one_of("Taxon"),
               names_to = "region", values_to = "value") %>%
  mutate(value = round(value, 2))
df_heat$Taxon = factor(df_heat$Taxon, levels = sort(sig_taxa))

lo = floor(min(df_heat$value))
up = ceiling(max(df_heat$value))
mid = (lo + up)/2
p_heat = df_heat %>%
  ggplot(aes(x = region, y = Taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(region, Taxon, label = value), color = "black", size = 7) +
  labs(x = NULL, y = NULL, 
       title = "Log fold changes for globally significant taxa") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))
p_heat


##################################################################################################################
# 12/17/2024
# Now do this for Archaea: CA film + S3 + AS and CA film + AS

## load data
setwd("/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast")

pswd <- qza_to_phyloseq(features="caqiime2b/arc2analysis/16S-table-nocmb-rarefied-3000_filtered.qza",
                        tree="caqiime2b/arc2analysis/rooted-16S-tree-filteredSVs.qza", 
                        taxonomy = "output_arc2/qiime2/input/taxonomy.qza", 
                        metadata= "metadata/metadata_a_12122024.txt")

sample_variables(pswd)
metadata_a <- read.delim("metadata/metadata_a_12122024.txt")

# Subset for film
pswd_film <- subset_samples(pswd, Experiment == "BP_Film") # subset only film experiment
print(pswd_filmsub)

# Subset for CA film + S3 + AS and CA film + AS conditions only (all timepoints because replicates are limited)
pswd_filmsub <- subset_samples(pswd, Testcond == "CA film + S3 + AS" | Testcond == "CA film + AS")
print(pswd_filmsub)
unique(pswd_filmsub@sam_data$Testcond)


# Run ancombc function using the phyloseq object
# Since I have small sample sizes, I set "conserve = TRUE" and lowered the p value cutoff for significance 
# from 0.05 to 0.01
# But I can loosen these restraints if we want
out = ancombc(data = pswd_filmsub, tax_level = "Genus", 
              formula = "Testcond", 
              p_adj_method = "holm",
              conserve = TRUE,
              alpha = 0.01,
              verbose = TRUE)

res = out$res
res_global = out$res_global


## ANCOMBC primary result
# LFC
tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "CA film + S3 + AS")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# SE
tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Test statistic
tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# P-values
tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Adjusted p-values
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Differentially abundant taxa
tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")

cadiff <- tab_diff[tab_diff$`CA film + S3 + AS` == TRUE,]

tab_lfc_filt <- tab_lfc[tab_lfc$Taxon %in% cadiff$Taxon,]
colnames(tab_lfc_filt) <- c("Taxon", "Intercept LFC", "CA film + S3 + AS LFC")

tab_se_filt <- tab_se[tab_se$Taxon %in% cadiff$Taxon,]
colnames(tab_se_filt) <- c("Taxon", "Intercept SE", "CA film + S3 + AS SE")

tab_w_filt <- tab_w[tab_w$Taxon %in% cadiff$Taxon,]
colnames(tab_w_filt) <- c("Taxon", "Intercept W", "CA film + S3 + AS W")

tab_p_filt <- tab_p[tab_p$Taxon %in% cadiff$Taxon,]
colnames(tab_p_filt) <- c("Taxon", "Intercept P", "CA film + S3 + AS P")

tab_q_filt <- tab_q[tab_q$Taxon %in% cadiff$Taxon,]
colnames(tab_q_filt) <- c("Taxon", "Intercept Q", "CA film + S3 + AS Q")

tab_diff_filt <- tab_diff[tab_diff$Taxon %in% cadiff$Taxon,]
colnames(tab_diff_filt) <- c("Taxon", "Intercept diff", "CA film + S3 + AS diff")

#put all data frames into list
df_list <- list(tab_lfc_filt, tab_se_filt, tab_w_filt, tab_p_filt, tab_q_filt, tab_diff_filt)

#merge all data frames in list
res_filt <- df_list %>% reduce(full_join, by='Taxon')

# save list
writexl::write_xlsx(res_filt, "/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast/ANCOMBC/CAfilm_S3_AS_archaea_ancombc.xlsx")




### Visualization of differentially abundant taxa heatmap
tab_diff = res$diff_abn
tab_diff %>% datatable(caption = "Differentially Abundant Taxa 
                       from the Global Test Result")

sig_taxa = tab_diff %>%
  dplyr::filter(`TestcondCA film + S3 + AS` == TRUE) %>%
  .$taxon

tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "CA film + S3 + AS vs. CA film + AS")
colnames(tab_lfc) = col_name
tab_lfc <- tab_lfc[,-2]

df_CAS3AS = tab_lfc %>%
  filter(Taxon %in% sig_taxa)

df_heat = df_CAS3AS %>%
  pivot_longer(cols = -one_of("Taxon"),
               names_to = "region", values_to = "value") %>%
  mutate(value = round(value, 2))
df_heat$Taxon = factor(df_heat$Taxon, levels = sort(sig_taxa))

lo = floor(min(df_heat$value))
up = ceiling(max(df_heat$value))
mid = (lo + up)/2
p_heat = df_heat %>%
  ggplot(aes(x = region, y = Taxon, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(region, Taxon, label = value), color = "black", size = 7) +
  labs(x = NULL, y = NULL, 
       title = "Log fold changes for globally significant taxa") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20))
p_heat


##################################################################################################################
# 12/17/2024
# Now do this for Archaea: PHA film + AS and PHA film + G1 + AS

## load data
setwd("/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast")

pswd <- qza_to_phyloseq(features="caqiime2b/arc2analysis/16S-table-nocmb-rarefied-3000_filtered.qza",
                        tree="caqiime2b/arc2analysis/rooted-16S-tree-filteredSVs.qza", 
                        taxonomy = "output_arc2/qiime2/input/taxonomy.qza", 
                        metadata= "metadata/metadata_a_12122024.txt")

sample_variables(pswd)
metadata_a <- read.delim("metadata/metadata_a_12122024.txt")

# Subset for film
pswd_film <- subset_samples(pswd, Experiment == "BP_Film") # subset only film experiment
print(pswd_filmsub)

# Subset for PHA film + G1 + AS and PHA film + AS conditions only (all timepoints because replicates are limited)
pswd_filmsub <- subset_samples(pswd, Testcond == "PHA film + G1 + AS" | Testcond == "PHA film + AS")
print(pswd_filmsub)
unique(pswd_filmsub@sam_data$Testcond)


# Run ancombc function using the phyloseq object
# Since I have small sample sizes, I set "conserve = TRUE" and lowered the p value cutoff for significance 
# from 0.05 to 0.01
# But I can loosen these restraints if we want
out = ancombc(data = pswd_filmsub, tax_level = "Genus", 
              formula = "Testcond", 
              p_adj_method = "holm",
              conserve = TRUE,
              alpha = 0.01,
              verbose = TRUE)

## NOTE: a lot of warnings produced. There is only one sample in the PHA film + AS condition so it is 
# likely unreliable
# Were any warnings produced in the above 3 runs?
# #1 no
# #2 no
# #3 no
# so warnings were only produced for this comparison

res = out$res
res_global = out$res_global


## ANCOMBC primary result
# LFC
tab_lfc = res$lfc
col_name = c("Taxon", "Intercept", "PHA film + G1 + AS")
colnames(tab_lfc) = col_name
tab_lfc %>% 
  datatable(caption = "Log Fold Changes from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# SE
tab_se = res$se
colnames(tab_se) = col_name
tab_se %>% 
  datatable(caption = "SEs from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Test statistic
tab_w = res$W
colnames(tab_w) = col_name
tab_w %>% 
  datatable(caption = "Test Statistics from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# P-values
tab_p = res$p_val
colnames(tab_p) = col_name
tab_p %>% 
  datatable(caption = "P-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Adjusted p-values
tab_q = res$q
colnames(tab_q) = col_name
tab_q %>% 
  datatable(caption = "Adjusted p-values from the Primary Result") %>%
  formatRound(col_name[-1], digits = 2)

# Differentially abundant taxa
tab_diff = res$diff_abn
colnames(tab_diff) = col_name
tab_diff %>% 
  datatable(caption = "Differentially Abundant Taxa from the Primary Result")

cadiff <- tab_diff[tab_diff$`PHA film + G1 + AS` == TRUE,]

tab_lfc_filt <- tab_lfc[tab_lfc$Taxon %in% cadiff$Taxon,]
colnames(tab_lfc_filt) <- c("Taxon", "Intercept LFC", "PHA film + G1 + AS LFC")

tab_se_filt <- tab_se[tab_se$Taxon %in% cadiff$Taxon,]
colnames(tab_se_filt) <- c("Taxon", "Intercept SE", "PHA film + G1 + AS SE")

tab_w_filt <- tab_w[tab_w$Taxon %in% cadiff$Taxon,]
colnames(tab_w_filt) <- c("Taxon", "Intercept W", "PHA film + G1 + AS W")

tab_p_filt <- tab_p[tab_p$Taxon %in% cadiff$Taxon,]
colnames(tab_p_filt) <- c("Taxon", "Intercept P", "PHA film + G1 + AS P")

tab_q_filt <- tab_q[tab_q$Taxon %in% cadiff$Taxon,]
colnames(tab_q_filt) <- c("Taxon", "Intercept Q", "PHA film + G1 + AS Q")

tab_diff_filt <- tab_diff[tab_diff$Taxon %in% cadiff$Taxon,]
colnames(tab_diff_filt) <- c("Taxon", "Intercept diff", "PHA film + G1 + AS diff")

#put all data frames into list
df_list <- list(tab_lfc_filt, tab_se_filt, tab_w_filt, tab_p_filt, tab_q_filt, tab_diff_filt)

#merge all data frames in list
res_filt <- df_list %>% reduce(full_join, by='Taxon')

# save list
writexl::write_xlsx(res_filt, "/Users/colleenahern/Documents/Magda_BPs_experiment/ampliseq_fungiplast/ANCOMBC/PHAfilm_G1_AS_archaea_ancombc.xlsx")





























# Bias-corrected abundances
samp_frac = out$samp_frac
# Replace NA with 0
samp_frac[is.na(samp_frac)] = 0 
# Add pesudo-count (1) to avoid taking the log of 0
log_obs_abn = log(out$feature_table + 1)
# Adjust the log observed abundances
log_corr_abn = t(t(log_obs_abn) - samp_frac)
# Show the first 6 samples
round(log_corr_abn[, 1:6], 2) %>% 
  datatable(caption = "Bias-corrected log observed abundances")

# Visualization for age
df_lfc = data.frame(res$lfc[, -1] * res$diff_abn[, -1], check.names = FALSE) %>%
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
df_se = data.frame(res$se[, -1] * res$diff_abn[, -1], check.names = FALSE) %>% 
  mutate(taxon_id = res$diff_abn$taxon) %>%
  dplyr::select(taxon_id, everything())
colnames(df_se)[-1] = paste0(colnames(df_se)[-1], "SE")

df_fig_age = df_lfc %>% 
  dplyr::left_join(df_se, by = "taxon_id") %>%
  dplyr::transmute(taxon_id, age, ageSE) %>%
  dplyr::filter(age != 0) %>% 
  dplyr::arrange(desc(age)) %>%
  dplyr::mutate(direct = ifelse(age > 0, "Positive LFC", "Negative LFC"))
df_fig_age$taxon_id = factor(df_fig_age$taxon_id, levels = df_fig_age$taxon_id)
df_fig_age$direct = factor(df_fig_age$direct, 
                           levels = c("Positive LFC", "Negative LFC"))

p_age = ggplot(data = df_fig_age, 
               aes(x = taxon_id, y = age, fill = direct, color = direct)) + 
  geom_bar(stat = "identity", width = 0.7, 
           position = position_dodge(width = 0.4)) +
  geom_errorbar(aes(ymin = age - ageSE, ymax = age + ageSE), width = 0.2,
                position = position_dodge(0.05), color = "black") + 
  labs(x = NULL, y = "Log fold change", 
       title = "Log fold changes as one unit increase of age") + 
  scale_fill_discrete(name = NULL) +
  scale_color_discrete(name = NULL) +
  theme_bw() + 
  theme(plot.title = element_text(hjust = 0.5),
        panel.grid.minor.y = element_blank(),
        axis.text.x = element_text(angle = 60, hjust = 1))
p_age

# Visualization for BMI
df_fig_bmi = df_lfc %>% 
  filter(bmioverweight != 0 | bmilean != 0) %>%
  transmute(taxon_id, 
            `Overweight vs. Obese` = round(bmioverweight, 2),
            `Lean vs. Obese` = round(bmilean, 2)) %>%
  pivot_longer(cols = `Overweight vs. Obese`:`Lean vs. Obese`, 
               names_to = "group", values_to = "value") %>%
  arrange(taxon_id)
lo = floor(min(df_fig_bmi$value))
up = ceiling(max(df_fig_bmi$value))
mid = (lo + up)/2
p_bmi = df_fig_bmi %>%
  ggplot(aes(x = group, y = taxon_id, fill = value)) + 
  geom_tile(color = "black") +
  scale_fill_gradient2(low = "blue", high = "red", mid = "white", 
                       na.value = "white", midpoint = mid, limit = c(lo, up),
                       name = NULL) +
  geom_text(aes(group, taxon_id, label = value), color = "black", size = 4) +
  labs(x = NULL, y = NULL, title = "Log fold changes as compared to obese subjects") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))
p_bmi

# ANCOMBC global test result
tab_w = res_global[, c("taxon", "W")]
tab_w %>% datatable(caption = "Test Statistics 
                    from the Global Test Result") %>%
  formatRound(c("W"), digits = 2)

# P-values
tab_p = res_global[, c("taxon", "p_val")]
tab_p %>% datatable(caption = "P-values 
                    from the Global Test Result") %>%
  formatRound(c("p_val"), digits = 2)

# Adjusted p-values
tab_q = res_global[, c("taxon", "q_val")]
tab_q %>% datatable(caption = "Adjusted p-values 
                    from the Global Test Result") %>%
  formatRound(c("q_val"), digits = 2)




















# Run ancombc function using the tse object
tse = mia::makeTreeSummarizedExperimentFromPhyloseq(pseq)

out = ancombc(data = tse, assay_name = "counts", tax_level = "Family", 
              formula = "age + region + bmi", 
              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
              group = "bmi", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
              n_cl = 1, verbose = TRUE)

res = out$res
res_global = out$res_global

# Run ancombc function by directly providing the abundance and metadata
abundance_data = microbiome::abundances(pseq)
aggregate_data = microbiome::abundances(microbiome::aggregate_taxa(pseq, "Family"))
meta_data = microbiome::meta(pseq)

out = ancombc(data = abundance_data, aggregate_data = aggregate_data, 
              meta_data = meta_data, formula = "age + region + bmi", 
              p_adj_method = "holm", prv_cut = 0.10, lib_cut = 1000, 
              group = "bmi", struc_zero = TRUE, neg_lb = TRUE, tol = 1e-5, 
              max_iter = 100, conserve = TRUE, alpha = 0.05, global = TRUE,
              n_cl = 1, verbose = TRUE)

res = out$res
res_global = out$res_global