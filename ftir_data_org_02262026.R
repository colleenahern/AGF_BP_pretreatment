# C. Ahern
# 05/26/2026
# Making graph of Magda's FTIR data for supplemental manuscript

# Source - https://stackoverflow.com/a/74254982
# Posted by Rodrigo Zepeda, modified by community. See post 'Timeline' for change history
# Retrieved 2026-02-26, License - CC BY-SA 4.0

all_files <- list.files(path = "/Users/colleenahern/Documents/Magda_BPs_experiment/FTIR_data", pattern = "*.CSV", full.names = TRUE)
all_files


setwd("~")

for (filepath in all_files) {
  # Generate variable name from filename (e.g., "file1")
  var_name <- sub(pattern = ".CSV", replacement = "", filepath)
  var_name <- sub(pattern = "/Users/colleenahern/Documents/Magda_BPs_experiment/FTIR_data/", replacement = "", var_name)
  
  # Read the file
  data <- read_csv(filepath, col_names = FALSE)
  
  # Assign the data to a variable name in the global environment
  assign(var_name, data, envir = .GlobalEnv)
}

all_varnames <- gsub(".CSV|/Users/colleenahern/Documents/Magda_BPs_experiment/FTIR_data/", "", all_files)
all_varnames <- as.data.frame(all_varnames)
all_varnames



for (f in 1:nrow(all_varnames)) {
  varname <- as.character(all_varnames[f, ])
  obj <- get(varname)
  names(obj) <- c(
    paste0(varname, "_wl"),
    paste0(varname, "_tm")
  )
  assign(varname, obj)
}

all_varnames <- as.data.frame(all_varnames)

matching_dfs <- all_varnames[grepl("PHA|58-60days-rep-3", all_varnames$all_varnames),]
matching_dfs
ca_ftir <- do.call(
  cbind,
  lapply(matching_dfs, get)
)

ca_ftir[,grepl("wl", names(ca_ftir))]
ca_ftir$AllEqual <- apply(ca_ftir[,grepl("wl", names(ca_ftir))], 1, function(row) all(row == row[1]))
all(ca_ftir$AllEqual) == TRUE
ca_ftir2 <- ca_ftir[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32)]

matching_dfs <- all_varnames[grepl("PLA|pla", all_varnames$all_varnames),]
matching_dfs
pla_ftir <- do.call(
  cbind,
  lapply(matching_dfs, get)
)
pla_ftir$AllEqual <- apply(pla_ftir[,grepl("wl", names(pla_ftir))], 1, function(row) all(row == row[1]))
all(pla_ftir$AllEqual) == TRUE
pla_ftir2 <- pla_ftir[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36)]

matching_dfs <- all_varnames[grepl("PHA|pha", all_varnames$all_varnames),]
matching_dfs
pha_ftir <- do.call(
  cbind,
  lapply(matching_dfs, get)
)
pha_ftir$AllEqual <- apply(pha_ftir[,grepl("wl", names(pha_ftir))], 1, function(row) all(row == row[1]))
all(pha_ftir$AllEqual) == TRUE
pha_ftir2 <- pha_ftir[,c(1,2,4,6,8,10,12,14,16,18,20,22,24,26,28,30,32,34,36,38,40)]

# plot CA FTIR
ca_ftir3 <- ca_ftir2
ca_ftir3$CA_58_60days_avg_tm <- rowMeans(ca_ftir3[, c("CA-58-60days-rep1_tm", "CA-58-60days-rep2_tm", "58-60days-rep-3_tm")])
ca_ftir3$CA_60A_4days_AF_avg_tm <- rowMeans(ca_ftir3[, c("CA_60A-4days-AF-rep1_tm", "CA_60A-4days-AF-rep2_tm", "CA_60A-4days-AF-rep3_tm")])
ca_ftir3$CA_61_60days_avg_tm <- rowMeans(ca_ftir3[, c("CA_61-60days-rep1_tm", "CA_61-60days-rep2_tm")])
ca_ftir3$CA_77_60days_avg_tm <- rowMeans(ca_ftir3[, c("CA_77-60days-rep1_tm", "CA_77-60days-rep2_tm")])
ca_ftir3$CA_65A_4days_avg_tm <- rowMeans(ca_ftir3[, c("CA-65A-4-days-rep1_tm", "CA-65A-4-days-rep2_tm")])
ca_ftir3$CA_control_avg_tm <- rowMeans(ca_ftir3[, c("CA-control_rep1_tm", "CA-control_rep2_tm", "CA-control_rep3_tm")])
ca_ftir4 <- ca_ftir3[,c(1,10,18:23)]
ca_ftir5 <- ca_ftir4
names(ca_ftir5) <- gsub("58-60days-rep-3_wl", "wl", names(ca_ftir5))
names(ca_ftir5) <- gsub("CA-12-60days", "CA_AS_Day60", names(ca_ftir5))
names(ca_ftir5) <- gsub("CA_58_60days", "CA_S3_AS_Day60_a", names(ca_ftir5))
names(ca_ftir5) <- gsub("CA_65A_4days", "tCA_S3_Day4", names(ca_ftir5))
names(ca_ftir5) <- gsub("CA_60A_4days_AF", "CA_S3_Day4", names(ca_ftir5))
names(ca_ftir5) <- gsub("CA_61_60days", "CA_S3_AS_Day60_b", names(ca_ftir5))
names(ca_ftir5) <- gsub("CA_77_60days", "tCA_AS_Day60", names(ca_ftir5))
# no tCA_S3_AS_Day60? ask magda

cols <- c("CA control"="#0072B2","CA + AS Day 60"="#E69F00","tCA + AS Day 60"="#56B4E9","CA + S3 + AS Day 60 (a)"="#D55E00", "CA + S3 + AS Day 60 (b)"="#D55E00","tCA + S3 Day 4"="#882255","CA + S3 Day 4"="#332288")
caft <- ggplot(ca_ftir5, aes(x = wl)) + geom_line(aes(y = CA_control_avg_tm, group = 1, color = "CA control")) + 
  geom_line(aes(y = CA_AS_Day60_tm, group = 1, color = "CA + AS Day 60")) + 
  geom_line(aes(y = tCA_AS_Day60_avg_tm, group = 1, color = "tCA + AS Day 60")) + 
  geom_line(aes(y = CA_S3_AS_Day60_a_avg_tm, group = 1, color = "CA + S3 + AS Day 60 (a)")) + 
  geom_line(aes(y = CA_S3_AS_Day60_b_avg_tm, group = 1, color = "CA + S3 + AS Day 60 (b)")) + 
  geom_line(aes(y = CA_S3_Day4_avg_tm, group = 1, color = "CA + S3 Day 4")) + 
  geom_line(aes(y = tCA_S3_Day4_avg_tm, group = 1, color = "tCA + S3 Day 4")) + 
  labs(
    #title = "CA FTIR",
    x = expression("Wavenumber (cm"^-1*")"),
    y = "Transmittance (%)") +
  ylim(0, 100) +
  xlim(3999,525) +
  # scale_x_reverse() +
  # scale_x_continuous(breaks = c(4,32,60), expand = c(0, 0)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face = "bold", vjust = -2))
caft

# plot PHA FTIR
pha_ftir3 <- pha_ftir2
pha_ftir3$PHA_45A_20days_avg_tm <- rowMeans(pha_ftir3[, c("PHA-45A-20-days-rep1_tm", "PHA-45A-20-days-rep2_tm")])
pha_ftir3$PHA_45A_4days_avg_tm <- rowMeans(pha_ftir3[, c("PHA-45A-4-days-rep1_tm", "PHA-45A-4-days-rep2_tm")])
pha_ftir3$PHA_47_60days_AF_avg_tm <- rowMeans(pha_ftir3[, c("PHA-46-60-days-rep1_tm", "pha-47-60days-rep2_tm", "pha-47-60days-rep3_tm")])
pha_ftir3$PHA_52_60days_avg_tm <- rowMeans(pha_ftir3[, c("PHA-52-60-days-rep1_tm", "PHA-52-60-days-rep2_tm")])
pha_ftir3$PHA_53A_20days_avg_tm <- rowMeans(pha_ftir3[, c("PHA-53A-20-days-rep1_tm", "PHA-53A-20-days-rep2_tm")])
pha_ftir3$PHA_53A_4days_avg_tm <- rowMeans(pha_ftir3[, c("PHA-53A-4-days-rep1_tm", "PHA-53A-4-days-rep2_tm")])
pha_ftir3$PHA_74_60days_avg_tm <- rowMeans(pha_ftir3[, c("PHA-74-60-days-rep1_tm", "PHA-74-60-days-rep2_tm")])
pha_ftir3$PHA_8_60days_avg_tm <- rowMeans(pha_ftir3[, c("PHA-8-60-days-rep1_tm", "PHA-8-60-days-rep2_tm")])
pha_ftir3$PHA_control_avg_tm <- rowMeans(pha_ftir3[, c("PHA-ctrl-rep1_tm", "PHA-ctrl-rep2_tm")])
pha_ftir4 <- pha_ftir3[,c(1,11,22:30)]
pha_ftir5 <- pha_ftir4

names(pha_ftir5) <- gsub("PHA-45A-20-days-rep1_wl", "wl", names(pha_ftir5))
names(pha_ftir5) <- gsub("PHA-53-60-days", "tPHA_G1_AS_Day60", names(pha_ftir5)) 
names(pha_ftir5) <- gsub("PHA_45A_20days", "PHA_G1_Day20", names(pha_ftir5)) # confirm with Magda, what is this?? does it have sludge
names(pha_ftir5) <- gsub("PHA_45A_4days", "PHA_G1_Day4", names(pha_ftir5)) 
names(pha_ftir5) <- gsub("PHA_47_60days_AF", "PHA_G1_AS_Day60", names(pha_ftir5)) 
names(pha_ftir5) <- gsub("PHA_52_60days", "tPHA_G1_AS_Day60", names(pha_ftir5))
names(pha_ftir5) <- gsub("PHA_53A_20days", "tPHA_G1_Day20", names(pha_ftir5)) # confirm with Magda, what is this?? does it have sludge
names(pha_ftir5) <- gsub("PHA_53A_4days", "tPHA_G1_Day4", names(pha_ftir5))
names(pha_ftir5) <- gsub("PHA_74_60days", "tPHA_AS_Day60", names(pha_ftir5))
names(pha_ftir5) <- gsub("PHA_8_60days", "PHA_AS_Day60", names(pha_ftir5))
names(pha_ftir5)
#PHA 52 and PHA 53 both tPHA_G1_AS_Day60, should i average?

cols <- c("PHA control"="#0072B2","PHA + AS Day 60"="#E69F00","tPHA + AS Day 60"="#56B4E9","PHA + G1 + AS Day 60"="#D55E00", "tPHA + G1 + AS Day 60 (a)"="#CC79A7","tPHA + G1 + AS Day 60 (b)"="#b66dff","tPHA + G1 Day 4"="#882255","PHA + G1 Day 4"="#332288","tPHA + G1 Day 20"="#A2DBB5","PHA + G1 Day 20"="#000000") 
phaft <- ggplot(pha_ftir5, aes(x = wl)) + geom_line(aes(y = PHA_control_avg_tm, group = 1, color = "PHA control")) + 
  geom_line(aes(y = PHA_AS_Day60_avg_tm, group = 1, color = "PHA + AS Day 60")) + 
  geom_line(aes(y = tPHA_AS_Day60_avg_tm, group = 1, color = "tPHA + AS Day 60")) + 
  geom_line(aes(y = PHA_G1_AS_Day60_avg_tm, group = 1, color = "PHA + G1 + AS Day 60")) + 
  geom_line(aes(y = `tPHA_G1_AS_Day60-rep1_tm`, group = 1, color = "tPHA + G1 + AS Day 60 (a)")) + 
  geom_line(aes(y = tPHA_G1_AS_Day60_avg_tm, group = 1, color = "tPHA + G1 + AS Day 60 (b)")) + 
  geom_line(aes(y = PHA_G1_Day4_avg_tm, group = 1, color = "PHA + G1 Day 4")) + 
  geom_line(aes(y = tPHA_G1_Day4_avg_tm, group = 1, color = "tPHA + G1 Day 4")) + 
  geom_line(aes(y = PHA_G1_Day20_avg_tm, group = 1, color = "PHA + G1 Day 20")) + 
  geom_line(aes(y = tPHA_G1_Day20_avg_tm, group = 1, color = "tPHA + G1 Day 20")) + 
  labs(
    #title = "PHA FTIR",
    x = expression("Wavenumber (cm"^-1*")"),
    y = "Transmittance (%)") +
  ylim(0, 100) +
  xlim(3999,525) +
  # sphale_x_reverse() +
  # sphale_x_continuous(breaks = c(4,32,60), expand = c(0, 0)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face = "bold", vjust = -2))
phaft

# plot PLA FTIR
pla_ftir3 <- pla_ftir2
pla_ftir3$PLA_30A_4days_avg_tm <- rowMeans(pla_ftir3[, c("30A-PLA-4days-rep1_tm", "30A-PLA-4days-rep2_tm")])
pla_ftir3$PLA_31_60days_avg_tm <- rowMeans(pla_ftir3[, c("31-pla-60days-rep1_tm", "pla-31-60-days-rep2_tm")])
pla_ftir3$PLA_31A_20days_avg_tm <- rowMeans(pla_ftir3[, c("PLA-31A--20-days-rep1_tm", "PLA-31A-20-days-rep2_tm")])
pla_ftir3$PLA_39_60days_avg_tm <- rowMeans(pla_ftir3[, c("pla-39-60-days-rep2_tm", "pla-39-60-days_tm")])
pla_ftir3$PLA_39A_20days_avg_tm <- rowMeans(pla_ftir3[, c("PLA-39A-20-days-rep1_tm", "PLA-39A-20-days-rep2_tm")])
pla_ftir3$PLA_5_60days_avg_tm <- rowMeans(pla_ftir3[, c("pla-5-60-days-rep1_tm", "pla-5-60-days-rep2_tm")])
pla_ftir3$PLA_72_60days_avg_tm <- rowMeans(pla_ftir3[, c("pla-72-60-days-rep1_tm", "pla-72-60-days-rep2_tm", "pla-72-60-days-rep3_tm")])
pla_ftir3$PLA_control_avg_tm <- rowMeans(pla_ftir3[, c("PLA-ctrl-rep1_tm", "PLA-ctrl-rep2_tm")])
pla_ftir4 <- pla_ftir3[,c(1,8,20:27)]
pla_ftir5 <- pla_ftir4

names(pla_ftir5) <- gsub("30A-PLA-4days-rep1_wl", "wl", names(pla_ftir5))
names(pla_ftir5) <- gsub("PLA-37A-4-days", "tPLA_S3_Day4", names(pla_ftir5)) 
names(pla_ftir5) <- gsub("PLA_30A_4days", "PLA_S3_Day4", names(pla_ftir5)) # confirm with Magda, what is this?? does it have sludge
names(pla_ftir5) <- gsub("PLA_31_60days", "PLA_S3_AS_Day60", names(pla_ftir5)) 
names(pla_ftir5) <- gsub("PLA_31A_20days", "PLA_S3_Day20", names(pla_ftir5)) # confirm with Magda, what is this?? does it have sludge
names(pla_ftir5) <- gsub("PLA_39_60days", "tPLA_S3_AS_Day60", names(pla_ftir5))
names(pla_ftir5) <- gsub("PLA_39A_20days", "tPLA_S3_Day20", names(pla_ftir5)) # confirm with Magda, what is this?? does it have sludge
names(pla_ftir5) <- gsub("PLA_5_60days", "PLA_AS_Day60", names(pla_ftir5))
names(pla_ftir5) <- gsub("PLA_72_60days", "tPLA_AS_Day60", names(pla_ftir5))
names(pla_ftir5)

cols <- c("PLA control"="#0072B2","PLA + AS Day 60"="#E69F00","tPLA + AS Day 60"="#56B4E9","PLA + S3 + AS Day 60"="#D55E00", "tPLA + S3 + AS Day 60"="#CC79A7","tPLA + S3 Day 4"="#882255","PLA + S3 Day 4"="#332288","tPLA + S3 Day 20"="#A2DBB5","PLA + S3 Day 20"="#000000") 
plaft <- ggplot(pla_ftir5, aes(x = wl)) + geom_line(aes(y = PLA_control_avg_tm, group = 1, color = "PLA control")) + 
  geom_line(aes(y = PLA_AS_Day60_avg_tm, group = 1, color = "PLA + AS Day 60")) + 
  geom_line(aes(y = tPLA_AS_Day60_avg_tm, group = 1, color = "tPLA + AS Day 60")) + 
  geom_line(aes(y = PLA_S3_AS_Day60_avg_tm, group = 1, color = "PLA + S3 + AS Day 60")) + 
  geom_line(aes(y = tPLA_S3_AS_Day60_avg_tm, group = 1, color = "tPLA + S3 + AS Day 60")) + 
  geom_line(aes(y = PLA_S3_Day4_avg_tm, group = 1, color = "PLA + S3 Day 4")) + 
  geom_line(aes(y = `tPLA_S3_Day4-rep1_tm`, group = 1, color = "tPLA + S3 Day 4")) + 
  geom_line(aes(y = PLA_S3_Day20_avg_tm, group = 1, color = "PLA + S3 Day 20")) + 
  geom_line(aes(y = tPLA_S3_Day20_avg_tm, group = 1, color = "tPLA + S3 Day 20")) + 
  labs(
    #title = "PLA FTIR",
    x = expression("Wavenumber (cm"^-1*")"),
    y = "Transmittance (%)") +
  ylim(0, 100) +
  xlim(3999,525) +
  # splale_x_reverse() +
  # splale_x_continuous(breaks = c(4,32,60), expand = c(0, 0)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face = "bold", vjust = -2))
plaft

ggsave(plot = caft, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/CA_FTIR_03012026.pdf", width = 10, height = 10, units = "in", dpi = 300)
ggsave(plot = caft, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/CA_FTIR_03012026.tiff", width = 10, height = 10, units = "in", dpi = 200, bg="white")

ggsave(plot = phaft, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/PHA_FTIR_03012026.pdf", width = 10, height = 10, units = "in", dpi = 300)
ggsave(plot = phaft, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/PHA_FTIR_03012026.tiff", width = 10, height = 10, units = "in", dpi = 200, bg="white")

ggsave(plot = plaft, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/PLA_FTIR_03012026.pdf", width = 10, height = 10, units = "in", dpi = 300)
ggsave(plot = plaft, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/PLA_FTIR_03012026.tiff", width = 10, height = 10, units = "in", dpi = 200, bg="white")
