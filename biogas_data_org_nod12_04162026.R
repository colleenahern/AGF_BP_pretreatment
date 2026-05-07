# Colleen Ahern
# 04/16/2026
# Based on biogas_data_org_01192026.R file
# Remaking Magda's biogas graphs in R for figures but this time removing Day 12 data for the PLA and CA samples because of bad S3 + CA data for that day
# only remove day 12 data for PLA and CA because PHA uses G1 not S3

library(ggplot2)
library(readr)
library(ggpubr)
library(cowplot)

# Data for PHA samples
data <- read_csv("/Users/colleenahern/Documents/Magda_BPs_experiment/biogas_data/BP_expII_gas_organized_forR_01192026.csv")
colnames(data) <- gsub("PLA_70C", "tPLA", colnames(data))
colnames(data) <- gsub("PHA_70C", "tPHA", colnames(data))
colnames(data) <- gsub("CA_70C", "tCA", colnames(data))

# PHA + G1 + AS
# Make biogas plot

palette.colors(palette = "Okabe-Ito")

cols <- c("AS"="#0072B2","PHA + AS"="#E69F00","tPHA + AS"="#56B4E9", "PHA + G1 + AS"="#D55E00","tPHA + G1 + AS"="#CC79A7", "G1 + AS"="#009E73")
phalin <- ggplot(data, aes(x = Days)) + geom_line(aes(y = averageAS, group = 1, color = "AS")) + 
  geom_point(aes(y = averageAS, color = "AS")) +
  geom_errorbar(aes(ymax = averageAS + sdAS, 
                    ymin = averageAS - sdAS, color = "AS")) +
  geom_line(aes(y = averagePHAAS, group = 1, color = "PHA + AS")) + 
  geom_point(aes(y = averagePHAAS, color = "PHA + AS")) + 
  geom_errorbar(aes(ymax = averagePHAAS + sdPHAAS, 
                    ymin = averagePHAAS - sdPHAAS, color = "PHA + AS")) +
  geom_line(aes(y = averagePHA70CAS, group = 1, color = "tPHA + AS")) + 
  geom_point(aes(y = averagePHA70CAS, color = "tPHA + AS")) + 
  geom_errorbar(aes(ymax = averagePHA70CAS + sdPHA70CAS, 
                    ymin = averagePHA70CAS - sdPHA70CAS, color = "tPHA + AS")) +
  geom_line(aes(y = averagePHAG1AS, group = 1, color = "PHA + G1 + AS")) + 
  geom_point(aes(y = averagePHAG1AS, color = "PHA + G1 + AS")) + 
  geom_errorbar(aes(ymax = averagePHAG1AS + sdPHAG1AS, 
                    ymin = averagePHAG1AS - sdPHAG1AS, color = "PHA + G1 + AS")) +
  geom_line(aes(y = averagePHA70CG1AS, group = 1, color = "tPHA + G1 + AS")) + 
  geom_point(aes(y = averagePHA70CG1AS, color = "tPHA + G1 + AS")) +
  geom_errorbar(aes(ymax = averagePHA70CG1AS + sdPHA70CG1AS, 
                    ymin = averagePHA70CG1AS - sdPHA70CG1AS, color = "tPHA + G1 + AS")) +
  geom_line(aes(y = averageG1AS, group = 1, color = "G1 + AS")) + 
  geom_point(aes(y = averageG1AS, color = "G1 + AS")) + 
  geom_errorbar(aes(ymax = averageG1AS + sdG1AS, 
                    ymin = averageG1AS - sdG1AS, color = "G1 + AS")) +
  labs(
    title = "PHA",
    x = "Day",
    y = "Cumulative Biogas Production (mL)") +
  ylim(0, 250) +
  scale_x_continuous(breaks = c(4,32,60), expand = c(0, 0)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols, breaks = c("AS", "PHA + AS", "PHA + G1 + AS", "G1 + AS", "tPHA + AS", "tPHA + G1 + AS")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face = "bold", vjust = -2))

phalin

# Make bar plot with significance
data_subphag1 <- data[13,grepl("PHA \\+ AS|PHA \\+ G1 \\+ AS|tPHA \\+ G1 \\+ AS", colnames(data))]
data_subphag1 <- cbind(data[13, c(39:41,2:4)], data_subphag1) # G1 + AS and AS
data_subphag1 <- t(data_subphag1)
data_subphag1 <- as.data.frame(data_subphag1)
data_subphag1$group <- rownames(data_subphag1)
data_subphag1$group <- substr(data_subphag1$group, 1, nchar(data_subphag1$group) - 2)

# Basic bar plot with mean ± SE
data_subphag1$group <- factor(data_subphag1$group, levels = c("PHA + G1 + AS", "PHA + AS", "G1 + AS", "AS", "tPHA + AS", "tPHA + G1 + AS"))

phabar <- ggbarplot(data_subphag1, x = "group", y = "V1",
               add = "mean_se",
               fill = "group",
               # xlab = "Condition",
               ylab = "Total Biogas Production (mL)") +
  labs(x = NULL,
       title = "PHA",) +
  ylim(0,400) +
  scale_fill_manual(values = c("AS"="#0072B2",
                               "PHA + AS"="#E69F00",
                               "tPHA + AS"="#56B4E9", 
                               "G1 + AS"="#009E73", 
                               "PHA + G1 + AS"="#D55E00",
                               "tPHA + G1 + AS"="#CC79A7")) +
  scale_x_discrete(labels = c("AS"="AS",
                              "PHA + AS"="PHA +\nAS",
                              "tPHA + AS"="tPHA +\nAS", 
                              "G1 + AS"="G1 +\nAS", 
                              "PHA + G1 + AS"="PHA +\nG1 +\nAS",
                              "tPHA + G1 + AS"="tPHA +\nG1 +\nAS")) +
  stat_compare_means(
  method = "t.test",
  method.args = list(formula = log(V1) ~ group),
  comparisons = list(c("PHA + G1 + AS", "PHA + AS"), c("PHA + G1 + AS", "G1 + AS"), c("PHA + G1 + AS", "AS"),
                     c("tPHA + G1 + AS", "tPHA + AS"), c("tPHA + G1 + AS", "AS"), c("tPHA + G1 + AS", "G1 + AS")),
  label = "p.signif",  # shows *, **, ***
  label.y = c(240,270,300,330,360,390) 
) +   
  theme_classic() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        # axis.text.x = element_blank(),
        legend.position = "none"#, plot.margin = margin(10, 10, 33, 10)
        )

phabar

# Data for CA and PLA samples (no Day 12)
datanod12 <- read_csv("/Users/colleenahern/Documents/Magda_BPs_experiment/biogas_data/BP_expII_gas_organized_forR_noD12_CAPLA_04162026.csv")
rem <- colnames(datanod12)[grepl("ind", colnames(datanod12))]
datanod12f <- datanod12[,!(colnames(datanod12) %in% rem)]

# CA + S3 + AS
# Make biogas plot

palette.colors(palette = "Okabe-Ito")

cols <- c("AS"="#0072B2","CA + AS"="#E69F00","tCA + AS"="#56B4E9", "S3 + AS"="#009E73", "CA + S3 + AS"="#D55E00","tCA + S3 + AS"="#CC79A7")
calin <- ggplot(datanod12f, aes(x = Day)) + geom_line(aes(y = averageAS, group = 1, color = "AS")) + 
  geom_point(aes(y = averageAS, color = "AS")) +
  geom_errorbar(aes(ymax = averageAS + sdAS, 
                    ymin = averageAS - sdAS, color = "AS")) +
  geom_line(aes(y = averageCAAS, group = 1, color = "CA + AS")) + 
  geom_point(aes(y = averageCAAS, color = "CA + AS")) + 
  geom_errorbar(aes(ymax = averageCAAS + sdCAAS, 
                    ymin = averageCAAS - sdCAAS, color = "CA + AS")) +
  geom_line(aes(y = averagetCAAS, group = 1, color = "tCA + AS")) + 
  geom_point(aes(y = averagetCAAS, color = "tCA + AS")) + 
  geom_errorbar(aes(ymax = averagetCAAS + sdtCAAS, 
                    ymin = averagetCAAS - sdtCAAS, color = "tCA + AS")) +
  geom_line(aes(y = averageS3AS, group = 1, color = "S3 + AS")) + 
  geom_point(aes(y = averageS3AS, color = "S3 + AS")) + 
  geom_errorbar(aes(ymax = averageS3AS + sdS3AS, 
                    ymin = averageS3AS - sdS3AS, color = "S3 + AS")) +
  geom_line(aes(y = averageCAS3AS123, group = 1, color = "CA + S3 + AS")) + 
  geom_point(aes(y = averageCAS3AS123, color = "CA + S3 + AS")) + 
  geom_errorbar(aes(ymax = averageCAS3AS123 + sdCAS3AS123, 
                    ymin = averageCAS3AS123 - sdCAS3AS123, color = "CA + S3 + AS")) +
  geom_line(aes(y = averagetCAS3AS, group = 1, color = "tCA + S3 + AS")) + 
  geom_point(aes(y = averagetCAS3AS, color = "tCA + S3 + AS")) +
  geom_errorbar(aes(ymax = averagetCAS3AS + sdtCAS3AS, 
                    ymin = averagetCAS3AS - sdtCAS3AS, color = "tCA + S3 + AS")) +
  labs(
    title = "CA",
    x = "Day",
    y = "Cumulative Biogas Production (mL)") +
  ylim(0, 250) +
  scale_x_continuous(breaks = c(4,32,60), expand = c(0, 0)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face = "bold", vjust = -2))

calin

# Make bar plot with significance

data_subcas3nod12 <- datanod12f[12,grepl("CA \\+ AS|CA \\+ S3 \\+ AS|tCA \\+ S3 \\+ AS", colnames(datanod12f))]
data_subcas3nod12 <- cbind(datanod12f[12, c(2,3,6:8)], data_subcas3nod12) # S3 + AS and AS
data_subcas3nod12 <- t(data_subcas3nod12)
data_subcas3nod12 <- as.data.frame(data_subcas3nod12)
data_subcas3nod12$group <- rownames(data_subcas3nod12)
data_subcas3nod12$group <- substr(data_subcas3nod12$group, 1, nchar(data_subcas3nod12$group) - 2)

# Basic bar plot with mean ± SE
data_subcas3nod12$group <- factor(data_subcas3nod12$group, levels = c("CA + S3 + AS", "CA + AS", "S3 + AS", "AS", "tCA + AS", "tCA + S3 + AS"))

cabar <- ggbarplot(data_subcas3nod12, x = "group", y = "V1",
               add = "mean_se",
               fill = "group",
               # xlab = "Condition",
               ylab = "Total Biogas Production (mL)") +
  labs(x = NULL,
       title = "CA",) +
  ylim(0,400) +
  scale_fill_manual(values = c("AS"="#0072B2",
                               "CA + AS"="#E69F00",
                               "tCA + AS"="#56B4E9", 
                               "S3 + AS"="#009E73", 
                               "CA + S3 + AS"="#D55E00",
                               "tCA + S3 + AS"="#CC79A7")) +
  scale_x_discrete(labels = c("AS"="AS",
                              "CA + AS"="CA +\nAS",
                              "tCA + AS"="tCA +\nAS", 
                              "S3 + AS"="S3 +\nAS", 
                              "CA + S3 + AS"="CA +\nS3 +\nAS",
                              "tCA + S3 + AS"="tCA +\nS3 +\nAS")) +
  stat_compare_means(
  method = "t.test",
  method.args = list(formula = log(V1) ~ group),
  comparisons = list(c("CA + S3 + AS", "CA + AS"), c("CA + S3 + AS", "S3 + AS"), c("CA + S3 + AS", "AS"), 
                     c("tCA + S3 + AS", "tCA + AS"), c("tCA + S3 + AS", "AS"), c("tCA + S3 + AS", "S3 + AS")),
  label = "p.signif",  # shows *, **, ***
  label.y = c(240,270,300,330,360,390) 
) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        # axis.text.x = element_blank(),
        legend.position = "none"#, plot.margin = margin(10, 10, 33, 10)
        )
cabar

# make plots for CA + S3 + AS 13 instead of CA + S3 + AS 123 as backup
cols <- c("AS"="#0072B2","CA + AS"="#E69F00","tCA + AS"="#56B4E9", "S3 + AS"="#009E73", "CA + S3 + AS"="#D55E00","tCA + S3 + AS"="#CC79A7")
calin13 <- ggplot(datanod12f, aes(x = Day)) + geom_line(aes(y = averageAS, group = 1, color = "AS")) + 
  geom_point(aes(y = averageAS, color = "AS")) +
  geom_errorbar(aes(ymax = averageAS + sdAS, 
                    ymin = averageAS - sdAS, color = "AS")) +
  geom_line(aes(y = averageCAAS, group = 1, color = "CA + AS")) + 
  geom_point(aes(y = averageCAAS, color = "CA + AS")) + 
  geom_errorbar(aes(ymax = averageCAAS + sdCAAS, 
                    ymin = averageCAAS - sdCAAS, color = "CA + AS")) +
  geom_line(aes(y = averagetCAAS, group = 1, color = "tCA + AS")) + 
  geom_point(aes(y = averagetCAAS, color = "tCA + AS")) + 
  geom_errorbar(aes(ymax = averagetCAAS + sdtCAAS, 
                    ymin = averagetCAAS - sdtCAAS, color = "tCA + AS")) +
  geom_line(aes(y = averageS3AS, group = 1, color = "S3 + AS")) + 
  geom_point(aes(y = averageS3AS, color = "S3 + AS")) + 
  geom_errorbar(aes(ymax = averageS3AS + sdS3AS, 
                    ymin = averageS3AS - sdS3AS, color = "S3 + AS")) +
  geom_line(aes(y = averageCAS3AS13, group = 1, color = "CA + S3 + AS")) + 
  geom_point(aes(y = averageCAS3AS13, color = "CA + S3 + AS")) + 
  geom_errorbar(aes(ymax = averageCAS3AS13 + sdCAS3AS13, 
                    ymin = averageCAS3AS13 - sdCAS3AS13, color = "CA + S3 + AS")) +
  geom_line(aes(y = averagetCAS3AS, group = 1, color = "tCA + S3 + AS")) + 
  geom_point(aes(y = averagetCAS3AS, color = "tCA + S3 + AS")) +
  geom_errorbar(aes(ymax = averagetCAS3AS + sdtCAS3AS, 
                    ymin = averagetCAS3AS - sdtCAS3AS, color = "tCA + S3 + AS")) +
  labs(
    title = "CA",
    x = "Day",
    y = "Cumulative Biogas Production (mL)") +
  ylim(0, 250) +
  scale_x_continuous(breaks = c(4,32,60), expand = c(0, 0)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face = "bold", vjust = -2))

calin13

# Make bar plot with significance

data_subcas3nod12 <- datanod12f[12,grepl("CA \\+ AS|CA \\+ S3 \\+ AS|tCA \\+ S3 \\+ AS", colnames(datanod12f))]
data_subcas3nod12 <- cbind(datanod12f[12, c(2,3,6:8)], data_subcas3nod12) # S3 + AS and AS
data_subcas3nod12 <- t(data_subcas3nod12)
data_subcas3nod12 <- as.data.frame(data_subcas3nod12)
data_subcas3nod12$group <- rownames(data_subcas3nod12)
data_subcas3nod12$group <- substr(data_subcas3nod12$group, 1, nchar(data_subcas3nod12$group) - 2)
data_subcas3nod12_13 <- data_subcas3nod12[-13,]

# Basic bar plot with mean ± SE
data_subcas3nod12_13$group <- factor(data_subcas3nod12_13$group, levels = c("CA + S3 + AS", "CA + AS", "S3 + AS", "AS", "tCA + AS", "tCA + S3 + AS"))

cabar13 <- ggbarplot(data_subcas3nod12_13, x = "group", y = "V1",
                   add = "mean_se",
                   fill = "group",
                   # xlab = "Condition",
                   ylab = "Total Biogas Production (mL)") +
  labs(x = NULL,
       title = "CA",) +
  ylim(0,400) +
  scale_fill_manual(values = c("AS"="#0072B2",
                               "CA + AS"="#E69F00",
                               "tCA + AS"="#56B4E9", 
                               "S3 + AS"="#009E73", 
                               "CA + S3 + AS"="#D55E00",
                               "tCA + S3 + AS"="#CC79A7")) +
  scale_x_discrete(labels = c("AS"="AS",
                              "CA + AS"="CA +\nAS",
                              "tCA + AS"="tCA +\nAS", 
                              "S3 + AS"="S3 +\nAS", 
                              "CA + S3 + AS"="CA +\nS3 +\nAS",
                              "tCA + S3 + AS"="tCA +\nS3 +\nAS")) +
  stat_compare_means(
    method = "t.test",
    method.args = list(formula = log(V1) ~ group),
    comparisons = list(c("CA + S3 + AS", "CA + AS"), c("CA + S3 + AS", "S3 + AS"), c("CA + S3 + AS", "AS"), 
                       c("tCA + S3 + AS", "tCA + AS"), c("tCA + S3 + AS", "AS"), c("tCA + S3 + AS", "S3 + AS")),
    label = "p.signif",  # shows *, **, ***
    label.y = c(240,270,300,330,360,390) 
  ) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        # axis.text.x = element_blank(),
        legend.position = "none"#, plot.margin = margin(10, 10, 33, 10)
  )
cabar13


# PLA + S3 + AS
# Make biogas plot

palette.colors(palette = "Okabe-Ito")

cols <- c("AS"="#0072B2","PLA + AS"="#E69F00","tPLA + AS"="#56B4E9", "S3 + AS"="#009E73", "PLA + S3 + AS"="#D55E00","tPLA + S3 + AS"="#CC79A7")
plalin <- ggplot(datanod12f, aes(x = Day)) + geom_line(aes(y = averageAS, group = 1, color = "AS")) + 
  geom_point(aes(y = averageAS, color = "AS")) +
  geom_errorbar(aes(ymax = averageAS + sdAS, 
                    ymin = averageAS - sdAS, color = "AS")) +
  geom_line(aes(y = averagePLAAS, group = 1, color = "PLA + AS")) + 
  geom_point(aes(y = averagePLAAS, color = "PLA + AS")) + 
  geom_errorbar(aes(ymax = averagePLAAS + sdPLAAS, 
                    ymin = averagePLAAS - sdPLAAS, color = "PLA + AS")) +
  geom_line(aes(y = averagetPLAAS, group = 1, color = "tPLA + AS")) + 
  geom_point(aes(y = averagetPLAAS, color = "tPLA + AS")) + 
  geom_errorbar(aes(ymax = averagetPLAAS + sdtPLAAS, 
                    ymin = averagetPLAAS - sdtPLAAS, color = "tPLA + AS")) +
  geom_line(aes(y = averageS3AS, group = 1, color = "S3 + AS")) + 
  geom_point(aes(y = averageS3AS, color = "S3 + AS")) + 
  geom_errorbar(aes(ymax = averageS3AS + sdS3AS, 
                    ymin = averageS3AS - sdS3AS, color = "S3 + AS")) +
  geom_line(aes(y = averagePLAS3AS, group = 1, color = "PLA + S3 + AS")) + 
  geom_point(aes(y = averagePLAS3AS, color = "PLA + S3 + AS")) + 
  geom_errorbar(aes(ymax = averagePLAS3AS + sdPLAS3AS, 
                    ymin = averagePLAS3AS - sdPLAS3AS, color = "PLA + S3 + AS")) +
  geom_line(aes(y = averagetPLAS3AS, group = 1, color = "tPLA + S3 + AS")) + 
  geom_point(aes(y = averagetPLAS3AS, color = "tPLA + S3 + AS")) +
  geom_errorbar(aes(ymax = averagetPLAS3AS + sdtPLAS3AS, 
                    ymin = averagetPLAS3AS - sdtPLAS3AS, color = "tPLA + S3 + AS")) +
  labs(
    title = "PLA",
    x = "Day",
    y = "Cumulative Biogas Production (mL)") +
  ylim(0, 250) +
  scale_x_continuous(breaks = c(4,32,60), expand = c(0, 0)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols) +
  theme( axis.text = element_text(size = 12),
         axis.title = element_text(size = 12),
         plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
         legend.text = element_text(size=12),
         legend.title = element_text(size=12, face = "bold", vjust = -2))

plalin

# Make bar plots with significance

data_subplas3nod12 <- datanod12f[12,grepl("PLA \\+ AS|PLA \\+ S3 \\+ AS|tPLA \\+ S3 \\+ AS", colnames(datanod12f))]
data_subplas3nod12 <- cbind(datanod12f[12, c(2,3,6:8)], data_subplas3nod12) # S3 + AS and AS
data_subplas3nod12 <- t(data_subplas3nod12)
data_subplas3nod12 <- as.data.frame(data_subplas3nod12)
data_subplas3nod12$group <- rownames(data_subplas3nod12)
data_subplas3nod12$group <- substr(data_subplas3nod12$group, 1, nchar(data_subplas3nod12$group) - 2)

# Basic bar plot with mean ± SE
data_subplas3nod12$group <- factor(data_subplas3nod12$group, levels = c("PLA + S3 + AS", "PLA + AS", "S3 + AS", "AS", "tPLA + AS", "tPLA + S3 + AS"))

plabar <- ggbarplot(data_subplas3nod12, x = "group", y = "V1",
               add = "mean_se",
               fill = "group",
               # xlab = "Condition",
               ylab = "Total Biogas Production (mL)") +
  labs(x = NULL,
       title = "PLA",) +
  ylim(0,400) +
  scale_fill_manual(values = c("AS"="#0072B2",
                               "PLA + AS"="#E69F00",
                               "tPLA + AS"="#56B4E9", 
                               "S3 + AS"="#009E73", 
                               "PLA + S3 + AS"="#D55E00",
                               "tPLA + S3 + AS"="#CC79A7")) +
  scale_x_discrete(labels = c("AS"="AS",
                              "PLA + AS"="PLA +\nAS",
                              "tPLA + AS"="tPLA +\nAS", 
                              "S3 + AS"="S3 +\nAS", 
                              "PLA + S3 + AS"="PLA +\nS3 +\nAS",
                              "tPLA + S3 + AS"="tPLA +\nS3 +\nAS")) +
  stat_compare_means(
  method = "t.test",
  method.args = list(formula = log(V1) ~ group),
  comparisons = list(c("PLA + S3 + AS", "PLA + AS"), c("PLA + S3 + AS", "S3 + AS"), c("PLA + S3 + AS", "AS"), 
                     c("tPLA + S3 + AS", "tPLA + AS"), c("tPLA + S3 + AS", "AS"), c("tPLA + S3 + AS", "S3 + AS")),
  label = "p.signif",  # shows *, **, ***
  label.y = c(240,270,300,330,360,390) 
) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        axis.text.x = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        # axis.text.x = element_blank(),
        legend.position = "none"#, plot.margin = margin(10, 10, 33, 10)
        )
plabar

# Combine all plots into one figure
finbiogas <- plot_grid(phalin, phabar, calin, cabar, plalin, plabar, nrow = 3, ncol = 2, labels = c("A)", "B)", "C)", "D)", "E)", "F)"))
finbiogas
ggsave(plot = finbiogas, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Biogas_fig_noD12CAPLA_123_04162026.pdf", width = 10, height = 10, units = "in", dpi = 300)
ggsave(plot = finbiogas, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Biogas_fig_noD12CAPLA_123_04162026.tiff", width = 10, height = 10, units = "in", dpi = 200, bg="white")

finbiogas13 <- plot_grid(phalin, phabar, calin13, cabar13, plalin, plabar, nrow = 3, ncol = 2, labels = c("A)", "B)", "C)", "D)", "E)", "F)"))
finbiogas13
ggsave(plot = finbiogas13, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Biogas_fig_noD12CAPLA_13_04162026.pdf", width = 10, height = 10, units = "in", dpi = 300)
ggsave(plot = finbiogas13, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Biogas_fig_noD12CAPLA_13_04162026.tiff", width = 10, height = 10, units = "in", dpi = 200, bg="white")


##########################################################################################################################
# S3 CA monoculture
s3ca <- read_csv("/Users/colleenahern/Documents/Magda_BPs_experiment/biogas_data/S3CA_mono_gas_org_01272026.csv")
colnames(s3ca) <- gsub(" \\+ ", "_", colnames(s3ca))
colnames(s3ca) <- gsub(" ", "_", colnames(s3ca))

cbf_palette <- c("#785EF0", "#DC267F", "#648FFF")
cols <- c("S3"="#785EF0","CA + S3"="#DC267F","Blank media"="#648FFF")

v1 <- ggplot(s3ca, aes(x = Days)) + geom_line(aes(y = S3_media_avg_cp, group = 1, color = "S3")) + 
  geom_point(aes(y = S3_media_avg_cp, color = "S3")) +
  geom_errorbar(aes(ymax = S3_media_avg_cp + S3_media_std_cp, 
                    ymin = S3_media_avg_cp - S3_media_std_cp, color = "S3")) +
  geom_line(aes(y = S3_media_CA_avg_cp, group = 1, color = "CA + S3")) + 
  geom_point(aes(y = S3_media_CA_avg_cp, color = "CA + S3")) + 
  geom_errorbar(aes(ymax = S3_media_CA_avg_cp + S3_media_CA_std_cp, 
                    ymin = S3_media_CA_avg_cp - S3_media_CA_std_cp, color = "CA + S3")) +
  geom_line(aes(y = media_avg_cp, group = 1, color = "Blank media")) + 
  geom_point(aes(y = media_avg_cp, color = "Blank media")) + 
  geom_errorbar(aes(ymax = media_avg_cp + media_std_cp, 
                    ymin = media_avg_cp - media_std_cp, color = "Blank media")) +
  labs(
    #title = "CA",
    x = "Day",
    y = "Cumulative Biogas Production (psi)") +
  ylim(0, 10) +
  scale_x_continuous(breaks = c(0,7,13), expand = c(0, 0)) +
  scale_y_continuous(limits = c(-.7, 10), breaks = seq(0, 10, by = 2)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols, breaks = c("CA + S3", "S3", "Blank media")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face = "bold", vjust = -2),
        plot.margin = margin(23, 0, 5, 5))

v1

# G1 PHA monoculture
g1pha <- read_csv("/Users/colleenahern/Documents/Magda_BPs_experiment/biogas_data/G1PHA_mono_gas_org_01272026.csv")
colnames(g1pha) <- gsub(" \\+ ", "_", colnames(g1pha))
colnames(g1pha) <- gsub(" ", "_", colnames(g1pha))

cbf_palette <- c("#785EF0", "#DC267F", "#648FFF")
cols <- c("G1"="#785EF0","PHA + G1"="#DC267F","Blank media"="#648FFF")

w1 <- ggplot(g1pha, aes(x = Day)) + geom_line(aes(y = G1_media_avg_cp, group = 1, color = "G1")) + 
  geom_point(aes(y = G1_media_avg_cp, color = "G1")) +
  geom_errorbar(aes(ymax = G1_media_avg_cp + G1_media_std_cp, 
                    ymin = G1_media_avg_cp - G1_media_std_cp, color = "G1")) +
  geom_line(aes(y = G1_PHA_media_avg_cp, group = 1, color = "PHA + G1")) + 
  geom_point(aes(y = G1_PHA_media_avg_cp, color = "PHA + G1")) + 
  geom_errorbar(aes(ymax = G1_PHA_media_avg_cp + G1_PHA_media_std_cp, 
                    ymin = G1_PHA_media_avg_cp - G1_PHA_media_std_cp, color = "PHA + G1")) +
  geom_line(aes(y = media_avg_cp, group = 1, color = "Blank media")) + 
  geom_point(aes(y = media_avg_cp, color = "Blank media")) + 
  geom_errorbar(aes(ymax = media_avg_cp + media_std_cp, 
                    ymin = media_avg_cp - media_std_cp, color = "Blank media")) +
  labs(
    #title = "CA",
    x = "Day",
    y = "Cumulative Biogas Production (psi)") +
  ylim(-.7, 10) +
  scale_x_continuous(breaks = c(0,4,8)) +
  scale_y_continuous(limits = c(-.7, 10), breaks = seq(0, 10, by = 2)) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols, breaks = c("PHA + G1", "G1", "Blank media")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face = "bold", vjust = -2),
        plot.margin = margin(23, 0, 5, 5))

w1

library(ggpubr)
deg_plot4 <- ggarrange(v1, w1, ncol=1, nrow=2, labels = c("A)", "B)"))
deg_plot4


library(cowplot)
dp <- plot_grid(deg_plot4, hm1$gtable, nrow = 1, ncol = 2, rel_widths = c(1, 1.5), labels = c("", "C)"))
dp

ggsave(plot = dp, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Monoculture_S3CA_deseqheatmap_fig_02052026.pdf", width = 12, height = 11, units = "in", dpi = 300)
# THIS TIFF IS DPI 200, NEED DPI 300 FOR PUBLICATION
# I am making this tiff just so I can add it to the shared word doc for Magda. Use the pdf when submitting to the journal
ggsave(plot = dp, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Monoculture_S3CA_deseqheatmap_fig_02052026.tiff", width = 12, height = 11, units = "in", dpi = 200, bg="white")

# 02/05/2026 bringing sem fig into here
library(magick)
mypic <- image_read("/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/SEM_fig_size12_02052026.png")
sem1 <- rasterGrob(mypic, interpolate=TRUE)
sem1 <- as_ggplot(sem1)
sem1 <- sem1 + theme(plot.margin = margin(5, 5, 5, 5))
sem1

ggsave(plot = sem1, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/SEM_ggplot_fig_size12_02052026.pdf", width = 10, height = 6, units = "in", dpi = 300)
ggsave(plot = sem1, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/SEM_ggplot_fig_size12_02052026.tiff", width = 10, height = 6, units = "in", dpi = 300, bg="white")


# 04/14/26
# plot CA + S3 + AS 1 to see discrepancies
data2 <- data[,16:17]
data2$Day <- c(4,6,9,12,16,21,24,27,31,34,41,47,60)
colnames(data2) <- c("V1", "V3", "Day")
cols <- c("S3 + AS 1" = "#DC267F", "S3 + AS 3" = "#648FFF")

d21  <- ggplot(data2, aes(x = Day)) + geom_line(aes(y = V1, color = "S3 + AS 1")) + 
  geom_point(aes(y = V1, color = "S3 + AS 1")) +
  geom_line(aes(y = V3, color = "S3 + AS 3")) + 
  geom_point(aes(y = V3, color = "S3 + AS 3")) +
  labs(
    x = "Day",
    y = "Cumulative Biogas Production (psi)") +
  ylim(0, 250) +
  scale_x_continuous(breaks = c(4,32,60), expand = expansion(mult = c(0, 0.03))) +
  theme_classic() +
  scale_color_manual(name = "Condition", values = cols, breaks = c("S3 + AS 1", "S3 + AS 3")) +
  theme(axis.text = element_text(size = 12),
        axis.title = element_text(size = 12),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        legend.text = element_text(size=12),
        legend.title = element_text(size=12, face = "bold", vjust = -2))


d21
