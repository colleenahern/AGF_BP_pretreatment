# Colleen Ahern
# 01/19/2026
# Remaking Magda's biogas graphs in R for figures

library(ggplot2)
library(readr)

data <- read_csv("/Users/colleenahern/Documents/Magda_BPs_experiment/biogas_data/BP_expII_gas_organized_forR_01192026.csv")
colnames(data) <- gsub("PLA_70C", "tPLA", colnames(data))
colnames(data) <- gsub("PHA_70C", "tPHA", colnames(data))
colnames(data) <- gsub("CA_70C", "tCA", colnames(data))

# T-tests to determine statistical significance of all comparisons - 04/07/2026
# G1-PHA pairing
data_sub <- data[13,grepl("PHA \\+ AS|PHA \\+ G1 \\+ AS|tPHA \\+ G1 \\+ AS", colnames(data))]
data_sub <- cbind(data[13, c(39:41,2:4)], data_sub) # G1 + AS and AS

data_sub <- data[13,grepl("CA \\+ AS|CA \\+ S3 \\+ AS|tCA \\+ S3 \\+ AS", colnames(data))]
data_sub <- cbind(data[13, c(16:17,2:4)], data_sub) # S3 + AS and AS

data_sub <- data[13,grepl("PLA \\+ AS|PLA \\+ S3 \\+ AS|tPLA \\+ S3 \\+ AS", colnames(data))]
data_sub <- cbind(data[13, c(16:17,2:4)], data_sub) # S3 + AS and AS
data_sub <- t(data_sub)

# making significance bar plots - 04/13/2026

# install.packages("ggplot2")
# install.packages("ggpubr")

library(ggplot2)
library(ggpubr)

# PHA + G1 + AS
data_subphag1 <- data[13,grepl("PHA \\+ AS|PHA \\+ G1 \\+ AS|tPHA \\+ G1 \\+ AS", colnames(data))]
data_subphag1 <- cbind(data[13, c(39:41,2:4)], data_subphag1) # G1 + AS and AS
data_subphag1 <- t(data_subphag1)
data_subphag1 <- as.data.frame(data_subphag1)
data_subphag1$group <- rownames(data_subphag1)
data_subphag1$group <- substr(data_subphag1$group, 1, nchar(data_subphag1$group) - 2)

# Basic bar plot with mean ± SE
data_subphag1$group <- factor(data_subphag1$group, levels = c("PHA + G1 + AS", "PHA + AS", "G1 + AS", "AS", "tPHA + AS", "tPHA + G1 + AS"))

p <- ggbarplot(data_subphag1, x = "group", y = "V1",
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
                               "tPHA + G1 + AS"="#CC79A7"))

# Add significance comparisons
p <- p + stat_compare_means(
  method = "t.test",
  method.args = list(formula = log(V1) ~ group),
  # comparisons = list(c("PHA + G1 + AS", "AS"), c("PHA + G1 + AS", "G1 + AS"), c("PHA + G1 + AS", "PHA + AS"), c("tPHA + G1 + AS", "G1 + AS"), c("tPHA + G1 + AS", "AS"), c("tPHA + G1 + AS", "tPHA + AS")),
  comparisons = list(c("PHA + G1 + AS", "PHA + AS"), c("PHA + G1 + AS", "G1 + AS"), c("PHA + G1 + AS", "AS"),
                     c("tPHA + G1 + AS", "tPHA + AS"), c("tPHA + G1 + AS", "AS"), c("tPHA + G1 + AS", "G1 + AS")),
  label = "p.signif",  # shows *, **, ***
  label.y = c(265,290,315,340,365,390) 
) +   
  theme_classic() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        axis.text.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(10, 10, 33, 10))

p

# CA + S3 + AS
data_subcas3 <- data[13,grepl("CA \\+ AS|CA \\+ S3 \\+ AS|tCA \\+ S3 \\+ AS", colnames(data))]
data_subcas3 <- cbind(data[13, c(16:17,2:4)], data_subcas3) # S3 + AS and AS
data_subcas3 <- t(data_subcas3)
data_subcas3 <- as.data.frame(data_subcas3)
data_subcas3$group <- rownames(data_subcas3)
data_subcas3$group <- substr(data_subcas3$group, 1, nchar(data_subcas3$group) - 2)

# Basic bar plot with mean ± SE
data_subcas3$group <- factor(data_subcas3$group, levels = c("CA + S3 + AS", "CA + AS", "S3 + AS", "AS", "tCA + AS", "tCA + S3 + AS"))

q <- ggbarplot(data_subcas3, x = "group", y = "V1",
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
                               "tCA + S3 + AS"="#CC79A7"))

# Add significance comparisons
q <- q + stat_compare_means(
  method = "t.test",
  method.args = list(formula = log(V1) ~ group),
  # comparisons = list(c("CA + S3 + AS", "AS"), c("CA + S3 + AS", "S3 + AS"), c("CA + S3 + AS", "CA + AS"), c("tCA + S3 + AS", "S3 + AS"), c("tCA + S3 + AS", "AS"), c("tCA + S3 + AS", "tCA + AS")),
  comparisons = list(c("CA + S3 + AS", "CA + AS"), c("CA + S3 + AS", "S3 + AS"), c("CA + S3 + AS", "AS"), 
                     c("tCA + S3 + AS", "tCA + AS"), c("tCA + S3 + AS", "AS"), c("tCA + S3 + AS", "S3 + AS")),
  label = "p.signif",  # shows *, **, ***
  label.y = c(265,290,315,340,365,390) 
) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        axis.text.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(10, 10, 33, 10))
q

# PLA + S3 + AS
data_subplas3 <- data[13,grepl("PLA \\+ AS|PLA \\+ S3 \\+ AS|tPLA \\+ S3 \\+ AS", colnames(data))]
data_subplas3 <- cbind(data[13, c(16:17,2:4)], data_subplas3) # S3 + AS and AS
data_subplas3 <- t(data_subplas3)
data_subplas3 <- as.data.frame(data_subplas3)
data_subplas3$group <- rownames(data_subplas3)
data_subplas3$group <- substr(data_subplas3$group, 1, nchar(data_subplas3$group) - 2)

# Basic bar plot with mean ± SE
data_subplas3$group <- factor(data_subplas3$group, levels = c("PLA + S3 + AS", "PLA + AS", "S3 + AS", "AS", "tPLA + AS", "tPLA + S3 + AS"))

r <- ggbarplot(data_subplas3, x = "group", y = "V1",
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
                               "tPLA + S3 + AS"="#CC79A7"))
 
# Add significance comparisons
r <- r + stat_compare_means(
    method = "t.test",
    method.args = list(formula = log(V1) ~ group),
    # comparisons = list(c("PLA + S3 + AS", "AS"), c("PLA + S3 + AS", "S3 + AS"), c("PLA + S3 + AS", "PLA + AS"), c("tPLA + S3 + AS", "S3 + AS"), c("tPLA + S3 + AS", "AS"), c("tPLA + S3 + AS", "tPLA + AS")),
    comparisons = list(c("PLA + S3 + AS", "PLA + AS"), c("PLA + S3 + AS", "S3 + AS"), c("PLA + S3 + AS", "AS"), 
                       c("tPLA + S3 + AS", "tPLA + AS"), c("tPLA + S3 + AS", "AS"), c("tPLA + S3 + AS", "S3 + AS")),
    label = "p.signif",  # shows *, **, ***
    label.y = c(265,290,315,340,365,390) 
  ) +
  theme_classic() +
  theme(axis.text.y = element_text(size = 12),
        axis.title.y = element_text(size = 12),
        # axis.text.x = element_text(angle = 45, hjust = 1),
        plot.title = element_text(size = 14, hjust = 0.5, face = "bold"),
        axis.text.x = element_blank(),
        legend.position = "none",
        plot.margin = margin(10, 10, 33, 10))
r

library(cowplot)
ss <- plot_grid(p, q, r, nrow = 2, ncol = 2, labels = c("A)", "B)", "C)"))
ss

# keep <- c("PHA + G1 + AS 1", "PHA + G1 + AS 2", "PHA + G1 + AS 3", "G1 + AS 1", "G1 + AS 2", "G1 + AS 3")
# keep <- c("PHA + G1 + AS 1", "PHA + G1 + AS 2", "PHA + G1 + AS 3", "AS 1", "AS 2", "AS 3")
# keep <- c("PHA + G1 + AS 1", "PHA + G1 + AS 2", "PHA + G1 + AS 3", "PHA + AS 1", "PHA + AS 2")
# 
# keep <- c("tPHA + G1 + AS 2", "tPHA + G1 + AS 3", "G1 + AS 1", "G1 + AS 2", "G1 + AS 3")
# keep <- c("tPHA + G1 + AS 2", "tPHA + G1 + AS 3", "AS 1", "AS 2", "AS 3")
# keep <- c("tPHA + G1 + AS 2", "tPHA + G1 + AS 3", "tPHA + AS 1", "tPHA + AS 2")

# keep <- c("CA + S3 + AS 1", "CA + S3 + AS 3", "S3 + AS 1", "S3 + AS 2")
# keep <- c("CA + S3 + AS 1", "CA + S3 + AS 2", "CA + S3 + AS 3", "AS 1", "AS 2", "AS 3")
# keep <- c("CA + S3 + AS 1", "CA + S3 + AS 2", "CA + S3 + AS 3", "CA + AS 1", "CA + AS 2")
# 
# keep <- c("tCA + S3 + AS 1", "tCA + S3 + AS 2", "tCA + S3 + AS 3", "S3 + AS 1", "S3 + AS 2", "S3 + AS 3")
# keep <- c("tCA + S3 + AS 1", "tCA + S3 + AS 2", "tCA + S3 + AS 3", "AS 1", "AS 2", "AS 3")
# keep <- c("tCA + S3 + AS 1", "tCA + S3 + AS 2", "tCA + S3 + AS 3", "tCA + AS 1", "tCA + AS 2", "tCA + AS 3")

keep <- c("PLA + S3 + AS 1", "PLA + S3 + AS 2", "PLA + S3 + AS 3", "S3 + AS 1", "S3 + AS 2", "S3 + AS 3")
keep <- c("PLA + S3 + AS 1", "PLA + S3 + AS 2", "PLA + S3 + AS 3", "AS 1", "AS 2", "AS 3")
keep <- c("PLA + S3 + AS 1", "PLA + S3 + AS 2", "PLA + S3 + AS 3", "PLA + AS 1", "PLA + AS 2", "PLA + AS 3")

keep <- c("tPLA + S3 + AS 1", "tPLA + S3 + AS 2", "tPLA + S3 + AS 3", "S3 + AS 1", "S3 + AS 2", "S3 + AS 3")
keep <- c("tPLA + S3 + AS 1", "tPLA + S3 + AS 2", "tPLA + S3 + AS 3", "AS 1", "AS 2", "AS 3")
keep <- c("tPLA + S3 + AS 1", "tPLA + S3 + AS 2", "tPLA + S3 + AS 3", "tPLA + AS 1", "tPLA + AS 2", "tPLA + AS 3")

data_subsub <- as.data.frame(data_sub[rownames(data_sub) %in% keep,])
data_subsub$group <- c(rep("A",3), rep("B",3))
# data_subsub$group <- c(rep("A",2), rep("B",3))
# data_subsub$group <- c(rep("A",3), rep("B",2))
# data_subsub$group <- c(rep("A",2), rep("B",2))

colnames(data_subsub) <- c("value", "group")
t.test(value ~ group, data = data_subsub)

# T-test results: not log transformed
# PHA + G1 + AS is statistically significantly more than G1 + AS; p-value = 0.01143
# PHA + G1 + AS is statistically significantly more than AS; p-value = 0.0329
# PHA + G1 + AS is statistically significantly more than PHA + AS; p-value = 0.03367

# tPHA + G1 + AS is statistically significantly more than G1 + AS; p-value = 0.007606
# tPHA + G1 + AS is statistically significantly more than AS; p-value = 0.002166
# tPHA + G1 + AS is NOT statistically significantly more than tPHA + AS

# CA + S3 + AS NOT statistically significantly more than S3 + AS
# CA + S3 + AS is statistically significantly more than AS
# CA + S3 + AS is statistically significantly more than CA + AS

# tCA + S3 + AS NOT statistically significantly more than S3 + AS
# tCA + S3 + AS is statistically significantly more than AS
# tCA + S3 + AS is statistically significantly more than tCA + AS

# PLA + S3 + AS NOT statistically significantly more than S3 + AS
# PLA + S3 + AS is statistically significantly more than AS
# PLA + S3 + AS is statistically significantly more than PLA + AS


palette.colors(palette = "Okabe-Ito")

cols <- c("AS"="#0072B2","PLA + AS"="#E69F00","tPLA + AS"="#56B4E9", "S3 + AS"="#009E73", "PLA + S3 + AS"="#D55E00","tPLA + S3 + AS"="#CC79A7")
p1 <- ggplot(data, aes(x = Days)) + geom_line(aes(y = averageAS, group = 1, color = "AS")) + 
  geom_point(aes(y = averageAS, color = "AS")) +
  geom_errorbar(aes(ymax = averageAS + sdAS, 
                    ymin = averageAS - sdAS, color = "AS")) +
  geom_line(aes(y = averagePLAAS, group = 1, color = "PLA + AS")) + 
  geom_point(aes(y = averagePLAAS, color = "PLA + AS")) + 
  geom_errorbar(aes(ymax = averagePLAAS + sdPLAAS, 
                    ymin = averagePLAAS - sdPLAAS, color = "PLA + AS")) +
  geom_line(aes(y = averagePLA70CAS, group = 1, color = "tPLA + AS")) + 
  geom_point(aes(y = averagePLA70CAS, color = "tPLA + AS")) + 
  geom_errorbar(aes(ymax = averagePLA70CAS + sdPLA70CAS, 
                    ymin = averagePLA70CAS - sdPLA70CAS, color = "tPLA + AS")) +
  geom_line(aes(y = averageS3AS, group = 1, color = "S3 + AS")) + 
  geom_point(aes(y = averageS3AS, color = "S3 + AS")) + 
  geom_errorbar(aes(ymax = averageS3AS + sdS3AS, 
                    ymin = averageS3AS - sdS3AS, color = "S3 + AS")) +
  geom_line(aes(y = averagePLAS3AS, group = 1, color = "PLA + S3 + AS")) + 
  geom_point(aes(y = averagePLAS3AS, color = "PLA + S3 + AS")) + 
  geom_errorbar(aes(ymax = averagePLAS3AS + sdPLAS3AS, 
                    ymin = averagePLAS3AS - sdPLAS3AS, color = "PLA + S3 + AS")) +
  geom_line(aes(y = averagePLA70CS3AS, group = 1, color = "tPLA + S3 + AS")) + 
  geom_point(aes(y = averagePLA70CS3AS, color = "tPLA + S3 + AS")) +
  geom_errorbar(aes(ymax = averagePLA70CS3AS + sdPLA70CS3AS, 
                    ymin = averagePLA70CS3AS - sdPLA70CS3AS, color = "tPLA + S3 + AS")) +
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
                                    
p1

cols <- c("AS"="#0072B2","PHA + AS"="#E69F00","tPHA + AS"="#56B4E9", "PHA + G1 + AS"="#D55E00","tPHA + G1 + AS"="#CC79A7", "G1 + AS"="#009E73")
q1 <- ggplot(data, aes(x = Days)) + geom_line(aes(y = averageAS, group = 1, color = "AS")) + 
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

q1

cols <- c("AS"="#0072B2","CA + AS"="#E69F00","tCA + AS"="#56B4E9", "S3 + AS"="#009E73", "CA + S3 + AS"="#D55E00","tCA + S3 + AS"="#CC79A7")
r1 <- ggplot(data, aes(x = Days)) + geom_line(aes(y = averageAS, group = 1, color = "AS")) + 
  geom_point(aes(y = averageAS, color = "AS")) +
  geom_errorbar(aes(ymax = averageAS + sdAS, 
                    ymin = averageAS - sdAS, color = "AS")) +
  geom_line(aes(y = averageCAAS, group = 1, color = "CA + AS")) + 
  geom_point(aes(y = averageCAAS, color = "CA + AS")) + 
  geom_errorbar(aes(ymax = averageCAAS + sdCAAS, 
                    ymin = averageCAAS - sdCAAS, color = "CA + AS")) +
  geom_line(aes(y = averageCA70CAS, group = 1, color = "tCA + AS")) + 
  geom_point(aes(y = averageCA70CAS, color = "tCA + AS")) + 
  geom_errorbar(aes(ymax = averageCA70CAS + sdCA70CAS, 
                    ymin = averageCA70CAS - sdCA70CAS, color = "tCA + AS")) +
  geom_line(aes(y = averageS3AS, group = 1, color = "S3 + AS")) + 
  geom_point(aes(y = averageS3AS, color = "S3 + AS")) + 
  geom_errorbar(aes(ymax = averageS3AS + sdS3AS, 
                    ymin = averageS3AS - sdS3AS, color = "S3 + AS")) +
  geom_line(aes(y = averageCAS3AS, group = 1, color = "CA + S3 + AS")) + 
  geom_point(aes(y = averageCAS3AS, color = "CA + S3 + AS")) + 
  geom_errorbar(aes(ymax = averageCAS3AS + sdCAS3AS, 
                    ymin = averageCAS3AS - sdCAS3AS, color = "CA + S3 + AS")) +
  geom_line(aes(y = averageCA70CS3AS, group = 1, color = "tCA + S3 + AS")) + 
  geom_point(aes(y = averageCA70CS3AS, color = "tCA + S3 + AS")) +
  geom_errorbar(aes(ymax = averageCA70CS3AS + sdCA70CS3AS, 
                    ymin = averageCA70CS3AS - sdCA70CS3AS, color = "tCA + S3 + AS")) +
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

r1

# 04/14/2026 adding sig plots to my growth curves
psps <- plot_grid(q1, p, r1, q, p1, r, nrow = 3, ncol = 2, labels = c("A)", "B)", "C)", "D)", "E)", "F)"))
psps
ggsave(plot = psps, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Biogas_fig_04142026.pdf", width = 12, height = 12, units = "in", dpi = 300)
ggsave(plot = psps, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Biogas_fig_04142026.tiff", width = 12, height = 12, units = "in", dpi = 200, bg="white")


# read in part d pic
# install.packages("magick")
library(magick)
mypic <- image_read("/Users/colleenahern/Documents/Magda_BPs_experiment/biogas_data/biodeg_pic02052026_b.png")
s1 <- rasterGrob(mypic, interpolate=TRUE)
s1 <- as_ggplot(s1)
s1 <- s1 + theme(plot.margin = margin(20, 0, 5, 0)) 
s1

# tp1 <- p1  + r1 + q1 + s1 + plot_layout(ncol = 2)
tp1 <- plot_grid(p1, r1, q1, s1, nrow = 2, ncol = 2, labels = c("A)", "B)", "C)", "D)"))
tp1

ggsave(plot = tp1, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Biogas_fig_02052026.pdf", width = 10, height = 10, units = "in", dpi = 300)
ggsave(plot = tp1, filename = "/Users/colleenahern/Documents/Magda_BPs_experiment/Figures/Biogas_fig_02052026.tiff", width = 10, height = 10, units = "in", dpi = 200, bg="white")


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


# Hello world!
# Hello hello :)
