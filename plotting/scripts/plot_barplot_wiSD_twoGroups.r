#!/usr/bin/env Rscript

suppressMessages({
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(gridExtra)
})


dta <- read.delim("in.table-for_grouped_barplot_with_SD", header = TRUE, sep = "\t")

dta_long <- dta %>%
  pivot_longer(cols = starts_with("Rep.."), names_to = "Replicate", values_to = "Expression")
dta_long$Individual <- factor(dta_long$Individual, levels=unique(dta_long$Individual))

# Compute mean and standard deviation for each Individual
dta_summary <- dta_long %>%
  group_by(Group, Individual) %>%
  summarise(Mean_Expression = mean(Expression), SD_Expression = sd(Expression), .groups = 'drop')


p1 <- ggplot(dta_summary, aes(x = Group, y = Mean_Expression, fill = Individual)) +
  geom_bar(stat = "identity", position = position_dodge(width = 0.9), width = 0.7) + 
  geom_errorbar(aes(ymin = Mean_Expression - SD_Expression, ymax = Mean_Expression + SD_Expression),
                width = 0.3, position = position_dodge(width = 0.9)) +  # Add SD error bars
  labs(y = "Caratenoids detected (mg/kg, FW)", x = "Individual") +
  theme_minimal() +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1),
    panel.grid = element_blank(),
    axis.line  = element_line(linewidth = 0.5, colour = 'black'),
    axis.ticks = element_line(linewidth = 0.5)
)

print (p1)
pdf(file= "in.table-for_grouped_barplot_with_SD.pdf", family= "Helvetica", width= 3.5, height= 3.5)
print (p1)
dev.off()

