#!/usr/bin/env Rscript

# Parse command-line arguments
args <- commandArgs(trailingOnly = TRUE)
# args <- c('data/65K_DEL-three_class', 'tt');

if (length(args) < 2) {
  stop("\nRscript this.r  in-group_col1_col2.tsv  output_prefix\n");
}

# Load necessary libraries
suppressMessages({
library(ggplot2)
library(dplyr)
library(reshape2)
});

data_file <- as.character(args[1]);  # First argument: Path to data file
out_pref  <- as.character(args[2]);

# Read the input data
data <- read.delim(data_file, header = TRUE, sep = "\t"); # First column name should be 'Group'.

# Reshape data for ggplot
data_long <- melt(data, id.vars = "Group", variable.name = "ColumnType", value.name = "Count");

# Calculate percentages
data_long <- data_long %>% 
  group_by(Group) %>% 
  mutate(Percentage = Count / sum(Count) * 100)

# Plot by input order.
data_long$Group <- factor(data_long$Group, levels=data$Group);
data_long$ColumnType <- factor(data_long$ColumnType, levels= colnames(data)[-1]);

# Plot the data
## Y axis as real sum count.
pdf(file=paste0(out_pref, "-barplot_count.pdf"))
p <- ggplot(data_long, aes(x = Group, y = Count, fill = ColumnType)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = "Group",
       y = "Count",
       fill = "ColumnType") +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(linewidth = 0.5, colour = 'black'),
    axis.ticks = element_line(linewidth = 0.5)
  )
print(p);
dev.off();

pdf(file=paste0(out_pref, "-barplot_percent.pdf"))
p <- ggplot(data_long, aes(x = Group, y = Percentage, fill = ColumnType)) +
  geom_bar(stat = "identity", position = "stack") +
  theme_minimal() +
  labs(x = "Group",
       y = "Percentage",
       fill = "ColumnType") +
  scale_y_continuous(labels = scales::percent_format(scale = 1)) +
  theme(
    panel.grid = element_blank(),
    axis.line  = element_line(linewidth = 0.5, colour = 'black'),
    axis.ticks = element_line(linewidth = 0.5)
  )
print(p)
dev.off();


