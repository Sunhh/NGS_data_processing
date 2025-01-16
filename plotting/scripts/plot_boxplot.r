#!/usr/bin/env Rscript
# Load necessary libraries
suppressMessages({ library(optparse) });

# Define command-line options
option_list <- list(
  make_option(c("-p", "--pheno_files"), type = "character", help = "Comma-separated list of phenotypic data files", metavar = "FILES"),
  make_option(c("-t", "--test_groups"), type = "character", help = ":-separated groups with columns separated with , ;"),
  make_option(c("-o", "--output_prefix"), type = "character", help = "Prefix of output files.")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

# opt <- list(
# pheno_files = 'FleshBrix_22HN_1,FleshBrix_22HN_2,FleshBrix_19YQ_1,FleshBrix_19YQ_2',
# test_groups = 'FleshBrix_22HN_1,FleshBrix_22HN_2:FleshBrix_19YQ_1,FleshBrix_19YQ_2',
# output_prefix = 'test_out'
# );

# Check for mandatory arguments
if (is.null(opt$pheno_files) || is.null(opt$test_groups) ) {
  print_help(opt_parser)
  stop("\nMandatory arguments missing! Provide phenotypic data files and group definitions.\n")
}

if (is.null(opt$output_prefix)) {
  opt$output_prefix <- 'test_out'
}

# Load necessary libraries which may take long time.
suppressMessages({
library(ggplot2)
library(dplyr)
library(reshape2)
library(stats)
});


# Load data files
phenotype_files <- strsplit(opt$pheno_files, ",")[[1]]
if (length(phenotype_files) == 0) {
  stop("No phenotypic data files provided!")
}

# Generate grouped columns.
groups <- strsplit(opt$test_groups, ":")[[1]]
grouped_columns <- lapply(groups, function(group) strsplit(group, ",")[[1]])

# Combine data from files
melted_data <- data.frame()
for (file1 in phenotype_files) {
  file1_data <- read.delim(file1, header = F, col.names = "Value")
  file1_data$Phenotype <- gsub(".txt$", "", basename(file1))
  melted_data <- rbind(melted_data, file1_data)
}

melted_data$Phenotype <- factor(melted_data$Phenotype, levels= phenotype_files)

# Calculate mean values and sample counts for each group
summary_values <- melted_data %>%
  group_by(Phenotype) %>%
  summarise(mean_phenotype = mean(Value, na.rm = TRUE),
            sample_count = n())

# Create violin plot
violin_plot <- ggplot(melted_data, aes(x = Phenotype, y = Value, fill = Phenotype)) +
  geom_violin(trim = FALSE) +
  geom_jitter(width = 0.2, alpha = 0.7) +
  geom_crossbar(data = summary_values, aes(x = as.factor(Phenotype), 
                                        y = mean_phenotype, 
                                        ymin = mean_phenotype, 
                                        ymax = mean_phenotype),
                width = 0.5, color = "red", linewidth = 0.8) +
  theme_minimal() +
  labs(x = "Phenotype data set", y = "Phenotypic value") +
  theme(
    panel.grid = element_blank(), # Remove grids.
    axis.line = element_line(linewidth = 0.5, colour = "black"), # Add X and Y axises.
    axis.ticks = element_line(linewidth = 0.5) # Add ticks on X and Y axises.
  )

box_plot <- ggplot(melted_data, aes(x = Phenotype, y = Value, fill = Phenotype)) +
  geom_boxplot(width = 0.1, position = position_dodge(width = 0.75), outlier.shape = NA) +
  theme_minimal() +
  labs(x = "Phenotype data set", y = "Phenotypic value") +
  theme(
    panel.grid = element_blank(), # Remove grids.
    axis.line = element_line(linewidth = 0.5, colour = "black"), # Add X and Y axises.
    axis.ticks = element_line(linewidth = 0.5) # Add ticks on X and Y axises.
  )

pdf(file= paste0(opt$output_prefix, "-violinplot.pdf"))
print(violin_plot)
dev.off()
pdf(file= paste0(opt$output_prefix, "-boxplot.pdf"))
print(box_plot)
dev.off()

# Compute and return P-values for grouped columns
p_values <- data.frame()

for (i in seq_along(grouped_columns)) {
  group1 <- grouped_columns[[i]]
  group1_data <- melted_data %>% filter(Phenotype %in% group1)
  if (length(unique(group1_data$Phenotype)) == 2) {
    p_value <- t.test(Value ~ Phenotype, data = group1_data)$p.value
    p_values <- rbind(p_values, data.frame(Group = paste(group1, collapse = ","), P_value = p_value))
  } else {
    warning(paste("Group", i, "does not have exactly two levels, skipping t-test."))
  }
}

# Save P-values to a file
write.table(p_values, paste0(opt$output_prefix, "-Pval.tsv"), row.names = FALSE, quote = FALSE, sep = "\t");

# Print and save summary values
write.table(summary_values, paste0(opt$output_prefix, "-summary.tsv"), row.names = FALSE, quote = FALSE, sep = "\t");

