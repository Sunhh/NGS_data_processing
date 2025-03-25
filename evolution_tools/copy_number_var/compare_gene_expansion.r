#!/usr/bin/env Rscript

# Load necessary libraries
suppressMessages({
if (!requireNamespace("optparse", quietly = TRUE)) install.packages("optparse")
if (!requireNamespace("dplyr", quietly = TRUE)) install.packages("dplyr")
library(optparse)
library(dplyr)
});

# Command-line options
option_list <- list(
  make_option(c("-a", "--inputA"), type = "character", help = "Input file A (gene abundance table)"),
  make_option(c("-b", "--inputB"), type = "character", help = "Input file B (individual grouping info)"),
  make_option(c("-f", "--groupF"), type = "character", help = "Group From"),
  make_option(c("-t", "--groupT"), type = "character", help = "Group To"),
  make_option(c("-o", "--output"), type = "character", help = "Output file")
)

opt <- parse_args(OptionParser(option_list = option_list))

# Load input data
gene_data <- read.table(opt$inputA, header = TRUE, sep = "\t")
group_data <- read.table(opt$inputB, header = FALSE, sep = "\t")

# Ensure matching individuals
common_individuals <- intersect(colnames(gene_data)[-1], group_data$V1)
gene_data <- gene_data %>% select(OGID, all_of(common_individuals))
group_data <- group_data %>% filter(V1 %in% common_individuals)

# Group individuals by F and T
group_data <- group_data %>% filter(V2 %in% c(opt$groupF, opt$groupT))
group_f <- group_data$V1[group_data$V2 == opt$groupF]
group_t <- group_data$V1[group_data$V2 == opt$groupT]
size_f <- length(group_f);
size_t <- length(group_t);
F_high <- paste0(opt$groupF, "_high");
T_high <- paste0(opt$groupT, "_high");

# Define a function to calculate the mode
Mode <- function(x) {
  ux <- unique(x)
  ux[which.max(tabulate(match(x, ux)))]
}

# Define function to test expansion/contraction
analyze_gene <- function(row) {
  
  # Calculate means for F and T groups
  ### mean_f <- mean(as.numeric(row[group_f]))
  ### mean_t <- mean(as.numeric(row[group_t]))
  # overall_mean <- mean(c(mean_f, mean_t))
  overall_mode <- Mode(c(row[group_f], row[group_t]));
  
  # Expansion analysis
  expand_f <- sum(as.numeric(row[group_f]) > overall_mode)
  expand_t <- sum(as.numeric(row[group_t]) > overall_mode)
  if ((expand_f + expand_t > 0) && ((size_f - expand_f) + (size_t - expand_t) > 0)) {
    expand_pval <- fisher.test(matrix(c(expand_f, size_f - expand_f, expand_t, size_t - expand_t), nrow = 2))$p.value
  } else {
    expand_pval <- NA;
  }
  
  # Contraction analysis
  contract_f <- sum(as.numeric(row[group_f]) < overall_mode)
  contract_t <- sum(as.numeric(row[group_t]) < overall_mode)
  if ((contract_f + contract_t > 0) && ((size_f - contract_f) + (size_t - contract_t) > 0)) {
    contract_pval <- fisher.test(matrix(c(contract_f, size_f - contract_f, contract_t, size_t - contract_t), nrow = 2))$p.value
  } else {
    contract_pval <- NA
  }
  
  # c(gene_name, mean_f, mean_t, expand_f, expand_t, expand_pval, contract_f, contract_t, contract_pval)
  c(overall_mode, expand_f, expand_t, expand_pval, contract_f, contract_t, contract_pval)
}

# Apply analysis for each gene
# results <- as.data.frame(t(apply(gene_data, 1, analyze_gene)))
# colnames(results) <- c("Gene", "Mean_F", "Mean_T", "Expand_F", "Expand_T", "Expand_Pval", "Contract_F", "Contract_T", "Contract_Pval")
results <- as.data.frame(t(apply(gene_data[,-1], 1, analyze_gene)))
colnames(results) <- c("Mode_all", "Expand_F", "Expand_T", "Expand_Pval", "Contract_F", "Contract_T", "Contract_Pval")
results <- data.frame(Gene= gene_data$OGID, results);
results <- results %>%
  mutate(Expand_Pval = as.numeric(Expand_Pval),
         Contract_Pval = as.numeric(Contract_Pval)) %>%
  mutate(Expand_FDR = p.adjust(Expand_Pval, method = "BH"),
         Contract_FDR = p.adjust(Contract_Pval, method = "BH"),
         Expansion_rawP = case_when(
           Expand_Pval < 0.05 & Expand_F/size_f > Expand_T/size_t ~ F_high,
           Expand_Pval < 0.05 & Expand_T/size_t > Expand_F/size_f ~ T_high,
           TRUE ~ "No"
         ),
         Contraction_rawP = case_when(
           Contract_Pval < 0.05 & Contract_F/size_f > Contract_T/size_t ~ T_high,
           Contract_Pval < 0.05 & Contract_T/size_t > Contract_F/size_f ~ F_high,
           TRUE ~ "No"
         ),
         Expansion = case_when(
           Expand_FDR < 0.05 & Expand_F/size_f > Expand_T/size_t ~ F_high,
           Expand_FDR < 0.05 & Expand_T/size_t > Expand_F/size_f ~ T_high,
           TRUE ~ "No"
         ),
         Contraction = case_when(
           Contract_FDR < 0.05 & Contract_F/size_f > Contract_T/size_t ~ T_high,
           Contract_FDR < 0.05 & Contract_T/size_t > Contract_F/size_f ~ F_high,
           TRUE ~ "No"
         )) %>%
  mutate(Size_F = size_f, Size_T = size_t ) %>%
  mutate(Final = case_when(
           Expansion == Contraction | Contraction == "No" ~ Expansion,
           Expansion == "No" ~ Contraction,
           Expand_FDR < 0.01 * Contract_FDR ~ Expansion,
           Contract_FDR < 0.01 * Expand_FDR ~ Contraction,
           TRUE ~ "No"
        ),
         Final_FDR = case_when(
           is.na(Expand_FDR) & is.na(Contract_FDR) ~ NA,
           is.na(Expand_FDR)   ~ Contract_FDR,
           is.na(Contract_FDR) ~ Expand_FDR,
           Expansion == Contraction & Expand_FDR < Contract_FDR  ~ Expand_FDR,
           Expansion == Contraction & Expand_FDR >= Contract_FDR ~ Contract_FDR,
           Contraction == "No" ~ Expand_FDR,
           Expansion   == "No" ~ Contract_FDR,
           Expand_FDR < Contract_FDR ~ Expand_FDR,
           Expand_FDR > Contract_FDR ~ Contract_FDR,
           TRUE ~ NA
        ))

# Select and save output
output_table <- results %>%
  select(Gene, Final, Expansion, Contraction, Expansion_rawP, Contraction_rawP, Final_FDR, Expand_FDR, Contract_FDR, Expand_Pval, Contract_Pval, Mode_all, Expand_F, Expand_T, Contract_F, Contract_T, Size_F, Size_T)

write.table(output_table, file = opt$output, sep = "\t", row.names = FALSE, quote = FALSE, col.names = TRUE)

