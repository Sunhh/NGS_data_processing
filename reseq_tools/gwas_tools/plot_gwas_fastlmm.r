#!/usr/bin/env Rscript
#r Guide 1: https://r-graph-gallery.com/101_Manhattan_plot.html
# Guide 2: https://cran.r-project.org/web/packages/qqman/vignettes/qqman.htmul
# 20250122: Change shapes of SVs.

suppressMessages({ library(optparse) });
# Define command-line options
option_list <- list(
  make_option("--fastlmm_res", type = "character", help = "fl.var.res; Result file from fastlmm", metavar = "FILES"),
  make_option("--pvalue_threshold_file", type = "character", help = "var.cuts file with cut1_signif and cut2_suggest in -log10;", metavar = "FILES"),
  make_option("--cut1_signif", type = "numeric", help = "cut1_signif value for alpha=0.05"),
  make_option("--cut2_suggest", type = "numeric", help = "cut2_suggest value for alpha=0.1"),
  make_option("--region", type = "character", help = "Chromosome region in the format '10:100-3000'"),
  make_option(c("-o", "--output_prefix"), default= "output", type = "character", help = "Prefix of output files.")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser);


# Check for mandatory arguments.
if (is.null(opt$fastlmm_res)) {
  print_help(opt_parser);
  stop("\nNeed --fastlmm_res\n");
}

tgt_chr <- NA;
tgt_s   <- NA;
tgt_e   <- NA;
if (!is.null(opt$region)) {
  match <- regexec("^([0-9]+):(\\d+)-(\\d+)$", opt$region)
  parsed <- regmatches(opt$region, match)
  if (length(parsed[[1]]) == 4) {
    # Extract chromosome, start, and end
    tgt_chr <- as.numeric(parsed[[1]][2])
    tgt_s   <- as.numeric(parsed[[1]][3])
    tgt_e   <- as.numeric(parsed[[1]][4])
  }
}
a1 <- read.table(opt$fastlmm_res, header=F, stringsAsFactors=F, sep= "\t", skip=1);
names(a1)[1:6] <- c('SNP', 'Chromosome', 'GeneticDistance', 'Position', 'Pvalue', 'Qvalue');
n1 <- sum(!is.na(a1$Pvalue))
if (!is.null(opt$pvalue_threshold_file)) {
  a2 <- read.table(opt$pvalue_threshold_file, header=T, stringsAsFactors=F, sep= "\t");
  cut1Red <- as.numeric(a2$cut1_signif);
  cut2Blue <- as.numeric(a2$cut2_suggest);
} else if (!is.null(opt$cut1_signif) || !is.null(opt$cut2_suggest)) {
  if (is.null(opt$cut2_suggest)) {
    stop("\n  --cut1_signif and --cut2_suggest need to be assigned together.\n");
  }
  cut1Red <- opt$cut1_signif
  cut2Blue <- opt$cut2_suggest
} else {
  pv  <- a1$Pvalue[!is.na(a1$Pvalue)];
  fdr <- p.adjust(pv, method= 'BH');
  if (sum(fdr <= 0.05) > 0) {
    cut1Red <- -log10(max(pv[fdr <= 0.01]));
  } else {
    cut1Red <- -log10(0.05/n1);
  }
  if (sum(fdr <= 0.1) > 0) {
    cut2Blue <- -log10(max(pv[fdr <= 0.1]));
  } else {
    cut2Blue <- -log10(0.1/n1);
  }
}
if (!is.na(tgt_chr)) {
  a1 <- a1[a1$Chromosome == tgt_chr & a1$Position >= tgt_s & a1$Position <= tgt_e, ]
}

# Calculate lambda: https://genometoolbox.blogspot.com/2014/08/how-to-calculate-genomic-inflation.html
### Explain lambda: https://parkinsonsroadmap.org/understanding-gwas-part-2-additional-insights-and-tips/#
### p_chi_sq <- -2 * log(a1$V6) # This is suggested by ChatGPT and it seems wrong!!!
if (is.na(tgt_chr)) {
  p_chi_sq <- qchisq(1-a1$Pvalue, 1)
  obs_median_chi_sq <- median(p_chi_sq)
  gwas_lambda <- obs_median_chi_sq / qchisq(0.5,1)
  write(gwas_lambda, file=paste0(opt$output_prefix, '.lambda'))
}


# k <- a1$V1 == 2
# plot(a1$V2[k]/1e6, -log10(a1$V6[k]), xlab= "Position (Mb)", ylab= "-log10(p-value)")
suppressMessages({
library(qqman)
library(showtext)
library(extrafont)
showtext_auto(FALSE)

});

SNP_data <- a1;
SNP_data$Type <- ifelse(grepl('v', SNP_data$SNP), 'sv', 'snp');
SNP_data <- SNP_data[order(SNP_data$Chromosome, SNP_data$Position), ];
SNP_data$cum_pos <- NA;
chrom_start <- 0;
unique_chromosomes <- unique(SNP_data$Chromosome);
for (chr in unique_chromosomes) {
    chr_indices <- which(SNP_data$Chromosome == chr)
    SNP_data$cum_pos[chr_indices] <- SNP_data$Position[chr_indices] + chrom_start
    chrom_start <- max(SNP_data$cum_pos[chr_indices]) + 1e6  # Add buffer for separation
}

# manhattan plot of P value.
# pdf(file= paste0(opt$output_prefix, ".P.pdf"), width=14, height=7, useDingbats = FALSE, family = "Arial")
cairo_pdf(file= paste0(opt$output_prefix, ".P.pdf"), width=14, height=7, family = "Arial")
# manhattan(a1, chr="V2", bp="V4", snp="V1", p="V5",
#   annotatePval = NULL, annotateTop= F,
#   col = c("blue4", "orange3"), cex= 0.6, 
#   suggestiveline=cut2Blue, genomewideline=cut1Red)
SNP_snps <- subset(SNP_data, Type == 'snp');
SNP_svs  <- subset(SNP_data, Type == 'sv');
chr_col_snp <- '#909295'
chr_col_sv  <- 'red';

if (nrow(SNP_snps) == 0) {
  if (length(unique(SNP_svs$Chromosome)) > 1) {
    chr_col_sv <- ifelse(SNP_svs$Chromosome %% 2 == 0, "blue", "red");
  }
  plot(
    SNP_svs$cum_pos, -log10(SNP_svs$Pvalue), 
    col = chr_col_sv,
    pch = 2, cex = 2, xlab = "Chromosome", ylab = "-log10(P)",
    ylim = c(0, max(-log10(SNP_data$Pvalue)) + 1), xaxt = 'n',
    bty = "l"
  );
} else {
  if (length(unique(SNP_snps$Chromosome)) > 1) {
    chr_col_snp <- ifelse(SNP_snps$Chromosome %% 2 == 0, "#909295", "#dcdddf");
  }
  if (length(unique(SNP_svs$Chromosome)) > 1) {
    chr_col_sv <- ifelse(SNP_svs$Chromosome %% 2 == 0, "blue", "red");
  }
  plot(
    SNP_snps$cum_pos, -log10(SNP_snps$Pvalue), 
    col = chr_col_snp,
    pch = 20, cex = 2, xlab = "Chromosome", ylab = "-log10(P)",
    ylim = c(0, max(-log10(SNP_data$Pvalue)) + 1), xaxt = 'n',
    bty = "l"
  );
  if (nrow(SNP_svs) > 0) {
    # Add SVs on top of SNPs, alternating colors by chromosome
    points(
      SNP_svs$cum_pos, -log10(SNP_svs$Pvalue), 
      col = chr_col_sv,
      pch = 2, cex = 2
    );
  }
}

if (is.na(tgt_chr)) {
  # Customize x-axis to display chromosome labels at the midpoint of each chromosome's range
  axis(1, at = tapply(SNP_data$cum_pos, SNP_data$Chromosome, mean), labels = unique(SNP_data$Chromosome))
  # title("Manhattan Plot")
} else {
  legend('topright', legend=c("SNP", "SV"), pch= c(20, 2))
  axis(1)
}
abline(h= cut1Red,  col= 'red',    lty= 2, lwd=3 );
abline(h= cut2Blue, col= 'green',  lty= 2, lwd=3 );
dev.off()

# Plot manhattan frame with blank plotting area.
# pdf(file= paste0(opt$output_prefix, ".P-blank.pdf"), width=14, height=7, useDingbats = FALSE, family = "Arial")
cairo_pdf(file= paste0(opt$output_prefix, ".P-blank.pdf"), width=14, height=7, family = "Arial")
if (nrow(SNP_snps) == 0) {
  plot(
    SNP_svs$cum_pos, -log10(SNP_svs$Pvalue), 
    type = 'n',
    xlab = "Chromosome", ylab = "-log10(P)",
    ylim = c(0, max(-log10(SNP_data$Pvalue)) + 1), xaxt = 'n',
    bty = "l"
  );
} else {
  plot(
    SNP_snps$cum_pos, -log10(SNP_snps$Pvalue), 
    type = 'n',
    xlab = "Chromosome", ylab = "-log10(P)",
    ylim = c(0, max(-log10(SNP_data$Pvalue)) + 1), xaxt = 'n',
    bty = "l"
  );
}
if (is.na(tgt_chr)) {
  # Customize x-axis to display chromosome labels at the midpoint of each chromosome's range
  axis(1, at = tapply(SNP_data$cum_pos, SNP_data$Chromosome, mean), labels = unique(SNP_data$Chromosome))
  # title("Manhattan Plot")
} else {
  legend('topright', legend=c("SNP", "SV"), pch= c(20, 2))
  axis(1)
}
dev.off()


# manhattan plot of Z-score.
#  png(file= paste0(opt$output_prefix, ".Z.png"), width=960, height=480)
#  a1z <- transform(a1, zscore= qnorm(V6/2, lower.tail= FALSE) )
#  manhattan(a1z, p= 'zscore', logp=F, ylab= "Z-score", chr="V1", bp="V2", snp="V3",
#    col = c("blue4", "orange3"), cex= 0.6)
#  dev.off()


# qq-plot
# png(file= paste0(opt$output_prefix, ".qq.png"))
if (is.na(tgt_chr)) {
  # pdf(file= paste0(opt$output_prefix, ".qq.pdf"), family = "Arial", useDingbats = FALSE)
  cairo_pdf(file= paste0(opt$output_prefix, ".qq.pdf"), family = "Arial")
  qq(a1$Pvalue)
  dev.off()
  # Plot a blank QQ plot.
  # pdf(file= paste0(opt$output_prefix, ".qq-blank.pdf"), family = "Arial", useDingbats = FALSE)
  cairo_pdf(file= paste0(opt$output_prefix, ".qq-blank.pdf"), family = "Arial")
  qq(a1$Pvalue, type= 'n')
  dev.off()
}

