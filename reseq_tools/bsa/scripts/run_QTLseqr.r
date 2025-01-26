#!/usr/bin/env Rscript

suppressMessages({ library(optparse) });

# Manually corrrect this variable to include all wanted chromosome names!!!
Chroms   <- 1:11;
# window_sizes <- c(1e6, 2e6, 500e3);
# window_size_IDs <- c('1M', '2M', '500K')

# Define command-line options
option_list <- list(
  make_option(c("-i", "--in_vcf_table"), type = "character", help = "GATK VCF-derived table with only high and low bulk samples", metavar = "FILES"),
  make_option("--high_bulk", type = "character", help = "Sample ID of high bulk"),
  make_option("--low_bulk", type = "character", help = "Sample ID of low bulk"),
  make_option("--indvN_high", type = "numeric", help = "Sequenced sample number in high bulk", default= 30),
  make_option("--indvN_low", type = "numeric", help = "Sequenced sample number in low bulk", default= 30),
  make_option("--window_size", type = "numeric", help = "Length of window size", default= 1e6),
  make_option("--window_name", type = "character", help = "Label of window size", default= "1M"),
  make_option("--plot_chr",    type= "character", help = "ID of one chromosome to be plotted"),
  make_option(c("-o", "--output_prefix"), type = "character", help = "Prefix of output files.", default= "out")
)

# Parse command-line arguments
opt_parser <- OptionParser(option_list = option_list)
opt <- parse_args(opt_parser)

if (is.null(opt$in_vcf_table) || is.null(opt$high_bulk) || is.null(opt$low_bulk)) {
  print_help(opt_parser)
  stop("\nMandatory arguments missing! Provide phenotypic data files and group definitions.\n")
}

# rawData <- 'vcf/use-2bulk-chrV.table'
# HighBulk <- 'hongchi'
# LowBulk  <- 'fenchi'

# bulk_size_High <- 20
# bulk_size_Low  <- 20


suppressMessages({
library(QTLseqr)
library(ggplot2)
library(dplyr)
});


#import data
df <- importFromGATK(
 file = opt$in_vcf_table,
 highBulk = opt$high_bulk,
 lowBulk = opt$low_bulk,
 chromList = Chroms
)

# Quality control.

### total read depth;
pdf(file= paste0(opt$output_prefix, '-depth_ttl.pdf'))
ggplot(data = df) +
 geom_histogram(aes(x = DP.HIGH + DP.LOW)) +
 xlim(0,80)
dev.off()
pdf(file= paste0(opt$output_prefix, '-depth_high.pdf'))
ggplot(data = df) +
 geom_histogram(aes(x = DP.HIGH)) +
 xlim(0,80)
dev.off()
pdf(file= paste0(opt$output_prefix, '-depth_low.pdf'))
ggplot(data = df) +
 geom_histogram(aes(x = DP.LOW)) +
 xlim(0,80)
dev.off()
### total reference allele frequency
pdf(file= paste0(opt$output_prefix, '-ref_AF_ttl.pdf'))
ggplot(data = df) + geom_histogram(aes(x = REF_FRQ))
dev.off()
### SNP-index for high bulk.
pdf(file= paste0(opt$output_prefix, '-SNPidx_high.pdf'))
ggplot(data = df) + geom_histogram(aes(x = SNPindex.HIGH))
dev.off()
### SNP-index for low bulk.
pdf(file= paste0(opt$output_prefix, '-SNPidx_low.pdf'))
ggplot(data = df) + geom_histogram(aes(x = SNPindex.LOW))
dev.off()

### Filter SNPs by depth
df_filt <- filterSNPs(
 SNPset = df,
 refAlleleFreq = 0.10,
 minTotalDepth = 5,
 maxTotalDepth = 60,
 depthDifference = 22,
 minSampleDepth = 3,
 minGQ = 30,
 verbose = TRUE
)

# Run BSA analysis.

  df_filt_ana <- df_filt # Duplicate data to avoid changing origin data.

  # Run analysis:
  ### runQTLseqAnalysis: performs Takagi et al type QTLseq analysis
  df_filt_ana <- runQTLseqAnalysis(df_filt_ana,
   windowSize = opt$window_size,
   popStruc = "F2",
   bulkSize = c(opt$indvN_high, opt$indvN_low),
   replications = 10000,
   intervals = c(95, 99),
   maxk = 200
  )
  ### runGprimeAnalysis: performs Magwene et al type G′ analysis
  df_filt_ana <- runGprimeAnalysis(df_filt_ana,
   windowSize = opt$window_size,
   outlierFilter = "deltaSNP",
   filterThreshold = 0.1,
   maxk = 200
  )
  
  pdf(file= paste0(opt$output_prefix, '-ana-', opt$window_name, '.pdf'), width= 14, height= 7)
  ### Plot G prime distribution: Is it close to log-normal distribution?
  plotGprimeDist(SNPset = df_filt_ana, outlierFilter = "Hampel") + xlim(0, 10)
  plotGprimeDist(SNPset = df_filt_ana, outlierFilter = "deltaSNP", filterThreshold = 0.1) + xlim(0,10)
  ### plot the SNP/window distribution along chromosome.
  p1 <- plotQTLStats(SNPset = df_filt_ana, var = "nSNPs")
  print(p1);
  ### plot the confidence intervals to identify QTL using G prime method.
  p2 <- plotQTLStats(SNPset = df_filt_ana, var = "Gprime", plotThreshold = TRUE, q = 0.01)
  print(p2)
  ### plot the confidence intervals to identify QTL using QTL seq method.
  p3 <- plotQTLStats(SNPset = df_filt_ana, var = "deltaSNP", plotIntervals = TRUE) + 
    theme(
      panel.grid.major = element_blank(), # Remove grids.
      panel.grid.minor = element_blank()  # Remove grids.
    )
  print(p3)
  if (!is.null(opt$plot_chr)) {
    ### Subset plot
    ##### QTLseqr
    p4 <- plotQTLStats(
     SNPset = df_filt_ana,
     var = "deltaSNP",
     plotIntervals = TRUE,
     plotThreshold = FALSE,
     q = 0.01,
     subset = opt$plot_chr
    ) +
    theme(
      panel.grid.major = element_blank(), # Remove grids.
      panel.grid.minor = element_blank()  # Remove grids.
    )
    print(p4)
    ##### Gprime
    p5 <- plotQTLStats(
     SNPset = df_filt_ana,
     var = "negLog10Pval",
     plotThreshold = TRUE,
     q = 0.01,
     subset = opt$plot_chr
    ) +
    theme(
      panel.grid.major = element_blank(), # Remove grids.
      panel.grid.minor = element_blank()  # Remove grids.
    )
    print(p5)
  }
  dev.off()
 

  # Extract QTL data.
  res_qtl_var  <- getSigRegions(SNPset = df_filt_ana, alpha = 0.01, method= 'QTLseq')
  res_qtl_var_merged <- bind_rows(res_qtl_var)
  res_qtl_loci <- getQTLTable(SNPset = df_filt_ana, method = "QTLseq", interval = 99, export = FALSE) # 6   1 23817744 25756406 1938662;
  res_gpr_var  <- getSigRegions(SNPset = df_filt_ana, alpha = 0.01, method= 'Gprime')
  res_gpr_var_merged <- bind_rows(res_gpr_var)
  res_gpr_loci <- getQTLTable(SNPset = df_filt_ana, method = "Gprime", interval = 99, export = FALSE) # 6   1 23817744 25756406 1938662;

  # Output tables.
  write.table(x= df_filt_ana,  file= paste0(opt$output_prefix, '-ana-', opt$window_name, '-data.tab'), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(x= res_qtl_loci, file= paste0(opt$output_prefix, '-ana-', opt$window_name, '-QTLseq_loci.tab'), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(x= res_qtl_var_merged, file= paste0(opt$output_prefix, '-ana-', opt$window_name, '-QTLseq_vars.tab'), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(x= res_gpr_loci, file= paste0(opt$output_prefix, '-ana-', opt$window_name, '-Gprime_loci.tab'), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)
  write.table(x= res_gpr_var_merged, file= paste0(opt$output_prefix, '-ana-', opt$window_name, '-Gprime_vars.tab'), append = FALSE, quote = FALSE, sep = "\t", row.names = FALSE, col.names = TRUE)

