#!/usr/bin/env Rscript
library(optparse);

# bg_total <- 34429
# bg_slct  <- 6983
# sub_total <- 15
# sub_slct  <- 9

option_list <- list(
  make_option("--bg_total",  type="double", default=34429),
  make_option("--bg_slct",   type="double", default=6983),
  make_option("--sub_total", type="double", default=15),
  make_option("--sub_slct",  type="double", default=9)
)

opt <- parse_args(OptionParser(option_list=option_list))

bg_total  <- opt$bg_total;
bg_slct   <- opt$bg_slct;
sub_total <- opt$sub_total;
sub_slct  <- opt$sub_slct;

# Compute
mat <- matrix(
  c(
    sub_slct,
    sub_total - sub_slct,
    bg_slct - sub_slct,
    bg_total - bg_slct - (sub_total - sub_slct)
  ),
  nrow = 2,
  byrow = TRUE
)

res <- fisher.test(mat, alternative = "greater");

# Enrichment ratio
enrich_ratio <- (sub_slct / sub_total) / (bg_slct / bg_total);

cat(sprintf("P value\t%.3g\n", res$p.value));
cat(sprintf("Enrichment ratio\t%.4f\n", enrich_ratio));
cat(sprintf("Odds ratio\t%.4f\n", res$estimate));

