#!/usr/bin/env Rscript
###!/home/sunhh/bin/R
# A script to perform rank inverse-normal transformation (INT) to deal with extremely biased phenotypic data.
## Could be useful for resistance data in natural populations.

args <- commandArgs(trailingOnly = TRUE)

if (length(args) < 2) {
  stop("Usage: Perform rank inverse-normal transformation.\n      Rscript this.r input_file output_file\n");
}

infile  <- as.character(args[1])
outfile <- as.character(args[2])
# infile  <- 'pheno-race1rmWild'
# outfile <- 'pheno-race1rmWildINT'

ph <- read.table(infile, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE)

# Rank-based inverse normal transformation (Blom's formula: (r-0.5)/n)
out_ph <- ph
for (i in 2:ncol(out_ph)) {
  # i <- 2
  ok <- is.finite(ph[,i])
  r <- rank(ph[ok,i], ties.method = "average")
  out_ph[,i] <- NA_real_
  out_ph[ok,i] <- qnorm( (r - 0.5) / length(r) )
  colnames(out_ph)[i] <- paste0(colnames(ph)[i],"INT")
}

# Write output (keep ID and INT phenotype)
write.table(out_ph, outfile, sep= "\t", quote= FALSE, row.names= FALSE)


