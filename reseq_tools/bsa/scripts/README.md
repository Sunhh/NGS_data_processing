# bsa/scripts — BSA analysis scripts

Two BSA analyses run on the two-bulk VCF produced in the parent [`../README.md`](../README.md)
(GATK → filtered `two_bulks.vcf.gz` → GATK `VariantsToTable`):

- **QTLseqr** (`run_QTLseqr.r`) — the path documented in the parent README (ΔSNP-index and
  G' with the R/QTLseqr package).
- **Custom G' pipeline** (the Perl + R scripts here) — a from-scratch ΔSNP-index / G' scan,
  useful when a parent genotype is available or QTLseqr does not fit.

| script | role |
|--------|------|
| `run_QTLseqr.r` | QTLseqr BSA (`-i table --high_bulk --low_bulk --indvN_high --indvN_low --window_size`) |
| `cnvt_var2tab.pl` | `VariantsToTable` output → `bsaTab` (allele-depth per bulk; `-highID -lowID`) |
| `slct_sites_forBsa.pl` | VCF → per-site ΔSNP-index with parent constraints (`-P1ID -P2ID -B1ID`) |
| `filter_vcfTab_forBsaParent.pl` | filter a vcfTab by per-column rules (`-col_noN -col_homo -col_diffPair -col_samePair`) |
| `ana_bsa_Gprime.pl` | `bsaTab` → windowed G' statistic (`-windSize`) |
| `get_pval_for_Gprime.R` | windowed G' → p-values / FDR (null distribution) |
| `slct_sites_by_windows.pl` | keep sites falling in a window/position list |
| `plot_pipeResult.R` | plot a windowed track (SNP-index / G' / p) per chromosome; optional min-count mask |

---

## QTLseqr path (see parent README for the VCF→table steps)
```sh
Rscript run_QTLseqr.r -i vcf/use-2bulk-chrV.table --high_bulk highBulk --low_bulk lowBulk \
  --indvN_high 20 --indvN_low 20 --window_size 2000000 --window_name 2M --plot_chr 6 -o output
```

## Custom G' pipeline
```sh
# 1. table -> bsaTab (allele depths for the high/low bulk)
perl cnvt_var2tab.pl input.VariantsToTable.table -highID F2BMut -lowID F2H > input.bsaTab
#    (with a parent, slct_sites_forBsa.pl -P1ID -P2ID -B1ID goes straight from the VCF to a SNP-index)
# 2. windowed G'
perl ana_bsa_Gprime.pl input.bsaTab -windSize 100000 > input.bsaGprime
# 3. p-values for G'
Rscript get_pval_for_Gprime.R input.bsaGprime > input.bsaGprime.pval
# 4. plot a column (e.g. Gprime) per chromosome
Rscript plot_pipeResult.R input.bsaGprime.pval chrList Gprime out.pdf point 10
#    args: <table> <chrList> <column> <out.pdf> [point|line|both=point] [minCnt=0]
```

_Hand-written — `gen_script_index.sh` leaves it alone._
