# gwas_tools — GWAS with FaST-LMM (prep, peaks, LD, plot)

Helpers around a **FaST-LMM** mixed-model GWAS: prepare genotype/phenotype, derive a
significance threshold, pull peaks, define their LD blocks, and plot. Population-structure
covariates come from [`../pca/`](../pca/); LD-decay analysis is in [`../LD_ana/`](../LD_ana/);
`../tassel/` has a taxa-list helper for the Tassel alternative.

Needs the external **FaST-LMM** and **plink** binaries, plus `R` (`Rscript`).

| script | role |
|--------|------|
| `cnvt_sv2snpVCF.pl` | recode structural variants as bi-allelic SNP-style VCF so SVs can be tested in the SNP GWAS |
| `perform_INT_pheno.r` | rank inverse-normal transform (INT) of a skewed phenotype |
| `parse_gec.pl` | GEC output (`snp_gec.sum`) → significance cutoffs (`.cuts`, effective-test corrected) |
| `get_gwasPeaks_fastlmm.pl` | call peaks from a FaST-LMM result by **distance** (`minuslogP_cut res peak_distance`) |
| `get_gwasPeaks_fastlmm_byLD.pl` | call peaks by **LD block** (uses `peak2LD_plink_fastlmm.pl`; joins nearby blocks) |
| `peak2LD_plink_fastlmm.pl` | define the plink LD block around each peak (`outPref peak vcf.gz minR2 flank`) |
| `plot_gwas_fastlmm.r` | Manhattan + QQ plot of a FaST-LMM result |

---

## Workflow
1. **Genotype**: a SNP table/VCF; add SVs with `cnvt_sv2snpVCF.pl` if testing them.
2. **Phenotype**: for skewed traits, `Rscript perform_INT_pheno.r pheno.tsv > pheno.int.tsv`.
3. **Covariates**: population-structure PCs from `../pca/` (top eigenvectors).
4. **Association**: run FaST-LMM (external) with the genotype, INT phenotype and PC covariates
   → `in_fastlmm.res`.
5. **Threshold**: run GEC (external) for the effective number of tests, then
   `perl parse_gec.pl snp_gec.sum > snp_gec.sum.cuts`.
6. **Peaks** (pick one):
```sh
perl get_gwasPeaks_fastlmm.pl        snp_gec.sum.cuts in_fastlmm.res 1000000 > res.peak
perl get_gwasPeaks_fastlmm_byLD.pl   out res.cuts in_fastlmm.res geno.vcf.gz 1e6 1e5   # LD-block peaks
```
   (both accept either a `-log10P` number or a `.cuts` file as the cutoff.)
7. **Peak LD** (for the distance caller): `perl peak2LD_plink_fastlmm.pl out res.peak geno.vcf.gz 0.8 1e6`
8. **Plot**: `Rscript plot_gwas_fastlmm.r in_fastlmm.res out_prefix`

_Hand-written — `gen_script_index.sh` leaves it alone._
