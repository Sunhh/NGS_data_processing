# pca — population-structure PCA (for GWAS covariates / structure plots)

Turn a SNP table into PCA input (EIGENSOFT/`smartpca` style), then tabulate and plot the
eigenvectors. The top PCs are used as covariates in [`../gwas_tools/`](../gwas_tools/) and
to visualise population structure.

| script | role |
|--------|------|
| `cnvt_snp2pca.pl` | SNP table → PCA input files (`-opref out`) |
| `cnvt_pcaEvec_to_tbl.pl` | `smartpca` `*.pca.evec` → a plain table (sample × PC) |
| `plot_EVs.R` | scatter the first eigenvectors, coloured by group |

---

## Workflow
```sh
perl cnvt_snp2pca.pl in.snp -opref set01           # -> set01.* PCA input
smartpca -p set01.par > set01.pca.evec             # external EIGENSOFT run
perl cnvt_pcaEvec_to_tbl.pl set01.pca.evec > set01.pca.evec.tbl
Rscript plot_EVs.R set01.pca.evec.tbl              # edit group colours at the top
```
`set01.pca.evec.tbl` (top PCs) feeds `../gwas_tools/` as structure covariates. See
`cmd_list_set01` for a fuller worked example (note: it hardcodes some `../../structure/`
paths from the original project).

_Hand-written — `gen_script_index.sh` leaves it alone._
