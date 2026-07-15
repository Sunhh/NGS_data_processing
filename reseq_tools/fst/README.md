# fst — population differentiation (FST) from a SNP table

Compute windowed FST between groups with `hierfstat`, starting from a genotype table
(`cols` format) and an individual→group mapping. `pipe_get_fst.pl` orchestrates the whole
run; the other scripts are its steps or downstream helpers.

Needs `R` with the **hierfstat** package (and `Rscript` on PATH, or `-exe_Rscript`).

| script | role |
|--------|------|
| `pipe_get_fst.pl` | **orchestrator:** SNP table + ind→group → per-site & per-window FST |
| `snpTbl_sepByWind.pl` | split the SNP table into sliding windows | 
| `cols2fstat.pl` (in `../cnvt_tools/`) | per-window genotype table → hierfstat input |
| `run_hierfstat.pl` | run `hierfstat` on a window (Nei / Weir-Cockerham FST); drives a temp R script |
| `join_fst_siteChrPos.pl` | merge the per-window results into per-site / per-window tables |
| `get_stat.pl` | quick summary stats of a windowed-FST file |
| `extract_top1_fst_site.R` | pull the top-FST sites from `*.fst.perSiteChrPos` |
| `extract_top1_vcfFst_wind.R` | pull the top-FST windows from a **vcftools** `*.windowed.weir.fst` (alternative FST source) |

---

## Inputs
- **`-snp_tbl`**: `cols`-format genotype table — row1 header, col0 = chr, col1 = pos, then
  one column per individual (e.g. produced by `../cnvt_tools/tab2cols.pl` from a VCF table).
- **`-ind2grp_list`**: two columns, `individual_ID`  `group_number`; individuals not listed
  are dropped, and each site is compared between the listed groups.

## Run
```sh
perl pipe_get_fst.pl \
  -snp_tbl in_snp.tbl \
  -ind2grp_list indiv_to_grpNum \
  -o_pref out_prefix \
  -wind_length 10000 -wind_step 10000 -ncpu 20
# optional site filters: -maxNmissR 0.5   -rmNegNeiFst   -rmNegWcFst
```
Steps run for you: `snpTbl_sepByWind.pl` (windows) → per window `../cnvt_tools/cols2fstat.pl`
+ `run_hierfstat.pl` → `join_fst_siteChrPos.pl`. Outputs:
- `out_prefix.fst.perSiteChrPos` — per-site FST (chr, pos, Ho, Hs, Ht, Dst, …, Fst, Fis, Dest)
- `out_prefix.fst.perWindLine` — per-window FST.

## Downstream
```sh
perl get_stat.pl out_prefix.fst.perWindLine            # summary stats
Rscript extract_top1_fst_site.R out_prefix.fst.perSiteChrPos   # top-FST sites
```
The per-window FST feeds the selective-sweep step (`../slct_sweep/`).

_Hand-written — `gen_script_index.sh` leaves it alone._
