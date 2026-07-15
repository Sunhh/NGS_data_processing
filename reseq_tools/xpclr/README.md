# xpclr — XP-CLR selection scan between two populations

Run XP-CLR from a genetic-map-annotated SNP table, window the scores, then reuse the
selective-sweep window caller and annotator from [`../slct_sweep/`](../slct_sweep/). Needs
the external **`XPCLR`** binary on PATH and `perl`.

| script | role |
|--------|------|
| `pipe_xpclr_fromSNPtbl.pl` | **orchestrator** (parallel, `-cpuN`): prepare input → XPCLR per chr → window → pick + annotate sweep regions |
| `pipe_xpclr_forWM97.pl` | legacy, hardcoded WM97 run (fixed window/params, stops at the windowed XPCLR table) — kept for reference |
| `set_GMpos_to_SNP.pl` | add a genetic-map position (cM) column to a SNP table (`marker_loc2GM  in.snp`) |
| `get_uniq_cM.pl` | collapse a marker→cM map to unique cM positions |
| `prepare_xpclr_input_wiGmP.pl` | build XPCLR per-population geno/snp input files (with cM) |
| `sep_run_xpclr.pl` | run `XPCLR` on one chromosome, splitting large inputs |
| `xpclr_wind_cmd_wiGmP.pl` | emit per-window XPCLR commands |
| `cluster_xpclrscore.pl` | window/average per-SNP XPCLR scores (`-wind_size -wind_step`) |
| `cumPos.pl` | cumulative genome position (for Manhattan plotting) |
| `chk_nonsyn.pl` | flag non-synonymous variants in candidate regions |
| `plot_manhattan/plot_xpclr_cum.R` | Manhattan plot of windowed XP-CLR |

---

## Inputs
- **set file** (`set01_Grp14_to_Grp13`): `sample_ID <TAB> group_ID`, **reference population
  first, object (selected) population after** (or use `-firstAsObjPop`).
- **`-in_snpTbl`**: SNP table with a genetic-map column — `chr, pos, cM, <per-sample genotype…>`.
  Add the `cM` column with `set_GMpos_to_SNP.pl marker_loc2GM in.snp` first if needed.
- **`-in_wind`**: window list (`chrID, WindS, WindE, WindL, BpCnt`).
- **`-in_annot`**: gene annotation (`chrID, start, end, strand, …`).

## Run
```sh
perl pipe_xpclr_fromSNPtbl.pl set01_Grp14_to_Grp13 out_dir \
  -in_snpTbl Acc131_mask.snp_addGmP -in_wind WM97_w10ks10k.wind \
  -in_annot genes.annot -lis_chrID2num chrID2num -cpuN 20
#  -firstAsObjPop : use the 1st population as the object/selected pop
#  -getCmdXPCLR   : print the XPCLR commands instead of running them
#  -chk_scripts   : only check helper scripts are found
```
Steps: `prepare_xpclr_input_wiGmP.pl` → per-chr `sep_run_xpclr.pl` (XPCLR) →
`cluster_xpclrscore.pl` (windowed scores) → `ColLink.pl` → then the shared
`../slct_sweep/slct_sweep_wind.pl` (top windows) and `../slct_sweep/ret_annot_by_loc.pl`
(annotation). Output in `out_dir/`: `xpclr_<windTag>`, `wind.xpclr.compare`,
`…med01_slct01`, `…med01_slct01.annot`.

## Plot
```sh
perl cumPos.pl chrLen xpclr_w10ks10k > xpclr_w10ks10k.jnChr
Rscript plot_manhattan/plot_xpclr_cum.R   # edit the paths at the top for your run
```

_Hand-written — `gen_script_index.sh` leaves it alone._
