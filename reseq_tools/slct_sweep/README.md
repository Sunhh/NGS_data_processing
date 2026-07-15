# slct_sweep — selective-sweep windows from population diversity

Detect candidate selective-sweep regions from windowed nucleotide diversity (π): reduce
the per-window, per-group π to a between-group ratio (ROD), pick the extreme windows, join
them into regions, and annotate. The windowed π comes from `../snpTbl_stats.pl`
(pi/theta/tajima_D per window); windowed FST from `../fst/` can be used the same way.

| script | role |
|--------|------|
| `get_mean.pl` | collapse a per-window, multi-group π/θ/D table to per-group means (`.avg`) |
| `rod_from_PIavg.pl` | two groups' π → ROD = π_high / (π_high + π_low) per window |
| `slct_sweep_wind.pl` | pick windows in the extreme `-qtCutoff` percentile and group neighbours → selected windows |
| `merge_wind_pos.pl` | merge adjacent selected windows into regions (`-dist2join`) |
| `rm_overlap_wind.pl` | collapse overlapping sliding windows |
| `wind_in_region_list.pl` | keep only windows falling in a region list |
| `ret_annot_by_loc.pl` | annotate selected regions with gene annotations (`-inLocLis -inAnnotLis`) |

---

## Workflow
1. Windowed π per group (from `../snpTbl_stats.pl -value_types 'pi,theta,tajima_D'`), then
   per-group means:
```sh
perl get_mean.pl in_w100ks10k.grpAll_grp4_grp3_grp2 > in.PIavg
```
2. Between-group ROD (paste the two groups' π columns first):
```sh
paste filt_w50ks5k.pi.gCC.PIavg filt_w50ks5k.pi.gCA.PIavg | deal_table.pl -column 0-5,11 \
  | perl rod_from_PIavg.pl > CC_CA.avgComp
```
3. Call sweep windows (extreme percentile; group neighbouring windows):
```sh
perl slct_sweep_wind.pl -inRatioFile CC_CA.avgComp -qtCutoff 0.01 \
  -slct_colN 5 -bpCnt_colN 4 -grpLen 5 -grpGood 3 > CC_CA.slct
#   -qtFromLow : take the LOW tail;  -joinNeighbor : always join adjacent selected windows
```
4. Merge / de-overlap the selected windows into regions:
```sh
perl merge_wind_pos.pl CC_CA.slct -dist2join 1 > CC_CA.slct.merged
perl rm_overlap_wind.pl CC_CA.slct > CC_CA.slct.noOverlap      # alternative
```
5. Annotate (and/or restrict to a region list):
```sh
perl ret_annot_by_loc.pl -inLocLis CC_CA.slct.merged -inAnnotLis genes.gff3.annot > CC_CA.slct.annot
perl wind_in_region_list.pl CC_CA.slct.merged in.PIavg > CC_CA.slct.inRegion
```

_Hand-written — `gen_script_index.sh` leaves it alone._
