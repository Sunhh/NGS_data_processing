# LD_ana — linkage-disequilibrium decay

Compute pairwise LD within groups over sliding windows (via **Haploview**) and bin the
results into an LD-decay curve.

Needs `java` + the **Haploview** jar.

| script | role |
|--------|------|
| `pipe_cnt_LD.pl` | **orchestrator**: SNP table + group list → per-pair LD counts over windows (Haploview) |
| `bin_LD_cnt.pl` | bin the LD counts by distance (`-bin_len 50`) into a decay table |

---

## Workflow
```sh
perl pipe_cnt_LD.pl -snp_tbl in_snp.tbl -grp_list indiv_grp.list \
  -o_pref out_prefix -chr_kln in_chr.fa.key_len_number
perl bin_LD_cnt.pl -bin_len 50 out_prefix.LD_cnt > out_prefix.LD_cnt_bin50bp
```
`pipe_cnt_LD.pl` selects each group's individuals, splits the SNP table into windows
(`../fst/snpTbl_sepByWind.pl` style), runs Haploview per window, and tallies r^2 by distance;
`bin_LD_cnt.pl` averages into fixed-bp bins for plotting the decay curve.

_Hand-written — `gen_script_index.sh` leaves it alone._
