# `evolution_tools/structure` — population structure with STRUCTURE

Infer population ancestry proportions from a SNP table using Pritchard's
**STRUCTURE** (admixture model), run over many replicates and a range of *K*,
then summarise across replicates with **CLUMPP** and report per-individual
membership (Q) matrices ready for bar-plotting.

External programs: **STRUCTURE** (a 32-bit `structure` binary and its
`mainparams` / `extraparams` templates are vendored in this folder) and
**CLUMPP** (for aligning replicate runs). The pipeline scripts wrap them and
handle the SNP-table <-> STRUCTURE-format conversions, replicate bookkeeping,
and *K*-selection.

Genotype input is a repo-standard `cols` SNP table (chr, pos, one column per
individual) or a `vcf.tab`; STRUCTURE itself wants a recoded integer matrix
(alleles -> `1 2 3 4`, missing -> `-9`), which the prep scripts produce.

---

## Scripts by role

| script | role |
|--------|------|
| **filter / relabel** | |
| `03.rm_Nmiss_maf.pl` | drop SNP sites by missing-`N` rate and minor-allele frequency -> `in.cols.filter` |
| `shrt_col0.pl` | shorten long individual IDs in column 0 (STRUCTURE truncates long labels) |
| **build STRUCTURE input** | |
| `prepare_structure_input.pl` | **modern (2019) entry point** — `cols` SNP table or `vcf.tab` -> STRUCTURE input, and can emit per-*K* `mainparams` + a run command list (`-task cnvt_SNP` / `cnvt_vcfTab` / `get_param_*`) |
| `get_structure_input.pl` | legacy recode (ATGC/IUPAC -> `1 2 3 4`, `-`->`-9`); called *internally* by the driver's default (non-`-isStructTbl`) route |
| `new_mainparams.pl` | write a `mainparams` per *K* from the `mainparams` template |
| `rand_small_position.pl` | draw a random subset of SNP positions for one replicate |
| **run** | |
| `00.run_Structure.pl` | **driver** — per replicate: subsample SNPs, build `Input.file`, write per-*K* `mainparams`, and print the `structure` run commands (one per *K*) to STDOUT |
| **collect / choose K** | |
| `mv_result_files.pl` | rename scattered `structure_<rep>/structure_K<k>` result files by replicate+*K* |
| `collect_rand_result.pl` | gather per-replicate result files into one output dir |
| `get_LnPD.pl` | extract Ln P(D) (estimated log-likelihood) across runs -> pick best *K* |
| `get_time_from_structScrn.pl` | parse STRUCTURE run time from its screen log (bookkeeping) |
| **align replicates (CLUMPP) + report** | |
| `order_structureIndv_byIndID.pl` | reorder a STRUCTURE indfile's individuals to a given ID order (to line up CLUMPP inputs) |
| `order_ClumppIndFileOut_byIndID.pl` | reorder CLUMPP output rows by a given ID order |
| `cnvt_clumppOut_to_tab.pl` | CLUMPP `ClumppIndFile.output` -> tidy Q-matrix table for plotting |
| `result_report.pl` | assemble the final per-individual membership report |

---

## Workflow

### 1. Filter and relabel the SNP table
```sh
perl 03.rm_Nmiss_maf.pl all.cols          # -> all.cols.filter (by N-miss rate + MAF)
perl shrt_col0.pl indvID.long > indvID.shrtC-3   # trim long individual IDs
```

### 2. Build the STRUCTURE input
Two routes to the same recoded matrix:

- **Modern route** — convert the SNP table yourself, then feed the driver with
  `-isStructTbl` so it uses the table as-is:
  ```sh
  perl prepare_structure_input.pl -task cnvt_SNP  -inFmt SNP    -infile_SNP all.cols.filter -outfile_struct struct.tbl
  # (or)                          -task cnvt_vcfTab -inFmt vcfTab -infile_SNP all.vcf.tab -outfile_struct struct.tbl
  ```
- **Legacy route** — hand the driver a raw genotype table and let it call
  `get_structure_input.pl` internally (no `-isStructTbl`).

`prepare_structure_input.pl` also has `-task get_param_*` helpers to emit per-*K*
`mainparams` and a `structure` command list without the driver, driven by the
vendored `mainparams` / `extraparams` and the `structure` binary
(`<REPO>/evolution_tools/structure/structure`).

### 3. Run STRUCTURE over replicates x K
The driver loops replicates `minR..maxR`; for each it subsamples `-snp_number`
SNPs (`-1` = use all), builds `Input.file`, writes a `mainparams` per *K*
(`minK..maxK`), and **prints one `structure` command per K to STDOUT** — collect
them and run (optionally in parallel):
```sh
perl 00.run_Structure.pl \
 -indv_txt    indvID.shrtC-3 \
 -geno_tbl    struct.tbl -isStructTbl \
 -snp_number  1000 \
 -minR 1 -maxR 20 -minK 2 -maxK 20 \
 -OutDir      ./struct_out  > cmd_list_runStruct
bash cmd_list_runStruct        # each: cd structure_<rep>/structure_K<k>; ./structure ...
```

### 4. Collect results and choose K
```sh
perl mv_result_files.pl     result.paths        # normalise result filenames
perl collect_rand_result.pl in_dir out_dir      # gather per-replicate results
perl get_LnPD.pl            struct_result_list   # Ln P(D) per K/rep -> best K
```

### 5. Align replicates with CLUMPP, then report
Run CLUMPP on the chosen *K*'s replicate Q-matrices, then tidy and report:
```sh
perl order_structureIndv_byIndID.pl   order_list P2_var_rep_3_f    > P2_var_rep_3_f.srt   # before CLUMPP
# ... run CLUMPP -> ClumppIndFile.output ...
perl order_ClumppIndFileOut_byIndID.pl order_list ClumppIndFile.output > ClumppIndFile.output.srt
perl cnvt_clumppOut_to_tab.pl -clumpp_out ClumppIndFile.output.srt -new_order ordered_ID_list > Q.tab
perl result_report.pl result.list individual.txt final_report
```
`Q.tab` is the per-individual ancestry-proportion matrix for structure bar-plots.

---

_Two input-prep routes coexist: `prepare_structure_input.pl` (2019, general
`cols`/`vcf.tab`) feeds the driver with `-isStructTbl`, while the driver's
default path recodes via `get_structure_input.pl`; `rand_small_position.pl` and
`new_mainparams.pl` are driver-internal helpers. The `structure` binary and
`mainparams`/`extraparams` are vendored 3rd-party STRUCTURE assets. Hand-written
— not the `gen_script_index.sh` autogen index._
