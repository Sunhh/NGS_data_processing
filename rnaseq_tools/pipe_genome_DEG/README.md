# Standard genome RNA-seq → DEG pipeline (self-contained)

Map clean RNA-seq reads to a **single reference genome**, quantify per-gene
expression as TPM, and call differentially expressed genes (DEGs).

**Design of the DEG step:** DESeq2 (or edgeR) is used **only to get FDRs** — no
TPM/fold-change enters the statistical model. The actual up / down / not-changed
**DEG call is made downstream in `DEG_byList_vClaude.pl`**, which merges per-group
mean TPM with those FDRs and applies **adjustable thresholds** (FDR default `0.05`,
TPM-mean fold change default `2`). This keeps the significance test and the
effect-size/threshold decision cleanly separated and tunable without re-running R.

This directory holds every script the workflow needs. The pipeline-specific scripts are
real files here; the general helpers (`DEG_byList*.pl`, `cnvt_cnt_to_normExpr.pl`, `get_meanTPM.pl`, `plot_*`,
`rmRRNA_in_fqFiles.pl`, `summary_ht2_log.pl`) are **symlinks** to their canonical
top-level `rnaseq_tools/` copies (one source of truth). External programs still
required: `hisat2` / `hisat2-build`, `samtools`, `featureCounts` (subread), `gffread`,
`R` with **DESeq2** (and **edgeR** if you use that method), and (step 0 only) `sortmerna`.

Scripts in this folder:

| script | role | step |
|--------|------|------|
| `rmRRNA_in_fqFiles.pl` | drop rRNA reads (SortMeRNA) | 0 (optional QC) |
| `runHisat2_with2pass.pl` (calls `fix_NHnum.pl`) | 2-pass HISAT2 mapping → `*_fixNH.bam` | 1 |
| `summary_ht2_log.pl` | tabulate HISAT2 log stats | 1 (optional) |
| `cnvt_cnt_to_normExpr.pl` | counts → **TPM and/or RPKM**; reads featureCounts output or a plain matrix+length; optional true library size | 2 |
| `DEG_byList.pl` | **FDR producer:** compareList (DESeq2/edgeR) → one `ds.*`/`et.*` FDR table; `-writeRcode` dumps the R and exits | 3 |
| `DEG_byList_vClaude.pl` | **DEG caller:** TPM + FDR table → up/down/not, adjustable `--fdr_cut`/`--fc` | 4 |
| `get_meanTPM.pl` | per-replicate TPM → per-group mean | optional |
| `cnvt_gene2group_val.pl`, `cnvt_gene2group_DEGlabel.pl` | roll gene values / labels up to orthogroups | 5 (optional) |
| `plot_heatmap_by_geneList.pl`, `plot_expr_heatmap.r` | expression heatmap for a gene list | optional |

---

## Try it

A tiny runnable dataset lives in [`example/`](example/): a featureCounts-style count
matrix plus `list.comparison` / `list.in_rd` / `synOG.grp`. `bash example/run_example.sh`
runs steps 2–5 (needs `Rscript` with DESeq2 + edgeR) and writes results to `example/_out/`.

---

## (0) Optional QC — remove rRNA reads
```sh
perl rmRRNA_in_fqFiles.pl -inFq in_R1.fq,in_R2.fq -outPref cleaned -db_sortmerna db.fasta,db.idx
```

## (1) Map reads with two-step HISAT2 (no GFF used)
- Build the HISAT2 database.
```sh
hisat2-build in_chr.fa ref_db
```
- Prepare an input read list `list.in_rd` (dataPrefix and Fq names must be unique; Fq may be gz).

| #sample | readGroup | library | dataPrefix | inFq1 | inFq2 | PL | PU | Others |
|---------|-----------|---------|------------|-------|-------|----|----|--------|
| S1      | S1\_Rep1  |S1\_Rep1 | S1\_Rep1   | path1 | NA    | NA | NA | NA     |
| S1      | S1\_Rep2  |S1\_Rep2 | S1\_Rep2   | path1 | path2 | NA | NA | NA     |

- Map. Result: `*_fixNH.bam` in `out_bams/`. ` --rna-strandness R ` fits reverse-strand strand-specific RNA-seq.
```sh
mkdir out_bams/
perl runHisat2_with2pass.pl \
 -cpuN         20 \
 -in_pref_list list.in_rd \
 -db_hisat2    ref_db \
 -wrk_dir      out_bams/ \
 -para_hisat2  ' -p 4 --dta --dta-cufflinks -q --phred33 --rna-strandness R '
```
- (optional) Summarise mapping rates from the HISAT2 logs.
```sh
perl summary_ht2_log.pl out_prefix out_bams/LOG/*.ht2.log > mapping_summary.tbl
```

## (2) Count on-gene reads with featureCounts → TPM
- GTF from GFF3.
```sh
gffread ref_protein.gff3 -T -o ref_protein.gtf
```
- Count (`-Q 0 -M` keeps multi-mapping reads, good for duplicated genes; `-s 2` = reverse strand).

All samples in `list.in_rd` into one matrix `out_cnt/joint-cnt.txt`, then rename columns:
```sh
mkdir -p out_cnt/sep/
cut -f 4 list.in_rd|tail -n +2|perl -e 'my @f=<>;chomp(@f);print join(" ", "featureCounts -Q 0 -M -T 8 -s 2 -a ref_protein.gtf -o out_cnt/joint-cnt.txt", map {"out_bams/${_}_fixNH.bam"} @f)."\n";' > cx2cnt
bash cx2cnt
perl -e 'while (<>) {s!out_bams/!!g; s!_fixNH.bam!!g; print;}' out_cnt/joint-cnt.txt > out_cnt/joint-cnt_rename.txt
```
Counts → TPM with the unified converter. It auto-detects the featureCounts format and
reads gene length from the `Length` column; output columns are the per-replicate
`<group>_Rep<N>` TPMs (each column sums to 1e6):
```sh
perl cnvt_cnt_to_normExpr.pl --cntFn out_cnt/joint-cnt_rename.txt --out_type tpm > out_cnt/joint-tpm.txt
```
The same script also produces RPKM and accepts other input shapes:
- `--out_type rpkm` (or `both` with `--rpkm_out rpkm.tsv`). RPKM needs a true library
  size; by default it uses each column's count sum, but with featureCounts `-M` that
  over-counts multi-mapping reads, so pass `--libFn <sample TAB total_mapped_reads>`
  for an accurate RPKM denominator. (TPM self-normalizes and ignores `--libFn`.)
- `--fmt matrix --lenFn <geneID TAB length>` for a plain count matrix that lacks the
  featureCounts annotation columns (a 3-column `id TAB trans_len TAB length` file also works).

## (3) Get FDRs (DESeq2/edgeR — significance only, no TPM)

The comparison groups are just the replicate column names minus their trailing
`_Rep<N>` (case-insensitive), e.g. `S1G1T1P10D_Rep1` → group `S1G1T1P10D`.

`DEG_byList.pl` handles many comparisons in one call, DESeq2 or edgeR per line.
- `DEG_byList.pl` needs a **plain count matrix**: `GeneID` + one `<group>_Rep<N>`
  column per sample (drop featureCounts' Chr/Start/End/Strand/Length):
```sh
cut -f1,7- out_cnt/joint-cnt_rename.txt > out_cnt/joint-cnt_matrix.txt
```
- Write `list.comparison` — **col1 = method**, **col2 = group1**, **col3 = group2**,
  **col4 = comparison name**. Groups may pool several names with `;`. This same file is
  reused as vClaude's `--pair` in step 4. The **col1 method keyword** picks the test:

  | col1 keyword | test | FDR column |
  |--------------|------|------------|
  | `DESeq2` | DESeq2 | `ds.<g1>_VS_<g2>` |
  | `edgeR` / `edgeR_glm` / `glm` | edgeR **GLM** (glmFit + glmLRT) — *default edgeR route* | `eg.<g1>_VS_<g2>` |
  | `exactTest` / `edgeR_classic` | edgeR **classic** exact test | `et.<g1>_VS_<g2>` |

| method | group1 | group2 | name |
|--------|--------|--------|------|
| DESeq2 | S1G1T1P10D | S1G1T1P18D | flesh_10D_18D |
| edgeR  | S1G1T1P18D | S1G1T1P26D | flesh_18D_26D |

```sh
perl DEG_byList.pl \
 -rdCntFn     out_cnt/joint-cnt_matrix.txt \
 -compareList list.comparison \
 -outFDRFn    out_cnt/joint-fdr.txt
```
This writes **one** FDR table: `GeneID` plus one column per comparison, named by method
(`ds.` DESeq2, `eg.` edgeR GLM, `et.` edgeR classic) — e.g. `ds.<group1>_VS_<group2>`.
Step 4 (vClaude) matches each pair's column regardless of the method prefix. (Add
`-exe_Rscript /path/to/Rscript` if `Rscript` is not on PATH.)

To inspect or run the exact DESeq2/edgeR recipe on its own (say, in another R
environment for tuning), add `-writeRcode <dir>`: `DEG_byList.pl` writes each
comparison's R script plus its count subset under `<dir>/<comparison>/` and exits
without running R (so `-outFDRFn` is not required). That is why no separate standalone
DESeq2/edgeR R script is shipped — the R code lives in `DEG_byList.pl` and is
extracted on demand:
```sh
perl DEG_byList.pl -rdCntFn out_cnt/joint-cnt_matrix.txt -compareList list.comparison -writeRcode Rcode_out/
```

## (4) Call DEGs in vClaude (adjustable thresholds)

`DEG_byList_vClaude.pl` merges per-group mean TPM with the FDR table and labels each
gene per comparison. Reuse `list.comparison` from step 3 as `--pair` (its col2/col3 are
group1/group2; the method column is ignored here).
```sh
perl DEG_byList_vClaude.pl \
 --tpm     out_cnt/joint-tpm.txt \
 --pair    list.comparison \
 --fdr     out_cnt/joint-fdr.txt \
 --fdr_cut 0.05 \
 --fc      2 \
 --out     DEG_calls.tsv
# --add_repTPM also appends each replicate's TPM after every comparison block
```
Output — per comparison `<name>`: `<name>.v1` (mean TPM group1), `<name>.v2` (mean TPM
group2), `<name>.R` (v2/v1, each mean floored to `--floor`, default 0.01), `<name>.FDR`,
and `<name>.DEG`:
- `up`   if `R > fc`   and `FDR < fdr_cut`
- `down` if `R < 1/fc` and `FDR < fdr_cut`
- `not`  otherwise (FDR is NA, `FDR >= fdr_cut`, or `1/fc <= R <= fc`)

Because vClaude emits every comparison's call in one combined table, no separate
merge + label step is needed.

## (5) Optional — roll up to orthogroups / heatmaps
```sh
perl get_meanTPM.pl out_cnt/joint-tpm.txt > joint-tpm.mean.tsv          # per-group mean TPM
perl cnvt_gene2group_val.pl synOG.grp joint-tpm.mean.tsv > og-tpm.tsv   # mean TPM per orthogroup
perl plot_heatmap_by_geneList.pl out_cnt/joint-tpm.txt my_genes.list out_heatmap
```
To roll the **DEG calls** up to orthogroups, feed the `.DEG` columns of `DEG_calls.tsv`
(a GeneID + up/down/not label table) to `cnvt_gene2group_DEGlabel.pl`; each orthogroup gets
its non-`not` labels (both `not` and the legacy `N` count as unchanged). See
`example/run_example.sh` for the exact extraction.

---

_Reorganised, self-contained standard pipeline. DESeq2/edgeR produce FDRs only; the DEG
decision (with adjustable FDR and TPM fold-change thresholds) is made by
`DEG_byList_vClaude.pl`. This folder is the canonical genome->DEG pipeline (it replaced
the retired `map_to_genome/`); its general helpers are symlinks to the same-named
top-level `rnaseq_tools/` tools, so there is a single real copy of each. Hand-written — `gen_script_index.sh` leaves it alone (no AUTOGEN
marker)._
