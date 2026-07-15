# Transcriptome RNA-seq → DEG pipeline (Salmon *or* kallisto)

Quantify expression by mapping clean reads to a **transcriptome / cDNA fasta** with a
pseudo-aligner, then call DEGs. You pick the quantifier — **Salmon or kallisto** — and the
rest of the chain is identical: `tximport` merges the per-sample quant to gene level, a
tximport-aware DESeq2 produces FDRs, and the up/down/not decision is made by
`DEG_byList_vClaude.pl` with adjustable thresholds — the **same DEG step as the genome
pipeline** (`../pipe_genome_DEG/`), driven by the same `list.comparison` file.

Because the transcript→gene map here collapses transcripts onto **orthogroups** (`OGID`),
the "gene" level already *is* the orthogroup — no separate roll-up step is needed.

External programs: `salmon` **or** `kallisto`; `R` with `tximport`, `DESeq2` (and `rhdf5`
for kallisto's `abundance.h5`). Reference: the tximport vignette
(`https://bioconductor.org/packages/release/bioc/vignettes/tximport/inst/doc/tximport.html`).

Scripts in this folder:

| script | role | step |
|--------|------|------|
| `cnvt_synOGgrp_to_trans2gene.pl` | orthogroup `.grp` → `transcript_id`→`gene_id(OGID)` map | 1 |
| `merge_tx_quant_to_gene.r` | tximport merge of per-sample quant → gene-level counts/TPM/lengths + `.RDS`; `--type salmon\|kallisto` | 3 |
| `run_deseq2_tximport.r` | tximport DESeq2, **FDR only**, one `ds.<g1>_VS_<g2>` column per comparison | 4 |
| `../DEG_byList_vClaude.pl` | **DEG caller** (shared): TPM + FDR → up/down/not, adjustable `--fdr_cut`/`--fc` | 5 |

Data layout: `db/` (sequences + index), `data/` (sample & comparison lists),
`salmon_out/` or `kallisto_out/` (per-sample quant), `result/` (merged tables + `.RDS`),
`DEGs/` (FDR + DEG tables).

---

## (1) Group & project cDNA sequences → transcript↔gene(OGID) map
Pool the cDNA of the genomes you compare (e.g. orthogroup members) into one fasta and
build the transcript→gene table.
```sh
cat db/97103.c.fa db/PI296341.c.fa > db/slct-synOG.c.fa
deal_table.pl data/update_t2t/t2t-synOG.grp -column 0,6,46 | awk -F"\t" '!($2=="" && $3=="")' > data/slct-synOG.grp
perl cnvt_synOGgrp_to_trans2gene.pl data/slct-synOG.grp > data/slct-synOG.trans2gene.tsv
```

## (2) Quantify each sample — pick ONE quantifier
The per-replicate sample name (the quant output dir) should be `<group>_Rep<N>` so the
downstream group means line up (e.g. `S1G1T1P10D_Rep1`). Read list `data/list.in_rd-fruit`
uses the `runHisat2_with2pass.pl` format (col4 = dataPrefix = sample name, col5 = reads).

### 2a. Salmon
```sh
mamba activate salmon
salmon index -t db/slct-synOG.c.fa -i db/synOG_index.salmon
mkdir -p salmon_out/
deal_table.pl -column 3,4 data/list.in_rd-fruit | tail -n +2 | perl -e 'while (<>) {chomp; my @a=split; print "salmon quant -i db/synOG_index.salmon -l A -p 8 -r $a[1] -o salmon_out/$a[0]\n";}' > cx1quant
run_cmd_in_batch.pl -cpuN 8 -wait_sec 1 cx1quant
mamba deactivate
```

### 2b. kallisto
```sh
mamba activate kallisto
kallisto index -i db/synOG_index.kallisto db/slct-synOG.c.fa
mkdir -p kallisto_out/
# single-end needs -l/-s fragment length mean/sd; paired-end: drop --single and give both fastqs
deal_table.pl -column 3,4 data/list.in_rd-fruit | tail -n +2 | perl -e 'while (<>) {chomp; my @a=split; print "kallisto quant -i db/synOG_index.kallisto -t 8 --single -l 200 -s 20 -o kallisto_out/$a[0] $a[1]\n";}' > cx1quant
run_cmd_in_batch.pl -cpuN 8 -wait_sec 1 cx1quant
mamba deactivate
```

## (3) Merge to gene level with tximport
List the per-sample quant dirs, then merge. Use the matching `--type`; everything after
this step is quantifier-independent.
```sh
mkdir -p result/
# Salmon:
ls -d salmon_out/*/   | sed 's#/$##' > result/list.sample_dir
Rscript merge_tx_quant_to_gene.r -s result/list.sample_dir -g data/slct-synOG.trans2gene.tsv -t salmon   -o result/comb
# kallisto (reads abundance.h5 by default; add --quant_file abundance.tsv to skip rhdf5):
# ls -d kallisto_out/*/ | sed 's#/$##' > result/list.sample_dir
# Rscript merge_tx_quant_to_gene.r -s result/list.sample_dir -g data/slct-synOG.trans2gene.tsv -t kallisto -o result/comb
```
Outputs: `result/comb.gene_counts.tsv`, `comb.gene_tpm.tsv`, `comb.gene_lengths.tsv`,
`comb.RDS` (the tximport object).

## (4) FDRs with tximport-aware DESeq2 (significance only)
`list.sample_meta-fruit` has columns `sample` and `group`. Reuse `list.comparison` (col1
method, ignored here; col2 group1; col3 group2; col4 name) — the same file feeds step 5.
```sh
mkdir -p DEGs/
Rscript run_deseq2_tximport.r -r result/comb.RDS -m data/list.sample_meta-fruit -c data/list.comparison -o DEGs/joint-fdr.txt
```
Output: `GeneID` + one `ds.<group1>_VS_<group2>` padj column per comparison (uses the
transcript-length offsets from Salmon/kallisto via `DESeqDataSetFromTximport`).

## (5) Call DEGs in vClaude (adjustable thresholds)
The shared caller merges the gene-level TPM with the FDR table (`--pair` = the same
`list.comparison`; its method column is ignored here).
```sh
perl ../DEG_byList_vClaude.pl --tpm result/comb.gene_tpm.tsv --pair data/list.comparison --fdr DEGs/joint-fdr.txt --fdr_cut 0.05 --fc 2 --out DEGs/DEG_calls.tsv
```
Per comparison `<name>`: `<name>.v1/.v2` (mean TPM), `.R` (fold change), `.FDR`, `.DEG`
(`up` if R>fc & FDR<fdr_cut; `down` if R<1/fc & FDR<fdr_cut; else `not`). These calls are
already at the orthogroup level.

---

_See [`example/`](example/): a tiny kallisto dataset + `run_example.sh` running steps 3→5
(needs `Rscript` with tximport + DESeq2). Salmon is identical with `--type salmon` on
`salmon quant` output dirs. Hand-written — `gen_script_index.sh` leaves it alone._
