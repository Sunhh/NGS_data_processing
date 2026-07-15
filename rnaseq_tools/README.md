# rnaseq_tools — RNA-seq expression & DEG analysis

Tools for quantifying gene expression from RNA-seq reads, calling differentially
expressed genes (DEGs), and downstream analyses (co-expression, heatmaps), plus a
specialised graft read-origin toolkit.

The collection grew across several projects, so it spans **two eras**: a modern,
R-driven pipeline (recommended) and an older, self-contained **Perl** count→DEG
toolkit. This page maps both and says which script to reach for.

---

## 1. Pick a route

| Goal | Pipeline | Where |
|------|----------|-------|
| Expression + DEG from a **reference genome** *(recommended)* | HISAT2 (2-pass) → featureCounts → TPM → DESeq2/edgeR FDR → DEG calls | [`pipe_genome_DEG/README.md`](pipe_genome_DEG/README.md) |
| Expression + DEG from a **transcriptome fasta** | Salmon **or** kallisto → tximport → DESeq2 | [`map_to_transcriptome/README.md`](map_to_transcriptome/README.md) |
| **Co-expression** modules & trait association | WGCNA (`run_wgcna.r` → `down_phenoAssoc.r`) | [`coexp/`](coexp/) |
| **Graft** read-origin detection (which parent genome a read maps to) | dual-map + read sorting | [`graft/`](graft/) and §4 |
| A **Perl-only** count→TPM→DEG toolkit (no featureCounts/DESeq2) | top-level scripts | §3 below |

For most new work use **`pipe_genome_DEG/`** — a self-contained, tested,
copy-pasteable walk-through (build DB → 2-pass HISAT2 → featureCounts → TPM/RPKM →
DESeq2 or edgeR FDRs → DEG calls with adjustable thresholds → orthogroup roll-up),
with a runnable `example/`. The rest of this page documents the shared utilities and
the legacy Perl toolkit.

---

## 2. The recommended genome route, in one glance

```
reads.fq ─(rmRRNA_in_fqFiles.pl, optional QC)─▶ clean.fq
   │  runHisat2_with2pass.pl  (+ fix_NHnum.pl to normalise NH tags)
   ▼
*_fixNH.bam ──featureCounts (-M -s 2)──▶ joint-cnt.txt
   ├─ cnvt_cnt_to_normExpr.pl ─▶ joint-tpm.txt   (TPM and/or RPKM)
   └─ DEG_byList.pl (DESeq2 / edgeR, FDR only) ─▶ joint-fdr.txt
              │  DEG_byList_vClaude.pl  (merge TPM + FDR; adjustable --fdr_cut / --fc)
              ▼
        DEG_calls.tsv ──cnvt_gene2group_DEGlabel.pl──▶ per-gene / -orthogroup up/down/not
```
Full commands + a runnable `example/`: [`pipe_genome_DEG/README.md`](pipe_genome_DEG/README.md).


---

## 3. Perl toolkit (top-level `*.pl`) — an alternative count→DEG chain

Use these when you want a featureCounts/DESeq2-free path, or the individual helpers.
General-purpose scripts only are listed here; project-specific ones are in §5.
Every script prints its own usage with `-h` or no arguments.

### 3a. Preprocess / QC
| script | does |
|--------|------|
| `rmRRNA_in_fqFiles.pl` | drop rRNA reads with SortMeRNA (`-inFq R1,R2 -outPref p -db_sortmerna fa,idx`) |
| `fix_SRAfqID.pl` | rewrite SRA `.fq` read IDs to a plain `/1 /2` form |
| `summary_ht2_log.pl` | tabulate HISAT2 `*.ht2.log` mapping stats |

### 3b. Count reads per gene (BAM → per-sample counts)
| script | does |
|--------|------|
| `cntRdInGene_wiHTSeq.pl` | wrapper over `htseq-count` (`-inBam -outCnt -gff`); default is strand-reverse, union, CDS/Parent |
| `cnt_rd_in_bam.pl` | count reads over a BED of gene models (no external tool) |
| `cnt_uniqMap_in_bam.pl` | per-read unique-vs-multi hit summary from a BAM |
| `compare_SensAnti.pl` | compare a sample's sense vs antisense counts |

### 3c. Combine per-sample counts → a matrix
| script | does |
|--------|------|
| `combine_samCnt.pl` | combine the count files listed in `file_list_to_combine` into one matrix (GeneID + 1 column/sample) |
| `join_samCnt.pl` | `paste`-style join of several count files (`-toJoinList`) |

### 3d. Normalise / express as TPM
| script | does |
|--------|------|
| `cnvt_cnt_to_normExpr.pl` | counts → **TPM and/or RPKM** (`--cntFn --lenFn --out_type tpm\|rpkm\|both`); reads featureCounts output or a plain matrix+length; `--libFn` for a true library size |
| `get_meanTPM.pl` | collapse per-replicate TPM to per-group mean |
| `get_DESeqNormCnt.pl` | apply DESeq2 size factors (embedded in the count table) to get normalised counts |

### 3e. Call DEGs
| script | engine / input | note |
|--------|----------------|------|
| `DEG_byList.pl` | **FDR producer:** `-rdCntFn -compareList -outFDRFn`; DESeq2 or edgeR (classic/GLM) per comparison → one `ds.`/`et.`/`eg.` FDR table (`-writeRcode` dumps the R) | pairs with `_vClaude` |
| `DEG_byList_vClaude.pl` | **DEG caller:** `--tpm --pair --fdr` → up/down/not per comparison, adjustable `--fdr_cut`/`--fc` | merges TPM + FDR |

### 3f. Downstream
| script | does |
|--------|------|
| `plot_heatmap_by_geneList.pl` | expression heatmap for a gene list (`<exp.txt> <gene_list> <out_pref>`; uses `plot_expr_heatmap.r`) |
| `fix_excelV.pl` | rewrite tiny/scientific values so Excel doesn't mangle them |
| `combine_DEGs.pl` | merge per-comparison DESeq2 `res-*` files (`meta_file [cols]`) into one wide table |
| `label_DEGs.pl` | label that merged table `N`/`U`/`D` per comparison (FDR<0.05 & |log2FC|>1); used by `map_to_transcriptome/` |

### Typical Perl-toolkit chain
```
BAM ─cntRdInGene_wiHTSeq.pl─▶ per-sample .cnt
    ─combine_samCnt.pl (list)─▶ count matrix
    ├─cnvt_cnt_to_normExpr.pl ─▶ TPM/RPKM ─get_meanTPM.pl─▶ group-mean TPM
    └─DEG_byList.pl (DESeq2/edgeR) ─▶ FDR table
                                    └─DEG_byList_vClaude.pl (+ TPM) ─▶ DEG calls
    ─▶ plot_heatmap_by_geneList.pl
```

---

## 4. Read-origin / graft utilities

`sep_reads_by_toRef.pl` and `extract_samAln_by_fq.pl` (top level) split reads by
which reference they align to. The full **graft** work-flow (map source reads to
both parent genomes, then decide each read's origin and transmit labels) lives in
[`graft/`](graft/): `simple_pipe_find_graft_rd_SE.pl` is the SE entry point;
`find_transmit_step1/2/3.pl` (and `stepX1/X2`) are the staged PE pipeline;
`get_alnBam_by_src2tgt_rdList.pl` extracts the chosen alignments.

---

## 5. Project-specific scripts (watermelon graft / pan-genome) — not general

These hardcode that project's sample-naming (`[SM](FL|FR|LV|RT|SD|ST)(F1|P1|P3)_rep\d+`,
`P[13]g_*`, `_toP1/_toP3` count suffixes), so they only work on that dataset:
`add_sizefactor.pl`, `get_MPCnt.pl`, `FET_getCnt.pl`, `DEGtool_withSizeFactor.pl`,
`combine_samCnt_by_Pref.pl`, `combine_samCnt_by_Pref_antiSense.pl` (the latter two are
the sense / antisense variants of the same per-prefix combiner). Reuse them only after
editing the naming regexes, or start from the general scripts in §3.

`draw_SNP.pl` (top level) is a per-chromosome SNP karyotype drawer, unrelated to the
expression pipeline.

---

## 6. Variants & duplicates (which to prefer)

- **DEG calling:** `DEG_byList.pl` (DESeq2/edgeR → FDR table) feeds `DEG_byList_vClaude.pl`
  (merges TPM + FDR → up/down/not, tunable thresholds) — they are a pair, not rivals.
  `DEGtool_withSizeFactor.pl` is project-specific.
- **Count combiner:** `combine_samCnt.pl` (general, list-driven) vs the project-specific
  `combine_samCnt_by_Pref*.pl`.
- `runHisat2_with2pass.pl` / `fix_NHnum.pl` now live under `pipe_genome_DEG/`; both the
  genome pipeline and the `graft/` scripts call them from there (the old `map_to_genome/`
  pipeline they belonged to has been retired, superseded by `pipe_genome_DEG/`).

_Maintenance note: this README is hand-written; `gen_script_index.sh` will NOT overwrite it
(only files carrying the AUTOGEN marker). For a flat per-script index of the whole repo see
`../SCRIPTS_INDEX.md`._
