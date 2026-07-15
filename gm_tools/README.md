# `gm_tools` — genetic-map (ABH) genotyping from population SNPs

Turn a population SNP table (parents + offspring, `vcf.tab` format) into **ABH
genotype calls** (`a` = parent-1 homozygous, `b` = parent-2 homozygous,
`h` = heterozygous, `-`/`u` = missing) for R/qtl-style linkage / QTL mapping,
then collapse per-site calls into per-individual genotype **blocks**.

## Main engine — `rnaseq_snp_tools.pl`
A multi-task tool that covers the whole ABH workflow (`-task ...`):

| task | does |
|------|------|
| `cnt_pop_allele` | per-site population allele/genotype counts from a `vcf.tab` (`-P1_colN`/`-P2_colN` parents, `-offs_colN`; `-infer_PG` infers a missing parent genotype). `-out_itayCnt <file>` also writes the compact `chr pos Total A_N B_N AF%` table |
| `get_abh` | SNP `vcf.tab` -> one ABH table (chr, pos, REF, P1, P2, then `a/b/h/-` per offspring). `-out_itayCsv <pref>` also writes one R/qtl-style CSV per chromosome (`<pref><chr>.csv`) |
| `blk_abh` | per-site ABH -> per-individual blocks, with window-genotype correction, double-crossover smoothing, and short-block joining (`-len_jnBlk`, `-window_geno`, `-smooth_geno`, …) |
| `rejoin_abhblk` | re-join adjacent blocks (`-abhblk_tab`, `-len_jnBlk`, `-min_isiteN2`) |

```sh
# SNP table -> ABH table (STDOUT) + per-chromosome R/qtl CSVs
perl rnaseq_snp_tools.pl -task get_abh -vcf_tab pop.vcf.tab -P1_colN 3 -P2_colN 4 \
     -out_itayCsv pop_ > pop.abh
# per-site allele-frequency summary
perl rnaseq_snp_tools.pl -task cnt_pop_allele -vcf_tab pop.vcf.tab -P1_colN 3 -P2_colN 4 \
     -out_itayCnt pop.afcnt > pop.allele_counts
# per-site ABH -> per-individual blocks
perl rnaseq_snp_tools.pl -task blk_abh -abh_tab pop.abh -window_geno -smooth_geno > pop.abh.blk
```

## ABH format helpers
| script | role |
|--------|------|
| `abh_to_qtlCsv.pl` | ABH table (`get_abh` output) -> R/qtl `.csv` (transposed; `-gm_id_list` to pick/order markers) |
| `cnvt_abhJn_to_abh.pl` | joined block table (`.abh.jn`: Offs_ID/Chr/Start/End/abh) -> per-site `.abh` |

## Genome-version lift
| script | role |
|--------|------|
| `tomato_loc_V2p4_to_V2p5.pl` | lift locus coordinates from tomato genome V2.4 AGP to V2.5 AGP |

---

_Hand-written. The former standalone converters `cnvt_snp_to_itayNeed.pl`,
`cvt_snp_to_itayNeed.pl` and `cnvt_snp_to_itayNeed_cnt.pl` were retired: their
outputs are now produced (byte-identical, smoke-tested) by
`rnaseq_snp_tools.pl -task get_abh -out_itayCsv` and `-task cnt_pop_allele
-out_itayCnt`. The verbatim duplicate `cnt_allele_withPop.pl` was dropped here
(the copy under `evolution_tools/vcf_tab/` is kept as the vcf.tab-utility home)._
