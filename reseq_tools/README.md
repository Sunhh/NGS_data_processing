# reseq_tools — resequencing / population-genomics toolkit

From short reads to variants to population-genetics results. Each **analysis pipeline** lives
in its own subdirectory with a workflow README; the loose top-level scripts are a **SNP-table
toolbox** (stats, filtering, counting, conversion) shared by those pipelines.

Genotypes flow as either a **`cols` SNP table** (`chr, pos, <per-sample allele…>`) or a
**`vcf.tab`** (`bcftools query`/VariantsToTable style); `cnvt_tools/` moves between formats.

---

## 1. Analysis pipelines (see each subdir's README)
| goal | pipeline | where |
|------|----------|-------|
| **Variant calling** (reads → SNP/InDel) | GATK3 GVCF + samtools/mao (config-driven) | [`gatk/`](gatk/) |
| **Format conversion** | `cols`↔`vcf.tab`↔ped/vcf/fasta/phylip/fstat/pca/… | [`cnvt_tools/`](cnvt_tools/) |
| **Population differentiation (FST)** | hierfstat, windowed | [`fst/`](fst/) |
| **Selective sweeps** | π-ratio / ROD windows | [`slct_sweep/`](slct_sweep/) |
| **XP-CLR selection scan** | XP-CLR → shared sweep caller | [`xpclr/`](xpclr/) |
| **GWAS** | FaST-LMM (+ PCA covariates, peaks, LD, plot) | [`gwas_tools/`](gwas_tools/), [`pca/`](pca/), [`LD_ana/`](LD_ana/), [`tassel/`](tassel/) |
| **BSA** (bulk-segregant) | QTLseqr + custom G′ | [`bsa/`](bsa/) |
| **Introgression / admixture** | group major-allele classification | [`detect_mix/`](detect_mix/) |
| **SNP phylogeny** | random-sampled datasets → NJ (MEGA) | [`phylo_tools/`](phylo_tools/) |
| **Variant-effect post-processing** | simplify SnpEff SV/InDel annotation | [`snpeff_data/`](snpeff_data/) |

Compiled helpers used by `gatk/`: [`mao_exe/`](mao_exe/), [`C_exe/`](C_exe/) (with `.c` sources).
`filter_tools/cnt_depth/`, `for_yb63/` are small/project-specific; `example_data/` holds demo files.

---

## 2. SV / read-based workflows (`scripts/`)
Count read pairs spanning a target region — genotypes large deletions/insertions (a whole-genome
reference is recommended so repeat-related pairs are not forced onto the junction):
```sh
perl scripts/cntRd_spanJunctionSite_inBam.pl -in_bam PE2allele.bam \
  -target_loc seq:S-E -positive_control_loc seq:nearbyS-nearbyE > PE2allele.read_counts
```
Genotype a list of long deletions across many BAMs → melt + matrix:
```sh
perl scripts/genotype_longDEL_byListBam.pl list.long_deletions list.sample_bam out_geno
perl scripts/cnvt_melt_to_matrix.pl out_geno-melt.tab > out_geno-mat.tab
```
Others: `scripts/cnvt_msaFa_to_sam.pl` (ClustalW MSA → SAM), `scripts/cntRdMismatch_inSam.pl`
(per-read mismatch counts). Example inputs in `example_data/`.

---

## 3. SNP-table toolbox (top-level scripts)
Every script prints its own usage with `-h` / no arguments.

**Per-site / per-sample stats**
- `basic_snp_infor_bySite.pl` (`cols` input) / `basic_snp_infor_bySite_inVcfTab.pl` (`vcf.tab`, parallel, more fields) — per-site missing/allele/homo/hete/MAF
- `snpTbl_stats.pl` — windowed π / θ / Tajima's D (uses `PopGenSunhh`)
- `cnt_maf_ratio.pl`, `cnt_homo_hete_ratio.pl`, `cnt_Nmiss_ratio_heteNotN.pl`, `cnt_NHH_byIndv.pl`, `cnt_diff_inTbl.pl`, `cnt_genotype_in_1col.pl`

**Pileup** — `extract_pileup.pl`, `cnt_pileup_depC.pl`, `cnt_genotype_inpileup.pl`

**BAM read counts / depth** — `cnt_mapRdN_inBam.pl`, `rdNum_in_bam.pl`

**Variant effect** — `SNP_effect.pl` (classify effects from a CDS GFF + chr fasta), `class_SNPeffect_tbl.pl`

**Filter / mask** — `mask_vcf_geno_byGQ.pl`, `mask_weiredSNP.pl`, `maskClose_in_1col.pl`,
`rm_adjacent_sites.pl`, `rm_Nmiss_sites.pl`, `rm_same_site.pl` (`-hete2N`), `slim_SNP_sites.pl`,
`get_set2_varOnlyHete.pl`, `get_set3_varWiIndel.pl`

**Region / window / site selection** — `SNP_in_region.pl`, `extract_sites_by_list.pl`,
`mk_wind_from_noNlis.pl`, `rand_site_wiWind.pl`

**Convert / edit** — `vcf_simplify_addRef.pl`, `add_SNP1col_to_basic.pl`, `cnvt_agp2chain.pl`
(AGP → UCSC chain), `rename_plink_map.pl`, `lcnt_to_represent_allele.pl`

**Plot** — `draw_SNP_dist.pl`, `extract_top1_vcfPI_wind.R`

_Hand-written — `gen_script_index.sh` leaves it alone. Per-directory `README.md` files give the
per-analysis walk-throughs; `../SCRIPTS_INDEX.md` is a flat index of every script in the repo._
