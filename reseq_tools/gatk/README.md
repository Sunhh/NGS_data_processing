# gatk — resequencing variant calling (GATK GVCF + samtools/mao)

A **config-driven** SNP/InDel calling system: every pipeline script reads a `-conf_file`
(tool paths, reference, dbSNP, parameters), so you edit one config and run. Two calling
paths share that config:
- **GATK path** (per-sample GVCF → joint genotyping),
- **samtools/mao path** (mpileup-based, using the small C helpers in `../mao_exe/` + `../C_exe/`).

> Note: this pipeline targets **GATK3** (`GenomeAnalysisTK.jar`, `-T HaplotypeCaller/GenotypeGVCFs`).
> For new projects a GATK4 / Sentieon workflow (as in [`../bsa/README.md`](../bsa/README.md)) is
> recommended; this directory is kept for the established watermelon-resequencing pipeline.

## Config
Copy and edit `pipe_gatk_conf` (or `conf_pipe_gatk`) — a `key <TAB> value` file read by every
pipe script. Key fields: `exe_java exe_bwa exe_samtools exe_bcftools`, `jar_picard jar_gatk`,
`ref_fasta`, `known_dbsnp dirDbsnp`, `para_bwa para_bqsr nct_hapCaller`, `dir_tmp`, and the mao
binaries `dirMao exe_combine2PileFiles exe_reSeqPrintRefChr exe_reSeqPrintSample`.

## Scripts
| script | role |
|--------|------|
| `pipe_gatk_inFqList.pl` | **main GATK entry**: fastq list → per-sample GVCF (uBAM → mark adapters → BWA-MEM → merge → dedup → BQSR → HaplotypeCaller `-ERC GVCF`) |
| `pipe_gatk.pl` | earlier single-config GATK pipeline (fastq → GVCF) |
| `pipe_gatk_singleAlnBam.pl` | GATK pipeline starting from one already-aligned BAM |
| `gatk_dvd_step8_combineGVCF_interval.pl` | large cohort: `CombineGVCFs` split by interval |
| `gatk_dvd_step9_gvcf2var.pl` | `GenotypeGVCFs`: combined GVCF → variant VCF |
| `gatk_dvd_step9_gvcf2var_fromGVCF.pl` | step-9 variant for an already single/VAR gVCF input |
| `CatVariants.pl` | concatenate the per-interval variant VCFs (GATK `CatVariants`) |
| `get_pass_vcf.pl` | keep only `PASS` records (`... \| perl get_pass_vcf.pl > out_PASS.vcf`) |
| `pipe_maoSNP.pl` | samtools/mao SNP calling from `*_dedup_pipe2.bam` (mpileup + mao C helpers) |
| `pipe_sbSNP.pl` | samtools + bcftools SNP calling (htslib workflow) |
| `est_depth_inBam.pl` | per-region depth stats from a BAM (`-inBam -inBed`) |
| `revert_alnBam_to_uBam.pl` | aligned BAM → unaligned BAM/fastq (for re-processing) |

Config example files: `conf_pipe_gatk`, `pipe_gatk_conf`; `pref_lmyPM` is a sample prefix list.

---

## GATK path
```sh
# 1. per-sample GVCFs from a fastq list (SAMPLE, READ_GROUP, LIBRARY, dataPrefix, fq1, fq2)
perl pipe_gatk_inFqList.pl -conf_file pipe_gatk_conf -in_pref_list list.in_rd \
  -out_dir out/ -wrk_dir wrk/
#    (pipe_gatk_singleAlnBam.pl to start from an aligned BAM instead)

# 2. joint genotyping (large cohorts: split by interval)
perl gatk_dvd_step8_combineGVCF_interval.pl  -conf_file pipe_gatk_conf ...   # CombineGVCFs
perl gatk_dvd_step9_gvcf2var.pl              -conf_file pipe_gatk_conf ...    # GenotypeGVCFs
perl CatVariants.pl sorted_GVCF_list                                         # merge intervals

# 3. keep PASS variants (after VariantFiltration)
gzip -cd cohort_filtV.vcf.gz | perl get_pass_vcf.pl > cohort_filtV_PASS.vcf
```

## samtools / mao path (alternative caller)
```sh
perl pipe_maoSNP.pl -conf_file pipe_gatk_conf -in_pref_list step5b_out_pref_list > cmd_list.sMao
# or the samtools+bcftools variant:
perl pipe_sbSNP.pl  -conf_file pipe_gatk_conf -in_pref_list step5b_out_pref_list > cmd_list.sSb
```
These emit command lists that run mpileup + the compiled helpers (`../mao_exe/combine2PileFiles`,
`reSeqPrintRefChr/Sample`; `../C_exe/` has `rmSameSite`, `maskClose_in_1col` with their C sources).

_Hand-written — `gen_script_index.sh` leaves it alone._
